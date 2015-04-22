module matrix_io
  use mpi
  use command_argument, only : argument, matrix_info
  use event_logger_m, only : add_event
  use eigenpairs_types, only : eigenpairs_types_union
  use processes, only : check_master, terminate
  implicit none

  type sparse_mat
    integer :: size, num_non_zeros
    double precision, allocatable :: value(:)
    integer, allocatable :: suffix(:, :)
  end type sparse_mat

  private
  public :: read_matrix_file, print_matrix, sparse_mat, print_eigenvectors

contains

  subroutine read_matrix_file(filename, info, matrix, ierr)
    character(*), intent(in) :: filename
    type(matrix_info), intent(in) :: info
    type(sparse_mat), intent(out) :: matrix
    integer, intent(out) :: ierr

    double precision :: time_start, time_start_part, time_end
    integer, parameter :: iunit = 8
    integer :: rows, cols, num_non_zeros

    time_start = mpi_wtime()
    time_start_part = time_start

    if (check_master()) then
      print '("start reading matrix file ", a)', trim(filename)
    end if

    matrix%size = info%rows
    matrix%num_non_zeros = info%entries

    allocate(matrix%suffix(2, info%entries), matrix%value(info%entries), stat = ierr)
    if (ierr /= 0) then
      return
    end if

    time_end = mpi_wtime()
    call add_event('read_matrix_file:allocate', time_end - time_start_part)
    time_start_part = time_end

    open(iunit, file=filename, status='old', iostat=ierr)
    if (ierr /= 0) then
      return
    end if
    ! read_matrix_file_header is added to skip comment lines
    call read_matrix_file_header(iunit, rows, cols, num_non_zeros)

    time_end = mpi_wtime()
    call add_event('read_matrix_file:header', time_end - time_start_part)
    time_start_part = time_end

    call read_matrix_file_value(0, iunit, info%rows, &
         info%entries, matrix%value, matrix%suffix)
    close(iunit)

    time_end = mpi_wtime()
    call add_event('read_matrix_file:value', time_end - time_start_part)
    call add_event('read_matrix_file', time_end - time_start)
  end subroutine read_matrix_file


  subroutine read_matrix_file_header(unit_num, rows, cols, num_non_zeros)
    integer, intent(in) :: unit_num

    integer :: ierr, rows, cols, num_non_zeros
    character(len=1024) :: chara_wrk

    ! Read the first line
    read(unit_num,'(a)') chara_wrk

    ! Plot the comment lines
    do while (.true.)
      read(unit_num,'(a)') chara_wrk
      if (index(chara_wrk, '%') /= 1) exit
    enddo

    read (chara_wrk, *, iostat = ierr) rows, cols, num_non_zeros
  end subroutine read_matrix_file_header


  subroutine read_matrix_file_value(verbose_level, unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
    integer, intent(in) :: verbose_level
    integer, intent(in) :: unit_num
    integer, intent(in) :: mat_size
    integer, intent(in) :: num_non_zeros
    double precision, allocatable, intent(inout) :: mat_value(:)
    integer, allocatable, intent(inout) :: mat_suffix(:,:)

    integer :: line_count, print_count
    logical :: debug_mode
    integer :: i, j, ierr
    double precision :: value_wrk

    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif

    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_value'
    if (debug_mode) write(*,'(a)')' ONLY REAL-SYMMETRIC MATRIX is supported'

    if (debug_mode) write(*,*)'size(mat_value)=',size(mat_value)

    if (debug_mode) write(*,*)'size(mat_suffix,1)=',size(mat_suffix,1)
    if (debug_mode) write(*,*)'size(mat_suffix,2)=',size(mat_suffix,2)

    print_count = 1
    do line_count = 1, num_non_zeros
      if (line_count > num_non_zeros / 10 * print_count .and. check_master()) then
         write (0, *) '[Event', mpi_wtime(), '] read matrix element number ', line_count
         print_count = print_count + 1
      end if

      read(unit_num, *, iostat=ierr) i, j, value_wrk
      if (ierr /= 0) then
        call terminate('read_matrix_file_value: invalid format of matrix value', ierr)
      endif

      ! if (debug_mode) write(*,*)'i,j, matrix_value=',i, j, value_wrk

      if (i < 1 .or. i > mat_size .or. j < 1 .or. j > mat_size) then
        call terminate('read_matrix_file_value: index of matrix out of range', ierr)
      endif

      mat_value(line_count) = value_wrk
      mat_suffix(1,line_count) = i
      mat_suffix(2,line_count) = j
    enddo

    ! if (num_non_zeros /= size(mat_value,1)) then
    ! endif
  end subroutine read_matrix_file_value


  subroutine print_matrix(name, mat, m, n)
    character(*), intent(in) :: name
    double precision, intent(in) :: mat(:, :)

    integer :: i, j, m, n
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    if (m < 0) then
       m = size(mat, 1)
    end if
    if (n < 0) then
       n = size(mat, 2)
    end if
    do j = 1, n
       do i = 1, m
          print ' (A, "(", I6, ",", I6, ")=", D30.18) ', name, i, j, mat(i, j)
       end do
    end do

    time_end = mpi_wtime()
    call add_event('print_matrix', time_end - time_start)
  end subroutine print_matrix


  subroutine print_eigenvectors(arg, eigenpairs)
    type(argument) :: arg
    type(eigenpairs_types_union) :: eigenpairs

    double precision :: work(arg%matrix_A_info%rows), time_start, time_end
    integer :: max_num_digits, len, i, digit, stat, err
    character(512) :: num_str, filename
    integer, parameter :: iunit = 10

    time_start = mpi_wtime()

    if (eigenpairs%type_number == 1) then
      stop '[Error] print_eigenvectors: printer for a local matrix not implemented yet'
    else if (eigenpairs%type_number == 2) then
      max_num_digits = int(log10(real(arg%printed_vecs_end))) + 1
      do i = arg%printed_vecs_start, arg%printed_vecs_end
        if (check_master()) then
          write (num_str, '(i0)') i
          len = len_trim(num_str)
          num_str(max_num_digits - len + 1 : max_num_digits) = num_str(1 : len)
          do digit = 1, max_num_digits - len
            num_str(digit:digit) = '0'
          end do

          filename = trim(arg%eigenvector_dir) // '/' // trim(num_str) // '.dat'
          open (iunit, file=trim(filename), status='replace', iostat=stat)
        end if

        call mpi_bcast(stat, 1, mpi_integer, 0, mpi_comm_world, err) ! Share stat for file opening
        if (stat /= 0) then
          if (check_master ()) print *, 'iostat: ', stat
          call terminate('print_eigenvectors: cannot open ' // trim(filename), stat)
        end if

        call mpi_barrier(mpi_comm_world, stat)
        call eigentest_pdlaprnt(arg%matrix_A_info%rows, 1, eigenpairs%blacs%Vectors, &
             1, i, eigenpairs%blacs%desc, 0, 0, '', iunit, work)

        if (check_master()) then
          close (iunit)
        end if
      end do
    end if

    time_end = mpi_wtime()
    call add_event('print_eigenvectors', time_end - time_start)
  end subroutine print_eigenvectors
end module matrix_io
