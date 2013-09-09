module matrix_io
  use command_argument, only : argument, matrix_info
  use eigenpairs_types, only : eigenpairs_types_union
  implicit none

  type sparse_mat
    integer :: size, num_non_zeros
    double precision, allocatable :: value(:)
    integer, allocatable :: suffix(:, :)
  end type sparse_mat

  private
  public :: read_matrix_file, print_matrix, sparse_mat, print_eigenvectors

contains

  subroutine read_matrix_file(filename, info, matrix)
    character(*), intent(in) :: filename
    type(matrix_info), intent(in) :: info
    type(sparse_mat), intent(out) :: matrix

    integer, parameter :: iunit = 8

    matrix%size = info%rows
    matrix%num_non_zeros = info%entries

    allocate(matrix%suffix(2, info%entries), matrix%value(info%entries))

    open(iunit, file = filename)

    ! read_matrix_file_header is added to skip comment lines
    call read_matrix_file_header(iunit)
    call read_matrix_file_value(0, iunit, info%rows, &
         info%entries, matrix%value, matrix%suffix)

    close(iunit)
  end subroutine read_matrix_file


  subroutine read_matrix_file_header(unit_num)
    integer, intent(in) :: unit_num

    integer :: ierr, mat_size, mat_size2, num_non_zeros
    character(len=1024) :: chara_wrk

    ! Read the first line
    read(unit_num,'(a)') chara_wrk

    ! Plot the comment lines
    do while (.true.)
      read(unit_num,'(a)') chara_wrk
      if (index(chara_wrk, '%') /= 1) exit
    enddo

    read (chara_wrk, *, iostat = ierr) mat_size, mat_size2, num_non_zeros
  end subroutine read_matrix_file_header


  subroutine read_matrix_file_value(verbose_level, unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
    integer, intent(in) :: verbose_level
    integer, intent(in) :: unit_num
    integer, intent(in) :: mat_size
    integer, intent(in) :: num_non_zeros
    double precision, allocatable, intent(inout) :: mat_value(:)
    integer, allocatable, intent(inout) :: mat_suffix(:,:)

    integer :: line_count
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

    do line_count = 1, num_non_zeros
      read(unit_num, *, iostat=ierr) i, j, value_wrk
      if (ierr /= 0) then
        stop '[Error] read_matrix_file_value: Invalid format of matrix value'
      endif

      ! if (debug_mode) write(*,*)'i,j, matrix_value=',i, j, value_wrk

      if (i < 1 .or. i > mat_size) then
        stop '[Error] read_matrix_file_value: Index of matrix out of range'
      endif

      if (j < 1 .or. j > mat_size) then
        stop '[Error] read_matrix_file_value: Index of matrix out of range'
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
  end subroutine print_matrix

  subroutine print_eigenvectors(arg, eigenpairs)
    type(argument) :: arg
    type(eigenpairs_types_union) :: eigenpairs

    double precision :: work(arg%matrix_A_info%rows)
    integer :: max_num_digits, len, i, j
    character(16) :: num_str

    if (eigenpairs%type_number == 1) then
      stop '[Error] print_eigenvectors: printer for' // &
           'a local matrix not implemented yet'
    else if (eigenpairs%type_number == 2) then
      max_num_digits = int(log10(real( &
           arg%printed_vecs_end - arg%printed_vecs_start + 1))) + 1
      do i = arg%printed_vecs_start, arg%printed_vecs_end
        write (num_str, '(i0)') i
        len = len_trim(num_str)
        num_str(max_num_digits - len + 1 : max_num_digits) = num_str(1 : len)
        num_str(1 : max_num_digits - len) = '0'
        call pdlawrite(trim(arg%eigenvector_dir) // '/' // &
             trim(num_str) // '.dat', &
             arg%matrix_A_info%rows, 1, eigenpairs%blacs%Vectors, &
             1, i, eigenpairs%blacs%desc, 0, 0, work)
      end do
    end if
  end subroutine print_eigenvectors
end module matrix_io
