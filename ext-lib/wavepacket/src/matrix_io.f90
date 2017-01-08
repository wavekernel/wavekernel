module wp_matrix_io_m
  use mpi
  use wp_event_logger_m, only : add_event
  use wp_processes_m, only : check_master, terminate
  use wp_global_variables_m
  implicit none

  type sparse_mat
    integer :: size, num_non_zeros
    real(8), allocatable :: value(:)
    integer, allocatable :: suffix(:, :)
  end type sparse_mat

  private
  public :: get_dimension, read_matrix_file, set_sparse_matrix_identity, print_matrix, sparse_mat, &
       destroy_sparse_mat, copy_sparse_matrix, sparse_matrix_to_diag, print_sparse_matrix

contains

  subroutine get_dimension(filename, dim)
    character(*), intent(in) :: filename
    integer, intent(out) :: dim

    integer, parameter :: iunit = 9
    integer :: ierr, rows, cols, num_non_zeros
    character(len=1024) :: line, filename_filled

    call fill_integer_format(filename, 1, filename_filled)
    open(iunit, file=filename_filled)
    do while (.true.)
      read(iunit,'(a)') line
      if (index(line, '%') /= 1) exit
    enddo
    read (line, *, iostat=ierr) rows, cols, num_non_zeros
    close(iunit)

    if (rows /= cols) then
      call terminate('get_dimension: not square matrix', 1)
    else
      dim = rows
    end if
  end subroutine get_dimension


  subroutine read_matrix_file(filename, step, matrix)
    character(*), intent(in) :: filename
    integer, intent(in) :: step
    type(sparse_mat), intent(out) :: matrix

    double precision :: time_start, time_start_part, time_end
    integer, parameter :: iunit = 8
    integer :: rows, cols, num_non_zeros, ierr
    character(len=1024) :: filename_filled

    time_start = mpi_wtime()
    time_start_part = time_start

    call fill_integer_format(filename, step, filename_filled)

    if (check_master()) then
      if (trim(filename) /= trim(filename_filled)) then
        write (0, '(A, F16.6, 4A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
             '] read_matrix_file() integer format filled: ', trim(filename), ' -> ', trim(filename_filled)
      end if
      write (0, '(A, F16.6, 2A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] read_matrix_file() filename: ', trim(filename_filled)
    end if

    open(iunit, file=filename_filled)
    call read_matrix_file_header(iunit, rows, cols, num_non_zeros)

    if (rows /= cols) then
      call terminate('read_matrix_file: matrix must be square', 1)
    end if

    matrix%size = rows
    matrix%num_non_zeros = num_non_zeros

    allocate(matrix%suffix(2, num_non_zeros), matrix%value(num_non_zeros), stat = ierr)
    if (ierr /= 0) then
      call terminate('read_matrix_file: allocation failed', ierr)
    end if

    time_end = mpi_wtime()
    call add_event('read_matrix_file:allocate', time_end - time_start_part)
    time_start_part = time_end

    call read_matrix_file_value(iunit, rows, num_non_zeros, matrix%value, matrix%suffix)
    close(iunit)

    time_end = mpi_wtime()
    call add_event('read_matrix_file:value', time_end - time_start_part)
    call add_event('read_matrix_file:total', time_end - time_start)
  end subroutine read_matrix_file


  subroutine read_matrix_file_skip_steps(iunit, step)
    integer, intent(in) :: iunit, step
    integer :: i, rows, cols, num_non_zeros, line_num
    character(len=1024) :: line

    do i = 1, step - 1  ! step == 1 means first step, no skip.
      call read_matrix_file_header(iunit, rows, cols, num_non_zeros)
      do line_num = 1, num_non_zeros
        read(iunit, '(A)') line
      end do
    end do
  end subroutine read_matrix_file_skip_steps


  subroutine read_matrix_file_header(unit_num, rows, cols, num_non_zeros)
    integer, intent(in) :: unit_num

    integer :: ierr, rows, cols, num_non_zeros
    character(len=1024) :: chara_wrk

    do while (.true.)
      read(unit_num,'(a)') chara_wrk
      if (index(chara_wrk, '%') /= 1) exit
    enddo

    read (chara_wrk, *, iostat = ierr) rows, cols, num_non_zeros
  end subroutine read_matrix_file_header


  subroutine read_matrix_file_value(unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
    integer, intent(in) :: unit_num
    integer, intent(in) :: mat_size
    integer, intent(in) :: num_non_zeros
    double precision, allocatable, intent(inout) :: mat_value(:)
    integer, allocatable, intent(inout) :: mat_suffix(:,:)

    integer :: line_count, print_count
    integer :: i, j, ierr
    double precision :: value_wrk

    print_count = 1
    do line_count = 1, num_non_zeros
      if (num_non_zeros > 10000 .and. &
           line_count > num_non_zeros / 10 * print_count .and. &
           check_master()) then
         write (0, '(A, F16.6, A, I0)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
              '] read matrix element number ', line_count
         print_count = print_count + 1
      end if

      read(unit_num, *, iostat=ierr) i, j, value_wrk
      if (ierr /= 0) then
        call terminate('read_matrix_file_value: invalid format of matrix value', ierr)
      endif

      if (i < 1 .or. i > mat_size .or. j < 1 .or. j > mat_size) then
        call terminate('read_matrix_file_value: index of matrix out of range', 1)
      endif

      mat_value(line_count) = value_wrk
      mat_suffix(1,line_count) = i
      mat_suffix(2,line_count) = j
    enddo
  end subroutine read_matrix_file_value


  subroutine set_sparse_matrix_identity(dim, matrix)
    integer, intent(in) :: dim
    type(sparse_mat), intent(out) :: matrix

    real(8) :: time_start, time_end
    integer :: i, ierr

    time_start = mpi_wtime()
     if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', time_start - g_wp_mpi_wtime_init, &
           '] set_sparse_matrix_identity()'
    end if

    matrix%size = dim
    matrix%num_non_zeros = dim
    allocate(matrix%suffix(2, dim), matrix%value(dim), stat = ierr)
    if (ierr /= 0) then
      call terminate('set_sparse_matrix_identity: allocation failed', ierr)
    end if

    do i = 1, dim
      matrix%suffix(1, i) = i
      matrix%suffix(2, i) = i
      matrix%value(i) = 1d0
    end do

    time_end = mpi_wtime()
    call add_event('set_sparse_matrix_identity', time_end - time_start)
  end subroutine set_sparse_matrix_identity


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


  subroutine destroy_sparse_mat(matrix)
    type(sparse_mat), intent(inout) :: matrix

    matrix%size = -1
    matrix%num_non_zeros = -1
    if (allocated(matrix%value)) then
      deallocate(matrix%value)
    end if
    if (allocated(matrix%suffix)) then
      deallocate(matrix%suffix)
    end if
  end subroutine destroy_sparse_mat


  subroutine fill_integer_format(str_in, step, str_out)
    character(*), intent(in) :: str_in
    integer, intent(in) :: step
    character(*), intent(out) :: str_out

    integer :: first, last, max_num_digits, len, digit
    character(len=1024) :: num_str

    first = index(str_in, '%0')
    last = index(str_in(first + 2 :), 'd')
    if (first == 0 .or. last == 0) then
      write(str_out, '(A)') str_in
    else
      read(str_in(first + 2 : first + last), *) max_num_digits
      write(num_str, '(I0)') step
      len = len_trim(num_str)
      num_str(max_num_digits - len + 1 : max_num_digits) = num_str(1 : len)
      do digit = 1, max_num_digits - len
        num_str(digit : digit) = '0'
      end do
      write(str_out, '(A, A, A)') str_in(: first - 1), trim(num_str), trim(str_in(first + last + 2 :))
    end if
  end subroutine fill_integer_format


  subroutine copy_sparse_matrix(matrix_in, matrix_out)
    type(sparse_mat), intent(in) :: matrix_in
    type(sparse_mat), intent(out) :: matrix_out

    matrix_out%size = matrix_in%size
    matrix_out%num_non_zeros = matrix_in%num_non_zeros
    if (allocated(matrix_out%value)) then
      deallocate(matrix_out%value)
    end if
    if (allocated(matrix_out%suffix)) then
      deallocate(matrix_out%suffix)
    end if
    allocate(matrix_out%value(size(matrix_in%value, 1)))
    allocate(matrix_out%suffix(size(matrix_in%suffix, 1), size(matrix_in%suffix, 2)))
    matrix_out%value(:) = matrix_in%value(:)
    matrix_out%suffix(:, :) = matrix_in%suffix(:, :)
  end subroutine copy_sparse_matrix


  subroutine sparse_matrix_to_diag(A, diag)
    type(sparse_mat), intent(in) :: A
    real(8), intent(out) :: diag(A%size)

    integer :: i, j

    do i = 1, A%num_non_zeros
      j = A%suffix(1, i)
      if (j == A%suffix(2, i)) then
        diag(j) = A%value(i)
      end if
    end do
  end subroutine sparse_matrix_to_diag


  subroutine print_sparse_matrix(name, matrix, iunit)
    character(len=*), intent(in) :: name
    type(sparse_mat), intent(in) :: matrix
    integer, intent(in) :: iunit

    integer :: i

    write(iunit, '(A)') '%%MatrixMarket matrix coordinate real symmetric'
    write(iunit, '(A)') '%'
    write(iunit, *) matrix%size, matrix%size, matrix%num_non_zeros
    do i = 1, matrix%num_non_zeros
      write(iunit, *) matrix%suffix(1, i), matrix%suffix(2, i), matrix%value(i)
    end do
  end subroutine print_sparse_matrix
end module wp_matrix_io_m
