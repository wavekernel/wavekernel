module m
  use processes
  implicit none

  integer, parameter :: desc_size = 9, desc_type_ = 1, context_ = 2, &
       rows_ = 3, cols_ = 4, block_row_ = 5, block_col_ = 6, &
       rsrc_ = 7, csrc_ = 8, local_rows_ = 9

contains

  integer function get_local_cols(proc, desc)
    type(process), intent(in) :: proc
    integer, intent(in) :: desc(desc_size)

    integer :: numroc

    get_local_cols = max(1, numroc(desc(cols_), desc(block_col_), &
         proc%my_proc_col, 0, proc%n_procs_col))
  end function get_local_cols


  subroutine setup_distributed_matrix(proc, rows, cols, desc, mat)
    type(process), intent(in) :: proc
    integer, intent(in) :: rows, cols
    integer, intent(out) :: desc(desc_size)
    real(8), intent(out), allocatable :: mat(:, :)

    integer :: numroc
    integer :: actual_block_size, local_rows, info

    actual_block_size = 128

    local_rows = max(1, numroc(rows, actual_block_size, &
         proc%my_proc_row, 0, proc%n_procs_row))

    call descinit(desc, rows, cols, actual_block_size, actual_block_size, &
         0, 0, proc%context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols(proc, desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('setup_distributed_matrix: allocation failed', info)
    end if

    mat(:, :) = 0.0d0
  end subroutine setup_distributed_matrix


  subroutine transpose(mat, desc_mat, uplo)
    include 'mpif.h'

    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size)
    character, intent(in), optional :: uplo

    integer :: dim, i, ierr
    double precision :: times(3)

    call mpi_barrier(mpi_comm_world, ierr)
    times(1) = mpi_wtime()

    dim = desc_mat(rows_)

    if (present(uplo)) then
      if (uplo == 'U') then ! Note: Not tested
        do i = 1, dim - 1
          call pdcopy(dim - i, mat, i, i + 1, desc_mat, dim, &
               mat, i + 1, i, desc_mat, 1)
        end do
      else if (uplo == 'L') then
        do i = 1, dim - 1
          call pdcopy(dim - i, mat, i + 1, i, desc_mat, 1, &
               mat, i, i + 1, desc_mat, dim)
        end do
      else
        call terminate("eigen_solver_eigenexa: uplo must be 'U' or 'L'", 1)
      end if
    end if
  end subroutine transpose
end module m


program main
  use m
  use processes
  use mpi
  implicit none

  real(8), allocatable :: mat(:, :)
  character(len=256) :: arg_str
  integer :: dim, ierr, desc(desc_size)
  type(process) :: proc
  real(8) :: times(3)

  if (command_argument_count() < 1) then
    stop 'Specify dimension'
  else
    call get_command_argument(1, arg_str)
    read (arg_str, *) dim
  end if

  call setup_distribution(proc)

  call mpi_barrier(mpi_comm_world, ierr)
  times(1) = mpi_wtime()

  call setup_distributed_matrix(proc, dim, dim, desc, mat)

  call mpi_barrier(mpi_comm_world, ierr)
  times(2) = mpi_wtime()

  call transpose(mat, desc, 'L')

  call mpi_barrier(mpi_comm_world, ierr)
  times(3) = mpi_wtime()

  if (check_master()) then
    print *, 'setup : ', times(2) - times(1)
    print *, 'trans : ', times(3) - times(2)
  end if

  call mpi_finalize(ierr)

end program main
