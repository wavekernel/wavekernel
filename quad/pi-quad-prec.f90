program main
  use mpi
  implicit none

  ! Select an appropriate kind number from the number of significant figures
  integer, parameter :: qp = selected_real_kind(30)
  real(qp) :: pi_divided, pi_sum
  character(len = 1) :: mode

  integer :: nprocs, myrank, ierr

  call get_command_argument(1, mode)

  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  if (myrank == 0) then
    print *, 'mode char: ', mode
    print *, 'quad precision kind: ', qp
  end if

  ! Compute (pi / nprocs), which becomes pi after summation
  pi_divided = 4.0_qp * atan(1.0_qp) / nprocs
  print *, nprocs, myrank, pi_divided

  select case (mode)
    case ('r')
      pi_sum = sum_with_reduce(pi_divided)
    case ('g')
      pi_sum = sum_with_gather(pi_divided, nprocs, myrank)
    case default
      stop 'unknown mode'
    end select

  if (myrank == 0) then
    print *, "sum:", pi_sum
  end if

  call mpi_finalize(ierr)

contains
  real(qp) function sum_with_reduce(x)
    real(qp), intent(in) :: x
    integer :: ierr

    call mpi_reduce(x, sum_with_reduce, 1, mpi_real16, &
         mpi_sum, 0, mpi_comm_world, ierr)

    if (ierr /= 0) then
      stop "mpi_reduce failed"
    end if
  end function sum_with_reduce

  real(qp) function sum_with_gather(x, nprocs, myrank)
    real(qp), intent(in) :: x
    integer, intent(in) :: nprocs, myrank
    real(qp), allocatable :: recv(:)
    integer :: ierr

    allocate(recv(nprocs))
    call mpi_gather(x, 1, mpi_real16, recv, 1, mpi_real16, &
         0, mpi_comm_world, ierr)
    if (ierr /= 0) then
      stop "mpi_gather failed"
    end if

    if (myrank == 0) then
      sum_with_gather = sum(recv)
    end if
  end function sum_with_gather
end program main
