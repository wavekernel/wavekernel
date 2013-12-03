module solver_eigenexa
  use eigen_libs

  use distribute_matrix, only : process
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use matrix_io, only : sparse_mat
  use processes, only : terminate
  use time, only : get_wclock_time

  implicit none

  private
  public :: setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa

contains

  subroutine setup_distributed_matrix_for_eigenexa(dim, desc_A, matrix_A, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) ::desc_A(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: nx, ny, context, info

    call eigen_init()
    call eigen_get_matdims(dim, nx, ny)
    context = eigen_get_blacs_context()

    call descinit(desc_A, dim, dim, 1, 1, 0, 0, context, nx, info)

    eigenpairs%type_number = 2
    call descinit(eigenpairs%blacs%desc, dim, dim, 1, 1, 0, 0, context, nx, info)
    if (info /= 0) then
      print '(a, i0)', 'info(descinit): ', info
      call terminate('[Error] setup_distributed_matrix_for_eigenexa: descinit failed')
    end if

    allocate(matrix_A(nx, ny), eigenpairs%blacs%Vectors(nx, ny), &
         eigenpairs%blacs%values(dim), stat = info)
    if (info /= 0) then
      print '(a, i0)', 'stat(allocate): ', info
      call terminate('[Error] setup_distributed_matrix_for_eigenexa: allocation failed')
    end if

    matrix_A(:, :) = 0.0d0
    eigenpairs%blacs%Vectors(:, :) = 0.0d0
  end subroutine setup_distributed_matrix_for_eigenexa


  subroutine eigen_solver_eigenexa(mat, dim, n_vec, eigenpairs)
    include 'mpif.h'

    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: dim, n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs

    integer :: nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128

    ! Time
    integer, parameter :: n_intervals = 1
    integer :: i
    double precision :: t_intervals(n_intervals)
    double precision :: t_init, t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'total'/)

    call get_wclock_time(t_init)

    if (dim /= n_vec) then
      call terminate('[Error] eigen_solver_eigenexa: current version of EigenExa does not support partial eigenvector computation')
    end if

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

    call eigen_sx(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_all_end - t_init

    if (my_rank == 0) then
      call MPI_Reduce(MPI_IN_PLACE, t_intervals, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      print '("elapsed time (sec)")'
      do i = 1, n_intervals
        print '("  ", a, " ", f12.2)', interval_names(i), t_intervals(i)
      end do
    else
      call MPI_Reduce(t_intervals, 0, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    end if
  end subroutine eigen_solver_eigenexa
end module solver_eigenexa
