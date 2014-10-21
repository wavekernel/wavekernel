module solver_eigenexa
  use eigen_libs

  use distribute_matrix, only : process
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use matrix_io, only : sparse_mat
  use processes, only : check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time

  implicit none

  private
  public :: setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa, eigen_solver_eigenk

contains

  subroutine setup_distributed_matrix_for_eigenexa(dim, desc_A, matrix_A, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) ::desc_A(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: nx, ny, context, info

    if (check_master()) then
      print '( "Creating 2 distributed matrices for EigenExa with &
           &M, N, MB, NB: ", I0, ", ", I0, ", ", I0, ", ", I0 )', &
           dim, dim, 1, 1
    end if

    call eigen_init()
    call eigen_get_matdims(dim, nx, ny)
    context = eigen_get_blacs_context()

    call descinit(desc_A, dim, dim, 1, 1, 0, 0, context, nx, info)

    eigenpairs%type_number = 2
    call descinit(eigenpairs%blacs%desc, dim, dim, 1, 1, 0, 0, context, nx, info)
    if (info /= 0) then
      print '(a, i0)', 'info(descinit): ', info
      call terminate('setup_distributed_matrix_for_eigenexa: descinit failed', info)
    end if

    allocate(matrix_A(nx, ny), eigenpairs%blacs%Vectors(nx, ny), &
         eigenpairs%blacs%values(dim), stat = info)
    if (info /= 0) then
      call terminate('setup_distributed_matrix_for_eigenexa: allocation failed', info)
    end if

    matrix_A(:, :) = 0.0d0
    eigenpairs%blacs%Vectors(:, :) = 0.0d0
  end subroutine setup_distributed_matrix_for_eigenexa


  ! uplo: takes the value of 'L' or 'U' or (not present).
  !       Specifies how the entries of the input matrix is stored.
  !       If not present, it is assumed that both of upper and lower
  !       triangular parts are stored.
  subroutine eigen_solver_eigenexa(mat, desc_mat, n_vec, eigenpairs, uplo)
    include 'mpif.h'

    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    integer :: dim, nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128

    ! Time
    integer, parameter :: n_intervals = 1
    integer :: i, base_count
    double precision :: t_intervals(n_intervals)
    double precision :: t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'total'/)

    call get_wall_clock_base_count(base_count)

    dim = desc_mat(rows_)
    if (dim /= n_vec) then
      call terminate('eigen_solver_eigenexa: current version of EigenExa does not support partial eigenvector computation', 1)
    end if

    ! Unlike usual ScaLAPACK routines, EigenExa requires both of upper and lower
    ! triangular parts of the input matrix. Therefore copying from the lower
    ! triangular parts to upper one (or inverse) is needed.
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

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

    call eigen_sx(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)

    call get_wall_clock_time(base_count, t_all_end)

    t_intervals(1) = t_all_end

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


  subroutine eigen_solver_eigenk(mat, desc_mat, n_vec, eigenpairs, uplo)
    include 'mpif.h'

    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    integer :: dim, nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128

    ! Time
    integer, parameter :: n_intervals = 1
    integer :: i, base_count
    double precision :: t_intervals(n_intervals)
    double precision :: t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'total'/)

    call get_wall_clock_base_count(base_count)

    dim = desc_mat(rows_)
    if (dim /= n_vec) then
      call terminate('eigen_solver_eigenexa: current version of EigenExa does not support partial eigenvector computation', 1)
    end if

    ! Unlike usual ScaLAPACK routines, EigenExa requires both of upper and lower
    ! triangular parts of the input matrix. Therefore copying from the lower
    ! triangular parts to upper one (or inverse) is needed.
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
        call terminate("eigen_solver_eigenk: uplo must be 'U' or 'L'", 1)
      end if
    end if

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

    call eigen_s(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)

    call get_wall_clock_time(base_count, t_all_end)

    t_intervals(1) = t_all_end

    if (my_rank == 0) then
      call MPI_Reduce(MPI_IN_PLACE, t_intervals, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      print '("elapsed time (sec)")'
      do i = 1, n_intervals
        print '("  ", a, " ", f12.2)', interval_names(i), t_intervals(i)
      end do
    else
      call MPI_Reduce(t_intervals, 0, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    end if
  end subroutine eigen_solver_eigenk
end module solver_eigenexa
