module solver_eigenexa
  use eigen_libs
  use mpi
  use distribute_matrix, only : process
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use event_logger_m, only : add_event
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
    double precision :: time_start, time_end

    time_start = mpi_wtime()

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

    time_end = mpi_wtime()
    call add_event('setup_distributed_matrix_for_eigenexa', time_end - time_start)
  end subroutine setup_distributed_matrix_for_eigenexa


  ! uplo: takes the value of 'L' or 'U' or (not present).
  !       Specifies how the entries of the input matrix is stored.
  !       If not present, it is assumed that both of upper and lower
  !       triangular parts are stored.
  subroutine eigen_solver_eigenexa(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    integer :: i, dim, nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

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

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenexa:transpose', time_end - time_start_part)
    time_start_part = time_end

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

    call eigen_sx(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenexa:eigen_sx', time_end - time_start_part)
    call add_event('eigen_solver_eigenexa:total', time_end - time_start_part)
  end subroutine eigen_solver_eigenexa


  subroutine eigen_solver_eigenk(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    integer :: i, dim, nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

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

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenk:transpose', time_end - time_start_part)
    time_start_part = time_end

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

    call eigen_s(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenk:eigen_s', time_end - time_start_part)
    call add_event('eigen_solver_eigenk:total', time_end - time_start)
  end subroutine eigen_solver_eigenk
end module solver_eigenexa
