module solver_main
  use descriptor_parameters
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use eigenpairs_types, only: eigenpairs_types_union, eigenpairs_blacs
  use generalized_to_standard, only : reduce_generalized, recovery_generalized
  use global_variables, only : g_block_size
  use processes, only : process, setup_distribution, print_map_of_grid_to_processes, &
       check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time

  implicit none

  private
  public :: eigen_solver

contains

  subroutine eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
    use solver_lapack, only : eigen_solver_lapack
    use solver_scalapack_all, only : eigen_solver_scalapack_all
    use solver_scalapack_select, only : eigen_solver_scalapack_select
    use solver_eigenexa, only : setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa, eigen_solver_eigenk
    use solver_elpa
    use solver_elpa_eigenexa
    include 'mpif.h'

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: n, desc_A(desc_size), desc_B(desc_size), desc_A_re(desc_size)
    integer :: ierr, mpierr
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs_tmp

    integer :: base_count
    double precision :: times(10)

    n = arg%matrix_A_info%rows
    if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
    if (arg%block_size > 0) then  ! Do not use the default block size.
      g_block_size = arg%block_size
    end if

    if (trim(arg%solver_type) /= 'lapack') then
      call setup_distribution(proc)
    end if

    select case (trim(arg%solver_type))
    case ('lapack')
      call eigen_solver_lapack(matrix_A, eigenpairs)
    case ('scalapack')
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
    case ('scalapack_select')
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
    case ('general_scalapack')
      call mpi_barrier(mpi_comm_world, mpierr)
      times(1) = mpi_wtime()
      call get_wall_clock_base_count(base_count)
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(2) = mpi_wtime()
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(3) = mpi_wtime()
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(4) = mpi_wtime()
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(5) = mpi_wtime()
      if (check_master()) then
        print *, 'general_scalapack elapsed time (printed in solver_main): '
        print *, 'general_scalapack setup_matrix           : ', times(2) - times(1)
        print *, 'general_scalapack reduce_generalized     : ', times(3) - times(2)
        print *, 'general_scalapack eigen_solver_scalapack : ', times(4) - times(3)
        print *, 'general_scalapack recovery_generalized   : ', times(5) - times(4)
        print *, 'general_scalapack total                  : ', times(5) - times(1)
      end if
    case ('general_scalapack_select')
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
      call recovery_generalized(n, arg%n_vec, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('eigenexa')
      call setup_distribution(proc)
      call setup_distributed_matrix_for_eigenexa(n, desc_A, matrix_A_dist, eigenpairs)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_eigenexa(matrix_A_dist, desc_A, arg%n_vec, eigenpairs)
    case ('general_scalapack_eigenexa')
      call mpi_barrier(mpi_comm_world, mpierr)
      times(1) = mpi_wtime()
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(2) = mpi_wtime()
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(3) = mpi_wtime()
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      deallocate(matrix_A_dist)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(4) = mpi_wtime()
      call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, arg%n_vec, eigenpairs_tmp, 'L')
      call mpi_barrier(mpi_comm_world, mpierr)
      times(5) = mpi_wtime()
      deallocate(matrix_A_redist)
      eigenpairs%type_number = 2
      allocate(eigenpairs%blacs%values(n), stat = ierr)
      if (ierr /= 0) then
        call terminate('eigen_solver, general_scalapack_eigenexa: allocation failed', ierr)
      end if
      eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
      eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
           eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
           eigenpairs_tmp%blacs%desc(context_))
            call mpi_barrier(mpi_comm_world, mpierr)
      times(6) = mpi_wtime()
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(7) = mpi_wtime()
      if (check_master()) then
        print *, 'general_scalapack_eigenexa elapsed time (printed in solver_main): '
        print *, 'general_scalapack_eigenexa setup_matrix          : ', times(2) - times(1)
        print *, 'general_scalapack_eigenexa reduce_generalized    : ', times(3) - times(2)
        print *, 'general_scalapack_eigenexa pdgemr2d_A            : ', times(4) - times(3)
        print *, 'general_scalapack_eigenexa eigen_solver_eigenexa : ', times(5) - times(4)
        print *, 'general_scalapack_eigenexa pdgemr2d_B            : ', times(6) - times(5)
        print *, 'general_scalapack_eigenexa recovery_generalized  : ', times(7) - times(6)
        print *, 'general_scalapack_eigenexa total                 : ', times(7) - times(1)
      end if
    case ('general_scalapack_eigenk')
      call mpi_barrier(mpi_comm_world, mpierr)
      times(1) = mpi_wtime()
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(2) = mpi_wtime()
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(3) = mpi_wtime()
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      deallocate(matrix_A_dist)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(4) = mpi_wtime()
      call eigen_solver_eigenk(matrix_A_redist, desc_A_re, arg%n_vec, eigenpairs_tmp, 'L')
      call mpi_barrier(mpi_comm_world, mpierr)
      times(5) = mpi_wtime()
      deallocate(matrix_A_redist)
      eigenpairs%type_number = 2
      allocate(eigenpairs%blacs%values(n), stat = ierr)
      if (ierr /= 0) then
        call terminate('eigen_solver, general_scalapack_eigenk: allocation failed', ierr)
      end if
      eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
      eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
           eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
           eigenpairs_tmp%blacs%desc(context_))
            call mpi_barrier(mpi_comm_world, mpierr)
      times(6) = mpi_wtime()
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call mpi_barrier(mpi_comm_world, mpierr)
      times(7) = mpi_wtime()
      if (check_master()) then
        print *, 'general_scalapack_eigenk elapsed time (printed in solver_main): '
        print *, 'general_scalapack_eigenk setup_matrix         : ', times(2) - times(1)
        print *, 'general_scalapack_eigenk reduce_generalized   : ', times(3) - times(2)
        print *, 'general_scalapack_eigenk pdgemr2d_A           : ', times(4) - times(3)
        print *, 'general_scalapack_eigenk eigen_solver_eigenk  : ', times(5) - times(4)
        print *, 'general_scalapack_eigenk pdgemr2d_B           : ', times(6) - times(5)
        print *, 'general_scalapack_eigenk recovery_generalized : ', times(7) - times(6)
        print *, 'general_scalapack_eigenk total                : ', times(7) - times(1)
      end if
    case ('general_elpa_scalapack')
      call solve_with_general_elpa_scalapck(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa1')
      call solve_with_general_elpa1(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa2')
      call solve_with_general_elpa2(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa_eigenexa')
      call solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa_eigenk')
      call solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    case default
      call terminate('eigen_solver: Unknown solver', 1)
    end select
  end subroutine eigen_solver
end module solver_main
