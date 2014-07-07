module solver_main
  use descriptor_parameters
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use eigenpairs_types, only: eigenpairs_types_union, eigenpairs_blacs
  use generalized_to_standard, only : reduce_generalized, recovery_generalized
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
    use solver_eigenexa, only : setup_distributed_matrix_for_eigenexa, &
         setup_distributed_matrix_for_general_eigenexa, eigen_solver_eigenexa

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: n, desc_A(desc_size), desc_B(desc_size), desc_A_re(desc_size)
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)

    integer :: base_count
    double precision :: times(10)

    n = arg%matrix_A_info%rows

    select case (trim(arg%solver_type))
    case ('lapack')
      call eigen_solver_lapack(matrix_A, eigenpairs)
    case ('scalapack_all')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
    case ('scalapack_select')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
    case ('general_scalapack_all')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('general_scalapack_select')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
      call recovery_generalized(n, arg%n_vec, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('eigenexa')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      call setup_distributed_matrix_for_eigenexa(n, desc_A, matrix_A_dist, eigenpairs)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_eigenexa(matrix_A_dist, desc_A, arg%n_vec, eigenpairs)
    case ('general_eigenexa')
      call get_wall_clock_base_count(base_count)
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      call get_wall_clock_time(base_count, times(1))
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call get_wall_clock_time(base_count, times(2))
      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs)
      call get_wall_clock_time(base_count, times(3))
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call get_wall_clock_time(base_count, times(4))
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call get_wall_clock_time(base_count, times(5))
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call get_wall_clock_time(base_count, times(6))
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      call get_wall_clock_time(base_count, times(7))
      call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, arg%n_vec, eigenpairs, 'L')
      call get_wall_clock_time(base_count, times(8))
      call pdgemr2d(n, n, matrix_B_dist, 1, 1, desc_B, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      call get_wall_clock_time(base_count, times(9))
      call recovery_generalized(n, n, matrix_A_redist, desc_A_re, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call get_wall_clock_time(base_count, times(10))
      if (check_master()) then
        print *, 'in solver_main: ', times
      end if
    case default
      stop '[Error] lib_eigen_solver: Unknown solver'
    end select
  end subroutine eigen_solver
end module solver_main
