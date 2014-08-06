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
    use solver_eigenexa, only : setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: n, desc_A(desc_size), desc_B(desc_size), desc_A_re(desc_size)
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs_tmp

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
      call get_wall_clock_base_count(base_count)
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
      call get_wall_clock_time(base_count, times(1))
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call get_wall_clock_time(base_count, times(2))
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
      call get_wall_clock_time(base_count, times(3))
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call get_wall_clock_time(base_count, times(4))
      if (check_master()) then
        print *, 'general_scalapack_all elapsed time (printed in solver_main): '
        print *, 'general_scalapack_all setup_matrix: ', times(1)
        print *, 'general_scalapack_all reduce_generalized: ', times(2)
        print *, 'general_scalapack_all eigen_solver_scalapack: ', times(3)
        print *, 'general_scalapack_all recovery_generalized: ', times(4)
      end if
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
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call get_wall_clock_time(base_count, times(1))
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call get_wall_clock_time(base_count, times(2))
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      call get_wall_clock_time(base_count, times(3))
      call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, arg%n_vec, eigenpairs_tmp, 'L')
      call get_wall_clock_time(base_count, times(4))
      eigenpairs%type_number = 2
      allocate(eigenpairs%blacs%values(n))
      eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
      eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
           eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
           eigenpairs_tmp%blacs%desc(context_))
      call get_wall_clock_time(base_count, times(5))
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call get_wall_clock_time(base_count, times(6))
      if (check_master()) then
        print *, 'general_eigenexa elapsed time (printed in solver_main): '
        print *, 'general_eigenexa setup_matrix: ', times(1)
        print *, 'general_eigenexa reduce_generalized: ', times(2)
        print *, 'general_eigenexa pdgemr2d_A: ', times(3)
        print *, 'general_eigenexa eigen_solver_eigenexa: ', times(4)
        print *, 'general_eigenexa pdgemr2d_B: ', times(5)
        print *, 'general_eigenexa recovery_generalized: ', times(6)
      end if
    case default
      stop '[Error] lib_eigen_solver: Unknown solver'
    end select
  end subroutine eigen_solver
end module solver_main
