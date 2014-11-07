module solver_elpa_eigenexa
  use global_variables, only : g_block_size
  use distribute_matrix, only : process, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use matrix_io, only : sparse_mat
  use processes, only : check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time

  implicit none

  !private

contains

  subroutine solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('', 1)
    end subroutine solve_with_general_elpa_eigenexa


    subroutine solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('', 1)
    end subroutine solve_with_general_elpa_eigenk
end module solver_elpa_eigenexa
