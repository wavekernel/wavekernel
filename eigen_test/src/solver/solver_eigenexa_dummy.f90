module solver_eigenexa
  use distribute_matrix, only : desc_size
  use processes, only : terminate
  use matrix_io, only : sparse_mat
  use eigenpairs_types, only : eigenpairs_types_union

  implicit none

  private
  public :: setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa

contains

  subroutine setup_distributed_matrix_for_eigenexa(dim, desc_A, matrix_A, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) ::desc_A(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs
  end subroutine setup_distributed_matrix_for_eigenexa


  subroutine eigen_solver_eigenexa(mat, dim, n_vec, eigenpairs)
    double precision, intent(in) :: mat(:, :)
    integer, intent(in) :: dim, n_vec
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('[Error] lib_eigen_solver: EigenExa is not supported in this build')
  end subroutine eigen_solver_eigenexa
end module solver_eigenexa
