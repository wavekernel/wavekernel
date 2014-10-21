module solver_eigenexa
  use distribute_matrix, only : desc_size
  use processes, only : terminate
  use matrix_io, only : sparse_mat
  use eigenpairs_types, only : eigenpairs_types_union

  implicit none

  private
  public :: setup_distributed_matrix_for_eigenexa, &
       setup_distributed_matrix_for_general_eigenexa, eigen_solver_eigenexa, eigen_solver_eigenk

contains

  subroutine setup_distributed_matrix_for_eigenexa(dim, desc_A, matrix_A, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) ::desc_A(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('lib_eigen_solver: EigenExa is not supported in this build', 1)
  end subroutine setup_distributed_matrix_for_eigenexa


  subroutine eigen_solver_eigenexa(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(in) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(out) :: eigenpairs
    character, intent(in), optional :: uplo

    call terminate('lib_eigen_solver: EigenExa is not supported in this build', 1)
  end subroutine eigen_solver_eigenexa


  subroutine setup_distributed_matrix_for_general_eigenexa( &
       dim, desc_A, matrix_A, desc_B, matrix_B, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) :: desc_A(desc_size), desc_B(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :), matrix_B(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('lib_eigen_solver: EigenExa is not supported in this build', 1)
  end subroutine setup_distributed_matrix_for_general_eigenexa


  subroutine setup_distributed_matrix_for_general_eigenk( &
       dim, desc_A, matrix_A, desc_B, matrix_B, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) :: desc_A(desc_size), desc_B(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :), matrix_B(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('lib_eigen_solver: EigenExa is not supported in this build', 1)
  end subroutine setup_distributed_matrix_for_general_eigenk
end module solver_eigenexa
