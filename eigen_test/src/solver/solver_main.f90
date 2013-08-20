module solver_main
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : create_dense_matrix, gather_matrix
  use eigenpairs_types, only: eigenpairs_types_union
  use processes, only : check_master
  implicit none

  private
  public :: lib_eigen_solver
  public :: lib_eigen_checker

contains
  subroutine eigen_checker_local(arg, matrix_A, eigenpairs, &
       res_norm_ave, res_norm_max, matrix_B)
    type(argument), intent(in) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(inout) :: eigenpairs

    double precision, intent(out) :: res_norm_ave, res_norm_max

    double precision, allocatable :: a(:,:)
    double precision, allocatable :: b(:,:)
    double precision, allocatable :: left(:), right(:) ! residual = left - right
    double precision, allocatable :: res_norm(:) ! residual norm

    integer :: j, dim

    dim = matrix_A%size

    allocate(left(dim), right(dim), res_norm(arg%n_check_vec))
    res_norm(:) = 0.0d0

    call create_dense_matrix(matrix_A, a)
    if (arg%is_generalized_problem) then
      call create_dense_matrix(matrix_B, b)
    end if

    do j = 1, arg%n_check_vec
      left(:) = matmul(a, eigenpairs%local%vectors(:, j))
      if (arg%is_generalized_problem) then
        right(:) = eigenpairs%local%values(j) * &
             matmul(b, eigenpairs%local%vectors(:, j))
      else
        right(:) = eigenpairs%local%values(j) * eigenpairs%local%vectors(:, j)
      end if
      res_norm(j) = sqrt(abs(dot_product(left - right, left - right) / &
           dot_product(eigenpairs%local%vectors(:, j), &
           eigenpairs%local%vectors(:, j))))
    enddo

    res_norm_max = maxval(res_norm)
    res_norm_ave = sum(res_norm) / dble(arg%n_vec)
  end subroutine eigen_checker_local


  subroutine lib_eigen_checker(arg, matrix_A, eigenpairs, &
       res_norm_ave, res_norm_max, matrix_B)
    implicit none

    type(argument), intent(in) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    ! residual norm average, max
    double precision, intent(out) :: res_norm_ave, res_norm_max

    logical :: is_master

    call check_master(arg%solver_type, is_master)

    if (eigenpairs%type_number == 2) then ! Temporal implementation
      eigenpairs%type_number = 1

      allocate(eigenpairs%local%values(matrix_A%size))
      eigenpairs%local%values(:) = eigenpairs%blacs%values(:)

      if (is_master) then
        allocate(eigenpairs%local%vectors(matrix_A%size, matrix_A%size))
        eigenpairs%local%vectors(:, :) = 0.0d0
      end if

      call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
           0, 0, eigenpairs%local%vectors)

      deallocate(eigenpairs%blacs%values, eigenpairs%blacs%Vectors)

      if (.not. is_master) then
        return
      end if
    end if

    if (eigenpairs%type_number == 1) then
      call eigen_checker_local(arg, matrix_A, eigenpairs, &
           res_norm_ave, res_norm_max, matrix_B)
    else
      print *, 'result checker for distributed output is not implemeted yet'
    end if
  end subroutine lib_eigen_checker


  subroutine lib_eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
    use command_argument, only : argument
    use solver_lapack, only : eigen_solver_lapack
    use solver_scalapack_all, only : eigen_solver_scalapack_all
    use solver_scalapack_select, only : eigen_solver_scalapack_select
    use solver_eigenexa, only : eigen_solver_eigenexa
    use matrix_io, only : sparse_mat
    use distribute_matrix, only : process, setup_distribution, &
         setup_distributed_matrix, copy_global_dense_matrix_to_local, &
         copy_global_sparse_matrix_to_local
    use eigenpairs_types, only : eigenpairs_types_union
    implicit none

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: n, desc_A(9), desc_B(9)
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :)

    n = arg%matrix_A_info%rows

    select case (trim(arg%solver_type))
    case ('lapack')
      call eigen_solver_lapack(matrix_A, eigenpairs)
    case ('scalapack_all')
      call setup_distribution(proc)
      call setup_distributed_matrix(proc, n, n, desc_A, matrix_A_dist)
      call copy_global_sparse_matrix_to_local(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
    case ('scalapack_select')
      call setup_distribution(proc)
      call setup_distributed_matrix(proc, n, n, desc_A, matrix_A_dist)
      call copy_global_sparse_matrix_to_local(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
    case ('general_scalapack_all')
      call setup_distribution(proc)
      call setup_distributed_matrix(proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix(proc, n, n, desc_B, matrix_B_dist)
      call copy_global_sparse_matrix_to_local(matrix_A, desc_A, matrix_A_dist)
      call copy_global_sparse_matrix_to_local(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('general_scalapack_select')
      call setup_distribution(proc)
      call setup_distributed_matrix(proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix(proc, n, n, desc_B, matrix_B_dist)
      call copy_global_sparse_matrix_to_local(matrix_A, desc_A, matrix_A_dist)
      call copy_global_sparse_matrix_to_local(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
      call recovery_generalized(n, arg%n_vec, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('eigenexa')
      stop '[Error] lib_eigen_solver: Eigen Exa is not supported yet'
      !call eigen_solver_eigenexa(matrix_A, arg%n_vec, eigenpairs)
    case default
      stop '[Error] lib_eigen_solver: Unknown solver'
    end select
  end subroutine lib_eigen_solver


  subroutine reduce_generalized(dim, A, desc_A, B, desc_B)
    integer, intent(in) :: dim, desc_A(9), desc_B(9)
    double precision, intent(inout) :: A(:, :), B(:, :)

    integer :: info
    double precision :: scale

    ! B = LL', overwritten to B
    call pdpotrf('L', dim, B, 1, 1, desc_B, info)
    if (info /= 0) then
      stop '[Error] reduce_generalized: pdpotrf failed'
    end if
    ! Reduction to standard problem by A <- L^(-1) * A * L'^(-1)
    call pdsygst(1, 'L', dim, A, 1, 1, desc_A, B, 1, 1, desc_B, scale, info)
    if (info /= 0) then
      stop '[Error] reduce_generalized: pdsygst failed'
    end if
  end subroutine reduce_generalized


  subroutine recovery_generalized(dim, n_vec, B, desc_B, Vectors, desc_Vectors)
    integer, intent(in) :: dim, n_vec, desc_B(9), desc_Vectors(9)
    double precision, intent(in) :: B(:, :)
    double precision, intent(inout) :: Vectors(:, :)

    integer :: info

    ! Recovery eigenvectors by V <- L'^(-1) * V
    call pdtrtrs('L', 'T', 'N', dim, n_vec, B, 1, 1, desc_B, &
         Vectors, 1, 1, desc_Vectors, info)
    if (info /= 0) then
      stop '[Error] reduce_generalized: pdtrtrs failed'
    end if
  end subroutine recovery_generalized
end module solver_main
