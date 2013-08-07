module solver_main
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : create_dense_matrix
  use eigenpairs_types, only: eigenpairs_types_union
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
    type(eigenpairs_types_union), intent(in) :: eigenpairs

    double precision, intent(out) :: res_norm_ave, res_norm_max

    double precision, allocatable :: a(:,:)
    double precision, pointer :: v(:,:)
    double precision, allocatable :: res(:) ! residual vector
    double precision, allocatable :: res_norm(:) ! residual norm

    integer :: j, dim, ierr

    if (present(matrix_B)) then
      print *, 'result checker for generalized eigenvalue problem is not implemeted yet'
      return
    end if

    dim = matrix_A%size

    allocate(res(dim), res_norm(arg%n_check_vec), a(dim, dim), stat = ierr)

    res_norm(:) = 0.0d0
    call create_dense_matrix(0, matrix_A, a)

    v => eigenpairs%local%vectors

    do j = 1, arg%n_check_vec
      res(:) = matmul(a, v(:, j)) - eigenpairs%local%values(j) * v(:, j)
      res_norm(j) = sqrt(abs(dot_product(res, res) / dot_product(v(:, j), v(:, j))))
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
    type(eigenpairs_types_union), intent(in) :: eigenpairs
    ! residual norm average, max
    double precision, intent(out) :: res_norm_ave, res_norm_max

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

    integer :: n, desc_A(9), desc_B(9), info
    double precision :: scale
    type(process) :: proc
    real(kind(1.d0)), allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :)

    n = arg%matrix_A_info%rows

    if ((arg%n_vec < 0) .or. (arg%n_vec > n)) then
      write(*,*) 'Error(lib_eigen_solver): n, n_vec=', n, arg%n_vec
      stop
    endif

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
      call pdpotrf('L', n, matrix_B_dist, 1, 1, desc_B, info)
      call pdsygst(1, 'L', n, matrix_A_dist, 1, 1, desc_A, matrix_B_dist, &
           1, 1, desc_B, scale, info)
      !call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
!      call pdtrtrs('L', 'T', 'N', n, n, matrix_A_dist, 1, 1, desc_A, &

    case ('eigenexa')
      stop 'Eigen Exa is not supported yet'
      !call eigen_solver_eigenexa(matrix_A, arg%n_vec, eigenpairs)
    case default
      write(*,*) 'Error(lib_eigen_solver): solver type=', trim(arg%solver_type)
      stop
    end select

  end subroutine lib_eigen_solver
end module solver_main
