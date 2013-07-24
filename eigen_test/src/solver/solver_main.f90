module solver_main
  use distribute_matrix, only : create_dense_matrix !(routine)
  implicit none

  private
  public :: lib_eigen_solver
  public :: lib_eigen_checker

contains
  subroutine lib_eigen_checker(a, n_vec, n_check_vec, eigen_level, v, rn_ave, rn_max)
    implicit none

    real(kind(1.d0)), intent(in) :: a(:,:)
    real(kind(1.d0)), intent(in) :: v(:,:)
    real(kind(1.d0)), intent(in) :: eigen_level(:)
    integer, intent(in) :: n_vec, n_check_vec
    real(kind(1.d0)), intent(out) :: rn_ave, rn_max      ! residual norm average, max

    real(kind(1.d0)), allocatable :: r(:)   ! residual vector (work array)
    real(kind(1.d0)), allocatable :: rn(:)  ! residual norm

    integer :: n, ierr, j

    n = size(a, 1)

    if (n_check_vec < 1 .or. n_check_vec > n_vec) then
      print *, 'ERROR(lib_eigen_checker): n, n_vec, n_check_vec = ', n, n_vec, n_check_vec
      stop
    endif

    allocate(r(n), stat=ierr)
    if (ierr /= 0) stop 'Alloc error (la_eigen_solver_check) for w'

    allocate(rn(n_check_vec), stat=ierr)
    if (ierr /= 0) stop 'Alloc error (la_eigen_solver_check) for rn'
    rn(:) = 0.0d0

    do j=1,n_check_vec
      r(:) = matmul(a, v(:, j)) - eigen_level(j) * v(:, j)  ! redisual vector
      rn(j) = sqrt(abs(dot_product(r, r))) / sqrt(abs(dot_product(v(:, j), v(:, j))))
    enddo

    rn_max = maxval(rn)
    rn_ave = sum(rn) / dble(n_vec)
  end subroutine lib_eigen_checker


  subroutine lib_eigen_solver(arg, matrix_A, eigenvalues, eigenvectors, matrix_B)
    use command_argument, only : argument
    use solver_lapack, only : eigen_solver_lapack
    use solver_scalapack_all, only : eigen_solver_scalapack_all
    use solver_scalapack_select, only : eigen_solver_scalapack_select
    use solver_eigenexa, only : eigen_solver_eigenexa
    use matrix_io, only : sparse_mat
    use distribute_matrix, only : conf_distribution, setup_distribution, &
         setup_distributed_matrix, copy_global_dense_matrix_to_local, &
         copy_global_sparse_matrix_to_local
    implicit none

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    real(kind(1.d0)), intent(out), allocatable :: eigenvalues(:), eigenvectors(:, :)

    integer :: n, desc(9)
    type(conf_distribution) :: conf
    real(kind(1.d0)), allocatable :: a(:, :), mat_dist(:, :)

    n = arg%matrix_A_info%rows

    if ((arg%n_vec < 0) .or. (arg%n_vec > n)) then
      write(*,*) 'Error(lib_eigen_solver): n, n_vec=', n, arg%n_vec
      stop
    endif

    allocate(eigenvalues(n))
    eigenvalues(:) = 0.0d0
    allocate(eigenvectors(n, n))
    eigenvectors(:, :) = 0.0d0

    select case (trim(arg%solver_type))
    case ('lapack')
      call create_dense_matrix(0, matrix_A, a)
      call eigen_solver_lapack(a, eigenvalues)
      eigenvectors = a
    case ('scalapack_all')
      if (arg%n_vec < n) then
        stop 'Error(lib_eigen_solver): this solver does not support partial eigenvalue computation'
      end if
      call setup_distribution(n, conf)
      call setup_distributed_matrix(conf, desc, mat_dist)
      call copy_global_sparse_matrix_to_local(matrix_A, desc, mat_dist)
      call eigen_solver_scalapack_all(conf, desc, mat_dist, eigenvalues, eigenvectors)
    case ('scalapack_select')
      call setup_distribution(n, conf)
      call setup_distributed_matrix(conf, desc, mat_dist)
      call copy_global_sparse_matrix_to_local(matrix_A, desc, mat_dist)
      call eigen_solver_scalapack_select(conf, desc, mat_dist, arg%n_vec, &
           eigenvalues, eigenvectors)
    case ('eigenexa')
      call eigen_solver_eigenexa(matrix_A, arg%n_vec, eigenvalues, eigenvectors)
    case default
      write(*,*) 'Error(lib_eigen_solver): solver type=', trim(arg%solver_type)
      stop
    end select

  end subroutine lib_eigen_solver
end module solver_main
