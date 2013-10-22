module solver_main
  use descriptor_parameters
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : process, create_dense_matrix, &
       setup_distributed_matrix, gather_matrix, copy_global_sparse_matrix_to_local
  use eigenpairs_types, only: eigenpairs_types_union, eigenpairs_blacs
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
    type(eigenpairs_types_union), intent(in) :: eigenpairs

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
    res_norm_ave = sum(res_norm) / dble(arg%n_check_vec)
  end subroutine eigen_checker_local


  subroutine eigen_checker_blacs_standard(arg, matrix_A, eigenpairs, &
       res_norm_ave, res_norm_max)
    include 'mpif.h'

    type(argument), intent(in) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(eigenpairs_blacs), intent(inout) :: eigenpairs
    ! residual norm average, max
    double precision, intent(out) :: res_norm_ave, res_norm_max

    type(process) :: proc
    integer :: dim, desc_Residual(desc_size), desc_A(desc_size)
    double precision, allocatable :: Residual(:, :), matrix_A_dist(:, :)
    ! ave_and_max is declared as array due to usage of bcast
    ! 3rd element is for the index of the max value (discarded currently)
    double precision :: norm, ave_and_max(3)
    integer :: j, owner_proc_col, ierr
    integer :: indxg2p ! ScaLAPACK function

    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    dim = arg%matrix_A_info%rows
    call setup_distributed_matrix(proc, dim, dim, desc_A, matrix_A_dist)
    call copy_global_sparse_matrix_to_local(matrix_A, desc_A, matrix_A_dist)
    call setup_distributed_matrix(proc, dim, arg%n_check_vec, desc_Residual, Residual)

    ! Residual <- A * Eigenvectors
    call pdsymm('L', 'L', dim, arg%n_check_vec, 1.0d0, &
         matrix_A_dist, 1, 1, desc_A, &
         eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
         0.0d0, Residual, 1, 1, desc_Residual)

    ! For each column j, Residual(*, j) <- Residual(*, j) - eigenvalues(j) * Eigenvectors(*, j)
    ! Then store the 2-norm of the residual column in the first row of the column
    do j = 1, arg%n_check_vec
      call pdaxpy(dim, -1.0d0 * eigenpairs%values(j), &
           eigenpairs%Vectors, 1, j, eigenpairs%desc, 1, &
           Residual, 1, j, desc_Residual, 1)
      call pdnrm2(dim, norm, Residual, 1, j, desc_Residual, 1)

      owner_proc_col = indxg2p(j, desc_Residual(block_col_), 0, &
           desc_Residual(csrc_), proc%n_procs_col)
      if (proc%my_proc_row == 0 .and. proc%my_proc_col == owner_proc_col) then
        call pdelset(Residual, 1, j, desc_Residual, norm)
      end if
    end do

    ! Although these pd* routines return result values only in the scope of focused subvector,
    ! processor rank 0 must be in the scope because the subvector is the first row of Residual here
    call pdasum(arg%n_check_vec, ave_and_max(1), &
         Residual, 1, 1, desc_Residual, desc_Residual(rows_))
    call pdamax(arg%n_check_vec, ave_and_max(2), ave_and_max(3), &
         Residual, 1, 1, desc_Residual, desc_Residual(rows_))
    print *, proc%my_rank, ave_and_max(2)
    call mpi_bcast(ave_and_max(1), 3, mpi_double_precision, 0, mpi_comm_world, ierr)

    res_norm_ave = ave_and_max(1) / dble(arg%n_check_vec)
    res_norm_max = ave_and_max(2)
  end subroutine eigen_checker_blacs_standard


  subroutine lib_eigen_checker(arg, matrix_A, eigenpairs, &
       res_norm_ave, res_norm_max, matrix_B)
    type(argument), intent(in) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    ! residual norm average, max
    double precision, intent(out) :: res_norm_ave, res_norm_max

    logical :: is_master

    call check_master(arg%solver_type, is_master)

    if (eigenpairs%type_number == 2 .and. .not. arg%is_generalized_problem) then
      call eigen_checker_blacs_standard(arg, matrix_A, eigenpairs%blacs, &
           res_norm_ave, res_norm_max)
      return
    end if

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
    use matrix_io, only : sparse_mat
    use distribute_matrix, only : process, setup_distribution, &
         setup_distributed_matrix, copy_global_dense_matrix_to_local, &
         copy_global_sparse_matrix_to_local
    use eigenpairs_types, only : eigenpairs_types_union

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
      print *, 'info(pdpotrf): ', info
      stop '[Error] reduce_generalized: pdpotrf failed'
    end if
    ! Reduction to standard problem by A <- L^(-1) * A * L'^(-1)
    call pdsygst(1, 'L', dim, A, 1, 1, desc_A, B, 1, 1, desc_B, scale, info)
    if (info /= 0) then
      print *, 'info(pdsygst): ', info
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
      print *, 'info(pdtrtrs): ', info
      stop '[Error] reduce_generalized: pdtrtrs failed'
    end if
  end subroutine recovery_generalized
end module solver_main
