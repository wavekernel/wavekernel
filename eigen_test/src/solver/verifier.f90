module verifier
  use command_argument, only : argument
  use descriptor_parameters
  use distribute_matrix, only : convert_sparse_matrix_to_dense, &
       setup_distributed_matrix, distribute_global_sparse_matrix
  use eigenpairs_types, only: eigenpairs_types_union, eigenpairs_local, &
       eigenpairs_blacs
  use global_variables, only : g_block_size
  use matrix_io, only : sparse_mat
  use processes, only : process, terminate
  implicit none

  private
  public :: eval_residual_norm, eval_orthogonality

contains

  subroutine eval_residual_norm_local(arg, matrix_A, eigenpairs, &
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

    integer :: j, dim, ierr

    dim = matrix_A%size

    allocate(left(dim), right(dim), res_norm(arg%n_check_vec), stat = ierr)
    if (ierr /= 0) then
      call terminate('eval_residual_norm_local: allocation failed', ierr)
    end if
    res_norm(:) = 0.0d0

    call convert_sparse_matrix_to_dense(matrix_A, a)
    if (arg%is_generalized_problem) then
      call convert_sparse_matrix_to_dense(matrix_B, b)
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
  end subroutine eval_residual_norm_local


  subroutine eval_residual_norm_blacs(arg, matrix_A, eigenpairs, &
       res_norm_ave, res_norm_max, matrix_B)
    include 'mpif.h'

    type(argument), intent(in) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(eigenpairs_blacs), intent(inout) :: eigenpairs
    type(sparse_mat), intent(in), optional :: matrix_B
    ! residual norm average, max
    double precision, intent(out) :: res_norm_ave, res_norm_max

    type(process) :: proc
    integer :: dim, desc_Residual(desc_size), desc_A(desc_size), desc_B(desc_size)
    double precision, allocatable :: Residual(:, :), matrix_A_dist(:, :), matrix_B_dist(:, :)
    ! ave_and_max is declared as array due to usage of bcast
    ! 3rd element is for the index of the max value (discarded currently)
    double precision :: res_norm, A_norm, ave_and_max(3)
    integer :: j, block_size, owner_proc_col, ierr
    ! ScaLAPACK functions
    integer :: indxg2p
    double precision :: pdlange

    if (arg%is_generalized_problem .and. .not. present(matrix_B)) then
      call terminate('eval_residual_norm_blacs: matrix_B is not provided', 1)
    end if

    if (trim(arg%solver_type) == 'eigenexa' .or. &
         trim(arg%solver_type) == 'general_eigenexa') then
      block_size = 1
    else if (arg%block_size > 0) then  ! Do not use the default block size
      block_size = arg%block_size
    else
      block_size = g_block_size
    end if

    ! call blacs_get(-1, 0, proc%context)
    ! Because context acquiring by blacs_get can fail in some environments,
    ! use context in the descriptor of eigenpairs instead.
    proc%context = eigenpairs%desc(context_)
    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    dim = arg%matrix_A_info%rows

    call setup_distributed_matrix('A', proc, dim, dim, desc_A, matrix_A_dist, &
         block_size = block_size)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    A_norm = pdlange('F', dim, dim, matrix_A_dist, 1, 1, desc_A, 0)
    print *, 'A_norm: ', A_norm

    if (arg%is_generalized_problem) then
      call setup_distributed_matrix('B', proc, dim, dim, &
           desc_B, matrix_B_dist, block_size = block_size)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
    end if

    call setup_distributed_matrix('Residual', proc, dim, arg%n_check_vec, &
         desc_Residual, Residual, block_size = block_size)

    if (arg%is_generalized_problem) then
      ! Residual <- B * Eigenvectors
      call pdsymm('L', 'L', dim, arg%n_check_vec, &
           1.0d0, matrix_B_dist, 1, 1, desc_B, &
           eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
           0.0d0, Residual, 1, 1, desc_Residual)
    else
      ! Residual <- Eigenvectors
      call pdlacpy('A', dim, arg%n_check_vec, &
           eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
           Residual, 1, 1, desc_Residual)
    end if

    ! For each column j, Residual(*, j) <- Residual(*, j) * -eigenvalues(j)
    do j = 1, arg%n_check_vec
      call pdscal(dim, -1.0d0 * eigenpairs%values(j), &
           Residual, 1, j, desc_Residual, 1)
    end do

    ! Residual <- Residual + A * Eigenvectors
    call pdsymm('L', 'L', dim, arg%n_check_vec, &
         1.0d0, matrix_A_dist, 1, 1, desc_A, &
         eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
         1.0d0, Residual, 1, 1, desc_Residual)

    ! Store 2-norm of each residual column in the first row of the column
    do j = 1, arg%n_check_vec
      call pdnrm2(dim, res_norm, Residual, 1, j, desc_Residual, 1)
      owner_proc_col = indxg2p(j, desc_Residual(block_col_), 0, &
           desc_Residual(csrc_), proc%n_procs_col)
      if (proc%my_proc_row == 0 .and. proc%my_proc_col == owner_proc_col) then
        call pdelset(Residual, 1, j, desc_Residual, res_norm)
      end if
    end do

    ! Although these pd* routines return result values only in the scope of focused subvector,
    ! processor rank 0 must be in the scope because the subvector is the first row of Residual here
    call pdasum(arg%n_check_vec, ave_and_max(1), &
         Residual, 1, 1, desc_Residual, desc_Residual(rows_))
    call pdamax(arg%n_check_vec, ave_and_max(2), ave_and_max(3), &
         Residual, 1, 1, desc_Residual, desc_Residual(rows_))
    call mpi_bcast(ave_and_max(1), 3, mpi_double_precision, 0, mpi_comm_world, ierr)

    res_norm_ave = ave_and_max(1) / A_norm / dble(arg%n_check_vec)
    res_norm_max = ave_and_max(2) / A_norm
  end subroutine eval_residual_norm_blacs


  subroutine eval_residual_norm(arg, matrix_A, eigenpairs, &
       res_norm_ave, res_norm_max, matrix_B)
    type(argument), intent(in) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    ! residual norm average, max
    double precision, intent(out) :: res_norm_ave, res_norm_max

    if (eigenpairs%type_number == 1) then
      call eval_residual_norm_local(arg, matrix_A, eigenpairs, &
           res_norm_ave, res_norm_max, matrix_B)
    else if (eigenpairs%type_number == 2) then
      if (arg%is_generalized_problem) then
        call eval_residual_norm_blacs(arg, matrix_A, eigenpairs%blacs, &
             res_norm_ave, res_norm_max, matrix_B)
      else
        call eval_residual_norm_blacs(arg, matrix_A, eigenpairs%blacs, &
             res_norm_ave, res_norm_max)
      end if
    else
      print '("[Warning] eval_residual_norm: output of this type is not supported yet")'
    end if
  end subroutine eval_residual_norm


  ! For each pair of vectors within the index range [index1, index2],
  ! computes cosine of the angle between them, and report the average and
  ! the maximum (with the index pair) in the cosine values.
  subroutine eval_orthogonality_local(index1, index2, eigenpairs, &
       cos_ave, cos_max, cos_max_index1, cos_max_index2)
    integer, intent(in) :: index1, index2
    type(eigenpairs_local), intent(in) :: eigenpairs
    double precision, intent(out) :: cos_ave, cos_max
    integer, intent(out) :: cos_max_index1, cos_max_index2

    integer :: dim, n, i, j, k, ierr
    double precision :: dot
    double precision, allocatable :: norms(:), coss(:)

    double precision :: ddot, dnrm2  ! Functions.

    dim = size(eigenpairs%vectors, 1)
    n = index2 - index1 + 1
    allocate(norms(n), coss(n * (n - 1) / 2), stat = ierr)
    if (ierr /= 0) then
      call terminate('eval_orthogonality_local: allocation failed', ierr)
    end if

    do i = 1, n
      norms(i) = dnrm2(dim, eigenpairs%vectors(index1 + i - 1, :), 1)
    end do

    ! Computes inner products for all index pairs.
    ! Conversion between an index pair and an entire pair is:
    ! (i, j) -> k = (j - 1) * (j - 2) / 2 + i
    ! k -> j = floor((3 + sqrt(8 * k - 6)) / 2)
    do j = 2, n
      do i = 1, j - 1
        k = (j - 1) * (j - 2) / 2 + i
        dot = ddot(dim, eigenpairs%vectors(i + index1 - 1, :), 1, &
        eigenpairs%vectors(j + index1 - 1, :), 1)
        coss(k) = abs(dot) / (norms(i) * norms(j))
      end do
    end do

    ! Finds max value and its index.
    k = maxloc(coss, dim = 1)
    cos_max = coss(k)
    cos_max_index2 = int((3.0d0 + sqrt(8.0d0 * dble(k) - 6.0d0)) / 2.0d0)
    cos_max_index1 = k - (cos_max_index2 - 1) * (cos_max_index2 - 2) / 2
    cos_max_index2 = cos_max_index2 + index1 - 1
    cos_max_index1 = cos_max_index1 + index1 - 1

    cos_ave = sum(coss) / dble(n * (n - 1) / 2)
  end subroutine eval_orthogonality_local


  ! Distributed parallel version of subroutine eval_orthogonality_local.
  ! Returns correct values only in the master processor (coordinate (0, 0)).
  subroutine eval_orthogonality_blacs(index1, index2, eigenpairs, &
       cos_ave, cos_max, cos_max_index1, cos_max_index2)
    include 'mpif.h'

    integer, intent(in) :: index1, index2
    type(eigenpairs_blacs), intent(in) :: eigenpairs
    double precision, intent(out) :: cos_ave, cos_max
    integer, intent(out) :: cos_max_index1, cos_max_index2

    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col, context  ! Process info.
    integer :: n, block_row, block_col, dim  ! Global matrix info.
    integer :: local_cols, local_rows  ! Local matrix info.
    integer :: diag, i, j, owner_proc_row, owner_proc_col, max_in_col_index, info, ierr
    double precision :: max_in_col, scale, sum_in_col
    double precision, allocatable :: InnerProducts(:, :), Cos_maxes(:), Cos_sums(:)
    integer, allocatable :: Cos_maxes_index(:)
    integer :: desc_InnerProducts(desc_size), desc_Cos_maxes(desc_size), desc_Cos_maxes_index(desc_size), desc_Cos_sums(desc_size)

    integer :: indxg2p, indxg2l, numroc, blacs_pnum  ! Functions.

    dim = eigenpairs%desc(rows_)  ! Dimension of the original problem.
    n = index2 - index1 + 1  ! The number of vectors to be checked.
    context = eigenpairs%desc(context_)

    call blacs_gridinfo(eigenpairs%desc(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    block_row = eigenpairs%desc(block_row_)
    block_col = eigenpairs%desc(block_col_)
    local_rows = max(1, numroc(n, block_row, my_proc_row, 0, n_procs_row))
    local_cols = max(1, numroc(n, block_col, my_proc_col, 0, n_procs_col))
    call descinit(desc_InnerProducts, n, n, block_row, block_col, 0, 0, context, local_rows, info)
    call descinit(desc_Cos_maxes, 1, n, 1, block_col, 0, 0, context, 1, info)
    call descinit(desc_Cos_maxes_index, 1, n, 1, block_col, 0, 0, context, 1, info)
    call descinit(desc_Cos_sums, 1, n, 1, block_col, 0, 0, context, 1, info)
    allocate(InnerProducts(local_rows, local_cols), stat = ierr)
    if (ierr /= 0) then
      call terminate('eval_orthogonality_blacs: allocation failed', ierr)
    end if
    if (my_proc_row == 0) then
      allocate(Cos_maxes(local_cols), Cos_maxes_index(local_cols), Cos_sums(local_cols), stat = ierr)
      if (ierr /= 0) then
        call terminate('eval_orthogonality_blacs: allocation failed', ierr)
      end if
    end if

    call pdgemm('T', 'N', n, n, dim, 1.0d0, &
         eigenpairs%Vectors, 1, index1, eigenpairs%desc, &
         eigenpairs%Vectors, 1, index1, eigenpairs%desc, &
         0.0d0, InnerProducts, 1, 1, desc_InnerProducts)

    ! Scales elements to obtain cosine values.
    do diag = 1, n
      owner_proc_row = indxg2p(diag, block_row, 0, 0, n_procs_row)
      owner_proc_col = indxg2p(diag, block_col, 0, 0, n_procs_col)
      if (my_proc_col == owner_proc_col .and. my_proc_row == owner_proc_row) then
        i = indxg2l(diag, block_row, 0, 0, n_procs_row)
        j = indxg2l(diag, block_col, 0, 0, n_procs_col)
        scale = 1.0d0 / sqrt(InnerProducts(i, j))  ! Inverse of norm.
      end if
      call mpi_bcast(scale, 1, mpi_double_precision, &
           blacs_pnum(context, owner_proc_row, owner_proc_col), mpi_comm_world, info)
      call pdscal(n, scale, InnerProducts, 1, diag, desc_InnerProducts, 1)
      call pdscal(n, scale, InnerProducts, diag, 1, desc_InnerProducts, n)
    end do

    do diag = 1, n - 1
      ! Finds the maximum absolute cosine value and its index in each column
      call pdamax(n - diag, max_in_col, max_in_col_index, InnerProducts, diag + 1, diag, desc_InnerProducts, 1)
      ! Computes sum of absolute cosine values in each column
      call pdasum(n - diag, sum_in_col, InnerProducts, diag + 1, diag, desc_InnerProducts, 1)
      owner_proc_col = indxg2p(diag, block_col, 0, 0, n_procs_col)
      if (my_proc_col == owner_proc_col .and. my_proc_row == 0) then
        ! This process owns the result of the pdamax/pdasum and
        ! the corresponding local array of Cos_maxes/Cos_maxes_index/Cos_sums simultaneously.
        j = indxg2l(diag, block_col, 0, 0, n_procs_col)
        Cos_maxes(j) = max_in_col
        Cos_maxes_index(j) = max_in_col_index
        Cos_sums(j) = sum_in_col
      end if
    end do

    ! Aggregates in the direction of row
    call pdamax(n - 1, cos_max, cos_max_index1, Cos_maxes, 1, 1, desc_Cos_maxes, 1)
    cos_max = abs(cos_max)
    call pielget('R', ' ', cos_max_index2, Cos_maxes_index, 1, cos_max_index1, desc_Cos_maxes_index)
    cos_max_index1 = cos_max_index1 + index1 - 1
    cos_max_index2 = cos_max_index2 + index1 - 1
    call pdasum(n - 1, cos_ave, Cos_sums, 1, 1, desc_Cos_sums, 1)
    cos_ave = cos_ave / dble(n * (n - 1) / 2)
  end subroutine eval_orthogonality_blacs


  subroutine eval_orthogonality(arg, eigenpairs, &
       cos_ave, cos_max, cos_max_index1, cos_max_index2)
    type(argument), intent(in) :: arg
    type(eigenpairs_types_union), intent(in) :: eigenpairs
    double precision, intent(out) :: cos_ave, cos_max
    integer, intent(out) :: cos_max_index1, cos_max_index2

    if (eigenpairs%type_number == 1) then
      call eval_orthogonality_local(arg%ortho_check_index_start, &
           arg%ortho_check_index_end, eigenpairs%local, &
           cos_ave, cos_max, cos_max_index1, cos_max_index2)
    else if (eigenpairs%type_number == 2) then
      call eval_orthogonality_blacs(arg%ortho_check_index_start, &
           arg%ortho_check_index_end, eigenpairs%blacs, &
           cos_ave, cos_max, cos_max_index1, cos_max_index2)
    else
      print '("[Warning] eval_orthogonality: orthogonality evaluator for output of this type is not implemeted yet")'
    end if
  end subroutine eval_orthogonality
end module verifier
