module wk_linear_algebra_m
  use mpi
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_event_logger_m
  use wk_processes_m
  use wk_global_variables_m
  use wk_matrix_io_m
  use wk_util_m
  implicit none

  private
  external :: zdotc, numroc, indxg2p, blacs_pnum
  complex(kind(0d0)) :: zdotc
  integer :: numroc, indxg2p, blacs_pnum

  public :: matvec_sd_z, matvec_dd_z, matvec_bd_z, matmul_sd_d, matmul_sb_d_diagonal, matmul_bsb_d, &
       elementwise_matmul_bb_d_inplace, &
       matvec_time_evolution_by_matrix_replace, &
       matvec_nearest_orthonormal_matrix, scale_columns_to_stochastic_matrix, &
       solve_gevp, normalize_vector, get_A_inner_product, get_A_sparse_inner_product, &
       get_ipratio, set_inv_sqrt, reduce_hamiltonian, &
       get_symmetricity, print_offdiag_norm, cutoff_vector, get_moment_for_each_column, &
       add_diag, sum_columns, sum_columns_block, maxloc_columns, maxloc_columns_block, &
       pdscal_elem_block, pdscal_column_block

contains

  ! Compuet y <- alpha * A * x + beta * y (A: symmetric).
  subroutine matvec_sd_z(trans, A, alpha, dv_x, beta, dv_y)
    character(len=*), intent(in) :: trans
    type(sparse_mat), intent(in) :: A
    complex(kind(0d0)), intent(in) :: alpha, beta, dv_x(A%size)
    complex(kind(0d0)), intent(inout) :: dv_y(A%size)

    integer, parameter :: block_size = 32
    integer :: n, i, j, rank, size, num_blocks, b, bl, n1, n2, ierr
    real(8) :: elem, wtime_start
    complex(kind(0d0)) :: dv_Ax(A%size), dv_Ax_recv(A%size)

    wtime_start = mpi_wtime()

    if (trans(1 : 1) == 'T') then
      call terminate('not implemented', 83)
    end if

    call check_nan_vector('matvec_sd_z input real', dreal(dv_x))
    call check_nan_vector('matvec_sd_z input imag', aimag(dv_x))

    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, size, ierr)
    num_blocks = (numroc(A%num_non_zeros, block_size, rank, 0, size) - 1) / block_size + 1

    dv_Ax(:) = kZero
    do b = 1, num_blocks
      n1 = block_size * (size * (b - 1) + rank) + 1
      n2 = min(n1 + block_size - 1, A%num_non_zeros)
      do n = n1, n2
        elem = A%value(n)
        i = A%suffix(1, n)
        j = A%suffix(2, n)
        dv_Ax(i) = dv_Ax(i) + elem * dv_x(j)
        if (i /= j) then
          dv_Ax(j) = dv_Ax(j) + elem * dv_x(i)
        end if
      end do
    end do

    dv_Ax_recv(:) = kZero
    call mpi_allreduce(dv_Ax, dv_Ax_recv, A%size, mpi_double_complex, mpi_sum, mpi_comm_world, ierr)
    if (beta == kZero) then
      dv_y(:) = kZero
    end if
    dv_y(:) = beta * dv_y(:) + dv_Ax_recv(:)

    call check_nan_vector('matvec_sd_z output real', dreal(dv_y))
    call check_nan_vector('matvec_sd_z output imag', aimag(dv_y))
    call add_event('matvec_sd_z', mpi_wtime() - wtime_start)
  end subroutine matvec_sd_z


  ! Compuet y <- alpha * A * x + beta * y.
  subroutine matvec_dd_z(trans, A, A_desc, alpha, dv_x, beta, dv_y)
    character(len=*), intent(in) :: trans
    real(8), intent(in) :: A(:, :)
    integer, intent(in) :: A_desc(desc_size)
    complex(kind(0d0)), intent(in) :: alpha, beta, dv_x(:)
    complex(kind(0d0)), intent(inout) :: dv_y(:)

    integer :: i, j, m, n, nprow, npcol, myrow, mycol, li, lj, pi, pj, ierr
    integer :: mb, nb, brow, bcol, bi, bj, ml, nl, li1, li2, lj1, lj2
    real(8) :: elem, wtime_start
    complex(kind(0d0)), allocatable :: dv_Ax(:), dv_Ax_recv(:)

    wtime_start = mpi_wtime()
    call check_nan_vector('matvec_dd_z input real', dreal(dv_x))
    call check_nan_vector('matvec_dd_z input imag', aimag(dv_x))
    call check_nan_matrix('matvec_dd_z input matrix', A)

    call blacs_gridinfo(A_desc(context_), nprow, npcol, myrow, mycol)
    if (trans(1 : 1) == 'T') then
      m = A_desc(cols_)
      n = A_desc(rows_)
    else
      m = A_desc(rows_)
      n = A_desc(cols_)
    end if
    allocate(dv_Ax(m), dv_Ax_recv(m))
    dv_Ax(:) = kZero
    dv_Ax_recv(:) = kZero
    if (trans(1 : 1) == 'T') then
      mb = A_desc(block_col_)
      nb = A_desc(block_row_)
      brow = (numroc(m, mb, mycol, A_desc(csrc_), npcol) - 1) / mb + 1
      bcol = (numroc(n, nb, myrow, A_desc(rsrc_), nprow) - 1) / nb + 1
      do bi = 1, brow
        do bj = 1, bcol
          i = mb * (npcol * (bi - 1) + mycol) + 1
          j = nb * (nprow * (bj - 1) + myrow) + 1
          ml = min(mb, m - i + 1)
          nl = min(nb, n - j + 1)
          li1 = mb * (bi - 1) + 1
          li2 = li1 + ml - 1
          lj1 = nb * (bj - 1) + 1
          lj2 = lj1 + nl - 1
          do li = li1, li2
            do lj = lj1, lj2
              dv_Ax(i + li - li1) = dv_Ax(i + li - li1) + A(lj, li) * dv_x(j + lj - lj1)
            end do
          end do
        end do
      end do
    else
      mb = A_desc(block_row_)
      nb = A_desc(block_col_)
      brow = (numroc(m, mb, myrow, A_desc(rsrc_), nprow) - 1) / mb + 1
      bcol = (numroc(n, nb, mycol, A_desc(csrc_), npcol) - 1) / nb + 1
      do bi = 1, brow
        do bj = 1, bcol
          i = mb * (nprow * (bi - 1) + myrow) + 1
          j = nb * (npcol * (bj - 1) + mycol) + 1
          ml = min(mb, m - i + 1)
          nl = min(nb, n - j + 1)
          li1 = mb * (bi - 1) + 1
          li2 = li1 + ml - 1
          lj1 = nb * (bj - 1) + 1
          lj2 = lj1 + nl - 1
          do li = li1, li2
            do lj = lj1, lj2
              dv_Ax(i + li - li1) = dv_Ax(i + li - li1) + A(li, lj) * dv_x(j + lj - lj1)
            end do
          end do
        end do
      end do
    end if
    call mpi_allreduce(dv_Ax, dv_Ax_recv, m, mpi_double_complex, mpi_sum, mpi_comm_world, ierr)

    if (beta == kZero) then
      dv_y(:) = kZero
    end if
    dv_y(:) = alpha * dv_Ax_recv(:) + beta * dv_y(:)

    call check_nan_vector('matvec_dd_z output real', dreal(dv_y))
    call check_nan_vector('matvec_dd_z output imag', aimag(dv_y))
    call add_event('matvec_dd_z', mpi_wtime() - wtime_start)
  end subroutine matvec_dd_z


  subroutine matvec_dd_z_comm(trans, A, A_desc, alpha, dv_x, beta, dv_y)
    character(len=*), intent(in) :: trans
    real(8), intent(in) :: A(:, :)
    integer, intent(in) :: A_desc(desc_size)
    complex(kind(0d0)), intent(in) :: alpha, beta, dv_x(:)
    complex(kind(0d0)), intent(inout) :: dv_y(:)

    integer :: i, j, m, n, nprow, npcol, myrow, mycol, li, lj, pi, pj, comm, ierr
    integer :: mb, nb, brow, bcol, bi, bj, ml, nl, li1, li2, lj1, lj2
    real(8) :: elem, wtime_start
    complex(kind(0d0)), allocatable :: dv_Ax(:), dv_Ax_recv(:)

    wtime_start = mpi_wtime()
    call check_nan_vector('matvec_dd_z input real', dreal(dv_x))
    call check_nan_vector('matvec_dd_z input imag', aimag(dv_x))
    call check_nan_matrix('matvec_dd_z input matrix', A)

    call blacs_gridinfo(A_desc(context_), nprow, npcol, myrow, mycol)
    if (trans(1 : 1) == 'T') then
      m = A_desc(cols_)
      n = A_desc(rows_)
    else
      m = A_desc(rows_)
      n = A_desc(cols_)
    end if
    allocate(dv_Ax(m), dv_Ax_recv(m))
    dv_Ax(:) = kZero
    dv_Ax_recv(:) = kZero
    if (trans(1 : 1) == 'T') then
      mb = A_desc(block_col_)
      nb = A_desc(block_row_)
      brow = (numroc(m, mb, mycol, A_desc(csrc_), npcol) - 1) / mb + 1
      bcol = (numroc(n, nb, myrow, A_desc(rsrc_), nprow) - 1) / nb + 1
      do bi = 1, brow
        do bj = 1, bcol
          i = mb * (npcol * (bi - 1) + mycol) + 1
          j = nb * (nprow * (bj - 1) + myrow) + 1
          ml = min(mb, m - i + 1)
          nl = min(nb, n - j + 1)
          li1 = mb * (bi - 1) + 1
          li2 = li1 + ml - 1
          lj1 = nb * (bj - 1) + 1
          lj2 = lj1 + nl - 1
          do li = li1, li2
            do lj = lj1, lj2
              dv_Ax(i + li - li1) = dv_Ax(i + li - li1) + A(lj, li) * dv_x(j + lj - lj1)
            end do
          end do
        end do
      end do
    else
      mb = A_desc(block_row_)
      nb = A_desc(block_col_)
      brow = (numroc(m, mb, myrow, A_desc(rsrc_), nprow) - 1) / mb + 1
      bcol = (numroc(n, nb, mycol, A_desc(csrc_), npcol) - 1) / nb + 1
      do bi = 1, brow
        do bj = 1, bcol
          i = mb * (nprow * (bi - 1) + myrow) + 1
          j = nb * (npcol * (bj - 1) + mycol) + 1
          ml = min(mb, m - i + 1)
          nl = min(nb, n - j + 1)
          li1 = mb * (bi - 1) + 1
          li2 = li1 + ml - 1
          lj1 = nb * (bj - 1) + 1
          lj2 = lj1 + nl - 1
          do li = li1, li2
            do lj = lj1, lj2
              dv_Ax(i + li - li1) = dv_Ax(i + li - li1) + A(li, lj) * dv_x(j + lj - lj1)
            end do
          end do
        end do
      end do
    end if
    call blacs_get(A_desc(context_), 10, comm)
    call mpi_allreduce(dv_Ax, dv_Ax_recv, m, mpi_double_complex, mpi_sum, comm, ierr)

    if (beta == kZero) then
      dv_y(:) = kZero
    end if
    dv_y(:) = alpha * dv_Ax_recv(:) + beta * dv_y(:)

    call check_nan_vector('matvec_dd_z output real', dreal(dv_y))
    call check_nan_vector('matvec_dd_z output imag', aimag(dv_y))
    call add_event('matvec_dd_z', mpi_wtime() - wtime_start)
  end subroutine matvec_dd_z_comm


  ! Compuet y <- alpha * op(A) * x + beta * y.
  ! A: block diagonal matrix.
  ! x, y: redundant vector.
  ! op(): if trans starts with 'T', transpose the matrix.
  subroutine matvec_bd_z(trans, A, alpha, dv_x, beta, dv_y)
    character(len=*), intent(in) :: trans
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: A
    complex(kind(0d0)), intent(in) :: alpha, beta, dv_x(:)
    complex(kind(0d0)), intent(inout) :: dv_y(:)

    integer :: l, g, i1, i2, j1, j2, i1_local, i2_local, j1_local, j2_local, my_rank_in_color, ierr
    complex(kind(0d0)), allocatable :: dv_x_work(:), dv_y_work(:)
    complex(kind(0d0)), allocatable :: dv_y_buf(:)
    ! Functions.
    integer :: numroc, indxl2g

    if (trans(1 : 1) == 'T') then
      allocate(dv_y_buf(A%n))
    else
      allocate(dv_y_buf(A%m))
    end if
    dv_y_buf(:) = kZero
    do l = 1, numroc(A%num_blocks, 1, A%my_color_index, 0, A%num_colors)  ! Iteration for local blocks.
      g = indxl2g(l, 1, A%my_color_index, 0, A%num_colors)
      i1 = A%block_to_row(g)
      i2 = A%block_to_row(g + 1) - 1
      j1 = A%block_to_col(g)
      j2 = A%block_to_col(g + 1) - 1
      i1_local = 1
      i2_local = i2 - i1 + 1
      j1_local = 1
      j2_local = j2 - j1 + 1
      if (trans(1 : 1) == 'T') then
        allocate(dv_x_work(i2_local), dv_y_work(j2_local))
        dv_x_work(:) = real(dv_x(i1 : i2))
        call matvec_dd_z_comm('T', A%local_matrices(l)%val, A%descs(:, l), &
             kOne, dv_x_work, kZero, dv_y_work)
        dv_y_buf(j1 : j2) = dv_y_work(1 : j2_local)
      else
        allocate(dv_x_work(j2_local), dv_y_work(i2_local))
        dv_x_work(:) = real(dv_x(j1 : j2))
        call matvec_dd_z_comm('N', A%local_matrices(l)%val, A%descs(:, l), &
             kOne, dv_x_work, kZero, dv_y_work)
        dv_y_buf(i1 : i2) = dv_y_work(1 : i2_local)
      end if
      deallocate(dv_x_work, dv_y_work)
    end do
    call mpi_comm_rank(A%my_comm, my_rank_in_color, ierr)
    if (my_rank_in_color > 0) then  ! Alleviate duplicated sum.
      dv_y_buf(:) = 0d0
    end if
    if (trans(1 : 1) == 'T') then
      call mpi_allreduce(dv_y_buf, dv_y, A%n, mpi_double_complex, mpi_sum, mpi_comm_world, ierr)
    else
      call mpi_allreduce(dv_y_buf, dv_y, A%m, mpi_double_complex, mpi_sum, mpi_comm_world, ierr)
    end if
  end subroutine matvec_bd_z


  ! Compute C <- alpha * A * B + beta * C.
  ! where A: symmetric sparse redundant.
  !       B, C: dense distributed and share descriptor.
  subroutine matmul_sd_d(A, B, B_desc, alpha, C, beta)
    type(sparse_mat), intent(in) :: A
    real(8), intent(in) :: B(:, :), alpha, beta
    integer, intent(in) :: B_desc(desc_size)
    real(8), intent(inout) :: C(:, :)

    integer :: m, n, A_desc(desc_size)
    real(8), allocatable :: A_dist(:, :)

    m = B_desc(rows_)
    n = B_desc(cols_)
    if (m /= A%size) then
      call terminate('matmul_sd_d: invalid matrix size', 1)
    end if
    call setup_distributed_matrix_real('A', m, m, A_desc, A_dist, .true.)
    call distribute_global_sparse_matrix_wk(A, A_desc, A_dist)
    call pdgemm('No', 'No', m, n, m, alpha, &
         A_dist, 1, 1, A_desc, &
         B, 1, 1, B_desc, &
         beta, &
         C, 1, 1, B_desc)
  end subroutine matmul_sd_d


  ! Compute C <- alpha * bd(A) * B + beta * C.
  ! where A: symmetric sparse matrix (redundant in all processes).
  !       B, C: block diagonal matrices with the same block structure.
  !       bd(): select only block diagonal part based on the block structure of B and C.
  subroutine matmul_sb_d_diagonal(A, B, alpha, C, beta)
    type(sparse_mat), intent(in) :: A
    real(8), intent(in) :: alpha, beta
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: B
    type(wk_distributed_block_diagonal_matrices_t), intent(inout) :: C

    integer :: l, g, col, index_offset, k, k1, k2, i, j
    type(sparse_mat) :: A_work
    integer, allocatable :: group_nonzero_indices(:)
    real(8) :: x, elem_B, elem_C
    ! Functions.
    integer :: indxl2g

    call check_eq_distributed_block_diagonal_matrices(B, C)
    ! Prepare sorted A.
    call copy_sparse_matrix(A, A_work)
    allocate(group_nonzero_indices(B%num_blocks * B%num_blocks + 1))
    call sort_sparse_matrix_by_group(B%num_blocks, B%block_to_row, A_work, group_nonzero_indices)
    do l = 1, numroc(B%num_blocks, 1, B%my_color_index, 0, B%num_colors)
      g = indxl2g(l, 1, B%my_color_index, 0, B%num_colors)
      index_offset = B%block_to_row(g) - 1
      ! Get corresponding block diagonal part of A.
      k1 = group_nonzero_indices(B%num_blocks * (g - 1) + g)
      k2 = group_nonzero_indices(B%num_blocks * (g - 1) + g + 1)
      ! C_g = alpha * bd(A)_g * B_g + beta * C_g.
      do k = k1, k2
        i = A_work%suffix(1, k)
        j = A_work%suffix(2, k)
        x = A_work%value(k)
        do col = 1, B%block_to_col(g + 1) - B%block_to_col(g)
          call pdelget('Self', ' ', elem_B, B%local_matrices(l)%val, j - index_offset, col)
          call pdelget('Self', ' ', elem_C, C%local_matrices(l)%val, i - index_offset, col)
          elem_C = alpha * x * elem_B + beta * elem_C
    !      C%local_matrices(l)%val(i - index_offset, col) = C%local_matrices(l)%val(i - index_offset, col) + &
    !           alpha * x * B%local_matrices(l)%val(j - index_offset, col)
        end do
      end do
    end do
  end subroutine matmul_sb_d_diagonal


  integer function get_rg(X, large_block_col, rg_cycle)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: X
    integer, intent(in) :: large_block_col, rg_cycle

    stop 'IMPLEMENT HERE 8'

    !get_rg = mod(X%num_procs * (large_block_col - 1) + X%my_rank + rg_cycle - 1, X%num_blocks) + 1
  end function get_rg


  ! Compute B <- alpha * X^T * A * Y + beta * B.
  ! where A: symmetric sparse matrix (redundant in all processes).
  !       X, Y: block diagonal matrices with the same block structure.
  !       B: distributed matrix in ScaLAPACK.
  subroutine matmul_bsb_d(X, A, Y, alpha, B, B_desc, beta)
    type(sparse_mat), intent(in) :: A
    real(8), intent(in) :: alpha, beta
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: X, Y
    real(8), intent(inout) :: B(:, :)
    integer, intent(in) :: B_desc(desc_size)

    integer :: large_block_col, cg, cg1, cg2, rg_cycle, rg
    integer :: send_edge_rank, recv_edge_rank
    logical :: is_active
    ! Functions.
    integer :: numroc

    stop 'IMPLEMENT HERE 9'

    ! Check block structure identity.
    !if (X%m /= Y%m .or. X%n /= Y%n) then
    !  call terminate('matmul_bsb_d: different matrix size', 1)
    !end if
    !if ((X%num_blocks /= Y%num_blocks) .or. &
    !     any(X%block_to_row(:) /= Y%block_to_row(:)) .or. &
    !     any(X%block_to_col(:) /= Y%block_to_col(:))) then
    !  call terminate('matmul_bsb_d: different block structure', 1)
    !end if
    !if (X%num_procs /= Y%num_procs .or. X%my_rank /= Y%my_rank) then
    !  call terminate('matmul_bsb_d: different process assignment', 1)
    !end if
    !
    !do large_block_col = 1, numroc(X%num_blocks)
    !  ! cg1 and cg2 represent a global column block index range (same among all processes).
    !  ! cg is a global column block index (different among all processes).
    !  cg1 = X%num_procs * (large_block_col - 1) + 1
    !  cg2 = X%num_procs * (large_block_col - 1) + X%num_procs
    !  cg = X%num_procs * (large_block_col - 1) + X%my_rank + 1
    !  is_active = cg1 <= cg .and. cg <= cg2
    !  recv_edge_rank = XXXXXXXXXXXXX
    !  ! In each cycle, X_rg^T * A_(rg, cg) * Y_cg is computed.
    !  ! Special treatment for rg_cycle = 1 (diagonal blocks).
    !  rg_cycle = 1
    !  rg = get_rg(X, large_block_col, rg_cycle)
    !  if (rg /= cg) then
    !    call terminate('matmul_bsb_d: index calculation bug', 1)
    !  end if
    !  stop
    !
    !  do rg_cycle = 2, X%num_blocks
    !    rg = get_rg(X, large_block_col, rg_cycle)
    !    send_edge_rank = XXXXXXXXXXXXX
    !  end do
    !end do
  end subroutine matmul_bsb_d


  ! Compute X <- A *. X,
  ! where '*.': elementwise mulitplication.
  !       A, X: block diagonal matrices with the same block structure.
  subroutine elementwise_matmul_bb_d_inplace(A, X)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: A
    type(wk_distributed_block_diagonal_matrices_t), intent(inout) :: X

    integer :: l
    ! Functions.
    integer :: numroc

    call check_eq_distributed_block_diagonal_matrices(A, X)
    do l = 1, numroc(A%num_blocks, 1, A%my_color_index, 0, A%num_colors)
      X%local_matrices(l)%val(:, :) = A%local_matrices(l)%val(:, :) * X%local_matrices(l)%val(:, :)
    end do
  end subroutine elementwise_matmul_bb_d_inplace


  subroutine matvec_time_evolution_by_matrix_replace(delta_t, &
       H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
       dv_psi_in, dv_psi_out)

    real(8), intent(in) :: delta_t
    type(sparse_mat), intent(in) :: H_sparse, S_sparse, H_sparse_prev, S_sparse_prev
    complex(kind(0d0)), intent(in) :: dv_psi_in(H_sparse%size)
    complex(kind(0d0)), intent(out) :: dv_psi_out(H_sparse%size)

    integer :: dim, desc(desc_size), lwork, liwork, ierr, i
    integer, allocatable :: iwork(:)
    real(8) :: eigenvalues(H_sparse%size), lwork_real
    real(8), allocatable :: H0(:, :), H1(:, :), Hdiff(:, :), Eigenvectors(:, :), work(:)
    complex(kind(0d0)) :: dv_psi_work(H_sparse%size)

    dim = H_sparse%size
    call setup_distributed_matrix_real('H0', dim, dim, desc, H0, .true.)
    call setup_distributed_matrix_real('H1', dim, dim, desc, H1, .true.)
    call setup_distributed_matrix_real('Hdiff', dim, dim, desc, Hdiff, .true.)
    call setup_distributed_matrix_real('V', dim, dim, desc, Eigenvectors, .true.)
    call distribute_global_sparse_matrix_wk(H_sparse_prev, desc, H0)
    call distribute_global_sparse_matrix_wk(H_sparse, desc, H1)
    Hdiff(:, :) = -1d0 * (H1(:, :) - H0(:, :)) * delta_t / 2.0
    call pdsyevd('V', 'L', dim, Hdiff, 1, 1, desc, eigenvalues, Eigenvectors, 1, 1, desc, &
         lwork_real, -1, liwork, 0, ierr)
    lwork = ceiling(lwork_real)
    allocate(work(lwork), iwork(liwork))
    call pdsyevd('V', 'L', dim, Hdiff, 1, 1, desc, eigenvalues, Eigenvectors, 1, 1, desc, &
         work, lwork, iwork, liwork, ierr)
    dv_psi_work(:) = kZero
    dv_psi_out(:) = kZero
    call matvec_dd_z('T', Eigenvectors, desc, kOne, dv_psi_in, kZero, dv_psi_work)
    do i = 1, dim
      dv_psi_work(i) = dv_psi_work(i) * exp(kImagUnit * eigenvalues(i))
    end do
    call matvec_dd_z('N', Eigenvectors, desc, kOne, dv_psi_work, kZero, dv_psi_out)
  end subroutine matvec_time_evolution_by_matrix_replace


  subroutine matvec_nearest_orthonormal_matrix(A, A_desc, dv_x, dv_y, B)
    real(8), intent(in) :: A(:, :)
    integer, intent(in) :: A_desc(desc_size)
    complex(kind(0d0)), intent(in) :: dv_x(:)
    complex(kind(0d0)), intent(out) :: dv_y(:)
    real(8), intent(out) :: B(:, :)

    integer :: dim, lwork, ierr
    real(8) :: A_work(size(A, 1), size(A, 2)), singular_values(A_desc(rows_))
    real(8) :: U(size(A, 1), size(A, 2)), VT(size(A, 1), size(A, 2)), lwork_real
    real(8), allocatable :: work(:)
    complex(kind(0d0)) :: dv_VTx(A_desc(rows_))

    dim = A_desc(rows_)
    A_work(:, :) = A(:, :)
    call pdgesvd('V', 'V', dim, dim, A_work, 1, 1, A_desc, singular_values, &
         U, 1, 1, A_desc, VT, 1, 1, A_desc, lwork_real, -1, ierr)
    lwork = ceiling(lwork_real)
    allocate(work(lwork))
    call pdgesvd('V', 'V', dim, dim, A_work, 1, 1, A_desc, singular_values, &
         U, 1, 1, A_desc, VT, 1, 1, A_desc, work, lwork, ierr)
    dv_VTx(:) = kZero
    dv_y(:) = kZero
    call matvec_dd_z('N', VT, A_desc, kOne, dv_x, kZero, dv_VTx)
    call matvec_dd_z('N', U, A_desc, kOne, dv_VTx, kZero, dv_y)
    call pdgemm('N', 'N', dim, dim, dim, 1d0, &
         U, 1, 1, A_desc, VT, 1, 1, A_desc, 0d0, B, 1, 1, A_desc)
  end subroutine matvec_nearest_orthonormal_matrix


  subroutine scale_columns_to_stochastic_matrix(A_desc, A)
    integer, intent(in) :: A_desc(desc_size)
    real(8), intent(inout) :: A(:, :)

    integer :: dim, i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc, ierr
    real(8) :: buf_column_sum(A_desc(cols_)), buf_column_sum_recv(A_desc(cols_))

    call blacs_gridinfo(A_desc(context_), nprow, npcol, myrow, mycol)
    dim = A_desc(cols_)
    buf_column_sum(:) = 0d0
    buf_column_sum_recv(:) = 0d0
    do i = 1, dim
      do j = 1, dim
        call infog2l(i, j, A_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
        if (myrow == rsrc .and. mycol == csrc) then
          buf_column_sum(j) = buf_column_sum(j) + A(i_local, j_local)
        end if
      end do
    end do
    call mpi_allreduce(buf_column_sum, buf_column_sum_recv, dim, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    do j = 1, dim
      if (abs(buf_column_sum_recv(j)) >= 1d-10) then
        call pdscal(dim, 1d0 / buf_column_sum_recv(j), A, 1, j, A_desc, 1)
      end if
    end do
    call check_nan_matrix('scale_columns_to_stochastic_matrix', A)
  end subroutine scale_columns_to_stochastic_matrix


  ! 最初に一般化固有値問題を解く部分.
  ! Complexity: O(m^3).
  subroutine solve_gevp(dim, origin, H, H_desc, S, S_desc, eigenvalues, Y_all, Y_all_desc)
    integer, intent(in) :: dim, origin
    integer, intent(in) :: H_desc(desc_size), S_desc(desc_size)
    integer, intent(in) :: Y_all_desc(desc_size)
    real(8), intent(in) :: H(:, :), S(:, :)
    real(8), intent(out) :: eigenvalues(dim)
    real(8), intent(out) :: Y_all(:, :)

    integer :: lwork_pdsyngst, lwork, lrwork, liwork, info, trilwmin
    integer :: nb, npg, npe, nqe
    integer :: H2_desc(desc_size), L_desc(desc_size)
    integer, allocatable :: iwork(:)
    real(8) :: scale
    !complex(kind(0d0)) :: elem
    real(8), allocatable :: H2(:, :), L(:, :), work_pdsyngst(:), work(:)
    !real(8) :: work_pzlaprnt(1000)

    if (check_master()) then
       write (0, '(A, F16.6, A, I0, A, I0)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] solve_gevp: workspace preparation start. dim, origin: ', dim, ', ', origin
    end if

    call setup_distributed_matrix_with_context_real('H2', H_desc(context_), dim, dim, H2_desc, H2, .true.)
    call setup_distributed_matrix_with_context_real('L', H_desc(context_), dim, dim, L_desc, L, .true.)
    call pdgemr2d(dim, dim, H, origin, origin, H_desc, H2, 1, 1, H2_desc, H_desc(context_))
    call pdgemr2d(dim, dim, S, origin, origin, S_desc, L, 1, 1, L_desc, S_desc(context_))

    ! Workspace for pdsyngst
    nb = H_desc(block_col_)
    npg = numroc(dim, nb, 0, 0, g_n_procs_row)
    lwork_pdsyngst = 3 * npg * nb + nb * nb
    ! Workspace for pdsyevd
    npe = numroc(dim, nb, g_my_proc_row, 0, g_n_procs_row)
    nqe = numroc(dim, nb, g_my_proc_col, 0, g_n_procs_col)
    trilwmin = 3 * dim + max(nb * (npe + 1), 3 * nb)
    lwork = max(1 + 6 * dim + 2 * npe * nqe, trilwmin) + 2 * dim
    ! Temporary implementation. Under some unknown condition, pdormtr in pdsyevd fails due to small lwork.
    lwork = (lwork + 1000) * 2
    liwork = 7 * dim + 8 * g_n_procs_col + 2
    allocate(work_pdsyngst(lwork_pdsyngst), work(lwork), iwork(liwork), stat=info)
    if (info /= 0) then
      call terminate('allocation failed in solve_gevp', info)
    end if

    ! find L' s.t. S = L' L'^\dagger, L <- L'.
    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] solve_gevp: Cholesky factorization start'
    end if
    call pdpotrf('Lower', dim, L, 1, 1, L_desc, info)
    if (info /= 0) then
      call terminate('pdpotrf failed', info)
    end if

    ! H' = L^{-1} H L^{-\dagger}, H <- H'
    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] solve_gevp: matrix reduction start'
    end if
    call pdsyngst(1, 'Lower', dim, &
         H2, 1, 1, H2_desc, &
         L, 1, 1, L_desc, &
         scale, work_pdsyngst, lwork_pdsyngst, info)
    if (info /= 0) then
      call terminate('pdsyngst failed', info)
    end if

    ! find (evs, Y) s.t. H' Y = Y diag(evs),
    ! eigenvalues <- evs, H <- Y
    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] solve_gevp: solving standard eigenvalue problem start'
    end if
    call pdsyevd('Vectors', 'Lower', dim, &
         H2, 1, 1, H2_desc, &
         eigenvalues, Y_all, 1, 1, Y_all_desc, &
         work, lwork, iwork, liwork, info)
    if (info /= 0) then
      call terminate('pdsyevd failed', info)
    end if

    ! Y' = L^{-\dagger} Y, H <- Y'
    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] solve_gevp: inverse transformation of eigenvectors start'
    end if
    call pdtrtrs('Lower', 'Trans', 'No', dim, dim, &
         L, 1, 1, L_desc, &
         Y_all, 1, 1, Y_all_desc, info)
    if (info /= 0) then
      call terminate('pdtrtrs failed', info)
    end if

    deallocate(work_pdsyngst, work, iwork, H2, L)
  end subroutine solve_gevp


  subroutine normalize_vector(dim, xs)
    integer, intent(in) :: dim
    complex(kind(0d0)), intent(inout) :: xs(:)
    real(8) :: norm
    real(8) :: dznrm2  ! Function.

    norm = dznrm2(dim, xs, 1)
    xs(1 : dim) = xs(1 : dim) / norm
  end subroutine normalize_vector


  ! Complexity: O(dim^2).
  ! z <- x^\dagger A y
  subroutine get_A_inner_product(dim, A, A_desc, dv_x, dv_y, z)
    integer, intent(in) :: dim, A_desc(desc_size)
    real(8), intent(in) :: A(:, :)
    complex(kind(0d0)), intent(in) :: dv_x(dim), dv_y(dim)
    complex(kind(0d0)), intent(out) :: z

    complex(kind(0d0)) :: dv_work(dim)

    dv_work(:) = kZero
    call matvec_dd_z('No', A, A_desc, kOne, dv_y, kZero, dv_work)
    z = dot_product(dv_x, dv_work)
  end subroutine get_A_inner_product


  ! z <- x^\dagger A y
  subroutine get_A_sparse_inner_product(dim, A_sparse, dv_x, dv_y, z)
    integer, intent(in) :: dim
    type(sparse_mat), intent(in) :: A_sparse
    complex(kind(0d0)), intent(in) :: dv_x(:), dv_y(:)
    complex(kind(0d0)), intent(out) :: z

    complex(kind(0d0)) :: dv_work(dim)

    dv_work(:) = kZero
    call matvec_sd_z('No', A_sparse, kOne, dv_y, kZero, dv_work)
    z = dot_product(dv_x, dv_work)
  end subroutine get_A_sparse_inner_product


  double precision function get_ipratio(dim, psi)
    integer, intent(in) :: dim
    complex(kind(0d0)), intent(in) :: psi(dim)

    integer :: i
    double precision :: sum_2nd_power, sum_4th_power

    sum_2nd_power = 0d0
    sum_4th_power = 0d0
    do i = 1, dim
      sum_2nd_power = sum_2nd_power + abs(psi(i)) ** 2d0
      sum_4th_power = sum_4th_power + abs(psi(i)) ** 4d0
    end do
    get_ipratio = sum_4th_power / (sum_2nd_power ** 2d0)
  end function get_ipratio


  ! S^{-1/2} ~ I - (1 / 2) (S - I) + (3 / 8) (S - I)^2. S is used as workspace after computation.
  subroutine set_inv_sqrt(S_desc, S, S_inv_sqrt_desc, S_inv_sqrt)
    integer, intent(in) :: S_desc(desc_size), S_inv_sqrt_desc(desc_size)
    complex(kind(0d0)), intent(inout) :: S(:, :)
    complex(kind(0d0)), intent(out) :: S_inv_sqrt(:, :)

    integer :: row, col, i, k, dim, nprow, npcol, myrow, mycol, iacol, nq0
    integer :: approx_order = 3
    real(8) :: norm_1, norm_f, error_1, error_f, expand_coef, k_
    real(8), allocatable :: work_pzlange(:)
    complex(kind(0d0)), allocatable :: work(:, :), work_acc1(:, :), work_acc2(:, :)
    complex(kind(0d0)) :: elem
    ! Functions.
    integer :: indxg2p, numroc
    real(8) :: pzlange

    dim = S_desc(rows_)
    row = size(S, 1)
    col = size(S, 2)

    call blacs_gridinfo(S_desc(context_), nprow, npcol, myrow, mycol)
    iacol = indxg2p(1, S_desc(block_col_), mycol, S_desc(csrc_), npcol)
    nq0 = numroc(S_desc(cols_), S_desc(block_col_), mycol, iacol, npcol)

    allocate(work(row, col), work_acc1(row, col), work_acc2(row, col), work_pzlange(nq0))

    ! work <- c (S - I)
    work(:, :) = S(:, :)
    do i = 1, dim
      call pzelget('Self', ' ', elem, work, i, i, S_desc)
      elem = elem - kOne
      call pzelset(work, i, i, S_desc, elem)
    end do

    ! 0-th order
    S_inv_sqrt(:, :) = kZero
    do i = 1, dim
      call pzelset(S_inv_sqrt, i, i, S_inv_sqrt_desc, kOne)
    end do
    ! k-th order
    expand_coef = 1d0
    do k = 1, approx_order
      k_ = real(k, kind(0d0))
      expand_coef = expand_coef * (0.5d0 / k_ - 1d0)
      ! k: odd  => work_acc1 <- (c S - c I) ** k
      ! k: even => work_acc2 <- (c S - c I) ** k
      if (k == 1) then
        work_acc1(:, :) = work(:, :)
      else if (mod(k, 2) == 0) then
        call pzhemm('Left', 'Upper', dim, dim, &
             kOne, work, 1, 1, S_desc, &
             work_acc1, 1, 1, S_desc, &
             kZero, work_acc2, 1, 1, S_desc)
      else
        call pzhemm('Left', 'Upper', dim, dim, &
             kOne, work, 1, 1, S_desc, &
             work_acc2, 1, 1, S_desc, &
             kZero, work_acc1, 1, 1, S_desc)
      end if
      if (mod(k, 2) == 0) then
        norm_1 = pzlange('1', dim, dim, work_acc2, 1, 1, S_desc, work_pzlange)
        norm_f = pzlange('Frobenius', dim, dim, work_acc2, 1, 1, S_desc, 0)
        call pzgeadd('No', dim, dim, &
             cmplx(expand_coef, 0d0, kind(0d0)), work_acc2, 1, 1, S_desc, &
             kOne, S_inv_sqrt, 1, 1, S_inv_sqrt_desc)
      else
        norm_1 = pzlange('1', dim, dim, work_acc1, 1, 1, S_desc, work_pzlange)
        norm_f = pzlange('Frobenius', dim, dim, work_acc1, 1, 1, S_desc, 0)
        call pzgeadd('No', dim, dim, &
             cmplx(expand_coef, 0d0, kind(0d0)), work_acc1, 1, 1, S_desc, &
             kOne, S_inv_sqrt, 1, 1, S_inv_sqrt_desc)
      end if
      if (check_master()) then
        write (0, '(A, F16.6, A, I0, A, F16.6, A, F16.6)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
             '] set_inv_sqrt: k is ', k, ', ||S - I||^k_1 is', norm_1, ', ||S - I||^k_F is', norm_f
      end if
    end do

    ! Error checking hereafter.
    norm_1 = pzlange('1', dim, dim, S, 1, 1, S_desc, work_pzlange)
    norm_f = pzlange('Frobenius', dim, dim, S, 1, 1, S_desc, 0)
    ! work <- S^{-1/2} S
    call pzhemm('Left', 'Upper', dim, dim, &
         kOne, S_inv_sqrt, 1, 1, S_inv_sqrt_desc, &
         S, 1, 1, S_desc, &
         kZero, work, 1, 1, S_desc)
    ! S <- S^{-1/2} (S^{-1/2} S). S (LHS) is used as workspace.
    call pzhemm('Left', 'Upper', dim, dim, &
         kOne, S_inv_sqrt, 1, 1, S_inv_sqrt_desc, &
         work, 1, 1, S_desc, &
         kZero, S, 1, 1, S_desc)
    ! S <- S^{-1/2} (S^{-1/2} S) - I
    do i = 1, dim
      call pzelget('Self', ' ', elem, S, i, i, S_desc)
      elem = elem - kOne
      call pzelset(S, i, i, S_desc, elem)
    end do
    ! error <- ||S^{-1/2} (S^{-1/2} S) - I|| / ||S||
    error_1 = pzlange('1', dim, dim, S, 1, 1, S_desc, work_pzlange) / norm_1
    error_f = pzlange('Frobenius', dim, dim, S, 1, 1, S_desc, 0) / norm_f
    if (check_master()) then
      write (0, '(A, F16.6, A, F16.6, A, F16.6)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] set_inv_sqrt: relative error in approximating inverse square root of S (1-norm, Frobenius norm) is', &
           error_1, ',', error_f
    end if
  end subroutine set_inv_sqrt


  ! S is used as workspace temporarily.
  subroutine reduce_hamiltonian(S_inv_sqrt_desc, S_inv_sqrt, H_desc, H, S_desc, S)
    integer, intent(in) :: S_inv_sqrt_desc(desc_size), H_desc(desc_size), S_desc(desc_size)
    complex(kind(0d0)), intent(in) :: S_inv_sqrt(:, :)
    complex(kind(0d0)), intent(inout) :: H(:, :)
    complex(kind(0d0)), intent(out) :: S(:, :)
    integer :: i, dim

    dim = H_desc(rows_)
    ! S <- S^{-1/2} H
    call pzhemm('Left', 'Upper', dim, dim, &
         kOne, S_inv_sqrt, 1, 1, S_inv_sqrt_desc, &
         H, 1, 1, H_desc, &
         kZero, S, 1, 1, S_desc)
    ! H <- (S^{-1/2} H) S^{-1/2}
    call pzhemm('Right', 'Upper', dim, dim, &
         kOne, S_inv_sqrt, 1, 1, S_inv_sqrt_desc, &
         S, 1, 1, S_desc, &
         kZero, H, 1, 1, H_desc)
    ! Set S to I.
    S(:, :) = 0d0
    do i = 1, dim
      call pzelset(S, i, i, S_desc, kOne)
    end do
  end subroutine reduce_hamiltonian


  ! Calculate ||A - A^T||_F.
  subroutine get_symmetricity(A_desc, A, symmetricity)
    integer, intent(in) :: A_desc(desc_size)
    real(8), intent(in) :: A(:, :)
    real(8), intent(out) :: symmetricity
    integer :: m, n, i, j
    real(8) :: elem1, elem2
    real(8), allocatable :: B(:, :)
    real(8) :: pdlange

    m = A_desc(rows_)
    n = A_desc(cols_)
    allocate(B(size(A, 1), size(A, 2)))
    B(:, :) = 0d0
    do i = 1, m
      do j = i + 1, n
        call pdelget('All', ' ', elem1, A, i, j, A_desc)
        call pdelget('All', ' ', elem2, A, j, i, A_desc)
        call pdelset(B, i, j, A_desc, elem1 - elem2)
      end do
    end do
    symmetricity = pdlange('Frobenius', m, n, B, 1, 1, A_desc, 0) * sqrt(2d0)
    deallocate(B)
  end subroutine get_symmetricity


  subroutine add_diag(X_desc, X, diag)
    integer, intent(in) :: X_desc(desc_size)
    complex(kind(0d0)) :: X(:, :)
    complex(kind(0d0)) :: diag, elem
    integer :: m, n, i, j

    m = X_desc(rows_)
    n = X_desc(cols_)
    do i = 1, min(m, n)
      call pzelget('Self', ' ', elem, X, i, i, X_desc)
      call pzelset(X, i, i, X_desc, elem + diag)
    end do
  end subroutine add_diag


  subroutine print_offdiag_norm(name, X, X_desc)
    character(len=*), intent(in) :: name
    complex(kind(0d0)), intent(in) :: X(:, :)
    integer, intent(in) :: X_desc(desc_size)

    complex(kind(0d0)), allocatable :: X_offdiag(:, :)
    real(8), allocatable :: pzlange_work(:)
    integer :: dim, nprow, npcol, myrow, mycol, iarow, iacol, mp0, nq0, i
    real(8) :: X_offdiag_norm_1, X_offdiag_norm_inf, X_offdiag_norm_frobenius
    integer :: indxg2p, numroc
    real(8) :: pzlange

    dim = min(X_desc(rows_), X_desc(cols_))
    call blacs_gridinfo(X_desc(context_), nprow, npcol, myrow, mycol)
    iarow = indxg2p(1, X_desc(block_row_), myrow, X_desc(rsrc_), nprow)
    iacol = indxg2p(1, X_desc(block_col_), mycol, X_desc(csrc_), npcol)
    mp0 = numroc(X_desc(rows_), X_desc(block_row_), myrow, iarow, nprow)
    nq0 = numroc(X_desc(cols_), X_desc(block_col_), mycol, iacol, npcol)

    allocate(X_offdiag(size(X, 1), size(X, 2)), pzlange_work(max(mp0, nq0)))
    X_offdiag(:, :) = X(:, :)
    do i = 1, dim
      call pzelset(X_offdiag, i, i, X_desc, kZero)
    end do
    X_offdiag_norm_1 = pzlange('1', X_desc(rows_), X_desc(cols_), X_offdiag, 1, 1, X_desc, pzlange_work)
    X_offdiag_norm_inf = pzlange('I', X_desc(rows_), X_desc(cols_), X_offdiag, 1, 1, X_desc, pzlange_work)
    X_offdiag_norm_frobenius = pzlange('F', X_desc(rows_), X_desc(cols_), X_offdiag, 1, 1, X_desc, 0)
    if (check_master()) then
      write (0, *) '[1-norm of offdiagonal part of ', name, '] ', X_offdiag_norm_1
      write (0, *) '[inf-norm of offdiagonal part of ', name, '] ', X_offdiag_norm_inf
      write (0, *) '[Frobenius norm of offdiagonal part of ', name, '] ', X_offdiag_norm_frobenius
    end if
    deallocate(X_offdiag, pzlange_work)
  end subroutine print_offdiag_norm


  subroutine cutoff_vector(dim, cutoff, vec_in, vec_out)
    integer, intent(in) :: dim
    real(8), intent(in) :: cutoff
    complex(kind(0d0)), intent(in) :: vec_in(dim)
    complex(kind(0d0)), intent(out) :: vec_out(dim)

    integer :: i, vec_indices(dim), cutoff_index
    real(8) :: norm2, acc, vec_abs2(dim)

    do i = 1, dim
      vec_abs2(i) = abs(vec_in(i)) ** 2d0
    end do
    norm2 = sum(vec_abs2(:))
    call comb_sort(dim, vec_abs2, vec_indices)
    acc = 0d0
    do i = 1, dim
      cutoff_index = i
      acc = acc + vec_abs2(i)
      if (acc / norm2 >= 1d0 - cutoff) then
        exit
      end if
    end do
    if (check_master()) then
      write (0, *) '[vector cutoff index] ', cutoff_index, '/', dim
      write (0, *) '[vector cutoff accumulated amplitude] ', acc / norm2
    end if

    vec_out(:) = vec_in(:)
    do i = cutoff_index + 1, dim
      vec_out(vec_indices(i)) = kZero
    end do
  end subroutine cutoff_vector


  subroutine get_moment_for_each_column(fs, X_desc, X, mean, dev)
    integer, intent(in) :: X_desc(desc_size)
    real(8), intent(in) :: fs(X_desc(rows_)), X(:, :)
    real(8), intent(out) :: mean(X_desc(cols_)), dev(X_desc(cols_))

    integer :: i, j, ierr
    real(8) :: column_buf(X_desc(rows_)), column(X_desc(rows_)), normalizer

    mean(:) = 0d0
    dev(:) = 0d0
    do j = 1, X_desc(cols_)
      column_buf(:) = 0d0
      column(:) = 0d0
      do i = 1, X_desc(rows_)
        call pdelget('Self', ' ', column_buf(i), X, i, j, X_desc)
        column_buf(i) = column_buf(i) ** 2d0
      end do
      call mpi_allreduce(column_buf, column, X_desc(rows_), mpi_real8, mpi_sum, mpi_comm_world, ierr)
      normalizer = sum(column)
      !print *, 'ZZZZZZnormal', normalizer
      column(:) = column(:) / normalizer
      do i = 1, X_desc(rows_)
        mean(j) = mean(j) + fs(i) * column(i)
      end do
      do i = 1, X_desc(rows_)
        dev(j) = dev(j) + (fs(i) - mean(j)) ** 2d0 * column(i)
      end do
      dev(j) = dsqrt(dev(j))
    end do

    if (check_master()) then
      do j = 1, X_desc(cols_)
        write (0, '(A, F16.6, A, I0, 2E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
             '] get_moment_for_each_column() : j, mean, dev : ', j, mean(j), dev(j)
      end do
    end if
  end subroutine get_moment_for_each_column


  ! sums(j) = \sum_i X(i, j)
  subroutine sum_columns(X, X_desc, sums)
    real(8), intent(in) :: X(:, :)
    integer, intent(in) :: X_desc(desc_size)
    real(8), intent(out) :: sums(X_desc(cols_))

    integer :: n, n_procs_row, n_procs_col, my_proc_row, my_proc_col, cols_local, c_l, c_g, ierr
    real(8) :: sums_buf(X_desc(cols_))
    ! Functions.
    integer :: numroc, indxl2g

    n = X_desc(cols_)
    call blacs_gridinfo(X_desc(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    sums_buf(:) = 0d0
    cols_local = numroc(n, X_desc(block_col_), my_proc_col, X_desc(csrc_), n_procs_col)
    do c_l = 1, cols_local
      c_g = indxl2g(c_l, X_desc(block_col_), my_proc_col, X_desc(csrc_), n_procs_col)
      sums_buf(c_g) = sums_buf(c_g) + sum(X(:, c_l))
    end do
    call mpi_allreduce(sums_buf, sums, n, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    if (ierr /= 0) then
      call terminate('sum_columns: mpi_allreduce failed', 1)
    end if
  end subroutine sum_columns


  ! sums(j) = \sum_i X(i, j) for block diagonal matrix X.
  subroutine sum_columns_block(X, sums)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: X
    real(8) :: sums(X%n)

    integer :: l, g, c1, c2, my_rank_in_color, ierr
    real(8) :: sums_buf(X%n)
    ! Functions.
    integer :: numroc, indxl2g

    sums_buf(:) = 0d0
    do l = 1, numroc(X%num_blocks, 1, X%my_color_index, 0, X%num_colors)
      g = indxl2g(l, 1, X%my_color_index, 0, X%num_colors)
      c1 = X%block_to_col(g)
      c2 = X%block_to_col(g + 1) - 1
      call sum_columns(X%local_matrices(l)%val, X%descs(:, l), sums_buf(c1 : c2))
    end do
    call mpi_comm_rank(X%my_comm, my_rank_in_color, ierr)
    if (my_rank_in_color > 0) then  ! Alleviate duplicated sum.
      sums_buf(:) = 0d0
    end if
    call mpi_allreduce(sums_buf, sums, X%n, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  end subroutine sum_columns_block


  ! maxlocs(j) = \maxloc_i X(i, j)
  subroutine maxloc_columns(X, X_desc, maxlocs)
    real(8), intent(in) :: X(:, :)
    integer, intent(in) :: X_desc(desc_size)
    integer, intent(out) :: maxlocs(X_desc(cols_))

    integer :: n, n_procs_row, n_procs_col, my_proc_row, my_proc_col, cols_local
    integer :: c_l, c_g, r_l, r_g, r_p, ierr1, ierr2
    real(8) :: val
    integer, allocatable :: maxlocs_local(:)  ! rank-1 array for column of X as a local matrix.
    integer :: maxlocs_global(X_desc(cols_))  ! rank-1 array for column of X as a distributed matrix.
    ! For below 2 buffers, dim 1: global processor row, dim 2: global matrix col.
    real(8), allocatable :: values_buf(:, :), values(:, :)
    integer, allocatable :: indices_buf(:, :), indices(:, :)
    ! Functions.
    integer :: numroc, indxl2g, indxg2p

    n = X_desc(cols_)
    call blacs_gridinfo(X_desc(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    cols_local = numroc(n, X_desc(block_col_), my_proc_col, X_desc(csrc_), n_procs_col)
    allocate(maxlocs_local(cols_local))
    allocate(values_buf(n_procs_row, n), indices_buf(n_procs_row, n))
    allocate(values(n_procs_row, n), indices(n_procs_row, n))
    values_buf(:, :) = 0d0
    indices_buf(:, :) = 0
    maxlocs_local = maxloc(X, 1)
    do c_l = 1, cols_local
      r_l = maxlocs_local(c_l)
      c_g = indxl2g(c_l, X_desc(block_col_), my_proc_col, X_desc(csrc_), n_procs_col)
      r_g = indxl2g(r_l, X_desc(block_row_), my_proc_row, X_desc(rsrc_), n_procs_row)
      r_p = indxg2p(r_g, X_desc(block_row_), my_proc_row, X_desc(rsrc_), n_procs_row)
      values_buf(r_p + 1, c_g) = X(r_l, c_l)
      indices_buf(r_p + 1, c_g) = r_g
    end do
    call mpi_allreduce(values_buf, values, n_procs_row * n, mpi_double_precision, mpi_sum, mpi_comm_world, ierr1)
    call mpi_allreduce(indices_buf, indices, n_procs_row * n, mpi_integer, mpi_sum, mpi_comm_world, ierr2)
    if (ierr1 /= 0 .or. ierr2 /= 0) then
      call terminate('maxloc_columns: mpi_allreduce failed', 1)
    end if
    maxlocs_global = maxloc(values_buf, 1)
    do c_g = 1, n
      maxlocs(c_g) = indices(maxlocs_global(c_g), c_g)
    end do
  end subroutine maxloc_columns


  ! maxlocs(j) = \maxloc_i X(i, j) for block diagonal matrix X.
  subroutine maxloc_columns_block(X, maxlocs)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: X
    integer :: maxlocs(X%n)

    integer :: l, g, c1, c2, my_rank_in_color, ierr
    integer :: maxlocs_buf(X%n)
    ! Functions.
    integer :: numroc, indxl2g

    maxlocs_buf(:) = 0
    do l = 1, numroc(X%num_blocks, 1, X%my_color_index, 0, X%num_colors)
      g = indxl2g(l, 1, X%my_color_index, 0, X%num_colors)
      c1 = X%block_to_col(g)
      c2 = X%block_to_col(g + 1) - 1
      call maxloc_columns(X%local_matrices(l)%val, X%descs(:, l), maxlocs_buf(c1 : c2))
    end do
    call mpi_comm_rank(X%my_comm, my_rank_in_color, ierr)
    if (my_rank_in_color > 0) then  ! Alleviate duplicated sum.
      maxlocs_buf(:) = 0
    end if
    call mpi_allreduce(maxlocs_buf, maxlocs, X%n, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  end subroutine maxloc_columns_block


  subroutine pdscal_elem_block(alpha, i, j, X)
    real(8), intent(in) :: alpha
    integer, intent(in) :: i, j
    type(wk_distributed_block_diagonal_matrices_t), intent(inout) :: X

    integer :: g, c, l, i_in_block, j_in_block
    ! Functions.
    integer :: indxg2p, indxg2l

    g = find_range_index(X%block_to_row, i)
    if (g /= find_range_index(X%block_to_col, j)) then  ! (i, j) is out of the blocks.
      return
    end if
    c = indxg2p(g, 1, X%my_color_index, 0, X%num_colors)
    if (X%my_color_index == c) then
      l = indxg2l(g, 1, X%my_color_index, 0, X%num_colors)
      i_in_block = i - X%block_to_row(g) + 1
      j_in_block = j - X%block_to_col(g) + 1
      call pdscal(1, alpha, X%local_matrices(l)%val, i_in_block, j_in_block, X%descs(:, l), 1)
    end if
  end subroutine pdscal_elem_block


  ! X(:, j) *= alpha for block diagonal matrix X.
  subroutine pdscal_column_block(alpha, j, X)
    real(8), intent(in) :: alpha
    integer, intent(in) :: j
    type(wk_distributed_block_diagonal_matrices_t), intent(inout) :: X

    integer :: g, l, c, j_in_block, m, ierr
    ! Functions.
    integer :: indxg2p, indxg2l

    g = find_range_index(X%block_to_col, j)
    c = indxg2p(g, 1, X%my_color_index, 0, X%num_colors)
    if (X%my_color_index == c) then
      l = indxg2l(g, 1, X%my_color_index, 0, X%num_colors)
      j_in_block = j - X%block_to_col(g) + 1
      m = X%block_to_row(g + 1) - X%block_to_row(g)
      call pdscal(m, alpha, X%local_matrices(l)%val, 1, j_in_block, X%descs(:, l), 1)
    end if
  end subroutine pdscal_column_block
end module wk_linear_algebra_m
