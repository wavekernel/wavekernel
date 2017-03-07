module wk_linear_algebra_m
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_event_logger_m
  use wk_processes_m
  use wk_global_variables_m
  use wk_util_m
  implicit none

  private
  external :: zdotc, numroc, indxg2p, blacs_pnum
  complex(kind(0d0)) :: zdotc
  integer :: numroc, indxg2p, blacs_pnum

  public :: matvec_sd_z, matvec_dd_z, matvec_time_evolution_by_matrix_replace, &
       matvec_nearest_orthonormal_matrix, scale_columns_to_stochastic_matrix, &
       solve_gevp, normalize_vector, get_A_inner_product, get_A_sparse_inner_product, &
       get_ipratio, set_inv_sqrt, reduce_hamiltonian, &
       get_symmetricity, print_offdiag_norm, cutoff_vector, get_moment_for_each_column, &
       add_diag

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


  subroutine matvec_time_evolution_by_matrix_replace(proc, delta_t, &
       H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
       dv_psi_in, dv_psi_out)
    type(wk_process_t), intent(in) :: proc
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
    call setup_distributed_matrix_real('H0', proc, dim, dim, desc, H0, .true.)
    call setup_distributed_matrix_real('H1', proc, dim, dim, desc, H1, .true.)
    call setup_distributed_matrix_real('Hdiff', proc, dim, dim, desc, Hdiff, .true.)
    call setup_distributed_matrix_real('V', proc, dim, dim, desc, Eigenvectors, .true.)
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
    buf_column_sum(:) = kZero
    buf_column_sum_recv(:) = kZero
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
  subroutine solve_gevp(dim, origin, proc, H, H_desc, S, S_desc, eigenvalues, Y_all, Y_all_desc)
    integer, intent(in) :: dim, origin
    integer, intent(in) :: H_desc(desc_size), S_desc(desc_size)
    integer, intent(in) :: Y_all_desc(desc_size)
    type(wk_process_t), intent(in) :: proc
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

    call setup_distributed_matrix_real('H2', proc, dim, dim, H2_desc, H2, .true.)
    call setup_distributed_matrix_real('L', proc, dim, dim, L_desc, L, .true.)
    call pdgemr2d(dim, dim, H, origin, origin, H_desc, H2, 1, 1, H2_desc, H_desc(context_))
    call pdgemr2d(dim, dim, S, origin, origin, S_desc, L, 1, 1, L_desc, S_desc(context_))

    ! Workspace for pdsyngst
    nb = H_desc(block_col_)
    npg = numroc(dim, nb, 0, 0, proc%n_procs_row)
    lwork_pdsyngst = 3 * npg * nb + nb * nb
    ! Workspace for pdsyevd
    npe = numroc(dim, nb, proc%my_proc_row, 0, proc%n_procs_row)
    nqe = numroc(dim, nb, proc%my_proc_col, 0, proc%n_procs_col)
    trilwmin = 3 * dim + max(nb * (npe + 1), 3 * nb)
    lwork = max(1 + 6 * dim + 2 * npe * nqe, trilwmin) + 2 * dim
    ! Temporary implementation. Under some unknown condition, pdormtr in pdsyevd fails due to small lwork.
    lwork = (lwork + 1000) * 2
    liwork = 7 * dim + 8 * proc%n_procs_col + 2
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
end module wk_linear_algebra_m
