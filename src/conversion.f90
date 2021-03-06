module wk_conversion_m
  use mpi
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_event_logger_m
  use wk_linear_algebra_m
  use wk_processes_m
  use wk_global_variables_m
  use wk_state_m
  implicit none

  private
  public :: alpha_to_lcao_coef, &
       change_basis_lcao_to_alpha, &
       change_basis_lcao_to_alpha_group_filter, &
       change_basis_lcao_diag_to_alpha, &
       change_basis_lcao_diag_to_alpha_group_filter, &
       lcao_coef_to_alpha

contains

  ! psi = Y alpha.
  subroutine alpha_to_lcao_coef(basis, dv_alpha, dv_psi)
    type(wk_basis_t), intent(in) :: basis
    complex(kind(0d0)), intent(in) :: dv_alpha(basis%num_basis)
    complex(kind(0d0)), intent(out) :: dv_psi(basis%dim)

    complex(kind(0d0)) :: dv_beta(basis%num_basis)

    if (basis%is_group_filter_mode) then
      dv_beta(:) = kZero
      dv_psi(:) = kZero
      call matvec_dd_z('No', basis%Y_all, basis%Y_all_desc, kOne, dv_alpha, kZero, dv_beta)
      call matvec_bd_z('No', basis%Y_local, kOne, dv_beta, kZero, dv_psi)
    else
      dv_psi(:) = kZero
      call matvec_dd_z('No', basis%Y_filtered, basis%Y_filtered_desc, kOne, dv_alpha, kZero, dv_psi)
    end if
  end subroutine alpha_to_lcao_coef


  subroutine change_basis_lcao_to_alpha(basis, A_lcao_sparse, A_alpha, A_alpha_desc)
    type(wk_basis_t), intent(in) :: basis
    type(sparse_mat), intent(in) :: A_lcao_sparse
    integer, intent(in) :: A_alpha_desc(desc_size)
    real(8), intent(out) :: A_alpha(:, :)

    if (basis%is_group_filter_mode) then
      call change_basis_lcao_to_alpha_group_filter(basis%filter_group_indices, basis%Y_local, &
           A_lcao_sparse, A_alpha, A_alpha_desc)
    else
      call change_basis_lcao_to_alpha_all(basis%Y_filtered, basis%Y_filtered_desc, &
           A_lcao_sparse, A_alpha, A_alpha_desc)
    end if
  end subroutine change_basis_lcao_to_alpha


  subroutine change_basis_lcao_to_alpha_all(Y_filtered, Y_filtered_desc, A_lcao_sparse, A_alpha, A_alpha_desc)
    real(8), intent(in) :: Y_filtered(:, :)
    type(sparse_mat), intent(in) :: A_lcao_sparse
    integer, intent(in) :: Y_filtered_desc(desc_size), A_alpha_desc(desc_size)
    real(8), intent(out) :: A_alpha(:, :)

    integer :: dim, num_filter, AY_desc(desc_size), A_lcao_desc(desc_size)
    real(8), allocatable :: A_lcao(:, :), AY(:, :)

    dim = Y_filtered_desc(rows_)
    num_filter = Y_filtered_desc(cols_)
    call setup_distributed_matrix_real('A', dim, dim, A_lcao_desc, A_lcao, .true.)
    call distribute_global_sparse_matrix_wk(A_lcao_sparse, A_lcao_desc, A_lcao)
    call setup_distributed_matrix_real('AY', dim, num_filter, AY_desc, AY)

    ! A' = A_lcao * Y, AY <- A'
    call pdgemm('No', 'No', dim, num_filter, dim, 1d0, &
         A_lcao, 1, 1, A_lcao_desc, &
         Y_filtered, 1, 1, Y_filtered_desc, &
         0d0, &
         AY, 1, 1, AY_desc)
    ! A'' = Y^\dagger A', A_alpha <- A''
    call pdgemm('Trans', 'No', num_filter, num_filter, dim, 1d0, &
         Y_filtered, 1, 1, Y_filtered_desc, &
         AY, 1, 1, AY_desc, &
         0d0, &
         A_alpha, 1, 1, A_alpha_desc)
    deallocate(A_lcao, AY)
  end subroutine change_basis_lcao_to_alpha_all


  subroutine change_basis_lcao_to_alpha_group_filter(filter_group_indices, Y_local, &
       A_lcao_sparse, A_alpha, A_alpha_desc, to_ignore_diag_block_)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: Y_local
    type(sparse_mat), intent(in) :: A_lcao_sparse
    real(8), intent(out) :: A_alpha(:, :)
    integer, intent(in) :: filter_group_indices(:, :), A_alpha_desc(desc_size)
    logical, intent(in), optional :: to_ignore_diag_block_

    call matmul_bsb_d(Y_local, A_lcao_sparse, Y_local, 1d0, A_alpha, A_alpha_desc, 0d0)
  !
  !  integer :: num_groups, dim, num_filter
  !  integer :: g, p, nprow, npcol, myp, myg, myprow, mypcol, np, prow, pcol
  !  integer :: i_lcao, j_lcao, m_right, n_right, g2, ierr
  !  complex(kind(0d0)), allocatable :: A_lcao_cmplx(:, :)
  !  integer :: A_lcao_desc(desc_size)
  !  real(8), allocatable :: A_lcao(:, :), A_alpha(:, :)  ! Distributed.
  !  real(8), allocatable :: A_lcao_local(:, :), Y_local_right(:, :), Y_local_left(:, :), AY(:, :), YAY(:, :)  ! Local.
  !  integer :: prow_tmp, pcol_tmp, j_lcao_tmp, j_alpha_tmp, m_right_tmp, n_right_tmp, i_alpha, j_alpha, m_left, n_left
  !  integer :: blacs_pnum
  !  logical :: is_group_member, to_ignore_diag_block
  !  real(8) :: wtime_start, wtime_end
  !  !real(8) :: work_pdlaprnt(10000), symmetricity
  !
  !  wtime_start = mpi_wtime()
  !
  !
  !  if (present(to_ignore_diag_block_)) then
  !    to_ignore_diag_block = to_ignore_diag_block_
  !  else
  !    to_ignore_diag_block = .false.  ! Default value.
  !  end if
  !
  !  num_groups = size(filter_group_indices, 2) - 1
  !  dim = filter_group_indices(1, num_groups + 1) - 1
  !  num_filter = filter_group_indices(2, num_groups + 1) - 1
  !  call blacs_gridinfo(proc%context, nprow, npcol, myprow, mypcol)
  !  np = nprow * npcol
  !  myp = blacs_pnum(proc%context, myprow, mypcol)
  !  myg = myp + 1
  !  is_group_member = myg <= num_groups
  !
  !  call setup_distributed_matrix_complex('A', proc, dim, dim, A_lcao_desc, A_lcao_cmplx, .true.)
  !  call distribute_global_sparse_matrix_wk(A_lcao_sparse, A_lcao_desc, A_lcao_cmplx)
  !
  !  allocate(A_lcao(size(A_lcao_cmplx, 1), size(A_lcao_cmplx, 2)), &
  !       A_alpha(size(A_alpha_cmplx, 1), size(A_alpha_cmplx, 2)))
  !  A_lcao(:, :) = real(A_lcao_cmplx(:, :), kind(0d0))
  !  A_alpha(:, :) = 1d20  ! For debug.
  !
  !  if (is_group_member) then
  !    j_lcao = filter_group_indices(1, myg)
  !    j_alpha = filter_group_indices(2, myg)
  !    m_right = filter_group_indices(1, myg + 1) - j_lcao
  !    n_right = filter_group_indices(2, myg + 1) - j_alpha
  !    allocate(Y_local_right(m_right, n_right))
  !    !print *, 'U', myp, ':', filter_group_indices
  !    !print *, 'Y', myp, ':', m_right, n_right, ',', size(Y_local_cmplx, 1), size(Y_local_cmplx, 2)
  !    Y_local_right(:, :) = real(Y_local_cmplx(1)%val(:, :), kind(0d0))
  !  else
  !    allocate(Y_local_right(0, 0))
  !  end if
  !
  !  wtime_end = mpi_wtime()
  !  call add_event('change_basis_lcao_diag_to_alpha_group_filter_real:setup', wtime_end - wtime_start)
  !  wtime_start = wtime_end
  !
  !  !print *, 'P100', myp
  !  !call mpi_barrier(mpi_comm_world, ierr)
  !  do g = 1, num_groups
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    !print *, 'group start', g
  !    p = g - 1
  !    call blacs_pcoord(proc%context, p, prow, pcol)
  !    i_lcao = filter_group_indices(1, g)
  !    i_alpha = filter_group_indices(2, g)
  !    m_left = filter_group_indices(1, g + 1) - i_lcao
  !    n_left = filter_group_indices(2, g + 1) - i_alpha
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    !print *, 'P', myp, ':', i_lcao, i_alpha, m_left, n_left
  !    if (g > 1) then
  !      deallocate(A_lcao_local, Y_local_left, AY, YAY)
  !    end if
  !    if (is_group_member) then
  !      allocate(A_lcao_local(m_left, m_right), Y_local_left(m_left, n_left), AY(m_left, n_right), YAY(n_left, n_right))
  !      !print *, 'GG', myg, g, ',', m_left, m_right
  !    else
  !      allocate(A_lcao_local(0, 0), Y_local_left(m_left, n_left), AY(0, 0), YAY(0, 0))
  !    end if
  !    A_lcao_local(:, :) = 1d30  ! For debug.
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    !print *, 'Palloc', myp
  !    ! Send part of A.
  !
  !    wtime_end = mpi_wtime()
  !    call add_event('change_basis_lcao_diag_to_alpha_group_filter_real:inner_loop1', wtime_end - wtime_start)
  !    wtime_start = wtime_end
  !
  !    do g2 = 1, num_groups
  !      if (to_ignore_diag_block .and. g2 == g) then
  !        ! Without separating this if statement,
  !        ! arguments of gather_matrix_real_part below will differ among processes.
  !        if (myg == g) then
  !          A_lcao_local(:, :) = 0d0
  !        end if
  !      else
  !        j_lcao_tmp = filter_group_indices(1, g2)
  !        j_alpha_tmp = filter_group_indices(2, g2)
  !        m_right_tmp = filter_group_indices(1, g2 + 1) - j_lcao_tmp
  !        n_right_tmp = filter_group_indices(2, g2 + 1) - j_alpha_tmp
  !        call blacs_pcoord(proc%context, g2 - 1, prow_tmp, pcol_tmp)
  !
  !        !print *, 'FF', g, g2, myg, ',', i_lcao, '-', i_lcao + m_left - 1, ',', j_lcao_tmp, '-', j_lcao_tmp + m_right_tmp - 1
  !        call gather_matrix_real_part(A_lcao, A_lcao_desc, i_lcao, j_lcao_tmp, m_left, m_right_tmp, &
  !             prow_tmp, pcol_tmp, A_lcao_local)
  !      end if
  !    end do
  !
  !    wtime_end = mpi_wtime()
  !    call add_event('change_basis_lcao_diag_to_alpha_group_filter_real:inner_loop2', wtime_end - wtime_start)
  !    wtime_start = wtime_end
  !
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    !print *, 'P0', myp
  !    ! Send Y_left.
  !    if (myp == p) then
  !      Y_local_left(:, :) = Y_local_right(:, :)
  !      !print *, 'Px send', myp, m_right, n_right, size(Y_local_left, 1), size(Y_local_left, 2)
  !      call dgebs2d(proc%context, 'All', ' ', m_right, n_right, Y_local_left, m_right)
  !    else
  !      !print *, 'Px recv', myp, m_left, n_left, size(Y_local_left, 1), size(Y_local_left, 2), ',', prow, pcol
  !      call dgebr2d(proc%context, 'All', ' ', m_left, n_left, Y_local_left, m_left, prow, pcol)
  !    end if
  !
  !    wtime_end = mpi_wtime()
  !    call add_event('change_basis_lcao_diag_to_alpha_group_filter_real:inner_loop3', wtime_end - wtime_start)
  !    wtime_start = wtime_end
  !
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    !print *, 'P1', myp
  !    ! Calculate YAY.
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    if (is_group_member) then
  !      call dgemm('No', 'No', m_left, n_right, m_right, &
  !           1d0, A_lcao_local, m_left, Y_local_right, m_right, &
  !           0d0, AY, m_left)
  !    end if
  !    !do g2 = 1, num_groups
  !    !  call mpi_barrier(mpi_comm_world, ierr)
  !    !  if (myp == g2 - 1) then
  !    !    print *, 'AY', myp, ':', AY
  !    !  end if
  !    !end do
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    if (is_group_member) then
  !      call dgemm('Trans', 'No', n_left, n_right, m_left, &
  !           1d0, Y_local_left, m_left, AY, m_left, &
  !           0d0, YAY, n_left)
  !    end if
  !    !do g2 = 1, num_groups
  !    !  call mpi_barrier(mpi_comm_world, ierr)
  !    !  if (myp == g2 - 1) then
  !    !    print *, 'YAY', myp, ':', YAY
  !    !  end if
  !    !end do
  !
  !    !call mpi_barrier(mpi_comm_world, ierr)
  !    !print *, 'P2', myp
  !
  !    wtime_end = mpi_wtime()
  !    call add_event('change_basis_lcao_diag_to_alpha_group_filter_real:inner_loop4', wtime_end - wtime_start)
  !    wtime_start = wtime_end
  !
  !    do g2 = 1, num_groups
  !      j_alpha_tmp = filter_group_indices(2, g2)
  !      n_right_tmp = filter_group_indices(2, g2 + 1) - j_alpha_tmp
  !      call blacs_pcoord(proc%context, g2 - 1, prow_tmp, pcol_tmp)
  !      !call mpi_barrier(mpi_comm_world, ierr)
  !      !print *, 'P3', myp, ':', g2, ',', n_left, n_right_tmp, ',', i_alpha, j_alpha_tmp, ',', prow_tmp, pcol_tmp
  !      call distribute_matrix_real_part(proc%context, n_left, n_right_tmp, i_alpha, j_alpha_tmp, &
  !           prow_tmp, pcol_tmp, n_left, YAY, A_alpha_desc, A_alpha)
  !      !call mpi_barrier(mpi_comm_world, ierr)
  !      !if (is_group_member .and. myp == g2 - 1) then
  !      !  print *, 'P4', myp, YAY
  !      !end if
  !    end do
  !
  !    wtime_end = mpi_wtime()
  !    call add_event('change_basis_lcao_diag_to_alpha_group_filter_real:inner_loop5', wtime_end - wtime_start)
  !    wtime_start = wtime_end
  !  end do
  !
  !  !call get_symmetricity(A_lcao_desc, A_lcao, symmetricity)
  !  !print *, 'A_lcao symmetricity', symmetricity
  !  !call get_symmetricity(A_alpha_desc, A_alpha, symmetricity)
  !  !print *, 'A_alpha symmetricity', symmetricity
  !
  !  !call pdlaprnt(dim, dim, A_lcao, 1, 1, A_lcao_desc, 0, 0, 'AAA', 6, work_pdlaprnt)
  !  !call pdlaprnt(num_filter, num_filter, A_alpha, 1, 1, A_alpha_desc, 0, 0, 'XXX', 6, work_pdlaprnt)
  !  A_alpha_cmplx(:, :) = A_alpha(:, :)
  !  deallocate(Y_local_right, A_lcao, A_alpha, A_lcao_local, Y_local_left, AY, YAY)
  end subroutine change_basis_lcao_to_alpha_group_filter


  subroutine change_basis_lcao_diag_to_alpha(basis, A_lcao_diag, A_alpha, A_alpha_desc)
    type(wk_basis_t), intent(in) :: basis
    real(8), intent(in) :: A_lcao_diag(:)
    integer, intent(in) :: A_alpha_desc(desc_size)
    real(8), intent(out) :: A_alpha(:, :)

    if (basis%is_group_filter_mode) then
      call change_basis_lcao_diag_to_alpha_group_filter(basis%filter_group_indices, basis%Y_local, &
           A_lcao_diag, A_alpha, A_alpha_desc)
    else
      call change_basis_lcao_diag_to_alpha_all(basis%Y_filtered, basis%Y_filtered_desc, &
       A_lcao_diag, A_alpha, A_alpha_desc)
    end if
  end subroutine change_basis_lcao_diag_to_alpha


  ! A_alpha <- Y^\dagger diag(A_lcao_diag) Y.
  ! Complexity: O(m n^2).
  subroutine change_basis_lcao_diag_to_alpha_all(Y_filtered, Y_filtered_desc, &
       A_lcao_diag, A_alpha, A_alpha_desc)
    real(8), intent(in) :: Y_filtered(:, :), A_lcao_diag(:)
    real(8), intent(out) :: A_alpha(:, :)
    integer, intent(in) :: Y_filtered_desc(desc_size), A_alpha_desc(desc_size)

    integer :: i, j, dim, num_filter, Y_dagger_desc(desc_size)
    real(8), allocatable :: Y_dagger(:, :)
    real(8) :: wtime_start, wtime_end

    wtime_start = mpi_wtime()

    dim = Y_filtered_desc(rows_)
    num_filter = Y_filtered_desc(cols_)
    call setup_distributed_matrix_real('', num_filter, dim, Y_dagger_desc, Y_dagger)

    wtime_end = mpi_wtime()
    call add_event('change_basis_lcao_diag_to_alpha:setup', wtime_end - wtime_start)
    wtime_start = wtime_end

    call pdtran(num_filter, dim, 1d0, &
         Y_filtered, 1, 1, Y_filtered_desc, &
         0d0, &
         Y_dagger, 1, 1, Y_dagger_desc)

    wtime_end = mpi_wtime()
    call add_event('change_basis_lcao_diag_to_alpha:pdtran', wtime_end - wtime_start)
    wtime_start = wtime_end

    do i = 1, dim
      call pdscal(num_filter, A_lcao_diag(i), Y_dagger, 1, i, Y_dagger_desc, 1)
    end do

    wtime_end = mpi_wtime()
    call add_event('change_basis_lcao_diag_to_alpha:pdscal', wtime_end - wtime_start)
    wtime_start = wtime_end

    call pdgemm('No', 'No', num_filter, num_filter, dim, 1d0, &
         Y_dagger, 1, 1, Y_dagger_desc, &
         Y_filtered, 1, 1, Y_filtered_desc, &
         0d0, &
         A_alpha, 1, 1, A_alpha_desc)

    wtime_end = mpi_wtime()
    call add_event('change_basis_lcao_diag_to_alpha:pdgemm', wtime_end - wtime_start)
  end subroutine change_basis_lcao_diag_to_alpha_all


  subroutine change_basis_lcao_diag_to_alpha_group_filter(filter_group_indices, Y_local, &
       A_lcao_diag, A_alpha, A_alpha_desc)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: Y_local
    real(8), intent(in) :: A_lcao_diag(:)
    real(8), intent(out) :: A_alpha(:, :)
    integer, intent(in) :: A_alpha_desc(desc_size), filter_group_indices(:, :)

    type(sparse_mat) :: A_lcao_sparse

    call diag_to_sparse_matrix(Y_local%m, A_lcao_diag, A_lcao_sparse)
    call matmul_bsb_d(Y_local, A_lcao_sparse, Y_local, 1d0, A_alpha, A_alpha_desc, 0d0)
  !
  !  integer :: num_groups, dim, num_filter
  !  integer :: g, p, nprow, npcol, myp, myg, myprow, mypcol, np, prow, pcol
  !  integer :: i_lcao, j_lcao, m_right, n_right, g2, i, ierr
  !  real(8), allocatable :: A_lcao_diag(:), A_alpha(:, :)  ! Distributed.
  !  real(8), allocatable :: Y_local_right(:, :), Y_local_left(:, :), YAY(:, :)  ! Local.
  !  integer :: prow_tmp, pcol_tmp, j_lcao_tmp, j_alpha_tmp, m_right_tmp, n_right_tmp, i_alpha, j_alpha, m_left, n_left
  !  integer :: blacs_pnum
  !  logical :: is_group_member
  !
  !  !call mpi_barrier(mpi_comm_world, ierr)
  !
  !  num_groups = size(filter_group_indices, 2) - 1
  !  dim = filter_group_indices(1, num_groups + 1) - 1
  !  num_filter = filter_group_indices(2, num_groups + 1) - 1
  !  call blacs_gridinfo(proc%context, nprow, npcol, myprow, mypcol)
  !  np = nprow * npcol
  !  myp = blacs_pnum(proc%context, myprow, mypcol)
  !  myg = myp + 1
  !  is_group_member = myg <= num_groups
  !
  !  allocate(A_lcao_diag(size(A_lcao_diag_cmplx)), &
  !       A_alpha(size(A_alpha_cmplx, 1), size(A_alpha_cmplx, 2)))
  !  A_lcao_diag(:) = real(A_lcao_diag_cmplx(:), kind(0d0))
  !  A_alpha(:, :) = 1d20  ! For debug.
  !
  !  if (is_group_member) then
  !    j_lcao = filter_group_indices(1, myg)
  !    j_alpha = filter_group_indices(2, myg)
  !    m_right = filter_group_indices(1, myg + 1) - j_lcao
  !    n_right = filter_group_indices(2, myg + 1) - j_alpha
  !    allocate(Y_local_right(m_right, n_right))
  !    Y_local_right(:, :) = real(Y_local_cmplx(1)%val(:, :), kind(0d0))
  !  else
  !    allocate(Y_local_right(0, 0))
  !  end if
  !
  !  do g = 1, num_groups
  !    p = g - 1
  !    call blacs_pcoord(proc%context, p, prow, pcol)
  !    i_lcao = filter_group_indices(1, g)
  !    i_alpha = filter_group_indices(2, g)
  !    m_left = filter_group_indices(1, g + 1) - i_lcao
  !    n_left = filter_group_indices(2, g + 1) - i_alpha
  !    if (g > 1) then
  !      deallocate(Y_local_left, YAY)
  !    end if
  !    if (is_group_member) then
  !      allocate(Y_local_left(m_left, n_left), YAY(n_left, n_right))
  !    else
  !      allocate(Y_local_left(m_left, n_left), YAY(0, 0))
  !    end if
  !
  !    ! Send Y_left.
  !    if (myp == p) then
  !      Y_local_left(:, :) = Y_local_right(:, :)
  !      call dgebs2d(proc%context, 'All', ' ', m_right, n_right, Y_local_left, m_right)
  !    else
  !      call dgebr2d(proc%context, 'All', ' ', m_left, n_left, Y_local_left, m_left, prow, pcol)
  !    end if
  !
  !    if (myg == g) then
  !      do i = 1, m_right
  !        Y_local_right(i, :) = Y_local_right(i, :) * A_lcao_diag(j_lcao + i - 1)
  !      end do
  !    else
  !      Y_local_right(:, :) = 0d0
  !    end if
  !
  !    if (is_group_member) then
  !      call dgemm('Trans', 'No', n_left, n_right, m_left, &
  !           1d0, Y_local_left, m_left, Y_local_right, m_left, &
  !           0d0, YAY, n_left)
  !    end if
  !
  !    do g2 = 1, num_groups
  !      j_alpha_tmp = filter_group_indices(2, g2)
  !      n_right_tmp = filter_group_indices(2, g2 + 1) - j_alpha_tmp
  !      call blacs_pcoord(proc%context, g2 - 1, prow_tmp, pcol_tmp)
  !      call mpi_barrier(mpi_comm_world, ierr)
  !      call distribute_matrix_real_part(proc%context, n_left, n_right_tmp, i_alpha, j_alpha_tmp, &
  !           prow_tmp, pcol_tmp, n_left, YAY, A_alpha_desc, A_alpha)
  !      call mpi_barrier(mpi_comm_world, ierr)
  !    end do
  !  end do
  !
  !  A_alpha_cmplx(:, :) = A_alpha(:, :)
  !  deallocate(Y_local_right, A_alpha, Y_local_left, YAY)
  end subroutine change_basis_lcao_diag_to_alpha_group_filter


  subroutine lcao_coef_to_alpha(S_sparse, basis, dv_psi, dv_alpha)
    type(sparse_mat), intent(in) :: S_sparse
    type(wk_basis_t), intent(in) :: basis
    complex(kind(0d0)), intent(in) :: dv_psi(basis%dim)
    complex(kind(0d0)), intent(out) :: dv_alpha(basis%num_basis)

    complex(kind(0d0)) :: dv_s_psi(basis%dim)

    ! c = Y^\dagger S \psi, \alpha <- c
    if (basis%is_group_filter_mode) then
      !call matvec_sd_z('No', S_sparse, kOne, dv_psi, kZero, dv_s_psi)
      !call matvec_bd_z('Trans', basis%Y_local, kOne, dv_s_psi, kZero, dv_alpha)
      call terminate('lcao_coef_to_alpha: group filter mode is not supported', 1)
    else
      dv_s_psi(:) = kZero
      call matvec_sd_z('No', S_sparse, kOne, dv_psi, kZero, dv_s_psi)
      call matvec_dd_z('Trans', basis%Y_filtered, basis%Y_filtered_desc, kOne, dv_s_psi, kZero, dv_alpha)
    end if
  end subroutine lcao_coef_to_alpha
end module wk_conversion_m
