module wk_initialization_m
  use mpi
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_global_variables_m
  use wk_processes_m
  use wk_atom_m
  use wk_charge_m
  use wk_conversion_m
  use wk_global_variables_m
  use wk_linear_algebra_m
  use wk_matrix_generation_m
  use wk_setting_m
  use wk_state_m
  use wk_util_m
  implicit none

  private
  public :: initialize, re_initialize_state, clear_offdiag_blocks, compute_energies

contains

  subroutine make_initial_psi_gauss(dim, dv_psi)
    integer, intent(in) :: dim
    complex(kind(0d0)), intent(out) :: dv_psi(:)

    real(8), parameter :: center = 0.5d0, sigma2 = 0.001d0, pi = atan(1.0d0) * 4.0d0
    integer :: i
    real(8) :: x

    do i = 1, dim
      x = index_to_coordinate(dim, i)
      dv_psi(i) = 1.0d0 / sqrt(2.0d0 * pi * sigma2) &
           * exp(- (x - center) ** 2.0d0 / (2.0d0 * sigma2))
    end do
  end subroutine make_initial_psi_gauss


  subroutine make_initial_psi_delta(dim, i, dv_psi)
    integer, intent(in) :: dim, i
    complex(kind(0d0)), intent(out) :: dv_psi(:)

    dv_psi(:) = kZero
    dv_psi(i) = kOne
  end subroutine make_initial_psi_delta


  ! Because group_id file describes (group -> [atom index]),
  ! atom_indices is needed to convert an atom index to orbital indices.
  subroutine make_localized_indices(setting, dim, group_id, atom_indices, localized_indices)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: dim, group_id(:, :), atom_indices(:)
    logical, intent(out) :: localized_indices(dim)
    integer :: i, j, atom

    localized_indices(:) = .false.
    if (trim(setting%init_type) == 'local_alpha_delta') then
      localized_indices(setting%localize_start : setting%localize_end) = .true.
    else if (trim(setting%init_type) == 'local_alpha_delta_group') then
      do i = setting%localize_start, setting%localize_end
        do j = 1, group_id(1, i)
          atom = group_id(j + 1, i)
          !print *, 'A', atom, atom_indices(atom), atom_indices(atom + 1) - 1
          localized_indices(atom_indices(atom) : atom_indices(atom + 1) - 1) = .true.
        end do
      end do
    end if
  end subroutine make_localized_indices


  subroutine print_localized_indices(dim, localized_indices)
    integer, intent(in) :: dim
    logical, intent(in) :: localized_indices(dim)
    integer :: i
    write (0, '(A, F16.6, A)', advance='no') &
         ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, '] initialize() : localized indices are '
    do i = 1, dim
      if (localized_indices(i)) then
        write(0, '(I0, " ")', advance='no') i
      end if
    end do
    write(0, *)
  end subroutine print_localized_indices


  integer function get_actual_alpha_delta_index(setting, dim, localized_indices, &
       alpha_delta_index)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: dim, alpha_delta_index
    logical, intent(in) :: localized_indices(dim)
    integer :: num_localized_indices, i

    ! alpha_delta_index == 0 means HOMO.
    if (trim(setting%init_type) == 'alpha_delta' .and. alpha_delta_index == 0) then
      get_actual_alpha_delta_index = dim / 2
    else if ((trim(setting%init_type) == 'local_alpha_delta' .or. &
         trim(setting%init_type) == 'local_alpha_delta_group') .and. &
         alpha_delta_index == 0) then
      ! Count the number of localized indices.
      num_localized_indices = 0
      do i = 1, dim
        if (localized_indices(i)) then
          num_localized_indices = num_localized_indices + 1
        end if
      end do
      ! "No orbital is localized" is considered "All orbitals are localized".
      if (num_localized_indices == 0) then
        num_localized_indices = dim
      end if
      ! "+ 1" below avoids that num_localized_indices = 1 leads to
      ! get_actual_alpha_delta_index = 0.
      get_actual_alpha_delta_index = (num_localized_indices + 1) / 2
    else
      get_actual_alpha_delta_index = alpha_delta_index
    end if
    if (check_master()) then
      write (0, '(A, F16.6, A, I0)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] initialize() : set index of initial eigenstate to ', get_actual_alpha_delta_index
    end if
  end function get_actual_alpha_delta_index


  ! Calculate errors:
  ! absolute_filter_error = ||psi (col_psi_2) - psi' (col_psi)|| and
  ! relative_filter_error = absolute_filter_error / ||psi||.
  ! Brakes psi (col_psi_2).
  subroutine get_filtering_errors(dim, dv_psi_orig, dv_psi_filtered, errors)
    integer, intent(in) :: dim
    complex(kind(0d0)), intent(in) :: dv_psi_orig(dim), dv_psi_filtered(dim)
    type(wk_error_t), intent(out) :: errors

    complex(kind(0d0)) :: dv_psi_diff(dim)
    real(8) :: psi_norm
    real(8) :: dznrm2  ! Function.

    call check_nan_vector('get_filtering_errors orig real', dreal(dv_psi_orig))
    call check_nan_vector('get_filtering_errors orig imag', aimag(dv_psi_orig))
    call check_nan_vector('get_filtering_errors filtered real', dreal(dv_psi_filtered))
    call check_nan_vector('get_filtering_errors filtered imag', aimag(dv_psi_filtered))
    dv_psi_diff(1 : dim) = dv_psi_orig(1 : dim) - dv_psi_filtered(1 : dim)
    errors%absolute = dznrm2(dim, dv_psi_diff, 1)
    psi_norm = dznrm2(dim, dv_psi_orig, 1)
    errors%relative = errors%absolute / psi_norm
  end subroutine get_filtering_errors


  ! Calculate errors:
  ! absolute_filter_error = ||alpha (col_filtered_work1) - alpha' (col_alpha)|| and
  ! relative_filter_error = absolute_filter_error / ||alpha||.
  ! Brakes alpha (col_filtered_work1).
  !subroutine get_re_initialize_errors(dim, filtered_vecs, filtered_vecs_desc, absolute_error, relative_error)
  !  integer, intent(in) :: dim
  !  complex(kind(0d0)), intent(inout) :: filtered_vecs(:, :)
  !  integer, intent(in) :: filtered_vecs_desc(desc_size)
  !  real(8), intent(out) :: absolute_error, relative_error
  !
  !  integer :: nprow, npcol, myrow, mycol, prow_owns_error, pcol_owns_error, pnum_owns_error, ierr
  !  real(8) :: alpha_norm
  !  integer :: indxg2p, blacs_pnum
  !
  !  call blacs_gridinfo(filtered_vecs_desc(context_), nprow, npcol, myrow, mycol)
  !  prow_owns_error = 0
  !  pcol_owns_error = indxg2p(col_filtered_work1, filtered_vecs_desc(block_col_), 0, filtered_vecs_desc(csrc_), npcol)
  !  pnum_owns_error = blacs_pnum(filtered_vecs_desc(context_), prow_owns_error, pcol_owns_error)
  !  call pdznrm2(dim, alpha_norm, filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1)
  !  call mpi_bcast(alpha_norm, 1, mpi_double_precision, pnum_owns_error, mpi_comm_world, ierr)
  !  call pzaxpy(dim, -kOne, filtered_vecs, 1, col_alpha, filtered_vecs_desc, 1, &
  !       filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1)
  !  call pdznrm2(dim, absolute_error, filtered_vecs, 1, col_filtered_work1, filtered_vecs_desc, 1)
  !  call mpi_bcast(absolute_error, 1, mpi_double_precision, pnum_owns_error, mpi_comm_world, ierr)
  !  relative_error = absolute_error / alpha_norm
  !end subroutine get_re_initialize_errors


  !subroutine reconcile_from_lcao_coef(num_filter, S_sparse, Y_filtered, Y_filtered_desc, &
  !     dv_psi, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size)
  !  real(8), intent(in) :: Y_filtered(:, :)
  !  type(sparse_mat), intent(in) :: S_sparse
  !  complex(kind(0d0)), intent(in) :: dv_psi(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  complex(kind(0d0)) :: dv_alpha(num_filter)
  !
  !  call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, dv_psi, dv_alpha_reconcile)
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha_reconcile, dv_psi_reconcile)
  !end subroutine reconcile_from_lcao_coef
  !
  !
  !subroutine reconcile_from_lcao_coef_cutoff(num_filter, cutoff, &
  !     S_sparse, Y_filtered, Y_filtered_desc, &
  !     dv_psi, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size)
  !  real(8), intent(in) :: cutoff, Y_filtered(:, :)
  !  type(sparse_mat), intent(in) :: S_sparse
  !  complex(kind(0d0)), intent(in) :: dv_psi(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  complex(kind(0d0)) :: dv_alpha(num_filter)
  !
  !  call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, dv_psi, dv_alpha)
  !  call cutoff_vector(num_filter, cutoff, dv_alpha, dv_alpha_reconcile)
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha_reconcile, dv_psi_reconcile)
  !end subroutine reconcile_from_lcao_coef_cutoff
  !
  !
  !subroutine reconcile_from_lcao_coef_suppress(num_filter, suppress_constant, &
  !     S_sparse, eigenvalues, Y_filtered, Y_filtered_desc, &
  !     dv_psi, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size)
  !  real(8), intent(in) :: suppress_constant, eigenvalues(:), Y_filtered(:, :)
  !  type(sparse_mat), intent(in) :: S_sparse
  !  complex(kind(0d0)), intent(in) :: dv_psi(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  integer :: i, min_abs_i, ierr
  !  complex(kind(0d0)) :: dv_alpha(num_filter), diff
  !  real(8) :: min_abs
  !
  !  call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, dv_psi, dv_alpha)
  !  i = 0
  !  min_abs = -1d0
  !  do i = 1, num_filter
  !    if (abs(dv_alpha(i)) > min_abs) then
  !      min_abs = abs(dv_alpha(i))
  !      min_abs_i = i
  !    end if
  !  end do
  !  do i = 1, num_filter
  !    diff = eigenvalues(i) - eigenvalues(min_abs_i)
  !    dv_alpha_reconcile(i) = dv_alpha(i) * exp(- abs(diff * suppress_constant) ** 2d0)
  !  end do
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha_reconcile, dv_psi_reconcile)
  !end subroutine reconcile_from_lcao_coef_suppress


  subroutine reconcile_from_alpha_matrix_suppress(t, suppress_constant, &
       dv_eigenvalues_prev, dv_eigenvalues, basis, YSY_filtered, YSY_filtered_desc, &
       dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: YSY_filtered_desc(desc_size)
    type(wk_basis_t), intent(in) :: basis
    real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
    real(8), intent(in) :: YSY_filtered(:, :)
    complex(kind(0d0)), intent(in) :: dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc, ks(basis%num_basis)
    real(8) :: dv_suppress_factor(basis%num_basis), suppress_factor_sum, dv_suppress_factor_copy(basis%num_basis)
    real(8) :: dv_suppress_factor2(basis%num_basis), suppress_factor_sum2
    complex(kind(0d0)) :: dv_evcoef_amplitude(basis%num_basis), dv_evcoef_amplitude_reconcile(basis%num_basis)
    real(8), allocatable :: ENE_suppress(:, :), YSY_filtered_suppress(:, :)
    real(8) :: ysy_ipratios(basis%num_basis), ysy_suppress_ipratios(basis%num_basis)
    real(8) :: average_pratio, weighted_average_pratio, average_pratio_suppress, weighted_average_pratio_suppress
    real(8) :: work(1000), energy_normalizer, diff
    integer, save :: count = 0
    character(len=50) :: filename
    character(len=6) :: count_str
    character(len=12) :: time_str
    integer, parameter :: iunit_YSY = 20
    real(8) :: dznrm2

    call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
    allocate(ENE_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
    allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
    ENE_suppress(:, :) = 0d0

    do j = 1, basis%num_basis
      do i = 1, basis%num_basis
        diff = dv_eigenvalues(i) - dv_eigenvalues_prev(j)
        dv_suppress_factor(i) = exp(- (suppress_constant * diff) ** 2d0)
      end do
      call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor', dv_suppress_factor)
      do i = 1, basis%num_basis
        call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
        if (myrow == rsrc .and. mycol == csrc) then
          YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
        end if
      end do
    end do

    !call get_ipratio_of_eigenstates(YSY_filtered, YSY_filtered_desc, ysy_ipratios)
    !call get_ipratio_of_eigenstates(YSY_filtered_suppress, YSY_filtered_desc, ysy_suppress_ipratios)
    !average_pratio = 0d0
    !weighted_average_pratio = 0d0
    !average_pratio_suppress = 0d0
    !weighted_average_pratio_suppress = 0d0
    !do j = 1, basis%num_basis
    !  average_pratio = average_pratio + 1d0 / ysy_ipratios(j)
    !  average_pratio_suppress = average_pratio_suppress + 1d0 / ysy_suppress_ipratios(j)
    !  weighted_average_pratio = weighted_average_pratio + (abs(dv_alpha(j)) ** 2d0) / ysy_ipratios(j)
    !  weighted_average_pratio_suppress = weighted_average_pratio_suppress + &
    !       (abs(dv_alpha(j)) ** 2d0) / ysy_suppress_ipratios(j)
    !end do
    !average_pratio = average_pratio / basis%num_basis
    !weighted_average_pratio = weighted_average_pratio
    !average_pratio_suppress = average_pratio_suppress / basis%num_basis
    !weighted_average_pratio_suppress = weighted_average_pratio_suppress
    !print *, 'YSYPRATIO1 average pratio', t, average_pratio
    !print *, 'YSYPRATIO2 weighted average pratio', t, weighted_average_pratio
    !print *, 'YSYPRATIO3 average pratio (suppressed)', t, average_pratio_suppress
    !print *, 'YSYPRATIO4 weighted average pratio (suppressed)', t, weighted_average_pratio_suppress

    !write(count_str, '(I6.6)') count
    !write(time_str, '(F0.4)') t * 2.418884326505e-5  ! kPsecPerAu.
    !filename = 'YSYsuppressed_' // count_str // '_' // trim(time_str) // 'ps.mtx'
    !if (check_master()) then
    !  open(iunit_YSY, file=trim(filename))
    !end if
    !call pdlaprnt(basis%num_basis, basis%num_basis, YSY_filtered_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'YSY', iunit_YSY, work)
    !if (check_master()) then
    !  close(iunit_YSY)
    !end if

    dv_alpha_reconcile(:) = kZero
    print *, '------------------ ZZZZstart', t * 2.418884326505d-5, 'ps ------------------'
    call matvec_dd_z('No', YSY_filtered_suppress, YSY_filtered_desc, kOne, dv_alpha, kZero, dv_alpha_reconcile)
    !print *, 'Energy correction X', dv_evcoef
    !print *, 'Energy correction Y', dv_evcoef_reconcile
    !
    !do j = 1, basis%num_basis
    !  dv_evcoef_amplitude(j) = conjg(dv_evcoef(j)) * dv_evcoef(j)
    !  do i = 1, basis%num_basis
    !    if (abs(dv_evcoef_reconcile(i)) < 1d-10) then  ! αの値ではなくsmoothed deltaの方で切り落としを定める.
    !      !print *, 'ZZZZ i: alphaprime0(zero), abs(alphaprime0)', i, ':', dv_evcoef_reconcile(i), abs(dv_evcoef_reconcile(i))
    !      dv_suppress_factor2(i) = kZero
    !    else
    !      dv_suppress_factor2(i) = exp(- 1d4 * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
    !    end if
    !  end do
    !  call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor2', dv_suppress_factor2)
    !  do i = 1, basis%num_basis
    !    call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
    !    if (myrow == rsrc .and. mycol == csrc) then
    !      ENE_suppress(i_local, j_local) = dv_suppress_factor2(i)
    !    end if
    !  end do
    !end do
    !call scale_columns_to_stochastic_matrix(YSY_filtered_desc, ENE_suppress)
    !call check_nan_matrix('ENE_suppress', ENE_suppress)
    !
    !filename = 'ENEsuppress_' // count_str // '_' // trim(time_str) // 'ps.mtx'
    !if (check_master()) then
    !  open(iunit_YSY, file=trim(filename))
    !end if
    !call pdlaprnt(basis%num_basis, basis%num_basis, ENE_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'ENEsuppress', &
    !     iunit_YSY, work)
    !if (check_master()) then
    !  close(iunit_YSY)
    !end if
    count = count + 1
    !
    !dv_evcoef_amplitude_reconcile(:) = kZero
    !call matvec_dd_z('No', ENE_suppress, YSY_filtered_desc, kOne, dv_evcoef_amplitude, kZero, dv_evcoef_amplitude_reconcile)
    !print *, 'Energy correction A', dreal(dv_evcoef_amplitude)
    !print *, 'Energy correction B', dreal(dv_evcoef_amplitude_reconcile)
    !do i = 1, basis%num_basis
    !  energy_normalizer = dsqrt(dreal(dv_evcoef_amplitude_reconcile(i))) / abs(dv_evcoef_reconcile(i))
    !  if (abs(dv_evcoef_reconcile(i)) < 1d-10) then
    !    !print *, 'Energy correction(1) i: alphaprime0, abs(alphaprime0), targetamp', i, ':', &
    !    !     dv_evcoef_reconcile(i), &
    !    !     abs(dv_evcoef_reconcile(i)), dsqrt(dreal(dv_evcoef_amplitude_reconcile(i)))
    !    print *, 'Energy correction(1) i: abs(alphaprime0), targetamp', i, ':', &
    !         abs(dv_evcoef_reconcile(i)), dsqrt(dreal(dv_evcoef_amplitude_reconcile(i)))
    !  else
    !    !print *, 'Energy correction(2) i: alphaprime0, abs(alphaprime0), targetamp, normalizer', i, ':', &
    !    !     dv_evcoef_reconcile(i), &
    !    !     abs(dv_evcoef_reconcile(i)), dsqrt(dreal(dv_evcoef_amplitude_reconcile(i))), energy_normalizer
    !    print *, 'Energy correction(2) i: abs(alphaprime0), targetamp, normalizer', i, ':', &
    !         abs(dv_evcoef_reconcile(i)), dsqrt(dreal(dv_evcoef_amplitude_reconcile(i))), energy_normalizer
    !    dv_evcoef_reconcile(i) = dv_evcoef_reconcile(i) * energy_normalizer
    !  end if
    !end do

    if (check_master()) then
      write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] reconcile_from_alpha_matrix_suppress() : evcoef norm before normalization ', &
           dznrm2(basis%num_basis, dv_alpha_reconcile, 1)
    end if
    call normalize_vector(basis%num_basis, dv_alpha_reconcile)
    call alpha_to_lcao_coef(basis, dv_alpha_reconcile, dv_psi_reconcile)
    deallocate(YSY_filtered_suppress, ENE_suppress)
  end subroutine reconcile_from_alpha_matrix_suppress


  !subroutine reconcile_from_alpha_matrix_suppress_orthogonal(num_filter, t, suppress_constant, &
  !     dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
  !     dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
  !  real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
  !  real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
  !  complex(kind(0d0)), intent(in) :: dv_alpha(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc
  !  real(8) :: dv_suppress_factor(num_filter)
  !  real(8), allocatable :: YSY_filtered_orthogonal(:, :), YSY_filtered_suppress(:, :)
  !  real(8) :: dznrm2
  !
  !  call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
  !  allocate(YSY_filtered_orthogonal(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !  allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !  YSY_filtered_orthogonal(:, :) = 0d0
  !
  !  do j = 1, num_filter
  !    do i = 1, num_filter
  !      dv_suppress_factor(i) = exp(- suppress_constant * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
  !    end do
  !    call check_nan_vector('reconcile_from_alpha_matrix_suppress_orthogonal dv_suppress_factor', dv_suppress_factor)
  !    do i = 1, num_filter
  !      call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
  !      if (myrow == rsrc .and. mycol == csrc) then
  !        YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
  !      end if
  !    end do
  !  end do
  !
  !  dv_alpha_reconcile(:) = kZero
  !  call matvec_nearest_orthonormal_matrix(YSY_filtered_suppress, YSY_filtered_desc, &
  !       dv_alpha, dv_alpha_reconcile, YSY_filtered_orthogonal)
  !
  !  if (check_master()) then
  !    write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
  !         '] reconcile_from_alpha_matrix_suppress_orthogonal() : evcoef norm before normalization ', &
  !         dznrm2(num_filter, dv_alpha_reconcile, 1)
  !  end if
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha_reconcile, dv_psi_reconcile)
  !  deallocate(YSY_filtered_suppress, YSY_filtered_orthogonal)
  !end subroutine reconcile_from_alpha_matrix_suppress_orthogonal


  !subroutine reconcile_from_alpha_matrix_suppress_adaptive(num_filter, t, suppress_constant, &
  !     is_first_in_multiple_initials, &
  !     dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
  !     dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
  !  logical, intent(in) :: is_first_in_multiple_initials
  !  real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
  !  real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
  !  complex(kind(0d0)), intent(in) :: dv_alpha(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc
  !  real(8) :: dv_suppress_factor(num_filter)
  !  real(8), allocatable, save :: YSY_filtered_suppress(:, :)
  !  type(wk_energy_t) :: energies
  !
  !  integer, save :: count = 0
  !  character(len=100) :: filename
  !  character(len=6) :: count_str
  !  character(len=12) :: time_str
  !  character(len=12) :: suppress_str
  !  integer, parameter :: iunit_YSY = 20
  !  real(8) :: energy_diff, work(1000)
  !  real(8), parameter :: relax_constant_for_suppression = 1d-10
  !  real(8) :: energy_mean(Y_filtered_desc(cols_)), energy_dev(Y_filtered_desc(cols_))
  !  ! Function.
  !  real(8) :: dznrm2
  !
  !  if (check_master()) then
  !    write (0, '(A, F16.6, A, F16.6, A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
  !         '] reconcile_from_alpha_matrix_suppress_adaptive() : start for time ', t, ' au, ', &
  !         t * kPsecPerAu, ' ps'
  !  end if
  !
  !  call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
  !
  !  if (is_first_in_multiple_initials) then
  !    if (allocated(YSY_filtered_suppress)) then
  !      deallocate(YSY_filtered_suppress)
  !    end if
  !    allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !
  !    call get_moment_for_each_column(dv_eigenvalues, YSY_filtered_desc, YSY_filtered, energy_mean, energy_dev)
  !    call check_nan_vector('reconcile_from_alpha_matrix_suppress_adaptive energy_mean', energy_mean)
  !    call check_nan_vector('reconcile_from_alpha_matrix_suppress_adaptive energy_dev', energy_dev)
  !
  !    do j = 1, num_filter
  !      do i = 1, num_filter
  !        energy_diff = dv_eigenvalues(i) - energy_mean(j)
  !        dv_suppress_factor(i) = exp(- suppress_constant * &
  !             (energy_diff / (energy_dev(j) + relax_constant_for_suppression)) ** 2d0)
  !        !if (dv_suppress_factor(i) > 1e-100) then
  !        !  print *, 'ZZZZZZ', j, i, suppress_constant, energy_mean(j), energy_diff, &
  !        !       energy_dev(j), dv_suppress_factor(i)
  !        !end if
  !      end do
  !      call check_nan_vector('reconcile_from_alpha_matrix_suppress_adaptive dv_suppress_factor', dv_suppress_factor)
  !      do i = 1, num_filter
  !        call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
  !        if (myrow == rsrc .and. mycol == csrc) then
  !          YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
  !        end if
  !      end do
  !    end do
  !
  !    !write(count_str, '(I6.6)') count
  !    !if (mod(count, 6) == 0) then
  !    !  write(time_str, '(F0.5)') t * 2.418884326505e-5  ! kPsecPerAu.
  !    !  write(suppress_str, '(F0.5)') suppress_constant
  !    !  filename = 'YSY_' // trim(count_str) // '_' // trim(suppress_str) // '_' //  trim(time_str) // 'ps.mtx'
  !    !  if (check_master()) then
  !    !    open(iunit_YSY, file=trim(filename))
  !    !  end if
  !    !  call pdlaprnt(num_filter, num_filter, YSY_filtered, 1, 1, YSY_filtered_desc, 0, 0, 'YSY', iunit_YSY, work)
  !    !  if (check_master()) then
  !    !    close(iunit_YSY)
  !    !  end if
  !    !  filename = 'YSYsuppressed_' // trim(count_str) // '_' // trim(suppress_str) // '_' // trim(time_str) // 'ps.mtx'
  !    !  if (check_master()) then
  !    !    open(iunit_YSY, file=trim(filename))
  !    !  end if
  !    !  call pdlaprnt(num_filter, num_filter, YSY_filtered_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'YSYf', iunit_YSY, work)
  !    !  if (check_master()) then
  !    !    close(iunit_YSY)
  !    !  end if
  !    !end if
  !    !count = count + 1
  !  end if
  !
  !  dv_alpha_reconcile(:) = kZero
  !  call matvec_dd_z('No', YSY_filtered_suppress, YSY_filtered_desc, kOne, dv_alpha, kZero, dv_alpha_reconcile)
  !
  !  if (check_master()) then
  !    write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
  !         '] reconcile_from_alpha_matrix_suppress() : evcoef norm before normalization ', &
  !         dznrm2(num_filter, dv_alpha_reconcile, 1)
  !  end if
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha_reconcile, dv_psi_reconcile)
  !end subroutine reconcile_from_alpha_matrix_suppress_adaptive



  !subroutine reconcile_from_alpha_matrix_suppress_select(num_filter, t, suppress_constant, &
  !     is_first_in_multiple_initials, &
  !     dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
  !     dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
  !  logical, intent(in) :: is_first_in_multiple_initials
  !  real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
  !  real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
  !  complex(kind(0d0)), intent(in) :: dv_alpha(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc, max_index_in_col, max_index_in_row
  !  real(8), allocatable, save :: YSY_filtered_select(:, :), YSY_filtered_orthogonal(:, :)
  !  integer, allocatable, save :: col_to_max_index_in_col(:), row_to_max_index_in_row(:)
  !  integer :: tmp_max_index(2), ierr
  !  real(8) :: col_vector_local(num_filter, 1), row_vector_local(1, num_filter), suppress_factor, norm
  !
  !  integer, save :: count = 0
  !  character(len=100) :: filename
  !  character(len=6) :: count_str
  !  character(len=12) :: time_str
  !  character(len=12) :: suppress_str
  !  integer, parameter :: iunit_YSY = 20
  !  real(8) :: energy_diff, energy_diff_prev, energy_diff_diff, work(1000)
  !  ! Function.
  !  integer :: blacs_pnum
  !  real(8) :: dznrm2
  !
  !  if (check_master()) then
  !    write (0, '(A, F16.6, A, F16.6, A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
  !         '] reconcile_from_alpha_matrix_suppress_select() : start for time ', t, ' au, ', &
  !         t * kPsecPerAu, ' ps'
  !  end if
  !
  !  call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
  !
  !  if (is_first_in_multiple_initials) then
  !    print *, 'ZZZZZZZX'
  !    if (allocated(YSY_filtered_select)) then
  !      deallocate(YSY_filtered_select)
  !    end if
  !    if (allocated(YSY_filtered_orthogonal)) then
  !      deallocate(YSY_filtered_orthogonal)
  !    end if
  !    allocate(YSY_filtered_select(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !    allocate(YSY_filtered_orthogonal(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !    YSY_filtered_select(:, :) = YSY_filtered(:, :)
  !    YSY_filtered_orthogonal(:, :) = 0d0
  !
  !    if (allocated(col_to_max_index_in_col)) then
  !      deallocate(col_to_max_index_in_col)
  !    end if
  !    if (allocated(row_to_max_index_in_row)) then
  !      deallocate(row_to_max_index_in_row)
  !    end if
  !    allocate(col_to_max_index_in_col(num_filter), row_to_max_index_in_row(num_filter))
  !
  !    do j = 1, num_filter
  !      call gather_matrix_real_part(YSY_filtered, YSY_filtered_desc, 1, j, num_filter, 1, &
  !           0, 0, col_vector_local)
  !      if (myrow == 0 .and. mycol == 0) then
  !        tmp_max_index = maxloc(abs(col_vector_local))
  !        col_to_max_index_in_col(j) = tmp_max_index(1)
  !      end if
  !    end do
  !    do i = 1, num_filter
  !      call gather_matrix_real_part(YSY_filtered, YSY_filtered_desc, i, 1, 1, num_filter, &
  !           0, 0, row_vector_local)
  !      if (myrow == 0 .and. mycol == 0) then
  !        tmp_max_index = maxloc(abs(row_vector_local))
  !        row_to_max_index_in_row(i) = tmp_max_index(2)
  !      end if
  !    end do
  !    call mpi_bcast(col_to_max_index_in_col, num_filter, mpi_integer, 0, mpi_comm_world, ierr)
  !    call mpi_bcast(row_to_max_index_in_row, num_filter, mpi_integer, 0, mpi_comm_world, ierr)
  !    print *, 'ZZZZZZcol_to_max_index_in_col ', col_to_max_index_in_col
  !    print *, 'ZZZZZZrow_to_max_index_in_row ', row_to_max_index_in_row
  !
  !    do j = 1, num_filter
  !      do i = 1, num_filter
  !        max_index_in_col = col_to_max_index_in_col(j)
  !        max_index_in_row = row_to_max_index_in_row(i)
  !        if (i /= max_index_in_col .or. j /= max_index_in_row) then
  !          call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
  !          if (myrow == rsrc .and. mycol == csrc) then
  !            energy_diff = dv_eigenvalues(i) - dv_eigenvalues(max_index_in_col)
  !            energy_diff_prev = dv_eigenvalues_prev(max_index_in_row) - dv_eigenvalues_prev(j)
  !            energy_diff_diff = abs(energy_diff - energy_diff_prev)
  !            !if (abs(i - max_index_in_col) > 3) then
  !            !  suppress_factor = 0d0
  !            !else
  !            if (i == max_index_in_col) then
  !              suppress_factor = 1d0 / abs(YSY_filtered_select(i_local, j_local))
  !            else
  !              if (energy_diff_diff < 1d-16 .or. suppress_constant / energy_diff_diff > 37d0) then
  !                suppress_factor = 0d0
  !              else
  !                suppress_factor = exp(- suppress_constant / energy_diff_diff)
  !              end if
  !            end if
  !            YSY_filtered_select(i_local, j_local) = YSY_filtered_select(i_local, j_local) * suppress_factor
  !          end if
  !        end if
  !      end do
  !      !call pdnrm2(num_filter, norm, YSY_filtered_select, 1, j, YSY_filtered_desc, 1)
  !      !call infog2l(1, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
  !      !call mpi_bcast(norm, 1, mpi_double_precision, blacs_pnum(YSY_filtered_desc(context_), rsrc, csrc), &
  !      !     mpi_comm_world, ierr)
  !      !print *, 'ZZZZZZZZYSY_filtered_select_norm', j, norm
  !      !if (norm >= 1d-16) then
  !      !  call pdscal(num_filter, 1d0 / norm, YSY_filtered_select, 1, j, YSY_filtered_desc, 1)
  !      !end if
  !    end do
  !  end if
  !
  !  !write(count_str, '(I6.6)') count
  !  !if (mod(count, 6) == 0) then
  !  !  write(time_str, '(F0.5)') t * 2.418884326505e-5  ! kPsecPerAu.
  !  !  write(suppress_str, '(F0.5)') suppress_constant
  !  !  filename = 'YSY' // trim(count_str) // '_' // trim(suppress_str) // '_' // trim(time_str) // 'ps.mtx'
  !  !  if (check_master()) then
  !  !    open(iunit_YSY, file=trim(filename))
  !  !  end if
  !  !  call pdlaprnt(num_filter, num_filter, YSY_filtered, 1, 1, YSY_filtered_desc, 0, 0, 'YSY', iunit_YSY, work)
  !  !  if (check_master()) then
  !  !    close(iunit_YSY)
  !  !  end if
  !  !  filename = 'YSYselect' // trim(count_str) // '_' // trim(suppress_str) // '_' // trim(time_str) // 'ps.mtx'
  !  !  if (check_master()) then
  !  !    open(iunit_YSY, file=trim(filename))
  !  !  end if
  !  !  call pdlaprnt(num_filter, num_filter, YSY_filtered_select, 1, 1, YSY_filtered_desc, 0, 0, 'YSYs', iunit_YSY, work)
  !  !  if (check_master()) then
  !  !    close(iunit_YSY)
  !  !  end if
  !  !end if
  !  !count = count + 1
  !  !!stop
  !
  !  dv_alpha_reconcile(:) = kZero
  !  !call matvec_nearest_orthonormal_matrix(YSY_filtered_select, YSY_filtered_desc, &
  !  !     dv_evcoef, dv_evcoef_reconcile, YSY_filtered_orthogonal)
  !  call matvec_dd_z('No', YSY_filtered_select, YSY_filtered_desc, kOne, dv_alpha, kZero, dv_alpha_reconcile)
  !
  !  print *, 'ZZZZTT', row_to_max_index_in_row
  !  do i = 1, num_filter
  !    j = row_to_max_index_in_row(i)
  !    if (abs(dv_alpha(j)) > 1d-16) then
  !      dv_alpha_reconcile(i) = abs(dv_alpha_reconcile(i))! * (dv_alpha(j) / abs(dv_alpha(j)))
  !    end if
  !  end do
  !
  !  if (check_master()) then
  !    write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
  !         '] reconcile_from_alpha_matrix_suppress_select() : evcoef norm before normalization ', &
  !         dznrm2(num_filter, dv_alpha_reconcile, 1)
  !  end if
  !  if (dznrm2(num_filter, dv_alpha_reconcile, 1) < 1d-16) then
  !    if (check_master()) then
  !      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
  !           '] reconcile_from_alpha_matrix_suppress_select() : evcoef norm missing'
  !    end if
  !    dv_alpha_reconcile(:) = dv_alpha(:)
  !  end if
  !
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha_reconcile, dv_psi_reconcile)
  !end subroutine reconcile_from_alpha_matrix_suppress_select


  subroutine set_initial_value(setting, dim, structure, group_id, &
       H_sparse, S_sparse, basis, &
       dv_psi, dv_alpha, errors)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: dim, group_id(:, :)
    type(sparse_mat), intent(in) :: H_sparse, S_sparse
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(in) :: basis
    complex(kind(0d0)), intent(out) :: dv_psi(:, :), dv_alpha(:, :)
    type(wk_error_t), intent(out) :: errors(:)

    integer :: nprow, npcol, myrow, mycol, ierr
    real(8) :: min_msd, eigenvalues(dim)
    complex(kind(0d0)) :: diag_tmp, msd_temp, mean_temp
    integer :: indxg2p
    real(8), allocatable :: Y_localized(:, :)
    integer :: Y_localized_desc(desc_size), i, actual_alpha_delta_index, i_min_msd
    ! i-th element is .true. <=> i-th element is localized.
    logical :: localized_indices(dim)

    if (setting%is_restart_mode) then
      call terminate('not implemented', 47)
      !do i = 1, dim
      !  call pzelset(full_vecs, i, col_psi2, full_vecs_desc, setting%restart_psi(i))
      !end do
      !! State in each time step is saved after time count is incremented (see main.f90),
      !! so the time for conversion between LCAO coefficient and eigenstate expansion
      !! is setting%restart_t - setting%delta_t, not setting%restart_t.
      !call reconcile_from_lcao_coef(setting, setting%restart_t - setting%delta_t, &
      ! S_sparse, Y_filtered, Y_filtered_desc, &
      ! full_vecs, full_vecs_desc, filtered_vecs, filtered_vecs_desc)
      !call get_filtering_errors(dim, full_vecs, full_vecs_desc, absolute_filter_error, relative_filter_error)
    else if (setting%init_type(1:5) == 'alpha') then
      ! 初期値を \alpha で決める.
      ! 時間発展を計算するのは \alpha を介してであり, \psi と Mulliken charge は
      ! \alpha から固有ベクトル Y を用いて変換することで得る.
      if (trim(setting%init_type) == 'alpha_delta') then
        actual_alpha_delta_index = get_actual_alpha_delta_index(setting, &
             dim, localized_indices, setting%alpha_delta_index)
        call make_initial_psi_delta(dim, actual_alpha_delta_index - setting%fst_filter + 1, dv_alpha(:, 1))
      else if (trim(setting%init_type) == 'alpha_delta_multiple') then
        do i = 1, setting%num_multiple_initials
          actual_alpha_delta_index = get_actual_alpha_delta_index(setting, &
               dim, localized_indices, setting%alpha_delta_multiple_indices(i))
          call make_initial_psi_delta(dim, actual_alpha_delta_index - setting%fst_filter + 1, dv_alpha(:, i))
        end do
      !else if (trim(setting%init_type) == 'alpha_gauss') then
      !  call make_initial_psi_gauss(filtered_vecs, filtered_vecs_desc, col_alpha)
      !else if (trim(setting%init_type) == 'alpha_file') then
      !  call terminate('initialize: not implemented yet', 1)
      !  !call read_vector(setting%num_filter, trim(setting%alpha_filename), alpha)
      !else if (trim(setting%init_type) == 'alpha_delta_low_msd') then
      !  i_min_msd = 0
      !  min_msd = 1d100
      !  do i = 1, setting%num_filter
      !    call pzelget('All', ' ', msd_temp, filtered_vecs, i, col_eigenstate_msd_total, filtered_vecs_desc)
      !    call pzelget('All', ' ', mean_temp, filtered_vecs, i, col_eigenstate_mean_x, filtered_vecs_desc)
      !    if (dble(mean_temp) < setting%alpha_delta_min_x .or. setting%alpha_delta_max_x < dble(mean_temp)) then
      !      cycle
      !    end if
      !    if (dble(msd_temp) < min_msd) then
      !      min_msd = dble(msd_temp)
      !      i_min_msd = i
      !    end if
      !  end do
      !  if (i_min_msd == 0) then
      !    call terminate('initialize: no eigenstate in the specified region', 1)
      !  else
      !    call make_initial_psi_delta(i_min_msd, filtered_vecs, filtered_vecs_desc, col_alpha)
      !  end if
      else
        stop 'unknown initialization method (fail in command line option parsing)'
      end if
      do i = 1, setting%num_multiple_initials
        ! \alpha のノルムはプログラムがうまく動いていれば保存されるはずなので規格化は最初の一度のみ行う.
        ! 1 =: \Psi^\dagger S \Psi = c^\dagger Y^\dagger S Y c = ||c||^2 = ||\alpha||^2
        call normalize_vector(setting%num_filter, dv_alpha(:, i))
        ! フィルターした後の状態の重ね合わせで \psi を決めているのでこの時点では
        ! \psi に対するエラーは存在しない.
        call alpha_to_lcao_coef(basis, dv_alpha(:, i), dv_psi(:, i))
      end do
      if (setting%to_multiply_phase_factor) then
        call terminate('not implemented', 49)
      !  ! LCAO -> (phase factor multiplication) -> LCAO -> alpha -> (normalize) ->
      !  ! alpha -> LCAO. Without this procedure, Mulliken charge error occurs.
      !  call multiply_phase_factors(num_atoms, atom_indices, atom_coordinates, &
      !       setting%phase_factor_coef, full_vecs, full_vecs_desc, col_psi2)
      !  call reconcile_from_lcao_coef(setting, 0d0, &
      !       S_sparse, Y_filtered, Y_filtered_desc, &
      !       full_vecs, full_vecs_desc, filtered_vecs, filtered_vecs_desc)
      !  call get_filtering_errors(dim, full_vecs, full_vecs_desc, absolute_filter_error, relative_filter_error)
      else
        do i = 1, setting%num_multiple_initials
          errors(i)%absolute = 0d0
          errors(i)%relative = 0d0
        end do
      end if
    else if (trim(setting%init_type) == 'local_alpha_delta' .or. &
      trim(setting%init_type) == 'local_alpha_delta_group') then
      call terminate('init type not supported in this version', 1)
      !call make_localized_indices(setting, dim, group_id, atom_indices, localized_indices)
      !if (check_master()) then
      !  call print_localized_indices(dim, localized_indices)
      !end if
      !call blacs_gridinfo(proc%context, nprow, npcol, myrow, mycol)
      !call setup_distributed_matrix_complex('Y', proc, dim, dim, &
      !     Y_localized_desc, Y_localized, .true.)
      !do i = 1, dim
      !  if (.not. localized_indices(i)) then
      !    if (myrow == indxg2p(i, H_desc(block_row_), 0, H_desc(rsrc_), nprow) .and. &
      !         mycol == indxg2p(i, H_desc(block_col_), 0, H_desc(csrc_), npcol)) then
      !      call pzelget('Self', ' ', diag_tmp, H, i, i, H_desc)
      !      call pzelset(H, i, i, H_desc, diag_tmp + setting%localize_potential_depth)
      !    end if
      !  end if
      !end do
      !call solve_gevp(dim, 1, proc, H, H_desc, S, S_desc, eigenvalues, Y_localized, Y_localized_desc)
      !do i = 1, dim
      !  if (.not. localized_indices(i)) then
      !    if (myrow == indxg2p(i, H_desc(block_row_), 0, H_desc(rsrc_), nprow) .and. &
      !         mycol == indxg2p(i, H_desc(block_col_), 0, H_desc(csrc_), npcol)) then
      !      call pzelget('Self', ' ', diag_tmp, H, i, i, H_desc)
      !      call pzelset(H, i, i, H_desc, diag_tmp - setting%localize_potential_depth)
      !    end if
      !  end if
      !end do
      !actual_alpha_delta_index = get_actual_alpha_delta_index(setting, &
      !     dim, localized_indices, setting%alpha_delta_index)
      !call pzcopy(dim, Y_localized, 1, actual_alpha_delta_index, Y_localized_desc, 1, &
      !     full_vecs, 1, col_psi2, full_vecs_desc, 1)
      !call reconcile_from_lcao_coef(setting, 0d0, &
      !     S, S_desc, Y_filtered, Y_filtered_desc, &
      !     full_vecs, full_vecs_desc, filtered_vecs, filtered_vecs_desc)
      !call get_filtering_errors(dim, full_vecs, full_vecs_desc, absolute_filter_error, relative_filter_error)
    else if (trim(setting%init_type) == 'lcao_file') then
      call terminate('not implemented', 49)
      !if (check_master()) then
      !  call read_vector(dim, trim(setting%lcao_filename), dv_psi)
      !end if
      !call mpi_bcast(dv_psi, dim, mpi_double_complex, g_wk_master_pnum, mpi_comm_world, ierr)
      !!if (setting%to_multiply_phase_factor) then
      !!  call multiply_phase_factors(dim, num_atoms, atom_indices, atom_coordinates, &
      !!       setting%phase_factor_coef, psi)
      !!end if
      !call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, dv_psi(:, 1), dv_alpha)
      !call normalize_vector(setting%num_filter, dv_alpha(:, 1))
      !call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_alpha(:, 1), dv_psi(:, 1))
    end if
  end subroutine set_initial_value


  subroutine compute_charges(dim, structure, S_sparse, dv_psi, dv_charge_on_basis, dv_charge_on_atoms, charge_moment)
    integer, intent(in) :: dim
    type(wk_structure_t), intent(in) :: structure
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(dim)
    real(8), intent(out) :: dv_charge_on_basis(:), dv_charge_on_atoms(:)
    type(wk_charge_moment_t), intent(out) :: charge_moment

    call get_mulliken_charges_on_basis(dim, S_sparse, dv_psi, dv_charge_on_basis)
    call get_mulliken_charges_on_atoms(dim, structure, S_sparse, dv_psi, dv_charge_on_atoms)
    call get_mulliken_charge_coordinate_moments(structure, dv_charge_on_atoms, charge_moment)
  end subroutine compute_charges


  subroutine compute_tightbinding_energy(num_filter, dv_alpha, eigenvalues, energies)
    integer, intent(in) :: num_filter
    real(8), intent(in) :: eigenvalues(num_filter)
    complex(kind(0d0)), intent(in) :: dv_alpha(num_filter)
    ! Because this subroutine does not determine all fields of energies,
    ! intent(in) is needed to save the undetermined fields.
    type(wk_energy_t), intent(inout) :: energies

    integer :: i
    real(8) :: dv_alpha_squared(num_filter), alpha_norm
    ! Function.
    real(8) :: dznrm2

    alpha_norm = dznrm2(num_filter, dv_alpha, 1)
    do i = 1, num_filter
      dv_alpha_squared(i) = (abs(dv_alpha(i)) / alpha_norm) ** 2d0
    end do
    energies%tightbinding = 0d0
    energies%tightbinding_deviation = 0d0
    do i = 1, num_filter
      energies%tightbinding = energies%tightbinding + dv_alpha_squared(i) * eigenvalues(i)
    end do
    do i = 1, num_filter
      energies%tightbinding_deviation = energies%tightbinding_deviation + &
           dv_alpha_squared(i) * (eigenvalues(i) - energies%tightbinding) ** 2d0
    end do
    energies%tightbinding_deviation = sqrt(energies%tightbinding_deviation)
  end subroutine compute_tightbinding_energy


  subroutine compute_energies(setting, structure, &
       H_sparse, S_sparse, basis, &
       dv_charge_on_basis, dv_charge_on_atoms, charge_factor, &
       dv_psi, dv_alpha, &
       dv_atom_perturb, &
       H1, H1_desc, energies)
    type(wk_setting_t), intent(in) :: setting
    type(wk_structure_t), intent(in) :: structure
    type(sparse_mat), intent(in) :: H_sparse, S_sparse
    type(wk_basis_t), intent(in) :: basis
    integer, intent(in) :: H1_desc(desc_size)
    type(wk_charge_factor_t), intent(in) :: charge_factor
    real(8), intent(in) :: dv_charge_on_basis(:), dv_charge_on_atoms(structure%num_atoms)
    complex(kind(0d0)), intent(in) :: dv_psi(:), dv_alpha(:)
    real(8), intent(inout) :: dv_atom_perturb(structure%num_atoms)
    real(8), intent(out) :: H1(:, :)
    type(wk_energy_t), intent(out) :: energies

    integer :: dim
    complex(kind(0d0)) :: energy_tmp
    type(sparse_mat) :: H1_lcao_sparse_charge_overlap
    ! Functions.
    complex(kind(0d0)) :: zdotc

    call compute_tightbinding_energy(basis%num_basis, dv_alpha, basis%dv_eigenvalues, energies)

    call make_H1(setting%h1_type, structure, &
         S_sparse, basis, &
         .true., &  ! H and H_desc are not referenced because is_init is true.
         setting%is_restart_mode, &
         0d0, setting%temperature, setting%delta_t, setting%perturb_interval, &
         dv_charge_on_basis, dv_charge_on_atoms, charge_factor, &
         dv_atom_perturb, &
         H1, H1_desc)
    if (trim(setting%filter_mode) == 'group') then
      H1(:, :) = H1(:, :) + basis%H1_base(:, :)  ! Sum of distributed matrices.
    end if

    if (trim(setting%h1_type) == 'charge_overlap') then
      call get_charge_overlap_energy(structure, charge_factor, dv_charge_on_atoms, energies%nonlinear)
    else
      call get_A_inner_product(basis%num_basis, H1, H1_desc, dv_alpha, dv_alpha, energy_tmp)
      energies%nonlinear = truncate_imag(energy_tmp)
      if (trim(setting%h1_type) == 'charge') then
        energies%nonlinear = energies%nonlinear / 2d0
      end if
    end if
    energies%total = energies%tightbinding + energies%nonlinear
  end subroutine compute_energies


  subroutine initialize(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state
    ! intent(out) parameters in state are:
    ! H1(:, :), dv_psi(dim), dv_alpha(Y_filtered_desc(cols_))
    ! dv_charge_on_basis(dim), dv_charge_on_atoms(structure%num_atoms)
    ! charge_moment, energies, errors

    integer :: i

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] initialize() : start, set initial value'
    end if

    call set_initial_value(setting, state%dim, state%structure, state%group_id, &
         state%H_sparse, state%S_sparse, state%basis, &
         state%dv_psi, state%dv_alpha, state%errors)

    if (setting%is_restart_mode) then
      if (check_master()) then
        write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
             '] initialize() : read restarting atom speed or atom perturb'
      end if
      if (trim(setting%h1_type) == 'maxwell') then
        !do i = 1, num_atoms
        !  call pzelset(full_vecs, i, col_atom_speed, full_vecs_desc, &
        !       cmplx(setting%restart_atom_speed(i), 0d0, kind(0d0)))
        !end do
      else if (trim(setting%h1_type) == 'harmonic') then
        !do i = 1, num_atoms
        !  call pzelset(full_vecs, i, col_atom_perturb, full_vecs_desc, &
        !       cmplx(setting%restart_atom_perturb(i), 0d0, kind(0d0)))
        !end do
      end if
    end if

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] initialize() : compute charges form initial value'
    end if

    do i = 1, setting%num_multiple_initials
      call compute_charges(state%dim, state%structure, state%S_sparse, &
           state%dv_psi(:, i), state%dv_charge_on_basis(:, i), state%dv_charge_on_atoms(:, i), &
           state%charge_moment(i))
    end do

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] initialize() : compute energies form initial value'
    end if

    do i = 1, setting%num_multiple_initials
      call compute_energies(setting, state%structure, &
           state%H_sparse, state%S_sparse, state%basis, &
           state%dv_charge_on_basis(:, i), state%dv_charge_on_atoms(:, i), state%charge_factor, &
           state%dv_psi(:, i), state%dv_alpha(:, i), &
           state%dv_atom_perturb, &
           state%H1, state%basis%H1_desc, state%energies(i))
    end do
  end subroutine initialize


  subroutine re_initialize_state(setting, dim, t, structure, charge_factor, is_first_in_multiple_initials, &
       H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
       basis, YSY_filtered, YSY_filtered_desc, &
       dv_eigenvalues_prev, H1, H1_desc, &
       dv_psi, dv_alpha, dv_psi_reconcile, dv_alpha_reconcile, &
       dv_charge_on_basis, dv_charge_on_atoms, &
       dv_atom_perturb, &
       charge_moment, energies, errors)
    type(wk_setting_t), intent(in) :: setting
    type(wk_structure_t), intent(in) :: structure
    type(wk_charge_factor_t), intent(in) :: charge_factor
    logical, intent(in) :: is_first_in_multiple_initials
    integer, intent(in) :: dim
    real(8), intent(in) :: t, dv_eigenvalues_prev(setting%num_filter)!, dv_eigenvalues(setting%num_filter)
    type(sparse_mat), intent(in) :: H_sparse, S_sparse, H_sparse_prev, S_sparse_prev
    type(wk_basis_t), intent(in) :: basis
    integer, intent(in) :: YSY_filtered_desc(desc_size), H1_desc(desc_size)
    real(8), intent(in) :: YSY_filtered(:, :)
    complex(kind(0d0)), intent(in) :: dv_psi(:), dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_psi_reconcile(:), dv_alpha_reconcile(:)
    real(8), intent(out) :: H1(:, :)
    real(8), intent(out) :: dv_charge_on_basis(dim), dv_charge_on_atoms(structure%num_atoms)
    real(8), intent(inout) :: dv_atom_perturb(structure%num_atoms)
    type(wk_charge_moment_t), intent(out) :: charge_moment
    type(wk_energy_t), intent(out) :: energies
    type(wk_error_t), intent(out) :: errors

    complex(kind(0d0)) :: dv_psi_evol(dim), dv_alpha_evol(setting%num_filter)
    if (.true.) then
      dv_psi_evol = dv_psi
      dv_alpha_evol = dv_alpha
    else
      call matvec_time_evolution_by_matrix_replace(setting%delta_t, &
           H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
           dv_psi, dv_psi_evol)
      call lcao_coef_to_alpha(S_sparse, basis, dv_psi_evol, dv_alpha_evol)
    end if

    !if (trim(setting%re_initialize_method) == 'minimize_lcao_error') then
    !  call reconcile_from_lcao_coef(setting%num_filter, &
    !       S_sparse, Y_filtered, Y_filtered_desc, &
    !       dv_psi_evol, dv_alpha_reconcile, dv_psi_reconcile)
    !else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_cutoff') then
    !  call reconcile_from_lcao_coef_cutoff(setting%num_filter, setting%vector_cutoff_residual, &
    !       S_sparse, Y_filtered, Y_filtered_desc, &
    !       dv_psi_evol, dv_alpha_reconcile, dv_psi_reconcile)
    !else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_suppress') then
    !  call reconcile_from_lcao_coef_suppress(setting%num_filter, setting%suppress_constant, &
    !       S_sparse, dv_eigenvalues, Y_filtered, Y_filtered_desc, &
    !       dv_psi_evol, dv_alpha_reconcile, dv_psi_reconcile)
    if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress') then
      call reconcile_from_alpha_matrix_suppress(t, setting%suppress_constant, &
           dv_eigenvalues_prev, basis%dv_eigenvalues, basis, YSY_filtered, YSY_filtered_desc, &
           dv_alpha_evol, dv_alpha_reconcile, dv_psi_reconcile)
    !else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_orthogonal') then
    !  call reconcile_from_alpha_matrix_suppress_orthogonal(setting%num_filter, t, setting%suppress_constant, &
    !       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
    !       dv_alpha_evol, dv_alpha_reconcile, dv_psi_reconcile)
    !else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_adaptive') then
    !  call reconcile_from_alpha_matrix_suppress_adaptive(setting%num_filter, t, setting%suppress_constant, &
    !       is_first_in_multiple_initials, &
    !       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
    !       dv_alpha_evol, dv_alpha_reconcile, dv_psi_reconcile)
    !else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_select') then
    !  call reconcile_from_alpha_matrix_suppress_select(setting%num_filter, t, setting%suppress_constant, &
    !       is_first_in_multiple_initials, &
    !       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
    !       dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_alpha_error') then
      dv_alpha_reconcile(:) = dv_alpha(:)
      call alpha_to_lcao_coef(basis, dv_alpha_reconcile, dv_psi_reconcile)
    else
      call terminate('re_initialize_state: unknown re-initialization method', 1)
    end if
    call get_filtering_errors(dim, dv_psi, dv_psi_reconcile, errors)
    ! backed up col_filtered_work1 is used here.
    !call get_re_initialize_errors(setting%num_filter, filtered_vecs, filtered_vecs_desc, &
    !     absolute_alpha_error, relative_alpha_error)

    if (check_master()) then
      write (0, '(A, F16.6, A, E26.16e3, A, E26.16e3, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] alpha error after re-initialization is ', &
           errors%absolute, ' (absolute) and ', errors%relative, ' (relative)'
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] re_initialize_from_lcao_coef() : compute charges form initial value'
    end if

    call compute_charges(dim, structure, S_sparse, dv_psi_reconcile, &
         dv_charge_on_basis, dv_charge_on_atoms, charge_moment)

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] re_initialize_from_lcao_coef() : compute energies form initial value'
    end if

    call compute_energies(setting, structure, &
       H_sparse, S_sparse, basis, &
       dv_charge_on_basis, dv_charge_on_atoms, charge_factor, &
       dv_psi_reconcile, dv_alpha_reconcile, &
       dv_atom_perturb, &
       H1, H1_desc, energies)
  end subroutine re_initialize_state


  subroutine clear_offdiag_blocks(filter_group_indices, X_sparse)
    integer, intent(in) :: filter_group_indices(:, :)
    type(sparse_mat), intent(inout) :: X_sparse

    integer :: num_groups, dim, g, i_start, j_start, j_end, i, j, k

    num_groups = size(filter_group_indices, 2) - 1
    !print *, filter_group_indices, num_groups
    dim = filter_group_indices(1, num_groups + 1) - 1
    do g = 1, num_groups - 1
      ! These indices are four corners of sub-diagonal blocks below the g-th diagonal block.
      i_start = filter_group_indices(1, g + 1)
      j_start = filter_group_indices(1, g)
      j_end = i_start - 1
      do k = 1, X_sparse%num_non_zeros
        i = X_sparse%suffix(1, k)
        j = X_sparse%suffix(2, k)
        if (i_start <= i .and. i <= dim .and. j_start <= j .and. j <= j_end) then
          X_sparse%value(k) = 0d0
        else if (i_start <= j .and. j <= dim .and. j_start <= i .and. i <= j_end) then
          X_sparse%value(k) = 0d0
        end if
      end do
    end do
  end subroutine clear_offdiag_blocks
end module wk_initialization_m
