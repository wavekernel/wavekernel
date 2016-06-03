module wp_initialization_m
  use mpi
  use wp_descriptor_parameters_m
  use wp_distribute_matrix_m
  use wp_processes_m
  use wp_atom_m
  use wp_charge_m
  use wp_conversion_m
  use wp_global_variables_m
  use wp_linear_algebra_m
  use wp_matrix_generation_m
  use wp_setting_m
  use wp_state_m
  use wp_util_m
  implicit none

  private
  public :: initialize, re_initialize_state, clear_offdiag_blocks

contains

  subroutine print_eigenvectors(Y, Y_desc)
    complex(kind(0d0)) :: Y(:, :)
    integer :: Y_desc(desc_size)

    real(8) :: work(Y_desc(rows_))
    real(8), allocatable :: Y_real(:, :)
    integer :: len, i, digit, stat, err
    character(512) :: num_str, filename
    integer, parameter :: iunit = 10, max_num_digits = 8

    allocate(Y_real(size(Y, 1), size(Y, 2)))
    Y_real(:, :) = real(Y(:, :), kind(0d0))

    do i = 1, Y_desc(cols_)
      if (check_master()) then
        write (num_str, '(i0)') i
        len = len_trim(num_str)
        num_str(max_num_digits - len + 1 : max_num_digits) = num_str(1 : len)
        do digit = 1, max_num_digits - len
          num_str(digit:digit) = '0'
        end do

        filename = trim(num_str) // '.dat'
        open (iunit, file=trim(filename), status='replace', iostat=stat)
      end if

      call mpi_bcast(stat, 1, mpi_integer, 0, mpi_comm_world, err) ! Share stat for file opening
      if (stat /= 0) then
        if (check_master ()) print *, 'iostat: ', stat
        call terminate('print_eigenvectors: cannot open ' // trim(filename), stat)
      end if

      call mpi_barrier(mpi_comm_world, stat)
      call eigentest_pdlaprnt(Y_desc(rows_), 1, Y_real, &
           1, i, Y_desc, 0, 0, iunit, work)

      if (check_master()) then
        close (iunit)
      end if
    end do
  end subroutine print_eigenvectors


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
    type(wp_setting_t), intent(in) :: setting
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
         ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, '] initialize() : localized indices are '
    do i = 1, dim
      if (localized_indices(i)) then
        write(0, '(I0, " ")', advance='no') i
      end if
    end do
    write(0, *)
  end subroutine print_localized_indices


  integer function get_actual_alpha_delta_index(setting, dim, localized_indices, &
       alpha_delta_index)
    type(wp_setting_t), intent(in) :: setting
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
      write (0, '(A, F16.6, A, I0)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
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
    type(wp_error_t), intent(out) :: errors

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


  subroutine reconcile_from_lcao_coef(num_filter, t, S_sparse, eigenvalues, Y_filtered, Y_filtered_desc, &
       dv_psi, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size)
    real(8), intent(in) :: t, eigenvalues(:), Y_filtered(:, :)
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    complex(kind(0d0)) :: dv_alpha(num_filter)

    call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, eigenvalues, t, dv_psi, dv_alpha_reconcile)
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
  end subroutine reconcile_from_lcao_coef


  subroutine reconcile_from_lcao_coef_cutoff(num_filter, t, cutoff, &
       S_sparse, eigenvalues, Y_filtered, Y_filtered_desc, &
       dv_psi, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size)
    real(8), intent(in) :: t, cutoff, eigenvalues(:), Y_filtered(:, :)
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    complex(kind(0d0)) :: dv_alpha(num_filter)

    call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, eigenvalues, t, dv_psi, dv_alpha)
    call cutoff_vector(num_filter, cutoff, dv_alpha, dv_alpha_reconcile)
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
  end subroutine reconcile_from_lcao_coef_cutoff


  subroutine reconcile_from_lcao_coef_suppress(num_filter, t, suppress_constant, &
       S_sparse, eigenvalues, Y_filtered, Y_filtered_desc, &
       dv_psi, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size)
    real(8), intent(in) :: t, suppress_constant, eigenvalues(:), Y_filtered(:, :)
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    integer :: i, min_abs_i, ierr
    complex(kind(0d0)) :: dv_alpha(num_filter), diff
    real(8) :: min_abs

    call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, eigenvalues, t, dv_psi, dv_alpha)
    i = 0
    min_abs = -1d0
    do i = 1, num_filter
      if (abs(dv_alpha(i)) > min_abs) then
        min_abs = abs(dv_alpha(i))
        min_abs_i = i
      end if
    end do
    do i = 1, num_filter
      diff = eigenvalues(i) - eigenvalues(min_abs_i)
      dv_alpha_reconcile(i) = dv_alpha(i) * exp(- abs(diff * suppress_constant) ** 2d0)
    end do
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
  end subroutine reconcile_from_lcao_coef_suppress


  subroutine reconcile_from_alpha_matrix_suppress_old(num_filter, t, suppress_constant, &
       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
       dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
    real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
    real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
    complex(kind(0d0)), intent(in) :: dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc, ks(num_filter)
    real(8) :: dv_suppress_factor(num_filter), suppress_factor_sum, dv_suppress_factor_copy(num_filter)
    complex(kind(0d0)) :: dv_evcoef(num_filter), dv_evcoef_reconcile(num_filter)
    real(8), allocatable :: YSY_filtered_suppress(:, :), B(:, :)
    real(8) :: work(1000)
    integer, save :: count = 0
    character(len=50) :: filename
    character(len=6) :: count_str
    character(len=12) :: time_str
    integer, parameter :: iunit_YSY = 20
    real(8) :: dznrm2
    print *, 'ZZZZZstart', count, t * 2.418884326505e-5
    print *, 'ZZZZZdiffeigenval', abs(dv_eigenvalues(:) - dv_eigenvalues_prev(:))
    call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
    allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
    allocate(B(size(YSY_filtered, 1), size(YSY_filtered, 2)))
    YSY_filtered_suppress(:, :) = 0d0
    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues_prev, t, dv_alpha, dv_evcoef)
    do j = 1, num_filter
      do i = 1, num_filter
        dv_suppress_factor(i) = exp(- suppress_constant * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
        !dv_suppress_factor(i) = exp(- suppress_constant * (dv_eigenvalues(i) - dv_eigenvalues(j)) ** 2d0)
      end do
      suppress_factor_sum = sum(dv_suppress_factor)
      dv_suppress_factor_copy(:) = dv_suppress_factor(:)
      if (check_master()) then
        call comb_sort(num_filter, dv_suppress_factor_copy, ks)
        write (0, '(A, F16.6, A, F16.6, I10, 4F16.6)') ' [Event', &
             mpi_wtime() - g_wp_mpi_wtime_init, &
             '] reconcile_from_alpha_matrix_suppress() : t, j, min, 2ndmax, max, sum of dv_suppress_factor ', &
             t, j, dv_suppress_factor_copy(num_filter), &
             dv_suppress_factor_copy(2), dv_suppress_factor_copy(1), suppress_factor_sum
        !if (j /= ks(1)) then
        !  write (0, '(A, F16.6, A, I10, I0)') ' [Event', &
        !       mpi_wtime() - g_wp_mpi_wtime_init, &
        !       '] reconcile_from_alpha_matrix_suppress() : swap occured between ', j, ks(1)
        !end if
      end if
      call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor', dv_suppress_factor)
      dv_suppress_factor(:) = dv_suppress_factor(:) !/ suppress_factor_sum
      print *, 'ZZZZZZZsuppressfactor j, suppress_factor_sum, all: ', j, suppress_factor_sum, ',', dv_suppress_factor(:)
      do i = 1, num_filter
        call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
        if (myrow == rsrc .and. mycol == csrc) then
          YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
          if (abs(YSY_filtered_suppress(i_local, j_local)) > 5d-3) then
            print *, 'ZZZZZZZYSYsuppressed j, i, YSYs, YSY, fac: ', j, ' -> ', i, YSY_filtered_suppress(i_local, j_local), ', ', &
                 YSY_filtered(i_local, j_local), dv_suppress_factor(i)
          end if
        end if
      end do
    end do

    write(count_str, '(I6.6)') count
    write(time_str, '(F0.4)') t * 2.418884326505e-5  ! kPsecPerAu.
    filename = 'YSY_' // count_str // '_' // trim(time_str) // 'ps.mtx'
    if (check_master()) then
      open(iunit_YSY, file=trim(filename))
    end if
    call pdlaprnt(num_filter, num_filter, YSY_filtered, 1, 1, YSY_filtered_desc, 0, 0, 'YSY', iunit_YSY, work)
    if (check_master()) then
      close(iunit_YSY)
    end if
    filename = 'YSYsuppress_' // count_str // '_' // trim(time_str) // 'ps.mtx'
    if (check_master()) then
      open(iunit_YSY, file=trim(filename))
    end if
    call pdlaprnt(num_filter, num_filter, YSY_filtered_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'YSYsuppress', &
         iunit_YSY, work)
    if (check_master()) then
      close(iunit_YSY)
    end if

    dv_evcoef_reconcile(:) = kZero
    !stop
    !print *, 'ZZZZ1', t * 2.418884326505e-5, 'ps', dv_evcoef
    !call matvec_nearest_orthonormal_matrix(YSY_filtered_suppress, YSY_filtered_desc, &
    !     dv_evcoef, dv_evcoef_reconcile, B)
    !filename = 'YSYB_' // count_str // '_' // trim(time_str) // 'ps.mtx'
    !if (check_master()) then
    !  open(iunit_YSY, file=trim(filename))
    !end if
    !call pdlaprnt(num_filter, num_filter, B, 1, 1, YSY_filtered_desc, 0, 0, 'YSYB', &
    !     iunit_YSY, work)
    !if (check_master()) then
    !  close(iunit_YSY)
    !end if
    count = count + 1
    call matvec_dd_z('No', YSY_filtered_suppress, YSY_filtered_desc, kOne, dv_evcoef, kZero, dv_evcoef_reconcile)
    !print *, 'ZZZZ2', t * 2.418884326505e-5, 'ps', dv_evcoef_reconcile
    !stop
    !if (t * 2.418884326505e-5 > 0.002) then
    !  stop
    !end if
    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues, -t, dv_evcoef_reconcile, dv_alpha_reconcile)
    print *, 'ZZZZZZalphanorm', dznrm2(num_filter, dv_alpha_reconcile, 1)
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
    deallocate(YSY_filtered_suppress)
  end subroutine reconcile_from_alpha_matrix_suppress_old


  subroutine reconcile_from_alpha_matrix_suppress(num_filter, t, suppress_constant, &
       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
       dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
    real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
    real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
    complex(kind(0d0)), intent(in) :: dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc, ks(num_filter)
    real(8) :: dv_suppress_factor(num_filter), suppress_factor_sum, dv_suppress_factor_copy(num_filter)
    real(8) :: dv_suppress_factor2(num_filter), suppress_factor_sum2
    complex(kind(0d0)) :: dv_evcoef(num_filter), dv_evcoef_reconcile(num_filter)
    complex(kind(0d0)) :: dv_evcoef_amplitude(num_filter), dv_evcoef_amplitude_reconcile(num_filter)
    real(8), allocatable :: ENE_suppress(:, :), YSY_filtered_suppress(:, :)
    real(8) :: work(1000), energy_normalizer
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

    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues_prev, t, dv_alpha, dv_evcoef)
    do j = 1, num_filter
      do i = 1, num_filter
        dv_suppress_factor(i) = exp(- suppress_constant * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
      end do
      call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor', dv_suppress_factor)
      do i = 1, num_filter
        call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
        if (myrow == rsrc .and. mycol == csrc) then
          YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
        end if
      end do
    end do

    !write(count_str, '(I6.6)') count
    !write(time_str, '(F0.4)') t * 2.418884326505e-5  ! kPsecPerAu.
    !filename = 'YSYsuppressed_' // count_str // '_' // trim(time_str) // 'ps.mtx'
    !if (check_master()) then
    !  open(iunit_YSY, file=trim(filename))
    !end if
    !call pdlaprnt(num_filter, num_filter, YSY_filtered_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'YSY', iunit_YSY, work)
    !if (check_master()) then
    !  close(iunit_YSY)
    !end if

    dv_evcoef_reconcile(:) = kZero
    print *, '------------------ ZZZZstart', t * 2.418884326505d-5, 'ps ------------------'
    call matvec_dd_z2('No', YSY_filtered_suppress, YSY_filtered_desc, kOne, dv_evcoef, kZero, dv_evcoef_reconcile)
    !print *, 'Energy correction X', dv_evcoef
    !print *, 'Energy correction Y', dv_evcoef_reconcile
    !
    !do j = 1, num_filter
    !  dv_evcoef_amplitude(j) = conjg(dv_evcoef(j)) * dv_evcoef(j)
    !  do i = 1, num_filter
    !    if (abs(dv_evcoef_reconcile(i)) < 1d-10) then  ! αの値ではなくsmoothed deltaの方で切り落としを定める.
    !      !print *, 'ZZZZ i: alphaprime0(zero), abs(alphaprime0)', i, ':', dv_evcoef_reconcile(i), abs(dv_evcoef_reconcile(i))
    !      dv_suppress_factor2(i) = kZero
    !    else
    !      dv_suppress_factor2(i) = exp(- 1d4 * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
    !    end if
    !  end do
    !  call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor2', dv_suppress_factor2)
    !  do i = 1, num_filter
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
    !call pdlaprnt(num_filter, num_filter, ENE_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'ENEsuppress', &
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
    !do i = 1, num_filter
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

    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues, -t, dv_evcoef_reconcile, dv_alpha_reconcile)
    if (check_master()) then
      write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] reconcile_from_alpha_matrix_suppress() : evcoef norm before normalization ', &
           dznrm2(num_filter, dv_alpha_reconcile, 1)
    end if
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
    deallocate(YSY_filtered_suppress, ENE_suppress)
  end subroutine reconcile_from_alpha_matrix_suppress


  subroutine reconcile_from_alpha_matrix_suppress_orthogonal(num_filter, t, suppress_constant, &
       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
       dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
    real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
    real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
    complex(kind(0d0)), intent(in) :: dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc
    real(8) :: dv_suppress_factor(num_filter)
    complex(kind(0d0)) :: dv_evcoef(num_filter), dv_evcoef_reconcile(num_filter)
    real(8), allocatable :: YSY_filtered_orthogonal(:, :), YSY_filtered_suppress(:, :)
    real(8) :: dznrm2

    call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
    allocate(YSY_filtered_orthogonal(size(YSY_filtered, 1), size(YSY_filtered, 2)))
    allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
    YSY_filtered_orthogonal(:, :) = 0d0

    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues_prev, t, dv_alpha, dv_evcoef)
    do j = 1, num_filter
      do i = 1, num_filter
        dv_suppress_factor(i) = exp(- suppress_constant * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
      end do
      call check_nan_vector('reconcile_from_alpha_matrix_suppress_orthogonal dv_suppress_factor', dv_suppress_factor)
      do i = 1, num_filter
        call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
        if (myrow == rsrc .and. mycol == csrc) then
          YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
        end if
      end do
    end do

    dv_evcoef_reconcile(:) = kZero
    call matvec_nearest_orthonormal_matrix(YSY_filtered_suppress, YSY_filtered_desc, &
         dv_evcoef, dv_evcoef_reconcile, YSY_filtered_orthogonal)

    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues, -t, dv_evcoef_reconcile, dv_alpha_reconcile)
    if (check_master()) then
      write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] reconcile_from_alpha_matrix_suppress_orthogonal() : evcoef norm before normalization ', &
           dznrm2(num_filter, dv_alpha_reconcile, 1)
    end if
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
    deallocate(YSY_filtered_suppress, YSY_filtered_orthogonal)
  end subroutine reconcile_from_alpha_matrix_suppress_orthogonal


  subroutine reconcile_from_alpha_matrix_suppress_adaptive(num_filter, t, suppress_constant, &
       dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
       dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
    integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
    real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
    real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
    complex(kind(0d0)), intent(in) :: dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)

    integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc
    real(8) :: dv_suppress_factor(num_filter)
    complex(kind(0d0)) :: dv_evcoef(num_filter), dv_evcoef_reconcile(num_filter)
    real(8), allocatable :: YSY_filtered_suppress(:, :)
    type(wp_energy_t) :: energies

    real(8) :: suppress_exponent
    real(8), parameter :: relax_constant_for_suppression = 1d-3
    ! Function.
    real(8) :: dznrm2

    call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
    allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))

    call compute_tightbinding_energy(num_filter, dv_alpha, dv_eigenvalues, energies)

    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues_prev, t, dv_alpha, dv_evcoef)
    do j = 1, num_filter
      do i = 1, num_filter
        suppress_exponent = suppress_constant / (energies%tightbinding_deviation + relax_constant_for_suppression)
        dv_suppress_factor(i) = exp(- suppress_exponent * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
      end do
      call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor', dv_suppress_factor)
      do i = 1, num_filter
        call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
        if (myrow == rsrc .and. mycol == csrc) then
          YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
        end if
      end do
    end do

    dv_evcoef_reconcile(:) = kZero
    call matvec_dd_z2('No', YSY_filtered_suppress, YSY_filtered_desc, kOne, dv_evcoef, kZero, dv_evcoef_reconcile)

    call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues, -t, dv_evcoef_reconcile, dv_alpha_reconcile)
    if (check_master()) then
      write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] reconcile_from_alpha_matrix_suppress() : evcoef norm before normalization ', &
           dznrm2(num_filter, dv_alpha_reconcile, 1)
    end if
    call normalize_vector(num_filter, dv_alpha_reconcile)
    call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
    deallocate(YSY_filtered_suppress)
  end subroutine reconcile_from_alpha_matrix_suppress_adaptive


  !subroutine reconcile_from_alpha_matrix_suppress(num_filter, t, suppress_constant, &
  !     dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
  !     dv_alpha, dv_alpha_reconcile, dv_psi_reconcile)
  !  integer, intent(in) :: num_filter, Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
  !  real(8), intent(in) :: t, suppress_constant, dv_eigenvalues_prev(:), dv_eigenvalues(:)
  !  real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :)
  !  complex(kind(0d0)), intent(in) :: dv_alpha(:)
  !  complex(kind(0d0)), intent(out) :: dv_alpha_reconcile(:), dv_psi_reconcile(:)
  !
  !  integer :: i, j, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc, ks(num_filter)
  !  real(8) :: dv_suppress_factor(num_filter), suppress_factor_sum, dv_suppress_factor_copy(num_filter)
  !  real(8) :: dv_suppress_factor2(num_filter), suppress_factor_sum2
  !  complex(kind(0d0)) :: dv_evcoef(num_filter), dv_evcoef_reconcile(num_filter)
  !  complex(kind(0d0)) :: dv_evcoef_amplitude(num_filter), dv_evcoef_amplitude_reconcile(num_filter)
  !  real(8), allocatable :: ENE_suppress(:, :), YSY_filtered_suppress(:, :)
  !  real(8) :: work(1000), energy_normalizer
  !  integer, save :: count = 0
  !  character(len=50) :: filename
  !  character(len=6) :: count_str
  !  character(len=12) :: time_str
  !  integer, parameter :: iunit_YSY = 20
  !  real(8) :: dznrm2
  !
  !  print *, 'ZZZZZstart', count, t * 2.418884326505e-5
  !
  !  call blacs_gridinfo(YSY_filtered_desc(context_), nprow, npcol, myrow, mycol)
  !  !allocate(ENE_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !  allocate(YSY_filtered_suppress(size(YSY_filtered, 1), size(YSY_filtered, 2)))
  !  !ENE_suppress(:, :) = 0d0
  !
  !  call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues_prev, t, dv_alpha, dv_evcoef)
  !  do j = 1, num_filter
  !    dv_evcoef_amplitude(j) = conjg(dv_evcoef(j)) * dv_evcoef(j)
  !    do i = 1, num_filter
  !      dv_suppress_factor(i) = exp(- suppress_constant * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
  !      !dv_suppress_factor2(i) = exp(- 1d4 * (dv_eigenvalues(i) - dv_eigenvalues_prev(j)) ** 2d0)
  !    end do
  !    !suppress_factor_sum = sum(dv_suppress_factor)
  !    !suppress_factor_sum2 = sum(dv_suppress_factor2)
  !    dv_suppress_factor(:) = dv_suppress_factor(:) !/ suppress_factor_sum
  !    !dv_suppress_factor2(:) = dv_suppress_factor2(:) / suppress_factor_sum2
  !    call check_nan_vector('reconcile_from_alpha_matrix_suppress dv_suppress_factor', dv_suppress_factor)
  !    do i = 1, num_filter
  !      call infog2l(i, j, YSY_filtered_desc, nprow, npcol, myrow, mycol, i_local, j_local, rsrc, csrc)
  !      if (myrow == rsrc .and. mycol == csrc) then
  !        YSY_filtered_suppress(i_local, j_local) = YSY_filtered(i_local, j_local) * dv_suppress_factor(i)
  !        !ENE_suppress(i_local, j_local) = dv_suppress_factor2(i)
  !      end if
  !    end do
  !  end do
  !  call scale_columns_to_stochastic_matrix(YSY_filtered_desc, YSY_filtered_suppress)
  !
  !  write(count_str, '(I6.6)') count
  !  write(time_str, '(F0.4)') t * 2.418884326505e-5  ! kPsecPerAu.
  !  filename = 'YSY_' // count_str // '_' // trim(time_str) // 'ps.mtx'
  !  if (check_master()) then
  !    open(iunit_YSY, file=trim(filename))
  !  end if
  !  call pdlaprnt(num_filter, num_filter, YSY_filtered, 1, 1, YSY_filtered_desc, 0, 0, 'YSY', iunit_YSY, work)
  !  if (check_master()) then
  !    close(iunit_YSY)
  !  end if
  !  !filename = 'ENEsuppress_' // count_str // '_' // trim(time_str) // 'ps.mtx'
  !  !if (check_master()) then
  !  !  open(iunit_YSY, file=trim(filename))
  !  !end if
  !  !call pdlaprnt(num_filter, num_filter, ENE_suppress, 1, 1, YSY_filtered_desc, 0, 0, 'ENEsuppress', &
  !  !     iunit_YSY, work)
  !  !if (check_master()) then
  !  !  close(iunit_YSY)
  !  !end if
  !
  !  count = count + 1
  !  dv_evcoef_reconcile(:) = kZero
  !  dv_evcoef_amplitude_reconcile(:) = kZero
  !  call matvec_dd_z('No', YSY_filtered_suppress, YSY_filtered_desc, kOne, dv_evcoef, kZero, dv_evcoef_reconcile)
  !  call matvec_dd_z('No', YSY_filtered_suppress, YSY_filtered_desc, &
  !       kOne, dv_evcoef_amplitude, kZero, dv_evcoef_amplitude_reconcile)
  !  print *, 'ZZZ1', dv_evcoef_amplitude
  !  print *, 'ZZZ2', dv_evcoef_amplitude_reconcile
  !  !call matvec_dd_z('No', ENE_suppress, YSY_filtered_desc, kOne, dv_evcoef_amplitude, kZero, dv_evcoef_amplitude_reconcile)
  !  do i = 1, num_filter
  !    if (abs(dv_evcoef_reconcile(i)) < 1d-16) then
  !      print *, 'ZZZZ i: alphaprime0(zero), abs(alphaprime0), targetamp(ignored)', i, ':', &
  !           dv_evcoef_reconcile(i), abs(dv_evcoef_reconcile(i)), dsqrt(dreal(dv_evcoef_amplitude_reconcile(i)))
  !      dv_evcoef_reconcile(i) = kZero
  !    else
  !      energy_normalizer = dsqrt(dreal(dv_evcoef_amplitude_reconcile(i))) / abs(dv_evcoef_reconcile(i))
  !      print *, 'ZZZZ i: alphaprime0, targetamp, normalizer', i, ':', dv_evcoef_reconcile(i), &
  !           dsqrt(dreal(dv_evcoef_amplitude_reconcile(i))), energy_normalizer
  !      dv_evcoef_reconcile(i) = dv_evcoef_reconcile(i) * energy_normalizer
  !    end if
  !  end do
  !
  !  call alpha_to_eigenvector_coef(num_filter, dv_eigenvalues, -t, dv_evcoef_reconcile, dv_alpha_reconcile)
  !  print *, 'ZZZZZZalpha norm before normalization', dznrm2(num_filter, dv_alpha_reconcile, 1)
  !  call normalize_vector(num_filter, dv_alpha_reconcile)
  !  call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
  !  !deallocate(ENE_suppress)
  !end subroutine reconcile_from_alpha_matrix_suppress


  subroutine set_initial_value(setting, proc, dim, structure, group_id, &
       H_sparse, S_sparse, Y_filtered, Y_filtered_desc, &
       dv_psi, dv_alpha, errors)
    type(wp_setting_t), intent(in) :: setting
    type(wp_process_t), intent(in) :: proc
    integer, intent(in) :: dim, group_id(:, :)
    type(sparse_mat), intent(in) :: H_sparse, S_sparse
    type(wp_structure_t), intent(in) :: structure
    integer, intent(in) :: Y_filtered_desc(desc_size)
    real(8), intent(in) :: Y_filtered(:, :)
    complex(kind(0d0)), intent(out) :: dv_psi(:), dv_alpha(:)
    type(wp_error_t), intent(out) :: errors

    integer :: nprow, npcol, myrow, mycol
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
        call make_initial_psi_delta(dim, actual_alpha_delta_index - setting%fst_filter + 1, dv_alpha)
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
      ! \alpha のノルムはプログラムがうまく動いていれば保存されるはずなので規格化は最初の一度のみ行う.
      ! 1 =: \Psi^\dagger S \Psi = c^\dagger Y^\dagger S Y c = ||c||^2 = ||\alpha||^2
      call normalize_vector(setting%num_filter, dv_alpha)
      ! フィルターした後の状態の重ね合わせで \psi を決めているのでこの時点では
      ! \psi に対するエラーは存在しない.
      call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, eigenvalues, 0d0, dv_alpha, dv_psi)
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
        errors%absolute = 0d0
        errors%relative = 0d0
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
      call terminate('initialize: not implemented yet', 2)
      !if (check_master()) then
      !  call read_vector_real(dim, trim(setting%lcao_filename), psi)
      !end if
      !if (setting%to_multiply_phase_factor) then
      !  call multiply_phase_factors(dim, num_atoms, atom_indices, atom_coordinates, &
      !       setting%phase_factor_coef, psi)
      !end if
      !call lcao_coef_to_alpha(S, S_desc, Y, Y_desc, 0d0, full_vecs, full_vecs_desc, &
      !     filtered_vecs, filtered_vecs_desc, col_psi, col_alpha)
      !call normalize_vector(setting%num_filter, filtered_vecs, filtered_vecs_desc, col_alpha)
      !call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, full_vecs, full_vecs_desc, &
      ! filtered_vecs, filtered_vecs_desc, 0d0, col_psi)
    end if
  end subroutine set_initial_value


  subroutine compute_charges(dim, structure, S_sparse, dv_psi, dv_charge_on_basis, dv_charge_on_atoms, charge_moment)
    integer, intent(in) :: dim
    type(wp_structure_t), intent(in) :: structure
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(dim)
    real(8), intent(out) :: dv_charge_on_basis(:), dv_charge_on_atoms(:)
    type(wp_charge_moment_t), intent(out) :: charge_moment

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
    type(wp_energy_t), intent(inout) :: energies

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


  subroutine compute_energies(setting, proc, structure, &
       H_sparse, S_sparse, Y_filtered, Y_filtered_desc, dv_charge_on_atoms, charge_factor, &
       filter_group_indices, Y_local, eigenvalues, dv_psi, dv_alpha, &
       dv_atom_perturb, &
       H1_base, H1, H1_desc, energies)
    type(wp_setting_t), intent(in) :: setting
    type(wp_process_t), intent(in) :: proc
    integer, intent(in) :: filter_group_indices(:, :)
    type(wp_structure_t), intent(in) :: structure
    type(sparse_mat), intent(in) :: H_sparse, S_sparse
    integer, intent(in) :: Y_filtered_desc(desc_size), H1_desc(desc_size)
    type(wp_charge_factor_t), intent(in) :: charge_factor
    real(8), intent(in) :: dv_charge_on_atoms(structure%num_atoms), eigenvalues(:)
    real(8), intent(in) :: Y_filtered(:, :), H1_base(:, :)
    complex(kind(0d0)), intent(in) :: dv_psi(:), dv_alpha(:)
    type(wp_local_matrix_t), intent(in) :: Y_local(:)
    real(8), intent(inout) :: dv_atom_perturb(structure%num_atoms)
    real(8), intent(out) :: H1(:, :)
    type(wp_energy_t), intent(out) :: energies

    integer :: dim, num_filter
    complex(kind(0d0)) :: energy_tmp, dv_h_psi(Y_filtered_desc(rows_)), dv_evcoef(Y_filtered_desc(cols_))
    ! Functions.
    complex(kind(0d0)) :: zdotc

    dim = Y_filtered_desc(rows_)
    num_filter = Y_filtered_desc(cols_)

    call compute_tightbinding_energy(num_filter, dv_alpha, eigenvalues, energies)

    call make_H1(proc, setting%h1_type, structure, &
         S_sparse, Y_filtered, Y_filtered_desc, &
         .true., &  ! H and H_desc are not referenced because is_init is true.
         setting%is_restart_mode, &
         trim(setting%filter_mode) == 'group', filter_group_indices, Y_local, &
         0d0, setting%temperature, setting%delta_t, setting%perturb_interval, &
         dv_charge_on_atoms, charge_factor, &
         dv_atom_perturb, &
         H1, H1_desc)
    if (trim(setting%filter_mode) == 'group') then
      H1(:, :) = H1(:, :) + H1_base(:, :)  ! Sum of distributed matrices.
    end if
    call alpha_to_eigenvector_coef(num_filter, eigenvalues, 0d0, dv_alpha, dv_evcoef)

    if (trim(setting%h1_type) == 'charge_overlap') then
      call get_charge_overlap_energy(structure, charge_factor, dv_charge_on_atoms, energies%nonlinear)
    else
      call get_A_inner_product(setting%num_filter, H1, H1_desc, dv_evcoef, dv_evcoef, energy_tmp)
      energies%nonlinear = truncate_imag(energy_tmp)
      if (trim(setting%h1_type) == 'charge') then
        energies%nonlinear = energies%nonlinear / 2d0
      end if
    end if
    energies%total = energies%tightbinding + energies%nonlinear
  end subroutine compute_energies


  subroutine initialize(setting, proc, state)
    type(wp_setting_t), intent(in) :: setting
    type(wp_process_t), intent(in) :: proc
    type(wp_state_t), intent(inout) :: state
    ! intent(out) parameters in state are:
    ! H1(:, :), dv_psi(dim), dv_alpha(Y_filtered_desc(cols_))
    ! dv_charge_on_basis(dim), dv_charge_on_atoms(structure%num_atoms)
    ! charge_moment, energies, errors

    integer :: i

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] initialize() : start, set initial value'
    end if

    call set_initial_value(setting, proc, state%dim, state%structure, state%group_id, &
         state%H_sparse, state%S_sparse, state%Y_filtered, state%Y_filtered_desc, &
         state%dv_psi, state%dv_alpha, state%errors)

    if (setting%is_restart_mode) then
      if (check_master()) then
        write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
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
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] initialize() : compute charges form initial value'
    end if

    call compute_charges(state%dim, state%structure, state%S_sparse, &
         state%dv_psi, state%dv_charge_on_basis, state%dv_charge_on_atoms, state%charge_moment)

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] initialize() : compute energies form initial value'
    end if

    call compute_energies(setting, proc, state%structure, &
         state%H_sparse, state%S_sparse, state%Y_filtered, state%Y_filtered_desc, &
         state%dv_charge_on_atoms, state%charge_factor, &
         state%filter_group_indices, state%Y_local, state%dv_eigenvalues, state%dv_psi, state%dv_alpha, &
         state%dv_atom_perturb, &
         state%H1_base, state%H1, state%H1_desc, state%energies)
  end subroutine initialize


  subroutine re_initialize_state(setting, proc, dim, t, structure, charge_factor, &
       H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, Y_filtered, Y_filtered_desc, &
       YSY_filtered, YSY_filtered_desc, &
       dv_eigenvalues_prev, dv_eigenvalues, filter_group_indices, Y_local, &
       H1_base, H1, H1_desc, dv_psi, dv_alpha, dv_psi_reconcile, dv_alpha_reconcile, &
       dv_charge_on_basis, dv_charge_on_atoms, &
       dv_atom_perturb, &
       charge_moment, energies, errors)
    type(wp_setting_t), intent(in) :: setting
    type(wp_process_t), intent(in) :: proc
    type(wp_structure_t), intent(in) :: structure
    type(wp_charge_factor_t), intent(in) :: charge_factor
    integer, intent(in) :: dim, filter_group_indices(:, :)
    real(8), intent(in) :: t, dv_eigenvalues_prev(setting%num_filter), dv_eigenvalues(setting%num_filter)
    type(sparse_mat), intent(in) :: H_sparse, S_sparse, H_sparse_prev, S_sparse_prev
    integer, intent(in) :: Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size), H1_desc(desc_size)
    real(8), intent(in) :: Y_filtered(:, :), YSY_filtered(:, :), H1_base(:, :)
    type(wp_local_matrix_t), intent(in) :: Y_local(:)
    complex(kind(0d0)), intent(in) :: dv_psi(:), dv_alpha(:)
    complex(kind(0d0)), intent(out) :: dv_psi_reconcile(:), dv_alpha_reconcile(:)
    real(8), intent(out) :: H1(:, :)
    real(8), intent(out) :: dv_charge_on_basis(dim), dv_charge_on_atoms(structure%num_atoms)
    real(8), intent(inout) :: dv_atom_perturb(structure%num_atoms)
    type(wp_charge_moment_t), intent(out) :: charge_moment
    type(wp_energy_t), intent(out) :: energies
    type(wp_error_t), intent(out) :: errors

    complex(kind(0d0)) :: dv_psi_evol(dim), dv_alpha_evol(setting%num_filter)
    if (.true.) then
      dv_psi_evol = dv_psi
      dv_alpha_evol = dv_alpha
    else
      call matvec_time_evolution_by_matrix_replace(proc, setting%delta_t, &
           H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
           dv_psi, dv_psi_evol)
      call lcao_coef_to_alpha(S_sparse, Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_psi_evol, &
           dv_alpha_evol)
    end if

    if (trim(setting%re_initialize_method) == 'minimize_lcao_error') then
      call reconcile_from_lcao_coef(setting%num_filter, t, &
           S_sparse, dv_eigenvalues, Y_filtered, Y_filtered_desc, &
           dv_psi_evol, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_cutoff') then
      call reconcile_from_lcao_coef_cutoff(setting%num_filter, t, setting%vector_cutoff_residual, &
           S_sparse, dv_eigenvalues, Y_filtered, Y_filtered_desc, &
           dv_psi_evol, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_suppress') then
      call reconcile_from_lcao_coef_suppress(setting%num_filter, t, setting%suppress_constant, &
           S_sparse, dv_eigenvalues, Y_filtered, Y_filtered_desc, &
           dv_psi_evol, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress') then
      call reconcile_from_alpha_matrix_suppress(setting%num_filter, t, setting%suppress_constant, &
           dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
           dv_alpha_evol, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_orthogonal') then
      call reconcile_from_alpha_matrix_suppress_orthogonal(setting%num_filter, t, setting%suppress_constant, &
           dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
           dv_alpha_evol, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_adaptive') then
      call reconcile_from_alpha_matrix_suppress_adaptive(setting%num_filter, t, setting%suppress_constant, &
           dv_eigenvalues_prev, dv_eigenvalues, Y_filtered, Y_filtered_desc, YSY_filtered, YSY_filtered_desc, &
           dv_alpha_evol, dv_alpha_reconcile, dv_psi_reconcile)
    else if (trim(setting%re_initialize_method) == 'minimize_alpha_error') then
      dv_alpha_reconcile(:) = dv_alpha(:)
      call alpha_to_lcao_coef(Y_filtered, Y_filtered_desc, dv_eigenvalues, t, dv_alpha_reconcile, dv_psi_reconcile)
    else
      call terminate('re_initialize_state: unknown re-initialization method', 1)
    end if
    call get_filtering_errors(dim, dv_psi, dv_psi_reconcile, errors)
    ! backed up col_filtered_work1 is used here.
    !call get_re_initialize_errors(setting%num_filter, filtered_vecs, filtered_vecs_desc, &
    !     absolute_alpha_error, relative_alpha_error)

    if (check_master()) then
      write (0, '(A, F16.6, A, E26.16e3, A, E26.16e3, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] alpha error after re-initialization is ', &
           errors%absolute, ' (absolute) and ', errors%relative, ' (relative)'
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] re_initialize_from_lcao_coef() : compute charges form initial value'
    end if

    call compute_charges(dim, structure, S_sparse, dv_psi_reconcile, &
         dv_charge_on_basis, dv_charge_on_atoms, charge_moment)

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] re_initialize_from_lcao_coef() : compute energies form initial value'
    end if

    call compute_energies(setting, proc, structure, &
       H_sparse, S_sparse, Y_filtered, Y_filtered_desc, dv_charge_on_atoms, charge_factor, &
       filter_group_indices, Y_local, dv_eigenvalues, dv_psi_reconcile, dv_alpha_reconcile, &
       dv_atom_perturb, &
       H1_base, H1, H1_desc, energies)
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
end module wp_initialization_m
