module wk_main_aux_m
  use mpi
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_processes_m
  use wk_fson_m
  use wk_fson_path_m
  use wk_fson_string_m
  use wk_fson_value_m
  use wk_event_logger_m
  use wk_atom_m
  use wk_charge_m
  use wk_conversion_m
  use wk_global_variables_m
  use wk_initialization_m
  use wk_linear_algebra_m
  use wk_matrix_generation_m
  use wk_output_m
  use wk_setting_m
  use wk_state_m
  use wk_time_evolution_m
  use wk_util_m
  implicit none

contains

  subroutine init_timers(wtime_total, wtime)
    real(8), intent(out) :: wtime_total, wtime
    integer :: ierr

    call mpi_barrier(mpi_comm_world, ierr)
    wtime_total = mpi_wtime()
    g_wk_mpi_wtime_init = wtime_total
    wtime = wtime_total
  end subroutine init_timers


  subroutine print_help_and_stop()
    integer :: ierr

    if (check_master()) then
      call print_help()
    end if
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)
    stop
  end subroutine print_help_and_stop


  subroutine add_timer_event(routine, event, wtime)
    character(*), intent(in) :: routine, event
    real(8), intent(inout) :: wtime

    real(8) :: wtime_end

    wtime_end = mpi_wtime()
    call add_event(routine // ':' // event, wtime_end - wtime)
    wtime = wtime_end
  end subroutine add_timer_event


  subroutine read_bcast_setting(setting, dim)
    type(wk_setting_t), intent(out) :: setting
    integer, intent(out) :: dim
    integer :: dim_overlap, num_atoms, num_groups, max_group_size, ierr

    if (check_master()) then
      call read_setting(setting)
      call inherit_setting(setting)

      call get_dimension(trim(setting%filename_hamiltonian), dim)
      if (.not. setting%is_overlap_ignored) then
        call get_dimension(trim(setting%filename_overlap), dim_overlap)
        if (dim /= dim_overlap) then
          call terminate('matrix sizes are inconsistent', ierr)
        end if
      end if

      if (trim(setting%filter_mode) == 'group') then
        call read_group_id_header(trim(setting%filter_group_filename), num_atoms, num_groups, max_group_size)
      end if

      call fill_filtering_setting(dim, num_groups, setting)
      call verify_setting(dim, setting)
      call print_setting(setting)
    end if
    call bcast_setting(g_wk_master_pnum, setting)
    call mpi_bcast(dim, 1, mpi_integer, g_wk_master_pnum, mpi_comm_world, ierr)
    call mpi_bcast(g_wk_block_size, 1, mpi_integer, g_wk_master_pnum, mpi_comm_world, ierr)
  end subroutine read_bcast_setting


  subroutine allocate_dv_vectors(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state

    integer :: n, a, m

    n = setting%num_filter
    a = state%structure%num_atoms
    m = setting%num_multiple_initials
    allocate(state%dv_alpha(n, m), state%dv_psi(state%dim, m), state%dv_alpha_next(n, m), state%dv_psi_next(state%dim, m))
    allocate(state%dv_alpha_reconcile(n, m), state%dv_psi_reconcile(state%dim, m))
    allocate(state%dv_atom_perturb(a), state%dv_atom_speed(a))
    allocate(state%dv_charge_on_basis(state%dim, m), state%dv_charge_on_atoms(a, m))
    allocate(state%charge_moment(m), state%energies(m), state%errors(m))
    allocate(state%basis%dv_eigenvalues(n), state%basis%dv_ipratios(n))
    allocate(state%basis%eigenstate_mean(3, n), state%basis%eigenstate_msd(4, n))
    allocate(state%basis%eigenstate_ipratios(setting%num_filter))
    allocate(state%fsons(m))
  end subroutine allocate_dv_vectors


  subroutine setup_distributed_matrices(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state

    integer :: num_groups, group, m_local, n_local, ierr

    if (check_master()) then
      call print_proc()
    end if

    state%basis%is_group_filter_mode = (trim(setting%filter_mode) == 'group')
    state%basis%dim = state%dim
    if (state%basis%is_group_filter_mode) then
      num_groups = size(state%basis%filter_group_id, 2)
      state%basis%num_basis = setting%num_filter * num_groups
      call filter_group_id_to_indices(state%basis%dim, setting%num_filter, state%structure, &
           state%basis%filter_group_id, state%basis%filter_group_indices)
      call initialize_distributed_block_diagonal_matrices(mpi_comm_world, state%basis%dim, state%basis%num_basis, &
           num_groups, state%basis%filter_group_indices(1, 1 : num_groups + 1), &
           state%basis%filter_group_indices(2, 1 : num_groups + 1), state%basis%Y_local)
    else
      state%basis%num_basis = setting%num_filter
    end if

    if (state%basis%is_group_filter_mode) then
      call setup_distributed_matrix_real('H1_base', &
           state%basis%num_basis, state%basis%num_basis, state%basis%H1_desc, state%basis%H1_base, .true.)
      call setup_distributed_matrix_real('S1', &
           state%basis%num_basis, state%basis%num_basis, state%basis%H1_desc, state%basis%S1, .true.)
      call setup_distributed_matrix_real('Y_filtered', &
           state%basis%num_basis, state%basis%num_basis, state%basis%Y_filtered_desc, state%basis%Y_filtered)
    else
      call setup_distributed_matrix_real('Y_all', state%dim, state%dim, &
           state%basis%Y_all_desc, state%basis%Y_all, .true.)
      call setup_distributed_matrix_real('Y_filtered', &
           state%dim, state%basis%num_basis, state%basis%Y_filtered_desc, state%basis%Y_filtered)
    end if
    call setup_distributed_matrix_real('H1', &
         state%basis%num_basis, state%basis%num_basis, state%basis%H1_desc, state%H1, .true.)
    call setup_distributed_matrix_complex('A', &
         state%basis%num_basis, state%basis%num_basis, state%basis%H1_desc, state%A, .true.)
    call setup_distributed_matrix_real('YSY_filtered', &
         state%basis%num_basis, state%basis%num_basis, state%YSY_filtered_desc, state%YSY_filtered)
  end subroutine setup_distributed_matrices


  ! read_matrix_file は allocation も同時に行う
  subroutine read_bcast_matrix_files(dim, setting, step, H_sparse, S_sparse)
    integer, intent(in) :: dim, step
    type(wk_setting_t), intent(in) :: setting
    type(sparse_mat), intent(inout) :: H_sparse, S_sparse

    real(8) :: wtime

    wtime = mpi_wtime()

    call destroy_sparse_mat(H_sparse)
    call destroy_sparse_mat(S_sparse)
    call add_timer_event('read_bcast_matrix_files', 'destroy_sparse_matrices', wtime)
    if (check_master()) then
      call read_matrix_file(trim(setting%filename_hamiltonian), step, H_sparse)
      call add_timer_event('read_bcast_matrix_files', 'read_matrix_file_H', wtime)
      if (setting%is_overlap_ignored) then
        call set_sparse_matrix_identity(dim, S_sparse)
        call add_timer_event('read_bcast_matrix_files', 'set_sparse_matrix_identity', wtime)
      else
        call read_matrix_file(trim(setting%filename_overlap), step, S_sparse)
        call add_timer_event('read_bcast_matrix_files', 'read_matrix_file_S', wtime)
      end if
    end if
    call bcast_sparse_matrix(g_wk_master_pnum, H_sparse)
    call bcast_sparse_matrix(g_wk_master_pnum, S_sparse)
    call add_timer_event('read_bcast_matrix_files', 'bcast_sparse_matrices', wtime)
  end subroutine read_bcast_matrix_files


  ! dv_psi is not referenced when is_initial is true.
  subroutine add_perturbation_for_zero_damp_charge(setting, is_initial, &
       dim, structure, S_sparse, basis, dv_psi, H_sparse)
    type(wk_setting_t), intent(in) :: setting
    logical, intent(in) :: is_initial
    integer, intent(in) :: dim
    type(sparse_mat), intent(in) :: S_sparse
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(in) :: basis
    complex(kind(0d0)), intent(in) :: dv_psi(:, :)

    type(sparse_mat), intent(inout) :: H_sparse

    integer :: i, j, k, ierr
    complex(kind(0d0)) :: psi(dim)
    real(8) :: diag, diag_acc, H1_lcao_diag(dim), dv_charge_on_atoms(structure%num_atoms), eigenvector(dim, 1)

    if (basis%is_group_filter_mode) then
      stop 'IMPLEMENT HERE (not needed unless zero_damp_charge is used)'
    end if

    if (setting%num_multiple_initials > 1) then
      call terminate('num_multiple_initials must be 1 in zero_damp_charge_* mode', 1)
    end if

    if (is_initial) then  ! In initiailization, perturbation source is HOMO eigenvector.
      call gather_matrix_real_part(basis%Y_filtered, basis%Y_filtered_desc, &
           1, setting%alpha_delta_index - setting%fst_filter + 1, dim, 1, 0, 0, eigenvector)
      if (check_master()) then
        psi(:) = cmplx(eigenvector(:, 1), 0d0, kind(0d0))
      end if
      call mpi_bcast(psi, dim, mpi_double_complex, g_wk_master_pnum, mpi_comm_world, ierr)
    else
      psi(:) = dv_psi(:, 1)
    end if

    ! Distributed redundant computation.
    if (trim(setting%h1_type) == 'zero_damp_charge_base') then
      do k = 1, H_sparse%num_non_zeros
        i = H_sparse%suffix(1, k)
        j = H_sparse%suffix(2, k)
        if (i == j) then
          H_sparse%value(k) = H_sparse%value(k) + setting%charge_factor_common * (abs(psi(i)) ** 2d0)
        end if
      end do
    else if (trim(setting%h1_type) == 'zero_damp_charge_atom') then
      diag_acc = 0d0
      call get_mulliken_charges_on_atoms(dim, structure, S_sparse, psi, dv_charge_on_atoms)
      do i = 1, structure%num_atoms
        diag = setting%charge_factor_common * dv_charge_on_atoms(i)  ! Positive value.
        print *, 'ZZZZcharge', i, diag
        do j = structure%atom_indices(i), structure%atom_indices(i + 1) - 1
          H1_lcao_diag(j) = diag / (structure%atom_indices(i + 1) - structure%atom_indices(i))
          diag_acc = diag_acc + H1_lcao_diag(j)
        end do
      end do
      print *, 'ZZZZchargediagacc', diag_acc, setting%charge_factor_common
      do k = 1, H_sparse%num_non_zeros
        i = H_sparse%suffix(1, k)
        j = H_sparse%suffix(2, k)
        if (i == j) then
          H_sparse%value(k) = H_sparse%value(k) + H1_lcao_diag(i)
        end if
      end do
    end if
  end subroutine add_perturbation_for_zero_damp_charge


  subroutine set_aux_matrices(setting, &
       to_use_precomputed_eigenpairs, eigenvalues, desc_eigenvectors, eigenvectors, state)
    type(wk_setting_t), intent(in) :: setting
    logical, intent(in) :: to_use_precomputed_eigenpairs
    integer, intent(in) :: desc_eigenvectors(desc_size)
    real(8), intent(in) :: eigenvalues(:), eigenvectors(:, :)
    type(wk_state_t), intent(inout) :: state

    integer :: end_filter, ierr
    real(8), allocatable :: eigenstate_charges(:, :)
    real(8) :: wtime

    wtime = mpi_wtime()

    if (.not. allocated(state%group_id)) then
      call make_dummy_group_id(state%structure%num_atoms, state%group_id)
    end if

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] set_aux_matrices() start '
    end if

    ! Needed only for the first step.
    call copy_sparse_matrix(state%H_sparse, state%H_sparse_prev)
    call copy_sparse_matrix(state%S_sparse, state%S_sparse_prev)

    if (state%basis%is_group_filter_mode) then
      call clear_offdiag_blocks_of_overlap(setting, state%basis%filter_group_indices, state%S_sparse)
      call add_timer_event('set_aux_matrices', 'clear_offdiag_blocks_of_overlap', wtime)
    end if

    if (to_use_precomputed_eigenpairs) then
      if (state%basis%is_group_filter_mode) then
        call terminate('to_use_precomputed_eigenpairs cannot be used with group filter mode', 1)
      end if
      if (size(eigenvalues, 1) /= state%dim) then
        call terminate('wrong size of eigenvalues', 1)
      end if
      if (desc_eigenvectors(rows_) /= state%dim .or. desc_eigenvectors(cols_) /= state%dim) then
        call terminate('wrong size of eigenvectors', 1)
      end if
      end_filter = setting%fst_filter + setting%num_filter - 1
      state%basis%dv_eigenvalues(1 : setting%num_filter) = eigenvalues(setting%fst_filter : end_filter)
      call pdgemr2d(state%dim, setting%num_filter, &
           eigenvectors, 1, setting%fst_filter, desc_eigenvectors, &
           state%basis%Y_filtered, 1, 1, state%basis%Y_filtered_desc, &
           state%basis%Y_filtered_desc(context_))
    else
      call set_eigenpairs(state%dim, setting, state%basis%filter_group_id, state%structure, &
           state%H_sparse, state%S_sparse, state%basis, &
           state%H1_lcao_sparse_charge_overlap)
    end if
    call add_timer_event('set_aux_matrices', 'set_eigenpairs', wtime)

    if (trim(setting%h1_type(1 : 16)) == 'zero_damp_charge') then
      call add_perturbation_for_zero_damp_charge(setting, .true., &
           state%dim, state%structure, state%S_sparse, &
           state%basis, state%dv_psi, state%H_sparse)
      ! After adding perturbation to the Hamiltonian, calculate eigenpairs again.
      call set_eigenpairs(state%dim, setting, state%basis%filter_group_id, state%structure, &
           state%H_sparse, state%S_sparse, state%basis, &
           state%H1_lcao_sparse_charge_overlap)
      call add_timer_event('set_aux_matrices', 'set_eigenpairs_second_for_zero_damp_charge_*', wtime)
    end if

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] get_msd_of_eigenstates() start '
    end if

    call get_msd_of_eigenstates(state%structure, state%S_sparse, state%basis)
    call get_ipratio_of_eigenstates(state%basis)
    call add_timer_event('set_aux_matrices', 'get_msd_and_ipratio', wtime)
  end subroutine set_aux_matrices


  subroutine set_aux_matrices_for_multistep(setting, &
       to_use_precomputed_eigenpairs, eigenvalues, desc_eigenvectors, eigenvectors, state)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: desc_eigenvectors(desc_size)
    logical, intent(in) :: to_use_precomputed_eigenpairs
    real(8), intent(in) :: eigenvalues(:), eigenvectors(:, :)
    type(wk_state_t), intent(inout) :: state

    integer :: j, end_filter, S_desc(desc_size), SY_desc(desc_size)
    real(8), allocatable :: dv_eigenvalues_prev(:)
    real(8), allocatable :: S(:, :), SY(:, :)
    real(8) :: wtime
    character(len=50) :: filename
    character(len=6) :: count_str
    integer, save :: count = 1
    integer, parameter :: iunit = 90
    type(wk_distributed_block_diagonal_matrices_t) :: Y_local_prev

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] set_aux_matrices_for_multistep() start '
    end if

    wtime = mpi_wtime()

    if (state%basis%is_group_filter_mode) then
      call clear_offdiag_blocks_of_overlap(setting, state%basis%filter_group_indices, state%S_sparse)
      call add_timer_event('read_next_input_step_with_basis_replace', 'clear_offdiag_blocks_of_overlap', wtime)
    end if

    ! Save eigenvalues in the previous step.
    allocate(dv_eigenvalues_prev(size(state%basis%dv_eigenvalues, 1)))
    dv_eigenvalues_prev(:) = state%basis%dv_eigenvalues(:)
    ! Compute SY or backup Y_local in advance because YSY needs eigenvectors in two steps.
    ! TODO: copy the basis variable and abstract generating YSY.
    if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_orthogonal' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_adaptive' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_select' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_damp') then
      if (state%basis%is_group_filter_mode) then
        call copy_allocation_distributed_block_diagonal_matrices(state%basis%Y_local, Y_local_prev)
        call copy_distributed_block_diagonal_matrices(state%basis%Y_local, Y_local_prev)
      else
        call setup_distributed_matrix_real('SY', state%dim, setting%num_filter, SY_desc, SY)
        call matmul_sd_d(state%S_sparse, state%basis%Y_filtered, state%basis%Y_filtered_desc, 1d0, SY, 1d0)
      end if
    end if

    if (to_use_precomputed_eigenpairs) then
      if (state%basis%is_group_filter_mode) then
        call terminate('to_use_precomputed_eigenpairs cannot be used with group filter mode', 1)
      end if
      if (size(eigenvalues, 1) /= state%dim) then
        call terminate('wrong size of eigenvalues', 1)
      end if
      if (desc_eigenvectors(rows_) /= state%dim .or. desc_eigenvectors(cols_) /= state%dim) then
        call terminate('wrong size of eigenvectors', 1)
      end if
      end_filter = setting%fst_filter + setting%num_filter - 1
      state%basis%dv_eigenvalues(1 : setting%num_filter) = eigenvalues(setting%fst_filter : end_filter)
      call pdgemr2d(state%dim, setting%num_filter, &
           eigenvectors, 1, setting%fst_filter, desc_eigenvectors, &
           state%basis%Y_filtered, 1, 1, state%basis%Y_filtered_desc, &
           desc_eigenvectors(context_))
    else
      if (trim(setting%h1_type(1 : 16)) == 'zero_damp_charge') then
        call add_perturbation_for_zero_damp_charge(setting, .false., &
             state%dim, state%structure, state%S_sparse, &
             state%basis, state%dv_psi, state%H_sparse)
      end if
      call set_eigenpairs(state%dim, setting, state%basis%filter_group_id, state%structure, &
           state%H_sparse, state%S_sparse, state%basis, &
           state%H1_lcao_sparse_charge_overlap)
    end if
    call add_timer_event('read_next_input_step_with_basis_replace', 'set_eigenpairs', wtime)

    if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_orthogonal' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_adaptive' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_select' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_damp') then
      if (state%basis%is_group_filter_mode) then
        call matmul_bsb_d(state%basis%Y_local, state%S_sparse, Y_local_prev, 1d0, &
             state%YSY_filtered, state%YSY_filtered_desc, 0d0)
      else
        call pdgemm('Trans', 'No', setting%num_filter, setting%num_filter, state%dim, 1d0, &
             state%basis%Y_filtered, 1, 1, state%basis%Y_filtered_desc, &
             SY, 1, 1, SY_desc, &
             0d0, &
             state%YSY_filtered, 1, 1, state%YSY_filtered_desc)
      end if
    end if

    if (setting%to_calculate_eigenstate_moment_every_step) then
      call get_msd_of_eigenstates(state%structure, state%S_sparse, state%basis)
      call add_timer_event('read_next_input_step_with_basis_replace', 'get_msd_of_eigenstates', wtime)
    end if

    !call get_ipratio_of_eigenstates(Y_filtered, Y_filtered_desc, ipratios)
    call add_timer_event('read_next_input_step_with_basis_replace', 'get_ipratio_of_eigenstates (skipped)', wtime)
    ! Should output information of eigenstates to JSON output.

    ! Time step was incremented after the last calculation of alpha,
    ! so the time for conversion between LCAO coefficient and eigenstate expansion
    ! is t - setting%delta_t, not t.
    do j = 1, setting%num_multiple_initials
      call re_initialize_state(setting, state%dim, state%t - setting%delta_t, state%structure, &
           state%charge_factor, j == 1, &
           state%H_sparse, state%S_sparse, state%H_sparse_prev, state%S_sparse_prev, &
           state%basis, state%YSY_filtered, state%YSY_filtered_desc, &
           dv_eigenvalues_prev, state%H1, state%basis%H1_desc, &
           state%dv_psi(:, j), state%dv_alpha(:, j), &
           state%dv_psi_reconcile(:, j), state%dv_alpha_reconcile(:, j), &
           state%dv_charge_on_basis(:, j), state%dv_charge_on_atoms(:, j), &
           state%dv_atom_perturb, &
           state%charge_moment(j), state%energies(j), state%errors(j))
    end do
    call add_timer_event('read_next_input_step_with_basis_replace', 're_initialize_from_lcao_coef', wtime)

    !if (.false.) then
    !  write(count_str, '(I6.6)') count
    !  filename = 'S_' // count_str // '.mtx'
    !  open(iunit, file=trim(filename))
    !  call print_sparse_matrix('S', state%S_sparse, iunit)
    !  close(iunit)
    !  filename = 'H0_' // count_str // '.mtx'
    !  open(iunit, file=trim(filename))
    !  call print_sparse_matrix('H0', state%H_sparse, iunit)
    !  close(iunit)
    !  filename = 'H1_' // count_str // '.mtx'
    !  open(iunit, file=trim(filename))
    !  call print_sparse_matrix('H1', state%H1_lcao_sparse_charge_overlap, iunit)
    !  close(iunit)
    !  filename = 'structure_' // count_str // '.txt'
    !  open(iunit, file=trim(filename))
    !  call print_structure(state%structure, iunit)
    !  close(iunit)
    !  count = count + 1
    !end if

    if (check_master()) then
      do j = 1, setting%num_multiple_initials
        write (0, '(A, F16.6, A, E26.16e3, A, E26.16e3, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
             '] psi error after re-initialization is ', state%errors(j)%absolute, &
             ' (absolute) and ', state%errors(j)%relative, ' (relative)'
      end do
    end if
  end subroutine set_aux_matrices_for_multistep


  subroutine read_bcast_structure(dim, setting, step, structure)
    integer, intent(in) :: dim, step
    type(wk_setting_t), intent(in) :: setting
    type(wk_structure_t), intent(out) :: structure

    if (allocated(structure%atom_indices)) then
      deallocate(structure%atom_indices)
    end if
    if (allocated(structure%atom_coordinates)) then
      deallocate(structure%atom_coordinates)
    end if
    if (allocated(structure%atom_elements)) then
      deallocate(structure%atom_elements)
    end if

    if (check_master()) then
      if (setting%is_atom_indices_enabled) then
        call read_structure(setting%atom_indices_filename, step, structure)
      else
        call make_dummy_structure(dim, structure)
      end if
    end if
    call bcast_structure(g_wk_master_pnum, structure)
    if (structure%atom_indices(structure%num_atoms + 1) - 1 /= dim) then
      stop 'invalid atom index file'
    else if (setting%is_restart_mode .and. trim(setting%h1_type) == 'maxwell') then
      if (size(setting%restart_atom_speed) /= structure%num_atoms) then
        stop 'dimension of the atom_speed in the restart file is not equal to num_atoms'
      end if
    else if (setting%is_restart_mode .and. trim(setting%h1_type) == 'harmonic') then
      if (size(setting%restart_atom_perturb) /= structure%num_atoms) then
        stop 'dimension of the atom_perturb in the restart file is not equal to num_atoms'
      end if
    end if
  end subroutine read_bcast_structure


  subroutine read_bcast_group_id(setting, num_atoms, group_id)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: num_atoms
    integer, allocatable, intent(out) :: group_id(:, :)

    if (check_master()) then
      if (setting%is_group_id_used) then
        call read_group_id(trim(setting%group_id_filename), group_id)
      else
        call make_dummy_group_id(num_atoms, group_id)
      end if
    end if
    call bcast_group_id(g_wk_master_pnum, group_id)
  end subroutine read_bcast_group_id


  subroutine set_eigenpairs(dim, setting, filter_group_id, structure, H_sparse, S_sparse, basis, &
       H1_lcao_sparse_charge_overlap)
    integer, intent(in) :: dim
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: filter_group_id(:, :)
    type(wk_structure_t), intent(in) :: structure
    type(sparse_mat), intent(in) :: H_sparse, S_sparse, H1_lcao_sparse_charge_overlap
    type(wk_basis_t), intent(inout) :: basis  ! value of last step is needed for offdiag norm calculation.

    integer :: i, j, num_groups, atom_min, atom_max, index_min, index_max, dim_sub, Y_sub_desc(desc_size)
    integer :: homo_level, filter_lowest_level, base_index_of_group, sum_dim_sub, my_rank_in_color
    integer :: H_desc(desc_size), S_desc(desc_size)
    real(8) :: elem
    real(8) :: dv_eigenvalues(dim), wtime
    integer :: YSY_desc(desc_size), l, i1, i2, j1, j2, m, n, lwork, liwork, g, p, ierr
    real(8), allocatable :: H(:, :), S(:, :), Y_sub(:, :), Y_filtered_last(:, :), SY(:, :), YSY(:, :), w(:)
    real(8), allocatable :: H1_lcao_charge_overlap(:, :), work(:)
    real(8), allocatable :: dv_block_eigenvalues_full_buf(:), dv_block_eigenvalues_filtered_buf(:)
    integer, allocatable :: filter_group_indices(:, :), iwork(:)
    type(sparse_mat) :: H_sparse_work

    ! For harmonic_for_nn_exciton.
    real(8) :: dv_atom_perturb(structure%num_atoms / 3 * 2), t_last_refresh, H1_lcao_diag(dim)

    ! Functions.
    integer :: numroc, indxg2p, indxg2l, indxl2g

    wtime = mpi_wtime()

    if (.not. basis%is_group_filter_mode) then
      call setup_distributed_matrix_real('H', dim, dim, H_desc, H, .true.)
      call setup_distributed_matrix_real('S', dim, dim, S_desc, S, .true.)
      call distribute_global_sparse_matrix_wk(H_sparse, H_desc, H)
      call distribute_global_sparse_matrix_wk(S_sparse, S_desc, S)

      if (trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
        ! Add perturbation term for Hamiltonian to get non-degenerated eigenstates.
        ! Now refreshing atom perturbation is forced. It may cause problem when multiple input mode.
        call aux_make_H1_harmonic_for_nn_exciton(dv_atom_perturb, setting%temperature, .true., &
             0d0, t_last_refresh, structure%num_atoms / 3, H1_lcao_diag)
        do i = 1, dim
          call pdelget('Self', ' ', elem, H, i, i, H_desc)
          call pdelset(H, i, i, H_desc, elem + H1_lcao_diag(i))
        end do
      else if (setting%h1_type == 'charge_overlap' .and. allocated(H1_lcao_sparse_charge_overlap%value)) then
        if (check_master()) then
          write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
               '] set_eigenpairs() : add charge overlap type Hamiltonian perturbation of previous step'
        end if
        call setup_distributed_matrix_real('H1_lcao_charge_overlap', dim, dim, &
             H_desc, H1_lcao_charge_overlap, .true.)
        call distribute_global_sparse_matrix_wk(H1_lcao_sparse_charge_overlap, H_desc, H1_lcao_charge_overlap)
        H(:, :) = H(:, :) + H1_lcao_charge_overlap(:, :)  ! distributed sum.
      end if

      allocate(Y_filtered_last(size(basis%Y_filtered, 1), size(basis%Y_filtered, 2)), &
           SY(size(basis%Y_filtered, 1), size(basis%Y_filtered, 2)))
      call setup_distributed_matrix_real('YSY', setting%num_filter, setting%num_filter, YSY_desc, YSY)
      !Y_filtered_last(:, :) = Y_filtered(:, :)
      call add_timer_event('set_eigenpairs', 'prepare_matrices_for_offdiag_norm_calculation', wtime)
      call solve_gevp(dim, 1, H, H_desc, S, S_desc, dv_eigenvalues, basis%Y_all, basis%Y_all_desc)
      call add_timer_event('set_eigenpairs', 'solve_gevp_in_filter_all', wtime)
      ! Filter eigenvalues.
      basis%dv_eigenvalues(1 : setting%num_filter) = &
           dv_eigenvalues(setting%fst_filter : setting%fst_filter + setting%num_filter - 1)
      ! Filter eigenvectors.
      call pdgemr2d(dim, setting%num_filter, &
           basis%Y_all, 1, setting%fst_filter, basis%Y_all_desc, &
           basis%Y_filtered, 1, 1, basis%Y_filtered_desc, &
           basis%Y_all_desc(context_))
      call add_timer_event('set_eigenpairs', 'copy_filtering_eigenvectors_in_filter_all', wtime)
    else
      call copy_sparse_matrix(H_sparse, H_sparse_work)
      if (trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
        ! Add perturbation term for Hamiltonian to get non-degenerated eigenstates.
        ! Now refreshing atom perturbation is forced. It may cause problem when multiple input mode.
        call aux_make_H1_harmonic_for_nn_exciton(dv_atom_perturb, setting%temperature, .true., &
             0d0, t_last_refresh, structure%num_atoms / 3, H1_lcao_diag)
        call add_diag_to_sparse_matrix(dim, H1_lcao_diag, H_sparse_work)
      else if (setting%h1_type == 'charge_overlap' .and. allocated(H1_lcao_sparse_charge_overlap%value)) then
        call terminate('set_eigenpairs: not implemented (charge_overlap with group filtering mode)', 1)
      end if
      ! Initialization for mpi_allreduce.
      allocate(dv_block_eigenvalues_full_buf(basis%dim), dv_block_eigenvalues_filtered_buf(basis%num_basis))
      dv_block_eigenvalues_full_buf(:) = 0d0
      dv_block_eigenvalues_filtered_buf(:) = 0d0
      do l = 1, numroc(basis%Y_local%num_blocks, 1, basis%Y_local%my_color_index, 0, basis%Y_local%num_colors)
        g = indxl2g(l, 1, basis%Y_local%my_color_index, 0, basis%Y_local%num_colors)
        i1 = basis%filter_group_indices(1, g)
        i2 = basis%filter_group_indices(1, g + 1)
        j1 = basis%filter_group_indices(2, g)
        j2 = basis%filter_group_indices(2, g + 1)
        m = i2 - i1
        n = j2 - j1
        allocate(w(m))
        call setup_distributed_matrix_with_context_real('H_group', basis%Y_local%my_blacs_context, m, m, &
             H_desc, H, .true.)
        call setup_distributed_matrix_with_context_real('S_group', basis%Y_local%my_blacs_context, m, m, &
             S_desc, S, .true.)
        call setup_distributed_matrix_with_context_real('Y_group', basis%Y_local%my_blacs_context, m, m, &
             Y_sub_desc, Y_sub, .true.)
        call distribute_global_sparse_matrix_wk_part(H_sparse_work, i1, i1, m, m, H_desc, H)
        call distribute_global_sparse_matrix_wk_part(S_sparse, i1, i1, m, m, S_desc, S)
        call solve_gevp(m, 1, H, H_desc, S, S_desc, w, Y_sub, Y_sub_desc)
        dv_block_eigenvalues_full_buf(i1 : i2 - 1) = w(1 : m)
        dv_block_eigenvalues_filtered_buf(j1 : j2 - 1) = &
             w(setting%fst_filter : setting%fst_filter + setting%num_filter - 1)
        call pdgemr2d(m, n, &
             Y_sub, 1, setting%fst_filter, Y_sub_desc, &
             basis%Y_local%local_matrices(l), 1, 1, basis%Y_local%descs(:, l), &
             basis%Y_local%descs(context_, l))
        deallocate(w, H, S, Y_sub)
      end do
      call mpi_comm_rank(basis%Y_local%my_comm, my_rank_in_color, ierr)
      if (my_rank_in_color > 0) then  ! Alleviate duplicated sum.
        dv_block_eigenvalues_full_buf(:) = 0d0
        dv_block_eigenvalues_filtered_buf(:) = 0d0
      end if
      call mpi_allreduce(dv_block_eigenvalues_full_buf, basis%Y_local%dv_block_eigenvalues_full, &
           basis%dim, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(dv_block_eigenvalues_filtered_buf, basis%Y_local%dv_block_eigenvalues_filtered, &
           basis%num_basis, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
      call change_basis_lcao_to_alpha_group_filter(filter_group_indices, basis%Y_local, &
           H_sparse_work, basis%H1_base, basis%H1_desc, .true.)
      call change_basis_lcao_to_alpha_group_filter(filter_group_indices, basis%Y_local, &
           S_sparse, basis%S1, basis%H1_desc, .true.)
      ! Filtering not needed after solve_gevp.
      call solve_gevp(basis%num_basis, 1, basis%H1_base, basis%H1_desc, basis%S1, basis%H1_desc, &
           basis%dv_eigenvalues, basis%Y_filtered, basis%Y_filtered_desc)
    end if

    deallocate(H, S)
    if (allocated(H1_lcao_charge_overlap)) then
      deallocate(H1_lcao_charge_overlap)
    end if
  end subroutine set_eigenpairs


  !subroutine set_local_eigenvectors(proc, Y_filtered_desc, Y_filtered, filter_group_indices, Y_local)
  !  integer, intent(in) :: Y_filtered_desc(desc_size), filter_group_indices(:, :)
  !  type(wk_process_t), intent(in) :: proc
  !  real(8), intent(in) :: Y_filtered(:, :)
  !  type(wk_local_matrix_t), allocatable, intent(out) :: Y_local(:)
  !  integer :: g, num_groups, i, j, m, n, p, nprow, npcol, myp, myprow, mypcol, np, prow, pcol
  !  integer :: blacs_pnum
  !  real(8) :: wtime
  !
  !  wtime = mpi_wtime()
  !
  !  num_groups = size(filter_group_indices, 2) - 1
  !  call blacs_gridinfo(proc%context, nprow, npcol, myprow, mypcol)
  !  np = nprow * npcol
  !  if (num_groups > np) then
  !    call terminate('the number of processes must be the number of groups or more', 1)
  !  end if
  !  myp = blacs_pnum(proc%context, myprow, mypcol)
  !  do g = 1, num_groups
  !    p = g - 1  ! Destination process.
  !    call blacs_pcoord(proc%context, p, prow, pcol)
  !    i = filter_group_indices(1, g)
  !    j = filter_group_indices(2, g)
  !    m = filter_group_indices(1, g + 1) - i
  !    n = filter_group_indices(2, g + 1) - j
  !    if (p == myp) then
  !      if (allocated(Y_local)) then
  !        deallocate(Y_local(1)%val)
  !        deallocate(Y_local)
  !      end if
  !      allocate(Y_local(1)%val(m, n))  ! Temporary implementation:
  !    end if
  !    call gather_matrix_real_part(Y_filtered, Y_filtered_desc, i, j, m, n, prow, pcol, Y_local(1)%val)
  !  end do
  !
  !  call add_timer_event('set_local_eigenvectors', 'set_local_eigenvectors', wtime)
  !end subroutine set_local_eigenvectors


  !subroutine set_H1_base(filter_group_indices, Y_local, H_sparse, H1_base, H1_desc)
  !  integer, intent(in) :: filter_group_indices(:, :), H1_desc(desc_size)
  !  type(sparse_mat), intent(in) :: H_sparse
  !  type(wk_distributed_block_diagonal_matrices_t), intent(in) :: Y_local
  !  real(8), intent(out) :: H1_base(:, :)
  !
  !  call change_basis_lcao_to_alpha_group_filter(filter_group_indices, Y_local, &
  !       H_sparse, H1_base, H1_desc, .true.)
  !end subroutine set_H1_base


  subroutine clear_offdiag_blocks_of_overlap(setting, filter_group_indices, S_sparse)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: filter_group_indices(:, :)
    type(sparse_mat) :: S_sparse

    if (trim(setting%filter_mode) == 'group') then
      call clear_offdiag_blocks(filter_group_indices, S_sparse)
    end if
  end subroutine clear_offdiag_blocks_of_overlap


  subroutine prepare_json(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state
    integer :: i, j, master_prow, master_pcol
    type(fson_value), pointer :: split_files_metadata_elem

    integer, parameter :: iunit_header = 21

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] prepare_json(): start'
    end if

    do j = 1, setting%num_multiple_initials
      ! Create whole output fson object.
      state%fsons(j)%output => fson_value_create()
      state%fsons(j)%output%value_type = TYPE_OBJECT
      call add_setting_json(setting, j, state%fsons(j)%output)
      if (check_master()) then
        call add_condition_json(state%dim, setting, state%structure, state%errors(j), &
             state%group_id, state%basis%filter_group_id, &
             state%basis%dv_eigenvalues, state%basis%eigenstate_mean, state%basis%eigenstate_msd, &
             state%basis%dv_ipratios, &
             state%fsons(j)%output)
        open(iunit_header, file=trim(add_postfix_to_filename(setting%output_filename, '_header')))
        call fson_print(iunit_header, state%fsons(j)%output)
        close(iunit_header)
      end if

      state%fsons(j)%states => fson_value_create()
      call fson_set_name('states', state%fsons(j)%states)
      state%fsons(j)%states%value_type = TYPE_ARRAY
      state%structures => fson_value_create()
      call fson_set_name('structures', state%structures)
      state%structures%value_type = TYPE_ARRAY

      if (setting%is_output_split) then
        state%fsons(j)%split_files_metadata => fson_value_create()
        call fson_set_name('split_files_metadata', state%fsons(j)%split_files_metadata)
        state%fsons(j)%split_files_metadata%value_type = TYPE_ARRAY
        if (setting%is_restart_mode) then
          call terminate('prepare_json: restart mode is not supported now', 1)
          !do i = 1, size(setting%restart_split_files_metadata)
          !  split_files_metadata_elem => fson_value_create()
          !  call fson_set_as_string(trim(setting%restart_split_files_metadata(i)), split_files_metadata_elem)
          !  call fson_value_add(state%fsons(j)%split_files_metadata, split_files_metadata_elem)
          !end do
        end if
      end if
    end do
  end subroutine prepare_json


  !subroutine read_next_input_step_with_basis_replace(dim, setting, proc, group_id, filter_group_id, t, structure, &
  !     H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
  !     Y_all_desc, Y_all, Y_filtered_desc, Y_filtered, YSY_filtered_desc, YSY_filtered, &
  !     H1_desc, H1, H1_base, &
  !     filter_group_indices, Y_local, charge_factor, &
  !     dv_psi, dv_alpha, dv_psi_reconcile, dv_alpha_reconcile, &
  !     dv_charge_on_basis, dv_charge_on_atoms, dv_atom_perturb, &
  !     dv_eigenvalues, charge_moment, energies, errors, eigenstate_ipratios, &
  !     eigenstate_mean, eigenstate_msd, &
  !     to_use_precomputed_eigenpairs, eigenvalues, desc_eigenvectors, eigenvectors, &
  !     H1_lcao_sparse_charge_overlap)
  !  integer, intent(in) :: dim, group_id(:, :), filter_group_id(:, :)
  !  type(wk_process_t), intent(in) :: proc
  !  type(wk_setting_t), intent(in) :: setting
  !  real(8), intent(in) :: t
  !  type(wk_structure_t), intent(in) :: structure
  !  type(wk_charge_factor_t), intent(in) :: charge_factor
  !  integer, allocatable, intent(inout) :: filter_group_indices(:, :)
  !  type(sparse_mat), intent(inout) :: H_sparse  ! Modified when h1_type == zero_damp_charge_*
  !  type(sparse_mat), intent(in) :: S_sparse, H_sparse_prev, S_sparse_prev, H1_lcao_sparse_charge_overlap
  !  integer, intent(in) :: Y_all_desc(desc_size), Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
  !  integer, intent(in) :: H1_desc(desc_size), desc_eigenvectors(desc_size)
  !  ! value of last step is needed for offdiag norm calculation.
  !  real(8), intent(inout) :: dv_eigenvalues(:), Y_filtered(:, :)
  !  real(8), intent(out) :: Y_all(:, :), YSY_filtered(:, :), H1(:, :), H1_base(:, :)
  !  real(8), intent(out) :: dv_charge_on_basis(:, :), dv_charge_on_atoms(:, :)
  !  real(8), intent(inout) :: dv_atom_perturb(:)
  !  complex(kind(0d0)), intent(in) :: dv_psi(:, :), dv_alpha(:, :)
  !  complex(kind(0d0)), intent(out) :: dv_psi_reconcile(:, :), dv_alpha_reconcile(:, :)
  !  type(wk_local_matrix_t), allocatable, intent(out) :: Y_local(:)
  !  type(wk_charge_moment_t), intent(out) :: charge_moment(:)
  !  type(wk_energy_t), intent(out) :: energies(:)
  !  type(wk_error_t), intent(out) :: errors(:)
  !  real(8), intent(out) :: eigenstate_ipratios(:)
  !  real(8), intent(out) :: eigenstate_mean(:, :), eigenstate_msd(:, :)
  !  logical, intent(in) :: to_use_precomputed_eigenpairs
  !  real(8), intent(in) :: eigenvalues(:), eigenvectors(:, :)
  !
  !  integer :: end_filter, S_desc(desc_size), SY_desc(desc_size), j
  !  real(8) :: dv_eigenvalues_prev(setting%num_filter)
  !  real(8), allocatable :: eigenstate_charges(:, :), S(:, :), SY(:, :)
  !  real(8) :: wtime
  !
  !end subroutine read_next_input_step_with_basis_replace


  !subroutine read_next_input_step_matrix(setting, &
  !     S_multistep_sparse, &
  !     filter_group_indices)
  !  integer, intent(in) :: filter_group_indices(:, :)
  !  type(wk_setting_t), intent(in) :: setting
  !  type(sparse_mat), intent(inout) :: S_multistep_sparse
  !
  !  call clear_offdiag_blocks_of_overlap(setting, filter_group_indices, S_multistep_sparse)
  !end subroutine read_next_input_step_matrix


  subroutine post_process_after_matrix_replace(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state

    state%dv_psi(1 : state%dim, :) = state%dv_psi_reconcile(1 : state%dim, :)
    state%dv_alpha(1 : setting%num_filter, :) = state%dv_alpha_reconcile(1 : setting%num_filter, :)
    state%input_step = state%input_step + 1
    state%t_last_replace = state%t
  end subroutine post_process_after_matrix_replace


  subroutine step_forward_sparse(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state

    real(8) :: norm
    real(8), allocatable :: Sw(:, :)
    complex(kind(0d0)), allocatable :: SP1(:, :), SP2(:, :), psi(:, :), SP1psi(:, :)
    integer, allocatable :: SP2IPIV(:)
    complex(kind(0d0)), allocatable :: Spsi(:)
    integer :: j, psi_desc(desc_size), desc(desc_size), ierr
    complex(kind(0d0)) :: dot
    ! Functions.
    complex(kind(0d0)) :: zdotc

    call make_matrices_for_sparse(setting%delta_t, setting%multistep_input_read_interval, &
         state%t, state%t_last_replace, &
         state%H_sparse, state%S_sparse, state%H_sparse_prev, state%S_sparse_prev, &
         SP1, SP2, SP2IPIV, Sw, desc)
    allocate(Spsi(state%dim))
    do j = 1, setting%num_multiple_initials
      call matvec_dd_z('No', Sw, desc, kOne, state%dv_psi(:, j), kZero, Spsi)
      dot = zdotc(state%dim, state%dv_psi(:, j), 1, Spsi, 1)
      norm = sqrt(truncate_imag(dot))
    end do
    call setup_distributed_matrix_complex('psi', state%dim, setting%num_multiple_initials, &
         psi_desc, psi)
    call setup_distributed_matrix_complex('SP1psi', state%dim, setting%num_multiple_initials, &
         psi_desc, SP1psi)
    call distribute_matrix_complex(g_context, 0, 0, state%dim, state%dv_psi, psi_desc, psi)

    call pzgemm('No', 'No', state%dim, setting%num_multiple_initials, state%dim, kOne, &
         SP1, 1, 1, desc, &
         psi, 1, 1, psi_desc, &
         kZero, &
         SP1psi, 1, 1, psi_desc)
    call pzgetrs('No', state%dim, setting%num_multiple_initials, SP2, 1, 1, desc, SP2IPIV, &
         SP1psi, 1, 1, psi_desc, ierr)
    call gather_matrix_complex_with_pzgemr2d(g_context, psi_desc, SP1psi, &
         0, 0, state%dim, state%dv_psi_next)

    do j = 1, setting%num_multiple_initials
      call matvec_dd_z('No', Sw, desc, kOne, state%dv_psi_next(:, j), kZero, Spsi)
      dot = zdotc(state%dim, state%dv_psi_next(:, j), 1, Spsi, 1)
      norm = sqrt(truncate_imag(dot))
      if (norm > 1d-16) then
        state%dv_psi_next(:, j) = state%dv_psi_next(:, j) / norm
      end if
    end do
    call mpi_bcast(state%dv_psi_next, state%dim * setting%num_multiple_initials, mpi_double_complex, &
         g_wk_master_pnum, mpi_comm_world, ierr)
  end subroutine step_forward_sparse


  subroutine make_matrix_step_forward(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state

    integer :: j, i, k
    real(8) :: wtime, norm, diff_eigenvalue, damp_factor, alpha_next_norm, amplitude_after_normalize
    integer :: ierr

    real(8) :: dznrm2

    wtime = mpi_wtime()

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] make_matrix_step_forward() start '
    end if

    if (trim(setting%h1_type) == 'zero') then
      ! Skip time evolution calculation.
      do j = 1, setting%num_multiple_initials
        do i = 1, setting%num_filter
          state%dv_alpha_next(i, j) = exp(- kImagUnit * state%basis%dv_eigenvalues(i) * setting%delta_t) * state%dv_alpha(i, j)
        end do
      end do
      call add_timer_event('make_matrix_step_forward', 'step_forward_linear', wtime)
    else if (trim(setting%h1_type(1 : 9)) == 'zero_damp') then
      ! Skip time evolution calculation.
      do j = 1, setting%num_multiple_initials
        if (trim(setting%init_type) == 'alpha_delta') then
          k = setting%num_filter !setting%alpha_delta_index - setting%fst_filter + 1
        else if (trim(setting%init_type) == 'alpha_delta_multiple') then
          k = setting%num_filter !setting%alpha_delta_multiple_indices(j) - setting%fst_filter + 1
        else
          stop 'zero_damp must be used with alpha_delta or alpha_delta_multiple'
        end if
        do i = 1, setting%num_filter
          diff_eigenvalue = state%basis%dv_eigenvalues(k) - state%basis%dv_eigenvalues(i)
          if (trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_damp') then
            damp_factor = 0d0
          else
            damp_factor = setting%eigenstate_damp_constant * (diff_eigenvalue ** 2d0)
          end if
          !print *, 'ZZZZZ zero_damp ', j, k, i, ': ', &
          !     state%dv_eigenvalues(k), state%dv_eigenvalues(i), damp_factor * setting%delta_t
          state%dv_alpha_next(i, j) = &
               exp(- kImagUnit * state%basis%dv_eigenvalues(i) * setting%delta_t) * &
               exp(- damp_factor * setting%delta_t) * &
               state%dv_alpha(i, j)
        end do
        alpha_next_norm = dznrm2(setting%num_filter, state%dv_alpha_next(:, j), 1)
        state%dv_alpha_next(:, j) = state%dv_alpha_next(:, j) / alpha_next_norm
        amplitude_after_normalize = abs(state%dv_alpha_next(k, j)) / abs(state%dv_alpha(k, j))
        if (check_master()) then
          write (0, '(A, F16.6, A, 2E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
               '] make_matrix_step_forward() t, t + dt ', state%t, state%t + setting%delta_t
          write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
               '] make_matrix_step_forward() alpha_next_norm ', alpha_next_norm
          write (0, '(A, F16.6, A, E26.16e3)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
               '] make_matrix_step_forward()  amplitude_after_normalize', amplitude_after_normalize
        end if
      end do
    else if (trim(setting%h1_type) == 'zero_sparse' .and. .not. state%basis%is_group_filter_mode) then
      call step_forward_sparse(setting, state)
    else
      ! dv_charge_on_atoms must not be referenced when setting%num_multiple_initials > 1.
      call make_H1(trim(setting%h1_type), state%structure, &
           state%S_sparse, state%basis, .false., setting%is_restart_mode, &
           state%t, setting%temperature, setting%delta_t, setting%perturb_interval, &
           state%dv_charge_on_basis(:, 1), state%dv_charge_on_atoms(:, 1), state%charge_factor, &
           state%dv_atom_perturb, &
           state%H1, state%basis%H1_desc)
      call add_timer_event('make_matrix_step_forward', 'make_matrix_H1_not_multistep', wtime)

      if (state%basis%is_group_filter_mode) then
        state%H1(:, :) = state%H1(:, :) + state%basis%H1_base(:, :)  ! Sum of distributed matrices.
      end if

      call make_A(state%H1, state%basis%H1_desc, state%basis%dv_eigenvalues, state%t, &
           setting%amplitude_print_threshold, setting%amplitude_print_interval, &
           setting%delta_t, setting%fst_filter, &
           state%A, state%basis%H1_desc)
      call add_timer_event('make_matrix_step_forward', 'make_matrix_A', wtime)

      ! In group filtering mode, S is ignored in this implementation.
      ! Note that step_forward destroys A.
      do j = 1, setting%num_multiple_initials
        call step_forward(setting%time_evolution_mode, state%A, state%basis%H1_desc, setting%delta_t, &
             state%dv_alpha(:, j), state%dv_alpha_next(:, j))
        call add_timer_event('make_matrix_step_forward', 'step_forward', wtime)
      end do
    end if
  end subroutine make_matrix_step_forward


  subroutine make_matrices_for_sparse(delta_t, multistep_input_read_interval, t, t_last_replace, &
       H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, &
       SP1, SP2, SP2IPIV, Sw, desc)
    real(8) :: delta_t, t, multistep_input_read_interval, t_last_replace, t_next, weight
    type(sparse_mat) :: H_sparse, S_sparse, H_sparse_prev, S_sparse_prev
    complex(kind(0d0)), allocatable :: SP1(:, :), SP2(:, :)
    integer, allocatable :: SP2IPIV(:), SpIPIV(:)
    real(8), allocatable :: H0(:, :), H1(:, :), S0(:, :), S1(:, :), Sw(:, :)
    complex(kind(0d0)), allocatable :: Hp(:, :), Sp(:, :)
    integer :: desc(desc_size), ierr, dim, m, n

    dim = H_sparse%size
    if (multistep_input_read_interval > 0d0) then
      weight = (t - t_last_replace) / multistep_input_read_interval
    else
      weight = 0d0
    end if
    !print *, 'ZZZZZZZZZZZZZzweight', weight
    call setup_distributed_matrix_real('H0', dim, dim, desc, H0, .true.)
    call setup_distributed_matrix_real('H1', dim, dim, desc, H1, .true.)
    call setup_distributed_matrix_real('S0', dim, dim, desc, S0, .true.)
    call setup_distributed_matrix_real('S1', dim, dim, desc, S1, .true.)
    call distribute_global_sparse_matrix_wk(H_sparse_prev, desc, H0)
    call distribute_global_sparse_matrix_wk(H_sparse, desc, H1)
    call distribute_global_sparse_matrix_wk(S_sparse_prev, desc, S0)
    call distribute_global_sparse_matrix_wk(S_sparse, desc, S1)
    m = size(H0, 1)
    n = size(H0, 2)
    allocate(Hp(m, n), Sp(m, n), SP1(m, n), SP2(m, n), Sw(m, n), &
         SpIPIV(dim), SP2IPIV(dim))
    Hp(:, :) = dcmplx(weight * H1(:, :) + (1d0 - weight) * H0(:, :), 0d0)
    Sw(:, :) = weight * S1(:, :) + (1d0 - weight) * S0(:, :)
    Sp(:, :) = dcmplx(Sw, 0d0)

    call pzgetrf(dim, dim, Sp, 1, 1, desc, SpIPIV, ierr)
    call pzgetrs('No', dim, dim, Sp, 1, 1, desc, SpIPIV, Hp, 1, 1, desc, ierr)
    SP1(:, :) = -kImagUnit * delta_t * 0.5d0 * Hp(:, :)
    call add_diag(desc, SP1, kOne)

    SP2(:, :) = kImagUnit * delta_t * 0.5d0 * Hp(:, :)
    call add_diag(desc, SP2, kOne)
    call pzgetrf(dim, dim, SP2, 1, 1, desc, SP2IPIV, ierr)
  end subroutine make_matrices_for_sparse


  subroutine step_forward_post_process(setting, state)
    type(wk_setting_t), intent(in) :: setting
    type(wk_state_t), intent(inout) :: state

    integer :: num_filter, j
    real(8) :: wtime
    complex(kind(0d0)) :: energy_tmp

    wtime = mpi_wtime()

    num_filter = state%basis%num_basis

    do j = 1, setting%num_multiple_initials
      if (trim(setting%h1_type) == 'zero_sparse' .and. trim(setting%filter_mode) /= 'group') then
        call lcao_coef_to_alpha(state%S_sparse, state%basis, &
             state%dv_psi_next(:, j), state%dv_alpha_next(:, j))
        state%dv_psi(:, j) = state%dv_psi_next(:, j)
      else
        call alpha_to_lcao_coef(state%basis, state%dv_alpha_next(:, j), state%dv_psi(:, j))
        call add_timer_event('step_forward_post_process', 'alpha_to_lcao_coef', wtime)
      end if

      call get_mulliken_charges_on_basis(state%dim, state%S_sparse, state%dv_psi(:, j), state%dv_charge_on_basis(:, j))
      call add_timer_event('step_forward_post_process', 'get_mulliken_charges_on_basis', wtime)

      call get_mulliken_charges_on_atoms(state%dim, state%structure, state%S_sparse, &
           state%dv_psi(:, j), state%dv_charge_on_atoms(:, j))
      call add_timer_event('step_forward_post_process', 'get_mulliken_charges_on_atoms', wtime)

      call get_mulliken_charge_coordinate_moments(state%structure, state%dv_charge_on_atoms(:, j), &
           state%charge_moment(j))
      call add_timer_event('step_forward_post_process', 'get_mulliken_charge_coordinate_moments', wtime)

      call compute_energies(setting, state%structure, &
           state%H_sparse, state%S_sparse, state%basis, &
           state%dv_charge_on_basis(:, j), state%dv_charge_on_atoms(:, j), state%charge_factor, &
           state%dv_psi(:, j), state%dv_alpha(:, j), &
           state%dv_atom_perturb, &
           state%H1, state%basis%H1_desc, state%energies(j))

      state%dv_alpha(:, j) = state%dv_alpha_next(:, j)
      call add_timer_event('step_forward_post_process', 'move_alpha_between_time_steps', wtime)
    end do
  end subroutine step_forward_post_process


  subroutine get_simulation_ranges(states, min_time, max_time, min_step_num, max_step_num, &
       min_input_step, max_input_step)
    type(fson_value), pointer, intent(in) :: states
    real(8), intent(out) :: min_time, max_time
    integer, intent(out) :: min_step_num, max_step_num, min_input_step, max_input_step
    type(fson_value), pointer :: state, state_elem

    min_time = huge(min_time)
    max_time = -huge(max_time)
    min_step_num = huge(min_step_num)
    max_step_num = -huge(max_step_num)
    min_input_step = huge(min_input_step)
    max_input_step = -huge(max_input_step)

    state => states%children
    do while (associated(state))
      state_elem => fson_value_get(state, 'time')
      min_time = min(min_time, state_elem%value_real)
      max_time = max(max_time, state_elem%value_real)
      state_elem => fson_value_get(state, 'step_num')
      min_step_num = min(min_step_num, state_elem%value_integer)
      max_step_num = max(max_step_num, state_elem%value_integer)
      state_elem => fson_value_get(state, 'input_step')
      min_input_step = min(min_input_step, state_elem%value_integer)
      max_input_step = max(max_input_step, state_elem%value_integer)
      state => state%next
    end do
  end subroutine get_simulation_ranges


  subroutine save_state(setting, is_after_matrix_replace, state)
    type(wk_setting_t), intent(in) :: setting
    logical, intent(in) :: is_after_matrix_replace
    type(wk_state_t), intent(inout) :: state

    real(8) :: wtime
    integer :: master_prow, master_pcol, states_count, input_step_backup, j
    real(8), allocatable :: atom_coordinates_backup(:, :)
    real(8) :: t_backup
    type(fson_value), pointer :: last_structure

    wtime = mpi_wtime()

    call add_timer_event('save_state', 'gather_matrix_to_save_state', wtime)

    if (check_master()) then
      do j = 1, setting%num_multiple_initials
        call add_state_json(state%dim, setting%num_filter, state%structure, &
             state%i, state%t, state%input_step, &
             state%energies(j), state%dv_alpha(:, j), state%dv_psi(:, j), &
             state%dv_charge_on_basis(:, j), state%dv_charge_on_atoms(:, j), state%charge_moment(j), &
             setting%h1_type, state%dv_atom_speed, state%dv_atom_perturb, &
             is_after_matrix_replace, state%H_sparse, state%fsons(j)%states)
      end do
      state%total_state_count = state%total_state_count + 1
      if (setting%is_output_split) then
        states_count = fson_value_count(state%fsons(1)%states)  ! states_count is same on all initials.
        ! 'i + setting%output_interval' in the next condition expression is
        ! the next step that satisfies 'mod(i, setting%output_interval) == 0'.
        if (states_count >= setting%num_steps_per_output_split .or. &
             setting%delta_t * (state%i + setting%output_interval) >= setting%limit_t) then
          do j = 1, setting%num_multiple_initials
            call output_split_states(setting, state%total_state_count, j, state%fsons(j)%states, &
                 state%structures, state%fsons(j)%split_files_metadata)
          end do
          ! Reset json array of states.
          ! Back up the last structure.
          allocate(atom_coordinates_backup(3, state%structure%num_atoms))
          last_structure => fson_value_get(state%structures, fson_value_count(state%structures))
          call fson_path_get(last_structure, 'time', t_backup)
          call fson_path_get(last_structure, 'input_step', input_step_backup)
          call fson_path_get(last_structure, 'coordinates_x', atom_coordinates_backup(1, :))
          call fson_path_get(last_structure, 'coordinates_y', atom_coordinates_backup(2, :))
          call fson_path_get(last_structure, 'coordinates_z', atom_coordinates_backup(3, :))
          ! Reset jsons.
          do j = 1, setting%num_multiple_initials
            call fson_destroy(state%fsons(j)%states)
            state%fsons(j)%states => fson_value_create()
            call fson_set_name('states', state%fsons(j)%states)
            state%fsons(j)%states%value_type = TYPE_ARRAY
          end do
          call fson_destroy(state%structures)
          state%structures => fson_value_create()
          call fson_set_name('structures', state%structures)
          state%structures%value_type = TYPE_ARRAY
          ! Write the last atomic structure to make the split file independent from others.
          call add_structure_json(setting, state)
        end if
      end if
    end if
    call add_timer_event('save_state', 'add_and_print_state', wtime)
  end subroutine save_state


  subroutine output_split_states(setting, total_state_count, multiple_initial_index, states, &
       structures, split_files_metadata)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: total_state_count, multiple_initial_index
    type(fson_value), pointer, intent(in) :: states, structures
    type(fson_value), pointer, intent(inout) :: split_files_metadata

    integer, parameter :: iunit_split = 12, iunit_split_binary = 18
    real(8) :: min_time, max_time
    integer :: min_step_num, max_step_num, min_input_step, max_input_step
    type(fson_value), pointer :: metadata, metadata_elem, whole
    real(8) :: wtime
    integer :: states_count, binary_num_last_value
    character(len=1024) :: filename, binary_filename, binary_filename_basename  ! Output splitting.
    integer :: basename_start_index

    wtime = mpi_wtime()

    states_count = fson_value_count(states)
    if (setting%num_multiple_initials == 1) then
      filename = add_numbers_to_filename(setting%output_filename, &
           total_state_count - states_count + 1, total_state_count)
    else
      filename = add_numbers_to_filename_multiple_initial(setting%output_filename, multiple_initial_index, &
           total_state_count - states_count + 1, total_state_count)
    end if
    call get_simulation_ranges(states, min_time, max_time, min_step_num, max_step_num, &
         min_input_step, max_input_step)
    metadata => fson_value_create()
    metadata%value_type = TYPE_OBJECt
    ! metadata element: filename.
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_OBJECT
    call fson_set_name('filename', metadata_elem)
    call fson_set_as_string(trim(remove_directory_from_filename(filename)), metadata_elem)
    call fson_value_add(metadata, metadata_elem)
    ! metadata element min_time:
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_REAL
    call fson_set_name('min_time', metadata_elem)
    metadata_elem%value_real = min_time
    call fson_value_add(metadata, metadata_elem)
    ! metadata element max_time:
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_REAL
    call fson_set_name('max_time', metadata_elem)
    metadata_elem%value_real = max_time
    call fson_value_add(metadata, metadata_elem)
    ! metadata element min_step_num:
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_INTEGER
    call fson_set_name('min_step_num', metadata_elem)
    metadata_elem%value_integer = min_step_num
    call fson_value_add(metadata, metadata_elem)
    ! metadata element max_step_num:
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_INTEGER
    call fson_set_name('max_step_num', metadata_elem)
    metadata_elem%value_integer = max_step_num
    call fson_value_add(metadata, metadata_elem)
    ! metadata element min_input_step:
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_INTEGER
    call fson_set_name('min_input_step', metadata_elem)
    metadata_elem%value_integer = min_input_step
    call fson_value_add(metadata, metadata_elem)
    ! metadata element max_input_step:
    metadata_elem => fson_value_create()
    metadata_elem%value_type = TYPE_INTEGER
    call fson_set_name('max_input_step', metadata_elem)
    metadata_elem%value_integer = max_input_step
    call fson_value_add(metadata, metadata_elem)
    ! Add metadata of current split file to the accumulation array.
    call fson_value_add(split_files_metadata, metadata)

    whole => fson_value_create()
    whole%value_type = TYPE_OBJECT
    call fson_value_add(whole, states)
    call fson_value_add(whole, structures)

    call add_timer_event('output_split_states', 'prepare_split_file_metadata', wtime)

    open(iunit_split, file=trim(filename))
    if (setting%is_binary_output_mode) then
      binary_filename = trim(filename) // '.dat'
      basename_start_index = index(binary_filename, '/', back=.true.) + 1
      binary_filename_basename = trim(binary_filename(basename_start_index : ))
      open(iunit_split_binary, file=trim(binary_filename), form='unformatted', access='stream')
      binary_num_last_value = 0
      call fson_print(iunit_split, whole, .false., .true., 0, &
           iunit_split_binary, trim(binary_filename_basename), binary_num_last_value)
      close(iunit_split_binary)
    else
      call fson_print(iunit_split, whole)
    end if
    close(iunit_split)

    call add_timer_event('output_split_states', 'fson_print', wtime)
  end subroutine output_split_states


  subroutine output_fson_and_destroy(setting, multiple_initial_index, output, split_files_metadata, &
       states, structures, wtime_total)
    type(wk_setting_t), intent(in) :: setting
    integer, intent(in) :: multiple_initial_index
    type(fson_value), pointer, intent(inout) :: output, split_files_metadata, states, structures
    real(8), intent(inout) :: wtime_total

    integer :: iunit_master
    character(len=1024) :: actual_output_filename
    real(8) :: wtime

    wtime = mpi_wtime()

    if (check_master()) then
      if (setting%is_output_split) then
        call fson_value_add(output, split_files_metadata)
      else
        call fson_value_add(output, states)
        call fson_value_add(output, structures)
      end if

      call add_timer_event('main', 'write_states_to_fson', wtime)
      call add_timer_event('main', 'total', wtime_total)

      if (setting%num_multiple_initials == 1) then
        actual_output_filename = setting%output_filename
      else
        actual_output_filename = add_number_to_filename(setting%output_filename, multiple_initial_index)
      end if

      if (trim(actual_output_filename) /= '-') then
        iunit_master = 11
        open(iunit_master, file=actual_output_filename)
      else
        iunit_master = 6
      end if
      call fson_events_add(output)
      call fson_print(iunit_master, output)
      if (trim(actual_output_filename) /= '-') then
        close(iunit_master)
      end if
    end if

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] main loop end'
    end if

    call fson_destroy(output)
  end subroutine output_fson_and_destroy
end module wk_main_aux_m
