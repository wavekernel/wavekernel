program main
  use wp_main_aux_m
  implicit none

  type(wp_setting_t) :: setting
  type(wp_process_t) :: proc
  type(wp_state_t) :: state
  ! Dummy variables velow are not used from this main.f90.
  real(8) :: dummy_eigenvalues(1), dummy_eigenvectors(1, 1)
  integer :: dummy_desc_eigenvectors(desc_size)

  ! Functions.
  integer :: iargc

  call mpi_init(state%ierr)
  call init_timers(state%wtime_total, state%wtime)

  if (iargc() == 0) then
    call print_help_and_stop()
  end if

  ! Setting, restart file and dimension of the problem are broadcast from the master node.
  ! Following variables are also the same:
  ! num_atoms, atom_indices, atom_coordinates, atom_elements, H_sparse, S_sparse, group_id
  call read_bcast_setting(setting, state%dim)
  call add_timer_event('main', 'read_bcast_setting', state%wtime)

  call setup_distribution(proc)
  call setup_distributed_matrices(state%dim, setting, proc, state)
  call add_timer_event('main', 'setup_distributed_matrices', state%wtime)

  call read_bcast_matrix_files(state%dim, setting, 1, state%H_sparse, state%S_sparse)
  call add_timer_event('main', 'read_bcast_matrix_files', state%wtime)

  call read_bcast_structure(state%dim, setting, 1, state%structure)
  ! Group ids are fixed through the whole simulation.
  if (setting%is_group_id_used) then
    call read_bcast_group_id(setting%group_id_filename, state%group_id)
  end if
  if (trim(setting%filter_mode) == 'group') then
    call read_bcast_group_id(setting%filter_group_filename, state%filter_group_id)
  end if
  call add_timer_event('main', 'read_bcast_atom_indices_and_coordinates', state%wtime)

  call allocate_dv_vectors(setting, state)

  state%charge_factor%charge_factor_common = setting%charge_factor_common
  state%charge_factor%charge_factor_H = setting%charge_factor_H
  state%charge_factor%charge_factor_C = setting%charge_factor_C

  call set_aux_matrices(state%dim, setting, proc, state, &
       .false., dummy_eigenvalues, dummy_desc_eigenvectors, dummy_eigenvectors)
  call add_timer_event('main', 'set_aux_matrices', state%wtime)

  call initialize(setting, proc, state)
  call add_timer_event('main', 'initialize', state%wtime)

  call prepare_json(setting, proc, state)
  ! add_structure_json() must be called after both of coordinates reading and prepare_json().
  call add_structure_json(0d0, 1, state)
  call add_timer_event('main', 'prepare_json', state%wtime)
  if (check_master()) then
    write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
         '] main loop start'
  end if

  if (setting%is_restart_mode) then
    state%i = setting%restart_step_num
    state%t = setting%restart_t
    state%t_last_replace = state%t
    state%total_state_count = setting%restart_total_states_count
    state%input_step = setting%restart_input_step
  else
    state%i = 0
    state%t = 0d0
    state%t_last_replace = state%t
    setting%restart_t = 0d0  ! For simulation time report.
    state%total_state_count = 0
    state%input_step = 1  ! Valid only in multiple step input mode.
  end if
  state%print_count = 1
  do
    if (check_master()) then
      if (state%t > (setting%limit_t - setting%restart_t) / 10d0 * dble(state%print_count) .or. &
           state%i == 1 .or. state%i == 10) then
        write (0, '(A, F16.6, A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
             '] simulation time ', state%t, ' a.u. start'
        state%print_count = state%print_count + 1
      end if
    end if

    if (mod(state%i, setting%output_interval) == 0) then
      call save_state(setting, .false., state)
      call add_timer_event('main', 'save_state', state%wtime)
    end if

    if (setting%delta_t * (state%i + 1) >= setting%limit_t) then
      exit  ! Exit before vain computation.
    end if

    ! Read next input step if multiple step input mode.
    ! Note that group_id does not change.
    if (setting%is_multistep_input_mode .and. &
         setting%multistep_input_read_interval * state%input_step <= state%t .and. &
         state%t - setting%delta_t < setting%multistep_input_read_interval * state%input_step) then

      if (check_master()) then
        write (0, '(A, F16.6, A, I0, A, F16.6)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
             '] read next input step ', state%input_step + 1, ' at simulation time ', state%t
      end if
      ! The step to be read is 'input_step + 1', not 'input_step + 2' because XYZ information is not interpolated.
      call read_bcast_structure(state%dim, setting, state%input_step + 1, state%structure)
      if (setting%to_replace_basis) then
        call copy_sparse_matrix(state%H_sparse, state%H_sparse_prev)
        call copy_sparse_matrix(state%S_sparse, state%S_sparse_prev)
        call read_bcast_matrix_files(state%dim, setting, state%input_step + 1, state%H_sparse, state%S_sparse)
      else
        call read_bcast_matrix_files(state%dim, setting, state%input_step + 1, &
             state%H_multistep_sparse, state%S_multistep_sparse)
      end if

      call set_aux_matrices_for_multistep(setting, proc, &
           .false., dummy_eigenvalues, dummy_desc_eigenvectors, dummy_eigenvectors, state)
      state%dv_psi(1 : state%dim) = state%dv_psi_reconcile(1 : state%dim)
      state%dv_alpha(1 : setting%num_filter) = state%dv_alpha_reconcile(1 : setting%num_filter)
      state%input_step = state%input_step + 1
      state%t_last_replace = state%t

      call add_structure_json(state%t, state%input_step, state)

      ! re-save state after matrix replacement.
      call save_state(setting, .true., state)
      call add_timer_event('main', 'save_state', state%wtime)
    end if

    call make_matrix_step_forward(setting, proc, state)

    call step_forward_post_process(setting, state)

    state%i = state%i + 1
    state%t = setting%delta_t * state%i
  end do

  call output_fson_and_destroy(setting, state%fsons(1)%output, state%fsons(1)%split_files_metadata, &
       state%fsons(1)%states, state%structures, state%wtime_total)
  call mpi_finalize(state%ierr)
end program main
