!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
module M_wavepacket
  use M_config  ! Global variables: 'config'.
  use M_group_id_setting  ! Global variables: 'num_groups', 'num_group_mem', 'group_mem'.
  use M_ext_matrix_data  ! Global variables: 'matrix_data'.
  use wp_main_aux_m
  use mpi
  implicit none

  type wp_state_t
    type(process) :: proc
    integer :: dim, num_atoms
    integer :: i, total_state_count, input_step, print_count
    real(8) :: t, wtime_total, wtime
    real(8) :: charge_coordinate_mean(3), charge_coordinate_msd(4)
    real(8) :: tightbinding_energy, nonlinear_term_energy, total_energy
    real(8) :: absolute_filter_error, relative_filter_error
    integer, allocatable :: atom_indices(:), group_id(:, :), filter_group_id(:, :), filter_group_indices(:, :)
    complex(kind(0d0)), allocatable :: full_vecs_local(:, :), filtered_vecs_local(:, :)
    real(8), allocatable :: atom_coordinates(:, :)
    character, allocatable :: atom_elements(:)
    type(fson_value), pointer :: output, states, split_files_metadata
    ! Distributed matrices.
    type(wp_local_matrix_t), allocatable :: Y_local(:)  ! Not ScaLAPACK-like manner.
    complex(kind(0d0)), allocatable :: H(:, :), S(:, :), Y(:, :), S_inv_sqrt(:, :)  ! m x m
    ! Following matrices for multistep perturbation mode share descriptor with H, S and S_inv_sqrt.
    complex(kind(0d0)), allocatable :: H_multistep(:, :), S_multistep(:, :), S_inv_sqrt_multistep(:, :)
    complex(kind(0d0)), allocatable :: H_multistep1(:, :), S_multistep1(:, :), S_inv_sqrt_multistep1(:, :)  ! m x m
    complex(kind(0d0)), allocatable :: H_multistep2(:, :), S_multistep2(:, :), S_inv_sqrt_multistep2(:, :)  ! m x m
    complex(kind(0d0)), allocatable :: Y_filtered(:, :)  ! m x n
    complex(kind(0d0)), allocatable :: H1(:, :), H1_base(:, :), A(:, :)  ! n x n
    complex(kind(0d0)), allocatable :: H1_multistep(:, :), H1_multistep1(:, :), H1_multistep2(:, :)
    complex(kind(0d0)), allocatable :: full_vecs(:, :), filtered_vecs(:, :)
    integer :: H_desc(desc_size), S_desc(desc_size), Y_desc(desc_size), S_inv_sqrt_desc(desc_size), &
         Y_filtered_desc(desc_size), H1_desc(desc_size), A_desc(desc_size), &
         full_vecs_desc(desc_size), filtered_vecs_desc(desc_size)
  end type wp_state_t

  private
  public :: wp_state_t, wavepacket_init, wavepacket_replace_matrix, wavepacket_main

contains

  subroutine read_atom_indices_and_coordinates_from_ELSES_config(dim, &
       num_atoms, atom_indices, atom_coordinates, atom_elements)
    integer, intent(in) :: dim
    integer, intent(out) :: num_atoms
    integer, allocatable, intent(out) :: atom_indices(:)
    real(8), allocatable, intent(out) :: atom_coordinates(:, :)
    character, allocatable, intent(out) :: atom_elements(:)
    integer :: i, atom_valence

    if (allocated(atom_indices)) then
      deallocate(atom_indices)
    end if
    if (allocated(atom_coordinates)) then
      deallocate(atom_coordinates)
    end if
    if (allocated(atom_elements)) then
      deallocate(atom_elements)
    end if

    num_atoms = config%system%structure%natom
    allocate(atom_indices(num_atoms + 1), atom_coordinates(3, num_atoms), atom_elements(num_atoms))

    atom_indices(1) = 1
    do i = 1, num_atoms
      atom_coordinates(1, i) = config%system%structure%vatom(i)%position(1)
      atom_coordinates(2, i) = config%system%structure%vatom(i)%position(2)
      atom_coordinates(3, i) = config%system%structure%vatom(i)%position(3)
      atom_elements(i) = config%system%structure%vatom(i)%element%name
      if (atom_elements(i) == 'H') then
        atom_valence = 1
      else if (atom_elements(i) == 'C') then
        atom_valence = 4
      else
        stop 'unknown atom name'
      end if
      atom_indices(i + 1) = atom_indices(i) + atom_valence
    end do

    ! bcast of atom data is not needed because the procedure above is executed on whole the nodes redundantly.
  end subroutine read_atom_indices_and_coordinates_from_ELSES_config


  subroutine convert_group_elses_to_wp(num_groups_, num_group_mem_, group_mem_, group_id)
    ! Using underbars, identifiers below are distinguished from the global variables
    ! in the module M_group_id_setting.
    integer, intent(in) :: num_groups_, num_group_mem_(:), group_mem_(:, :)
    integer, allocatable, intent(out) :: group_id(:, :)

    integer :: max_group_size

    print *, num_groups_, num_group_mem_, group_mem_

    if (allocated(group_id)) then
      deallocate(group_id)
    end if

    max_group_size = maxval(num_group_mem_, 1)
    allocate(group_id(max_group_size + 1, num_groups))
    group_id(1, :) = num_group_mem_(:)
    group_id(2 :, :) = group_mem_(:, :)
  end subroutine convert_group_elses_to_wp


  subroutine convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data_, wp_sparse)
    type(sparce_matrix_data_real_type) :: matrix_data_
    type(sparse_mat), intent(out) :: wp_sparse

    wp_sparse%size = matrix_data_%matrix_size
    wp_sparse%num_non_zeros = matrix_data_%num_of_non_zero_elements
    allocate(wp_sparse%suffix(2, wp_sparse%num_non_zeros), wp_sparse%value(wp_sparse%num_non_zeros))
    wp_sparse%suffix(:, :) = matrix_data_%element_index(:, :)
    wp_sparse%value(:) = matrix_data_%element_data(:)
  end subroutine convert_sparse_matrix_data_real_type_to_wp_sparse


  subroutine copy_settings_from_elses_config_xml(setting)
    type(wp_setting_t), intent(out) :: setting

    ! Required settings.
    setting%delta_t = config%calc%wave_packet%delta_t
    setting%limit_t = config%calc%wave_packet%limit_t
    setting%multistep_input_read_interval = config%calc%wave_packet%replace_t

    ! Optional settings.
    ! Copy perturbation (H1) settings.
    if (trim(config%calc%wave_packet%h1_type) == '') then  ! Set default value.
      setting%h1_type = 'zero'
    else
      setting%h1_type = trim(config%calc%wave_packet%h1_type)
    end if
    if (config%calc%wave_packet%charge_factor_common < 0d0) then  ! Set default value.
      setting%charge_factor_common = 0d0
    else
      setting%charge_factor_common = config%calc%wave_packet%charge_factor_common
    end if
    if (config%calc%wave_packet%charge_factor_H >= 0d0) then
      setting%charge_factor_H = config%calc%wave_packet%charge_factor_H
    end if
    if (config%calc%wave_packet%charge_factor_C >= 0d0) then
      setting%charge_factor_C = config%calc%wave_packet%charge_factor_C
    end if
    setting%temperature = config%calc%wave_packet%temperature
    setting%perturb_interval = config%calc%wave_packet%perturb_interval

    ! Copy filtering settings.
    if (config%calc%wave_packet%filter_mode == '') then  ! Set default value.
      setting%filter_mode = 'all'
    else
      setting%filter_mode = config%calc%wave_packet%filter_mode
    end if
    if (config%calc%wave_packet%fst_filter <= 0 .or. config%calc%wave_packet%end_filter <= 0) then
      ! Set default value.
      ! If fst_filter == 0, fill_filtering_setting will work.
      setting%fst_filter = 0
    else
      setting%fst_filter = config%calc%wave_packet%fst_filter
      setting%num_filter = config%calc%wave_packet%end_filter - setting%fst_filter + 1
    end if
    setting%filter_group_filename = config%calc%wave_packet%filter_group_filename
    if (trim(setting%filter_mode) == 'group') then
      setting%num_filter = config%calc%wave_packet%num_group_filter_from_homo
    end if
    ! Copy initializatin settings.
    if (trim(config%calc%wave_packet%init_type) == '') then  ! Set default value.
      setting%init_type = 'alpha_delta'
    else
      setting%init_type = trim(config%calc%wave_packet%init_type)
    end if
    if (config%calc%wave_packet%alpha_delta_index <= 0) then  ! Set default value.
      setting%alpha_delta_index = setting%num_filter - setting%fst_filter + 1
    else
      setting%alpha_delta_index = config%calc%wave_packet%alpha_delta_index
    end if
    setting%to_multiply_phase_factor = config%calc%wave_packet%to_multiply_phase_factor
    setting%phase_factor_coef = config%calc%wave_packet%phase_factor_coef
    setting%localize_potential_depth = config%calc%wave_packet%localize_potential_depth
    setting%alpha_delta_min_x = config%calc%wave_packet%alpha_delta_min_x
    setting%alpha_delta_max_x = config%calc%wave_packet%alpha_delta_max_x
    setting%localize_start = config%calc%wave_packet%localize_start
    setting%localize_end = config%calc%wave_packet%localize_end
    ! Copy time evolution settings.
    if (trim(config%calc%wave_packet%time_evolution_mode) == '') then  ! Set default value.
      setting%time_evolution_mode = 'crank_nicolson'
    else
      setting%time_evolution_mode = trim(config%calc%wave_packet%time_evolution_mode)
    end if
    ! Copy message settings.
    if (config%calc%wave_packet%amplitude_print_threshold >= 0d0) then
      setting%amplitude_print_threshold = config%calc%wave_packet%amplitude_print_threshold
    end if
    if (config%calc%wave_packet%amplitude_print_interval >= 0d0) then
      setting%amplitude_print_interval = config%calc%wave_packet%amplitude_print_interval
    end if

    if (trim(config%calc%wave_packet%output_filename) == '') then  ! Set default value.
      setting%output_filename = 'out.json'
    else
      setting%output_filename = trim(config%calc%wave_packet%output_filename)
    end if
    setting%is_binary_output_mode = config%calc%wave_packet%is_binary_output_mode
    setting%is_output_split = config%calc%wave_packet%is_output_split
    if (config%calc%wave_packet%output_interval > 0) then
      setting%output_interval = config%calc%wave_packet%output_interval
    end if
    if (config%calc%wave_packet%num_steps_per_output_split > 0) then
      setting%num_steps_per_output_split = config%calc%wave_packet%num_steps_per_output_split
    end if
    ! Copy overlap matrix settings.
    setting%is_overlap_ignored = config%calc%wave_packet%is_overlap_ignored
    ! Settings that automatically determined when called from ELSES.
    setting%is_atom_indices_enabled = .true.
    setting%is_group_id_used = num_groups > 0 .and. allocated(num_group_mem) .and. allocated(group_mem)
    setting%is_multistep_input_mode = .true.
    setting%to_replace_basis = .true.
  end subroutine copy_settings_from_elses_config_xml


  subroutine wavepacket_init(setting, state)
    type(wp_setting_t), intent(inout) :: setting
    type(wp_state_t), intent(inout) :: state
    type(sparse_mat) :: H_sparse, S_sparse, H_multistep_sparse, S_multistep_sparse
    integer :: ierr

    call init_timers(state%wtime_total, state%wtime)

    state%dim = matrix_data(1)%matrix_size
    call copy_settings_from_elses_config_xml(setting)
    call fill_filtering_setting(state%dim, num_groups, setting)
    call verify_setting(state%dim, setting)
    call print_setting(setting)

    call setup_distributed_matrices(state%dim, setting, state%proc, &
         state%H_desc, state%H, state%S_desc, state%S, state%Y_desc, state%Y, &
         state%Y_filtered_desc, state%Y_filtered, state%S_inv_sqrt_desc, state%S_inv_sqrt, &
         state%H1_desc, state%H1, state%H1_base, state%A_desc, state%A, &
         state%H_multistep, state%S_multistep, state%S_inv_sqrt_multistep, &
         state%H_multistep1, state%S_multistep1, state%S_inv_sqrt_multistep1, &
         state%H_multistep2, state%S_multistep2, state%S_inv_sqrt_multistep2, &
         state%H1_multistep, state%H1_multistep1, state%H1_multistep2, &
         state%full_vecs_desc, state%full_vecs, state%filtered_vecs_desc, state%filtered_vecs, &
         state%full_vecs_local, state%filtered_vecs_local)
    call add_timer_event('main', 'setup_distributed_matrices', state%wtime)

    ! step = 1
    call convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data(1), H_sparse)
    call convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data(2), S_sparse)
    call add_timer_event('main', 'set_matrices', state%wtime)

    call read_atom_indices_and_coordinates_from_ELSES_config(state%dim, &
         state%num_atoms, state%atom_indices, state%atom_coordinates, state%atom_elements)
    ! Group ids are fixed through the whole simulation.
    if (setting%is_group_id_used) then
      call convert_group_elses_to_wp(num_groups, num_group_mem, group_mem, state%group_id)
    end if
    if (trim(setting%filter_mode) == 'group') then
      call convert_group_elses_to_wp(num_groups, num_group_mem, group_mem, state%filter_group_id)
    end if
    call add_timer_event('main', 'set_atom_indices_and_coordinates_and_group_id', state%wtime)

    if (trim(setting%h1_type) == 'multistep') then
      ! step = 2
      call convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data(1), H_multistep_sparse)
      call convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data(2), S_multistep_sparse)
    end if

    call set_aux_matrices(state%dim, setting, state%proc, state%num_atoms, &
         H_sparse, S_sparse, H_multistep_sparse, S_multistep_sparse, &
         state%H_desc, state%S_desc, state%H, state%S, state%S_inv_sqrt_desc, state%S_inv_sqrt, &
         state%H_multistep1, state%H_multistep2, &
         state%S_multistep1, state%S_multistep2, &
         state%S_inv_sqrt_multistep1, state%S_inv_sqrt_multistep2, &
         state%H1_multistep1, state%H1_multistep2, &
         state%full_vecs_desc, state%full_vecs, state%filtered_vecs_desc, state%filtered_vecs, &
         state%Y_desc, state%Y, state%Y_filtered_desc, state%Y_filtered, state%Y_local, &
         state%atom_indices, state%atom_coordinates, state%filter_group_id, state%filter_group_indices, &
         state%H1_base, state%H1_desc)
    call add_timer_event('main', 'set_aux_matrices', state%wtime)

    call initialize(setting, state%proc, state%dim, &
         state%num_atoms, state%atom_indices, state%atom_coordinates, state%atom_elements, state%group_id, &
         state%H, state%H_desc, state%S, state%S_desc, &
         state%full_vecs, state%full_vecs_desc, state%filtered_vecs, state%filtered_vecs_desc, &
         state%Y_filtered, state%Y_filtered_desc, &
         state%filter_group_indices, state%Y_local, &
         state%H1_base, state%H1, state%H1_desc, &
         state%charge_coordinate_mean, state%charge_coordinate_msd, &
         state%tightbinding_energy, state%nonlinear_term_energy, state%total_energy, &
         state%absolute_filter_error, state%relative_filter_error, state%A_desc)
    call add_timer_event('main', 'initialize', state%wtime)

    call prepare_json(state%dim, state%num_atoms, setting, state%proc, &
         state%absolute_filter_error, state%relative_filter_error, &
         state%filtered_vecs_desc, state%filtered_vecs, state%filtered_vecs_local, &
         state%output, state%states, state%split_files_metadata)
    call add_timer_event('main', 'prepare_json', state%wtime)

    if (check_master()) then
      write (0, '(A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] main loop start'
    end if

    if (setting%is_restart_mode) then
      state%i = setting%restart_step_num
      state%t = setting%restart_t
      state%total_state_count = setting%restart_total_states_count
      state%input_step = setting%restart_input_step
    else
      state%i = 0
      state%t = 0d0
      setting%restart_t = 0d0  ! For simulation time report.
      state%total_state_count = 0
      state%input_step = 1  ! Valid only in multiple step input mode.
    end if
    state%print_count = 1
  end subroutine wavepacket_init


  ! Read next input step if multiple step input mode.
  ! Note that group_id does not change.
  subroutine wavepacket_replace_matrix(setting, state)
    type(wp_setting_t), intent(in) :: setting
    type(wp_state_t), intent(inout) :: state
    type(sparse_mat) :: H_sparse, S_sparse, H_multistep_sparse, S_multistep_sparse

    if (check_master()) then
      write (0, '(A, F16.6, A, I0, A, F16.6)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
           '] read next input step ', state%input_step + 1, ' at simulation time ', state%t
    end if
    ! The step to be read is 'input_step + 1', not 'input_step + 2' because XYZ information is not interpolated.
    call read_atom_indices_and_coordinates_from_ELSES_config(state%dim, &
         state%num_atoms, state%atom_indices, state%atom_coordinates, state%atom_elements)
    if (setting%to_replace_basis) then
      ! step = input_step + 1 (call set_sample_matrices(dim, setting, input_step + 1, H_sparse, S_sparse))
      call convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data(1), H_sparse)
      call convert_sparse_matrix_data_real_type_to_wp_sparse(matrix_data(2), S_sparse)

      if (trim(setting%h1_type) == 'multistep') then
        stop 'matrix interpolation is not supported when called from ELSES'
      end if
    else
      stop 'basis replace mode must be used when called from ELSES'
    end if

    call set_aux_matrices_for_multistep(state%dim, state%num_atoms, setting, state%proc, state%filter_group_id, state%t, &
         state%atom_indices, state%atom_coordinates, state%atom_elements, &
         H_sparse, state%H_desc, state%H, &
         S_sparse, state%S_desc, state%S, &
         state%S_inv_sqrt_desc, state%S_inv_sqrt, state%Y_desc, state%Y, state%Y_filtered_desc, state%Y_filtered, &
         H_multistep_sparse, S_multistep_sparse, &
         state%H_multistep, state%H_multistep1, state%H_multistep2, &
         state%S_multistep, state%S_multistep1, state%S_multistep2, &
         state%S_inv_sqrt_multistep, state%S_inv_sqrt_multistep1, state%S_inv_sqrt_multistep2, &
         state%H1_desc, state%H1, state%H1_base, &
         state%H1_multistep, state%H1_multistep1, state%H1_multistep2, &
         state%full_vecs_desc, state%full_vecs, state%filtered_vecs_desc, state%filtered_vecs, &
         state%filter_group_indices, state%Y_local, &
         state%charge_coordinate_mean, state%charge_coordinate_msd, &
         state%tightbinding_energy, state%nonlinear_term_energy, state%total_energy, &
         state%absolute_filter_error, state%relative_filter_error)
    state%input_step = state%input_step + 1

    ! re-save state after matrix replacement.
    call save_state(state%dim, state%num_atoms, setting, state%i, state%t, &
         state%tightbinding_energy, state%nonlinear_term_energy, state%total_energy, &
         state%charge_coordinate_mean, state%charge_coordinate_msd, &
         state%full_vecs_desc, state%full_vecs, state%filtered_vecs_desc, state%filtered_vecs, &
         state%full_vecs_local, state%filtered_vecs_local, &
         state%split_files_metadata, state%total_state_count, state%input_step, .true., state%states)
  end subroutine wavepacket_replace_matrix


  subroutine wavepacket_main(setting, state)
    type(wp_setting_t), intent(in) :: setting
    type(wp_state_t), intent(inout) :: state
    integer :: ierr

    do
      ! Print progress report.
      if (check_master()) then
        if (state%t > (setting%limit_t - setting%restart_t) / 10d0 * dble(state%print_count) .or. &
             state%i == 1 .or. state%i == 10) then
          write (0, '(A, F16.6, A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
               '] simulation time ', state%t, ' a.u. start'
          state%print_count = state%print_count + 1
        end if
      end if

      ! Output for files.
      if (mod(state%i, setting%output_interval) == 0) then
        call save_state(state%dim, state%num_atoms, setting, state%i, state%t, &
             state%tightbinding_energy, state%nonlinear_term_energy, state%total_energy, &
             state%charge_coordinate_mean, state%charge_coordinate_msd, &
             state%full_vecs_desc, state%full_vecs, state%filtered_vecs_desc, state%filtered_vecs, &
             state%full_vecs_local, state%filtered_vecs_local, &
             state%split_files_metadata, state%total_state_count, state%input_step, .false., state%states)
        call add_timer_event('main', 'save_state', state%wtime)
      end if

      ! Exit loop before vain computation.
      if (setting%delta_t * (state%i + 1) >= setting%limit_t) then
        write (0, '(A, F16.6, A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
             '] wavepacket_main: next step simulation time ', setting%delta_t * (state%i + 1), &
             ' a.u., finish wavepacket calculation'
        exit
      end if

      ! Exit loop if matrix replacement occurs in the next step.
      if (setting%is_multistep_input_mode .and. &
           setting%multistep_input_read_interval * state%input_step <= state%t  .and. &
           state%t - setting%delta_t < setting%multistep_input_read_interval * state%input_step) then
         write (0, '(A, F16.6, A, F16.6, A, F16.6, A)') ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, &
              '] wavepacket_main: current step simulation time ', setting%delta_t * state%i, &
              ', next matrix replacement occurs at ', &
              setting%multistep_input_read_interval * state%input_step, &
              ', pause wavepacket calculation'
        exit
      end if

      call make_matrix_step_forward(setting, state%proc, state%input_step, state%t, &
           state%num_atoms, state%atom_indices, state%atom_elements, state%filter_group_indices, &
           state%Y_filtered_desc, state%Y_filtered, state%Y_local, &
           state%H_desc, state%H_multistep1, &
           state%H1_desc, state%H1, state%H1_base, state%H1_multistep, state%H1_multistep1, state%H1_multistep2,  &
           state%S_desc, state%S, state%A_desc, state%A, &
           state%filtered_vecs_desc, state%filtered_vecs, &
           state%full_vecs_desc, state%full_vecs)

      call step_forward_post_process(state%dim, setting, state%t, &
           state%num_atoms, state%atom_indices, state%atom_coordinates, state%atom_elements, &
           state%H_desc, state%H, state%S_desc, state%S, state%Y_filtered_desc, state%Y_filtered, &
           state%H1_desc, state%H1, &
           state%full_vecs_desc, state%full_vecs, state%filtered_vecs_desc, state%filtered_vecs, &
           state%charge_coordinate_mean, state%charge_coordinate_msd, &
           state%tightbinding_energy, state%nonlinear_term_energy, state%total_energy)

      state%i = state%i + 1
      state%t = setting%delta_t * state%i
    end do
  end subroutine wavepacket_main
end module M_wavepacket
