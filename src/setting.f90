module wk_setting_m
  use wk_fson_m
  use wk_fson_value_m
  use wk_fson_path_m
  use wk_global_variables_m
  implicit none

  private
  public :: wk_setting_t, read_alpha_delta_multiple_indices, read_setting, inherit_setting, fill_filtering_setting, &
       verify_setting, print_setting, bcast_setting, print_help

  type wk_setting_t
    ! Temperature: commandline input is in Kelvin and internal value is in atomic unit.
    real(8) :: delta_t = 0.1d0, limit_t = 1.0d0, &
         charge_factor_common, charge_factor_H = -1d0, charge_factor_C = -1d0, &  ! Negative means unset.
         phase_factor_coef, localize_potential_depth = 10d0, &
         temperature, restart_t, &
         alpha_delta_min_x, alpha_delta_max_x, perturb_interval, &
         amplitude_print_threshold = -1d0, amplitude_print_interval, multistep_input_read_interval, &
         vector_cutoff_residual = 0d0, suppress_constant = 0d0, eigenstate_damp_constant = 100d0
    ! Special value of alpha_delta_index: 0 -> HOMO
    integer :: alpha_delta_index = 1, num_steps_per_output_split = 100
    integer, allocatable :: alpha_delta_multiple_indices(:)
    integer :: fst_filter = 0, num_filter  ! fst_filter = 0 means not being set.
    ! These values have different meaning due to init_type.
    integer :: localize_start, localize_end
    integer :: output_interval = 1, restart_step_num, restart_total_states_count, restart_input_step
    integer :: num_multiple_initials = 1
    logical :: is_atom_indices_enabled = .false., to_multiply_phase_factor = .false., &
         is_output_split = .false., &
         is_group_id_used = .false., is_overlap_ignored = .false., is_restart_mode = .false., &
         is_binary_output_mode = .false., is_multistep_input_mode = .false., &
         to_calculate_eigenstate_moment_every_step = .true.
    character(len=1024) :: &
         filename_hamiltonian = '', &
         filename_overlap = '', &
         output_filename = 'out.json', &
         h1_type = 'zero', &
         time_evolution_mode = 'crank_nicolson', &
         init_type = 'alpha_delta', &
         atom_indices_filename = '', &
         alpha_filename = '', &
         lcao_filename = '', &
         group_id_filename = '', &
         restart_filename = '', &
         filter_mode = 'all', &
         filter_group_filename = '', &
         re_initialize_method = 'minimize_lcao_error'
    complex(kind(0d0)), allocatable :: restart_psi(:)
    character(len=1024), allocatable :: restart_split_files_metadata(:)
    real(8), allocatable :: restart_atom_speed(:), restart_atom_perturb(:)
  end type wk_setting_t

contains

  logical function is_digit(c)
    character, intent(in) :: c

    is_digit = c == '0' .or. c == '1' .or. c == '2' .or. c == '3' .or. c == '4' .or. c == '5' .or. &
         c == '6' .or. c == '7' .or. c == '8' .or. c == '9'
  end function is_digit


  subroutine read_charge_factor_for_atoms(argv, charge_factor_common, charge_factor_H, charge_factor_C)
    character(len=1024), intent(inout) :: argv
    real(8), intent(out) :: charge_factor_common, charge_factor_H, charge_factor_C
    integer :: i, j1, j2
    character(len=1024) :: atom_str, charge_factor_str

    if (is_digit(argv(1 : 1))) then
      read(argv, *) charge_factor_common
    else
      i = 1
      do while (i > 0)
        j1 = index(argv(i :), ':')
        if (j1 == 0) then
          stop 'read_charge_factor_for_atoms: invalid charge factor setting format'
        end if
        atom_str = argv(i : i + j1 - 2)
        j2 = index(argv(i :), ',')
        if (j2 > 0) then
          charge_factor_str = argv(i + j1 : i + j2 - 2)
          i = i + j2
        else
          charge_factor_str = argv(i + j1 :)
          i = -1
        end if
        if (trim(atom_str) == 'X') then
          read(charge_factor_str, *) charge_factor_common
        else if (trim(atom_str) == 'H') then
          read(charge_factor_str, *) charge_factor_H
        else if (trim(atom_str) == 'C') then
          read(charge_factor_str, *) charge_factor_C
        else if (trim(atom_str) == 'O') then  ! Temporary implementation: ignore specified value.
        else
          stop 'read_charge_factor_for_atoms: unknown atom'
        end if
      end do
    end if
  end subroutine read_charge_factor_for_atoms


  subroutine read_alpha_delta_multiple_indices(argv, alpha_delta_multiple_indices, num_multiple_initials)
    character(len=*), intent(in) :: argv
    integer, allocatable, intent(out) :: alpha_delta_multiple_indices(:)
    integer, intent(out) :: num_multiple_initials

    integer :: i, pos_rel, pos_abs, buf(10000), j

    pos_abs = 1
    i = 0
    do
      i = i + 1
      pos_rel = index(argv(pos_abs :), ',')
      if (pos_rel == 0) then
        read (argv(pos_abs :), *) buf(i)
        exit
      else
        read (argv(pos_abs : pos_abs + pos_rel - 1), *) buf(i)
        pos_abs = pos_abs + pos_rel
      end if
    end do

    if (allocated(alpha_delta_multiple_indices)) then
      stop
    end if
    allocate(alpha_delta_multiple_indices(i))
    alpha_delta_multiple_indices(1 : i) = buf(1 : i)
    num_multiple_initials = i
  end subroutine read_alpha_delta_multiple_indices


  subroutine read_setting(setting)
    type(wk_setting_t), intent(out) :: setting

    integer :: index_arg, i, j
    character(len=1024) :: argv
    logical :: is_init_type_specified

    is_init_type_specified = .false.
    index_arg = 1
    do while (index_arg <= command_argument_count())  ! Read options.
      call getarg(index_arg, argv)
      if (argv(1:1) == '-') then
        select case (trim(argv(2:)))
        case ('e')
          call getarg(index_arg + 1, setting%time_evolution_mode)
          index_arg = index_arg + 1
        case ('a')
          call getarg(index_arg + 1, setting%atom_indices_filename)
          setting%is_atom_indices_enabled = .true.
          index_arg = index_arg + 1
        case ('i')
          call getarg(index_arg + 1, setting%init_type)
          if (trim(setting%init_type) == 'alpha_delta') then
            call getarg(index_arg + 2, argv)
            if (trim(argv) == 'homo') then
              setting%alpha_delta_index = 0
            else
              read(argv, *) setting%alpha_delta_index
            end if
            index_arg = index_arg + 2
          else if (trim(setting%init_type) == 'alpha_delta_multiple') then
            call getarg(index_arg + 2, argv)
            call read_alpha_delta_multiple_indices(argv, setting%alpha_delta_multiple_indices, &
                 setting%num_multiple_initials)
            index_arg = index_arg + 2
          else if (trim(setting%init_type) == 'alpha_gauss') then
            index_arg = index_arg + 1
          else if (trim(setting%init_type) == 'alpha_file') then
            call getarg(index_arg + 2, setting%alpha_filename)
            index_arg = index_arg + 2
          else if (trim(setting%init_type) == 'lcao_file') then
            call getarg(index_arg + 2, setting%lcao_filename)
            index_arg = index_arg + 2
          else if (trim(setting%init_type) == 'local_alpha_delta' .or. &
               trim(setting%init_type) == 'local_alpha_delta_group') then
            call getarg(index_arg + 2, argv)
            i = index(argv, ',')
            if (i == 0) then
              stop 'missing comma in initial value localization specification (1)'
            else
              j = index(argv(i + 1 :), ',')
              if (j == 0) then
                stop 'missing comma in initial value localization specification (2)'
              else
                read(argv(1 : i - 1), *) setting%localize_start
                read(argv(i + 1 : i + j - 1), *) setting%localize_end
                if (trim(argv(i + j + 1 :)) == 'homo') then
                  setting%alpha_delta_index = 0
                else
                  read(argv(i + j + 1 :), *) setting%alpha_delta_index
                end if
              end if
            end if
            index_arg = index_arg + 2
          else if (trim(setting%init_type) == 'alpha_delta_low_msd') then
            call getarg(index_arg + 2, argv)
            i = index(argv, ',')
            if (i == 0) then
              stop 'missing comma in X value range specification'
            else
              read(argv(: i - 1), *) setting%alpha_delta_min_x
              read(argv(i + 1 :), *) setting%alpha_delta_max_x
              setting%alpha_delta_min_x = setting%alpha_delta_min_x * kAuPerAngstrom
              setting%alpha_delta_max_x = setting%alpha_delta_max_x * kAuPerAngstrom
            end if
            index_arg = index_arg + 2
          else
            stop 'unknown initialization method'
          end if
          is_init_type_specified = .true.
        case ('h')
          call getarg(index_arg + 1, setting%h1_type)
          if (trim(setting%h1_type) == 'diag' .or. trim(setting%h1_type) == 'zero' .or. &
               trim(setting%h1_type) == 'zero_sparse') then
            index_arg = index_arg + 1
          else if (trim(setting%h1_type) == 'zero_damp') then
            call getarg(index_arg + 2, argv)
            read(argv, *) setting%eigenstate_damp_constant
            index_arg = index_arg + 2
          else if (trim(setting%h1_type) == 'charge' .or. trim(setting%h1_type) == 'charge_overlap') then
            call getarg(index_arg + 2, argv)
            call read_charge_factor_for_atoms(argv, setting%charge_factor_common, &
                 setting%charge_factor_H, setting%charge_factor_C)
            index_arg = index_arg + 2
          else if (trim(setting%h1_type) == 'zero_damp_charge_base' .or. &
            trim(setting%h1_type) == 'zero_damp_charge_atom') then
            call getarg(index_arg + 2, argv)
            read(argv, *) setting%eigenstate_damp_constant
            call getarg(index_arg + 3, argv)
            call read_charge_factor_for_atoms(argv, setting%charge_factor_common, &
                 setting%charge_factor_H, setting%charge_factor_C)
            index_arg = index_arg + 3
          else if (trim(setting%h1_type) == 'maxwell') then
            call getarg(index_arg + 2, argv)
            read(argv, *) setting%temperature
            setting%temperature = setting%temperature * kAuPerKelvin
            index_arg = index_arg + 2
          else if (trim(setting%h1_type) == 'harmonic' .or. trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
            call getarg(index_arg + 2, argv)
            i = index(argv, ',')
            if (i == 0) then
              stop 'missing comma in harmonic oscillator type perturbation specification'
            else
              read(argv(1 : i - 1), *) setting%temperature
              read(argv(i + 1 :), *) setting%perturb_interval
              setting%temperature = setting%temperature * kAuPerKelvin
            end if
            index_arg = index_arg + 2
          else
            stop 'unknown H1 function type'
          end if
        case ('f')
          call getarg(index_arg + 1, setting%filter_mode)
          call getarg(index_arg + 2, argv)
          if (trim(setting%filter_mode) == 'all') then  ! format: all <fst>,<end>
            i = index(argv, ',')
            if (i == 0) then
              stop 'missing comma in bulk eigenstate filtering specification'
            else
              read(argv(1 : i - 1), *) setting%fst_filter
              ! Commandline argument is <start>,<end> but internal data is <start> and <num>.
              read(argv(i + 1 :), *) setting%num_filter
              setting%num_filter = setting%num_filter - setting%fst_filter + 1
              if (setting%fst_filter <= 0 .or. setting%num_filter <= 0) then
                stop 'invalid eigenstate filtering specification'
              end if
            end if
          else if (trim(setting%filter_mode) == 'group') then  ! format: group <filename>,<start>,<end>
            i = index(argv, ',')
            if (i == 0) then
              stop 'missing comma in group eigenstate filtering specification (1)'
            else
              j = index(argv(i + 1 :), ',')
              if (j == 0) then
                stop 'missing comma in group eigenstate filtering specification (2)'
              else
                read(argv(1 : i - 1), '(A)') setting%filter_group_filename
                read(argv(i + 1 : i + j - 1), *) setting%fst_filter
                read(argv(i + j + 1 :), *) setting%num_filter  ! Temporary substitution.
                setting%num_filter = setting%num_filter - setting%fst_filter + 1
              end if
            end if
          else
            stop "filter mode must be 'all' or 'group'"
          end if
          index_arg = index_arg + 2
        case ('d')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%delta_t
          if (setting%delta_t <= 0.0d0) then
            stop 'invalid delta_t'
          end if
          index_arg = index_arg + 1
        case ('t')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%limit_t
          if (setting%limit_t <= 0.0d0) then
            stop 'invalid limit_t'
          end if
          index_arg = index_arg + 1
        case ('o')
          call getarg(index_arg + 1, setting%output_filename)
          index_arg = index_arg + 1
        case ('k')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%phase_factor_coef
          setting%to_multiply_phase_factor = .true.
          index_arg = index_arg + 1
        case ('s')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%num_steps_per_output_split
          setting%is_output_split = .true.
          index_arg = index_arg + 1
        case ('g')
          call getarg(index_arg + 1, setting%group_id_filename)
          setting%is_group_id_used = .true.
          index_arg = index_arg + 1
        case ('r')
          call getarg(index_arg + 1, setting%restart_filename)
          setting%is_restart_mode = .true.
          index_arg = index_arg + 1
        case ('b')
          setting%is_binary_output_mode = .true.
        case ('m')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%multistep_input_read_interval
          if (setting%multistep_input_read_interval <= 0.0d0) then
            stop 'invalid multistep_input_read_interval'
          end if
          setting%is_multistep_input_mode = .true.
          index_arg = index_arg + 1
        case ('-block-size')
          call getarg(index_arg + 1, argv)
          read(argv, *) g_wk_block_size
          index_arg = index_arg + 1
        case ('-skip-out')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%output_interval
          index_arg = index_arg + 1
        case ('-local-depth')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%localize_potential_depth
          index_arg = index_arg + 1
        case ('-ignore-overlap')
          setting%is_overlap_ignored = .true.
        case ('-print-amp')
          call getarg(index_arg + 1, argv)
          i = index(argv, ',')
          if (i == 0) then
            stop 'missing comma in amplitude printing setting'
          else
            read(argv(1 : i - 1), *) setting%amplitude_print_threshold
            read(argv(i + 1 :), *) setting%amplitude_print_interval
          end if
          index_arg = index_arg + 1
        case ('-re-init')
          call getarg(index_arg + 1, setting%re_initialize_method)
          index_arg = index_arg + 1
        case ('-cutoff-resid')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%vector_cutoff_residual
          index_arg = index_arg + 1
        case ('-suppress-const')
          call getarg(index_arg + 1, argv)
          read(argv, *) setting%suppress_constant
          index_arg = index_arg + 1
        case default
          stop 'unknown option'
        end select
        index_arg = index_arg + 1  ! Common operation.
      else
        exit
      end if
    end do

    if (.not. setting%is_restart_mode) then
      if (setting%is_overlap_ignored) then
        if (command_argument_count() < index_arg) then
          stop 'specify paths to hamiltonian matrix'
        else
          call getarg(index_arg, setting%filename_hamiltonian)
        end if
      else
        if (command_argument_count() < index_arg + 1) then
          stop 'specify paths to hamiltonian and overlap matrices'
        else
          call getarg(index_arg, setting%filename_hamiltonian)
          call getarg(index_arg + 1, setting%filename_overlap)
        end if
      end if
    end if

    if (.not. is_init_type_specified) then
      write(0, *) '[Warning] initial value type is not specified'
    end if
  end subroutine read_setting


  subroutine inherit_setting(setting)
    type(wk_setting_t), intent(inout) :: setting
    type(fson_value), pointer :: restart_fson, setting_fson, split_files_metadata_fson
    type(fson_value), pointer :: atom_speed_fson, atom_perturb_fson
    type(fson_value), pointer :: last_state_fson, psi_fson, psi_real_fson, psi_imag_fson, value_fson
    integer :: iunit = 17
    integer :: i
    real(8) :: real_value, imag_value

    if (.not. setting%is_restart_mode) then
      allocate(setting%restart_psi(0), setting%restart_split_files_metadata(0))
      allocate(setting%restart_atom_speed(0))
      allocate(setting%restart_atom_perturb(0))
      return
    end if

    if (.not. setting%is_restart_mode .or. trim(setting%restart_filename) == '') then
      stop 'restart setting is not valid'
    end if
    restart_fson => fson_parse(trim(setting%restart_filename), iunit)
    setting_fson => fson_value_get(restart_fson, 'setting')
    split_files_metadata_fson => fson_value_get(restart_fson, 'split_files_metadata')
    last_state_fson => fson_value_get(restart_fson, 'states')
    last_state_fson => fson_value_get(last_state_fson, 0)
    psi_fson => fson_value_get(last_state_fson, 'psi')
    psi_real_fson => fson_value_get(psi_fson, 'real')
    psi_imag_fson => fson_value_get(psi_fson, 'imag')

    ! Read init type.
    call fson_path_get(setting_fson, 'init_type', setting%init_type)
    if (trim(setting%init_type) == 'alpha_delta') then
      call fson_path_get(setting_fson, 'alpha_delta_index', setting%alpha_delta_index)
    else if (trim(setting%init_type) == 'alpha_file') then
      call fson_path_get(setting_fson, 'alpha_filename', setting%alpha_filename)
    else if (trim(setting%init_type) == 'lcao_file') then
      call fson_path_get(setting_fson, 'lcao_filename', setting%lcao_filename)
    else if (trim(setting%init_type) == 'local_alpha_delta' .or. &
         trim(setting%init_type) == 'local_alpha_delta_group') then
      call fson_path_get(setting_fson, 'alpha_delta_index', setting%alpha_delta_index)
      call fson_path_get(setting_fson, 'localize_start', setting%localize_start)
      call fson_path_get(setting_fson, 'localize_end', setting%localize_end)
      call fson_path_get(setting_fson, 'localize_potential_depth', setting%localize_potential_depth)
    end if
    ! Read h1 type.  ! TODO: support zero_damp_charge_base etc.
    call fson_path_get(setting_fson, 'h1_type', setting%h1_type)
    if (trim(setting%h1_type) == 'charge' .or. trim(setting%h1_type) == 'charge_overlap') then
      call fson_path_get(setting_fson, 'charge_factor_common', setting%charge_factor_common)
      call fson_path_get(setting_fson, 'charge_factor_H', setting%charge_factor_H)
      call fson_path_get(setting_fson, 'charge_factor_C', setting%charge_factor_C)
    else if (trim(setting%h1_type) == 'maxwell') then
      call fson_path_get(setting_fson, 'temperature', setting%temperature)
      setting%temperature = setting%temperature * kAuPerKelvin
    else if (trim(setting%h1_type) == 'harmonic' .or. trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
      call fson_path_get(setting_fson, 'temperature', setting%temperature)
      setting%temperature = setting%temperature * kAuPerKelvin
      call fson_path_get(setting_fson, 'perturb_interval', setting%perturb_interval)
    end if
    call fson_path_get(setting_fson, 'time_evolution_mode', setting%time_evolution_mode)
    call fson_path_get(setting_fson, 'is_atom_indices_enabled', setting%is_atom_indices_enabled)
    if (setting%is_atom_indices_enabled) then
          call fson_path_get(setting_fson, 'atom_indices_filename', setting%atom_indices_filename)
    end if
    call fson_path_get(setting_fson, 'fst_filter', setting%fst_filter)
    call fson_path_get(setting_fson, 'end_filter', setting%num_filter)  ! Temporary substitution.
    setting%num_filter = setting%num_filter - setting%fst_filter + 1
    call fson_path_get(setting_fson, 'delta_t', setting%delta_t)
    ! limit_t should not be inherited.
    call fson_path_get(setting_fson, 'output_filename', setting%output_filename)
    call fson_path_get(setting_fson, 'to_multiply_phase_factor', setting%to_multiply_phase_factor)
    if (setting%to_multiply_phase_factor) then
      call fson_path_get(setting_fson, 'phase_factor_coef', setting%phase_factor_coef)
    end if
    call fson_path_get(setting_fson, 'is_output_split', setting%is_output_split)
    if (setting%is_output_split) then
      call fson_path_get(setting_fson, 'num_steps_per_output_split', setting%num_steps_per_output_split)
    end if
    call fson_path_get(setting_fson, 'is_group_id_used', setting%is_group_id_used)
    if (setting%is_group_id_used) then
      call fson_path_get(setting_fson, 'group_id_filename', setting%group_id_filename)
    end if
    call fson_path_get(setting_fson, 'output_interval', setting%output_interval)
    call fson_path_get(setting_fson, 'is_binary_output_mode', setting%is_binary_output_mode)
    call fson_path_get(setting_fson, 'is_overlap_ignored', setting%is_overlap_ignored)  !
    call fson_path_get(setting_fson, 'filename_hamiltonian', setting%filename_hamiltonian)
    if (.not. setting%is_overlap_ignored) then
      call fson_path_get(setting_fson, 'filename_overlap', setting%filename_overlap)
    end if
    call fson_path_get(setting_fson, 'is_multistep_input_mode', setting%is_multistep_input_mode)
    if (setting%is_multistep_input_mode) then
      call fson_path_get(setting_fson, 'multistep_input_read_interval', setting%multistep_input_read_interval)
    end if

    allocate(setting%restart_psi(fson_value_count(psi_real_fson)), &
         setting%restart_split_files_metadata(fson_value_count(split_files_metadata_fson)))
    do i = 1, fson_value_count(psi_real_fson)
      value_fson => fson_value_get(psi_real_fson, i)
      call fson_path_get(value_fson, '', real_value)
      value_fson => fson_value_get(psi_imag_fson, i)
      call fson_path_get(value_fson, '', imag_value)
      setting%restart_psi(i) = cmplx(real_value, imag_value, kind(0d0))
    end do
    do i = 1, fson_value_count(split_files_metadata_fson)
      value_fson => fson_value_get(split_files_metadata_fson, i)
      call fson_path_get(value_fson, '', setting%restart_split_files_metadata(i))
    end do

    if (trim(setting%h1_type) == 'zero') then  ! Nothing to inherit.
    else if (trim(setting%h1_type) == 'maxwell') then
      atom_speed_fson => fson_value_get(last_state_fson, 'atom_speed')
      allocate(setting%restart_atom_speed(fson_value_count(atom_speed_fson)))
      do i = 1, fson_value_count(atom_speed_fson)
        value_fson => fson_value_get(atom_speed_fson, i)
        call fson_path_get(value_fson, '', setting%restart_atom_speed(i))
      end do
    else if (trim(setting%h1_type) == 'harmonic' .or. trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
      atom_perturb_fson => fson_value_get(last_state_fson, 'atom_perturb')
      allocate(setting%restart_atom_perturb(fson_value_count(atom_perturb_fson)))
      do i = 1, fson_value_count(atom_perturb_fson)
        value_fson => fson_value_get(atom_perturb_fson, i)
        call fson_path_get(value_fson, '', setting%restart_atom_perturb(i))
      end do
    else
      stop 'restart is available only for -h zero, -h maxwell and -h harmonic'
    end if

    call fson_path_get(last_state_fson, 'step_num', setting%restart_step_num)
    call fson_path_get(last_state_fson, 'time', setting%restart_t)
    call fson_path_get(setting_fson, 'restart_total_states_count', setting%restart_total_states_count)
    call fson_path_get(setting_fson, 'restart_input_step', setting%restart_input_step)
    call fson_path_get(setting_fson, 're_initialize_method', setting%re_initialize_method)
    call fson_value_destroy(restart_fson)
  end subroutine inherit_setting


  subroutine fill_filtering_setting(dim, num_groups, setting)
    integer, intent(in) :: dim, num_groups
    type(wk_setting_t), intent(inout) :: setting

    if (trim(setting%filter_mode) == 'all' .and. setting%fst_filter == 0) then  ! -f option was not specified.
      setting%fst_filter = 1
      setting%num_filter = dim
    else if (trim(setting%filter_mode) == 'group') then
      setting%fst_filter = 1
      setting%num_filter = dim / num_groups
    end if
  end subroutine fill_filtering_setting


  subroutine verify_setting(dim, setting)
    integer, intent(in) :: dim
    type(wk_setting_t), intent(in) :: setting

    if (setting%time_evolution_mode /= 'crank_nicolson' .and. &
         setting%time_evolution_mode /= 'taylor1' .and. &
         setting%time_evolution_mode /= 'taylor2' .and. &
         setting%time_evolution_mode /= 'taylor3') then
      print *, 'a',trim(setting%time_evolution_mode),'b'
      stop 'invalid time evolution scheme'
    end if

    if (setting%to_multiply_phase_factor .and. &
         .not. setting%is_atom_indices_enabled) then
      stop '-k option must be used with -a option'
    end if

    if (setting%is_output_split .and. trim(setting%output_filename) == '-') then
      stop 'output to stdout cannot be split'
    end if

    !if (setting%fst_filter < 1 .or. setting%fst_filter + setting%num_filter - 1 > dim) then
    !  stop 'eigenvector filtering is out of index'
    !end if

    if (trim(setting%init_type) == 'local_alpha_delta' .and. &
         setting%to_multiply_phase_factor) then
      stop '-i local_alpha_delta and -k must not be set simultaneously'
    end if

    if (setting%alpha_delta_index /= 0 .and. setting%alpha_delta_index < 0) then
      stop 'index of eigenstate for initial value must be positive integer or "homo"'
    end if

    if (setting%localize_start > setting%localize_end) then
      stop 'localized indices are empty'
    end if

    if (trim(setting%group_id_filename) /= '' .and. .not. setting%is_atom_indices_enabled) then
      stop 'group id file (-g) must be used with XYZ file (-a)'
    end if

    if (trim(setting%init_type) == 'local_alpha_delta_group' .and. &
         trim(setting%group_id_filename) == '') then
      stop 'when -i local_alpha_delta_group is used, group id file (-g) also must be specified'
    end if

    if (setting%is_restart_mode) then
      if (size(setting%restart_psi) /= dim) then
        stop 'dimension of the psi in the restart file is not equal to the matrix dimension'
      end if
    end if

    if (setting%is_binary_output_mode .and. .not. setting%is_output_split) then
      stop 'binary output mode must be used with split output mode'
    end if

    if (setting%is_output_split .and. trim(setting%output_filename) == '-') then
      stop 'split output mode cannot be used with standard output mode'
    end if

    !if (setting%delta_t >= setting%multistep_input_read_interval) then
    !  stop 'time step size must be smaller than the interval of multistep input'
    !end if

    if (.not. &
         (trim(setting%re_initialize_method) == 'minimize_lcao_error' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_cutoff' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_suppress' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_orthogonal' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_adaptive' .or. &
         trim(setting%re_initialize_method) == 'minimize_lcao_error_matrix_suppress_select' .or. &
         trim(setting%re_initialize_method) == 'minimize_alpha_error')) then
      stop 'unknown state re-initialization method'
    end if

    if (setting%num_multiple_initials > 1 .and. .not. &
         (trim(setting%h1_type) == 'zero' .or. &
         trim(setting%h1_type) == 'zero_damp' .or. &
         trim(setting%h1_type) == 'zero_sparse' .or. &
         trim(setting%h1_type) == 'maxwell' .or. &
         trim(setting%h1_type) == 'harmonic' .or. &
         trim(setting%h1_type) == 'harmonic_for_nn_exciton')) then
      stop 'specified h1 type does not support multiple initials'
    end if

    if (setting%num_multiple_initials > 1 .and. .not. setting%is_output_split) then
      stop 'multiple initials mode must be used with output split'
    end if

    if (trim(setting%h1_type) == 'zero_damp' .and. setting%eigenstate_damp_constant < 0d0) then
      stop 'eigenstate damp constant must be positive'
    end if
  end subroutine verify_setting


  subroutine print_setting(setting)
    type(wk_setting_t), intent(in) :: setting

    integer :: index_arg
    character(len=1024) :: argv

    print *, 'version: ', kVersion

    do index_arg = 0, command_argument_count()
      call getarg(index_arg, argv)
      if (index_arg < command_argument_count()) then
        write(6, '(2A)', advance='no') trim(argv), ' '
      else
        write(6, '(A)') trim(argv)
      end if
    end do

    print *, 'h1_type: ', trim(setting%h1_type)
    if (trim(setting%h1_type) == 'zero_damp') then
      print '(A, E26.16e3)', ' eigenstate_damp_constant: ', setting%eigenstate_damp_constant
    else if (trim(setting%h1_type(1 : 6)) == 'charge') then
      print '(A, E26.16e3)', ' charge_factor_common: ', setting%charge_factor_common
      print '(A, E26.16e3)', ' charge_factor_H: ', setting%charge_factor_H
      print '(A, E26.16e3)', ' charge_factor_C: ', setting%charge_factor_C
    else if (trim(setting%h1_type) == 'zero_damp_charge_base' .or. &
      trim(setting%h1_type) == 'zero_damp_charge_atom') then
      print '(A, E26.16e3)', ' eigenstate_damp_constant: ', setting%eigenstate_damp_constant
      print '(A, E26.16e3)', ' charge_factor_common: ', setting%charge_factor_common
      print '(A, E26.16e3)', ' charge_factor_H: ', setting%charge_factor_H
      print '(A, E26.16e3)', ' charge_factor_C: ', setting%charge_factor_C
    else if (trim(setting%h1_type) == 'maxwell') then
      print '(A, E26.16e3)', ' temperature: ', setting%temperature / kAuPerKelvin
      print '(A, E26.16e3)', ' kAccelRatio: ', kAccelRatio
    else if (trim(setting%h1_type) == 'harmonic') then
      print '(A, E26.16e3)', ' temperature: ', setting%temperature / kAuPerKelvin
      print '(A, E26.16e3)', ' perturb_interval: ', setting%perturb_interval
    else if (trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
      print '(A, E26.16e3)', ' temperature: ', setting%temperature / kAuPerKelvin
      print '(A, E26.16e3)', ' perturb_interval: ', setting%perturb_interval
    end if
    print *, 'eigenstate_filtering: ', setting%fst_filter, '-', &
         setting%fst_filter + setting%num_filter - 1
    print *, 'init_type: ', trim(setting%init_type)
    if (trim(setting%init_type) == 'alpha_delta') then
      print *, 'alpha_delta_index: ', setting%alpha_delta_index
    else if (trim(setting%init_type) == 'alpha_file') then
      print *, 'alpha_filename: ', trim(setting%alpha_filename)
    else if (trim(setting%init_type) == 'lcao_file') then
      print *, 'lcao_filename: ', trim(setting%lcao_filename)
    end if
    print *, 'filename_hamiltonian: ', trim(setting%filename_hamiltonian)
    print *, 'filename_overlap ', trim(setting%filename_overlap)
    print '(A, E26.16e3)', ' delta_t: ', setting%delta_t
    print '(A, E26.16e3)', ' limit_t: ', setting%limit_t
    if (setting%is_group_id_used) then
      print *, 'is_group_id_used: true'
      print *, 'group_id_filename: ', trim(setting%group_id_filename)
    else
      print *, 'is_group_id_used: false'
    end if
    if (setting%is_overlap_ignored) then
      print *, 'is_overlap_ignored: true'
    else
      print *, 'is_overlap_ignored: false'
    end if
    print *, 're_initialize_method: ', trim(setting%re_initialize_method)
    if (setting%re_initialize_method == 'minimize_lcao_error_cutoff') then
      print *, 'vector_cutoff_residual: ', setting%vector_cutoff_residual
    else if (setting%re_initialize_method == 'minimize_lcao_error_suppress') then
      print *, 'vector_suppress_constant: ', setting%suppress_constant
    else if (setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress') then
      print *, 'matrix_suppress_constant: ', setting%suppress_constant
    else if (setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress_orthogonal') then
      print *, 'matrix_suppress_constant: ', setting%suppress_constant
    else if (setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress_adaptive') then
      print *, 'adaptive_matrix_suppress_constant: ', setting%suppress_constant
    else if (setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress_select') then
      print *, 'select_matrix_suppress_constant: ', setting%suppress_constant
    end if
    print *, 'output_filename: ', trim(setting%output_filename)
    if (setting%is_atom_indices_enabled) then
      print *, 'is_atom_indices_enabled: true'
      print *, 'atom_indices_filename: ', trim(setting%atom_indices_filename)
    else
      print *, 'is_atom_indices_enabled: false'
    end if
    if (setting%to_multiply_phase_factor) then
      print *, 'to_multiply_phase_factor: true'
      print '(A, E26.16e3)', ' phase_factor_coef: ', setting%phase_factor_coef
    else
      print *, 'to_multiply_phase_factor: false'
    end if
    if (setting%is_output_split) then
      print *, 'is_output_split: true'
      print '(A, I0)', ' num_steps_per_output_split: ', setting%num_steps_per_output_split
    else
      print *, 'is_output_split: false'
    end if
    if (setting%is_binary_output_mode) then
      print *, 'is_binary_output_mode: true'
    else
      print *, 'is_binary_output_mode: false'
    end if
    if (setting%amplitude_print_threshold >= 0d0) then
      print *, 'to_print_amplitude: true'
      print *, 'amplitude_print_threshold: ', setting%amplitude_print_threshold
      print *, 'amplitude_print_interval: ', setting%amplitude_print_interval
    else
      print *, 'to_print_amplitude: false'
    end if
    print *, 'num_multiple_initials: ', setting%num_multiple_initials
    if (allocated(setting%alpha_delta_multiple_indices)) then
      print *, 'alpha_delta_multiple_indices: ', setting%alpha_delta_multiple_indices(:)
    end if
  end subroutine print_setting


  subroutine bcast_setting(root, setting)
    use mpi
    integer, intent(in) :: root
    type(wk_setting_t), intent(inout) :: setting
    integer, parameter :: num_real = 18, num_integer = 17, num_logical = 8, num_character = 16
    real(8) :: buf_real(num_real)
    integer :: buf_integer(num_integer), my_rank, ierr, size_psi, size_split_files_metadata, size_atom_speed, size_atom_perturb
    logical :: buf_logical(num_logical)
    character(len=1024), allocatable :: buf_character(:)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (my_rank == root) then
      buf_real(1) = setting%delta_t
      buf_real(2) = setting%limit_t
      buf_real(3) = setting%charge_factor_common
      buf_real(4) = setting%phase_factor_coef
      buf_real(5) = setting%localize_potential_depth
      buf_real(6) = setting%temperature
      buf_real(7) = setting%restart_t
      buf_real(8) = setting%alpha_delta_min_x
      buf_real(9) = setting%alpha_delta_max_x
      buf_real(10) = setting%perturb_interval
      buf_real(11) = setting%amplitude_print_threshold
      buf_real(12) = setting%amplitude_print_interval
      buf_real(13) = setting%multistep_input_read_interval
      buf_real(14) = setting%charge_factor_H
      buf_real(15) = setting%charge_factor_C
      buf_real(16) = setting%vector_cutoff_residual
      buf_real(17) = setting%suppress_constant
      buf_real(18) = setting%eigenstate_damp_constant
      buf_integer(1) = setting%alpha_delta_index
      buf_integer(2) = setting%num_steps_per_output_split
      buf_integer(3) = setting%fst_filter
      buf_integer(4) = setting%num_filter
      buf_integer(5) = setting%localize_start
      buf_integer(6) = setting%localize_end
      buf_integer(7) = setting%output_interval
      buf_integer(8) = setting%restart_step_num
      buf_integer(9) = setting%restart_total_states_count
      buf_integer(10) = setting%restart_input_step
      if (setting%is_restart_mode) then
        buf_integer(11) = size(setting%restart_psi)
        buf_integer(12) = size(setting%restart_split_files_metadata)
        if (trim(setting%h1_type) == 'maxwell') then
          buf_integer(13) = size(setting%restart_atom_speed)
        else
          buf_integer(13) = 0
        end if
        if (trim(setting%h1_type) == 'harmonic' .or. trim(setting%h1_type) == 'harmonic_for_nn_exciton') then
          buf_integer(14) = size(setting%restart_atom_perturb)
        else
          buf_integer(14) = 0
        end if
      else
        buf_integer(11 : 14) = 0
      end if
      buf_integer(15) = -1  ! Dummy value.
      buf_integer(16) = setting%num_multiple_initials
      if (allocated(setting%alpha_delta_multiple_indices)) then
        buf_integer(17) = size(setting%alpha_delta_multiple_indices, 1)
      else
        buf_integer(17) = 0
      end if
      buf_logical(1) = setting%is_atom_indices_enabled
      buf_logical(2) = setting%to_multiply_phase_factor
      buf_logical(3) = setting%is_output_split
      buf_logical(4) = .true.  ! Dummy value.
      buf_logical(5) = setting%is_group_id_used
      buf_logical(6) = setting%is_overlap_ignored
      buf_logical(7) = setting%is_restart_mode
      buf_logical(8) = setting%is_multistep_input_mode
    end if
    call mpi_bcast(buf_real, num_real, mpi_double_precision, root, mpi_comm_world, ierr)
    call mpi_bcast(buf_integer, num_integer, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(buf_logical, num_logical, mpi_logical, root, mpi_comm_world, ierr)
    size_psi = buf_integer(11)
    size_split_files_metadata = buf_integer(12)
    size_atom_speed = buf_integer(13)
    size_atom_perturb = buf_integer(14)
    allocate(buf_character(num_character + size_split_files_metadata))
    if (my_rank == root) then
      buf_character(1) = setting%filename_hamiltonian
      buf_character(2) = setting%filename_overlap
      buf_character(3) = setting%output_filename
      buf_character(4) = setting%h1_type
      buf_character(5) = setting%time_evolution_mode
      buf_character(6) = setting%init_type
      buf_character(7) = setting%atom_indices_filename
      buf_character(8) = setting%alpha_filename
      buf_character(9) = setting%lcao_filename
      buf_character(10) = ''  ! Dummy value.
      buf_character(11) = ''
      buf_character(12) = setting%group_id_filename
      buf_character(13) = setting%restart_filename
      buf_character(14) = setting%filter_mode
      buf_character(15) = setting%filter_group_filename
      buf_character(16) = setting%re_initialize_method
      buf_character(num_character + 1 : num_character + size_split_files_metadata) = &
           setting%restart_split_files_metadata(1 : size_split_files_metadata)
    else
      allocate(setting%restart_psi(size_psi), setting%restart_split_files_metadata(size_split_files_metadata), &
           setting%restart_atom_speed(size_atom_speed), setting%restart_atom_perturb(size_atom_perturb))
      if (buf_integer(17) > 0) then
        allocate(setting%alpha_delta_multiple_indices(buf_integer(17)))
      end if
    end if
    call mpi_bcast(setting%restart_psi, size_psi, mpi_double_complex, root, mpi_comm_world, ierr)
    call mpi_bcast(buf_character, 1024 * (num_character + size_split_files_metadata), mpi_character, &
         root, mpi_comm_world, ierr)
    call mpi_bcast(setting%restart_atom_speed, size_atom_speed, mpi_double_precision, root, mpi_comm_world, ierr)
    call mpi_bcast(setting%restart_atom_perturb, size_atom_perturb, mpi_double_precision, root, mpi_comm_world, ierr)
    if (buf_integer(17) > 0) then
      call mpi_bcast(setting%alpha_delta_multiple_indices, buf_integer(17), mpi_integer, root, mpi_comm_world, ierr)
    end if
    if (my_rank /= root) then
      setting%delta_t = buf_real(1)
      setting%limit_t = buf_real(2)
      setting%charge_factor_common = buf_real(3)
      setting%phase_factor_coef = buf_real(4)
      setting%localize_potential_depth = buf_real(5)
      setting%temperature = buf_real(6)
      setting%restart_t = buf_real(7)
      setting%alpha_delta_min_x = buf_real(8)
      setting%alpha_delta_max_x = buf_real(9)
      setting%perturb_interval = buf_real(10)
      setting%amplitude_print_threshold = buf_real(11)
      setting%amplitude_print_interval = buf_real(12)
      setting%multistep_input_read_interval = buf_real(13)
      setting%charge_factor_H = buf_real(14)
      setting%charge_factor_C = buf_real(15)
      setting%vector_cutoff_residual = buf_real(16)
      setting%suppress_constant = buf_real(17)
      setting%eigenstate_damp_constant = buf_real(18)
      setting%alpha_delta_index = buf_integer(1)
      setting%num_steps_per_output_split = buf_integer(2)
      setting%fst_filter = buf_integer(3)
      setting%num_filter = buf_integer(4)
      setting%localize_start = buf_integer(5)
      setting%localize_end = buf_integer(6)
      setting%output_interval = buf_integer(7)
      setting%restart_step_num = buf_integer(8)
      setting%restart_total_states_count = buf_integer(9)
      setting%restart_input_step = buf_integer(10)
      ! buf_integer(15) is not used.
      setting%num_multiple_initials = buf_integer(16)
      setting%is_atom_indices_enabled = buf_logical(1)
      setting%to_multiply_phase_factor = buf_logical(2)
      setting%is_output_split = buf_logical(3)
      ! buf_logical(4) is not used.
      setting%is_group_id_used = buf_logical(5)
      setting%is_overlap_ignored = buf_logical(6)
      setting%is_restart_mode = buf_logical(7)
      setting%is_multistep_input_mode = buf_logical(8)
      setting%filename_hamiltonian  = buf_character(1)
      setting%filename_overlap = buf_character(2)
      setting%output_filename = buf_character(3)
      setting%h1_type = buf_character(4)
      setting%time_evolution_mode = buf_character(5)
      setting%init_type  = buf_character(6)
      setting%atom_indices_filename = buf_character(7)
      setting%alpha_filename = buf_character(8)
      setting%lcao_filename = buf_character(9)
      ! buf_character(10) and (11) are not used.
      setting%group_id_filename = buf_character(12)
      setting%restart_filename = buf_character(13)
      setting%filter_mode = buf_character(14)
      setting%filter_group_filename = buf_character(15)
      setting%re_initialize_method = buf_character(16)
      setting%restart_split_files_metadata(1 : size_split_files_metadata) = &
           buf_character(num_character + 1 : num_character + size_split_files_metadata)
    end if
  end subroutine bcast_setting


  subroutine print_help()
    print *, 'Usage: ./wavepacket [OPTION]... FILE_H FILE_S'
    print *, '       FILE_H/S is path to Hamiltonian/overlap matrix in MatrixMarket format'
    print *, 'Options are:'
    print *, '  -a <file>  Specify xyz file name'
    print *, '  -i [alpha_gauss|alpha_delta <n>|alpha_file <file>|lcao_file <file>]'
    print *, '  (cont.) -i [local_alpha_delta <i>,<j>,<k>|local_alpha_delta_group <i>,<j>,<k>]', &
         'Select state initialization method (Default: alpha_delta 1)'
    print *, '  -h [zero|zero_damp <float>|diag|charge <float>|charge_overlap <float>|maxwell <float>'
    print *, '      zero_damp_charge_base <float> <float>|zero_damp_charge_atom <float> <float>]'
    print *, '     Select how to generate H1, perturbation term (Default: zero)'
    print *, '  -e [crank_nicolson|taylor1|taylor2|taylor3]  Select time evolution mode (Default: crank_nicolson)'
    print *, '  -f <m>,<n>  Enable eigenstate filtering using from m-th to n-th (inclusive) generalized eigenvectors'
    print *, '  -d <float>  Specify time step width in atomic unit (Default: 0.1)'
    print *, '  -t <float>  Specify time simulation ends in atomic unit (Default: 1.0)'
    print *, '  -m <float>  Enable multiple input step mode and specify its time invertal in atomic unit (Default: disabled)'
    print *, '  -k <float>  Multiply phase factor to LCAO represented ', &
         'initial state and specify wavenumber (for x axis only). Must be used with -a'
    print *, '  -g <file>  Specify group id file, which is used in -i local_alpha_delta_group'
    print *, '  -o <file>  Specify output file name. "-" means stdout (Default: out.txt)'
    print *, '  -s <n>  Split output files per <n> times steps'
    print *, '  -r <file>  Restart calculation using <file>'
    print *, '  -b Output array of real values in binary mode'
    print *, '  --ignore-overlap  Use identity matrix for overlap matrix'
    print *, '  --no-replace-basis  With -m option, treat changes of matrices as perturbation without replacing eigenvectors'
    print *, '  --block-size <n>  Set default block size to <n> (Default: 64)'
    print *, '  --skip-out <n>  Set output printing interval to <n> (Default: 1)'
    print *, '  --local-depth <float>  Specify depth of potential well in atomic unit (Default: 10.0)'
    print *, '  --print-amp <threshold>,<interval>  ', &
         'Determine the threshold and interval to print high amplitude elements of A (Default: no print)'
    print *, '  --re-init [minimize_lcao_error | minimize_alpha_error]'
    print *, '  --cutoff-resid <float>  vector cutoff residual'
  end subroutine print_help
end module wk_setting_m
