module wp_output_m
  use mpi
  use wp_atom_m
  use wp_charge_m
  use wp_descriptor_parameters_m
  use wp_fson_m
  use wp_fson_value_m
  use wp_fson_string_m
  use wp_event_logger_m
  use wp_processes_m
  use wp_global_variables_m
  use wp_linear_algebra_m
  use wp_setting_m
  use wp_state_m
  use wp_util_m
  implicit none

  private
  public :: add_setting_json, add_condition_json, add_state_json, add_structure_json

contains

  subroutine add_setting_json(setting, proc, multiple_initial_index, output)
    type(wp_setting_t), intent(in) :: setting
    type(wp_process_t), intent(in) :: proc
    integer, intent(in) :: multiple_initial_index
    type(fson_value), pointer, intent(inout) :: output

    type(fson_value), pointer :: setting_in_fson, setting_elem
    type(fson_string), pointer :: str
    character(len=1024) :: argv
    integer :: index_arg

    setting_in_fson => fson_value_create()
    setting_in_fson%value_type = TYPE_OBJECT
    call fson_set_name('setting', setting_in_fson)

    ! Set version.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_OBJECT
    call fson_set_name('version', setting_elem)
    call fson_set_as_string(kVersion, setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set command.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_STRING
    call fson_set_name('command', setting_elem)
    str => fson_string_create()
    do index_arg = 0, command_argument_count()
      call getarg(index_arg, argv)
      call fson_string_append(str, trim(argv))
      if (index_arg < command_argument_count()) then
        call fson_string_append(str, ' ')
      end if
    end do
    setting_elem%value_string => str
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set h1_type.
    setting_elem => fson_value_create()
    call fson_set_name('h1_type', setting_elem)
    call fson_set_as_string(trim(setting%h1_type), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)
    if (trim(setting%h1_type(1 : 6)) == 'charge' .or. &
         trim(setting%h1_type) == 'zero_damp_charge_base' .or. &
         trim(setting%h1_type) == 'zero_damp_charge_atom') then
      ! charge_factor_common
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('charge_factor_common', setting_elem)
      setting_elem%value_real = setting%charge_factor_common
      call fson_value_add(setting_in_fson, setting_elem)
      ! charge_factor_H
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('charge_factor_H', setting_elem)
      setting_elem%value_real = setting%charge_factor_H
      call fson_value_add(setting_in_fson, setting_elem)
      ! charge_factor_C
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('charge_factor_C', setting_elem)
      setting_elem%value_real = setting%charge_factor_C
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%h1_type) == 'maxwell') then
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('temperature', setting_elem)
      setting_elem%value_real = setting%temperature / kAuPerKelvin
      call fson_value_add(setting_in_fson, setting_elem)
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('kAccelRatio', setting_elem)
      setting_elem%value_real = kAccelRatio
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%h1_type) == 'harmonic') then
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('temperature', setting_elem)
      setting_elem%value_real = setting%temperature / kAuPerKelvin
      call fson_value_add(setting_in_fson, setting_elem)
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('perturb_interval', setting_elem)
      setting_elem%value_real = setting%perturb_interval
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%h1_type(1 : 9)) == 'zero_damp') then
      setting_elem => fson_value_create()
      setting_elem%value_type = TYPE_REAL
      call fson_set_name('eigenstate_damp_constant', setting_elem)
      setting_elem%value_real = setting%eigenstate_damp_constant
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set setting%time_evolution_mode
    setting_elem => fson_value_create()
    call fson_set_name('time_evolution_mode', setting_elem)
    call fson_set_as_string(trim(setting%time_evolution_mode), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set eigenstate_filtering.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('fst_filter', setting_elem)
    setting_elem%value_integer = setting%fst_filter
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('end_filter', setting_elem)
    setting_elem%value_integer = setting%fst_filter + setting%num_filter - 1
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set init_type.
    setting_elem => fson_value_create()
    call fson_set_name('init_type', setting_elem)
    call fson_set_as_string(trim(setting%init_type), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)
    if (trim(setting%init_type) == 'alpha_delta') then
      setting_elem => fson_value_create()
      call fson_set_name('alpha_delta_index', setting_elem)
      setting_elem%value_type = TYPE_INTEGER
      setting_elem%value_integer = setting%alpha_delta_index
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%init_type) == 'alpha_delta_multiple') then
      setting_elem => fson_value_create()
      call fson_set_name('alpha_delta_multiple_index', setting_elem)
      setting_elem%value_type = TYPE_INTEGER
      setting_elem%value_integer = setting%alpha_delta_multiple_indices(multiple_initial_index)
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%init_type) == 'alpha_file') then
      setting_elem => fson_value_create()
      call fson_set_name('alpha_filename', setting_elem)
      call fson_set_as_string(trim(setting%alpha_filename), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%init_type) == 'lcao_file') then
      setting_elem => fson_value_create()
      call fson_set_name('lcao_filename', setting_elem)
      call fson_set_as_string(trim(setting%lcao_filename), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
    else if (trim(setting%init_type) == 'local_alpha_delta' .or. &
         trim(setting%init_type) == 'local_alpha_delta_group') then
      setting_elem => fson_value_create()
      call fson_set_name('alpha_delta_index', setting_elem)
      setting_elem%value_type = TYPE_INTEGER
      setting_elem%value_integer = setting%alpha_delta_index
      call fson_value_add(setting_in_fson, setting_elem)
      setting_elem => fson_value_create()
      call fson_set_name('localize_start', setting_elem)
      setting_elem%value_type = TYPE_INTEGER
      setting_elem%value_integer = setting%localize_start
      call fson_value_add(setting_in_fson, setting_elem)
      setting_elem => fson_value_create()
      call fson_set_name('localize_end', setting_elem)
      setting_elem%value_type = TYPE_INTEGER
      setting_elem%value_integer = setting%localize_end
      call fson_value_add(setting_in_fson, setting_elem)
      setting_elem => fson_value_create()
      call fson_set_name('localize_potential_depth', setting_elem)
      setting_elem%value_type = TYPE_REAL
      setting_elem%value_real = setting%localize_potential_depth
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set filenames.
    setting_elem => fson_value_create()
    call fson_set_name('filename_hamiltonian', setting_elem)
    call fson_set_as_string(trim(setting%filename_hamiltonian), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_overlap_ignored', setting_elem)
    setting_elem%value_logical = setting%is_overlap_ignored
    call fson_value_add(setting_in_fson, setting_elem)
    if (.not. setting%is_overlap_ignored) then
      setting_elem => fson_value_create()
      call fson_set_name('filename_overlap', setting_elem)
      call fson_set_as_string(trim(setting%filename_overlap), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
    end if
    setting_elem => fson_value_create()
    call fson_set_name('output_filename', setting_elem)
    call fson_set_as_string(trim(setting%output_filename), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set multistep input mode related settings.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_multistep_input_mode', setting_elem)
    setting_elem%value_logical = setting%is_multistep_input_mode
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%is_multistep_input_mode) then
      setting_elem => fson_value_create()
      call fson_set_name('multistep_input_read_interval', setting_elem)
      setting_elem%value_type = TYPE_REAL
      setting_elem%value_real = setting%multistep_input_read_interval
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set delta_t and limit_t.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_REAL
    call fson_set_name('delta_t', setting_elem)
    setting_elem%value_real = setting%delta_t
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_REAL
    call fson_set_name('limit_t', setting_elem)
    setting_elem%value_real = setting%limit_t
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set atom index related settings.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_atom_indices_enabled', setting_elem)
    setting_elem%value_logical = setting%is_atom_indices_enabled
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%is_atom_indices_enabled) then
      setting_elem => fson_value_create()
      call fson_set_name('atom_indices_filename', setting_elem)
      call fson_set_as_string(trim(setting%atom_indices_filename), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set phase factor related settings.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('to_multiply_phase_factor', setting_elem)
    setting_elem%value_logical = setting%to_multiply_phase_factor
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%to_multiply_phase_factor) then
      setting_elem => fson_value_create()
      call fson_set_name('phase_factor_coef', setting_elem)
      setting_elem%value_type = TYPE_REAL
      setting_elem%value_real = setting%phase_factor_coef
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set output file splitting, interval and binalize related settings.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_output_split', setting_elem)
    setting_elem%value_logical = setting%is_output_split
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%is_output_split) then
      setting_elem => fson_value_create()
      call fson_set_name('num_steps_per_output_split', setting_elem)
      setting_elem%value_type = TYPE_INTEGER
      setting_elem%value_integer = setting%num_steps_per_output_split
      call fson_value_add(setting_in_fson, setting_elem)
    end if
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('output_interval', setting_elem)
    setting_elem%value_integer = setting%output_interval
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_binary_output_mode', setting_elem)
    setting_elem%value_logical = setting%is_binary_output_mode
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set default block size.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('block_size', setting_elem)
    setting_elem%value_integer = g_wp_block_size
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set process information.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('num_mpi_processes', setting_elem)
    setting_elem%value_integer = proc%n_procs
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('num_mpi_processes_row', setting_elem)
    setting_elem%value_integer = proc%n_procs_row
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('num_mpi_processes_col', setting_elem)
    setting_elem%value_integer = proc%n_procs_col
    call fson_value_add(setting_in_fson, setting_elem)
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('num_omp_threads', setting_elem)
    setting_elem%value_integer = proc%n_omp_threads
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set group id related settings.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_group_id_used', setting_elem)
    setting_elem%value_logical = setting%is_group_id_used
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%is_group_id_used) then
      setting_elem => fson_value_create()
      call fson_set_name('group_id_filename', setting_elem)
      call fson_set_as_string(trim(setting%group_id_filename), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set settings for precomputed eigenpairs.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('to_use_precomputed_eigenpairs', setting_elem)
    setting_elem%value_logical = setting%to_use_precomputed_eigenpairs
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%to_use_precomputed_eigenpairs) then
      setting_elem => fson_value_create()
      call fson_set_name('eigenvalue_filename', setting_elem)
      call fson_set_as_string(trim(setting%eigenvalue_filename), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
      setting_elem => fson_value_create()
      call fson_set_name('eigenvector_dirname', setting_elem)
      call fson_set_as_string(trim(setting%eigenvector_dirname), setting_elem)
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    ! Set misc flags.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_LOGICAL
    call fson_set_name('is_reduction_mode', setting_elem)
    setting_elem%value_logical = setting%is_reduction_mode
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set re-initialize settings
    setting_elem => fson_value_create()
    call fson_set_name('re_initialize_method', setting_elem)
    call fson_set_as_string(trim(setting%re_initialize_method), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)
    if (setting%re_initialize_method == 'minimize_lcao_error_cutoff') then
      setting_elem => fson_value_create()
      call fson_set_name('vector_cutoff_residual', setting_elem)
      setting_elem%value_type = TYPE_REAL
      setting_elem%value_real = setting%vector_cutoff_residual
      call fson_value_add(setting_in_fson, setting_elem)
    else if (setting%re_initialize_method == 'minimize_lcao_error_suppress' .or. &
         setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress' .or. &
         setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress_orthogonal' .or. &
         setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress_adaptive' .or. &
         setting%re_initialize_method == 'minimize_lcao_error_matrix_suppress_select') then
      setting_elem => fson_value_create()
      call fson_set_name('suppress_constant', setting_elem)
      setting_elem%value_type = TYPE_REAL
      setting_elem%value_real = setting%suppress_constant
      call fson_value_add(setting_in_fson, setting_elem)
    end if

    call fson_value_add(output, setting_in_fson)
  end subroutine add_setting_json


  subroutine add_condition_json(dim, setting, structure, errors, group_id, filter_group_id, &
       eigenvalues, eigenstate_mean, eigenstate_msd, eigenstate_ipratio, output)
    integer, intent(in) :: dim
    type(wp_setting_t), intent(in) :: setting
    type(wp_structure_t), intent(in) :: structure
    type(wp_error_t), intent(in) :: errors
    integer, intent(in) :: group_id(:, :), filter_group_id(:, :)
    real(8), intent(in) :: eigenvalues(setting%num_filter)
    real(8), intent(in) :: eigenstate_mean(3, setting%num_filter), eigenstate_msd(4, setting%num_filter)
    real(8), intent(in) :: eigenstate_ipratio(setting%num_filter)
    type(fson_value), pointer, intent(inout) :: output

    integer :: i, num_groups, num_atoms_in_group, a, g

    type(fson_value), pointer :: condition, condition_elem, group_id_row, element

    condition => fson_value_create()
    condition%value_type = TYPE_OBJECT
    call fson_set_name('condition', condition)

    ! Set dim and num_atoms
    condition_elem => fson_value_create()
    call fson_set_name('dim', condition_elem)
    condition_elem%value_type = TYPE_INTEGER
    condition_elem%value_integer = dim
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_name('num_atoms', condition_elem)
    condition_elem%value_type = TYPE_INTEGER
    condition_elem%value_integer = structure%num_atoms
    call fson_value_add(condition, condition_elem)

    ! Set eigenvalues.
    allocate(condition_elem)
    call fson_set_as_real_array(setting%num_filter, eigenvalues, condition_elem)
    !condition_elem => array_to_json_real(setting%num_filter, eigenvalues)
    call fson_set_name('eigenvalues', condition_elem)
    call fson_value_add(condition, condition_elem)

    ! Set mean coordinates of eigenstates.
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_mean(1, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_mean_x', condition_elem)
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_mean(2, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_mean_y', condition_elem)
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_mean(3, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_mean_z', condition_elem)
    call fson_value_add(condition, condition_elem)

    ! Set MSD of eigenstates.
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_msd(1, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_msd_x', condition_elem)
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_msd(2, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_msd_y', condition_elem)
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_msd(3, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_msd_z', condition_elem)
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_msd(4, 1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_msd_total', condition_elem)
    call fson_value_add(condition, condition_elem)

    ! Set ipratio of eigenstates.
    condition_elem => fson_value_create()
    call fson_set_as_real_array(setting%num_filter, eigenstate_ipratio(1 : setting%num_filter), condition_elem)
    call fson_set_name('eigenstate_ipratio', condition_elem)
    call fson_value_add(condition, condition_elem)

    ! Set filter errors.
    condition_elem => fson_value_create()
    call fson_set_name('absolute_filter_error', condition_elem)
    condition_elem%value_type = TYPE_REAL
    condition_elem%value_real = errors%absolute
    call fson_value_add(condition, condition_elem)
    condition_elem => fson_value_create()
    call fson_set_name('relative_filter_error', condition_elem)
    condition_elem%value_type = TYPE_REAL
    condition_elem%value_real = errors%relative
    call fson_value_add(condition, condition_elem)

    ! Set atom elemets.
    condition_elem => fson_value_create()
    condition_elem%value_type = TYPE_ARRAY
    call fson_set_name('elements', condition_elem)
    do i = 1, structure%num_atoms
      element => fson_value_create()
      call fson_set_as_string(structure%atom_elements(i), element)
      call fson_value_add(condition_elem, element)
    end do
    call fson_value_add(condition, condition_elem)

    ! Set group id.
    if (setting%is_group_id_used) then
      call add_group_id_to_condition_json(group_id, condition)
    end if
    if (trim(setting%filter_mode) == 'group') then
      call add_group_id_to_condition_json(filter_group_id, condition)
    end if

    call fson_value_add(output, condition)
  end subroutine add_condition_json


  subroutine add_group_id_to_condition_json(group_id, condition)
    integer, intent(in) :: group_id(:, :)
    type(fson_value), pointer, intent(inout) :: condition

    integer :: num_groups, num_atoms_in_group, a, g
    type(fson_value), pointer :: condition_item, group_id_row, atom_number

    condition_item => fson_value_create()
    condition_item%value_type = TYPE_ARRAY
    call fson_set_name('group_id', condition_item)
    num_groups = size(group_id, 2)
    do g = 1, num_groups
      group_id_row => fson_value_create()
      group_id_row%value_type = TYPE_ARRAY
      num_atoms_in_group = group_id(1, g)
      do a = 1, num_atoms_in_group
        atom_number => fson_value_create()
        atom_number%value_type = TYPE_INTEGER
        atom_number%value_integer = group_id(a + 1, g)
        call fson_value_add(group_id_row, atom_number)
      end do
      call fson_value_add(condition_item, group_id_row)
    end do
    call fson_value_add(condition, condition_item)
  end subroutine add_group_id_to_condition_json


  subroutine add_state_json(dim, num_filter, structure, i, t, input_step, &
       energies, dv_alpha, dv_psi, dv_charges_on_basis, dv_charges_on_atoms, charge_moment, &
       h1_type, dv_atom_speed, dv_atom_perturb, is_after_matrix_replace, H_sparse, &
       states)
    integer, intent(in) :: dim, num_filter, i, input_step
    type(wp_structure_t), intent(in) :: structure
    type(wp_energy_t), intent(in) :: energies
    type(wp_charge_moment_t), intent(in) :: charge_moment
    double precision, intent(in) :: t
    complex(kind(0d0)), intent(in) :: dv_alpha(num_filter), dv_psi(dim)
    real(8), intent(in) :: dv_charges_on_basis(dim), dv_charges_on_atoms(structure%num_atoms)
    real(8), intent(in) :: dv_atom_speed(dim), dv_atom_perturb(dim)
    character(len=*), intent(in) :: h1_type
    logical, intent(in) :: is_after_matrix_replace
    type(sparse_mat), intent(in) :: H_sparse

    type(fson_value), pointer, intent(inout) :: states

    type(fson_value), pointer :: state, state_elem, state_elem_sub
    real(8) :: norm_alpha, wtime_start, wtime_end
    external :: dznrm2
    double precision :: dznrm2
    real(8) :: diag(dim)

    wtime_start = mpi_wtime()

    state => fson_value_create()
    state%value_type = TYPE_OBJECT

    ! Set time, energies.
    state_elem => fson_value_create()
    call fson_set_name('time', state_elem)
    state_elem%value_type = TYPE_REAL
    state_elem%value_real = t
    call fson_value_add(state, state_elem)
    state_elem => fson_value_create()
    call fson_set_name('step_num', state_elem)
    state_elem%value_type = TYPE_INTEGER
    state_elem%value_integer = i
    call fson_value_add(state, state_elem)
    state_elem => fson_value_create()
    call fson_set_name('TB_energy', state_elem)
    state_elem%value_type = TYPE_REAL
    state_elem%value_real = energies%tightbinding
    call fson_value_add(state, state_elem)
    state_elem => fson_value_create()
    call fson_set_name('NL_energy', state_elem)
    state_elem%value_type = TYPE_REAL
    state_elem%value_real = energies%nonlinear
    call fson_value_add(state, state_elem)
    state_elem => fson_value_create()
    call fson_set_name('total_energy', state_elem)
    state_elem%value_type = TYPE_REAL
    state_elem%value_real = energies%total
    call fson_value_add(state, state_elem)
    state_elem => fson_value_create()
    call fson_set_name('TB_energy_deviation', state_elem)
    state_elem%value_type = TYPE_REAL
    state_elem%value_real = energies%tightbinding_deviation
    call fson_value_add(state, state_elem)

    ! Set alpha (EVE coefficients).
    allocate(state_elem)
    call fson_set_as_cmplx_array(num_filter, dv_alpha, state_elem)
    call fson_set_name('alpha', state_elem)
    ! norm of alpha.
    state_elem_sub => fson_value_create()
    call fson_set_name('norm', state_elem_sub)
    state_elem_sub%value_type = TYPE_REAL
    norm_alpha = dznrm2(num_filter, dv_alpha, 1)
    if (abs(norm_alpha - 1d0) > 1d-5 .and. check_master()) then
      write(0, *) '[Warning] norm of alpha is distant from one: ', norm_alpha, ', in step: ', i
    end if
    state_elem_sub%value_real = norm_alpha
    call fson_value_add(state_elem, state_elem_sub)
    ! ipratio of alpha.
    state_elem_sub => fson_value_create()
    call fson_set_name('ipratio', state_elem_sub)
    state_elem_sub%value_type = TYPE_REAL
    state_elem_sub%value_real = get_ipratio(num_filter, dv_alpha)
    call fson_value_add(state_elem, state_elem_sub)
    call fson_value_add(state, state_elem)

    ! Set psi (LCAO coefficients).
    allocate(state_elem)
    call fson_set_as_cmplx_array(dim, dv_psi, state_elem)
    call fson_set_name('psi', state_elem)
    ! norm of psi.
    state_elem_sub => fson_value_create()
    call fson_set_name('norm', state_elem_sub)
    state_elem_sub%value_type = TYPE_REAL
    state_elem_sub%value_real = dznrm2(dim, dv_psi, 1)
    call fson_value_add(state_elem, state_elem_sub)
    ! ipratio of psi.
    state_elem_sub => fson_value_create()
    call fson_set_name('ipratio', state_elem_sub)
    state_elem_sub%value_type = TYPE_REAL
    state_elem_sub%value_real = get_ipratio(dim, dv_psi)
    call fson_value_add(state_elem, state_elem_sub)
    call fson_value_add(state, state_elem)

    if (trim(h1_type) == 'maxwell') then
      ! Set speed of atoms.
      allocate(state_elem)
      call fson_set_as_real_array(structure%num_atoms, dv_atom_speed(1:structure%num_atoms), state_elem)
      call fson_set_name('atom_speed', state_elem)
      call fson_value_add(state, state_elem)
    else if (trim(h1_type) == 'harmonic') then
      allocate(state_elem)
      call fson_set_as_real_array(structure%num_atoms, dv_atom_perturb(1:structure%num_atoms), state_elem)
      call fson_set_name('atom_perturb', state_elem)
      call fson_value_add(state, state_elem)
    else if (trim(h1_type) == 'harmonic_for_nn_exciton') then
      allocate(state_elem)
      call fson_set_as_real_array(structure%num_atoms / 3 * 2, &
           dv_atom_perturb(1 : structure%num_atoms / 3 * 2), state_elem)
      call fson_set_name('atom_perturb', state_elem)
      call fson_value_add(state, state_elem)
      call sparse_matrix_to_diag(H_sparse, diag)
      allocate(state_elem)
      call fson_set_as_real_array(dim, diag, state_elem)
      call fson_set_name('H_diag', state_elem)
      call fson_value_add(state, state_elem)
    end if

    ! Set Mulliken charges.
    allocate(state_elem)
    !state_elem => array_to_json_cmplx(dim, charges_on_basis)
    call fson_set_as_real_array(dim, dv_charges_on_basis(1:dim), state_elem)
    call fson_set_name('charges_on_basis', state_elem)
    call fson_value_add(state, state_elem)
    allocate(state_elem)
    !state_elem => array_to_json_cmplx(num_atoms, charges_on_atoms)
    call fson_set_as_real_array(structure%num_atoms, dv_charges_on_atoms(1:structure%num_atoms), state_elem)
    call fson_set_name('charges_on_atoms', state_elem)
    call fson_value_add(state, state_elem)

    ! Set statistical information of Mulliken charges.
    if (check_master() .and. mod(i, 100) == 0) then
      write (0, '(A, F16.6, A, E26.16e3)') ' [MSD x at time = ', t, '] ', charge_moment%msds(1)
      write (0, '(A, F16.6, A, E26.16e3)') ' [MSD y at time = ', t, '] ', charge_moment%msds(2)
      write (0, '(A, F16.6, A, E26.16e3)') ' [MSD z at time = ', t, '] ', charge_moment%msds(3)
      write (0, '(A, F16.6, A, E26.16e3)') ' [MSD total at time = ', t, '] ', charge_moment%msds(4)
    end if
    allocate(state_elem)
    !state_elem => array_to_json_real(3, charge_coordinate_mean)
    call fson_set_as_real_array(3, charge_moment%means, state_elem)
    call fson_set_name('charge_coordinate_mean', state_elem)
    call fson_value_add(state, state_elem)
    allocate(state_elem)
    !state_elem => array_to_json_real(4, charge_coordinate_msd)
    call fson_set_as_real_array(4, charge_moment%msds, state_elem)
    call fson_set_name('charge_coordinate_msd', state_elem)
    call fson_value_add(state, state_elem)

    ! Set input step.
    state_elem => fson_value_create()
    call fson_set_name('input_step', state_elem)
    state_elem%value_type = TYPE_INTEGER
    state_elem%value_integer = input_step
    call fson_value_add(state, state_elem)

    ! Set matrix replacement info.
    state_elem => fson_value_create()
    call fson_set_name('is_after_matrix_replace', state_elem)
    state_elem%value_type = TYPE_LOGICAL
    state_elem%value_logical = is_after_matrix_replace
    call fson_value_add(state, state_elem)

    call fson_value_add(states, state)

    wtime_end = mpi_wtime()
    call add_event('add_state_json', wtime_end - wtime_start)
  end subroutine add_state_json


  subroutine add_structure_json(setting, state)
    type(wp_setting_t), intent(in) :: setting
    type(wp_state_t), intent(inout) :: state

    integer :: i, coord, num_filter
    type(fson_value), pointer :: structure_fson, structure_item, element
    character(len=*), parameter :: coordinate_names = 'xyz'

    num_filter = state%Y_filtered_desc(cols_)

    structure_fson => fson_value_create()
    structure_fson%value_type = TYPE_OBJECT
    call fson_set_name('structure', structure_fson)

    ! Save time.
    structure_item => fson_value_create()
    structure_item%value_type = TYPE_REAL
    structure_item%value_real = state%t
    call fson_set_name('time', structure_item)
    call fson_value_add(structure_fson, structure_item)
    ! Save input_step.
    structure_item => fson_value_create()
    structure_item%value_type = TYPE_INTEGER
    structure_item%value_integer = state%input_step
    call fson_set_name('input_step', structure_item)
    call fson_value_add(structure_fson, structure_item)
    ! Save coordinates.
    do coord = 1, 3
      structure_item => fson_value_create()
      structure_item%value_type = TYPE_ARRAY
      call fson_set_as_real_array(state%structure%num_atoms, state%structure%atom_coordinates(coord, :), &
           structure_item)
      call fson_set_name('coordinates_' // coordinate_names(coord : coord), structure_item)
      call fson_value_add(structure_fson, structure_item)
    end do
    ! Save ipratios.
    structure_item => fson_value_create()
    call fson_set_as_real_array(num_filter, state%eigenstate_ipratios, structure_item)
    call fson_set_name('eigenstate_ipratio_on_groups', structure_item)
    call fson_value_add(structure_fson, structure_item)
    ! Save eigenvalues.
    structure_item => fson_value_create()
    call fson_set_as_real_array(num_filter, state%dv_eigenvalues, structure_item)
    call fson_set_name('eigenvalues', structure_item)
    call fson_value_add(structure_fson, structure_item)
    if (setting%to_calculate_eigenstate_moment_every_step) then
      ! Set mean coordinates of eigenstates.
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_mean(1, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_mean_x', structure_item)
      call fson_value_add(structure_fson, structure_item)
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_mean(2, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_mean_y', structure_item)
      call fson_value_add(structure_fson, structure_item)
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_mean(3, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_mean_z', structure_item)
      call fson_value_add(structure_fson, structure_item)
      ! Set MSD of eigenstates.
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_msd(1, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_msd_x', structure_item)
      call fson_value_add(structure_fson, structure_item)
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_msd(2, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_msd_y', structure_item)
      call fson_value_add(structure_fson, structure_item)
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_msd(3, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_msd_z', structure_item)
      call fson_value_add(structure_fson, structure_item)
      structure_item => fson_value_create()
      call fson_set_as_real_array(num_filter, state%eigenstate_msd(4, 1 : num_filter), structure_item)
      call fson_set_name('eigenstate_msd_total', structure_item)
      call fson_value_add(structure_fson, structure_item)
    end if
    ! Save all items.
    call fson_value_add(state%structures, structure_fson)
  end subroutine add_structure_json
end module wp_output_m
