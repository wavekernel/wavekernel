module wk_state_m
  use wk_atom_m
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_fson_m
  use wk_global_variables_m
  use wk_linear_algebra_m
  use wk_matrix_io_m
  use wk_util_m
  implicit none

  type wk_charge_moment_t
    real(8) :: means(3)
    real(8) :: msds(4)
  end type wk_charge_moment_t

  type wk_charge_factor_t
    real(8) :: charge_factor_common, charge_factor_H, charge_factor_C
  end type wk_charge_factor_t

  type wk_fson_values_for_state_t
    type(fson_value), pointer :: output, states, split_files_metadata
  end type wk_fson_values_for_state_t

  type wk_basis_t
    ! Common variables.
    logical :: is_group_filter_mode
    integer :: dim, num_basis
    real(8), allocatable :: dv_eigenvalues(:), dv_ipratios(:), eigenstate_mean(:, :), eigenstate_msd(:, :)
    real(8), allocatable :: eigenstate_ipratios(:)
    integer :: Y_filtered_desc(desc_size)
    real(8), allocatable :: Y_filtered(:, :)  ! m x n
    ! Variables for no filter mode.
    integer :: Y_all_desc(desc_size)
    real(8), allocatable :: Y_all(:, :)  ! m x m
    ! Variables for filter mode.
    real(8), allocatable :: H1_base(:, :), S1(:, :)  ! n x n
    type(wk_distributed_block_diagonal_matrices_t) :: Y_local
    integer, allocatable :: filter_group_id(:, :), filter_group_indices(:, :)
    integer :: H1_desc(desc_size)
  end type wk_basis_t

  type wk_state_t
    type(sparse_mat) :: H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, H1_lcao_sparse_charge_overlap
    type(wk_structure_t) :: structure
    type(wk_energy_t), allocatable :: energies(:)
    type(wk_error_t), allocatable :: errors(:)
    integer :: dim, i, ierr, input_step
    real(8) :: t, t_last_replace
    complex(kind(0d0)), allocatable :: dv_alpha(:, :), dv_psi(:, :), dv_alpha_next(:, :), dv_psi_next(:, :)
    complex(kind(0d0)), allocatable :: dv_alpha_reconcile(:, :), dv_psi_reconcile(:, :)
    real(8), allocatable :: dv_atom_perturb(:), dv_atom_speed(:)
    real(8), allocatable :: dv_charge_on_basis(:, :), dv_charge_on_atoms(:, :)
    type(wk_charge_moment_t), allocatable :: charge_moment(:)
    type(wk_charge_factor_t) :: charge_factor
    real(8) :: wtime_total, wtime
    integer, allocatable :: group_id(:, :)
    type(wk_fson_values_for_state_t), allocatable :: fsons(:)
    type(fson_value), pointer :: structures
    integer :: print_count, total_state_count
    ! MPI realated
    type(wk_basis_t) :: basis
    real(8), allocatable :: YSY_filtered(:, :)  ! n x n
    real(8), allocatable :: H1(:, :)  ! n x n
    complex(kind(0d0)), allocatable :: A(:, :)  ! n x n
    integer :: YSY_filtered_desc(desc_size)
  end type wk_state_t

contains

  subroutine check_psi_S_normal(state, is_psi_S_normal)
    type(wk_state_t), intent(in) :: state
    logical, intent(out) :: is_psi_S_normal
    complex(kind(0d0)) :: dv_work(state%dim)
    complex(kind(0d0)) :: diff
    complex(kind(0d0)) :: zdotc  ! Function.
    dv_work(:) = kZero
    call matvec_sd_z('No', state%S_sparse, kOne, state%dv_psi, kZero, dv_work)
    diff = zdotc(state%dim, state%dv_psi, 1, dv_work, 1) - kOne
    is_psi_S_normal = abs(diff) < 1d-10
    if (.not. is_psi_S_normal) then
      print *, '[Warn] \psi^\dagger S \psi - 1: ', diff
    end if
  end subroutine check_psi_S_normal
end module wk_state_m
