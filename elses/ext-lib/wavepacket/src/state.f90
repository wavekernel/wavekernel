module wp_state_m
  use wp_atom_m
  use wp_charge_m
  use wp_distribute_matrix_m
  use wp_fson_m
  use wp_linear_algebra_m
  use wp_matrix_io_m
  use wp_util_m
  implicit none

  type wp_fson_values_for_state_t
    type(fson_value), pointer :: output, states, split_files_metadata
  end type wp_fson_values_for_state_t

  type wp_state_t
    type(sparse_mat) :: H_sparse, S_sparse, H_sparse_prev, S_sparse_prev, H_multistep_sparse, S_multistep_sparse
    type(wp_structure_t) :: structure
    type(wp_energy_t), allocatable :: energies(:)
    type(wp_error_t), allocatable :: errors(:)
    integer :: dim, i, ierr, input_step
    real(8) :: t, t_last_replace
    complex(kind(0d0)), allocatable :: dv_alpha(:, :), dv_psi(:, :), dv_alpha_next(:, :), dv_psi_next(:, :)
    complex(kind(0d0)), allocatable :: dv_alpha_reconcile(:, :), dv_psi_reconcile(:, :)
    real(8), allocatable :: dv_atom_perturb(:), dv_atom_speed(:)
    real(8), allocatable :: dv_charge_on_basis(:, :), dv_charge_on_atoms(:, :)
    real(8), allocatable :: dv_eigenvalues(:), dv_ipratios(:), eigenstate_mean(:, :), eigenstate_msd(:, :)
    real(8), allocatable :: eigenstate_ipratios(:)
    type(wp_charge_moment_t), allocatable :: charge_moment(:)
    type(wp_charge_factor_t) :: charge_factor
    real(8) :: wtime_total, wtime
    integer, allocatable :: group_id(:, :), filter_group_id(:, :), filter_group_indices(:, :)
    type(wp_fson_values_for_state_t), allocatable :: fsons(:)
    type(fson_value), pointer :: structures
    integer :: print_count, total_state_count
    type(wp_local_matrix_t), allocatable :: Y_local(:)

    ! MPI realated
    real(8), allocatable :: Y(:, :)  ! m x m
    real(8), allocatable :: Y_filtered(:, :)  ! m x n
    real(8), allocatable :: YSY_filtered(:, :)  ! n x n
    real(8), allocatable :: H1(:, :), H1_base(:, :)  ! n x n
    complex(kind(0d0)), allocatable :: A(:, :)  ! n x n
    real(8), allocatable :: H1_multistep(:, :)
    integer :: Y_desc(desc_size), Y_filtered_desc(desc_size), YSY_filtered_desc(desc_size)
    integer :: H1_desc(desc_size), A_desc(desc_size)
  end type wp_state_t

contains

  subroutine check_psi_S_normal(state, is_psi_S_normal)
    type(wp_state_t), intent(in) :: state
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
end module wp_state_m
