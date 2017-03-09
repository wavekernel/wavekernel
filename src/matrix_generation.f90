module wk_matrix_generation_m
  use mpi
  use wk_atom_m
  use wk_charge_m
  use wk_descriptor_parameters_m
  use wk_distribute_matrix_m
  use wk_event_logger_m
  use wk_maxwell_distribution_m
  use wk_processes_m
  use wk_conversion_m
  use wk_global_variables_m
  use wk_linear_algebra_m
  use wk_state_m
  use wk_util_m
  implicit none

  private
  public :: make_H1, make_A, aux_make_H1_charge_with_overlap, aux_make_H1_harmonic_for_nn_exciton

contains

  !! LACO の世界で対角の H1 を作り, Y^\dagger H1 Y とすることで
  !! alpha の世界における対応する H1 を作る. これは対角にならないことに注意.
  !! eve stands for eigenvalue expansion.
  !! Complexity: O(m n).
  !subroutine make_H1_diag(proc, Y, Y_desc, is_group_filter_mode, filter_group_indices, Y_local, &
  !     H1_alpha, H1_alpha_desc)
  !  type(wk_process_t), intent(in) :: proc
  !  real(8), intent(in) :: Y(:, :)
  !  type(wk_local_matrix_t), intent(in) :: Y_local(:)
  !  integer, intent(in) :: Y_desc(desc_size), H1_alpha_desc(desc_size), filter_group_indices(:, :)
  !  logical, intent(in) :: is_group_filter_mode
  !  real(8), intent(out) :: H1_alpha(:, :)
  !
  !  integer :: i, dim
  !  real(8) :: H1_lcao_diag(Y_desc(rows_))
  !
  !  dim = Y_desc(rows_)
  !  do i = 1, dim
  !    H1_lcao_diag(i) = 0.05d0
  !  end do
  !
  !  if (is_group_filter_mode) then
  !    call terminate('not implemented', 68)
  !    !call change_basis_lcao_diag_to_alpha_group_filter_real(proc, filter_group_indices, Y_local, &
  !    !     H1_lcao_diag, H1_alpha, H1_alpha_desc)
  !  else
  !    call change_basis_lcao_diag_to_alpha(proc, Y, Y_desc, H1_lcao_diag, H1_alpha, H1_alpha_desc)
  !  end if
  !end subroutine make_H1_diag


  ! make_H1_diag のゼロ行列版.
  ! 変換してもゼロなのは分かっているので基底変換の計算は不要.
  ! Complexity: O(n^2).
  subroutine make_H1_zero(H1_alpha)
    real(8), intent(out) :: H1_alpha(:, :)

    H1_alpha(:, :) = 0.0d0
  end subroutine make_H1_zero


  !! Complexity: O(m n^2).
  !subroutine make_H1_charge_on_atoms(proc, structure, Y, Y_desc, &
  !     is_group_filter_mode, filter_group_indices, Y_local, &
  !     charge_on_atoms, charge_factor, &
  !     H1_alpha, H1_alpha_desc)
  !  type(wk_process_t), intent(in) :: proc
  !  type(wk_structure_t), intent(in) :: structure
  !  integer, intent(in) :: filter_group_indices(:, :)
  !  integer, intent(in) :: Y_desc(desc_size), H1_alpha_desc(desc_size)
  !  real(8), intent(in) :: Y(:, :)
  !  type(wk_local_matrix_t), intent(in) :: Y_local(:)
  !  real(8) :: charge_on_atoms(structure%num_atoms)
  !  logical, intent(in) :: is_group_filter_mode
  !  type(wk_charge_factor_t), intent(in) :: charge_factor
  !  real(8), intent(out) :: H1_alpha(:, :)
  !
  !  integer :: i, j, dim, num_filter
  !  !integer :: atom_valence
  !  real(8) :: diag, H1_lcao_diag(Y_desc(rows_)), wtime_start, wtime_end
  !
  !  wtime_start = mpi_wtime()
  !
  !  dim = Y_desc(rows_)
  !  num_filter = Y_desc(cols_)
  !
  !  do i = 1, structure%num_atoms
  !    diag = - get_charge_factor(i, structure, charge_factor) * charge_on_atoms(i)
  !    ! TODO: Should atom_valence be considered? The lines comment out above are for this consideration.
  !    !diag = -charge_factor * abs(charge_on_atom_i) / dble(atom_valence)
  !
  !    do j = structure%atom_indices(i), structure%atom_indices(i + 1) - 1
  !      H1_lcao_diag(j) = diag
  !    end do
  !  end do
  !
  !  wtime_end = mpi_wtime()
  !  call add_event('make_H1_charge_on_atoms:set_diag', wtime_end - wtime_start)
  !  wtime_start = wtime_end
  !
  !  if (is_group_filter_mode) then
  !    call terminate('not implemented', 49)
  !    !call change_basis_lcao_diag_to_alpha_group_filter_real(proc, filter_group_indices, Y_local, &
  !    !     H1_lcao_diag, H1_alpha, H1_alpha_desc)
  !  else
  !    call change_basis_lcao_diag_to_alpha(proc, Y, Y_desc, H1_lcao_diag, H1_alpha, H1_alpha_desc)
  !  end if
  !
  !  wtime_end = mpi_wtime()
  !  call add_event('make_H1_charge_on_atoms:change_basis', wtime_end - wtime_start)
  !end subroutine make_H1_charge_on_atoms


  subroutine aux_make_H1_charge_with_overlap(structure, S_sparse, &
       charge_on_basis, charge_on_atoms, charge_factor, H1_lcao_sparse)
    type(wk_structure_t), intent(in) :: structure
    type(sparse_mat), intent(in) :: S_sparse
    real(8), intent(in) :: charge_on_basis(:), charge_on_atoms(structure%num_atoms)
    type(wk_charge_factor_t), intent(in) :: charge_factor
    type(sparse_mat), intent(inout) :: H1_lcao_sparse

    integer :: i, j, k, atom_i, atom_j, ierr
    real(8) :: overlap_factor, overlap_factor_i, overlap_factor_j, wtime_start, wtime_end
    integer :: S_desc(desc_size)

    if (allocated(H1_lcao_sparse%value)) then
      deallocate(H1_lcao_sparse%value)
    end if
    if (allocated(H1_lcao_sparse%suffix)) then
      deallocate(H1_lcao_sparse%suffix)
    end if
    H1_lcao_sparse%size = S_sparse%size
    H1_lcao_sparse%num_non_zeros = S_sparse%num_non_zeros
    allocate(H1_lcao_sparse%value(S_sparse%num_non_zeros), H1_lcao_sparse%suffix(2, S_sparse%num_non_zeros))
    H1_lcao_sparse%value(:) = S_sparse%value(:)
    H1_lcao_sparse%suffix(:, :) = S_sparse%suffix(:, :)

    !wtime_end = mpi_wtime()
    !call add_event('make_H1_charge_with_overlap:copy_overlap', wtime_end - wtime_start)
    !wtime_start = wtime_end

    do k = 1, H1_lcao_sparse%num_non_zeros
      i = H1_lcao_sparse%suffix(1, k)
      j = H1_lcao_sparse%suffix(2, k)
      call lcao_index_to_atom_index(structure, i, atom_i)
      call lcao_index_to_atom_index(structure, j, atom_j)
      ! For diagonal terms (i == j), sum two -0.5d0 -> -1d0
      if (.false.) then
        overlap_factor_i = - 0.5d0 * charge_on_basis(i) * get_charge_factor(atom_i, structure, charge_factor)
        overlap_factor_j = - 0.5d0 * charge_on_basis(j) * get_charge_factor(atom_j, structure, charge_factor)
      else
        overlap_factor_i = - 0.5d0 * charge_on_atoms(atom_i) * get_charge_factor(atom_i, structure, charge_factor)
        overlap_factor_j = - 0.5d0 * charge_on_atoms(atom_j) * get_charge_factor(atom_j, structure, charge_factor)
      end if
      overlap_factor = overlap_factor_i + overlap_factor_j
      H1_lcao_sparse%value(k) = overlap_factor * H1_lcao_sparse%value(k)
    end do
  end subroutine aux_make_H1_charge_with_overlap


  subroutine make_H1_charge_with_overlap(structure, S_sparse, basis, &
       charge_on_basis, charge_on_atoms, charge_factor, &
       H1_alpha, H1_alpha_desc)
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(in) :: basis
    type(sparse_mat), intent(in) :: S_sparse
    integer, intent(in) :: H1_alpha_desc(desc_size)
    real(8), intent(in) :: charge_on_basis(:), charge_on_atoms(structure%num_atoms)
    type(wk_charge_factor_t), intent(in) :: charge_factor
    real(8), intent(out) :: H1_alpha(:, :)

    integer :: i, j, k, dim, atom_i, atom_j, i1, i2, j1, j2
    integer :: nprow, npcol, myrow, mycol, i_local, j_local, ierr
    real(8) :: overlap_factor, overlap_factor_i, overlap_factor_j, wtime_start, wtime_end
    type(sparse_mat) :: H1_lcao_sparse

    wtime_start = mpi_wtime()

    call aux_make_H1_charge_with_overlap(structure, S_sparse, &
         charge_on_basis, charge_on_atoms, charge_factor, H1_lcao_sparse)

    wtime_end = mpi_wtime()
    call add_event('make_H1_charge_with_overlap:multiply_overlap_by_charge', wtime_end - wtime_start)
    wtime_start = wtime_end

    call change_basis_lcao_to_alpha(basis, H1_lcao_sparse, H1_alpha, H1_alpha_desc)

    wtime_end = mpi_wtime()
    call add_event('make_H1_charge_with_overlap:change_basis', wtime_end - wtime_start)
  end subroutine make_H1_charge_with_overlap


  !subroutine init_speed_with_maxwell(proc, num_atoms, atom_elements, &
  !     temperature, full_vecs, full_vecs_desc)
  !  type(wk_process_t), intent(in) :: proc
  !  integer, intent(in) :: num_atoms
  !  character, intent(in) :: atom_elements(num_atoms)
  !  complex(kind(0d0)), intent(out) :: full_vecs(:, :)
  !  integer, intent(in) :: full_vecs_desc(desc_size)
  !  real(8), intent(in) :: temperature
  !
  !  integer :: atom, nprow, npcol, myprow, mypcol, prow_own, pcol_own, i, j
  !  real(8) :: mass, scale
  !  integer :: indxg2p, indxg2l
  !
  !  call blacs_gridinfo(proc%context, nprow, npcol, myprow, mypcol)
  !  pcol_own = indxg2p(col_atom_speed, full_vecs_desc(block_col_), 0, full_vecs_desc(csrc_), npcol)
  !  do atom = 1, num_atoms
  !    prow_own = indxg2p(atom, full_vecs_desc(block_row_), 0, full_vecs_desc(rsrc_), nprow)
  !    if (myprow == prow_own .and. mypcol == pcol_own) then
  !      if (atom_elements(atom) == 'H') then
  !        mass = kMassHydrogen
  !      else if (atom_elements(atom) == 'C') then
  !        mass = kMassCarbon
  !      else
  !        call terminate('make_H1_maxwell: unknown atom element', 1)
  !      end if
  !      scale = dsqrt(kBoltzmann * temperature / mass)
  !      i = indxg2l(atom, full_vecs_desc(block_row_), 0, 0, nprow)
  !      j = indxg2l(col_atom_speed, full_vecs_desc(block_col_), 0, 0, npcol)
  !      full_vecs(i, j) = cmplx(maxwell_rvs(scale), 0d0, kind(0d0))
  !    end if
  !  end do
  !end subroutine init_speed_with_maxwell
  !
  !
  !subroutine make_H1_maxwell(proc, num_atoms, atom_indices, atom_elements, &
  !     Y, Y_desc, is_group_filter_mode, filter_group_indices, Y_local, &
  !     is_init, is_restart_mode, temperature, delta_time, &
  !     full_vecs, full_vecs_desc, H1_alpha, H1_alpha_desc)
  !  type(wk_process_t), intent(in) :: proc
  !  integer, intent(in) :: num_atoms, atom_indices(num_atoms + 1), filter_group_indices(:, :)
  !  character, intent(in) :: atom_elements(num_atoms)
  !  complex(kind(0d0)), intent(inout) :: full_vecs(:, :)
  !  complex(kind(0d0)), intent(in) :: Y(:, :)
  !  type(wk_local_matrix_t), intent(in) :: Y_local(:)
  !  integer, intent(in) :: Y_desc(desc_size), full_vecs_desc(desc_size), H1_alpha_desc(desc_size)
  !  logical, intent(in) :: is_init, is_restart_mode, is_group_filter_mode
  !  real(8), intent(in) :: temperature, delta_time
  !  complex(kind(0d0)), intent(out) :: H1_alpha(:, :)
  !
  !  integer :: dim, atom, nprow, npcol, myprow, mypcol, prow_own, pcol_own, i, j
  !  real(8) :: mass, rand, cumulation, scale, speed, energy, wtime_start, wtime_end
  !  complex(kind(0d0)), allocatable :: H1_lcao_diag(:)
  !  logical :: accelerate
  !  integer :: indxg2p, indxg2l
  !
  !  wtime_start = mpi_wtime()
  !  call blacs_gridinfo(proc%context, nprow, npcol, myprow, mypcol)
  !
  !  if (is_init .and. .not. is_restart_mode) then
  !    call init_speed_with_maxwell(proc, num_atoms, atom_elements, &
  !         temperature, full_vecs, full_vecs_desc)
  !    wtime_end = mpi_wtime()
  !    call add_event('make_H1_maxwell:init_speed_with_maxwell', wtime_end - wtime_start)
  !    wtime_start = wtime_end
  !  end if
  !
  !  dim = Y_desc(rows_)
  !  allocate(H1_lcao_diag(dim))
  !  H1_lcao_diag(:) = kZero
  !  pcol_own = indxg2p(col_atom_speed, full_vecs_desc(block_col_), 0, full_vecs_desc(csrc_), npcol)
  !  do atom = 1, num_atoms
  !    prow_own = indxg2p(atom, full_vecs_desc(block_row_), 0, full_vecs_desc(rsrc_), nprow)
  !    if (myprow == prow_own .and. mypcol == pcol_own) then
  !      if (atom_elements(atom) == 'H') then
  !        mass = kMassHydrogen
  !      else if (atom_elements(atom) == 'C') then
  !        mass = kMassCarbon
  !      else
  !        call terminate('make_H1_maxwell: unknown atom element', 1)
  !      end if
  !      scale = dsqrt(kBoltzmann * temperature / mass)
  !      i = indxg2l(atom, full_vecs_desc(block_row_), 0, 0, nprow)
  !      j = indxg2l(col_atom_speed, full_vecs_desc(block_col_), 0, 0, npcol)
  !      speed = real(full_vecs(i, j), kind(0d0))
  !      if (.not. is_init) then
  !        cumulation = maxwell_cdf(scale, speed)
  !        call random_number(rand)
  !        accelerate = rand >= cumulation
  !        if (accelerate) then
  !          speed = speed * (1.0 + delta_time * kAccelRatio)
  !        else
  !          speed = max(0d0, speed * (1.0 - delta_time * kAccelRatio))
  !        end if
  !        full_vecs(i, j) = speed
  !      end if
  !      energy = 0.5d0 * mass * (speed ** 2d0)
  !      H1_lcao_diag(atom_indices(atom) : atom_indices(atom + 1) - 1) = &
  !           cmplx(energy, 0d0, kind(0d0))
  !    end if
  !  end do
  !
  !  wtime_end = mpi_wtime()
  !  call add_event('make_H1_maxwell:set_diag', wtime_end - wtime_start)
  !  wtime_start = wtime_end
  !
  !  call zgsum2d(proc%context, 'Row', ' ', dim, 1, H1_lcao_diag, dim, 0, pcol_own)
  !  if (myprow == 0 .and. mypcol == pcol_own) then
  !    call zgebs2d(proc%context, 'All', ' ', dim, 1, H1_lcao_diag, dim)
  !  else
  !    call zgebr2d(proc%context, 'All', ' ', dim, 1, H1_lcao_diag, dim, 0, pcol_own)
  !  end if
  !
  !  wtime_end = mpi_wtime()
  !  call add_event('make_H1_maxwell:sum_and_broadcast', wtime_end - wtime_start)
  !  wtime_start = wtime_end
  !
  !  if (is_group_filter_mode) then
  !    call change_basis_lcao_diag_to_alpha_group_filter_real(proc, filter_group_indices, Y_local, &
  !         H1_lcao_diag, H1_alpha, H1_alpha_desc)
  !  else
  !    call change_basis_lcao_diag_to_alpha_real(proc, Y, Y_desc, H1_lcao_diag, H1_alpha, H1_alpha_desc)
  !  end if
  !
  !  wtime_end = mpi_wtime()
  !  call add_event('make_H1_maxwell:change_basis', wtime_end - wtime_start)
  !  wtime_start = wtime_end
  !end subroutine make_H1_maxwell
  !

  subroutine make_H1_harmonic(structure, basis, &
       is_init, is_restart_mode, &
       t, temperature, delta_time, perturb_interval, &
       dv_atom_perturb, &
       H1_alpha, H1_alpha_desc)
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(in) :: basis
    logical, intent(in) :: is_init, is_restart_mode
    integer, intent(in) :: H1_alpha_desc(desc_size)
    real(8), intent(in) :: t, temperature, delta_time, perturb_interval
    real(8), intent(inout) :: dv_atom_perturb(:)
    real(8), intent(out) :: H1_alpha(:, :)

    integer :: atom, i, j, ierr
    integer :: seed = 10  ! Static variable.
    real(8) :: phase, energy, wtime_start, wtime_end, random_val(structure%num_atoms), val
    real(8) :: H1_lcao_diag(basis%dim)
    logical :: to_refresh

    wtime_start = mpi_wtime()
    ! Temporary implementation. Perturbation term is not inherited in restart_mode now.
    to_refresh = (is_init .and. .not. is_restart_mode) .or. &
         t < floor((t + delta_time) / perturb_interval) * perturb_interval - 1d-10

    H1_lcao_diag(:) = 0d0

    if (check_master()) then
      do atom = 1, structure%num_atoms
        call wk_random_number(seed, val)  ! 0 <= phase < 1
        random_val(atom) = val
      end do
    end if
    call mpi_bcast(random_val, structure%num_atoms, mpi_double_precision, g_wk_master_pnum, mpi_comm_world, ierr)

    do atom = 1, structure%num_atoms
      if (to_refresh) then
        energy = temperature * (sin(random_val(atom) * 2d0 * kPi) ** 2d0 - 0.5)
        dv_atom_perturb(atom) = energy
      end if
      H1_lcao_diag(structure%atom_indices(atom) : structure%atom_indices(atom + 1) - 1) = dv_atom_perturb(atom)
    end do

    wtime_end = mpi_wtime()
    call add_event('make_H1_harmonic:set_diag', wtime_end - wtime_start)
    wtime_start = wtime_end

    call change_basis_lcao_diag_to_alpha(basis, H1_lcao_diag, H1_alpha, H1_alpha_desc)

    wtime_end = mpi_wtime()
    call add_event('make_H1_harmonic:change_basis', wtime_end - wtime_start)
  end subroutine make_H1_harmonic


  subroutine aux_make_H1_harmonic_for_nn_exciton(dv_atom_perturb, temperature, to_refresh, &
       t, t_last_refresh, num_sites, H1_lcao_diag)
    integer, intent(in) :: num_sites
    real(8), intent(in) :: t, temperature
    real(8), intent(inout) :: dv_atom_perturb(:), t_last_refresh
    real(8), intent(out) :: H1_lcao_diag(:)
    logical, intent(in) :: to_refresh

    integer :: site, e_or_h, ierr
    integer :: seed = 10  ! Static variable.
    real(8) :: energy, random_val(num_sites * 2), val

    H1_lcao_diag(:) = 0d0

    if (to_refresh) then
      t_last_refresh = t
      if (check_master()) then
        do site = 1, num_sites
          do e_or_h = 1, 2
            call wk_random_number(seed, val)  ! 0 <= phase < 1
            random_val((site - 1) * 2 + e_or_h) = val
          end do
        end do
      end if
      call mpi_bcast(random_val, num_sites * 2, mpi_double_precision, &
           g_wk_master_pnum, mpi_comm_world, ierr)

      do site = 1, num_sites
        do e_or_h = 1, 2
          energy = temperature * (sin(random_val((site - 1) * 2 + e_or_h) * 2d0 * kPi) ** 2d0 - 0.5)
          dv_atom_perturb((site - 1) * 2 + e_or_h) = energy
        end do
      end do
    end if

    do site = 1, num_sites
      H1_lcao_diag((site - 1) * 3 + 1) = dv_atom_perturb((site - 1) * 2 + 1) + dv_atom_perturb((site - 1) * 2 + 2)
      if (site < num_sites) then
        H1_lcao_diag((site - 1) * 3 + 2) = dv_atom_perturb((site - 1) * 2 + 1) + dv_atom_perturb(site * 2 + 2)
        H1_lcao_diag((site - 1) * 3 + 3) = dv_atom_perturb(site * 2 + 1) + dv_atom_perturb((site - 1) * 2 + 2)
      end if
    end do
  end subroutine aux_make_H1_harmonic_for_nn_exciton


  subroutine make_H1_harmonic_for_nn_exciton(structure, basis, &
       is_init, is_restart_mode, &
       t, temperature, delta_time, perturb_interval, &
       dv_atom_perturb, &
       H1_alpha, H1_alpha_desc)
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(in) :: basis
    logical, intent(in) :: is_init, is_restart_mode
    integer, intent(in) :: H1_alpha_desc(desc_size)
    real(8), intent(in) :: t, temperature, delta_time, perturb_interval
    real(8), intent(inout) :: dv_atom_perturb(:)
    real(8), intent(out) :: H1_alpha(:, :)

    real(8) :: wtime_start, wtime_end
    real(8) :: H1_lcao_diag(basis%dim)
    logical :: to_refresh
    real(8), save :: t_last_refresh = -1d100

    wtime_start = mpi_wtime()
    ! Temporary implementation. Perturbation term is not inherited in restart_mode now.
    to_refresh = (is_init .and. .not. is_restart_mode) .or. &
         (t - t_last_refresh >= perturb_interval)

    call aux_make_H1_harmonic_for_nn_exciton(dv_atom_perturb, temperature, to_refresh, &
         t, t_last_refresh, structure%num_atoms / 3, H1_lcao_diag)

    wtime_end = mpi_wtime()
    call add_event('make_H1_harmonic_for_nn_exciton:set_diag', wtime_end - wtime_start)
    wtime_start = wtime_end

    call change_basis_lcao_diag_to_alpha(basis, H1_lcao_diag, H1_alpha, H1_alpha_desc)

    wtime_end = mpi_wtime()
    call add_event('make_H1_harmonic_for_nn_exciton:change_basis', wtime_end - wtime_start)
  end subroutine make_H1_harmonic_for_nn_exciton


  subroutine make_H1(h1_type, structure, &
       S_sparse, basis, is_init, is_restart_mode, &
       t, temperature, delta_time, perturb_interval, &
       charge_on_basis, charge_on_atoms, charge_factor, &
       dv_atom_perturb, &
       H1_alpha, H1_alpha_desc)
    character(*), intent(in) :: h1_type
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(in) :: basis
    type(sparse_mat), intent(in) :: S_sparse
    integer, intent(in) :: H1_alpha_desc(desc_size)
    real(8), intent(in) :: t, temperature, delta_time, perturb_interval
    type(wk_charge_factor_t), intent(in) :: charge_factor
    real(8), intent(in) :: charge_on_basis(:), charge_on_atoms(structure%num_atoms)
    logical, intent(in) :: is_init, is_restart_mode
    real(8), intent(inout) :: dv_atom_perturb(:)
    real(8), intent(out) :: H1_alpha(:, :)

    !if (trim(h1_type) == 'diag') then
    !  call make_H1_diag(Y, Y_desc, is_group_filter_mode, filter_group_indices, Y_local, &
    !       H1_alpha, H1_alpha_desc)
    if (trim(h1_type) == 'zero' .or. trim(h1_type) == 'zero_damp' .or. trim(h1_type) == 'zero_sparse' .or. &
         trim(h1_type) == 'zero_damp_charge_base' .or. trim(h1_type) == 'zero_damp_charge_atom') then
      call make_H1_zero(H1_alpha)
    !else if (trim(h1_type) == 'charge') then
    !  call make_H1_charge_on_atoms(structure, Y, Y_desc, &
    !       is_group_filter_mode, filter_group_indices, Y_local, &
    !       charge_on_atoms, charge_factor, &
    !       H1_alpha, H1_alpha_desc)
    else if (trim(h1_type) == 'charge_overlap') then
      call make_H1_charge_with_overlap(structure, S_sparse, basis, &
           charge_on_basis, charge_on_atoms, charge_factor, &
           H1_alpha, H1_alpha_desc)
    !else if (trim(h1_type) == 'maxwell') then
    !  call make_H1_maxwell(structure, Y, Y_desc, &
    !       is_group_filter_mode, filter_group_indices, Y_local, &
    !       is_init, is_restart_mode, temperature, delta_time, &
    !       H1_alpha, H1_alpha_desc)
    else if (trim(h1_type) == 'harmonic') then
      call make_H1_harmonic(structure, basis, &
           is_init, is_restart_mode, &
           t, temperature, delta_time, perturb_interval, &
           dv_atom_perturb, &
           H1_alpha, H1_alpha_desc)
    else if (trim(h1_type) == 'harmonic_for_nn_exciton') then
      call make_H1_harmonic_for_nn_exciton(structure, basis, &
           is_init, is_restart_mode, &
           t, temperature, delta_time, perturb_interval, &
           dv_atom_perturb, &
           H1_alpha, H1_alpha_desc)
    else
      stop 'unknown H1 function type (fail in command line option parsing)'
    end if
  end subroutine make_H1


  ! ノート (21) 式.
  ! (A')_{kj} = (H1)_{kj} * \exp{i (\lambda_k - \lambda_j) t}, A <- A'
  ! Complexity: O(n^2).
  ! amplitude_print_threshold, amplitude_print_interval, delta_time and fst_filter do not affect simulation.
  subroutine make_A(H1_eve, H1_eve_desc, eigenvalues, t, &
       amplitude_print_threshold, amplitude_print_interval, delta_time, fst_filter, A, A_desc)
    real(8), intent(in) :: H1_eve(:, :)
    complex(kind(0d0)), intent(out) :: A(:, :)
    integer, intent(in) :: H1_eve_desc(desc_size), A_desc(desc_size), fst_filter
    real(8), intent(in) :: eigenvalues(H1_eve_desc(rows_)), t
    real(8), intent(in) :: amplitude_print_threshold, amplitude_print_interval, delta_time

    integer :: i, j, num_filter
    complex(kind(0d0)) :: elem
    integer :: rsrc, csrc, myrow, mycol, nprow, npcol, myrank
    integer :: lrindx, lcindx, ierr
    integer :: max_amp_row, max_amp_col
    real(8) :: amp, max_amp

    num_filter = H1_eve_desc(rows_)
    A(:, :) = dcmplx(H1_eve(:, :), 0d0)
    do i = 1, num_filter
      call pzelget('Self', ' ', elem, A, i, i, A_desc)
      call pzelset(A, i, i, A_desc, elem + dcmplx(eigenvalues(i), 0d0))
    end do

    ! Print information of large amplitude elements.
    ! Skip if amplitude_print_threshold is negative.
    if (amplitude_print_threshold >= 0d0) then
      if (t == 0d0 .or. &
           t < floor((t + delta_time) / amplitude_print_interval) * amplitude_print_interval - 1d-10) then
        call print_offdiag_norm('A', A, A_desc)
        call mpi_barrier(mpi_comm_world, ierr)
        if (check_master()) then
          write (0, '(A, F16.6, A, F16.6)') ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
               '] make_A() : Print large amplitude elements in t = ', t
        end if
        max_amp_row = 0
        max_amp_col = 0
        max_amp = -1d0
        do j = 1, num_filter
          do i = 1, j  ! Only upper triangle.
            call infog2l(i, j, A_desc, &
                 nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc)
            if (rsrc == myrow .and. csrc == mycol) then
              amp = abs(A(lrindx, lcindx))
              if (amp >= amplitude_print_threshold) then
                print *, '[Large amplitude element in t = ', t, '] A(', &
                     i + fst_filter - 1, ', ', j + fst_filter - 1, ') ', A(lrindx, lcindx)
              end if
              if (amp > max_amp) then
                max_amp = amp
                max_amp_row = i
                max_amp_col = j
              end if
            end if
          end do
        end do
        if (max_amp_row /= 0 .and. max_amp_col /= 0) then
          print *, '[Max amplitude in t = ', t, 'on processes] proc ', myrank, ':  A(', &
               max_amp_row + fst_filter - 1, ', ', max_amp_col + fst_filter - 1, ') ', max_amp
        end if
      end if
    end if
  end subroutine make_A
end module wk_matrix_generation_m
