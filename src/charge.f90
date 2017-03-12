module wk_charge_m
  use mpi
  use wk_atom_m
  use wk_descriptor_parameters_m
  use wk_event_logger_m
  use wk_linear_algebra_m
  use wk_matrix_io_m
  use wk_global_variables_m
  use wk_state_m
  use wk_util_m
  implicit none

  private
  external :: blacs_pnum
  integer :: blacs_pnum

  public :: get_mulliken_charges_on_basis, get_mulliken_charges_on_atoms, &
       get_eigenstate_charges_on_groups, get_charge_overlap_energy, &
       get_mulliken_charge_coordinate_moments, get_msd_of_eigenstates, get_ipratio_of_eigenstates, get_charge_factor

contains

  ! Complexity: O(m^2).
  ! full_work <- S psi
  ! charge(col_charges_on_basis)_k = real(conjg(psi_k) * full_work_k)
  subroutine get_mulliken_charges_on_basis(dim, S_sparse, dv_psi, dv_charge_on_basis)
    integer, intent(in) :: dim
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(dim)
    real(8), intent(out) :: dv_charge_on_basis(dim)

    integer :: k
    complex(kind(0d0)) :: dv_work(dim)
    real(8) :: charge_k, sum_mulliken_charges
    real(8) :: wtime_start, wtime_end

    wtime_start = mpi_wtime()

    dv_work(:) = kZero
    call matvec_sd_z('No', S_sparse, kOne, dv_psi, kZero, dv_work)

    wtime_end = mpi_wtime()
    call add_event('get_mulliken_charges_on_basis:matvec_sd_z', wtime_end - wtime_start)
    wtime_start = wtime_end

    sum_mulliken_charges = 0d0  ! Must be 1 after summation. Valid only in master.
    do k = 1, dim
      charge_k = dble(dconjg(dv_psi(k)) * dv_work(k))
      dv_charge_on_basis(k) = charge_k
      sum_mulliken_charges = sum_mulliken_charges + charge_k
    end do

    wtime_end = mpi_wtime()
    call add_event('get_mulliken_charges_on_basis:sum_charges', wtime_end - wtime_start)

    if (check_master()) then
      if (abs(sum_mulliken_charges - 1d0) > charge_sum_error_tol) then
        write(0, *) 'Warn: sum of mulliken charges on LCAO basis is distant from one'
        write(0, '(A, E26.16e3)') ' sum_mulliken_charges: ', sum_mulliken_charges
      end if
    end if
  end subroutine get_mulliken_charges_on_basis


  ! Complexity: O(m^2).
  subroutine get_mulliken_charges_on_atoms(dim, structure, S_sparse, dv_psi, dv_charge_on_atoms)
    integer, intent(in) :: dim
    type(wk_structure_t) :: structure
    type(sparse_mat), intent(in) :: S_sparse
    complex(kind(0d0)), intent(in) :: dv_psi(dim)
    real(8), intent(out) :: dv_charge_on_atoms(structure%num_atoms)

    integer :: a, atom_index_start, atom_index_end, k
    complex(kind(0d0)) :: dv_work(dim)
    real(8) :: wtime_start, wtime_end, sum_mulliken_charges, charge_a

    wtime_start = mpi_wtime()

    dv_work(:) = kZero
    call matvec_sd_z('No', S_sparse, kOne, dv_psi, kZero, dv_work)

    wtime_end = mpi_wtime()
    call add_event('get_mulliken_charges_on_atoms:matvec_sd_z', wtime_end - wtime_start)
    wtime_start = wtime_end

    sum_mulliken_charges = 0d0  ! Must be 1 after summation.
    do a = 1, structure%num_atoms
      charge_a = 0d0
      atom_index_start = structure%atom_indices(a)
      atom_index_end = structure%atom_indices(a + 1) - 1
      do k = atom_index_start, atom_index_end
        charge_a = charge_a + dble(dconjg(dv_psi(k)) * dv_work(k))
        call check_nan_scalar('get_mulliken_charges_on_atoms', charge_a)
      end do
      dv_charge_on_atoms(a) = charge_a
      sum_mulliken_charges = sum_mulliken_charges + charge_a
    end do

    wtime_end = mpi_wtime()
    call add_event('get_mulliken_charges_on_atoms:sum_charges', wtime_end - wtime_start)

    if (check_master()) then
      if (abs(sum_mulliken_charges - 1d0) > charge_sum_error_tol) then
        write(0, *) 'Warn: sum of mulliken charges on atoms is distant from one'
        write(0, '(A, E26.16e3)') ' sum_mulliken_charges: ', sum_mulliken_charges
      end if
    end if
  end subroutine get_mulliken_charges_on_atoms


  ! charges(group, k): \sum_{i \in (atom \in group)} y_{i, k}^2 / \sum_{i} y_{i, k}^2
  ! ipratios(k): \sum_{g} charges(g, k)^2 (ipratio for sqrt{\sum_{i \in (atom \in group(g))} y_{i, k}^2})
  subroutine get_eigenstate_charges_on_groups(basis, structure, group_id)
    type(wk_basis_t), intent(inout) :: basis
    type(wk_structure_t), intent(in) :: structure
    integer, intent(in) :: group_id(:, :)

    integer :: nprow, npcol, myrow, mycol, myrank
    integer :: lrindx, lcindx, rsrc, csrc
    integer :: k, g, ai, a, atom_index_start, atom_index_end, i, ierr
    real(8), allocatable :: charges_buf(:, :), ipratios_buf(:), charges(:, :)
    real(8) :: charges_sum
    integer :: indxl2g, blacs_pnum  ! Function.

    if (basis%is_group_filter_mode) then
      stop 'IMPLEMENT HERE (get_eigenstate_charges_on_groups)'
    else
      allocate(charges(size(group_id, 2), basis%Y_filtered_desc(cols_)))
      allocate(charges_buf(size(group_id, 2), basis%Y_filtered_desc(cols_)))
      allocate(ipratios_buf(basis%Y_filtered_desc(cols_)))
      call blacs_gridinfo(basis%Y_filtered_desc(context_), nprow, npcol, myrow, mycol)
      charges_buf(:, :) = 0d0
      do k = 1, basis%Y_filtered_desc(cols_)
        do g = 1, size(group_id, 2)
          do ai = 1, group_id(1, g)
            a = group_id(ai + 1, g)
            atom_index_start = structure%atom_indices(a)
            atom_index_end = structure%atom_indices(a + 1) - 1
            do i = atom_index_start, atom_index_end
              call infog2l(i, k, basis%Y_filtered_desc, nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc)
              if (myrow == rsrc .and. mycol == csrc) then
                charges_buf(g, k) = charges_buf(g, k) + basis%Y_filtered(lrindx, lcindx) ** 2d0
              end if
            end do
          end do
        end do
      end do
      charges(:, :) = 0d0
      call mpi_allreduce(charges_buf, charges, size(group_id, 2) * basis%Y_filtered_desc(cols_), &
           mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
      charges_buf(:, :) = 0d0
      ipratios_buf(:) = 0d0
      call mpi_comm_rank(mpi_comm_world, myrank, ierr)
      do k = 1, basis%Y_filtered_desc(cols_)
        if (myrank == mod(k, nprow * npcol)) then
          charges_sum = sum(charges(:, k))
          charges_buf(:, k) = charges(:, k) / charges_sum
          ipratios_buf(k) = sum(charges_buf(:, k) ** 2d0)
        end if
      end do
      charges(:, :) = 0d0
      basis%eigenstate_ipratios(:) = 0d0
      call mpi_allreduce(charges_buf, charges, size(group_id, 2) * basis%Y_filtered_desc(cols_), &
           mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(ipratios_buf, basis%eigenstate_ipratios, basis%Y_filtered_desc(cols_), &
           mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    end if
  end subroutine get_eigenstate_charges_on_groups


  ! Definition of energy from the nonlinear term when h1_type is "charge_overlap".
  ! Called after calling get_mulliken_charges_on_atoms.
  subroutine get_charge_overlap_energy(structure, charge_factor, dv_charge_on_atoms, energy)
    type(wk_structure_t), intent(in) :: structure
    type(wk_charge_factor_t), intent(in) :: charge_factor
    real(8), intent(in) :: dv_charge_on_atoms(structure%num_atoms)
    real(8), intent(out) :: energy

    integer :: i
    real(8) :: charge, factor, sum_charge_squared

    sum_charge_squared = 0d0
    do i = 1, structure%num_atoms
      charge = dv_charge_on_atoms(i)
      factor = get_charge_factor(i, structure, charge_factor)
      sum_charge_squared = sum_charge_squared + factor * (charge ** 2d0)
    end do
    energy = - 0.5d0 * sum_charge_squared
  end subroutine get_charge_overlap_energy


  real(8) function wrap_around_center(coord, unit, center)
    real(8), intent(in) :: coord, unit, center

    real(8) :: c

    c = coord
    do while (c < center - unit / 2d0)  ! Move from left outside to unitcell.
      c = c + unit
    end do
    do while (center + unit / 2d0 <= c)  ! Move from right outside to unitcell.
      c = c - unit
    end do
    wrap_around_center = c
  end function wrap_around_center


  subroutine get_mulliken_charge_coordinate_moments(structure, dv_charge_on_atoms, charge_moment)
    type(wk_structure_t), intent(in) :: structure
    real(8), intent(in) :: dv_charge_on_atoms(structure%num_atoms)
    type(wk_charge_moment_t), intent(out) :: charge_moment

    real(8) :: normalizer, ratio, coord, unit, center, coord_old
    integer :: i, j, max_charge_atom_index
    real(8) :: atom_coordinates_copy(size(structure%atom_coordinates, 1), size(structure%atom_coordinates, 2))
    logical :: is_debug = .false.

    max_charge_atom_index = maxloc(dv_charge_on_atoms, 1)  ! The i-th atom is considered as a new center.
    do j = 1, 3
      if (structure%periodic_xyz(j)) then  ! Periodic boundary condition is imposed.
        center = structure%atom_coordinates(j, max_charge_atom_index)
        unit = structure%unitcell_xyz(j)
        ! New unitcell range is center - unit / 2 <= x < center + unit / 2
        do i = 1, structure%num_atoms
          coord = structure%atom_coordinates(j, i)
          atom_coordinates_copy(j, i) = wrap_around_center(coord, unit, center)
        end do
      else  ! Non-periodic boundary condition is imposed.
        atom_coordinates_copy(j, :) = structure%atom_coordinates(j, :)
      end if
    end do

    ! Output moved atoms for debugging.
    if (is_debug .and. check_master()) then
      do i = 1, structure%num_atoms
        do j = 1, 3
          coord_old = structure%atom_coordinates(j, i)
          coord = atom_coordinates_copy(j, i)
          if (coord_old /= coord) then
            print *, 'get_mulliken_charge_coordinate_moments: move atom', i, 'coord', j, 'from', &
                 coord_old, '->', coord, '(center, unit:', structure%atom_coordinates(j, max_charge_atom_index), &
                 ',', structure%unitcell_xyz(j), ')'
          end if
        end do
      end do
    end if

    normalizer = sum(dv_charge_on_atoms(:))
    charge_moment%means(:) = 0d0
    charge_moment%msds(:) = 0d0
    do i = 1, structure%num_atoms
      ratio = dv_charge_on_atoms(i) / normalizer
      do j = 1, 3
        charge_moment%means(j) = charge_moment%means(j) + atom_coordinates_copy(j, i) * ratio
        charge_moment%msds(j) = charge_moment%msds(j) + (atom_coordinates_copy(j, i) ** 2d0) * ratio
      end do
    end do
    do j = 1, 3
      charge_moment%msds(j) = charge_moment%msds(j) - charge_moment%means(j) ** 2d0
    end do
    charge_moment%msds(4) = charge_moment%msds(1) + charge_moment%msds(2) + charge_moment%msds(3)
  end subroutine get_mulliken_charge_coordinate_moments


  subroutine multiply_coordinates_to_charge_matrix(structure, axis, order, X_desc, X)
    type(wk_structure_t), intent(in) :: structure
    integer, intent(in) :: axis  ! 1 / 2 / 3 is x / y / z, respectively.
    integer, intent(in) :: order  ! Order of moment to compute (1 or 2).
    integer, intent(in) :: X_desc(desc_size)
    real(8), intent(inout) :: X(:, :)

    integer :: a, i, j, maxlocs(X_desc(cols_))
    real(8) :: center, unit, coord, coord_original, elem

    ! Warning: searching basis with max charge, not atom.
    !          This is incompatible with get_mulliken_charge_coordinate_moments.
    call maxloc_columns(X, X_desc, maxlocs)
    do j = 1, X_desc(cols_)
      if (structure%periodic_xyz(axis)) then  ! Periodic boundary condition is imposed.
        unit = structure%unitcell_xyz(axis)
        center = structure%atom_coordinates(axis, base_index_to_atom(structure, maxlocs(j)))
      end if
      do a = 1, structure%num_atoms
        coord_original = structure%atom_coordinates(axis, a)
        if (structure%periodic_xyz(axis)) then
          coord = wrap_around_center(coord_original, unit, center)
        else
          coord = coord_original
        end if
        do i = structure%atom_indices(a), structure%atom_indices(a + 1) - 1
          call pdelget('Self', ' ', elem, X, i, j, X_desc)
          call pdelset(X, i, j, X_desc, elem * (coord ** order))
        end do
      end do
    end do
  end subroutine multiply_coordinates_to_charge_matrix


  subroutine get_msd_of_eigenstates(structure, S_sparse, basis)
    type(sparse_mat), intent(in) :: S_sparse
    type(wk_structure_t), intent(in) :: structure
    type(wk_basis_t), intent(inout) :: basis

    real(8), allocatable :: dv_psi_local(:), dv_psi(:), Y_work(:, :), Y_work_coord(:, :), sums(:)
    real(8) :: dv_charge_on_atoms(structure%num_atoms)
    integer :: dim, num_filter, i, j, print_count, ierr, m_local, n_local, axis
    type(wk_charge_moment_t) ::charge_moment

    if (basis%is_group_filter_mode) then
      stop 'IMPLEMENT HERE (get_msd_of_eigenstates)'
    else
      dim = basis%Y_filtered_desc(rows_)
      num_filter = basis%Y_filtered_desc(cols_)
      m_local = size(basis%Y_filtered, 1)
      n_local = size(basis%Y_filtered, 2)
      allocate(Y_work(m_local, n_local), Y_work_coord(m_local, n_local))
      allocate(sums(num_filter))
      ! Make matrices SY and Y .* SY where '.*' means element-wise multiplication.
      call matmul_sd_d(S_sparse, basis%Y_filtered, basis%Y_filtered_desc, 1d0, Y_work, 0d0)
      Y_work(:, :) = basis%Y_filtered(:, :) * Y_work(:, :)
      ! Normalize charge distributions.
      call sum_columns(Y_work, basis%Y_filtered_desc, sums)
      do j = 1, num_filter
        call pdscal(dim, 1d0 / sums(j), Y_work, 1, j, basis%Y_filtered_desc, 1)
      end do
      do axis = 1, 3
        ! Get mean.
        Y_work_coord(:, :) = Y_work(:, :)
        call multiply_coordinates_to_charge_matrix(structure, axis, 1, basis%Y_filtered_desc, Y_work_coord)
        call sum_columns(Y_work_coord, basis%Y_filtered_desc, sums)
        basis%eigenstate_mean(axis, :) = sums(:)
        ! Get MSD.
        Y_work_coord(:, :) = Y_work(:, :)
        call multiply_coordinates_to_charge_matrix(structure, axis, 2, basis%Y_filtered_desc, Y_work_coord)
        call sum_columns(Y_work_coord, basis%Y_filtered_desc, sums)
        basis%eigenstate_msd(axis, :) = sums(:)
        ! Use formula for varianve <(x - <x>)^2> = <x^2> - <x>^2.
        basis%eigenstate_msd(axis, :) = basis%eigenstate_msd(axis, :) - (basis%eigenstate_mean(axis, :) ** 2d0)
      end do
      do j = 1, num_filter
        basis%eigenstate_msd(4, j) = sum(basis%eigenstate_msd(1 : 3, j))
      end do
    end if
  end subroutine get_msd_of_eigenstates


  subroutine get_ipratio_of_eigenstates(basis)
    type(wk_basis_t), intent(inout) :: basis

    integer :: dim, num_filter, j, m_local, n_local, l, g, c1, c2, ierr
    real(8), allocatable :: ipratios_buf(:)
    real(8), allocatable :: Y_work_power4(:, :), Y_work_power2(:, :), sum_power4(:), sum_power2(:)
    ! Functions.
    integer :: numroc, indxl2g

    if (basis%is_group_filter_mode) then
      allocate(ipratios_buf(basis%num_basis))  ! basis%num_basis == basis%Y_local%n.
      ipratios_buf(:) = 0d0
      do l = 1, numroc(basis%Y_local%num_blocks, 1, basis%Y_local%my_rank, 0, basis%Y_local%num_procs)
        g = indxl2g(l, 1, basis%Y_local%my_rank, 0, basis%Y_local%num_procs)
        c1 = basis%Y_local%block_to_col(g)
        c2 = basis%Y_local%block_to_col(g + 1)
        ipratios_buf(c1 : c2 - 1) = sum(basis%Y_local%local_matrices(l)%val(:, :) ** 4d0, 1) / &
             (sum(basis%Y_local%local_matrices(l)%val(:, :) ** 2d0, 1) ** 2d0)
      end do
      call mpi_allreduce(ipratios_buf, basis%dv_ipratios, basis%num_basis, mpi_double_precision, mpi_sum, &
           mpi_comm_world, ierr)
    else
      dim = basis%Y_filtered_desc(rows_)
      num_filter = basis%Y_filtered_desc(cols_)
      m_local = size(basis%Y_filtered, 1)
      n_local = size(basis%Y_filtered, 2)
      allocate(Y_work_power4(m_local, n_local), Y_work_power2(m_local, n_local))
      allocate(sum_power4(num_filter), sum_power2(num_filter))
      Y_work_power4(:, :) = basis%Y_filtered(:, :) ** 4d0
      Y_work_power2(:, :) = basis%Y_filtered(:, :) ** 2d0
      call sum_columns(Y_work_power4, basis%Y_filtered_desc, sum_power4)
      call sum_columns(Y_work_power2, basis%Y_filtered_desc, sum_power2)
      basis%dv_ipratios(1 : num_filter) = sum_power4(1 : num_filter) / (sum_power2(1 : num_filter) ** 2d0)
    end if
  end subroutine get_ipratio_of_eigenstates


  real(8) function get_charge_factor(atom_i, structure, charge_factor)
    integer, intent(in) :: atom_i
    type(wk_structure_t), intent(in) :: structure
    type(wk_charge_factor_t), intent(in) :: charge_factor

    if (structure%atom_elements(atom_i) == 'H' .and. charge_factor%charge_factor_H >= 0d0) then
      get_charge_factor = charge_factor%charge_factor_H
    else if (structure%atom_elements(atom_i) == 'C' .and. charge_factor%charge_factor_C >= 0d0) then
      get_charge_factor = charge_factor%charge_factor_C
    else
      get_charge_factor = charge_factor%charge_factor_common
    end if
  end function get_charge_factor
end module wk_charge_m
