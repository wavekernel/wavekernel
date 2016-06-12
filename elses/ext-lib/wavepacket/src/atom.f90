module wp_atom_m
  use mpi
  use wp_descriptor_parameters_m
  use wp_global_variables_m
  implicit none

  private
  public :: wp_structure_t, read_structure, multiply_phase_factors, &
       make_dummy_structure, read_group_id, read_group_id_header, make_dummy_group_id, print_group_id, &
       bcast_structure, bcast_group_id, lcao_index_to_atom_index

  type wp_structure_t
    integer :: num_atoms
    integer, allocatable :: atom_indices(:)
    real(8), allocatable :: atom_coordinates(:, :)
    character, allocatable :: atom_elements(:)
  end type wp_structure_t

contains

  subroutine read_structure_skip_steps(iunit, step)
    integer, intent(in) :: iunit, step
    integer :: i, j, num_atoms

    do i = 1, step - 1
      read(iunit, *) num_atoms
      read(iunit, *)  ! Skip a line.
      do j = 1, num_atoms
        read(iunit, *)
      end do
    end do
  end subroutine read_structure_skip_steps


  ! 原子のインデックスと LCAO 係数のインデックスの対応を xyz ファイルから読む.
  ! a 番目の原子のインデックスは,
  ! atom_indices(a) から atom_indices(a + 1) - 1 までとなる.
  ! 例えば xyz ファイルに C H C (H, C の基底数はそれぞれ 1, 4)
  ! という順番で原子が記述されていたとすると, それぞれの LCAO 係数でのインデックスは
  ! 1 - 4, 5, 6 - 9 であるから,
  ! atom_indices = (/1, 5, 6, 10/) となる.

  ! read_atom_indices
  subroutine read_structure(filename, step, structure)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: step
    type(wp_structure_t), intent(out) :: structure

    integer, parameter :: iunit = 12
    integer :: i, n, atom_valence
    double precision :: coordinate(3)  ! x, y, z.

    open(iunit, file=filename)

    call read_structure_skip_steps(iunit, step)

    read(iunit, *) n
    structure%num_atoms = n
    allocate(structure%atom_indices(n + 1), structure%atom_coordinates(3, n), structure%atom_elements(n))
    structure%atom_indices(1) = 1

    read(iunit, *)  ! Skip a line.
    do i = 1, n
      read(iunit, *) structure%atom_elements(i), coordinate(1 : 3)
      structure%atom_coordinates(1 : 3, i) = coordinate(1 : 3) * kAuPerAngstrom
      if (structure%atom_elements(i) == 'H') then
        atom_valence = 1
      else if (structure%atom_elements(i) == 'C' .or. structure%atom_elements(i) == 'O') then
        atom_valence = 4
      else
        print *, 'unknown atom name!!! ', structure%atom_elements(i)
        stop
      end if
      structure%atom_indices(i + 1) = structure%atom_indices(i) + atom_valence
    end do

    close(iunit)
  end subroutine read_structure


  subroutine multiply_phase_factors(structure, phase_factor_coef, dv_psi_in, dv_psi_out)
    type(wp_structure_t), intent(in) :: structure
    real(8), intent(in) :: phase_factor_coef
    complex(kind(0d0)), intent(in) :: dv_psi_in(:)
    complex(kind(0d0)), intent(out) :: dv_psi_out(:)

    integer :: atom, atom_index_start, atom_index_end
    complex(kind(0d0)) :: phase_factor

    do atom = 1, structure%num_atoms
      atom_index_start = structure%atom_indices(atom)
      atom_index_end = structure%atom_indices(atom + 1) - 1
      phase_factor = exp(kImagUnit * phase_factor_coef * structure%atom_coordinates(1, atom))
      dv_psi_out(atom_index_start : atom_index_end) = phase_factor_coef * dv_psi_in(atom_index_start : atom_index_end)
    end do
  end subroutine multiply_phase_factors

  ! make_dummy_atom_indices_and_coordinates
  subroutine make_dummy_structure(dim, structure)
    integer, intent(in) :: dim
    type(wp_structure_t), intent(out) :: structure

    integer :: i

    structure%num_atoms = dim
    allocate(structure%atom_indices(dim + 1), structure%atom_coordinates(3, dim), structure%atom_elements(dim))
    do i = 1, dim + 1
      structure%atom_indices(i) = i
    end do
    do i = 1, dim
      structure%atom_coordinates(:, i) = (/dble(i), 0d0, 0d0/)
      structure%atom_elements(i) = 'H'
    end do
  end subroutine make_dummy_structure


  ! Read mapping: (group -> [orbital index]).
  subroutine read_group_id(filename, group_id)
    character(len=*), intent(in) :: filename
    integer, allocatable, intent(out) :: group_id(:, :)
    integer, parameter :: iunit = 16
    integer :: num_atoms, num_groups, max_group_size
    integer :: i, atom, group
    character(len=1024) :: line
    logical :: is_setting_numbers_read
    open(iunit, file=filename, status='old')
    is_setting_numbers_read = .false.
    do while (.not. is_setting_numbers_read)
      read(iunit, '(A)') line
      if (line(1 : 1) == '#') then
        cycle
      else
        read(line, *) num_atoms, num_groups, max_group_size
        ! The first row of each group is the number of atoms in the group.
        allocate(group_id(max_group_size + 1, num_groups))
        group_id(:, :) = 0
        is_setting_numbers_read = .true.
      end if
    end do
    i = 1
    do while (i <= num_atoms)
      read(iunit, '(A)') line
      if (line(1 : 1) == '#') then
        cycle
      else
        read(line, *) atom, group
        group_id(1, group) = group_id(1, group) + 1
        group_id(group_id(1, group) + 1, group) = atom
        i = i + 1
      end if
    end do
    close(iunit)
  end subroutine read_group_id


  subroutine read_group_id_header(filename, num_atoms, num_groups, max_group_size)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: num_atoms, num_groups, max_group_size
    integer, parameter :: iunit = 18
    character(len=1024) :: line
    logical :: is_setting_numbers_read = .false.
    open(iunit, file=filename, status='old')
    do while (.not. is_setting_numbers_read)
      read(iunit, '(A)') line
      if (line(1 : 1) == '#') then
        cycle
      else
        read(line, *) num_atoms, num_groups, max_group_size
        is_setting_numbers_read = .true.
      end if
    end do
    close(iunit)
  end subroutine read_group_id_header


  subroutine make_dummy_group_id(num_atoms, group_id)
    integer, intent(in) :: num_atoms
    integer, allocatable, intent(out) :: group_id(:, :)
    integer :: i

    allocate(group_id(num_atoms + 1, 1))
    group_id(1, 1) = num_atoms
    do i = 1, num_atoms
      group_id(i + 1, 1) = i
    end do
  end subroutine make_dummy_group_id


  subroutine print_group_id(group_id)
    integer, intent(in) :: group_id(:, :)
    integer :: num_groups, i, j

    num_groups = size(group_id, 2)
    do i = 1, num_groups
      write(*, '(I0, "(", I0, "): ")', advance='no') i, group_id(1, i)
      do j = 1, group_id(1, i)
        write(*, '(I0, " ")', advance='no') group_id(j + 1, i)
      end do
      print *
    end do
  end subroutine print_group_id


  ! bcast_atom_indices
  subroutine bcast_structure(root, structure)
    integer, intent(in) :: root
    type(wp_structure_t), intent(inout) :: structure
    integer :: my_rank, ierr, n

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    call mpi_bcast(structure%num_atoms, 1, mpi_integer, root, mpi_comm_world, ierr)
    n = structure%num_atoms
    if (my_rank /= root) then
      allocate(structure%atom_indices(n + 1), structure%atom_coordinates(3, n), structure%atom_elements(n))
    end if
    call mpi_bcast(structure%atom_indices, n + 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(structure%atom_coordinates, 3 * n, mpi_double_precision, root, mpi_comm_world, ierr)
    call mpi_bcast(structure%atom_elements, n, mpi_character, root, mpi_comm_world, ierr)
  end subroutine bcast_structure


  subroutine bcast_group_id(root, group_id)
    integer, intent(in) :: root
    integer, allocatable, intent(inout) :: group_id(:, :)
    integer :: max_group_size, num_groups, my_rank, ierr

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (my_rank == root) then
      max_group_size = size(group_id, 1) - 1
      num_groups = size(group_id, 2)
    end if
    call mpi_bcast(max_group_size, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(num_groups, 1, mpi_integer, root, mpi_comm_world, ierr)
    if (my_rank /= root) then
      allocate(group_id(max_group_size + 1, num_groups))
    end if
    call mpi_bcast(group_id, (max_group_size + 1) * num_groups, mpi_integer, root, mpi_comm_world, ierr)
  end subroutine bcast_group_id


  subroutine lcao_index_to_atom_index(structure, lcao_index, atom_index)
    type(wp_structure_t), intent(in) :: structure
    integer, intent(in) :: lcao_index
    integer, intent(out) :: atom_index

    integer :: a_min, a, a_max

    a_min = 1
    a_max = structure%num_atoms
    do while (.true.)
      a = (a_min + a_max) / 2
      if (structure%atom_indices(a) <= lcao_index .and. lcao_index < structure%atom_indices(a + 1)) then
        exit
      else if (lcao_index < structure%atom_indices(a)) then
        a_max = min(a, a_max - 1)
      else
        a_min = max(a, a_min + 1)
      end if
    end do
    atom_index = a
  end subroutine lcao_index_to_atom_index
end module wp_atom_m
