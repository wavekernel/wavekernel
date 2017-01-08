!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! module for NRL atom loop 
!
module MNrlSystem
  implicit none

  private

  public TNrlAtomPair
  public TNrlAtomIterator

  public Initialize
  public Release
  public Reset
  public GetNextAtomPair

  interface Initialize
    module procedure NrlInitAtomIterator
  end interface

  interface Release
    module procedure NrlReleaseAtomIterator
  end interface

  interface Reset
    module procedure NrlResetAtomIterator
  end interface

  interface GetNextAtomPair
    module Procedure NrlGetNextAtomPair
  end interface

  type TNrlAtomIterator
    integer, pointer :: neighbor_list(:, :)
    integer, pointer :: neighbor_count(:)
    integer, pointer :: neighbor_atom(:, :)
    integer, pointer :: sub_to_global(:)
    integer, pointer :: species(:)
    integer, pointer :: orbit_count(:)
    integer :: atom1
    integer :: atom2
    logical :: end_of_loop
  end type

  !!
  ! @param atom1   global index of atom1
  ! @param atom2   global index of atom2
  ! @param atom_n2 index of atom2 as a neighbor of atom1
  !
  type TNrlAtomPair
    integer :: atom1
    integer :: atom2
    integer :: atom_n2
    integer :: species1
    integer :: species2
    integer :: orbit_count1
    integer :: orbit_count2
    integer :: row_to_lm1(16)
    integer :: row_to_lm2(16)
    integer :: il_size1
    integer :: il_size2
    integer :: row_to_il1(16)
    integer :: row_to_il2(16)
    real(8) :: distance
    real(8) :: l, m, n
  end type

  integer, parameter :: NRL_SUPPORTED_ROW_TO_LM(16) = &
      (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15/)

  ! Note that this il is the index of angular momentum, 
  ! not angular momentum itself.
  integer, parameter :: NRL_SUPPORTED_IL_SIZE = 3
  integer, parameter :: NRL_SUPPORTED_ROW_TO_IL(16) = &
      (/1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4/)

contains

  !!
  ! Initialize routine for TNRLAtomIterator
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlInitAtomIterator(io_atom_iterator)
    type(TNrlAtomIterator), intent(inout) :: io_atom_iterator
    io_atom_iterator%atom1 = 0
    nullify(io_atom_iterator%species)
    nullify(io_atom_iterator%neighbor_count)
    nullify(io_atom_iterator%neighbor_list)
    nullify(io_atom_iterator%sub_to_global)
    nullify(io_atom_iterator%orbit_count)
    nullify(io_atom_iterator%species)
  end subroutine

  !!
  ! Finalize routine for TNRLAtomIterator
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlReleaseAtomIterator(io_atom_iterator)
    type(TNrlAtomIterator), intent(inout) :: io_atom_iterator
    nullify(io_atom_iterator%species)
    nullify(io_atom_iterator%neighbor_count)
    nullify(io_atom_iterator%neighbor_list)
    nullify(io_atom_iterator%sub_to_global)
    nullify(io_atom_iterator%orbit_count)
    nullify(io_atom_iterator%species)
  end subroutine

  !!
  ! Resets the counter of loop
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlResetAtomIterator(io_atom_iterator, atom1)
    use MNrlVariables, only: js4jsv, jsv4jsd, njsd, jsv4jsd, jsei, nval
    type(TNrlAtomIterator), intent(inout) :: io_atom_iterator
    integer, intent(in) :: atom1
    io_atom_iterator%atom1 = atom1
    io_atom_iterator%atom2 = 1
    io_atom_iterator%sub_to_global  => js4jsv
    io_atom_iterator%neighbor_list  => jsv4jsd
    io_atom_iterator%neighbor_count => njsd(:, 1)
    io_atom_iterator%neighbor_atom  => jsv4jsd
    io_atom_iterator%species        => jsei
    io_atom_iterator%orbit_count    => nval
    io_atom_iterator%end_of_loop    =  .false.
  end subroutine

  !!
  ! Update information of atomic pair such as distance
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine UPdateDistance(io_atom_iterator, io_pair)
!   use MNrlVariables, only: iperiodic, tx, ty, tz, ax, ay, az
    use MNrlVariables, only: i_pbc_x, i_pbc_y, i_pbc_z, tx, ty, tz, ax, ay, az
    type(TNrlAtomIterator), intent(inout) :: io_atom_iterator
    type(TNrlAtomPair), intent(inout) :: io_pair
    integer :: i
    real(8) :: dvec(3)
    dvec(1) = tx(io_pair%atom2) - tx(io_pair%atom1)
    dvec(2) = ty(io_pair%atom2) - ty(io_pair%atom1)
    dvec(3) = tz(io_pair%atom2) - tz(io_pair%atom1)
    
!   if (iperiodic == 1) then
!     do i = 1, 3
!       dvec(i) = modulo(dvec(i) + 0.5d0, 1.0d0) - 0.5d0
!     end do
!   end if

    if (i_pbc_x == 1) dvec(1) = modulo(dvec(1) + 0.5d0, 1.0d0) - 0.5d0
    if (i_pbc_y == 1) dvec(2) = modulo(dvec(2) + 0.5d0, 1.0d0) - 0.5d0
    if (i_pbc_z == 1) dvec(3) = modulo(dvec(3) + 0.5d0, 1.0d0) - 0.5d0
 
    dvec(1) = dvec(1) * ax
    dvec(2) = dvec(2) * ay
    dvec(3) = dvec(3) * az 

    io_pair%distance = sqrt(dvec(1) * dvec(1) + dvec(2) * dvec(2) + dvec(3) * dvec(3))
    if (io_pair%distance > 1d-12) then
      io_pair%l = dvec(1) / io_pair%distance
      io_pair%m = dvec(2) / io_pair%distance
      io_pair%n = dvec(3) / io_pair%distance
    else
      io_pair%l = 0d0
      io_pair%m = 0d0
      io_pair%n = 1d0
    end if
    !print "(3F12.4)", io_pair%l, io_pair%m, io_pair%n
  end subroutine
    
  !!
  ! Increments the counter of atom loop and return the atom pair. 
  ! Note that this function returns not only off-site pairs but also
  ! the on-site pair.
  !
  logical function NrlGetNextAtomPair(io_atom_iterator, opair) result(oresult)
    type(TNrlAtomIterator), intent(inout) :: io_atom_iterator
    type(TNrlAtomPair), intent(out) :: opair

    integer :: atom_s1, atom_s2, atom_n2, species1, species2
    if (io_atom_iterator%end_of_loop) then
      oresult = .false.
      return
    end if

    atom_s1 = io_atom_iterator%atom1
    atom_n2 = io_atom_iterator%atom2
    atom_s2 = io_atom_iterator%neighbor_atom(atom_n2, atom_s1)

    opair%atom1 = io_atom_iterator%sub_to_global(atom_s1)
    opair%atom2 = io_atom_iterator%sub_to_global(atom_s2)
    opair%atom_n2 = atom_n2
    species1 = io_atom_iterator%species(opair%atom1)
    species2 = io_atom_iterator%species(opair%atom2)
    opair%species1 = species1
    opair%species2 = species2
    opair%orbit_count1 = io_atom_iterator%orbit_count(species1)
    opair%orbit_count2 = io_atom_iterator%orbit_count(species2)
    opair%row_to_lm1(1:opair%orbit_count1) = &
        NRL_SUPPORTED_ROW_TO_LM(1:opair%orbit_count1)
    opair%row_to_lm2(1:opair%orbit_count2) = &
        NRL_SUPPORTED_ROW_TO_LM(1:opair%orbit_count2)
    opair%il_size1 = NRL_SUPPORTED_IL_SIZE
    opair%il_size2 = NRL_SUPPORTED_IL_SIZE
    opair%row_to_il1(1:opair%orbit_count1) = &
        NRL_SUPPORTED_ROW_TO_IL(1:opair%orbit_count1)
    opair%row_to_il2(1:opair%orbit_count2) = &
        NRL_SUPPORTED_ROW_TO_IL(1:opair%orbit_count2)

    call UpdateDistance(io_atom_iterator, opair)

    call Increment(io_atom_iterator)
    oresult = .true.
  end function

  !!
  ! Increments atom loop count.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine Increment(io_atom_iterator)
    type(TNrlAtomIterator), intent(inout), target :: io_atom_iterator
    integer, pointer :: atom2

    if (io_atom_iterator%end_of_loop) return

    atom2 => io_atom_iterator%atom2
    atom2 = atom2 + 1

    if (atom2 > io_atom_iterator%neighbor_count(io_atom_iterator%atom1)) then
      io_atom_iterator%end_of_loop = .true.
    end if  
  end subroutine
end module
