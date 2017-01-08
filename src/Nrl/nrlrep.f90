!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! Module for repulsive term.
!
! Note that Nrl is only header name representing name-space.
! Original NRL scheme does not have repulsive term.
!
module MNrlRepulsive
  use MNrlSystem
  implicit none

  private

  public TNrlRepulsiveBuilder
  public NrlInitialize
  public NrlRelease
  public NrlBuildRepulsiveEnergyOnAtom
  public NrlBuildRepulsiveForceOnAtom

  interface NrlInitialize
    module procedure InitRepulsiveBuilder
  end interface 

  interface NrlRelease
    module procedure ReleaseRepulsiveBuilder
  end interface

  interface NrlBuildRepulsiveEnergyOnAtom
    module procedure BuildRepulsiveEnergyOnAtom
  end interface

  interface NrlBuildRepulsiveForceOnAtom
    module procedure BuildRepulsiveForceOnAtom 
  end interface

  ! Lattice constant of Pt in Bohr unit
  real(8), parameter :: LATTICE_CONST_PT = 7.41566621376


  type TNrlRepulsiveBuilder
    real(8) :: a
    real(8) :: b
    real(8) :: gamma
    real(8) :: cut_off
  end type 

contains

  !!
  ! Initialize procedure for TNrlRepulsiveBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitRepulsiveBuilder(io_builder)
    type(TNrlRepulsiveBuilder), intent(inout) :: io_builder
    write(*,*)' InitRepulsiveBuilder'
    io_builder%a       = 0.5d0
    io_builder%b       = 1.0d0 / (LATTICE_CONST_PT / sqrt(2d0) * 0.02d0)
    io_builder%gamma   = 3.0d0
    io_builder%cut_off = LATTICE_CONST_PT / sqrt(2d0) * 0.9d0
  end subroutine

  !!
  ! Finalize procedure for TNrlRepulsiveBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseRepulsiveBuilder(io_builder)
    type(TNrlRepulsiveBuilder), intent(inout) :: io_builder

  end subroutine

  !!
  ! Get Repulsive energy of given distance
  !
  function RepulsivePairEnergy(io_builder, i_r) result(o_energy)
    real(8) :: o_energy
    type(TNrlRepulsiveBuilder), intent(inout) :: io_builder
    real(8), intent(in) :: i_r
    real(8) :: r_c, a, b, gamma

    r_c   = io_builder%cut_off
    a     = io_builder%a
    b     = io_builder%b
    gamma = io_builder%gamma

    if(i_r > r_c) then
      o_energy = 0d0
    else
      o_energy = a * (exp((b * (r_c - i_r))**gamma) - 1)
    end if
  end function

  !!
  ! Get Repulsive energy of given distance
  !
  function RepulsivePairDerivative(io_builder, i_r) result (o_deriv)
    real(8) :: o_deriv
    type(TNrlRepulsiveBuilder), intent(inout) :: io_builder
    real(8), intent(in) :: i_r
    real(8) :: r_c, a, b, gamma

    r_c   = io_builder%cut_off
    a     = io_builder%a
    b     = io_builder%b
    gamma = io_builder%gamma

    if(i_r > r_c) then
      o_deriv = 0d0
    else
      o_deriv = - a * b * gamma * (b *(r_c - i_r))**(gamma - 1)  &
                    * exp((b * (r_c - i_r))**gamma)
    end if
  end function

  !!
  ! Calculates repulsive energy of system.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildRepulsiveEnergyOnAtom(io_builder, i_atom, o_energy, o_err)
    type(TNrlRepulsiveBuilder), intent(inout) :: io_builder
    integer, intent(in) :: i_atom
    real(8), intent(out) :: o_energy
    integer, intent(out) :: o_err

    type(TNrlAtomIterator) :: the_atom_iterator
    type(TNrlAtomPair) :: the_pair

    o_energy = 0d0
    call Reset(the_atom_iterator, i_atom)
    do while(GetNextAtomPair(the_atom_iterator, the_pair))
      if (the_pair%atom1 == the_pair%atom2) cycle
      o_energy = o_energy &
                 + RepulsivePairEnergy(io_builder, the_pair%distance)
    end do
    o_err = 0
  end subroutine

  !!
  ! Calculates repulsive force on atom
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildRepulsiveForceOnAtom(io_builder, i_atom, o_force, o_err)
    type(TNrlRepulsiveBuilder), intent(inout) :: io_builder
    integer, intent(in) :: i_atom
    real(8), intent(out) :: o_force(3)
    integer, intent(out) :: o_err
    
    type(TNrlAtomIterator) :: the_atom_iterator
    type(TNrlAtomPair) :: the_pair
    real(8) :: the_d

    o_force(1:3) = 0d0
    call Reset(the_atom_iterator, i_atom)
    do while(GetNextAtomPair(the_atom_iterator, the_pair))
      if (the_pair%atom1 == the_pair%atom2) cycle
      the_d = RepulsivePairDerivative(io_builder, the_pair%distance)
      ! note the sign of followings. (l, m, n) = (r2 - r1) / |r2 - r1|
      o_force(1) = o_force(1) + the_d * the_pair%l 
      o_force(2) = o_force(2) + the_d * the_pair%m
      o_force(3) = o_force(3) + the_d * the_pair%n
    end do

    o_err = 0
  end subroutine

end module
