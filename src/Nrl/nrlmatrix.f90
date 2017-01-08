!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! Director of Hopping and on-site builder module.
!
module MNrlMatrix
  use MNrlPairMatrix
  use MNrlOnSite
  use MNrlHopping
  use MNrlSystem

  implicit none
  logical, parameter :: suffix_check_in_MNrlMatrix= .true.

  private

  ! public types
  public TNrlMatrixBuilder

  ! public procedures for TMatrixBuilder
  public NrlInitialize
  public NrlRelease
  public NrlBuildHamiltonianMatrix
  public NrlBuildOverlapMatrix
  public NrlBuildHamiltonianForceOnAtom
  public NrlBuildOverlapForceOnAtom

  type TNrlMatrixBuilder
    type(TNrlOnSiteBuilder) :: onsite
    type(TNrlHopBuilder)    :: hop
    type(TNrlOnSiteBuilder) :: overlap_onsite
    type(TNrlHopBuilder)    :: overlap
    real(8)                 :: cut_off
  end type

  interface NrlInitialize
    module procedure InitMatrixBuilder
  end interface

  interface NrlRelease
    module procedure ReleaseMatrixBuilder
  end interface

  interface NrlBuildHamiltonianMatrix
    module procedure BuildHamiltonianMatrix
  end interface

  interface NrlBuildOverlapMatrix
    module procedure BuildOverlapMatrix
  end interface

  interface NrlBuildHamiltonianForceOnAtom
    module procedure BuildHamiltonianForceOnAtom
  end interface

  interface NrlBuildOverlapForceOnAtom
    module procedure BuildOverlapForceOnAtom
  end interface

contains

  !!
  ! Initializes TNrlMatrixBuilder object.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitMatrixBuilder(io_builder)
    type(TNrlMatrixBuilder), intent(inout) :: io_builder

    io_builder%cut_off = 0d0
    call NrlInitialize(io_builder%onsite)
    call NrlInitialize(io_builder%hop)
    call NrlInitialize(io_builder%overlap_onsite)
    call NrlInitialize(io_builder%overlap)
  end subroutine

  !!
  ! Finalizes TNrlMatrixBuilder object.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseMatrixBuilder(io_builder)
    type(TNrlMatrixBuilder), intent(inout) :: io_builder

    call NrlRelease(io_builder%onsite)
    call NrlRelease(io_builder%hop)
    call NrlRelease(io_builder%overlap_onsite)
    call NrlRelease(io_builder%overlap)
  end subroutine


  !!
  ! Adds local matrix element to global matrix.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AddPairMatrix(io_global_matrix, i_pair, i_matrix)
    real(8) :: io_global_matrix(:,:,:,:)
    type(TNrlAtomPair), intent(in) :: i_pair
    type(TNrlPairMatrix), intent(inout) :: i_matrix

    integer :: atom1, atom2, row1, row2, i, j

    atom1 = i_pair%atom1
    atom2 = i_pair%atom_n2
    row1 = i_matrix%col_count
    row2 = i_matrix%row_count

    io_global_matrix(1:row2, 1:row1, atom2, atom1) = i_matrix%e(1:row2, 1:row1)
  end subroutine

  !!
  ! Makes a Hamiltonian matrix from tight-binding parameters.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildHamiltonianMatrix(io_builder, io_global_matrix, o_err) 
    use MNrlVariables, only: noav
    type(TNrlMatrixBuilder), intent(inout) :: io_builder
    real(8), intent(inout) :: io_global_matrix(:,:,:,:)
    integer, intent(out)   :: o_err
    type(TNrlAtomIterator) :: the_atom_iterator
    type(TNrlAtomPair)     :: the_pair
    type(TNrlPairMatrix)   :: the_matrix
    integer :: jsv2

    do jsv2 = 1, noav
      call Reset(the_atom_iterator, jsv2)
      do while(GetNextAtomPair(the_atom_iterator, the_pair))
        if (the_pair%atom1 == the_pair%atom2) then
          call NrlGetMatrix(io_builder%onsite, the_pair%atom1, the_matrix, o_err)
          if(o_err /= 0) return
        else
!          the_matrix%e(:,:) = 0d0
          call NrlGetMatrix(io_builder%hop, the_pair, the_matrix)
        end if
        call AddPairMatrix(io_global_matrix, the_pair, the_matrix)
!       write(63, *) "atom: ", the_pair%atom1, the_pair%atom2
!       call WriteTo(imatrix, 63)        
      end do
    end do
    o_err = 0
  end subroutine

  !!
  ! Makes a overlap matrix from tight-binding parameters.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildOverlapMatrix(io_builder, io_global_matrix, o_err) 
    use MNrlVariables, only: noav
    use MNrlSystem
    type(TNrlMatrixBuilder), intent(inout) :: io_builder
    real(8), intent(inout) :: io_global_matrix(:,:,:,:)
    integer, intent(out)   :: o_err
    type(TNrlAtomIterator) :: the_atom_iterator
    type(TNrlAtomPair)     :: the_pair
    type(TNrlPairMatrix)   :: the_matrix
    integer :: jsv2

    do jsv2 = 1, noav
      call Reset(the_atom_iterator, jsv2)
      do while(GetNextAtomPair(the_atom_iterator, the_pair))
        if (the_pair%atom1 == the_pair%atom2) then
          call NrlGetMatrix(io_builder%overlap_onsite, the_pair%atom1, &
                            the_matrix, o_err)
          if(o_err /= 0) return
        else
!          the_matrix%e(:,:) = 0d0
          call NrlGetMatrix(io_builder%overlap, the_pair, the_matrix)
        end if
        call AddPairMatrix(io_global_matrix, the_pair, the_matrix)
      end do
    end do
    o_err = 0
  end subroutine

  !!
  ! Calculates sum (B_jj (d H_jj)/(d r_i))  and adds to force.
  ! i, j are the indexes of atom.
  !           
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AddDerivativeToOnSiteForce(i_density_matrix, i_deriv, &
      i_pair, io_force)
    real(8), intent(in) :: i_density_matrix(:,:,:,:)
    type(TNrlPairForceMatrix), intent(in) :: i_deriv
    type(TNrlAtomPair) :: i_pair
    real(8), intent(inout) :: io_force(:)

    integer :: atom, i, j

    atom = i_pair%atom2
    do i = 1, i_pair%orbit_count2
      do j = 1, i_pair%orbit_count2
        ! Note that index "one" represents on-site but this is not 
        ! a guaranteed behavior. 
        ! TODO: check the onsite index.
        if (suffix_check_in_MNrlMatrix) then
          if (i > size(i_deriv%e,3)) then
            write(*,*)'ERROR(AddDerivativeToOffSiteForce):i=',i 
            stop
          endif   
          if (j > size(i_deriv%e,2)) then
            write(*,*)'ERROR(AddDerivativeToOffSiteForce):j=',j 
            stop
          endif   
        endif  
        io_force(1:3) = io_force(1:3) - &
          i_density_matrix(j, i, 1, atom) * i_deriv%e(1:3, j, i)
      end do
    end do
  end subroutine

  !!
  ! Calculates sum (B_ij d(H_ij)/d(r_i)) and adds to force.
  ! i, j are the indexes of atom.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AddDerivativeToOffSiteForce(i_density_matrix, i_deriv, &
      i_pair, io_force)
    real(8), intent(in) :: i_density_matrix(:,:,:,:)
    type(TNrlPairForceMatrix), intent(in) :: i_deriv
    type(TNrlAtomPair) :: i_pair
    real(8), intent(inout) :: io_force(:)

    integer :: atom1, atom2, i, j

    atom1 = i_pair%atom1
    atom2 = i_pair%atom_n2    
    do i = 1, i_pair%orbit_count1
      do j = 1, i_pair%orbit_count2
        ! two represents the both effect of V_ij and V_ji
        if (suffix_check_in_MNrlMatrix) then
          if (i > size(i_deriv%e,3)) then
            write(*,*)'ERROR(AddDerivativeToOffSiteForce):i=',i 
            stop
          endif   
          if (j > size(i_deriv%e,2)) then
            write(*,*)'ERROR(AddDerivativeToOffSiteForce):j=',j 
            stop
          endif   
        endif  
        io_force(1:3) = io_force(1:3) - &
          2.0d0 * i_density_matrix(j, i, atom2, atom1) * i_deriv%e(1:3, j, i)
      end do
    end do
  end subroutine

  !!
  ! Calculates Hamiltonian part.in the force on an atom.
  !     
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildHamiltonianForceOnAtom(io_builder, i_density_matrix, &
                       i_atom, o_force, o_err)
    type(TNrlMatrixBuilder), intent(inout) :: io_builder
    real(8), intent(in) :: i_density_matrix(:,:,:,:)
    integer, intent(in) :: i_atom
    real(8), intent(out) :: o_force(:)
    integer, intent(out) :: o_err
    type(TNrlPairForceMatrix) :: the_deriv

    type(TNrlAtomIterator) :: the_atom_iterator
    type(TNrlAtomPair) :: the_pair
    
    o_force(1:3) = 0d0

    call Reset(the_atom_iterator, i_atom)
    do while(GetNextAtomPair(the_atom_iterator, the_pair))
      ! on-site
      call NrlGetDerivative(io_builder%onsite, the_pair, the_deriv)
      call AddDerivativeToOnSiteForce(i_density_matrix, the_deriv, &
                                        the_pair, o_force)

      if (the_pair%atom1 /= the_pair%atom2) then

        ! off-site
        call NrlGetDerivative(io_builder%hop, the_pair, the_deriv)
        call AddDerivativeToOffSiteForce(i_density_matrix, the_deriv, &
                                         the_pair, o_force)
      end if
    end do
    ! print "(A,I6,3G12.4)", "HOPPING FORCE", the_pair%atom1, o_force(1:3)
    o_err = 0
  end subroutine

  !!
  ! Calculates overlap part in the force on an atom.
  ! the sign of the contribution to the force by the sum S_ij P_ij is 
  ! negative where
  !     S_ij:  overlap matrix 
  !     P_ij:  sum_n {|c_i> f_n e_n <c_j|} 
  !     f_n: occupation   
  !     e_n: eigen energy value.
  !     
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildOverlapForceOnAtom(io_builder, i_density_matrix, &
                       i_atom, o_force, o_err)
    type(TNrlMatrixBuilder), intent(inout) :: io_builder
    real(8), intent(in) :: i_density_matrix(:,:,:,:)
    integer, intent(in) :: i_atom
    real(8), intent(out) :: o_force(:)
    integer, intent(out) :: o_err
    type(TNrlPairForceMatrix) :: the_deriv

    type(TNrlAtomIterator) :: the_atom_iterator
    type(TNrlAtomPair) :: the_pair
    
    o_force(1:3) = 0d0

    ! Off-site calculation
    call Reset(the_atom_iterator, i_atom)
    do while(GetNextAtomPair(the_atom_iterator, the_pair))
      if (the_pair%atom1 /= the_pair%atom2) then
        call NrlGetDerivative(io_builder%overlap, the_pair, the_deriv)
        call AddDerivativeToOffSiteForce(i_density_matrix, the_deriv, &
                                         the_pair, o_force)
      end if
    end do
    o_force(1:3) = - o_force(1:3)
    ! print "(A,I6,3G12.4)", "OVERLAP FORCE", the_pair%atom1, o_force(1:3)
    o_err = 0
  end subroutine
end module
