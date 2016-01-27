!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! module for TNrlPairMatrix
!
module MNrlPairMatrix
  implicit none

  private

  public TNrlPairMatrix
  public TNrlPairForceMatrix

  ! public procedures for TNrlPairMatrix
  public NrlSetIdentity
  public NrlWriteTo
  public NrlSetSize

  !!
  ! Stores the interaction between two atoms.
  !
  type TNrlPairMatrix
    integer :: row_count
    integer :: col_count
    real(8) :: e(16, 16)
  end type

  !!
  ! Stores the xyz vector between two atoms.
  !
  type TNrlPairForceMatrix
    integer :: row_count
    integer :: col_count
    real(8) :: e(3, 16, 16)
  end type

  interface NrlSetIdentity
    module procedure SetIdentity
  end interface

  interface NrlWriteTo
    module procedure WriteTo
  end interface

  interface NrlSetSize
    module procedure SetSizeOfPairMatrix
    module procedure SetSizeOfPairForceMatrix
  end interface

contains

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetIdentity(iomatrix, isize)
    type(TNrlPairMatrix), intent(inout) :: iomatrix
    integer, intent(in) :: isize
    integer :: i

    iomatrix%row_count = isize
    iomatrix%col_count = isize

    iomatrix%e(1:isize, 1:isize) = 0d0
    do i = 1, isize
      iomatrix%e(i, i) = 1d0
    end do
  end subroutine

  !!
  ! Set size of TNrlPairMatrix
  ! Currently this routine only set two size.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetSizeOfPairMatrix(iomatrix, irow, icol)
    type(TNrlPairMatrix), intent(inout) :: iomatrix
    integer, intent(in) :: irow
    integer, intent(in) :: icol

    iomatrix%row_count = irow
    iomatrix%col_count = icol
  end subroutine

  !!
  ! Set size of TNrlPairForceMatrix
  ! Currently this routine only set two size.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetSizeOfPairForceMatrix(iomatrix, irow, icol)
    type(TNrlPairForceMatrix), intent(inout) :: iomatrix
    integer, intent(in) :: irow
    integer, intent(in) :: icol

    iomatrix%row_count = irow
    iomatrix%col_count = icol
  end subroutine

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine WriteTo(imatrix, iunit)
    type(TNrlPairMatrix), intent(in) :: imatrix
    integer, intent(in) :: iunit
    integer :: i, j

    do i = 1, imatrix%col_count
      do j = 1, imatrix%row_count
        if (abs(imatrix%e(j, i)) > 1e-12) then 
          write(iunit, "(4I6,F12.4)") i, j, imatrix%e(j, i)
        end if
      end do
    end do
  end subroutine
end module
