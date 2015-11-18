!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
!! Copyright (C) ELSES. 2007-2015 all rights reserved
module M_00_v_info
!
  implicit none
  character(len=32), parameter :: version_info='ELSES version 0.05.03'
!
  private
  public :: elses_00_version_info 
  public :: version_info 
!
  contains
!
  subroutine elses_00_version_info
!
   write(*,'(a)') trim(version_info)
!
  end subroutine elses_00_version_info
!
end module M_00_v_info
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
