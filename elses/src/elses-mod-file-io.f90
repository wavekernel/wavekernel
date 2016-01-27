!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!zzz  @@@@ elses-mod-file_io.f90 @@@@@
!zzz  @@@@@ 2007/12/29 @@@@@
!2007!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!1229: T.Hoshi; Prepared (v.0.0.0a-2007_12_29_hoshi)
!        (Copied from a file written by Susumu Yamamoto)
!2008!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!0101T.Hoshi; Module name is changed into 'elses_mod_file_io'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Find vacant unit
!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
module elses_mod_file_io

   private
   public :: vacant_unit

   contains

   integer function vacant_unit()
   integer, parameter :: MAX_UNIT=99
   integer            :: IUNIT
   logical            :: not_found

   do IUNIT=1, MAX_UNIT
      inquire(IUNIT,OPENED=not_found)
      if(.not. not_found) exit
   end do
   if(not_found) then
      write(*,'("file_io%vacant_unit: all UNIT is busy.")')
      stop
   end if
!  write(*,'("file_io%vacant_unit: IUNIT=",I3)') IUNIT
   vacant_unit=IUNIT

   end function vacant_unit

end module elses_mod_file_io

