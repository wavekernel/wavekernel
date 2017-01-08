!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
module elses_mod_vis_main
!
   integer :: i_init, i_cohp
   integer :: n_interval
   integer :: iunit_for_rasmol_script
   character(len=256) :: snap_for_rasmol
!
   integer :: iunit_for_extended
!
end module elses_mod_vis_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_vis_init(n_interval2,snap_for_rasmol2)
!
   use elses_mod_md_dat,     only : itemdmx
   use elses_mod_ctrl,     only : i_verbose
   use elses_mod_vis_main, only : i_init, i_cohp, n_interval, &
&                                 iunit_for_rasmol_script, &
&                                 snap_for_rasmol, &
&                                 iunit_for_extended
!
!
   implicit none
   integer, intent(in) :: n_interval2
   character(len=*), intent(in) :: snap_for_rasmol2
!
   i_init=1
   i_cohp=1
   iunit_for_rasmol_script=0
   iunit_for_extended=0
!
!
   if (i_verbose >= 1) then 
      write(*,*)'@@ LSES_VIS_INIT'
      write(*,*)' setting: n_interval2=',n_interval2
   endif   
!
   n_interval=n_interval2
   snap_for_rasmol=trim(snap_for_rasmol2)
!
   if (i_verbose >= 1) then 
      write(*,*)' snap_for_rasmol=',trim(snap_for_rasmol)
   endif   
!
!  if ((n_interval < 0) .or. (n_interval > itemdmx)) then
   if (n_interval < 0) then
      write(*,*)'ERROR(LSES_VIS_INIT)'
      write(*,*)'  n_interval=',n_interval
      stop
   endif   
!
end subroutine elses_vis_init
!
