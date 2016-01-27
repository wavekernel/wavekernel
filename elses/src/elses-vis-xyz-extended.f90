!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!zzz  @@@@ elses-vis-exteded.f90 @@@@@
!zzz  @@@@@ 2008/01/03 @@@@@
!ccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0103T.Hoshi: Prepared (NT07E-122p29)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_vis_xyz_extended(i_extended)
!    Store the structure data in XYZ format with bond analysis
!        ----> 'structure-bond-data.txt'
!
   use elses_mod_sel_sys,    only : c_system
   use elses_mod_phys_const, only : angst
   use elses_mod_md_dat,   only : itemd
!  use elses_mod_sim_cell, only : iperiodic, noa, nos, ax, ay, az
   use elses_mod_sim_cell, only : noa, nos, ax, ay, az
   use elses_mod_tx,       only : tx, ty, tz, jsei
   use elses_mod_txp,      only : txp, typ, tzp
   use elses_mod_ctrl,     only : i_verbose
   use elses_mod_file_io,  only : vacant_unit
   use elses_mod_vis_main, only : i_cohp, iunit_for_extended
   use elses_mod_elem_name, only : elem_name
!
   implicit none
   integer :: iunit, n_bond
   character(len=32) :: fname, c_name
   integer           :: js, nss
   real(8)           :: txd, tyd, tzd, ddd
   integer           :: i_extended
!
!
   fname='structure-bond-data.txt'
!
   i_extended=1
!
   if (i_verbose >= 1) then 
      write(*,*)'@@ LSES_VIS_XYZ_EXTENDED:itemd=',itemd
      write(*,*)' i_extended',i_extended
   endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Prepare the batch file and setting file
!     at the first time
!  
   if (iunit_for_extended == 0) then
     iunit=vacant_unit()
     iunit_for_extended=iunit
     if (i_verbose >= 1) then 
       write(*,*)' get vacant unit:iunit_for_extended'
       write(*,*)'     =', iunit_for_extended
     endif   
     open(iunit,file=fname,status='replace')
     close(iunit)
   endif   
!
!
   iunit=iunit_for_extended
   open(iunit,file=fname,position='append',status='old')
!
   write(iunit, *) noa
   write(iunit, *) 'mdstep=',itemd
   ddd=0.0d0
   do js=1,noa
     nss=jsei(js)
     if ((nss <= 0) .or. (nss > nos)) then
        write(*,*)'ERROR(elses-vis-rasmol)'
        write(*,*)'js, nss=',js,nss
        stop
     endif   
     c_name=elem_name(nss)
     txd=tx(js)*ax*angst
     tyd=ty(js)*ay*angst
     tzd=tz(js)*az*angst
     write(iunit, "(a4,3f23.15)") trim(c_name),txd,tyd,tzd
   enddo   
!
   close(iunit)
!
   if (i_extended == 1) then
     n_bond=-1
     call elses_ana_icohp(0,0.0d0,fname,n_bond)
!     ---> count the number of bonds
!
     open(iunit,file=fname,position='append',status='old')
       write(iunit,*)n_bond
     close(iunit)
!
     call elses_ana_icohp(0,0.0d0,fname,n_bond)
!     ---> add the bond analysis
   endif  
!
!
end subroutine elses_vis_xyz_extended



