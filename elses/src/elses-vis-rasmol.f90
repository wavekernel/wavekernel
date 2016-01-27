!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
subroutine elses_vis_rasmol(step_count, snapname_header)
!
   use elses_mod_sel_sys,    only : c_system
   use elses_mod_phys_const, only : angst
!  use elses_mod_md_dat,   only : itemd
!  use elses_mod_sim_cell, only : iperiodic, noa, nos, ax, ay, az
   use elses_mod_sim_cell, only : noa, nos, ax, ay, az
   use elses_mod_tx,       only : tx, ty, tz, jsei
   use elses_mod_txp,      only : txp, typ, tzp
   use elses_mod_ctrl,     only : i_verbose
   use elses_mod_file_io
   use elses_mod_vis_main, only : i_cohp, &
&                                 iunit_for_rasmol_script, &
&                                 snap_for_rasmol
   use elses_mod_elem_name, only : elem_name
!
   implicit none
   integer :: iunit, iunit2
   character(len=32) :: fname, c_name
   character(len=32) :: snapname
   character(len=*) :: snapname_header
   integer           :: step_count
   integer           :: js, nss
   real(8)           :: txd, tyd, tzd, ddd
   real(8)           :: d_cri
   integer           :: n_bond
   integer           :: ierr
!
!
   if (i_verbose >= 1) then 
      write(*,*)'@@ LSES_VIS_RASMOL:step_count=',step_count
      write(*,*)'snapname_header=',trim(snapname_header)
   endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Prepare the batch file and setting file
!     at the first time
!  
   if (iunit_for_rasmol_script == 0) then
     iunit=vacant_unit()
     open (iunit,file='rasmol-setting.txt',status='replace') 
       write(iunit,*) 'select all'
       write(iunit,*) 'color green'
       write(iunit,*) 'spacefill 150'
       write(iunit,*) 'wireframe 100'
       write(iunit,*) '#pause'
     close (iunit)
     iunit=vacant_unit()
     iunit_for_rasmol_script=iunit
     if (i_verbose >= 1) then 
       write(*,*)' get vacant unit:iunit_for_rasmol_script'
       write(*,*)'     =', iunit_for_rasmol_script
     endif   
     iunit2=iunit_for_rasmol_script
     open (iunit2,file='rasmol-batch.txt',status='replace') 
        write(iunit2,*) 'echo batch mode'
     close(iunit2)
   endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the file unit
!
   iunit=vacant_unit()
   if (i_verbose >= 1) then 
     write(*,*)' get vacant unit:iunit=',iunit
   endif   
!
   write(snapname, '(i10.10)') step_count
   fname=trim(snapname_header)//trim(snapname)//'.pdb'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Writing the script file
!
   iunit2=iunit_for_rasmol_script
   if (i_verbose >= 1) then 
      write(*,*)'.. writing data into batch file:iunit2=',iunit2
   endif   
   open (iunit2,file='rasmol-batch.txt',position='append',status='unknown') 
      write(iunit2,*) 'zap'
      write(iunit2,*) 'set specular on'
      write(iunit2,*) 'background white'
      write(iunit2,*) 'load pdb '//trim(snapname_header)//trim(snapname)//'.pdb'
      write(iunit2,*) 'echo load '//trim(snapname_header)//trim(snapname)//'.pdb'
      write(iunit2,*) 'script rasmol-setting.txt'
      write(iunit2,*) 'write gif '//trim(snapname_header)//trim(snapname)//'.gif'
   close(iunit2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   if ((nos <= 0) .or. (nos > noa)) then
        write(*,*)'ERROR(elses-vis-rasmol)'
        write(*,*)'  nos=',nos
        stop
     endif   
!
   if (i_verbose >= 1) then 
       write(*,*)' iunit, fname=',iunit, fname
   endif   
!
   open(iunit,file=fname,status='replace')
!
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
!    write(*,*) js,trim(c_name),js,txd,tyd,tzd,ddd,ddd
!     write(iunit, "(I5,1X,A4,I4,3F8.3,2F6.2)") &
!&      js,trim(c_name),js,txd,tyd,tzd,ddd,ddd
     write(iunit, "('ATOM  ',I5,1X,A4,' MOL','  ',I4,'    ',3F8.3,2F6.2)", iostat=ierr) &
&      js,trim(c_name),nss,txd,tyd,tzd,ddd,ddd
     if (ierr /=0) then
       write(*,*)'ERROR(elses_vis_rasmol):write command'
       stop
     endif   
   enddo   
!
   close(iunit)
!
   n_bond=0
!    ----> initialize n_bond
   d_cri=-10.0d0
!    ----> critical energy for bond (in eV)
!  call elses_ana_icohp(1,d_cri,fname,n_bond)
!
   open(iunit,file=fname,position='append',status='old')
     write(iunit, "(a3)")'END'
   close(iunit)
!
end subroutine elses_vis_rasmol


