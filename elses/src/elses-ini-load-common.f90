!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
module M_ini_load_common
!
  use M_io_dst_write_log, only : log_unit !(unchanged)
  private
  public load_xml_config
  public ini_load_common
  public ini_misc_set
!
  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!@@ Read the infomation from the argument and the configuration file.
!
!! Copyright (C) ELSES. 2007-2015 all rights reserved
  subroutine load_xml_config
!
   use M_config
   use M_qm_domain,        only : i_verbose     ! (unchaged)
   use M_qm_domain,        only : c_system, nos, S_is_I ! (CHANGED!!)
   use M_qm_domain,        only : csc_ite_converged ! (CHANGED!!)
   use M_qm_domain,        only : csc_dq_converged  ! (CHANGED!!)
   use elses_mod_sim_cell, only : noa
   use elses_xml_config
   use M_io_dst_system_load,  only : dst_structure_load_from_txt !(routine)
   use M_md_dst,              only : set_dst_final               !(routine)
   use M_sax_data_sync,       only : struc_data_sync             !(routine)
   use M_sax_data_sync_matom, only : struc_data_sync_matom       !(routine)
!
   implicit none
   logical, parameter :: dst_micro_treat_huge_file=.true.
!  logical, parameter :: dst_micro_treat_huge_file=.false.
!
!  integer :: i, noa2, nos_def, imode ! unused
!
!
!   config%option%filename = "config.xml"
!     ----> dummy setting
!
!   call elses_get_arg
!        -----> Get the argument 
!
   if( config%option%filename == "" ) then
      write(*,'(a)') "# Error! : no config file is given"
      write(*,'(a)') "# Usage of this program"
      write(*,'(a)') "#   elses sample/C60/config.xml"
      stop
   end if
!
!  i_verbose=config%option%verbose
!
   if (i_verbose >= 1) then
!    write(*,*) '@@ elses_xml_load_config.. is done'
!    write(*,*)'  Config file name =', trim(config%option%filename)
!    write(*,*)'  The verbose level =',i_verbose
!    write(*,*) '...call config_load'
     if (log_unit > 0) then
       write(log_unit,*) '@@ elses_xml_load_config.. is done'
       write(log_unit,*)'  Config file name =', trim(config%option%filename)
       write(log_unit,*)'  The verbose level =',i_verbose
       write(log_unit,*) '...call config_load'
     endif
   endif  
!
   call config_load( config, config%option%filename )
!
!  stop 'Stop manually (before data sync)'
!
   if (config%system%structure%use_vatom) then
     call struc_data_sync  
!   ----> Syncronize the atomic data amond nodes, if needed.
   endif  
!
   if (config%system%structure%use_matom) then
     call struc_data_sync_matom  
!   ----> Syncronize the atomic data amond nodes, if needed.
   endif  
!
!  stop 'Stop manually (after data sync)'
!
   if (i_verbose >= 1) then
!    write(*,*) '...end config_load'
     if (log_unit > 0) then
       write(log_unit,*) '...end config_load'
     endif
   endif  
!
   if (dst_micro_treat_huge_file) then
     if (config%system%structure%natom == 0) then
       call dst_structure_load_from_txt
     endif
   endif
!
   if (i_verbose >= 1) then
     if (log_unit > 0) then
       write(log_unit,*) '...ends: dst_structure_load_from_txt'
     endif
   endif
!
   noa=config%system%structure%natom
!
   if ((noa <= 0) .or. (dble(noa) >= 0.9d0*dble(huge(0)))) then
     write(*,*) 'ERROR?(load_xml_structure:# atom is zero or too large !?)'
     write(*,*) '  noa=',noa
     stop
   endif   
!
   nos=config%system%structure%nelement
!
   if ((nos <= 0) .or. (nos >=100)) then
     write(*,*) 'ERROR(elses_xml_load_config):nos=',nos
     write(*,*) '#ELSES: Stop by error: This may be due to mismatch'
     write(*,*) '#  between the configuration XML file'
     write(*,*) '#  and the structure XML file'
     write(*,*) '# The element tag is required for each atom species'
     write(*,*) '#  in the configuration file.'
     stop
   endif   
!
   c_system=trim(config%system%structure%velement(1)%quantum%type)
!
   if (c_system == 'Xu.1992') c_system='C_Xu'
   if (c_system == 'Kwon.1994') c_system='Si_Kwon'
   if (c_system == 'NRL') c_system='NRL_leg'
!
   S_is_I = .true.
   if (c_system == 'geno')     S_is_I = .false.
   if (c_system == 'NRL_leg')  S_is_I = .false.
!
   csc_ite_converged = -1     ! dummy value
   csc_dq_converged  = -1.0d0 ! dummy value
!
   if (i_verbose >= 1) then
!     write(*,*) '  nos   =',nos
!     write(*,*) '  noa   =',noa
!     write(*,*) 'c_system=',c_system
!     write(*,*) 'S_is_I  =',S_is_I
      if (log_unit > 0) then
        write(log_unit,*) '  nos   =',nos
        write(log_unit,*) '  noa   =',noa
        write(log_unit,*) 'c_system=',c_system
        write(log_unit,*) 'S_is_I  =',S_is_I
      endif
   endif   
!
   if (trim(config%calc%solver%scheme) == 'scheme_default') then
       config%calc%solver%scheme = 'eigen'
!
!    if ((c_system == 'geno') .or. (c_system == 'NRL_leg')) then
!      config%calc%solver%scheme = 'eigen'
!    else
!      config%calc%solver%scheme = 'krylov'
!
     if (i_verbose >= 1) then
!      write(*,*) '  solver_scheme =',trim(config%calc%solver%scheme)
       if (log_unit > 0) write(log_unit,*) '  solver_scheme =',trim(config%calc%solver%scheme)
     endif   
!
   endif   
!
  end subroutine load_xml_config
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!@@ Initial load
!
!! Copyright (C) ELSES. 2007-2015 all rights reserved
  subroutine ini_load_common
!
!   Module variables used: unchanged 
!     -->   i_verbose
!   Module variables used: allcoated (in slaved routines)
!       ---> nval(nos), elem_name(nos),
!            awt(nos), val_elec_atm(nos)
!             tx(noa),   ty(noa),   tz(noa)
!            txp(noa),  typ(noa),  tzp(noa)
!           velx(noa), vely(noa), velz(noa)
!          foi(noa,3), foiold(noa,3),
!          iflag(noa), amm(noa)
!
     use M_config, only: config !(unchanged)(only config%calc%mode)
     use elses_mod_ctrl,     only : i_verbose               !(unchanged)
     use M_md_motion,        only : elses_md_motion_ini_set !(routine)
     use M_md_velocity_routines, only : allocate_velocity   !(routine)
!
     implicit none
!
     if (i_verbose >=1) then
       if (log_unit > 0) then
         write(log_unit,'(a)')'@@ ini_load_common_alloc'
       endif
     endif   
!
     call load_element_alloc
!       ---> Allocate   : nval(nos), elem_name(nos),
!                         awt(nos), val_elec_atm(nos)
!
     if (config%system%structure%use_vatom) then
       call load_structure_alloc
     else  
       if (i_verbose >=1) then
        if (log_unit > 0) then
           write(log_unit,*)'INFO:load_structure_alloc is skipped:use_vatom=', config%system%structure%use_vatom
        endif
       endif   
     endif  
!       ---> Allocate   :
!               tx(noa),   ty(noa),   tz(noa)
!              txp(noa),  typ(noa),  tzp(noa)
!            foi(noa,3), foiold(noa,3),
!            iflag(noa), amm(noa)
!
     if (config%calc%mode == "matrix_generation") then
       config%calc%solver%scheme="krylov"
       if (i_verbose >=1) then
        if (log_unit > 0) then
           write(log_unit,*)'INFO:matrix generation mode: solver mode is set to ', trim(config%calc%solver%scheme)
        endif
       endif   
     endif   
!
     if ((config%calc%mode == "dynamics") .or. (config%calc%mode == "conversion")) then
       call allocate_velocity 
!       ---> Allocate : velx(noa), vely(noa), velz(noa)
     endif   
!
     call set_unitcell_info
!       ---> Set values :ax, ay, az, i_pbc_x, i_pbc_y, i_pbc_z
!
     call elses_md_motion_ini_set
!       ---> Set values : dtmd, itemdmx
!
     call set_structure_data
!       ---> Set values:
!           itemdorg, 
!           r_hb_per_atom, amq, thb, vhb, 
!           txp(noa), typ(noa), tzp(noa), jsei(noa), iflag(noa)
!           tx(noa), ty(noa), tz(noa)
!           velx(noa), vely(noa), velz(noa),
!           foi(noa,3), 
!           ntlimit, r_mem_liit,
!           noak_min_def, nreclc_def
!
     call load_xml_conditions
!       ---> Set values:
!           tempk0, temp_for_electron,
!           r_mem_limit, ntlimit
!
     if (i_verbose >=1) then
       if (log_unit > 0) then
         write(log_unit,*)'.. ended : ini_load_common_alloc'
       endif
     endif   
!
  end subroutine ini_load_common
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine load_element_alloc
!
!
!   Module variables used: unchanged 
!     -->   i_verbose, c_system, nos
!   Module variables used: changed (in slaved routines)
!     -->   nval(nos), elem_name(nos), awt(nos)
!           val_elec_atm(nos)
!
    use M_qm_domain,        only : i_verbose, c_system, nos
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!
     implicit none
     integer imode
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!  Verbose message and trivial checkings (non-essential)

     if (i_verbose >= 1) then
       if (log_unit > 0) then
         write(log_unit,*) '@@ elses_xml_load_element_geno'
         write(log_unit,*) '      nos =',nos
         write(log_unit,*) ' c_system =',c_system
       endif 
     endif  
!
     if ((nos <= 0) .or. (nos >=100)) then
       write(*,*) 'ERROR(elses_xml_load_element_geno)'
       write(*,*) '   nos=',nos
       write(*,*) 'nos should be already set properly'
       stop
     endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!  Allocation
!
     imode=0
     call elses_alloc_nval(nos,imode)
!       ----> alocate nval(nos)
!
     imode=0
     call elses_alloc_elem_name(nos,imode)
!       ----> alocate elem_name(nos)
!
     call elses_alloc_awt(nos)
!       ----> alocate awt(nos)
!
     call elses_alloc_val_elec(nos)
!       ----> alocate call val_elec_atm(nos)
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!
!
  end subroutine load_element_alloc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine load_structure_alloc
!
!
!   Module variables used: unchanged 
!     -->   i_verbose, noa
!   Module variables used: allocated (in slaved routines)
!     -->     tx(noa),   ty(noa),   tz(noa)
!     -->    txp(noa),  typ(noa),  tzp(noa)
!     -->   velx(noa), vely(noa), velz(noa)
!     -->    foi(noa,3)
!     -->    foiold(noa,3)
!     -->    iflag(noa), amm(noa)
!
     use elses_mod_ctrl,     only : i_verbose 
     use elses_mod_sim_cell, only : noa
     use M_md_alloc,         only : array_alloc_for_md !(routine)
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!
     implicit none
     integer noa2
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!  Verbose message and trivial checkings (non-essential)

     if (i_verbose >= 1) then
       if (log_unit > 0) then
         write(log_unit,*)'@@ load_structure_alloc'
         write(log_unit,*)'      noa =',noa
       endif
     endif  
!
     if ((noa <= 0) .or. (noa >= huge(0))) then
       write(*,*) 'ERROR(load_structure_alloc)'
       write(*,*) '   noa=',noa
       write(*,*) '# atoms is zero or too large ?!'
       stop
     endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!  Allocation
!
     noa2=noa
     call array_alloc_for_md(noa2)
!    call elses_alloc_ini_txp(noa2)
!      -->allocation : tx(noa),   ty(noa),   tz(noa)
!                     txp(noa),  typ(noa),  tzp(noa)
!                    velx(noa), vely(noa), velz(noa)
!                   foi(noa,3),
!                foiold(noa,3),
!                   iflag(noa)
!
     call elses_alloc_amm(noa2)
!      -->allocation : amm(noa)
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!
!
  end subroutine load_structure_alloc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set strcuture data
!
  !! Copyright (C) ELSES. 2007-2015 all rights reserved
  subroutine set_structure_data
    use M_config
    use elses_mod_ctrl,     only : i_verbose           !(unchanged)
    use elses_mod_tx,       only : tx, ty, tz, jsei    !(CHANGED)
    use elses_mod_txp,      only : txp, typ, tzp       !(CHANGED)
    use elses_mod_md_dat,   only : dtmd                !(unchagend)   
    use elses_mod_sim_cell, only : noa,ax,ay,az, i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
    use elses_mod_vel,      only : velx, vely, velz    !(CHANGED) 
    use elses_mod_foi,      only : foi                 !(CHANGED)
    use elses_mod_md_dat,   only : itemdorg            !(CHANGED)
    use elses_mod_thermo,   only : amq, thb,thbold,vhb,vhbold !(CHANGED)
    use elses_mod_iflag,    only : iflag                      !(CHANGED)
    use elses_param_ctl_kr, only : nreclc_def, noak_min_def  !(CHANGED) 
    use M_md_motion,        only : elses_md_motion_ini_set           !(routine)
    use M_ini_load_vatom,   only : set_structure_data_from_vatom     !(routine)
    use M_ini_load_vatom,   only : set_atm_position_from_tx          !(routine)
    use M_ini_load_matom,   only : set_structure_data_from_matom     !(routine)
!
    implicit none
    integer :: j, k
    type(atom_type), pointer :: atom
    real(8) :: r_mem_limit
    real(8) :: ntlimit
    integer :: noak_min_tmp, nrecl_tmp
    real(8) :: r_hb_mass_per_atom
    real(8) :: dtmd_wrk
!      
    integer :: log_unit
!    
    if (i_verbose >=1 ) then
       write(6,*) '@@ set_structure_data'
    endif  
!
    log_unit=config%calc%distributed%log_unit
!
    itemdorg = config%system%structure%mdstep
!
    r_hb_mass_per_atom= config%system%structure%heatbath%massperatom
!   write(*,*)' r_hb_mass_per_atom=',r_hb_mass_per_atom
!
    amq=1.0d0/r_hb_mass_per_atom/dble(noa)
!
    thb = config%system%structure%heatbath%position
    vhb = config%system%structure%heatbath%velocity
!
    if (i_verbose >=1 ) then
      if (log_unit > 0 ) then
        write(log_unit,*) ' i_pbc_x  =',i_pbc_x
        write(log_unit,*) ' i_pbc_y  =',i_pbc_y
        write(log_unit,*) ' i_pbc_z  =',i_pbc_z
        write(log_unit,*) ' thb, vhb =', thb, vhb
        write(log_unit,*) ' r_hb_mass_per_atom  =', r_hb_mass_per_atom
       endif 
    endif  
!
    dtmd_wrk=dtmd
    if (dtmd_wrk < 1.0d-10) then
      if (i_verbose >=1 ) then
        if (log_unit > 0 ) then
          write(log_unit,*) 'INFO: No meaningful value of dtmd is found'
          write(log_unit,*) 'INFO: The velocity data in the input file will be ignored, if exist'
        endif  
      endif  
      dtmd_wrk=0.0d0
    endif
!   
    if (config%system%structure%use_matom) then
      call set_structure_data_from_matom
    endif  
!
    if (config%system%structure%use_vatom) then
      call set_structure_data_from_vatom
      call set_atm_position_from_tx
    endif  
!
    ntlimit = config%calc%limit%time ! sec
    r_mem_limit = config%calc%limit%memory  ! byte
!
!   call elses_md_motion_ini_set
!         -----> set itemdmx and dtmd 
!
    noak_min_tmp = config%calc%solver%projection
    nrecl_tmp    = config%calc%solver%dimension
!
    if (i_verbose >= 1 ) then
      if (log_unit > 0 ) then
        write(log_unit,*)'solver scheme=',config%calc%solver%scheme
        write(log_unit,*)'noak_min_tmp =',noak_min_tmp
        write(log_unit,*)'nrecl_tmp    =',nrecl_tmp
      endif
    endif  
!
    noak_min_def=noak_min_tmp
    nreclc_def=nrecl_tmp
!
    if (nreclc_def <= 0) then
      write(6,*)'ERROR:nreclc_def=',nreclc_def
      write(6,*)'   .. at ELSES_XML_LOAD_STR_EXP'
      stop
    endif   
!
  end subroutine set_structure_data
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set the unit cell information 
!      -->   ax, ay, az, i_pbc_x, i_pbc_y, i_pbc_z
!
  subroutine set_unitcell_info
    use M_qm_domain,        only : i_verbose !(unchanged)
    use M_qm_domain, only : ax,ay,az, &               !(CHANGED)
&                           i_pbc_x, i_pbc_y, i_pbc_z !(CHANGED)
    use M_config !(unchanged)
    use M_lib_math_func, only : get_determinant !(routines)
!
    implicit none
    real(8) :: cell_vectors(3,3)
    real(8) :: ddd, det
    integer :: i,j
!
    if (i_verbose >= 1) then
      write(*,*)'@@ set_unitcell_info'
    endif   
!
!  @ Define the cell vectors 
!
    cell_vectors(:,1)= config%system%structure%unitcell%vectorA(:) ! in a.u.
    cell_vectors(:,2)= config%system%structure%unitcell%vectorB(:) ! in a.u.
    cell_vectors(:,3)= config%system%structure%unitcell%vectorC(:) ! in a.u.
!
    if (i_verbose >= 1) then
      write(*,*)'Info : The unit cell vectors in the input file are as follows'
      write(*,'(a,3F20.10,a)')'  vector a=(',cell_vectors(:,1),')'
      write(*,'(a,3F20.10,a)')'  vector b=(',cell_vectors(:,2),')'
      write(*,'(a,3F20.10,a)')'  vector c=(',cell_vectors(:,3),')'
    endif  
!  
!  @ Check the orthogonality
!
    do j=1,3
      do i=1,3
        if (i == j) cycle
        ddd=dot_product(cell_vectors(:,i),cell_vectors(:,j))
        if (i_verbose >= 1) then
           write(*,*)'deviation from the orthogonality (i,j)',i,j,dabs(ddd)
        endif   
        if (dabs(ddd) >= 1.0d-10) then
          write(*,*)'ERROR(set_unitcell_info): Unsupported unit cell vectors'
          stop
        endif   
      enddo
    enddo  
!
!  @ Check the determinant
!
    call get_determinant(cell_vectors,det)
    if (i_verbose >= 1) then
      write(*,*)' det=',det
    endif  
    if ( det < 1.0d-10 ) then
      write(*,*)'ERROR(set_unitcell_info): Unsupported unit cell vectors'
      write(*,*)' det=',det
      stop
    endif   
!
!  @ Set (ax, ay, az) 
    ax = dsqrt(dot_product(cell_vectors(:,1),cell_vectors(:,1)))
    ay = dsqrt(dot_product(cell_vectors(:,2),cell_vectors(:,2)))
    az = dsqrt(dot_product(cell_vectors(:,3),cell_vectors(:,3)))
    config%system%structure%unitcell%vectorA(:) = 0
    config%system%structure%unitcell%vectorB(:) = 0
    config%system%structure%unitcell%vectorC(:) = 0
    config%system%structure%unitcell%vectorA(1) = ax
    config%system%structure%unitcell%vectorB(2) = ay
    config%system%structure%unitcell%vectorC(3) = az
    cell_vectors(:,1)= config%system%structure%unitcell%vectorA(:) ! in a.u.
    cell_vectors(:,2)= config%system%structure%unitcell%vectorB(:) ! in a.u.
    cell_vectors(:,3)= config%system%structure%unitcell%vectorC(:) ! in a.u.
!
    if (i_verbose >= 1) then
      write(*,*)'Info : The unit cell vectors are redefined as follows'
      write(*,'(a,3F20.10,a)')'  vector a=(',cell_vectors(:,1),')'
      write(*,'(a,3F20.10,a)')'  vector b=(',cell_vectors(:,2),')'
      write(*,'(a,3F20.10,a)')'  vector c=(',cell_vectors(:,3),')'
    endif  

!    
!  @ Check the volume
!
    if (i_verbose >= 1) then
      write(*,*)' ax*ay*ax =',ax*ay*az
    endif  
!
    if ( det < 1.0d-10 ) then
      write(*,*)'ERROR(set_unitcell_info): Unsupported unit cell vectors'
      write(*,*)'     det  =',det
      stop
    endif   
!
    ddd=dabs(ax*ay*az/det)-1.0d0
    if ( ddd > 1.0d-10 ) then
      write(*,*)'ERROR(set_unitcell_info): Unsupported unit cell vectors'
      write(*,*)' ax*ay*ax =',ax*ay*az
      write(*,*)'     det  =',det
      stop
    endif   
!
!  @ Set (i_pbc_x, i_pbc_y, i_pbc_z)
!
    i_pbc_x=0
    i_pbc_y=0
    i_pbc_z=0
    if( config%system%boundary%periodic_x ) i_pbc_x=1
    if( config%system%boundary%periodic_y ) i_pbc_y=1
    if( config%system%boundary%periodic_z ) i_pbc_z=1
!
    write(*,*)'ax =',ax
    write(*,*)'ay =',ay
    write(*,*)'az =',az
!   stop
!
  end subroutine set_unitcell_info
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Load various conditions from XML file
!
  subroutine load_xml_conditions
!
!   Module variables used: unchanged 
!     -->   i_verbose, noa
!   Module variables used: changed (in slaved routines)
!     -->    tempk0, temp_for_electron,
!            r_mem_limit, ntlimit
!
     use elses_mod_ctrl,     only : i_verbose 
     use elses_mod_sim_cell, only : noa
!
!cccccccccccccccccccccccccccccccccccccccccccccc
!
     implicit none
!
     if (i_verbose >= 1 ) then
       write(*,*)'@@ load_xml_conditions'
     endif   
!
!    call elses_xml_load_temperature
!      ---> set atomic temperature : tempk0
!
     call elses_xml_load_elec_temperature
!      ---> set electronic temperature : temp_for_electron
!
     call elses_xml_load_ctrl
!      ---> set values : r_mem_limit, ntlimit
!
!    call elses_xml_vis
!      ---> sepecify the visualization tool
!
  end subroutine load_xml_conditions
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Initial misc settings 
!     including ones for compatibility for older code
!
!   Module variables used: unchanged 
!     -->   i_verbose, c_system
!   Module variables used: changed (in slaved routines)
!     --> amm, j2js, j2ja, js2j, n_tot_base
!         ict4h, ict4l 
!
  subroutine ini_misc_set
!
     use M_config                 !(unchanged)
     use M_qm_domain,        only : i_verbose, c_system
!
     if (i_verbose >= 1 ) then
       write(*,*)'@@ ini_set_compat'
     endif   
!
     if (.not. config%system%structure%use_matom) then
       call elses_set_mass2
!         ---> Set values : amm(noa)
     endif  
!
     if (c_system == 'geno') then
       call ini_misc_set_geno
!        --> Set j2js, j2ja, js2j, n_tot_base, ict4h, ict4l
     else
       call elses_init_compat2
!        --> Set j2js, j2ja, js2j, n_tot_base, ict4h, ict4l
!                intf, dbx, dby, dbz, idngl, r_base
     endif
!
  end subroutine ini_misc_set
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Initial misc settings for geno
!
!  Module variables unchanged
!    -->  noa, i_verbose, c_system, nval, jsei
!  Module variables changed
!    -->  j2js, j2ja, js2j, n_tot_base, ict4h, ict4l
!
  subroutine ini_misc_set_geno
!
    use M_config
!
    use elses_mod_sim_cell, only : noa
    use M_qm_domain,        only : i_verbose, c_system, nval
    use M_qm_domain,        only : atm_element
!   use elses_mod_tx,       only : jsei
    use elses_mod_orb2,     only : j2js,j2ja,js2j,n_tot_base
    use elses_mod_multi,    only : ict4h, ict4l
!
    implicit none
    integer js,ja,iii,nss,nval2,j,ierr
    integer nvl_tmp
!
    if (i_verbose >= 1) then
      write(*,*)'@@ init_compat_geno'
    endif  
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Error checking
!
    if (c_system /= 'geno') then
      write(6,*)'ERROR:elses_init_compat_geno'
      write(6,*)'c_system=',c_system
      stop
    endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Count up the total number of bases 
!
    iii=0
    do js=1,noa
!     nss=jsei(js)
      nss=atm_element(js)
      if ((nss <= 0) .or. (nss > 100)) then
        write(*,*)'ERROR(ini_misc_set_geno):nss=',nss
        stop
      endif   
      nval2=nval(nss)
      if ((nval2 <= 0) .or. (nval2 > 10)) then
        write(*,*)'ERROR(ini_misc_set_geno):nval2=',nss
        stop
      endif   
      iii=iii+nval2
      if (dble(iii) > 0.9d0*dble(huge(0))) then
        write(*,*)'ERROR(ini_misc_set_geno):iii=',iii
        stop
      endif   
    enddo   
    n_tot_base=iii
!
    if (dble(n_tot_base) > 0.9d0*dble(huge(0))) then
      write(*,*)'ERROR(ini_misc_set_geno):n_tot_base=',n_tot_base
      stop
    endif   
!
    if ( (n_tot_base  <= 0) .or. (dble(n_tot_base) >= 1000.0d0*dble(noa)) ) then
      write(*,*)'ERROR:n_tot_base=',n_tot_base
      stop
    endif
!
    if (i_verbose >= 1) then
      write(*,*)' total number of bases=',n_tot_base
    endif  
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Set nvl_tmp as the maximum orbital number
!
    nvl_tmp=maxval(nval)
!
    if ((nvl_tmp <= 0) .or. (nvl_tmp >= 100)) then
      write(6,*)'ERROR:nvl_tmp=',nvl_tmp
      stop
    endif
!   
   if (config%system%structure%use_matom) then
     if (log_unit > 0) then
       write(log_unit,*) 'INFO:(use_matom):j2js is not allocated'
       write(log_unit,*) 'INFO:(use_matom):j2ja is not allocated'
       write(log_unit,*) 'INFO:(use_matom):js2j is not allocated'
     endif  
     return
   endif
!   
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Allocate arrays: j2js,j2ja,js2j
!
    allocate (j2js(n_tot_base),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'allocation error!(J2JS):ierr=',ierr
      stop
    endif
!
    allocate (j2ja(n_tot_base),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'allocation error!(J2JA):ierr=',ierr
      stop
    endif
!
!   allocate (js2j(nvl_tmp,n_tot_base),stat=ierr)
    allocate (js2j(nvl_tmp,noa),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'allocation error!(JS2J):ierr=',ierr
      stop
    endif
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Setting ICT4H, ICT4L 
!
    ict4h=1
    ict4l=1
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Setting j2js, j2ja, js2j
!
    j2js(:)=0
    j2ja(:)=0
    js2j(:,:)=0
    j=0
    do js=1,noa
!     nss=jsei(js)
      nss=atm_element(js)
      if ((nss <= 0) .or. (nss > 100)) then
        write(*,*)'ERROR(ini_misc_set_geno):nss=',nss
        stop
      endif   
      nval2=nval(nss)
      do ja=1,nval2
        j=j+1 
        j2js(j)=js
        j2ja(j)=ja
        js2j(ja,js)=j
      enddo
    enddo
!
    if (i_verbose >= 1) then
       write(*,*)'...ended: ini_misc_set_geno'
    endif   
!
  end subroutine ini_misc_set_geno
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module M_ini_load_common

