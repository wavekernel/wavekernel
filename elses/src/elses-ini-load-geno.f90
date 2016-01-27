!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_ini_load_geno
!
  use M_qm_domain,      only : i_verbose
!
  private
  public ini_load_geno
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine ini_load_geno
!
!   Module procedure
!
    use M_config,               only : config !(unchanged)
    use M_ini_load_common, &
         only:  load_xml_config, ini_load_common, ini_misc_set
    use M_md_cell_change, only : elses_md_cell_change
    use M_lib_element_database, only : make_db_element_name, show_db_element_name
    use M_qm_flexible_cutoff,   only : flexible_cutoff_ini !(routine)
    use M_md_output,            only : setting_for_main_output !(routine)
    use M_xml_compat_chk,       only : check_xml_compat !(routine)
    use M_lib_quad_prec,        only : quad_prec_test_for_pi !(routine)
    use M_lib_phys_const,       only : test_convert_unit     !(routine)
    use M_check_xml_element_info, only : check_xml_element_info !(routine)
    use M_group_id_setting,       only : ini_group_id           !(routine)
    use M_lib_read_w_split,       only : test_read_w_split      !(routine)

    ! 
    implicit none
    logical, parameter :: cell_change_for_next_step = .false.
    logical, parameter :: show_db_element_name_is_called = .false.
    integer            :: log_unit
    logical            :: quad_prec_test_verbose
    logical            :: call_test_read_w_split
!
    logical, parameter :: test_convert_unit_is_called = .true.
!
    log_unit=config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      write(*,*)'@@ ini_load_geno'
    endif
!   
    call make_db_element_name
!       ---> Make element-name database
!
    if (show_db_element_name_is_called) then
      call show_db_element_name
    endif  
!       ---> Show element-name database (optional)
!
    if (test_convert_unit_is_called) then
      call test_convert_unit(log_unit)
    endif  
!       ---> Test 'convert unit' routine (optional)
!
    if (log_unit > 0) then
      if (i_verbose >= 1) then
        quad_prec_test_verbose = .true.
        call_test_read_w_split      = .true.
      else
        quad_prec_test_verbose = .false.  
        call_test_read_w_split      = .false.
      endif   
    endif  
!
    if (call_test_read_w_split) then
      call test_read_w_split(log_unit)
    endif
!
    if (quad_prec_test_verbose) then
      call quad_prec_test_for_pi(log_unit, quad_prec_test_verbose)
      write(log_unit,'(a)') 'INFO:Quad precision test ... OK'
    endif  
!
    call load_xml_config
!       ---> Set values : i_verbose, c_system, nos
!
    call check_xml_compat
!       ---> Check the XML file setting for compatibility
!
    call setting_for_main_output
!       ---> Set the main output file
!
    call ini_load_common
!      ---> Allocate variables and load various information 
!             in config. and structure XML files
!            Ex. atom position, velocity, cell lenghts 
!
    call load_element_info
!      ---> Set value  : nval(nos), elem_name(nos), 
!                        awt(nos), val_elec_atm(nos), r_cut_book
!                      : parameters for Hamiltonian
!
    call check_xml_element_info
!      ---> Check the element info
!
    call set_tot_elec_num
!      ---> Set value  : val_elec_tot
!                       (total electron number)
!
    call ini_misc_set
!       ---> Initial misc. setting 
!       ---> Set values (mainly):
!                j2js, j2ja, js2j, n_tot_base
!
    call ini_group_id
!       ---> Initial setting for group id
!
    call for_conversion_mode
!       ---> Routines for conversion mode
!             This routine is meaningful, only if calc%mode='conversion'
!
    call elses_md_cell_change(cell_change_for_next_step)
!       ---> Optional routines for possible cell change
!
!   call flexible_cutoff_ini
!       ---> Optional initial routine for the flexible cutoff
!
    call chk_input_data
!       ---> Check the consistency of input data
!            No module variable is modified
!
  end subroutine ini_load_geno
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Routines for conversion mode
!      This routine is meaningful, only if calc%mode='conversion'
! 
  subroutine for_conversion_mode
    use M_config,           only: config    !(unchanged)
    use elses_mod_ctrl,     only: i_verbose !(unchanged)
    use elses_mod_md_dat,   only: itemd     !(CHANGED)
    use M_md_save_struct,   only: elses_md_save_struct !(routine)
    use M_md_dst,           only: set_dst_final !(routines)
    implicit none
    logical :: final_iteration
!
!
    if (i_verbose >= 1) then
      write(*,*)'@@ for_conversion_mode'
    endif
!
    if ( config%calc%mode /= 'conversion' ) then
      if (i_verbose >= 1) write(*,*)' ... is ignored'
      return
    endif   
!
    write(*,*)'The file conversion mode is started'
!
    itemd=1
    final_iteration=.true.
    config%output%position%interval=1
    call elses_md_save_struct(final_iteration)
!
    call set_dst_final
!
    write(*,*)'The file conversion mode is ended successfully'
    stop
!
  end subroutine for_conversion_mode
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set values for following element information
!           nvl, nval(nos), elem_name(nos), awt(nos)
!           val_elec_atm(nos), r_cut_book
! 
  subroutine load_element_info
!
!   module variables: unchanged 
!     -->   i_verbose, c_system, nos, au_mass, angst
!
!   module variables: changed 
!     -->   nvl, nval(nos), elem_name(nos), awt(nos)
!           val_elec_atm(nos), r_cut_book
!          ( They are already allocated )
!
   use elses_mod_phys_const, only : au_mass, angst
   use M_config  ! CHANGED ONLY IN config%system%cutoff_radius
                 !                 config%calc%interaction_range%cutoff_rest
                 !                 config%calc%interaction_range%cutoff_rest_non_vdW
   use M_qm_domain,      only : c_system, nval, nos, &
        i_pbc_x, i_pbc_y, i_pbc_z, ax, ay, az
   use elses_mod_sel_sys,    only : r_cut_book
!
   use elses_mod_orb1,         only : nvl
   use elses_mod_elem_name,    only : elem_name
   use elses_mod_mass,         only : awt
   use elses_mod_val_elec,     only : val_elec_atm
!
   use elses_xml_element
!
   use M_qm_geno_Huckel_atom_params, only : setAtomParameters
   implicit none
   real(8) ddd
   integer ierr, nss
   character(len=32) :: c_name
   integer           :: log_unit, lu
   real(8) :: rcut_cell_max
!  real(8) :: work_value
!
!
   log_unit = config%calc%distributed%log_unit
   lu       = config%calc%distributed%log_unit
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
! Message and trivial checking
!
   if (i_verbose >= 1) then
      if (lu > 0) write(lu,*)'@@ load_element_info:nos=',nos
   endif   
!
   if ((nos <= 0) .or. (nos >= 100)) then
     write(*,*)'ERROR(load_element_info)'
     write(*,*)'nos=',nos
     stop
   endif   
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
! Set elem_name(nos)
!
   do nss=1,nos
     elem_name(nss)=trim(config%system%structure%velement(nss)%name)
     if (i_verbose >= 1) then
       if (lu > 0) write(lu,*)'nss, elem_name=',nss,elem_name(nss)
     endif   
   enddo
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
! Set rcut_cell_max as possible maximum of interaction radius
!
   rcut_cell_max=huge(1.0d0)
   if( i_pbc_x == 1 ) rcut_cell_max = min(rcut_cell_max,ax/2.01d0)
   if( i_pbc_y == 1 ) rcut_cell_max = min(rcut_cell_max,ay/2.01d0)
   if( i_pbc_z == 1 ) rcut_cell_max = min(rcut_cell_max,az/2.01d0)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
! Set values for each case
!
   ierr=1
!
   if (c_system /= 'geno') then
     if (trim(config%calc%solver%scheme) == 'eigen') then
       if ( (c_system == 'C_Xu') .or. (c_system == 'Si_Kwon') ) then
         if (i_verbose >= 1) then
           if (lu > 0) write(lu,*)'call setAtomParameters for the compatibility to GENO'
         endif   
         call setAtomParameters(NOS,elem_name,awt,nval,val_elec_atm)
       endif    
     endif
   else
       if (i_verbose >= 1) then
         if (lu > 0) write(lu,*)'call setAtomParameters for GENO'
       endif   
      call setAtomParameters(NOS,elem_name,awt,nval,val_elec_atm)
      r_cut_book = config%system%cutoff_radius
      ddd=r_cut_book
      if( i_pbc_x == 1 ) then 
        if (ax/2.01d0 < r_cut_book) then ! 0.01 is a tolerance factor
          r_cut_book=min(r_cut_book,ax/2.01d0)
!         write(*,*)'Inconsistent parameter: ax, r_cut_book=', ax, r_cut_book
!         write(*,*)'Tentatively, AX should be larger than 2 x (r_cut_book) in PBC'
!         stop
        endif   
      endif   
      if( i_pbc_y == 1 ) then 
        if (ay/2.01d0 < r_cut_book) then 
          r_cut_book=min(r_cut_book,ay/2.01d0)
!         write(*,*)'Inconsistent parameter: ay, r_cut_book=', ay, r_cut_book
!         write(*,*)'Tentatively, AY should be larger than 2 x (r_cut_book) in PBC'
!         stop
        endif   
      endif   
      if( i_pbc_z == 1 ) then 
        if (az/2.01d0 < r_cut_book) then 
          r_cut_book=min(r_cut_book,az/2.01d0)
!         write(*,*)'Inconsistent parameter: az, r_cut_book=', az, r_cut_book
!         write(*,*)'Tentatively, AZ should be larger than 2 x (r_cut_book) in PBC'
!         stop
        endif   
      endif
      if (lu > 0) then
        write(lu,'("load_element_info: Are boundaries periodic? ",3L1)') i_pbc_x == 1, i_pbc_y == 1, i_pbc_z == 1
        write(lu,'("load_element_info: r_cut_book for geno =",ES21.14)') r_cut_book
      endif  
!
      config%system%cutoff_radius=r_cut_book
!      
      If (dabs(r_cut_book-ddd) .gt. 1.0d-12) then
!       write(*,*)'Inconsistent parameter: r_cut_book=', r_cut_book
!       stop
        if (lu > 0) then
          write(lu,*)'INFO:r_cut_book is modified by the sizes of the periodic simulation cell'
          write(lu,*)'INFO-XML:changed:config%system%cutoff_radius '
        endif  
      endif   
!
      if (r_cut_book .lt. 2.0d0/angst) then
        write(*,*)' ABORT: r_cut_book seems to be too short ( < 2 A )'
        stop
      endif   
!
      if (config%calc%interaction_range%cutoff_rest < 0.0d0) then
        if (config%calc%interaction_range%cutoff_rest_cellmax) then
         config%calc%interaction_range%cutoff_rest = rcut_cell_max
        else
         config%calc%interaction_range%cutoff_rest = huge(1.0d0)
        endif
        if (lu > 0) then
          write(lu,*)'INFO:interaction_range%cutoff_rest         =', & 
&                          config%calc%interaction_range%cutoff_rest
        endif
      endif
!
      if (config%calc%interaction_range%cutoff_rest_non_vdW < 0.0d0 ) then
        config%calc%interaction_range%cutoff_rest_non_vdW  = config%calc%interaction_range%cutoff_rest
        if (lu > 0) then
          write(lu,*)'INFO:interaction_range%cutoff_rest_non_vdW =', & 
&                          config%calc%interaction_range%cutoff_rest_non_vdw
        endif
      endif
!  
      ierr=0
      do nss=1,nos
         if (awt(nss) < 1.0d-10) then
            write(*,*)'ERROR:load_element_info'
            write(*,*)'awt(nss)=',awt(nss)
            stop
         endif
!         nval(nss)=9
!         val_elec_atm(nss)=config%system%structure%velement(nss)%classic%charge
      enddo

!!$     do nss=1,nos
!!$       ierr=0
!!$       call elses_nrl_ini(ddd)
!!$!        ---> awt(nss) is set in elses_nrl_ini
!!$       r_cut_book=ddd
!!$       if (awt(nss) < 1.0d-10) then
!!$         write(*,*)'ERROR:load_element_info'
!!$         write(*,*)'awt(nss)=',awt(nss)
!!$         stop
!!$       endif
!!$       nval(nss)=9
!!$       val_elec_atm(nss)=config%system%structure%velement(nss)%classic%charge
!!$     enddo  
   endif  
!
   if (c_system == 'C_Xu') then
     nss=1
     ierr=0
     nval(nss)=4
     c_name=elem_name(nss)
     if (c_name /= 'C') then
       write(*,*)'ERROR:load_element_info'
       write(*,*)'elem_name(nss)=',elem_name(nss)
       stop
     endif
     awt(nss)=12.01d0
     val_elec_atm(nss)=4.0d0
     r_cut_book=2.6d0/angst
   endif  
!
   if (c_system == 'Si_Kwon') then
     nss=1
     ierr=0
     nval(nss)=4
     c_name=elem_name(nss)
     if (c_name /= 'Si') then
       write(*,*)'ERROR:load_element_info'
       write(*,*)'elem_name(nss)=',elem_name(nss)
       stop
     endif
     awt(nss)=28.0855d0
     val_elec_atm(nss)=4.0d0
     r_cut_book=4.16d0/angst
   endif  
!
   if (c_system == 'NRL_leg') then
     nss=1
     ierr=0
     call elses_nrl_ini(ddd)
!        ---> awt(nss) is set in elses_nrl_ini
     if (config%system%cutoff_radius < 0.99d0*huge(1d0)) then
        r_cut_book=config%system%cutoff_radius
     else   
        r_cut_book=ddd
     endif   
     if (awt(nss) < 1.0d-10) then
       write(*,*)'ERROR:load_element_info'
       write(*,*)'awt(nss)=',awt(nss)
       stop
     endif
     nval(nss)=9
     val_elec_atm(nss)=config%system%structure%velement(nss)%classic%charge
   endif  
!
   config%system%cutoff_radius=r_cut_book
!
   if (config%calc%use_integer_elec_num) then
     do nss=1,nos
      if (log_unit > 0) write(log_unit,'(a,i5,f40.20)') &
&           'INFO:integer_elec_num corrction (befor) : val_elec_atom=',nss, val_elec_atm(nss)
      val_elec_atm(nss)=anint(val_elec_atm(nss)) 
      if (log_unit > 0) write(log_unit,'(a,i5,f40.20)') &
&           'INFO:integer_elec_num corrction (after) : val_elec_atom=',nss, val_elec_atm(nss)
     enddo   
   endif
!
   if (ierr == 1) then
     write(*,*)'ERROR(load_element_info)'
     write(*,*)'c_system =',c_system
     write(*,*)'Stop'
     stop
   endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!
   nvl=maxval(nval)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!
   if (i_verbose >= 1) then
     write(*,*)'       nos,nvl   =',nos,nvl
     write(*,*)'   r_cut_book    =',r_cut_book
     do nss=1,nos
       write(*,*)'           nss   =',nss
       write(*,*)'elem_name(nss)   =',elem_name(nss)
       write(*,*)'      awt(nss)   =',awt(nss)
       write(*,*)'val_elec_atm(nss)=',val_elec_atm(nss)
     enddo
   endif  
!
  end subroutine load_element_info
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine set_tot_elec_num
!
!   module variables: changed
!     -->   val_elec_tot     : total electron number
!           tot_elec_eig_leg : the same quantity
!
!   Note:  'tot_elec_eig_leg' is also set as the total elenctron nubmer,
!            just for compatibility to older code
!
    use M_config
    use M_qm_domain,         only : nos, atm_element
    use elses_mod_sim_cell,  only : noa
!   use elses_mod_tx,        only : jsei
    use elses_mod_val_elec,  only : val_elec_atm, val_elec_tot
!
    use elses_mod_eig_leg,   only :  tot_elec_eig_leg
!
    implicit none
    integer js, nss
    real(8) ddsum
    logical, parameter :: debug_mode = .true.
    integer           :: log_unit
!
   log_unit=config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      write(*,*)' @@ set_tot_elec_num:nos=',nos
      do nss=1,nos
        write(*,*)'nss,val_elec_atm=',nss,val_elec_atm(nss)
      enddo   
    endif   
!
    ddsum=0.0d0
    do js=1,noa
!     nss=jsei(js)
      nss=atm_element(js)
      if (debug_mode) then
        if ((nss < 1) .or. (nss > size(val_elec_atm,1))) then
          write(*,*)'ERROR(set_tot_elec_num):js,nss=',js, nss 
          stop
        endif   
      endif  
      ddsum=ddsum+val_elec_atm(nss)
    enddo   
    val_elec_tot=ddsum
!
    if (config%calc%use_integer_elec_num) then
      if (log_unit > 0) write(log_unit,'(a,f40.20)') &
&           'INFO:integer_elec_num corrction (befor) : val_elec_tot=', val_elec_tot
      val_elec_tot=anint(val_elec_tot)
      if (log_unit > 0) write(log_unit,'(a,f40.20)') &
&           'INFO:integer_elec_num corrction (after) : val_elec_tot=', val_elec_tot
    endif
!
    if (val_elec_tot <= 1.0d-10) then
       write(*,*)'ERROR:set_tot_elec_num'
       write(*,*)'val_elec_tot=',val_elec_tot
       write(*,*)'stop'
       stop
    endif
!   
    if (i_verbose >= 1) then
       write(*,*)' total number of elec.: val_elec_tot=',val_elec_tot
    endif  
!
    tot_elec_eig_leg=val_elec_tot
    if (i_verbose >= 1) then
       write(*,*)' Note: tot_elec_eig_leg is also set'
       write(*,*)'           for compatibility to old code'
    endif  
!
  end subroutine set_tot_elec_num
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!@@ Check the consistency of input data
!
!   No module data is modifed in this routine.
!
  subroutine chk_input_data
    use elses_mod_sel_sys,   only : r_cut_book
    use M_qm_domain, only : ax, ay, az, i_pbc_x, i_pbc_y, i_pbc_z
!
!
!
    if (i_verbose >= 1) then
       write(*,*)' @@ chk_input_data'
    endif  
!
    if (r_cut_book <= 1.0d-10) then
      write(*,*)'ERROR(chk_input_data)'
      write(*,*)' r_cut_book=',r_cut_book
      stop
    endif   
!
    if (i_pbc_x == 1) then
       if (ax <= 2.0d0*r_cut_book) then
         write(*,*)'ERROR(chk_input_data):ax=',ax
         write(*,*)'      r_cut_book        =',r_cut_book
         write(*,*)' Too small simulation cell is not supported'
         stop
       endif   
    endif   
!
    if (i_pbc_y == 1) then
       if (ay <= 2.0d0*r_cut_book) then
         write(*,*)'ERROR(chk_input_data):ay=',ay
         write(*,*)'      r_cut_book        =',r_cut_book
         write(*,*)' Too small simulation cell is not supported'
         stop
       endif   
    endif   
!
    if (i_pbc_z == 1) then
       if (az <= 2.0d0*r_cut_book) then
         write(*,*)'ERROR(chk_input_data):az=',az
         write(*,*)'      r_cut_book        =',r_cut_book
         write(*,*)' Too small simulation cell is not supported'
         stop
       endif   
    endif   
!
    if (i_verbose >= 1) then
       write(*,*)' ... ended successfully; chk_input_data'
    endif  
!
  end subroutine chk_input_data
!
end module M_ini_load_geno

