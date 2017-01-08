!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dstm_ini_setting
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_wall_clock_time,  only : get_system_clock_time !(routine)
   implicit none
   logical, allocatable :: CalledFirstFlag(:)
!
   private
   public :: qm_domain_ini_setting_dstm
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine qm_domain_ini_setting_dstm
!       called only at the first MD step
!
     use M_config,           only : config !(unchanged)
     use M_qm_domain,        only : DOUBLE_PRECISION !(parameter)
!
!    use M_qm_domain,        only : jsv4jsd, njsd !(CHANGED)
!    use elses_mod_ene,      only : enea          !(CHANGED)
     use elses_mod_sel_sys,  only : r_cut_book !(unchanged)
     use elses_mod_sim_cell, only : noa ! (unchanged)
     use elses_mod_noav,     only : noao,noab,noas,noab0,nncut ! (CHANGED)
     use M_qm_domain,        only : noav, nval, c_system !(unchanged)
!    use M_qm_domain,        only : dhij, dsij, ddhij, ddsij, dbij, dpij !(CHANGED)
!    use elses_arr_dbij2,    only : dbij2  !(CHANGED)
     use elses_mod_js4jsv,   only : jsv4js, js4jsv !(CHANGED)
     use M_qm_domain,        only : atm_force, atm_force_tb0, atm_force_csc !(CHANGED)
     use M_qm_domain,        only : get_length_of_list, qm_domain_allocate_csc !(routine)
     use M_md_set_bk_list_w_cell, only : make_bk_list_w_cell !(routine)
     use M_qm_alloc_charge,    only : allocate_charge_arrays !(routine)
     use M_qm_alloc_charge,    only : allocate_csc_arrays    !(routine)
!
     implicit none
     integer ierr
     integer nval_max
     real(DOUBLE_PRECISION) :: cutoff_radius
!    integer length_of_list
     logical :: non_orderN_memory
     real(DOUBLE_PRECISION) :: memory_size  ! Estimated memory size in GB
!    real(DOUBLE_PRECISION) :: ddd1
     integer :: noa_def, noav_def, noas_def, noao_def, nvl_def
     integer :: js, jsv
     integer :: n_csc_loop
     logical :: CalledFirst
!    logical, parameter :: use_jsv4jsd = .false.
!
     n_csc_loop = config%calc%genoOption%CSC_max_loop_count
!
     memory_size=0.0d0
     nvl_def=maxval(nval)
!
     if (i_verbose >= 1) then
       write(*,*)'@@ qm_domain_ini_setting_dstm'
       write(*,*)'  solver scheme=',trim(config%calc%solver%scheme)
     endif
!
     non_orderN_memory = .false. 
     if (trim(config%calc%solver%scheme) == 'eigen') then
        non_orderN_memory = .true. 
     endif   
!
!    call qm_domain_ini_compat
!
     nval_max=maxval(nval)
     nncut=2
!
     cutoff_radius=r_cut_book
!
     if (.not. allocated(CalledFirstFlag)) then
       CalledFirst = .true.
       allocate(CalledFirstFlag(1),stat=ierr)
       if( ierr /= 0) stop 'Alloc. error for CalledFirstFlag'
     else
       CalledFirst = .false.
     endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    if (CalledFirst) then
!      nncut_def=1
!      allocate (njsd(noav,0:nncut_def),stat=ierr)
!      if( ierr /= 0) stop 'Alloc. error for njsd'
!    endif  
!
!    njsd(:,0)=1 ! zero-th shell for the 'self'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocate the booking list with the estimation of the required array size
!     only at the first iteration
!
!    if (CalledFirst) then
!      njsd(:,1)=0
!      imode=1
!      call make_bk_list_w_cell(imode,cutoff_radius,njsd(:,1))
!      call set_booking_list_w_cell(imode,cutoff_radius,njsd(:,1))
!      length_of_list=maxval(njsd(:,1))
!
!      if (i_verbose >= 1) then
!        write(*,*)'INFO:length of list =',length_of_list
!        if (non_orderN_memory) then
!          write(*,*)'INFO: NOAO=NOAV (Non-order-N memory for eigen solver)'
!        else  
!          write(*,*)'INFO: NOAO is set to be min (noav, 2*length_of_list) at the first MD step'
!        endif  
!      endif
!
!      if (non_orderN_memory) then
!        noao  = noav
!      else  
!        noao  = min(length_of_list*2,noav)
!      endif  
!
!      if (use_jsv4jsd) then
!        allocate (jsv4jsd(noao,noav),stat=ierr)
!        if (ierr /= 0) stop 'Error in alloc. : jsv4jsd'
!      else
!        write(*,*)'INFO:jsv4jsd is NOT allocated'
!      endif
!
!    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set the booking list ---> njsd(:,1), jsv4jsd(:,:)
!
!    njsd(:,1)=0
!    if (use_jsv4jsd) then
!      jsv4jsd(:,:)=0
!      imode=2
!      call make_bk_list_w_cell(imode,cutoff_radius,njsd(:,1),jsv4jsd)
!      call set_booking_list_w_cell(imode,cutoff_radius,njsd(:,1),jsv4jsd)
!    else
!      imode=1
!      call make_bk_list_w_cell(imode,cutoff_radius,njsd(:,1))
!    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if (.not. CalledFirst) then 
       if (i_verbose >= 1) then
         write(*,*)'...return without allocation'
       endif   
       return
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     noao=0   ! TEMPORAL SETTING (T. Hoshi, 01.Oct.2011)
!
     noas  = noao
     noab  = noao
     noab0 = noas
!
     noa_def =noa
     noav_def=noav
     noao_def=noao
     noas_def=noas
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocations of arrays
!
     allocate (js4jsv(noav),stat=ierr)
     if (ierr /= 0) stop 'Error in alloc. : js4jsv'
     do js=1,noav
       js4jsv(js)=js
     enddo   
!
     allocate (jsv4js(noav),stat=ierr)
     if (ierr /= 0) stop 'Error in alloc. : jsv4js'
     do js=1,noav
       jsv4js(js)=js
     enddo   
!
!    allocate (enea(noa_def,2),stat=ierr)
!    if (ierr /= 0) stop 'Error in alloc. : enea'
!    enea(:,:)=0.0d0
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     do jsv=1,noav_def
       js=jsv
       jsv4js(js)=jsv
       js4jsv(jsv)=js
     enddo   
!     
     if (config%system%structure%use_vatom) then
       allocate (atm_force(3,noav),stat=ierr)
       if( ierr .ne. 0) then
          write(6,*)'alloc. error!(atom_force):ierr=',ierr
          stop
       endif
       if (i_verbose >= 1) then
         if (log_unit > 0) then
           write(log_unit, *)'INFO:Alloc. of atm_force     : size [GB] =', & 
&                                8.0d0*3.0d0*dble(noav)/1.0d9
         endif   
       endif   
     endif  
!
     if (trim(config%calc%calc_force_mode) /= 'off') then
!
       allocate (atm_force_tb0(3,noav),stat=ierr)
       if( ierr .ne. 0) then
         write(6,*)'alloc. error!(atom_force_tb0):ierr=',ierr
         stop
       endif
       if (i_verbose >= 1) then
         if (log_unit > 0) then
           write(log_unit, *)'INFO:Alloc. of atm_force_tb0 : size [GB] =', & 
&                                8.0d0*3.0d0*dble(noav)/1.0d9
         endif   
       endif   
!
       allocate (atm_force_csc(3,noav),stat=ierr)
       if( ierr .ne. 0) then
         write(6,*)'alloc. error!(atm_force_csc):ierr=',ierr
         stop
       endif
       if (i_verbose >= 1) then
         if (log_unit > 0) then
           write(log_unit, *)'INFO:Alloc. of atm_force_csc : size [GB] =', & 
&                                8.0d0*3.0d0*dble(noav)/1.0d9
         endif   
       endif   
!
     endif  
!
     if (c_system == "geno") then
       if (config%system%structure%use_vatom) then
         call allocate_charge_arrays(nval_max, noav)
       else
         if (log_unit > 0) then
           write(log_unit, *)'INFO:allocate_charge_arrays is skipped (no use_vatom)'
         endif   
       endif  
       if (n_csc_loop > 0) call allocate_csc_arrays(nval_max, noas, noav)
     endif
!  
   end subroutine qm_domain_ini_setting_dstm
!
!
end module M_qm_dstm_ini_setting


