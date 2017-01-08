!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_engine_dst
!
   use M_qm_domain,            only : i_verbose, c_system, &
        ham_tb0, ham_csc, dhij, ddhij, dham_tb0, jsv4jsd, DOUBLE_PRECISION
   use M_io_dst_write_log,     only : log_file_is_set, log_unit !(unchanged)
!
   private
   public qm_engine_dst
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_engine_dst
!
!   module variables unchanged
!      --> 
!   module variables changed
!      --> 
!
     use elses_mod_md_dat,     only : itemd, itemdmx
!    use elses_mod_md_dat,     only : e_kin
!    use elses_mod_foi,        only : foi
     use elses_mod_ene,        only : etb, ecc
     use elses_mod_ene,        only : ecsc !(CHANGED)
     use elses_mod_md_dat,     only : itemdorg
!
     use M_config, only : config !(unchanged except below)
!                       ! CHANGED : config%system%structure%mdstep
     use M_qm_solver_geno_main, &
&         only : qm_solver_geno_main, qm_solver_ortho_main
     use M_qm_domain, &
&         only : qm_domain_setting, generate_foi, mulliken, &
&         atm_force, atm_force_csc
     use M_qm_geno, &
&         only : set_hamiltonian_and_overlap_geno, set_rest_geno, &
&         set_qm_force_geno, renew_charge, rms_delta_q, &
&         initialize_charge
     use M_qm_geno_csc, &
&         only : set_CSC_parameters_geno, set_hamiltonian_csc, &
&         set_atm_force_csc
     use M_qm_geno_output, only : qm_geno_output
     use M_qm_nrl, only : qm_nrl_set_ham, qm_nrl_set_force
     use M_qm_engine_nrl, only : qm_engine_nrl
     use M_qm_engine_csc,     only : qm_engine_csc
     use M_qm_engine_csc_tb0, only : qm_engine_csc_tb0
!    use M_md_velocity_routines, only : calc_kinetic_energy
     use M_qm_solver_gkrylov_dst,  only : qm_solver_gkrylov_dst  !(routine)
!
     use M_qm_domain_dst, only : qm_domain_setting_dst !(routine)
     use M_qm_domain_dst, only : global_dens_mat !(CHANGED)
     use M_qm_domain_dst, only : global_ham_mat  !(CHANGED)
     use M_qm_domain_dst, only : overlap_is_on   !(CHANGED)
!
     use M_qm_engine_csc, only : tune_mix_ratio !(routine)
!
!    use M_qm_geno_dst,   only : set_ham_tb0_and_overlap_dst !(routine)
!
!    use M_qm_dst_global_ham_mat, only : set_mat_tb0_global
!
     use M_md_dst,        only : set_dst_final
     use M_qm_geno_CSC_dst, only : set_atm_force_csc_dst !(routine)
!
     use M_qm_geno_output, only :  qm_calc_ecsc !(routine)
     use M_wall_clock_time, only : get_system_clock_time !(routine)
!
!    use M_md_dst_cell,     only :  initial_cell_booking !(routine)
!
     use M_cohp_dstm_plot,  only : prep_cohp_dst !(routine)
!
     use M_qm_domain,       only : noav
!
     use M_qm_dst_matrix_gene, only : qm_dst_mat_gene !(routine)
!
     implicit none
     real(DOUBLE_PRECISION) :: dq, dq_1, x, mix_r
     integer :: ierr, i_csc_loop, jsd1, jsv2, ja1, ja2, n_csc_loop
     character(len=8) :: CSC_METHOD
     real(DOUBLE_PRECISION) :: kinetic_energy
     character(len=32) :: scheme_mode
!
     real(DOUBLE_PRECISION), allocatable :: dq_hist(:),mix_r_hist(:)
!         ----> history for dq and mix_r
     integer  :: k, mix_mode
     real(DOUBLE_PRECISION) :: value_of_ecsc
     real(DOUBLE_PRECISION) :: time_wrk, time_wrk_previous ! work variable for measuring the time
!
     CSC_METHOD = config%calc%genoOption%CSC_method
     n_csc_loop = config%calc%genoOption%CSC_max_loop_count
     mix_mode   = config%calc%genoOption%CSC_mode_for_tuning_mixing_ratio
!
     scheme_mode=trim(config%calc%solver%scheme)
!
     itemd=itemd+1
     config%system%structure%mdstep = itemdorg+itemd-1
!
     if (i_verbose >= 1) then
       if (log_unit > 0) then
         write(log_unit,*)'@@ elses_qm_engine-dst:itemd=',itemd
       endif  
     endif   
!
     call get_system_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
!    call elses_calc_kin_ene
!    call calc_kinetic_energy(kinetic_energy)
!    e_kin=kinetic_energy
!    if (i_verbose >= 1) then
!      write(*,*)'  kinetic energy [au]=',e_kin
!    endif   
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:calc_kinetic_energy =',time_wrk-time_wrk_previous
     if (log_file_is_set) write(log_unit,'(a,f20.10)') & 
&                         'TIME:qm_engine_dst:calc_kinetic_energy =',time_wrk-time_wrk_previous
     time_wrk_previous=time_wrk
!
     if (config%calc%distributed%global_dens_mat) then
       global_dens_mat =.true.
     else
       global_dens_mat =.false.
     endif
!
     if (config%calc%distributed%global_ham_mat) then
       global_ham_mat  =.true.
     else
       global_ham_mat  =.false.
     endif
!
     overlap_is_on   =.true.
!
     call get_system_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
     call qm_domain_setting_dst
!     
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst =',time_wrk-time_wrk_previous
     if (log_file_is_set) write(log_unit,'(a,f20.10)') & 
&                         'TIME:qm_engine_dst:qm_domain_setting_dst =',time_wrk-time_wrk_previous
     time_wrk_previous=time_wrk
!
!    write(*,*)'noav=',noav
     call prep_cohp_dst
!
!    if (c_system /= "geno") then
!      write(*,*)'ERROR(qm_engine_dst):c_system=',c_system
!      stop
!    endif   
!
     if (CSC_METHOD /= "ELSTNER") then
       write(*,*)'ERROR(qm_engine_dst):CSC_METHOD=',CSC_METHOD
       stop
     endif   
!
     if (n_csc_loop > 0) then 
       allocate (dq_hist(n_csc_loop), stat=ierr)
       if (ierr /= 0) stop 'Abort:alloc. error'
       dq_hist(:)=0.0d0
!
       allocate (mix_r_hist(n_csc_loop), stat=ierr)
       if (ierr /= 0) stop 'Abort:alloc. error'
       mix_r_hist(:)=0.0d0
!
       mix_r = config%calc%genoOption%CSC_charge_mixing_ratio ; dq=1d0
     endif  
!
     if (i_verbose >=1) then 
       if (log_unit > 0) then
         write(log_unit,'(a,i5)')'@@ qm_engine_dst : mode_for_tuning_mixing_ratio = ',mix_mode
       endif  
     endif  
!
     if (allocated(atm_force)) then 
       atm_force(:,:)=0.0d0
     else
       if (i_verbose >=1) then 
         if (log_unit > 0) then
           write(log_unit,'(a)')'INFO:Zero clear of atm_force is skiiped (not allocated)'
         endif  
       endif  
     endif
  
!    enea(:,:)=0.0d0
     ecc=0.0d0
     etb=0.0d0
!
!
     if (c_system == "geno") then
       if (itemd == 1) then 
         if (config%system%structure%use_vatom) then  
           call initialize_charge
         else
           if (i_verbose >= 1) then
             if (log_unit > 0) then
               write(log_unit,*)'INFO:initialize_charge is skipped (no use_vatom)' 
             endif   
           endif   
         endif
  
       endif   
     endif  
!
     call get_system_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
!    call set_rest_geno 
       ! initialize atm_force & set repulsive force 
       ! & calculate repulsive energy
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:set_rest_geno =',time_wrk-time_wrk_previous
     if (log_file_is_set) write(log_unit,'(a,f20.10)') & 
&                         'TIME:qm_engine_dst:set_rest_geno =',time_wrk-time_wrk_previous
     time_wrk_previous=time_wrk
!
!    call set_hamiltonian_and_overlap_geno
!     dhij(:,:,:,:)  = ham_tb0(:,:,:,:)
!    ddhij(:,:,:,:,:)=dham_tb0(:,:,:,:,:)
!
!    call set_ham_tb0_and_overlap_dst
!    call set_mat_tb0_global
!
!    call set_dst_final
!    stop 'Stop manually'
!
!    if (n_csc_loop > 0) then 
!      call set_CSC_parameters_geno
!    endif
!
     if (n_csc_loop == 0) then 
!      scheme_mode='ekrgl'
!
       call get_system_clock_time(time_wrk)
       time_wrk_previous=time_wrk
!
       if (config%calc%mode == "matrix_generation") then
         call qm_dst_mat_gene(scheme_mode)
         return
       endif
!
       call qm_solver_gkrylov_dst(scheme_mode)
!
       call get_system_clock_time(time_wrk)
!      write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_solver_gkrylov_dst =',time_wrk-time_wrk_previous
       if (log_file_is_set) write(log_unit,'(a,f20.10)') & 
&                           'TIME:qm_engine_dst:qm_solver_gkrylov_dst =',time_wrk-time_wrk_previous
       time_wrk_previous=time_wrk
!
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CSC loop
     if (n_csc_loop > 0) then 
       do i_csc_loop=1, n_csc_loop
!        
!        scheme_mode='ekrgl'
         call qm_solver_gkrylov_dst(scheme_mode)
!
         dq_1= dq
         dq  = rms_delta_q()
         dq_hist(i_csc_loop)=dq
         if( dq < config%calc%genoOption%CSC_charge_convergence ) then
            if(i_verbose > 0) then
              write(*,'(a,i10,E13.6)') 'INFO:CSC loop is converged:loop num., dq=',i_csc_loop, dq
            endif  
            exit
         endif   
!
         if (mix_mode /= 0) then
           call tune_mix_ratio(dq_hist, mix_r_hist,i_csc_loop, mix_r,mix_mode)
           mix_r_hist(i_csc_loop)=mix_r
           if(i_verbose >= 0) write(*,*)'mix_r_hist=', mix_r_hist(i_csc_loop)
         endif  
!
         call renew_charge(mix_r)
          ! renew e_num_on_basis, e_num_on_atom  (charge mixing)
!          
         if(i_verbose > 0) write(*,'("dq=",ES10.2,", mix_r=",ES10.2)') dq,mix_r
!         
       end do
     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if (allocated(dq_hist)) then
       deallocate (dq_hist, stat=ierr)
       if (ierr /= 0) stop 'Abort:dealloc. error'
     endif  
!
     if (allocated(mix_r_hist)) then
       deallocate (mix_r_hist, stat=ierr)
       if (ierr /= 0) stop 'Abort:dealloc. error'
     endif  
!
     if (i_verbose > 0 .and. n_csc_loop >0) then 
        write(*,'("mix_r=",ES10.2,", i_csc_loop=",I8)')mix_r, i_csc_loop
       if(i_csc_loop > n_csc_loop .and. n_csc_loop > 0) then
         stop "elses_qm_engine: !!!ERROR!!! charge did't converge"
       end if
     endif   
!
     call get_system_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
     call calc_total_force
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:calc_total_force =',time_wrk-time_wrk_previous
     if (log_file_is_set) write(log_unit,'(a,f20.10)') & 
&                         'TIME:qm_engine_dst:calc_total_force =',time_wrk-time_wrk_previous
     time_wrk_previous=time_wrk
!
     call qm_geno_output
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_geno_output   =',time_wrk-time_wrk_previous
     if (log_file_is_set) write(log_unit,'(a,f20.10)') & 
&                         'TIME:qm_engine_dst:qm_geno_output   =',time_wrk-time_wrk_previous
     time_wrk_previous=time_wrk
!
     if ( n_csc_loop > 0 ) then
       call qm_calc_ecsc(value_of_ecsc)
       ecsc=value_of_ecsc
       if (i_verbose >= 1) then
         write(6,*) 'n_csc_loop=',n_csc_loop
         write(6,*) 'Ecsc as 0.5*sum_{alpha,beta}gamma*dq_{alpha}*dq_{beta}=', value_of_ecsc
         write(6,*) '--> set as ECSC'
       endif   
     endif  
!
   end subroutine qm_engine_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine calc_total_force
!
!
     use M_config,           only : config !(unchanged)
     use elses_mod_js4jsv,   only : js4jsv !(unchanged)
     use elses_mod_foi,      only : foi    !(CHANGED)
     use M_qm_domain,        only : atm_force     !(CHANGED)
     use M_qm_domain,        only : atm_force_tb0 !(unchanged)
     use M_qm_domain,        only : atm_force_csc !(unchanged)
     use M_qm_domain,        only : noav !(unchanged)
!
     implicit none
     integer :: jsv,js
     logical :: error_flag
!
     if (i_verbose >= 1) then
       if (log_unit > 0) write(log_unit,'(a)') '@@ generate_foi'
     endif   
!
     if (trim(config%calc%calc_force_mode) == 'off') return
!
     error_flag = .false.
     if (.not. allocated(atm_force)) error_flag = .true. 
     if (.not. allocated(atm_force_tb0)) error_flag = .true. 
     if (.not. allocated(atm_force_csc)) error_flag = .true. 
     if (.not. allocated(foi)) error_flag = .true. 
!
     if (error_flag) then
       write(*,*)'ERROR(calc_total_force):error_flag=', error_flag
       stop 
     endif   
!
     atm_force(:,:)=atm_force(:,:)+atm_force_tb0(:,:)+atm_force_csc(:,:)
!
     foi(:,:)=0.0d0
!
     do jsv=1,noav
       js=js4jsv(jsv)
       foi(js,1:3)=atm_force(1:3,jsv)
     enddo
!   
    if (i_verbose >=1) then
      do jsv=1,min(noav,3)
       if (log_unit > 0) write(log_unit,'(a,i8,3f20.10)')'  foi=',jsv,atm_force(1:3,jsv)
      enddo
    endif
! 
   end subroutine calc_total_force
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_engine_dst
