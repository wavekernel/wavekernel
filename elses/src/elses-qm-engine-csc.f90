!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
module M_qm_engine_csc
!
   use M_qm_domain,            only : i_verbose, c_system, &
        ham_tb0, ham_csc, dhij, ddhij, dham_tb0, jsv4jsd, DOUBLE_PRECISION
!
   private
   public qm_engine_csc
   public tune_mix_ratio
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_engine_csc
!
!   module variables unchanged
!      --> 
!   module variables changed
!      --> 
!
!    use elses_mod_md_dat,     only : itemd, itemdmx
!    use elses_mod_foi,        only : foi
!    use elses_mod_ene,        only : etb, ecc, enea
!
!
!    Module procedures
     use M_config,         only : config   ! (unchanged)
     use M_qm_domain,      only : noav     ! (unchanged)
     use M_lib_phys_const, only : ev4au    ! (unchanged)
     use M_qm_solver_geno_main, &
          only : qm_solver_geno_main, qm_solver_ortho_main
     use M_qm_domain, &
          only : qm_domain_setting, generate_foi, mulliken, &
          atm_force, atm_force_csc
     use M_qm_geno, &
          only : set_hamiltonian_and_overlap_geno, set_rest_geno, &
          set_qm_force_geno, renew_charge, rms_delta_q, &
          initialize_charge
     use M_qm_geno_csc, &
          only : set_CSC_parameters_geno, set_hamiltonian_csc, &
          set_atm_force_csc
     use M_qm_geno_output, only : qm_geno_output, qm_calc_ecsc ! (routine)
!    use M_qm_nrl, only : qm_nrl_set_ham, qm_nrl_set_force
!    use M_qm_engine_nrl, only : qm_engine_nrl
!
     use M_qm_domain,     only : csc_ite_converged, csc_dq_converged
     use M_qm_geno_CSC_basic, only : gamma_csc_func_plot !(subroutine)
!
     implicit none
     real(DOUBLE_PRECISION)              :: dq, dq_1, mix_r
     real(DOUBLE_PRECISION), allocatable :: dq_hist(:),mix_r_hist(:)
     integer  :: k, ierr, i_csc_loop, n_csc_loop, mix_mode
     character(len=8) :: CSC_METHOD
!         ----> history for dq and mix_r
     real(DOUBLE_PRECISION)              :: mix_r0, mix_r_wrk
     integer  :: lu, step_count
     logical  :: plot_ecsc, plot_ecsc_atom ! Plot E_csc in the CSC calc. 
     real(DOUBLE_PRECISION)              :: value_of_ecsc
     real(DOUBLE_PRECISION), allocatable :: atom_csc_energy(:)
!
     lu = config%calc%distributed%log_unit
!
     CSC_METHOD = config%calc%genoOption%CSC_method
     n_csc_loop = config%calc%genoOption%CSC_max_loop_count
     mix_mode   = config%calc%genoOption%CSC_mode_for_tuning_mixing_ratio
     step_count = config%system%structure%mdstep
!
     plot_ecsc = .false.
     plot_ecsc_atom = .false.
!
     if (n_csc_loop < 1) then
       write(*,*)'ERROR(qm_engine_csc):n_csc_loop=',n_csc_loop
       stop
     endif
!
     if (n_csc_loop /= 0) then
       if (i_verbose >= 5)  plot_ecsc      = .true.
       if (i_verbose >= 7)  plot_ecsc_atom = .true.
     endif
!
     if (step_count == 0) then
       if (i_verbose >= 1) call gamma_csc_func_plot ! Plot gamma CSC function for demonstration
     endif
!
     if (lu > 0) write(lu,'(a,i5)')'@@ qm_engine_csc : mode_for_tuning_mixing_ratio = ',mix_mode
!
     allocate (dq_hist(n_csc_loop), stat=ierr)
     if (ierr /= 0) stop 'Abort:alloc. error'
     dq_hist(:)=0.0d0
!
     allocate (mix_r_hist(n_csc_loop), stat=ierr)
     if (ierr /= 0) stop 'Abort:alloc. error'
     mix_r_hist(:)=0.0d0
!
     mix_r = config%calc%genoOption%CSC_charge_mixing_ratio ; dq=1d0
!
     do i_csc_loop=1, n_csc_loop
         call set_CSC_parameters_geno
         call set_hamiltonian_csc
             ! Simple checking routine of Ham_{\rm csc}
         if ( i_verbose >= 100 ) call write_ham_csc
!$omp parallel do
         do k=1, size(dhij,4)
            dhij(:,:,:,k)  = ham_tb0(:,:,:,k)  + ham_csc(:,:,:,k)
            !           Generate dhij as H = H_{TB0} + H_{CSC}
            ddhij(:,:,:,:,k)=dham_tb0(:,:,:,:,k)!+dham_csc(:,:,:,:,k) !
            ! Differentiation of 2nd order (rho^2) term H_{CSC}
            ! is a little bit more complicated
            ! than that of H_{TB0} term (rho^1).
            ! Therefore, we must prepare a specific subroutine
            ! for each term.
         end do
        call qm_solver_geno_main ! input dhij(Hamiltonian), dsij(overlap)
        call mulliken           ! calculate e_num_on_basis, e_num_on_atom,
                                !to be used in the next csc iteration
        dq_1= dq
        dq  = rms_delta_q()
        dq_hist(i_csc_loop)=dq
        if ( dq < config%calc%genoOption%CSC_charge_convergence ) then
          if (lu > 0) write(lu,'(a,i10,E13.6)') 'INFO:CSC loop is converged:loop num., dq=',i_csc_loop, dq
          csc_ite_converged = i_csc_loop
          csc_dq_converged  = dq
          exit
        endif   
!
!       mix_r = max(mix_r, config%calc%genoOption%CSC_charge_mixing_ratio)
!       mix_r = min(mix_r, 0.5d0)
        ! clip CSC_charge_mixing_ratio
!
        if (mix_mode /= 0) then
           call tune_mix_ratio(dq_hist, mix_r_hist,i_csc_loop, mix_r,mix_mode)
           mix_r_hist(i_csc_loop)=mix_r
           if (lu > 0) write(lu,*)'mix_r_hist=', mix_r_hist(i_csc_loop)
        endif  
!
        if (i_csc_loop == 1) then
          mix_r0=1.0d0
          call renew_charge(mix_r0)
          mix_r_wrk=mix_r0
          ! renew e_num_on_basis, e_num_on_atom  (charge mixing)
        else
          call renew_charge(mix_r)
          mix_r_wrk=mix_r
          ! renew e_num_on_basis, e_num_on_atom  (charge mixing)
        endif
!
       value_of_ecsc=0.0d0
       if (plot_ecsc) then
         call qm_calc_ecsc(value_of_ecsc,plot_ecsc_atom)
         value_of_ecsc=value_of_ecsc/dble(noav)*ev4au
       endif
!
        if (lu > 0) then
          if(i_verbose > 0) then 
            write(lu,'("INFO-CSC:step, csc_step, dq, mix_r, Ecsc[eV/atom] =",2I10,ES10.2,ES10.2,F30.20)') & 
&                            step_count, i_csc_loop, dq, mix_r_wrk, value_of_ecsc
          endif
        endif
!
     end do

     deallocate (dq_hist, stat=ierr)
     if (ierr /= 0) stop 'Abort:dealloc. error'
!
     deallocate (mix_r_hist, stat=ierr)
     if (ierr /= 0) stop 'Abort:dealloc. error'

       if(i_verbose > 0 .and. n_csc_loop >0) then
          if (lu > 0) write(lu,'("mix_r=",ES10.2,", i_csc_loop=",I8)') mix_r, i_csc_loop
       endif

       if(i_csc_loop > n_csc_loop .and. n_csc_loop > 0) then
          stop "elses_qm_engine: !!!ERROR!!! charge did't converge"
       end if

!
     contains
       
       subroutine write_ham_csc
         integer::  jsd1, jsv2, ja1, ja2
         write(*,'("elses_qm_engine: Printing Ham_{\rm csc}")') 
         do jsv2=1, size(dhij,4)
            do jsd1=1, size(dhij,3)
               write(*,'(" jsv1=",I8," jsv2=",I8)') jsv4jsd(jsd1,jsv2), jsv2
               do ja1=1, size(dhij,1)
                  write(*,'(ES23.14$)') (ham_csc(ja1,ja2,jsd1,jsv2), ja2=1, size(dhij,2))
                  write(*,'()')
               end do
               write(*,'()')
            end do
         end do
       end subroutine write_ham_csc

   end subroutine qm_engine_csc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine tune_mix_ratio(dq_hist, mix_r_hist,i_csc_loop, mix_r_value,i_level)
!
     use M_config,         only : config   ! (unchanged)
     implicit none
     real(DOUBLE_PRECISION), intent(in)    :: dq_hist(:),mix_r_hist(:)
     integer,                intent(in)    :: i_csc_loop
     real(DOUBLE_PRECISION), intent(inout) :: mix_r_value
     real(DOUBLE_PRECISION) :: dq_min, dq_value
     real(DOUBLE_PRECISION) :: dq_ratio, dq_ratio_min
     integer                :: ite_dq_min
     integer                :: ite, ite_same_mix_r
     integer                :: i_level
     integer                :: lu
!
       lu = config%calc%distributed%log_unit
       dq_value    = dq_hist(i_csc_loop)
       dq_min      = minval(dq_hist(1:i_csc_loop))
       ite_dq_min  = minloc(dq_hist(1:i_csc_loop),1)
!
       if (i_csc_loop == 1) return
!
       if (lu > 0) then
         if (i_verbose >=1) then 
           write(lu,'(a,i10,2f20.10)')'@@ tune mix ratio:level, dq, dq_min=',i_level,dq_value, dq_min
         endif
       endif
!
       do ite=2, i_csc_loop
         if ( dq_hist(ite-1) < 1.0d-14 ) then
           write(*,*) '..dq is small enough:return'
           return
         endif   
         dq_ratio=dq_hist(ite)/dq_hist(ite-1)
         if (lu > 0) then
           if (i_verbose >=1) then 
             write(lu,'(a,i10,2f20.10)')' ite, dq, dq_ratio=',ite, dq_hist(ite), dq_ratio
           endif
         endif  
         if (ite == 2) dq_ratio_min=dq_ratio 
         if (dq_ratio < dq_ratio_min) dq_ratio_min=dq_ratio 
       enddo
       dq_ratio=dq_hist(i_csc_loop)/dq_hist(i_csc_loop-1)
       dq_ratio_min=min(1.0d0, dq_ratio_min)
!
       if (i_verbose >=1) then 
         if (lu > 0) then
           write(lu,*)' dq_ratio, dq_ratio_min*1.2=',dq_ratio, dq_ratio_min*1.2d0
         endif
       endif
!
!
       ite_same_mix_r=i_csc_loop-1
       do ite=1,i_csc_loop-2
!        write(*,*)'ite, mix_r_value, mix_r_hist=',ite, mix_r_value, mix_r_hist(i_csc_loop-ite)
         if ( dabs(mix_r_value-mix_r_hist(i_csc_loop-ite)) > 1.0d-10) then 
           ite_same_mix_r=ite
           exit
         endif  
       enddo
!         ----> get ite_same_mix_r as the iteration number for the same mix_r value
!
       if (lu > 0) then
         if (i_verbose >=1) then 
           write(lu,'(a,3i10)')'  ite, ite_dq_min, ite_same_mix_r =',i_csc_loop, ite_dq_min, ite_same_mix_r
         endif
       endif  
!
       if (ite_same_mix_r < 3) return
!         ----> The mix_r value will not be changed, if it was changed recently.
!   
       select case (i_level)
          case (1)
            if ((i_csc_loop >= ite_dq_min+3) .and. (dq_value > 0.9d0*dq_min)) then
              mix_r_value=mix_r_value/2.0d0
             if (lu > 0) then
                if (i_verbose >=1) then 
                  write(lu,*)'new_mix_r =',mix_r_value
                endif
              endif  
            endif   
          case (2)
            if ((ite_same_mix_r >= 3) .and. ( dq_ratio > dq_ratio_min * 1.2d0 )) then 
              mix_r_value=mix_r_value/2.0d0
              if (lu > 0) then
                if (i_verbose >=1) then 
                  write(lu,*)'new_mix_r (2) =',mix_r_value
                endif
             endif  
            else
              if ((i_csc_loop >= ite_dq_min+3) .and. (dq_value > 0.9d0*dq_min)) then
                mix_r_value=mix_r_value/2.0d0
                if (lu > 0) then
                  if (i_verbose >=1) then 
                    write(lu,*)'new_mix_r (1) =',mix_r_value
                  endif
                endif  
              endif   
            endif   
          case default
            stop'ERROR(modify_mix_ratio):wrong i_level'
       end select
!
   end subroutine tune_mix_ratio
!
end module M_qm_engine_csc
