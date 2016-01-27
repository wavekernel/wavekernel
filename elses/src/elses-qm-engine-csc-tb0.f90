!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_engine_csc_tb0
!
   use M_qm_domain,            only : i_verbose, c_system, &
&       ham_tb0, ham_csc, dhij, ddhij, dham_tb0, jsv4jsd, DOUBLE_PRECISION
!
!  private
   public qm_engine_csc_tb0
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_engine_csc_tb0
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
     use M_config, only : config  !(unchanged)
     use M_qm_solver_geno_main,  only : qm_solver_geno_main !(routine)
     use M_qm_domain,  only :  mulliken
     use M_qm_geno, &
          only : set_hamiltonian_and_overlap_geno, set_rest_geno, &
          set_qm_force_geno, renew_charge, rms_delta_q
!    use M_qm_geno_csc, &
!         only : set_CSC_parameters_geno, set_hamiltonian_csc, &
!         set_atm_force_csc
!    use M_qm_geno_output, only : qm_geno_output
!    use M_qm_nrl, only : qm_nrl_set_ham, qm_nrl_set_force
!    use M_qm_engine_nrl, only : qm_engine_nrl
!
     implicit none
     real(DOUBLE_PRECISION) :: dq, dq_1, x, mix_r
     integer :: ierr, i_csc_loop, jsd1, jsv2, ja1, ja2, n_csc_loop
     character(len=8) :: CSC_METHOD
!
     if (i_verbose > 0) then
       write(*,*)'@@ elses_qm_engine_geno_csc_tb0'
       write(*,*)'WARNING : This routine is NOT well tested'
     endif   

     CSC_METHOD = config%calc%genoOption%CSC_method
     n_csc_loop = config%calc%genoOption%CSC_max_loop_count
!
     if( CSC_METHOD /="TB0") then
         stop "elses_qm_engine_csc_tb0: !!!ERROR!!! unexpected CSC_METHOD"
     endif   
!
     if( n_csc_loop <= 0 ) then
         stop "ERROR:elses_qm_engine_csc_tb0"
     endif
!   
     mix_r = config%calc%genoOption%CSC_charge_mixing_ratio ; dq=1d0
     if(CSC_METHOD=="TB0" .and. n_csc_loop>0 )&
          call renew_charge(mix_r)
     !  If above two lines are disabled, Hamiltonian matrix does not change
     ! until i_csc_loop=2, in TB0-CSC mode.
     do i_csc_loop=1, n_csc_loop
        select case (CSC_METHOD)
        case("TB0")
           call set_hamiltonian_and_overlap_geno
            dhij(:,:,:,:)  = ham_tb0(:,:,:,:)
           ddhij(:,:,:,:,:)=dham_tb0(:,:,:,:,:)
        case default
           stop "elses_qm_engine: !!!ERROR!!! unexpected CSC_METHOD"
        end select
        call qm_solver_geno_main ! input dhij(Hamiltonian), dsij(overlap)
        call mulliken           ! calculate e_num_on_basis, e_num_on_atom,
                                !to be used in the next csc iteration
        dq_1= dq
        dq  = rms_delta_q()
        if( dq < config%calc%genoOption%CSC_charge_convergence ) then
           if (i_verbose >= 1) then
             write(*,*) 'INFO:CSC loop is converged:loop num., dq=',i_csc_loop, dq
           endif
           exit
        endif   
!
        mix_r = max(mix_r, config%calc%genoOption%CSC_charge_mixing_ratio)
        mix_r = min(mix_r, 0.5d0)
        ! clip CSC_charge_mixing_ratio
        call renew_charge(mix_r)
        ! renew e_num_on_basis, e_num_on_atom  (charge mixing)
          
        if(i_verbose > 0) write(*,'("dq=",ES10.2,", mix_r=",ES10.2)') dq,mix_r
       
     end do
!
     if(i_verbose > 0 .and. n_csc_loop >0) write(*,'("mix_r=",ES10.2,", i_csc_loop=",I8)')&
            mix_r, i_csc_loop
!
     if(i_csc_loop > n_csc_loop .and. n_csc_loop > 0) then
        stop "elses_qm_engine: !!!ERROR!!! charge did't converge"
     end if
!
   end subroutine qm_engine_csc_tb0
!
end module M_qm_engine_csc_tb0


