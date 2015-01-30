!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
module M_qm_engine
!
   use M_qm_domain,            only : i_verbose, c_system, &
        ham_tb0, ham_csc, dhij, ddhij, dham_tb0, jsv4jsd, DOUBLE_PRECISION
!
   private
   public elses_qm_engine
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine elses_qm_engine
!
!   module variables unchanged
!      --> 
!   module variables changed
!      --> 
!
     use elses_mod_md_dat,     only : itemd, itemdmx, e_kin
     use elses_mod_foi,        only : foi
     use elses_mod_ene,        only : etb, ecc, enea
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
     use M_md_velocity_routines, only : calc_kinetic_energy
     use M_qm_population, only : initial_guess_for_charge  !(routine)
     use M_qm_population, only : save_population_after_convergence !(routine)
     use M_qm_set_s_as_i, only : set_s_mat_as_i !(routine) 
!
     use M_qm_population_ortho, only : save_population_ortho !(routine)
     use M_qm_output_icohp,     only : qm_output_icohp ! (subroutine)
     use M_qm_output_matrices,  only : output_levels_matrices !(subroutine)
     use M_output_eigenstates,  only : output_eigenstates     !(subroutine)
!
     use M_ext_test_qd,         only : test_for_quantum_dynamics !(subroutine)
!
     implicit none
!    real(DOUBLE_PRECISION) :: dq, dq_1, x, mix_r
     integer :: ierr, n_csc_loop
     character(len=8) :: CSC_METHOD
     real(DOUBLE_PRECISION) :: kinetic_energy
     logical          :: result_flag
     logical, parameter  :: population_in_ortho_cases = .true.
     integer          :: imode_icohp
     logical :: matrix_generation_and_stop 
!    logical, parameter  :: matrix_generation_and_stop = .false.  ! for quick hack
!
     CSC_METHOD = config%calc%genoOption%CSC_method
     n_csc_loop = config%calc%genoOption%CSC_max_loop_count
!
     itemd=itemd+1
     config%system%structure%mdstep = itemdorg+itemd-1
!
     if (i_verbose >= 1) then
       write(*,*)'@@ elses_qm_engine:itemd=',itemd
     endif   
!
     if (config%calc%mode == 'matrix_generation') then
       matrix_generation_and_stop = .true.  
     else
       matrix_generation_and_stop = .false.  
     endif   
! S.Y Sep 21, 2013 two cases "calc%mode = 'matrix_generation'" and "calc%mode = 'write_Hamiltonian'" 
!      should be merged into single case. But before merging, compatibility check should be done.
     if(config%calc%mode == "optimization" .and. &
          & config%calc%optimization%scheme=="none" )  matrix_generation_and_stop = .true.
!
!
!    call elses_calc_kin_ene
     call calc_kinetic_energy(kinetic_energy)
     e_kin=kinetic_energy
     if (i_verbose >= 1) then
       write(*,*)'  kinetic energy [au]=',e_kin
     endif   
!
     call qm_domain_setting
!     
     foi(:,:)=0.0d0
     enea(:,:)=0.0d0
     ecc=0.0d0
     etb=0.0d0
!
     ierr=1
!
     if (trim(config%calc%solver%scheme) == "rest_part_only") then
       if (i_verbose >= 1) write(*,*)'INFO:scheme=',trim(config%calc%solver%scheme)
       if (c_system /= "geno") then
         write(*,*)'ERROR(rest_part_only mode):c_system=',c_system
         stop
       endif   
       call set_rest_geno 
       call generate_foi
       call qm_geno_output
       return
     endif   
!
     if (c_system == "geno") then
       ierr=0
       call initial_guess_for_charge
!      if (itemd == 1) call initialize_charge
       call set_rest_geno 
       ! initialize atm_force & set repulsive force 
       ! & calculate repulsive energy
       call set_hamiltonian_and_overlap_geno
       call set_s_mat_as_i(result_flag) ! S is set to be I artificialy'
!
       if (result_flag) then
         write(*,*)'INFO:S is set to be I artificialy, because of the set_S_as_I'
       endif   
        dhij(:,:,:,:)  = ham_tb0(:,:,:,:)
       ddhij(:,:,:,:,:)=dham_tb0(:,:,:,:,:)
!
       if ( matrix_generation_and_stop ) then
         write(*,*)'INFO:matrix_generation_and_stop is selected'
         call output_levels_matrices
!           ---> Plot the levels and matrices (dhij, dsij)
         call output_eigenstates('output_basis_information.txt', .true.)
!        call test_for_quantum_dynamics ! ADDED FOR EXPERIMENTAL QD CALCULATION (2015.1.17)
         write(*,'(a)')'.... ELSES ended in the matrix_generation mode.'
         stop
       endif   
!
       if( n_csc_loop ==0 ) then
         call qm_solver_geno_main ! input dhij(Hamiltonian), dsij(overlap)
         call mulliken            ! calculate new e_num_on_atom, e_num_on_basis
       endif  
!
       if (n_csc_loop > 0) then
         select case (CSC_METHOD)
           case("ELSTNER")
              call qm_engine_csc
           case("TB0")
              call qm_engine_csc_tb0
           case default
              stop "elses_qm_engine: !!!ERROR!!! unexpected CSC_METHOD"
         end select     
       endif  
!
       call save_population_after_convergence
!        -------> set config%system%structure%vatom(:)%population
!                 set config%system%structure%vatom(:)%population_guess
!
!!!!!!!! SY Oct 11, 2008
!!!!!!!! Force Calculation
       call set_qm_force_geno                       ! atm_force = F_{REP} + F_{TB}
       if(CSC_METHOD == "ELSTNER" .and. n_csc_loop>0) then
          call set_atm_force_csc
          atm_force(:,:) = atm_force(:,:) + atm_force_csc(:,:) ! (F_{REP} + F_{TB}) + F_{CSC}
          if(i_verbose > 0) then
             write(*,'("max atm_force_csc(1,:) = ",ES15.8)') maxval(atm_force_csc(1,:))
             write(*,'("max atm_force_csc(2,:) = ",ES15.8)') maxval(atm_force_csc(2,:))
             write(*,'("max atm_force_csc(3,:) = ",ES15.8)') maxval(atm_force_csc(3,:))
          end if
       end if
       call generate_foi
       call qm_geno_output
!      call output_levels_matrices
    end if
!
     if ((c_system == 'C_Xu') .or. (c_system == 'Si_Kwon') .or. (c_system == 'GaAs.Molteni.1993')) then
       ierr=0
       call elses_es_rest
       call elses_set_hami
       if ( matrix_generation_and_stop ) then
         write(*,*)'INFO:matrix_generation_and_stop is selected'
         call output_levels_matrices
!           ---> Plot the levels and matrices (dhij, dsij)
         write(*,'(a)')'.... ELSES ended in the matrix_generation mode.'
         stop
       endif   
       call qm_solver_ortho_main 
       if (population_in_ortho_cases) call save_population_ortho
       call elses_qm_force
       imode_icohp=1
       call qm_output_icohp(imode_icohp)
       if (itemd .eq. itemdmx) then
          call elses_plot_hami
       end if
     end if   
!
     if (c_system == 'NRL_leg') then
       ierr=0
       call qm_engine_nrl
!      call qm_nrl_set_ham
!      call qm_solver_geno_main ! input dhij(Hamiltonian), dsij(overlap)
!      call qm_nrl_set_force
!      call qm_geno_output
     endif   
!
     if (ierr == 1) then
       write(*,*)'ERROR(elses_qm_engine)'
       write(*,*)'ierr=',ierr
       write(*,*)'STOP'
       stop
     end if   
!
   end subroutine elses_qm_engine
!
end module M_qm_engine




