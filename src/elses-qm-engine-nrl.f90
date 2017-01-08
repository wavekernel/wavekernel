!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_engine_nrl
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   implicit none
!
   private
   public qm_engine_nrl
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_engine_nrl
!
!
     use elses_mod_md_dat,      only : itemd
!
!
!    Module procedures
     use M_config,              only : config
     use M_qm_solver_geno_main, only : qm_solver_geno_main 
     use M_qm_geno_output,      only : qm_geno_output
     use M_qm_nrl,              only : qm_nrl_set_ham, qm_nrl_set_force
     use M_qm_output_matrices,  only : output_levels_matrices !(subroutine)
!
     implicit none
     real(DOUBLE_PRECISION):: dq, dq_1, x, mix_r
     integer :: ierr, i_csc_loop, jsd1, jsv2, ja1, ja2, n_csc_loop
     character(len=8) :: CSC_METHOD
!
     CSC_METHOD = config%calc%genoOption%CSC_method
     n_csc_loop = config%calc%genoOption%CSC_max_loop_count
!
     if (i_verbose >= 1) then
       write(*,*)'@@ elses_qm_engine_nrl:itemd=',itemd
     endif   
!
     call qm_nrl_set_ham
     call qm_solver_geno_main ! input dhij(Hamiltonian), dsij(overlap)
     call qm_nrl_set_force
     call qm_geno_output
     call output_levels_matrices
!
   end subroutine qm_engine_nrl
!
end module M_qm_engine_nrl





