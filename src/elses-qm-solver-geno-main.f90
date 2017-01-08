!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_solver_geno_main
!
  implicit none
!  
  private
!
! Public routines
  public qm_solver_geno_main
  public qm_solver_ortho_main
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine qm_solver_geno_main
!
    use M_qm_domain,          only : i_verbose          !(unchanged)
    use M_config,             only : config             !(CHANGED:config%calc%solver%scheme)
    use M_qm_solver_eig_geno, only : qm_solver_eig_geno !(routine)
    use M_qm_solver_gscocg,   only : qm_solver_gscocg   !(routine)
    use M_qm_solver_krgl,     only : qm_solver_krgl     !(routine)
    use M_la_krgl_main  ,     only : reset_dsij         !(routine)
    implicit none
    character(len=32) :: scheme_mode
!    
!   if (config%calc%solver%scheme == 'GSCOCG') then
!     config%calc%solver%scheme='gscocg'
!   endif   
!
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_solver_geno_main'
    endif   
!
!   call reset_dsij
!
    scheme_mode=trim(config%calc%solver%scheme)
    if (scheme_mode == 'scheme_default') then 
      if (i_verbose >= 1) write(*,*)' default solver is chosen'
      scheme_mode='eigen'
    endif   
!
    if (scheme_mode == 'gKrylov')   scheme_mode='ekrgl'
    if (scheme_mode == 'gKrylov_A') scheme_mode='ekrgl'
    if (scheme_mode == 'gkrylov_a') scheme_mode='ekrgl'
    if (scheme_mode == 'GSCOCG')  scheme_mode='gscocg'
!
    config%calc%solver%scheme=scheme_mode
!
    if (i_verbose >= 1) then
      write(*,*)' solver     =',config%calc%solver%scheme
!     write(*,*)' projection =',config%calc%solver%projection
!     write(*,*)' dimension  =',config%calc%solver%dimension
    endif   
!
    select case (trim(config%calc%solver%scheme))
     case ('gscocg')
       call qm_solver_gscocg
     case ('krgl', 'ekrgl' )
       call qm_solver_krgl(scheme_mode)
     case ('eigen', 'eigen_mpi')
       call qm_solver_eig_geno
     case default
      write(*,*) 'ERROR: No solver is chosen'
      stop
    end select    
!
!
  end subroutine qm_solver_geno_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine qm_solver_ortho_main
!
    use M_qm_domain,          only : i_verbose          !(unchanged)
    use M_config,             only : config             !(unchanged)
    use elses_mod_md_dat,     only : itemd              !(unchanged)
    implicit none
    character(len=32) :: scheme_mode
    integer :: i_eig_main_init
    integer :: imode
!    
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_solver_ortho_main'
    endif   
!
    scheme_mode=trim(config%calc%solver%scheme)
    if (scheme_mode == 'scheme_default') then 
      if (i_verbose >= 1) write(*,*)' default solver is chosen'
      scheme_mode='krylov'
    endif   
!
    if (scheme_mode == 'Krylov subspace') scheme_mode='krylov'
    if (scheme_mode == 'krylov subspace') scheme_mode='krylov'
!
    config%calc%solver%scheme=scheme_mode
!
    if (i_verbose >= 1) then
      write(*,*)' solver     =',config%calc%solver%scheme
    endif   
!
    select case (scheme_mode)
     case ('eigen')
       i_eig_main_init=0
       if (itemd == 1) i_eig_main_init=1
       imode=0
       call elses_eig_main(imode, i_eig_main_init)
     case ('krylov')
       call elses_kr_main
     case default
      write(*,*) 'ERROR: No solver is chosen'
      stop
   end select    
!
!
  end subroutine qm_solver_ortho_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_solver_geno_main

