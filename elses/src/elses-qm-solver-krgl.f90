!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_solver_krgl
!
!
  implicit none
!  
  private
!
! Public routines
  public qm_solver_krgl
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine qm_solver_krgl(scheme_mode)
!
    use elses_mod_orb2,     only : js2j ! (unchanged)
    use M_qm_domain,        only : i_verbose, noav, atm_element, nval, &
&                                    total_electron_number, chemical_potential ! (unchanged)
    use elses_mod_noav,     only : noav ! (unchanged)
    use M_qm_projection,    only : proj_init_end, proj_get_mat_size, &
&                                  proj_get_list, convert_dm_and_edm ! (routines)
    use M_la_krgl_main,     only : krgl_alloc, krgl_main, set_chemical_potential, &
&                                  plot_ldos_in_krgl,  reset_dsij, sym_dbij_and_dpij ! (routines)
    use M_la_gkrylov_hist,  only : check_proj_list, alloc_for_hist, set_s_inv_e_j_for_hist ! (routines)
    use elses_mod_md_dat,    only : itemd, itemdmx
    use M_la_gkrylov_output, only : calc_eigen_in_gkrylov
!
    implicit none
    character(len=32), intent(in) :: scheme_mode
!
    logical, parameter :: flag_for_sym_dbij_and_dpij = .false.
!
    integer :: atm_index, orb_index, mat_size, num_energy_points, j_src
    integer :: num_atom_proj
    integer :: nss2, nval2, m, n, ierr, i, k
    integer :: prc_index
    integer :: imode
    integer :: kr_dim_max
    integer :: i_show
!
    real(8) :: elec_num, Nelec
!
!   real(8) :: sigma
!   real(8) :: elec_num,miu,miu1,miu2,bieps,Nelec,Eelec,nf,nf1,nf2,Fermi,tau,dk,temp !mended by Teng
!
    integer, allocatable :: jsv4jsk(:)
    integer, allocatable :: jsk4jsv(:)
    integer, allocatable :: jjkset(:)
    real(8), allocatable :: b(:)
!
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_solver_krgl:scheme mode=',trim(scheme_mode)
    endif  
!
!
    imode=1
    call proj_init_end(imode)
!
    call check_proj_list
    call alloc_for_hist
!
    imode=1
    call set_s_inv_e_j_for_hist(imode)
!
    imode=1
    call krgl_alloc(imode)
!    
    do prc_index=1,2
!
      if (i_verbose >= 1) then
        write(*,*)' prc_index = ',prc_index
      endif  
!
      Nelec=0.0d0 
!
!$omp  parallel default(shared) &
!$omp& private (atm_index, orb_index, i_show) &
!$omp& private (elec_num, mat_size, num_atom_proj, nss2, nval2, ierr, j_src) &
!$omp& private (jsv4jsk, jsk4jsv, jjkset, b) &
!$omp& reduction (+ : Nelec)
!$omp  do schedule(static)
      do atm_index=1,noav
!
        i_show=0
        if (atm_index <= 2) then
          i_show=1
        endif
!   
        call get_elec_num_for_atom(atm_index, elec_num)
        Nelec=Nelec+elec_num ! summing up of electron numbers
!
        call proj_get_mat_size(atm_index, mat_size, num_atom_proj)
        if (i_show >= 1) then
          write(*,*)'atm_index     =',atm_index
          write(*,*)'mat_size      =',mat_size
          write(*,*)'num_atom_proj =',num_atom_proj
        endif  
        nss2=atm_element(atm_index)
        nval2=nval(nss2)
!
        if (i_show >= 1) then
          write(*,*)'nss2          =',nss2
          write(*,*)'nval2         =',nval2
        endif  
!
        allocate (jsv4jsk(num_atom_proj),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        allocate (jsk4jsv(noav),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        allocate (jjkset(num_atom_proj),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        allocate (b(mat_size),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
        call proj_get_list(atm_index, jsv4jsk, jsk4jsv, jjkset)
!
!       call get_interac_list_num_proj_atom(atm_index,m_int)
!           ---> get number of atom in the interaction list : m_int
!
        do orb_index=1,nval2
          if (i_show >= 1) then
            write(*,*)"atom=",atm_index,"orbit=",orb_index 
          endif  
          j_src=js2j(orb_index,atm_index)
          b(:)=0.0d0
          b(orb_index)=1.0d0
          call krgl_main(j_src,b,jsv4jsk,jsk4jsv,jjkset,prc_index,scheme_mode)
        enddo   
!
        deallocate(jsv4jsk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate(jsk4jsv,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate(jjkset,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate(b, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      enddo ! loop end for atm_index 
!
!$omp end do
!$omp end parallel
!
!
    if (prc_index == 1) then ! (calculation of chemical potential)
      write(*,*)'Nelec=',Nelec
      if ( dabs( total_electron_number - Nelec) >= 1.0d-10 ) then
         write(*,*)'ERROR in Nelec : Nelec=',Nelec
         stop
      endif   
      call set_chemical_potential(Nelec)
!     chemical_potential=1.0d20
!     write(*,*)'Chem pot is set to be a huge value (TEST CASE):',chemical_potential
    endif   
!

!
    if ( ( itemd == itemdmx ) .and. (prc_index == 2)) then ! (calculation of LDOS)
      call plot_ldos_in_krgl
!     call calc_eigen_in_gkrylov
    endif   
!
    enddo ! loop end for prc_index
!
    if (i_verbose >= 1) write(*,*)'INFO: flag_for_sym_dbij_and_dpij=',flag_for_sym_dbij_and_dpij
    if (flag_for_sym_dbij_and_dpij) then
      call sym_dbij_and_dpij
    else
      if (i_verbose >= 1) write(*,'(a)')'INFO: sym_dbij_and_dpij IS SKIPPED'
    endif  
!
!   write(*,*)'Stop by T. Hoshi'
!   stop
!
    imode=2
    call krgl_alloc(imode)
!
    imode=2
    call proj_init_end(imode)
!
  end subroutine qm_solver_krgl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ subroutine for obtaining the valence elec. num. for the given atom
!
  subroutine get_elec_num_for_atom(atm_index,elec_num)
!
    use elses_mod_val_elec,     only : val_elec_atm ! (unchanged)
    use M_qm_domain,            only : atm_element  ! (unchanged)
!
    implicit none
    integer,          intent(in)  :: atm_index
    real(8),          intent(out) :: elec_num
    integer :: nss, size1, size2
!
    size1=size(atm_element,1)
    size2=size(val_elec_atm,1)
!
    if ((atm_index <= 0) .or. (atm_index > size1)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'atm_index=',atm_index
      stop
    endif   
!
    nss=atm_element(atm_index)
!
    if ((nss <= 0) .or. (nss > size2)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'nss=',nss
      stop
    endif   
!
    elec_num=val_elec_atm(nss)
!
    if ((elec_num <= 1.0d-10) .or. (elec_num >= 1000.0d0)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'  elec_num=',elec_num
      stop
    endif   
!
  end subroutine get_elec_num_for_atom
!
end module M_qm_solver_krgl

