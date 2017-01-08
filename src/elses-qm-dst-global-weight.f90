!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_global_weight
!
!
  implicit none
!  
  private
!
! Public routines
  public set_global_weight_mat
  public set_chem_pot_by_global_weight_mat
  public plot_ldos_by_global_weight_mat
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Store the global weight matrix
!
  subroutine set_global_weight_mat(dst_atm_list, len_dst_atm_list, & 
&                                       kr_dim_dst, wt_kr_dst, eig_kr_dst)
!
    use M_la_krgl_main,     only : eig_kr, wt_kr !(CHANGED)
    use M_qm_domain,         only : i_verbose, noav, atm_element, nval !(unchanged)
    use M_lib_mpi_wrapper,   only : mpi_wrapper_allreduce_r3 !(routine)
!
    implicit none
    integer,          intent(in)  :: dst_atm_list(:)
    integer,          intent(in)  :: len_dst_atm_list(:)
    integer,          intent(in)  :: kr_dim_dst(:,:)
    real(8),          intent(in)  :: wt_kr_dst(:,:,:)
    real(8),          intent(in)  :: eig_kr_dst(:,:,:)
!
    integer :: dst_atm_index, orb_index, atm_index
    integer :: imode
!
    imode=1
    call alloc_global_weight_mat(imode)
!
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_gkrylov_mpi_global_weight_mat'
    endif
!  
    wt_kr(:,:,:)=0.0d0
    eig_kr(:,:,:)=0.0d0
!
    do dst_atm_index=1,len_dst_atm_list(1)
      atm_index=dst_atm_list(dst_atm_index)
      do orb_index=1, nval(atm_element(atm_index))
        wt_kr (:, orb_index, atm_index) = wt_kr_dst (:, orb_index, dst_atm_index)
        eig_kr(:, orb_index, atm_index) = eig_kr_dst(:, orb_index, dst_atm_index)
      enddo   
    enddo  
!
    call mpi_wrapper_allreduce_r3(wt_kr)
    call mpi_wrapper_allreduce_r3(eig_kr)
!
  end subroutine set_global_weight_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the chemical potential by the global matrices (wt_kr, eig_kr)
!
  subroutine set_chem_pot_by_global_weight_mat(rNelec)
!
    use M_la_krgl_main,      only : set_chemical_potential !(routine)
    implicit none
    real(8), intent(in) :: rNelec
!
    call set_chemical_potential(rNelec)
    return
!
  end subroutine set_chem_pot_by_global_weight_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Allocate or deallocate the global weight matices (wt_kr, eig_kr)
!
  subroutine alloc_global_weight_mat(imode)
!
    use M_la_krgl_main,     only : eig_kr, wt_kr !(CHANGED)
    use M_qm_domain, only : i_verbose, noav, nval !(unchanged)
    use M_config,  only : config !(unchaged)
    implicit none
    integer,          intent(in)  :: imode
    integer                       :: nval_max, ierr
    integer                       :: kr_dim_max_input
!
    kr_dim_max_input=config%calc%solver%dimension
!
    if (i_verbose >= 1) then
      write(*,*)'@@ alloc_global_weight_mat:imode=',imode
      write(*,*)'   kr_dim_max_input=',kr_dim_max_input
    endif
!   
    nval_max=maxval(nval)
    if (i_verbose >= 1) then
      write(*,*)'   noav    =',noav
      write(*,*)'   nval_max=',nval_max
    endif
!
    if (imode == 1) then
!
      if (.not. allocated(wt_kr)) then 
        allocate (wt_kr(kr_dim_max_input,nval_max,noav),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (wt_kr)'
      endif  
      wt_kr(:,:,:)=0.0d0
!
      if (.not. allocated(eig_kr)) then 
        allocate (eig_kr(kr_dim_max_input,nval_max,noav),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (eig_kr)'
      endif  
      eig_kr(:,:,:)=0.0d0
!
    endif  
!
    if (imode == 2) then
      deallocate(wt_kr,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
      deallocate (eig_kr,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
    endif  
!
  end subroutine alloc_global_weight_mat
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot LDOS by the global weight matrices
!
  subroutine plot_ldos_by_global_weight_mat
!
    use M_la_krgl_main,     only : plot_ldos_in_krgl !(routine)
    implicit none
!
    call plot_ldos_in_krgl
!
  end subroutine plot_ldos_by_global_weight_mat
!    
end module M_qm_dst_global_weight
