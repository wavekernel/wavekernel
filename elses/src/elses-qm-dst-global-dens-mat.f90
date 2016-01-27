!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_global_dens_mat
!
!
  implicit none
!  
  private
!
! Public routines
  public set_global_density_mat
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Store the global density matrix
!
  subroutine set_global_density_mat(dst_atm_list, len_dst_atm_list, & 
&                                    dens_mat_dst, e_dens_mat_dst)
!
    use M_qm_domain,         only : dbij, dpij !(CHANGED)
    use M_qm_domain,         only : dhij       !(unchanged)
    use M_qm_domain,         only : i_verbose, noav, nval !(unchanged)
    use M_lib_mpi_wrapper,   only : mpi_wrapper_allreduce_r4 !(routine)
    use elses_mod_sel_sys,   only : r_cut_book !(unchanged)
    use M_la_krgl_main,      only : sym_dbij_and_dpij  ! (routine)
    use M_qm_domain,         only : get_length_of_list !(routine)
!
    implicit none
    integer,          intent(in)  :: dst_atm_list(:)
    integer,          intent(in)  :: len_dst_atm_list(:)
    real(8),          intent(in)  ::   dens_mat_dst(:,:,:,:)
    real(8),          intent(in)  :: e_dens_mat_dst(:,:,:,:)
!
    integer :: dst_atm_index, orb_index, atm_index
    real(8) :: cutoff_radius
    integer :: length_of_list
    integer :: alloc_mat_size, mat_size1, mat_size2, mat_size3
    integer :: ierr
!
    if (i_verbose >= 1) then
      write(*,*)'@@ set_global_density_mat'
    endif
!  
    cutoff_radius=r_cut_book
!   call get_length_of_list(cutoff_radius,length_of_list)
    alloc_mat_size=size(dhij,3)
!   alloc_mat_size=0
!   alloc_mat_size=minval(length_of_list*2, noav)
!      ----> two is a tolarence factor
!
    mat_size1=size(dens_mat_dst,1)
    mat_size2=size(dens_mat_dst,2)
    mat_size3=size(dens_mat_dst,3)
!
    if (alloc_mat_size < mat_size3) then
      write(*,*)'ERROR in set_global_density_mat'
      write(*,*)' alloc_mat_size =',alloc_mat_size
      write(*,*)' mat_size3      =',mat_size3
      stop
    endif
!   
    if (allocated(dbij)) then
      deallocate(dbij, stat=ierr)
      if ( ierr /= 0 ) stop 'ERROR in dealloc. dbij'
    endif   
!
    if (allocated(dpij)) then
      deallocate(dpij, stat=ierr)
      if ( ierr /= 0 ) stop 'ERROR in dealloc. dpij'
    endif   
!
    allocate (dbij(mat_size1, mat_size2, alloc_mat_size, noav), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (dbij)'
    dbij(:,:,:,:)=0.0d0
!
    allocate (dpij(mat_size1, mat_size2, alloc_mat_size, noav), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (dpij)'
    dpij(:,:,:,:)=0.0d0
!
!$omp  parallel default(shared) &
!$omp& private (dst_atm_index, atm_index)
!$omp  do schedule(static)
    do dst_atm_index=1,len_dst_atm_list(1)
      atm_index=dst_atm_list(dst_atm_index)
      dbij(:, :, 1:mat_size3, atm_index) =   dens_mat_dst(:, :, 1:mat_size3, dst_atm_index)
      dpij(:, :, 1:mat_size3, atm_index) = e_dens_mat_dst(:, :, 1:mat_size3, dst_atm_index)
    enddo  
!$omp end do
!$omp end parallel
!
    call mpi_wrapper_allreduce_r4(dbij)
    call mpi_wrapper_allreduce_r4(dpij)
!
    call sym_dbij_and_dpij
!
  end subroutine set_global_density_mat
!
!    
end module M_qm_dst_global_dens_mat
