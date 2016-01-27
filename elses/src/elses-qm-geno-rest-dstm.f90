!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_rest_dstm
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
!
! Public routines
  public :: calc_geno_rest_dstm
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine calc_geno_rest_dstm(atm_index_seed, jsv4jsk, num_atom_proj, &
&                 booking_list_dstm, booking_list_dstm_len, atm_energy_rep, atm_force)
!
    use M_config             ! unchanged
    use M_qm_domain,       only : noav !(unchanged)
    use M_qm_geno_Huckel,  only : RepulsiveEnergy !(routine)
    use M_md_get_distance, only : get_vector_for_atom_pair !(routine)
    use M_qm_geno_rest_kernel, only : calc_geno_rest_kernel !(routine)
!
!
    implicit none
    integer, intent(in)     :: jsv4jsk(:)
    integer, intent(in)     :: num_atom_proj
    integer, intent(in)     :: atm_index_seed
    integer, intent(in)     :: booking_list_dstm(:,:)
    integer, intent(in)     :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(out)    :: atm_energy_rep
    real(DOUBLE_PRECISION), intent(inout), optional  :: atm_force(:,:)
!
    integer            :: ierr
    integer            :: jsv1, jsv2, jsk1, jsk2, jsd1
    integer            :: atm_index1, atm_index2
    real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz, ddsum
!   real(DOUBLE_PRECISION) :: repulsive_force(3)
!
    integer, allocatable   :: local_booking_list(:)
    integer                :: len_local_booking_list
    real(DOUBLE_PRECISION) :: w, w_vdW
!
    jsk2=1
    jsv2=atm_index_seed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (jsv4jsk(jsk2) /= atm_index_seed) then
      write(*,'(a,i20)')'ERROR(cals_geno_rest_dstm):atm_index_seed=',atm_index_seed
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (config%calc%interaction_range%cutoff_rest_cellmax) then
      len_local_booking_list = noav
    else
      len_local_booking_list = booking_list_dstm_len(jsk2)
    endif   
!
    allocate (local_booking_list(len_local_booking_list), stat=ierr)
    if (ierr /= 0) then
      write(*,*)' Alloc. Error in calc_geno_rest_dstm (local_booking_list)'
      stop
    endif   
!
    if (config%calc%interaction_range%cutoff_rest_cellmax) then
      do jsv1=1, noav 
        local_booking_list(jsv1)=jsv1
      enddo  
    else
      do jsd1=1, len_local_booking_list
        jsk1 = booking_list_dstm(jsd1,jsk2)
        jsv1 = jsv4jsk(jsk1)
        local_booking_list(jsd1) = jsv1
      enddo
    endif   
!
    if (present(atm_force)) then
      call calc_geno_rest_kernel(atm_index_seed, local_booking_list, len_local_booking_list, w, w_vdw, atm_force) 
    else
      call calc_geno_rest_kernel(atm_index_seed, local_booking_list, len_local_booking_list, w, w_vdw)
    endif   
!
    deallocate (local_booking_list, stat=ierr)
    if (ierr /= 0) then
      write(*,*)' Dealloc. Error in calc_geno_rest_dstm (local_booking_list)'
      stop
    endif   
!
    atm_energy_rep = w
!
  end subroutine calc_geno_rest_dstm
!
end module M_qm_geno_rest_dstm


