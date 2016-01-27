!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_rest_kernel
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
!
! Public routines
  public :: calc_geno_rest_kernel
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Non order-N rest-part calculation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine calc_geno_rest_kernel(atm_index_seed, booking_list, len_booking_list, atm_rest_energy, atm_vdW_energy, atm_force)
!
    use M_qm_geno_Huckel,  only : RepulsiveEnergy !(routine)
    use M_md_get_distance, only : get_vector_for_atom_pair !(routine)
!
    implicit none
    integer, intent(in)     :: atm_index_seed
    integer, intent(in)     :: booking_list(:)
    integer, intent(in)     :: len_booking_list
    real(DOUBLE_PRECISION), intent(out)    :: atm_rest_energy
    real(DOUBLE_PRECISION), intent(out)    :: atm_vdW_energy
    real(DOUBLE_PRECISION), intent(inout), optional  :: atm_force(:,:)
!
    integer            :: jsv1, jsv2, jsk1, jsk2, jsd1
    integer            :: atm_index1, atm_index2
    real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz, ddsum, ddsum2
    real(DOUBLE_PRECISION) :: repulsive_force(3)
    real(DOUBLE_PRECISION) :: E_atompair_vdW_only
!
    jsv2=atm_index_seed
!
    ddsum=0.0d0
    ddsum2=0.0d0
!
    do jsd1=1, len_booking_list
      jsv1=booking_list(jsd1) 
!
      if (jsv1 == jsv2) cycle
!
      call get_vector_for_atom_pair(jsv1, jsv2, dvecx, dvecy, dvecz)
!          ----> (dvecx, dvecy, dvecz) : Vector r = R_2 - R_1 
!
      ddsum  = ddsum  + RepulsiveEnergy(jsv1, jsv2, dvecx, dvecy, dvecz, repulsive_force, E_atompair_vdW_only)
      ddsum2 = ddsum2 + E_atompair_vdW_only
      if (present(atm_force)) then
        atm_force(:,jsv2)=atm_force(:,jsv2)-repulsive_force(:)
        atm_force(:,jsv1)=atm_force(:,jsv1)+repulsive_force(:)
      endif  
    enddo
    atm_rest_energy = ddsum
    atm_vdW_energy  = ddsum2
!
  end subroutine calc_geno_rest_kernel
!
end module M_qm_geno_rest_kernel



