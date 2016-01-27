!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_dstm_get_nna_distance
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_wall_clock_time,  only : get_system_clock_time !(routine)
   implicit none
   logical :: global_dens_mat
   logical :: global_ham_mat 
   logical :: overlap_is_on  
!
   private
!
! Public module variable
   public  :: get_nna_distance_dstm
!  public  :: global_ham_mat 
!  public  :: overlap_is_on  
!
! Public routine
!  public :: qm_domain_setting_dst
!  public :: proj_init_end_dst
!  public :: set_projection_dst_cell
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_nna_distance_dstm(atm_index, orb_index, jjkset, jsv4jsk, & 
&                                  booking_list_dstm, booking_list_dstm_len, mat_dstm, atm_index2, nna_distance)
!
    use M_md_get_distance,         only : get_distance !(routine)
    implicit none
    integer,                intent(in)  :: atm_index, orb_index
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: jsv4jsk(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_dstm(:,:,:,:)
    integer,                intent(out) :: atm_index2
    real(DOUBLE_PRECISION), intent(out) :: nna_distance
!
    logical, parameter :: debug_mode = .true.
!
    integer :: jsk2, jsk1, jsd, jsv2, jsv1
    integer :: num_atom_proj
    integer :: jjkset1, jjkset2, jjk1, jjk2
!
    real(DOUBLE_PRECISION)  :: w
    integer                 :: jsv_min
    real(DOUBLE_PRECISION)  :: dist_min
!
!
    num_atom_proj=size(mat_dstm, 4)
!
    dist_min=1.0d100
    jsv_min =0
    do jsk2=1, num_atom_proj
      jjkset2=jjkset(jsk2)
      jsv2=jsv4jsk(jsk2)
      if (jsv2 /= atm_index) cycle
      do jsd=1, booking_list_dstm_len(jsk2)
        jsk1=booking_list_dstm(jsd, jsk2)
        if (debug_mode) then
          if ( ( jsk1 <= 0 ) .or. (jsk1 > num_atom_proj)) then
            stop 'ERROR(matvec_mul_dstm)'
          endif
        endif
        jsv1=jsv4jsk(jsk1)
        if (jsv1 == jsv2) cycle
        call get_distance(jsv1, jsv2, w) 
        if (w < dist_min) then
          jsv_min =jsv1   
          dist_min=w
        endif   
      enddo
    enddo
!
    atm_index2   = jsv_min
    nna_distance = w
!
  end subroutine get_nna_distance_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_dstm_get_nna_distance
