!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_matvec_dst
!
  use M_qm_domain
  implicit none
!
  private
  public :: matvec_mul_dst_r
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with H (REAL VARIABLE)
!   mat_kind='H' : (vect_out) = H (vect_in)
!   mat_kind='S' : (vect_out) = S (vect_in)
!
  subroutine matvec_mul_dst_r(vect_in,vect_out,mat_kind,jsv4jsk,jsk4jsv,jjkset)
!
    use M_qm_domain,   only : noav, nval, njsd, jsv4jsd, atm_element  ! (unchanged)
    use M_qm_dst_proj_cell,  only : dst_atm_list_rev, len_dst_atm_list      ! (unchanged)
    use M_qm_geno_dst, only : ham_tot_dst, overlap_dst                ! (unchanged)
    implicit none
    real(8),          intent(in)  :: vect_in(:)
    real(8),          intent(out) :: vect_out(:)
    integer,          intent(in)  :: jsv4jsk(:)
    integer,          intent(in)  :: jsk4jsv(:)
    integer,          intent(in)  :: jjkset(:)
    character(len=*),  intent(in)  :: mat_kind
    integer :: m, noak, noav2
    integer :: ierr, ict4h
    integer :: jsk2, jsv2, jjkset2, jsd1, jsv1, jsk1 
    integer :: jjkset1, ja2, jjk2, ja1, jjk1, dst_atm_index
    real(8) :: dbigd
!
    ict4h=1
    m=size(vect_in,1)
    noav2=size(jsk4jsv,1)
    noak=size(jsv4jsk,1)
!
    if (.not. allocated(dst_atm_list_rev)) then
      stop 'ERROR(matvec_mul_dst_r):dst_atm_list_rev is not allocated'
    endif
!
    if (m /= size(vect_out,1)) stop 'Abort:Size mismatch of m'
    if (noav2 /= noav) stop 'Abort:Size mismatch of noav'
    if (noak /= size(jjkset,1)) stop 'Abort:Size mismatch of noak'
    if ((mat_kind /='H') .and. (mat_kind /='S')) stop 'Abort:Wrong mode'
!
    vect_out(:)=0.0d0
!
    if (mat_kind == 'H') then
      do jsk2=1,noak
        jsv2=jsv4jsk(jsk2)
        dst_atm_index=dst_atm_list_rev(jsv2)
        if ((dst_atm_index <=0) .or. (dst_atm_index > len_dst_atm_list(2))) then
          write(*,*)'ERROR(matvec_mul_proj):dst_atm_index=',dst_atm_index
          stop
        endif   
        jjkset2=jjkset(jsk2)
        do jsd1=1,njsd(jsv2,ict4h)
          jsv1=jsv4jsd(jsd1,jsv2)
          jsk1=jsk4jsv(jsv1)
          if ((jsk1 < 0) .or. (jsk1 > noak)) then
            write(*,*)'ERROR(matvec_mul_proj)'
            write(*,*)'jsv1, jsk1=',jsv1,jsk1
            stop
          endif   
          if (jsk1 .ne. 0) then
            jjkset1=jjkset(jsk1)
            do ja2=1,nval(atm_element(jsv2))
              jjk2=jjkset2+ja2
              do ja1=1,nval(atm_element(jsv1))
                jjk1=jjkset1+ja1
                dbigd=ham_tot_dst(ja1,ja2,jsd1,dst_atm_index)
                vect_out(jjk2)=vect_out(jjk2)+dbigd*vect_in(jjk1)
              enddo
            enddo
          endif  
        enddo   
      enddo   
    endif
!
    if (mat_kind == 'S') then
      do jsk2=1,noak
        jsv2=jsv4jsk(jsk2)
        dst_atm_index=dst_atm_list_rev(jsv2)
        if ((dst_atm_index <=0) .or. (dst_atm_index > len_dst_atm_list(2))) then
          write(*,*)'ERROR(matvec_mul_proj):dst_atm_index=',dst_atm_index
          stop
        endif   
        jjkset2=jjkset(jsk2)
        do jsd1=1,njsd(jsv2,ict4h)
          jsv1=jsv4jsd(jsd1,jsv2)
          jsk1=jsk4jsv(jsv1)
          if ((jsk1 < 0) .or. (jsk1 > noak)) then
            write(*,*)'ERROR(matvec_mul_proj)'
            write(*,*)'jsv1, jsk1=',jsv1,jsk1
            stop
          endif   
          if (jsk1 .ne. 0) then
            jjkset1=jjkset(jsk1)
            do ja2=1,nval(atm_element(jsv2))
              jjk2=jjkset2+ja2
              do ja1=1,nval(atm_element(jsv1))
                jjk1=jjkset1+ja1
                dbigd=overlap_dst(ja1,ja2,jsd1,dst_atm_index)
                vect_out(jjk2)=vect_out(jjk2)+dbigd*vect_in(jjk1)
              enddo
            enddo
          endif  
        enddo   
      enddo   
    endif
!
  end subroutine matvec_mul_dst_r
!
!
end module M_la_matvec_dst

