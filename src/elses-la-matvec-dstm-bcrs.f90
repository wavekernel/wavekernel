!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_matvec_bcrs
!
  use M_config  !(unchanged)
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
  public :: switch_to_bcrs
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with DSTM matrix : A (= H or S)
!      (vect_out): = A (vect_in)
!
  subroutine switch_to_bcrs(vect_in, vect_out, jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, mat_dstm)
!
    implicit none
    real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
    real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: jsv4jsk(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_dstm(:,:,:,:)
    integer :: n_block
    logical, parameter                  :: use_bcrs_oc= .false.
!
    if (config%calc%distributed%mat_vec_const_num_orbital) then
      n_block=config%calc%distributed%mat_vec_max_num_orbital
      if (.not. use_bcrs_oc) then
        call matvec_mul_bcrs(vect_in, vect_out, jjkset, booking_list_dstm, booking_list_dstm_len, mat_dstm, n_block) 
      else
        call matvec_mul_bcrs_oc(vect_in, vect_out, jjkset, booking_list_dstm, booking_list_dstm_len, mat_dstm, n_block) 
      endif
    else
      stop 'STOP switch_to_bcrs'
    endif
!
  end subroutine switch_to_bcrs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with DSTM matrix : A (= H or S)
!      (vect_out): = A (vect_in)
!
  subroutine matvec_mul_bcrs(vect_in, vect_out, jjkset, booking_list_dstm, booking_list_dstm_len, mat_dstm, nb)
!
    implicit none
    real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
    real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_dstm(:,:,:,:)
    integer,                intent(in)  :: nb ! block size
!
    integer :: num_atom_proj, n
    integer :: jsd, jsk1, jsk2, ja1, ja2
    integer :: jjk1, jjk2, jjkset1, jjkset2
    integer :: i_sta, i_end, j_sta, j_end
    real(DOUBLE_PRECISION)  :: mat_value
    real(DOUBLE_PRECISION), parameter :: alpha=1.0d0
    real(DOUBLE_PRECISION), parameter :: beta=0.0d0
    character(len=*), parameter :: mat_vec_mode='blas'
!   character(len=*), parameter :: mat_vec_mode='matmul'
!   character(len=*), parameter :: mat_vec_mode='simple'
!   character(len=*), parameter :: mat_vec_mode='loop_order_changed'

!
    num_atom_proj=size(mat_dstm, 4)
    vect_out(:)=0.0d0
!
    if (mat_vec_mode == 'blas') then
      do jsk2=1, num_atom_proj
        i_sta=(jsk2-1)*nb+1
        i_end=i_sta-1+nb
        do jsd=1, booking_list_dstm_len(jsk2)
         jsk1=booking_list_dstm(jsd, jsk2)
         j_sta=jjkset(jsk1)+1
         j_end=j_sta-1+nb
         call dgemv('N', nb, nb, alpha, mat_dstm(:,:,jsd,jsk2), nb, vect_in(j_sta:j_end), 1, beta, vect_out(i_sta:i_end),1)
        enddo
      enddo
      return
    endif
!
    if (mat_vec_mode == 'matmul') then
      do jsk2=1, num_atom_proj
        i_sta=(jsk2-1)*nb+1
        i_end=i_sta-1+nb
        do jsd=1, booking_list_dstm_len(jsk2)
         jsk1=booking_list_dstm(jsd, jsk2)
         j_sta=jjkset(jsk1)+1
         j_end=j_sta-1+nb
         vect_out(i_sta:i_end)=vect_out(i_sta:i_end) &
&            +matmul(mat_dstm(:,:,jsd,jsk2),vect_in(j_sta:j_end))
        enddo
      enddo
      return
    endif
!
    if (mat_vec_mode == 'simple') then
      do jsk2=1, num_atom_proj
        jjkset2=(jsk2-1)*nb+1
        do jsd=1, booking_list_dstm_len(jsk2)
          jsk1=booking_list_dstm(jsd, jsk2)
          jjkset1=jjkset(jsk1)
          do ja2=1,nb
            jjk2=jjkset2+ja2
            do ja1=1,nb
              jjk1=jjkset1+ja1
              mat_value=mat_dstm(ja1,ja2,jsd,jsk2)
              vect_out(jjk2)=vect_out(jjk2)+mat_value*vect_in(jjk1)
            enddo
          enddo
        enddo
      enddo
      return
    endif
!
    if (mat_vec_mode == 'loop_order_changed') then
      do jsk2=1, num_atom_proj
        jjkset2=(jsk2-1)*nb+1
        do ja2=1,nb
          jjk2=jjkset2+ja2
          do ja1=1,nb
            do jsd=1, booking_list_dstm_len(jsk2)
              jsk1=booking_list_dstm(jsd, jsk2)
              jjkset1=jjkset(jsk1)
              jjk1=jjkset1+ja1
              mat_value=mat_dstm(ja1,ja2,jsd,jsk2)
              vect_out(jjk2)=vect_out(jjk2)+mat_value*vect_in(jjk1)
            enddo
          enddo
        enddo
      enddo
      return
    endif
!
  end subroutine matvec_mul_bcrs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with DSTM matrix : A (= H or S)
!      (vect_out): = A (vect_in)
!
  subroutine matvec_mul_bcrs_oc(vect_in, vect_out, jjkset, booking_list_dstm, booking_list_dstm_len, mat_dstm, nb)
!
    implicit none
    real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
    real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_dstm(:,:,:,:)
    integer,                intent(in)  :: nb ! block size
!
    real(DOUBLE_PRECISION), allocatable :: mat_dstm_oc(:,:,:,:)
!
    integer :: num_atom_proj
    integer :: jsk2, ja2, ja1, ierr
!
    num_atom_proj=size(mat_dstm, 4)
    allocate (mat_dstm_oc(num_atom_proj, nb, nb, num_atom_proj), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error. mat_dstm_oc'
    endif
!
    do jsk2=1, num_atom_proj
      do ja2=1,nb
        do ja1=1,nb
          mat_dstm_oc(1:booking_list_dstm_len(jsk2),ja1,ja2,jsk2) =   mat_dstm(ja1,ja2,1:booking_list_dstm_len(jsk2),jsk2)
        enddo
      enddo
    enddo
!
    call matvec_mul_bcrs_oc_core(vect_in, vect_out, jjkset, booking_list_dstm, booking_list_dstm_len, mat_dstm_oc, nb)
!
    deallocate(mat_dstm_oc, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Delloc. error. mat_dstm_oc'
    endif

  end subroutine matvec_mul_bcrs_oc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with DSTM matrix : A (= H or S)
!      (vect_out): = A (vect_in)
!
  subroutine matvec_mul_bcrs_oc_core(vect_in, vect_out, jjkset, booking_list_dstm, booking_list_dstm_len, mat_dstm_oc, nb)
!
    implicit none
    real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
    real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_dstm_oc(:,:,:,:)
    integer,                intent(in)  :: nb ! block size
!
    integer :: num_atom_proj, n
    integer :: jsd, jsk1, jsk2, ja1, ja2
    integer :: jjk1, jjk2, jjkset1, jjkset2
    real(DOUBLE_PRECISION)  :: s
    character(len=*), parameter :: mat_vec_mode ='simple'
!
    num_atom_proj=size(mat_dstm_oc, 4)
    vect_out(:)=0.0d0
!
    if (mat_vec_mode == 'simple') then
      do jsk2=1, num_atom_proj
        jjkset2=(jsk2-1)*nb+1
        do ja2=1,nb
          jjk2=jjkset2+ja2
          s=0.0d0
          do ja1=1,nb
            do jsd=1, booking_list_dstm_len(jsk2)
              jsk1=booking_list_dstm(jsd, jsk2)
              jjkset1=jjkset(jsk1)
              jjk1=jjkset1+ja1
              s=s+mat_dstm_oc(jsd,ja1,ja2,jsk2)*vect_in(jjk1)
            enddo
          enddo
          vect_out(jjk2)=s
        enddo
      enddo
      return
    endif
!
  end subroutine matvec_mul_bcrs_oc_core
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_la_matvec_bcrs

