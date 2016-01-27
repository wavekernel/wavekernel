!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_gscocg_tru_res
!
  implicit none
  integer m
  integer nz
!
  real(8), allocatable :: val(:)
  integer, allocatable :: col_ind(:)
  integer, allocatable :: row_ptr(:)
  integer, allocatable :: diag_ptr(:)
!
  real(8), allocatable :: val_s(:)
  integer, allocatable :: col_ind_s(:)
  integer, allocatable :: row_ptr_s(:)
  integer, allocatable :: diag_ptr_s(:)
!
  integer, allocatable :: flag_for_init(:)
!
  private
  public :: calc_tru_res
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_tru_res(b,x,z,log_res)
!    ---> Calculate the tru resisial r = b - ( z S - H ) x
!          output log_res = log_10 ( |r/b| )
!
    use M_la_matvec_io, only : matvec_mul
    implicit none
    complex(8),       intent(in)  :: b(:)
    complex(8),       intent(in)  :: x(:)
    complex(8),       intent(in)  :: z
    real(8),          intent(out) :: log_res
!
    complex(8), allocatable       :: sax(:), sbx(:)
    integer                       :: m, m2, i, ierr
    real(8)                       :: hnor, hes, tbs, tbs2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get matrix size : m
!
    m=size(b,1)
    m2=size(x,1)
    if (m /= m2) stop 'Abort in Tru_Res:Parameter Mismatch:m,m2'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix allocation
!
    allocate (sax(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (sbx(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call matvec_mul(x,sax,'H')
!      ---> Mat-vec multiplication : sax= H x
!
    call matvec_mul(x,sbx,'S')
!      ---> Mat-vec multiplication : sbx= S x
!
    sax(:)=b(:)- (z*sbx(:)-sax(:))
!      ---> r = b - (z S - H) x
!
    hes=0.0d0
    hnor=0.0d0
    do i=1,m
      tbs=abs(sax(i))
      hes=hes+tbs*tbs
      tbs2=abs(b(i))
      hnor=hnor+tbs2*tbs2
    enddo
!
    hnor=dlog10(hnor)/2.0d0
    log_res=dlog10(hes)/2.0d0-hnor
!
    return
!
!
  end subroutine calc_tru_res
!
end module M_la_gscocg_tru_res
