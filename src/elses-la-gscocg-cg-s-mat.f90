!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_gscocg_cg_s_mat
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
  public :: cg_s_mat_proj
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cg_s_mat_proj(b,x,eps,kend,jsv4jsk,jsk4jsv,jjkset)
!    ---> Solve linear eq.: S x = b 
!          convergence criteria : 
!             log_10 | r / b | < EPS
!      kend : (in input ) the maximum iteration number
!           : (in output) the last iteration number that is executed
!
!   use M_la_matvec_io, only : matvec_mul
    use M_qm_projection, only : matvec_mul_proj ! (routine)
!
    implicit none
    complex(8),       intent(in)  :: b(:)
    complex(8),       intent(out) :: x(:)
    real(8),          intent(in)  :: eps
    integer,        intent(inout) :: kend
    integer,           intent(in) :: jsv4jsk(:)
    integer,           intent(in) :: jsk4jsv(:)
    integer,           intent(in) :: jjkset(:)
!
    complex(8), allocatable       :: r(:), p(:), ap(:)
    complex(8)                    :: alpha, beta, rho0, rho1, tbs
    real(8)                       :: hg, hal, pap,hnor
    integer                       :: m, m2, i, ierr
    integer                       :: kk, kend_in
!
!
    kend_in=kend
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get matrix size : m
!
    m=size(b,1)
    m2=size(x,1)
    if (m /= m2) stop 'Parameter Mismatch:m,m2'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix allocation
!
    allocate (r(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (p(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (ap(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iniital value
!
    hnor = 0.0d0
    rho0 = dcmplx(0.0d0,0.0d0)
    beta = dcmplx(0.0d0,0.0d0)
!
    do i = 1, m
      x(i)  = dcmplx(0.0d0,0.0d0)
      r (i) = b(i)
      p (i) = dcmplx(0.0d0,0.0d0)
      rho0  = rho0 + conjg(r(i))*r(i)
      hnor  = hnor + conjg(b(i))*b(i)
    enddo
    hnor=dlog10(hnor)/2.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main CG loop
!
    do kk = 0, kend_in
!
     hal = 0.0d0
     do i = 1, m
       p(i) = r(i) + beta*p(i)
     enddo
!
     call matvec_mul_proj(p,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!        ---> Mat-vec multiplication : ap = S p
!
     pap = 0.0d0
     do i = 1, m
       pap  = pap + conjg(p(i))*ap(i)
     enddo
!
     hal = rho0
     hg = dlog10(hal)/2.0d0 - hnor
!
     if(hg .le. eps) then
!        --  If converged, 
!           the true residual is calculated and exit---
!
         call matvec_mul_proj(x,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!          ---> Mat-vec multiplication : ap = S x
!
        hal = dcmplx(0.0d0,0.0d0) 
        do i = 1, m
          tbs = b(i)-ap(i)
          hal = hal + conjg(tbs)*tbs
        enddo  
        write(50,*) kk,  dlog10(hal)/2.0d0 - hnor
!
        kend=kk
        return
      end if
!
      alpha = rho0/pap
      rho1 = 0.0d0
      do i = 1, m
        x(i) = x(i) + alpha*p(i)
        r(i) = r(i) - alpha*ap(i)
        rho1 = rho1 + conjg(r(i))*r(i)
      enddo  
      beta = rho1/rho0
      rho0 = rho1
!
    enddo
!
    return
!
!
  end subroutine cg_s_mat_proj
!
end module M_la_gscocg_cg_s_mat




