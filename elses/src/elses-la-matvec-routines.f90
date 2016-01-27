!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_matvec_routines
!
  use M_qm_domain
  implicit none
!
  private
  public :: matvec_mul_proj_r
  public :: cg_s_mat_proj_r2
  public :: cg_s_mat_proj_r
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
  subroutine matvec_mul_proj_r(vect_in,vect_out,mat_kind,jsv4jsk,jsk4jsv,jjkset)
!
    use M_qm_domain,     only : noav, nval, dhij, dsij, njsd, jsv4jsd, atm_element  ! (unchanged)
    use M_qm_geno_dst,   only : ham_tot_dst, overlap_dst                ! (unchanged)
    use M_la_matvec_dst, only : matvec_mul_dst_r !(routine)
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
    integer :: jjkset1, ja2, jjk2, ja1, jjk1
    logical, parameter :: use_dst_matrices = .true.
    real(8) :: dbigd
!
    if (use_dst_matrices) then 
!
      if (mat_kind =='H') then
        if ( allocated(ham_tot_dst) ) then 
          call matvec_mul_dst_r(vect_in,vect_out,mat_kind,jsv4jsk,jsk4jsv,jjkset)
          return
        endif  
        if ( .not. allocated(dhij) ) stop 'ERROR(matvec_mul_proj_r):dhij is not allocated'
      endif
!
      if (mat_kind =='S') then
        if ( allocated(overlap_dst) ) then 
          call matvec_mul_dst_r(vect_in,vect_out,mat_kind,jsv4jsk,jsk4jsv,jjkset)
          return
        endif  
        if ( .not. allocated(dsij) ) stop 'ERROR(matvec_mul_proj_r):dsij is not allocated'
      endif
!
    endif  
!
!   write(*,*)'INFO(matvec_mul_proj_r):global matrices are used'
!
    ict4h=1
    m=size(vect_in,1)
    noav2=size(jsk4jsv,1)
    noak=size(jsv4jsk,1)
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
                dbigd=dhij(ja1,ja2,jsd1,jsv2)
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
                dbigd=dsij(ja1,ja2,jsd1,jsv2)
                vect_out(jjk2)=vect_out(jjk2)+dbigd*vect_in(jjk1)
              enddo
            enddo
          endif  
        enddo   
      enddo   
    endif
!
  end subroutine matvec_mul_proj_r
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cg_s_mat_proj_r2(b,x,eps,kend,jsv4jsk,jsk4jsv,jjkset)
!    ---> Solve linear eq.: S x = b (REAL VARIALE)
!          convergence criteria : 
!             log_10 | r / b | < EPS
!      kend : (in input ) the maximum iteration number
!           : (in output) the last iteration number that is executed
!
!   use M_la_matvec_io, only : matvec_mul
    use M_qm_projection, only : matvec_mul_proj ! (routine)
!
    implicit none
    real(8),       intent(in)     :: b(:)
    real(8),       intent(inout)  :: x(:)
    real(8),       intent(inout)  :: eps
    integer,        intent(inout) :: kend
    integer,           intent(in) :: jsv4jsk(:)
    integer,           intent(in) :: jsk4jsv(:)
    integer,           intent(in) :: jjkset(:)
!
    real(8), allocatable          :: r(:), p(:), ap(:)
    real(8)                       :: alpha, beta, rho0, rho1, tbs
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
    beta = 0.0d0
!
    call matvec_mul_proj_r(x,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!        ---> Mat-vec multiplication : ap = S x
!
    r(:) = b(:)-ap(:)  ! r_0 = b - A x_0
    p(:) = r(:)        ! p_0 = r_0
!
    rho0 = dot_product(r(:),r(:))
    hnor = dot_product(b(:),b(:))
    hnor=dlog10(hnor)/2.0d0
!
    kk=-1
    hal = rho0
    hg = dlog10(hal)/2.0d0 - hnor
    if(hg .le. eps) then
      eps= dlog10(hal)/2.0d0 - hnor
      kend=kk
      return
    endif    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main CG loop
!
    do kk = 0, kend_in
!
      hal = rho0
      hg = dlog10(hal)/2.0d0 - hnor
!
      if(hg .le. eps) then
!        --  If converged, the true residual is calculated and exit---
!
         call matvec_mul_proj_r(x,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!          ---> Mat-vec multiplication : ap = S x
!
        r(:) = b(:)-ap(:)
        rho0 = dot_product(r(:),r(:))
        hal = rho0 + 1.0d-100
        eps= dlog10(hal)/2.0d0 - hnor
        kend=kk
        return
      end if
!
      p(:) = r(:) + beta*p(:)
!
      call matvec_mul_proj_r(p,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!        ---> Mat-vec multiplication : ap = S p
      pap = dot_product(p(:),ap(:))
!
      alpha = rho0/pap
      x(:) = x(:) +alpha*p(:)
      r(:) = r(:) -alpha*ap(:)
      rho1=dot_product(r(:),r(:))
      beta = rho1/rho0
      rho0 = rho1
!
    enddo
!
    return
!
!
  end subroutine cg_s_mat_proj_r2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cg_s_mat_proj_r(b,x,eps,kend,jsv4jsk,jsk4jsv,jjkset)
!    ---> Solve linear eq.: S x = b (REAL VARIALE)
!          convergence criteria : 
!             log_10 | r / b | < EPS
!      kend : (in input ) the maximum iteration number
!           : (in output) the last iteration number that is executed
!
!   use M_la_matvec_io, only : matvec_mul
    use M_qm_projection, only : matvec_mul_proj ! (routine)
!
    implicit none
    real(8),       intent(in)  :: b(:)
    real(8),       intent(out) :: x(:)
    real(8),       intent(inout)  :: eps
    integer,        intent(inout) :: kend
    integer,           intent(in) :: jsv4jsk(:)
    integer,           intent(in) :: jsk4jsv(:)
    integer,           intent(in) :: jjkset(:)
!
    real(8), allocatable          :: r(:), p(:), ap(:)
    real(8)                       :: alpha, beta, rho0, rho1, tbs
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
      rho0  = rho0 + r(i)*r(i)
      hnor  = hnor + b(i)*b(i)
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
     call matvec_mul_proj_r(p,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!        ---> Mat-vec multiplication : ap = S p
!
     pap = 0.0d0
     do i = 1, m
       pap  = pap + p(i)*ap(i)
     enddo
!
     hal = rho0
     hg = dlog10(hal)/2.0d0 - hnor
!
     if(hg .le. eps) then
!        --  If converged, 
!           the true residual is calculated and exit---
!
         call matvec_mul_proj_r(x,ap,'S',jsv4jsk,jsk4jsv,jjkset)
!          ---> Mat-vec multiplication : ap = S x
!
        hal = 0.0d0
        do i = 1, m
          tbs = b(i)-ap(i)
          hal = hal + tbs*tbs
        enddo  
        hal=hal+1.0d-100
!
        write(50,*) kk,  dlog10(hal)/2.0d0 - hnor
!
        eps= dlog10(hal)/2.0d0 - hnor
        kend=kk
        return
      end if
!
      alpha = rho0/pap
      rho1 = 0.0d0
      do i = 1, m
        x(i) = x(i) + alpha*p(i)
        r(i) = r(i) - alpha*ap(i)
        rho1 = rho1 + r(i)*r(i)
      enddo  
      beta = rho1/rho0
      rho0 = rho1
!
    enddo
!
    return
!
!
  end subroutine cg_s_mat_proj_r
!
end module M_la_matvec_routines
