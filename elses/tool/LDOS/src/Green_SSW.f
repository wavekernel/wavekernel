!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
************************************************************************
      module sCOCG
*     apply_cocg, apply_cg, Gjk
************************************************************************
      use common_constants
      use hamiltonian, only: apply_hamiltonian
      implicit none

      contains

************************************************************************
      subroutine apply_cocg(nd,n,z,b,v,alpha,beta,br,yr,y,yk,res,rnorm,
     $     xnorm,interval,rvecs)
*     On entering
*     nd: dimension of the vectors
*     n : Maximum number of CG step
*     z : reference energy (complex)
*     b : solve (z-H)x=b. Here after "A" denotes  (z-H)
*     y(nd,m) : Adjoint vectors.
*     yk(m)   : If y is absent, adjoint vectors are
*     {e_i}i=yk(1),...,yk(n). All element of e_i vector equals to zero,
*     except yk(i)-th element equals to one.
*     res:      acceptable accuracy of the approximate solution
*     interval: interval of printing ||r_k|| on the screen.
*     On exit
*     b : preserved
*     n : actual number of iteration
*     v(:,1)  = Real      part of x (solution)
*     v(:,2)  = Imaginary part of x (solution)
*     v(:,3:4)is used for work area and its contents are destroyed.
*     alpha(k)=alpha_{k-1}
*     beta (k)= beta_{k-1}
*     If br is present,       br(k)  = b        r_{k-1}
*     If y or yk is present, yr(i,k) = y(:,i)^T r_{k-1}
*     If rnorm is present, rnorm(k)  = ||r_{k-1}||
*     If xnorm is present, xnorm     = ||x_{m-1}||
*     If rvecs is present, rvecs(:,1)=\bm{r}_{m-2},
*                          rvecs(:,2)=\bm{r}_{m-1}
*                          rvecs(:,3)=\bm{r}_{m}
*     Reference:
*     R.Takayama, T.Hoshi, T.Sogabe, S.-L.Zhang, T.Fujiwara,
*     PRB 73, 165108
*     Recurrance equation:
*     x_0=p_{-1}=0, r_0=b, \alpha_{-1}=1, \beta_{-1}=0
*     x_n=x_{n-1}+\alpha_{n-1}p_{n-1}
*     r_n=r_{n-1}-\alpha_{n-1}p_{n-1}
*     p_n=r_{n}  + \beta_{n-1}p_{n-1}
*     \alpha_{n-1}=(r^T_{n-1}r_{n-1})/(p^T_{n-1}Ap_{n-1})
*     \beta_{n-1} =(r^T_{n}r_{n})/(r^T_{n-1}r{n-1})
*     Extension (May 28, 2007)
*     If y(:,m) is given, a dot product of y(:,i) and r(:) is
*     calculated on each step and stored into yr(:).
*     "yr" is used for later calculation of some matrix elements of
*     resolvent (Green's function).
*     m=nproj
************************************************************************
      use parallel_operation, only: init_workshare, get_NPROC, pdprod
      implicit none
      integer      :: nd, n, yk(:), interval
      real(knd)    :: b(nd), v(:,:), y(:,:), res, rnorm(:), xnorm
c                                    y(nd,nproj)
      complex(knd) :: z, alpha(:), beta(:), br(n), yr(:,:), rvecs(:,:)
c                                                  yr(nproj,n)
      intent(in)   :: nd, z, b, y, yk, res, interval
      intent(out)  :: v, alpha, beta, br, yr, rnorm, xnorm, rvecs
      intent(inout):: n
      optional  :: alpha, beta, br, yr, y, yk, res, rnorm, rvecs
      optional  :: xnorm, interval
*     local variables
      complex(knd) :: norm_r1, norm_r0, Anormp, csum, alp, bet
      complex(knd),allocatable,dimension(:)  :: x, p, r, Ap, csum_p
      complex(knd),allocatable,dimension(:,:):: csum_q
      real(knd) :: tmp,zr,zi
      integer :: np,id,lb,ub, i,j,k, nproj
      logical :: failed
*     FROM BLAS
      complex(knd) :: zdotu
      real(knd)    :: dznrm2
*     FOR DEBUG (check non-trivial relation)
c!!!      complex(knd) :: Ar(nd)

c     initialization & parameter checking

      zr=real(z,knd)
      zi=imag(z)

      if(present(y).and.present(yk))
     $     stop 'apply_cocg: y and yk are exclusive with each other.'
      if(size(v,1).lt.nd)stop'apply_cocg: length of v is insufficient.'
      if(size(v,2).lt.4 )stop'apply_cocg: work area is insufficient.'

      failed = .false.
      
      allocate(x(nd))
      allocate(p(nd))
      allocate(r(nd))
      allocate(Ap(nd))
      call get_NPROC(i)
      if(present(br)) allocate(csum_p(i))
      if(present(yr)) then
         nproj=size(yr,1)
         if(present(y))then
            if(size(y,2).ne.nproj) stop
     $           'apply_cocg: size mismatch (y,yr).'
            if(size(y,1).ne.nd) stop
     $           'apply_cocg: size mismatch first dim. of y.'
            if(size(yr,2).ne.n) stop
     $           'apply_cocg: size mismatch second dim. of yr.'
            allocate(csum_q(i,nproj))
         else if(present(yk)) then
            if(size(yk).ne.nproj) stop
     $           'apply_cocg: size mismatch (yk,yr).'
         else
            stop 'apply_cocg: no indication about y'
         end if
      else if(present(y) .or. present(yk))then
         stop 'apply_cocg: no y but yr exist'
      end if
c     end of initialization & parameter checking

*-begin-Conjugate-Orthogonal-Conjugate-Gradient block
*     initialize
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(nd,np,id,lb,ub)
      x(lb:ub)=0d0
      r(lb:ub)=b(lb:ub)
      p(lb:ub)=b(lb:ub)
      if(present(rvecs))then
         ! for strict compiler which complaints about the
         ! usage of non-initialized arrays or variables.
         rvecs(lb:ub,1)=0d0
         rvecs(lb:ub,2)=0d0
         rvecs(lb:ub,3)=0d0
      end if
C$OMP END PARALLEL
      norm_r0=zdotu(nd,r,1,r,1)
      do i=1, n
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         v(lb:ub,1)=real(p(lb:ub),knd)
         v(lb:ub,2)=imag(p(lb:ub))
C$OMP END PARALLEL
         call apply_hamiltonian(nd,v(:,1),v(:,3)) !Ap=Hp (Re)
         call apply_hamiltonian(nd,v(:,2),v(:,4)) !Ap=Hp (Im)
C$OMP PARALLEL PRIVATE(np,id,lb,ub,csum)
         call init_workshare(nd,np,id,lb,ub)
         Ap(lb:ub)=z*p(lb:ub)-cmplx(v(lb:ub,3),v(lb:ub,4),knd)
                                !Ap=(z-H)p_{i-1}
         if(present(br)) then
            csum=cmplx(0d0,0d0,knd)
            do j=lb, ub
               csum=csum+b(j)*r(j)
            end do
            csum_p(id)=csum
         end if
         if(present(y))then
            do j=1, nproj
               csum=cmplx(0d0,0d0,knd)
               do k=lb, ub
                  csum=csum+y(k,j)*r(k)
               end do
               csum_q(id,j)=csum
            end do
         end if
C$OMP END PARALLEL
         if(present(br))
     $        br(i)=sum(csum_p)  ! br(i)=b.r_{i-1}
         if(present(y))then
            do j=1, nproj
               yr(j,i)=sum(csum_q(:,j)) !yr(j,i)=y(:,j).r_{i-1}
            end do
         else if(present(yk))then
            do j=1, nproj
               yr(j,i)=r(yk(j))
            end do
         end if
         Anormp=zdotu(nd,p,1,Ap,1)
         tmp=abs(Anormp)
c     ERROR handling
         if (tmp.le.tiny(tmp)) then
            failed=.true.
         else if(tmp.lt.1d0) then
            if(huge(tmp)*tmp .lt. abs(norm_r0)) failed=.true.
         end if
         if(failed)then
            write(*,'("apply_cocg: CANNOT CONTINUE(1)")')
            write(*,'("|r0|=",ES14.7,", |pAp|=",ES14.7)')
     $           sqrt(abs(norm_r0)),sqrt(tmp)
            exit
         end if
c
         alp=norm_r0/Anormp     ! (alpha_{i-1})
         if(present(alpha)) alpha(i)=alp ! (alpha_{i-1})
         if(present(rnorm).or.present(interval)) then
            tmp=dznrm2(nd,r,1) ! dznorm2 = sqrt(r.r)
            if(present(rnorm)) rnorm(i)=tmp
            if(present(interval))then
               if(mod(i,interval)==0) then
                  write(*,'("apply_cocg:iter=",I4,
     $                 ", ||r_k||=",ES21.14)') i, tmp
               end if
            end if
         end if
* for output ||r_n||
c         write(40,*) i, dznrm2(nd,r,1), 'COCG'
c REDUCTION does not properly work in ifort 10
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         if(present(rvecs))then
            rvecs(lb:ub,1)=rvecs(lb:ub,2)
            rvecs(lb:ub,2)=rvecs(lb:ub,3)
            rvecs(lb:ub,3)=r    (lb:ub)
         end if
         x(lb:ub)=x(lb:ub)+alp* p(lb:ub)
                                !x_{i}=x_{i-1}+\alpha_{i-1}p_{i-1}
         r(lb:ub)=r(lb:ub)-alp*Ap(lb:ub)
                                !r_{i}=r_{i-1}-\alpha_{i-1}Ap_{i-1}
C$OMP END PARALLEL
         norm_r1=zdotu(nd,r,1,r,1)
         tmp    =abs(norm_r0)
*     TERMINAL CONDITION
         if(present(res))then
            if( sqrt(tmp) .lt. res ) exit
         end if
c     ERROR handling
         if(tmp.le.tiny(tmp)) then
            failed=.true.
         else if(tmp .lt. 1d0) then
            if(huge(tmp)*tmp .lt. abs(norm_r1)) failed=.true.
         end if
         if(failed)then
            write(*,'("apply_cocg: CANNOT CONTINUE(2)")')
            write(*,'("r0=",ES14.7,", r1=",ES14.7)')
     $           sqrt(tmp),sqrt(abs(norm_r1))
            exit
         end if
c
         bet=norm_r1/norm_r0    !  (beta_{i-1})
         if(present(beta)) beta(i)=bet !  (beta_{i-1})
         norm_r0=norm_r1
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         p(lb:ub)=r(lb:ub)+bet*p(lb:ub)
                                !p_{i}=r_{i}+\beta_{i-1}p_{i-1}
C$OMP END PARALLEL
C     DEBUG (check non-trivial relation)
c!!!         v(:,1)=real(r(:),knd)
c!!!         v(:,2)=imag(r(:))
c!!!         call apply_hamiltonian(nd,v(:,1),v(:,3))
c!!!         call apply_hamiltonian(nd,v(:,2),v(:,4))
c!!!         Ar(:)=z*r(:)-cmplx(v(:,3),v(:,4),knd)
c!!!         write(30,'("("ES21.13,",",ES21.13,")",
c!!!     $        " (",ES21.13,",",ES21.13,")")')
c!!!     $        norm_r1/(zdotu(nd,r,1,Ar,1)-(bet/alp*norm_r1)),
c!!!     $        alp
c     \alpha_i=\frac{r_i^T r_i}
c     {r_i^T A r_i - \frac{\beta_{i-1}}{\alpha_{i-1}} r_i^T r_i}
c     alp=\alpha_{i-1}, bet=\beta_{i-1}
C     Another check (r_{i}^T r_{i-1}=0: OK)
c!!!         write(30,'("("ES21.13,",",ES21.13,")")') zdotu(nd,r,1,Ar,1)
c!!!         Ar=r
C     END DEBUG
      end do
*-end-Conjugate-Orthogonal-Conjugate-Gradient block
*
      write(*,'("apply_cocg: sqrt(r1.r1) =",ES21.14,
     $     ", iteration=",I6)') sqrt(abs(norm_r1)), i
*
*     Check actual residue
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(nd,np,id,lb,ub)
      v(lb:ub,1)=real(x(lb:ub),knd)
      v(lb:ub,2)=imag(x(lb:ub))
C$OMP END PARALLEL
      call apply_hamiltonian(nd, v(:,1), v(:,3)) ! Ap=Ax
      call apply_hamiltonian(nd, v(:,2), v(:,4)) ! Ap=Ax
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(nd,np,id,lb,ub)
      v(lb:ub,3)=zr*v(lb:ub,1)-zi*v(lb:ub,2)-v(lb:ub,3)-b(lb:ub)
      v(lb:ub,4)=zr*v(lb:ub,2)+zi*v(lb:ub,1)-v(lb:ub,4)
                                ! Ap=(z-H)x-b
C$OMP END PARALLEL
      write(*,'("apply_cocg: ||(e-A)x-b||=",ES21.14)')
     $        sqrt(pdprod(v(:,3),v(:,3))+pdprod(v(:,4),v(:,4)))

      if(present(xnorm))xnorm=pdprod(v(:,1),v(:,1))
     $                       +pdprod(v(:,2),v(:,2))
      n=i-1                     ! n=i-1
      deallocate(x)
      deallocate(p)
      deallocate(r)
      deallocate(Ap)
      if(allocated(csum_p)) deallocate(csum_p)
      if(allocated(csum_q)) deallocate(csum_q)
      end subroutine apply_cocg

************************************************************************
      subroutine apply_cg(nd,n,e,b,v,alpha,beta,br,yr,y,yk,res,rnorm,
     $     xnorm,interval,rvecs)
*     On entering
*     nd: dimension of the vectors
*     n : Maximum number of CG step
*     e : reference energy
*     b : solve (e-H)x=b. Here after "A" denotes  (e-H)
*     y(nd,m) : Adjoint vectors.
*     yk(m)   : If y is absent, adjoint vectors are
*     {e_i}i=yk(1),...,yk(n). All element of e_i vector equals to zero,
*     except yk(i)-th element equals to one.
*     res:      acceptable accuracy of the approximate solution
*     interval: interval of printing ||r_k|| on the screen.
*     On exit
*     b : preserved
*     n : actual number of iteration
*     v(:,1)  = x (solution)
*     v(:,2:4)is used for work area and its contents are destroyed.
*     alpha(k)=alpha_{k-1}
*     beta (k)= beta_{k-1}
*     If br is present,       br(k)  = b        r_{k-1}
*     If y or yk is present, yr(i,k) = y(:,i)^T r_{k-1}
*     If rnorm is present, rnorm(k)  = ||r_{k-1}||
*     If xnorm is present, xnorm     = ||x_{m-1}||
*     If rvecs is present, rvecs(:,1)=\bm{r}_{m-2},
*                          rvecs(:,2)=\bm{r}_{m-1}
*                          rvecs(:,3)=\bm{r}_{m}
*     Reference:
*     R.Takayama, T.Hoshi, T.Sogabe, S.-L.Zhang, T.Fujiwara,
*     PRB 73, 165108
*     Recurrance equation:
*     x_0=p_{-1}=0, r_0=b, \alpha_{-1}=1, \beta_{-1}=0
*     x_n=x_{n-1}+\alpha_{n-1}p_{n-1}
*     r_n=r_{n-1}-\alpha_{n-1}p_{n-1}
*     p_n=r_{n}  + \beta_{n-1}p_{n-1}
*     \alpha_{n-1}=(r^T_{n-1}r_{n-1})/(p^T_{n-1}Ap_{n-1})
*     \beta_{n-1} =(r^T_{n}r_{n})/(r^T_{n-1}r{n-1})
*     Extension (May 28, 2007)
*     If y(:,m) is given, a dot product of y(:,i) and r(:) is
*     calculated on each step and stored into yr(:).
*     "yr" is used for later calculation of some matrix elements of
*     resolvent (Green's function).
*     m=nproj
************************************************************************
      use parallel_operation, only: init_workshare, get_NPROC, pdprod
      implicit none
      integer      :: nd, n, yk(:), interval
      real(knd)    :: b(nd), v(:,:), y(:,:), res, rnorm(:), xnorm
c                                    y(nd,nproj)
      real(knd)    :: e, alpha(:), beta(:), br(n), yr(:,:), rvecs(:,:)
c                                                  yr(nproj,n)
      intent(in)   :: nd, e, b, y, yk, res, interval
      intent(out)  :: v, alpha, beta, br, yr, rnorm, xnorm, rvecs
      intent(inout):: n
      optional  :: alpha, beta, br, yr, y, yk, res, rnorm, rvecs
      optional  :: xnorm, interval
*     local variables
      real(knd) :: norm_r1, norm_r0, Anormp, tmp, alp, bet
      real(knd),allocatable,dimension(:)   :: sum_p
      real(knd),allocatable,dimension(:,:) :: sum_q
      integer :: np,id,lb,ub, i,j, nproj
      logical :: failed

* for output ||r_n||
c     write(40,*)
c     initialization & parameter checking
      if(present(y).and.present(yk))
     $     stop 'apply_cg: y and yk are exclusive with each other.'
      if(size(v,1).lt.nd)stop'apply_cg: length of v is insufficient.'
      if(size(v,2).lt.4 )stop'apply_cg: work area is insufficient.'

      failed = .false.

      call get_NPROC(i)
      if(present(br)) allocate(sum_p(i))
      if(present(yr)) then
         nproj=size(yr,1)
         if(present(y))then
            if(size(y,2).ne.nproj) stop
     $           'apply_cg: size mismatch (y,yr).'
            if(size(y,1).ne.nd) stop
     $           'apply_cg: size mismatch first dim. of y.'
            if(size(yr,2).ne.n) stop
     $           'apply_cg: size mismatch second dim. of yr.'
            allocate(sum_q(i,nproj))
         else if(present(yk)) then
            if(size(yk).ne.nproj) stop
     $           'apply_cg: size mismatch (yk,yr).'
         else
            stop 'apply_cg: no indication about y'
         end if
      else if(present(y) .or. present(yk))then
         stop 'apply_cg: no y but yr exist'
      end if
c     end of initialization & parameter checking

*-begin-Conjugate-Gradient block
*     initialize
c     x->v(:,1),p->v(:,2),Ap->,v(:,3),r->v(:,4)
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(nd,np,id,lb,ub)
      v(lb:ub,1)=0d0
      v(lb:ub,2)=b(lb:ub) !p=b
      v(lb:ub,4)=b(lb:ub) !r=b
      if(present(rvecs))then
         rvecs(lb:ub,1)=0d0
         rvecs(lb:ub,2)=0d0
         rvecs(lb:ub,3)=0d0
      end if
C$OMP END PARALLEL
      norm_r0=pdprod(v(:,4),v(:,4)) ! r_{0} r_{0}
      do i=1, n
         call apply_hamiltonian(nd,v(:,2),v(:,3)) !Ap=Hp_{i-1}
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         v(lb:ub,3)=e*v(lb:ub,2)-v(lb:ub,3) !Ap=(e-H)p_{i-1}
         if(present(br))
     $         sum_p(id)  =dot_product(b(lb:ub),  v(lb:ub,4))
         if(present(y))then
            do j=1, nproj
               sum_q(id,j)=dot_product(y(lb:ub,j),v(lb:ub,4))
            end do
         end if
C$OMP END PARALLEL
         if(present(br))
     $        br(i)=sum(sum_p)  ! br(i)=b.r_{i-1}
         if(present(y))then
            do j=1, nproj
               yr(j,i)=sum(sum_q(:,j)) !yr(j,i)=y(:,j).r_{i-1}
            end do
         else if(present(yk))then
            do j=1, nproj
               yr(j,i)=v(yk(j),4) ! picking up yk(j)-th elements of r
            end do
         end if
         Anormp=pdprod(v(:,2),v(:,3)) ! Anormp=p.Ap
         tmp=abs(Anormp)
c     ERROR handling
         if (tmp.le.tiny(tmp)) then
            failed=.true.
         else if(tmp.lt.1d0) then
            if(huge(tmp)*tmp .lt. norm_r0) failed=.true.
         end if
         if(failed)then
            write(*,'("apply_cg: CANNOT CONTINUE(1)")')
            write(*,'("|r0|=",ES14.7,", |pAp|=",ES14.7)')
     $           sqrt(norm_r0),sqrt(tmp)
            exit
         end if
c
         alp=norm_r0/Anormp     ! (alpha_{i-1})
         if(present(alpha)) alpha(i)=alp ! (alpha_{i-1})
         if(present(rnorm).or.present(interval)) then
            tmp=sqrt(norm_r0)
            if(present(rnorm)) rnorm(i)=tmp
            if(present(interval))then
               if(mod(i,interval)==0) then
                  write(*,'("apply_cg:iter=",I4,
     $                 ", ||r_k||=",ES21.14)') i, tmp
               end if
            end if
         end if
* for output ||r_n||
c        write(40,*) i, sqrt(norm_r0), 'CG'
c REDUCTION does not properly work in ifort 10
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         if(present(rvecs))then
            rvecs(lb:ub,1)=rvecs(lb:ub,2) !r_{i-3}
            rvecs(lb:ub,2)=rvecs(lb:ub,3) !r_{i-2}
            rvecs(lb:ub,3)=v    (lb:ub,4) !r_{i-1}
         end if
         v(lb:ub,1)=v(lb:ub,1)+alp*v(lb:ub,2)
                                !x_{i}=x_{i-1}+alpha*p_{i-1}
         v(lb:ub,4)=v(lb:ub,4)-alp*v(lb:ub,3)
                                !r_{i}=r_{i-1}-alpha*Ap_{i-1}
C$OMP END PARALLEL
         norm_r1=pdprod(v(:,4),v(:,4))    !r_i.r_i
         tmp    =norm_r0
*     TERMINAL CONDITION
         if(present(res))then
            if( sqrt(tmp) .lt. res ) exit
         end if
c     ERROR handling
         if(tmp.le.tiny(tmp)) then
            failed=.true.
         else if(tmp .lt. 1d0) then
            if(huge(tmp)*tmp .lt. norm_r1) failed=.true.
         end if
         if(failed)then
            write(*,'("apply_cg: CANNOT CONTINUE(2)")')
            write(*,'("r0=",ES14.7,", r1=",ES14.7)')
     $           sqrt(tmp),sqrt(norm_r1)
            exit
         end if
c
         bet=norm_r1/norm_r0    !  (beta_{i-1})
         if(present(beta)) beta(i)=bet !  (beta_{i-1})
         norm_r0=norm_r1
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         v(lb:ub,2)=v(lb:ub,4)+bet*v(lb:ub,2) ! p=r+beta*p
C$OMP END PARALLEL
      end do
*-end-Conjugate-Gradient block
*
      write(*,'("apply_cg: sqrt(r1.r1) =",ES21.14,
     $     ", iteration=",I6)') sqrt(norm_r1), i
*
*     Check actual residue
      call apply_hamiltonian(nd, v(:,1), v(:,3)) ! Ap=Ax
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(nd,np,id,lb,ub)
      v(lb:ub,3)=e*v(lb:ub,1)-v(lb:ub,3)-b(lb:ub) ! Ap=(e-H)x-b
C$OMP END PARALLEL
      write(*,'("apply_cg: ||(e-A)x-b||=",ES21.14)')
     $        sqrt(pdprod(v(:,3),v(:,3)))

      if(present(xnorm)) xnorm=pdprod(v(:,1),v(:,1))
      n=i-1                     ! n=i-1
      end subroutine apply_cg

************************************************************************
      function Gjk(sigma,m,j,k,a,b,yr,rnorm,er,interval,res)
*     calculate g_{jk}=y_j^T (sigma+Eref-H)^{-1} b_k
*     on entering
*     sigma: shift (z-Eref)
*     j,k:   indeces of Green function
*     a:     alpha_{i-1}=a(i)
*     b:     beta_{i-1} =b(i)
*     yr:    yr(j,i,k)=y_j \cdot r^{(i-1)}_{k}
*           (initial vector = b_k, residue in (i-1)-th iteration.)
*     er:     the error of approximate solution.
*     interval: interval of printing ||r_k|| or |y^T r_k| on the screen.
*     res:   terminal condition
************************************************************************
      implicit none
      complex(knd) :: gjk, sigma




c     complex version
      complex(knd) :: a(:), b(:), yr(:,:,:), alpha_m2, beta_m2

c
      real(knd)    :: er, res, rnorm(:)
      integer      :: m,j,k, interval
      intent(in)   :: sigma,m,j,k,a,b,yr,rnorm,interval,res
      intent(out)  :: er
      optional     :: rnorm, er, interval, res
*     local variables
c     alpha_m2 and beta_m2 are also local variables, but
c     located above for easy transformation to real version
      real(knd)    :: err
      complex(knd) :: xsigma, psigma
                                !x^\sigma_{i-1}=xsigma
                                !p^\sigma_{i-1}=psigma
      complex(knd) :: alphasigma, betasigma
                                !\alpha^\sigma_{i-1}=alphasigma
                                ! \beta^\sigma_{i-1}= betasigma
      complex(knd) :: pi, pi_m1, pi_m2, c, d
                                !pi_{i}  =pi
                                !pi_{i-1}=pi_m1
                                !pi_{i-2}=pi_m2
      integer      :: i
      character(len=17) :: err_string
      if(m.lt.2)         stop 'Gjk: too small Kryrov space'
      if(m.gt.size(a,1)) stop 'Gjk: too large m, or too small a'
      if(m.gt.size(b,1)) stop 'Gjk: too large m, or too small b'
      if(present(interval).and.(.not.present(rnorm)))then
         write(*,'("Gjk: Projected Error is printed out.")')
      end if

      if(present(rnorm))then
         err_string=", ||r^\sigma_k||="
      else
         err_string=", |y.r^\sigma_k|="
      end if

      pi_m1=cmplx(1d0,0,knd)    !pi_{-1}
      pi   =pi_m1               !pi_{0}
      beta_m2 =0d0              !beta_{-1}=0d0
      alpha_m2=1d0              !alpha_{-1}=1d0
      xsigma=cmplx(0d0,0d0,knd) !x_0
      betasigma=0d0
      psigma   =0d0             !arbitrary
      do i=1, m
c     Recurrence equations
         pi_m2=pi_m1
         pi_m1=pi
         psigma=yr(j,i,k)/pi_m1+betasigma*psigma
c        y^T p^\sigma_{i-1}= y^T ( r^\sigma_{i-1} /pi_{i-1}
c                             +\beta^\sigma_{i-2}p^\sigma_{i-2} )
         c=beta_m2*a(i)/alpha_m2
c        c=\frac{\alpha_{i-1}}{\alpha_{i-2}} \beta_{i-2}
         pi=(1d0+a(i)*sigma+c)*pi_m1-c*pi_m2 ! pi=pi_{i+1}
         d=pi_m1/pi             !d =\frac{\pi_{i-1}}{\pi_{i}}
         alphasigma=  d*a(i)    !alpha^\sigma_{i-1} a(i)=alpha_{i-1}
         betasigma =d*d*b(i)    ! beta^\sigma_{i-1} b(i)= beta_{i-1}
         xsigma=xsigma        +alphasigma*psigma
c        x^\sigma_{i}=x^\sigma_{i-1}+\alpha^\sigma_{i-1}p^\sigma_{i-1}
c         write(*,*) alphasigma
         alpha_m2=a(i)
         beta_m2 =b(i)
c   end of  Recurrence equations

c     Output estimated errors of Green's function to STDOUT
         if(present(rnorm))then
            err=abs(rnorm(i)/pi_m1)
            if(present(res))then
               if(err.lt.res) exit
            end if
         else
            err=abs(yr(j,i,k)/pi_m1)
            if(present(res))then
               if(err.lt.res) exit
            end if
         end if
         if(present(interval))then
            if(mod(i,interval)==0 .or. i==m)then
               write(*,'("Gjk:\sigma=(",ES21.14,",",ES21.14,
     $              "), iter=",I4,A17,ES21.14)')
     $              sigma,i,err_string,err
            end if
         end if
      end do

      Gjk=xsigma
c      if( err .gt. 1d-6 )
c     $     write(*,'("Gjk: large error at sigma=(",
c     $     ES21.14,","ES21.14,"). Err=",ES21.14)')
c     $     sigma, err
      if(present(er)) then
c     if "er" is present, return the error of approximate solution
         if(present(rnorm)) then
c     if "rnorm" is present, the error is estimated as ||r^\sigma_k||
            er = abs(rnorm(m)/pi_m1)
         else
c     if "rnorm" is not present, the error is estimated as
c     |y^T r^\sigma_k|
            er = abs(yr(j,m,k)/pi_m1)
         end if
      end if

      end function Gjk


************************************************************************
      subroutine switch_seed(nd,m,j,k,z,sigma,alpha,beta,yr,
     $     res,rnorm,rvecs,yk,y,interval)
*     Seed SWitching (Set new reference energy and alpha, beta, y^T r )
*     Ref. T. Sogabe, T. Hoshi, S.-L. Zhang, and T. Fujiwara,
*     arXiv:math/0602652v1 (2006).
*
*     From the eq.(13) of PRB 73, 165108 (2006)
*     \bm{r}_{n}
*     =(1+\alpha_{n-1}(\frac{\beta_{n-2}}{\alpha{n-2}}-A))\bm{r}_{n-1}
*     -\alpha_{n-1}\frac{\beta_{n-2}}{\alpha{n-2}}\bm{r}_{n-2}.
*     Because \bm{r}_{i}'s are orthogonal (for the "inner-product"
*     (\bm{a},\bm{b})=\bm{a}^T \bm{b}) to each other,
*     \bm{r}^T (1+\alpha_{n-1}(\frac{\beta_{n-2}}{\alpha{n-2}}-A))
*     \bm{r}_{n-1}=0. Consequently,
*     \alpha_{n-1}=\frac{\bm{r}_{n-1}^T \bm{r}_{n-1}}
*                       {\bm{r}_{n-1}^T A \bm{r}_{n-1}
*                             -\frac{\beta_{n-2}}{\alpha_{n-2}}
*                              \bm{r}_{n-1}^T \bm{r}_{n-1}},
*     \beta_{n-1}=\frac{\bm{r}_{n}^T   \bm{r}_{n}}
*                      {\bm{r}_{n-1}^T \bm{r}_{n-1}},
*     \bm{r}_{0}=\bm{b}, \alpha_{-1}=1, \beta{-1}=0,
*     \alpha_{0}=\frac{\bm{b}^T \bm{b}}{\bm{b}^T A \bm{b}},
*     \bm{r}_{1}=\bm{b}-\alpha_{0} A \bm{b}.
*        These three recurrence equations and five initial conditions
*     are the complete set of equations needed to determine \bm{r}'s.
*
*     nd    : dimension of Hamiltonian
*     m     : length of alpha, beta, r
*     j     :
*     k     :
*     yk    : input vector = (0,...,0,1[yk(k)-th column],0,...,0)^T
*           adjoint vector = (0,...,0,1[yk(j)-th column],0,...,0)
*     y(:,:): alternative of yk. y(:,k) gives k-th input vec.
*     sigma : energy shift
*     alpha : \alpha_{i-1}
*     beta  :  \beta_{i-1}
*     yr(j,i,k) : \bm{y}_j^T \bm{r}_{i-1} (\bm{r} : k-th input vec.)
*     rnorm : ||\bm{r}_{i-1}||
*     rvecs : rvecs(:,1)=\bm{r}_{m-2}
*           : rvecs(:,2)=\bm{r}_{m-1}
*           : rvecs(:,3)=\bm{r}_{m}
*     On exit,
*     z     : new seed (reference energy) = z+sigma
*     m     : new length of alpha, beta, r
*     alpha : new alpha
*     beta  : new beta
*     rvecs : new r (vectors in last three iterations)
*     yr    : new \bm{y}^T \bm{r}
************************************************************************
      use parallel_operation, only: init_workshare
      implicit none





c     complex version
      complex(knd) :: z, sigma, alpha(:), beta(:), yr(:,:,:), rvecs(:,:)

c
      real(knd)    :: rnorm(:), res, y(:,:) ! (nd,nproj)
      integer      :: nd, m, j, k,  yk(:), interval
      intent(in)   :: nd, j, k, sigma, res, yk, y, interval
      intent(inout):: z, m, alpha, beta, yr, rnorm, rvecs
      optional     :: rnorm, yk, y, interval
*     Margin PARAMETER
      real(knd),parameter :: RMARGIN=1d-1 ! RMARGIN=0.7 for paper
*     local variables
c     x^\sigma_{i-1}=xsigma
c     p^\sigma_{i-1}=psigma
c     \alpha^\sigma_{i-1}=alphasigma
c     \beta^\sigma_{i-1} =betasigma
c     pi_{i}  =pi
c     pi_{i-1}=pi_m1
c     pi_{i-2}=pi_m2
c     complex version
      complex(knd) :: norm_r1, norm_r0, Anormp
      complex(knd) :: xsigma, psigma
      complex(knd) :: alphasigma, betasigma
      complex(knd) :: pi, pi_m1, pi_m2, c, d
      complex(knd) :: alpha_m2, beta_m2
*     local arrays
      complex(knd) :: Ar
      real(knd),allocatable    :: tvec(:,:)

      real(knd)    :: err, tmp
      integer      :: i, ij, nproj, np, id, lb, ub
      logical      :: failed
*     local arrays
      allocatable :: Ar(:)






*     functions
      complex(knd) :: zdotu
      real(knd)    :: dznrm2

      write(*,'("switch_seed: COMPLEX_VERSION")')
      allocate(tvec(nd,4))

      allocate(Ar(nd))

      nproj=size(yr,1)
      if(.not.(present(yk) .or. present(y)))
     $     stop 'switch_seed: "y" or "yk" is needed'
      if(present(yk) .and. present(y))then
         write(*,'("yk is present :",L1)') present(yk)
         write(*,'("y  is present :",L1)') present(y)
         stop 'switch_seed: "yk" and "y" are exclusive.'
      end if
      if(present(yk))then
         if(size(yk,1) .ne. nproj)then
            write(*,*) size(yk,1), nproj
            stop 'switch_seed: size mismatch (yk,yr).'
         end if
      else
         if(size(y,2) .ne. nproj)then
            write(*,*) size(y,2), nproj
            stop 'switch_seed: size mismatch  (y,yr).'
         end if
      end if         

c      write(*,*) 'm',m,size(rnorm,1)
      
c     we assume j=k=1
      if(j .ne. 1 .or. k .ne. 1 )then
         write(*,'("Switch_Seed: j or k may be wrong. (j,k)=",2I16)')j,k
      end if
      if(m.le.2) stop
     $     'switch_seed: Too small CGSTEP. Replace appropreate seed (ref
     $erence energy) may solve this problem.'

      pi_m1=cmplx(1d0,0,knd)    !pi_{-1}
      pi   =pi_m1               !pi_{0}
      beta_m2 =0d0              !beta_{-1}=0d0
      alpha_m2=1d0              !alpha_{-1}=1d0
      xsigma=cmplx(0d0,0d0,knd) !x_0
      betasigma=0d0
      psigma   =0d0             !arbitrary
* for output ||r_n||
c        write(40,*)
      do i=1, m
c     Recurrence equations
         pi_m2=pi_m1
         pi_m1=pi
         yr(j,i,k)=yr(j,i,k)/pi_m1
c     debug
c         write(*,*)'i=',i,',rnorm=',rnorm(i),', pi=',abs(pi_m1)
         if(present(rnorm))rnorm(i)=rnorm(i)/abs(pi_m1)
         psigma=yr(j,i,k)+betasigma*psigma
c        p^\sigma_{i-1}=    r^\sigma_{i-1} /pi_{i-1}
c                      +\beta^\sigma_{i-2}p^\sigma_{i-2}

         c=beta_m2*alpha(i)/alpha_m2
c        c=\frac{\alpha_{i-1}}{\alpha_{i-2}} \beta_{i-2}
         pi=(1d0+alpha(i)*sigma+c)*pi_m1-c*pi_m2 ! pi=pi_{i+1}
cdebug Dec 4, 2007
c         write(*,*) pi_m1,pi
* for output ||r_n||
c        write(40,*) i, rnorm(i), 'Seed_switch:1~M'
         d=pi_m1/pi             !d =\frac{\pi_{i-1}}{\pi_{i}}
         alphasigma=  d*alpha(i)!alpha^\sigma_{i-1} a(i)=alpha_{i-1}
         betasigma =d*d*beta(i) ! beta^\sigma_{i-1} b(i)= beta_{i-1}
         xsigma=xsigma       +alphasigma*psigma
c        x^\sigma_{i}=x^\sigma_{i-1}+\alpha^\sigma_{i-1}p^\sigma_{i-1}
         alpha_m2=alpha(i)
         beta_m2 =beta (i)
         alpha(i)=alphasigma
         beta (i)=betasigma
c   end of  Recurrence equations
      end do
* for output ||r_n||
c     write(40,*)
      failed=.false.
      z=z+sigma
c     debug
c      write(*,*) 'r_m_2=',dznrm2(nd,rvecs(:,1),1)
c      write(*,*) 'r_m_1=',dznrm2(nd,rvecs(:,2),1)
c      write(*,*) 'r_m='  ,dznrm2(nd,rvecs(:,3),1)
c      write(*,*) 'xxxx',dznrm2(nd,rvecs(:,2),1),abs(pi_m1)
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(nd,np,id,lb,ub)
c     rvecs(:,1)=\bm{r}_{m-2}
c     rvecs(:,2)=\bm{r}_{m-1}
      rvecs(lb:ub,1)=rvecs(lb:ub,1)/pi_m2
      rvecs(lb:ub,2)=rvecs(lb:ub,2)/pi_m1
C$OMP END PARALLEL
c     debug
c      write(*,*) 'yyyy',dznrm2(nd,rvecs(:,2),1),abs(pi_m1)



      norm_r0=zdotu(nd,rvecs(:,2),1,rvecs(:,2),1) !r_{m-1}^T r_{m-1}

*-Conjugate-Orthogonal-Conjugate-Gradient block
      do i=m, size(yr,2)



C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         tvec(lb:ub,1)=real(rvecs(lb:ub,2),knd)
                                !     real part of \bm{r}_{i-1}
         tvec(lb:ub,2)=imag(rvecs(lb:ub,2))    
                                !imaginary part of \bm{r}_{i-1}
C$OMP END PARALLEL
         call apply_hamiltonian(nd,tvec(:,1),tvec(:,3))!tvec= H r_{i-1}
         call apply_hamiltonian(nd,tvec(:,2),tvec(:,4))

C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)



         Ar(lb:ub)=z*rvecs(lb:ub,2)
     $            -cmplx(tvec(lb:ub,3),tvec(lb:ub,4),knd)

C$OMP END PARALLEL
         c=beta(i-1)/alpha(i-1) ! \beta_{i-2}/\alpha_{i-2}
c         write(*,*) alpha(i),beta(i)




         alpha(i)=norm_r0/(zdotu(nd,rvecs(:,2),1,Ar,1)-c*norm_r0)
         tmp=dznrm2(nd,rvecs(:,2),1)

* for output ||r_n||
c        write(40,*) i, tmp, 'Seed switch:M~'
c     debug
c         if(i.eq.m) write(*,*) 'rnorm',rnorm(i)
c
c     renew rnorm if it is available
         if(present(rnorm))rnorm(i)=tmp
c
c        write(*,*) 'i=',i,'||r||=',tmp
c     debug
c         if(i.eq.m) write(*,*) 'rnorm',rnorm(i)
         if(present(yk))then
            do ij=1, nproj
               yr(ij,i,k)=rvecs(yk(ij),2) ! yr(j,i,k)=y^T r_{i-1}
            end do
         else
            do ij=1, nproj



c               yr(ij,i,k)=zdotu(nd,y(:,ij),1,rvecs(:,2),1) ! bad (y:real)
               call dgemm('N','N',2,1,nd,1d0,rvecs(:,2),2,
     $              y(:,ij),nd,0d0,yr(ij,i,k),2)

            end do
         end if
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         rvecs(lb:ub,3)=rvecs(lb:ub,2)-alpha(i)*Ar(lb:ub)
     $        +c*alpha(i)*(rvecs(lb:ub,2)-rvecs(lb:ub,1))
              !r_{i}=(1+(c-z+H)*\alpha_{i-1})\bm{r}_{i-1}
              !              -c*\alpha_{i-1} \bm{r}_{i-2}
C$OMP END PARALLEL



         norm_r1=zdotu(nd,rvecs(:,3),1,rvecs(:,3),1)

         beta(i)=norm_r1/norm_r0
c         write(*,*) alpha(i),beta(i)
         norm_r0=norm_r1
* TERMINAL CONDITION
         if(tmp.lt.res*RMARGIN) exit
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
         call init_workshare(nd,np,id,lb,ub)
         rvecs(lb:ub,1)=rvecs(lb:ub,2)
         rvecs(lb:ub,2)=rvecs(lb:ub,3)
C$OMP END PARALLEL
      end do
*-end-Conjugate-Orthogonal-Conjugate-Gradient block
*
      write(*,'("switch_seed: sqrt(r.r) =",ES21.14,
     $     ", iteration=",I6)') sqrt(abs(norm_r0)), i
*
      m=size(yr,2)
      m=min(i,m)


      deallocate(tvec)

      deallocate(Ar)
      end subroutine switch_seed


      end module sCOCG

************************************************************************
      program Green
*     Green's function calculator
************************************************************************
      use common_constants
      use control
      use hamiltonian, only: read_hamiltonian, Ham
      use sparse_matrix, only: apply_matrix
      use parallel_operation, only: init_workshare, pdprod
      use sCOCG
      implicit none
      integer :: CGSTEP, NE1, NE2
      type(control_parameters) :: input

      real(knd), allocatable :: b(:,:), v(:,:), rnorm(:)
      real(knd) :: eps, err, dtmp




      complex(knd),allocatable :: rvecs(:,:)







c     complex version
      complex(knd) :: eref, alpha, beta, yr

      allocatable  :: alpha(:), beta(:), yr(:,:,:)
      integer, allocatable :: yk(:)
      integer   :: ND, j

*     reading parameters from file
      call read_parameters("input.txt",input)
      NE1=1
      NE2=1                     ! G(1:NE1, 1:NE2)
      CGSTEP=input%nitemax

*     reading Hamiltonian from file
c      call read_hamiltonian("Hamiltonian.dat")
                                ! binary data (Private format)
      call read_hamiltonian("Hamiltonian.txt") 
                   ! ascii data (See Dr. Takayama's documents)
*
      call set_list_of_orbitals_to_be_shown(input)
*     allocating prerequisites
      ND=Ham%dim
      allocate(b(ND,NE2))
      allocate(v(ND,4))
      allocate(alpha(CGSTEP))
      allocate(beta (CGSTEP))
      allocate(yr(NE1,CGSTEP,NE2))
      allocate(rnorm(CGSTEP))

      allocate(rvecs(ND,3))

*     end allocating prerequisites

      eps =input%eps_ref
      eps =sqrt(eps)            ! to keep backward compatibility
                                ! (With Takayama version)

c     complex version
      !eref=cmplx(input%Esample,input%gamma,knd)      ! for SSW
      if(input%target_i /= 0) then
         do j=1, input%number_of_orbitals_to_be_shown
            eref=cmplx(input%Esample,input%gamma,knd) ! for SSW
            call set_yk_and_b(input%list_of_orbitals_to_be_shown(j))

            call apply_cocg(ND,CGSTEP,eref,b,v,alpha,beta,
     $           yr=yr(:,:,1),yk=yk,res=eps,rnorm=rnorm,interval=1,
     $           rvecs=rvecs)




            write(*,*) "Eref=",eref
            call write_green(j/=1) ! Append_mode=(j/=1)
         end do
      else
         call set_random_b

         call apply_cocg(ND,CGSTEP,eref,b,v,alpha,beta,
     $        br=yr(:,:,1),res=eps,rnorm=rnorm,rvecs=rvecs)




         write(*,*) "Eref=",eref
         call write_green(.false.)
      end if
c     end complex version


      contains

      subroutine set_yk_and_b(j)
      implicit none
      integer, intent(in) :: j
      integer ::  np,id,lb,ub
c     LDOS calculation Compatible to Takayama
      if(.not. allocated(yk)) allocate(yk(NE1))
      yk(1)=j
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(ND,np,id,lb,ub)
      b(lb:ub,NE2)=0d0
C$OMP END PARALLEL
      b(j,NE2)=1d0
      end subroutine set_yk_and_b

      subroutine set_random_b
      implicit none
      integer ::  np,id,lb,ub
c     Total DOS calc.
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(ND,np,id,lb,ub)
      call random_number(b(lb:ub,NE2)) ! it does not work on PGF90
      b(lb:ub,NE2)=b(lb:ub,NE2)-5d-1
      !because random_number is not thread safe on PGF90
C$OMP END PARALLEL
      dtmp=pdprod(b(:,NE2),b(:,NE2))
C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(ND,np,id,lb,ub)
      b(lb:ub,NE2)=b(lb:ub,NE2)/dsqrt(dtmp)
C$OMP END PARALLEL
      end subroutine set_random_b

      subroutine write_green(append_mode)
*     write g_{jk}(\omega) to the file (for all \omega)
      use common_constants
      use file_io
      logical, intent(in) :: append_mode
      integer :: IUNIT, IUNT2
      character(len=*), parameter :: FNAME="Green.dat"
      character(len=*), parameter :: FNAM2="LDOS.dat"
      character(len=8) :: MODE_STRING
      integer :: i, nmesh, j, nfinished
      real(knd)    :: gamma, denergy, emin
      real(knd)    :: pi, integrated_dos, dos
      complex(knd) :: shift, omega, resolvent

      complex(knd),allocatable :: Gval(:)
      real   (knd),allocatable :: Gerr(:)
      logical :: lfinished
      complex(knd) :: sigma_max
      real(knd)    :: max_rnorm
c     debug
c      integer     ,allocatable :: NSSW(:)
c      integer                  :: NofSSW
c      NofSSW=0

      pi=atan(1d0)*4

      nmesh   = input%nenemax
      denergy = input%denergy
      emin    = input%enmin
      gamma   = input%gamma

      if(append_mode)then
         MODE_STRING='OLD'
      else
         MODE_STRING='REPLACE'
      end if

      IUNIT=vacant_unit()
      open(IUNIT,FILE=FNAME,STATUS=trim(MODE_STRING),POSITION='APPEND')
      IUNT2=vacant_unit()
      open(IUNT2,FILE=FNAM2,STATUS=trim(MODE_STRING),POSITION='APPEND')

      j=1
      integrated_dos=0d0        ! initialize integrated_dos


      allocate(Gval(nmesh))
      allocate(Gerr(nmesh))
c     debug
c     allocate(NSSW(nmesh))
c
      Gerr=huge(1d0)
      lfinished=.false.
      do while (.not. lfinished)
         lfinished=.true.
         max_rnorm=0d0
         nfinished=0
C$OMP PARALLEL DO PRIVATE(omega,shift,err,resolvent),
C$OMP$ REDUCTION(+:nfinished)
         do i=1, nmesh
            if(Gerr(i).lt.eps) cycle
            nfinished=nfinished+1
            omega=cmplx(emin+(i-1)*denergy,gamma,knd)
            shift=omega-eref
            resolvent=Gjk(shift,CGSTEP,j,j,alpha,beta,yr,
     $           rnorm=rnorm,er=err,interval=1,res=eps)
            Gval(i)=resolvent
            Gerr(i)=err
c     debug
c           NSSW(i)=NofSSW
c
c     renew max_rnorm if the latest rnorm is greater than max_rnorm
            if(err.gt.max_rnorm)then
               sigma_max=shift
               max_rnorm=err
            end if
            if(CGSTEP.lt.size(yr,2) .and. err.gt.eps) then
      !if CGSTEP reaches maximum, discontinue regardless of precision.
               lfinished=.false.
            end if
         end do
         if(.not. lfinished) then
            write(*,'("Calling switch_seed...")')
            write(*,*)"Eref(old)=", eref
            write(*,*)"Sigma    =", sigma_max
            write(*,*)"Eref(new)=", eref+sigma_max
            write(*,*)"Max(err) =", max_rnorm
            write(*,*)"CGSTEP(old)=", CGSTEP
            write(*,*)"# of left mesh points", nfinished
            shift=sigma_max





            call switch_seed(nd,CGSTEP,1,1,eref,shift,alpha,beta,yr,
     $           eps,rnorm,rvecs,yk,interval=1)

            write(*,*)"CGSTEP(new)=", CGSTEP
                                !eref=eref+shift
c     debug
c           NofSSW=NofSSW+1
         endif
      end do
c     debug
      if(max_rnorm .ge. eps) then
         write(*,'("Warning!: Max(Err)=",ES21.13," > eps")') max_rnorm
         write(*,'(" Relax the precision constraint or"
     $        ,    " increase the maximum CGSTEP.")')
         write(*,'(" Choosing another seed (reference energy)"
     $        ,    " may also resolve the problem.")')
      end if

      do i=1, nmesh
         omega=cmplx(emin+(i-1)*denergy,gamma,knd)
         shift=omega-eref

         resolvent=Gval(i)
         err      =Gerr(i)
         write(IUNIT,'(3ES21.13)') dble(omega),
     $        dble(resolvent),imag(resolvent)

         dos= -imag(resolvent)/pi
         if(i==1 .or. i==nmesh)then
            integrated_dos=integrated_dos+dos/2
         else
            integrated_dos=integrated_dos+dos
         end if
         write(IUNT2,'(4ES21.13)') dble(omega), dos,
     $        max(err/pi,9d-99), integrated_dos*denergy
                  !max() here is a quick hack for gnuplot.
                  !gnuplot cannot read floting point without E or D,
                  !for example,  "1.0-108".
c !NSSW for debug
c         write(IUNT2,'(4ES21.13,I4)') dble(omega), dos,
c     $        max(err/pi,9d-99), integrated_dos*denergy,NSSW(i)
      end do
      write(IUNIT,*)
      write(IUNT2,*)
      close(IUNIT)
      close(IUNT2)
      end subroutine write_green
      end program Green
