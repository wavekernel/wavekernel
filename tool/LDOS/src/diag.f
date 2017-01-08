!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
      module diag
      use common_constants

      character,parameter :: WITH_VEC='V', VAL_ONLY='N', UP='U'

      private
      public :: diag_herm, diag_sym, diag_gen

      contains

*     wrapper routine
*     input: integer   :: ndim
*         complex(knd) :: mat(ndim,ndim) ! symmetric matrix
*     output: mat
*            real(knd) :: eig(ndim)
*     diagonalize hermitian matrix "mat", and returns eigenvalues(eig)
*     and eigenvectors(mat)
*     
      subroutine diag_herm(ndim,mat,eig)
      implicit none
      integer,intent(in)        :: ndim
      complex(knd),intent(inout):: mat(ndim,ndim)
      real(knd),intent(out)     :: eig(ndim)
      character :: JOBZ, UPLO
      integer   :: INFO, LDA, LWORK, NM
      complex(knd),allocatable::WK(:)
      real(knd),allocatable::RWORK(:)
      JOBZ=WITH_VEC
      UPLO=UP
      NM=ndim
      LDA=NM
      LWORK=max(1,2*NM-1)
      allocate(WK(LWORK))
      allocate(RWORK(max(1,3*NM-2)))
      call zheev(JOBZ,UPLO,NM,mat,LDA,EIG,WK,LWORK,RWORK,INFO)
      if(INFO.ne.0) stop 'diag_herm: INFO <> 0'
      deallocate(RWORK)
      deallocate(WK)
      end subroutine diag_herm

*     wrapper routine
*     input: integer   :: ndim
*            real(knd) :: mat(ndim,ndim) ! symmetric matrix
*     output: mat
*            real(knd) :: eig(ndim)
*     diagonalize matrix "mat", and returns eigenvalues(eig)
*     and eigenvectors(mat)
*     
      subroutine diag_sym(ndim,mat,eig)
      implicit none
      integer,intent(in)     :: ndim
      real(knd),intent(inout):: mat(ndim,ndim)
      real(knd),intent(out)  :: eig(ndim)
      character:: JOBZ, UPLO
      integer:: INFO, LDA, LWORK, NM
      real(knd),allocatable::WK(:)
      JOBZ=WITH_VEC
      UPLO=UP
      NM=ndim
      LDA=NM
!      LWORK=3*NM-1    ! temp patch for bugs in mkl10.0 
      LWORK=6*NM
      allocate(WK(LWORK))
      call dsyev(JOBZ,UPLO,NM,mat,LDA,EIG,WK,LWORK,INFO)
      if(INFO.ne.0) stop 'diag_sym: INFO <> 0'
      deallocate(WK)
      end subroutine diag_sym


      subroutine diag_gen(ndim,mat,eig)
      implicit none
      integer,intent(in)      :: ndim
      real(knd),intent(inout) :: mat(ndim,ndim)
      complex(knd),intent(out):: eig(ndim)
      integer   :: LDA,LDVL,LDVR,LWORK,INFO, i, j
      real(knd),allocatable,dimension(:)   :: WR, WI, WORK
      real(knd),allocatable,dimension(:,:) :: VL, VR
      character :: JOBVL, JOBVR
      JOBVL='N' ! not calculate left eigen vectors
      JOBVR='V' ! calculate right eigen vectors
      LDA =ndim
      LDVL=1
      LDVR=ndim
      LWORK=max(1,4*ndim)
      allocate(WR(ndim))
      allocate(WI(ndim))
      allocate(VL(LDVL,ndim))
      allocate(VR(LDVR,ndim))
      allocate(WORK(LWORK))
**********FROM LAPACK 3.0 SOURCE**********
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues.  Complex
*          conjugate pairs of eigenvalues appear consecutively
*          with the eigenvalue having the positive imaginary part
*          first.
*
*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
*          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*          after another in the columns of VL, in the same order
*          as their eigenvalues.
*          If JOBVL = 'N', VL is not referenced.
*          If the j-th eigenvalue is real, then u(j) = VL(:,j),
*          the j-th column of VL.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
*          u(j+1) = VL(:,j) - i*VL(:,j+1).
*******************************************
      call DGEEV( JOBVL, JOBVR, ndim, mat, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
      if(INFO.ne.0) stop 'diag_gen: INFO <> 0'
      do i=1, ndim
         eig(i)=cmplx(WR(i),WI(i),knd)
         do j=1, ndim
            mat(j,i)=VR(j,i)
         end do
      end do
      deallocate(WR)
      deallocate(WI)
      deallocate(VL)
      deallocate(VR)
      deallocate(WORK)
      end subroutine diag_gen

      end module diag
