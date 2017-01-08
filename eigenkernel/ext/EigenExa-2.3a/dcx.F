      subroutine eigen_dcx(n, d, e, lde, z, ldz, INFO, RET)
!$    use OMP_LIB
      use MPI
      use eigen_devel
      use eigen_libs
      use eigen_blacs, only : eigen_get_blacs_context
      use eigen_dc
      implicit NONE

      integer, intent(IN)    :: n, lde, ldz
      real(8), intent(INOUT) :: d(1:n), e(1:2*lde)
      real(8), intent(OUT)   :: z(1:ldz,*)
      integer, intent(OUT)   :: info
      real(8), intent(OUT)   :: ret

! Parameters BLACS array descritor(the position of entry tags), etc
      INTEGER, PARAMETER     :: BLOCK_CYCLIC_2D = 1
      INTEGER, PARAMETER     :: DLEN_  = 9
      INTEGER, PARAMETER     :: DTYPE_ = 1
      INTEGER, PARAMETER     :: CTXT_  = 2
      INTEGER, PARAMETER     :: M_     = 3
      INTEGER, PARAMETER     :: N_     = 4
      INTEGER, PARAMETER     :: MB_    = 5
      INTEGER, PARAMETER     :: NB_    = 6
      INTEGER, PARAMETER     :: RSRC_  = 7
      INTEGER, PARAMETER     :: CSRC_  = 8
      INTEGER, PARAMETER     :: LLD_   = 9

      integer                :: DESCZ( DLEN_ )
      integer                :: DESCW( DLEN_ )
 
      integer                :: i,j,k,nx,istat,NB,lddz,lddw
      integer                :: NP, NQ, NPCOL, NPROW
      integer                :: NPROCS, IAM, MYCOL, MYROW
      integer                :: ICTXT

      integer                :: world_size, my_rank, ierr
      integer                :: TRILWMIN, LWORK, LIWORK
 
      real(8), pointer       :: work(:)
      integer, pointer       :: iwork(:)

      logical, parameter     :: USE_MY_REDIST = .TRUE.

      external               :: DESCINIT
      external               :: DLAED6_INIT, MY_PDSXEDC
      external               :: dc_redist1, dc_redist2
#if defined(__INTEL_COMPILER)
      integer                :: MKL_MODE
      integer, external      :: MKL_GET_DYNAMIC
#endif
      integer, external      :: NUMROC

#if TIMER_PRINT
      real(8)  :: d1,d2
#endif

 
      flops = 0D0
      dgemm_time = 0D0
 
      call eigen_timer_reset( 1, 1, 1, 1 )

#if TIMER_PRINT>1
      d1 = eigen_get_wtime()
#endif

      NPROCS = TRD_nnod
      IAM    = TRD_inod-1

      NPROW = x_nnod
      NPCOL = y_nnod

      ICTXT = eigen_get_blacs_context( )

      MYROW = x_inod-1
      MYCOL = y_inod-1


! BLACS array registration
      NB = eigen_NB
      NB = MIN(NB, N)

      NP = NUMROC( n, NB, MYROW, 0, NPROW )
      NQ = NUMROC( n, NB, MYCOL, 0, NPCOL )
      lddz = (n-1)/NPROW+1
      lddz = ((lddz-1)/NB+1)*NB
      lddw = (n-1)/NPCOL+1
      lddw = ((lddw-1)/NB+1)*NB
!
      call DESCINIT( DESCZ, n, n, NB, NB, 0, 0, ICTXT, lddz, INFO )

! preparing working arrays
      nx     = (N-1)/NPCOL+1
      LWORK  = MAX(1+6*N+2*NP*(NQ+MAX(NQ,NB)), lddz*lddw, ldz*nx)
      LIWORK = 2+7*n+8*NPCOL
      allocate(work(lwork), iwork(liwork), stat=istat)
      if(istat.ne.0) then
           print*,"Memory exhausted"
           call flush(6)
           call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_OTHER, ierr )
      end if

! Somehow, Z must be nullified (Originally next loop is not required.)
!$OMP PARALLEL DO
!OCL XFILL
      do i=1, lddz*lddw
         z(i,1) = ZERO
      enddo
!$OMP END PARALLEL DO
#if TIMER_PRINT>1
      d2 = eigen_get_wtime()
      IF(IAM==0)print*,"before PDSTEDC",d2-d1
#endif

      call DLAED6_INIT()

#if defined(__INTEL_COMPILER)
      MKL_MODE = MKL_GET_DYNAMIC()
      call MKL_SET_DYNAMIC(0)
#endif

#if TIMER_PRINT>1
      d1 = eigen_get_wtime ()
#endif
      CALL MY_PDSXEDC('I',2, n, d(1), e(1), lde, z(1,1), 1, 1, DESCZ,
     $               WORK(1), LWORK, IWORK(1), LIWORK, INFO)
#if TIMER_PRINT>1
      d2 = eigen_get_wtime()
#endif

#if defined(__INTEL_COMPILER)
      call MKL_SET_DYNAMIC(MKL_MODE)
#endif

#if TIMER_PRINT>1
      IF(IAM==0)print*,"PDSTEDC",d2-d1
#endif

      IF ( NB == 1 .OR. TRD_nnod == 1 ) THEN
         if ( lddz /= ldz ) then
         k = MIN(lddz,ldz)
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
         do i=1, lddz*lddw
            work(i) = z(i,1)
         enddo
!$OMP ENDDO
!$OMP BARRIER
!$OMP DO
         do i=1,nx
            j=(i-1)*lddz
            z(1:k,i)=work(1+j:k+j)
         enddo
!$OMP ENDDO
         if ( lddz < ldz ) then
!$OMP BARRIER
!$OMP DO
         do i=nx,1,-1
            z(lddz+1:ldz,i) = ZERO
         enddo
!$OMP ENDDO
         endif
!$OMP END PARALLEL
         endif
         call bcast_dbl( d(1), n, 1, 1, TRD_COMM_WORLD )
      ELSE

      IF( USE_MY_REDIST )THEN
#if TIMER_PRINT>1
         d1 = eigen_get_wtime()
#endif
         call dc_redist1C( n, NB, z, work, lddz, iwork, liwork/2 )
#if TIMER_PRINT>1
         d2 = eigen_get_wtime()
         IF(IAM==0)print*,"MY-REDIST1",d2-d1
         d1 = eigen_get_wtime()
#endif
         call dc_redist2C( n, NB, work, lddz, z, ldz, iwork, liwork/2 )
#if TIMER_PRINT>1
         d2 = eigen_get_wtime()
         IF(IAM==0)print*,"MY-REDIST2",d2-d1
#endif
      ELSE
#if TIMER_PRINT>1
         d1 = eigen_get_wtime()
#endif
         call DESCINIT( DESCW, n, n, 1, 1, 0, 0, ICTXT, ldz, INFO )
         call PDGEMR2D( n, n, z, 1, 1, DESCZ, work, 1, 1, DESCW, ICTXT )
#if TIMER_PRINT>1
         d2 = eigen_get_wtime()
         IF(IAM==0)print*,"PDGEMR2D",d2-d1
#endif
      END IF ! if ( USE_MY_REDIST ) .. else ..

#if TIMER_PRINT>1
         d1 = eigen_get_wtime()
#endif
!$OMP PARALLEL PRIVATE(i)
!$    IF ( omp_get_thread_num()==0 ) THEN
         call bcast_dbl( d(1), n, 1, 1, TRD_COMM_WORLD )
!$    ENDIF
!$    IF ( omp_get_num_threads()==1 .OR.
!$   $     omp_get_thread_num()==1 ) THEN
      IF( .NOT.USE_MY_REDIST )THEN
         do i=1,ldz*nx
            z(i,1)=work(i)
         enddo
      ENDIF
!$    ENDIF
!$OMP END PARALLEL
#if TIMER_PRINT>1
         d2 = eigen_get_wtime()
         IF(IAM==0)print*,"RERE1",d2-d1
#endif
      END IF ! if ( NB == 1 ) .. else ..

      call eigen_timer_print( 'EigenExa (dcx)' )

! freeing working arrays
      deallocate(work)
      deallocate(iwork)

#if TIMER_PRINT>1
      IF(IAM==0)THEN
         print*,"PDGEMM", dgemm_time, 
     &          1D-9*flops/dgemm_time,"GFLOPS",flops
      ENDIF
#endif

      ret = flops

      return
      end subroutine  eigen_dcx

