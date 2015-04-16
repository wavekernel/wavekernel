      subroutine eigen_dc2(n, d, e, z, ldz, INFO, RET)
!$    use omp_lib    
      use MPI
      use eigen_devel
      use eigen_libs
      use eigen_blacs, only : eigen_get_blacs_context
      use eigen_dc
      implicit NONE

      integer, intent(in)           :: n, ldz
      real(8), intent(inout)        :: d(1:n), e(1:n)
      real(8), intent(out)          :: z(1:ldz,*)
      integer, intent(out)          :: info
      real(8), intent(out)          :: ret

! Parameters BLACS array descritor(the position of entry tags), etc
      INTEGER, PARAMETER            :: BLOCK_CYCLIC_2D = 1
      INTEGER, PARAMETER            :: DLEN_  = 9
      INTEGER, PARAMETER            :: DTYPE_ = 1
      INTEGER, PARAMETER            :: CTXT_  = 2
      INTEGER, PARAMETER            :: M_     = 3
      INTEGER, PARAMETER            :: N_     = 4
      INTEGER, PARAMETER            :: MB_    = 5
      INTEGER, PARAMETER            :: NB_    = 6
      INTEGER, PARAMETER            :: RSRC_  = 7
      INTEGER, PARAMETER            :: CSRC_  = 8
      INTEGER, PARAMETER            :: LLD_   = 9

      integer                       :: DESCZ( DLEN_ )
      integer                       :: DESCW( DLEN_ )

      integer                       :: i,j,k,nx,istat,NB,lddz,lddw
      integer                       :: NP, NQ, NPCOL, NPROW
      integer                       :: NPROCS, IAM, MYCOL, MYROW
      integer                       :: ICTXT
 
      integer                       :: world_size, my_rank, ierr
      integer                       :: TRILWMIN, LWORK, LIWORK
 
      real(8), pointer              :: work(:)
      integer, pointer              :: iwork(:)

      logical, parameter            :: USE_MY_REDIST = .TRUE.

#if defined(__INTEL_COMPILER)
      integer                       :: MKL_MODE
      integer, external             :: MKL_GET_DYNAMIC
#endif

      integer, external             :: NUMROC

#if TIMER_PRINT
      real(8)  :: d1,d2
      real(8)  :: dd1,dd2
#endif


#if TIMER_PRINT
      flops = 0D0
      dgemm_time = 0D0
      p_time0 = 0D0
      p_timer = 0D0
      p_time2 = 0D0
      p_time3 = 0D0
      p_times = 0D0
      p_timez = 0D0
#endif

      call eigen_timer_reset( 1, 0, 0, 0 )

#if TIMER_PRINT
      dd1 = MPI_Wtime()
#endif

      NPROCS = TRD_nnod
      IAM    = TRD_inod-1

      NPROW = x_nnod
      NPCOL = y_nnod

      ICTXT = eigen_get_blacs_context( )
!      ICTXT = BLACS_ICONTXT_FOR_EIGEN

!      CALL BLACS_GRIDINIT( ICTXT, 'Column-major', NPROW, NPCOL )

      MYROW = x_inod-1
      MYCOL = y_inod-1

!      call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )


!      d1 = MPI_Wtime()
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
      call DESCINIT( DESCZ, n, n, NB, NB, 0, 0, ICTXT, lddz, ierr )
!      d2 = MPI_Wtime()
!      IF(IAM==0)print*,"RERE0-1",d2-d1
 
!      d1 = MPI_Wtime()
! preparing working arrays
      nx     = (N-1)/NPCOL+1
      LWORK  = MAX(1+6*N+2*NP*(NQ+MAX(NQ,NB)), lddz*lddw, ldz*nx)
      LIWORK = 2+7*n+8*NPCOL
      allocate(work(lwork), iwork(liwork), stat=istat)
      if(istat.ne.0) then
         print*,"Memory exhausted"
         call flush(6)
         call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
      end if
!      d2 = MPI_Wtime()
!      IF(IAM==0)print*,"RERE0-2",d2-d1


!      d1 = MPI_Wtime()
! Somehow, Z must be nullified (Originally next loop is not required.)
!$OMP PARALLEL DO
      do i=1, lddz*lddw
         z(i,1) = ZERO
      end do
!$OMP END PARALLEL DO
!      d2 = MPI_Wtime()
!      IF(IAM==0)print*,"RERE0-3",d2-d1


#if defined(__INTEL_COMPILER)
      MKL_MODE = MKL_GET_DYNAMIC()
      call MKL_SET_DYNAMIC(0)
#endif
#if TIMER_PRINT>1
      d1 = MPI_Wtime()
#endif
      CALL MX_PDSTEDC( 'I', n, d(1), e(1), z(1,1), 1, 1, DESCZ,
     $              WORK(1), LWORK, IWORK(1), LIWORK, INFO )
#if TIMER_PRINT>1
      d2 = MPI_Wtime()
#endif
#if defined(__INTEL_COMPILER)
      call MKL_SET_DYNAMIC(MKL_MODE)
#endif

#if TIMER_PRINT>1
      IF(IAM==0)print*,"PDSTEDC",d2-d1
#endif


      IF(NB==1.OR.TRD_nnod==1)THEN
#if TIMER_PRINT>1
         d1 = MPI_Wtime()
#endif
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
         do i=1,nx
            z(lddz+1:ldz,i) = ZERO
         enddo
!$OMP ENDDO
         endif
!$OMP END PARALLEL
         endif
         call bcast_dbl( d(1), n, 1, 1, TRD_COMM_WORLD )
#if TIMER_PRINT>1
         d2 = MPI_Wtime()
         IF(IAM==0)print*,"RERE1",d2-d1
#endif
      ELSE

      IF( USE_MY_REDIST )THEN
#if TIMER_PRINT>1
         d1 = MPI_Wtime()
#endif
         call dc_redist1C( n, NB, z, work, lddz, iwork, liwork/2 )
#if TIMER_PRINT>1
         d2 = MPI_Wtime()
         IF(IAM==0)print*,"MY-REDIST1",d2-d1
         d1 = MPI_Wtime()
#endif
         call dc_redist2C( n, NB, work, lddz, z, ldz, iwork, liwork/2 )
#if TIMER_PRINT>1
         d2 = MPI_Wtime()
         IF(IAM==0)print*,"MY-REDIST2",d2-d1
#endif
      ELSE
#if TIMER_PRINT>1
         d1 = MPI_Wtime()
#endif
         call DESCINIT( DESCW, n, n, 1, 1, 0, 0, ICTXT, ldz, ierr )
         call PDGEMR2D( n, n, z, 1, 1, DESCZ, work, 1, 1, DESCW, ICTXT )
#if TIMER_PRINT>1
         d2 = MPI_Wtime()
         IF(IAM==0)print*,"PDGEMR2D",d2-d1
#endif
      ENDIF

!         d1 = MPI_Wtime()
!$OMP PARALLEL private(i)
!$    IF ( omp_get_thread_num()==0 ) THEN
      IF( .NOT.USE_MY_REDIST )THEN
         do i=1, ldz*nx
            z(i,1)=work(i)
         end do
      ENDIF
!$    ENDIF
!$    IF ( omp_get_num_threads()==1 .OR.
!$   $     omp_get_thread_num()==1 ) THEN
         call bcast_dbl( d(1), n, 1, 1, TRD_COMM_WORLD )
!$    ENDIF
!$OMP END PARALLEL
!         d2 = MPI_Wtime()
!         IF(IAM==0)print*,"RERE1",d2-d1
      END IF !

! freeing working arrays
      deallocate(work)
      deallocate(iwork)

 
! BLACS/PBLAS/SCALAPACK finalize
!      call BLACS_GRIDEXIT( ICTXT )

      call eigen_timer_print( 'EigenExa (D-and-C for tri-diag)' )

#if TIMER_PRINT
      dd1 = MPI_Wtime() - dd1
#endif

#if TIMER_PRINT>1
         IF(IAM==0)print*,"DIVIDE.",p_time0
         IF(IAM==0)print*,"PDLASRT",p_timer
         IF(IAM==0)print*,"PDLAED2",p_time2
         IF(IAM==0)print*,"PDLAED3",p_time3
         IF(IAM==0)print*,"PDLAEDZ",p_timez
         IF(IAM==0)print*,"PDLASET",p_times
         IF(IAM==0)print*,"PDGEMM", dgemm_time
#endif

      ret = flops ! dgemm_time ! flops/dgemm_time

      return
      end subroutine

