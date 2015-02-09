         subroutine eigen_sx( n, nvec, a, lda, w, z, ldz,
     &                        m_forward, m_backward, mode )
!$       use OMP_LIB
         use MPI
         use eigen_devel
         implicit NONE
*-
         integer, intent(in)    :: n, nvec, lda, ldz
         integer, intent(in), optional :: m_forward, m_backward
         character, intent(in), optional :: mode
         real(8), intent(inout) :: a(lda,*)
         real(8), intent(out)   :: w(*), z(ldz,*)
         real(8), pointer       :: d(:), e(:), e2(:)
*-
         real(8)                :: hs0, hs1, s0, s1
         real(8)                :: ret, ret_1, ret_2, ret_3
         integer                :: m0, nb, nme
         integer                :: i
         integer                :: my_rank, world_size, INFO

#define	DEBUG	0

         real(8)                :: SIGMA
*-
*-------------------------------------------------------------------
*-
         if ( TRD_COMM_WORLD == MPI_COMM_NULL ) then
            return
         end if
*-
*-
*-------------------------------------------------------------------
*-
         hs0 = eigen_get_wtime()
         ret_1 = 0D0; ret_2 = 0D0; ret_3 = 0D0
*-
*-------------------------------------------------------------------
*-
         SIGMA = ONE
!         if ( PRESENT( mode ) ) then
!         if ( mode == 'X' ) then
            call eigen_scaling( n, a(1,1), lda, SIGMA )
!         end if
!         end if
*-
*-------------------------------------------------------------------
*-
         nme = ((n-1)/2+1)*2
         allocate( d(1:n), e(1:nme*2), e2(1:nme*2) )
*-
         world_size = TRD_nnod
         my_rank    = TRD_inod-1
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         s0 = eigen_get_wtime( )
#endif

        if ( PRESENT( m_forward ) ) then
            m0 = m_forward
        else
            m0 = 48
        end if

            call eigen_prd( n, a(1,1), lda, d(1), e(1), nme, m0 )
            ret_1 = dble(n)**3*4/3

#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if ( my_rank == 0 ) then
            print*,"TRD-BLK", n, s1-s0, 1D-9*ret_1/(s1-s0), "GFLOPS"
            print*,"TRD-BLK-INFO",n," ",m0
         end if
         call flush(6)
#endif
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         s0 = eigen_get_wtime( )
#endif


         if ( PRESENT( mode ) ) then
         if ( mode == 'A' .OR. mode == 'X' ) then
            e2(0*nme+1:0*nme+N-1) = e(0*nme+2:0*nme+N)
            e2(0*nme+N) = 0
            e2(1*nme+1:1*nme+N-2) = e(1*nme+3:1*nme+N)
            e2(1*nme+N-1) = 0
            e2(1*nme+N) = 0
            w(1:n)=d(1:n)
            call eigen_dcx( n, w, e2, nme, z, ldz, INFO, ret_2 )
#if DEBUG
            e2(1:n)=w(1:n)
#endif
         end if
         if ( mode == 'N' ) then
            call eigen_bisect2( d, e, e(1+nme), w, n, 0 )
         end if
         if ( mode == 'X' ) then
            call eigen_bisect2( d, e, e(1+nme), w, n, 1 )
         end if
         else
            e2(0*nme+1:0*nme+N-1) = e(0*nme+2:0*nme+N)
            e2(0*nme+N) = 0
            e2(1*nme+1:1*nme+N-2) = e(1*nme+3:1*nme+N)
            e2(1*nme+N-1) = 0
            e2(1*nme+N) = 0
            w(1:n)=d(1:n)
            call eigen_dcx( n, w, e2, nme, z, ldz, INFO, ret_2 )
#if DEBUG
            e2(1:n)=w(1:n)
#endif
         end if


#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if ( my_rank == 0 ) then
            print*,"D&C", s1-s0, 1D-9*ret_2/(s1-s0), "GFLOPS"
#if DEBUG
         if ( mode == 'N' .OR. mode == 'X' ) then
            print*, "D & C=", e2(1:n)
            print*, "STURM=", w(1:n)
            print*, "DC-STURM=", ABS(e2(1:n)-w(1:n))
         end if
#endif
         end if
         call flush(6)
#endif
         if ( PRESENT( mode ) ) then
         if ( mode == 'N' ) then
            goto 99999
         end if
         end if
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         s0 = eigen_get_wtime( )
#endif

        if ( PRESENT( m_backward ) ) then
            m0 = m_backward
        else
            m0 = 128
        end if

            nb = MBAND

            call eigen_trbakwyx( n, nvec, a(1,1), lda, z(1,1), ldz,
     $                   e(1+nme), m0, nb )
            ret_3 = 2*dble(nvec)*dble(n)**2

#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if ( my_rank == 0 ) then
            print*,"TRDBAK", n, s1-s0, 1D-9*ret_3/(s1-s0), "GFLOPS"
         end if
         call flush(6)
#endif
*-
*-------------------------------------------------------------------
*-
         IF ( SIGMA /= ONE ) THEN
            SIGMA = ONE / SIGMA
            call dscal( n, w, SIGMA, 1 )
         END IF
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         hs1 = eigen_get_wtime ( )
         ret = ret_1 + ret_2 + ret_3
         if ( my_rank == 0 ) then
            print*,"Total", hs1-hs0, 1D-9*ret/(hs1-hs0), "GFLOPS"
         end if
         call flush(6)
#endif
*-
*-------------------------------------------------------------------
*-
99999    continue
*-
*-------------------------------------------------------------------
*-
         hs1 = eigen_get_wtime ( )
         ret = ret_1 + ret_2 + ret_3
         a(1, 1) = ret
         a(2, 1) = hs1-hs0
*-
*-------------------------------------------------------------------
*-
         deallocate( d, e, e2 )
*-
*-------------------------------------------------------------------
*-
         end subroutine  eigen_sx

