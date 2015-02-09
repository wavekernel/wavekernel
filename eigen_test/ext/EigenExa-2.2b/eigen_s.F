         subroutine eigen_s( n, nvec, a, lda, w, z, ldz,
     $                       m_forward, m_backward, mode )
!$       use OMP_LIB
         use MPI
         use eigen_devel
         implicit NONE
*-
         integer, intent(in)    :: n, nvec, lda, ldz
         real(8), intent(inout) :: a(lda,*)
         real(8), intent(out)   :: w(*), z(ldz,*)
         integer, intent(in), optional :: m_forward, m_backward
         character, intent(in), optional :: mode
*-
         real(8), pointer       :: d(:), e(:), e2(:)
*-
         real(8)                :: hs0, hs1, s0, s1
         real(8)                :: ret, ret_1, ret_2, ret_3
         integer                :: m0, nb, my_rank, world_size
         integer                :: INFO, ierr, i
*-
         real(8)                :: SIGMA
*-
*-------------------------------------------------------------------
*-
         if ( TRD_COMM_WORLD == MPI_COMM_NULL ) then
            return
         end if
*-
*-------------------------------------------------------------------
*-
         hs0 = eigen_get_wtime( )
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
         print*,"Scaling factor = ",SIGMA
*-
*-------------------------------------------------------------------
*-
         allocate( d(1:n), e(1:n), e2(1:n) )
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

            call eigen_trd( n, a(1,1), lda, d(1), e(1), m0 )
            ret_1 = dble(n)**3*4/3

#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if(my_rank==0)then
            print*,"TRD-BLK", n, s1-s0, 1D-9*ret_1/(s1-s0), "GFLOPS"
         endif
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
            e2(1:n) = e(1:n)
            w (1:n) = d(1:n)
            call eigen_dc2( n, w(1), e2(2), z(1,1), ldz, INFO, ret_2 )
         end if
         if ( mode == 'N' ) then
            call eigen_bisect( d(1), e(1), w(1), n, 0 )
         end if
         if ( mode == 'X' ) then
            call eigen_bisect( d(1), e(1), w(1), n, 1 )
         end if
         else
            e2(1:n) = e(1:n)
            w (1:n) = d(1:n)
            call eigen_dc2( n, w(1), e2(2), z(1,1), ldz, INFO, ret_2 )
         end if


#if TIMER_PRINT
         s1 = eigen_get_wtime( )
         if ( my_rank == 0 ) then
            print*,"D&C", n, s1-s0, 1D-9*ret_2/(s1-s0), 'GFLOPS'
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

            nb = 1
            call eigen_trbakwyx( n, nvec, a(1,1), lda, z(1,1), ldz,
     $                   e, m0, nb )
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
         IF( SIGMA /= ONE ) THEN
            SIGMA = ONE / SIGMA
            call dscal( n, w, SIGMA, 1 )
         END IF
*-
*-------------------------------------------------------------------
*-
#if TIMER_PRINT
         hs1 = eigen_get_wtime( )
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
         hs1 = eigen_get_wtime( )
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
         end subroutine  eigen_s

