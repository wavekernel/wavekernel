      module eigen_dc
      implicit NONE

       real(8) ::  flops, dgemm_time, dgemm_dummy(1:2)
       real(8) ::  p_time0 = 0D0, p_time2 = 0D0, p_time3 = 0D0
       real(8) ::  p_times = 0D0, p_timez = 0D0, p_timer = 0D0

       integer, parameter :: DC_NB    = 64
       integer, parameter :: BLK_SIZE = 64

       contains

       logical function get_ID( IWORK, I, N, NB, N1, ID, MATSIZ )
     $         result( RET )
       implicit NONE
       integer, intent(IN)   :: IWORK(*), I, N, NB
       integer, intent(OUT)  :: N1, ID, MATSIZ

       integer :: NBL, NBL1

       RET = .FALSE.

       IF( I.EQ.0 ) THEN
          NBL = IWORK( 2 )
          NBL1 = IWORK( 1 )
          IF( NBL1.EQ.0 ) RETURN
          ID = 1
          MATSIZ = MIN( N, NBL*NB )
       ELSE
          NBL = IWORK( I+2 ) - IWORK( I )
          NBL1 = NBL / 2
          IF( NBL1.EQ.0 ) RETURN
          ID = IWORK( I )*NB + 1
          MATSIZ = MIN( NB*NBL, N-ID+1 )
       END IF
       N1 = NBL1*NB

       RET = .TRUE.

       return
       end function

      end module eigen_dc

