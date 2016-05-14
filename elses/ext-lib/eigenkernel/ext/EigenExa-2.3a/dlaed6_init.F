      SUBROUTINE DLAED6_INIT( )
      use dlaed6_common
      implicit NONE
*
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(8), PARAMETER  :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                       THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0
*     ..
*     .. External Functions ..
      REAL(8), EXTERNAL   :: DLAMCH
*     ..
*     .. Local Scalars ..
      REAL(8)             :: BASE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC           :: INT, LOG
*     ..
*     .. Executable Statements ..
*
!$OMP CRITICAL
      IF( FIRST ) THEN
         EPS = DLAMCH( 'Epsilon' )
         BASE = DLAMCH( 'Base' )
         SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) /
     $            THREE ) )
         SMINV1 = ONE / SMALL1
         SMALL2 = SMALL1*SMALL1
         SMINV2 = SMINV1*SMINV1
         FIRST = .FALSE.
      END IF
!$OMP END CRITICAL
*
      RETURN
*
*     End of DLAED6_INIT
*
      END SUBROUTINE DLAED6_INIT

