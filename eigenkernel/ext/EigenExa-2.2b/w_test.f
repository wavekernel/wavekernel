       subroutine  w_test( n, w, mtype )
       use MPI
       use eigen_libs
       implicit NONE

       integer, intent(IN)    :: n, mtype
       real(8)                :: w(1:n)

       real(8), parameter     :: ZERO = 0.0D0
       real(8), parameter     :: ONE  = 1.0D0

       real(8)                :: PAI, EPS, EPS2, EPS4
       real(8)                :: ax, bx, x, y, z, theta, s, t
       real(8), pointer       :: ww(:)
       integer, pointer       :: iseed(:)

       integer                :: i, j, k, kk, l


          if ( mtype < 1 .OR. mtype > 10 ) then
          return
          end if

          if ( mtype == 1 .OR. mtype == 2 ) then
          print*,"-----------------------------------------------"
          print*,"*** Eigenvalue Error Test *** : SKIP"
          print*,"Because the eigenvalues are not given by user"
          print*,"or it is hard to solve them analytically."
          print*,"-----------------------------------------------"
          return
          end if

          PAI  = 4*ATAN(ONE)
          EPS  = machine_epsilon()
          EPS2 = SQRT(EPS)
          EPS4 = SQRT(EPS2)

          allocate ( ww(1:n) )
          call w_set( n, ww, mtype )

          do i = 1, n-1
          do j = i+1, n
             if ( ww(i) > ww(j) ) then
                x = ww(i); ww(i) = ww(j); ww(j) = x
             end if
          end do
          end do

          ax=0D0; k  = 1
          bx=0D0; kk = 1
          do i = 1, n
             x = ww(i)
             y = ABS(w(i)-x)
             if ( x == ZERO ) then
                z = ZERO ! we do not check the relative error in this case
             else
                z = y / ABS(x)
             end if
             if ( ax < z ) then
                ax = MAX(ax, z)
                k = i
             end if
             if ( bx < y ) then
                bx = MAX(bx, y)
                kk = i
             end if
          end do

          deallocate ( ww )


          print*,"-----------------------------------------------"

          print*, "max|w(i)-w(i).true|/|w.true|=", ax, w(k)
          if ( ax < EPS2 ) then
          print*, "*** Eigenvalue Relative Error *** : PASSED"
          else if ( ax < EPS4 ) then
          print*, "*** Eigenvalue Relative Error *** : CAUTION"
          else
          print*, "*** Eigenvalue Relative Error *** : FAILED"
          if ( ABS(w(k)) < EPS ) then
          print*, " |w| is too small, so it is not severe."
          end if
          end if

          print*, "max|w(i)-w(i).true|         =", bx, w(kk)
          if ( bx < EPS2 ) then
          print*, "*** Eigenvalue Absolute Error *** : PASSED"
          else if ( bx < EPS4 ) then
          print*, "*** Eigenvalue Absolute Error *** : CAUTION"
          else
          print*, "*** Eigenvalue Absolute Error *** : FAILED"
          if ( ax < EPS4 ) then
          print*, " Do not mind it. Relative error is small enough."
          else
          print*, " Both absolute test and relative test are failed."
          print*, " Since it is severe, and some bugs might be hidden."
          print*, " Please inform this result to the developper soon."
          end if
          end if

          print*,"-----------------------------------------------"


       return
       end subroutine  w_test

