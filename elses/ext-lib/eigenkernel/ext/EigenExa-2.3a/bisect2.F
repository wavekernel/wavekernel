       subroutine eigen_bisect2( d, e, f, w_out, n, mode )
*
       use MPI
       use eigen_devel
       implicit NONE

*
       integer, intent(IN)    ::  n
       real(8), intent(IN)    ::  d(1:n)
       real(8), intent(IN)    ::  e(1:n)
       real(8), intent(INOUT) ::  f(1:n)
       real(8), intent(INOUT) ::  w_out(1:n)
       real(8), pointer       ::  wk(:)
       real(8), pointer       ::  w(:)
       integer, intent(IN)    ::  mode
       integer                ::  nv
#if TIMER_PRINT
       real(8)                ::  timer
#endif
*
#if TIMER_PRINT
       timer = eigen_get_wtime( )
#endif
       nv = (n-1)/TRD_nnod+1
*
       allocate ( wk(1:nv*TRD_nnod+2) )
       allocate ( w(1:nv) )

       call bisect2_body( d, e, f, w_out, n, wk, w, mode )

       deallocate ( wk )
       deallocate ( w )
*
#if TIMER_PRINT
       timer = eigen_get_wtime( ) - timer
       if ( TRD_inod == 1 ) print*,"BISECT2 ", timer
#endif
       return
       contains
*
!------------------------------------------------------------------
!------------------------------------------------------------------
       subroutine bisect2_body( d, e, f, w_out, n, wk, w, mode )
!$     use OMP_LIB
       use eigen_devel
       use eigen_libs
       implicit NONE

#define	DEBUG	0

       integer, intent(IN)    ::  n
       real(8), intent(IN)    ::  d(1:n)
       real(8), intent(IN)    ::  e(1:n)
       real(8), intent(INOUT) ::  f(1:n)
       real(8), intent(INOUT) ::  w_out(1:n)
       real(8), intent(INOUT) ::  wk(*)
       real(8), intent(INOUT) ::  w(*)
       integer, intent(IN)    ::  mode

       real(8)                ::  u, v, x, t
       real(8)                ::  x_lower, x_upper
       real(8)                ::  m_epsilon, epsilon
#if DEBUG
       real(8)                ::  PAI, theta
#endif
       real(8)                ::  lb, ub, lb_, ub_
       integer                ::  i, j, k, s
       integer                ::  i_1,i_2,i_3, nv
       integer                ::  sturm_val
       integer                ::  itr
       integer, parameter     ::  ITRMAX = 128
*
*
       if ( n <= 0 ) return
       if ( n == 1 ) then
          w_out(1) = d(1)
          return
       end if
*
       m_epsilon = machine_epsilon( )
#if DEBUG
       PAI = 4*ATAN(ONE)
#endif
*
       x_lower = d(1)
       x_upper = d(1)
       u = ZERO
       v = ZERO
!$OMP PARALLEL DO PRIVATE(u, v, t)
!$OMP+            REDUCTION(MAX:x_upper) REDUCTION(MIN:x_lower)
       do i = 1, n
          u = ABS(e(i))
          v = ABS(f(i))
          t = u + v
          if ( i == n ) then
             u = ZERO
             v = ZERO
          else if ( i == n-1 ) then
             u = ABS(e(i+1))
             v = ZERO
          else
             u = ABS(e(i+1))
             v = ABS(f(i+2))
          end if
          t = t + u + v
          x_lower = MIN(x_lower, d(i)-t)
          x_upper = MAX(x_upper, d(i)+t)
       end do
!$OMP END PARALLEL DO
*
       t = ZERO
!$OMP PARALLEL DO REDUCTION(MAX:t)
       do  i = 1, n
          t = MAX(t, ABS(e(i))+ABS(f(i)))
       end do! j
!$OMP END PARALLEL DO
       epsilon = m_epsilon * t
*
       x  = (ABS(x_lower)+ABS(x_upper))*m_epsilon
       lb_ = (x_lower - x) - epsilon
       ub_ = (x_upper + x) + epsilon
*
       nv = (n-1)/TRD_nnod+1
       i_2 = nv*(TRD_inod-1)+1
       i_3 = MIN(n, i_2+nv-1)
*
       if ( mode == 0 ) then
!$OMP PARALLEL DO
          do i=1,n
             w_out(i) = (lb_*(n+1-i)+ub_*(i+1-1))/(n+1)
          end do
!$OMP END PARALLEL DO
       end if
*
       if ( i_2 <= i_3 ) then
!$OMP PARALLEL DO PRIVATE(s)
          do i=i_2,i_3
             s = sturm2( w_out(i), d, e, f, n, epsilon )
             w(i-i_2+1) = DBLE(s)
          end do
!$OMP END PARALLEL DO
       end if
       call allgather_dbl( w(1), wk(2), nv, 1, TRD_COMM_WORLD )
       s = sturm2( lb_, d, e, f, n, epsilon )
       k = sturm2( ub_, d, e, f, n, epsilon )
       wk(1  ) = DBLE(s)
       wk(n+2) = DBLE(k)
#if DEBUG
       print*,wk(1:n+2)
       print*,"LB_",lb_
       print*,"UB_",ub_
#endif
*
#if !DEBUG
!$OMP PARALLEL DO PRIVATE(i, j, k, s, sturm_val, t, x, lb, ub)
#endif
       do  i_1 = i_2, i_3
          j = i_1
*
          sturm_val = j
*
!  Find i such that "wk(i) < sturm_val <= wk(i+1)" holds
          s = INT(wk((j+1)))
          if ( s < sturm_val ) then
             k = n+1
             do i=j+1,n
                if ( INT(wk((i+1))) >= sturm_val ) then
                   k = i
                   exit
                end if
             end do
             if ( k == n+1 ) then
                lb = MIN(w_out(k-1), w_out(j))
                ub = ub_
#if DEBUG
          print*,"T1",lb,ub
#endif
             else
                lb = MIN(w_out(k-1), w_out(j))
                ub = MAX(w_out(k), w_out(j))
#if DEBUG
          print*,"T2",lb,ub
#endif
             end if
          else
! sturm_val <= s
             k = 0
             do i=j-1,1,-1
                if ( INT(wk(i+1)) < sturm_val ) then
                   k = i
                   exit
                end if
             end do
! w_out(k) < val <= w_out(k+1)
             if ( k == 0 ) then
                lb = lb_
                ub = MAX(w_out(k+1), w_out(j))
#if DEBUG
          print*,"T3",lb,ub
#endif
             else
                lb = MIN(w_out(k), w_out(j))
                ub = MAX(w_out(k+1), w_out(j))
#if DEBUG
          print*,"T4",lb,ub
#endif
             end if
          end if
#if DEBUG
          print*,j,"++++  :",lb,".",ub,",", w_out(j)
#endif
*
          if ( mode == 1 ) then
             x = MAX(lb, ub - epsilon)
             if ( x > lb ) then
                s = sturm2( x, d, e, f, n, epsilon )
                if ( s < sturm_val ) then
                   lb = x
#if DEBUG
          print*,j,"LB:",lb,".",s,sturm_val
#endif
                else
                   ub = x
#if DEBUG
          print*,j,"UB:",ub,".",s,sturm_val
#endif
                end if
             end if
             x = MIN(ub, lb + epsilon)
             if ( x < ub ) then
                s = sturm2( x, d, e, f, n, epsilon )
                if ( s < sturm_val ) then
                   lb = x
#if DEBUG
          print*,j,"LB:",lb,".",s,sturm_val
#endif
                else
                   ub = x
#if DEBUG
          print*,j,"UB:",ub,".",s,sturm_val
#endif
                end if
             end if
          end if
*-
          x = lb
          do itr=1,ITRMAX
             t = x
             x = (lb+ub)/2
             if ( x == t ) exit
             s = sturm2( x, d, e, f, n, epsilon )
             if ( s < sturm_val ) then
                lb = x
             else
                ub = x
             end if
#if DEBUG
          print*,j,":",lb,".",x,".",ub,s,sturm_val
#endif
          end do
*
#if DEBUG
!          theta = PAI*(2*(n-j)+1)/(2*n+1)
!          u = 5D-1/(1D0-COS(theta))
!*-
          print*,j,":",w_out(j),"->",x
#endif
*
          w(j-i_2+1) = x
*
       end do! i_1
#if !DEBUG
!$OMP END PARALLEL DO
#endif
*
*
       call allgather_dbl( w(1), wk(1), nv, 2, TRD_COMM_WORLD )
       w_out(1:n) = wk(1:n)
*
       return
       end subroutine bisect2_body
!------------------------------------------------------------------
!------------------------------------------------------------------
       integer function sturm2( x, d, e, f, n, epsilon )
     &         result( ret )
*
       implicit none
*
       integer, intent(in)    ::  n
       real(8), intent(in)    ::  x
       real(8), intent(in)    ::  d(1:n)
       real(8), intent(in)    ::  e(1:n)
       real(8), intent(in)    ::  f(1:n)
       real(8), intent(in)    ::  epsilon
*
! Assumption
!  e(1) = f(1) = f(2) = 0
*
           ret = sturm2_LDLT( x, d, e, f, n,  epsilon )
!           ret = sturm2_STRM( x, d, e, f, n,  epsilon )
!
       return
       end function sturm2
!------------------------------------------------------------------
!------------------------------------------------------------------
*
       integer function sturm2_LDLT( x, d, e, f, n, epsilon )
     &         result( ret )
*
       implicit none
*
       integer, intent(in)    ::  n
       real(8), intent(in)    ::  x
       real(8), intent(in)    ::  d(1:n)
       real(8), intent(in)    ::  e(1:n)
       real(8), intent(in)    ::  f(1:n)
       real(8), intent(in)    ::  epsilon

       integer                ::  i, j, k
       real(8)                ::  d0, e0, e1, f0, f1, g0, g1, t

#define  p(i,j)                   p_ i _ j

       real(8)                ::  p_1_1
       real(8)                ::  p_2_1, p_2_2
       real(8)                ::  p_3_1, p_3_2, p_3_3
       real(8)                ::  p_4_1, p_4_2, p_4_3, p_4_4

       logical                ::  trail
*
*
!
!  e(1) = f(1) = f(2) = 0
!
!===============================================================
!
!  We take a diagonal-neibor pivoting strategy
!  in order to stabilize the numerical process to count 
!  the number of negative signs for the diagonal matrix
!  resulted in the LDL decomposition
!
       j = 0
!
       p(1,1) = 1
       p(2,1) = 0; p(2,2) = 1
       p(3,1) = 0; p(3,2) = 0; p(3,3) = 1
       p(4,1) = 0; p(4,2) = 0; p(4,3) = 0; p(4,4) = 1

       trail = .FALSE.
!OCL UXSIMD
       do i = 1, n

       if ( trail ) then

          e0=p(2,1)
          if ( ABS(e0) < epsilon ) e0 = epsilon
          f0=p(3,1)
          g0=p(4,1)
          f1=p(3,2)
          g1=p(4,2)

          p(1,1) = p(3,3) - 2*f0*f1/e0
          p(2,1) = p(3,2) - (g0*f1+g1*f0)/e0
          p(2,2) = p(3,3) - 2*g0*g1/e0
          p(3,1) = f(i-1)
          p(3,2) = e(i-1)
          p(3,3) = d(i-1)-x
          p(4,1) = 0
          p(4,2) = f(i)
          p(4,3) = e(i)
          p(4,4) = d(i)-x

          if ( ABS(e0) > epsilon ) then
             j = j + 1
          else
             j = j + 2
          end if
          trail = .FALSE.

       else

          if ( ABS(p(1,1)) < ABS(p(2,2)) ) then
             t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
             t = p(3,1); p(3,1) = p(3,2); p(3,2) = t
             p(4,1) = p(4,2); p(4,2) = 0 ! since p(4,1) was set to 0
          end if

       if ( p(1,1) /= 0 ) then
!
! If p(1,1)==0, then
!  we can do block LDL decomposition for a 2x2 submatrix.
!
          d0=p(1,1)
          e0=p(2,1)
          f0=p(3,1)
          g0=p(4,1)

          p(1,1) = p(2,2) - e0*e0/d0
          p(2,1) = p(3,2) - e0*f0/d0
          p(2,2) = p(3,3) - f0*f0/d0
          p(3,1) = p(4,2) - e0*g0/d0
          p(3,2) = p(4,3) - f0*g0/d0
          p(3,3) = p(4,4) - g0*g0/d0
          p(4,1) = 0
          p(4,2) = f(i)
          p(4,3) = e(i)
          p(4,4) = d(i)-x

          if ( d0 < 0 ) j = j + 1
          trail = .FALSE.

       else ! p(1,1) ==0 && p(2,2) == 0

          trail = .TRUE.

       end if
       end if

       end do

       if ( trail ) then

          e0=p(2,1)
!
! When p(2,1) is zero, we need another vipoting,
! however present version reduce such cost by making p(2,1) element
! as a small number (machine epsilon).
!
          if ( ABS(e0) < epsilon ) e0 = epsilon
          f0=p(3,1)
          g0=p(4,1)
          f1=p(3,2)
          g1=p(4,2)

          p(1,1) = p(3,3) - 2*f0*f1/e0
          p(2,1) = p(3,2) - (g0*f1+g1*f0)/e0
          p(2,2) = p(3,3) - 2*g0*g1/e0
          p(3,1) = f(n)
          p(3,2) = e(n)
          p(3,3) = d(n)-x

          if ( ABS(e0) > epsilon ) then
             j = j + 1
          else
             j = j + 2
          end if
          trail = .FALSE.

       else

          if ( ABS(p(1,1)) < ABS(p(2,2)) ) then
             t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
             t = p(3,1); p(3,1) = p(3,2); p(3,2) = t
             p(4,1) = p(4,2); p(4,2) = 0 ! since p(4,1) was set to 0
          end if

       if ( p(1,1) /= 0 ) then

          d0=p(1,1)
          e0=p(2,1)
          f0=p(3,1)
          g0=p(4,1)

          p(1,1) = p(2,2) - e0*e0/d0
          p(2,1) = p(3,2) - e0*f0/d0
          p(2,2) = p(3,3) - f0*f0/d0
          p(3,1) = p(4,2) - e0*g0/d0
          p(3,2) = p(4,3) - f0*g0/d0
          p(3,3) = p(4,4) - g0*g0/d0

          if ( d0 < 0 ) j = j + 1
          trail = .FALSE.

       else ! p(1,1) ==0 && p(2,2) == 0

          trail = .TRUE.

       end if
       end if

       if ( trail ) then

          e0=p(2,1)
          if ( ABS(e0) < epsilon ) e0 = epsilon
          f0=p(3,1)
          g0=p(4,1)
          f1=p(3,2)
          g1=p(4,2)

          p(1,1) = p(3,3) - 2*f0*f1/e0
          p(2,1) = p(3,2) - (g0*f1+g1*f0)/e0
          p(2,2) = p(3,3) - 2*g0*g1/e0

          if ( ABS(e0) > epsilon ) then
             j = j + 1
          else
             j = j + 2
          end if
          trail = .FALSE.

       else

          if ( ABS(p(1,1)) < ABS(p(2,2)) ) then
             t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
             t = p(3,1); p(3,1) = p(3,2); p(3,2) = t
          end if

       if ( p(1,1) /= 0 ) then

          d0=p(1,1)
          e0=p(2,1)
          f0=p(3,1)

          p(1,1) = p(2,2) - e0*e0/d0
          p(2,1) = p(3,2) - e0*f0/d0
          p(2,2) = p(3,3) - f0*f0/d0

          if ( d0 < 0 ) j = j + 1
          trail = .FALSE.

       else ! p(1,1) ==0 && p(2,2) == 0

          trail = .TRUE.

       end if
       end if

       if ( trail ) then

          e0=p(2,1)
          if ( ABS(e0) < epsilon ) e0 = epsilon
          f0=p(3,1)
          f1=p(3,2)

          p(1,1) = p(3,3) - 2*f0*f1/e0

          if ( ABS(e0) > epsilon ) then
             j = j + 1
          else
             j = j + 2
          end if
          trail = .FALSE.

       else

          if ( ABS(p(1,1)) < ABS(p(2,2)) ) then
             t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
          end if

       if ( p(1,1) /= 0 ) then

          d0=p(1,1)
          e0=p(2,1)

          p(1,1) = p(2,2) - e0*e0/d0

          if ( d0 < 0 ) j = j + 1
          trail = .FALSE.

       else ! p(1,1) ==0 && p(2,2) == 0

          trail = .TRUE.

       end if
       end if

       if ( trail ) then

          e0 = p(2,1)
          if ( ABS(e0) < epsilon ) e0 = epsilon

          if ( ABS(e0) > epsilon ) then
             j = j + 1
          else
             j = j + 2
          end if

       else

          d0=p(1,1)

          if ( d0 < 0 ) j = j + 1

       end if

#undef  p

             ret =  j
*
          return
       end function sturm2_LDLT
!------------------------------------------------------------------
!------------------------------------------------------------------
       integer function sturm2_STRM( x, d, e, f, n, epsilon )
     &         result( ret )
*
       implicit none
*
       integer, intent(in)    ::  n
       real(8), intent(in)    ::  x
       real(8), intent(in)    ::  d(1:n)
       real(8), intent(in)    ::  e(1:n)
       real(8), intent(in)    ::  f(1:n)
       real(8), intent(in)    ::  epsilon

       integer                ::  i, j
       real(8)                ::  H1,H2,H3,H4
       real(8)                ::  G2
       real(8)                ::  F1
       real(8)                ::  d0, d1, e0, f0, t
*
*
!
!  e(1) = f(1) = f(2) = 0
!
          H1 = ONE; H2 = ONE; H3 = ONE
          G2 = ZERO; F1 = ZERO
*
#if DEBUG
          print*,"D0=",d(1)-x,d(1),x
          print*,"E0=",e(1)
#endif
*
       j = 0
       t = ONE
!OCL UXSIMD
       do i = 1, n
          d0 = d(i)-x
          e0 = e(i-1+1)
          f0 = f(i-2+2)

          H4 = f0 * (e0*G2 - f0*F1)
          G2 =      (e0*H2 - f0*G2)
          H4 = H4 + (d0*H3 - e0*G2)
          F1 =      (d0*H2 - (f0**2)*H1)

#if DEBUG
          print*,i,"===>>>",x, H4, j
#endif

          if ( ABS(H4)  > epsilon * 8*(ABS(d0)+ABS(e0)) ) then
             if ( ABS(H4) > 1D0 ) then
                t = ONE / H4
             else if ( H4 > ZERO ) then
                t = ONE
             else
                t = MONE
             end if
             F1 = F1 * t
             G2 = G2 * t
             H1 = H2 * t
             H2 = H3 * t
             H3 = H4 * t
             if ( H2 < 0 ) j = j + 1
#if DEBUG
          print*,i,">>>",x,H2, j
#endif
          else
             if ( H4 > ZERO ) then
                t = ONE
             else
                t = MONE
             end if
             H1 = H2 * t
             H2 = H3 * t
!             H1 = H2
!             H2 = H3
             H3 = epsilon
             if ( H1 < 0 ) j = j + 1
#if DEBUG
          print*,i,">>>",x,H2,H1,j
#endif
          end if
       end do
          ret = j
*
          return
       end function sturm2_STRM
!------------------------------------------------------------------
!------------------------------------------------------------------
       end subroutine  eigen_bisect2

