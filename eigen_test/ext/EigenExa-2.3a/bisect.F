       subroutine eigen_bisect( d, e, w_out, n, mode )
*
       use MPI
       use eigen_devel
       implicit NONE

*
       integer, intent(IN)    ::  n
       real(8), intent(IN)    ::  d(1:n)
       real(8), intent(IN)    ::  e(1:n)
       real(8), intent(INOUT) ::  w_out(1:n)
       real(8), pointer       ::  wk(:)
       real(8), pointer       ::  w(:)
       integer, intent(IN)    ::  mode
       integer                ::  nv

#if TIMER_PRINT>1
       real(8)                ::  timer
#endif
*
#if TIMER_PRINT>1
       timer = eigen_get_wtime( )
#endif
*
       nv = (n-1)/TRD_nnod+1
*
       allocate ( wk(1:nv*TRD_nnod+2) )
       allocate ( w(1:nv) )

       call bisect_body( d, e, w_out, n, wk, w, mode )

       deallocate ( wk )
       deallocate ( w )
*
#if TIMER_PRINT>1
       timer = eigen_get_wtime( ) - timer
       if ( TRD_inod == 1 ) print*, "BISECT ", timer
#endif
       return
       contains
*
       subroutine bisect_body( d, e, w_out, n, wk, w, mode )
!$     use OMP_LIB
       use eigen_devel
       use eigen_libs
       implicit NONE

#define	DEBUG	0

       integer, intent(IN)    ::  n
       real(8), intent(IN)    ::  d(1:n)
       real(8), intent(IN)    ::  e(1:n)
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
       integer, parameter     ::  ITR_MAX = 128
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
!$OMP PARALLEL DO PRIVATE(u, v, t)
!$OMP+            REDUCTION(MAX:x_upper) REDUCTION(MIN:x_lower)
       do i = 1, n
          u = ABS(e(i))
          t = u
          if ( i == n ) then
             u = ZERO
          else
             u = ABS(e(i+1))
          end if
          t = t + u
          x_lower = MIN(x_lower, d(i)-t)
          x_upper = MAX(x_upper, d(i)+t)
       end do
!$OMP END PARALLEL DO
*
       t = ZERO
!$OMP PARALLEL DO REDUCTION(MAX:t)
       do  i = 1, n
          t = MAX(t, ABS(e(i)) )
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
             s = sturm( w_out(i), d, e, n, m_epsilon )
             w(i-i_2+1) = DBLE(s)
          end do
!$OMP END PARALLEL DO
       end if
       call allgather_dbl( w(1), wk(2), nv, 1, TRD_COMM_WORLD )
       s = sturm( lb_, d, e, n, m_epsilon )
       k = sturm( ub_, d, e, n, m_epsilon )
       wk(1  ) = DBLE(s)
       wk(n+2) = DBLE(k)
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
                lb = w_out(k-1)
                ub = ub_
             else
                lb = w_out(k-1)
                ub = w_out(k  )
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
                ub = w_out(1)
             else
                lb = w_out(k)
                ub = w_out(k+1)
             end if
          end if
*
          if ( mode == 1 ) then
             x = MAX(lb, ub - epsilon)
             if ( x > lb ) then
                s = sturm( x, d, e, n, m_epsilon )
                if ( s < sturm_val ) then
                   lb = x
                else
                   ub = x
                end if
             end if
             x = MIN(ub, lb + epsilon)
             if ( x < ub ) then
                s = sturm( x, d, e, n, m_epsilon )
                if ( s < sturm_val ) then
                   lb = x
                else
                   ub = x
                end if
             end if
          end if
*
          x = lb
          do itr=1,ITR_MAX
             t = x
             x = (lb+ub)/2
             if ( x == t ) exit
             s = sturm( x, d, e, n, m_epsilon )
             if ( s < sturm_val ) then
                lb = x
             else
                ub = x
             end if
          end do
*
#if DEBUG
          theta = PAI*(2*(n-j)+1)/(2*n+1)
          u = 5D-1/(1D0-COS(theta))
*-
          print*,j,":",w_out(j),"->",x,"/",u,ABS(x-u)/ABS(x)
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
       end subroutine bisect_body
!------------------------------------------------------------------
       integer function sturm( x, d, e, n, epsilon )
     &         result( ret )
*
       implicit none
*
       integer, intent(in)    ::  n
       real(8), intent(in)    ::  x
       real(8), intent(in)    ::  d(1:n)
       real(8), intent(in)    ::  e(1:n)
       real(8), intent(in)    ::  epsilon
*
       integer                ::  i, j
       real(8)                ::  t

       real(8)                ::  H2,H3,H4, d0,e0
*
          H2 = ONE; H3 = ONE
*
       j = 0
!OCL UXSIMD
       do i = 1, n
          d0 = (d(i)-x)
          e0 = e(i-1+1)
          H4 = (d0*H3 - (e0**2)*H2)
          if ( ABS(H4)  > epsilon * (ABS(d0)+ABS(e0)) ) then
              t = ONE / H4
              H2 = H3 * t
              H3 = ONE
              if ( H2 < 0 ) j = j + 1
          else
              H2 = H3
              H3 = epsilon
          end if
       end do
          ret = j
*
          return
       end function sturm
!------------------------------------------------------------------

       end subroutine  eigen_bisect

