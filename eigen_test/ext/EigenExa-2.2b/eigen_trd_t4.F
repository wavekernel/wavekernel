       subroutine  eigen_trd_u(
     &               a, nm,
     &               u_x, u_y, nv,
     &               u_t, v_t, i, beta, e)
!$     use OMP_LIB
       use eigen_devel
       use eigen_libs
       implicit NONE

       integer, intent(in)    ::  nm, nv, i
       real(8), intent(inout) ::  a(1:nm)
       real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) ::  u_t(*), v_t(*)
       real(8), intent(out)   ::  beta
       real(8), intent(out)   ::  e(*)

       real(8)                ::  anorm2, a_n, g_n
       real(8)                ::  tt(4), ss(4), t, s

       integer                ::  x_owner_nod, x_pos
       integer                ::  y_owner_nod, y_pos
       integer                ::  j_1, j_2, j_3
       integer                ::  jj_1, jj_2, jj_3
       integer                ::  L


       L = i-1

       x_owner_nod = eigen_owner_node (L, x_nnod, x_inod)
       y_owner_nod = eigen_owner_node (i, y_nnod, y_inod)

       x_pos       = eigen_translate_g2l(L, x_nnod, x_inod)
       y_pos       = eigen_translate_g2l(i, y_nnod, y_inod)

       j_2         = eigen_loop_start(1, x_nnod, x_inod)
       j_3         = eigen_loop_end  (L, x_nnod, x_inod)
!
! u=...
!
       if ( y_owner_nod == y_inod ) then
          anorm2 = ZERO
          do  j_1=j_2,j_3
             t = a(j_1)
             anorm2   =  anorm2 + t**2
             u_x(j_1) =  t
          end do ! j_1
          u_x(x_pos+1) =  anorm2
       end if

       call bcast_dbl(u_x(1), x_pos+1, y_owner_nod, y_COMM_WORLD)

       anorm2 = u_x(x_pos+1)
       u_x(x_pos+1) = ZERO

       if ( x_owner_nod == x_inod ) then
          a_n =  u_x(x_pos)
       else
          a_n =  ZERO
       end if

       tt(1) =  anorm2
       tt(2) =  a_n
       call reduce_dbl(tt(1), ss(1), 2, 1, x_COMM_WORLD)
       anorm2 =  tt(1)
       a_n    =  tt(2)

       if ( anorm2 /= ZERO ) then

          g_n   = -SIGN(SQRT(anorm2), a_n)
          beta  =  anorm2 - a_n * g_n
          e (i) =  g_n

          if ( x_owner_nod == x_inod ) then
             u_x(x_pos) =  a_n - g_n
             if ( y_owner_nod == y_inod ) then
                a(x_pos) =  u_x(x_pos)
             end if
          else
             if ( j_3 < x_pos ) then
                u_x(j_3+1:x_pos) = ZERO
             end if
          end if

          call datacast_dbl(u_y(1), u_x(1), u_t(1), v_t(1), x_pos)

       else

          beta  =  ONE
          e (i) =  ZERO

       end if

       return
       end subroutine  eigen_trd_u

