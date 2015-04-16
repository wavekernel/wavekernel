       subroutine  eigen_trd_u(
     &               a, nm,
     &               u_x, u_y, nv,
     &               u_t, v_t, i, kk, v_x, v_y)
!$     use OMP_LIB
       use eigen_devel
       use eigen_libs
       implicit NONE

       integer, intent(in)    ::  nm, nv, i, kk
       real(8), intent(inout) ::  a(1:nm)
       real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) ::  u_t(*), v_t(*)
       real(8), intent(inout) ::  v_x(*)
       real(8), intent(out)   ::  v_y(*)

       real(8)                ::  anorm2, a_n, g_n
       real(8)                ::  tt(4), ss(4), t, s

       integer                ::  x_owner_nod, x_pos
       integer                ::  y_owner_nod, y_pos
       integer                ::  j_1, j_2, j_3
       integer                ::  jj_1, jj_2, jj_3
       integer                ::  L, n


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
          u_x(j_2:j_3) = a(j_2:j_3)
       end if

       call bcast_dbl(u_x(1), x_pos, y_owner_nod, 1, y_COMM_WORLD)

       if ( kk == 0 ) then
       x_pos       = eigen_translate_g2l(L, x_nnod, x_inod)
       call datacast_dbl(u_y(1), u_x(1), u_t(1), v_t(1), x_pos, 1)
       else
       x_pos       = eigen_translate_g2l(i, x_nnod, x_inod)
       call datacast_dbl2(u_y(1), v_y(1), u_x(1), v_x(1),
     &                    u_t(1), v_t(1), x_pos, 1)
       end if

       call eigen_vector_zeropad_x( u_x(1), L )
       call eigen_vector_zeropad_y( u_y(1), L )
       if ( kk /= 0 ) then
       call eigen_vector_zeropad_y( v_y(1), i )
       end if


       return
       end subroutine  eigen_trd_u

