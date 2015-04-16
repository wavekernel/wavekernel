       subroutine eigen_trd_vo(
     &               u_x, v_x, v_y,
     &               ux, vx, nv,
     &               u_t, v_t, beta, u_n, v_n,
     &               i, i_base, m0, kk)
!$     use OMP_LIB!
       use eigen_devel
       use eigen_libs
       implicit none
!
       integer, intent(in)    ::  nv, i, i_base, m0, kk
       real(8), intent(in)    ::  u_x(*)
       real(8), intent(out)   ::  v_x(*), v_y(*)
       real(8), intent(in)    ::  ux(1:nv,*)
       real(8), intent(in)    ::  vx(1:nv,*)
       real(8), intent(inout) ::  u_t(*)
       real(8), intent(inout) ::  v_t(*)
       real(8), intent(in)    ::  beta
       real(8), intent(in)    ::  u_n
       real(8), intent(inout) ::  v_n
!
       integer                ::  j, k, L, n, LL
       integer                ::  k_1, k_2, k_3
!
       integer                ::  i_1, i_2, i_3, i_4
       integer                ::  j_1, j_2, j_3, j_4
       integer                ::  l_1, l_2, l_3, l_4
       integer                ::  jj_1, jj_2, jj_3, jj_4
       integer                ::  LX
!
#if __IBM_REGISTER_VARS
       integer, parameter     ::  VTOL = 102400
#else
       integer, parameter     ::  VTOL = 1024
#endif
!
       include 'CSTAB.h'
!
       real(8)                ::  w0
       real(8)                ::  u0, v0
       real(8)                ::  u1, v1
       real(8)                ::  u2, v2
       real(8)                ::  ux0, vx0
       real(8)                ::  ux1, vx1
       real(8)                ::  ux2, vx2
!
       real(8)                ::  alpha, prod_uv
       integer                ::  x_pos, x_owner_nod
       integer                ::  y_pos, y_owner_nod



          k_1 = i - i_base
          k_2 = m0

          L = i-1
          n  = eigen_translate_g2l(L, x_nnod,x_inod)

!$OMP MASTER
          call eigen_vector_zeropad_x( v_x(1), L )
!$OMP END MASTER

!$OMP BARRIER

          prod_uv = u_t(1)

          if ( k_2 <= k_1 ) then

             if ( beta /= ZERO ) then

             alpha = prod_uv/(2*beta)
!$OMP DO
             do j_1=1,n
                v_x(j_1) = (v_x(j_1)-alpha*u_x(j_1))/beta
             end do! j_1
!$OMP ENDDO

             end if

          else
!
! v=v-(UV+VU)u
!
             l_4 = MOD(k_2-k_1, 3)+k_1+1
             LX  = 64 ! L1_LSIZE*L1_WAY/16

             LL = (n-1)/y_nnod+1
             LL = ((LL-1)/2+1)*2

             alpha = ZERO
             do l_1=k_1+1,k_2
                j = l_1-k_1
                alpha = alpha + u_t(2*(j+0)-1+1)*u_t(2*(j+0)-0+1)
             enddo
             prod_uv = prod_uv - 2*alpha

             if ( n > VTOL ) then
                jj_2 = 1+LL*(y_inod-1)
                jj_3 = MIN(n, LL*y_inod)
             else
                jj_2 = 1
                jj_3 = n
             endif

!$OMP DO
             do jj_1=jj_2,jj_3,LX

                j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

                if(l_4-1==k_1+1)then
                   l_1 = k_1+1                           ! 0

                   j = l_1-k_1

                   u0 = u_t(2*(j+0)-1+1)
                   v0 = u_t(2*(j+0)-0+1)
                   do j_1=j_2,j_3
                      w0 = v_x(j_1)
                      ux0 = ux(j_1,l_1+0)
                      vx0 = vx(j_1,l_1+0)
                      w0 = w0
     &                   -ux0*u0
     &                   -vx0*v0
                      v_x(j_1) = w0
                   end do! j_1
                end if
                if(l_4-2==k_1+1)then
                   l_1 = k_1+1                           ! 1

                   j = l_1-k_1

                   u0 = u_t(2*(j+0)-1+1)
                   v0 = u_t(2*(j+0)-0+1)
                   u1 = u_t(2*(j+1)-1+1)
                   v1 = u_t(2*(j+1)-0+1)
                   do j_1=j_2,j_3
                      w0 = v_x(j_1)
                      ux0 = ux(j_1,l_1+0)
                      vx0 = vx(j_1,l_1+0)
                      w0 = w0
     &                   -ux0*u0
     &                   -vx0*v0
                      ux1 = ux(j_1,l_1+1)
                      vx1 = vx(j_1,l_1+1)
                      w0 = w0
     &                   -ux1*u1
     &                   -vx1*v1
                      v_x(j_1) = w0
                   end do! j_1
                end if
                do l_1=l_4,k_2,3                  ! 2

                   j = l_1-k_1

                   u0 = u_t(2*(j+0)-1+1)
                   v0 = u_t(2*(j+0)-0+1)
                   u1 = u_t(2*(j+1)-1+1)
                   v1 = u_t(2*(j+1)-0+1)
                   u2 = u_t(2*(j+2)-1+1)
                   v2 = u_t(2*(j+2)-0+1)
                   do j_1=j_2,j_3
                      w0 = v_x(j_1)
                      ux0 = ux(j_1,l_1+0)
                      vx0 = vx(j_1,l_1+0)
                      w0 = w0
     &                   -ux0*u0
     &                   -vx0*v0
                      ux1 = ux(j_1,l_1+1)
                      vx1 = vx(j_1,l_1+1)
                      w0 = w0
     &                   -ux1*u1
     &                   -vx1*v1
                      ux2 = ux(j_1,l_1+2)
                      vx2 = vx(j_1,l_1+2)
                      w0 = w0
     &                   -ux2*u2
     &                   -vx2*v2
                      v_x(j_1) = w0
                   end do! j_1
                end do! l_1

             end do! jj_1
!$OMP ENDDO

             if ( beta /= ZERO ) then

             alpha = prod_uv/(2*beta)

!$OMP DO
             do j_1=jj_2,jj_3
                v_x(j_1) = (v_x(j_1)-alpha*u_x(j_1))/beta
             end do! j_1
!$OMP ENDDO

             end if

!$OMP MASTER
             if ( n > VTOL ) then
                call allgather_dbl(v_x(jj_2), v_t, LL, 1, y_COMM_WORLD)
                j_3 = eigen_loop_end(L, x_nnod, x_inod)
                v_x(1:j_3) = v_t(1:j_3)
             end if
!$OMP END MASTER

          end if

!$OMP BARRIER

!$OMP MASTER

          v_n    = (v_n-alpha*u_n)/beta

          x_owner_nod = eigen_owner_node   (L, x_nnod,x_inod)
          x_pos       = eigen_translate_g2l(L, x_nnod,x_inod)
          if ( x_inod == x_owner_nod ) then
             v_x(x_pos) = v_n
          end if

          if ( kk == 0 ) then
          call datacast_dbl(v_y(1), v_x(1), u_t(1), v_t(1), x_pos, 2)
          end if

          call eigen_vector_zeropad_x( v_x(1), L )
          if ( kk == 0 ) then
          call eigen_vector_zeropad_y( v_y(1), L )
          end if

!$OMP END MASTER


      return
      end subroutine eigen_trd_vo

