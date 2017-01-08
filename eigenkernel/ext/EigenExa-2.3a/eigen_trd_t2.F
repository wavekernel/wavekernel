       subroutine  eigen_trd_au(
     &                a, nm,
     &                u_x, u_y, v_x, uy, vy, nv,
     &                u_t, v_t, d_t,
     &                i, i_base, m,
     &                w, e, beta, u_n, v_n)
!$     use omp_lib
       use MPI
       use eigen_devel
       use eigen_libs
       use eigen_house
       implicit NONE

       integer, intent(in)    ::  nm, nv, i, i_base, m
       real(8), intent(inout) ::  a(1:nm,*)
       real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) ::  v_x(1:nv)
       real(8), intent(in)    ::  uy(1:nv,*), vy(1:nv,*)
       real(8), intent(inout) ::  u_t(*), v_t(*)
       real(8), intent(in)    ::  d_t(*)
       real(8), intent(inout) ::  w(*), e(*), beta
       real(8), intent(out)   ::  u_n, v_n

       integer                ::  n, blk_0
       integer                ::  x_pos, y_pos, x_root, y_root
       integer                ::  n1, n2
       integer                ::  i_1, i_2, i_3
       integer                ::  j_1, j_2, j_3
       integer                ::  k_1, k_2, k_3, k_4
       integer                ::  L, j
       integer                ::  ierr

       real(8) :: anorm2, a_n, g_n
       real(8) :: prod_uv, prod_ua, d_n, w0, u0, v0, ux0, vx0

       real(8) :: d1, d2

       integer                ::  local_size, local_rank


          local_rank = 0
          local_size = 1
!$        local_rank = omp_get_thread_num()
!$        local_size = omp_get_num_threads()

          L = i-1

          i_2 = eigen_loop_start(1, y_nnod,y_inod)
          i_3 = eigen_loop_end  (L, y_nnod,y_inod)
          j_2 = eigen_loop_start(1, x_nnod,x_inod)
          j_3 = eigen_loop_end  (L, x_nnod,x_inod)

          y_root = eigen_owner_node (L, y_nnod,y_inod)
          y_pos  = eigen_translate_g2l(i, y_nnod,y_inod)

          x_root = eigen_owner_node (L, x_nnod,x_inod)
          x_pos  = eigen_translate_g2l(L, x_nnod,x_inod)

          k_1   = i-i_base
          k_2   = m
          n     = i_3
          blk_0 = k_1-k_2


       if ( local_size > 1 ) then

          n1 = offset1+nv*local_rank
          n2 = offset2+nv*local_rank

          do j_1=1,x_pos
             u0_z(j_1+n1) = ZERO
          enddo
          do j_1=1,y_pos
             v0_z(j_1+n2) = ZERO
          enddo

          call  eigen_trd_au_body1(
     &             a, nm,
     &             u_x, u_y, u0_z(1+n1), v0_z(1+n2),
     &             1, n, x_pos, y_pos, nv, blk_0
     &             ,local_rank, local_size
     &             )

!$OMP BARRIER

          call  eigen_trd_au_body2(
     &             v_x, u0_z(1+offset1),
     &             v_t, v0_z(1+offset2),
     &             nv, x_pos, y_pos
     &             ,local_rank, local_size
     &             )

       else

          do j_1=1,x_pos
             v_x(j_1) = ZERO
          enddo
          do j_1=1,y_pos
             v_t(j_1) = ZERO
          enddo

          call  eigen_trd_au_body1(
     &             a, nm,
     &             u_x, u_y, v_x, v_t,
     &             1, n, x_pos, y_pos, nv, blk_0
     &             ,local_rank, local_size
     &             )

       endif

!$OMP BARRIER

!$OMP MASTER

          anorm2  = ZERO ! norm(u0, 2)^2
          a_n     = ZERO ! u0_{L} = a_{L,i}
          prod_uv = ZERO ! u^tAu
          d_n     = ZERO ! a_{L,L}
          prod_ua = ZERO ! u0^t A_{:,L}


          ! u = u0 - s * e_L, s = - sign(u0_{L})*norm(u0,2)
          ! -> 
          ! u^t A u = (u0-se)^tA(u0-se)
          !         = u0^t A u0 - 2sa_{:,L}^tu0 + s^2a_{L,L}
          !         = u0^t A u0 - 2*g_n*prod_ua + g_n^2*d_n
          ! prod_uv = prod_uv + g_n * (g_n * d_n - 2 * prod_ua)

!DIR$ IVDEP
          do i_1=1,y_pos
             prod_uv = prod_uv + v_t(i_1)*u_y(i_1)
          end do

          if ( x_root == x_inod ) then
             a_n = u_x(x_pos)
          endif

          if ( y_root == y_inod ) then
             i_1  = eigen_translate_g2l(L, y_nnod,y_inod)
!DIR$ IVDEP
             do  j_1=1,x_pos
                w0  = u_x(j_1)
                anorm2  = anorm2  + w0**2
                prod_uv = prod_uv + w0 * v_x(j_1)
                prod_ua = prod_ua + w0 * a(j_1,i_1)
             end do
          if ( x_root == x_inod ) then
             d_n = d_t(i_1)
             prod_ua = prod_ua + d_n * u_x(x_pos)
          end if
          else
!DIR$ IVDEP
             do j_1=1,x_pos
                w0  = u_x(j_1)
                anorm2  = anorm2  + w0**2
                prod_uv = prod_uv + w0 * v_x(j_1)
             end do
          endif


       !  :: v_x += diag(A) * u_y   if diag(<L)
          call  eigen_trd_au_body4(
     &                u_x, u_y, v_x,
     &                u_t, v_t, d_t,
     &                L-1, x_pos, y_pos, nv, prod_uv
     &                )


          v_t(y_pos+1) = anorm2
          v_t(y_pos+2) = a_n
          v_t(y_pos+3) = prod_uv
          v_t(y_pos+4) = d_n
          v_t(y_pos+5) = prod_ua
          if ( TRD_nnod > 1 ) then
             call reduce_dbl(v_t, u_t, y_pos+5, 1, x_COMM_WORLD)
          end if
          anorm2  =  v_t(y_pos+1)
          a_n     =  v_t(y_pos+2)
          prod_uv =  v_t(y_pos+3)
          d_n     =  v_t(y_pos+4)
          prod_ua =  v_t(y_pos+5)

          if ( anorm2 /= ZERO ) then

          g_n     = -SIGN(SQRT(anorm2), a_n) ! s
          beta    =  anorm2 - a_n * g_n      ! norm(u)^2 / 2
          e (i)   =  g_n
          u_n     =  a_n - g_n

          else

          g_n     =  ZERO
          beta    =  ONE
          e (i)   =  ZERO
          u_n     =  ZERO

          end if

          if ( x_root == x_inod ) then
             u_x(x_pos) =  u_n
             i_1 = eigen_owner_node (i, y_nnod,y_inod)
             if ( i_1 == y_inod ) then
                w(x_pos) =  u_x(x_pos)
             end if
          end if
!
          if ( y_root == y_inod ) then
             i_1  = eigen_translate_g2l(L, y_nnod,y_inod)
             prod_uv = prod_uv + d_n * u_y(i_1)**2
             u_y(i_1) =  u_n

            ! v = Au = A(u0-se) = v0 - g_n * a_{:,L}
!DIR$ IVDEP
             do  j_1=1,x_pos
                v_x(j_1) = v_x(j_1) - g_n * a(j_1,i_1)
             end do

          if ( x_root == x_inod ) then
             j_1  = eigen_translate_g2l(L, x_nnod,x_inod)
                v_x(j_1) = v_x(j_1) + d_t(i_1) * u_y(i_1)
          end if
          end if

          if ( y_root == y_inod ) then
             i_1  = eigen_translate_g2l(L, y_nnod,y_inod)
             v_n  = v_t(i_1)
          else
             v_n  = ZERO
          end if

       !  :: v_x += v_t   if diag
          call  eigen_trd_au_body3(
     &                u_x, u_y, v_x,
     &                u_t, v_t, d_t,
     &                L-0, x_pos, y_pos, nv
     &                )


!$OMP END MASTER

!$OMP BARRIER

          k_1 = i - i_base
          k_2 = m
          k_4 = k_2-k_1
          k_3 = 1+k_4*2

          j = k_3+3+x_pos-k_1

       if ( k_2 > k_1 ) then
!$OMP DO
          do L=k_1+1,k_2
             u0 = ZERO
             v0 = ZERO
             do i_1=i_2,i_3
                w0 = u_y(i_1)
                u0 = u0 + uy(i_1, L) * w0
                v0 = v0 + vy(i_1, L) * w0
             enddo
             ux0 = ZERO
             vx0 = ZERO
             if ( y_root == y_inod ) then
                ux0 = uy(i_3, L)
                vx0 = vy(i_3, L)
             end if
             u_t(1+(L-k_1-1)*2+1) = v0
             u_t(1+(L-k_1-1)*2+2) = u0
             u_t(L+j+0*k_4)       = ux0
             u_t(L+j+1*k_4)       = vx0
          enddo
!$OMP END DO
       end if

!$OMP BARRIER

!$OMP MASTER

          u_t(1)     = prod_uv
          u_t(k_3+1) = prod_ua
          u_t(k_3+2) = d_n
          u_t(k_3+3) = v_n

          if ( TRD_nnod > 1 ) then
             u_t((k_3+3)+1:(k_3+3)+x_pos) = v_x(1:x_pos)
             call reduce_dbl(u_t, v_t, k_3+3+x_pos+2*(k_2-k_1),
     &                       2, y_COMM_WORLD)
             v_x(1:x_pos) = u_t((k_3+3)+1:(k_3+3)+x_pos)
          end if
          prod_uv = u_t(1)
          prod_ua = u_t(k_3+1)
          d_n     = u_t(k_3+2)
          v_n     = u_t(k_3+3)

          prod_uv = prod_uv + g_n * (g_n * d_n - 2 * prod_ua)
          u_t(1)  = prod_uv

          v_n     = v_n + d_n * u_n
          do L=k_1+1,k_2
             ux0 = u_t(L+j+0*k_4)
             vx0 = u_t(L+j+1*k_4)
             u0 = u_t(1+(L-k_1-1)*2+1)
             v0 = u_t(1+(L-k_1-1)*2+2)
             v_n = v_n - ux0*u0 - vx0*v0
          end do

!$OMP END MASTER


       return
       end subroutine  eigen_trd_au

       subroutine  eigen_trd_au_body1(
     &                a, nm,
     &                u_x, u_y, u_t, v_t,
     &                n1, n2, x_pos, y_pos, nv, blk_0
     &                ,local_rank, local_size
     &                )
!$     use omp_lib
       use eigen_devel
       use eigen_libs
       implicit NONE

       integer, intent(in)    :: nm, nv, n1, n2
       real(8), intent(inout) :: a(1:nm,*)
       real(8), intent(inout) :: u_x(1:nv), u_y(1:nv)
       real(8), intent(inout) :: u_t(1:nv), v_t(1:nv)
       integer, intent(in)    :: blk_0
       integer                :: x_pos, y_pos, x_root, y_root
       integer                :: local_size, local_rank

       integer                :: i_0
       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: k_1, k_2, k_3, k_4
       integer                :: l_1, l_2, l_3, l_4
       integer                :: i, j, k

       real(8)                :: v0, u0
       real(8)                :: a0_0
       real(8)                :: a0_1
       real(8)                :: a0_2
       real(8)                :: a0_3
       real(8)                :: a0_4
       real(8)                :: a0_5
       real(8)                :: a0_6
       real(8)                :: a0_7
       real(8)                :: w0_(0:7)
       real(8)                :: v_t0, u_y0
       real(8)                :: v_t1, u_y1
       real(8)                :: v_t2, u_y2
       real(8)                :: v_t3, u_y3
       real(8)                :: v_t4, u_y4
       real(8)                :: v_t5, u_y5
       real(8)                :: v_t6, u_y6
       real(8)                :: v_t7, u_y7

       integer                ::   LX, LY
       integer                ::   ii_1, ii_2, ii_3, ii_4, ii_5
       integer                ::   jj_1, jj_2, jj_3, jj_4, jj_5
       integer                ::   kk_1, kk_2, kk_3, kk_4, kk_5

#define                            DO_UNROLL   6
!#define                            DO_UNROLL   5
!#define                            DO_UNROLL   4
       integer, parameter     ::   ULEN   = 16*16*16
       integer                ::   VLEN
       integer                ::   DO_ITR = 1
!#       integer, parameter     ::   VLEN   = 16*9
#define                            K_LOOP   1


          VLEN = 16*15

          i_2 = n1
          i_3 = n2
!
! v:= Au
!
          if ( blk_0 == 0 ) then
             do i_1=i_2+local_rank,i_3,local_size
                j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
                j_3 = eigen_loop_end  (j, x_nnod,x_inod)
                j   = j+y_nnod*(DO_UNROLL*2-1)
                j_4 = eigen_loop_end  (j, x_nnod,x_inod)
                do j_1=j_3+1,MIN(j_4,nm)
                   a(j_1,i_1) = ZERO
                end do! j_1
             end do! i_1
!$OMP BARRIER
          end if
!
!
!          LX = 32*1000; LY = 32*1
          LX = 16*30*6; LY = 32*1
!
          l_2 = i_2; l_3 = i_3
          l_1 = l_2; l_4 = l_3

             k_2 = 1
             k   = eigen_translate_l2g(l_4, y_nnod,y_inod)
             k_3 = eigen_loop_end  (k-1, x_nnod,x_inod)
#if K_LOOP
             do k_1=k_2,k_3,LX; k_4 = MIN(k_3,k_1+LX-1)
#else
             k_1=k_2; k_4 = MIN(k_3,k_1+LX-1)
#endif

                j    = eigen_translate_l2g(k_1, x_nnod,x_inod)
                ii_2 = eigen_loop_start(j,y_nnod,y_inod)
                ii_2 = MAX(l_1,ii_2)
                ii_3 = l_4
#if K_LOOP
                if ( ii_2 > ii_3 ) CYCLE
#else
                if ( ii_2 > ii_3 ) goto 99999
#endif
                ii_4 = MOD(ii_3-ii_2+1,DO_UNROLL*2)+ii_2


!-------------------------------------------------------------
             do i_1=ii_2+local_rank,ii_4-1,local_size
                   j    = eigen_translate_l2g(i_1, y_nnod,y_inod)
                   j    = j+(1-1)*y_nnod
                   jj_2 = k_1
                   jj_3 = eigen_loop_end  (j-1, x_nnod,x_inod)
                   jj_3 = MIN(k_4,jj_3)
                   if ( jj_2 > jj_3 ) cycle
!=============================================================
                do kk_1=jj_2,jj_3,ULEN; kk_4=MIN(kk_1+ULEN-1,jj_3)

                   w0_(0) = v_t(i_1+0)
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=kk_1,kk_4
!                   do j_1=jj_2,jj_3

                      u0 = u_x(j_1+0)

                      a0_0 = a(j_1+0,i_1+0)
                      w0_(0) = w0_(0)
     &                            + (a0_0*u0)

                   end do! j_1
                   v_t(i_1+0) = w0_(0)

                   u_y0 = u_y(i_1+0)
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
                   do j_1=kk_1,kk_4
!                   do j_1=jj_2,jj_3

                      v0 = ZERO

                      a0_0 = a(j_1+0,i_1+0)
                      v0 = v0
     &                            + (a0_0*u_y0)
                      u_t(j_1+0) = v0 + u_t(j_1+0)

                   end do! j_1

                end do! kk_1
!=============================================================
             end do! i_1
!-------------------------------------------------------------


!-------------------------------------------------------------
             do i_0 = ii_4+local_rank*DO_UNROLL*DO_ITR,
     &                   ii_3,
     &                   local_size*DO_UNROLL*DO_ITR

                   j    = eigen_translate_l2g(i_0, y_nnod,y_inod)
                   j    = j+(DO_UNROLL*DO_ITR-1)*y_nnod
                   jj_2 = k_1
                   jj_3 = eigen_loop_end  (j-1, x_nnod,x_inod)
                   jj_3 = MIN(k_4,jj_3)
                   if ( jj_2 > jj_3 ) cycle
!=============================================================
             do kk_1=jj_2,jj_3,VLEN; kk_4=MIN(kk_1+VLEN-1,jj_3)

                do i_1 = i_0,i_0+DO_UNROLL*DO_ITR-1,DO_UNROLL

                   w0_(0) = v_t(i_1+0)
                   w0_(1) = v_t(i_1+1)
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
!DIR$ UNROLL(8)
!OCL LOOP_NOFUSION
!OCL UNROLL(8)
!OCL SIMD
                   do j_1=kk_1,kk_4

                      u0 = u_x(j_1+0)

                      a0_0 = a(j_1+0,i_1+0)
                      a0_1 = a(j_1+0,i_1+1)

                      w0_(0) = w0_(0)
     &                            + (a0_0*u0)
                      w0_(1) = w0_(1)
     &                            + (a0_1*u0)

                   end do

                   v_t(i_1+0) = w0_(0)
                   v_t(i_1+1) = w0_(1)

#                  if DO_UNROLL>=3
                   w0_(2) = v_t(i_1+2)
#                  endif
#                  if DO_UNROLL>=4
                   w0_(3) = v_t(i_1+3)
#                  endif
#                  if DO_UNROLL>=5
                   w0_(4) = v_t(i_1+4)
#                  endif
#                  if DO_UNROLL>=6
                   w0_(5) = v_t(i_1+5)
#                  endif
#                  if DO_UNROLL>=7
                   w0_(6) = v_t(i_1+6)
#                  endif
#                  if DO_UNROLL>=8
                   w0_(7) = v_t(i_1+7)
#                  endif

                   u_y0 = u_y(i_1+0)
                   u_y1 = u_y(i_1+1)
#               if DO_UNROLL>=3
                   u_y2 = u_y(i_1+2)
#               endif
#               if DO_UNROLL>=4
                   u_y3 = u_y(i_1+3)
#               endif
#               if DO_UNROLL>=5
                   u_y4 = u_y(i_1+4)
#               endif
#               if DO_UNROLL>=6
                   u_y5 = u_y(i_1+5)
#               endif
#               if DO_UNROLL>=7
                   u_y6 = u_y(i_1+6)
#               endif
#               if DO_UNROLL>=8
                   u_y7 = u_y(i_1+7)
#               endif

!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
!DIR$ UNROLL(4)
!OCL LOOP_NOFUSION
!OCL UNROLL(4)
!OCL SIMD
                   do j_1=kk_1,kk_4

                      u0 = u_x(j_1+0)
                      v0 = ZERO

                      a0_0 = a(j_1+0,i_1+0)
                      a0_1 = a(j_1+0,i_1+1)

                      v0 = v0
     &                            + (a0_0*u_y0)
     &                            + (a0_1*u_y1)

#                  if DO_UNROLL>=3
                      a0_2 = a(j_1+0,i_1+2)
#                  endif
#                  if DO_UNROLL>=4
                      a0_3 = a(j_1+0,i_1+3)
#                  endif
#                  if DO_UNROLL>=3
                      w0_(2) = w0_(2)
     &                            + (a0_2*u0)
#                  endif
#                  if DO_UNROLL>=4
                      w0_(3) = w0_(3)
     &                            + (a0_3*u0)
#                  endif
#                  if DO_UNROLL>=3
                      v0 = v0
     &                            + (a0_2*u_y2)
#                  endif
#                  if DO_UNROLL>=4
     &                            + (a0_3*u_y3)
#                  endif

#                  if DO_UNROLL>=5
                      a0_4 = a(j_1+0,i_1+4)
#                  endif
#                  if DO_UNROLL>=6
                      a0_5 = a(j_1+0,i_1+5)
#                  endif
#                  if DO_UNROLL>=5
                      w0_(4) = w0_(4)
     &                            + (a0_4*u0)
#                  endif
#                  if DO_UNROLL>=6
                      w0_(5) = w0_(5)
     &                            + (a0_5*u0)
#                  endif
#                  if DO_UNROLL>=5
                      v0 = v0
     &                            + (a0_4*u_y4)
#                  endif
#                  if DO_UNROLL>=6
     &                            + (a0_5*u_y5)
#                  endif

#                  if DO_UNROLL>=7
                      a0_6 = a(j_1+0,i_1+6)
#                  endif
#                  if DO_UNROLL>=8
                      a0_7 = a(j_1+0,i_1+7)
#                  endif
#                  if DO_UNROLL>=7
                      w0_(6) = w0_(6)
     &                            + (a0_6*u0)
#                  endif
#                  if DO_UNROLL>=8
                      w0_(7) = w0_(7)
     &                            + (a0_7*u0)
#                  endif
#                  if DO_UNROLL>=7
                      v0 = v0
     &                            + (a0_6*u_y6)
#                  endif
#                  if DO_UNROLL>=8
     &                            + (a0_7*u_y7)
#                  endif

                      u_t(j_1+0) = v0 + u_t(j_1+0)

                   end do! j_1

#                  if DO_UNROLL>=3
                   v_t(i_1+2) = w0_(2)
#                  endif
#                  if DO_UNROLL>=4
                   v_t(i_1+3) = w0_(3)
#                  endif
#                  if DO_UNROLL>=5
                   v_t(i_1+4) = w0_(4)
#                  endif
#                  if DO_UNROLL>=6
                   v_t(i_1+5) = w0_(5)
#                  endif
#                  if DO_UNROLL>=7
                   v_t(i_1+6) = w0_(6)
#                  endif
#                  if DO_UNROLL>=8
                   v_t(i_1+7) = w0_(7)
#                  endif

                end do! i_1

             end do! kk_1
!=============================================================
             end do! i_1
!-------------------------------------------------------------
99999           continue
#if K_LOOP
             end do! k_1
#endif


       return
       end subroutine  eigen_trd_au_body1

       subroutine  eigen_trd_au_body2(
     &                u_t, u_z, v_t, v_z, nv,
     &                x_pos, y_pos
     &                ,local_rank, local_size
     &                )
       use eigen_devel
       implicit NONE

       integer, intent(in)    :: nv, x_pos, y_pos
       real(8), intent(out)   :: u_t(1:nv),   v_t(1:nv)
       real(8), intent(in)    :: u_z(1:nv,*), v_z(1:nv,*)

       integer                :: local_rank, local_size

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: jj_1, jj_2, jj_3, jj_4
       integer                :: i, j, k

       integer, parameter     :: LX = 256

       integer                :: ll_size, ll_rank

       if ( MOD(local_size,2) == 0 ) then
          ll_size = local_size/2
       else
          ll_size = local_size
       end if
          ll_rank = MOD( local_rank, ll_size )


       if ( ll_size == local_size .OR.
     &     (local_rank < ll_size) ) then

          jj_1 = x_pos
          jj_2 = MAX(512,(jj_1-1)/ll_size+1)
          jj_3 =    (jj_2*(ll_rank+0)     )+1
          jj_4 = MIN(jj_2*(ll_rank+1),jj_1)

          do jj_1=jj_3,jj_4,LX
             j_3=jj_1; j_4=MIN(jj_1+LX-1,jj_4)

         if ( MOD(local_size,2) == 1 ) then
            do j_1=j_3,j_4
               u_t(j_1) = u_z(j_1,1)
            end do
         else
            do j_1=j_3,j_4
               u_t(j_1) = u_z(j_1,1)+u_z(j_1,2)
            end do
         end if
         if ( local_size > 2 ) then
         do j=3-MOD(local_size,2),local_size,2
            do j_1=j_3,j_4
               u_t(j_1) = u_t(j_1)+u_z(j_1,j+0)+u_z(j_1,j+1)
            end do
         end do
         end if

          end do

       endif

       if ( ll_size == local_size .OR.
     &     (local_rank >= ll_size .AND. local_rank < 2*ll_size) ) then

          jj_1 = y_pos
          jj_2 = MAX(512,(jj_1-1)/ll_size+1)
          jj_3 =    (jj_2*(ll_rank+0)     )+1
          jj_4 = MIN(jj_2*(ll_rank+1),jj_1)

          do jj_1=jj_3,jj_4,LX
             j_3=jj_1; j_4=MIN(jj_1+LX-1,jj_4)
           if ( local_size == 4 ) then
                do j_1=j_3,j_4
                  v_t(j_1) = v_z(j_1,1)+v_z(j_1,2)+v_z(j_1,3)+v_z(j_1,4)
                end do
           else
             if ( MOD(local_size,2) == 1 ) then
                do j_1=j_3,j_4
                   v_t(j_1) = v_z(j_1,1)
                end do
             else
                do j_1=j_3,j_4
                   v_t(j_1) = ZERO
                end do
             end if
             do j=MOD(local_size,2)+1,local_size,2
                do j_1=j_3,j_4
                   v_t(j_1) = v_t(j_1)+v_z(j_1,j+0)+v_z(j_1,j+1)
                end do
             end do
           end if
          end do

       end if


       return
       end subroutine  eigen_trd_au_body2

       subroutine  eigen_trd_au_body3(
     &                u_x, u_y, v_x,
     &                u_t,v_t,d_t,
     &                n, x_pos, y_pos, nv
     &                )
       use eigen_devel
       use eigen_libs
       implicit NONE

       integer, intent(in)    :: nv, n
       real(8), intent(in)    :: u_x(1:nv),u_y(1:nv)
       real(8), intent(inout) :: v_x(1:nv)
       real(8), intent(in)    :: u_t(*),v_t(*)
       real(8), intent(in)    :: d_t(1:nv)
       integer                :: x_pos, y_pos, x_root, y_root

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: k_1, k_2, k_3, k_4
       integer                :: i, j, k
       integer                :: nx, nm1, nm2


       if ( diag_0 > 0 ) then
          i_2 = diag_0
          i_3 = eigen_loop_end  (n, y_nnod, y_inod)
       if ( i_2 <= i_3 ) then
          j_2 = diag_1
          j_3 = eigen_loop_end  (n, x_nnod, x_inod)
       if ( j_2 <= j_3 ) then
          nm1 = y_nnod/n_common
          nm2 = x_nnod/n_common
          nx = (i_3-i_2)/nm2+1

!         do  k_1 = k_2, k_3
!            i_1 = i_2 + k_1 * nm2
!            j_1 = j_2 + k_1 * nm1
!            v_x(j_1) = v_x(j_1) + d_t(i_1) * u_y(i_1)
!         end do ! k_1

          if ( nm2 == 1 ) then
          if ( nm1 == 1 ) then
                call eigen_trd_au_body3_sub11(
     &               v_x(diag_1), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               nx, nm1, nm2)
          else
                call eigen_trd_au_body3_subX1(
     &               v_x(diag_1), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               nx, nm1, nm2)
          end if
          else
                call eigen_trd_au_body3_subXX(
     &               v_x(diag_1), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               nx, nm1, nm2)

          end if

       end if
       end if
       end if


       return
       end subroutine  eigen_trd_au_body3

       subroutine eigen_trd_au_body3_subXX(v_x,v_t,d_t,u_y, n,nm1,nm2)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(nm1,*)
       real(8), intent(in)    :: v_t(nm2,*)
       real(8), intent(in)    :: d_t(nm2,*)
       real(8), intent(in)    :: u_y(nm2,*)

       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             v_x(1,i) = v_x(1,i)+v_t(1,i)
          end do! i


       return
       end subroutine  eigen_trd_au_body3_subXX

       subroutine eigen_trd_au_body3_subX1(v_x,v_t,d_t,u_y, n,nm1,nm2)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(nm1,*)
       real(8), intent(in)    :: v_t(*)
       real(8), intent(in)    :: d_t(*)
       real(8), intent(in)    :: u_y(*)

       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             v_x(1,i) = v_x(1,i)+v_t(i)
          end do! i


       return
       end subroutine  eigen_trd_au_body3_subX1

       subroutine eigen_trd_au_body3_sub11(v_x,v_t,d_t,u_y, n,nm1,nm2)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(1:n)
       real(8), intent(in)    :: v_t(1:n)
       real(8), intent(in)    :: d_t(1:n)
       real(8), intent(in)    :: u_y(1:n)

       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             v_x(i) = v_x(i)+v_t(i)
          end do! i


       return
       end subroutine  eigen_trd_au_body3_sub11

       subroutine  eigen_trd_au_body4(
     &                u_x, u_y, v_x,
     &                u_t,v_t,d_t,
     &                n, x_pos, y_pos, nv, prod_uv
     &                )
       use eigen_devel
       use eigen_libs
       implicit NONE

       integer, intent(in)    :: nv, n
       real(8), intent(in)    :: u_x(1:nv),u_y(1:nv)
       real(8), intent(inout) :: v_x(1:nv)
       real(8), intent(in)    :: u_t(*),v_t(*)
       real(8), intent(in)    :: d_t(1:nv)
       real(8), intent(inout) :: prod_uv
       integer                :: x_pos, y_pos, x_root, y_root

       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       integer                :: k_1, k_2, k_3, k_4
       integer                :: i, j, k
       integer                :: nx, nm1, nm2

       real(8) :: t


       if ( diag_0 > 0 ) then
          i_2 = diag_0
          i_3 = eigen_loop_end  (n, y_nnod, y_inod)
       if ( i_2 <= i_3 ) then
          j_2 = diag_1
          j_3 = eigen_loop_end  (n, x_nnod, x_inod)
       if ( j_2 <= j_3 ) then
          nm1 = y_nnod/n_common
          nm2 = x_nnod/n_common
          nx = (i_3-i_2)/nm2+1

!         do  k_1 = k_2, k_3
!            i_1 = i_2 + k_1 * nm2
!            j_1 = j_2 + k_1 * nm1
!            v_x(j_1) = v_x(j_1) + d_t(i_1) * u_y(i_1)
!         end do ! k_1

          if ( nm2 == 1 ) then
          if ( nm1 == 1 ) then
                call eigen_trd_au_body4_sub11(
     &               v_x(diag_1), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               nx, nm1, nm2, prod_uv)
          else
                call eigen_trd_au_body4_subX1(
     &               v_x(diag_1), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               nx, nm1, nm2, prod_uv)
          end if
          else
                call eigen_trd_au_body4_subXX(
     &               v_x(diag_1), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &               nx, nm1, nm2, prod_uv)

          end if

       end if
       end if
       end if


       return
       end subroutine  eigen_trd_au_body4

       subroutine eigen_trd_au_body4_subXX(v_x,v_t,d_t,u_y, n,nm1,nm2,
     &                                     prod_uv)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(nm1,*)
       real(8), intent(in)    :: v_t(nm2,*)
       real(8), intent(in)    :: d_t(nm2,*)
       real(8), intent(in)    :: u_y(nm2,*)
       real(8), intent(inout) :: prod_uv

       real(8)                :: s
       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             s = d_t(1,i) * u_y(1,i)
             v_x(1,i) = v_x(1,i) + s
             prod_uv  = prod_uv  + s * u_y(1,i)
          end do! i


       return
       end subroutine  eigen_trd_au_body4_subXX

       subroutine eigen_trd_au_body4_subX1(v_x,v_t,d_t,u_y, n,nm1,nm2,
     &                                     prod_uv)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(nm1,*)
       real(8), intent(in)    :: v_t(*)
       real(8), intent(in)    :: d_t(*)
       real(8), intent(in)    :: u_y(*)
       real(8), intent(inout) :: prod_uv

       real(8)                :: s
       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             s = d_t(i) * u_y(i)
             v_x(1,i) = v_x(1,i) + s
             prod_uv  = prod_uv  + s * u_y(i)
          end do! i


       return
       end subroutine  eigen_trd_au_body4_subX1

       subroutine eigen_trd_au_body4_sub11(v_x,v_t,d_t,u_y, n,nm1,nm2,
     &                                     prod_uv)
       implicit NONE

       integer, intent(in)    :: n, nm1, nm2
       real(8), intent(inout) :: v_x(1:n)
       real(8), intent(in)    :: v_t(1:n)
       real(8), intent(in)    :: d_t(1:n)
       real(8), intent(in)    :: u_y(1:n)
       real(8), intent(inout) :: prod_uv

       real(8)                :: s
       integer                :: i


!DIR$ VECTOR ALWAYS
          do i=1,n
             s = d_t(i) * u_y(i)
             v_x(i)  = v_x(i) + s
             prod_uv = prod_uv + s * u_y(i)
          end do! i


       return
       end subroutine  eigen_trd_au_body4_sub11


