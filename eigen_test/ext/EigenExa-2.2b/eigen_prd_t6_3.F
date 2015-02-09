      subroutine eigen_prd_vo(
     &              ux, vx, vy, nv,
     &              u_t, v_t, c,
     &              i, i_base, m0)
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      use eigen_house, only : tsave
      implicit NONE
!
      integer, intent(in)    ::  nv, i, i_base, m0
      real(8), intent(inout) ::  ux(1:nv, *)
      real(8), intent(inout) ::  vx(1:nv, *), vy(1:nv, *)
      real(8), intent(inout) ::  u_t(4,*)
      real(8), intent(inout) ::  v_t(*)
!
      real(8), intent(inout) ::  c(1:2,1:2)
!
      real(8)                ::  g_g(1:2,1:2)
!
      integer                ::  j, k, L, n, LL
      integer                ::  k_1, k_2
!
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  l_1, l_2, l_3, l_4
      integer                ::  jj_1, jj_2, jj_3, jj_4
      integer                ::  kk_1, kk_2, kk_3, kk_4

      integer, parameter     ::  LX   = 1024
      integer, parameter     ::  LX2  = 32

#if __IBM_REGISTER_VARS
      integer, parameter     ::  VTOL = 1024000
#else
      integer, parameter     ::  VTOL = 1024
#endif
!
      include 'CSTAB.h'
!
      real(8)                ::  w0, w1, w2, w3
      real(8)                ::  u0_0, v0_0
      real(8)                ::  u1_0, v1_0
      real(8)                ::  u0_1, v0_1
      real(8)                ::  u1_1, v1_1
      real(8)                ::  uu_x
      real(8)                ::  ux0, vx0
      real(8)                ::  ux1, vx1

      real(8)                ::  s11, s21, s12, s22
      real(8)                ::  t11, t21, t12, t22
      real(8)                ::  c11, c21, c12, c22
      real(8)                ::  u0, u1, v0, v1, u12


      k_1 = i - i_base
      k_2 = m0

      L = i - MBAND
      n   = eigen_translate_g2l(L, x_nnod, x_inod)


!$OMP MASTER
! FOR attention to unexpected overflow or NAN
      j_3 = eigen_loop_end     (L, x_nnod, x_inod)
      if ( j_3 < n ) then
         do i_1=1,MBAND
            vx(j_3+1:n,k_1-i_1+1) = ZERO ! in case
         end do
      end if
!$OMP END MASTER


      l_4 = MOD(k_2-k_1,2)+k_1+1

      if ( k_2 <= k_1 ) then

         jj_2 = 1
         jj_3 = eigen_loop_end(L, x_nnod,x_inod)

      else

         LL = (n-1)/y_nnod+1
         LL = ((LL-1)/2+1)*2

         if ( n > VTOL ) then
            jj_2 = 1+LL*(y_inod-1)
            jj_3 = MIN(n, LL*y_inod)
         else
            jj_2 = 1
            jj_3 = n
         endif

      end if


      if ( k_2 <= k_1 ) then

!$OMP MASTER
         u_t(1:4, 1:2)=ZERO

         do jj_1=jj_2,jj_3,LX
            u0_0 = ZERO 
            u0_1 = ZERO 
            u1_0 = ZERO 
            u1_1 = ZERO 
            uu_x = ZERO
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1,jj_3)
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
            do j_1=j_2,j_3
               u0 = ux(j_1, k_1-0)
               u1 = ux(j_1, k_1-1)
               v0 = vx(j_1, k_1-0)
               v1 = vx(j_1, k_1-1)
               u0_0 = u0_0 + u1 * v1
               u0_1 = u0_1 + u1 * v0
               u1_0 = u1_0 + u0 * v1
               u1_1 = u1_1 + u0 * v0
               uu_x = uu_x + u0 * u1
            end do
         j = 1
            u_t(1, j+0) = u_t(1, j+0) + u0_0
            u_t(2, j+0) = u_t(2, j+0) + u0_1
            u_t(3, j+0) = u_t(3, j+0) + u1_0
            u_t(4, j+0) = u_t(4, j+0) + u1_1
            u_t(1, j+1) = u_t(1, j+1) + uu_x
         end do

         call reduce_dbl(u_t, v_t, 5, 4, x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER

         g_g(1, 1) =  u_t(1, 1)
         g_g(1, 2) = (u_t(2, 1) + u_t(3, 1))/2
         g_g(2, 2) =  u_t(4, 1)
         g_g(2, 1) =  u_t(1, 2)

      else

!
! v=v-(UV+VU)u
!

         l_2 = k_2-k_1+2

!$OMP MASTER
         u_t(1:4, 1:l_2) = ZERO

         do jj_1=jj_2,jj_3,LX
            u0_0 = ZERO 
            u0_1 = ZERO 
            u1_0 = ZERO 
            u1_1 = ZERO 
            uu_x = ZERO
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1,jj_3)
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
            do j_1=j_2,j_3
               u0 = ux(j_1, k_1-0)
               u1 = ux(j_1, k_1-1)
               v0 = vx(j_1, k_1-0)
               v1 = vx(j_1, k_1-1)
               u0_0 = u0_0 + u1 * v1
               u0_1 = u0_1 + u1 * v0
               u1_0 = u1_0 + u0 * v1
               u1_1 = u1_1 + u0 * v0
               uu_x = uu_x + u0 * u1
            end do
         j = 1
            u_t(1, j+0) = u_t(1, j+0) + u0_0
            u_t(2, j+0) = u_t(2, j+0) + u0_1
            u_t(3, j+0) = u_t(3, j+0) + u1_0
            u_t(4, j+0) = u_t(4, j+0) + u1_1
            u_t(1, j+1) = u_t(1, j+1) + uu_x
         end do
!$OMP END MASTER
!$OMP BARRIER

         do jj_1=jj_2,jj_3,LX
            j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

!$OMP DO
            do l_1=k_1+1,l_4-1

               u0_0 = ZERO
               v0_0 = ZERO
               u1_0 = ZERO
               v1_0 = ZERO
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
               do j_1=j_2,j_3
                  w0 = ux(j_1, k_1-0)
                  u0_0 = u0_0 + vx(j_1, l_1+0) * w0
                  v0_0 = v0_0 + ux(j_1, l_1+0) * w0
                  w1 = ux(j_1, k_1-1)
                  u1_0 = u1_0 + vx(j_1, l_1+0) * w1
                  v1_0 = v1_0 + ux(j_1, l_1+0) * w1
               end do! j_1
            j = l_1-k_1 +2
               u_t(1, j+0) = u_t(1, j+0) + u0_0
               u_t(2, j+0) = u_t(2, j+0) + v0_0
               u_t(3, j+0) = u_t(3, j+0) + u1_0
               u_t(4, j+0) = u_t(4, j+0) + v1_0

            end do! l_1
!$OMP END DO NOWAIT
!$OMP DO
            do l_1=l_4,k_2,2                      ! 1

               u0_0 = ZERO
               v0_0 = ZERO
               u1_0 = ZERO
               v1_0 = ZERO
               u0_1 = ZERO
               v0_1 = ZERO
               u1_1 = ZERO
               v1_1 = ZERO
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
               do j_1=j_2,j_3
                  w0 = ux(j_1, k_1-0)
                  u0_0 = u0_0 + vx(j_1, l_1+0) * w0
                  v0_0 = v0_0 + ux(j_1, l_1+0) * w0
                  u0_1 = u0_1 + vx(j_1, l_1+1) * w0
                  v0_1 = v0_1 + ux(j_1, l_1+1) * w0
                  w1 = ux(j_1, k_1-1)
                  u1_0 = u1_0 + vx(j_1, l_1+0) * w1
                  v1_0 = v1_0 + ux(j_1, l_1+0) * w1
                  u1_1 = u1_1 + vx(j_1, l_1+1) * w1
                  v1_1 = v1_1 + ux(j_1, l_1+1) * w1
               end do! j_1
            j = l_1-k_1 +2
               u_t(1, j+0) = u_t(1, j+0) + u0_0
               u_t(2, j+0) = u_t(2, j+0) + v0_0
               u_t(3, j+0) = u_t(3, j+0) + u1_0
               u_t(4, j+0) = u_t(4, j+0) + v1_0
               u_t(1, j+1) = u_t(1, j+1) + u0_1
               u_t(2, j+1) = u_t(2, j+1) + v0_1
               u_t(3, j+1) = u_t(3, j+1) + u1_1
               u_t(4, j+1) = u_t(4, j+1) + v1_1

            end do! l_1
!$OMP END DO

         end do! jj_1

!$OMP BARRIER

!$OMP MASTER
         l_2 = (k_2-k_1+2)*2*MBAND
         if ( n > VTOL ) then
            call reduce_dbl(u_t, v_t, l_2, 5, y_COMM_WORLD)
         end if
         call reduce_dbl(u_t, v_t, l_2, 4, x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER


         g_g(1, 1) =  u_t(1, 1)
         g_g(1, 2) = (u_t(2, 1) + u_t(3, 1))/2
         g_g(2, 2) =  u_t(4, 1)
         g_g(2, 1) =  u_t(1, 2)

         c11 = ZERO
         c21 = ZERO
         c12 = ZERO
         c22 = ZERO

!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
         do l_1=k_1+1,k_2
            j = l_1-k_1 +2

            u0_0 = u_t(1, j+0)
            v0_0 = u_t(2, j+0)
            u1_0 = u_t(3, j+0)
            v1_0 = u_t(4, j+0)

            c11 = c11 + v1_0 * u1_0
            c21 = c21 + v1_0 * u0_0
            c12 = c12 + v0_0 * u1_0
            c22 = c22 + v0_0 * u0_0
         end do

         g_g(1, 1) = g_g(1, 1) - 2 * c11
         g_g(1, 2) = g_g(1, 2) -    (c12 + c21)
         g_g(2, 2) = g_g(2, 2) - 2 * c22

!$OMP DO
         do jj_1=jj_2,jj_3,LX2
            j_2 = jj_1; j_3 = MIN(jj_1+LX2-1, jj_3)

            do l_1=k_1+1,l_4-1

            j = l_1-k_1 +2
               u0_0 = u_t(1, j+0)
               v0_0 = u_t(2, j+0)
               u1_0 = u_t(3, j+0)
               v1_0 = u_t(4, j+0)

!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
               do j_1=j_2,j_3
                  ux0 = ux(j_1, l_1+0)
                  vx0 = vx(j_1, l_1+0)
                  w0 =
     &               - ux0 * u0_0
     &               - vx0 * v0_0
                  w1 =
     &               - ux0 * u1_0
     &               - vx0 * v1_0
                  vx(j_1, k_1-0) = vx(j_1, k_1-0) + w0
                  vx(j_1, k_1-1) = vx(j_1, k_1-1) + w1
               end do! j_1
            end do! l_1
            do l_1=l_4,k_2,2

            j = l_1-k_1 +2
               u0_0 = u_t(1, j+0)
               v0_0 = u_t(2, j+0)
               u1_0 = u_t(3, j+0)
               v1_0 = u_t(4, j+0)
               u0_1 = u_t(1, j+1)
               v0_1 = u_t(2, j+1)
               u1_1 = u_t(3, j+1)
               v1_1 = u_t(4, j+1)

!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
               do j_1=j_2,j_3
                  ux0 = ux(j_1, l_1+0)
                  vx0 = vx(j_1, l_1+0)
                  ux1 = ux(j_1, l_1+1)
                  vx1 = vx(j_1, l_1+1)
                  w0 =
     &               - ux0 * u0_0
     &               - vx0 * v0_0
     &               - ux1 * u0_1
     &               - vx1 * v0_1
                  w1 =
     &               - ux0 * u1_0
     &               - vx0 * v1_0
     &               - ux1 * u1_1
     &               - vx1 * v1_1
                  vx(j_1, k_1-0) = vx(j_1, k_1-0) + w0
                  vx(j_1, k_1-1) = vx(j_1, k_1-1) + w1
               end do! j_1
            end do! l_1

         end do! jj_1
!$OMP END DO

      end if

!$OMP BARRIER

!
! C:=[b0 -b0*b1(u0,u1)]
!    [    b1          ]
      c11 =  c(1, 1)
      c21 =  ZERO
      c12 = -c(2, 2) * c(1, 1) * g_g(2, 1) ! (u1, u0)
      c22 =  c(2, 2)

!
      s11 =  g_g(1, 1) ! (u0, u0)_A
      s21 =  g_g(1, 2) ! (u1, u0)_A
      s12 =  g_g(1, 2) ! (u1, u0)_A
      s22 =  g_g(2, 2) ! (u1, u1)_A

!
! TX:=SX*C^T=[s11 s12][c11    ]
!            [s21 s22][c12 c22]
      t11 =  s11 * c11 + s12 * c12
      t21 =  s21 * c11 + s22 * c12
      t12 =              s12 * c22
      t22 =              s22 * c22

!
! S:=C*TX   =[c11 c12][t11 t12]
!            [    c22][t21 t22]
      s11 =  c11 * t11 + c12 * t21
      s21 =              c22 * t21
      s12 =  c11 * t12 + c12 * t22
      s22 =              c22 * t22

!
! S:=diag(S)/2+upper(S)
      s11 = (s11      )/2
      s12 = (s21 + s12)/2
      s21 =  ZERO
      s22 = (s22      )/2

!
! V:=VC^T=[v0 v1][c11    ]
!                [c12 c22]
!$OMP DO
      do jj_1=jj_2,jj_3,LX2
         j_2 = jj_1; j_3 = MIN(jj_1+LX2-1, jj_3)
!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
         do j_1=j_2,j_3
            v0 = vx(j_1, k_1-1)
            v1 = vx(j_1, k_1-0)
            v0             =      v0 * c11 + v1 * c12
            v1             =                 v1 * c22
!
! V:=V-US=[v0 v1]-[u0 u1][s11 s12]
!                        [    s22]
            u0 = ux(j_1, k_1-1)
            u1 = ux(j_1, k_1-0)
            vx(j_1, k_1-1) = v0 - u0 * s11 - u1 * s12
            vx(j_1, k_1-0) = v1            - u1 * s22
         end do! j_1
      end do! jj_1
!$OMP END DO


!$OMP MASTER

      if ( k_2 > k_1 ) then
      if ( n > VTOL ) then
         j_3 = eigen_loop_end(L, x_nnod, x_inod)
         call allgather_dbl(vx(jj_2, k_1-0), v_t, LL, y_COMM_WORLD)
         vx(1:j_3, k_1-0) = v_t(1:j_3)
         call allgather_dbl(vx(jj_2, k_1-1), v_t, LL, y_COMM_WORLD)
         vx(1:j_3, k_1-1) = v_t(1:j_3)
      end if
      end if

      call datacast_dbl2(vy(1, k_1-0), vy(1, k_1-1),
     &                   vx(1, k_1-0), vx(1, k_1-1), u_t, v_t, n)


! FOR attention to unexpected overflow or NAN
      j_3 = eigen_loop_end     (L, x_nnod, x_inod)
      n   = eigen_translate_g2l(L, x_nnod, x_inod)
      if ( j_3 < n ) then
         vx(j_3+1:n, k_1-0) = ZERO ! in case
         vx(j_3+1:n, k_1-1) = ZERO ! in case
      end if

      j_3 = eigen_loop_end     (L, y_nnod, y_inod)
      n   = eigen_translate_g2l(L, y_nnod, y_inod)
      if ( j_3 < n ) then
         vy(j_3+1:n, k_1-0) = ZERO ! in case
         vy(j_3+1:n, k_1-1) = ZERO ! in case
      end if

!$OMP END MASTER


      return
      contains

      subroutine OMP_DO_SUM( x, y )
      implicit NONE
      real(8), intent(IN)    :: x
      real(8), intent(OUT)   :: y

!$    if ( omp_get_num_threads() == 1 ) then
         y = x
!$    else
!$OMP BARRIER
!$       tsave(omp_get_thread_num(), 1) = x
!$OMP BARRIER
!$OMP MASTER
!$       y = SUM( tsave(0:omp_get_num_threads()-1, 1) )
!$OMP END MASTER
!$OMP BARRIER
!$    end if

      return
      end subroutine OMP_DO_SUM

      end subroutine  eigen_prd_vo

