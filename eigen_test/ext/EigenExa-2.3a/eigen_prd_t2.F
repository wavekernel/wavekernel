      subroutine  eigen_prd_au(
     &               a, nm,
     &               u_x, u_y, v_x, nv,
     &               u_t, v_t, d_t,
     &               i, i_base, m)

!$    use OMP_LIB
      use MPI
      use eigen_devel
      use eigen_libs
      use eigen_house
      implicit NONE

      integer, intent(in)    :: nm, nv, i, i_base, m
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(inout) :: u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(inout) :: v_x(1:nv,*)
      real(8), intent(inout) :: u_t(*), v_t(1:nv,*)
      real(8), intent(in)    :: d_t(*)

      integer                :: n, blk_0
      integer                :: x_pos, y_pos
      integer                :: n1, n2, n3, n4
      integer                :: i_1, i_2, i_3
      integer                :: j_1, j_2, j_3
      integer                :: k_1, k_2, k_3
      integer                :: L

      integer                :: local_size, local_rank
      integer                :: ll_size, ll_rank

      real(8) :: d1, d2


      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      L = i-2

      i_2 = eigen_loop_start(1, y_nnod,y_inod)
      i_3 = eigen_loop_end  (L, y_nnod,y_inod)
      j_2 = eigen_loop_start(1, x_nnod,x_inod)
      j_3 = eigen_loop_end  (L, x_nnod,x_inod)

      y_pos  = eigen_translate_g2l(i, y_nnod,y_inod)
      x_pos  = eigen_translate_g2l(L, x_nnod,x_inod)

      k_1   = i-i_base
      k_2   = m
      n     = i_3
      blk_0 = k_1-k_2


      n1 = offset1+nv*local_rank
      n2 = offset2+nv*local_rank
      n3 = offset3+nv*local_rank
      n4 = offset4+nv*local_rank

      do j_1=1,x_pos
         u0_z(j_1+n1) = ZERO
         u1_z(j_1+n3) = ZERO
      end do! j_1
      do j_1=1,y_pos
         v0_z(j_1+n2) = ZERO
         v1_z(j_1+n4) = ZERO
      end do! j_1

         call  eigen_prd_au_body1(
     &            a, nm,
     &            u_x(1,1), u_y(1,1),
     &            u_x(1,2), u_y(1,2),
     &            u0_z(1+n1), v0_z(1+n2),
     &            u1_z(1+n3), v1_z(1+n4),
     &            1, n, x_pos, y_pos, nv, blk_0
     &            ,local_rank, local_size
     &            )

!$OMP BARRIER

      call  eigen_prd_au_body2(
     &         v_x(1,1), u0_z(1+offset1), u1_z(1+offset3),
     &         v_t(1,1), v0_z(1+offset2), v1_z(1+offset4),
     &         nv, x_pos, y_pos
     &         ,local_rank, local_size
     &         )

!$OMP BARRIER

!$OMP MASTER

      if ( TRD_nnod > 1 ) then
         u_t(1:y_pos) = v_t(1:y_pos,1)
         u_t(1+y_pos:2*y_pos) = v_t(1:y_pos,2)
         call reduce_dbl(u_t, v_t, 2*y_pos, 1, x_COMM_WORLD)
         v_t(1:y_pos,1) = u_t(1:y_pos)
         v_t(1:y_pos,2) = u_t(1+y_pos:2*y_pos)
      end if

!$OMP END MASTER

!$OMP BARRIER

!$    if ( local_size == 1 .OR. local_rank == 1 ) then
      call  eigen_prd_au_body3(
     &         u_x, u_y, v_x,
     &         u_t, v_t, d_t,
     &         1, n, x_pos, y_pos, nv
     &         )
!$    end if

!$OMP BARRIER

!$OMP MASTER

      if ( TRD_nnod > 1 ) then
         v_t(1:x_pos,1) = v_x(1:x_pos,1)
         v_t(1+x_pos:2*x_pos,1) = v_x(1:x_pos,2)
         call reduce_dbl(v_t, u_t, 2*x_pos, 2, y_COMM_WORLD)
         v_x(1:x_pos,1) = v_t(1:x_pos,1)
         v_x(1:x_pos,2) = v_t(1+x_pos:2*x_pos,1)
      end if

!$OMP END MASTER

!$OMP BARRIER

      return
      end subroutine  eigen_prd_au

      subroutine  eigen_prd_au_body1(
     &               a, nm,
     &               u0_x, u0_y, u1_x, u1_y,
     &               u0_t, v0_t, u1_t, v1_t,
     &               n1, n2, x_pos, y_pos, nv, blk_0
     &               ,local_rank, local_size
     &               )
      use eigen_devel
      use eigen_libs
      implicit NONE

      integer, intent(in)    :: nm, nv, n1, n2
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(inout) :: u0_x(1:nv), u0_y(1:nv)
      real(8), intent(inout) :: u1_x(1:nv), u1_y(1:nv)
      real(8), intent(inout) :: u0_t(1:nv), v0_t(1:nv)
      real(8), intent(inout) :: u1_t(1:nv), v1_t(1:nv)
      integer, intent(in)    :: blk_0
      integer                :: x_pos, y_pos
      integer                :: local_size, local_rank

      integer                :: i_0, i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: k_1, k_2, k_3, k_4
      integer                :: l_1, l_2, l_3, l_4
      integer                :: i, j, k

      real(8)                :: v0, u0
      real(8)                :: v1, u1
      real(8)                :: w0_(0:8), w1_(0:8)
!//#define w0_(x)  w0_ x
!//#define w1_(x)  w1_ x
      real(8)                :: w0_0, w1_0
      real(8)                :: w0_1, w1_1
      real(8)                :: w0_2, w1_2
      real(8)                :: w0_3, w1_3
      real(8)                :: w0_4, w1_4
      real(8)                :: w0_5, w1_5
      real(8)                :: w0_6, w1_6
      real(8)                :: w0_7, w1_7
      real(8)                :: a0_0
      real(8)                :: a0_1
      real(8)                :: a0_2
      real(8)                :: a0_3
      real(8)                :: a0_4
      real(8)                :: a0_5
      real(8)                :: a0_6
      real(8)                :: a0_7
      real(8)                :: v0_t0, u0_y0
      real(8)                :: v1_t0, u1_y0
      real(8)                :: v0_t1, u0_y1
      real(8)                :: v1_t1, u1_y1
      real(8)                :: v0_t2, u0_y2
      real(8)                :: v1_t2, u1_y2
      real(8)                :: v0_t3, u0_y3
      real(8)                :: v1_t3, u1_y3
      real(8)                :: v0_t4, u0_y4
      real(8)                :: v1_t4, u1_y4
      real(8)                :: v0_t5, u0_y5
      real(8)                :: v1_t5, u1_y5
      real(8)                :: v0_t6, u0_y6
      real(8)                :: v1_t6, u1_y6
      real(8)                :: v0_t7, u0_y7
      real(8)                :: v1_t7, u1_y7

      integer                ::   LX, LY
      integer                ::   ii_1, ii_2, ii_3, ii_4, ii_5
      integer                ::   jj_1, jj_2, jj_3, jj_4, jj_5
      integer                ::   kk_1, kk_2, kk_3, kk_4, kk_5

#if __FUJITSU
! for FX10 DO_UNROLL=7 is the best
   #define         DO_UNROLL       7
#else
#if __IBM_REGISTER_VARS
! for BG/Q DO_UNROLL=? is the best
#   define         DO_UNROLL       3
#else
#   define         DO_UNROLL       5
#endif
#endif

! for FX10 DO_ITR=1 is the best
#define         DO_ITR          1

! for FX10 VLEN=(16*30) is the best
#if __FUJITSU
#   define         VLEN            (16*30)
#else
#if __IBM_REGISTER_VARS
#   define         VLEN            (16*30)
#else
#   define         VLEN            (16*30)
#endif
#endif

#if __FUJITSU
#   define         LX_ITR           6
#else
#   define         LX_ITR           6
#endif

      i_2 = n1
      i_3 = n2

!
! v:= Au
!
      if ( blk_0 == 0 ) then
         do i_1=i_2+local_rank,i_3,local_size
            j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
            j_3 = eigen_loop_end  (j, x_nnod,x_inod)
            j   = j+y_nnod*(DO_UNROLL*DO_ITR-1)
            j_4 = eigen_loop_end  (j, x_nnod,x_inod)
            do j_1=j_3+1,MIN(j_4,nm)
               a(j_1,i_1) = ZERO
            end do! j_1
         end do! i_1
!$OMP BARRIER
      end if
!
!
      LX = VLEN*LX_ITR; LY = 32*1
!
      l_2 = i_2; l_3 = i_3
      l_1 = l_2; l_4 = l_3

      k_2 = 1
      k   = eigen_translate_l2g(l_4, y_nnod,y_inod)
      k_3 = eigen_loop_end  (k-1, x_nnod,x_inod)
      do k_1=k_2,k_3,LX; k_4 = MIN(k_3,k_1+LX-1)

         j    = eigen_translate_l2g(k_1, x_nnod,x_inod)
         ii_2 = eigen_loop_start(j,y_nnod,y_inod)
         ii_2 = MAX(l_1,ii_2)
         ii_3 = l_4
         if ( ii_2 > ii_3 ) cycle
         ii_4 = MOD(ii_3-ii_2+1,DO_UNROLL*DO_ITR)+ii_2

         do i_1=ii_2+local_rank,ii_4-1,local_size

            j    = eigen_translate_l2g(i_1, y_nnod,y_inod)
            j    = j+(1-1)*y_nnod
            jj_2 = k_1
            jj_3 = eigen_loop_end  (j-1, x_nnod,x_inod)
            jj_3 = MIN(k_4,jj_3)
            if ( jj_2 > jj_3 ) cycle

            do kk_1=jj_2,jj_3,VLEN; kk_4=MIN(kk_1+VLEN-1,jj_3)

                  w0_(0) = v0_t(i_1+0)
                  w1_(0) = v1_t(i_1+0)

!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
               do j_1=kk_1,kk_4
!               do j_1=jj_2,jj_3

                  u0 = u0_x(j_1+0)
                  u1 = u1_x(j_1+0)

                  a0_0 = a(j_1+0,i_1+0)
                  w0_(0) = w0_(0)
     &                        + (a0_0*u0)
                  w1_(0) = w1_(0)
     &                        + (a0_0*u1)

               end do! j_1

                  v0_t(i_1+0) = w0_(0)
                  v1_t(i_1+0) = w1_(0)

#if defined(__INTEL_COMPILER)
               call MM_PREFETCH(a(kk_1,i_1),1)
#endif

               u0_y0 = u0_y(i_1+0)
               u1_y0 = u1_y(i_1+0)

!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
               do j_1=kk_1,kk_4
!               do j_1=jj_2,jj_3

                  v0 = ZERO
                  v1 = ZERO

                  a0_0 = a(j_1+0,i_1+0)
                  v0 = v0
     &                        + (a0_0*u0_y0)

                  v1 = v1
     &                        + (a0_0*u1_y0)

                  u0_t(j_1+0) = u0_t(j_1+0) + v0
                  u1_t(j_1+0) = u1_t(j_1+0) + v1

               end do! j_1

            enddo ! kk_1

         end do! i_1

         do i_0=ii_4+local_rank*DO_UNROLL*DO_ITR,
     &          ii_3,
     6          local_size*DO_UNROLL*DO_ITR

            j    = eigen_translate_l2g(i_0, y_nnod,y_inod)
            j    = j+(DO_UNROLL*DO_ITR-1)*y_nnod
            jj_2 = k_1
            jj_3 = eigen_loop_end  (j-1, x_nnod,x_inod)
            jj_3 = MIN(k_4,jj_3)
            if ( jj_2 > jj_3 ) cycle

            do kk_1=jj_2,jj_3,VLEN; kk_4=MIN(kk_1+VLEN-1,jj_3)

            do i_1 = i_0,i_0+DO_UNROLL*DO_ITR-1,DO_UNROLL

                  w0_(0) = v0_t(i_1+0)
                  w1_(0) = v1_t(i_1+0)
                  w0_(1) = v0_t(i_1+1)
                  w1_(1) = v1_t(i_1+1)

! Next two loops might be fused by compiler, however, it is intended
! that the first loop performs as a cache prefetcher.
! Please care compiler options or directives not to be fused.

!IBM* ASSERT(NODEPS)
!IBM* INDEPENDENT
!OCL LOOP_NOFUSION
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
               do j_1=kk_1,kk_4
!               do j_1=jj_2,jj_3

                  u0 = u0_x(j_1+0)
                  u1 = u1_x(j_1+0)

                  a0_0 = a(j_1+0,i_1+0)
                  a0_1 = a(j_1+0,i_1+1)

                  w0_(0) = w0_(0)
     &                        + (a0_0*u0)
                  w1_(0) = w1_(0)
     &                        + (a0_0*u1)

                  w0_(1) = w0_(1)
     &                        + (a0_1*u0)
                  w1_(1) = w1_(1)
     &                        + (a0_1*u1)

               end do! j_1

               v0_t(i_1+0) = w0_(0)
               v1_t(i_1+0) = w1_(0)
               v0_t(i_1+1) = w0_(1)
               v1_t(i_1+1) = w1_(1)

#if defined(__INTEL_COMPILER)
               call MM_PREFETCH(a(kk_1,i_1),1)
#endif
#if __IBM_REGISTER_VARS
!IBM* PREFETCH_FOR_LOAD(a(kk_1,i_1))
#endif

#if DO_UNROLL>=3
               w0_(2) = v0_t(i_1+2)
               w1_(2) = v1_t(i_1+2)
#endif
#if DO_UNROLL>=4
               w0_(3) = v0_t(i_1+3)
               w1_(3) = v1_t(i_1+3)
#endif
#if DO_UNROLL>=5
               w0_(4) = v0_t(i_1+4)
               w1_(4) = v1_t(i_1+4)
#endif
#if DO_UNROLL>=6
               w0_(5) = v0_t(i_1+5)
               w1_(5) = v1_t(i_1+5)
#endif
#if DO_UNROLL>=7
               w0_(6) = v0_t(i_1+6)
               w1_(6) = v1_t(i_1+6)
#endif
#if DO_UNROLL>=8
               w0_(7) = v0_t(i_1+7)
               w1_(7) = v1_t(i_1+7)
#endif

               u0_y0 = u0_y(i_1+0)
               u0_y1 = u0_y(i_1+1)
#if DO_UNROLL>=3
               u0_y2 = u0_y(i_1+2)
#endif
#if DO_UNROLL>=4
               u0_y3 = u0_y(i_1+3)
#endif
#if DO_UNROLL>=5
               u0_y4 = u0_y(i_1+4)
#endif
#if DO_UNROLL>=6
               u0_y5 = u0_y(i_1+5)
#endif
#if DO_UNROLL>=7
               u0_y6 = u0_y(i_1+6)
#endif
#if DO_UNROLL>=8
               u0_y7 = u0_y(i_1+7)
#endif
               u1_y0 = u1_y(i_1+0)
               u1_y1 = u1_y(i_1+1)
#if DO_UNROLL>=3
               u1_y2 = u1_y(i_1+2)
#endif
#if DO_UNROLL>=4
               u1_y3 = u1_y(i_1+3)
#endif
#if DO_UNROLL>=5
               u1_y4 = u1_y(i_1+4)
#endif
#if DO_UNROLL>=6
               u1_y5 = u1_y(i_1+5)
#endif
#if DO_UNROLL>=7
               u1_y6 = u1_y(i_1+6)
#endif
#if DO_UNROLL>=8
               u1_y7 = u1_y(i_1+7)
#endif

!IBM* ASSERT(NODEPS)
!OCL LOOP_NOFUSION
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
               do j_1=kk_1,kk_4
!               do j_1=jj_2,jj_3

                  u0 = u0_x(j_1+0)
                  u1 = u1_x(j_1+0)

                  a0_0 = a(j_1+0,i_1+0)
                  a0_1 = a(j_1+0,i_1+1)

                  v0 = ZERO
                  v1 = ZERO

                  v0 = v0
     &                        + (a0_0*u0_y0)
     &                        + (a0_1*u0_y1)
                  v1 = v1
     &                        + (a0_0*u1_y0)
     &                        + (a0_1*u1_y1)

#if DO_UNROLL>=3
                  a0_2 = a(j_1+0,i_1+2)
#endif
#if DO_UNROLL>=4
                  a0_3 = a(j_1+0,i_1+3)
#endif

#if DO_UNROLL>=3
                  w0_(2) = w0_(2)
     &                        + (a0_2*u0)
                  w1_(2) = w1_(2)
     &                        + (a0_2*u1)
#endif
#if DO_UNROLL>=4
                  w0_(3) = w0_(3)
     &                        + (a0_3*u0)
                  w1_(3) = w1_(3)
     &                        + (a0_3*u1)
#endif

#if DO_UNROLL>=3
                  v0 = v0
     &                        + (a0_2*u0_y2)
#endif
#if DO_UNROLL>=4
     &                        + (a0_3*u0_y3)
#endif
#if DO_UNROLL>=3
                  v1 = v1
     &                        + (a0_2*u1_y2)
#endif
#if DO_UNROLL>=4
     &                        + (a0_3*u1_y3)
#endif

#if DO_UNROLL>=5
                  a0_4 = a(j_1+0,i_1+4)
#endif
#if DO_UNROLL>=6
                  a0_5 = a(j_1+0,i_1+5)
#endif

#if DO_UNROLL>=5
                  w0_(4) = w0_(4)
     &                        + (a0_4*u0)
                  w1_(4) = w1_(4)
     &                        + (a0_4*u1)
#endif
#if DO_UNROLL>=6
                  w0_(5) = w0_(5)
     &                        + (a0_5*u0)
                  w1_(5) = w1_(5)
     &                        + (a0_5*u1)
#endif

#if DO_UNROLL>=5
                  v0 = v0
     &                        + (a0_4*u0_y4)
#endif
#if DO_UNROLL>=6
     &                        + (a0_5*u0_y5)
#endif
#if DO_UNROLL>=5
                  v1 = v1
     &                        + (a0_4*u1_y4)
#endif
#if DO_UNROLL>=6
     &                        + (a0_5*u1_y5)
#endif

#if DO_UNROLL>=7
                  a0_6 = a(j_1+0,i_1+6)
#endif
#if DO_UNROLL>=8
                  a0_7 = a(j_1+0,i_1+7)
#endif
#if DO_UNROLL>=7
                  w0_(6) = w0_(6)
     &                        + (a0_6*u0)
                  w1_(6) = w1_(6)
     &                        + (a0_6*u1)
#endif
#if DO_UNROLL>=8
                  w0_(7) = w0_(7)
     &                        + (a0_7*u0)
                  w1_(7) = w1_(7)
     &                        + (a0_7*u1)
#endif
                  v0 = v0
#if DO_UNROLL>=7
     &                        + (a0_6*u0_y6)
#endif
#if DO_UNROLL>=8
     &                        + (a0_7*u0_y7)
#endif
                  v1 = v1
#if DO_UNROLL>=7
     &                        + (a0_6*u1_y6)
#endif
#if DO_UNROLL>=8
     &                        + (a0_7*u1_y7)
#endif

                  u0_t(j_1+0) = u0_t(j_1+0) + v0
                  u1_t(j_1+0) = u1_t(j_1+0) + v1

               end do! j_1

#if DO_UNROLL>=3
               v0_t(i_1+2) = w0_(2)
               v1_t(i_1+2) = w1_(2)
#endif
#if DO_UNROLL>=4
               v0_t(i_1+3) = w0_(3)
               v1_t(i_1+3) = w1_(3)
#endif
#if DO_UNROLL>=5
               v0_t(i_1+4) = w0_(4)
               v1_t(i_1+4) = w1_(4)
#endif
#if DO_UNROLL>=6
               v0_t(i_1+5) = w0_(5)
               v1_t(i_1+5) = w1_(5)
#endif
#if DO_UNROLL>=7
               v0_t(i_1+6) = w0_(6)
               v1_t(i_1+6) = w1_(6)
#endif
#if DO_UNROLL>=8
               v0_t(i_1+7) = w0_(7)
               v1_t(i_1+7) = w1_(7)
#endif

           end do! i_1
           end do! kk_1

         end do! i_0

      end do! k_1


      return
      end subroutine  eigen_prd_au_body1

      subroutine  eigen_prd_au_body2(
     &               u_t, u0_z, u1_z, v_t, v0_z, v1_z, nv,
     &               x_pos, y_pos
     &               ,local_rank, local_size
     &               )
      implicit NONE

      integer, intent(in)    :: nv, x_pos, y_pos
      real(8), intent(out)   :: u_t(1:nv,*),  v_t(1:nv,*)
      real(8), intent(in)    :: u0_z(1:nv,*), v0_z(1:nv,*)
      real(8), intent(in)    :: u1_z(1:nv,*), v1_z(1:nv,*)

      integer                :: local_rank, local_size

      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: jj_1, jj_2, jj_3, jj_4
      integer                :: i, j, k

      integer, parameter     :: LX = 1024

      integer                :: ll_size, ll_rank

      ll_rank =    local_rank
      ll_size =    local_size
      if ( local_size > 1 ) then
         ll_rank =    ll_rank - 1
         ll_size =    ll_size - 1
      end if
      if ( ll_rank < 0 ) return

      jj_1 = x_pos
      jj_2 = (jj_1-1)/ll_size+1
      jj_3 =    (jj_2*(ll_rank+0)     )+1
      jj_4 = MIN(jj_2*(ll_rank+1),jj_1)

      do jj_1=jj_3,jj_4,LX
         j_3=jj_1; j_4=MIN(jj_1+LX-1,jj_4)
         if ( MOD(local_size,2) == 1 ) then
            do j_1=j_3,j_4
               u_t(j_1,1) = u0_z(j_1,1)
               u_t(j_1,2) = u1_z(j_1,1)
            end do
         else
            do j_1=j_3,j_4
               u_t(j_1,1) = u0_z(j_1,1)+u0_z(j_1,2)
               u_t(j_1,2) = u1_z(j_1,1)+u1_z(j_1,2)
            end do
         end if
         if ( local_size > 2 ) then
         do j=3-MOD(local_size,2),local_size,2
            do j_1=j_3,j_4
               u_t(j_1,1) = u_t(j_1,1)+u0_z(j_1,j+0)+u0_z(j_1,j+1)
               u_t(j_1,2) = u_t(j_1,2)+u1_z(j_1,j+0)+u1_z(j_1,j+1)
            end do
         end do
         end if
      end do

      jj_1 = y_pos
      jj_2 = (jj_1-1)/ll_size+1
      jj_3 =    (jj_2*(ll_rank+0)     )+1
      jj_4 = MIN(jj_2*(ll_rank+1),jj_1)

      do jj_1=jj_3,jj_4,LX
         j_3=jj_1; j_4=MIN(jj_1+LX-1,jj_4)
         if ( MOD(local_size,2) == 1 ) then
            do j_1=j_3,j_4
               v_t(j_1,1) = v0_z(j_1,1)
               v_t(j_1,2) = v1_z(j_1,1)
            end do
         else
            do j_1=j_3,j_4
               v_t(j_1,1) = v0_z(j_1,1)+v0_z(j_1,2)
               v_t(j_1,2) = v1_z(j_1,1)+v1_z(j_1,2)
            end do
         end if
         if ( local_size > 2 ) then
         do j=3-MOD(local_size,2),local_size,2
            do j_1=j_3,j_4
               v_t(j_1,1) = v_t(j_1,1)+v0_z(j_1,j+0)+v0_z(j_1,j+1)
               v_t(j_1,2) = v_t(j_1,2)+v1_z(j_1,j+0)+v1_z(j_1,j+1)
            end do
         end do
         end if
      end do

      return
      end subroutine  eigen_prd_au_body2

      subroutine  eigen_prd_au_body3(
     &               u_x, u_y, v_x,
     &               u_t, v_t, d_t,
     &               n1, n2, x_pos, y_pos, nv
     &               )
      use eigen_devel
      use eigen_libs
      implicit NONE

      integer, intent(in)    :: nv, n1, n2
      real(8), intent(in)    :: u_x(1:nv,*),u_y(1:nv,*)
      real(8), intent(inout) :: v_x(1:nv,*)
      real(8), intent(in)    :: u_t(1:nv,*),v_t(1:nv,*)
      real(8), intent(in)    :: d_t(1:nv)
      integer                :: x_pos, y_pos

      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: i, j, k
      integer                :: nm1, nm2


      i_2 = n1
      i_3 = n2

      if ( diag_0 > 0 ) then

         j = eigen_translate_l2g(diag_0, y_nnod,y_inod)
         j = eigen_translate_g2l(j, x_nnod,x_inod)
         if ( j > nv ) return

         nm1 = y_nnod/n_common
         nm2 = x_nnod/n_common

         do i=1,MBAND

            if ( nm2 == 1 ) then
            if ( nm1 == 1 ) then
               call eigen_prd_au_body3_sub11(
     &            v_x(j,i), v_t(diag_0,i), d_t(diag_0), u_y(diag_0,i),
     &            (i_3-diag_0)/(x_nnod/n_common)+1, nm1, nm2)
            else
               call eigen_prd_au_body3_subX1(
     &            v_x(j,i), v_t(diag_0,i), d_t(diag_0), u_y(diag_0,i),
     &            (i_3-diag_0)/(x_nnod/n_common)+1, nm1, nm2)
            end if
            else
               call eigen_prd_au_body3_subXX(
     &            v_x(j,i), v_t(diag_0,i), d_t(diag_0), u_y(diag_0,i),
     &            (i_3-diag_0)/(x_nnod/n_common)+1, nm1, nm2)
            end if

          end do! i

       end if


       return
       end subroutine  eigen_prd_au_body3

      subroutine eigen_prd_au_body3_subXX(v_x,v_t,d_t,u_y, n,nm1,nm2)
      implicit NONE

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: v_x(nm1,*)
      real(8), intent(in)    :: v_t(nm2,*)
      real(8), intent(in)    :: d_t(nm2,*)
      real(8), intent(in)    :: u_y(nm2,*)

      integer                :: i


*soption unroll(4)
!DIR$ VECTOR ALWAYS
      do i=1,n
         v_x(1,i) = v_x(1,i)+v_t(1,i)+d_t(1,i)*u_y(1,i)
      end do! i


      return
      end subroutine  eigen_prd_au_body3_subXX

      subroutine eigen_prd_au_body3_subX1(v_x,v_t,d_t,u_y, n,nm1,nm2)
      implicit NONE

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: v_x(nm1,*)
      real(8), intent(in)    :: v_t(*)
      real(8), intent(in)    :: d_t(*)
      real(8), intent(in)    :: u_y(*)

      integer                :: i


*soption unroll(4)
!DIR$ VECTOR ALWAYS
      do i=1,n
         v_x(1,i) = v_x(1,i)+v_t(i)+d_t(i)*u_y(i)
      end do! i


      return
      end subroutine  eigen_prd_au_body3_subX1

      subroutine eigen_prd_au_body3_sub11(v_x,v_t,d_t,u_y, n,nm1,nm2)
      implicit NONE

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: v_x(1:n)
      real(8), intent(in)    :: v_t(1:n)
      real(8), intent(in)    :: d_t(1:n)
      real(8), intent(in)    :: u_y(1:n)

      integer                :: i


*soption unroll(4)
!DIR$ VECTOR ALWAYS
      do i=1,n
         v_x(i) = v_x(i)+v_t(i)+d_t(i)*u_y(i)
      end do! i


      return
      end subroutine  eigen_prd_au_body3_sub11

