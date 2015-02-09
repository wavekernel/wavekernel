      subroutine eigen_prd_ux(
     &              a, nm,
     &              u_x, u_y, nv,
     &              u_t, v_t, i, c, e, ne)
!$    use OMP_LIB
      use MPI
      use eigen_devel
      use eigen_libs
      use eigen_house, only : tsave
      implicit NONE

      integer, intent(in)    :: nm, nv, ne, i
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(out)   :: u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(out)   :: u_t(*), v_t(*)

      real(8), intent(out)   :: c(MBAND,MBAND), e(1:ne,*)

      real(8)                :: g_n(MBAND), a_n, anorm2, beta
      integer, parameter     :: LWK = MBAND*(MBAND+1)
      real(8)                :: t, tt(LWK), s, ss(LWK)

!      real(8), save :: tsave(0:15, 0:MBAND)
      save :: tt, ss, anorm2

      integer                :: x_owner_nod, x_pos
      integer                :: y_owner_nod, y_pos
      integer                :: i_1, i_2, i_3, i_4, i_5
      integer                :: j_1, j_2, j_3, j_4, j_5
      integer                :: k_1, k_2, k_3, k_4, k_5
      integer                :: jj_1, jj_2, jj_3
      integer                :: L, LL

      real(8)                :: t11, t12, t22, tol
      real(8)                :: a1, a2
      real(8)                :: aL(MBAND,MBAND), aX(MBAND, MBAND)

!      integer, parameter     :: VLEN = 2048
      integer, parameter     :: VLEN = 512

      real(8) :: d1,d2


      L  = i-MBAND
      LL = L-MBAND

!$OMP MASTER

      c(1:MBAND, 1:MBAND) = ZERO

#if _DEBUG_
      d1=eigen_get_wtime()
#endif
      i_2 = MBAND
      i_3 = MAX(1, 1-LL)
#if 0
      do i_1=i_2,i_3,-1
!      do i_1=MBAND,1,-1
!      if ( LL + i_1 >= 1  ) then
         y_owner_nod = eigen_owner_node   (L+i_1, y_nnod, y_inod)
         x_pos       = eigen_translate_g2l(i-1,   x_nnod, x_inod)
         call bcast_dbl(a(1,i_1), x_pos, y_owner_nod, y_COMM_WORLD)
!      end if
      end do
#else
      i_1 = i_3
      y_owner_nod = eigen_owner_node   (L+i_1, y_nnod, y_inod)
      x_pos       = eigen_translate_g2l(i-1,   x_nnod, x_inod)
      call bcastw_dbl( a(1,i_1), x_pos, y_owner_nod,
     $                nm, i_2-i_3+1, u_t, y_COMM_WORLD )
#endif
#if _DEBUG_
      d2=eigen_get_wtime()
#if _DEBUG2_
      rt_ptr(1)=rt_ptr(1)+1
      rt_timer(rt_ptr(1),1)=d2-d1

      rt_ptr(2)=rt_ptr(2)+1
#endif
#endif

!$OMP END MASTER

!$OMP BARRIER

      i_2 = MBAND
      i_3 = MAX(1, 1-LL)
!OCL NOFP_RELAXED
      do i_1=i_2,i_3,-1
!      do i_1=MBAND,1,-1
!      if ( LL + i_1 >= 1  ) then

         j_2 = eigen_loop_start(1,      x_nnod, x_inod)
         j_3 = eigen_loop_end  (LL+i_1, x_nnod, x_inod)

!$       if ( omp_get_num_threads() /= 1 ) then
!$          tsave(omp_get_thread_num(), 0:i_1) = ZERO
!$       else
!$OMP MASTER
         tt(1:i_1) = ZERO
!$OMP END MASTER
!$       end if

!$OMP BARRIER

         do jj_1=j_2,j_3,VLEN; jj_2=jj_1; jj_3=MIN(jj_1+VLEN-1,j_3)

             t = ZERO
!$OMP DO
!IBM* ASSERT(NODEPS)
             do j_1=jj_2,jj_3
                s = a(j_1, i_1)
                u_x(j_1, i_1) = s
                t = t + s * s
             end do! j_1
!$OMP END DO
!-
!$        if ( omp_get_num_threads() /= 1 ) then
!$           tsave(omp_get_thread_num(), 0)
!$   +        = tsave(omp_get_thread_num(), 0) + t
!$        else
             tt(1) = tt(1) + t
!$        end if
!-
          do k_1=1,i_1-1
             t = ZERO
!$OMP DO
!IBM* ASSERT(NODEPS)
             do j_1=jj_2,jj_3
                t = t + u_x(j_1, i_1)*a(j_1, k_1)
             enddo! j_1
!$OMP END DO
!-
!$        if ( omp_get_num_threads() /= 1 ) then
!$           tsave(omp_get_thread_num(), k_1)
!$   +        = tsave(omp_get_thread_num(), k_1) + t
!$        else
             tt(k_1+1) = tt(k_1+1) + t
!$        end if
!-
          end do! k_1

         end do! jj_1

!$OMP BARRIER

!$        if ( omp_get_num_threads() /= 1 ) then
!$OMP MASTER
!$           t = ZERO
!$           do j_1=0,omp_get_num_threads()-1
!$              t = t + tsave(j_1, 0)
!$           end do
!$           tt(1) = t
!$        do k_1=1,i_1-1
!$           t = ZERO
!$           do j_1=0,omp_get_num_threads()-1
!$              t = t + tsave(j_1, k_1)
!$           end do ! j_1
!$           tt(k_1+1) = t
!$        end do ! k_1
!$OMP END MASTER
!$        end if

!$OMP BARRIER

!$OMP MASTER
         if ( i_1 == MBAND ) then
            do k_1=MBAND,1,-1
            if ( LL + k_1 >= 1  ) then
               do j_1=1,MBAND
                  x_pos    = eigen_translate_g2l(LL+j_1, x_nnod, x_inod)
                  x_owner_nod = eigen_owner_node(LL+j_1, x_nnod, x_inod)
                  if ( x_inod == x_owner_nod ) then
                     t = a(x_pos, k_1)
                  else
                     t = ZERO
                  endif
                  tt(j_1+MBAND*k_1) = t
               end do! j_1
            else
               tt(1+MBAND*k_1:MBAND+MBAND*k_1) = ZERO
            endif
            end do! k_1
#if _DEBUG_
      d1=eigen_get_wtime()
#endif
            call reduce_dbl(tt, ss, LWK, 1, x_COMM_WORLD)
#if _DEBUG_
      d2=eigen_get_wtime()
#if _DEBUG2_
      rt_timer(rt_ptr(2),2)=rt_timer(rt_ptr(2),2)+d2-d1
#endif
#endif
            do k_1=1,MBAND
               do j_1=1,MBAND
                  aL(j_1, k_1) = tt(j_1+MBAND*k_1)
               end do! j_1
            end do! k_1
         else
#if _DEBUG_
      d1=eigen_get_wtime()
#endif
            call reduce_dbl(tt, ss, i_1, 1, x_COMM_WORLD)
#if _DEBUG_
      d2=eigen_get_wtime()
#if _DEBUG2_
      rt_timer(rt_ptr(2),2)=rt_timer(rt_ptr(2),2)+d2-d1
#endif
#endif
         end if

         anorm2 = tt(1)
         do k_1=1,i_1-1
            aX(k_1, i_1) = tt(k_1+1)
         end do! k_1
!$OMP END MASTER

! // barrier to share 'anorm2'
!$OMP BARRIER

         if ( anorm2 > ZERO ) then

!$OMP MASTER
            a_n  =  aL(i_1,i_1)
            g_n(i_1)  = -SIGN(SQRT(anorm2), a_n)
            beta =  anorm2 - a_n * g_n(i_1)
            a_n  =  a_n - g_n(i_1)

            x_owner_nod = eigen_owner_node   (LL+i_1, x_nnod, x_inod)
            x_pos       = eigen_translate_g2l(LL+i_1, x_nnod, x_inod)
            if ( x_inod == x_owner_nod ) then
               u_x(x_pos, i_1) =  a_n
            end if

            do k_1=1,i_1-1
               tt(k_1) = (aX(k_1, i_1) - g_n(i_1) * aL(i_1, k_1))/ beta
            end do
!$OMP END MASTER

! // barrier to share 'tt()'
!$OMP BARRIER

!OCL NOFP_CONTRACT
            do k_1=1,i_1-1
               t = tt(k_1)
!$OMP MASTER
!IBM* ASSERT(NODEPS)
               do j_1=1,i_1-1
                  aL(j_1, k_1) = aL(j_1, k_1) - t * aL (j_1, i_1)
               end do! j_1
!$OMP END MASTER
!$OMP DO
!IBM* ASSERT(NODEPS)
             do j_1=j_2,j_3
                a (j_1, k_1) = a (j_1, k_1) - t * u_x(j_1, i_1)
             end do! j_1
!$OMP ENDDO
            end do! k_1

!$OMP MASTER
            if ( x_inod == x_owner_nod ) then
               a(x_pos, i_1) = a_n
            end if
!$OMP END MASTER

         else

!$OMP MASTER
            g_n(i_1)  = ZERO
            beta = ONE
!$OMP END MASTER

         end if

!$OMP MASTER
         e(i-MBAND+i_1, MBAND) = g_n(i_1)
         c(i_1, i_1) = ONE / beta
!$OMP END MASTER

!      end if
      end do! i_1


!$OMP MASTER

      y_owner_nod = eigen_owner_node(i, y_nnod, y_inod)
      if ( y_owner_nod == y_inod ) then
      do i_1=1,MBAND
         if ( i+i_1 > MBAND+1 ) then
           x_owner_nod = eigen_owner_node(i-MBAND-1+i_1, x_nnod, x_inod)
            if ( x_owner_nod == x_inod ) then
              x_pos = eigen_translate_g2l(i-MBAND-1+i_1, x_nnod, x_inod)
               e(i-MBAND+i_1, 1) = a(x_pos, i_1)
            end if
         end if
      end do! i_1
      end if


      i_2 = MBAND
      i_3 = MAX(1, 1-LL)
      do i_1=i_3,i_2
!      do i_1=1,MBAND
!      if ( LL + i_1 >= 1  ) then
         j_2 = eigen_loop_start(LL+i_1+1, x_nnod, x_inod)
         j_3 = eigen_loop_end  (i-1,      x_nnod, x_inod)
         u_x(j_2:j_3,i_1) = ZERO

         j_2 = eigen_loop_start(LL+i_1+1, x_nnod, x_inod)
         j_3 = eigen_loop_end  (LL+i_1+MBAND-1, x_nnod, x_inod)
         a  (j_2:j_3,i_1) = ZERO

!      end if
      end do! i_1


#if _DEBUG_
      d1=eigen_get_wtime()
#endif
      x_pos = eigen_translate_g2l(i-1, x_nnod, x_inod)
      if ( MOD(MBAND,2)==1 ) then; i_1 = 1
         call datacast_dbl (u_y(1, i_1), u_x(1, i_1),   u_t, v_t, x_pos)
      end if
      do i_1=MOD(MBAND,2)+1,MBAND,2
         call datacast_dbl2(u_y(1, i_1), u_y(1, i_1+1),
     &                      u_x(1, i_1), u_x(1, i_1+1), u_t, v_t, x_pos)
      end do! i_1
#if _DEBUG_
      d2=eigen_get_wtime()
#if _DEBUG2_
      rt_ptr(3)=rt_ptr(3)+1
      rt_timer(rt_ptr(3),3)=d2-d1
#endif
#endif

!$OMP END MASTER

!
! [u1, u2] = [a1, a2]|1 0| + [e1, e2]|a 0| - [u1, u2]|0 b|
!                    |0 1|           |0 c|           |0 0|
!
!            a = g_n(1)
!            c = g_n(2)
!            b = S1(1,2)
!
! [u1, u2]|1 b| = [a1, a2]|1 0| + [e1, e2]|a 0|
!         |  1|           |0 1|           |0 c|
!
! [u1, u2] = [a1, a2]|1 -b| + [e1, e2]|a c-ab| = A * S1 + E * S2
!                    |0  1|           |0 c   |
!
      return
      end subroutine eigen_prd_ux

