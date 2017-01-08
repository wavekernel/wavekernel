       subroutine eigen_trd_local_2update(
     &               w, nm,
     &               ux, uy, vx, vy, nv,
     &               i_base, i, ix)
!$     use OMP_LIB
       use eigen_devel
       use eigen_libs
       implicit none
!
       integer, intent(in)    ::  nm, nv
       real(8), intent(inout) ::  w(1:nm,*)
       real(8), intent(in)    ::  ux(1:nv,*),uy(1:nv,*)
       real(8), intent(in)    ::  vx(1:nv,*),vy(1:nv,*)
       integer, intent(in)    ::  i_base, i, ix
!
       integer                ::  k_1
       integer                ::  j, k, l
!
       integer                ::  i_1, i_2, i_3, i_4
       integer                ::  j_1, j_2, j_3, j_4
       integer                ::  l_1, l_2, l_3, l_4
       integer                ::  jj_1, jj_2, jj_3, jj_4
       integer                ::  LX
!
       include 'CSTAB.h'
!
       real(8)                :: u_x, v_x
       real(8)                :: uy0, vy0
       real(8)                :: uy1, vy1
       real(8)                :: uy2, vy2
       real(8)                :: uy3, vy3
       real(8)                :: w0
       real(8)                :: w1
       real(8)                :: w2
       real(8)                :: w3

          if ( i - i_base <= 1 ) return
!
          k_1 = 1
!
          LX  = L1_LSIZE*L1_WAY/16

          i_2 = eigen_loop_start(i_base+1, y_nnod,y_inod)
          i_3 = eigen_loop_end  (i-1,      y_nnod,y_inod)
          i_4 = MOD(i_3-i_2+1, 4)+i_2
          if ( i_2 > i_3 ) return
!
          L = ix - 1
          jj_2 = eigen_loop_start(1, x_nnod,x_inod)
          jj_3 = eigen_loop_end  (L, x_nnod,x_inod)

          do jj_1=jj_2,jj_3,LX
             j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)
             do i_1=i_2,i_4-1                         ! 0
                j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
                l_1 = j-i_base
                uy0 = uy(i_1+0,k_1)
                vy0 = vy(i_1+0,k_1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                do j_1=j_2,j_3
                   u_x = ux(j_1,k_1)
                   v_x = vx(j_1,k_1)
                   w0 = w(j_1,l_1+0*y_nnod)
                   w0 = w0
     &                    -(u_x*vy0)
                   w0 = w0
     &                    -(v_x*uy0)
                   w(j_1,l_1+0*y_nnod) = w0
                end do! j_1
             end do! l_1
             do i_1=i_4,i_3,4                     ! 3
                j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
                l_1 = j-i_base
                uy0 = uy(i_1+0,k_1)
                vy0 = vy(i_1+0,k_1)
                uy1 = uy(i_1+1,k_1)
                vy1 = vy(i_1+1,k_1)
                uy2 = uy(i_1+2,k_1)
                vy2 = vy(i_1+2,k_1)
                uy3 = uy(i_1+3,k_1)
                vy3 = vy(i_1+3,k_1)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                do j_1=j_2,j_3
                   u_x = ux(j_1,k_1)
                   v_x = vx(j_1,k_1)
                   w0 = w(j_1,l_1+0*y_nnod)
                   w1 = w(j_1,l_1+1*y_nnod)
                   w0 = w0
     &                    -(u_x*vy0)
                   w1 = w1
     &                    -(u_x*vy1)
                   w0 = w0
     &                    -(v_x*uy0)
                   w1 = w1
     &                    -(v_x*uy1)
                   w(j_1,l_1+0*y_nnod) = w0
                   w(j_1,l_1+1*y_nnod) = w1
                   w2 = w(j_1,l_1+2*y_nnod)
                   w3 = w(j_1,l_1+3*y_nnod)
                   w2 = w2
     &                    -(u_x*vy2)
                   w3 = w3
     &                    -(u_x*vy3)
                   w2 = w2
     &                    -(v_x*uy2)
                   w3 = w3
     &                    -(v_x*uy3)
                   w(j_1,l_1+2*y_nnod) = w2
                   w(j_1,l_1+3*y_nnod) = w3
                end do! j_1
             end do! l_1
          end do! jj_1
!
       return
       end subroutine  eigen_trd_local_2update

