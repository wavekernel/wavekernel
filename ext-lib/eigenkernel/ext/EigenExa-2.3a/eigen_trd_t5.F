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

       integer                ::  local_size, ll_size
       integer                ::  local_rank, ll_rank


          if ( i - i_base <= 1 ) return

          local_rank = 0
          local_size = 1
!$        local_rank = omp_get_thread_num()
!$        local_size = omp_get_num_threads()

!$        if ( local_size > 1 ) then
!$           ll_size  = local_size - 1
!$           ll_rank  = local_rank - 1
!$        else
             ll_size  = local_size
             ll_rank  = local_rank
!$        endif
!
          k_1 = 1
!
!          LX  = L1_LSIZE*L1_WAY/16
          LX  = 64

          i_2 = eigen_loop_start(i_base+1, y_nnod,y_inod)
          i_3 = eigen_loop_end  (i-1,      y_nnod,y_inod)
          if ( i_2 > i_3 ) return
!
          L = ix - 1
          jj_2 = eigen_loop_start(1, x_nnod,x_inod)
          jj_3 = eigen_loop_end  (L, x_nnod,x_inod)

          do jj_1=jj_2+LX*ll_rank,jj_3,LX*ll_size
             j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

             do i_1=i_2,i_3
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
                   w0 = u_x*vy0 + v_x*uy0
                   w(j_1,l_1) = w(j_1,l_1) - w0
                end do! j_1

             end do! l_1

          end do! jj_1
!
       return
       end subroutine  eigen_trd_local_2update

