      subroutine eigen_prd_load(a, w, nm,
     &              d_t,
     &              u_x, u_y, v_x, v_y, nv,
     &              m0, i_base, i_block)
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      implicit NONE

      integer, intent(in)    :: nm, nv, m0, i_base, i_block
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(out)   :: w(1:nm,*)
      real(8), intent(out)   :: u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(out)   :: v_x(1:nv,*), v_y(1:nv,*)
      real(8), intent(out)   :: d_t(*)

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3, k_4
      integer                ::  i, j, L


      i_2 = eigen_loop_start(i_base+1,  y_nnod,y_inod)
      i_3 = eigen_loop_end  (i_base+m0, y_nnod,y_inod)
      j_2 = eigen_loop_start(1,         x_nnod,x_inod)
      j_3 = eigen_loop_end  (i_base+m0, x_nnod,x_inod)

      do i_1=i_2,i_3
         j = eigen_translate_l2g(i_1, y_nnod,y_inod)
         do j_1=j_2,j_3
            w(j_1,j-i_base) = a(j_1,i_1)
         end do! j_1
      end do! i_1

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = eigen_loop_end  (i_base+m0-1, y_nnod,y_inod)
         i_4 = x_nnod/n_common
         if ( y_nnod == x_nnod ) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            do i_1=i_2,i_3,i_4
               d_t(i_1)   = a(i_1,i_1)
               a(i_1,i_1) = ZERO
            end do! i_1
         else
            do i_1=i_2,i_3,i_4
               j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
               j_1 = eigen_translate_g2l(j, x_nnod,x_inod)
               d_t(i_1)   = a(j_1,i_1)
               a(j_1,i_1) = ZERO
            end do! i_1
         end if
      end if

      i = i_base+m0
      L = i - 2
      y_pos  = eigen_translate_g2l(i, y_nnod,y_inod)
      x_pos  = eigen_translate_g2l(i, x_nnod,x_inod)

      k_2 = m0
      k_3 = MAX(1,3*(2-i_block))
      do k_1=k_3,k_2
         do j_1=1,x_pos
            u_x(j_1,k_1) = ZERO
            v_x(j_1,k_1) = ZERO
         end do! j_1
         do j_1=1,y_pos
            u_y(j_1,k_1) = ZERO
            v_y(j_1,k_1) = ZERO
         end do! j_1
      end do! k_1

      return
      end subroutine  eigen_prd_load


      subroutine eigen_prd_store(a, w, nm,
     &              d_t,
     &              m0, i_base)
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      implicit NONE

      integer, intent(in)    :: nm, m0, i_base
      real(8), intent(out)   :: a(1:nm,*)
      real(8), intent(in)    :: w(1:nm,*), d_t(*)

      integer                ::  i, i_1, i_2, i_3, i_4
      integer                ::  j, j_1, j_2, j_3, j_4


      i_2 = eigen_loop_start(i_base+1,  y_nnod,y_inod)
      i_3 = eigen_loop_end  (i_base+m0, y_nnod,y_inod)
      do i_1=i_2,i_3
         j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
         j_3 = eigen_loop_end  (j, x_nnod,x_inod)
         do j_1=1,j_3
            a(j_1,i_1) = w(j_1,j-i_base)
         end do! j_1
      end do! i_1

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = eigen_loop_end  (i_base, y_nnod,y_inod)
         i_4 = x_nnod/n_common
         if ( i_2 <= i_3 ) then
         if ( y_nnod == x_nnod ) then
            do i_1=i_2,i_3,i_4
               j_1 = i_1
               a(j_1,i_1) = d_t(i_1)
            end do! i_1
         else
            do i_1=i_2,i_3,i_4
!              j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
!              j_1 = eigen_translate_g2l(j, x_nnod,x_inod)
               j   = (i_1-1)*y_nnod+y_inod
               j_1 = (j-1)/x_nnod+1
               a(j_1,i_1) = d_t(i_1)
            end do! i_1
         end if
         end if
      end if


      return
      end subroutine  eigen_prd_store

