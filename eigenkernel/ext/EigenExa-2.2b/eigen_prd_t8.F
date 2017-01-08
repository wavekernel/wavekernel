      subroutine eigen_prd_init(a, nm, n,
     &              d_out, e_out, ne,
     &              u_t, v_t, nv)
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      use eigen_house
      implicit NONE

      integer, intent(in)    :: nm, n, ne, nv
      real(8), intent(inout) :: a(1:nm, *)
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne, *)
      real(8), intent(in)    :: u_t(*), v_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j

      integer                :: local_size, local_rank

      include 'CSTAB.h'


      d_out(1:n)    = ZERO
      e_out(1:n,1)  = ZERO
      e_out(1:n,2)  = ZERO

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = eigen_loop_end  (n, y_nnod, y_inod)
         if ( i_2 <= i_3 ) then
            i_4 = x_nnod/n_common
            j_2 = diag_1
            j_3 = eigen_loop_end  (n, x_nnod, x_inod)
            j_4 = y_nnod/n_common
            k_2 = 0
            k_3 = (i_3-i_2)/i_4
            do  k_1 = k_2, k_3
               i_1 = i_2 + k_1 * i_4
               j_1 = j_2 + k_1 * j_4
               j   = (i_1-1)*y_nnod+y_inod
               d_out(j) = a(j_1, i_1)
            end do ! k_1
         end if
      end if

      i_2 = eigen_loop_start(1, y_nnod,y_inod)
      i_3 = eigen_loop_end  (n, y_nnod,y_inod)
      do i_1=i_2,i_3
         j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
         j_2 = eigen_loop_start(j+1, x_nnod,x_inod)
         if ( j <= n ) then
            a(j_2:nm, i_1) = ZERO
         else
            a(1:nm, i_1) = ZERO
         end if
      end do! i_1

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      allocate(u0_z(nv*local_size+n_columns),
     &         v0_z(nv*local_size+n_columns))
      allocate(u1_z(nv*local_size+n_columns),
     &         v1_z(nv*local_size+n_columns))
      call CSTAB_adjust_base(u0_z(1), u_t(1), offset1)
      call CSTAB_adjust_base(v0_z(1), v_t(1), offset2)
      call CSTAB_adjust_base(u1_z(1), u_t(1), offset3)
      call CSTAB_adjust_base(v1_z(1), v_t(1), offset4)
      k_1 =  (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &           +(L1_LSIZE)
     &           +(L2_LSIZE/8)
      offset1 = offset1 + k_1 * 1
      offset2 = offset2 + k_1 * 2
      offset3 = offset3 + k_1 * 3
      offset4 = offset4 + k_1 * 4
      call CSTAB_round_offset(offset1)
      call CSTAB_round_offset(offset2)
      call CSTAB_round_offset(offset3)
      call CSTAB_round_offset(offset4)

      allocate ( tsave(0:local_size-1, 0:MBAND) )

      return
      end subroutine  eigen_prd_init


      subroutine eigen_prd_final(a, nm, n, d_out, e_out, ne, u_t)
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      use eigen_house
      implicit NONE

      integer, intent(in)    :: nm, ne, n
      real(8), intent(inout) :: a(1:nm, *)
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne, *)
      real(8), intent(out)   :: u_t(*)

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  i, j, k, L
      real(8)                ::  t


      if ( n >= 2 ) then
         do i = 2+MOD(n,2), 2, -1
            L = i-1
            y_owner_nod = eigen_owner_node (i, y_nnod,y_inod)
            x_owner_nod = eigen_owner_node (L, x_nnod,x_inod)
            if ( y_owner_nod == y_inod
     &         .AND. x_owner_nod == x_inod ) then
               i_1 = eigen_translate_g2l(i, y_nnod,y_inod)
               j_1 = eigen_translate_g2l(L, x_nnod,x_inod)
               e_out(i,1)  =  a(j_1,i_1)
               a(j_1,i_1)  =  ZERO
            end if
            L = i-2; if ( L < 1 ) cycle
            y_owner_nod = eigen_owner_node (i, y_nnod,y_inod)
            x_owner_nod = eigen_owner_node (L, x_nnod,x_inod)
            if ( y_owner_nod == y_inod
     &         .AND. x_owner_nod == x_inod ) then
               i_1 = eigen_translate_g2l(i, y_nnod,y_inod)
               j_1 = eigen_translate_g2l(L, x_nnod,x_inod)
               e_out(i,2)  =  a(j_1,i_1)
               a(j_1,i_1)  =  ZERO
            end if
         end do
      end if

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = eigen_loop_end  (n, y_nnod,y_inod)
         i_4 = x_nnod/n_common
         do i_1=i_2,i_3,i_4
            j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
            j_1 = eigen_translate_g2l(j, x_nnod,x_inod)
            t          = d_out(j)
            d_out(j)   = a(j_1,i_1)
            a(j_1,i_1) = t
         end do! i_1
      end if


      call reduce_dbl(d_out(1),   u_t(1), n, 1, x_COMM_WORLD)
      call reduce_dbl(d_out(1),   u_t(1), n, 1, y_COMM_WORLD)
      if ( MOD(n,2) /= 0 ) then
         e_out(1,1) = e_out(3,2)
      end if
      call reduce_dbl(e_out(1,1), u_t(1), n, 1, x_COMM_WORLD)
      call reduce_dbl(e_out(1,1), u_t(1), n, 1, y_COMM_WORLD)
      if ( MOD(n,2) /= 0 ) then
         e_out(3,2) = e_out(1,1)
         e_out(1,1) = ZERO
      end if


      deallocate(u0_z, v0_z)
      deallocate(u1_z, v1_z)

      deallocate ( tsave )


      return
      end subroutine  eigen_prd_final

