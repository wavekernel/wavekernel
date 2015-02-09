c======================================================================c
      subroutine eigen_prd_2update(
     &      ar, nma,
     &      ur, uyr, vr, vyr, nmv,
     &      m_size, i_base)
c======================================================================c
!$    use OMP_LIB
      use eigen_devel
      use eigen_libs
      implicit NONE
c
      real(8), intent(inout) :: ar (nma, *)
      real(8), intent(in)    :: ur (nmv, *)
      real(8), intent(in)    :: uyr(nmv, *)
      real(8), intent(in)    :: vr (nmv, *)
      real(8), intent(in)    :: vyr(nmv, *)
      integer, intent(in)    :: nma, nmv, m_size, i_base
c
      integer ::  i, j, k, k0, k1, k2, L, m, n

      integer ::  i_1,i_2,i_3,i_4
      integer ::  j_1,j_2,j_3,j_4

      integer ::  ii_1,ii_2,ii_3,ii_4
      integer ::  jj_1,jj_2,jj_3,jj_4

      integer ::  blk_size1, blk_size2
      integer ::  ii_step
c
      integer ::  local_rank, local_size
c
      integer, parameter :: BLAS_CHUNK = 7*8
      integer, parameter :: BLAS_NVECT = 256*3
c
      intrinsic :: MIN, MAX
      external  :: DGEMM
c
c======================================================================c
c
      if ( i_base <= 0 ) return
c
c======================================================================c
c
      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()
c
      m  = m_size                               ! preserve the argument variable 
      n  = eigen_loop_end  (i_base, y_nnod, y_inod) ! local matrix size
      k0 = eigen_translate_l2g(n, y_nnod, y_inod)  ! translation to the global index

      jj_2 = 1                                  ! beggining of loop
      jj_3 = eigen_loop_end  (k0, x_nnod, x_inod)   ! end of loop

      ii_step = 0
      do j_1 = jj_2, jj_3, BLAS_NVECT
         j_4 = MIN(j_1+BLAS_NVECT-1, jj_3)      ! [j_1:j_4] available on this iteration 

         k1   = eigen_translate_l2g(j_1, x_nnod, x_inod) ! translation to the global index
         ii_2 = eigen_loop_start(k1,     y_nnod, y_inod)  ! beggining of loop
         ii_2 = MAX(1, ii_2)                          !  ** should be .GE. 1
         ii_3 = n                                     ! end of loop

         do i_1 = ii_2, ii_3, BLAS_CHUNK
            i_4 = MIN(i_1+BLAS_CHUNK-1, ii_3)   ! [i_1:i_4] available on this iteration

            k2  = eigen_translate_l2g(i_4, y_nnod, y_inod) ! translation to the global index
            j_3 = eigen_loop_end  (k2, x_nnod, x_inod) ! end of loop
            j_3 = MIN(j_4, j_3)

            i_2 = i_1; i_3 = i_4
            j_2 = j_1

            blk_size1 = j_3-j_2+1
            blk_size2 = i_3-i_2+1

            if ( blk_size1 > 0 .AND. blk_size2 > 0 ) then

               if ( MOD(ii_step, local_size) == local_rank ) then

                  call DGEMM('N','T',
     &                   blk_size1, blk_size2, m,
     &                   MONE, ur (j_1, 1),   nmv,
     &                         vyr(i_1, 1),   nmv,
     &                   ONE,  ar (j_1, i_1), nma)
                  call DGEMM('N','T',
     &                   blk_size1, blk_size2, m,
     &                   MONE, vr (j_1, 1),   nmv,
     &                         uyr(i_1, 1),   nmv,
     &                   ONE,  ar (j_1, i_1), nma)

               end if

               ii_step = ii_step+1

            end if

         end do ! i_1

      end do ! j_1
c-
      return
      end subroutine eigen_prd_2update
c======================================================================c

