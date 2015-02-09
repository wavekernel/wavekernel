       subroutine matrix_adjust_s( n, a, b, nm )
       implicit NONE
       integer, intent(IN)    :: n, nm
       real(8), intent(IN)    :: a(1:n, *)
       real(8), intent(OUT)   :: b(1:nm, *)


          call matrix_adjust( n, a, n, b, nm )


       return
       end subroutine matrix_adjust_s
!---
       subroutine matrix_adjust( n, a, nma, b, nmb )
       use eigen_libs
       implicit NONE
!--
       integer, intent(IN)    :: n, nma, nmb
       real(8), intent(IN)    :: a(1:nma, *)
       real(8), intent(OUT)   :: b(1:nmb, *)

       integer                :: nnod, x_nnod, y_nnod
       integer                :: inod, x_inod, y_inod

       integer                :: i, i_1, i_2, i_3
       integer                :: j, j_1, j_2, j_3


          call eigen_get_procs( nnod, x_nnod, y_nnod )
          call eigen_get_id   ( inod, x_inod, y_inod )

          ! x_nnod = The number of rows
          ! y_nnod = The number of columns

          j_2 = eigen_loop_start( 1, x_nnod, x_inod )
          j_3 = eigen_loop_end  ( n, x_nnod, x_inod )
          i_2 = eigen_loop_start( 1, y_nnod, y_inod )
          i_3 = eigen_loop_end  ( n, y_nnod, y_inod )

          do i_1 = i_2, i_3
             i = eigen_translate_l2g( i_1, y_nnod, y_inod )
             do j_1 = j_2, j_3
                j = eigen_translate_l2g( j_1, x_nnod, x_inod )
                b(j_1, i_1) = a(j, i)
             end do
          end do


       return
       end subroutine matrix_adjust

