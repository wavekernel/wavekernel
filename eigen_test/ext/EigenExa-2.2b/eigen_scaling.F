      subroutine eigen_scaling( n, a, lda, x )
!$    use OMP_LIB
      use MPI
      use eigen_libs
      use eigen_devel
      implicit NONE

      integer, intent(IN)    :: lda, n
      real(8), intent(INOUT) :: a(lda, *)
      real(8), intent(OUT)   :: x

      real(8)                :: SAFMIN, EPS, SMLNUM, BIGNUM
      real(8)                :: RMIN, RMAX, ANRM, SIGMA
      real(8)                :: t
      real(8), external      :: DLAMCH

      integer                :: ERR_MPI
      integer                :: i, i_1, i_2, i_3
      integer                :: j, j_1, j_2, j_3


         SAFMIN = DLAMCH( 'Safe minimum' )
         EPS    = DLAMCH( 'Precision' )
         SMLNUM = SAFMIN / EPS
         BIGNUM = ONE / SMLNUM
         RMIN = SQRT( SMLNUM )
         RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*-
         i_2 = eigen_loop_start   ( 1, y_nnod, y_inod )
         i_3 = eigen_loop_end     ( n, y_nnod, y_inod )

         ANRM = ZERO
         do i_1=i_2,i_3

            i   = eigen_translate_l2g( i_1, y_nnod, y_inod )
            j   = i

            j_2 = eigen_loop_start   ( 1, x_nnod, x_inod )
            j_3 = eigen_loop_end     ( j, x_nnod, x_inod )

            do j_1=j_2,j_3
               ANRM = MAX( ANRM, ABS( a(j_1, i_1) ) )
            end do

         end do

         t = ANRM
         call MPI_Allreduce( t, ANRM,
     &                       1, MPI_DOUBLE_PRECISION,
     &                       MPI_MAX, TRD_COMM_WORLD, ERR_MPI )
*-
!
!   SIGMA is calcurated in order to avoid Under/over-flow
!
         SIGMA = ONE
         IF( ANRM /= ZERO .AND. ANRM < RMIN ) THEN
            SIGMA = RMIN / ANRM
         ELSE IF( ANRM > RMAX ) THEN
            SIGMA = RMAX / ANRM
         END IF
*-
         x = SIGMA
         if ( x == ONE ) return
*-
         i_2 = eigen_loop_start   ( 1, y_nnod, y_inod )
         i_3 = eigen_loop_end     ( n, y_nnod, y_inod )

         do i_1=i_2,i_3

            i   = eigen_translate_l2g( i_1, y_nnod, y_inod )
            j   = i

            j_2 = eigen_loop_start   ( 1, x_nnod, x_inod )
            j_3 = eigen_loop_end     ( j, x_nnod, x_inod )

            do j_1=j_2,j_3
               a(j_1, i_1) = a(j_1, i_1) * x
            end do

         end do


      return
      end subroutine  eigen_scaling

