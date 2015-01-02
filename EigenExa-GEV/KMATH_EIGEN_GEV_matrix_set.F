      subroutine matrix_set_a(n, a)
      implicit double precision(a-h,o-z),integer(i-n)
      real(8)   :: a(n,n)

      one=1.0D0
      do j = 1,n
         do i = 1, n
            if( i == j ) then
!               a(i,j) = ( dble( n-i+1 ) ) / dble( n )+one /
!     &                       ( dble( i+j )-one)
                a(i,j) = 2d0
            else
!               a(i,j) = one / ( dble( i+j )-one )
                a(i,j) = 1d0/ABS(DBLE(i-j))**2
            end if 
         end do
      end do

      return
      end

      subroutine matrix_set_b(n, a)
      implicit double precision(a-h,o-z),integer(i-n)
      real(8)   :: a(n,n)

      do j = 1,n
         do i = 1, n
!            a(i,j) = dble( n-abs( i-j )) + n
            if ( i == j ) then
               a(i,j) = 1d2
            else
               a(i,j) = 1d0/ABS(DBLE(i-j))
            end if
         end do
      end do

      return
      end
