       subroutine ev_test(n, a, nma, w, z, nmz)
!$     use OMP_LIB
       use MPI
       use eigen_libs
       implicit NONE

       integer, intent(in)    :: n, nma, nmz
       real(8), intent(in)    :: a(1:nma,*)
       real(8), intent(in)    :: w(1:n)
       real(8), intent(in)    :: z(1:nmz,*)

       real(8), pointer       :: b(:)
       real(8), pointer       :: c(:,:), d(:,:), e(:,:)

       real(8)                :: m_epsilon

       integer                :: i, i_1, i_2, i_3, i_4
       integer                :: j, j_1, j_2, j_3, j_4
       integer                :: i0, k, ii_1, jj_1
       integer                :: ierr, k_4, i_m, i_mx, larray
       integer, parameter     :: BLK = 3
       real(8)                :: a0, a1, a2, a3
       real(8)                :: z0, z1, z2, z3
       real(8)                :: s0, s1, s2, s3
       real(8)                :: t0, t1, t2, t3
       real(8)                :: ss0, ss1, ss2, ss3
       real(8)                :: tt0, tt1, tt2, tt3
       real(8)                :: s_0, s_1, s_2, s_3
       real(8)                :: t_0, t_1, t_2, t_3
       real(8)                :: b0, bb0, b_0, ra, rr, rz, r
       real(8)                :: f, ff(1:BLK)

       integer                :: TRD_COMM_WORLD, x_comm, y_comm
       integer                :: nnod, x_nnod, y_nnod
       integer                :: inod, x_inod, y_inod
       integer                :: x_owner_nod, y_owner_nod
*
       integer                :: local_size
       real(8), pointer       :: b0_(:)
       real(8), pointer       :: rz_(:)
*
       local_size = 1
!$OMP PARALLEL
!$OMP MASTER
!$     local_size = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
*
       allocate( b0_(0:local_size) )
       allocate( rz_(0:local_size) )
*
       m_epsilon = machine_epsilon()
*
       rz=0D0; rz_=0D0
       rr=0D0
       r =0D0
               b0_=0D0

       call eigen_get_comm ( TRD_COMM_WORLD, x_comm, y_comm )
       call eigen_get_procs( nnod, x_nnod, y_nnod )
       call eigen_get_id   ( inod, x_inod, y_inod )

       j_2 = eigen_loop_start   ( 1, x_nnod, x_inod )
       j_3 = eigen_loop_end     ( n, x_nnod, x_inod )
       j_4 = eigen_translate_g2l( n, x_nnod, x_inod )
       i_2 = eigen_loop_start   ( 1, y_nnod, y_inod )
       i_3 = eigen_loop_end     ( n, y_nnod, y_inod )
       i_4 = eigen_translate_g2l( n, y_nnod, y_inod )
       k_4 = MAX(i_4, j_4)+8

       larray = n
       allocate(b(1:larray))

       larray = j_4
       allocate(c(1:larray,1:BLK))

       larray = MAX(n+1, 3*k_4, x_nnod, y_nnod)
       allocate(d(1:larray,1:BLK), e(1:larray,1:BLK))

       b(1:n)=0d0

       ra =0.0D+00
       if ( i_2 <= i_3 ) then
!$OMP PARALLEL DO PRIVATE(i_1,j_1,t0,a0)
          do i_1 = i_2, i_3
             t0 = 0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1)
                t0 = t0 + ABS(a0)           ! u = |A|_1
             end do
             d(i_1,1) = t0
          end do
!$OMP END PARALLEL DO
          call MPI_Allreduce( d(1,1), e(1,1), (i_3-i_2+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, x_comm, ierr)
!$OMP PARALLEL DO REDUCTION(MAX:ra)
          do i_1 = i_2, i_3
             ra = MAX(ra, e(i_1,1))
          end do
!$OMP END PARALLEL DO
       end if

       if ( i_2 <= i_3 ) then
          i_m = MOD(i_3-i_2+1,4)+i_2
       else
          i_m = i_3+1
       end if

       i_mx = MOD(n,BLK)+1

       do i=1,i_mx-1,1

          x_owner_nod = eigen_owner_node( i, x_nnod, x_inod )
          y_owner_nod = eigen_owner_node( i, y_nnod, y_inod )
          jj_1        = eigen_translate_g2l( i, x_nnod, x_inod )
          ii_1        = eigen_translate_g2l( i, y_nnod, y_inod )

          if ( y_owner_nod == y_inod ) then
             do j_1 = j_2, j_3
                c(j_1,1) = z(j_1, ii_1)
             end do
             ff(1) = 1.0D+00
          else
             ff(1) = 0.0D+00
          end if
          call MPI_Bcast( c(1,1), j_4, MPI_DOUBLE_PRECISION,
     &                   y_owner_nod-1, y_comm, ierr)

          if ( i_2 <= i_3 ) then
          do i_1 = i_2, i_m-1
             t0 = 0.0D+00
             s0 = 0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1)
                b0 = c(j_1,1)
                t0 = t0 + b0 * a0           ! t = Az(i)
                s0 = s0 + b0 * z(j_1, i_1)  ! s = Zz(i)
             end do! k
             d(1+(i_1-1)*2,1) = t0
             d(2+(i_1-1)*2,1) = s0
          end do

!$OMP PARALLEL DO PRIVATE(i_1,j_1,a0,a1,a2,a3,b0,
!$OMP+                    t0,s0,t1,s1,t2,s2,t3,s3)
          do i_1 = i_m, i_3, 4
             t0 = 0.0D+00
             t1 = 0.0D+00
             t2 = 0.0D+00
             t3 = 0.0D+00
             s0 = 0.0D+00
             s1 = 0.0D+00
             s2 = 0.0D+00
             s3 = 0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1+0)
                a1 = a(j_1,i_1+1)
                a2 = a(j_1,i_1+2)
                a3 = a(j_1,i_1+3)
                b0 = c(j_1,1)
                t0 = t0 + b0 * a0           ! t = Az(i)
                t1 = t1 + b0 * a1           ! t = Az(i)
                t2 = t2 + b0 * a2           ! t = Az(i)
                t3 = t3 + b0 * a3           ! t = Az(i)
                s0 = s0 + b0 * z(j_1, i_1+0)  ! s = Zz(i)
                s1 = s1 + b0 * z(j_1, i_1+1)  ! s = Zz(i)
                s2 = s2 + b0 * z(j_1, i_1+2)  ! s = Zz(i)
                s3 = s3 + b0 * z(j_1, i_1+3)  ! s = Zz(i)
             end do! k
             d(1+(i_1-1)*2,1) = t0
             d(2+(i_1-1)*2,1) = s0
             d(1+(i_1  )*2,1) = t1
             d(2+(i_1  )*2,1) = s1
             d(1+(i_1+1)*2,1) = t2
             d(2+(i_1+1)*2,1) = s2
             d(1+(i_1+2)*2,1) = t3
             d(2+(i_1+2)*2,1) = s3
          end do
!$OMP END PARALLEL DO
          else

             d(1:2*j_4,1) = 0D0
             e(1:2*j_4,1) = 0D0

          end if

          if ( i_2 <= i_3 ) then
          call MPI_Allreduce( d(1,1), e(1,1), 2*(i_3-i_2+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, x_comm, ierr)
          end if

          call datacast_dbl(d(1,1), c(1,1),
     &                      d(1+k_4,1), d(1+2*k_4,1), j_4)

          x_owner_nod = eigen_owner_node( i, x_nnod, x_inod )
          ii_1        = eigen_translate_g2l( i, y_nnod, y_inod )
          if ( x_owner_nod == x_inod ) then
             b0_ =0.0D+00
             j = 0

          if ( i_2 <= i_3 ) then
             f = ff(1)
!$OMP PARALLEL DO PRIVATE(i_1,j,s0,t0)
             do i_1 = i_2, i_3
!$              j = omp_get_thread_num()
                t0 = e(1+(i_1-1)*2,1)
                s0 = e(2+(i_1-1)*2,1)
                t0 = t0 - w(i) * d(i_1,1)     ! t = Az(i) - wz(i)
                if ( i_1 == ii_1 ) s0 = s0 - f ! s = Zz(i) - 1(i)
!                s0 = s0 - ff(1)*((n-IABS(i_1-ii_1))/n)
                b0_(j) = b0_(j) + ABS(t0) ! b(i)=|Az(i)-wz(i)|_1
                rz_(j) = rz_(j) + s0**2   ! rz  =|ZZ-I|_F
             end do
!$OMP END PARALLEL DO
          end if
             b0 =0.0D+00
             do i_1=0,local_size-1
                b0 = b0 + b0_(i_1)
             end do
             b(i) = b0
          else
             b(i) =0.0D+00
          end if

       end do
*>
       do i=i_mx,n,BLK

       do i0=0,BLK-1
          x_owner_nod = eigen_owner_node( i+i0, x_nnod, x_inod )
          y_owner_nod = eigen_owner_node( i+i0, y_nnod, y_inod )
          jj_1        = eigen_translate_g2l( i+i0, x_nnod, x_inod )
          ii_1        = eigen_translate_g2l( i+i0, y_nnod, y_inod )
          if ( y_owner_nod == y_inod ) then
             do j_1 = j_2, j_3
                c(j_1,1+i0) = z(j_1, ii_1)
             end do
             ff(1+i0) = 1D0
          else
             ff(1+i0) =0.0D+00
          end if
          call MPI_Bcast( c(1,1+i0), j_4, MPI_DOUBLE_PRECISION,
     &                   y_owner_nod-1, y_comm, ierr)
       end do

          if ( i_2 <= i_3 ) then
          do i_1 = i_2, i_m-1
             t0 = 0.0D+00
             s0 = 0.0D+00
             tt0 = 0.0D+00
             ss0 = 0.0D+00
             t_0 = 0.0D+00
             s_0 = 0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1)
                z0 = z(j_1,i_1)
                b0  = c(j_1,1)
                bb0 = c(j_1,2)
                b_0 = c(j_1,3)
                t0  = t0  + b0  * a0  ! t = Az(i)
                s0  = s0  + b0  * z0  ! s = Zz(i)
                tt0 = tt0 + bb0 * a0  ! t = Az(i)
                ss0 = ss0 + bb0 * z0  ! s = Zz(i)
                t_0 = t_0 + b_0 * a0  ! t = Az(i)
                s_0 = s_0 + b_0 * z0  ! s = Zz(i)
             end do! k
             d(1+(i_1-1)*2,1) = t0
             d(2+(i_1-1)*2,1) = s0
             d(1+(i_1-1)*2,2) = tt0
             d(2+(i_1-1)*2,2) = ss0
             d(1+(i_1-1)*2,3) = t_0
             d(2+(i_1-1)*2,3) = s_0
          end do

!$OMP PARALLEL DO PRIVATE(i_1,j_1,jj_1,a0,a1,a2,a3,
!$OMP+                    z0,z1,z2,z3,b0,bb0,b_0,
!$OMP+                    t0,s0,t1,s1,t2,s2,t3,s3,
!$OMP+                    tt0,ss0,tt1,ss1,tt2,ss2,tt3,ss3,
!$OMP+                    t_0,s_0,t_1,s_1,t_2,s_2,t_3,s_3)
          do i_1 = i_m, i_3, 4
             t0 =0.0D+00
             s0 =0.0D+00
             t1 =0.0D+00
             s1 =0.0D+00
             t2 =0.0D+00
             s2 =0.0D+00
             t3 =0.0D+00
             s3 =0.0D+00
             tt0 =0.0D+00
             ss0 =0.0D+00
             tt1 =0.0D+00
             ss1 =0.0D+00
             tt2 =0.0D+00
             ss2 =0.0D+00
             tt3 =0.0D+00
             ss3 =0.0D+00
             t_0 =0.0D+00
             s_0 =0.0D+00
             t_1 =0.0D+00
             s_1 =0.0D+00
             t_2 =0.0D+00
             s_2 =0.0D+00
             t_3 =0.0D+00
             s_3 =0.0D+00
             do j_1 = j_2, j_3
                a0 = a(j_1,i_1+0)
                a1 = a(j_1,i_1+1)
                a2 = a(j_1,i_1+2)
                a3 = a(j_1,i_1+3)
                z0 = z(j_1,i_1+0)
                z1 = z(j_1,i_1+1)
                z2 = z(j_1,i_1+2)
                z3 = z(j_1,i_1+3)
                b0  = c(j_1,1)
                bb0 = c(j_1,2)
                b_0 = c(j_1,3)
                t0  = t0  + b0  * a0  ! t = Az(i)
                tt0 = tt0 + bb0 * a0  ! t = Az(i)
                t_0 = t_0 + b_0 * a0  ! t = Az(i)
                t1  = t1  + b0  * a1  ! t = Az(i)
                tt1 = tt1 + bb0 * a1  ! t = Az(i)
                t_1 = t_1 + b_0 * a1  ! t = Az(i)
                t2  = t2  + b0  * a2  ! t = Az(i)
                tt2 = tt2 + bb0 * a2  ! t = Az(i)
                t_2 = t_2 + b_0 * a2  ! t = Az(i)
                t3  = t3  + b0  * a3  ! t = Az(i)
                tt3 = tt3 + bb0 * a3  ! t = Az(i)
                t_3 = t_3 + b_0 * a3  ! t = Az(i)
                s0  = s0  + b0  * z0  ! s = Zz(i)
                ss0 = ss0 + bb0 * z0  ! s = Zz(i)
                s_0 = s_0 + b_0 * z0  ! s = Zz(i)
                s1  = s1  + b0  * z1  ! s = Zz(i)
                ss1 = ss1 + bb0 * z1  ! s = Zz(i)
                s_1 = s_1 + b_0 * z1  ! s = Zz(i)
                s2  = s2  + b0  * z2  ! s = Zz(i)
                ss2 = ss2 + bb0 * z2  ! s = Zz(i)
                s_2 = s_2 + b_0 * z2  ! s = Zz(i)
                s3  = s3  + b0  * z3  ! s = Zz(i)
                ss3 = ss3 + bb0 * z3  ! s = Zz(i)
                s_3 = s_3 + b_0 * z3  ! s = Zz(i)
             end do! k
             d(1+(i_1-1)*2,1) = t0
             d(2+(i_1-1)*2,1) = s0
             d(1+(i_1  )*2,1) = t1
             d(2+(i_1  )*2,1) = s1
             d(1+(i_1+1)*2,1) = t2
             d(2+(i_1+1)*2,1) = s2
             d(1+(i_1+2)*2,1) = t3
             d(2+(i_1+2)*2,1) = s3
             d(1+(i_1-1)*2,2) = tt0
             d(2+(i_1-1)*2,2) = ss0
             d(1+(i_1  )*2,2) = tt1
             d(2+(i_1  )*2,2) = ss1
             d(1+(i_1+1)*2,2) = tt2
             d(2+(i_1+1)*2,2) = ss2
             d(1+(i_1+2)*2,2) = tt3
             d(2+(i_1+2)*2,2) = ss3
             d(1+(i_1-1)*2,3) = t_0
             d(2+(i_1-1)*2,3) = s_0
             d(1+(i_1  )*2,3) = t_1
             d(2+(i_1  )*2,3) = s_1
             d(1+(i_1+1)*2,3) = t_2
             d(2+(i_1+1)*2,3) = s_2
             d(1+(i_1+2)*2,3) = t_3
             d(2+(i_1+2)*2,3) = s_3
          end do
!$OMP END PARALLEL DO
          else

             d(1:2*j_4,1) = 0D0
             e(1:2*j_4,1) = 0D0

          end if

          if ( i_2 <= i_3 ) then

       do i0=0,BLK-1
          call MPI_Allreduce( d(1,1+i0), e(1,1+i0), 2*(i_3-i_2+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, x_comm, ierr)
       end do

          end if

       do i0=0,BLK-1
          call datacast_dbl(d(1,1+i0), c(1,1+i0),
     &                      d(1+k_4,1+i0), d(1+2*k_4,1+i0), j_4)
       end do

       do i0=0,BLK-1
          x_owner_nod = eigen_owner_node( i+i0, x_nnod, x_inod )
          ii_1        = eigen_translate_g2l( i+i0, y_nnod, y_inod )
          if ( x_owner_nod == x_inod ) then
             b0_ =0.0D+00
             j = 0
          if ( i_2 <= i_3 ) then
             f = ff(1+i0)
!$OMP PARALLEL DO PRIVATE(i_1,j,s0,t0)
             do i_1 = i_2, i_3
!$              j = omp_get_thread_num()
                t0 = e(1+(i_1-1)*2,1+i0)
                s0 = e(2+(i_1-1)*2,1+i0)
                t0 = t0 - w(i+i0) * d(i_1,1+i0)     ! t = Az(i) - wz(i)
                if ( i_1 == ii_1 ) s0 = s0 - f ! s = Zz(i) - 1(i)
!                s0 = s0 - ff(1+i0)*((n-IABS(i_1-ii_1))/n)
                b0_(j) = b0_(j) + ABS(t0) ! b(i)=|Az(i)-wz(i)|_1
                rz_(j) = rz_(j) + s0**2   ! rz  =|ZZ-I|_F
             end do
!$OMP END PARALLEL DO
          end if
             b0 =0.0D+00
             do i_1=0,local_size-1
                b0 = b0 + b0_(i_1)
             end do
             b(i+i0) = b0
          else
             b(i+i0) =0.0D+00
          end if
       end do

       end do
*>
       rz = 0.0D0
       do i_1=0,local_size-1
          rz = rz + rz_(i_1)
       end do
*>
       d(1:n,1)=b(1:n)
       d(n+1,1)=rz
       call MPI_Allreduce( d(1,1), e(1,1), n+1, MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, TRD_COMM_WORLD, ierr )
       b(1:n)=e(1:n,1)
       rz    =e(n+1,1)

       r=ra
       call MPI_Allreduce( ra, r, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, TRD_COMM_WORLD, ierr )
       ra=r

       rr=0d0
       do i=1,n
          if ( b(i) > rr ) then
             rr=b(i); k=i
          end if
       end do! i

       if ( x_inod == 1 .AND. y_inod == 1 ) then
!          print*, "max|Ax-wx|=", rr, k
          print*, "|A|_{1}=", ra
          print*, "epsilon=", m_epsilon
          rr =  rr/(n*m_epsilon*ra)
          print*, "max|Ax-wx|_{1}/Ne|A|_{1}=", rr, k
          if ( rr < 2.5D+00 ) then
             print*,"*** Residual Error Test ***   : PASSED"
          else
             print*,"***=============================******"
             print*,"*** Residual Error Test ***   : FAILED"
             print*,"***=============================******"
          end if
          print*, "|ZZ-I|_{F}=", SQRT(rz)
          if ( rz < 1.0D-20 ) then
             print*,"*** Orthogonality  Test ***   : PASSED"
          else
             print*,"***=============================******"
             print*,"*** Orthogonality  Test ***   : FAILED"
             print*,"***=============================******"
          end if
       endif
*
       deallocate( b  )
       deallocate( c )
       deallocate( d )
       deallocate( e )
*
       deallocate( b0_ )
       deallocate( rz_ )
*
       return
       end

