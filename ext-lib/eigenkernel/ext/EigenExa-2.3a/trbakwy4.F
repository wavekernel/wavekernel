       subroutine eigen_trbakwyx(n, nvec, a, nma0, z, nmz0, e, m0, iblk)
!$     use OMP_LIB
       use MPI
       use eigen_devel
       use eigen_libs
       use eigen_trbak
       implicit NONE

       integer, intent(in)    :: n, nvec, nma0, nmz0, m0, iblk
       real(8), intent(in)    :: a(*)
       real(8), intent(out)   :: z(*)
       real(8), intent(in)    :: e(*)

       real(8) , pointer      :: w(:), d(:)
       real(8) , pointer      :: v1(:), ss1(:)
       real(8) , pointer      :: v2(:), ss2(:)
       real(8) , pointer      :: v3(:), ss3(:)
       real(8) , pointer      :: tt(:)

       integer                :: nma, nmz, m
       integer                :: na, nm, kx
       integer                :: i
       integer                :: i_v1, i_v2, i_v3
       integer                :: i_s1, i_s2, i_s3
       integer                :: i_t
       integer                :: i_stat, ierr

       real(8) :: d1, d2
       real(8) :: s, dd(10)

       include 'CSTAB.h'
*
       nma  = nma0
       nmz  = nmz0
       m    = MIN(nsm, m0)
       if ( m < 1 ) m = 1
*
       call eigen_timer_reset( 2, 4, 0, 1 )
*
       na   = (n-1)/y_nnod+1
       na   = na  +MOD(na  -1,2)
       call CSTAB_get_optdim(nma,9,16*4,16*6,nm)
*
       allocate(
     &          d(1:n),
     &          v1(1:MAX(nm*m,n)+n_columns),
     &          ss1(1:na*m+ns0+n_columns),
     &          v2(1:MAX(nm*m,n)+n_columns),
     &          ss2(1:na*m+ns0+n_columns),
     &          v3(1:MAX(nm*m,n)+n_columns),
     &          ss3(1:na*m+ns0+n_columns),
     &          tt(1:na*m+ns0+n_columns),
     &          stat=i_stat)
       if ( i_stat /= 0 ) then
         if ( TRD_inod == 1 ) print*, "Memory allocation error."
         call MPI_Abort( TRD_COMM_WORLD, MPI_ERR_OTHER, ierr)
       end if
*
       d = ZERO
       v1 = ZERO
       ss1 = ZERO
       v2 = ZERO
       ss2 = ZERO
       v3 = ZERO
       ss3 = ZERO
       tt = ZERO
*
       call CSTAB_adjust_base(v1(1), z(1), i_v1)
       call CSTAB_adjust_base(ss1(1), z(1), i_s1)
       call CSTAB_adjust_base(v2(1), z(1), i_v2)
       call CSTAB_adjust_base(ss2(1), z(1), i_s2)
       call CSTAB_adjust_base(v3(1), z(1), i_v3)
       call CSTAB_adjust_base(ss3(1), z(1), i_s3)
       call CSTAB_adjust_base(tt(1), z(1), i_t)
          kx = (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &           +(L1_LSIZE)
     &           +(L2_LSIZE/8)
       i_v1 = i_v1 + kx*3
       i_s1 = i_s1 + kx*1
       i_v2 = i_v2 + kx*5
       i_s2 = i_s2 + kx*1
       i_v3 = i_v3 + kx*7
       i_s3 = i_s3 + kx*1
       i_t = i_t + kx*1
       call CSTAB_round_offset(i_v1)
       call CSTAB_round_offset(i_s1)
       call CSTAB_round_offset(i_v2)
       call CSTAB_round_offset(i_s2)
       call CSTAB_round_offset(i_v3)
       call CSTAB_round_offset(i_s3)
       call CSTAB_round_offset(i_t)
*
       allocate( TRBK_buf(1:m*nm) )
*
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
*-
       TRBK_TIME_REDUC_OVERHEAD_x = Reduce_Overhead_x
*-
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP PARALLEL
       call trbakwy_body( n, nvec,
     $      a(1), nma,
     $      z(1), nmz,
     $      e(1), d(1),
     $      v1(1+i_v1), v2(1+i_v2), v3(1+i_v3), nm, m,
     $      ss1(1+i_s1), ss2(1+i_s2), ss3(1+i_s3),
     $      tt(1+i_t), iblk, na )
!$OMP END PARALLEL
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
*
       deallocate( TRBK_buf )
*
       deallocate( d )
       deallocate( v1 )
       deallocate( ss1 )
       deallocate( v2 )
       deallocate( ss2 )
       deallocate( v3 )
       deallocate( ss3 )
       deallocate( tt )
*
       call eigen_timer_print( 'EigenExa (Back-transformation)' )
*
       return
       end
*
*
*
       subroutine trbakwy_body(n, nvec, a, nma, z, nmz, e, d,
     &             v1, v2, v3, nm, m, ss1, ss2, ss3, tt, iblk, nss)
!$     use OMP_LIB
       use MPI
       use eigen_devel
       use eigen_libs
       use eigen_TRBAK
       implicit NONE

       integer, intent(in)    :: n, nvec, nma, nmz
       integer, intent(in)    :: nm, m, iblk, nss
       real(8), intent(in)    :: a(1:nma,*)
       real(8), intent(inout) :: z(1:nmz,*)
       real(8), intent(in)    :: e(1:n)
       real(8), intent(inout) :: d(1:n)

       real(8), intent(out)   :: v1(1:nm,*), ss1(*)
       real(8), intent(out)   :: v2(1:nm,*), ss2(*)
       real(8), intent(out)   :: v3(1:nm,*), ss3(*)
       real(8), intent(out)   :: tt(*)
*
       integer                :: nodes(0:nsm-1)

       real(8), pointer       :: wk(:)
       integer                :: nx, lwk, ierr

       integer                :: i, L, mode
       integer                :: x_pos, x_root
       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       real(8)                :: s0, s1, s2, s3
       real(8)                :: d0, d1, d2, ds, de, dcom, dx, dy, dz

*
       dx = ZERO; dy = ZERO; dz = ZERO

!$OMP MASTER
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
#if TIMER_PRINT
       d1 = eigen_get_wtime()
#endif
       dcom = ZERO
*-
#if AT_BCAST_OVERLAP
!$     DO_OVERLAP_BCAST_LEVEL  = 2
!$     TRBK_TIME_COUNTER  = 0
!$     TRBK_TIME_INTERVAL = 1
!$     TRBK_TIME_NEXT     = 2
!$     TRBK_SWITCHED      = 0
!$     TRBK_TIME_TRBK1    = ZERO
!$     TRBK_TIME_TRBK1_   = ZERO
!$     TRBK_TIME_TRBK1x   = ZERO
!$     TRBK_TIME_TRBK1x_  = ZERO
!$     TRBK_TIME_TRBK2    = ZERO
!$     TRBK_TIME_BCAST    = ZERO
!$     TRBK_TIME_REDUC    = ONE
!$     TRBK_TIME_FR       = ZERO
!===
!$     call sync_other_than_master_init( TRBK_LOCK, TRBK_MASK )
!===
#else
!$     DO_OVERLAP_BCAST_LEVEL  = 0
#endif
*-
       lwk=MAX(((m-1)/y_nnod+1)*((n-1)/x_nnod+1),n)
       allocate(wk(lwk))
!$OMP END MASTER
!$OMP BARRIER
*
       i_2 = eigen_loop_start(1+iblk, y_nnod,y_inod)
       i_3 = eigen_loop_end  (n,      y_nnod,y_inod)
*
!$OMP MASTER
       d(1:n) = ZERO
       do i_1=i_2,i_3
!          i = eigen_translate_l2g(i_1, y_nnod,y_inod)
          i = (i_1-1)*y_nnod+y_inod
          l = i-iblk
!          x_root = eigen_owner_node(l, x_nnod,x_inod)
          x_root = MOD(l-1,x_nnod)+1
          if ( x_root == x_inod ) then
!             j_1 = eigen_translate_g2l(l, x_nnod,x_inod)
             j_1  = (l-1)/x_nnod+1
             d(i) = a(j_1,i_1)
          end if
       end do! i_1
       call reduce_dbl( d, wk, n, 1, TRD_COMM_WORLD )
       do i=1+iblk,n
          if ( e(i)*d(i) == ZERO ) then
             s0 = ZERO
          else
             s0 = (ONE/d(i))/e(i)
          end if
          d(i) = s0
       end do! i
!$OMP END MASTER
*
       i_2 = eigen_loop_start(1,    y_nnod,y_inod)
       i_3 = eigen_loop_end  (nvec, y_nnod,y_inod)
*
*----------------------------------------------------------------
*
       nx = MIN(MOD(n-(1+iblk)+1,m)+(1+iblk)-1,n)
*
       do i=(1+iblk),nx

!$OMP BARRIER
*
          if ( e(i) == ZERO ) cycle
*
          j_2 = eigen_loop_start(1,      x_nnod,x_inod)
          j_3 = eigen_loop_end  (i-iblk, x_nnod,x_inod)
          i_4=MOD(i_3-i_2+1,4)+i_2
*
!$OMP MASTER
          ds = eigen_get_wtime()
          nodes(0) = eigen_owner_node(i, y_nnod, y_inod)
          if ( nodes(0) == y_inod ) then
             i_1 = eigen_translate_g2l(i, y_nnod, y_inod)
             do j_1=j_2,j_3
                v1(j_1,1) = a(j_1,i_1)
             end do! j_1
          end if
          call bcast_dbl( v1(j_2,1), j_3-j_2+1,
     &                    nodes(0), 1, y_COMM_WORLD )

          de = eigen_get_wtime()
          dcom = dcom + (de - ds)
!$OMP END MASTER
*
!$OMP BARRIER
*
!$OMP MASTER
          call barrier( y_COMM_WORLD )
          call barrier( x_COMM_WORLD )
!$OMP END MASTER
!$    if ( omp_get_num_threads() == 1 .OR.
!$   +     omp_get_thread_num() == 1 ) then
          if ( i_4 == i_2 + 1 ) then
             i_1 = i_2
             s0 = ZERO
             do j_1=j_2,j_3
                s0=s0+v1(j_1,1)*z(j_1,i_1+0)
             end do! j_1
             ss1(i_1+0)=s0
          end if
          if ( i_4 == i_2 + 2 ) then
             i_1 = i_2
             s0 = ZERO
             s1 = ZERO
             do j_1=j_2,j_3
                s0=s0+v1(j_1,1)*z(j_1,i_1+0)
                s1=s1+v1(j_1,1)*z(j_1,i_1+1)
             end do! j_1
             ss1(i_1+0)=s0
             ss1(i_1+1)=s1
          end if
          if ( i_4 == i_2 + 3 ) then
             i_1 = i_2
             s0 = ZERO
             s1 = ZERO
             s2 = ZERO
             do j_1=j_2,j_3
                s0=s0+v1(j_1,1)*z(j_1,i_1+0)
                s1=s1+v1(j_1,1)*z(j_1,i_1+1)
                s2=s2+v1(j_1,1)*z(j_1,i_1+2)
             end do! j_1
             ss1(i_1+0)=s0
             ss1(i_1+1)=s1
             ss1(i_1+2)=s2
          end if
!$        end if
!$OMP DO SCHEDULE(DYNAMIC)
          do i_1=i_4,i_3,4
             s0 = ZERO
             s1 = ZERO
             s2 = ZERO
             s3 = ZERO
             do j_1=j_2,j_3
                s0=s0+v1(j_1,1)*z(j_1,i_1+0)
                s1=s1+v1(j_1,1)*z(j_1,i_1+1)
                s2=s2+v1(j_1,1)*z(j_1,i_1+2)
                s3=s3+v1(j_1,1)*z(j_1,i_1+3)
             end do! j_1
             ss1(i_1+0)=s0
             ss1(i_1+1)=s1
             ss1(i_1+2)=s2
             ss1(i_1+3)=s3
          end do! i_1
!$OMP END DO
*
!$OMP BARRIER
*
!$OMP MASTER
          ds = eigen_get_wtime()
          call reduce_dbl(ss1(i_2),tt, i_3-i_2+1, 2, x_COMM_WORLD)
          de = eigen_get_wtime()
          dcom = dcom + (de - ds)

          s0 = d(i)
          do i_1=i_2,i_3
             ss1(i_1) = ss1(i_1) * s0
          end do! i_1
!$OMP END MASTER

!$OMP BARRIER

          j_2 = eigen_loop_start(1,      x_nnod,x_inod)
          j_3 = eigen_loop_end  (i-iblk, x_nnod,x_inod)
          i_4=MOD(i_3-i_2+1,4)+i_2

!$OMP MASTER
          if ( i_4 == i_2 + 1 ) then
             i_1 = i_2
             s0 = ss1(i_1+0)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
             end do! j_1
          end if
          if ( i_4 == i_2 + 2 ) then
             i_1 = i_2
             s0 = ss1(i_1+0)
             s1 = ss1(i_1+1)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
             end do! j_1
          end if
          if ( i_4 == i_2 + 3 ) then
             i_1 = i_2
             s0 = ss1(i_1+0)
             s1 = ss1(i_1+1)
             s2 = ss1(i_1+2)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
                z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v1(j_1,1)
             end do! j_1
          end if
!$OMP END MASTER
!$OMP DO
          do i_1=i_4,i_3,4
             s0 = ss1(i_1+0)
             s1 = ss1(i_1+1)
             s2 = ss1(i_1+2)
             s3 = ss1(i_1+3)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
                z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v1(j_1,1)
                z(j_1,i_1+3) = z(j_1,i_1+3) + s3 * v1(j_1,1)
             end do! j_1
          end do! i_1
!$OMP ENDDO
*
!$OMP BARRIER
*
       end do
*
*----------------------------------------------------------------
*
       d2 = eigen_get_wtime()
       d0 = d2 -d1

#define	USE_BCASTW	1

       mode = 0; i = nx+1
*
!$OMP MASTER
          if ( nx+1 <= n ) then
             ds = eigen_get_wtime()
             call trbakwy_datacast(i_3, m, i, a(1,1), nma,
     $                          v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     $                          wk(1),
#endif
     $                          iblk)
             de = eigen_get_wtime()
             dcom = dcom + (de - ds)
          end if
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
          if ( nx+1+m <= n ) then
             ds = eigen_get_wtime()
             call trbakwy_datacast(i_3, m, i+m, a(1,1), nma,
     $                          v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     $                          wk(1),
#endif
     $                          iblk)
             de = eigen_get_wtime()
             dcom = dcom + (de - ds)
          end if
!$OMP END MASTER
*
       do i=nx+1, n, m
*
          if ( mode == 0 ) then
             call trbakwy_body2(i_3, z(1,1), nmz,
     $            d(1), e(1), v1(1,1), nm, m,
     $            i, ss1(1), tt, nss, iblk,
     $            dcom, dx, dy, dz)
          endif
          if ( mode == 1 ) then
             call trbakwy_body2(i_3, z(1,1), nmz,
     $            d(1), e(1), v2(1,1), nm, m,
     $            i, ss2(1), tt, nss, iblk,
     $            dcom, dx, dy, dz)
          endif
          if ( mode == 2 ) then
             call trbakwy_body2(i_3, z(1,1), nmz,
     $            d(1), e(1), v3(1,1), nm, m,
     $            i, ss3(1), tt, nss, iblk,
     $            dcom, dx, dy, dz)
          endif
*
!$OMP MASTER
          if ( i+2*m <= n ) then
          ds = eigen_get_wtime()
          if ( mode == 0 ) then
             call trbakwy_datacast(i_3, m, i+2*m, a(1,1), nma,
     $                             v3(1,1), nm, ss3(1),
#if !USE_BCASTW
     $                             wk(1),
#endif
     $                             iblk)
          end if
          if ( mode == 1 ) then
             call trbakwy_datacast(i_3, m, i+2*m, a(1,1), nma,
     $                             v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     $                             wk(1),
#endif
     $                             iblk)
          end if
          if ( mode == 2 ) then
             call trbakwy_datacast(i_3, m, i+2*m, a(1,1), nma,
     $                             v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     $                             wk(1),
#endif
     $                             iblk)
          end if
          de = eigen_get_wtime()
          dcom = dcom + (de - ds)
          end if
!$OMP END MASTER
*
          mode = MOD( mode + 1, 3 )
*
       end do
*
*----------------------------------------------------------------
*
!$     if ( omp_get_num_threads() > 1 ) then
!$OMP BARRIER
!$        if ( omp_get_thread_num() == 1 ) then
!$           ss1(1) = dx; ss1(2) = dy
!$        endif
!$OMP BARRIER
!$        if ( omp_get_thread_num() == 0 ) then
!$           dx = ss1(1); dy = ss1(2)
!$        endif
!$     endif
!$OMP BARRIER
*
#if AT_BCAST_OVERLAP
!===
!$OMP MASTER
!$     call omp_destroy_lock( TRBK_LOCK )
!$OMP END MASTER
!===
#endif
!$OMP MASTER
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
       deallocate(wk)

#if TIMER_PRINT>1
       d2 = eigen_get_wtime()

       if ( 1 == TRD_inod ) then
          print*,"TRBAK=",(d2-d1)
          print*,"COMM=",dcom
          print*,"   ",d2-d1,(2d0*nvec*n*n)/(d2-d1)*1d-9,"GFLOPS"
          if ( dx > 0 )
     &    print*,"   ",dx,(1d0*nvec*n*n)/(dx)*1d-9,"GFLOPS"
          if ( dy > 0 )
     &    print*,"   ",dy,(1d0*nvec*n*n)/(dy)*1d-9,"GFLOPS"
          print*,"   ",d0,dz
       end if
#endif
!$OMP END MASTER
*
*
       return

       end subroutine

       subroutine trbakwy_datacast(i_3, m, i, a, nma, v, nm,
     $                             ss,
#if !USE_BCASTW
     $                             wk,
#endif
     $                             iblk)
!$     use OMP_LIB
       use MPI
       use eigen_devel
       use eigen_libs
       use eigen_TRBAK
       implicit NONE

       integer, intent(in)   :: i_3, m, i, nma, nm, iblk
       real(8), intent(in)   :: a(nma,*)
       real(8), intent(out)  :: ss(*)
#if !USE_BCASTW
       real(8), intent(out)  :: wk(*)
#endif
       real(8), intent(out)  :: v(nm,*)

       integer :: j, iy, k0
       integer :: i_1
       integer :: j_1, j_2, j_3, j_4, j_5
       integer :: nodes(0:nsm-1)

!$     real(8) :: ds, de
       real(8) :: du1,du2
       integer :: req(2), IERR

*
!$      ds = omp_get_wtime()

        ss(1:i_3*m+ns0) = ZERO

#if USE_BCASTW
        j_3 = MIN(eigen_loop_end (i+m-1-iblk, x_nnod, x_inod), nm)
!-
      if ( MOD( m, y_nnod ) == 0 ) then
!-
        do j=0,m-1
           if ( y_inod == eigen_owner_node(i+j, y_nnod, y_inod) ) then
              iy  = eigen_owner_node(i+j, y_nnod, y_inod)
              i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
              k0  = (j/y_nnod)*j_3
              v(1+k0:j_3+k0, 1) = a(1:j_3, i_1)
           end if
        end do

        k0  = (m/y_nnod)*j_3
        call allgather_dbl(v, TRBK_buf, k0, 1, y_COMM_WORLD)

        do j=0,m-1
           iy = eigen_owner_node(i+j, y_nnod, y_inod)
           k0 = ((iy-1)*(m/y_nnod)+(j/y_nnod))*j_3
           v(1:j_3, j+1) = TRBK_buf(1+k0:j_3+k0)
           j_2 = eigen_loop_start(i+j,  x_nnod,x_inod)
           if ( j_2 <= j_3 ) v(j_2:j_3, j+1) = ZERO
        end do
!-
      else
!-
        do j=0,m-1
           if ( y_inod == eigen_owner_node(i+j, y_nnod, y_inod) ) then
              iy  = eigen_owner_node(i+j, y_nnod, y_inod)
              i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
              v(1:j_3, j+1) = a(1:j_3, i_1)
              j_2 = eigen_loop_start(i+j,  x_nnod,x_inod)
              if ( j_2 <= j_3 ) v(j_2:j_3, j+1) = ZERO
           end if
        end do
        do j=0,m-1,y_nnod
           iy  = eigen_owner_node(i+j, y_nnod, y_inod)
           k0  = MIN(m-j, y_nnod)
           call bcastw_dbl(v(1,j+1), j_3, iy,
     &                     nm, k0, TRBK_buf, 2, y_COMM_WORLD)
        end do
!-
      end if
#else
!-
        IF ( m > y_nnod .AND. y_nnod > 1 ) THEN

          do j=0,m-1
             nodes(j) = eigen_owner_node(i+j, y_nnod, y_inod)
          enddo

          do iy=1,y_nnod

             j_2 = eigen_loop_start(1,         x_nnod,x_inod)
             j_3 = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)

             k0=0
             do j=0,m-1
             if ( nodes(j) == iy ) then
                i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
                if ( y_inod == iy ) then
                   do j_1=j_2,j_3
                      wk(k0+j_1) = a(j_1, i_1)
                   end do! k
                endif
                k0=k0+(j_3-j_2+1)
             endif
             enddo

             call bcast_dbl(wk, k0, iy, 2, y_COMM_WORLD)

             k0=0
             do j=0,m-1
             if ( nodes(j) == iy ) then
                do j_1=j_2,j_3
                   v(j_1, j+1) = wk(k0+j_1)
                end do! k
                k0=k0+(j_3-j_2+1)
             endif
             enddo

             do j=0,m-1
             if ( nodes(j) == iy ) then
                j_4 = eigen_loop_start(i+j,       x_nnod,x_inod)
                j_5 = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
!OCL XFILL
                do j_1=j_4,j_5
                   v(j_1, j+1) = ZERO
                end do
             endif
             enddo

          enddo

        ELSE

          do j=0,m-1
             nodes(j) = eigen_owner_node(i+j, y_nnod, y_inod)
             if ( nodes(j) == y_inod ) then
                i_1 = eigen_translate_g2l(i+j, y_nnod, y_inod)
                j_2 = eigen_loop_start(1,         x_nnod,x_inod)
                j_3 = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
                do j_1=j_2,j_3
                   v(j_1, j+1) = a(j_1, i_1)
                end do! k
                j_2 = eigen_loop_start(i+j,       x_nnod,x_inod)
                j_3 = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
!CDIR NOVECTOR
                do j_1=j_2,j_3
                   v(j_1, j+1) = ZERO
                end do
             end if
          end do
          j_2 = eigen_loop_start(1,          x_nnod,x_inod)
          j_3 = eigen_loop_end  (i+m-1-iblk, x_nnod,x_inod)
          do j=0,m-1
             call bcast_dbl(v(1,j+1), j_3-j_2+1,
     &                      nodes(j), 2, y_COMM_WORLD)
          enddo

        ENDIF
!-
#endif

!$      de = omp_get_wtime()
!$      TRBK_TIME_BCAST = de - ds
*
       return
       end subroutine
