       module eigen_libs
!--
          use eigen_blacs, only :  eigen_get_blacs_context
          implicit NONE

          character(32), private ::  CodeNAME = 'EigenExa'
          type version_t
             integer ::  Major_Version  ! Major version
             integer ::  Minor_Version  ! Minor version
             integer ::  Patch_Level    ! Patchlevel 0=none,1=a,2=b,...
             character(32) :: date      ! Release date
             character(32) :: vcode     ! Version code name
          end type  version_t
          type(version_t), private  :: Eigen_Version
     &       = version_t ( 2, 3, 1, 'April 14, 2015', 'Wakamurasaki' )
!--
          integer, parameter ::  eigen_NB    = 64
!--
       interface
!--
          subroutine eigen_sx( n, nvec, a, lda, w, z, ldz, times,
     &                         m_forward, m_backward, mode )
          integer, intent(in)    :: n, nvec, lda, ldz
            ! lda and ldz are the leading dimensions corresponding to
            !  a and z, respectively. These values must be equal or larger
            !  than one obtained by calling eigen_get_matdims
          integer, intent(in), optional :: m_forward, m_backward
            ! default values are
            !   m_forward = 48, and m_backward = 128.
          character, intent(in), optional :: mode
            ! = 'A' : eigenvalues and corresponding eigenvectors (default)
            ! = 'N' : eigenvalues only
            ! = 'X' : mode 'A' + improve accuracy of eigenvalues
          real(8), intent(inout) :: a(1:lda, *)
          real(8), intent(out)   :: w(*), z(1:ldz, *)
          real(8), intent(out)   :: times(8)
          end subroutine  eigen_sx
!--
          subroutine eigen_s( n, nvec, a, lda, w, z, ldz, times,
     &                         m_forward, m_backward, mode )
          integer, intent(in)    :: n, nvec, lda, ldz
            ! lda and ldz are the leading dimensions corresponding to
            !  a and z, respectively. These values must be equal or larger
            !  than one obtained by calling eigen_get_matdims
          integer, intent(in), optional :: m_forward, m_backward
            ! default values are
            !   m_forward = 48, and m_backward = 128.
          character, intent(in), optional :: mode
            ! = 'A' : eigenvalues and corresponding eigenvectors (default)
            ! = 'N' : eigenvalues only
            ! = 'X' : mode 'A' + improve accuracy of eigenvalues
          real(8), intent(inout) :: a(1:lda, *)
          real(8), intent(out)   :: w(*), z(1:ldz, *)
          real(8), intent(out)   :: times(8)
          end subroutine  eigen_s
!--
       end interface

!          private ::  eigen_loop_calc

       contains
!--
!--
          subroutine eigen_get_version( version, date, vcode )
          use eigen_devel
          implicit NONE
          integer, intent(out) :: version
          character*(*), intent(out), optional :: date
          character*(*), intent(out), optional :: vcode

             version = Eigen_Version%Major_Version * 100
     &               + Eigen_Version%Minor_Version * 10
     &               + Eigen_Version%Patch_Level

             if ( PRESENT(date) ) then
                date = Eigen_Version%date
             end if

             if ( PRESENT(vcode) ) then
                vcode = Eigen_Version%vcode
             end if

          return
          end subroutine eigen_get_version
!--
!--
          subroutine eigen_show_version( )
          use eigen_devel
          implicit NONE
          character*(256) :: version
          character*(1  ) :: patchlevel

             patchlevel = ' '
             select case ( Eigen_Version%Patch_Level )
                case (0); patchlevel = ' '
                case (1); patchlevel = 'a'
                case (2); patchlevel = 'b'
                case (3); patchlevel = 'c'
                case (4); patchlevel = 'd'
                case (5); patchlevel = 'e'
                case (6); patchlevel = 'f'
                case (7); patchlevel = 'g'
                case (8); patchlevel = 'h'
                case (9); patchlevel = 'i'
                case default ; patchlevel = '*'
             end select

             write( version, '(I1,A,I1,A)' )
     &              Eigen_Version%Major_Version,
     &                 '.',Eigen_Version%Minor_Version, patchlevel

             if ( TRD_inod == 1 ) then
                print*, "## EigenExa version (", trim(version),
     &                  ") / (", trim(Eigen_Version%date),
     &                  ") / (", trim(Eigen_Version%vcode), ")"
             end if

          return
          end subroutine eigen_show_version
!--
!--
          subroutine eigen_initialized( flag )
          use eigen_devel
          implicit NONE
*-
          logical, intent(out) :: flag

             call eigen_get_initialized( flag )

          return
          end subroutine eigen_initialized
!--
!--
          subroutine eigen_init( comm, order,
     &                           scalapack_context, gridmap )
*-
!$        use OMP_LIB
          use MPI
          use eigen_devel
          use eigen_blacs, only : eigen_blacs_init
          implicit NONE
*-
          integer, intent(in), optional ::  comm
          character*(*), intent(in), optional ::  order
          ! scalapack_context option will be supported in future version
          integer, intent(in), optional ::  scalapack_context
          ! gridmap option will be supported in future version
          integer, intent(in), optional ::  gridmap(:,:)
*-
          logical :: flag
          integer ::  n1, n2, n3, i, j, k, ierr
          integer :: topo_type, cart_dim, dims(2), coords(2)
          logical :: periods(2)
          integer :: old_grp, new_grp

          integer, pointer :: tmpgrid(:,:)
          integer          :: comm0
          integer          :: group0, group1, kk0(1), kk1(1)

          integer          :: local_size
!$        integer          :: th0(2), th1(2)

          real(8)          :: s1, s2, d1, d2, buff(4)
          real(8), pointer :: buff1(:), buff2(:)

          character*1 :: GRID_major = 'C'


             call eigen_timer_reset( 0, 0, 0, 0 )

             call eigen_initialized( flag )
          if ( flag ) then
          if ( TRD_inod == 1 ) then
             print*,"*************"
             print*,"** CAUTION **"
             print*,"*************"
             print*,"You are going to initialize EigenExa,"
             print*,"while EigenExa was not freed at last call."
             print*,"EigenExa restarts again by itself."
          end if
             call eigen_free( )
          end if

!---- Comm setup ---

          if ( PRESENT(comm) ) then
             comm0 = comm
          else
             comm0 = MPI_COMM_WORLD
          end if
          call MPI_Comm_Test_Inter( comm0, flag, ierr )
          if ( flag .OR. ierr /= MPI_SUCCESS ) then
             print*,"*************"
             print*,"** CAUTION **"
             print*,"*************"
             print*,"You are going to initialize EigenExa with"
             print*,"an invalid communicator."
             print*,"EigenExa terminates this run."
             call sleep(1)
             call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_OTHER, ierr )
          endif
          if ( comm0 == MPI_COMM_NULL ) then
             TRD_COMM_WORLD = MPI_COMM_NULL
             x_COMM_WORLD   = MPI_COMM_NULL
             y_COMM_WORLD   = MPI_COMM_NULL
             z_COMM_WORLD   = MPI_COMM_NULL
             w_COMM_WORLD   = MPI_COMM_NULL
             TRD_nnod = 0; TRD_inod = 0
             x_nnod   = 0; x_inod   = 0
             y_nnod   = 0; y_inod   = 0
             z_inod   = 0; z_nnod   = 0
             w_inod   = 0; w_nnod   = 0
          else
             call MPI_Comm_dup( comm0, TRD_COMM_WORLD, ierr )
          end if

!---- OpenMP thread check ---
!$        if ( TRD_COMM_WORLD /= MPI_COMM_NULL ) then
!$           call MPI_Query_thread( MPI_THREAD_MODE, ierr )
!$           local_size = 1
!$OMP PARALLEL
!$OMP MASTER
!$           local_size = omp_get_num_threads( )
!$OMP END MASTER
!$OMP END PARALLEL
!$           th0(1) = local_size
!$           th0(2) = -th0(1)
!$           call MPI_Allreduce( th0, th1, 2, MPI_INTEGER, MPI_MAX,
!$   &                           TRD_COMM_WORLD, ierr )
!$           j = th1(1) + th1(2)
!$        if ( j /= 0 ) then
!$           call MPI_Barrier( TRD_COMM_WORLD, ierr )
!$           call flush( 6 )
!$        if ( TRD_inod == 1 ) then
!$           print*,"*************"
!$           print*,"** CAUTION **"
!$           print*,"*************"
!$         print*,"EigenExa supports only homogeneous thread setting!"
!$           print*,"EigenExa terminates this run."
!$           print*,"*************"
!$        end if
!$        do j=1,2
!$           call MPI_Barrier( TRD_COMM_WORLD, ierr )
!$           call flush( 6 )
!$        end do
!$           call MPI_Barrier( TRD_COMM_WORLD, ierr )
!$           call sleep(1)
!$           call MPI_Abort( MPI_COMM_WORLD, MPI_ERR_OTHER, ierr )
!$        end if
!$        allocate ( tmp_save(0:local_size-1) )
!$        end if
!--------------------

!---- Cartesian check ---

          if ( TRD_COMM_WORLD /= MPI_COMM_NULL ) then

             call MPI_Comm_size( TRD_COMM_WORLD, TRD_nnod, ierr )
             call MPI_Comm_rank( TRD_COMM_WORLD, TRD_inod, ierr )
             TRD_inod = TRD_inod+1

             call MPI_Topo_test( TRD_COMM_WORLD, topo_type, ierr )
             if ( topo_type == MPI_CART ) then
                call MPI_Cartdim_get( TRD_COMM_WORLD, cart_dim, ierr )
             else
                cart_dim = 1
             end if

!---- 2D Cartesian setup ---

          if ( cart_dim == 2 ) then

             call MPI_Cart_get( TRD_COMM_WORLD, cart_dim,
     &                          dims, periods, coords, ierr )
             x_nnod = dims(1)
             y_nnod = dims(2)
             x_inod = coords(1) +1 ! fixed
             y_inod = coords(2) +1 ! fixed

             GRID_major = 'R'
          if ( PRESENT(order) ) then
          if ( x_inod == 1 .AND. y_inod == 1 ) then
          if ( order(1:1) == 'C' .or. order(1:1) == 'c' ) then
             print*,"*************"
             print*,"** CAUTION **"
             print*,"*************"
             print*,"The MPI_CART you specified is based on R-major,"
             print*,"but you also specified C-major option."
             print*,"EigenExa solve this conflict by taking R-major."
          end if
          end if
          end if

          else ! if ( cart_dim == 2 ) then

!---- Setup 2D process map ---

             x_nnod = INT(SQRT(DBLE(TRD_nnod)))
             i = 1 ! minimum factor, x_nnod must be multiple of k
             if ( MOD(TRD_nnod, i) == 0 ) then
                k = i
             else
                k = 1
             end if
             do
                if ( x_nnod <= k ) exit
                if ( MOD(x_nnod, k) == 0 .AND.
     &               MOD(TRD_nnod, x_nnod) == 0 ) exit
                x_nnod = x_nnod-1
             end do!!
             y_nnod = TRD_nnod/x_nnod

             if ( PRESENT(order) ) then
                GRID_major = order(1:1)
             else
                GRID_major = 'C'
             end if

             if ( GRID_major == 'R' .or. GRID_major == 'r' ) then
                GRID_major = 'R'
             else
                GRID_major = 'C'
             endif

             call eigen_set_grid_major( GRID_major )

             if ( GRID_major == 'R' ) then
! row-major
                x_inod =    (TRD_inod-1)/y_nnod +1
                y_inod = MOD(TRD_inod-1, y_nnod)+1
             else
! column-major
! ** EigenK adopts column-major in default
!    The process ordering on the Earth Simulator was done in the way
!    such as incrmental from internal processors to external nodes.
!    We want to make the processes, which possess a distributed vector,
!    close physically on the network connection.
!    These are the historical reason to adopt column-major.
                x_inod = MOD(TRD_inod-1, x_nnod)+1
                y_inod =    (TRD_inod-1)/x_nnod +1
             end if

          end if ! if ( cart_dim == 2 ) then

!---- x-Comm & y-Comm setup ---

             call MPI_Comm_split( TRD_COMM_WORLD, y_inod, x_inod,
     &            x_COMM_WORLD, ierr )
             call MPI_Comm_split( TRD_COMM_WORLD, x_inod, y_inod,
     &            y_COMM_WORLD, ierr )

!---- Parameter setup for datacast_dbl (row-column vector redistribution) ---

             n1 = MAX(x_nnod, y_nnod)
             n2 = MIN(x_nnod, y_nnod)
             do
                if ( n1 == n2 ) then
                   n_common = n1
                   exit
                end if
                n3 = n1-n2
                n1 = MAX(n2, n3)
                n2 = MIN(n2, n3)
             end do!!

             if ( ASSOCIATED( p0_ ) ) then
                deallocate( p0_ )
             end if
             allocate ( p0_(1:MAX(x_nnod, y_nnod)) )

             if ( ASSOCIATED( q0_ ) ) then
                deallocate( q0_ )
             end if
             allocate ( q0_(1:MAX(x_nnod, y_nnod)) )

             p0_(:) = -1
             q0_(:) = -1
             do i = 1, x_nnod
             if ( MOD(i-1, n_common) == MOD(y_inod-1, n_common) ) then
                n1 = y_inod-i
                if ( n1 >= 0 ) then
                   do j = 1, x_nnod
                      k = +n1+(j-1)*y_nnod
                      if ( MOD(k, x_nnod) == 0 ) then
                         p0_(i) = k/x_nnod
                         q0_(i) = (j-1)
                         exit
                      end if
                   end do! j
                else
                   do j = 1, y_nnod
                      k = -n1+(j-1)*x_nnod
                      if ( MOD(k, y_nnod) == 0 ) then
                         q0_(i) = k/y_nnod
                         p0_(i) = (j-1)
                         exit
                      end if
                   end do! j
                end if
             end if
             end do! i
             p0_(:) = p0_(:)+1
             q0_(:) = q0_(:)+1


             diag_0 = 0
             diag_1 = 0
             do i = 1, y_nnod/n_common
                j = (i-1)*y_nnod+y_inod
                k = MOD(j-1, x_nnod)+1
                if ( k == x_inod ) then
                   diag_0 = i
                   diag_1 = (j-1)/x_nnod+1
                   exit
                end if
             end do! i_1

!---- extra-Comm setup for datacast_dbl ---

             if ( n_common > 1 ) then
                w_inod =     (x_inod-1)/ n_common +1
                z_inod = MOD((x_inod-1), n_common)+1
                w_nnod = x_nnod/n_common
                z_nnod = n_common
                call MPI_Comm_split( x_COMM_WORLD, w_inod-1, z_inod-1,
     &            z_COMM_WORLD, ierr )
                call MPI_Comm_split( x_COMM_WORLD, z_inod-1, w_inod-1,
     &            w_COMM_WORLD, ierr )
             else
                z_COMM_WORLD = MPI_COMM_SELF
                w_COMM_WORLD = x_COMM_WORLD
                z_inod = 1; z_nnod = 1
                w_inod = x_inod; w_nnod = x_nnod
             end if

!---- Exceptional setup ---

          else !  if ( TRD_COMM_WORLD /= MPI_COMM_NULL )

             x_COMM_WORLD = MPI_COMM_NULL
             y_COMM_WORLD = MPI_COMM_NULL
             z_COMM_WORLD = MPI_COMM_NULL
             w_COMM_WORLD = MPI_COMM_NULL
             x_nnod = 0; x_inod = 0
             y_nnod = 0; y_inod = 0
             z_inod = 0; z_nnod = 0
             w_inod = 0; w_nnod = 0

             GRID_major = 'C'
             call eigen_set_grid_major( GRID_major )

          end if !  if ( TRD_COMM_WORLD /= MPI_COMM_NULL )

!---- BLACS setup ---

#if TIMER_PRINT>1
          if ( TRD_inod <= 1 ) then
             print*,"GRID major ",GRID_major," is specified."
          end if
#endif

          if ( TRD_COMM_WORLD == MPI_COMM_NULL ) then
             comm0  = MPI_COMM_SELF
             x_nnod = 1; x_inod = 1
             y_nnod = 1; y_inod = 1
          else
             comm0  = TRD_COMM_WORLD
          end if

             call Eigen_BLACS_Init( comm0, x_nnod, y_nnod, GRID_Major )

          if ( TRD_COMM_WORLD == MPI_COMM_NULL ) then
             x_nnod = 0; x_inod = 0
             y_nnod = 0; y_inod = 0
          end if

!--------------------

!---- Measure the overhead of collective communications ---

          allocate( buff1(2048), buff2(2048) )

!--------------------

          call MPI_Barrier(TRD_COMM_WORLD,ierr)

          s1 = ZERO
          do  i = 1, 3
             call MPI_Barrier(x_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call MPI_Barrier(x_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s1 = s1 + (d2-d1)
          end do

          call MPI_Barrier(TRD_COMM_WORLD,ierr)

          s2 = ZERO
          do  i = 1, 3
             call MPI_Barrier(y_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call MPI_Barrier(y_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s2 = s2 + (d2-d1)
          end do

          call MPI_Barrier(TRD_COMM_WORLD,ierr)

          buff(1) = s1 / 2
          buff(2) = s2 / 2
          call MPI_Allreduce( buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, TRD_COMM_WORLD, ierr )
          Barrier_Overhead_x = buff(3) / TRD_nnod
          Barrier_Overhead_y = buff(4) / TRD_nnod

!--------------------

          s1 = ZERO
          do  i = 1, 3
             call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call reduce_dbl(buff1(1),buff2(1),2048,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s1 = s1 + (d2-d1)
          end do

          s2 = ZERO
          do  i = 1, 3
             call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call reduce_dbl(buff1(1),buff2(1),2048,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s2 = s2 + (d2-d1)
          end do

          buff(1) = s1 / 2
          buff(2) = s2 / 2
          call MPI_Allreduce( buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, TRD_COMM_WORLD, ierr )
          Reduce_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
          Reduce_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y

          s1 = ZERO
          do  i = 1, 3
             call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
             call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s1 = s1 + (d2-d1)
          end do

          s2 = ZERO
          do  i = 1, 3
             call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
             call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s2 = s2 + (d2-d1)
          end do

          buff(1) = s1 / 2
          buff(2) = s2 / 2
          call MPI_Allreduce( buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, TRD_COMM_WORLD, ierr )
          Reduce_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
     &                      - Reduce_Overhead_x
          Reduce_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y
     &                      - Reduce_Overhead_y

!--------------------

          s1 = ZERO
          do  i = 1, 3
             call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call bcast_dbl(buff1(1),2048,1,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s1 = s1 + (d2-d1)
          end do

          s2 = ZERO
          do  i = 1, 3
             call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call bcast_dbl(buff1(1),2048,1,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s2 = s2 + (d2-d1)
          end do

          buff(1) = s1 / 2
          buff(2) = s2 / 2
          call MPI_Allreduce( buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, TRD_COMM_WORLD, ierr )
          Bcast_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
          Bcast_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y

          s1 = ZERO
          do  i = 1, 3
             call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
             call bcast_dbl(buff2(1),1024,1,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s1 = s1 + (d2-d1)
          end do

          s2 = ZERO
          do  i = 1, 3
             call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d1 = eigen_get_wtime()
             call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
             call bcast_dbl(buff2(1),1024,1,0,y_COMM_WORLD)
             call MPI_Barrier(TRD_COMM_WORLD,ierr)
             d2 = eigen_get_wtime()
             if ( i /= 1 ) s2 = s2 + (d2-d1)
          end do

          buff(1) = s1 / 2
          buff(2) = s2 / 2
          call MPI_Allreduce( buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, TRD_COMM_WORLD, ierr )
          Bcast_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
     &                     - Bcast_Overhead_x
          Bcast_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y
     &                     - Bcast_Overhead_y

!--------------------

          if ( Reduce_Overhead_x < ZERO ) Reduce_Overhead_x = ZERO
          if ( Reduce_Overhead_y < ZERO ) Reduce_Overhead_y = ZERO
          if ( Bcast_Overhead_x  < ZERO ) Bcast_Overhead_x  = ZERO
          if ( Bcast_Overhead_y  < ZERO ) Bcast_Overhead_y  = ZERO

          deallocate( buff1, buff2 )

!--------------------

             call eigen_set_initialized( )

!--------------------

          return
          end subroutine  eigen_init
!--
!--
          subroutine eigen_free( flag )
          use MPI
          use eigen_devel
          use eigen_blacs, only : eigen_blacs_exit
          implicit NONE
          integer, intent(in), optional ::  flag
*-
          integer                ::  ierr
          logical                ::  local_flag


          call eigen_get_initialized( local_flag )

          if ( .NOT. local_flag ) then
             return
          end if

          if ( TRD_COMM_WORLD /= MPI_COMM_NULL ) then

             call MPI_Comm_free( x_COMM_WORLD, ierr )
             call MPI_Comm_free( y_COMM_WORLD, ierr )

             if ( n_common > 1 ) then
             call MPI_Comm_free( z_COMM_WORLD, ierr )
             call MPI_Comm_free( w_COMM_WORLD, ierr )
             end if

             call MPI_Comm_free( TRD_COMM_WORLD, ierr )

             deallocate ( p0_)
             deallocate ( q0_)

             if ( PRESENT(flag) ) then
                if ( flag == 1 ) then
                   call eigen_timer_print( 'EigenExa(finalized)' )
                end if
             end if
             call eigen_timer_reset( 0, 0, 0, 0 )

!$           deallocate( tmp_save )

          end if

!             call Eigen_BLACS_Exit( )

             TRD_COMM_WORLD = MPI_COMM_WORLD

             call eigen_unset_initialized( )


          return
          end subroutine  eigen_free
!--
!--
          subroutine eigen_get_matdims( n, nx, ny )
          use eigen_devel
          implicit NONE
          integer, intent(in)    :: n
          integer, intent(out)   :: nx, ny
*-
          integer :: NPROW, NPCOL, NB
          integer :: n1, nm, nmz, nmw, nn, larray


             NPROW = x_nnod
             NPCOL = y_nnod

             n1 = ((n-1)/NPROW+1)
             call CSTAB_get_optdim( n1, 6, 16*4, 16*4*2, nm )

             NB  = eigen_NB

             nmz = ((n-1)/NPROW+1)
             nmz = ((nmz-1)/NB+1)*NB+1
             nn  = nmz
             nmz = (n-1)/NB+1
             nmz = ((nmz-1)/NPROW+1)*NB
! Fix on version 2.2b
! to avoid unexpected SIGSEGV,
! use the maximum of nn and nmz.
             nmz = MAX(nn, nmz)

             nmw = ((n-1)/NPCOL+1)
             nmw = ((nmw-1)/NB+1)*NB+1
             nn  = nmw
             nmw = (n-1)/NB+1
             nmw = ((nmw-1)/NPCOL+1)*NB
! Fix on version 2.2b
! to avoid unexpected SIGSEGV,
! use the maximum of nn and nmz.
             nmw = MAX(nn, nmw)

             larray = MAX(nmz, nm)*nmw

             nx = nm
             ny = (larray-1)/nm+1


          return
          end subroutine  eigen_get_matdims
!--
!--
          integer function eigen_memory_internal( n, lda, ldz,
     &                     m1_opt, m0_opt )
     &                     result( byte )
!$        use OMP_LIB
          use MPI
          use eigen_devel
          implicit NONE
          integer, intent(in)    :: n, lda, ldz
          integer, intent(in), optional :: m1_opt, m0_opt
*-
          integer :: m1, m0
          integer :: l_uv_t
          integer :: MYROW, MYCOL, NPROW, NPCOL, NP, NQ, NB
          integer :: lddz, lddw, LWORK, LIWORK
          integer :: m, na, nx, nv, nm, ierr
          integer :: local_rank, local_size

          integer, parameter     :: nm_max_L1 = 16*4
          integer, parameter     :: nm_max_L2 = 16*6

          integer :: byte_tridi
          integer :: byte_dcx
          integer :: byte_trbak

          integer, external :: NUMROC

          include 'CSTAB.h'


             local_rank = 0
             local_size = 1
!$           local_rank = omp_get_thread_num()
!$           local_size = omp_get_num_threads()

          if ( PRESENT(m1_opt) ) then
             m1 = m1_opt
          else
             m1 = 48
          end if

          if ( PRESENT(m0_opt) ) then
             m0 = m0_opt
          else
             m0 = 128
          end if

!  ===== for TRIDI ====

             nx = (n-1)/x_nnod+1 +2
             call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv)

             l_uv_t = MAX(2*nv,
     &                    2*(n+MAX(3,y_nnod/n_common)),
     &                    2*8*m1,
     &                    nv+3+4*m1)

             LWORK = nm*m1+n_columns
     &             + 2*(l_uv_t+n_columns)
     &             + (nv+n_columns)
     &             + 2*(nv*2*m1+2*n_columns)
     &             + 4*(nv*local_size+n_columns)

             byte_tridi = LWORK * 8

!  ===== for DCX ====

             NB = eigen_NB
             NB = MIN(NB, N)

             NPROW = x_nnod
             NPCOL = y_nnod
             MYROW = x_inod-1
             MYCOL = y_inod-1
             NP = NUMROC( n, NB, MYROW, 0, NPROW )
             NQ = NUMROC( n, NB, MYCOL, 0, NPCOL )

             lddz = (n-1)/NPROW+1
             lddz = ((lddz-1)/NB+1)*NB
             lddw = (n-1)/NPCOL+1
             lddw = ((lddw-1)/NB+1)*NB

             nx     = (N-1)/NPCOL+1
             LWORK  = MAX(1+6*N+2*NP*(NQ+MAX(NQ,NB)), lddz*lddw, ldz*nx)
     &              + n
             LIWORK = 2+7*n+8*NPCOL

             byte_dcx = LWORK * 8 + LIWORK * 4

!  ===== for TRBAK ====

             m    = MIN(nsm, m0)
             if ( m < 1 ) m = 1

             na   = (n-1)/y_nnod+1
             na   = na  +MOD(na  -1,2)
             call CSTAB_get_optdim(lda, 9, 16*4, 16*6, nm)

             LWORK = n
     &             + 3*(MAX(nm*m,n)+n_columns)
     &             + 4*(na*m+ns0+n_columns)
     &             + m*nm
             LIWORK = MAX(((m-1)/y_nnod+1)*((n-1)/x_nnod+1),n)

             byte_trbak = LWORK * 8 + LIWORK * 4


             na = MAX(byte_tridi, byte_dcx, byte_trbak)
             call MPI_Allreduce( na, nx, 1, MPI_INTEGER,
     &                           MPI_MAX, TRD_COMM_WORLD, ierr )
             byte = nx


          return
          end function  eigen_memory_internal
!--
!--
          subroutine eigen_get_procs( procs, x_procs, y_procs )
          use eigen_devel
          implicit NONE
          integer, intent(out)   ::  procs, x_procs, y_procs


             procs   = TRD_nnod
             x_procs = x_nnod
             y_procs = y_nnod


          return
          end subroutine  eigen_get_procs
!--
!--
          subroutine eigen_get_id( id, x_id, y_id )
          use eigen_devel
          implicit NONE
          integer, intent(out)   ::  id, x_id, y_id


             id   = TRD_inod
             x_id = x_inod
             y_id = y_inod


          return
          end subroutine  eigen_get_id
!--
!--
          subroutine eigen_get_comm( comm, x_comm, y_comm )
          use eigen_devel
          implicit NONE
          integer, intent(out)   ::  comm, x_comm, y_comm


             comm   = TRD_COMM_WORLD
             x_comm = x_COMM_WORLD
             y_comm = y_COMM_WORLD


          return
          end subroutine  eigen_get_comm
!--
!--
          integer function eigen_loop_start( istart, nnod, inod )
     &            result( ret )
          implicit NONE
          integer, intent(in)    ::  istart, nnod, inod


!             ret = eigen_loop_calc( istart, nnod, inod, 1 )
             ret = (istart + nnod - 1 - inod) / nnod + 1

          return
          end function  eigen_loop_start
!--
!--
          integer function eigen_loop_end( iend, nnod, inod )
     &            result( ret )
          implicit NONE
          integer, intent(in)    ::  iend, nnod, inod

!             ret   = eigen_loop_calc( iend, nnod, inod, 0 )
             ret = (iend + nnod - 0 - inod) / nnod + 0

          return
          end function  eigen_loop_end
!--
!--
          integer function eigen_translate_l2g( ictr, nnod, inod )
     &            result( ret )
          implicit NONE
          integer, intent(in)    ::  ictr, nnod, inod

             ret = (ictr-1)*nnod+inod

          return
          end function  eigen_translate_l2g
!--
          integer function eigen_translate_g2l( ictr, nnod, inod )
          implicit NONE
          integer, intent(in)    ::  ictr, nnod, inod

             eigen_translate_g2l = (ictr-1)/nnod+1

          return
          end function  eigen_translate_g2l
!--
!--
          integer function eigen_owner_node( ictr, nnod, inod )
     &            result( ret )
          implicit NONE
          integer, intent(in)    ::  ictr, nnod, inod

             ret = MOD(ictr-1, nnod)+1

          return
          end function  eigen_owner_node
!--
!--
!          integer function eigen_loop_calc( ictr, nnod, inod, d )
!     &            result( ret )
!          implicit NONE
!          integer, intent(in)    ::  ictr, nnod, inod, d
!
!             ret = (ictr+nnod-d-inod)/nnod+d
!
!          return
!          end function  eigen_loop_calc
!--
!--
          subroutine  eigen_vector_zeropad_x( v_x, L )
          use eigen_devel
          implicit NONE

          real(8), intent(inout) :: v_x(*)
          integer, intent(in)    :: L
          integer                :: j_3, n

          j_3 = eigen_loop_end     (L, x_nnod, x_inod)
          n   = eigen_translate_g2l(L, x_nnod, x_inod)
          if ( j_3 < n ) then
!OCL XFILL
             v_x(j_3+1:n) = ZERO ! in case
          end if

          end subroutine  eigen_vector_zeropad_x
!--
!--
          subroutine  eigen_vector_zeropad_y( v_y, L )
          use eigen_devel
          implicit NONE

          real(8), intent(inout) :: v_y(*)
          integer, intent(in)    :: L
          integer                :: j_3, n

          j_3 = eigen_loop_end     (L, y_nnod, y_inod)
          n   = eigen_translate_g2l(L, y_nnod, y_inod)
          if ( j_3 < n ) then
!OCL XFILL
             v_y(j_3+1:n) = ZERO ! in case
          end if

          end subroutine  eigen_vector_zeropad_y
!--
*
* NOTICE : THIS SOURCE IS LIMITTED FOR THE MICROPROCESSORS
*          WHICH ADOPT THE IEEE754 FLOATING POINT ARITHMETIC FORMAT
*          AND OPERATIONS
*
       real(8) function machine_epsilon( ) result( r )
       implicit NONE
       real(8) :: eps
       data eps /z'3CB0000000000000'/

          r = eps

       return
       end function  machine_epsilon
!--
       end module eigen_libs
