      use MPI
      use eigen_libs
      implicit double precision (a-h,o-z)

      real(8), pointer ::
     &           a(:,:),z(:,:),w(:)

      real(8), parameter :: ONE = 1D0
      real(8) :: PAI
!
      integer :: x_nnod, y_nnod, x_inod, y_inod
      integer :: comm, comm1, comm2
      integer :: orig_group, group1, group2 
      integer, allocatable :: ranks1(:), ranks2(:)
!
      call mpi_init_thread( MPI_THREAD_MULTIPLE, i, ierr )
      call mpi_comm_rank( mpi_comm_world, i_inod, ierr )
      call mpi_comm_size( mpi_comm_world, i_nnod, ierr )
*
      allocate(ranks1(i_nnod/2), ranks2(i_nnod/2))
      do i=0, i_nnod/2 - 1
         ranks1(i+1) = 2 * i
         ranks2(i+1) = 2 * i + 1
      end do
      print*, "ranks1=", ranks1
      print*, "ranks2=", ranks2

      call mpi_comm_group(mpi_comm_world, orig_group, ierr);
      call mpi_group_incl(orig_group, i_nnod / 2, ranks1, group1, ierr);
      call mpi_group_incl(orig_group, i_nnod / 2, ranks2, group2, ierr);
      
      call mpi_comm_create(mpi_comm_world, group1, comm1, ierr);
      call mpi_comm_create(mpi_comm_world, group2, comm2, ierr);

      if ( mod(i_inod,2) == 0 ) then
         MPI_COMM_EIGEN = comm1
      else
         MPI_COMM_EIGEN = comm2 !mpi_comm_null
      endif
*     
!     call eigen_init( MPI_COMM_WORLD, order='R' )
!     call Eigen_BLACS_Init( MPI_COMM_WORLD, x_nnod, y_nnod, 'R' )
      call eigen_init( MPI_COMM_EIGEN, order='R' )
!     call eigen_init( MPI_COMM_WORLD, order='C' )
!     call eigen_init( MPI_COMM_WORLD )
!     call eigen_init( order='R' )
!     call eigen_init( )
      
      call eigen_get_procs( nnod, x_nnod, y_nnod )
      call eigen_get_id   ( inod, x_inod, y_inod )

      if ( i_inod == 0 ) then
         open(10,file='IN')
      end if
*
      m = 48; nall = 2; mtype = 2

      PAI = 4*ATAN(ONE)
!
      DO

         if ( i_inod == 0 ) then
            read(10,*) n0, m, nall, mtype
         end if

         call mpi_bcast(n0,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
         if ( n0 < 1 .or. n0 > 400000 ) EXIT
         call mpi_bcast(m,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
         call mpi_bcast(nall,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
         call mpi_bcast(mtype,1,mpi_integer,0,MPI_COMM_WORLD,ierr)

!         n = INT(n0*(DBLE(i_nnod)**(1D0/3)))
         n = n0
         nall = 6
         if ( n0 > 1000 ) nall = 2

         call eigen_get_matdims( n, nm, ny )

         allocate(
     &               a(nm, ny),
     &               z(nm, ny),
     &               w(n),
     &               stat=istat)
         if ( istat /= 0 ) then
            print*,"Memory exhausted"
            call flush(6)
            EXIT
         end if
*
         call mat_set( n, a(1,1), nm, mtype )
*-
         i =  eigen_memory_internal( n, nm, nm, m, 128 )

         call MPI_BARRIER( MPI_COMM_WORLD, ierr )
         d1 = MPI_WTIME( )

!!      if (icolor == 1) then
         call eigen_sx( n, n, a, nm, w, z, nm,
     &                  m_forward=m, m_backward=128 )
!!      endif

         call MPI_BARRIER( MPI_COMM_WORLD, ierr )
         d2 = MPI_WTIME( )
*-
!!      if (icolor == 1) then
         if ( inod == 0 ) then
            print*,"Matrix dimension = ",N
            print*,"Internally required memory = ",i," [Byte]"
            print*,"Elapsed time = ",d2-d1," [sec]"
         end if
!!      endif
*-
!!      if (icolor == 1) then
         if ( mtype == 0 ) then
         if ( inod == 0 ) then
            ax=0D0
            do i = 1, n
               j = n-i
               theta = PAI*(2*j+1)/(2*n+1)
               x = 5D-1/(1D0-COS(theta))
               ax = MAX(ax, ABS(w(i)-x)/ABS(x))
            end do
            print*, "max|w(i)-w(i).true|/|w|=", ax
         end if
         end if
         call flush( 6 )
*-
         if ( IAND(nall, 4) /= 0 ) then
            call mat_set( n, a(1,1), nm, mtype )
            call ev_test( n, a(1,1), nm, w(1), z(1,1), nm, t )
         end if
!!      endif
*
         deallocate( a )
         deallocate( z )
         deallocate( w )

      END DO

      if ( i_inod == 0 ) then
         close(10)
      end if
*
      call eigen_free( )
*
      call MPI_Finalize( ierr )
      end


