      program EigenExa_benchmark
!$    use OMP_LIB
      use MPI
      use eigen_libs
      implicit NONE

      real(8), pointer   :: a(:,:),z(:,:),w(:)
      real(8), parameter :: ONE = 1D0
      real(8)            :: PAI, EPS, EPS2, EPS4
      real(8)            :: d1, d2, flops, time
!
      integer            :: nnod, x_nnod, y_nnod
      integer            :: inod, x_inod, y_inod
      integer            :: i, i_inod, i_nnod, ierror, istat
      integer            :: n, m, mb, nall, mtype, msolver
      integer            :: nm, nx, ny, imem
      integer            :: new_comm, color, key, n_comms
      integer            :: new_dim_x, new_dim_y, dims(2), coords(2)
      logical            :: periods(2), reorder
      character*1        :: mode = ' ', grid = ' '
      character*256      :: input_file = ' '
      logical            :: check_accuracy = .TRUE.
      character*1024     :: in_buffer
!
      integer            :: nargc
      character*256      :: argv
      real(8)            :: times(8)


*-
*----------------------------------------------------------------------
*-
      call MPI_Init_thread( MPI_THREAD_MULTIPLE, i, ierror )
      call MPI_Comm_rank( MPI_COMM_WORLD, i_inod, ierror )
      call MPI_Comm_size( MPI_COMM_WORLD, i_nnod, ierror )
*-
*----------------------------------------------------------------------
*-
      nargc = command_argument_count( )
      i = 0
      do
         i = i+1
         call get_command_argument( i, argv )
         if ( len_trim(argv) == 0 ) exit
         argv = trim(argv)

         if ( argv(1:2) == '-h' ) then
         if ( i_inod == 0 ) then

      print*,""
      print*," eigenexa_benchmark [options]"
      print*,""
      print*,"options:"
      print*,""
      print*," -h              displays this help and exit"
      print*," -L              displays the list of test matrices"
      print*," -f input_file   uses input_file"
      print*,"                   default is ./IN"
      print*," -c or -n        check or non-check accuracy."
      print*,"                   default is check."
      print*," -g mode         sets the process grid as follows"
      print*,"    R, r           MPI_COMM_WORLD row-major mode"
      print*,"    C, c           MPI_COMM_WORLD column-major mode"
      print*,"    A, a           MPI_COMM_SELF (embarrasingly parallel)"
      print*,"    1, 2,... 9     splitted MPI_COMM_WORLD"
      print*,"                      with the color=MOD(RANK,{number})"
      print*," -x dimX dimY    sets the cartecian shape (dimX, dimY)"
      print*,"                      dimX <= dimY must be hold."
      print*,""

         end if
         goto 99999
         end if

         if ( argv(1:2) == '-L' ) then
         if ( i_inod == 0 ) then

            print*,""
            print*," eigenexa_benchmark [options]"
            print*,""
            print*,"The list of test matrices:"
            print*,""
            call print_mat_name_list()
            print*,""

         end if
         goto 99999
         end if

         if ( argv(1:2) == '-c' ) then
            check_accuracy = .TRUE.
            cycle
         end if

         if ( argv(1:2) == '-n' ) then
            check_accuracy = .FALSE.
            cycle
         end if

         if ( argv(1:2) == '-f' ) then
            i = i+1
            call get_command_argument( i, argv )
            if ( len_trim(argv) == 0 ) exit
            input_file = trim(argv)
            cycle
         end if

         if ( argv(1:2) == '-g' ) then
            i = i+1
            call get_command_argument( i, argv )
            if ( len_trim(argv) == 0 ) exit
            grid = argv(1:1)
            if ( grid == 'x' ) grid = ' '
            cycle
         end if

         if ( argv(1:2) == '-x' ) then
            i = i+1
            call get_command_argument( i, argv )
            if ( len_trim(argv) == 0 ) exit
            read( argv, * ) new_dim_x
            i = i+1
            call get_command_argument( i, argv )
            if ( len_trim(argv) == 0 ) exit
            read( argv, * ) new_dim_y
            grid = 'x'
            cycle
         end if

         print*,"Command option '", argv(1:2), "' is inavailable. ",
     &     "This benchmark code skips it, and continues to proceed."

      end do
*
! User can specify the process grid shape by -g option
!
      n_comms = 1
      color   = 0
      select case ( grid )
      case ( 'R', 'r' )
         call eigen_init( order='r' )
      case ( 'C', 'c' )
         call eigen_init( order='c' )
      case ( 'A', 'a' )
         n_comms = i_nnod
         color   = i_inod
         call eigen_init( MPI_COMM_SELF )
      case ( '1', '2', '3', '4', '5', '6', '7', '8', '9' )
         read( grid, * ) n_comms
         color = MOD(i_inod, n_comms)
         key   = i_inod / n_comms
         call MPI_Comm_split( MPI_COMM_WORLD,
     &                        color, key,
     &                        new_comm, ierror )
         if ( color /= 0 ) new_comm = MPI_COMM_NULL
         call eigen_init( new_comm )
      case ( 'x' )
         if ( new_dim_x * new_dim_y /= i_nnod ) then
         if ( i_inod == 0 ) then
            print*,"Illegal dimensions are specified."
         end if
            call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
         end if
         if ( new_dim_x > new_dim_y ) then
         if ( i_inod == 0 ) then
            print*,"This process map is not supported."
            print*,"Px should be smaller than Py"
         end if
            call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
         end if
         dims(1:2)    = ( / new_dim_x, new_dim_y / )
         periods(1:2) = .FALSE.
         reorder      = .FALSE.
         call MPI_Cart_create( MPI_COMM_WORLD,
     &                         2, dims, periods, reorder,
     &                         new_comm, ierror )
         call eigen_init( new_comm )
      case default
         call eigen_init( )
      end select
*
! -f option specifies the input file, default is 'IN'
!
      if ( input_file(1:1) == ' ' ) then
         input_file = 'IN'
      end if
*-
*----------------------------------------------------------------------
*-
      if ( i_inod == 0 ) then
         open( unit=10, file=input_file, status='old', recl=1024,
     &             iostat=ierror )
         if ( ierror /= 0 ) then
            print*,"Can not open input file"
            call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
         end if
         print*,"INPUT FILE=",trim(input_file)
      end if
*-
*----------------------------------------------------------------------
*-
      call eigen_get_procs( nnod, x_nnod, y_nnod )
      call eigen_get_id   ( inod, x_inod, y_inod )
*-
*----------------------------------------------------------------------
*-
      PAI = 4*ATAN(ONE)
      EPS = machine_epsilon()
      EPS2 = SQRT(EPS)
      EPS4 = SQRT(EPS2)
*-
*----------------------------------------------------------------------
*-
      DO

         if ( i_inod == 0 ) then
!
! Input file format
!
! N bx by mode matrix solver
!
! N      : matrix dimension
! bx     : block width for the forward transformation
! by     : block width for the backward transformation
! mode   : eigensolver mode { 0 : only eigenvalues }
!                           { 1 : eigenvalues and corresponding eigenvectors}
!                           { 2 : mode 1 + accuracy improvement for eigenvalue}
! matrix : test matrix { 11 types, 0 ... 10 }
! solver : { 0 : eigen_sx, new algorithm, faster on the K }
!          { 1 : eigen_s,  conventional algorithm }
!
! if a line starts from '!', the line is treated as a comment
!
            do
               read( unit=10, fmt='(a1024)',
     &                err=10000, end=10000 ) in_buffer
               if ( in_buffer(1:1) /= '!' ) EXIT
            end do
            read( in_buffer, fmt=*,
     &                err=10000, end=10000 )
     &                n, m, mb, nall, mtype, msolver
         end if
         goto 20000
10000    continue
         n = -1
20000    continue

         call MPI_Bcast( n, 1, MPI_INTEGER,
     &                  0, MPI_COMM_WORLD, ierror )
         if ( n <= 0 ) exit
         call MPI_Bcast( m, 1, MPI_INTEGER,
     &                  0, MPI_COMM_WORLD, ierror )
         call MPI_Bcast( mb, 1, MPI_INTEGER,
     &                  0, MPI_COMM_WORLD, ierror )
         call MPI_Bcast( nall, 1, MPI_INTEGER,
     &                  0, MPI_COMM_WORLD, ierror )
         call MPI_Bcast( mtype, 1, MPI_INTEGER,
     &                  0, MPI_COMM_WORLD, ierror )
         call MPI_Bcast( msolver, 1, MPI_INTEGER,
     &                  0, MPI_COMM_WORLD, ierror )

         mode = 'A'
         select case (nall)
         case (0)
            mode = 'N' ! only eigenvalues, no eigenvector
         case (1)
            mode = 'A' ! all the eigenpairs
         case (2)
            mode = 'X' ! mode 'A' + accuracy improvement
         end select
*-
*----------------------------------------------------------------------
*-
         if ( new_comm /= MPI_COMM_NULL ) then
*
         call eigen_get_matdims( n, nm, ny )
         imem =  eigen_memory_internal( n, nm, nm, m, 128 )
*
         allocate(
     &               a(nm, ny),
     &               z(nm, ny),
     &               w(n),
     &               stat=istat )
         if ( istat /= 0 ) then
            print*,"Memory exhausted"
            call flush(6)
            EXIT
         end if
*
         call mat_set( n, a(1,1), nm, mtype )
*
         end if
*-
*----------------------------------------------------------------------
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
         d1 = MPI_Wtime( )
*
         if ( new_comm /= MPI_COMM_NULL ) then
*
         if ( msolver == 0 ) then
         call eigen_sx( n, n, a, nm, w, z, nm, times,
     &                  m_forward=m, m_backward=mb, mode=mode )
         else
         call eigen_s ( n, n, a, nm, w, z, nm, times,
     &                  m_forward=m, m_backward=mb, mode=mode )
         end if
         flops = a(1, 1)
*
!         do i=1,n
!            print*,w(i)
!         enddo
*
         end if
*
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
         d2 = MPI_Wtime( )
         time = d2 -d1
*-
*----------------------------------------------------------------------
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
*-
         if ( i_inod == 0 ) then
         print*,"======================================================"
*-
            if ( msolver == 0 ) then
            print*,"Solver = eigen_sx / via penta-diagonal format"
            else
            print*,"Solver = eigen_s  / via tri-diagonal format"
            end if
            print*,"Block width = ", m, "/", mb
            print*,"NUM.OF.PROCESS=",nnod,"(",x_nnod,y_nnod,")"
!$OMP PARALLEL
!$OMP MASTER
!$          print*,"NUM.OF.THREADS=", omp_get_num_threads( )
!$OMP END MASTER
!$OMP END PARALLEL
            print*,"Matrix dimension = ",N
            call print_mat_name( mtype )
            print*,"Internally required memory = ",imem," [Byte]"
            select case (nall)
            case (0)
               print*, "mode 'N' :: only eigenvalues, no eigenvector"
            case (1)
               print*, "mode 'A' :: all the eigenpairs"
            case (2)
               print*, "mode 'X' :: mode 'A' + accuracy improvement"
            end select
            print*,"Elapsed time = ",time," [sec]"
            print*,"FLOP = ",flops
            print*,"Performance = ",(flops/time)*1D-9," [GFLOPS]"
*-
         if ( check_accuracy ) then
            call w_test( n, w, mtype )
         end if
*-
         end if
*-
         if ( check_accuracy ) then
         if ( mode == 'X' .OR. mode == 'A' ) then
         do i=0,n_comms-1
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
*-
         if ( color == i ) then
         if ( new_comm /= MPI_COMM_NULL ) then
            if ( n_comms > 1 ) then
               if ( inod == 1 ) then
                  print*,"> Group ",i
               end if
            end if
            call mat_set( n, a(1,1), nm, mtype )
            call ev_test( n, a(1,1), nm, w(1), z(1,1), nm )
         end if
         end if
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
*-
         end do
         end if
         end if
*-
         if ( i_inod == 0 ) then
         print*,"======================================================"
         print*,""
         end if
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
*-
*----------------------------------------------------------------------
*-
         if ( new_comm /= MPI_COMM_NULL ) then
         deallocate( a, z, w )
         end if
*-
         call MPI_Barrier( MPI_COMM_WORLD, ierror )
*-
      end do
*-
*----------------------------------------------------------------------
*-
      call MPI_Barrier( MPI_COMM_WORLD, ierror )
*-
*----------------------------------------------------------------------
*-
      if ( i_inod == 0 ) then
         close( unit=10 )
      end if
*-
*----------------------------------------------------------------------
*-
      call eigen_free( )
*-
*----------------------------------------------------------------------
*-
      if ( i_inod == 0 ) then
         print*,"Benchmark completed"
      end if
*-
*----------------------------------------------------------------------
*-
99999 continue
*-
*----------------------------------------------------------------------
*-
      call MPI_Finalize( ierror )
*-
*----------------------------------------------------------------------
*-
      end program EigenExa_benchmark
