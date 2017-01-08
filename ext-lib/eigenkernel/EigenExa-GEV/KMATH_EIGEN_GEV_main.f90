program main
  use MPI
  use eigen_libs
  use distributed_read_m
  implicit NONE

  logical,parameter :: check_flag = .true., read_from_file = .true.
  real(8),allocatable :: a(:),b(:),z(:),w(:)
  integer :: ierr,inod,nnod
  integer :: n
  integer :: NPROW,NPCOL,nx,nm,nmw, ifl,nfl
  integer :: larray_a,larray_b,larray_z
  real(8) :: s_time,e_time
  integer :: lda,ldb,ldz
  integer :: desc(9), info
  character(len=1024) :: filename_matA, filename_matB

  integer :: numroc,iceil
  external :: numroc,iceil

  call MPI_Init_thread( MPI_THREAD_MULTIPLE, nx, ierr )
  call MPI_Comm_rank( MPI_COMM_WORLD, inod, ierr )
  call MPI_Comm_size( MPI_COMM_WORLD, nnod, ierr )

  call MPI_Barrier( MPI_COMM_WORLD, ierr )

  call eigen_init( )
  call eigen_get_procs( ierr, NPROW, NPCOL )

  if (read_from_file) then
    call get_command_argument(1, filename_matA)
    call get_command_argument(2, filename_matB)
    call distributed_read_header(trim(filename_matA), n)
  else
    n = 30
  end if
  call eigen_get_matdims( n, nm, nmw )

  lda = nm
  ldb = nm
  ldz = nm

  larray_a = lda*nmw
  larray_b = ldb*nmw
  larray_z = ldz*nmw

  allocate(a(larray_a), b(larray_b), z(larray_z), w(n))

  s_time = MPI_Wtime( )

  if (read_from_file) then
    call descinit(desc, n, n, 1, 1, 0, 0, eigen_get_blacs_context(), lda, info)
    call distributed_read(trim(filename_matA), desc, a)
    call distributed_read(trim(filename_matB), desc, b)
  else
    call mat_set( n, a, lda, 2 )
    !call mat_set( n, b, ldb, 10 )
    call mat_set( n, b, ldb, 3 )
  end if
  e_time = MPI_Wtime( )
  if(inod==0)then
    write(*,*)"Matset time:",e_time-s_time
  end if

  call MPI_Barrier( MPI_COMM_WORLD, ierr )
  s_time = MPI_Wtime( )

  call KMATH_EIGEN_GEV( n, a, lda, b, ldb, w, z, ldz )

  call MPI_Barrier( MPI_COMM_WORLD, ierr )
  e_time = MPI_Wtime( )

  if(inod==0)then
    write(*,*)"Solver time:",e_time-s_time
  end if

  if(check_flag) then
    s_time = MPI_Wtime( )
    if (read_from_file) then
      call distributed_read(trim(filename_matA), desc, a)
      call distributed_read(trim(filename_matB), desc, b)
    else
      call mat_set( n, a, lda, 2 )
      !call mat_set( n, b, ldb, 10)
      call mat_set( n, b, ldb, 3)
    end if

    e_time = MPI_Wtime( )
    if(inod==0)then
      write(*,*)"Matset time:",e_time-s_time
    end if
    ifl=0
    nfl=0
    s_time = MPI_Wtime( )
    call KMATH_EIGEN_GEV_check(n,a,lda,b,ldb,w,z,ldz,ifl,nfl)
    e_time = MPI_Wtime( )
    if(inod==0)then
      write(*,*)"Checking time:",e_time-s_time
    end if
  end if

  deallocate(a,b,z,w)

  call MPI_Finalize( ierr )
end program main
