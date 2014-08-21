program test_real2
  use ELPA1
  use ELPA2
  !#ifdef WITH_OPENMP
  !  use test_util
  !#endif
  implicit none
  include 'mpif.h'

  !-------------------------------------------------------------------------------
  ! Please set system size parameters below!
  ! na:   System size
  ! nev:  Number of eigenvectors to be calculated
  ! nblk: Blocking factor in block cyclic distribution
  !-------------------------------------------------------------------------------

  integer :: nblk
  integer :: na, nev

  !-------------------------------------------------------------------------------
  !  Local Variables

  integer :: np_rows, np_cols, na_rows, na_cols
  integer :: myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
  integer :: i, mpierr, context, sc_desc(9), info, nprow, npcol, ierr

  integer, external :: numroc

  double precision :: err, errmax
  double precision, allocatable :: a(:,:), b(:, :), z(:,:), tmp1(:,:), tmp2(:,:), &
       as(:,:), bs(:, :), ev(:)

  integer :: iseed(4096) ! Random seed, size should be sufficient for every generator
  integer :: status
  !#ifdef WITH_OPENMP
  !      integer :: omp_get_max_threads,  required_mpi_thread_level, &
  !      provided_mpi_thread_level
  !#endif
  !-------------------------------------------------------------------------------
  !  Parse command line argumnents, if given
  character*16 :: arg1, arg2, arg3
  integer, parameter   :: error_unit = 6
  logical :: success
  double precision :: times(11)

  nblk = 16
  na = 200
  nev = 10

  if (COMMAND_ARGUMENT_COUNT() == 3) then
    call GET_COMMAND_ARGUMENT(1, arg1)
    call GET_COMMAND_ARGUMENT(2, arg2)
    call GET_COMMAND_ARGUMENT(3, arg3)
    read(arg1, *) na
    read(arg2, *) nev
    read(arg3, *) nblk
  endif

  !-------------------------------------------------------------------------------
  !  MPI Initialization

  !#ifndef WITH_OPENMP
  call mpi_init(mpierr)
  !#else
  !   required_mpi_thread_level = MPI_THREAD_MULTIPLE
  !   call mpi_init_thread(required_mpi_thread_level,     &
  !                        provided_mpi_thread_level, mpierr)
  !
  !   if (required_mpi_thread_level .ne. provided_mpi_thread_level) then
  !      write(error_unit,*) "MPI ERROR: MPI_THREAD_MULTIPLE is not provided &
  !      & on this system"
  !      write(error_unit,*) "           ", &
  !      mpi_thread_level_name(provided_mpi_thread_level), " is available"
  !      call EXIT(1)
  !      stop 1
  !   endif
  !
  !#endif

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(1) = mpi_wtime()

  call mpi_comm_rank(mpi_comm_world, myid, mpierr)
  call mpi_comm_size(mpi_comm_world, nprocs, mpierr)

  status = 0
  !#ifdef WITH_OPENMP
  !   if (myid .eq. 0) then
  !      print *,"Threaded version of test program"
  !      print *,"Using ",omp_get_max_threads()," threads"
  !      print *," "
  !   endif
  !#endif

  ! Selection of number of processor rows/columns
  ! We try to set up the grid square-like, i.e. start the search for possible
  ! divisors of nprocs with a number next to the square root of nprocs
  ! and decrement it until a divisor is found.

  do np_cols = nint(sqrt(real(nprocs))), 2, -1
    if (mod(nprocs, np_cols) == 0) exit
  enddo
  ! at the end of the above loop, nprocs is always divisible by np_cols
  np_rows = nprocs/np_cols
  if (myid == 0) then
    print *
    print '(a)', 'Standard eigenvalue problem - REAL version'
    print *
    print '(3(a,i0))', 'Matrix size=', na, &
         ', Number of eigenvectors=', nev, ', Block size=', nblk
    print '(3(a,i0))', 'Number of processor rows=', np_rows, &
         ', cols=', np_cols, ', total=', nprocs
    print *
  endif

  !-------------------------------------------------------------------------------
  ! Set up BLACS context and MPI communicators
  !
  ! The BLACS context is only necessary for using Scalapack.
  !
  ! For ELPA, the MPI communicators along rows/cols are sufficient,
  ! and the grid setup may be done in an arbitrary way as long as it is
  ! consistent (i.e. 0<=my_prow<np_rows, 0<=my_pcol<np_cols and every
  ! process has a unique (my_prow,my_pcol) pair).

  context = mpi_comm_world
  call BLACS_Gridinit( context, 'C', np_rows, np_cols )
  call BLACS_Gridinfo( context, nprow, npcol, my_prow, my_pcol )

  if (myid == 0) then
    print '(a)','| Past BLACS_Gridinfo.'
  end if

  ! All ELPA routines need MPI communicators for communicating within
  ! rows or columns of processes, these are set in get_elpa_row_col_comms.

  call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
       mpi_comm_rows, mpi_comm_cols)

  if (myid == 0) then
    print '(a)','| Past split communicator setup for rows and columns.'
  end if

  ! Determine the necessary size of the distributed matrices,
  ! we use the Scalapack tools routine NUMROC for that.

  na_rows = numroc(na, nblk, my_prow, 0, np_rows)
  na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

  ! Set up a scalapack descriptor for the checks below.
  ! For ELPA the following restrictions hold:
  ! - block sizes in both directions must be identical (args 4+5)
  ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

  call descinit(sc_desc, na, na, nblk, nblk, 0, 0, context, na_rows, info)

  if (myid==0) then
    print '(a)','| Past scalapack descriptor setup.'
  end if

  !-------------------------------------------------------------------------------
  ! Allocate matrices and set up a test matrix for the eigenvalue problem

  allocate(a(na_rows, na_cols), b(na_rows, na_cols), z(na_rows, na_cols), &
       as(na_rows, na_cols), bs(na_rows, na_cols), ev(na), stat = ierr)
  if (ierr /= 0) then
    write(error_unit, *) 'memory allocation failed'
    call mpi_abort(mpi_comm_world, mpierr)
  end if

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(2) = mpi_wtime()

  ! For getting a symmetric test matrix A we get a random matrix Z
  ! and calculate A = Z + Z**T

  ! We want different random numbers on every process
  ! (otherways A might get rank deficient):
  iseed(:) = myid
  call RANDOM_SEED(put = iseed)

  ! Use z as a temporary matrix to generate a
  call RANDOM_NUMBER(z)
  a(:, :) = z(:, :)
  if (myid == 0) then
    print '(a)', '| Random matrix block has been set up. &
         & (only processor 0 confirms this step)'
  end if
  ! A = A + Z**T
  call pdtran(na, na, 1.d0, z, 1, 1, sc_desc, 1.d0, a, 1, 1, sc_desc)
  if (myid == 0) then
    print '(a)','| Random matrix has been symmetrized.'
  end if
  ! Save original matrix A for later accuracy checks
  as = a

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(3) = mpi_wtime()

  ! Set B.
  do i = 1, na
    if (i > 1) then
      call pdelset(b, i, i - 1, sc_desc, -1.0d0)
    end if
    call pdelset(b, i, i, sc_desc, 4.0d0)
    if (i < na) then
      call pdelset(b, i, i + 1, sc_desc, -1.0d0)
    end if
  end do
  ! Save original matrix B for later accuracy checks
  bs = b

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(4) = mpi_wtime()

  if (myid == 0) then
    print '(a)', '| cholesky_real start'
  end if
  ! Return of cholesky_real is stored in the upper triangle.
  call cholesky_real(na, b, na_rows, nblk, mpi_comm_rows, mpi_comm_cols, success)
  if (.not. success) then
    if (myid == 0) then
      print *, 'cholesky_real failed'
    end if
    stop
  end if

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(5) = mpi_wtime()

  call invert_trm_real(na, b, na_rows, nblk, mpi_comm_rows, mpi_comm_cols, success)
  ! invert_trm_real always returns fail
  !if (.not. success) then
  !  if (myid == 0) then
  !    print *, 'invert_trm_real failed'
  !  end if
  !  stop
  !end if

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(6) = mpi_wtime()

  ! Reduce A as U^-T A U^-1
  ! A <- U^-T A
  ! This operation can be done as below:
  ! call pdtrmm('Left', 'Upper', 'Trans', 'No_unit', na, na, 1.0d0, &
  !      b, 1, 1, sc_desc, a, 1, 1, sc_desc)
  ! but it is slow. Instead use mult_at_b_real.
  call mult_at_b_real('Upper', 'Full', na, na, &
       b, na_rows, as, na_rows, &
       nblk, mpi_comm_rows, mpi_comm_cols, a, na_rows)

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(7) = mpi_wtime()

  ! A <- A U^-1
  call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', na, na, 1.0d0, &
       b, 1, 1, sc_desc, a, 1, 1, sc_desc)

  !-------------------------------------------------------------------------------
  ! Calculate eigenvalues/eigenvectors

  if (myid == 0) then
    print '(a)', '| Entering two-stage ELPA solver ... '
    print *
  end if

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(8) = mpi_wtime()

  success = solve_evp_real_2stage(na, nev, a, na_rows, ev, z, &
       na_rows, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_world)
  if (.not. success) then
    write(error_unit, *) "solve_evp_real_2stage produced &
         & an error                  ! Aborting..."
    call MPI_ABORT(mpi_comm_world, mpierr)
  endif

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(9) = mpi_wtime()

  ! Z <- U^-1 Z
  call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', na, na, 1.0d0, &
       b, 1, 1, sc_desc, z, 1, 1, sc_desc)

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(10) = mpi_wtime()

  if (myid == 0) then
    print '(a)','| Two-step ELPA solver complete.'
    print *, 'init            : ', times(2) - times(1)
    print *, 'set_A           : ', times(3) - times(2)
    print *, 'set_B           : ', times(4) - times(3)
    print *, 'cholesky_real   : ', times(5) - times(4)
    print *, 'invert_trm_real : ', times(6) - times(5)
    print *, 'pdtrmm_A_left   : ', times(7) - times(6)
    print *, 'pdtrmm_A_right  : ', times(8) - times(7)
    print *, 'solve_evp       : ', times(9) - times(8)
    print *, 'pdtrmm_EVs      : ', times(10) - times(9)
    print *, 'total           : ', times(10) - times(1)
    print *
    print *, 'solve_evp_real_2stage transform_to_tridi :', time_evp_fwd
    print *, 'solve_evp_real_2stage solve_tridi        :', time_evp_solve
    print *, 'solve_evp_real_2stage transform_back_EVs :', time_evp_back
    print *, 'solve_evp_real_2stage total              :', &
         time_evp_back+time_evp_solve+time_evp_fwd
  end if

  !-------------------------------------------------------------------------------
  ! Test correctness of result (using plain scalapack routines)

  deallocate(a)
  allocate(tmp1(na_rows, na_cols))

  ! 1. Residual (maximum of || A*Zi - Zi*EVi ||)
  ! tmp1 = A * Z
  call pdgemm('N', 'N', na, nev, na,1.d0, as, 1, 1, sc_desc, &
       z, 1, 1, sc_desc, 0.d0, tmp1, 1, 1, sc_desc)
  deallocate(as)
  allocate(tmp2(na_rows, na_cols))
  ! tmp2 = B * Zi * EVi
  call pdsymm('Left', 'Upper', na, nev, 1.0d0, bs, 1, 1, sc_desc, &
       z, 1, 1, sc_desc, 0.0d0, tmp2, 1, 1, sc_desc)
  do i=1,nev
    call pdscal(na, ev(i), tmp2, 1, i, sc_desc, 1)
  enddo
  ! tmp1 = A*Zi - B*Zi*EVi
  tmp1(:, :) =  tmp1(:, :) - tmp2(:, :)

  ! Get maximum norm of columns of tmp1
  errmax = 0.0d0
  do i=1, nev
    call pdnrm2(na, err, tmp1, 1, i, sc_desc, 1)
    errmax = max(errmax, err)
  enddo

  ! Get maximum error norm over all processors
  err = errmax
  call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
  if (myid==0) then
    print *
    print *,'Error Residual     :',errmax
  end if

  if (errmax > 5e-12) then
    status = 1
  endif

  !! 2. Eigenvector orthogonality
  !
  !! tmp1 = Z**T * Z
  !tmp1 = 0
  !call pdgemm('T','N',nev,nev,na,1.d0,z,1,1,sc_desc, &
  !     z,1,1,sc_desc,0.d0,tmp1,1,1,sc_desc)
  !! Initialize tmp2 to unit matrix
  !tmp2 = 0
  !call pdlaset('A',nev,nev,0.d0,1.d0,tmp2,1,1,sc_desc)
  !
  !! tmp1 = Z**T * Z - Unit Matrix
  !tmp1(:,:) =  tmp1(:,:) - tmp2(:,:)
  !
  !! Get maximum error (max abs value in tmp1)
  !err = maxval(abs(tmp1))
  !call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
  !if(myid==0) print *,'Error Orthogonality:',errmax
  !
  !if (errmax .gt. 5e-12) then
  !  status = 1
  !endif

  call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
  times(11) = mpi_wtime()
  if (myid == 0) then
    print *, 'whole_elapse_time : ', times(11) - times(1)
  end if

  deallocate(z, tmp1, tmp2, ev)
  call blacs_gridexit(context)
  call mpi_finalize(mpierr)
  call exit(status)
end program test_real2
