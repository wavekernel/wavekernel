module solver_elpa
  use ELPA1
  use ELPA2

  use global_variables, only : g_block_size
  use distribute_matrix, only : process, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use matrix_io, only : sparse_mat
  use processes, only : check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time

  implicit none
  private
  public :: solve_with_general_elpa_scalapack, solve_with_general_elpa1, solve_with_general_elpa2

contains

  subroutine solve_with_general_elpa_scalapack(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         pdsyevd_lwork, pdsyevd_liwork, pdsyevd_trilwmin, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success
    integer, allocatable :: pdsyevd_iwork(:)
    double precision :: times(10)
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), pdsyevd_work(:)
    integer :: numroc
    include 'mpif.h'

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(1) = mpi_wtime()
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)
    call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)

    max_block_size = min(n / np_rows, n / np_cols)
    block_size = min(max_block_size, g_block_size)
    na_rows = numroc(n, block_size, my_prow, 0, np_rows)
    na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
    call descinit(sc_desc, n, n, block_size, block_size, 0, 0, proc%context, na_rows, info)
    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa_scalapack: allocation failed', ierr)
    end if
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(2) = mpi_wtime()

    ! Return of cholesky_real is stored in the upper triangle.
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
    if (.not. success) then
      call terminate('solver_main, general_elpa1: cholesky_real failed', ierr)
    end if

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(3) = mpi_wtime()

    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
    ! invert_trm_real always returns fail
    !if (.not. success) then
    !  call terminate('solver_main, general_elpa1: invert_trm_real failed', 1)
    !end if

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(4) = mpi_wtime()

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    ! This operation can be done as below:
    ! call pdtrmm('Left', 'Upper', 'Trans', 'No_unit', na, na, 1.0d0, &
    !      b, 1, 1, sc_desc, a, 1, 1, sc_desc)
    ! but it is slow. Instead use mult_at_b_real.
    call mult_at_b_real('Upper', 'Full', n, n, &
         matrix_B_dist, na_rows, matrix_A2_dist, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols, matrix_A_dist, na_rows)

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(5) = mpi_wtime()

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(6) = mpi_wtime()

    !success = solve_evp_real(n, n, matrix_A_dist, na_rows, &
    !     eigenpairs%blacs%values, eigenpairs%blacs%Vectors, na_rows, &
    !     block_size, mpi_comm_rows, mpi_comm_cols)
    pdsyevd_trilwmin = 3 * n + max(block_size * (na_rows + 1), 3 * block_size)
    pdsyevd_lwork = max(1 + 6 * n + 2 * na_rows * na_cols, pdsyevd_trilwmin) + 2 * n
    pdsyevd_liwork = 7 * n + 8 * np_cols + 2
    allocate(pdsyevd_work(pdsyevd_lwork), pdsyevd_iwork(pdsyevd_liwork), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa_scalapack: allocation failed', ierr)
    end if
    call pdsyevd('V', 'Upper', n, matrix_A_dist, 1, 1, desc_A, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         pdsyevd_work, pdsyevd_lwork, pdsyevd_iwork, pdsyevd_liwork, info)
    if (info /= 0) then
      call terminate('solver_main, general_elpa_scalapack: pdsyevd failed', 1)
    endif

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(7) = mpi_wtime()

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    eigenpairs%type_number = 2

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(8) = mpi_wtime()

    if (check_master()) then
      print *, 'solve_with_general_elpa_scalapack init            : ', times(2) - times(1)
      print *, 'solve_with_general_elpa_scalapack cholesky_real   : ', times(3) - times(2)
      print *, 'solve_with_general_elpa_scalapack invert_trm_real : ', times(4) - times(3)
      print *, 'solve_with_general_elpa_scalapack pdtrmm_A_left   : ', times(5) - times(4)
      print *, 'solve_with_general_elpa_scalapack pdtrmm_A_right  : ', times(6) - times(5)
      print *, 'solve_with_general_elpa_scalapack pdsyevd         : ', times(7) - times(6)
      print *, 'solve_with_general_elpa_scalapack pdtrmm_EVs      : ', times(8) - times(7)
      print *, 'solve_with_general_elpa_scalapack total           : ', times(8) - times(1)
    end if
  end subroutine solve_with_general_elpa_scalapack

  subroutine solve_with_general_elpa1(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success

    double precision :: times(10)
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :)
    integer :: numroc
    include 'mpif.h'

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(1) = mpi_wtime()
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)
    call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)

    max_block_size = min(n / np_rows, n / np_cols)
    block_size = min(max_block_size, g_block_size)
    na_rows = numroc(n, block_size, my_prow, 0, np_rows)
    na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
    call descinit(sc_desc, n, n, block_size, block_size, 0, 0, proc%context, na_rows, info)
    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa1: allocation failed', ierr)
    end if
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(2) = mpi_wtime()

    ! Return of cholesky_real is stored in the upper triangle.
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
    if (.not. success) then
      call terminate('solver_main, general_elpa1: cholesky_real failed', ierr)
    end if

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(3) = mpi_wtime()

    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
    ! invert_trm_real always returns fail
    !if (.not. success) then
    !  call terminate('solver_main, general_elpa1: invert_trm_real failed', 1)
    !end if

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(4) = mpi_wtime()

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    ! This operation can be done as below:
    ! call pdtrmm('Left', 'Upper', 'Trans', 'No_unit', na, na, 1.0d0, &
    !      b, 1, 1, sc_desc, a, 1, 1, sc_desc)
    ! but it is slow. Instead use mult_at_b_real.
    call mult_at_b_real('Upper', 'Full', n, n, &
         matrix_B_dist, na_rows, matrix_A2_dist, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols, matrix_A_dist, na_rows)

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(5) = mpi_wtime()

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(6) = mpi_wtime()

    success = solve_evp_real(n, n, matrix_A_dist, na_rows, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols)
    if (.not. success) then
      call terminate('solver_main, general_elpa1: solve_evp_real failed', 1)
    endif

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(7) = mpi_wtime()

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    eigenpairs%type_number = 2

    call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
    times(8) = mpi_wtime()

    if (check_master()) then
      print *, 'solve_with_general_elpa1 init            : ', times(2) - times(1)
      print *, 'solve_with_general_elpa1 cholesky_real   : ', times(3) - times(2)
      print *, 'solve_with_general_elpa1 invert_trm_real : ', times(4) - times(3)
      print *, 'solve_with_general_elpa1 pdtrmm_A_left   : ', times(5) - times(4)
      print *, 'solve_with_general_elpa1 pdtrmm_A_right  : ', times(6) - times(5)
      print *, 'solve_with_general_elpa1 solve_evp_real  : ', times(7) - times(6)
      print *, 'solve_with_general_elpa1 pdtrmm_EVs      : ', times(8) - times(7)
      print *, 'solve_with_general_elpa1 total           : ', times(8) - times(1)
      print *
      print *, 'solve_evp_real transform_to_tridi :', time_evp_fwd
      print *, 'solve_evp_real solve_tridi        :', time_evp_solve
      print *, 'solve_evp_real transform_back_EVs :', time_evp_back
      print *, 'solve_evp_real total              :', &
           time_evp_fwd + time_evp_solve + time_evp_back
    end if
  end subroutine solve_with_general_elpa1


  subroutine solve_with_general_elpa2(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs
    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success
    double precision :: times(10)
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :)
    integer :: numroc
    include 'mpif.h'

    call mpi_barrier(mpi_comm_world, mpierr)
    times(1) = mpi_wtime()
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)
    call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)

    max_block_size = min(n / np_rows, n / np_cols)
    block_size = min(max_block_size, g_block_size)
    na_rows = numroc(n, block_size, my_prow, 0, np_rows)
    na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
    call descinit(sc_desc, n, n, block_size, block_size, 0, 0, proc%context, na_rows, info)
    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa2: allocation failed', ierr)
    end if
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    call mpi_barrier(mpi_comm_world, mpierr)
    times(2) = mpi_wtime()

    ! Return of cholesky_real is stored in the upper triangle.
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
    if (.not. success) then
      call terminate('solver_main, general_elpa2: cholesky_real failed', ierr)
    end if

    call mpi_barrier(mpi_comm_world, mpierr)
    times(3) = mpi_wtime()

    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
    ! invert_trm_real always returns fail
    !if (.not. success) then
    !  call terminate('solver_main, general_elpa2: invert_trm_real failed', 1)
    !end if

    call mpi_barrier(mpi_comm_world, mpierr)
    times(4) = mpi_wtime()

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    ! This operation can be done as below:
    ! call pdtrmm('Left', 'Upper', 'Trans', 'No_unit', na, na, 1.0d0, &
    !      b, 1, 1, sc_desc, a, 1, 1, sc_desc)
    ! but it is slow. Instead use mult_at_b_real.
    call mult_at_b_real('Upper', 'Full', n, n, &
         matrix_B_dist, na_rows, matrix_A2_dist, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols, matrix_A_dist, na_rows)

    call mpi_barrier(mpi_comm_world, mpierr)
    times(5) = mpi_wtime()

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    call mpi_barrier(mpi_comm_world, mpierr)
    times(6) = mpi_wtime()

    success = solve_evp_real_2stage(n, n, matrix_A_dist, na_rows, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols, mpi_comm_world)
    if (.not. success) then
      call terminate('solver_main, general_elpa2: solve_evp_real failed', 1)
    endif

    call mpi_barrier(mpi_comm_world, mpierr)
    times(7) = mpi_wtime()

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    eigenpairs%type_number = 2

    call mpi_barrier(mpi_comm_world, mpierr)
    times(8) = mpi_wtime()

    if (check_master()) then
      print *, 'solve_with_general_elpa2 init                  : ', times(2) - times(1)
      print *, 'solve_with_general_elpa2 cholesky_real         : ', times(3) - times(2)
      print *, 'solve_with_general_elpa2 invert_trm_real       : ', times(4) - times(3)
      print *, 'solve_with_general_elpa2 pdtrmm_A_left         : ', times(5) - times(4)
      print *, 'solve_with_general_elpa2 pdtrmm_A_right        : ', times(6) - times(5)
      print *, 'solve_with_general_elpa2 solve_evp_real_2stage : ', times(7) - times(6)
      print *, 'solve_with_general_elpa2 pdtrmm_EVs            : ', times(8) - times(7)
      print *, 'solve_with_general_elpa2 total                 : ', times(8) - times(1)
      print *
      print *, 'solve_evp_real_2stage transform_to_tridi :', time_evp_fwd
      print *, 'solve_evp_real_2stage solve_tridi        :', time_evp_solve
      print *, 'solve_evp_real_2stage transform_back_EVs :', time_evp_back
      print *, 'solve_evp_real_2stage total              :', &
           time_evp_fwd + time_evp_solve + time_evp_back
    end if
  end subroutine solve_with_general_elpa2
end module solver_elpa
