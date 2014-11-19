module solver_elpa_eigenexa
  use eigen_libs
  use ELPA1
  use ELPA2

  use global_variables, only : g_block_size
  use distribute_matrix, only : process, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use matrix_io, only : sparse_mat
  use processes, only : check_master, setup_distribution, terminate
  use solver_eigenexa
  use time, only : get_wall_clock_base_count, get_wall_clock_time

  implicit none

  !private

contains

  subroutine solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs
    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         desc_A_re(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success
    type(eigenpairs_types_union) :: eigenpairs_tmp

    double precision :: times(10)
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    integer :: numroc
    include 'mpif.h'

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(1) = mpi_wtime()
      call mpi_comm_rank(mpi_comm_world, myid, mpierr)
      call BLACS_Gridinfo( mpi_comm_world, np_rows, np_cols, my_prow, my_pcol )
      call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
           mpi_comm_rows, mpi_comm_cols)

      max_block_size = min(n / np_rows, n / np_cols)
      block_size = min(max_block_size, g_block_size)
      na_rows = numroc(n, block_size, my_prow, 0, np_rows)
      na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
      call descinit(sc_desc, n, n, block_size, block_size, 0, 0, mpi_comm_world, na_rows, info)
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(2) = mpi_wtime()

      ! Return of cholesky_real is stored in the upper triangle.
      call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
      if (.not. success) then
        call terminate('solver_main, general_elpa_eigenexa: cholesky_real failed', 1)
      end if

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(3) = mpi_wtime()

      call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
      ! invert_trm_real always returns fail
      !if (.not. success) then
      !  if (myid == 0) then
      !    print *, 'invert_trm_real failed'
      !  end if
      !  stop
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
      deallocate(matrix_A2_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(5) = mpi_wtime()

      ! A <- A U^-1
      call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
           matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(6) = mpi_wtime()

      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      deallocate(matrix_A_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(7) = mpi_wtime()

      call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(8) = mpi_wtime()

      deallocate(matrix_A_redist)
      eigenpairs%type_number = 2
      allocate(eigenpairs%blacs%values(n), stat = ierr)
      if (ierr /= 0) then
        call terminate('eigen_solver, general_elpa_eigenexa: allocation failed', ierr)
      end if
      eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
      eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
           eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
           eigenpairs_tmp%blacs%desc(context_))

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(9) = mpi_wtime()

      ! Z <- U^-1 Z
      call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
           matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(10) = mpi_wtime()

      if (check_master()) then
        print *, 'init                  : ', times(2) - times(1)
        print *, 'cholesky_real         : ', times(3) - times(2)
        print *, 'invert_trm_real       : ', times(4) - times(3)
        print *, 'pdtrmm_A_left         : ', times(5) - times(4)
        print *, 'pdtrmm_A_right        : ', times(6) - times(5)
        print *, 'pdgemr2d_A            : ', times(7) - times(6)
        print *, 'eigen_solver_eigenexa : ', times(8) - times(7)
        print *, 'pdgemr2d_B            : ', times(9) - times(8)
        print *, 'pdtrmm_EVs            : ', times(10) - times(9)
        print *, 'total                 : ', times(10) - times(1)
      end if
    end subroutine solve_with_general_elpa_eigenexa


    subroutine solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs
    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         desc_A_re(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success
    type(eigenpairs_types_union) :: eigenpairs_tmp
    double precision :: times(10)
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    integer :: numroc
    include 'mpif.h'

          call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(1) = mpi_wtime()
      call mpi_comm_rank(mpi_comm_world, myid, mpierr)
      call BLACS_Gridinfo( mpi_comm_world, np_rows, np_cols, my_prow, my_pcol )
      call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
           mpi_comm_rows, mpi_comm_cols)

      max_block_size = min(n / np_rows, n / np_cols)
      block_size = min(max_block_size, g_block_size)
      na_rows = numroc(n, block_size, my_prow, 0, np_rows)
      na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
      call descinit(sc_desc, n, n, block_size, block_size, 0, 0, mpi_comm_world, na_rows, info)
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(2) = mpi_wtime()

      ! Return of cholesky_real is stored in the upper triangle.
      call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
      if (.not. success) then
        call terminate('solver_main, general_elpa_eigenk: cholesky_real failed', 1)
      end if

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(3) = mpi_wtime()

      call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
      ! invert_trm_real always returns fail
      !if (.not. success) then
      !  if (myid == 0) then
      !    print *, 'invert_trm_real failed'
      !  end if
      !  stop
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
      deallocate(matrix_A2_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(5) = mpi_wtime()

      ! A <- A U^-1
      call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
           matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(6) = mpi_wtime()

      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      deallocate(matrix_A_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(7) = mpi_wtime()

      call eigen_solver_eigenk(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(8) = mpi_wtime()

      deallocate(matrix_A_redist)
      eigenpairs%type_number = 2
      allocate(eigenpairs%blacs%values(n), stat = ierr)
      if (ierr /= 0) then
        call terminate('eigen_solver, general_elpa_eigenk: allocation failed', ierr)
      end if
      eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
      eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
           eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
           eigenpairs_tmp%blacs%desc(context_))

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(9) = mpi_wtime()

      ! Z <- U^-1 Z
      call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
           matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(10) = mpi_wtime()

      if (check_master()) then
        print *, 'init                : ', times(2) - times(1)
        print *, 'cholesky_real       : ', times(3) - times(2)
        print *, 'invert_trm_real     : ', times(4) - times(3)
        print *, 'pdtrmm_A_left       : ', times(5) - times(4)
        print *, 'pdtrmm_A_right      : ', times(6) - times(5)
        print *, 'pdgemr2d_A          : ', times(7) - times(6)
        print *, 'eigen_solver_eigenk : ', times(8) - times(7)
        print *, 'pdgemr2d_B          : ', times(9) - times(8)
        print *, 'pdtrmm_EVs          : ', times(10) - times(9)
        print *, 'total               : ', times(10) - times(1)
      end if
    end subroutine solve_with_general_elpa_eigenk
end module solver_elpa_eigenexa
