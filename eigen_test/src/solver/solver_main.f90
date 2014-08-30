module solver_main
  use descriptor_parameters
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use eigenpairs_types, only: eigenpairs_types_union, eigenpairs_blacs
  use generalized_to_standard, only : reduce_generalized, recovery_generalized
  use processes, only : process, setup_distribution, print_map_of_grid_to_processes, &
       check_master, terminate
  use time, only : get_wall_clock_base_count, get_wall_clock_time
  implicit none

  private
  public :: eigen_solver

contains

  subroutine eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
    use ELPA1
    use ELPA2
    use solver_lapack, only : eigen_solver_lapack
    use solver_scalapack_all, only : eigen_solver_scalapack_all
    use solver_scalapack_select, only : eigen_solver_scalapack_select
    use solver_eigenexa, only : setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa
    include 'mpif.h'

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: n, desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), desc_A_re(desc_size), &
         myid, nblk, np_rows, np_cols, nprow, npcol, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size)
    integer :: ierr, mpierr, info
    logical :: success
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs_tmp

    integer :: numroc

    integer :: base_count
    double precision :: times(10)

    n = arg%matrix_A_info%rows

    select case (trim(arg%solver_type))
    case ('lapack')
      call eigen_solver_lapack(matrix_A, eigenpairs)
    case ('scalapack_all')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
    case ('scalapack_select')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
    case ('general_scalapack_all')
      call get_wall_clock_base_count(base_count)
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call get_wall_clock_time(base_count, times(1))
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call get_wall_clock_time(base_count, times(2))
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
      call get_wall_clock_time(base_count, times(3))
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call get_wall_clock_time(base_count, times(4))
      if (check_master()) then
        print *, 'general_scalapack_all elapsed time (printed in solver_main): '
        print *, 'general_scalapack_all setup_matrix: ', times(1)
        print *, 'general_scalapack_all reduce_generalized: ', times(2)
        print *, 'general_scalapack_all eigen_solver_scalapack: ', times(3)
        print *, 'general_scalapack_all recovery_generalized: ', times(4)
      end if
    case ('general_scalapack_select')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      if (arg%block_size <= 0) then  ! Default block size
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      else
        call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist, arg%block_size)
        call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist, arg%block_size)
      end if
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
      call recovery_generalized(n, arg%n_vec, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('eigenexa')
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      call setup_distributed_matrix_for_eigenexa(n, desc_A, matrix_A_dist, eigenpairs)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_eigenexa(matrix_A_dist, desc_A, arg%n_vec, eigenpairs)
    case ('general_eigenexa')
      call get_wall_clock_base_count(base_count)
      call setup_distribution(proc)
      if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call get_wall_clock_time(base_count, times(1))
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call get_wall_clock_time(base_count, times(2))
      call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
      deallocate(matrix_A_dist)
      call get_wall_clock_time(base_count, times(3))
      call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, arg%n_vec, eigenpairs_tmp, 'L')
      call get_wall_clock_time(base_count, times(4))
      deallocate(matrix_A_redist)
      eigenpairs%type_number = 2
      allocate(eigenpairs%blacs%values(n))
      eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
      eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
           eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
           eigenpairs_tmp%blacs%desc(context_))
      call get_wall_clock_time(base_count, times(5))
      call recovery_generalized(n, n, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
      call get_wall_clock_time(base_count, times(6))
      if (check_master()) then
        print *, 'general_eigenexa elapsed time (printed in solver_main): '
        print *, 'general_eigenexa setup_matrix: ', times(1)
        print *, 'general_eigenexa reduce_generalized: ', times(2)
        print *, 'general_eigenexa pdgemr2d_A: ', times(3)
        print *, 'general_eigenexa eigen_solver_eigenexa: ', times(4)
        print *, 'general_eigenexa pdgemr2d_B: ', times(5)
        print *, 'general_eigenexa recovery_generalized: ', times(6)
      end if
    case ('general_elpa')
      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(1) = mpi_wtime()

      nblk = 64
      call setup_distribution(proc)
      call mpi_comm_rank(mpi_comm_world, myid, mpierr)
      !context = mpi_comm_world
      !call BLACS_Gridinit( mpi_comm_world, 'C', np_rows, np_cols )
      call BLACS_Gridinfo( mpi_comm_world, np_rows, np_cols, my_prow, my_pcol )
      call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
           mpi_comm_rows, mpi_comm_cols)
      na_rows = numroc(n, nblk, my_prow, 0, np_rows)
      na_cols = numroc(n, nblk, my_pcol, 0, np_cols)
      call descinit(sc_desc, n, n, nblk, nblk, 0, 0, mpi_comm_world, na_rows, info)
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call setup_distributed_matrix('Eigenvectors', proc, n, n, &
           eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
      allocate(eigenpairs%blacs%values(n))
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(2) = mpi_wtime()

      ! Return of cholesky_real is stored in the upper triangle.
      call cholesky_real(n, matrix_B_dist, na_rows, nblk, mpi_comm_rows, mpi_comm_cols, success)
      if (.not. success) then
        if (myid == 0) then
          print *, 'cholesky_real failed'
        end if
        stop
      end if

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(3) = mpi_wtime()

      call invert_trm_real(n, matrix_B_dist, na_rows, nblk, mpi_comm_rows, mpi_comm_cols, success)
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
           nblk, mpi_comm_rows, mpi_comm_cols, matrix_A_dist, na_rows)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(5) = mpi_wtime()

      ! A <- A U^-1
      call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
           matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(6) = mpi_wtime()

      success = solve_evp_real(n, n, matrix_A_dist, na_rows, &
           eigenpairs%blacs%values, eigenpairs%blacs%Vectors, na_rows, &
           nblk, mpi_comm_rows, mpi_comm_cols)
      if (.not. success) then
        write(6, *) "solve_evp_real_2stage produced &
             & an error                  ! Aborting..."
        call mpi_abort(mpi_comm_world, mpierr)
      endif

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(7) = mpi_wtime()

      ! Z <- U^-1 Z
      call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
           matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

      eigenpairs%type_number = 2

      call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
      times(8) = mpi_wtime()

      if (myid == 0) then
        print '(a)','| ELPA solver complete.'
        print *, 'init            : ', times(2) - times(1)
        print *, 'cholesky_real   : ', times(3) - times(2)
        print *, 'invert_trm_real : ', times(4) - times(3)
        print *, 'pdtrmm_A_left   : ', times(5) - times(4)
        print *, 'pdtrmm_A_right  : ', times(6) - times(5)
        print *, 'solve_evp       : ', times(7) - times(6)
        print *, 'pdtrmm_EVs      : ', times(8) - times(7)
        print *, 'total           : ', times(8) - times(1)
        print *
        print *, 'solve_evp_real_2stage transform_to_tridi :', time_evp_fwd
        print *, 'solve_evp_real_2stage solve_tridi        :', time_evp_solve
        print *, 'solve_evp_real_2stage transform_back_EVs :', time_evp_back
        print *, 'solve_evp_real_2stage total              :', &
             time_evp_back+time_evp_solve+time_evp_fwd
      end if
    case default
      stop '[Error] lib_eigen_solver: Unknown solver'
    end select
  end subroutine eigen_solver
end module solver_main
