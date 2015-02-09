program eigen_test
  use mpi
  use fson
  use fson_value_m
  use fson_string_m
  use distribute_matrix
  use event_logger_m
  use solver_main, only : eigen_solver
  use command_argument, only : argument, required_memory, &
       read_command_argument, print_command_argument, fson_setting_add
  use matrix_io, only : sparse_mat, read_matrix_file, print_eigenvectors
  use processes, only : get_num_procs, check_master
  use eigenpairs_types, only : eigenpairs_types_union
  use verifier, only : eval_residual_norm, eval_orthogonality
  implicit none

  type(argument) :: arg
  type(sparse_mat) :: matrix_A, matrix_B
  type(eigenpairs_types_union) :: eigenpairs
  double precision :: A_norm, rn_ave, rn_max, orthogonality
  integer :: num_mpi_procs, num_omp_procs, j, ierr
  integer, parameter :: iunit = 10
  double precision :: time_start, time_start_part, time_end
  type(fson_value), pointer :: output

  call mpi_init(ierr)
  if (ierr /= 0) then
    write (0, *) '[Error] eigen_test: mpi_init failed, error code is ', ierr
    stop
  end if

  time_start = mpi_wtime()
  time_start_part = time_start

  call read_command_argument(arg)

  if (check_master()) then
    print '("---------- Eigen Test start ----------")'
    print '("----- Configurations -----")'
    call print_command_argument(arg)
    print '("approximate required memory per process (Mbytes): ", f10.1)', &
         required_memory(arg) / real(2 ** 20)
    call get_num_procs(num_mpi_procs, num_omp_procs)
    print '("MPI processes: ", i0)', num_mpi_procs
    print '("OpenMP threads per process (may be inaccurate): ", i0)', num_omp_procs
  end if

  time_end = mpi_wtime()
  call add_event('main:read_command_argument', time_end - time_start_part)
  time_start_part = time_end

  if (check_master()) then
    call read_matrix_file(arg%matrix_A_filename, arg%matrix_A_info, matrix_A)
    if (arg%is_generalized_problem) then
      call read_matrix_file(arg%matrix_B_filename, arg%matrix_B_info, matrix_B)
    end if
    output => fson_value_create()
    output%value_type = TYPE_OBJECT
    call fson_setting_add(arg, output)
  end if

  if (arg%is_dry_run) then
    if (check_master()) print '(/, "dry run mode, exit")'
    call mpi_finalize(ierr)
    stop
  end if

  call bcast_sparse_matrix(0, arg%matrix_A_info, matrix_A)
  if (arg%is_generalized_problem) then
    call bcast_sparse_matrix(0, arg%matrix_B_info, matrix_B)
  end if

  time_end = mpi_wtime()
  call add_event('main:read_matrix_file', time_end - time_start_part)
  time_start_part = time_end

  if (check_master()) print '(/, "----- Solver Call -----")'
  if (arg%is_generalized_problem) then
    call eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
  else
    call eigen_solver(arg, matrix_A, eigenpairs)
  end if

  time_end = mpi_wtime()
  call add_event('main:eigen_solver', time_end - time_start_part)
  time_start_part = time_end

  ! Print eigenvalues and eigenvectors if required
  if (check_master()) then
    open(iunit, file=arg%output_filename, status='unknown')
    do j=1,arg%n_vec
      if (eigenpairs%type_number == 1) then
        write (iunit, '(i10, f20.12)') j, eigenpairs%local%values(j)
      else if (eigenpairs%type_number == 2) then
        write (iunit, '(i10, f20.12)') j, eigenpairs%blacs%values(j)
      end if
    enddo
    close(iunit)
  end if

  if (arg%printed_vecs_start /= 0) then
    call print_eigenvectors(arg, eigenpairs)
  end if

  time_end = mpi_wtime()
  call add_event('main:print_eigenpairs', time_end - time_start_part)
  time_start_part = time_end

  if (arg%n_check_vec /= 0) then
    if (check_master()) print '(/, "----- Checker Call -----")'
    if (arg%is_generalized_problem) then
      call eval_residual_norm(arg, matrix_A, eigenpairs, &
           A_norm, rn_ave, rn_max, matrix_B)
    else
      call eval_residual_norm(arg, matrix_A, eigenpairs, &
           A_norm, rn_ave, rn_max)
    end if

    if (check_master()) then
      print '("A norm: ", e15.8)', A_norm
      print '("residual norm (average): ", e15.8)', rn_ave
      print '("residual norm (max):     ", e15.8)', rn_max
    end if
  end if

  time_end = mpi_wtime()
  call add_event('main:eval_residual_norm', time_end - time_start_part)
  time_start_part = time_end

  if (arg%ortho_check_index_start /= 0) then
    if (arg%is_generalized_problem) then
      call eval_orthogonality(arg, eigenpairs, orthogonality, matrix_B)
    else
      call eval_orthogonality(arg, eigenpairs, orthogonality)
    end if
    if (check_master()) then
      print '("orthogonality criterion: ", e15.8)', orthogonality
    end if
  end if

  time_end = mpi_wtime()
  call add_event('main:eval_orthogonality', time_end - time_start_part)
  call add_event('main:total', time_end - time_start)

  if (check_master()) then
    call fson_events_add(output)
    open(iunit, file=trim(arg%log_filename), status='unknown')
    call fson_print(iunit, output)
    close(iunit)
  end if

  call mpi_finalize(ierr)
end program eigen_test
