program eigen_test
  use solver_main, only : eigen_solver, eigen_checker
  use command_argument, only : argument, required_memory, &
       read_command_argument, print_command_argument
  use matrix_io, only : sparse_mat, read_matrix_file, print_eigenvectors
  use distribute_matrix, only : create_dense_matrix !(routine)
  use time, only : data_and_time_wrapper !(routine)
  use processes, only : get_num_procs, check_master
  use eigenpairs_types, only : eigenpairs_types_union
  implicit none

  include 'mpif.h'

  type(argument) :: arg
  type(sparse_mat) :: matrix_A, matrix_B
  type(eigenpairs_types_union) :: eigenpairs
  double precision :: rn_ave, rn_max
  integer :: num_mpi_procs, num_omp_procs, j, ierr
  integer, parameter :: iunit = 10
  logical :: is_master

  call mpi_init(ierr)
  if (ierr /= 0) then
    write (0, *) '[Error] eigen_test: mpi_init failed, error code is ', ierr
    stop
  end if

  call read_command_argument(arg)

  is_master = check_master()

  if (is_master) then
    print '("----- Configurations -----")'
    call print_command_argument(arg)
    print '("approximate required memory per process (Mbytes): ", f10.1)', &
         required_memory(arg) / real(2 ** 20)
    call get_num_procs(num_mpi_procs, num_omp_procs)
    print '("MPI processes: ", i0)', num_mpi_procs
    print '("OpenMP threads per process: ", i0)', num_omp_procs
  end if

  call read_matrix_file(arg%matrix_A_filename, arg%matrix_A_info, matrix_A)
  if (arg%is_generalized_problem) then
    call read_matrix_file(arg%matrix_B_filename, arg%matrix_B_info, matrix_B)
  end if

  if (is_master) then
    print '(/, "----- Solver Call -----")'
  end if

  if (arg%is_generalized_problem) then
    call eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
  else
    call eigen_solver(arg, matrix_A, eigenpairs)
  end if

  if (arg%printed_vecs_start /= 0) then
    call print_eigenvectors(arg, eigenpairs)
  end if

  if (arg%n_check_vec /= 0) then
    if (is_master) then
      print '(/, "----- Checker Call -----")'
    end if
    if (arg%is_generalized_problem) then
      call eigen_checker(arg, matrix_A, eigenpairs, &
           rn_ave, rn_max, matrix_B)
    else
      call eigen_checker(arg, matrix_A, eigenpairs, &
           rn_ave, rn_max)
    end if

    if (is_master) then
      print '("residual norm (average): ", e15.8)', rn_ave
      print '("residual norm (max):     ", e15.8)', rn_max
    end if
  endif

  call mpi_finalize(ierr)
  if (.not. is_master) then
    stop
  end if

  open(iunit, file=arg%output_filename, status='unknown')
  do j=1,arg%n_vec
    if (eigenpairs%type_number == 1) then
      write (iunit, '(i10, f20.12)') j, eigenpairs%local%values(j)
    else if (eigenpairs%type_number == 2) then
      write (iunit, '(i10, f20.12)') j, eigenpairs%blacs%values(j)
    end if
  enddo
  close(iunit)
end program eigen_test
