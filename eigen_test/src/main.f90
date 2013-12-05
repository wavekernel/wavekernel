program eigen_test
  use solver_main, only : eigen_solver
  use command_argument, only : argument, required_memory, &
       read_command_argument, print_command_argument
  use matrix_io, only : sparse_mat, read_matrix_file, print_eigenvectors
  use time, only : get_wclock_time, data_and_time_wrapper
  use processes, only : get_num_procs, check_master
  use eigenpairs_types, only : eigenpairs_types_union
  use verifier, only : eigen_checker
  implicit none

  include 'mpif.h'

  type(argument) :: arg
  type(sparse_mat) :: matrix_A, matrix_B
  type(eigenpairs_types_union) :: eigenpairs
  double precision :: rn_ave, rn_max
  integer :: num_mpi_procs, num_omp_procs, j, ierr
  integer, parameter :: iunit = 10
  logical :: is_master
  double precision :: t_start, t_end

  call get_wclock_time(t_start)

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
    print '("OpenMP threads per process (may be inaccurate): ", i0)', num_omp_procs
  end if

  call read_matrix_file(arg%matrix_A_filename, arg%matrix_A_info, matrix_A)
  if (arg%is_generalized_problem) then
    call read_matrix_file(arg%matrix_B_filename, arg%matrix_B_info, matrix_B)
  end if

  if (is_master) print '(/, "----- Solver Call -----")'

  if (arg%is_generalized_problem) then
    call eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
  else
    call eigen_solver(arg, matrix_A, eigenpairs)
  end if

  ! Print eigenvalues and eigenvectors if required
  if (is_master) then
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

  if (arg%n_check_vec /= 0) then
    if (is_master) print '(/, "----- Checker Call -----")'
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

  if (is_master) then
    call get_wclock_time(t_end, t_start)
    print '("whole execution time (sec): ", f12.2)', t_end
  end if

  call mpi_finalize(ierr)
end program eigen_test
