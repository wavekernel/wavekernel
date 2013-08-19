program eigen_test
  use solver_main, only : lib_eigen_solver, lib_eigen_checker !(routine)
  use command_argument, only : argument, &
       read_command_argument, print_command_argument
  use matrix_io, only : sparse_mat, read_matrix_file !(type, rouine)
  use distribute_matrix, only : create_dense_matrix !(routine)
  use time, only : data_and_time_wrapper !(routine)
  use processes, only : check_master !(routine)
  use eigenpairs_types, only : eigenpairs_types_union
  implicit none

  type(argument) :: arg

  type(sparse_mat) :: matrix_A, matrix_B
  double precision, allocatable :: mat(:,:)
  type(eigenpairs_types_union) :: eigenpairs
  double precision :: rn_ave, rn_max

  integer :: j, iunit, ierr

  double precision, allocatable :: mat_value(:,:)
  integer, allocatable :: mat_suffix(:,:)

  character(len=256) :: chara_data_time

  logical :: is_master

  call read_command_argument(arg)

  call read_matrix_file(arg%matrix_A_filename, arg%matrix_A_info, matrix_A)
  if (arg%is_generalized_problem) then
    call read_matrix_file(arg%matrix_B_filename, arg%matrix_B_info, matrix_B)
  end if

  call check_master(arg%solver_type, is_master)

  if (is_master) then
    call print_command_argument(arg)
  end if

  if (arg%is_generalized_problem) then
    call lib_eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
  else
    call lib_eigen_solver(arg, matrix_A, eigenpairs)
  end if

  if (arg%n_check_vec /= 0) then
    if (arg%is_generalized_problem) then
      call lib_eigen_checker(arg, matrix_A, eigenpairs, &
           rn_ave, rn_max, matrix_B)
    else
      call lib_eigen_checker(arg, matrix_A, eigenpairs, &
           rn_ave, rn_max)
    end if

    if (is_master) then
      write(*,'(a,2e20.8)')' check residual norm: ave, max = ', rn_ave, rn_max
    end if
  endif

  call mpi_finalize(ierr)
  if (.not. is_master) then
    stop
  end if

  iunit=70
  open (iunit, file=arg%output_filename, status='unknown')
  do j=1,arg%n_vec
    if (eigenpairs%type_number == 1) then
      write (iunit, '(i10, f20.12)') j, eigenpairs%local%values(j)
    else if (eigenpairs%type_number == 2) then
      write (iunit, '(i10, f20.12)') j, eigenpairs%blacs%values(j)
    end if
  enddo

  write(*,'(a)') '...the program ends'
  write(*,'(a)') '--------------------------------------'
end program eigen_test
