program eigen_test
  use solver_main, only : lib_eigen_solver, lib_eigen_checker !(routine)
  use command_argument, only : argument, &
       read_command_argument, print_command_argument
  use matrix_io, only : sparse_mat, read_matrix_file !(type, rouine)
  use distribute_matrix, only : create_dense_matrix !(routine)
  use time, only : data_and_time_wrapper !(routine)
  use processes, only : check_master !(routine)
  implicit none

  type(argument) :: arg

  type(sparse_mat) :: matrix_A, matrix_B
  real(kind(1.d0)), allocatable :: mat(:,:)
  real(kind(1.d0)), allocatable :: eigenvalues(:), eigenvectors(:, :)
  real(kind(1.d0)) :: rn_ave, rn_max

  integer :: j, iunit

  real(kind(1.d0)), allocatable :: mat_value(:,:)
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
    call lib_eigen_solver(arg, matrix_A, eigenvalues, eigenvectors, matrix_B)
  else
    call lib_eigen_solver(arg, matrix_A, eigenvalues, eigenvectors)
  end if

  if (.not. is_master) then
     stop
  end if

  if (arg%n_check_vec /= 0) then
    !call create_dense_matrix(verbose_level, mat_in, mat) todo: integrate into below
    if (arg%is_generalized_problem) then
      !call lib_eigen_checker(arg, matrix_A, eigenvalues, eigenvectors, &
      !     rn_ave, rn_max, matrix_B)
    else
      !call lib_eigen_checker(arg, matrix_A, eigenvalues, eigenvectors, &
      !     rn_ave, rn_max)
    end if
    write(*,'(a,2e20.8)')' check residual norm: ave, max = ', rn_ave, rn_max
  endif

  iunit=70
  open (iunit, file=arg%output_filename, status='unknown')
  do j=1,arg%n_vec
    write(iunit,'(i10,f20.12)') j, eigenvalues(j)
  enddo

  write(*,'(a)') '...the program ends'
  write(*,'(a)') '--------------------------------------'
end program eigen_test
