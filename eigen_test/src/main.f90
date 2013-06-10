program eigen_test
!
  use solver_main, only : lib_eigen_solver, lib_eigen_checker !(routine)
  use command_argument, only : read_command_argument !(routine)
  use read_matrix, only : read_matrix_file !(rouine)
  use distribute_matrix, only : create_dense_matrix !(routine)
  use time, only : data_and_time_wrapper !(routine)
  use processes, only : check_master !(routine)
!
  implicit none
  integer :: verbose_level
!
  real(kind(1.d0)), allocatable :: mat(:,:)
  real(kind(1.d0)), allocatable :: mat_bak(:,:)
  real(kind(1.d0)), allocatable :: eig(:)
  real(kind(1.d0)) :: rn_ave, rn_max
!
  integer :: matrix_size, num_non_zeros
!
  integer :: n_vec = -1, n_check_vec = -1
  integer :: j, iunit
  integer, parameter :: num_filename = 2
  character(len=256) :: matrix_filename(num_filename)
!
  character(len=256) :: eqn_type
  character(len=256) :: matrix_type
  character(len=256) :: solver_type
!
  character(len=*), parameter :: filename='output_eigen_value.txt'
!
  real(kind(1.d0)),    allocatable :: mat_value(:,:)
  integer,             allocatable :: mat_suffix(:,:)
!
  character(len=256) :: chara_data_time
!
  logical :: is_master
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  call read_command_argument(verbose_level, matrix_filename, eqn_type, matrix_type, solver_type, n_vec, n_check_vec)

  call read_matrix_file(verbose_level, matrix_filename(1), matrix_type, matrix_size, num_non_zeros, mat_value, mat_suffix)

  if (n_vec == -1) then ! unspecified in command line arguments
    n_vec = matrix_size
  end if

  if (n_check_vec == -1) then
    n_check_vec = n_vec
  end if

  call check_master(solver_type, is_master)

  if (is_master) then
     print *, '--------------------------------------'
     print *, 'Eigen test'
     print *, '--------------------------------------'
     call data_and_time_wrapper(chara_data_time)
     print *, '  ', trim(chara_data_time)
     print *, '  verbose level        = ', verbose_level
     print *, '  equation type        = ', trim(eqn_type)
     print *, '  matrix type          = ', trim(matrix_type)
     print *, '  solver type          = ', trim(solver_type)
     print *, '  matrix filename      = ', trim(matrix_filename(1))
     print *, '  required eigenvalues = ', n_vec
     print *, '  checked eigenvalues  = ', n_check_vec
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!
!
  call create_dense_matrix(verbose_level, matrix_size, num_non_zeros, mat_value, mat_suffix, mat, mat_bak, eig)
!
  call lib_eigen_solver(mat, eig, solver_type, n_vec)

  if (.not. is_master) then
     stop
  end if

!
  if (n_check_vec /= 0) then
    call lib_eigen_checker(mat_bak, mat, eig, n_vec, n_check_vec, rn_ave, rn_max)
    write(*,'(a,2e20.8)')' check residual norm: ave, max = ', rn_ave, rn_max
  endif
!
  write(*,*)'Plot the eigen values into the file :',trim(filename)
  iunit=70
  open (iunit, file=filename, status='unknown')
  do j=1,n_vec
    write(iunit,'(i10,f20.12)') j, eig(j)
  enddo
!
  write(*,'(a)') '...the program ends'
  write(*,'(a)') '--------------------------------------'
!
end program eigen_test
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
