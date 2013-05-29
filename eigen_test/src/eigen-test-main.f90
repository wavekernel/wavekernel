program eigen_test
!
  use M_lib_eigen_routines,   only : lib_eigen_solver, lib_eigen_checker !(routine)
  use M_command_argument,     only : read_command_argument !(routine)
  use M_read_matrix_file,     only : read_matrix_file !(rouine)
  use M_create_dense_matrix,  only : create_dense_matrix !(routine)
  use M_wrap_data_and_time,   only : data_and_time_wrapper !(routine)
  use M_check_master,         only : check_master !(routine)
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
  integer :: n_vec
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
  call read_command_argument(verbose_level, matrix_filename, eqn_type, matrix_type, solver_type)

  call check_master(solver_type, is_master)

  if (is_master) then
     write(*,'(a)') '--------------------------------------'
     write(*,'(a)') 'Eigen test'
     write(*,'(a)') '--------------------------------------'
     call data_and_time_wrapper(chara_data_time)
     write(*,*) '  ',trim(chara_data_time)

     write(*,*)'  verbose level    = ',verbose_level
     write(*,*)'  equation type    = ',trim(eqn_type)
     write(*,*)'  matrix type      = ',trim(matrix_type)
     write(*,*)'  solver type      = ',trim(solver_type)
     write(*,*)'  matrix filename  = ',trim(matrix_filename(1))
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!
  call read_matrix_file(verbose_level, matrix_filename(1), matrix_type, matrix_size, num_non_zeros, mat_value, mat_suffix)
!
  call create_dense_matrix(verbose_level, matrix_size, num_non_zeros, mat_value, mat_suffix, mat, mat_bak, eig)
!
  n_vec=matrix_size
!
  call lib_eigen_solver(mat, eig, solver_type, n_vec)
!
  if (.not. is_master) then
     stop
  end if
!
  if (n_vec /=0) then
    call lib_eigen_checker(mat_bak, mat, eig, n_vec, rn_ave, rn_max)
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
