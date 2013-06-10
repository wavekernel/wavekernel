module solver_main
  use distribute_matrix, only : create_dense_matrix !(routine)

  implicit none

  private
  public :: lib_eigen_solver
  public :: lib_eigen_checker

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lib_eigen_checker(a, n_vec, n_check_vec, eigen_level, v, rn_ave, rn_max)
    implicit none
    real(kind(1.d0)),  intent(in) :: a(:,:)
    real(kind(1.d0)),  intent(in) :: v(:,:)
    real(kind(1.d0)),  intent(in) :: eigen_level(:)
    integer,           intent(in) :: n_vec, n_check_vec
    real(kind(1.d0)), intent(out) :: rn_ave, rn_max      ! residual norm average, max

    real(kind(1.d0)), allocatable :: r(:)   ! residual vector (work array)
    real(kind(1.d0)), allocatable :: rn(:)  ! residual norm

    integer :: n, ierr, j

    n=size(a,1)

    if ((n_check_vec < 1) .or. (n_check_vec > n_vec)) then
      print *, 'ERROR(lib_eigen_checker): n, n_vec, n_check_vec = ', n, n_vec, n_check_vec
      stop
    endif

    allocate(r(n), stat=ierr)
    if (ierr /= 0) stop 'Alloc error (la_eigen_solver_check) for w'

    allocate(rn(n_check_vec), stat=ierr)
    if (ierr /= 0) stop 'Alloc error (la_eigen_solver_check) for rn'
    rn(:)=0.0d0

    do j=1,n_check_vec
      r(:)=matmul(a,v(:,j))-eigen_level(j)*v(:,j)  ! redisual vector
      rn(j)=sqrt(abs(dot_product(r,r)))/sqrt(abs(dot_product(v(:,j),v(:,j))))
    enddo

    rn_max=maxval(rn)
    rn_ave=sum(rn)/dble(n_vec)

  end subroutine lib_eigen_checker

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lib_eigen_solver(mat, solver_type, n_vec, eigenvalues, eigenvectors)

   use solver_lapack, only : eigen_solver_lapack
   use solver_scalapack_all, only : eigen_solver_scalapack_all
   use solver_scalapack_select, only : eigen_solver_scalapack_select
   use matrix_io, only : sparse_mat
   implicit none
   type(sparse_mat), intent(in) :: mat
   character(len=*), intent(in) :: solver_type
   integer, intent(in) :: n_vec
   real(kind(1.d0)), intent(out), allocatable :: eigenvalues(:), eigenvectors(:, :)

   integer :: n
   real(kind(1.d0)), allocatable :: a(:, :)

   n = mat%size

   if ((n_vec < 0) .or. (n_vec > n)) then
     write(*,*) 'Error(lib_eigen_solver):n, n_vec=',n,n_vec
     stop
   endif

  call create_dense_matrix(0, mat, a)

  allocate(eigenvalues(n))
  eigenvalues(:) = 0.0d0
  allocate(eigenvectors(n, n))
  eigenvectors(:, :) = 0.0d0

  select case (trim(solver_type))
  case ('lapack')
    call eigen_solver_lapack(a, eigenvalues)
  case ('scalapack_all')
    call eigen_solver_scalapack_all(a, eigenvalues)
  case ('scalapack_select')
    call eigen_solver_scalapack_select(a, eigenvalues)
  case default
    write(*,*) 'Error(lib_eigen_solver):solver type=',trim(solver_type)
    stop
  end select

  eigenvectors = a

  end subroutine lib_eigen_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module solver_main
