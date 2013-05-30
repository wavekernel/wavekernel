module M_lib_eigen_routines
!
  implicit none
!
  private
  public :: lib_eigen_solver
  public :: lib_eigen_checker
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine lib_eigen_checker(a, v, eigen_level, n_vec, rn_ave, rn_max)
!
    implicit none
    real(kind(1.d0)),  intent(in) :: a(:,:)
    real(kind(1.d0)),  intent(in) :: v(:,:)
    real(kind(1.d0)),  intent(in) :: eigen_level(:)
    integer,           intent(in) :: n_vec
    real(kind(1.d0)), intent(out) :: rn_ave, rn_max      ! residual norm average, max
!
    real(kind(1.d0)), allocatable :: r(:)   ! residual vector (work array)
    real(kind(1.d0)), allocatable :: rn(:)  ! residual norm
!
    integer :: n, ierr, j
!
    n=size(a,1)
!
    if ((n_vec < 1 ) .or. (n_vec > n)) then
      write(*,*)'ERROR(lib_eigen_checker):n,n_vec=',n,n_vec
      stop
    endif
!
    allocate(r(n), stat=ierr)
    if (ierr /= 0) stop 'Alloc error (la_eigen_solver_check) for w'
!
    allocate(rn(n_vec), stat=ierr)
    if (ierr /= 0) stop 'Alloc error (la_eigen_solver_check) for rn'
    rn(:)=0.0d0
!
    do j=1,n_vec
      r(:)=matmul(a,v(:,j))-eigen_level(j)*v(:,j)  ! redisual vector
      rn(j)=sqrt(abs(dot_product(r,r)))/sqrt(abs(dot_product(v(:,j),v(:,j))))
    enddo
!
    rn_max=maxval(rn)
    rn_ave=sum(rn)/dble(n_vec)
!
  end subroutine lib_eigen_checker
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lib_eigen_solver(a,eigen_level,solver_type,n_vec)
!
   use M_lib_eigen_solver_lapack, only : eigen_solver_lapack
   use M_lib_eigen_solver_scalapack_all, only : eigen_solver_scalapack_all
   use M_lib_eigen_solver_scalapack_select, only : eigen_solver_scalapack_select
   implicit none
   real(kind(1.d0)), intent(inout) :: a(:,:)
   real(kind(1.d0)), intent(out)   :: eigen_level(:)
   character(len=*)      , intent(in)    :: solver_type
   integer               , intent(in)    :: n_vec
   integer :: n
!
   if (size(a,1) /= size(a,2)) then
     write(*,*)'ERROR(ext_eigen_solver_wrapper):size=',size(a,1), size(a,2)
     stop
   endif
!
   if (size(a,1) /= size(eigen_level,1)) then
     write(*,*)'ERROR(ext_eigen_solver_wrapper):size=',size(a,1), size(eigen_level,1)
     stop
   endif
!
   n=size(a,1)
!
   if ((n_vec < 0) .or. (n_vec > n)) then
     write(*,*) 'Error(lib_eigen_solver):n, n_vec=',n,n_vec
     stop
   endif
!
   select case (trim(solver_type))
     case ('lapack')
       call eigen_solver_lapack(a,eigen_level)
    case ('scalapack_all')
       call eigen_solver_scalapack_all(a,eigen_level)
    case ('scalapack_select')
       call eigen_solver_scalapack_select(a,eigen_level)
     case default
       write(*,*) 'Error(lib_eigen_solver):solver type=',trim(solver_type)
       stop
   end select
!
  end subroutine lib_eigen_solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_eigen_routines
!
