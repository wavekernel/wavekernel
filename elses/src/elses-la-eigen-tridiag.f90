!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_eig_tridiag
!
  implicit none
  integer m
  integer nz
!
  private
  public calc_eig_tridiag
!
  contains
!
  subroutine calc_eig_tridiag(mat,eig_val,eig_vec,ierr)
!
    implicit none
    real(8), intent(in)  :: mat(:,:)
    real(8), intent(out) :: eig_val(:)
    real(8), intent(out) :: eig_vec(:,:)
    integer, intent(out) :: ierr
    real(8), allocatable :: a(:), b(:)
    real(8), allocatable :: work(:) 
    real(8), allocatable :: eig_vec_wrk(:,:)
    integer              :: n, ldz, work_size, j
!
    n=size(mat,1)
!
    allocate(a(n), stat=ierr)
    if (ierr /=0) then
      stop 'Alloc error in calc_eig_tridiag'
    endif   
!
    allocate(b(n-1), stat=ierr)
    if (ierr /=0) then
      stop 'Alloc error in calc_eig_tridiag'
    endif   
!
    allocate(eig_vec_wrk(n,n), stat=ierr)
    if (ierr /=0) then
      stop 'Alloc error in calc_eig_tridiag'
    endif   
!
    ldz=n
    work_size=max(1, 2*n-2)
    allocate(work(work_size), stat=ierr)
    if (ierr /=0) then
      stop 'Alloc error in calc_eig_tridiag'
    endif   
!
    do j=1,n
      a(j) = mat(j,j)
    enddo   
!
    do j=1,n-1
      b(j) = ( mat(j,j+1) + mat(j+1,j) ) / 2.0d0
    enddo   
!
    call dstev('V', n, a, b, eig_vec_wrk, ldz, work ,ierr)
!
    eig_val(:)   = a(:)
    eig_vec(:,:) = eig_vec_wrk(:,:)
!
  end subroutine calc_eig_tridiag
!
end module M_la_eig_tridiag



