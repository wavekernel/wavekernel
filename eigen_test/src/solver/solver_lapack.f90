!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_lapack
  implicit none

  private
  public :: eigen_solver_lapack

contains
  subroutine eigen_solver_lapack(mat,eigen_level)
   use time, only : get_wclock_time !(routine)
   implicit none

   real(kind(1.d0)), intent(inout) :: mat(:,:)   ! ( n x n ) matrix
   real(kind(1.d0)), intent(out) :: eigen_level(:)

   integer :: ierr, info
   integer :: n, lda, lwork

   real(kind(1.d0)), allocatable :: work(:)

   real(kind(1.d0)) :: time_origin, elapse_time

   n = size(mat, 1)
   lda = n
   lwork = n * n  ! Note: (lwork > 3*n-1 ) should be satisfied.

   allocate(work(lwork), stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR(eigen_solver_lapack): Alloc. error work'
     stop
   endif

   call get_wclock_time(time_origin)

   call dsyev("V", "U", n, mat, lda, eigen_level, work, lwork, info)

   call get_wclock_time(elapse_time, time_origin)

   write(*,'(a,i10, f20.10)')  ' solver result (LAPACK, DSYEV) : matrix size, time (sec) =', n, elapse_time

   deallocate(work, stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR(eigen_solver_lapack): Dealloc. error work'
     stop
   endif
  end subroutine eigen_solver_lapack
end module solver_lapack
