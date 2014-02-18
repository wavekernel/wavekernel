!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_lapack
  implicit none

  private
  public :: eigen_solver_lapack

contains

  subroutine eigen_solver_lapack(mat, eigenpairs)
   use time, only : get_wall_clock_base_count, get_wall_clock_time
   use matrix_io, only : sparse_mat
   use distribute_matrix, only : convert_sparse_matrix_to_dense
   use eigenpairs_types, only : eigenpairs_types_union

   type(sparse_mat), target, intent(in) :: mat
   type(eigenpairs_types_union), intent(out) :: eigenpairs

   integer :: n, lda, lwork, info

   double precision, allocatable :: work(:)

   integer :: base_count
   double precision :: elapse_time

   eigenpairs%type_number = 1

   call convert_sparse_matrix_to_dense(mat, eigenpairs%local%vectors)

   n = mat%size
   lda = n
   lwork = n * n  ! Note: (lwork > 3*n-1 ) should be satisfied.

   allocate(eigenpairs%local%values(n), work(lwork))

   call get_wall_clock_base_count(base_count)

   call dsyev("V", "U", n, eigenpairs%local%vectors, lda, &
        eigenpairs%local%values, work, lwork, info)

   call get_wall_clock_time(base_count, elapse_time)

   write(*,'(a,i10, f20.10)')  ' solver result (LAPACK, DSYEV) : matrix size, time (sec) =', n, elapse_time
  end subroutine eigen_solver_lapack
end module solver_lapack
