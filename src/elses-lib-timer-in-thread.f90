!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_timer_in_thread 
!
  implicit none
  real(8), allocatable :: matvec_timer_in_thread(:,:)
!
  private
!
  public :: matvec_timer_in_thread
  public :: init_timer_in_thread
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_timer_in_thread
    implicit none
    integer :: n_threads
    integer :: ierr
!
    call wrap_get_num_thread(n_threads)
    write(*,*)'INFO:init_timer_in_thread:n_threads=', n_threads
!
    if (allocated(matvec_timer_in_thread)) deallocate(matvec_timer_in_thread)
!
    ierr=0
    if (n_threads < 1) ierr=1
    if (n_threads > 1024) ierr=1
!
    if (ierr == 1) then
      write(*,*)'ERROR(init_timer_in_thread):n_threads=',n_threads
      stop
    endif
!
    allocate(matvec_timer_in_thread(n_threads,2), stat=ierr)
    if (ierr /=0) then
      write(*,*)'Alloc. Error.matvec_timer_in_thread'
      stop
    endif
    matvec_timer_in_thread(:,:)=0.0d0
!
  end subroutine init_timer_in_thread
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wrap_get_num_thread(n_threads)
    implicit none
    integer, intent(out) :: n_threads
    integer :: ipe, npe
    integer :: omp_get_thread_num
    integer :: omp_get_num_threads
!
!$omp parallel default(shared) &
!$omp& private(ipe,npe)
    ipe=1
    npe=1
!$  ipe=omp_get_thread_num()+1
!$  npe=omp_get_num_threads()
    if (ipe == 1) n_threads=npe
!$omp end parallel 
!
  end subroutine wrap_get_num_thread
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_timer_in_thread 
