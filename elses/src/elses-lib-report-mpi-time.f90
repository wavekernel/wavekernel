!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_rep_mpi_time
!
  implicit none
  real(8), allocatable :: previous_allreduce_time(:) ! DIFFERENT VALUES AMONG NODES
  real(8), allocatable :: previous_barrier_time(:)   ! DIFFERENT VALUES AMONG NODES
!
  private
  public :: get_lap_of_allreduce_time
  public :: get_lap_of_barrier_time
  public :: record_allreduce_barrier_time
!
  contains
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_lap_of_allreduce_time(lap_time)
!
    use M_lib_dst_info,    only : myrank, mpi_is_active
    use M_lib_mpi_wrapper, only : sum_up_wtime         !(unchanged)
    use M_lib_mpi_wrapper, only : total_allreduce_time !(unchanged)
    implicit none
    real(8), intent(out) :: lap_time
    integer :: ierr
!
    if (.not. mpi_is_active) then 
      lap_time=0.0d0 
      return
    endif  
!
    if (.not. sum_up_wtime) then 
      lap_time=0.0d0 
      return
    endif  
!
    if (.not. allocated(previous_allreduce_time)) then
      allocate(previous_allreduce_time(1),stat=ierr) 
      if (ierr /= 0) then
        stop 'Alloc. error (get_lap_of_allreduce_time)'
      endif   
      previous_allreduce_time(1)=0.0d0
    endif   
!
    lap_time=total_allreduce_time(myrank)-previous_allreduce_time(1)
    previous_allreduce_time(1)=total_allreduce_time(myrank)
!
  end subroutine get_lap_of_allreduce_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_lap_of_barrier_time(lap_time)
!
    use M_lib_dst_info,    only : myrank, mpi_is_active
    use M_lib_mpi_wrapper, only : sum_up_wtime         !(unchanged)
    use M_lib_mpi_wrapper, only : total_barrier_time !(unchanged)
    implicit none
    real(8), intent(out) :: lap_time
    integer :: ierr
!
    if (.not. mpi_is_active) then 
      lap_time=0.0d0 
      return
    endif  
!
    if (.not. sum_up_wtime) then 
      lap_time=0.0d0 
      return
    endif  
!
    if (.not. allocated(previous_barrier_time)) then
      allocate(previous_barrier_time(1),stat=ierr) 
      if (ierr /= 0) then
        stop 'Alloc. error (get_lap_of_barrier_time)'
      endif   
      previous_barrier_time(1)=0.0d0
    endif   
!
    lap_time=total_barrier_time(myrank)-previous_barrier_time(1)
    previous_barrier_time(1)=total_barrier_time(myrank)
!
  end subroutine get_lap_of_barrier_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine record_allreduce_barrier_time(filename)
!
    use M_lib_dst_info,    only : myrank, nprocs, mpi_is_active    !(unchanged)
    use M_lib_mpi_wrapper, only : sum_up_wtime             !(unchanged)
    use M_lib_mpi_wrapper, only : total_allreduce_time     !(unchanged)
    use M_lib_mpi_wrapper, only : total_barrier_time       !(unchanged)
    use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_r1 !(routine)
    use elses_mod_file_io, only : vacant_unit              !(function)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: j, iunit
!
    if (.not. mpi_is_active) return
    if (.not. sum_up_wtime) return
    if (.not. allocated(total_allreduce_time)) return 
    if (.not. allocated(total_barrier_time)) return 
!
    call mpi_wrapper_allreduce_r1(total_allreduce_time)
    call mpi_wrapper_allreduce_r1(total_barrier_time)
!
    iunit=vacant_unit()
!
    open(iunit, file=filename , status='unknown')
!
    if (myrank == 0) then
      do j=0,nprocs-1
      write(iunit,'(a,i10,2f20.10)')  'node_id, allreduce_time, barrier_time= ', & 
&                                  j, total_allreduce_time(j), total_barrier_time(j) 
      enddo   
    endif  
!
  end subroutine record_allreduce_barrier_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_rep_mpi_time
