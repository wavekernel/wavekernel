!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_mpi_wrapper !(DUMMY ROUTINES)
!
  use M_lib_dst_info, only  : log_unit  
  implicit none
  logical, parameter :: sum_up_wtime      = .false.
!
  real(kind(1d0)), allocatable :: total_allreduce_time(:)
  real(kind(1d0)), allocatable :: total_barrier_time(:)
!
! real(kind(1d0)), allocatable :: matvec_timer_in_thread(:) 
!
  real(kind(1d0)), allocatable :: time_origin(:)  ! used in get_wall_clock_time
  integer                      :: count_previous  ! used in get_wall_clock_time
!
  private
  public :: mpi_wrapper_initial      !(DUMMY ROUTINE)
  public :: mpi_wrapper_final        !(DUMMY ROUTINE)
!
  public :: mpi_wrapper_allreduce_r1 !(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_r2 !(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_r3 !(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_r4 !(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_r5 !(DUMMY ROUTINE)
!
  public :: mpi_wrapper_allreduce_i2 !(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_i1 !(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_i1c!(DUMMY ROUTINE)
  public :: mpi_wrapper_allreduce_r0 !(DUMMY ROUTINE)
!
  public :: mpi_wrapper_allreduce_minmax_r1 !(DUMMY ROUTINE)
!
  public :: mpi_wrapper_barrier_time
  public :: mpi_wrapper_allreduce_chk
!
  public :: sum_up_wtime
  public :: total_allreduce_time
  public :: total_barrier_time
!
  public :: mpi_wrapper_wtime
  public :: get_wall_clock_time 
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_wtime(time_data)
    implicit none
    real(8), intent(out) :: time_data
    integer :: count, rate, max
!
    call system_clock(count, rate)
    time_data=dble(count)/dble(rate)
!
  end subroutine mpi_wrapper_wtime
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
   subroutine get_wall_clock_time(elapse_time)
!      ----> Get the elapse wall-clock time in sec.
!             from the 'time origin'
!
     implicit none
     real(kind=8), intent(out) :: elapse_time
     integer :: count, rate, max
     integer :: ierr
!
     call system_clock(count, rate, max)
!
     if ( .not. allocated(time_origin))then
        allocate(time_origin(1), stat=ierr)
        time_origin(1)=dble(count)/dble(rate)
        count_previous=count
     else
        if (count < count_previous) then
           time_origin(1)=time_origin(1)-dble(max)/dble(rate)
        endif
     endif
!
     elapse_time=dble(count)/dble(rate)-time_origin(1)
     count_previous=count
!
   end subroutine get_wall_clock_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_initial(myrank, nprocs, mpi_is_active, true_routine)
!
    implicit none
    integer, intent(out) :: myrank
    integer, intent(out) :: nprocs
    logical, intent(out) :: mpi_is_active
    logical, intent(out) :: true_routine
!
    write(*,*)'INFO-MPI:NO MPI routine is used.'
!
!   if (log_unit > 0 )then
!     write(log_unit,*)'INFO-MPI:mpi_wrapper_initial(dummy routine)'
!   endif
!
    true_routine = .false.
    myrank=0
    nprocs=1
    mpi_is_active = .false.
!
!    
  end subroutine mpi_wrapper_initial
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_final
!
    implicit none
!
    if (log_unit > 0 )then
      write(log_unit,'(a)') '..... ELSES ended successfully (without MPI)'
      write(*,'(a)')        '..... ELSES ended successfully (without MPI)'
    endif   
!
    return
!
  end subroutine mpi_wrapper_final
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r5(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:,:,:,:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_r5
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r4(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:,:,:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_r4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r3(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:,:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_r3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r2(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_r2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r1(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_r1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_i2(wrk_array)
!
    implicit none
    integer, intent(inout) :: wrk_array(:,:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_i2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_i1(wrk_array)
!
    implicit none
    integer, intent(inout) :: wrk_array(:)
!
    return
!
  end subroutine mpi_wrapper_allreduce_i1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_i1c(array_in,array_out) 
!
    implicit none
    integer, intent(in)  :: array_in(:)
    integer, intent(out) :: array_out(:)
!
    array_out(:)=array_in(:)
!
  end subroutine mpi_wrapper_allreduce_i1c
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r0(wrk_variable)
!
    implicit none
    real(8), intent(inout) :: wrk_variable
!
    return
!
  end subroutine mpi_wrapper_allreduce_r0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_minmax_r1(wrk_array,mode_chara)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:)
    character(len=*), intent(in) :: mode_chara
!
    return
!
  end subroutine mpi_wrapper_allreduce_minmax_r1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_barrier_time(time_check)  
!
    use M_config,    only : config !(unchanged)
    implicit none
    real(8), intent(out) :: time_check
!
    time_check=0.0d0
    if (config%calc%distributed%use_mpi_barrier) then
      time_check=0.0d0
    endif
!
  end subroutine mpi_wrapper_barrier_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_chk(time_check)  
!
    implicit none
    real(8), intent(out) :: time_check
!
    time_check=0.0d0
!
  end subroutine mpi_wrapper_allreduce_chk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_mpi_wrapper



