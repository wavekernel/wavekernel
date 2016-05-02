!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_wall_clock_time
!
!
   implicit none
   real(kind=8) :: time_origin
   integer      :: count_previous
   logical      :: initialized=.false.
!
   private
!
   public :: get_system_clock_time
   public :: get_elapse_wall_clock_time
!  public :: measure_clock_time_period
!  public :: show_reset_period_of_clock
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine get_system_clock_time(elapse_time)
!      ----> Get the system_clock time in sec.
!
     use M_lib_mpi_wrapper, only : mpi_wrapper_wtime
     implicit none
     real(kind=8), intent(out) :: elapse_time
!
     call mpi_wrapper_wtime(elapse_time)
!
   end subroutine get_system_clock_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine get_elapse_wall_clock_time(elapse_time)
!      ----> Get the elapse wall-clock time in sec.
!             from the 'time origin'
!
     use M_lib_mpi_wrapper, only : get_wall_clock_time
     implicit none
     real(kind=8), intent(out) :: elapse_time
!
     call get_wall_clock_time(elapse_time)
!
   end subroutine get_elapse_wall_clock_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_wall_clock_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine tclock(elapse_time)
   use  M_wall_clock_time, only : get_system_clock_time
   implicit none
   real(kind=8), intent(out) :: elapse_time
!
   call get_system_clock_time(elapse_time)
!
end subroutine tclock


