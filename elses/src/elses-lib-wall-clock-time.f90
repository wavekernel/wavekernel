!================================================================
! ELSES version 0.05
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
   public :: measure_clock_time_period
   public :: show_reset_period_of_clock
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
     implicit none
     real(kind=8), intent(out) :: elapse_time
     integer :: count, rate
!
     call system_clock(count, rate)
     elapse_time=dble(count)/dble(rate)
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
     implicit none
     real(kind=8), intent(out) :: elapse_time
     integer :: count, rate, max
!
     call system_clock(count, rate, max)
!
     if ( .not. initialized )then
        time_origin=dble(count)/dble(rate)
        initialized=.true.
     else
        if (count < count_previous) then
           time_origin=time_origin-dble(max)/dble(rate)
        endif
     endif
!
     elapse_time=dble(count)/dble(rate)-time_origin
     count_previous=count
!
   end subroutine get_elapse_wall_clock_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine measure_clock_time_period(time_present, time_previous, time_period)
!      ----> Measure the clock-time period 
!             with the correction for the possible reset of the system clock
!
     implicit none
     real(kind=8), intent(in)  :: time_present, time_previous
     real(kind=8), intent(out) :: time_period
     integer :: count, rate, max
!
     time_period = time_present - time_previous
!
     if (time_period < 0) then
       call system_clock(count, rate, max)
       time_period = time_period + dble(max)/dble(rate)
     endif   
!
   end subroutine measure_clock_time_period
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine show_reset_period_of_clock(time_reset_period)
!      ----> Show the reset period of the system clock
     implicit none
     real(kind=8), intent(out) :: time_reset_period
     integer :: count, rate, max
!
     call system_clock(count, rate, max)
     time_reset_period = dble(max)/dble(rate)
!
   end subroutine show_reset_period_of_clock
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


