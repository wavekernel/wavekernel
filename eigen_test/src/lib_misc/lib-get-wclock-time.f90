module M_get_clock_time
  implicit none
  private 
  public :: get_wclock_time
!
contains
!
  subroutine get_wclock_time(time_present, time_origin)
    implicit none
    real(kind(1.0d0)), intent(out) :: time_present
    real(kind(1.0d0)), optional    :: time_origin
!
    real(kind(1.0d0))              :: time_origin_wrk
    integer :: count, rate, max
!
    if (present(time_origin)) then
      time_origin_wrk=time_origin
    else
      time_origin_wrk=0.0d0
    endif
!   
    call system_clock(count, rate, max)
!
    time_present=dble(count)/dble(rate)-time_origin_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Correct the present time, when the time count was initialized
!
    if (time_present < 0) then
      time_present = time_present + dble(max)/dble(rate) 
    endif   
!
!
  end subroutine get_wclock_time
!
end module M_get_clock_time



