module time
  implicit none

  private
  public :: get_wclock_time, data_and_time_wrapper

contains
  subroutine get_wclock_time(time_present, time_origin)
    double precision, intent(out) :: time_present
    double precision, optional :: time_origin

    double precision :: time_origin_wrk
    integer :: count, rate, max

    if (present(time_origin)) then
      time_origin_wrk=time_origin
    else
      time_origin_wrk=0.0d0
    endif

    call system_clock(count, rate, max)

    time_present=dble(count)/dble(rate)-time_origin_wrk

! @ Correct the present time, when the time count was initialized
    if (time_present < 0) then
      time_present = time_present + dble(max)/dble(rate)
    endif
  end subroutine get_wclock_time


  subroutine data_and_time_wrapper(chara_wrk)
    character(len=256), intent(out) :: chara_wrk
    character(len=8) :: chara_date
    character(len=10) :: chara_time

    call date_and_time(chara_date, chara_time)
    chara_wrk='Date: '//chara_date(1:4)//' '//chara_date(5:6)//' '//chara_date(7:8)//'; '// &
&                             'Time: '//chara_time(1:2)//' '//chara_time(3:4)//' '//chara_time(5:6)
  end subroutine data_and_time_wrapper
end module time
