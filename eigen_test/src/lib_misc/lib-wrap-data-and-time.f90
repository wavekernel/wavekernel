module M_wrap_data_and_time
  implicit none
  private 
  public :: data_and_time_wrapper
!
contains
!
  subroutine data_and_time_wrapper(chara_wrk)
    implicit none
    character(len=256), intent(out) :: chara_wrk
    character(len=8)  :: chara_date
    character(len=10) :: chara_time
!
    call date_and_time(chara_date, chara_time)
    chara_wrk='Date: '//chara_date(1:4)//' '//chara_date(5:6)//' '//chara_date(7:8)//'; '// &
&                             'Time: '//chara_time(1:2)//' '//chara_time(3:4)//' '//chara_time(5:6)
!
  end subroutine data_and_time_wrapper
!
end module M_wrap_data_and_time



