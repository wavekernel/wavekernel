module time
  implicit none

  private
  public :: get_wall_clock_base_count, get_wall_clock_time

contains

  subroutine get_wall_clock_base_count(base_count)
    integer, intent(out) :: base_count

    integer :: rate, max  ! Dummy variables

    call system_clock(base_count, rate, max)
  end subroutine get_wall_clock_base_count


  subroutine get_wall_clock_time(base_count, time)
    integer, intent(in) :: base_count
    double precision, intent(out) :: time

    integer :: count, rate, max, diff

    call system_clock(count, rate, max)

    diff = count - base_count
    if (diff < 0) then
      diff = max + diff
    end if

    time = dble(diff) / dble(rate)
  end subroutine get_wall_clock_time
end module time
