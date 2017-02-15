module wp_maxwell_distribution_m
  implicit none
  integer, parameter :: num_max_stock = 2 * 100
  real(8), parameter :: pi = 3.1415926535897932
  real(8) :: xs(num_max_stock)
  integer :: num_stock = 0

contains

  subroutine init_stock()
    integer :: i
    real(8) :: x, y
    call random_number(xs)
    do i = 1, num_max_stock / 2
      x = 1d0 - xs(2 * i - 1)
      y = 1d0 - xs(2 * i)
      xs(2 * i - 1) = dsqrt(-2d0 * dlog(x)) * dcos(2d0 * pi * y)
      xs(2 * i) = dsqrt(-2d0 * dlog(x)) * dsin(2d0 * pi * y)
    end do
    num_stock = num_max_stock
  end subroutine init_stock


  real(8) function get_normal()
    if (num_stock == 0) then
      call init_stock()
    end if
    get_normal = xs(num_max_stock - num_stock + 1)
    num_stock = num_stock - 1
  end function get_normal


  real(8) function maxwell_rvs(scale)
    real(8), intent(in) :: scale
    real(8) :: v1, v2, v3
    v1 = get_normal()
    v2 = get_normal()
    v3 = get_normal()
    maxwell_rvs = dsqrt(v1 ** 2d0 + v2 ** 2d0 + v3 ** 2d0) * scale
  end function maxwell_rvs


  real(8) function maxwell_cdf(scale, x)
    real(8) :: scale, x
    maxwell_cdf = erfunc(x / dsqrt(2d0) / scale) - dsqrt(2d0 / pi) * &
         x * dexp(-1d0 * x ** 2d0 / (2d0 * scale ** 2d0)) / scale
  end function maxwell_cdf


  ! http://www.pa.msu.edu/~duxbury/phy201_f00/worksheet2_f00/node2.html
  real(8) function erfunc(x)
    real(8) :: x
    real(8) :: d, x2n1, factorial
    integer :: n, sign
    integer, parameter :: max_iter = 200

    if (dabs(x) > 5d0) then
      erfunc = 1d0
    end if

    factorial = 1d0
    x2n1 = x
    sign = 1
    erfunc = x
    do n = 1, max_iter
      factorial = factorial * n  ! n!
      x2n1 = x2n1 * (x ** 2d0)  ! x ** (2 * n + 1)
      sign = sign * (-1)  ! (-1) ** n
      d = sign * x2n1 / factorial / (2 * n + 1)
      erfunc = erfunc + d
      if (dabs(d) < 1e-16) then
        erfunc = erfunc * 2d0 / dsqrt(pi)
        return
      end if
    end do
  end function erfunc
end module wp_maxwell_distribution_m
