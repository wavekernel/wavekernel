!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_cubic_harmonics
!
  implicit none
!
!
  private
  public :: cubic_harmonics_func_mod
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ modified cubic harmonics func　for y = r^l Q(x/r, y/r, z,r)
!     x, y, z: 
!     k: orbital index
!
  function cubic_harmonics_func_mod (x, y, z, k) result(result_value)
!            
    implicit none
    real(8), intent(in) :: x, y, z
    integer, intent(in) :: k
    real(8), parameter :: eps=1.0d-10
    real(8), parameter :: pi=.314159265358979323D1
    integer :: ierr
    real(8) :: result_value
    real(8) :: r, a
!
    r = dsqrt(x*x + y*y + z*z)
    a = 1.0d0/dsqrt(4.0d0*pi)
!
    select case(k)
      case (1) ! s type 
        result_value = a 
      case (2) ! px type
        result_value = a * x * dsqrt(3.0d0)
      case (3) ! py type
        result_value = a * y * dsqrt(3.0d0)
      case (4) ! pz type
        result_value = a * z * dsqrt(3.0d0)
      case (5) ! xy type
        result_value = a* x*y * dsqrt(15.0d0)
      case (6) ! yz type
        result_value = a* y*z * dsqrt(15.0d0)
      case (7) ! zx type
        result_value = a* z*x * dsqrt(15.0d0)
      case (8) ! x^2-y^2 type
        result_value = a* (x*x-y*y) * dsqrt(15.0d0/4.0d0)
      case (9) ! 3z^2-r^2 type
        result_value = a* (3*z*z-r*r) * dsqrt(5.0d0/4.0d0)
      case default
        write(*,*)'ERROR(cubic_harmonics_func):unsupported:k=',k
        stop
    end select
!
  end function cubic_harmonics_func_mod
!
end module M_lib_cubic_harmonics



