!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_sto_radial_func
!
  implicit none
!
!
  private
  public :: sto_radial_func_mod
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ modified STO radial function : R(r)/r^l in the double-zeta form
!
  function sto_radial_func_mod(r_in, z1, z2, c1_in, c2_in, n_in, k) result(y)
!      k : orbital index 
!            
    implicit none
    real(8), intent(in) :: r_in, z1, z2, c1_in, c2_in
    integer, intent(in) :: n_in, k
    real(8), parameter :: eps=1.0d-14
    real(8), parameter :: pi=.314159265358979323D1
    logical :: double_zeta
    integer :: ierr
    real(8) :: y
    real(8) :: r, n, l, zeta, zeta2, c1, c2
    real(8) :: f, f2
    real(8) :: a1, a2
!
    n=dble(n_in)
    zeta =z1
    zeta2=z2
    c1=c1_in
    c2=c2_in
!
    select case(k)
      case (1)         ! s type
        l=0.0d0 
      case (2,3,4)     ! p type
        l=1.0d0
      case (5,6,7,8,9) ! d type
        l=2.0d0
      case default
        write(*,*)'ERROR(sto_radial_func_mod):k=', k
    end select
!
    if (n_in < 1) then
      write(*,*)'ERROR(sto_radial_func):n_in = ', n_in
      stop
    endif
!
    if (r_in < -eps) then
      write(*,*)'ERROR(sto_radial_func):r_in = ', r_in
      stop
    endif
!
    if (r_in < 0.0d0) then 
      r=0.0d0
    else
      r=r_in
    endif
!
    if (dabs(c2) < eps) then
      double_zeta = .false.
    else
      double_zeta = .true.
    endif
!
    if (.not. double_zeta) then
      f=1.0d0
      c1=1.0d0
      c2=0.0d0
      zeta2=1.0d0 ! dummy value
    else
      f2=c1*c1+c2*c2+2.0d0*c1*c2*(4*zeta*zeta2/(zeta+zeta2)**2.0d0)**(n+0.5d0)
      if (f2 < eps) then
        write(*,*)'ERROR(sto_radial_func):f2=', f2
        stop
      endif
      f=dsqrt(f2)
    endif
!
    a1=dsqrt( (2.0d0*zeta)**(2.0d0*n+1.0d0)  / factorial_func(2*n_in) )
    a2=dsqrt( (2.0d0*zeta2)**(2.0d0*n+1.0d0) / factorial_func(2*n_in) )
!
    y= r**(n-l-1.0d0) / f * ( c1*a1*dexp(-zeta*r) + c2*a2*dexp(-zeta2*r) )
!
  end function sto_radial_func_mod
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ factorial function ; y(n) = n!
!
  function factorial_func(n) result(y)
    implicit none
    integer, intent(in) :: n
    real(8)             :: y
    integer             :: k
!
    if (n < 0) then
      write(*,*)' ERROR(factorial_func):n=',n
      stop
    endif
!
    y=1.0d0
    if (n > 1) then
     do k=1,n
       y=y*dble(k)
     enddo
    endif
!
  end function factorial_func
!
end module M_lib_sto_radial_func




