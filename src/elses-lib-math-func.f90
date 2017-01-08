!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_math_func
!
  implicit none
!
!
  private
  public :: Fermi_Dirac_Func
  public :: get_determinant
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ Fermi Dirac function
!       f(x) = 1 / (1 + exp(x) )
!
  function Fermi_Dirac_Func (x)
!
! Note: This routine was written so as to avoid 
!          an floating point exception in exp(x) 
!           with a large value of x
!            
    implicit none
    real(8) :: Fermi_Dirac_Func
    real(8), intent(in) :: x
!
    if(x .gt. 100.d0) then
        Fermi_Dirac_Func=0.0d0
    elseif(x .lt.-100.d0) then
        Fermi_Dirac_Func=1.0d0
    else
        Fermi_Dirac_Func=1.0d0/(dexp(x)+1.0d0)
    endif
!
  end function Fermi_Dirac_Func
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ Calculation of the determinant for matrices
!     (now only for 3x3 matrices)
!       
  subroutine get_determinant(a,det)
    implicit none
    real(8), intent(in)  :: a(:,:)
    real(8), intent(out) :: det
!
    if ((size(a,1) /= 3) .or. (size(a,2) /= 3)) then
      write(*,*)'ERROR(get_determinant)',size(A,1)
      stop
    endif   
!
    det = a(1,1)*a(2,2)*a(3,3) +a(2,1)*a(3,2)*a(1,3) +a(3,1)*a(1,2)*a(2,3) & 
&        -a(1,1)*a(3,2)*a(2,3) -a(2,1)*a(1,2)*a(3,3) -a(3,1)*a(2,2)*a(1,3) 
!
  end subroutine get_determinant
!
end module M_lib_math_func


