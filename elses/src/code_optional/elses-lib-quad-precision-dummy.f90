!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_quad_prec
!
  implicit none
  private
  public :: quad_prec_test_for_pi
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@  DUMMY routine 
!     This routine should be used, 
!          only when the compiler does not support the quad precision format.
!
  subroutine quad_prec_test_for_pi(log_unit_input, quad_prec_test_verbose)
    implicit none
    integer,     intent(in)  :: log_unit_input
    logical,     intent(in)  :: quad_prec_test_verbose
    integer                  :: log_unit
!    
    if (quad_prec_test_verbose) then
      log_unit = log_unit_input
    else  
      log_unit = 0
    endif  
!
!
    if (log_unit > 0) then
      write(log_unit,'(a)')'@@ DUMMY routine of quad precision test'
    endif  
!
  end subroutine quad_prec_test_for_pi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_quad_prec


