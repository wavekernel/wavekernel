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
!@@  Quad precision test for pi (= 4.0 * atan(1.0))
!
  subroutine quad_prec_test_for_pi(log_unit_input, quad_prec_test_verbose)
    implicit none
    integer,     intent(in)  :: log_unit_input
    logical,     intent(in)  :: quad_prec_test_verbose
!    
    integer, parameter  :: DOUBLE_PRECISION=kind(1d0)
    integer,  parameter :: QUAD_PRECISION = selected_real_kind(p=30)  
    real(QUAD_PRECISION)   :: value_quad
    real(DOUBLE_PRECISION) :: value_double
    character(len=256)     :: chara_wrk
    integer                :: log_unit
!
    if (quad_prec_test_verbose) then
      log_unit = log_unit_input
    else  
      log_unit = 0
    endif  
!
    value_double = 4.0d0              * atan(1.0d0)
    value_quad   = 4.0_QUAD_PRECISION * atan(1.0_QUAD_PRECISION)
!
    if (log_unit > 0) then
      write(log_unit,'(a)')'@@ Quad precision test for pi (= 4.0 * atan(1.0))  '
      write(log_unit,'(a,f45.40)') ' Double precision : pi = ',value_double
      write(log_unit,'(a,f45.40)') ' Quad   precision : pi = ',value_quad
      write(log_unit,'(a)')       ' Reference  value : pi =    3.141592653589793238462643383279502884197169399375'
      write(log_unit,'(a)')       '     (eye guide)              123456789012345678901234567890123456789012345678'
    endif  
!
    write(chara_wrk,'(f45.40)') value_quad
!
    if (index(chara_wrk, '3.141592653589793238462643383279502') == 0) then
      write(*,'(a)')'@@ Quad precision test for pi (= 4.0 * atan(1.0))  '
      write(*,'(a,f45.40)') ' Double precision : pi = ',value_double
      write(*,'(a,f45.40)') ' Quad   precision=: pi = ',value_quad
      write(*,'(a)')       ' Reference  value : pi =    3.141592653589793238462643383279502884197169399375'
      write(*,'(a)')       '     (eye guide)              123456789012345678901234567890123456789012345678'
      write(*,*)' ERROR!(Quad precision test)'
      stop
    endif   
!
    
!
  end subroutine quad_prec_test_for_pi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_quad_prec


