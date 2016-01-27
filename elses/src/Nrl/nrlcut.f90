!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module MNrlCutOff
  private

  public NrlCutOffFunction
  public NrlCutOffDerivative


contains

  !!
  ! a cut off function for hopping in NRL scheme.
  !
  real(8) function NRLCutOffFunction(r, rc, delta) result (oresult)
    implicit none
    real(8), intent(in) :: r, rc, delta

    oresult = 1d0/(1d0 + exp((r - rc) / delta))
  end function

  !!
  ! the derivative of a cut-off function in NRL scheme
  !
  real(8) function NRLCutOffDerivative(r, rc, delta) result (oresult)
    implicit none
    real(8), intent(in) :: r, rc, delta

    real(8) :: the_exp

    the_exp = exp((r - rc) / delta)
    oresult = - the_exp / (delta * (1 + the_exp)**2)
  end function

end module
