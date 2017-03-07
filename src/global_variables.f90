module wk_global_variables_m
  implicit none

  integer :: g_wk_block_size = 64, g_wk_master_pnum = 0
  real(8) :: g_wk_mpi_wtime_init = 0d0

  character(*), parameter :: kVersion = '20160624'
  real(8), parameter :: kAuPerAngstrom = 1.8897259885789d0
  real(8), parameter :: kMassCarbon = 2.1874661e4  ! [a.u. (mass)]
  real(8), parameter :: kMassHydrogen = 1.8371526e3  ! [a.u. (mass)]
  real(8), parameter :: kMassOxygen = 2.9156944e4  ! [a.u. (mass)]
  real(8), parameter :: kBoltzmann = 1d0  ! Boltzmann's constant [a.u.]
  real(8), parameter :: kJoulePerEh = 4.35974434e-18  ! [J / a.u. (energy)]
  real(8), parameter :: kBoltzmannInSI = 1.3806488e-23  ! [J / K]
  real(8), parameter :: kAuPerKelvin = kBoltzmannInSI / kJoulePerEh  ! Temperature
  real(8), parameter :: kPsecPerAu = 2.418884326505e-5  ! Time
  real(8), parameter :: kAccelRatio = 1d-3
  real(8), parameter :: kPi = 3.141592653589793d0
  real(8), parameter :: charge_sum_error_tol = 1d-6
  complex(kind(0d0)), parameter :: kOne = (1d0, 0d0), kZero = (0d0, 0d0), &
       kImagUnit = (0d0, 1d0)
end module wk_global_variables_m
