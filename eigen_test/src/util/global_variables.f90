module global_variables
  implicit none

  public :: g_block_size, g_version, g_mpi_wtime_init
  integer :: g_block_size = 128
  character(*), parameter :: g_version = '20150625'
  real(8) :: g_mpi_wtime_init = 0d0
end module global_variables
