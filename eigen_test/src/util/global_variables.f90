module global_variables
  implicit none

  public :: g_block_size, g_version
  integer :: g_block_size = 128
  character(*), parameter :: g_version = '20150421'
end module global_variables
