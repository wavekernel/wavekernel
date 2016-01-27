!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_dst_info
!
  implicit none
  logical :: mpi_is_active
  integer :: myrank   
  integer :: nprocs   
  integer :: log_unit 
  integer :: output_unit
! character(len=*), parameter :: output_filename="Output.txt"
  logical :: root_node
  logical :: true_mpi_wrapper
!
  private
  public :: mpi_is_active
  public :: myrank   
  public :: nprocs   
  public :: log_unit
  public :: output_unit
! public :: output_filename
  public :: root_node
  public :: true_mpi_wrapper
!
!   
end module M_lib_dst_info

