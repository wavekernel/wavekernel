!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
program elses_main
   use M_md_main,   only : elses_md_main_01 !(routines)
   use M_00_v_info, only : elses_00_version_info !(routines)
   use M_band_calc, only : elses_band_calculation !(routines)
   use M_options                                     !(uchanged)
   use M_md_dst,    only : set_dst_initial, set_dst_final !(routines)
!  use M_lib_omp,   only : show_omp_info !(routines)
   use M_lib_omp,   only : show_mpi_omp_info !(routines)
   use M_io_dst_write_log, only : set_dst_write_log_file !(routine)
!
   implicit none
!
!  write(6,*)'@@ ELSES Project ' 
!
   call elses_process_options
!
!  call elses_00_version_info
!     ----> Display the version infomation (non-essential)
!
   call set_dst_initial
!
   call set_dst_write_log_file
!     ----> Set logfile in DST calculation
!
   call show_mpi_omp_info
!     ----> show MPI and OpenMP information
!
   select case(config%option%functionality)
      case("band":"band@") ! '@' is the previous character of 'A' in ascii code table.
         call elses_band_calculation
      case default
         call elses_md_main_01
!     ----> Main routine for MD
!
   end select
!
   call set_dst_final
!
!  write(6,*)'  ... ELSES ended successfully.'
!
 end program elses_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
