!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_eig_solver_center   ! DUMMY routines
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
! use eigen_test
! use global_variables  
  implicit none
!
  private
  public :: eig_solver_center, set_density_matrix_mpi, set_pratio_mpi
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ eigen solver center
!     imode = 1 : standard eigen-value problem with real-symmetrix matrix    (A y = e y )
!     imode = 2 : generalized eigen-value problem with real-symmetrix matrix (A y = e B y)
!
!         mat_a      : input  : input matrix A
!                      output : eigen vectors ( A(:,k) is the k-th eigen vector)
!         eig_levels : input  : empty
!                      output : eigen values e ( e(1) =< e(2) =< ...)
!         mat_b      : input  : input matrix B
!                      output : not preserved
!
  subroutine eig_solver_center(imode, log_unit, SEP_solver, GS_transformation, blocksize, level_low_high, &
&                eig_levels, desc_eigenvectors, eigenvectors)
!
   use elses_mod_md_dat, only : final_iteration
   use M_config, only : config
   use M_ext_matrix_data   
   use M_lib_mpi_wrapper
!  use wp_setting_m
!  use wp_processes_m   
!  use wp_main_aux_m
   use M_wavepacket  ! For testing wavepacket_main_ext().
   implicit none

   integer, intent(in) :: imode
   character(len=*),       intent(in)               :: SEP_solver
   character(len=*),       intent(in)               :: GS_transformation
   integer,                intent(in)               :: blocksize
   integer,                intent(in)               :: level_low_high(2)
   integer,                intent(in)               :: log_unit
   real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
   integer,                intent(out)              :: desc_eigenvectors(9)
   real(DOUBLE_PRECISION), intent(out), allocatable :: eigenvectors(:, :)
!
   
   write(*,*)'@@ eig_solver_center : dummy routine'
!
   write(*,*)'  SEP_solver        =', trim(SEP_solver)
   write(*,*)'  GS_transformation =', trim(GS_transformation)
   write(*,*)'  blocksize         =', blocksize
   write(*,*)'  lowest level      =', level_low_high(1)
   write(*,*)'  highest level     =', level_low_high(2)

   write(*,*)'ERROR:DUMMY routine is called(eig_solver_center)'
   write(*,*)'The external library called Eigen Engine may be required.'
   stop
!
!
  end subroutine eig_solver_center

  subroutine set_pratio_mpi()
    use elses_arr_eig_leg
    write(*,*)'ERROR:DUMMY routine is called(set_pratio_mpi)'
    write(*,*)'The external library called Eigen Engine may be required.'
  end subroutine set_pratio_mpi  

  subroutine set_density_matrix_mpi()
    use M_qm_domain
    use elses_arr_eig_leg
    use elses_mod_orb2
    use elses_mod_js4jsv
    use elses_mod_multi
    use M_lib_mpi_wrapper
    write(*,*)'ERROR:DUMMY routine is called(set_density_matrix_mpi)'
    write(*,*)'The external library called Eigen Engine may be required.'
  end subroutine set_density_matrix_mpi
!
!
end module M_eig_solver_center

