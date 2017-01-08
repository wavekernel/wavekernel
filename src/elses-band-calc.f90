!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_band_calc
  use M_qm_domain
  use M_wall_clock_time
  use M_qm_geno, only : set_hamiltonian_and_overlap_geno, &
       write_hamiltonian_and_overlap_geno, &
       set_rest_geno, write_Ecc_perAtom_geno
       
  implicit none
  real(DOUBLE_PRECISION) :: elapse_time
  real(DOUBLE_PRECISION),allocatable :: EccPerAtom(:)
  character(len=*),parameter :: outputFileName="OverlapAndHamiltonian.txt"
  character(len=*),parameter :: outputFileNameEcc="EccPerAtom.txt"
  contains
    subroutine elses_band_calculation
      use M_config
      use M_ini_load_geno
      logical :: BinaryIO = .false.

      if( index(config%option%functionality,"Bin") > 0 ) BinaryIO = .true.
      write(*,*) 'BAND:', config%option%functionality

      call get_system_clock_time(elapse_time)
      write(*,'("*** Entering BAND-CALC Preparation (",ES15.5," sec) ***")') elapse_time
      call ini_load_geno
      ! call elses_ini_set_velocity ! initialization for MD
      call get_system_clock_time(elapse_time)
      call qm_initialization
      allocate(EccPerAtom(noav))
      call set_rest_geno(EccPerAtom)
      call set_hamiltonian_and_overlap_geno
      ! parameter check
      write(*,'("Cluster radius (in atomic unit)=",ES14.5)')config%system%cutoff_radius
      if( i_pbc_x /= 1 .or. i_pbc_y /= 1 .or. i_pbc_z /= 1 ) &
           stop '!!!ERROR!!! To prepare Hamiltonian for band-calculation,&
           & periodic boundary condition should be applied'
      if( config%system%cutoff_radius > ax/2.01d0 ) then
         write(*,'("ax=",ES14.5)')ax
         stop '!!!ERROR!!! Too large cutoff_radius, or, too small system size along x-axis'
      end if
      if( config%system%cutoff_radius > ay/2.01d0 ) then
         write(*,'("ay=",ES14.5)')ay
         stop '!!!ERROR!!! Too large cutoff_radius, or, too small system size along y-axis'
      end if
      if( config%system%cutoff_radius > az/2.01d0 ) then
         write(*,'("az=",ES14.5)')az
         stop '!!!ERROR!!! Too large cutoff_radius, or, too small system size along z-axis'
      end if
      
      call write_hamiltonian_and_overlap_geno(outputFileName,BinaryIO)
      call write_Ecc_perAtom_geno(EccPerAtom,outputFileNameEcc,BinaryIO)
      
      write(*,'("*** Quitting BAND-CALC Preparation (",ES15.5," sec) ***")') elapse_time
      write(*,'("Now overlap & Hamiltonian matrix are written in ",A,".")') &
           outputFileName
      
      contains
        subroutine qm_initialization
          use M_qm_geno, only:initialize_charge
          call qm_domain_setting
          call initialize_charge
        end subroutine qm_initialization
    end subroutine elses_band_calculation
end module M_band_calc
