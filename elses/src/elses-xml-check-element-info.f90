!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_check_xml_element_info
!
  use M_config,           only: config    !(unchanged)
!
  implicit none
  integer  :: i_verbose, log_unit
!
  private
  public check_xml_element_info
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine check_xml_element_info
!
    use M_qm_geno_Huckel_atom_params, only: AtomParameter
    implicit none
    integer :: lu
    integer :: i, ierr
!
    i_verbose=config%option%verbose
    lu  = config%calc%distributed%log_unit
!
    if (i_verbose > 0) then
      if (lu>0) write(lu,*) '@@ check_xml_element_info'
    endif   
!
    if (config%calc%genoOption%vanDerWaals) then
      do i=1, size(AtomParameter,1)
        ierr=0
        if (AtomParameter(i)%ionization_energy < 0.0d0) ierr=1
        if (AtomParameter(i)%dipole_polarizability < 0.0d0) ierr=1
        if (AtomParameter(i)%vdW_radius < 0.0d0) ierr=1
        if (ierr == 1) then
          write(*,*)'ERROR(check_xml_element_info ):i=',i
          write(*,*)'ionization_energy     =', AtomParameter(i)%ionization_energy
          write(*,*)'dipole_polarizability =', AtomParameter(i)%dipole_polarizability
          write(*,*)'vdW_radius            =', AtomParameter(i)%vdW_radius
          stop
        endif   
      enddo   
    endif   
!
  end subroutine check_xml_element_info
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module M_check_xml_element_info




