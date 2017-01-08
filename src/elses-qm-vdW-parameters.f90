module M_qm_vdW_parameters
  use M_config, only: config
  use flib_dom
  use M_xml_get_value_with_unit, only : xml_get_value_w_unit        ! routine
  use M_xml_get_value_with_unit, only : xml_get_value_w_unit_detail ! routine
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  private
  public set_vdW_parameters_db
  public set_vdW_parameters_geno
  public show_vdW_parameters_db
!
  contains
!  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Show the vdW parameters from the daba base
!
  subroutine show_vdW_parameters_db
    use elses_mod_file_io,        only : vacant_unit                 !(function)
    use M_lib_element_database,   only : get_element_name            ! routine
    use M_lib_element_database,   only : get_atomic_num              ! routine
    use M_lib_phys_const,         only : convert_unit                ! routine
    use M_lib_db_ion_ene,         only : get_ion_ene_db_unit         ! routine
    use M_lib_db_polarization,    only : get_polarization_db_unit    ! routien
    use M_lib_db_covalent_radius, only : get_covalent_radius_db_unit ! routine
    implicit none
    character(len=8)       :: element_name
    real(DOUBLE_PRECISION) :: ionization_energy
    real(DOUBLE_PRECISION) :: dipole_polarizability
    real(DOUBLE_PRECISION) :: covalent_radius
    character(len=256)     :: unit_name1, unit_name2, unit_name3
    integer                :: i_verbose, lu, ierr, iunit
    integer                :: atomic_num, atomic_num_tmp
!
    integer, parameter          :: atomic_num_max = 86
    character(len=*), parameter :: filename_wrk='output_vdW_parameter_db.txt'
!
    i_verbose = config%option%verbose
    lu        = config%calc%distributed%log_unit
    iunit     = vacant_unit()
!
    open(iunit, file=filename_wrk, status='unknown')
!
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*)'@@ show_vdW_parameters_db:output file=', trim(filename_wrk)
    endif  
!
    do atomic_num=1, atomic_num_max
      call get_element_name(atomic_num, element_name) 
      call get_atomic_num(trim(element_name), atomic_num_tmp)
      if (atomic_num /= atomic_num_tmp) then
        write(*,*)'ERROR(show_vdW_parameters_db):atomic_num=', atomic_num, atomic_num_tmp
      endif   
      call get_ion_ene_db_unit(atomic_num, ionization_energy, unit_name1)
      call get_polarization_db_unit(atomic_num, dipole_polarizability, unit_name2)
      call get_covalent_radius_db_unit(atomic_num, covalent_radius, unit_name3)
      if (atomic_num == 1) then
        write(iunit,'(a)')   '# ELSES: Data base for vdW parameters'
        write(iunit,'(a)')   '#   Ref. CRC Handbook of Chemistry and Physics, 94th Edition'
        write(iunit,'(a)')   '#     ed. William M. Haynes, CRC Press, Boca Raton (2013); ISBN:978-1466571143'
        write(iunit,'(a,a)') '#   The unit of ionization_energy     : ', trim(unit_name1)
        write(iunit,'(a,a)') '#   The unit of dipole_polarizability : 10^{-24} ', trim(unit_name2)
        write(iunit,'(a,a)') '#   The unit of covalent_radius       : ', trim(unit_name3)
        write(iunit,'(a)')   '#'
        write(iunit,'(a)')   '#     Note: The ionization energy of At is a dummy value (-1.0),' 
        write(iunit,'(a)')   '#           since the value is missing in the book. '
        write(iunit,'(a)')   '#'
        write(iunit,'(a)')   '#                  ion_ene         polariz    cov_rad'
      endif
      write(iunit,'(i5, a8, f15.7, f15.7, f10.4)') atomic_num, trim(element_name), & 
&           ionization_energy, dipole_polarizability*1d24, covalent_radius
    enddo
!   
    close(iunit)
!
!   stop 'STOP manually'
!
  end subroutine show_vdW_parameters_db
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set vdW parameters from the daba base
!
  subroutine set_vdW_parameters_db(element_name, ionization_energy, dipole_polarizability, covalent_radius)
    use M_lib_element_database,   only : get_atomic_num              ! routine
    use M_lib_phys_const,         only : convert_unit                ! routine
    use M_lib_db_ion_ene,         only : get_ion_ene_db_unit         ! routine
    use M_lib_db_polarization,    only : get_polarization_db_unit    ! routien
    use M_lib_db_covalent_radius, only : get_covalent_radius_db_unit ! routine
    implicit none
!   type(fnode), pointer     :: node_input
!   type(fnode), pointer     :: node_wrk
    character(len=*),       intent(in)  :: element_name
    real(DOUBLE_PRECISION), intent(out) :: ionization_energy
    real(DOUBLE_PRECISION), intent(out) :: dipole_polarizability
    real(DOUBLE_PRECISION), intent(out) :: covalent_radius
    character(len=256)       :: unit_name
!   character(len=256)       :: value_chara, unit_name
!   character(len=256)       :: comment
    real(DOUBLE_PRECISION)   :: value_tmp
    integer                  :: i_verbose, lu, ierr
    integer                  :: atomic_num
    integer, parameter       :: atomic_num_max = 84
!
    i_verbose = config%option%verbose
    lu        = config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*)'INFO-vdW:element name=',trim(element_name)
    endif  
!
    call get_atomic_num(trim(element_name), atomic_num)
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'INFO-vdW:atomic_num=', atomic_num
    endif  
!
    if (atomic_num > atomic_num_max) then
      write(*,*)'ERROR(set_vdW_parameters_db):unsuppoted atomic number: atomic_num = ', atomic_num
      stop
    endif   
!
    unit_name =''                    ! dummy value
    ionization_energy     = -1.0d0   ! dummy value
    dipole_polarizability = -1.0d0   ! dummy value
    covalent_radius       = -1.0d0   ! dummy value
!
    call get_ion_ene_db_unit(atomic_num, value_tmp, unit_name)
    ionization_energy = value_tmp*convert_unit('from',trim(unit_name))  ! value in au
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'INFO-vdW:ionization_energy unit =', value_tmp, trim(unit_name)
      if (lu > 0) write(lu,*) 'INFO-vdW:ionization_energy [au] =', ionization_energy
    endif  
!
    call get_polarization_db_unit(atomic_num, value_tmp, unit_name)
    dipole_polarizability = value_tmp*convert_unit('from',trim(unit_name))  ! value in au 
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'INFO-vdW:dipole_polarizability unit =', value_tmp, trim(unit_name)
      if (lu > 0) write(lu,*) 'INFO-vdW:dipole_polarizability [au] =', dipole_polarizability
    endif  
!
    call get_covalent_radius_db_unit(atomic_num, value_tmp, unit_name)
    covalent_radius = value_tmp*convert_unit('from',trim(unit_name))  ! value in au
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'INFO-vdW:covalent_radius unit =', value_tmp, trim(unit_name)
      if (lu > 0) write(lu,*) 'INFO-vdW:covalent_radius [au] =', covalent_radius
    endif  
!
  end subroutine set_vdW_parameters_db
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set vdW parameters from the element XML files
!
  subroutine set_vdW_parameters_geno(node_input, ionization_energy, dipole_polarizability, covalent_radius)
    implicit none
    type(fnode), pointer     :: node_input
    type(fnode), pointer     :: node_wrk
    real(DOUBLE_PRECISION), intent(out) :: ionization_energy
    real(DOUBLE_PRECISION), intent(out) :: dipole_polarizability
    real(DOUBLE_PRECISION), intent(out) :: covalent_radius
    character(len=256)       :: value_chara, unit_name
    character(len=256)       :: comment
    real(DOUBLE_PRECISION)   :: value, value_org
!   logical, parameter       :: debug_mode = .true.
    integer                  :: i_verbose, lu, ierr
!
    i_verbose = config%option%verbose
    lu        = config%calc%distributed%log_unit
!
    ionization_energy     = -1.0d0   ! dummy value
    dipole_polarizability = -1.0d0   ! dummy value
    covalent_radius       = -1.0d0   ! dummy value
!
    comment=''
    unit_name=''
!
!   call xml_get_value_w_unit(node_input, 'ionization_energy', value, comment)
!   write(*,*)' 0:ionization energy [au] =', value
!
    call xml_get_value_w_unit_detail(node_input, 'ionization_energy', value, value_org, unit_name, comment)
    if (i_verbose >= 1) then
      if (lu > 0) then 
        write(lu,*)'INFO-XML:vdW:ionization energy, unit name =', value_org, ' ', trim(unit_name)
        write(lu,*)'INFO-XML:vdW:ionization energy [au]       =', value
      endif  
    endif   
    ionization_energy = value
!
    call xml_get_value_w_unit_detail(node_input, 'dipole_polarizability', value, value_org, unit_name, comment)
    if (i_verbose >= 1) then
      if (lu > 0) then 
        write(lu,*)'INFO-XML:vdW:dipole_polarizability, unit name =', value_org, ' ', trim(unit_name)
        write(lu,*)'INFO-XML:vdW:dipole_polarizability [au]       =', value
      endif  
    endif   
    dipole_polarizability = value
!
    call xml_get_value_w_unit_detail(node_input, 'covalent_radius', value, value_org, unit_name, comment)
    if (i_verbose >= 1) then
      if (lu > 0) then 
        write(lu,*)'INFO-XML:vdW:covalent_radius, unit name =', value_org, ' ', trim(unit_name)
        write(lu,*)'INFO-XML:vdW:covalent_radius [au]       =', value
      endif  
    endif   
    covalent_radius = value
!
    ierr=0
    if (ionization_energy < 0.0d0) ierr=1
    if (dipole_polarizability < 0.0d0) ierr=1
    if (covalent_radius < 0.0d0) ierr=1
    if (ierr == 1) then
      write(*,*)'ERROR!(set_vdW_parameters_geno)'
      write(*,*)' ionization energy =', ionization_energy
      write(*,*)' dipole_polarizability =', dipole_polarizability
      write(*,*)' covalent_radius =', covalent_radius
      stop
    endif   
!
  end subroutine set_vdW_parameters_geno
!  
end module M_qm_vdW_parameters
