!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module M_qm_population_ortho
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_lib_dst_info,     only : log_unit !(unchanged)
!
  implicit none
!
  private
!
! Public routines
!
  public :: save_population_ortho
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Save population in the orthogonal-system cases
!
  subroutine save_population_ortho
!
    use M_config  !(CHANGED : config%system%structure%vatom(:)%population)
    use M_qm_domain, only : c_system !(unchanged)
    use M_qm_domain, only : noav !(unchanged)
    use M_qm_domain, only : dbij, jsv4jsd !(unchanged)
!   use M_qm_domain, only : e_num_on_atom, previous_e_num_on_atom !(CHANGED)
    use M_structure, only : atom_type !(type)
    implicit none
    type(atom_type), pointer :: atom
    real(DOUBLE_PRECISION)   :: ddd
    integer :: atm_index, orb_index
    integer :: jsd1
    integer :: i_show
!
    i_show=3
!
!
    if ((c_system == 'C_Xu') .or. (c_system == 'Si_Kwon') .or. (c_system == 'GaAs.Molteni.1993')) then
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,*) '@@ save_population_ortho'
      endif  
    else
      write(*,*)'Error(save_population_ortho):c_system=',c_system
      stop
    endif   
!
    if (.not. allocated(dbij)) then
      write(*,*)'Error(save_population_ortho):dbij is not allocated'
      stop
    endif
!

!
    do atm_index=1, noav
      atom => config%system%structure%vatom(atm_index)
      atom%population_set       = .true.
      jsd1=1
      if (jsv4jsd(jsd1,atm_index) /= atm_index) then
        write(*,*)'ERROR(save_population_ortho):atm_index, jsv4jsd=',atm_index, jsv4jsd(jsd1,atm_index )
        stop
      endif   
      ddd=0.0d0
      do orb_index=1,size(dbij,1)
        ddd=ddd+2.0d0*dbij(orb_index, orb_index, jsd1, atm_index) 
      enddo  
      atom%population=ddd
      if (atm_index <= i_show) then
        if (log_unit > 0) write(log_unit,*) ' save popu=', atm_index, atom%population
      endif
      if ((atom%population < -100.d00) .or. (atom%population > 100.d0)) then
        write(*,*)'Error(population)=', atm_index, atom%population
        stop
      endif
    enddo
!   stop 'stop manually'
!
  end subroutine save_population_ortho
!
end module M_qm_population_ortho


