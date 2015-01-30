module M_ext_atom_info
!
  use M_qm_domain ,  only : i_verbose    !(unchanged)
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  private
!
  public get_num_of_atoms
  public get_basis_info
  public get_name_and_positon_for_atom
!
  contains
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Get the number of atoms
!
  subroutine get_num_of_atoms(n_atom)
    use M_qm_domain ,  only : noav  !(unchanged)
    implicit none
    integer,  intent(out) :: n_atom
!
    n_atom=noav
!
  end subroutine get_num_of_atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Get the number of bases
!
  subroutine get_basis_info(last_element_for_atom)
    use M_qm_domain ,  only : noav, nval, atm_element
    implicit none
    integer,  intent(inout) :: last_element_for_atom(:)
    integer :: jj, jsv
!
    jj=0
    do jsv=1,noav
      jj=jj+nval(atm_element(jsv))
      last_element_for_atom(jsv)=jj
    enddo
!
  end subroutine get_basis_info
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Get the element name and position for atom
!    The position is given in atomic unit.
!
  subroutine get_name_and_positon_for_atom(jsv, elem_name, position)
    use M_qm_domain ,  only : noav, nval, atm_element, atm_position, element_name, ax, ay, az
    implicit none
    integer,  intent(in)                :: jsv
    character(len=*), intent(inout)     :: elem_name
    real(DOUBLE_PRECISION), intent(out) :: position(3)
!
    position(1)=atm_position(1,jsv)*ax
    position(2)=atm_position(2,jsv)*ay
    position(3)=atm_position(3,jsv)*az
!
!   write(*,*)' atom_element =', atm_element(jsv)
!   write(*,*)' element_name =', trim(element_name(atm_element(jsv)))
!
    elem_name=trim(element_name(atm_element(jsv)))
!    
  end subroutine get_name_and_positon_for_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_ext_atom_info



