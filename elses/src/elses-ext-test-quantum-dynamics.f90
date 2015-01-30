module M_ext_test_qd
!
  use M_qm_domain ,  only : i_verbose    !(unchanged)
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  private
!
  public test_for_quantum_dynamics
!
  contains
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ test routine for quantum dynamics
!
  subroutine test_for_quantum_dynamics
!
    use M_qm_domain ,  only : i_pbc_x, i_pbc_y, i_pbc_z  !(unchanged)
!
    use M_ext_matrix_data, only: set_matrix_data  !(routine)
    use M_ext_matrix_data, only: plot_matrix_data !(routine)
    use M_ext_atom_info,   only: get_num_of_atoms !(routine)
    use M_ext_atom_info,   only: get_basis_info   !(routine)
    use M_ext_atom_info,   only: get_name_and_positon_for_atom !(routine)
!
    implicit none
    integer :: ierr
    integer :: n_atom
    integer, allocatable :: last_element_for_atom(:)
!
    integer :: j
    character(len=8)       :: element_name
    real(DOUBLE_PRECISION) :: position(3)
!
    if (i_verbose >= 1) then
      write(*,*)'@@ test_for_quantum_dynamics'
    endif 
!
    ierr=i_pbc_x*i_pbc_y*i_pbc_z
    if (ierr /=0) then
      write(*,*)'ERROR(test_for_quantum_dynamics):Periodic case is not supported'
      write(*,*)'ERROR(test_for_quantum_dynamics): for quantum dynamics calculation now'
      stop
    endif
!
    call get_num_of_atoms(n_atom)
    if (i_verbose >= 1) then
      write(*,*)' number of atoms =', n_atom
    endif 
!
    allocate(last_element_for_atom(n_atom), stat=ierr)
!
    call get_basis_info(last_element_for_atom)
!
    do j=1,n_atom
      call get_name_and_positon_for_atom(j, element_name, position)
      write(*,*) ' atom information : j=', j
      write(*,*)'    element name              =', j, trim(element_name)
      write(*,*)'    position (au)             =', j, position(:)
      write(*,*)'    last element for the atom =', j, last_element_for_atom(j)
    enddo
!
    call set_matrix_data
!     --> Set matrix data array for H and S 
!          as sparse matrix form
!
    call plot_matrix_data
!     --> Plot the matrix data of H and S 
!             as the MatrixMarket style format
!
  end subroutine test_for_quantum_dynamics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_ext_test_qd


