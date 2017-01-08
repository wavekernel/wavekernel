!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_output_atom_charge_ortho

  use M_qm_domain,        only : i_verbose !(unchanged)
  implicit none
! character(len=*), parameter :: filename="output_atom_charge.txt"
!
  public :: output_atom_charge_ortho
!
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Main routine
!
  subroutine output_atom_charge_ortho(unit_num, flag_for_init, filename_wrk)
!
    use M_qm_domain, only : c_system !(unchanged)
    use M_qm_domain, only : noav, njsd, atm_element, nval, jsv4jsd, element_name !(unchanged)
    use M_qm_domain, only : dbij   !(unchanged)
    use M_config,    only : config !(unchanged)
!
    implicit none
    integer,             intent(in) :: unit_num
    logical,             intent(in) :: flag_for_init
    character(len=*),    intent(in) :: filename_wrk

    integer :: atm_index, orb_index
    integer :: ict, ierr, nval_wrk, nval_max
    integer :: step_count
    integer :: jsd1
!
    real(kind=8), allocatable ::  d_e_num_on_atom(:)
    real(kind=8)              ::  e_num_on_basis_wrk(9)  ! for s, p, d bases
    real(kind=8) :: d_sum
    real(kind=8) :: e_neutral_charge
!
    ict=1
!
    if (c_system /= 'C_Xu') then
      if (i_verbose >= 1) write(*,*)'@@ output_atom_charge_ortho.. is skipped:c_system=',c_system
      return
    endif   
!
    if (i_verbose >= 1) write(*,*)'@@ output_atom_charge_ortho:c_system=',c_system
!
    nval_max=maxval(nval)
    step_count=config%system%structure%mdstep
!
    e_neutral_charge=4.0d0
!
    allocate (d_e_num_on_atom(noav), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error'
      stop
    endif   
!
    if (flag_for_init) then
       open(unit_num, file=filename_wrk)
    else
       open(unit_num, file=filename_wrk, position="append") 
    endif
!
    write(unit_num,*) noav
    write(unit_num,*) " mdstep=", step_count
!
    do atm_index=1,noav
      jsd1=1 
      if (jsv4jsd(jsd1,atm_index) /= atm_index) then
        write(*,*)'ERROR(save_population_ortho):atm_index, jsv4jsd=',atm_index, jsv4jsd(jsd1,atm_index ) 
        stop
      endif
      e_num_on_basis_wrk(:)=0.0d0
      do orb_index=1,size(dbij,1)
        e_num_on_basis_wrk(orb_index)=2.0d0*dbij(orb_index, orb_index, jsd1, atm_index)
      enddo   
      d_e_num_on_atom(atm_index)=sum(e_num_on_basis_wrk)-e_neutral_charge
      write(unit_num,'(i7,a6,f10.5,f15.10,5f10.5,i10)') &
&            atm_index, trim(element_name(atm_element(atm_index))), &
&            e_neutral_charge,  d_e_num_on_atom(atm_index), &
&            e_num_on_basis_wrk(1), sum(e_num_on_basis_wrk(2:4)), &
&            e_num_on_basis_wrk(2:4),  njsd(atm_index,ict)
    enddo
!
    close(unit_num)
!
    write(*,'(a,2i10,f10.5)')'INFO:Max of elec_num ( n - n_ini )=', &
&         step_count, maxloc(d_e_num_on_atom), maxval(d_e_num_on_atom)
    write(*,'(a,2i10,f10.5)')'INFO:Min of elec_num ( n - n_ini )=', &
&         step_count, minloc(d_e_num_on_atom), minval(d_e_num_on_atom)
!
    d_sum=sum( d_e_num_on_atom(:) )
!
    if (i_verbose >= 1) then
      write(*,'(a,f30.20)')     '  sum of d_e_num_on_atom (should be zero) =',d_sum 
    endif   
!
    deallocate (d_e_num_on_atom, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Dealloc. error'
      stop
    endif   
!
  end subroutine output_atom_charge_ortho
!
end module M_output_atom_charge_ortho

