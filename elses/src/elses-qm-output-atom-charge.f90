!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_output_atom_charge

  use M_qm_domain,        only : i_verbose !(unchanged)
  implicit none
  integer :: unit_num
! character(len=*), parameter :: filename="output_atom_charge.txt"
!
  public :: output_atom_charge
!
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Main routine
!
  subroutine output_atom_charge
!
    use elses_mod_phys_const, only : ev4au !(parameter)
    use M_qm_domain,    only : c_system !(unchanged)
    use  M_qm_domain,   only : mulliken !(routine)
    use  M_qm_domain,   only : noav, njsd, atm_element, nval, &
&                        e_num_on_atom,  e_num_on_basis, element_name !(unchanged)
    use M_qm_geno_Huckel_atom_params, only : GetInitialENum !(routine)
    use M_config,    only : config !(unchanged)
    use M_qm_geno_output, only : qm_calc_ecsc_list !(routine)
!
    implicit none
    integer :: atm_index
    integer :: ict, ierr, nval_wrk, nval_max
    integer :: step_count
    logical :: flag_for_init
    real(kind=8), allocatable ::  initial_e_num(:)
    real(kind=8), allocatable ::  d_e_num_on_atom(:)
    real(kind=8)              ::  e_num_on_basis_wrk(9)  ! for s, p, d bases
    real(kind=8), allocatable ::  atom_csc_energy_wrk(:,:)
    real(kind=8) :: d_sum
    character(len=100)        ::  filename_wrk
    integer :: lu, n_csc_loop
!
    ict=1
    lu = config%calc%distributed%log_unit
    n_csc_loop = config%calc%genoOption%CSC_max_loop_count
!
    if (c_system /= 'geno') then
      if (i_verbose >= 1) then 
        if (lu > 0) write(lu,*)'@@ output_atom_charge.. is skipped'
        return
      endif
    endif   
!
    nval_max=maxval(nval)
    step_count=config%system%structure%mdstep
!
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*) '@@ output_atom_charge'
    endif   
!
    if ((nval_max > 9) .or. (nval_max <1)) then
      write(*,*)'ERROR:nval_max=',nval_max
      stop
    endif   
!
    call output_setting(unit_num, flag_for_init,filename_wrk)
!
    if (unit_num == -1) then
      if (i_verbose >= 1) then
        if (lu > 0) write(lu,*) '....no file is generated for atom charge.'
      endif   
      return
    endif  
!
!   call mulliken
!     ----> calculate e_num_on_atom(:),  e_num_on_basis(:)
!
    allocate (initial_e_num(noav), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error'
      stop
    endif   
    initial_e_num(:)=0.0d0
!
    allocate (d_e_num_on_atom(noav), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error'
      stop
    endif   
    d_e_num_on_atom(:)=0.0d0
!
    allocate (atom_csc_energy_wrk(noav,2), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error'
      stop
    endif   
    atom_csc_energy_wrk(:,:)=0.0d0
!
    call GetInitialENum(initial_e_num)
!     ----> calculate the initial charge of the atoms
!
    d_e_num_on_atom(:)=e_num_on_atom(:)-initial_e_num(:)
!     ----> deviation from the initial charge
!
    if (n_csc_loop > 0) then
      if (i_verbose >= 1) then
        if (lu > 0) write(lu,'(a)') 'INFO-CSC:call qm_calc_ecsc_list'
      endif
      call qm_calc_ecsc_list(atom_csc_energy_wrk)
      atom_csc_energy_wrk(:,:)=atom_csc_energy_wrk(:,:)*ev4au
    endif
!
    if (lu > 0) then 
      write(lu,'(a,2i10,f10.5)')'INFO:Max of elec_num ( n - n_ini )=', &
&         step_count, maxloc(d_e_num_on_atom), maxval(d_e_num_on_atom)
      write(lu,'(a,2i10,f10.5)')'INFO:Min of elec_num ( n - n_ini )=', &
&         step_count, minloc(d_e_num_on_atom), minval(d_e_num_on_atom)
    endif
!
    d_sum=sum( d_e_num_on_atom(:) )
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,'(a)') 'INFO:The following quantities should be zero (within the bisection method accuracy)'
      if (lu > 0) write(lu,'(a,f30.20)') '  sum of d_e_num_on_atom (should be zero) =',d_sum 
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
    do atm_index=1,noav
       e_num_on_basis_wrk(:)=0.0d0
       nval_wrk=nval(atm_element(atm_index))
       e_num_on_basis_wrk(1:nval_wrk)=e_num_on_basis(1:nval_wrk,atm_index)
       if (nval_max <= 4) then
       write(unit_num,'(i7,a6,f10.5,f15.10,5f10.5,2i10,2f17.8)') &
&            atm_index, trim(element_name(atm_element(atm_index))), &
&            initial_e_num(atm_index),  d_e_num_on_atom(atm_index), &
&            e_num_on_basis_wrk(1), sum(e_num_on_basis_wrk(2:4)), &
&            e_num_on_basis_wrk(2:4),  njsd(atm_index,ict), step_count, &
&            atom_csc_energy_wrk(atm_index,1:2)
       else
       write(unit_num,'(i7,a6,f10.5,f15.10,11f10.5,2i10,2f17.8)') &
&            atm_index, trim(element_name(atm_element(atm_index))), &
&            initial_e_num(atm_index),  d_e_num_on_atom(atm_index), &
&            e_num_on_basis_wrk(1), sum(e_num_on_basis_wrk(2:4)), sum(e_num_on_basis_wrk(5:9)), &
&            e_num_on_basis_wrk(2:9), njsd(atm_index,ict), step_count, &
&            atom_csc_energy_wrk(atm_index,1:2)
       endif   
    enddo
!
    close(unit_num)
!
    deallocate (initial_e_num, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Dealloc. error'
      stop
    endif   
!
    deallocate (d_e_num_on_atom, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Dealloc. error'
      stop
    endif   
!
    deallocate (atom_csc_energy_wrk, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Dealloc. error(atom_csc_energy_wrk)'
      stop
    endif   
!
  end subroutine output_atom_charge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set the unit number (iunit), if the output file will be generated.
!  Otherwise iunit is set to be -1.
!
  subroutine output_setting(iunit,flag_for_init,filename_wrk)
!
    use M_config, only : config                     !(unchanged)
    use elses_mod_file_io,      only : vacant_unit  !(function)
    use elses_mod_md_dat,       only : itemd        !(unchanged)
    implicit none
    integer, intent(out) :: iunit
    logical, intent(out) :: flag_for_init
    character(len=*), intent(out) :: filename_wrk
!
    iunit=-1
    if (itemd == 1) then
      flag_for_init = .true.
    else
      flag_for_init = .false.
    endif  
!
    if ( .not. config%output%atom_charge%set ) return
!
    if (mod(itemd-1,config%output%atom_charge%interval) /= 0 ) return
!
    filename_wrk=trim(config%output%atom_charge%filename)
!
    iunit=vacant_unit()
!
  end subroutine output_setting
!
end module M_output_atom_charge
