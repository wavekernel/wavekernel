!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
module M_output_atom_energy
!
! Module variables : changed 
!   --> ETB (in qm_geno_output)
!
   private
!
! Public routines
   public output_atom_energy
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Output for atom energy
!
!
  subroutine output_atom_energy
!
!
    use M_qm_domain                              !(unchanged)
    use M_config,    only : config               !(unchanged)
    use elses_mod_phys_const, only : ev4au       !(parameter)
    use M_qm_geno_output,     only : qm_calc_etb_atom, qm_calc_ecsc_atom !(routines)
    use M_qm_geno,            only : calc_e_rep_atom                     !(routines)
    use M_md_velocity_routines, only : calc_kin_ene_atom          !(routines)
    implicit none
    integer      :: i
    integer, parameter  :: ict4h=1
    integer, parameter  :: number_of_items = 4
    real(kind=8) :: value_of_etb
    integer      :: jsv2, nss2, nval2, ja2, jsd1
    integer      :: jsv1, nss1, nval1, ja1, ierr, iplot
    real(kind=8) :: ddsum, ddd1, ddd2, ddd3
    real(kind=8), allocatable :: atom_energy(:,:)
    real(kind=8), allocatable :: atom_energy_sum(:,:) 
    integer                   :: unit_num
    integer                   :: step_count
    logical                   :: flag_for_init
    character(len=100)        :: filename_wrk
!
    step_count=config%system%structure%mdstep
!
    call output_setting(unit_num, flag_for_init, filename_wrk)
!
    if (unit_num == -1) then
      if (i_verbose >= 1) then
        write(*,*)'No file is generated for atom energy.'
      endif   
      return
    endif
!
    if (i_verbose >= 1) then
      write(*,*)'@@ output_atom_energy '
    endif
!
    allocate (atom_energy(noav, number_of_items), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. Error (output_atom_energy)'
      stop
    endif   
    atom_energy(:,:)=0.0d0
!
    allocate (atom_energy_sum(nos, number_of_items), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. Error (output_atom_energy)'
      stop
    endif   
    atom_energy_sum(:,:)=0.0d0
!
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      call qm_calc_etb_atom(jsv2, atom_energy(jsv2,1))
      call qm_calc_ecsc_atom(jsv2, atom_energy(jsv2,2))
      call calc_e_rep_atom(jsv2, atom_energy(jsv2,3))
      call calc_kin_ene_atom(jsv2, atom_energy(jsv2,4))
    enddo  
!
    if (flag_for_init) then
       open(unit_num, file=filename_wrk)
    else
       open(unit_num, file=filename_wrk, position="append") 
    endif
!
!   write(unit_num,*) noav
!   write(unit_num,*) " mdstep=", step_count
!
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      write(unit_num,'(2i10,a,a2,a,i10, 4f20.10)') step_count, jsv2, ' ', trim(element_name(nss2)), & 
&              ' ', njsd(jsv2,ict4h), atom_energy(jsv2,1:4)
    enddo  
!
    do i=1,number_of_items
      do jsv2=1,noav
        nss2=atm_element(jsv2)
        nval2=nval(nss2)
        atom_energy_sum(nss2,i)=atom_energy_sum(nss2,i)+atom_energy(jsv2,i)
      enddo  
    enddo  
!
    do nss2=1,nos
      write(*,'(a,i10,a,a2,a,4f20.10)')'atom species, atom_energy_sum =',nss2, & 
&         ' ', trim(element_name(nss2)), ' ', atom_energy_sum(nss2,1:4)
    enddo
!
    close(unit_num)
!
  end subroutine output_atom_energy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ Set the unit number (iunit), if the output file will be generated.
!     Otherwise iunit is set to be -1.
!
  subroutine output_setting(iunit,flag_for_init,filename_wrk)
!
    use M_config, only : config                     !(unchanged)
    use elses_mod_file_io,      only : vacant_unit  !(function)
    use elses_mod_md_dat,       only : itemd        !(unchanged)
    implicit none
    integer,          intent(out) :: iunit
    logical,          intent(out) :: flag_for_init
    character(len=*), intent(out) :: filename_wrk
!
    iunit=-1
    if (itemd == 1) then
      flag_for_init = .true.
    else
      flag_for_init = .false.
    endif  
!
    if ( .not. config%output%atom_energy%set ) return
!
    if (mod(itemd-1,config%output%atom_energy%interval) /= 0 ) return
!
    filename_wrk=trim(config%output%atom_energy%filename)
!   write(*,*)'filename (atom_energy) =',filename_wrk
!
    iunit=vacant_unit()
!
  end subroutine output_setting
!
!
end module M_output_atom_energy



