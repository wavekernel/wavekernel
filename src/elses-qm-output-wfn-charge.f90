!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_output_wfn_charge
!
! Module variables : changed 
!   --> ETB (in qm_geno_output)
!
  use M_qm_domain
!
   private
   public output_for_wfn_charge

   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine output_for_wfn_charge
!
    use M_config,               only : config  !(unchanged)
    use M_qm_domain,            only : e_num_on_basis !(unchanged)
    use M_qm_domain,            only : noav, nval, atm_element !(unchanged)
!   use elses_mod_md_dat,       only : final_iteration !(unchanged)
!   use M_io_ctrl_output_eigen, only : set_filename_for_save_wfn ! (routine)
!   use M_output_eigenstates,   only : output_eigenstates
!   use M_output_participation, only : calc_participation
!   use M_qm_solver_eig_geno,   only : plot_spectrum_from_eigen  ! (routine)
    use elses_mod_file_io,      only : vacant_unit               ! (function)
    use M_output_basis,         only : write_basis_info          ! (routine)
    use elses_mod_md_dat,       only : final_iteration             ! (unchanged)
    use elses_mod_md_dat,       only : itemd                       ! (unchanged)
    implicit none
    character(len=100) :: filename_wfn_chg
    character(len=100) :: filename_wrk
    character(len=100) :: chara_wrk
    integer            :: i_v, lu, iunit, lenf
    logical            :: save_now
    integer            :: jsv, ja, j
!
    i_v = config%option%verbose
    lu  = config%calc%distributed%log_unit
    filename_wfn_chg = trim(adjustl(config%output%wfn_charge%filename))
!
!   stop 'STOP=STOP'
!
    if (lu > 0) write(lu,*) '@@@ output_for_wavefunction_charge'
    if (lu > 0) write(lu,*) '  config%output%wfn_charge%set      = ', config%output%wfn_charge%set
    if (lu > 0) write(lu,*) '  config%output%wfn_charge%interval = ', config%output%wfn_charge%interval
!
    if (.not. allocated(e_num_on_basis)) then
      if (i_v >= 1) then
        if (lu > 0) write(lu,*) '@@@ output_for_wavefunction_charge...is skipped.'
      endif
      return
    endif
!
!   stop 'STOP=STOP2'
!
    save_now = .false.
    if (.not. config%output%wfn_charge%set) then 
      save_now = .false.
      if (i_v >= 1) then
        if (lu > 0) write(lu,*) '@@@ output_for_wavefunction_charge...is skipped.'
      endif
      return
    endif
!
!   stop 'STOP=STOP3'
!
    if (config%output%wfn_charge%interval < 0 ) then
      write(*,*)' ERROR! config%output%wfn_charge%interval :', config%output%wfn_charge%interval
      stop
    endif
!
!   stop 'STOP=STOP4'
!
    if (config%output%wfn_charge%interval == 0) then
      filename_wrk=trim(adjustl(filename_wfn_chg)) 
      if ( trim(adjustl(filename_wfn_chg)) == '') then
        save_now = .false.
      else
        if (final_iteration) save_now = .true.
      endif
    else
!     write(*,*)'interval = ', config%output%wfn_charge%interval
      if ( mod(itemd-1,config%output%wfn_charge%interval) ==  0 ) save_now = .true.
      write(chara_wrk, '(i8.8)') config%system%structure%mdstep
      lenf=len_trim(filename_wfn_chg)
      if (filename_wfn_chg(lenf-3:lenf) /= '.txt') then
        write(*,*)'ERROR(output_for_wfn_charge):incorret filename:',trim(adjustl(filename_wfn_chg))
        stop
      endif
      filename_wrk=trim(adjustl(filename_wfn_chg(1:lenf-4)))//'_'//trim(chara_wrk)//'.txt'
    endif
!
    if (.not. save_now) then
      if (i_v >= 1) then
        if (lu > 0) write(lu,*) '@@@ output_for_wavefunction_charge...is skipped.'
      endif
      return
    endif

    if (i_v >= 1) then
      if (lu > 0) write(lu,*) '@@@ output_for_wavefunction_charge:file=', trim(filename_wrk)
    endif
!
    iunit=vacant_unit()
    open(iunit, file=filename_wrk, status='unknown')
!
    call write_basis_info(iunit)
!
    write(iunit,'(a)') '# LCAO coefficients'
    write(iunit,'(a)') ' 1 '
!
    j=0
    do jsv=1,noav
      do ja=1,nval(atm_element(jsv))
        j=j+1
        write(iunit,'(I10,F25.15,a,I10)') j, e_num_on_basis(ja,jsv), '   k= 1'
      enddo
    enddo
!
    close(iunit)
!
    if (i_v >= 1) then
      if (lu > 0) write(lu,*) '  filename_wfn_chg = ',filename_wfn_chg
    endif
!
!   stop 'STOP Manually 1'
!
  end subroutine output_for_wfn_charge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_output_wfn_charge


