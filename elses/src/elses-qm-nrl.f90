!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_nrl
!
  use elses_mod_ctrl,       only : i_verbose
  implicit none
!
  private
  public qm_nrl_ini
  public qm_nrl_init_alloc
  public qm_nrl_set_ham
  public qm_nrl_set_force
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Read the parameters
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine qm_nrl_ini(rcut_ham)
    use M_config
    use param_ham_nrl,    only: rcut_ham_nrl
    use elses_mod_noav,    only: noav
    use MNrlAdaptor,      only: NRLReadTBParameters
    implicit none
    integer :: i
    real(8) :: rNelec, rcut_ham
    character(64) :: filename 
!
    write(6,*)'@@ LSES_NRL_INI'
!
    if (noav .le. 0) then
      write(6,*)'Note(LSES_NRL_INI):NOAV=',NOAV
    endif
!
    do i = 1, config%system%structure%nelement
       filename = config%system%structure%velement(i)%filename
       call NRLReadTBParameters(filename, rNelec)
    end do
!     ----> rNelec is dummy.
!
    rcut_ham=rcut_ham_nrl
    write(6,*)'rcut_ham=',rcut_ham
!
  end subroutine qm_nrl_ini
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Initial Alloc.for LSES for NRL-leg case
!        ----> allocation of dsij, dpij
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine qm_nrl_init_alloc
    use elses_mod_sim_cell, only : noa
    use elses_mod_noav,     only : noav,noao,noab,noas,noab0,nncut
    use elses_mod_orb1,     only : nvl
!       NOTE: ALL the above quantities are input parameters.
!
    implicit none
    integer noa_def,noav_def,noao_def,noas_def,nvl_def,nncut_def,imode
    real*8 ddd
!
    write(*,*)'@@ LSES_INIT_ALLOC_NRL_LEG:NOA=',noa
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set parameters for matrix size
!
    write(*,*)'   NOAV =',noav
    write(*,*)'   NOAO =',noao
    write(*,*)'   NOAS =',noas
    write(*,*)'   NOAB =',noab
    write(*,*)'   NOAB0=',noab0
    write(*,*)'   NVL  =',nvl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Check the input parameters
!
    ddd=dble(noav)*dble(noao)*dble(noas)*dble(noab)*dble(noab0)*dble(nvl)
!
    if (dabs(ddd) .lt. 1.0d-10) then
      write(6,*)'ERROR!:LSES_INIT_ALLOC2B'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Set working parameters 
!
    noa_def =noa
    noav_def=noav
    noao_def=noao
    noas_def=noas
    nncut_def=nncut
    nvl_def=nvl
!
!
    imode=0
    call elses_alloc_dsij_dpij(nvl_def,noas_def,noav_def,imode)
!        ---->  alloc : dsij(nvl,nvl,noas,noav)
!                     : dpij(nvl,nvl,noas,noav)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Initial setting
!
    call elses_set_js4jsv
!     ---->  set : js4jsv(noav), jsv4js(noa)
!
  end subroutine qm_nrl_init_alloc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Build Hamiltonian for 'NRL_LEG' case
!
   !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine qm_nrl_set_ham
!
    use elses_mod_sel_sys,    only : c_system
!
    use MNrlAdaptor, only : NRLBuildMatrices
    use elses_arr_dhij,       only : dhij
    use elses_arr_dsij,       only : dsij
!
    use elses_mod_tx,       only : jsei
    use elses_mod_js4jsv,   only : js4jsv 
    use elses_mod_orb1,     only : nval
    use elses_mod_jsv4jsd,  only : jsv4jsd,njsd
    use elses_mod_noav,     only : noav
    use elses_mod_orb1,     only : nvl
    use elses_mod_multi,    only : ict4h
    implicit none
    integer :: ishow, jsv2, njsd2, js2, nss2, nval2, jsd1, jsv1, js1
    integer :: nss1, nval1, ja2, ja1
!
    ishow=5
    write(6,*)'@@ LSES_SET_HAMI_NRL_LEG'
    write(6,*)'   --> c_system=',c_system
    if (c_system .ne. "NRL_leg") then
      write(6,*)'  : unsupported:c_system=',c_system
      stop
    endif   
!
    call NRLBuildMatrices(dhij, dsij)
!
    if (i_verbose .ge. 10) then
      do jsv2=1,noav
        njsd2=njsd(jsv2,ict4h)
        js2=js4jsv(jsv2)
        nss2=jsei(js2)
        nval2=nval(nss2)
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
          write(6,*)' js2, jsv2   =',js2,jsv2
          write(6,*)'  njsd2      =',njsd2
          write(6,*)' nss2, nval2 =',nss2,nval2
          write(6,*)' js1, jsv1   =',js1,jsv1
          write(6,*)' nss1, nval1 =',nss1,nval1
          do ja2=1,nval2
            do ja1=1,nval1
               if (jsv2 .le. ishow) then
                 write(6,*)'js1,js2,ja1,ja2,H=', &
&                  js1,js2,ja1,ja2,dhij(ja1,ja2,jsd1,jsv2)
               endif
            enddo
          enddo  
        enddo   
      enddo   
    endif  
!
    write(6,*)'.... is ended without error :LSES_SET_HAMI_NRL_LEG'
!
  end subroutine qm_nrl_set_ham
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Build Hamiltonian for 'NRL_LEG' case
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine qm_nrl_set_force
!
    use elses_mod_sel_sys,    only : c_system
!
    use MNrlAdaptor, only : NrlBuildForces
    use elses_arr_dbij,       only : dbij
    use elses_arr_dpij,       only : dpij
    use elses_mod_foi,        only : foi
    use elses_mod_noav,       only : noav
!
    implicit none
    integer :: jsv, ishow
!
    ishow=5
    write(6,*)'@@ LSES_NRL_SET_FORC'
    write(6,*)'   --> c_system=',c_system
!
    if (c_system .ne. "NRL_leg") then
      write(6,*)'  : unsupported:c_system=',c_system
      stop
    endif   
!
    call NrlBuildForces(foi, dbij, dpij)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ NRL Extra : additional repulsive part,
!         ( ignored, unless the specific input file exists )
!
    call qm_nrl_set_force_extra
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    do jsv=1,noav
      if (jsv .le. ishow) then
        write(6,*) "FORCE", jsv, foi(jsv,1)
        write(6,*) "FORCE", jsv, foi(jsv,2)
        write(6,*) "FORCE", jsv, foi(jsv,3)
      endif
    enddo
!
    write(6,*)'.... is ended without error :LSES_NRL_SET_FORC'
!
  end subroutine qm_nrl_set_force
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ NRL extra: additional repulsive
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine qm_nrl_set_force_extra
!
    use MNrlAdaptor, only : NrlBuildRepulsiveEnergy,  NrlAddRepulsiveForce
    use elses_mod_foi,        only : foi
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_ene,        only : ecc
    use elses_mod_elem_name,  only : elem_name
!
    implicit none
    real(8) :: ecc_extra
!
    character(len=64) :: file_name
    logical :: file_exist
    integer :: iunit, nss
    character(len=64) :: chara1
!
    nss=1
    if (elem_name(nss) /= 'Pt') then
      if (i_verbose >= 1) then
        write(6,*)'@@ NRL EXTRAt... is ignored'
      endif  
      return
    endif
!
    file_name='input_nrl_extra.txt'
    inquire (file=trim(file_name), exist=file_exist)
    if (file_exist .eqv. .false.) then 
      if (i_verbose >= 1) then
        write(6,*)'@@ NRL EXTRAt... is ignored'
        stop
      endif   
    endif  
!
    iunit=vacant_unit()
    open(iunit, file=file_name,status='old')
    read(iunit,*) chara1
    if (trim(chara1) .ne. 'Pt') then
      if (i_verbose >= 1) then
        write(6,*)'@@ NRL EXTRAt... is ignored'
        stop
      endif   
    endif   
!
    if (i_verbose >= 1) then
      write(6,*)'@@ NRL EXTRA: Repulsive part:',chara1
    endif  
!
    ecc_extra=0.0d0
    call NrlBuildRepulsiveEnergy(ecc_extra)
    if (i_verbose >= 1) then
      write(6,*)'NRL EXTRA:ecc_extra=',ecc_extra
    endif  
    call NrlAddRepulsiveForce(foi)
!
    ecc=ecc+ecc_extra
!
    close(iunit)
!
  end subroutine qm_nrl_set_force_extra
!
end module M_qm_nrl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elses_nrl_ini(rcut_ham)
!
  use M_qm_nrl, only : qm_nrl_ini ! (routine)
  implicit none
  real(8)   :: rcut_ham
!
  call qm_nrl_ini(rcut_ham)
!
end subroutine elses_nrl_ini
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elses_init_alloc_nrl_leg
  use M_qm_nrl, only : qm_nrl_init_alloc ! (routine)
!
  call qm_nrl_init_alloc
!
end subroutine elses_init_alloc_nrl_leg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elses_set_hami_nrl_leg
  use M_qm_nrl, only : qm_nrl_set_ham ! (routine)
!
  call qm_nrl_set_ham
!
end subroutine elses_set_hami_nrl_leg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elses_nrl_set_forc
  use M_qm_nrl, only : qm_nrl_set_force ! (routine)
!
  call qm_nrl_set_force
!
end subroutine elses_nrl_set_forc

