!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-qm-main.f @@@@@
czzz  @@@@@ 2008/07/23 @@@@@
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0109: Prepared. (NT07A,p01) 
c       based on elses-kr-ccvd.f 
c    : 'elses_kr_core' is prepared. Formerly named as 'elses_kr_main'
c0206: Supported : Silicon system (i_system=2) (NT07A,p21)
c0209: c_system is introduced instead ob i_system=2 (NT07A,p21)
c    : In elses_engin, 
c         'call elses_write_molden_dat' is prepared.
c0309: Renamed as 'elses-es-main.f' (NT07Ap29)
c    : 'call elses_kr_core' is renamed as 'elses_kr_main' (NT07Ap29)
c    : Subroutine 'elses_kr_main' is moved into elses-kr-main.f
c0326: In 'elses_set_hami', c_system="NRL_leg" supported.
c         (NT07A-118,p.39)
c0513: File name is renamed as elses-qm-main.f
c    : In 'elses_engine', switch for 'elses_engine_nrl_leg' is added.
c       (NT07B-119,p.6)
c0606: In elses_engin, 
c         'call elses_write_molden_dat' is now commented out.
c0725: GaAs case is supported in 'elses_es_rest' and 'elses_set_hami'
c0731: A modi. on elses_engine.
c0812: 'elses_qm_force' is parepared; NT07C-120p11
c1225: 'elses_plot_hami' is called; NT07E-122p17
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0103T.Hoshi; 'call elses_set_hami_01(0)' NT07E-122p29
c0105T.Hoshi; 'i_eig_main_init' is introduced for 'elses-engine'
c                 (NT07E-122p31, v000a-wrk05)
c0213T.Hoshi; In 'elses-qm-force', 
c        'call elses_force_xu' is replaced by 'call elses_qm_force_01'
c0213T.Hoshi; In 'elses-qm-force', 'c_system == Si_kwon' is supported
c                 (NT07E-122p31, v0.0.12a-wrk01)
c0401T.Hoshi; Delete unused entry of 'iperiodic' in 'use' command
c0402T.Hoshi; In elses_engine, (NT08A-123,p31)
c              elses_write_molden_dat is called, in the verbose mode
c0723T.Hoshi; Simple bufix in elses_engine, 
c               call elses_kr_main(imode)  ----> call elses_kr_main
c               call elses_qm_force(imode) ----> call elses_qm_force
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Main routine for 'rest' part of energy function
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_es_rest
      use elses_mod_sel_sys, only : c_system
      implicit none
      integer ierr
c
      write(6,*)'@@ LSES_ES_REST:c_system=',c_system
c
      ierr=1
c
      if (c_system .eq. "C_Xu") then 
        call elses_set_rest_01
        ierr=0
      endif  
c
      if (c_system .eq. "Si_Kwon") then 
        call elses_set_rest_01
        ierr=0
      endif  
c
      if (c_system .eq. "GaAs.Molteni.1993") then 
        call elses_set_rest_gaas_Molteni
        ierr=0
      endif  
c
      if (ierr .eq. 1) then
        write(6,*)'ERROR!(LSES_ES_REST):ierr=',ierr
        stop
      endif
c   
      end subroutine elses_es_rest
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Main routine for buliling hamiltonian
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_hami
      use elses_mod_md_dat,     only : itemd !(unchaged)
      use elses_mod_sel_sys,    only : c_system !(unchaged) 
      use elses_mod_ctrl,       only : i_verbose !(unchanged)
      use elses_mod_noav,     only : noav,noas !(unchanged)
      use elses_mod_orb1,     only : nvl !(unchanged)
      implicit none
      integer ierr
      integer nvl_def, noas_def, noav_def, imode
c
      if (i_verbose .ge. 1) then
        write(6,*)'@@ LSES_SET_HAMI:c_system=',c_system
      endif  
c
      noav_def=noav
      noas_def=noas
      nvl_def=nvl
      imode=1
c
      if (itemd .eq. 1) then
        call elses_alloc_dfij(nvl_def, noas_def, noav_def,imode)
c          ---->  alloc : dfij(nvl,nvl,noas,noav)
      endif  

      ierr=1
c
      if (c_system .eq. "C_Xu") then 
        call elses_set_hami_01(0)
c       call elses_set_hami_01_new
        ierr=0
      endif  
c
      if (c_system .eq. "Si_Kwon") then 
        call elses_set_hami_01(0)
        ierr=0
      endif  
c
      if (c_system .eq. "NRL_leg") then 
        call elses_set_hami_nrl_leg
        ierr=0
      endif  
c
      if (c_system .eq. "GaAs.Molteni.1993") then 
        call elses_set_hami_gaas_Molteni
        ierr=0
      endif  
c
      if (ierr .eq. 1) then
        write(6,*)'ERROR!(LSES_ES_REST):ierr=',ierr
        stop
      endif
c   
      end subroutine elses_set_hami
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Main routine for force calculation
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_qm_force
      use elses_mod_sel_sys, only : c_system
      use elses_mod_ctrl,       only : i_verbose
      implicit none
      integer ierr
c
      if (i_verbose >= 1) then
        write(6,*)'@@ LSES_QM_FORCE:c_system=',c_system
      endif  
c
      ierr=1
c
      if (c_system .eq. "C_Xu") then 
        call elses_qm_force_01
        ierr=0
      endif  
c
      if (c_system .eq. "Si_Kwon") then 
        call elses_qm_force_01
        ierr=0
      endif  
c
      if (c_system .eq. "GaAs.Molteni.1993") then 
        call elses_force_molteni
        ierr=0
      endif  
c
      if (ierr .eq. 1) then
        write(6,*)'ERROR!(LSES_ES_REST):ierr=',ierr
        stop
      endif
c   
      end subroutine elses_qm_force
c      
