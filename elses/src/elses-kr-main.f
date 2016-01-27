!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-kr-main.f @@@@@
czzz  @@@@@ 2007/08/12 @@@@@
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0109: Prepared. (NT07A,p29) only for 'elses_kr_main'
c0812: Now 'call elses_force_xu' is commented out; NT07C-120p11
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Main routine for LSES-KR
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_kr_main
      use elses_mod_md_dat,   only : itemd
      use elses_param_ctl_kr, only : iKnstr !(CHANGED!)
      use elses_mod_ene,      only : etb
      use M_config,           only : config !(unchanged)
c
      implicit none
      integer imode, iselect
      integer nrec_def, nvl_def, noav_def, noav_kr_def
      integer njjt_mx2_def, nreclc_mx2_def
      real*8  etbtmp
c
      write(*,*)'@@ LSES_KR_MAIN:ITEMD',itemd
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call elses_set_param_ctl_kr
c
c     nrec_def=nrecmx
c     noav_def=noav
c     nvl_def=nvl
c
      imode=1
      call elses_alloc_kry_glo(imode)
c
c
      imode=1
      call elses_set_rcut_kry(imode)
c
      imode=1
      call elses_alloc_jsv4jsk_str(imode)
      call elses_set_jsv4jsk_str(imode)
c
      iKnstr=0
c
      if (config%calc%solver%mode_for_large_memory .ge. 1) then
        write(*,*)'INFO'
c       write(*,*)'INFO: iKnstr is set to be one:  mode_for_large_memory=', 
c    +              config%calc%solver%mode_for_large_memory
        iKnstr=1
      endif   
c
      imode=1
c
      if ((iKnstr .eq. 1) .and. (itemd .eq. 1)) then 
         imode=1
         call elses_alloc_kry_rKnstr(imode)
c            ---> rKnstr is allocated, when itemd=0
      endif
c
      imode=1
      iselect=0
      call elses_set_kry(imode,iselect)
c
      call elses_check_rn
c
      call elses_kr_set_chempot
c
      call elses_cal_rho
c
      imode=2
      iselect=2
      call elses_set_kry(imode,iselect)
c
      imode=2
      call elses_alloc_kry_glo(imode)
      call elses_alloc_jsv4jsk_str(imode)
c
      call elses_sym_dbij
c
      etbtmp=0.0d0
      call elses_enetb(etbtmp)
      write(6,*)'ETB=,',etbtmp
      etb=etbtmp
c
!     call elses_force_xu
c
      end subroutine elses_kr_main
