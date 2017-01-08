!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-mod-02.f @@@@@
czzz  @@@@@ 2008/01/08 @@@@@
cccc2006cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0722: Prepared. (NT116p31): Modules for LSES engine
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0324: 'target' is added for the compatibility of NRL-part code.
c          (NT07A-118,p39)
c        ---> js4jsv, jsv4jsd, njsd 
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0108: elses_mod_elec_cond is added (NT07E-122,p39)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      module elses_mod_js4jsv
        integer, target, allocatable :: js4jsv(:)
        integer,         allocatable :: jsv4js(:)
      end module
c
      module elses_mod_jsv4jsd
        integer, target, allocatable :: jsv4jsd(:,:)
        integer, target, allocatable :: njsd(:,:)
      end module
c
      module elses_mod_noav
        integer  noav
        integer  noao
        integer  noab
        integer  noas
        integer  noab0
        integer  nncut
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Module for conditions for electronic system
c
      module elses_mod_elec_cond
        real*8 temp_for_electron
c         ---> temperature (level-broadening) parameter 
c                  for electron system 
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_ene
        real*8 etb, ecc, ecsc
        real*8, allocatable :: enea(:,:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Module for multi-solver scheme 
c
      module elses_mod_multi
        integer  ict4h, ict4l
        integer, allocatable :: intf(:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of matrices related to atom freedom 
c
c          jsv4js(noa)
c          js4sjsv(noav)
c          jsv4jsd(noao,noav)
c             njsd(noav,0:nncut)
c
c          enea(noa,2)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_atm
     +       (noa_def,noav_def,noao_def,nncut_def,imode)
c
      use elses_mod_js4jsv,  only : js4jsv,jsv4js
      use elses_mod_jsv4jsd, only : jsv4jsd,njsd
      use elses_mod_ene,     only : etb, ecc, enea
      implicit none
      integer imode, ierr
      integer noa_def, noav_def, noao_def, nncut_def
c
      if (imode .ne. 1) goto 1999
c
      write(6,*)'@@ LSES_ALLOC_ATM:allocation of matrices'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate (js4jsv(noav_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(JS4JSV):ierr=',ierr
        write(6,*)' noav_def=',noav_def
        stop
      endif
      js4jsv(:)=0
c
      allocate (jsv4js(noa_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(JSV4JS):ierr=',ierr
        stop
      endif
      jsv4js(:)=0
c
      allocate (jsv4jsd(noao_def,noav_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(JSV4JSD):ierr=',ierr
        stop
      endif
      jsv4jsd(:,:)=0
c
      allocate (njsd(noav_def,0:nncut_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(NJSD):ierr=',ierr
        stop
      endif
      njsd(:,:)=0
c
      allocate (enea(noa_def,2),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(JS4JSV):ierr=',ierr
        stop
      endif
      enea(:,:)=0.0d0
c
      write(6,*)'..is succesfully done'
c
 1999 continue
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of INTF
c
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_intf(noa_def, imode)
c
      use elses_mod_multi, only : intf
      implicit none
      integer noa_def, imode, ierr
c
      write(6,*)'@@ LSES_ALLOC_INTF'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate (intf(noa_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(JS4JSV):ierr=',ierr
        write(6,*)' noa_def=',noa_def
        stop
      endif
      intf(:)=0
c
      write(6,*)'..is succesfully done'
c
      end
