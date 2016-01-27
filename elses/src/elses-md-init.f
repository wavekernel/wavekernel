!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-md-init.f @@@@@
czzz  @@@@@ 2009/08/02 @@@@@
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0104: Prepared. (NT07Ap01) 
c       based on elses-kr-ccvd.f on md-ccvd/work09
c    : Prepared : 'elses_init_atm_file30' 
c            based on INIT30D in zcarbon.f on stable-001v12
c0105: Prepared:'elses_set_mass(awt)'
c            based on 'elses_mass_ccvd' in elses-kr-ccvd.f
c0106: Prepared:'elses_vel_ini' (NT07A,p01)
c0109: 'itemdmx' is set manually in 'elses_init_atm_file30'
c0115: In 'elses_init_atm_file30', 'amq=amq_atom*dble(noa)' is added.
c0131: Bug fix (NT07A,p13):   amq=amq_atom/dble(noa)
c0202: Fix for compatibility (NT07A,p13):   awt=12.01
c0205: PT41/work02 starts (NT07A,p21-)
c    : In 'elses_init_atm_file30', the followings are given (NT07A,p21)
c        i_system, r_cut_book
c0209: In 'elses_init_atm_file30', c_system is introduced (NT07A,p21)
c    : In 'elses_init_atm_file30', small systems are supported
c        ----> following lines added after reading 'noa' from file30
c               if (noao .gt. noav) noao=noav
c               if (noas .gt. noav) noas=noav
c               if (noab .gt. noav) noab=noav
c0322: In 'elses_init_atm_file30', "NRL_Au_only" is supported.
c        (NT07A,p37)
c        ----> file53 is used as input file, instead of file30
c              file30 is used only for indicating c_system
c0324: In 'elses_init_atm_file30' and 'elses_set_mass', 
c       awt is replaced by awt0, so as to avoid the conflict of name
c             (NT07A-118,p.39)
c    : awt(1:nos) is introduced, defined in 'elses_init_atm_file30'
c             (NT07A-118,p.39)
c0326: "NRL_Au_only" is renamed "NRL_leg" (leg means legacy)
c    : In 'elses_init_atm_file30',neig0 is introduced 
c         for eigen-state solver
c0426: : In 'elses_init_atm_file30', 'call elses_xml_test' is added.
c0430: : In 'elses_init_atm_file30', 'call elses_xml_test2' is added.
c0502: : In XML case, 'elses_init_atm_xml' is called, 
c               instead of elses_init_atm_file30
c0511: In 'elses_init_atm_file30', 'call elses_xml_test2' is deleted.
c    : A modi on elses_set_mass
c0517: In 'elses_init_atm_file30', (NT07B-119,p15)
c         Valence electron number per atom is prepared.
c    : Prepared; subroutine 'elses_set_val_elec_tot' (NT07B-119,p15)
c         ----> set the total valence electron number
c0529: In 'elses_md_init', 'call elses_vel_ini' is back. 
c             (NT07B-119,p.22)
c0606: Prepaed 'elses_vel_ini_chk' 
c        ---> checking the initial velocity
c               and adjust to ensure zero total momentum
c0607: A bugfix in 'elses_vel_ini_chk' 
c    : 'elses_test_ran_num' is now commented out with '!dm' mark
c    : Several write commands are commented out with  '!dm' mark
c0722: In 'elses_vel_ini_chk', elses_vel_ini is called
c         when no velocity is set in the input file.
c0731: In 'elses_set_js4jsv', an error checking  is added.
c0823: A modi. in 'elses_set_val_elec_tot' (NT07C-120-p26)
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0104T.Hoshi; Unused routine are deleted;
c         'elses_init_atm_file30', 'elses_init_compat'
c0724T.Hoshi; In elses_md_init, for 'geno' case,
c    'call elses_init_main' is replaced by 'call ini_load_geno'
cccc2009cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0530T.Hoshi(geno); In 'elses_vel_ini',
c       'use M_lib_random_num' is added NT09B-131 
c0802T.Hoshi(geno); In 'elses_md_init',
c       elses_ini_set_velocity is called, instead of elses_vel_ini_chk
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Initial Alloc.for LSES 
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_init_alloc2b
      use elses_mod_sim_cell, only : noa
      use elses_mod_noav,     only : noav,noao,noab,noas,noab0,nncut
      use elses_mod_orb1,     only : nvl
c        NOTE: ALL the above quantities are input parameters.
c
      implicit none
      integer noa_def,noav_def,noao_def,noas_def,nvl_def,nncut_def,imode
      real*8 ddd
c
      write(*,*)'@@ LSES_INIT_ALLOC2:NOA=',noa
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Set parameters for matrix size
c
      write(*,*)'   NOAV =',noav
      write(*,*)'   NOAO =',noao
      write(*,*)'   NOAS =',noas
      write(*,*)'   NOAB =',noab
      write(*,*)'   NOAB0=',noab0
      write(*,*)'   NVL  =',nvl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Check the input parameters
c
      ddd=dble(noav)*dble(noao)*dble(noas)
     +    *dble(noab)*dble(noab0)*dble(nvl)
c
      if (dabs(ddd) .lt. 1.0d-10) then
        write(6,*)'ERROR!:LSES_INIT_ALLOC2B'
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Set working parameters 
c
      noa_def =noa
      noav_def=noav
      noao_def=noao
      noas_def=noas
      nncut_def=nncut
      nvl_def=nvl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Allocations
c
      write(*,*)' Alloc: JSV4JSD, NJSD'
      imode=1
      call elses_alloc_atm(noa_def,noav_def,noao_def,nncut_def,imode)
c        ---->  alloc : js4jsv(noav)
c                     : jsv4js(noa)
c                     : jsv4jsd(noao,noav)
c                     : njsd(noav,0:nncut)
c
      call elses_alloc_mat1(nvl_def,noas_def,noav_def,imode)
c        ---->  alloc : dbij(nvl,nvl,noas,noav)
c                     : dhij(nvl,nvl,noas,noav)
c
c     call elses_alloc_dfij(nvl_def, noas_def, noav_def,imode)
c        ---->  alloc : dfij(nvl,nvl,noas,noav)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Initial setting
c
      call elses_set_js4jsv
c        ---->  set : js4jsv(noav), jsv4js(noa)
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Setting JS4JSV and JSV4JS
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_js4jsv
      use elses_mod_sim_cell, only : noa
      use elses_mod_noav,     only : noav
      use elses_mod_js4jsv,     only : js4jsv,jsv4js
c
      implicit none
!     integer noa_def,noav_def,noao_def,nncut_def,imode ! unused
      integer js,jsv
c
      write(*,*)'@@ LSES_SET_JS4JSV'
c
      if (noa .ne. noav) then
        write(*,*)'ERROR!:NOA,NOAV=',noa,noav
        stop
      endif
c
      if (allocated(js4jsv) .eqv. .false.) then
         write(*,*)'ERROR:JS4JSV is not allocated'
         stop
      endif   
c
      if (allocated(jsv4js) .eqv. .false.) then
         write(*,*)'ERROR:JSV4JS is not allocated'
         stop
      endif   
c
      if (noa .eq. noav) then
        write(*,*)' Set JSV4JS and JS4JSV for equivalent'
        do js=1,noa
         jsv=js
         js4jsv(jsv)=js
         jsv4js(js)=jsv
        enddo 
      endif   
c
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Set the total electron number : val_elec_tot
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_val_elec_tot
      use elses_mod_ctrl,     only : i_verbose
      use elses_mod_val_elec, only : val_elec_atm, val_elec_tot
      use elses_mod_sim_cell, only : noa,nos
      use elses_mod_tx,       only : jsei
      implicit none
      integer js,nss
      real(8) ddsum
!
      if ((nos <= 0) .or. (nos >=100)) then
         write(*,*)'ERROR:NOS='
         stop
      endif   
!
      if (i_verbose >=1) then
        write(*,*)'@@ elses_set_val_elec_tot:nos=',nos
      endif  
!
      if (allocated(val_elec_atm) .eqv. .false.) then
         write(*,*)'ERROR:val_elec_atom is not allocated'
         stop
      endif   
!
      if ((nos <= 0) .or. (nos >=100)) then
         write(*,*)'ERROR:NOS='
         stop
      endif   
!
      do nss=1,nos
        if (val_elec_atm(nss) <= 0.0d-10) then
           write(*,*)'EERROR:val_elec_atm=',nss,val_elec_atm(nss)
           stop
        endif  
        if (i_verbose >=1) then
           write(*,*)'val_elec_atm=',nss,val_elec_atm(nss)
        endif  
      enddo

!
      ddsum=0.0d0
      do js=1,noa
        nss=jsei(js)
        ddsum=ddsum+val_elec_atm(nss)
      enddo
      val_elec_tot=ddsum
!
      if (i_verbose >=1) then
        write(*,*)'val_elec_tot=',val_elec_tot
      endif  
!
      end subroutine elses_set_val_elec_tot
