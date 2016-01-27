!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-kr-sub.f @@@@@
czzz  @@@@@ 2009/06/30 @@@@@
cccc2006cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0727: Prepared: (NT116p33-)
c0801: Continued.
c0830: A modi.:Now '...all atoms are into list_atom' is not plotted.
c    : A modi.
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0109: Renamed elses-kr-sub.f, instead of elses-kr-eng.f
c    : 'elses_bkmatv_kr2' is copied elses-kr-ccvd.f
c0114: (continued) NT07A,p3
c0206: Prepared 'elses_es_bkmat' based on 'elses_bkmatv_kr2' 
c       --> rcut_bk is determined from i_system (NT07A,p21)
c0209: c_system is used instead of i_system (NT07A,p21)
c    : 'elses_set_rcut_kry' suppors smaller systems 
c        -->  'if ((noak .gt. noak_c) .or. (noak .eq. noav))'
c    : 'elses_set_param_ctl_kry' suppors smaller systems 
c          ---> if (noak_min_def .gt. noa) noak_min_def=noa
c0607: Several write commands are commented out with '!dm' mark.
c0805: In 'elses_es_bkmat', i_verbose is introduced. (NT07B-119p69)
c0809: A mod.
c0823: In 'elses_set_param_ctl_kr', (NT07C-120,p26)
c          rNelec is now set from val_elec_tot
c0910: In 'elses_set_param_ctl_kr', (NT07D-121,p07)
c        Krylov dimension is doubly defined; nreclc_def, nrec_kr_def
c0911: In 'elses_set_param_ctl_kr', (NT07D-121,p07)
c        noak_min_def, nreclc_def are now assumed to be already set.
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0108: Trivial rrror checking is added in 'elses_enetb'
c        (NT07E-122p39, v000a-wrk05)
c0122T.Hoshi; In 'elses_set_rcut_kry', r_cut_book is refered
c               for booking
c0128T.Hoshi; In 'elses_set_rcut_kry', noalst=noa (NT07E-112,p55)
c0129T.Hoshi; Written on v.0.0.10 (NT07E-112,p57)
c0401T.Hoshi; i_pbc_x, i_pbc_y, i_pbc_z are introduced. 
!             (NT08A123-p.25) in subroutine 
!              'elses_set_rcut_kry' 'elses_set_jsv4jsk_str'
!              'elses_es_bkmat' 'elses_bkmatv_kr2'
!0630T.Hoshi; In 'elses_set_rcut_kry', (NT08C-125,p35)
!              The following line is commented out
!   ---> if (c_system .eq. "C_Xu") akwon=4.16d0/angst
!0803T.Hoshi; in  'elses_kr_set_chempot' (NT08D-126p41)
!               temp_for_electron is now refered.
!           ; in  'elses_kr_set_param_ctrl_kr' (NT08D-126p41)
!               ' kr_green_ldos_ini' is called. 
!0812T.Hoshi; 'implicit none' is introduced in
!               'elses_kr_set_chempot' 'elses_sym_dbij'
!              Several variables are declared explicitly.
cccc2009cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0629T.Hoshi; Unused routines are deleted :
!                 elses_bkmatv_kr2 
!     ; In elses_es_bkmat, 'get_max_size_of_booking_list' is called
!           (log-2009-06-29b)
!0630T.Hoshi; In 'elses_cal_rho', (log-2009-06-29c-NT09C-132)
!               an error checking routine is modified
!           ; A modi on elses_cal_rho
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Set the parameters in PARAM_CTL_KR
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_param_ctl_kr
      use M_qm_domain, only : nval
      use elses_mod_phys_const, only : ev4au
c     use elses_mod_sim_cell,   only : noa
      use elses_mod_ctrl,       only : i_verbose
      use elses_mod_noav    ,   only : noav
      use elses_mod_val_elec,   only : val_elec_tot
      use elses_param_ctl_kr,   only : noav_kr_def, noak_min_def, 
     +                 nreclc_def,nrec_kr_def,nval_kr_def,nkene_def,
     +                 ddemin_def,ddemax_def,rNelec_kr_def
      implicit none
c
      if ( i_verbose >= 1 ) then
        write(6,*)'@@ SET_PARAM_CTL_KR'
      endif
c
!     noak_min_def=200
!
      if (noak_min_def .gt. noav) then 
         noak_min_def=noav
         write(6,*)'Message in SET_PARAM_CTL_KR'
         write(6,*)'....forced setting: noak_min_def=',noak_min_def
      endif   
!
      if (noak_min_def <= 0) then
        write(6,*)'ERROR:noak_min_def=',noak_min_def
        write(6,*)'   .. at SET_PARAM_CTL_KR'
        stop
      endif   
!
c
!     nreclc_def=30
!
!     if ((nreclc_def <= 0) .or. (dble(nreclc_def) >= dble(noav)*1000.0d0)) then
!       write(6,*)'ERROR:nreclc_def=',nreclc_def
!       write(6,*)'  noav, noav*1000=', noav, noav*1000
!       write(6,*)'   .. at SET_PARAM_CTL_KR'
!       stop
!     endif   
!
      noav_kr_def=noav
!     rNelec_kr_def=4.0d0*dble(noav)
      rNelec_kr_def=val_elec_tot
!
      if (rNelec_kr_def <= 1.0d-10) then
        write(6,*)'ERROR:# elec. : rNelec_kr_def=',rNelec_kr_def
        stop
      endif   
!
c     nrec_kr_def=180
c     nrec_kr_def=100
      nrec_kr_def=nreclc_def
c
!     nval_kr_def=10
      nval_kr_def=maxval(nval)
c
      nkene_def=26
      ddemin_def=-15.0d0/ev4au
      ddemax_def=10.0d0/ev4au
c
      if ( i_verbose >= 1 ) then
        write(6,*)'# atoms : noav_kr_def=',noav_kr_def
        write(6,*)'# elec. : rNelec_kr_def=',rNelec_kr_def
        write(6,*)' ( from val_elec_tot )'
c
        write(6,*)'# atoms in projection radius'
        write(6,*)'  : noak_min_def=',noak_min_def
c
        write(6,*)'KR subspace dimension'
        write(6,*)'  : nreclc_def  =',nreclc_def
c
        write(6,*)'KR subspace dimension for allocation'
        write(6,*)'  : nrec_kr_def  =',nrec_kr_def
c
        write(6,*)'# energy points for  G(e)'
        write(6,*)'  : nkene_def   =',nkene_def
c
        write(6,*)'Bounds of energy points for G(e)'
        write(6,*)'  : ddemin [eV]  =',ddemin_def*ev4au
        write(6,*)'  : ddemax [eV]  =',ddemax_def*ev4au
c
      endif
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Check the parameters in PARAM_CTL_KR
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_chk_param_ctl_kr
      use elses_mod_sim_cell, only : noa
      use elses_param_ctl_kr, only : noav_kr_def, noak_min_def, 
     +                 nreclc_def,nrec_kr_def,nval_kr_def,nkene_def,
     +                 ddemin_def,ddemax_def,rNelec_kr_def
c
      implicit none
      integer ierr
c
      ierr=0
      write(6,*)'@@ CHK_PARAM_CTL_KR'
c
      write(6,*)'# total atoms'
      write(6,*)'  :  noav_kr_def=',noav_kr_def
c
      write(6,*)'# atoms in projection radius'
      write(6,*)'  : noak_min_def=',noak_min_def
c
      write(6,*)'KR subspace dimension'
      write(6,*)'  : nreclc_def  =',nreclc_def
c
      write(6,*)'KR subspace dimension for allocation'
      write(6,*)'  : nrec_kr_def  =',nrec_kr_def
c
      write(6,*)'# energy points for  G(e)'
      write(6,*)'  : nkene_def   =',nkene_def
c
      write(6,*)'Bounds of energy points for G(e)'
      write(6,*)'  : ddemin [eV]  =',ddemin_def
      write(6,*)'  : ddemax [eV]  =',ddemax_def
c
      if (noav_kr_def .le.   0) ierr=1
      if (noav_kr_def .gt. noa) ierr=1
      if (noak_min_def .le.   0) ierr=1
      if (noak_min_def .gt. noa) ierr=1
      if (nrec_kr_def .le. 0) ierr=1
      if (nrec_kr_def .gt. 1000) ierr=1
      if (nreclc_def .le. 0) ierr=1
      if (nreclc_def .gt. nrec_kr_def) ierr=1
      if (nkene_def .le.    0) ierr=1
      if (nkene_def .gt. 1000) ierr=1
      if (ddemin_def .gt. ddemax_def) ierr=1
      if (dabs(ddemin_def-ddemax_def) .lt. 1.0d-10) ierr=1
c
      if (ierr .ne. 0) then
        write(6,*)'... checking fails:ierr=',ierr
        stop
      else   
        write(6,*)'... checking succeeds:ierr=',ierr
      endif   
c
 9999 continue
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Find rcut in the Krylov subspace -3-
c        ---> Search within booked atom
c
c      OUTPUTs: 
c        noak_min(1:noav_kr) 
c              : Required # atoms inside projection
c                ----> set before the main loop
c        rcut_kry(1:noav_kr)
c              : Redius for projection
c                ----> set in the main loop
c        noak_str(1:noav_kr)
c              : Actual # atoms inside projection
c                 ----> set in the main loop
c
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_rcut_kry(imode)
      use elses_mod_phys_const, only : angst
      use elses_mod_ctrl,     only : i_verbose
      use elses_mod_sel_sys,  only : c_system, r_cut_book
!     use elses_mod_sim_cell, only : noa, iperiodic, ax, ay, az
      use elses_mod_sim_cell, only : noa, ax, ay, az,
     +                                 i_pbc_x, i_pbc_y,i_pbc_z 
      use elses_mod_noav,     only : noav
      use elses_mod_tx,       only : tx, ty, tz
      use elses_mod_js4jsv,   only : js4jsv
      use elses_param_ctl_kr, only : noav_kr_def,noak_min_def 
      use elses_arr_kry_glo,  only : rcut_kry,noak_min,noak_str
      use M_qm_proj_suggest_radius, only : suggest_proj_radius_wrap !(routine)
c
c     include 'zconst.f'
c     include 'zconst-VR.f'
c
      implicit none
      integer imode, ierr, ipe, npe
      integer omp_get_thread_num
      integer omp_get_num_threads
      integer,allocatable :: list_atom(:)
c
      integer noalst,noav_kr
      integer jsv4seed,js4seed, ishow2, noak_c,iilst,niilst
      integer jsv,js, itry,ntrymax,jsk,noak
c
      real*8  akwon,acutkry,acutkry0,acutkry_tmp
      real*8  dd,dvecx,dvecy,dvecz,dxc,dyc,dzc,ddd
c
      integer noakmax, noakmin
      real*8  ddsum
c
      logical, parameter :: compat_mode = .false. 
c                 swtich for the exact compatibitily to the legacy code 
c
      noav_kr=noav_kr_def
      if (i_verbose >= 1) then
        write(6,*)'@@ elses_set_rcut_kry:noav_kr,imode=',noav_kr,imode
        write(6,*)'   ... NOTE:non-order-N routine !!!'
      endif  
c
c     akwon=4.16d0/angst
c
      call suggest_proj_radius_wrap(akwon)
      akwon=akwon/2.0d0
c
      if (compat_mode) akwon=r_cut_book
c
!     if (c_system .eq. "C_Xu") akwon=4.16d0/angst
c       NOTE:Non-essential; Just for the exact 
c             compatibility of the result of the old carbon calc.
c
      noalst=noa
c     noalst=min(4096,noa)
c     noalst=40960
c
      if (i_verbose >= 1) then
        write(6,*)'  noalst=',noalst
        write(6,*)'  cutoff for pre-booking(au)=',akwon
        write(6,*)'  cutoff for pre-booking(A )=',akwon*angst
      endif  
c
      if (akwon <= 1.0d-10) then
        write(6,*)'ERROR!(SET_RCUT_KRY3):r_cut_book=',r_cut_book
        stop
      endif
c
      if ((noav_kr .le. 0) .or. (noav_kr .gt. 10000000)) then
        write(6,*)'ERROR!(SET_RCUT_KRY3):NOAV_KR=',noav_kr
        stop
      endif
c   
      if (imode .eq. 0) then
        acutkry=10.0d0*ax
        write(6,*)' acutkry=',acutkry
        ntrymax=1
        rcut_kry(:)=acutkry
        noak_min(:)=0
      endif   
c
      if (imode .eq. 1) then
        ntrymax=1000
        rcut_kry(:)=akwon
        noak_min(:)=noak_min_def
      endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Main (parallel) loop
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe,ierr,ishow2)
!$omp& private(jsv4seed,js4seed,noak_c,iilst,niilst)
!$omp& private(jsv,js,itry,jsk,noak)
!$omp& private(list_atom)
!$omp& private(acutkry0,acutkry_tmp,ddd)
!$omp& private(dvecx,dvecy,dvecz,dd,dxc,dyc,dzc)
      ipe=0
      npe=0
c     ipe=omp_get_thread_num()+1
c     npe=omp_get_num_threads()
c     write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
      do 1000 jsv4seed=1,noav_kr
       ishow2=0
c      ishow2=1
       if (jsv4seed .le. 20) ishow2=1
c
       allocate (list_atom(noalst),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(set_rcut_kry_2)'
         write(6,*)'listr_atom:ierr=',ierr
         stop
       endif
c      write(6,*)'list_ atom is allocated:ipe,npe=',ipe,npe
c
       acutkry0=rcut_kry(jsv4seed)  
       noak_c=noak_min(jsv4seed)
       acutkry_tmp=akwon*10.0d0
c         ---> This value is non-essential
c           (If this value is insufficient, the pre-booking 
c              is done for all atoms)
c
       iilst=0
       do jsv=1,noav_kr
         js=js4jsv(jsv)
         js4seed=js4jsv(jsv4seed)
         if (jsv .ne. jsv4seed) then
          dvecx=tx(js)-tx(js4seed)
          dvecy=ty(js)-ty(js4seed)
          dvecz=tz(js)-tz(js4seed)

          if (i_pbc_x == 1)dvecx = modulo(dvecx+0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_y == 1)dvecy = modulo(dvecy+0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_z == 1)dvecz = modulo(dvecz+0.5d0, 1.0d0) - 0.5d0

!         if (iperiodic .eq. 1) then
!           dd=1.0d0
!           dvecx=dvecx+0.5d0*dd
!           dvecx=dmod(dvecx+dd,dd)-0.5d0*dd
!           dvecy=dvecy+0.5d0*dd
!           dvecy=dmod(dvecy+dd,dd)-0.5d0*dd
!           dvecz=dvecz+0.5d0*dd
!           dvecz=dmod(dvecz+dd,dd)-0.5d0*dd
!         endif

          dxc=dvecx*ax
          dyc=dvecy*ay
          dzc=dvecz*az
          ddd=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
c               ----> ddd : distance in a.u.
          if (ddd .lt. acutkry_tmp) then
             iilst=iilst+1
             if (iilst .gt. noalst) then
                write(6,*)'ERROR!(set_rcut_kry_2)'
                write(6,*)'iilst,noalst=',iilst,noalst
                write(6,*)'Please enlarge the value of noalst'
                stop
             endif   
             list_atom(iilst)=jsv
          endif   
         endif
       enddo
       niilst=iilst
       if (niilst .lt. noak_c) then
c        write(6,*)'Warning!(set_rcut_kry_3)'
c        write(6,*)'jsv4seed,niilst,noak_c=',jsv4seed,niilst,noak_c
c        write(6,*)'acutkry_tmp=',acutkry_tmp
c        write(6,*)'tx(seed)=',tx(js4seed)
c        write(6,*)'ty(seed)=',ty(js4seed)
c        write(6,*)'tz(seed)=',tz(js4seed)
c        if (noav_kr .gt. noalst) then
c          write(6,*)'ERROR!(set_rcut_kry_2)'
c          write(6,*)'     NOAV_kr,NOALST=',noav_kr,noalst
c          stop
c        endif   
         do iilst=1,noav_kr
            jsv=iilst
            list_atom(iilst)=jsv
         enddo   
         niilst=noav_kr
c        write(6,*)'...all atoms are into list_atom:jsv=',jsv4seed
c        write(6,*)'  niilst=',niilst
       endif   
c        ----> pre-booking is done.
c
       do 1100 itry=1,ntrymax  
         if (imode .eq. 0) acutkry_tmp=acutkry0
         if (imode .ne. 0) then 
            acutkry_tmp=akwon*(2.0d0+dble(itry-1)/10.0d0)
         endif   
         jsk=1
         do iilst=1,niilst
           jsv=list_atom(iilst) 
           if ((jsv .le. 0) .or. (jsv .gt. noav_kr)) then
             write(6,*)'ERROR!(set_rcut_kry_2)'
             write(6,*)'jsv4seed,iilst,jsv=',jsv4seed,iilst,jsv
             stop
           endif   
           js=js4jsv(jsv)
           if (jsv .ne. jsv4seed) then
            dvecx=tx(js)-tx(js4seed)
            dvecy=ty(js)-ty(js4seed)
            dvecz=tz(js)-tz(js4seed)

            if (i_pbc_x == 1)dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
            if (i_pbc_y == 1)dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
            if (i_pbc_z == 1)dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0

!            if (iperiodic .eq. 1) then
!              dd=1.0d0
!              dvecx=dvecx+0.5d0*dd
!              dvecx=dmod(dvecx+dd,dd)-0.5d0*dd
!              dvecy=dvecy+0.5d0*dd
!              dvecy=dmod(dvecy+dd,dd)-0.5d0*dd
!              dvecz=dvecz+0.5d0*dd
!              dvecz=dmod(dvecz+dd,dd)-0.5d0*dd
!            endif

             dxc=dvecx*ax
             dyc=dvecy*ay
             dzc=dvecz*az
             ddd=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
c                ----> ddd : distance in a.u.
             if (ddd .lt. acutkry_tmp) then
               jsk=jsk+1
               if (jsk .gt. noav_kr) then
                 write(6,*)'ERROR!(SET_RCUT_KRY):JSK=',jsk
                 stop
               endif
             endif   
          endif  
        enddo   
        noak=jsk
        if (ishow2 .eq. 1) then
          write(6,*)'jsv,itry,noak,cut=',jsv,itry,noak,acutkry_tmp
        endif   
        if ((noak .gt. noak_c) .or. (noak .eq. noav)) then
          noak_str(jsv4seed)=noak
          rcut_kry(jsv4seed)=acutkry_tmp
          if (ishow2 .eq. 1) then
!dm          write(6,9003)jsv4seed,itry,noak,acutkry_tmp
          endif   
c         if (niilst .eq. noav_kr) then
c           write(6,9003)jsv4seed,itry,noak,acutkry_tmp
c         endif   
          goto 1101
        endif   
        if (itry .eq. ntrymax) then
          if (imode .ne. 0) then 
            write(6,*)'ERROR!(SET_RCUT_KRY):IMODE=',IMODE
            write(6,*)' NOAK,NOAK_MIN=',noak,noak_c
            write(6,*)'    ITRY=',itry
            write(6,*)' NTRYMAX=',ntrymax
            stop
          endif   
          noak_str(jsv4seed)=noak
        endif   
 1100  continue  
 1101  continue  
       deallocate (list_atom)
 1000 continue  
!$omp end do
!$omp end parallel 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Plot the data (optional)
c
      noakmax=0
      noakmin=noav_kr
      ddsum=0.0d0
      do jsv=1,noav_kr
        ishow2=0
        ishow2=1
        if (jsv .le. 3) ishow2=1
        noak=noak_str(jsv)
        ddd=rcut_kry(jsv)
c       if (ishow2 .eq. 1) then
c         write(6,9002)jsv,noak,ddd
c       endif  
        if (noak .gt. noakmax) noakmax=noak
        if (noak .lt. noakmin) noakmin=noak
        ddsum=ddsum+dble(noak)
      enddo
      ddsum=ddsum/dble(noav_kr)
      write(6,9001)noav_kr,noakmax,noakmin,ddsum
 9001 format('noav_kr,noak_max,min,ave=',3I10,F13.3)
 9002 format('jsv,noak,rcut=',2I10,F10.5)
 9003 format('jsv,itry,noak,rcut=',3I10,F15.5)
c
 9999 continue
c

c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Find rcut in the Krylov subspace -4-
c        ---> Search within booked atom (Non-order-N calc.)
c
c      OUTPUTs: 
c        jsv4jsk_str(1:noak,noav)
c
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_jsv4jsk_str(imode)
      use elses_mod_phys_const, only : angst
!     use elses_mod_sim_cell, only : noa, iperiodic, ax, ay, az
      use elses_mod_sim_cell, only : noa, ax, ay, az,
     +                                 i_pbc_x, i_pbc_y,i_pbc_z 
      use elses_mod_tx,       only : tx, ty, tz
      use elses_mod_js4jsv,   only : js4jsv
      use elses_param_ctl_kr, only : noav_kr_def
      use elses_arr_kry_glo,  only : rcut_kry,noak_str,jsv4jsk_str
c
      implicit none
      integer imode
      integer noav_kr, noak_tmp, jsv4seed, js4seed, jsk, jsv, js
      integer noak
      real*8  acutkry_tmp,dvecx,dvecy,dvecz,dxc,dyc,dzc,dd,ddd
c
      noav_kr=noav_kr_def
      write(6,*)'@@ LSES_SET_JSV4JSK'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Check that jsv4jsk is already allocated ?
c
      if (allocated(jsv4jsk_str) .eqv. .false.) then
        write(6,*)'ERROR!:Not yet allocated: JSV4JSK'
        stop
      else  
        write(6,*)' JSV4JSK is already allocated; OK'
      endif 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do jsv4seed=1,noav_kr
        js4seed=js4jsv(jsv4seed)
        acutkry_tmp=rcut_kry(jsv4seed)
        noak_tmp=noak_str(jsv4seed)
c
        jsk=1
        jsv4jsk_str(jsk,jsv4seed)=jsv4seed
c          ----> (jsk=1) corresponds to (jsv=jsv4seed)
c
c
        do jsv=1,noav_kr
           js=js4jsv(jsv)
           if (js .ne. js4seed) then
            dvecx=tx(js)-tx(js4seed)
            dvecy=ty(js)-ty(js4seed)
            dvecz=tz(js)-tz(js4seed)

            if (i_pbc_x == 1)dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
            if (i_pbc_y == 1)dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
            if (i_pbc_z == 1)dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0

!           if (iperiodic .eq. 1) then
!             dd=1.0d0
!             dvecx=dvecx+0.5d0*dd
!             dvecx=dmod(dvecx+dd,dd)-0.5d0*dd
!             dvecy=dvecy+0.5d0*dd
!             dvecy=dmod(dvecy+dd,dd)-0.5d0*dd
!             dvecz=dvecz+0.5d0*dd
!             dvecz=dmod(dvecz+dd,dd)-0.5d0*dd
!           endif

            dxc=dvecx*ax
            dyc=dvecy*ay
            dzc=dvecz*az
            ddd=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
c               ----> ddd : distance in a.u.
            if (ddd .lt. acutkry_tmp) then
              jsk=jsk+1
              if (jsk .gt. noav_kr) then
                write(6,*)'ERROR!(SETKRY):JSK=',jsk
                stop
              endif
              jsv4jsk_str(jsk,jsv4seed)=jsv
            endif   
          endif  
        enddo   
        noak=jsk
        if (noak .ne. noak_tmp) then
           write(6,*)'ERROR!(SET_JSV4JSK_STR):NOAK=',noak,noak_tmp
           stop
        endif   
c
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 9999 continue
      write(6,*)'... ended succesfully: LSES_SET_JSV4JSK'
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Check the Residual norm 
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_check_rn
      use elses_mod_tx,       only : jsei
      use elses_mod_md_dat,   only : itemd
      use elses_mod_js4jsv,   only : js4jsv 
      use elses_mod_orb1,     only : nval
      use elses_param_ctl_kr, only : noav_kr_def
      use elses_arr_kry_glo,  only : rRNave, rRNdia
c
c     use arr_kry_glo
c     include 'zconst.f'
c     include 'zconst-VR.f'
c
      implicit none
      integer ishow2, ipe, npe, ierr
      integer iibase, iii1, iii2, iii3, iii4
      integer noav_kr, jsv, js, nss, nval2, ja
      real*8 tb, ptime
      real*8 ddsum, ddd, ddd1, ddave, ddmax1, ddmax2
c
      write(6,*)'@@ CHECK_RN_MP'
c
      call tclock(ptime)
      tb=ptime
c
      noav_kr=noav_kr_def
      ddsum=0.0d0
      iibase=0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(ishow2,ierr)
!$omp& private(jsv,js,nss,nval2,ja)
!$omp& private(ddd,ddd1)
!$omp& reduction(+ : ddsum)
!$omp& reduction(+ : iibase)
!$omp do schedule(static)
      do jsv=1,noav_kr
         ishow2=0
         ishow2=1
c        if (jsv .le. 4) ishow2=1
         js=js4jsv(jsv)
         nss=jsei(js)
         nval2=nval(nss)
         iibase=iibase+nval2
         if (nval2 .eq. 0) then
           write(6,*)'ERROR!(check_RN):nval2=',nval2
           stop
         endif   
         do ja=1,nval2
           ddd=rRNave(ja,jsv)
           ddd1=rRNdia(ja,jsv)
c          write(6,8002)jsv,ja,ddd,ddd1
           ddsum=ddsum+ddd
         enddo   
      enddo   
!$omp end do
!$omp end parallel 
c
      if (iibase .eq. 0) then
        write(6,*)'ERROR!(check_RN):iibase=',iibase
        stop
      endif   
      ddave=ddsum/dble(iibase)
      write(6,*)'RN_AVE=',itemd,ddave
c
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-0310 ',ptime-tb
      tb=ptime
c
      ddmax1=0.0d0
      ddmax2=0.0d0
      iii1=0
      iii2=0
      iii3=0
      iii4=0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(ishow2,ierr)
!$omp& private(jsv,js,nss,nval2,ja)
!$omp& private(ddd,ddd1)
!$omp& reduction(+ : ddsum)
!$omp& reduction(max : ddmax1,ddmax2)
!$omp& reduction(+ : iii1,iii2,iii3,iii4)
!$omp do schedule(static)
      do jsv=1,noav_kr
         js=js4jsv(jsv)
         nss=jsei(js)
         nval2=nval(nss)
         do ja=1,nval2
           ddd=rRNave(ja,jsv)
           ddd1=rRNdia(ja,jsv)
           ddmax1=max(ddmax1,ddd)
           ddmax2=max(ddmax2,ddd1)
           if (ddd  .gt. ddave*1.1d0)  iii1=iii1+1
           if (ddd  .gt. ddave*1.5d0)  iii2=iii2+1
           if (ddd  .gt. ddave*2.0d0)  iii3=iii3+1
           if (ddd  .gt. ddave*10.0d0) iii4=iii4+1
c          if (ddd .gt. ddave*2.0d0) then
c             write(6,8003)jsv,ja,ddd,ddave
c          endif  
         enddo   
      enddo   
!$omp end do
!$omp end parallel 
c
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-0320 ',ptime-tb
      tb=ptime
c
      write(6,8001)itemd,ddave,ddmax1,iibase,iii1,iii2,iii3,iii4
      write(6,*)'Max RN_DIA=',itemd,ddmax2
      if (ddmax2 .gt. 1.0d-10) then
        write(6,*)'Too large Max RN_DIA ?',itemd,ddmax2
      endif   
c
 8001 format('RN_AVE,sta=',I10,2F25.20,5I10)
 8002 format('RN_BASE,RN_DIA=',2I6,2F30.20)
 8003 format('RN_BASE,_AVE,_DIA=',2I6,3F30.20)
 8004 format('RN_AVE =',I6,F30.20)
c
c     write(6,*)'Stop manually'
c     stop
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Set the chemical potential
c       The total enectron number is given by rNelec
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_kr_set_chempot
      use elses_mod_phys_const, only : ev4au
      use elses_mod_tx,         only : jsei
      use elses_mod_md_dat,     only : itemd
      use elses_mod_js4jsv,     only : js4jsv 
      use elses_mod_orb1,       only : nval
      use elses_param_ctl_kr,   only : noav_kr_def, rNelec_kr_def
      use elses_arr_kry_glo,    only : rxmu, rE, rR,rFocc, nrecg
      use elses_mod_elec_cond,  only : temp_for_electron
c
c     use arr_kry_glo
c
c     include 'zconst.f'
c
      implicit none
      integer ipe,npe
      integer omp_get_thread_num
      integer omp_get_num_threads
      integer noav_kr
      real*8  rNelec
      real*8  tb, ptime
c
      integer iloop, nloopmax, irec
      integer jsv, js, nss, nval2, ja, nrecl
      real*8  xtemp, xbeta, ddd1, ddd2, ddemin, ddemax
      real*8  xmu0, xmumn, xmumx, xidos0, xidosE, tdossum, tensum
      real*8  rRd, rEd, yyy, xexp, xqinit, err
c

c
      noav_kr=noav_kr_def
      rNelec=rNelec_kr_def
c
      write(6,*)'@@ LSES_KR_SET_CHEMPOT:N_elec=',rNelec
c
      nloopmax=100
c     xtemp=0.0002d0
c     xtemp=0.002d0
c     xtemp=1.0d0/10.0d0
c     xtemp=1.0d0/20.0d0
c     xtemp=1.0d0/50.0d0
c     xtemp=1.0d0/100.0d0
c     xtemp=1.0d0/200.0d0
c     xtemp=1.0d0/2000.0d0
      xtemp=temp_for_electron
c        ----> Temperature in Fermi Distribution
      xbeta=1.0d0/xtemp
c
      write(6,*)'xtemp [au,eV]=',xtemp,xtemp*ev4au
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Trivial checking
c
      if (rNelec .lt. 0.99d0) then
        write(6,*)'@@ ERROR!:SET_KRY_CHEM_POT:Nelec=',rNelec
        stop
      endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call tclock(ptime)
      tb=ptime
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Find the min and max of energy levels
c
      jsv=1
      ja=1
      nrecl=nrecg(ja,jsv)
      ddemin=rE(1,ja,jsv)
      ddemax=rE(nrecl,ja,jsv)
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv,js,nss,nval2,ja,nrecl)
!$omp& private(ddd1,ddd2)
!$omp& reduction(max : ddemax)
!$omp& reduction(min : ddemin)
       ipe=0
       npe=0
c      ipe=omp_get_thread_num()+1
c      npe=omp_get_num_threads()
c      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
      do jsv=1,noav_kr
        js=js4jsv(jsv)
        nss=jsei(js)
        nval2=nval(nss)
        do ja=1,nval2
          nrecl=nrecg(ja,jsv)
          ddd1=rE(1,ja,jsv)
          ddd2=rE(nrecl,ja,jsv)
          ddemin=min(ddemin,ddd1)
          ddemax=max(ddemax,ddd2)
c         write(6,*)'ddemin,ddemax=',ddemin,ddemax
        enddo
      enddo
!$omp end do
!$omp end parallel 
c
      write(6,*)' ddemin [au,eV]=',ddemin,ddemin*ev4au
      write(6,*)' ddemax [au,eV]=',ddemax,ddemax*ev4au
      if (dabs(ddemin-ddemax) .lt. 1.0d-10) then
        write(6,*)'ERROR!:SET_KRY_CHEM_POT'
        stop
      endif
 1104 format('rEE=',3i5,f30.20)
c     stop
c  
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-410 ',ptime-tb
      tb=ptime
c
 1199 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ The bisection method
c
      xmu0=ddemin
      xmumn=ddemin
      xmumx=ddemax
c
      do 1000 iloop=1,nloopmax
        xidos0= 0.0d0
        xidosE= 0.0d0
        tdossum=0.0d0
        tensum =0.0d0
c
c      OMP-NOTE : Shared scalors : xmu0,xbeta
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv,js,nss,nval2,ja,irec,nrecl)
!$omp& private(rRd,rEd,yyy,xexp)
!$omp& reduction(+ : xidos0)
!$omp& reduction(+ : xidosE)
!$omp& reduction(+ : tdossum)
!$omp& reduction(+ : tensum)
       ipe=0
       npe=0
c      ipe=omp_get_thread_num()+1
c      npe=omp_get_num_threads()
c      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
        do jsv=1,noav_kr
          js=js4jsv(jsv)
          nss=jsei(js)
          nval2=nval(nss)
          do ja=1,nval2
            nrecl=nrecg(ja,jsv)
            do irec=nrecl,1,-1
              rRd=rR(irec,ja,jsv)
              rEd=rE(irec,ja,jsv)
              yyy=xbeta*(rEd-xmu0)
              if(yyy .gt. 100.d0) then
                xexp=0.0d0
              elseif(yyy.lt.-100.d0) then
                xexp=1.0d0
              else
                xexp=1.0d0/(dexp(yyy)+1.0d0)
              endif
              tdossum=tdossum +2.0d0*rRd
              tensum =tensum  +2.0d0*rRd*rEd
              xidos0=xidos0   +2.0d0*xexp*rRd
              xidosE=xidosE   +2.0d0*xexp*rRd*rEd
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel 
c
c
        xqinit=rNelec
        if (dabs(xqinit) .le. 1.0d-10) then
          write(6,*)'ERROR!(CHEM_POT):xqinit=',xqinit
          stop
        endif   
        err=(xqinit-xidos0)/xqinit
        if (iloop .le. 3) write(6,*)'Bisec.',iloop,xmu0,err
        if(abs(err) .le. 1.0d-12) then
          goto 1001
        endif   
c
        if(err .gt.  0.0d0) then
          xmumn=xmu0 
          xmumx=xmumx
          xmu0=(xmumn+xmumx)*0.5d0
        else
          xmumn=xmumn
          xmumx=xmu0 
          xmu0=(xmumn+xmumx)*0.5d0
        endif
c       
 1000  continue
       write(6,*)'ERROR!:chem_pot is not converged!!'
       stop
 1001  continue
       write(6,1091)itemd,iloop,xmu0,xmu0*ev4au
 1091  format('chem_pot [au,eV]=',2I10,2F10.5)
c      stop
c
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-420 ',ptime-tb
      tb=ptime
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Check the values and rFocc(nrecl,ja,jsv) are given.
c
       xidos0= 0.0d0
       xidosE= 0.0d0
       tdossum=0.0d0
       tensum =0.0d0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv,js,nss,nval2,ja,irec,nrecl)
!$omp& private(rRd,rEd,yyy,xexp)
!$omp& reduction(+ : xidos0)
!$omp& reduction(+ : xidosE)
!$omp& reduction(+ : tdossum)
!$omp& reduction(+ : tensum)
!$omp do schedule(static)
       do jsv=1,noav_kr
         js=js4jsv(jsv)
         nss=jsei(js)
         nval2=nval(nss)
         do ja=1,nval2
           nrecl=nrecg(ja,jsv)
           do irec=1,nrecl
             rRd=rR(irec,ja,jsv)
             rEd=rE(irec,ja,jsv)
             yyy=xbeta*(rEd-xmu0)
             if(yyy .gt. 100.d0) then
               xexp=0.0d0
             elseif(yyy.lt.-100.d0) then
               xexp=1.0d0
             else
               xexp=1.0d0/(exp(yyy)+1.0d0)
             endif
             rFocc(irec,ja,jsv)=xexp
             tdossum=tdossum +2.0d0*rRd
             tensum =tensum  +2.0d0*rRd*rEd
             xidos0=xidos0   +2.0d0*xexp*rRd
             xidosE=xidosE   +2.0d0*xexp*rRd*rEd
           enddo
         enddo
       enddo
!$omp end do
!$omp end parallel 
c
       write(6,*)'tdossum=',tdossum
       write(6,*)'tensum =',tensum
       write(6,*)' N_elec         =',xidos0
       write(6,*)' E_bs [au]      =',xidosE
       write(6,*)' E_bs [eV/atom] =',xidosE/dble(noav_kr)*ev4au
c
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-430 ',ptime-tb
      tb=ptime
c
       rxmu(:,:)=xmu0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
 9999 continue
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Calculate density matrix within Krylov space
c       < K_n | rho | j > = \sum_al c_{n al} f(E_al) < w_al | j >
c
c      NOTE: c_{n al} is re-calculated from the tridiagonal matrix.
c
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_cal_rho
      use elses_mod_ctrl,    only : i_verbose
      use elses_mod_tx,      only : jsei
      use elses_mod_js4jsv,  only : js4jsv 
      use elses_mod_orb1,    only : nval
      use elses_arr_kry_glo, only : nrecg, rA, rB, rRho, 
     +                             rW, rFocc

      use elses_mod_phys_const, only : ev4au

      use elses_mod_md_dat,     only : itemd


      use elses_param_ctl_kr,   only : noav_kr_def, rNelec_kr_def
      use elses_arr_kry_glo,    only : rxmu, rE, rR,rFocc, nrecg
c
c
c     use arr_kry_glo
c     include 'zconst.f'
c
      implicit none
      integer ipe, npe
      integer ierr, ishow2
      integer noav_kr, rNelec
      integer ishow
c
      integer jsv, js, nss, nval2, ja, nrecl, irec, irec1
      real*8  rWd, xexp, ddd, rFoccd, ddsum

      real*8,  allocatable :: rAtmp(:)
      real*8,  allocatable :: rBtmp(:)
      real*8,  allocatable :: rGe(:,:)
      real*8,  allocatable :: work_lapack(:,:)
c
c
      ishow=5
      noav_kr=noav_kr_def
      rNelec=rNelec_kr_def
      write(6,*)'@@ CAL_RHO_K_MP:rNelec=',rNelec
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv,js,ja,nss,nval2,irec,irec1,nrecl,ierr)
!$omp& private(rAtmp,rBtmp,rGe,work_lapack)
!$omp& private(rWd,xexp,ddd)
!$omp do schedule(static)
      do 1000 jsv=1,noav_kr
       js=js4jsv(jsv)
       nss=jsei(js)
       nval2=nval(nss)
       do 2000 ja=1,nval2
c 
c
          nrecl=nrecg(ja,jsv)
          if ((nrecl .le. 0) .or. (nrecl .gt. 10000000)) then
            write(6,*)'ERROR!(CAL_RHO_K_MP):nrecl=',nrecl
            stop
          endif
c
c         write(6,*)'nrecl=',jsv,ja,nrecl
c
          allocate(rAtmp(nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(CALRHO)rAtmp:ierr=',ierr
            write(6,*)'size:nrecl=',nrecl
            stop
          endif
          rAtmp(:)=0.0d0
c
          allocate(rBtmp(nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(CALRHO)rBtmp:ierr=',ierr
            stop
          endif
          rBtmp(:)=0.0d0
c
          allocate(rGe(nrecl,nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(CALRHO)rGe:ierr=',ierr
            stop
          endif
          rGe(:,:)=0.0d0
c
          allocate(work_lapack(1,2*nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(CALRHO)work_lapack:ierr=',ierr
            stop
          endif
          work_lapack(:,:)=0.0d0
c
        do irec=1,nrecl
           rAtmp(irec)=rA(irec,ja,jsv)
           rBtmp(irec)=rB(irec,ja,jsv)
c          write(6,*)'A,B=',irec,rAtmp(irec),rBtmp(irec)
        enddo   
c
        call DSTEV('V',nrecl,rAtmp,rBtmp,rGe,nrecl,work_lapack,ierr)
        if(ierr.ne.0) then
           write(6,*) 'in DSTEC, INFO=',ierr
           stop
        endif
c
        rRho(:,ja,jsv)=0.0d0
        do irec=1,nrecl
          rWd=rW(irec,ja,jsv)
          xexp=rFocc(irec,ja,jsv)
c            xexp :  FD function
          if ((xexp .gt. 1.0001d0) .or. (xexp .lt. -1.0d-5)) then
             write(6,*) 'ERROR!(CAL_RHO_K_MP)xexp=',xexp
             write(6,*) 'irec=',irec
             write(6,*) ' jsv=',jsv
             write(6,*) '  ja=',ja
             write(6,*) '  rE=',rE(irec,ja,jsv)
             stop
          endif   
c         write(6,*)'xexp,rWd=',irec,xexp,rWd
          do irec1=1,nrecl
            ddd=rGe(irec1,irec)*xexp*rWd
            rRho(irec1,ja,jsv)=rRho(irec1,ja,jsv)+ddd
          enddo
c
        enddo
c           < K_n | rho | j > = \sum_al c_{n al} f(E_al) < w_al | j >
c
        deallocate(rAtmp)
        deallocate(rBtmp)
        deallocate(rGe)
        deallocate(work_lapack)
c
 2000  continue
 1000 continue
!$omp end do
!$omp end parallel 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Check the electron number : sum_i < i | rho | i > 
c
      ddsum=0.0d0
      do jsv=1,noav_kr
        ishow2=0
        if (jsv .le. 2) ishow2=1
        js=js4jsv(jsv)
        nss=jsei(js)
        nval2=nval(nss)
        do ja=1,nval2
c         jj=js2j(ja,js)
c         idngl2=idngl(jj)
          nrecl=nrecg(ja,jsv)
c         xmu0=rxmu(ja,jsv)
c         write(6,*)'idngl2,xmu0=',idngl2,xmu0*ev4au
c         if (idngl2 .eq. 1) then
c           write(6,*)'  terminated basis'
c         endif   
          do irec=1,nrecl
            ddd=rRho(irec,ja,jsv)
c           rFoccd=rFocc(irec,ja,jsv)
            if (irec .eq. 1) then
               ddsum=ddsum+ddd
               if (i_verbose >= 1) then
                 if (jsv <= ishow) then 
                   write(6,3001)irec,ja,jsv,ddd
                 endif  
               endif
            endif   
          enddo
        enddo
      enddo  
 3001 format('rRho =',3I6,F20.10) 
c   
      ddsum=ddsum*2.0d0
      write(6,*)'N_elec (from rRho)=',ddsum
      if (dabs(ddsum-rNelec)/rNelec .gt. 1.0d-10) then
        write(6,*)'ERROR!(CAL_RHO):N_elec (defined)  =',rNelec
        write(6,*)'               :N_elec (from rRho)=',ddsum
        stop
      endif   
c
 9999 continue
c     stop
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Symmetrize dbij 
c      DBIJ2 is used as temporaly array
c          
c      dbij_asym (J) = \sum_{I,a,b} 
c           | < Ia | rho^(Jb) | Jb > -  < Jb | rho^(Ia) | Ia > |^2
c       J : js, ( not jsv )
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_sym_dbij
      use elses_mod_tx,       only : jsei
      use elses_mod_js4jsv,   only : js4jsv 
      use elses_mod_jsv4jsd,  only : jsv4jsd, njsd
      use elses_mod_orb1,     only : nval
      use elses_arr_dbij,     only : dbij
      use elses_arr_dbij2,    only : dbij2
      use elses_mod_multi,    only : ict4h, ict4l
      use elses_mod_orb2,     only : js2j
c
      use elses_mod_md_dat,   only : itemd
      use elses_param_ctl_kr, only : noav_kr_def
c
c     use arr_kry_glo
c     use arr_kry_hst
c     include 'zconst.f'
c     include 'zconst-VR.f'
c
      implicit none
      real*8  ptime, tb
      integer ipe,npe
      integer omp_get_thread_num
      integer omp_get_num_threads
c
      integer ishow2, iii, noav_kr
      integer jsv2, jsd2, jsd2b, js2, njsd2, nss2, nval2, ja2, ig2
      integer jsv1, jsd1,        js1, njsd1, nss1, nval1, ja1, ig1
      real*8  dsum1, dsum2,dlijd1,dlijd2,ddd,derr,ddd1,ddd2
      integer jsv2b
c
      noav_kr=noav_kr_def
      ishow2=0
      write(6,*)'@@ LSES_SYM_DBIJ : noav_kr=',noav_kr
c
c
      call tclock(ptime)
      tb=ptime
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    DBIJ is backuped as DBIJ2
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv2,js2,ja2,nss2,njsd2,nval2)
!$omp& private(jsv1,js1,ja1,nss1,njsd1,nval1,jsd1)
!$omp do schedule(static)
      do jsv2=1,noav_kr
        njsd2=njsd(jsv2,ict4h)
        js2=js4jsv(jsv2)
        nss2=jsei(js2)
        nval2=nval(nss2)
        dbij2(:,:,:,jsv2)=0.0d0
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
          do ja2=1,nval2
            do ja1=1,nval1
              dbij2(ja1,ja2,jsd1,jsv2)=dbij(ja1,ja2,jsd1,jsv2) 
            enddo
          enddo  
        enddo   
      enddo   
!$omp end do
!$omp end parallel 
c
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-0200 ',ptime-tb
      tb=ptime
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    DBIJ is symmetrized. 
c
      dsum1=0.0d0
      dsum2=0.0d0
      iii=0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv2,js2,ja2,nss2,njsd2,nval2,jsd2,jsv2b)
!$omp& private(jsv1,js1,ja1,nss1,njsd1,nval1,jsd1)
!$omp& private(ig1,ig2)
!$omp& private(dlijd1,dlijd2,ddd,derr)
!$omp& reduction(+ : iii)
!$omp& reduction(+ : dsum1)
!$omp& reduction(+ : dsum2)
!$omp do schedule(static)
      do jsv2=1,noav_kr
        njsd2=njsd(jsv2,ict4h)
        js2=js4jsv(jsv2)
        nss2=jsei(js2)
        nval2=nval(nss2)
c       ddsum=0.0d0
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          njsd1=njsd(jsv1,ict4h)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
c
c         jsd2=jsd4jsvf1(jsv2,jsv1)
          njsd1=njsd(jsv1,ict4h)
          do 133 jsd2=1,njsd1
           jsv2b=jsv4jsd(jsd2,jsv1)
           if (jsv2b .eq. jsv2) goto 134
  133     continue
          write(6,*)'ERROR!(JSD4JSVF1):JSV1,JSV2=',JSV1,JSV2
          stop
  134     continue
c
          if ((jsd2 .le. 0) .or. (jsd2 .gt. njsd1)) then
            write(6,*)'ERROR!(SYM_DBIJ):JSD2=',jsd2
            write(6,*)'jsd1,jsv2=',jsd1,jsv2
            write(6,*)'jsd2,jsv1=',jsd2,jsv1
            stop
          endif   
          jsv2b=jsv4jsd(jsd2,jsv1)
          if (jsv2b .ne. jsv2) then
            write(6,*)'ERROR!(SYM_DBIJ)'
            write(6,*)'jsd1,jsv2=',jsd1,jsv2
            write(6,*)'jsd2,jsv1=',jsd2,jsv1
            write(6,*)'    jsv2b=',jsv2b
            stop
          endif   
          do ja2=1,nval2
             ig2=js2j(ja2,js2)
            do ja1=1,nval1
              ig1=js2j(ja1,js1)
              dlijd1=dbij2(ja1,ja2,jsd1,jsv2)
              dlijd2=dbij2(ja2,ja1,jsd2,jsv1)
              ddd=(dlijd1+dlijd2)/2.0d0
              derr=dabs(dlijd1-dlijd2)
              dsum1=dsum1+derr
              dsum2=dsum2+derr/(dabs(ddd)+1.0d-6)
              iii=iii+1
c             if (ig1 .lt. ig2) then
c               dbij(ja1,ja2,jsd1,jsv2)=dlijd2
c             endif 
              dbij(ja1,ja2,jsd1,jsv2)=ddd
c             if (jsv2 .le. 3) then
c                write(6,129)ja1,ja2,jsv1,jsv2,ddd,derr
c             endif   
            enddo
          enddo  
        enddo   
      enddo   
!$omp end do
!$omp end parallel 
c
  129 format('DBIJ-SYM=',4i10,2f30.20)
c
      call tclock(ptime)
      write(6,*)' @@@OMP-TIME-0210 ',ptime-tb
      tb=ptime
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (iii .eq. 0) then
        write(6,*)'ERROR!(DBIJ-SYM):iii=',iii
        stop
      endif   
      write(6,139)itemd,iii,dsum1/dble(iii),dsum2/dble(iii)
  139 format('DBIJ-ERR=',2i15,2f30.20)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    DBIJ2 is re-defined as ( DBIJ - DBIJ_SYM ) (optional)
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv2,js2,ja2,nss2,njsd2,nval2)
!$omp& private(jsv1,js1,ja1,nss1,njsd1,nval1,jsd1)
!$omp& private(ddd1,ddd2)
!$omp do schedule(static)
      do jsv2=1,noav_kr
        njsd2=njsd(jsv2,ict4h)
        js2=js4jsv(jsv2)
        nss2=jsei(js2)
        nval2=nval(nss2)
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
          do ja2=1,nval2
            do ja1=1,nval1
              ddd1=dbij(ja1,ja2,jsd1,jsv2) 
              ddd2=dbij2(ja1,ja2,jsd1,jsv2) 
              dbij2(ja1,ja2,jsd1,jsv2)=ddd2-ddd1
            enddo
          enddo  
        enddo   
      enddo   
!$omp end do
!$omp end parallel 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    DBIJ2 is plotted (optional)
c
      do jsv2=1,noav_kr
        njsd2=njsd(jsv2,ict4h)
        js2=js4jsv(jsv2)
        nss2=jsei(js2)
        nval2=nval(nss2)
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
          do ja2=1,nval2
            do ja1=1,nval1
              ddd1=dbij(ja1,ja2,jsd1,jsv2) 
              ddd2=dbij2(ja1,ja2,jsd1,jsv2) 
              derr=dabs(ddd2/(ddd1+1.0d-10))
              if (jsv2 .le. ishow2) then
                write(6,137)ja1,ja2,js1,js2,ddd1,ddd2,derr
              endif   
            enddo
          enddo  
        enddo   
      enddo   
c
  137 format('DBIJ-ASYM-ATOM=',4I10,3F20.10)
c   
 9999 continue
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ TB energy
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_enetb(etbtmp)
      use elses_mod_ctrl,       only : i_verbose
      use elses_mod_sim_cell,   only : noa
      use elses_mod_tx,       only : jsei
      use elses_mod_js4jsv,   only : js4jsv 
      use elses_mod_jsv4jsd,  only : jsv4jsd, njsd
      use elses_mod_orb1,     only : nval
      use elses_arr_dbij,     only : dbij
      use elses_arr_dhij,     only : dhij
      use elses_mod_multi,    only : ict4h, ict4l
      use elses_mod_orb2,     only : js2j
c
      use elses_param_ctl_kr, only : noav_kr_def
c
      implicit none
      integer ipe, npe
      real*8  etbtmp

      integer noav_kr
      integer jsv2,      js2, ja2, njsd2, nss, nval2, jg
      integer jsv1,jsd1, js1, ja1,             nval1, ig
      real*8  dsum, dsum1, dlijd, ahijd
c
      noav_kr=noav_kr_def
      if (i_verbose >= 1) then
        write(6,*)'@@ ENETB_KR:noav=',noav_kr
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Trivial error checking (optional)
c
      if ((noav_kr <= 0) .or. (noav_kr > noa)) then
        write(6,*)'ERROR(ENETB_KR):noav=',noav_kr
        stop
      endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dsum=0.0D0
      dsum1=0.0D0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(jsv2,js2,ja2,njsd2,nss,nval2)
!$omp& private(jsv1,js1,ja1,jsd1,nval1)
!$omp& private(ig,jg)
!$omp& private(dlijd,ahijd)
!$omp& reduction(+ : dsum)
!$omp& reduction(+ : dsum1)
!$omp do schedule(static)
      do 1000 jsv2=1,noav_kr
         js2=js4jsv(jsv2)
c        if ((js2 .le. 0) .or. (js2 .gt. noa)) then
c          write(6,*)'error!(enetb):jsv2,js2=',jsv2,js2
c          stop
c        endif
         njsd2=njsd(jsv2,ict4h)
         nss=jsei(js2)
         nval2=nval(nss)
      do 1100 jsd1=1,njsd2
         jsv1=jsv4jsd(jsd1,jsv2)
         if ((jsv1 .le. 0) .or. (jsv1 .gt. noav_kr)) then
           write(6,*)'error!(enetb):jsv1,noav=',jsv1,noav_kr
           stop
         endif
         js1=js4jsv(jsv1)
c        if ((js1 .le. 0) .or. (js1 .gt. noa)) then
c          write(6,*)'error!(enetb):jsv1,js1=',jsv1,js1
c          stop
c        endif
         nss=jsei(js1)
         nval1=nval(nss)
      do ja2=1,nval2
         jg=js2j(ja2,js2)
         do ja1=1,nval1
           ig=js2j(ja1,js1)
           dlijd=dbij(ja1,ja2,jsd1,jsv2)
           ahijd=dhij(ja1,ja2,jsd1,jsv2)
           dsum=dsum+2.0d0*dlijd*ahijd
           if (ig .eq. jg) then
             dsum1=dsum1+2.0d0*dlijd
c            write(6,*)'ig,jg=',ig,jg
             if (js1 .ne. js2) then
               write(6,*)'error!:js1,js2=',js1,js2
               stop
             endif
           endif
         enddo
      enddo   
 1100 continue
 1000 continue
!$omp end do
!$omp end parallel 
c
      etbtmp=dsum
c
      WRITE(6,*)'ETB,N_elec  =',DSUM,DSUM1
c            ---> energy and electron number
c
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     @@ Booking the neighboring list 
c
c         Output : jsv4jsd(noao,noav)
c                : njsd(noav,0:nncut)
c
c         non-order-N calculation!!!
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_es_bkmat
      use elses_mod_ctrl,       only : i_verbose
      use elses_mod_sel_sys,    only : c_system, r_cut_book
      use elses_mod_noav,       only : noav,noao,nncut
!     use elses_mod_sim_cell,   only : iperiodic, noa, ax, ay, az 
      use elses_mod_sim_cell,   only : noa, ax, ay, az,
     +                                 i_pbc_x, i_pbc_y,i_pbc_z 
      use elses_mod_tx,         only : tx, ty, tz
      use elses_mod_jsv4jsd,    only : jsv4jsd, njsd
      use elses_mod_js4jsv,     only : js4jsv,jsv4js
      use M_md_booking,         only : get_max_size_of_booking_list
c
      implicit none
      real*8  rcut_bk
c
      integer ishow,ict
      real*8  rcut1
c
      integer ipe,npe
      integer js1,js2,jsv1,jsv2,js1d,js1b,jsd
      real*8  ddd,dvecx,dvecy,dvecz,dxc,dyc,dzc,dd,w
      integer max_size
c
      ishow=5
c
      ICT=1
c
      write(6,*)'@@@ LSES_ES_BKMAT'
      write(6,*)'   noav,noao,nncut=',noav,noao,nncut
      write(6,*)'c_system=',c_system
      write(6,*)'r_cut_book =',r_cut_book
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Get the maximum size of the booking list
c
      max_size=0
      call get_max_size_of_booking_list(r_cut_book,max_size)
      write(6,*)' max_size for booking list=',max_size
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Check the parapmeters
c
      if (nncut .ne. 2) then 
        write(6,*)'ERROR!(BKMATV)nncut=',nncut
        stop
      endif
c
      rcut_bk=r_cut_book
      if (rcut_bk .le. 1.0d0) then 
        write(6,*)'ERROR!(BKMATV)rcut_bk=',rcut_bk
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      RCUT1=RCUT_BK*1.001d0
c         radius  for cutoff
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(dd,dvecx,dvecy,dvecz,w)
!$omp& private(jsv1,jsv2,js1,js2,js1d,js1b,jsd)
c     ipe=omp_get_thread_num()+1
c     npe=omp_get_num_threads()
c     write(6,*)'ipe,npe=',ipe,npe
!$omp  do schedule(static)
      DO 200 JSV2=1,NOAV
c
        jsv4jsd(:,jsv2)=0
c         
        JS2=JS4JSV(JSV2)
        JSV4JSD(1,JSV2)=JSV2
        NJSD(JSV2,0)=1
c          ----> zero-th shell 
c
        JS1D=1
        JS1B=1
        ICT=1
c   
        DO 220 JSV1=1,NOAV
          JS1=JS4JSV(JSV1)
          dd=1.0d0
          dvecx=tx(js2)-tx(js1)
          dvecy=ty(js2)-ty(js1)
          dvecz=tz(js2)-tz(js1)

          if (i_pbc_x == 1)dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
          if (i_pbc_y == 1)dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
          if (i_pbc_z == 1)dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0

!         if (iperiodic .eq. 1) then
!           dvecx=dvecx+0.5d0*dd
!           dvecx=dmod(dvecx+dd,dd)-0.5d0*dd
!           dvecy=dvecy+0.5d0*dd
!           dvecy=dmod(dvecy+dd,dd)-0.5d0*dd
!           dvecz=dvecz+0.5d0*dd
!           dvecz=dmod(dvecz+dd,dd)-0.5d0*dd
!         endif

          dvecx=dvecx*ax
          dvecy=dvecy*ay
          dvecz=dvecz*az
          w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
c
c         write(6,*)'jsv2,jsv1,w,rcut=',jsv2,jsv1,w,rcut1
c
          IF (JSV1 .EQ. JSV2) GOTO 220
          IF (W .GT. RCUT1) GOTO 220
          JS1D=JS1D+1
c
          IF (JS1D .GT. NOAO) THEN
             WRITE(6,*)'ERROR!(BKMATV2):JS1D,NOAO=',JS1D,NOAO
             STOP
          ENDIF   
          JSV4JSD(JS1D,JSV2)=JSV1
  220   CONTINUE
  221   CONTINUE
        NJSD(JSV2,ICT)=JS1D
c
  210 CONTINUE
  200 CONTINUE
!$omp end do
!$omp end parallel 
c     
      WRITE(6,*)'loop200 is ended'
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Plot the result (optional)
c
      if (i_verbose .ge. 1) then
        do jsv2=1,noav
          if (jsv2 .le. ishow) then
            ict=1
            write(6,*)'jsv,njsd=',jsv2,njsd(jsv2,ict)
          endif   
        enddo
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 9999 CONTINUE
c
      END
