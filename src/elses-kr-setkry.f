!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-kr-setkry.f @@@@@
czzz  @@@@@ 2007/06/07 @@@@@
c2006cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0728: Prepared: elses_set_kry (NT116p33) 
c                 based on setkry_b in export/md-nrl-002v01
c0830: A modi.:Now 'Terminated SEED=' is not plotted for saving disk
c2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0109: Becomes a member of PT41, without any modifcation.
c0607: Several write commands are commented out with '!dm' mark.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_kry(imode,iselect)
      use elses_mod_phys_const, only : ev4au
c     use elses_mod_sim_cell,   only : iperiodic, ax, ay, az
      use elses_mod_tx,         only : jsei
      use elses_mod_jsv4jsd,    only : jsv4jsd,njsd
      use elses_mod_js4jsv,     only : js4jsv
      use elses_mod_orb1,       only : nval
      use elses_mod_orb2,       only : js2j, idngl
      use elses_arr_dbij,       only : dbij
      use elses_arr_dhij,       only : dhij
      use elses_mod_multi,      only : ict4h, ict4l, intf
      use elses_param_ctl_kr, only : iKnstr, noav_kr_def, nreclc_def,
     +                  nkene_def,ddemin_def,ddemax_def
      use elses_arr_kry_glo,    only : rA, rB, rE, rR, rW, rRho, 
     +          rRNave, rRNdia, nrecg, rKnstr, noak_str, jsv4jsk_str
c
c     use arr_kry_glo
c     use arr_kry_glo_b
c     include 'zconst.f'
c     include 'zconst-VR.f'
c
      implicit none
c
      integer, allocatable :: jsv4jsk(:)
      integer, allocatable :: jsk4jsv(:)
      integer, allocatable :: jjkset(:)
      real*8,  allocatable :: rU1(:)
      real*8,  allocatable :: rU2(:)
      real*8,  allocatable :: w0xl(:)
      real*8,  allocatable :: rAtmp(:)
      real*8,  allocatable :: rBtmp(:)
      real*8,  allocatable :: acwrk1(:)
      real*8,  allocatable :: acwrk2(:)
      real*8,  allocatable :: rGe(:,:)
      real*8,  allocatable :: work_lapack(:,:)
      real*8,  allocatable :: rKntmp(:,:)
      real*8,  allocatable :: rLntmp(:,:)
c
      complex*16,  allocatable :: cZ_k(:)
      complex*16,  allocatable :: cGn1tmp(:,:)
      complex*16,  allocatable :: cRNtmp(:,:)
c   
      integer, allocatable :: jjk4jjt(:)
      integer, allocatable :: jjt4jjk(:)
      integer, allocatable :: jsk4jjk(:)
      integer, allocatable :: ja4jjk(:)
c
      integer, allocatable :: jsv4jrun(:)
c
      integer  imode, iselect, ierr
      integer  iRN_calc, ishow2
c
      real*8   tb, ptime
      integer  ipe,npe
      integer  mr, nvalmx
      integer  omp_get_thread_num
      integer  omp_get_num_threads
c
      integer  nreclc, noav_kr,nrecl, irec
      integer  jj, jrun, njrun, nreclt, iibase, iii
      integer  kene, nkene, ial
      real*8   ddemin, ddemax, de, ene_gamma
c
      integer  noak, jsv4seed, js4seed, ja4seed
      integer  jjk4seed, jj4seed, idngl2
      integer  jsv,  js,        nss,         ja
      integer  jsv1, js1, jsd1, nss1, nval1, ja1
      integer  jsv2, js2,       nss2, nval2, ja2
      integer  jsk, jsk1, jsk2, ntot, njsd2 
      integer  jjt, njjt, jjk, jjk1, jjk2, jjkd, jjkset1, jjkset2
      real*8   ddd, ddsum, dsum1, dbigd, dlijd, acwrk2d
      complex*16  cddsum, cdd1, cdd2
c
      write(6,*)'@@ LSES_SET_KRY:imode,iKnstr=',imode,iKnstr
      if (iKnstr .eq. 1) write(6,*)'  NOTE:K_n is stored!!!!'
c
      mr=0
      nvalmx=10
c        : just used for trivial checking 
c
      nreclc =nreclc_def
      noav_kr=noav_kr_def
c
      write(6,*)'  nreclc=',nreclc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate (jsv4jrun(noav_kr),stat=ierr)
      if(ierr .ne. 0) then
        write(6,*)'alloc. error!(SETKRY)jsv4jsk:ierr=',ierr
        stop
      endif
      jsv4jrun(:)=0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Check that jsv4jsk is already allocated ?
c
      if (.not. allocated(jsv4jsk_str)) then
        write(6,*)'ERROR!:Not yet allocated: JSV4JSK'
        stop
      endif 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (iKnstr .eq. 1) then
        if (.not. allocated(rKnstr)) then
          write(6,*)'ERROR(SETKRY):rKnstr is not allocated!!'
          stop
        else
          write(6,*)'...rKnstr is already allocated. OK!'
        endif 
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Krylov subspace is calculated for ALL seeds
c
      if ((iselect .eq. 0) .or. (iselect .eq. 2)) then
        jj=0 
        do jsv=1,noav_kr
         jj=jj+1
         jsv4jrun(jj)=jsv
        enddo
        njrun=jj
        if (iselect .eq. 0) nrecg(:,:)=nreclc
        if (noav_kr .ne. njrun) then
          write(6,*)'ERROR!(SETKRY)njrun,noav_kr=',njrun,noav_kr
          stop
        endif   
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Krylov subspace is calculated for SELECTED seeds
c
      if (iselect .eq. 1) then
        jj=0 
        nreclt=0
        iibase=0
        do jsv=1,noav_kr
         js=js4jsv(jsv)
         nss=jsei(js)
         nval2=nval(nss)
         iii=0
         do ja=1,nval2
           nrecl=nrecg(ja,jsv)
           if (nrecl .gt. nreclc) iii=iii+1
         enddo  
         if (iii .ne. 0) then
           jj=jj+1
           jsv4jrun(jj)=jsv
           iibase=iibase+nval2
           do ja=1,nval2
             nrecl=nrecg(ja,jsv)
             nreclt=nreclt+nrecl
           enddo  
         endif  
        enddo
        njrun=jj
        if (njrun .eq. 0) then
          write(6,*)'...skipped (SETKRY):njrun=',njrun
          goto 9999
        endif   
        if (iibase .eq. 0) then
          write(6,*)'ERROR!(SETKRY)njrun,iibase=',njrun,iibase
          stop
        endif   
        ddd=dble(nreclt)/dble(iibase)
        write(6,*)'njrun=',njrun
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      if (imode .eq. 2) then
         write(6,*)'DBIJ is cleared'
         dbij(:,:,:,:)=0.0d0
      endif   
c
      iRN_calc=0
      if (imode .eq. 1) iRN_calc=1
c
      ddemin=ddemin_def
      ddemax=ddemax_def
      nkene=nkene_def
      if (nkene .eq. 0) then
         write(6,*)'ERROR!(SETKRY):NKENE=',nkene
         stop
      endif   
c
      de=(ddemax-ddemin)/dble(nkene-1)
      ene_gamma=de*3.0d0
c
      if (iRN_calc .eq. 1) then
         write(6,*)'Residual Norm will be calculated'
         write(6,*)'  ddemin[eV]=',ddemin*ev4au
         write(6,*)'  ddemax[eV]=',ddemax*ev4au
         write(6,*)'      de[eV]=',de*ev4au
         write(6,*)'   nkene    =',nkene
         write(6,*)'   gamma[eV]=',ene_gamma*ev4au
      endif   
c
      call tclock(ptime)
      tb=ptime
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(ishow2,ierr)
!$omp& private(irec,nrecl)
!$omp& private(jrun)
!$omp& private(jsv4seed,js4seed,ja4seed,jjk4seed,jj4seed,noak,ntot)
!$omp& private(ddd,ddsum)
!$omp& private(cdd1,cdd2,cddsum)
!$omp& private(dbigd,acwrk2d)
!$omp& private(jsv4jsk,jsk4jsv,jjkset)
!$omp& private(jjk4jjt,jjt4jjk,jsk4jjk,ja4jjk)
!$omp& private(rU1,rU2,acwrk1,acwrk2,w0xl)
!$omp& private(rAtmp,rBtmp,rGe,work_lapack)
!$omp& private(rKntmp,rLntmp,cZ_k,cGn1tmp,cRNtmp)
!$omp& private(jsk,jsv,jj,js,idngl2,jjk,kene,ial)
!$omp& private(jsk1,jsk2,jsv1,jsv2)
!$omp& private(js1,js2,nss,nss1,nss2,nval1,nval2,njsd2,jjkset1,jjkset2)
!$omp& private(jsd1,jjk1,jjk2,jjkd,jjt,njjt)
!$omp& private(ja,ja1,ja2)
       ipe=0
       npe=0
c      ipe=omp_get_thread_num()+1
c      npe=omp_get_num_threads()
c      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
      do 1000 jrun=1,njrun
c
       jsv4seed=jsv4jrun(jrun)
       if ((jsv4seed .le. 0) .or. (jsv4seed .gt. noav_kr)) then
         write(6,*)'ERROR!(SETKRY):jrun,jsv4seed=',jrun,jsv4seed
         stop
       endif   
c
c      acutkry_tmp=rcut_kry(jsv4seed)
c      if (acutkry_tmp .le. 1.0d-6) then
c        write(6,*)'ERROR!(SETKRY):acutkry_tmp=',acutkry_tmp
c        stop
c      endif   
c
       noak=noak_str(jsv4seed)
       if ((noak .le. 0) .or. (noak .gt. noav_kr)) then
         write(6,*)'ERROR!(SETKRY):noak=',noak
         stop
       endif   
c
       ishow2=0
       js4seed=js4jsv(jsv4seed)
       if (jsv4seed .le. 3) ishow2=1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Initial clearance (within loop 1000)
c
       allocate (jsv4jsk(noav_kr),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)jsv4jsk:ierr=',ierr
         stop
       endif
       jsv4jsk(:)=0
c
       allocate (jsk4jsv(noav_kr),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)jsk4jsv:ierr=',ierr
         stop
       endif
       jsk4jsv(:)=0
c
       allocate (jjkset(noav_kr),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)jjkset:ierr=',ierr
         stop
       endif
       jjkset(:)=0
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Booking atoms for Krylov subspace (within loop 1000)
c          ( jsk <--> jsv )
c
       jsv4jsk(1:noak)=jsv4jsk_str(1:noak,jsv4seed)
       do jsk=1, noak
          jsv=jsv4jsk(jsk)
          if ((jsv .le. 0) .or. (jsv .gt. noav_kr)) then
            write(6,*)'ERROR!(SETKRY):jsk,jsv=',jsk,jsv
            stop
          endif   
          jsk4jsv(jsv)=jsk
       enddo   
c
       if (ishow2 .eq. 1) then
          write(6,*)'SEED,NOAK=',jsv4seed,noak
       endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Determine the matrix size ( within loop 1000 )
c         ----> jjkset(1:noak), ntot
c          ( Trivial checkings are also done )
c            The elements for the jsk-th atom are
c              jjk = jjkset(jsk)+1, ..... , jjkset(jsk) + nval2
c
c
       jj=0
       do 1100 jsk=1,noak
          jjkset(jsk)=jj
          jsv=jsv4jsk(jsk)
          if ((jsv .le. 0) .or. (jsv .gt. noav_kr)) then
            write(6,*)'ERROR!(SETKRY:1100):JSK,JSV=',jsk,jsv
            stop
          endif   
          js=js4jsv(jsv)
c         if ((js .le. 0) .or. (js .gt. noa)) then
c           write(6,*)'ERROR!(SETKRY:1100):JS,JSV=',js,jsv
c           stop
c         endif   
          nss=jsei(js)
c         if ((nss .le. 0) .or. (nss .gt. nos)) then
c           write(6,*)'ERROR!(SETKRY:1100):JS,NSS=',js,nss
c           stop
c         endif   
          nval2=nval(nss)
          if ((nval2 .le. 0) .or. (nval2 .gt. nvalmx)) then
            write(6,*)'ERROR!(SETKRY:1100):JS,NVAL2=',js,nval2
            stop
          endif   
          jj=jj+nval2
 1100  continue
       ntot=jj
       if (ishow2 .eq. 1) then
          write(6,*)'SEED,NOAK,NTOT=',jsv4seed,noak,ntot
       endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Count the number of the target bases |i > ( within loop 1000 )
c
       jjt=0
       jsk2=1
       jsv2=jsv4jsk(jsk2)
       if (jsv2 .ne. jsv4seed) then
         write(6,*)'ERROR!(SETKRY):jsv2,jsv4seed=',jsv2,jsv4seed
         stop
       endif   
       if(intf(jsv2).eq.1) then
          njsd2=njsd(jsv2,ict4l)
          else
          njsd2=njsd(jsv2,ict4h)
       endif
       do jsd1=1,njsd2
         jsv1=jsv4jsd(jsd1,jsv2)
         jsk1=jsk4jsv(jsv1)
         if ((jsk1 .eq. 0) .or. (jsk1 .gt. noak)) then
           write(6,*)'jsk1,jsv1=',jsk1,jsv1
           write(6,*)'jsk2,jsv2=',jsk2,jsv2
           write(6,*)'njsd2    =',njsd2
           write(6,*)'     jsd1=',jsd1
           write(6,*)'ERROR!(SETKRY:0):jjt'
           stop
         endif   
         js1=js4jsv(jsv1)
         nss1=jsei(js1)
         nval1=nval(nss1)
         jjt=jjt+nval1
       enddo  
       njjt=jjt
c      write(6,*)'njsd2,njjt=',njsd2,njjt
       if ((njjt .le. 0) .or. (njjt .gt. ntot)) then
         write(6,*)'Error!(SETKRY):njjt,ntot=',njjt,ntot
         stop
       endif   
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Allocate the vectors ( within loop 1000 )
c
       if ((ntot .le. 0) .or. (ntot .gt. nvalmx*noav_kr)) then
         write(6,*)'ERROR?!(SETKRY):ntot=',ntot
         stop
       endif   
c
       if ((njjt .le. 0) .or. (njjt .gt. ntot)) then
         write(6,*)'ERROR?!(SETKRY):njjt=',njjt
         stop
       endif   
c
       allocate(rU1(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)rU1:ierr=',ierr
         stop
       endif
       rU1(:)=0.0d0
c
       allocate(rU2(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)rU2:ierr=',ierr
         stop
       endif
       rU2(:)=0.0d0
c
       allocate(acwrk1(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)acwrk1:ierr=',ierr
         stop
       endif
c
       allocate(acwrk2(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)acwrk2:ierr=',ierr
         stop
       endif
c
       allocate(w0xl(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)w0xl:ierr=',ierr
         stop
       endif
       w0xl(:)=0.0d0
c
       allocate(jjk4jjt(njjt),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)jjk4jjt:ierr=',ierr
         stop
       endif
c
       allocate(jjt4jjk(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)jjt4jjk:ierr=',ierr
         stop
       endif
c
       allocate(jsk4jjk(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)jsk4jjk:ierr=',ierr
         stop
       endif
c
       allocate(ja4jjk(ntot),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(SETKRY)ja4jjk:ierr=',ierr
         stop
       endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Define jsk4jjk(1:ntot),ja4jjk(1:ntot) (within loop 1000)
c
       jjk=0
       jj=0
       do jsk=1,noak
          jjkset2=jjkset(jsk)
          jsv=jsv4jsk(jsk)
          js=js4jsv(jsv)
          nss=jsei(js)
          nval2=nval(nss)
          do ja=1,nval2
           jjk=jjkset2+ja
           jj=jj+1
           if (jj .gt. ntot) then
             write(6,*)'ERROR(SETKRY:4):jj,ntot=',jj,ntot
             stop
           endif
           if (jjk .ne. jj) then
             write(6,*)'ERROR(SETKRY:4):jjk,jj=',jjk,jj
             stop
           endif
           jsk4jjk(jj)=jsk
           ja4jjk(jj)=ja
          enddo
       enddo   
       if (ntot .ne. jj) then
         write(6,*)'ERROR(SETKRY:3):jj,ntot=',jj,ntot
         stop
       endif   
c      if (ishow2 .eq. 1) then
c        do jjk=1,ntot
c          jsk=jsk4jjk(jjk)
c          ja=ja4jjk(jjk)
c          write(6,*)'jjk=',jjk,jsk,ja
c        enddo  
c     endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Determine the target bases |i > ( within loop 1000 )
c         for storing < i | K_n > and < i | L_n >
c
       jsk2=1
       jsv2=jsv4jsk(jsk2)
       if (jsv2 .ne. jsv4seed) then
         write(6,*)'ERROR!(SETKRY):jsv2,jsv4seed=',jsv2,jsv4seed
         stop
       endif   
c
       jjt=0
       if(intf(jsv2).eq.1) then
          njsd2=njsd(jsv2,ict4l)
          else
          njsd2=njsd(jsv2,ict4h)
       endif
c     write(6,*)'njsd2=',njsd2
       do jsd1=1,njsd2
         jsv1=jsv4jsd(jsd1,jsv2)
         jsk1=jsk4jsv(jsv1)
         if ((jsk1 .eq. 0) .or. (jsk1 .gt. noak)) then
           write(6,*)'jsk1,jsv1=',jsk1,jsv1
           write(6,*)'jsk2,jsv2=',jsk2,jsv2
           write(6,*)'njsd2    =',njsd2
           write(6,*)'     jsd1=',jsd1
           write(6,*)'ERROR!(SETKRY):jjt'
           stop
         else  
           js1=js4jsv(jsv1)
           jjkset1=jjkset(jsk1)
           nss1=jsei(js1)
           nval1=nval(nss1)
           do ja1=1,nval1
             jjt=jjt+1
             jjk=jjkset1+ja1
             if (jjt .gt. njjt) then
              write(6,*)'ERROR(SETKRY:98):jjt,njjt=',jjt,njjt
              write(6,*)'jsv4seed,njsd2=',jsv4seed,njsd2
              stop
             endif
             jjk4jjt(jjt)=jjk
             jjt4jjk(jjk)=jjt
c            write(6,*)'book:jjt,jjk=',jjt,jjk
           enddo   
         endif   
       enddo  
       if (jjt .ne. njjt) then
         write(6,*)'ERROR(SETKRY):njjt=',jjt,njjt
         stop
       endif
c
c      if (ishow2 .eq. 1) then
c         do jjt=1,njjt
c           jjk=jjk4jjt(jjt)
c           jsk=jsk4jjk(jjk)
c           ja=ja4jjk(jjk)
c           write(6,*)'target:',jjt,jjk,jsk,ja
c         enddo  
c     endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Lanczos loop start ( within loop 1000)
c
       do 2000 ja4seed=1,nval2
c         jj=js4jsv(jsv4seed)
          jjk4seed=ja4seed
          jj4seed=js2j(ja4seed,js4seed)
          idngl2=idngl(jj4seed)
          nrecl=nrecg(ja4seed,jsv4seed)
          if (ishow2 .eq. 1) then
!dm          write(6,*)'IPE,SEED Num. Idngl=',ipe,jj4seed,idngl2
          endif  
          if (idngl2 .eq. 1) then
             write(6,*)'TERMINATED BASE!!:',jj4seed
          endif   
c
          rU2(:)=0.0d0
          rU1(:)=0.0d0
          rU2(jjk4seed)=1.0d0
c           ----> Set the seed vector : | K_1 > = |j >
c
          if ((nrecl .le. 0) .or. (nrecl .gt. 1010)) then
            write(6,*)'ERROR?!(SETKRY):nrecl=',nrecl
            stop
          endif   
c
          if ((nkene .le. 0) .or. (nkene .gt. 1010)) then
            write(6,*)'ERROR?!(SETKRY):nkene=',nrecl
            stop
          endif   
c
          if ((njjt .le. 0) .or. (njjt .gt. ntot)) then
            write(6,*)'ERROR?!(SETKRY):njjt=',njjt
            stop
          endif   
c
          allocate(rAtmp(nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)rAtmp:ierr=',ierr
            write(6,*)'size:nrecl=',nrecl
            stop
          endif
          rAtmp(:)=10.0d0
c            ----> Dummy values are essential for terminated case !!
c
          allocate(rBtmp(nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)rBtmp:ierr=',ierr
            stop
          endif
          rBtmp(:)=0.0d0
c            ----> Dummy values are essential for terminated case !!
c
          allocate(rKntmp(njjt,nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)rKntmp:ierr=',ierr
            stop
          endif
          rKntmp(:,:)=0.0d0
c
          allocate(rLntmp(njjt,nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)rLntmp:ierr=',ierr
            stop
          endif
          rLntmp(:,:)=0.0d0
c
          allocate(cGn1tmp(nrecl,nkene),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)cGn1tmp:ierr=',ierr
            stop
          endif
          cGn1tmp(:,:)=0.0d0
c
          allocate(cRNtmp(njjt,nkene),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)cGn1tmp:ierr=',ierr
            stop
          endif
          cRNtmp(:,:)=0.0d0
c
          allocate(cZ_k(nkene),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)cZ_k:ierr=',ierr
            stop
          endif
          cZ_k(:)=0.0d0
c
          if (iKnstr .eq. 1) then
             if (imode .eq. 2) goto 3001 
          endif   
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          do 3000 irec=1,nrecl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ LANCZOS : Matrix-vector multiplication (in loop 3000)
c                          | K_n > : rU2(ntot)   [input]
c                        H | K_n > : w0xl(ntot)  [output]
c        a_n =   < K_n | H | K_n > : rAtmp(irec) [output]
c
            acwrk1(:)=rU2(:)
            acwrk2(:)=0.0d0
            do jsk2=1,noak
              jsv2=jsv4jsk(jsk2)
              js2=js4jsv(jsv2)
              nss2=jsei(js2)
              nval2=nval(nss2)
              if(intf(jsv2).eq.1) then
                njsd2=njsd(jsv2,ict4l)
              else
                njsd2=njsd(jsv2,ict4h)
              endif
              jjkset2=jjkset(jsk2)
              do jsd1=1,njsd2
                jsv1=jsv4jsd(jsd1,jsv2)
                jsk1=jsk4jsv(jsv1)
c               if (jsk1 .eq. 0) then
c                 write(6,*)'jsk1,jsv1=',jsk1,jsv1
c                 write(6,*)'jsk2,jsv2=',jsk2,jsv2
c                 write(6,*)'njsd2    =',njsd2
c                 write(6,*)'     jsd1=',jsd1
c                 write(6,*)'Is this modified Krylov?'
c                 write(6,*)'Stop manually'
c                 stop
c               endif   
                if (jsk1 .ne. 0) then
                  js1=js4jsv(jsv1)
                  jjkset1=jjkset(jsk1)
                  nss1=jsei(js1)
                  nval1=nval(nss1)
                  do ja2=1,nval2
                   jjk2=jjkset2+ja2
                   do ja1=1,nval1
                    jjk1=jjkset1+ja1
c                   dbigd=dlij(ja1,ja2,jsd1,jsv2)
c                     ---> H + 2 eta rho_{PT}
                    dbigd=dhij(ja1,ja2,jsd1,jsv2)
c                     ---> H 
                    acwrk2d=acwrk2(jjk2)
                    acwrk2(jjk2)=acwrk2d+dbigd*acwrk1(jjk1)
                   enddo
                  enddo
                endif  
              enddo   
            enddo   
            w0xl(:)=acwrk2(:)
c
            ddd=dot_product(rU2(1:ntot),w0xl(1:ntot))
c           write(6,*)'ddd=',ddd
            rAtmp(irec)=ddd
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ LANCZOS : Store <i | K_n > and < i | L_n > = < i | H | K_n >
c
            do jjt=1,njjt
              jjk=jjk4jjt(jjt)
              if ((jjk .le. 0) .or. (jjk .gt. ntot)) then
                write(6,*)'ERROR!(SETKRY):store K_n, L_n'
                write(6,*)'jjt,jjk,ntot=',jjt,jjk,ntot
                stop
              endif
              rKntmp(jjt,irec)=rU2(jjk)
              rLntmp(jjt,irec)=w0xl(jjk)
            enddo   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ LANCZOS : Generate non-normalized vector | L_n >
c          for n = 2,3...
c            | L_n > = H | K_n > - a_n | K_n > - b_{n-1} | K_{n-1} >    
c          for n = 1
c            | L_n > = H | K_n > - a_n | K_n > 
c
c              | K_n > : rU2(ntot)   [input]
c          | K_{n-1} > : rU1(ntot)   [input]
c            H | K_n > : w0xl(ntot)  [input,will be modified]
c              | L_n > : w0xl(ntot)  [output]
c
           w0xl(:)=w0xl(:)-rAtmp(irec)*rU2(:)
           if (irec .ge. 2) then
             w0xl(:)=w0xl(:)-rBtmp(irec-1)*rU1(:)
           endif   
c
           if (irec .eq. nrecl) goto 3000
c              ---> In the last Lanczos loop, b_n is not needed.
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ LANCZOS : Generate coefficient b_n 
c                       | L_n > : w0xl(ntot)  [input]
c           b_n = sqrt{ < L_n | L_n > }  : rBtmp [output]
c
            ddd=dot_product(w0xl(1:ntot),w0xl(1:ntot))
            if (dabs(ddd) .lt. 1.0d-15) then
              nrecl=irec
              nrecg(ja4seed,jsv4seed)=nrecl
c             write(6,*)'Terminated SEED=',ja4seed,jsv4seed,nrecl
              goto 3001
            endif   
            rBtmp(irec)=dsqrt(ddd)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ LANCZOS : Update the vector for n --> n+1
c         | K_{n+1} > = (1/b_n) | L_n > : rU2(ntot)  [output]
c
            ddd=1.0d0/rBtmp(irec)
            rU1(:)=rU2(:)
            rU2(:)=w0xl(:)*ddd
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 3000     continue
 3001     continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ AFTER-LANCZOS : Allocation 
c        nrecl is changed in the terminated case.
c
          allocate(rGe(nrecl,nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)rGe:ierr=',ierr
            stop
          endif
          rGe(:,:)=0.0d0
c
          allocate(work_lapack(1,2*nrecl),stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'alloc. error!(SETKRY)work_lapack:ierr=',ierr
            stop
          endif
          work_lapack(:,:)=0.0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ AFTER-LANCZOS : iKnstr=1
c
          if (iKnstr .eq. 1) then
c            ierr=0
c            if (njjt .gt. njjt_mx2) ierr=1
c            if (nreclc .gt. nreclc_mx2) ierr=1
c            if (ierr .eq. 1) then
c               write(6,*)'ERROR!(SETKRY_B):after Lanczos'
c               write(6,*)'   njjt,njjt_mx2=',njjt,njjt_mx2
c               write(6,*)'   nreclc,nreclc_mx2=', nreclc,nreclc_mx2
c               stop
c            endif   
             if (imode .eq. 1) then
               do irec=1,nrecl
                 do jjt=1,njjt
                   rKnstr(jjt,irec,ja4seed,jsv4seed)=rKntmp(jjt,irec)
                 enddo
               enddo  
             endif  
             if (imode .eq. 2) then
               do irec=1,nrecl
                 ddsum=0.0d0 
                 do jjt=1,njjt
                   ddd=rKnstr(jjt,irec,ja4seed,jsv4seed)
                   rKntmp(jjt,irec)=ddd
                   ddsum=ddsum+dabs(ddd)
                 enddo
                 if (ddsum .le. 1.0d-10) then
                    write(6,*)'ERROR!(SETKRY_B):after Lanczos'
                    write(6,*)'ddsum in rKnstr=',ddsum
                    write(6,*)'jrun,irec=',jrun,irec
                    stop
                 endif   
               enddo  
             endif  
          endif  
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ AFTER-LANCZOS : Density matrix (IMODE=2)
c
          if (imode .eq. 2) then
            jsk2=1
            jsv2=jsv4jsk(jsk2)
            if (jsv2 .ne. jsv4seed) then
               write(6,*)'ERROR!(SETKRY):jsv2,jsv4seed=',jsv2,jsv4seed
               stop
            endif   
            if(intf(jsv2).eq.1) then
              njsd2=njsd(jsv2,ict4l)
            else
              njsd2=njsd(jsv2,ict4h)
            endif
            do jsd1=1,njsd2
              jsv1=jsv4jsd(jsd1,jsv2)
              jsk1=jsk4jsv(jsv1)
              js1=js4jsv(jsv1)
              jjkset1=jjkset(jsk1)
              nss1=jsei(js1)
              nval1=nval(nss1)
              do ja1=1,nval1
                 jjk=jjkset1+ja1
                 jjt=jjt4jjk(jjk)
                 if ((jjt .le. 0) .or. (jjt .gt. njjt)) then
                   write(6,*)'ERROR(SETKRY):AF:jjt=',jjt
                   stop
                 endif
                 jjkd=jjk4jjt(jjt)
                 if (jjk .ne. jjkd) then
                   write(6,*)'ERROR(SETKRY):AF:jjk=',jjk,jjkd
                   stop
                 endif
                 ddsum=0.0d0
                 do irec=1,nrecl
                   ddd=rKntmp(jjt,irec)*rRho(irec,ja4seed,jsv4seed)
                   ddsum=ddsum+ddd 
                 enddo   
c                if (jsv4seed .le. 2) then
c                  write(6,8003)ja1,ja4seed,jsv1,jsv4seed,ddsum
c                endif   
                 dbij(ja1,ja4seed,jsd1,jsv4seed)=ddsum
              enddo   
            enddo   
          endif
c
 8003     format('dbij=',4I6,F13.5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ AFTER-LANCZOS : Diagonalize the tridiagonal matrix (IMODE=1)
c          using the LAPACK
c             NOTE: rAtmp is modied in DETEV !!
c
          if (imode .eq. 1) then
c
           rA(:,ja4seed,jsv4seed)=10.0d0
           rB(:,ja4seed,jsv4seed)=0.0d0
c              High dummy values for rA is essential !
c
           do irec=1,nrecl
             rA(irec,ja4seed,jsv4seed)=rAtmp(irec)
             rB(irec,ja4seed,jsv4seed)=rBtmp(irec)
           enddo   
c
c
           call DSTEV('V',nrecl,rAtmp,rBtmp,rGe,nrecl,work_lapack,ierr)
           if(ierr.ne.0) then
             write(6,*) 'in DSTEV, INFO=',ierr
             write(6,*) 'jsv4seed=',jsv4seed
             stop
           endif
c
           rW(:,ja4seed,jsv4seed)=0.0d0
           rE(:,ja4seed,jsv4seed)=0.0d0
           rR(:,ja4seed,jsv4seed)=0.0d0
c
           ddsum=0.0d0
           do irec=1,nrecl
             ddd=rGe(1,irec)**2
             rW(irec,ja4seed,jsv4seed)=rGe(1,irec)
c                   c(1,al) = < K_1 |  w_al >
             rE(irec,ja4seed,jsv4seed)=rAtmp(irec)
             rR(irec,ja4seed,jsv4seed)=ddd
             ddsum=ddsum+ddd
c            if (jsv4seed .eq. 770) then
c              write(6,*)jsv4seed,irec,rAtmp(irec)*ev4au,ddd
c            endif
           enddo   
c          if (jsv4seed .eq. 771) stop
c          write(6,*) 'ddsum=',ddsum
           if (dabs(ddsum-1.0d0) .gt. 1.0d-10) then
             write(6,*) 'ERROR!(SETKRY:After Lapack)ddsum=',ddsum
             stop
           endif   
c
         endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    @@ AFTER-LANCZOS : RN calc
c        rGe(irec,ial) : c_{n al}  eigen vectors
c      rAtmp(irec)     : e_{al}    eigen values
c
          if (iRN_calc .eq. 1) then
c
c
            if (njjt*nkene .eq. 0) then
              write(6,*)'ERROR!(SETKRY)RN_calc'
              write(6,*)'njjt,nkene=',njjt,nkene
              stop
            endif
c   
            do kene=1,nkene
               ddd=(ddemax-ddemin)/dble(nkene-1)
               cZ_k(kene)=ddemin+dble(kene-1)*ddd
     +                     +dcmplx(0.0d0,ene_gamma)
c                   ----> energy mesh z_k = e_k + i gamma
c              write(6,*)'z_k[eV]=',kene,cZ_k(kene)*ev4au
            enddo
c           stop
c   
c           do irec=1,nrecl
c              ddd=rAtmp(irec)
c              write(6,*)'al, E_al[eV]=',irec,ddd*ev4au
c              if ((ddd .le. ddemin) .or. (ddd .gt. ddemax)) then
c               write(6,*)'note :out of range : e_al [eV]=',ddd*ev4au
c              endif   
c           enddo  
c           stop
c
            cGn1tmp(:,:)=dcmplx(0.0d0,0.0d0)
            do irec=1,nrecl
              do kene=1,nkene
                cddsum=0.0d0
                do ial=1,nrecl
                   cdd1=cZ_k(kene)-rAtmp(ial)
c                    -->    ( z_k - e_al )
                   cdd2=rGe(irec,ial)*rGe(1,ial)
c                    -->    c_{n,ial} x c_{1,ial}
                   cddsum=cddsum+cdd2/cdd1
                enddo
                cGn1tmp(irec,kene)=cddsum
c                 -->    Green fn < K_n | G(z) | K_ 1 >
              enddo  
            enddo   
c           write(6,*)'stop manually'
c           stop
c
            do kene=1,nkene
              cRNtmp(:,kene)=dcmplx(rKntmp(:,1),0.0d0)
            enddo 
c
            do jjt=1,njjt
              do kene=1,nkene
                cddsum=0.0d0
                do irec=1,nrecl
                  cdd1=cZ_k(kene)*rKntmp(jjt,irec)
     +                          - rLntmp(jjt,irec)
                  cdd2=cGn1tmp(irec,kene)
                  cddsum=cddsum+cdd1*cdd2
c                 if (cdabs(cdd1) .le. 1.0d-10) then
c                    write(6,*)'cdd1=',irec,kene,cdd1
c                 endif   
c                 if (cdabs(cdd2) .le. 1.0d-10) then
c                    write(6,*)'cdd2=',irec,kene,cdd2
c                 endif   
c                 write(6,*)'jjt,irec,kene=',jjt,irec,kene
c                 write(6,*)'cddsum       =',cddsum
c                 write(6,*)'cdabs(cddsum)=',cdabs(cddsum)
c                 if (cdabs(cddsum) .le. 1.0d-10) then
c                    write(6,*)'cddsum=',irec,kene,cddsum
c                 endif   
                enddo
                cRNtmp(jjt,kene)=cRNtmp(jjt,kene)-cddsum
c               write(6,*)'cRNs  =',cRNtmp(jjt,kene)
              enddo 
            enddo   
c            -->  < i | r(z_k) > 
c             = < i | { I -  (z-H) G(z) } | K_1 >
c             = < i | K_1 >- < i | (z-H) G(z) | K_1 >
c             = < i | K_1 >
c               - sum_n < i | (z-H) | K_n > < K_n | G(z) | K_1 >
c
c           stop
c
c           do jjt=1,njjt
c             jjk =jjk4jjt(jjt) 
c             jsk1=jsk4jjk(jjk)
c             jsv1=jsv4jsk(jsk1)
c             ja1 =ja4jjk(jjk)
c             if ((ja1 .eq. ja4seed) .and. (jsv1 .eq. js4seed)) then
c               ddsum=0.0d0 
c               do kene=1,nkene
c                cdd1=cRNtmp(jjt,kene)
c                ddd=cdabs(cdd1)
c                write(6,8101)jjt,kene,dreal(cZ_k(kene))*ev4au,
c    +                dreal(cdd1),dimag(cdd1),ddsum
c                ddsum=ddsum+ddd**2
c               enddo
c             endif  
c           enddo
c  
            ddsum=0.0d0 
            jjt=ja4seed
            do kene=1,nkene
              cdd1=cRNtmp(jjt,kene)
              ddd=cdabs(cdd1)
              ddsum=ddsum+ddd**2
            enddo
            ddsum=ddsum/dble(nkene)
            rRNdia(ja4seed,jsv4seed)=ddsum
c
            ddsum=0.0d0
            do jjt=1,njjt
              do kene=1,nkene
                 cdd1=cRNtmp(jjt,kene)
                 ddd=cdabs(cdd1)
                 ddsum=ddsum+ddd**2
c                if (ddd .gt. 10.0d0) then
c                  write(6,8101)jjt,kene,dreal(cZ_k(kene))*ev4au,
c    +                 dreal(cdd1),dimag(cdd1),ddsum
c                endif  
              enddo
            enddo   
c           write(6,*)'ddsum =',ja4seed,jsv4seed,ddsum
            ddsum=ddsum/dble(nkene)/dble(njjt)
            rRNave(ja4seed,jsv4seed)=ddsum
c
c           write(6,*)'rRNave=',ja4seed,jsv4seed,ddsum
c           stop
c8101       format('cRNtmp =',2I6,4F30.20) 
c
          endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
          deallocate(rAtmp,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)rAtmp:ierr=',ierr
            stop
          endif
c
          deallocate(rBtmp,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)rBtmp:ierr=',ierr
            stop
          endif
c
          deallocate(rGe,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)rGe:ierr=',ierr
            stop
          endif
c
          deallocate(work_lapack,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)work_lapack:ierr=',ierr
            stop
          endif
c
          deallocate(rKntmp,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)rKntmp:ierr=',ierr
            stop
          endif
c
          deallocate(rLntmp,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)rLntmp:ierr=',ierr
            stop
          endif
c
          deallocate(cGn1tmp,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)cGn1tmp:ierr=',ierr
            stop
          endif
c
          deallocate(cRNtmp,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)cRNtmp:ierr=',ierr
            stop
          endif
c
          deallocate(cZ_k,stat=ierr)
          if(ierr .ne. 0) then
            write(6,*)'dealloc. error!(SETKRY)cZ_k:ierr=',ierr
            stop
          endif
c
 2000  continue
c
       deallocate(jsv4jsk)
       deallocate(jsk4jsv)
       deallocate(jjkset)
       deallocate(rU1)
       deallocate(rU2)
       deallocate(acwrk1)
       deallocate(acwrk2)
       deallocate(w0xl)
       deallocate(jjk4jjt)
       deallocate(jjt4jjk)
       deallocate(jsk4jjk)
       deallocate(ja4jjk)
c
 1000 continue
!$omp end do
!$omp end parallel 
c
      call tclock(ptime)
      if (mr .eq. 0) write(6,*)' @@@OMP-TIME-0100 ',ptime-tb
      tb=ptime
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Check the electron number (IMODE=2)  
c
      if (imode .eq. 2) then
c
        DSUM1=0.0D0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(ishow2,ierr)
!$omp& private(jsv2,jsd1,jsv1,js2,nss2,nval2,ja2,ja1)
!$omp& private(dlijd)
!$omp& reduction(+ : dsum1)
!$omp do schedule(static)
        DO 10 JSV2=1,NOAV_KR
           JSD1=1
           JSV1=JSV4JSD(JSD1,JSV2)
           IF (JSV1 .NE. JSV2) THEN
             WRITE(6,*)'ERROR!(ENETB)'
             WRITE(6,*)'JSV1,JSV2=',JSV1,JSV2
             STOP
           ENDIF
           JS2=JS4JSV(JSV2)
           NSS2=JSEI(JS2)
           NVAL2=NVAL(NSS2)
        DO 10 JA2=1,NVAL2
           JA1=JA2
           DLIJD=DBIJ(JA1,JA2,JSD1,JSV2)
           DSUM1=DSUM1+2.0D0*DLIJD
   10   CONTINUE
!$omp end do
!$omp end parallel 
c
        WRITE(6,*)'N (from dbij)=',DSUM1
c
      endif   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      deallocate(jsv4jrun)
c
      write(6,*)'...SETKRY is ended'
 9999 continue
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
