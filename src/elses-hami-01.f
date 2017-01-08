!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-hami-01.f @@@@@
czzz  @@@@@ 2009/05/08 @@@@@
cccc2006cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0724: Prepared: (NT116p33)
c0726: Prepared: elses_set_dhij based on setdhjimp2 on pt40/work04
c0802: A modi.
c    : A modi.: Negative value ( enea(ioa,2) < 0 ) is reported.
c0803: A modi.
c0922: Error report for negative value is modified. (NT117p27)
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0109: Become a member of pt41 without any modification
c    : A modi. 
c0206: Renamed into 'elses-hami-01.f' from 'elses-hami-xu.f'
c    : Supoorted : Silicon system with Kwon's Hamiltonian (i_system=2)
c        ---> 'elses_set_hami_01' ( based on 'elses_set_dhij' )
c0209: c_system is used instead of i_system 
c    : NOTE: c_system="C_Xu" :Carbon
c                C. H. Xu, C. Z. Wang and K. M. Ho
c                 J.Phys: Condens. Matter 4, (1992) 6047-6054
c            c_system="Si_Kwon" :Silicon
c                I. Kwon et al. PRB49, 7242 (1994).
c    : In 'elses_set_hami_01', 
c       error checking is added for too small system ( ax < 2 r_cut)
c0309: Prepared : Silicon system with Kwon's Hamiltonian 
c        ---> 'elses_set_rest_01' (but not checked)
c0513: 'elses_set_rest_01' for Si case is checked 
c             ephi0 should be one for Si case (NT07B-119,p7)
c    : 'elses_set_hami_01' for Si case is checked (NT07B-119,p9)
c0607: Several write commands are commented out with !dm' mark.
c0828: Prepared; elses_set_hami_01_new (NT07C-120,p43)
c        in which r_base is used and dfij is redined.
c      ----> experimental code, not used at the present moment.
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0103T.Hoshi: Partial interation is supported (NT07E-122,p.29)
c0108T.Hoshi: Modi. on write commands (NT07E-112p39,v000a-wrk05)
c0215T.Hoshi: 'elses_qm_force_01' prepared (NT07E-122p68)
c               (NT07E-122p68, v0.00.12a-wrk02)
c               as common force-part routine for C and Si
c0401T.Hoshi: function dvec is replaced by dvec2
c            in which iperiodic is not used. 
c                   (NT08A-123p23, 0.00.14a-wrk05) 
!           : i_pbc_x, i_pbc_y, i_pbc_z are introduced.
!              in subroutine 'elses_set_hami_01' 
!                    'elses_set_rest_01', 'elses_qm_force_01'
!                    'elses_set_hami_01_new' 
c0606T.Hoshi: A bugfix in 'elses_qm_force_01', (NT08B-124,p76)
c               'goto 130' is commented out.
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0508T.Hoshi: Error message is modified (NT09A-130p77)
c              in subroutine 'elses_set_rest_01'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     @@ Construction of  the DHIJ, DFIJ matrix
c
c       i_partial  = 0 for DHJI and DFIJ as full interatcion
c       i_partial /= 0  for DHJI as partial interatcion
c              =1 only for (ss sigma)
c              =2 only for (sp sigma) 
c              =3 only for (pp sigma) 
c              =4 only for (pp pi)
c              =5 for all but (pp pi)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_hami_01(i_partial)
c
      use elses_mod_ctrl,    only : i_verbose
      use elses_mod_sel_sys, only : c_system
      use elses_mod_phys_const,only : ev4au,angst
!     use elses_mod_sim_cell,  only : iperiodic,ax,ay,az
      use elses_mod_sim_cell,  only : ax,ay,az,
     +                                i_pbc_x, i_pbc_y, i_pbc_z
      use elses_mod_tx,        only : tx,ty,tz,jsei
      use elses_mod_js4jsv,    only : js4jsv,jsv4js
      use elses_mod_jsv4jsd,   only : jsv4jsd,njsd
      use elses_mod_noav,      only : noav
      use elses_arr_dhij,      only : dhij
      use elses_arr_dfij,      only : dfij
      use elses_mod_orb1,      only : nvl, nval
      use elses_mod_orb2,      only : js2j,dbx, dby, dbz, idngl
      use elses_mod_multi,     only : ict4h
      use elses_arr_dhij_cohp, only : dhij_cohp
c
      implicit none
      integer imode, ipe, npe, ierr
      real*8  eunphys
      real*8  dnal0,rnn0,rcc,es0,eps,ep0,esp3a,esp3b
      real*8  qc0,qc1,qc2,qc3,r_cut_tail
      integer jsv2,js2,njsd2,nss2,nval2
      integer jsd1,jsv1,js1,nss1,nval1,jg,ig,ja1,ja2,isym
      real*8  ahij4d 
      real*8  dxc,dyc,dzc,dd,drr,rnn
      real*8  dha,dna,rca,rat1,rat2,rat3,rat4
      real*8  fac1,fac2,fac3,fac4,dargexp,dexpon
      real*8  dddx,potij,dphidr
      real*8  dvss0,dvsp0,dvpp0,dvpp1
      real*8  fss0, fsp0, fpp0, fpp1
      real*8  dbx1,dby1,dbz1,dbx2,dby2,dbz2,ad1,ad2
      real*8  dbx1b,dby1b,dbz1b,dbx2b,dby2b,dbz2b
      real*8  app0,app1,aaa,aaa0,aaax,aaay,aaaz
c
c     common /cdfij/ dfij(3,nvl,nvl,noas,noav0)
c        ----> Matrix for force calc.
c
      real*8 dhal(4),dnal(4),rcal(4)
      real*8 dkwondd(4,2)
c       ----> work array for TB parameters
c
      integer i_partial
c
      ierr=1
      imode=0
      if (c_system .eq. "C_Xu") imode=1
      if (c_system .eq. "Si_Kwon") imode=2
c
      if (i_verbose >= 1) then
         WRITE(6,*)'@@ LSES_SETDHIJ:imode=',imode
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Checking allocation status
c       for partial interaction matrix
c
      if (i_partial /= 0) then
        if (allocated(dhij_cohp) .eqv. .false.) then
          write(*,*)'ERROR:DHIJ_COHP is not allocated'
          stop
        endif
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Paremeter for compatibility
c
      if (ict4h .ne. 1) then
        write(6,*)'ERROR!(SET_DHIJ):ict4h=',ict4h
        stop
      endif   
      eunphys=5.0d0
c
      if (i_verbose >= 1) then
        write(6,*)'   ict4h=',ict4h
        write(6,*)' EUNPHYS=',eunphys
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Paremeter for Xu's Hamiltonian (if imode=1)
c
      if (imode .eq. 1) then 
c
        ierr=0
        DNAL0=2.0D0
c           ---> n
        RNN0=1.536329D0
c           ---> r_0
c
c       DHAL(:)=0.0d0
        DHAL(1)=-5.0D0
        DHAL(2)= 4.7D0
        DHAL(3)= 5.5D0
        DHAL(4)=-1.55D0
c           ---> V_{ss sigma} etc.
c
        DNAL(1)=6.5D0
        DNAL(2)=6.5D0
        DNAL(3)=6.5D0
        DNAL(4)=6.5D0
c           ---> n_c (common)
c
        RCAL(1)=2.18D0
        RCAL(2)=2.18D0
        RCAL(3)=2.18D0
        RCAL(4)=2.18D0
c           ---> r_c (common)
c
        RCC=2.6D0/angst
c           ---> r_m (cut-off distance) in a.u.
c
        ES0=-2.99d0/EV4AU
c           ---> E_s
        EPS=(3.71d0+2.99d0)/EV4AU
        EP0=ES0+EPS
c           ---> E_p
        ESP3A=0.25D0*ES0+0.75D0*EP0
        ESP3B=-0.25D0*(EP0-ES0)
c
        qc0= 6.7392620074314d-3
        qc1=-8.1885359517898d-2
        qc2= 0.1932365259144d0
        qc3= 0.3542874332380d0
c          ---> c_0, c_1, c_2, c_3
        r_cut_tail=2.45d0
c          ---> r_1
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Paremeter for Kwon's Hamiltonian (if imode=2)
c
c
      if (imode .eq. 2) then 
c
        ierr=0
        DNAL0=2.0D0
c           ---> n
        RNN0=2.360352d0
c           ---> r_0 [in A]
c
        DHAL(1)=-2.038D0
        DHAL(2)= 1.745D0
        DHAL(3)= 2.75D0
        DHAL(4)=-1.075D0
c           ---> V_{ss sigma} etc.
c
        DNAL(1)=9.5D0
        DNAL(2)=8.5D0
        DNAL(3)=7.5D0
        DNAL(4)=7.5D0
c           ---> n_c (common)
c
        RCAL(1)=3.4D0
        RCAL(2)=3.55D0
        RCAL(3)=3.7D0
        RCAL(4)=3.7D0
c           ---> r_c (common)
c
        RCC=4.16D0/angst
c           ---> r_m (cut-off distance) in a.u.
c
        ES0=-5.25d0/EV4AU
c           ---> E_s
        EPS=6.45/EV4AU
c           ---> E_p - E_s
        EP0=ES0+EPS
c           ---> E_p
        ESP3A=0.25D0*ES0+0.75D0*EP0
        ESP3B=-0.25D0*(EP0-ES0)
c
        qc0=0.0d0
        qc1=0.0d0
        qc2=0.0d0
        qc3=0.0d0
c          ---> c_0, c_1, c_2, c_3
c
        r_cut_tail=4.16d0+1.0d10
c          ---> r_1 (in Angstrom)
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Error checking is added for too small system ( ax < 2 r_cut)
c       
      if (ax .lt. 2.0d0*rcc) then
        write(6,*)'ERROR!(LSES_SET_HAMI_01):ax,rcc=',ax,rcc
        write(6,*)'This version dose not support too small system'
        stop
      endif   
c
      if (ay .lt. 2.0d0*rcc) then
        write(6,*)'ERROR!(LSES_SET_HAMI_01):ay,rcc=',ay,rcc
        write(6,*)'This version dose not support too small system'
        stop
      endif   
c
      if (az .lt. 2.0d0*rcc) then
        write(6,*)'ERROR!(LSES_SET_HAMI_01):az,rcc=',az,rcc
        write(6,*)'This version dose not support too small system'
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ierr .eq. 1) then
        write(6,*)'ERROR!(SET_DHIJ):ierr=',ierr
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(dkwondd,eunphys)
!$omp& private(ahij4d)
!$omp& private(dxc,dyc,dzc,dd,drr)
!$omp& private(rnn,dha,dna,rca,rat1,rat2,rat3,rat4)
!$omp& private(fac1,fac2,fac3,fac4,dargexp,dexpon)
!$omp& private(dvss0,dvsp0,dvpp0,dvpp1)
!$omp& private(fss0,fsp0,fpp0,fpp1)
!$omp& private(dbx1,dby1,dbz1,dbx2,dby2,dbz2)
!$omp& private(ad1,ad2,dbx1b,dby1b,dbz1b,dbx2b,dby2b,dbz2b)
!$omp& private(app0,app1,aaa,aaa0,aaax,aaay,aaaz)
!$omp& private(dddx,potij,dphidr)
!$omp& private(jsv1,jsv2,js1,js2,njsd2,nss1,nss2,nval1,nval2)
!$omp& private(ja1,ja2,jsd1,ig,jg,isym)
c     ipe=omp_get_thread_num()+1
c     npe=omp_get_num_threads()
c     write(6,*)'ipe,npe=',ipe,npe
!$omp  do schedule(static)
      do 1000 jsv2=1,noav
        js2=js4jsv(jsv2)
        njsd2=njsd(jsv2,ict4h)
        nss2=jsei(js2)
        nval2=nval(nss2)
        do 1100 jsd1=1,njsd2
c
          if (i_partial == 0) then
            dhij(:,:,jsd1,jsv2)=0.0d0
            dfij(:,:,:,jsd1,jsv2)=0.0d0
          endif  
c
          if (i_partial /= 0) then
            dhij_cohp(:,:,jsd1,jsv2)=0.0d0
          endif  
c
          jsv1=jsv4jsd(jsd1,jsv2)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
c
          if (js1 .eq. js2) then
             do ja2=1,nval2
               jg=js2j(ja2,js2)
               do ja1=1,nval1
                 ig=js2j(ja1,js1)
                 if (ig .eq. jg) then
                   ahij4d=esp3a
                   if (idngl(ig) .eq. 1) ahij4d=eunphys
                 else  
                   ahij4d=esp3b
                   if (idngl(ig) .eq. 1) ahij4d=0.0d0
                   if (idngl(jg) .eq. 1) ahij4d=0.0d0
                 endif
                 dhij(ja1,ja2,jsd1,jsv2)=ahij4d
               enddo
             enddo
             goto 1100
          endif
c           ---> on-site terms
c
          dxc=tx(js2)-tx(js1)
          dyc=ty(js2)-ty(js1)
          dzc=tz(js2)-tz(js1)

          if (i_pbc_x == 1) dxc = modulo(dxc + 0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_y == 1) dyc = modulo(dyc + 0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_z == 1) dzc = modulo(dzc + 0.5d0, 1.0d0) - 0.5d0

!         if (iperiodic .eq. 1) then
!           dd=1.0d0
!           dxc=dxc+0.5d0*dd
!           dxc=dmod(dxc+dd,dd)-0.5d0*dd
!           dyc=dyc+0.5d0*dd
!           dyc=dmod(dyc+dd,dd)-0.5d0*dd
!           dzc=dzc+0.5d0*dd
!           dzc=dmod(dzc+dd,dd)-0.5d0*dd
!         endif

          dxc=dxc*ax
          dyc=dyc*ay
          dzc=dzc*az
          drr=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
c           ---> Distance | R_1 - R_2 | in a.u.
          if (drr .ge. rcc) goto 1100
          if (drr .le. 1.0d-10) then
             write(6,*)'ERROR!(SETDHIJ):JS1,JS2,R=',js1,js2,drr
             stop
          endif
c
          dxc=dxc/drr
          dyc=dyc/drr
          dzc=dzc/drr
c           ---> Direction vector : R_1 - R_2

          rnn=drr*0.529177d0
c            ---> Distance in [A]
c
          do isym=1,4
            dha=dhal(isym)
            dna=dnal(isym)
            rca=rcal(isym)
c
            rat1=rnn0/rnn
            rat2=rnn/rca
            rat3=rnn0/rca
            rat4=dnal0/rnn
c
            fac1=rat1**dnal0
            fac2=rat2**dna
            fac3=rat3**dna
            fac4=1.d0+dna*fac2
c
            dargexp=dnal0*(-fac2+fac3)
            dexpon=dexp(dargexp)
            dkwondd(isym,1)=dha*fac1*dexpon
            dkwondd(isym,2)=-dha*fac1*dexpon*rat4*fac4
            if (rnn .gt. r_cut_tail) then
              dddx=rnn-r_cut_tail
              potij=qc0+qc1*dddx+qc2*dddx*dddx
     +            +qc3*dddx*dddx*dddx
              dphidr=qc1+2.0d0*qc2*dddx+3.0d0*qc3*dddx*dddx
              dkwondd(isym,1)=dha*potij
              dkwondd(isym,2)=dha*dphidr
            endif   
          enddo
c
          dvss0=dkwondd(1,1)/ev4au
          dvsp0=dkwondd(2,1)/ev4au
          dvpp0=dkwondd(3,1)/ev4au
          dvpp1=dkwondd(4,1)/ev4au
c            ---> Slator-Koster parameters in au
c
          fss0=dkwondd(1,2)/ev4au*angst
          fsp0=dkwondd(2,2)/ev4au*angst
          fpp0=dkwondd(3,2)/ev4au*angst
          fpp1=dkwondd(4,2)/ev4au*angst
c            ---> Derivative of Slator-Koster parameters in au
c
          if (i_partial == 1) then
c           dvss0: unchanged
            dvsp0=0.0d0
            dvpp0=0.0d0
            dvpp1=0.0d0
          endif   
c
          if (i_partial == 2) then
            dvss0=0.0d0
c           dvsp0:unchanged
            dvpp0=0.0d0
            dvpp1=0.0d0
          endif   
c
          if (i_partial == 3) then
            dvss0=0.0d0
            dvsp0=0.0d0
c           dvpp0:unchanged
            dvpp1=0.0d0
          endif   
c
          if (i_partial == 4) then
            dvss0=0.0d0
            dvsp0=0.0d0
            dvpp0=0.0d0
c           dvpp1:unchanged
          endif   
c
          if (i_partial == 5) then
c           dvss0:unchanged
c           dvsp0:unchanged
c           dvpp0:unchanged
            dvpp1=0.0d0
          endif   
c
          do ja2=1,nval2
            jg=js2j(ja2,js2)
          do ja1=1,nval1
            ig=js2j(ja1,js1)
c
            dbx1=dbx(ig)
            dby1=dby(ig)
            dbz1=dbz(ig)
            dbx2=dbx(jg)
            dby2=dby(jg)
            dbz2=dbz(jg)
c
ccc Inner products
c
            ad2=dbx2*dxc+dby2*dyc+dbz2*dzc
c              inner product ( a_2 | d )
            ad1=dbx1*dxc+dby1*dyc+dbz1*dzc
c              inner product  ( a_1 | d )
c
ccc Vector :  a' = a - (ad) d
c
            dbx1b=dbx1-ad1*dxc
            dby1b=dby1-ad1*dyc
            dbz1b=dbz1-ad1*dzc
c
            dbx2b=dbx2-ad2*dxc
            dby2b=dby2-ad2*dyc
            dbz2b=dbz2-ad2*dzc
c
ccc < p_1 | p_2 > parts 
c
            app0=ad1*ad2
c             double inner product : ( a_1 | d ) ( a_2 | d ) 
            app1=dbx1b*dbx2b+dby1b*dby2b+dbz1b*dbz2b
c             inner product : ( a'_1 | a'_2 )
c
            aaa=dvss0+dsqrt(3.0d0)*(ad2-ad1)*dvsp0
     +            +3.0d0*app0*dvpp0+3.0d0*app1*dvpp1
            aaa=0.25d0*aaa
c
c
            aaa0=fss0+dsqrt(3.0d0)*(ad2-ad1)*fsp0
     +            +3.0d0*app0*fpp0+3.0d0*app1*fpp1
            aaa0=0.25d0*aaa0
c
            aaax=aaa0*(-dxc)*2.0d0
            aaay=aaa0*(-dyc)*2.0d0
            aaaz=aaa0*(-dzc)*2.0d0
c              the factor 2 is that for (ij) and (ji)
c
            if (idngl(ig)+idngl(jg) .ne. 0) then
               aaa=0.0d0
               aaax=0.0d0
               aaay=0.0d0
               aaaz=0.0d0
            endif  
c
c
            ahij4d=aaa
            if ( i_partial ==0 ) then
              dhij(ja1,ja2,jsd1,jsv2)=ahij4d
              dfij(1,ja1,ja2,jsd1,jsv2)=aaax
              dfij(2,ja1,ja2,jsd1,jsv2)=aaay
              dfij(3,ja1,ja2,jsd1,jsv2)=aaaz
            else
              dhij_cohp(ja1,ja2,jsd1,jsv2)=ahij4d
            endif  
c
!           if (jsv2 >= 1) then
!              write(6,*)'dhij=',ja1,ja2,jsv1,jsv2,ahij4d
!              write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaax
!              write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaay
!              write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaaz
!            endif  
c
          enddo
         enddo
 1100   continue
 1000 continue
!$omp end do
!$omp end parallel 
c
c
      if  (i_verbose >= 1 ) then
        write(6,*)' ... suscessfully ended (LSES_SET_DHIJ)'
      endif  
c
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Rest part in order-N cost
c     Note that the length unit is Angstrom, not a.u., in this routine.
c 
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_rest_01
      use elses_mod_ctrl,    only : i_verbose
      use elses_mod_sel_sys, only : c_system
      use elses_mod_phys_const, only : angst,ev4au
!     use elses_mod_sim_cell,   only : iperiodic, noa, ax, ay, az 
      use elses_mod_sim_cell,   only : noa, ax, ay, az, 
     +                                 i_pbc_x, i_pbc_y, i_pbc_z
      use elses_mod_tx,         only : tx, ty, tz
      use elses_mod_foi,        only : foi
      use elses_mod_jsv4jsd,    only : jsv4jsd,njsd
      use elses_mod_noav,       only : noav,noao,nncut
      use elses_mod_ene,        only : ecc, enea
      use elses_mod_md_dat,     only : itemd
c
      implicit none
      integer ierr,i_demo_mode,ishow, imode
      integer j,js,jsv,ioa,joa,joad,noa1,njsd2,nbkedd
      real*8  am,amc,rc,r0,r1,pc1,pc2,pc3,pc4
      real*8  qc0,qc1,qc2,qc3,ephi0,rcut
      real*8  axang, ayang, azang
      real*8  rd,rd2,poti,potij,dphidr,potij0,dphidr0,dddx
      real*8  erep,rxi,ryi,rzi
      real*8  rxx,ryy,rzz,dvec2
      real*8  f0,f1
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     Working arrays
c
      integer, allocatable :: l4oij(:,:)
      integer, allocatable :: nbke(:,:)
      real*8,  allocatable :: force(:,:),f1d(:)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      imode=0
      if (c_system .eq. "C_Xu") imode=1
      if (c_system .eq. "Si_Kwon") imode=2
      ierr=1
c
      if (imode .eq. 1) then
        ierr=0
        am = 3.30304d0
c         ----> m   in Table 1
        amc= 8.6655d0
c         ----> m_c in Table 1
        rc = 2.1052d0
c         ----> d_c in Table 1
        r0 = 1.64d0
c         ----> d_0 in Table 1
        r1 = 2.57d0
c         ----> d_1 in Table 1
        pc1= 0.5721151498619d0
        pc2=-1.7896349903996d-3
        pc3= 2.3539221516757d-5
        pc4=-1.24251169551587d-7
c         ----> C_1, C_2, C_3, C_4 for f(x) in Table 2
        qc0= 2.2504290109d-8
        qc1=-1.4408640561d-6
        qc2= 2.1043303374d-5
        qc3= 6.6024390226d-5
c         ----> C_0,C_1, C_2, C_3 for 
c             tail of phi(x) in Table 2
        ephi0=8.18555d0
c         ----> phi_0 in Table 1
        rcut=2.6d0
c         ----> cutoff radious c
      endif   
c
      if (imode .eq. 2) then
        ierr=0
        am = 6.8755d0
c         ----> m   in Table 1
        amc= 13.017d0
c         ----> m_c in Table 1
        rc = 3.66995d0
c         ----> r_c in Table 1
        r0 = 2.360352d0
c         ----> r_0 in Table 1
        pc1= 2.1604385d0
        pc2=-0.1384393d0
        pc3= 5.8398423d-3
        pc4=-8.0263577d-5
c         ----> C_1, C_2, C_3, C_4 for f(x) in Table 2
        qc0= 0.0d0
        qc1= 0.0d0
        qc2= 0.0d0
        qc3= 0.0d0
c         ----> dummy 
        ephi0=1.0d0
c         ----> dummy 
        rcut=4.16d0
c         ----> cutoff radious c
        r1 = rcut*1.00001d0
c         ----> dummy 
      endif   
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      if (ierr .eq. 1) then
         write(6,*) 'ERROR(LSES_SET_REST_01)'
         write(6,*) ' ierr,c_system=',ierr,c_system
         stop
      endif    
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     a.u. -> Angstrom
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      axang=ax*angst
      ayang=ay*angst
      azang=az*angst
c
      if  (i_verbose >= 1 ) then
        WRITE(6,*) '@@ LSES_SET_REST_01:AX=',AXANG,AYANG,AZANG
      endif  
c
      ierr=0
      if (axang .le. 2.0d0*rcut) ierr=1
      if (ayang .le. 2.0d0*rcut) ierr=1
      if (azang .le. 2.0d0*rcut) ierr=1
      if (ierr .eq. 1) then
       write(6,*) 'ERROR!(REPULC):rcut [A]=',rcut
       write(6,*) 'ERROR!(REPULC):ax   [A]=',axang
       write(6,*) 'ERROR!(REPULC):ay   [A]=',ayang
       write(6,*) 'ERROR!(REPULC):az   [A]=',azang
       write(6,*) 'Too small cell size for cutoff radius?!'
       stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Demonstration only
c
c     i_demo_mode=0
c     if (i_demo_mode .eq. 1) then
c       write(6,*) '...demonstration start'
c       do j=0,2000
c         rd=1.3d0+dble(j)/1000.0d0
c         potij=ephi0*(r0/rd)**am
c    +               *dexp(am*(-(rd/rc)**amc+(r0/rc)**amc))
c         dphidr=(1.d0+amc*(rd/rc)**amc)*(am/rd)*potij
c         potij0=potij
c         dphidr0=dphidr
c         if (rd .gt. r1) then
c           dddx=rd-r1 
c           potij=qc0+qc1*dddx+qc2*dddx*dddx
c    +          +qc3*dddx*dddx*dddx
c           dphidr=-(qc1+2.0d0*qc2*dddx+3.0d0*qc3*dddx*dddx)
c         endif   
c         if (rd .gt. rcut) then 
c           potij=0.0d0
c           dphidr=0.0d0
c         endif   
c         write(6,1011)j,rd,potij,potij0,dphidr,dphidr0
c       enddo
c       write(6,*) '...demonstration end'
c       stop
c     endif   
c1011 format('potij=',i10,f10.6,4f30.20)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ishow=5
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Checking the consistency 
c
      if  (nncut .ne. 2) then
        write(6,*)'alloc. error!(REPULC)NNCUT=',nncut
        stop
      endif   
c
      if  (noav .ne. noa) then
        write(6,*)'alloc. error!(REPULC)noa,noav=',noa,noav
        stop
      endif   
c
      if (allocated(njsd) .eqv. .false.) then
        write(6,*)'ERROR!:Not yet allocated: NJSD'
        stop
      endif  
c
      if (allocated(jsv4jsd) .eqv. .false.) then
        write(6,*)'ERROR!:Not yet allocated: JSV4JSD'
        stop
      endif  
c
      if (allocated(enea) .eqv. .false.) then
        write(6,*)'ERROR!:Not yet allocated: ENEA'
        stop
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ Allocation and setting working arrays
c
c
      allocate (l4oij(noao,noa),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(L4OIJ):ierr=',ierr
        stop
      endif
      l4oij(:,:)=0
c
      allocate (force(noa,3),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(force):ierr=',ierr
        stop
      endif
      force(:,:)=0.0d0
c
      allocate (f1d(noa),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(force):ierr=',ierr
        stop
      endif
      f1d(:)=0.0d0
c
      allocate (nbke(noav,0:nncut),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(NBKE):ierr=',ierr
        stop
      endif
      nbke(:,:)=0
c
      do js=1,noa
        jsv=js
        njsd2=njsd(jsv,1)
        nbke(jsv,1)=njsd2
        if ((njsd2 .le. 0) .or. (njsd2 .gt. noao)) then
          write(6,*)'Error!(REPULC)js,njsd2=',js,njsd2
          stop
        endif   
        l4oij(1:njsd2,jsv)=jsv4jsd(1:njsd2,jsv)
      enddo  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     initialization of the repulsive energy E_{rep}
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      erep=0.d0
c
      do 310 ioa = 1, noa
c
         NOA1=NBKE(IOA,1)
c
         rxi = tx(ioa)*axang
         ryi = ty(ioa)*ayang
         rzi = tz(ioa)*azang

C
         IF (IOA .LE. ISHOW) THEN
            NBKEDD=NBKE(IOA,1)
            if  (i_verbose >= 1 ) then
              WRITE(6,*)'IOA,NBKE=',IOA,NBKEDD
            endif  
         ENDIF
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    initialization of the repulsive pair potential sum_{j}phi(r_{ij})
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         poti=0.d0
         do 340 JOAD= 1, NOA1
            FORCE(JOAD,1)=0.D0
            FORCE(JOAD,2)=0.D0
            FORCE(JOAD,3)=0.D0
340      continue
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    calculation of distance between joa and ioa ion
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         DO 320 JOAD=1,NBKE(IOA,1)
           JOA=L4OIJ(JOAD,IOA)
c
           IF (IOA .EQ. JOA) GOTO 320
c
           IF ((JOA .LE. 0) .OR. (JOA .GT. NOA)) THEN
              WRITE(6,*)'ERROR!(ZREPUL2)'
              WRITE(6,*)'JOA,JOAD,IOA,NBKE=',JOA,JOAD,IOA,NBKE(IOA,1)
              STOP
           ENDIF
c
c          RXX=DVEC(TX(IOA),TX(JOA),1.0D0,IERR)*AXANG
c          RYY=DVEC(TY(IOA),TY(JOA),1.0D0,IERR)*AYANG
c          RZZ=DVEC(TZ(IOA),TZ(JOA),1.0D0,IERR)*AZANG
c
           RXX=DVEC2(TX(IOA),TX(JOA),1.0D0,IERR,i_pbc_x)*AXANG
           RYY=DVEC2(TY(IOA),TY(JOA),1.0D0,IERR,i_pbc_y)*AYANG
           RZZ=DVEC2(TZ(IOA),TZ(JOA),1.0D0,IERR,i_pbc_z)*AZANG
c
c
           rd2=rxx*rxx+ryy*ryy+rzz*rzz
           if(rd2.gt.rcut*rcut) goto 320
           rd =dsqrt(rd2)
C
c          if (ioa .le. ishow) then
c            write(*,*) 'rd[A]=',ioa,joa,rd
c          endif
c
           if(rd .lt. 1.0d-10) then
             write(*,*) 'ERROR!(REPULC):ioa,joa= ',ioa,joa
             write(*,*) 'The two atoms are too near?!'
             write(*,*) 'tx(ioa)=',tx(ioa)
             write(*,*) 'ty(ioa)=',ty(ioa)
             write(*,*) 'tz(ioa)=',tz(ioa)
             write(*,*) 'tx(joa)=',tx(joa)
             write(*,*) 'ty(joa)=',ty(joa)
             write(*,*) 'tz(joa)=',tz(joa)
             write(*,*) ' rxx   =',rxx
             write(*,*) ' ryy   =',ryy
             write(*,*) ' rzz   =',rzz
             stop
           endif   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Calculation of the repulsive pair potential(from joa to ioa)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           potij=ephi0*(r0/rd)**am
     +                *dexp(am*(-(rd/rc)**amc+(r0/rc)**amc))
           dphidr=(1.d0+amc*(rd/rc)**amc)*(am/rd)*potij
c
           if (rd .gt. r1) then
             dddx=rd-r1 
             potij=qc0+qc1*dddx+qc2*dddx*dddx
     +           +qc3*dddx*dddx*dddx
             dphidr=-(qc1+2.0d0*qc2*dddx+3.0d0*qc3*dddx*dddx)
           endif   
c
           poti=poti+potij
C
c          write(*,*) 'rd,potij= ',ioa,joa,rd,potij
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Repulsive Force (force_(from joa to ioa))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           force(joad,1)=dphidr*(rxx/rd)
           force(joad,2)=dphidr*(ryy/rd)
           force(joad,3)=dphidr*(rzz/rd)
c
320      continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Kwon's function f(X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         f0=(((pc4*poti+pc3)*poti+pc2)*poti+pc1)*poti
         f1=((4.d0*pc4*poti+3.d0*pc3)*poti+2.d0*pc2)*poti+pc1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Sum of Repulsive Energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         erep=erep+f0
         enea(ioa,2)=f0/ev4au
         f1d(ioa)=f1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Sum of Repulsive Force (sum_(joa) force_(from joa to ioa))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
         do 330 joad=1, noa1
           joa=l4oij(joad,ioa)
           if(joa.eq.ioa) goto 330
           foi(ioa,1)=foi(ioa,1)+f1*force(joad,1)*(angst/ev4au)
           foi(ioa,2)=foi(ioa,2)+f1*force(joad,2)*(angst/ev4au)
           foi(ioa,3)=foi(ioa,3)+f1*force(joad,3)*(angst/ev4au)
330      continue

100      continue
310   continue
c
        if  (i_verbose >= 1 ) then
          write(*,*) 'Erup [eV,a.u.]= ', erep, erep/ev4au
        endif  
        ecc=erep/ev4au
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Error report, if negative values
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  

        do ioa=1,noa
c         if ((ioa .eq. 41) .or. (ioa .eq. 381)) then
c           write(6,*) 'enea:',itemd,ioa,enea(ioa,2)
c         endif   
c
          if (enea(ioa,2) .le. -1.0d-10) then
            write(6,*) 'ERROR!?:Negative repul:',itemd,ioa,enea(ioa,2)
          endif   
        enddo   
c
        if (ecc .le. -1.0d10) then
          write(*,*) 'ERROR!:Negative value of repulsive energy!!'
          write(*,*) 'Several atoms are unphysically near?!'
          stop
        endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     initialization of the repulsive energy E_{rep}
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      do 410 ioa = 1, noa
c
         NOA1=NBKE(IOA,1)
c
         rxi = tx(ioa)*axang
         ryi = ty(ioa)*ayang
         rzi = tz(ioa)*azang

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    initialization of the repulsive pair potential sum_{j}phi(r_{ij})
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         poti=0.d0
         do 440 joad = 1, noa1
            force(joad,1)=0.d0
            force(joad,2)=0.d0
            force(joad,3)=0.d0
440      continue
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    calculation of distance between joa and ioa ion
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         DO 420 JOAD=1,NOA1
           JOA=L4OIJ(JOAD,IOA)
           IF (IOA .EQ. JOA) GOTO 420
c
           IF ((JOA .LE. 0) .OR. (JOA .GT. NOA)) THEN
              WRITE(6,*)'ERROR!(ZREPUL2)'
              WRITE(6,*)'JOA,JOAD,IOA,NBKE=',JOA,JOAD,IOA,NBKE(IOA,1)
           ENDIF
c
c          RXX=DVEC(TX(IOA),TX(JOA),1.0D0,IERR)*AXANG
c          RYY=DVEC(TY(IOA),TY(JOA),1.0D0,IERR)*AYANG
c          RZZ=DVEC(TZ(IOA),TZ(JOA),1.0D0,IERR)*AZANG
c
           RXX=DVEC2(TX(IOA),TX(JOA),1.0D0,IERR,i_pbc_x)*AXANG
           RYY=DVEC2(TY(IOA),TY(JOA),1.0D0,IERR,i_pbc_y)*AYANG
           RZZ=DVEC2(TZ(IOA),TZ(JOA),1.0D0,IERR,i_pbc_z)*AZANG
c
           rd2=rxx*rxx+ryy*ryy+rzz*rzz
           if(rd2.gt.rcut*rcut) goto 420
           rd =dsqrt(rd2)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Calculation of the repulsive pair potential(from joa to ioa)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           potij=ephi0*(r0/rd)**am*
     +                 dexp(am*(-(rd/rc)**amc+(r0/rc)**amc))
           dphidr=(1.d0+amc*(rd/rc)**amc)*(am/rd)*potij
c
           if (rd .gt. r1) then
             dddx=rd-r1 
             potij=qc0+qc1*dddx+qc2*dddx*dddx
     +           +qc3*dddx*dddx*dddx
             dphidr=-(qc1+2.0d0*qc2*dddx+3.0d0*qc3*dddx*dddx)
           endif   
c
           poti=poti+potij
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Repulsive Force (force_(from joa to ioa))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           force(joad,1)=dphidr*(rxx/rd)
           force(joad,2)=dphidr*(ryy/rd)
           force(joad,3)=dphidr*(rzz/rd)
c
420      continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Kwon's function f(X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         f0=(((pc4*poti+pc3)*poti+pc2)*poti+pc1)*poti
         f1=((4.d0*pc4*poti+3.d0*pc3)*poti+2.d0*pc2)*poti+pc1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    Sum of Repulsive Force (sum_(joa) force_(from joa to ioa))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         do 430 joad=1, noa1
           joa=l4oij(joad,ioa)
           if(joa.eq.ioa) goto 430
           foi(ioa,1)=foi(ioa,1)+f1d(joa)*force(joad,1)*(angst/ev4au)
           foi(ioa,2)=foi(ioa,2)+f1d(joa)*force(joad,2)*(angst/ev4au) 
           foi(ioa,3)=foi(ioa,3)+f1d(joa)*force(joad,3)*(angst/ev4au)
430      continue

400      continue
410   continue

c       write(*,*) 'Erup [eV,a.u.]= ', erep, erep/ev4au
c       ecc=erep/ev4au
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      deallocate(l4oij,stat=ierr)
      if( ierr .ne. 0) then 
        write(6,*)'dealloc. error!(L4OIJ):ierr=',ierr
        stop
      endif   
c
      deallocate(force,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'dealloc. error!(force):ierr=',ierr
        stop
      endif
c
      deallocate(nbke,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'dealloc. error!(NBKE):ierr=',ierr
        stop
      endif
c
      deallocate(f1d,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'dealloc. error!(F1D):ierr=',ierr
        stop
      endif
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ The vector (X2-X1) in the periodicity 
c        (periodicity) = DD
c           assumed 0< X1, X2 < DD
c         --> result [-0.5*DD, 0.5*DD]
c
c    Be careful for the boundary position !!!
c        ( | X2-X1 | = DD ) ---> IERR=1
c
      REAL*8 FUNCTION DVEC2(X2,X1,DD,IERR,I_PBC)
!     use elses_mod_sim_cell,   only : iperiodic
c
      implicit none
      real*8 x2,x1,dd
      integer ierr, i_pbc
      real*8 eps
c
      EPS=1.0D-6
c
      IF (DD .LT. EPS) THEN
         WRITE(6,*)'Error!!(DD =< 0)'
         STOP
      ENDIF   
c
      IERR=0
c
      DVEC2=X2-X1+0.5D0*DD
      DVEC2=DMOD(DVEC2+DD,DD)-0.5D0*DD
c
      IF (I_PBC .EQ. 0) THEN
        DVEC2=X2-X1
        GOTO 9999
      ENDIF
c
      IF (DABS(DABS(DVEC2)-DD/2.0D0) .LE. EPS) IERR=1
c
c
 9999 CONTINUE
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_qm_force_01
c   
      use elses_mod_ctrl,     only : i_verbose
      use elses_mod_sel_sys, only : c_system
      use elses_mod_phys_const,only : ev4au,angst
!     use elses_mod_sim_cell,  only : noa, iperiodic,ax,ay,az
      use elses_mod_sim_cell,  only : noa, ax,ay,az,
     +                                i_pbc_x, i_pbc_y, i_pbc_z
      use elses_mod_tx,        only : tx,ty,tz,jsei
      use elses_mod_js4jsv,    only : js4jsv,jsv4js
      use elses_mod_jsv4jsd,   only : jsv4jsd,njsd
      use elses_mod_foi,       only : foi
      use elses_mod_noav,      only : noav
      use elses_arr_dhij,      only : dhij
      use elses_arr_dbij,      only : dbij
      use elses_arr_dfij,      only : dfij
      use elses_mod_orb1,      only : nvl, nval
      use elses_mod_orb2,      only : js2j,j2js, dbx, dby, dbz, idngl
      use elses_mod_multi,     only : ict4h
c
      implicit none
      integer imode, ipe, npe
      real*8  detb, dens2
      integer js,nss,ja
c
      real*8  dnal0,rnn0,rcc,es0,eps,ep0,esp3a,esp3b
      real*8  qc0,qc1,qc2,qc3,r_cut_tail
c
      integer jsd2,jsv2,js2,njsd2,nss2,nval2,js2b
      integer jsd1,jsv1,js1,nss1,njsd1, nval1,jsv1b,js1b
      integer jg,ig,ja1,ja2,isym,i,j
      real*8  ddalij, ahij4dd
      real*8  dvecx,dvecy,dvecz,ddd,dddr,ddd0
      real*8  dxc,dyc,dzc,dd,drr,rnn
      real*8  dha,dna,rca,rat1,rat2,rat3,rat4
      real*8  fac1,fac2,fac3,fac4,dargexp,dexpon
      real*8  dddx, dddy, dddz, potij,dphidr
      real*8  dvss0,dvsp0,dvpp0,dvpp1
      real*8  fss0, fsp0, fpp0, fpp1
      real*8  dbx1,dby1,dbz1,dbx2,dby2,dbz2,ad1,ad2
      real*8  dbx1b,dby1b,dbz1b,dbx2b,dby2b,dbz2b
      real*8  app0,app1,aaa,aaa0,aaax,aaay,aaaz
c
      real*8  dddx1, dddy1, dddz1
      real*8  dbx0, dby0 , dbz0, pinner
      real*8  fdx1, fdy1 , fdz1, fdx2, fdy2 , fdz2
      real*8  acc0, acc3, acc4
      real*8  dddx2, dddy2, dddz2
      real*8  dddx3, dddy3, dddz3
      real*8  dddx4, dddy4, dddz4
c
      real*8 dhal(4),dnal(4),rcal(4)
      real*8 dkwondd(4,2)
c       ----> work array for TB parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      imode=0
      if (c_system .eq. "C_Xu") imode=1
      if (c_system .eq. "Si_Kwon") imode=2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (i_verbose >= 1) then
         WRITE(6,*)'@@@ ELSES_QM_FORCE_01:imode=',imode
      endif  
c
      if (imode == 0) then
        WRITE(6,*)'ERROR:ELSES_QM_FORCE_01'
        WRITE(6,*)'      not supported imode=',imode
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Paremeter for Xu's TB (if imode=1)
c
      if (imode .eq. 1) then 
        DNAL0=2.0D0
c           ---> n
        RNN0=1.536329D0
c           ---> r_0
c
        DHAL(1)=-5.0D0
        DHAL(2)= 4.7D0
        DHAL(3)= 5.5D0
        DHAL(4)=-1.55D0
c           ---> V_{ss sigma} etc.
c
        DNAL(1)=6.5D0
        DNAL(2)=6.5D0
        DNAL(3)=6.5D0
        DNAL(4)=6.5D0
c           ---> n_c (common)
c
        RCAL(1)=2.18D0
        RCAL(2)=2.18D0
        RCAL(3)=2.18D0
        RCAL(4)=2.18D0
c           ---> r_c (common)
c
c       EUNPHYS=ETERMI
c
        RCC=2.6D0/angst
c           ---> r_m (cut-off distance) in a.u.
c
        ES0=-2.99d0/EV4AU
c           ---> E_s
        EPS=(3.71d0+2.99d0)/EV4AU
        EP0=ES0+EPS
c           ---> E_p
        ESP3A=0.25D0*ES0+0.75D0*EP0
        ESP3B=-0.25D0*(EP0-ES0)
c
        qc0= 6.7392620074314d-3
        qc1=-8.1885359517898d-2
        qc2= 0.1932365259144d0
        qc3= 0.3542874332380d0
c          ---> c_0, c_1, c_2, c_3
        r_cut_tail=2.45d0
c          ---> r_1
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Paremeter for Kwon's Hamiltonian (if imode=2)
c
c
      if (imode .eq. 2) then 
c
        DNAL0=2.0D0
c           ---> n
        RNN0=2.360352d0
c           ---> r_0 [in A]
c
        DHAL(1)=-2.038D0
        DHAL(2)= 1.745D0
        DHAL(3)= 2.75D0
        DHAL(4)=-1.075D0
c           ---> V_{ss sigma} etc.
c
        DNAL(1)=9.5D0
        DNAL(2)=8.5D0
        DNAL(3)=7.5D0
        DNAL(4)=7.5D0
c           ---> n_c (common)
c
        RCAL(1)=3.4D0
        RCAL(2)=3.55D0
        RCAL(3)=3.7D0
        RCAL(4)=3.7D0
c           ---> r_c (common)
c
        RCC=4.16D0/angst
c           ---> r_m (cut-off distance) in a.u.
c
        ES0=-5.25d0/EV4AU
c           ---> E_s
        EPS=6.45/EV4AU
c           ---> E_p - E_s
        EP0=ES0+EPS
c           ---> E_p
        ESP3A=0.25D0*ES0+0.75D0*EP0
        ESP3B=-0.25D0*(EP0-ES0)
c
        qc0=0.0d0
        qc1=0.0d0
        qc2=0.0d0
        qc3=0.0d0
c          ---> c_0, c_1, c_2, c_3
c
        r_cut_tail=4.16d0+1.0d10
c          ---> r_1 (in Angstrom)
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      DETB=0.0D0
      DENS2=0.0D0
c
      if (noa .gt. 5) then
        do js=1,5
          if (i_verbose >= 1) then
            WRITE(6,*)'FORCE JS(X)=',JS,FOI(JS,1)
            WRITE(6,*)'FORCE JS(Y)=',JS,FOI(JS,2)
            WRITE(6,*)'FORCE JS(Z)=',JS,FOI(JS,3)
          endif   
        enddo
      endif  
c
c     stop
c
      ddd0=0.0d0
!$omp  parallel 
!$omp& default(shared)
!$omp& private(ipe,npe)
!$omp& private(dkwondd,ahij4dd)
!$omp& private(dvecx,dvecy,dvecz,dd)
!$omp& private(dxc,dyc,dzc,ddd,drr,rnn)
!$omp& private(dha,dna,rca,rat1,rat2,rat3,rat4,fac1,fac2,fac3,fac4)
!$omp& private(dargexp,dexpon)
!$omp& private(ddalij,dbx1,dby1,dbz1,dbx2,dby2,dbz2)
!$omp& private(dbx0,dby0,dbz0)
!$omp& private(pinner,fdx1,fdy1,fdz1,fdx2,fdy2,fdz2,acc0)
!$omp& private(dddx,dddy,dddz)
!$omp& private(dddx1,dddy1,dddz1)
!$omp& private(dddx2,dddy2,dddz2)
!$omp& private(dddx3,dddy3,dddz3)
!$omp& private(dddx4,dddy4,dddz4)
!$omp& private(ad1,ad2,dbx1b,dby1b,dbz1b,dbx2b,dby2b,dbz2b)
!$omp& private(app0,app1,acc3,acc4)
!$omp& private(dvss0,dvsp0,dvpp0,dvpp1)
!$omp& private(dddr,potij,dphidr)
!$omp& private(js,jsv1,jsv1b,jsv2,js1,js2,js1b,js2b,jsd1,jsd2)
!$omp& private(ja,ja1,ja2,i,j)
!$omp& private(nss,nss1,nss2,nval1,nval2)
!$omp& private(njsd1,njsd2,isym)
!$omp& reduction(+ : ddd0)
      ipe=0
      npe=0
c     ipe=omp_get_thread_num()+1
c     npe=omp_get_num_threads()
c     write(6,*)'ipe,npe=',ipe,npe
!$omp  do schedule(static)
      DO 100 JSV1=1,NOAV
         JS=JS4JSV(JSV1)
c        WRITE(6,*)'JSV1,JS=',JSV1,JS
         IF ((JS .LE. 0) .OR. (JS .GT. NOA)) THEN
           WRITE(6,*)'ERROR!(FORTSTV):JSV1,JS=',JSV1,JS
           STOP
         ENDIF
         JS1=JS
         DDDX=0.0D0
         DDDY=0.0D0
         DDDZ=0.0D0
         NSS=JSEI(JS)
         NVAL1=NVAL(NSS)
c
       NJSD1=NJSD(JSV1,ICT4H)
       DO 110 JSD2=1,NJSD1
          JSV2=JSV4JSD(JSD2,JSV1)
c         WRITE(6,*)'JSV1,JSV2=',JSV1,JSV2
          IF (JSV1 .EQ. JSV2) GOTO 110
          IF ((JSV2 .LE. 0) .OR. (JSV2 .GT. NOAV)) THEN
            WRITE(6,*)'ERROR!(FORTSTV)'
            WRITE(6,*)'JSV2,JSD2,JSV1=',JSV2,JSD2,JSV1
            STOP
          ENDIF
          JS2=JS4JSV(JSV2)
          IF ((JS2 .LE. 0) .OR. (JS2 .GT. NOA)) THEN
            WRITE(6,*)'ERROR!(FORTSTV):JS2,JSV2=',JS2,JSV2
            STOP
          ENDIF
          NSS2=JSEI(JS2)
          NVAL2=NVAL(NSS2)
c
          njsd2=njsd(jsv2,ict4h)
          do 133 jsd1=1,njsd2
           jsv1b=jsv4jsd(jsd1,jsv2)
           if (jsv1b .eq. jsv1) goto 134
  133     continue
          write(6,*)'ERROR!(FORSTVMP2):JSV1,JSV2=',JSV1,JSV2
          stop
  134     continue
c            ---> jsd1 is generated.
c
c         write(6,*)'JSD1=',JSD1
c
          dvecx=tx(js2)-tx(js1)
          dvecy=ty(js2)-ty(js1)
          dvecz=tz(js2)-tz(js1)

          if (i_pbc_x == 1) 
     +         dvecx = modulo(dvecx + 0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_y == 1) 
     +         dvecy = modulo(dvecy + 0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_z == 1) 
     +         dvecz = modulo(dvecz + 0.5d0, 1.0d0) - 0.5d0

!         if (iperiodic .eq. 1) then
!           dd=1.0d0
!           dvecx=dvecx+0.5d0*dd
!           dvecx=dmod(dvecx+dd,dd)-0.5d0*dd
!           dvecy=dvecy+0.5d0*dd
!           dvecy=dmod(dvecy+dd,dd)-0.5d0*dd
!           dvecz=dvecz+0.5d0*dd
!           dvecz=dmod(dvecz+dd,dd)-0.5d0*dd
!         endif
c
          dxc=dvecx*ax
          dyc=dvecy*ay
          dzc=dvecz*az
c
          DDD=DSQRT(DXC*DXC+DYC*DYC+DZC*DZC)
          DRR=DDD
c         write(6,*)'DRR=',DRR
          if (ddd .le. 1.0d-10) then
             write(6,*)'ERROR!(FORTSTVMP)'
             write(6,*)'JS1,JS2,DDD=',js1,js2,ddd
             stop
          endif   
          IF (DRR .GT. RCC) GOTO 110
c
c
          DXC=DXC/DDD
          DYC=DYC/DDD
          DZC=DZC/DDD
c
          RNN=DRR*0.529177D0
c            ---> distance in [A]
          do isym=1,4
            DHA=DHAL(ISYM)
            DNA=DNAL(ISYM)
            RCA=RCAL(ISYM)
c
            RAT1=RNN0/RNN
            RAT2=RNN/RCA
            RAT3=RNN0/RCA
            RAT4=DNAL0/RNN
c
            FAC1=RAT1**DNAL0
            FAC2=RAT2**DNA
            FAC3=RAT3**DNA
            FAC4=1.d0+DNA*FAC2
c
            DARGEXP=DNAL0*(-FAC2+FAC3)
            DEXPON=DEXP(DARGEXP)
            dkwondd(isym,1)=dha*fac1*dexpon
            dkwondd(isym,2)=-dha*fac1*dexpon*rat4*fac4
            if (rnn .gt. r_cut_tail) then
              dddr=rnn-r_cut_tail
              potij=qc0+qc1*dddr+qc2*dddr*dddr
     +            +qc3*dddr*dddr*dddr
              dphidr=qc1+2.0d0*qc2*dddr+3.0d0*qc3*dddr*dddr
              dkwondd(isym,1)=dha*potij
              dkwondd(isym,2)=dha*dphidr
            endif   
          enddo
c
c
          DVSS0=DKWONDD(1,1)/EV4AU
          DVSP0=DKWONDD(2,1)/EV4AU
          DVPP0=DKWONDD(3,1)/EV4AU
          DVPP1=DKWONDD(4,1)/EV4AU
c
c         DVSS0=DKWON2(1,DRR,0)/EV4AU
c         DVSP0=DKWON2(2,DRR,0)/EV4AU
c         DVPP0=DKWON2(3,DRR,0)/EV4AU
c         DVPP1=DKWON2(4,DRR,0)/EV4AU
c
c         write(6,*)'DVSS0=',DVSS0
c         write(6,*)'DVSP0=',DVSP0
c         write(6,*)'DVPP0=',DVPP0
c         write(6,*)'DVPP1=',DVPP1
c
c         write(6,*)'stop manually'
c         stop
c
c         WRITE(6,*)'  goto 120'
          DO 120 JA=1,NVAL1
            JA1=JA 
            I=JS2J(JA,JS)
c           WRITE(6,*)' loop 120:JS,JS2=',JS,JS2
            JS1B=J2JS(I)
            IF (JS1 .NE. JS1B) THEN
             WRITE(6,*)'ERROR!(FORTSTV):JS1,JS1B=',JS1,JS1B
             STOP
            ENDIF
c
          DO 130 JA2=1,NVAL2
c            WRITE(6,*)' loop 130:JS,JS2=',JS,JS2
             J=JS2J(JA2,JS2)
             JS2B=J2JS(J)
             IF (JS2 .NE. JS2B) THEN
              WRITE(6,*)'ERROR!(FORTSTV):JS2,JS2B=',JS2,JS2B
              STOP
             ENDIF
c
cc  Force for  V'(r) part 
c
c            WRITE(6,*)' DDD0=',DDD0
c
             DDALIJ=DBIJ(JA2,JA1,JSD2,JSV1)
c            ahij4dd=AHIJ4(I,J,0)
             ahij4dd=dhij(ja2,ja1,jsd2,jsv1)
c              equivalent to ahij(i,j)=ahij(j,i)
             DDD0=DDD0+2.0D0*DDALIJ*AHIJ4DD
c
c            WRITE(6,*)' DBIJ,DHIJ=',DDALIJ,AHIJ4DD
c
c            IF (I .EQ. 1) THEN
c             WRITE(6,*)'I,J=',I,J,DDALIJ,AHIJ4DD
c             WRITE(6,*)'I,J=',I,J,JS1,JS2,JA1,JA2,DDALIJ,AHIJ4(I,J,2)
c             WRITE(6,*)'I,J=',I,J,JS1,JS2,JA1,JA2,DDALIJ,AHIJ4(I,J,3)
c            ENDIF
c
             DBX1=DBX(I)
             DBY1=DBY(I)
             DBZ1=DBZ(I)
             DBX2=DBX(J)
             DBY2=DBY(J)
             DBZ2=DBZ(J)
c
c            JSD1=JSD4JSV(JSV1,JSV2)
c
             DDDX1=2.0D0*DDALIJ*DFIJ(1,JA1,JA2,JSD1,JSV2)
             DDDY1=2.0D0*DDALIJ*DFIJ(2,JA1,JA2,JSD1,JSV2)
             DDDZ1=2.0D0*DDALIJ*DFIJ(3,JA1,JA2,JSD1,JSV2)
c
c!!!             if (jsv2 >= 1) then
c!!!               write(6,*)'dbij=',ja1,ja2,jsv1,jsv2,ddalij
c!!!               write(6,*)'dhij=',ja1,ja2,jsv1,jsv2,ahij4dd
c!!!               write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,
c!!!     +                           dfij(1,ja1,ja2,jsd1,jsv2)
c!!!               write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,
c!!!     +                           dfij(2,ja1,ja2,jsd1,jsv2)
c!!!               write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,
c!!!     +                           dfij(3,ja1,ja2,jsd1,jsv2)
c!!!             endif  
c
c            DDDX1=2.0D0*DDALIJ*AHIJ4(I,J,1)
c            DDDY1=2.0D0*DDALIJ*AHIJ4(I,J,2)
c            DDDZ1=2.0D0*DDALIJ*AHIJ4(I,J,3)
c
c
c            DDDX1=0.0D0
c            DDDY1=0.0D0
c            DDDZ1=0.0D0
c
c            IF (I .EQ. 1) THEN
c             WRITE(6,*)'DDDX1=',I,J,DDDX1
c            ENDIF
c
             DDDX=DDDX+DDDX1
             DDDY=DDDY+DDDY1
             DDDZ=DDDZ+DDDZ1
c            goto 130
c
cc  Force for  V_{sp sigma} part 
c     : including < s^* | H | p > case (JA1=5)
c
             DBX0=DBX1-DBX2
             DBY0=DBY1-DBY2
             DBZ0=DBZ1-DBZ2
c
c            IF (JA1 .EQ. 5) THEN
c              DBX0=DBX1
c              DBY0=DBY1
c              DBZ0=DBZ1
c            ENDIF
c
c            IF (JA2 .EQ. 5) THEN
c              DBX0=DBX2
c              DBY0=DBY2
c              DBZ0=DBZ2
c            ENDIF
c
             PINNER=DXC*DBX0+DYC*DBY0+DZC*DBZ0
c
             FDX1=0.0D0
             FDY1=0.0D0
             FDZ1=0.0D0
             FDX2=0.0D0
             FDY2=0.0D0
             FDZ2=0.0D0
c
             FDX1=DBX0/DDD
             FDY1=DBY0/DDD
             FDZ1=DBZ0/DDD
c
             FDX2=PINNER*(-DXC)/DDD
             FDY2=PINNER*(-DYC)/DDD
             FDZ2=PINNER*(-DZC)/DDD
c
             ACC0=2.0D0*DDALIJ*DSQRT(3.0D0)/4.0D0*DVSP0
c               * factor 2.0 is the para-spin factor
c
c            IF (JA1 .EQ. 5) THEN
c              DVSP02=DRSPS*DKWON2(2,DRR,0)/EV4AU
c              ACC0=2.0D0*DDALIJ*DSQRT(3.0D0)/2.0D0*DVSP02
c            ENDIF
c
c            IF (JA2 .EQ. 5) THEN
c              DVSP02=-DRSPS*DKWON2(2,DRR,0)/EV4AU
c              ACC0=2.0D0*DDALIJ*DSQRT(3.0D0)/2.0D0*DVSP02
c            ENDIF
c
             DDDX2=ACC0*2.0D0*(FDX1+FDX2)
             DDDY2=ACC0*2.0D0*(FDY1+FDY2)
             DDDZ2=ACC0*2.0D0*(FDZ1+FDZ2)
c               * factor 2.0 is that for (ij) and (ji) 
c
             DDDX=DDDX+DDDX2
             DDDY=DDDY+DDDY2
             DDDZ=DDDZ+DDDZ2
c
c            WRITE(6,*)' ACC0=',ACC0
c            WRITE(6,*)'DDDX2=',DDDX2
c            WRITE(6,*)'DDDY2=',DDDY2
c            WRITE(6,*)'DDDZ2=',DDDZ2
c
c            IF (I .EQ. 2) THEN
c             WRITE(6,*)'DDDX2=',I,J,DDDX2
c            ENDIF
c
c            IF (JA1 .EQ. 5) GOTO 130
c            IF (JA2 .EQ. 5) GOTO 130
c
cc  Force for  V_{pp sigma} part
c
             AD2=DBX2*DXC+DBY2*DYC+DBZ2*DZC
c              inner product ( a_2 | r^ )
             AD1=DBX1*DXC+DBY1*DYC+DBZ1*DZC
c              inner product  ( a_1 | r^ )
c
c
             DBX1B=DBX1-AD1*DXC
             DBY1B=DBY1-AD1*DYC
             DBZ1B=DBZ1-AD1*DZC
c              -->  Vector :  a1' = a1 - (a1|r^) r^
             DBX2B=DBX2-AD2*DXC
             DBY2B=DBY2-AD2*DYC
             DBZ2B=DBZ2-AD2*DZC
c              -->  Vector :  a2' = a2 - (a2|r^) r^
c
             APP0=AD1*AD2
c              -->  double inner product 
c                   : ( a_1 | r^ ) ( a_2 | r^ ) 
             APP1=DBX1B*DBX2B+DBY1B*DBY2B+DBZ1B*DBZ2B
c              -->  inner product : ( a'_1 | a'_2 )
c
             ACC3=2.0D0*DDALIJ*3.0D0/4.0D0*DVPP0
c               * factor 2.0 is the para-spin factor
             DDDX3=-ACC3*2.0D0*(AD2*DBX1B+AD1*DBX2B)/DDD
             DDDY3=-ACC3*2.0D0*(AD2*DBY1B+AD1*DBY2B)/DDD
             DDDZ3=-ACC3*2.0D0*(AD2*DBZ1B+AD1*DBZ2B)/DDD
c               * factor 2.0 is that for (ij) and (ji) 
c
c            WRITE(6,*)' APP0=',APP0
c            WRITE(6,*)' APP1=',APP1
c            WRITE(6,*)' ACC3=',ACC3
c            WRITE(6,*)'DDDX3=',DDDX3
c            WRITE(6,*)'DDDY3=',DDDY3
c            WRITE(6,*)'DDDZ3=',DDDZ3
c
c
c            DDDX3=0.0D0
c            DDDY3=0.0D0
c            DDDZ3=0.0D0
c
             DDDX=DDDX+DDDX3
             DDDY=DDDY+DDDY3
             DDDZ=DDDZ+DDDZ3
c
c            IF (I .EQ. 1) THEN
c             W=DRR
c             WRITE(6,*)'ACC3,VPP0=',ACC3,DVPP0,DKWON2(3,W,0)/EV4AU
c             WRITE(6,*)'AD1=',AD1,AD2
c             WRITE(6,*)'DDD=',DDD
c             WRITE(6,*)'DDDX3=',I,J,DDDX3
c             WRITE(6,*)'DDDY3=',I,J,DDDY3
c             WRITE(6,*)'DDDZ3=',I,J,DDDZ3
c            ENDIF
c
             DDDX4=(-AD1)*DBX2B+(-AD2)*DBX1B
             DDDY4=(-AD1)*DBY2B+(-AD2)*DBY1B
             DDDZ4=(-AD1)*DBZ2B+(-AD2)*DBZ1B
c
             DDDX4=-DDDX4/DDD
             DDDY4=-DDDY4/DDD
             DDDZ4=-DDDZ4/DDD
c
             ACC4=2.0D0*DDALIJ*3.0D0/4.0D0*DVPP1
c               * factor 2.0 is the para-spin factor
c
c            DDDX4=0.0D0
c            DDDY4=0.0D0
c            DDDZ4=0.0D0
c
             DDDX=DDDX+DDDX4*ACC4*2.0D0
             DDDY=DDDY+DDDY4*ACC4*2.0D0
             DDDZ=DDDZ+DDDZ4*ACC4*2.0D0
c               * factor 2.0 is that for (ij) and (ji) 
c
c            WRITE(6,*)' ACC4=',ACC4
c            WRITE(6,*)'DDDX4=',DDDX4
c            WRITE(6,*)'DDDY4=',DDDY4
c            WRITE(6,*)'DDDZ4=',DDDZ4
c
c            IF (I .EQ. 1) THEN
c             WRITE(6,*)'DDDX4=',I,J,DDDX4
c            ENDIF
c
c
  130     CONTINUE
  120     CONTINUE
c
c
  110  CONTINUE
c
       FOI(JS,1)=FOI(JS,1)-DDDX
       FOI(JS,2)=FOI(JS,2)-DDDY
       FOI(JS,3)=FOI(JS,3)-DDDZ
c
  100 CONTINUE
!$omp end do
!$omp end parallel 
c
  199 CONTINUE
c
      WRITE(6,*)'ETB(off site only)=',DDD0
      DO 200 JS=1,5
         if (i_verbose >= 1) then
           WRITE(6,*)'FORCE JS(X)=',JS,FOI(JS,1)
           WRITE(6,*)'FORCE JS(Y)=',JS,FOI(JS,2)
           WRITE(6,*)'FORCE JS(Z)=',JS,FOI(JS,3)
         endif  
  200 CONTINUE
c
c     IF (IREPONLY .EQ. 1) THEN
c       WRITE(6,*)'REPULSIVE-only!!!'
c     ENDIF
c
c     DO 300 III=1,3
c     DO 300 JS=1,NOA
c        FOI(JS,III)=FEC(JS,III)+FC(JS,III)
c        IF (IREPONLY .EQ. 1) THEN
c          FOI(JS,III)=FC(JS,III)
c        ENDIF
c 300 CONTINUE
c
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     @@ Construction of  the DHIJ, DFIJ matrix
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_hami_01_new
c
      use elses_mod_sel_sys, only : c_system
      use elses_mod_phys_const,only : ev4au,angst
!     use elses_mod_sim_cell,  only : iperiodic,ax,ay,az
      use elses_mod_sim_cell,  only : ax,ay,az,
     +                                i_pbc_x, i_pbc_y, i_pbc_z
      use elses_mod_tx,        only : tx,ty,tz,jsei
      use elses_mod_js4jsv,    only : js4jsv,jsv4js
      use elses_mod_jsv4jsd,   only : jsv4jsd,njsd
      use elses_mod_noav,      only : noav
      use elses_arr_dhij,      only : dhij
      use elses_arr_dfij,      only : dfij
      use elses_mod_orb1,      only : nvl, nval
      use elses_mod_orb2,      only : js2j,dbx, dby, dbz, idngl
      use elses_mod_multi,     only : ict4h
      use elses_mod_r_base,    only : r_base
c
      implicit none
      integer imode, ipe, npe, ierr
      real*8  eunphys
      real*8  dnal0,rnn0,rcc,es0,eps,ep0,esp3a,esp3b
      real*8  qc0,qc1,qc2,qc3,r_cut_tail
      integer jsv2,js2,njsd2,nss2,nval2
      integer jsd1,jsv1,js1,nss1,nval1,jg,ig,ja1,ja2,isym
      real*8  ahij4d 
      real*8  dxc,dyc,dzc,dd,drr,rnn
      real*8  dha,dna,rca,rat1,rat2,rat3,rat4
      real*8  fac1,fac2,fac3,fac4,dargexp,dexpon
      real*8  dddx,potij,dphidr
      real*8  dvss0,dvsp0,dvpp0,dvpp1
      real*8  fss0, fsp0, fpp0, fpp1
      real*8  dbx1,dby1,dbz1,dbx2,dby2,dbz2,ad1,ad2
      real*8  dbx1b,dby1b,dbz1b,dbx2b,dby2b,dbz2b
      real*8  app0,app1,aaa,aaa0,aaax,aaay,aaaz
      real*8  bbbi, bbbj, bbbx, bbby, bbbz
      real*8  dvspi, dvspj, fspi, fspj
      real*8  dbs1,dbs2
c
c     common /cdfij/ dfij(3,nvl,nvl,noas,noav0)
c        ----> Matrix for force calc.
c
      real*8 dhal(4),dnal(4),rcal(4)
      real*8 dkwondd(4,2)
c       ----> work array for TB parameters
c
      ierr=1
      imode=0
      if (c_system .eq. "C_Xu") imode=1
      if (c_system .eq. "Si_Kwon") imode=2
c
      WRITE(6,*)'@@ LSES_SET_HAMI_01_NEW:imode=',imode
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Paremeter for compatibility
c
      if (ict4h .ne. 1) then
        write(6,*)'ERROR!(SET_DHIJ):ict4h=',ict4h
        stop
      endif   
      eunphys=5.0d0
c
      write(6,*)'   ict4h=',ict4h
      write(6,*)' EUNPHYS=',eunphys
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Paremeter for Xu's Hamiltonian (if imode=1)
c
      if (imode .eq. 1) then 
c
        ierr=0
        DNAL0=2.0D0
c           ---> n
        RNN0=1.536329D0
c           ---> r_0
c
c       DHAL(:)=0.0d0
        DHAL(1)=-5.0D0
        DHAL(2)= 4.7D0
        DHAL(3)= 5.5D0
        DHAL(4)=-1.55D0
c           ---> V_{ss sigma} etc.
c
        DNAL(1)=6.5D0
        DNAL(2)=6.5D0
        DNAL(3)=6.5D0
        DNAL(4)=6.5D0
c           ---> n_c (common)
c
        RCAL(1)=2.18D0
        RCAL(2)=2.18D0
        RCAL(3)=2.18D0
        RCAL(4)=2.18D0
c           ---> r_c (common)
c
        RCC=2.6D0/angst
c           ---> r_m (cut-off distance) in a.u.
c
        ES0=-2.99d0/EV4AU
c           ---> E_s
        EPS=(3.71d0+2.99d0)/EV4AU
        EP0=ES0+EPS
c           ---> E_p
        ESP3A=0.25D0*ES0+0.75D0*EP0
        ESP3B=-0.25D0*(EP0-ES0)
c
        qc0= 6.7392620074314d-3
        qc1=-8.1885359517898d-2
        qc2= 0.1932365259144d0
        qc3= 0.3542874332380d0
c          ---> c_0, c_1, c_2, c_3
        r_cut_tail=2.45d0
c          ---> r_1
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Paremeter for Kwon's Hamiltonian (if imode=2)
c
c
      if (imode .eq. 2) then 
c
        ierr=0
        DNAL0=2.0D0
c           ---> n
        RNN0=2.360352d0
c           ---> r_0 [in A]
c
        DHAL(1)=-2.038D0
        DHAL(2)= 1.745D0
        DHAL(3)= 2.75D0
        DHAL(4)=-1.075D0
c           ---> V_{ss sigma} etc.
c
        DNAL(1)=9.5D0
        DNAL(2)=8.5D0
        DNAL(3)=7.5D0
        DNAL(4)=7.5D0
c           ---> n_c (common)
c
        RCAL(1)=3.4D0
        RCAL(2)=3.55D0
        RCAL(3)=3.7D0
        RCAL(4)=3.7D0
c           ---> r_c (common)
c
        RCC=4.16D0/angst
c           ---> r_m (cut-off distance) in a.u.
c
        ES0=-5.25d0/EV4AU
c           ---> E_s
        EPS=6.45/EV4AU
c           ---> E_p - E_s
        EP0=ES0+EPS
c           ---> E_p
        ESP3A=0.25D0*ES0+0.75D0*EP0
        ESP3B=-0.25D0*(EP0-ES0)
c
        qc0=0.0d0
        qc1=0.0d0
        qc2=0.0d0
        qc3=0.0d0
c          ---> c_0, c_1, c_2, c_3
c
        r_cut_tail=4.16d0+1.0d10
c          ---> r_1 (in Angstrom)
      endif  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Error checking is added for too small system ( ax < 2 r_cut)
c       
      if (ax .lt. 2.0d0*rcc) then
        write(6,*)'ERROR!(LSES_SET_HAMI_01):ax,rcc=',ax,rcc
        write(6,*)'This version dose not support too small system'
        stop
      endif   
c
      if (ay .lt. 2.0d0*rcc) then
        write(6,*)'ERROR!(LSES_SET_HAMI_01):ay,rcc=',ay,rcc
        write(6,*)'This version dose not support too small system'
        stop
      endif   
c
      if (az .lt. 2.0d0*rcc) then
        write(6,*)'ERROR!(LSES_SET_HAMI_01):az,rcc=',az,rcc
        write(6,*)'This version dose not support too small system'
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (ierr .eq. 1) then
        write(6,*)'ERROR!(SET_DHIJ):ierr=',ierr
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 1000 jsv2=1,noav
        js2=js4jsv(jsv2)
        njsd2=njsd(jsv2,ict4h)
        nss2=jsei(js2)
        nval2=nval(nss2)
        do 1100 jsd1=1,njsd2
          dhij(:,:,jsd1,jsv2)=0.0d0
          dfij(:,:,:,jsd1,jsv2)=0.0d0
c
          jsv1=jsv4jsd(jsd1,jsv2)
          js1=js4jsv(jsv1)
          nss1=jsei(js1)
          nval1=nval(nss1)
c
          if (js1 .eq. js2) then
             do ja2=1,nval2
               jg=js2j(ja2,js2)
               do ja1=1,nval1
                 ig=js2j(ja1,js1)
                 if (ig .eq. jg) then
                   ahij4d=esp3a
                   if (idngl(ig) .eq. 1) ahij4d=eunphys
                 else  
                   ahij4d=esp3b
                   if (idngl(ig) .eq. 1) ahij4d=0.0d0
                   if (idngl(jg) .eq. 1) ahij4d=0.0d0
                 endif
                 dhij(ja1,ja2,jsd1,jsv2)=ahij4d
               enddo
             enddo
             goto 1100
          endif
c           ---> on-site terms
c
          dxc=tx(js2)-tx(js1)
          dyc=ty(js2)-ty(js1)
          dzc=tz(js2)-tz(js1)

          if (i_pbc_x == 1) dxc = modulo(dxc + 0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_y == 1) dyc = modulo(dyc + 0.5d0, 1.0d0) - 0.5d0
          if (i_pbc_z == 1) dzc = modulo(dzc + 0.5d0, 1.0d0) - 0.5d0

!         if (iperiodic .eq. 1) then
!           dd=1.0d0
!           dxc=dxc+0.5d0*dd
!           dxc=dmod(dxc+dd,dd)-0.5d0*dd
!           dyc=dyc+0.5d0*dd
!           dyc=dmod(dyc+dd,dd)-0.5d0*dd
!           dzc=dzc+0.5d0*dd
!           dzc=dmod(dzc+dd,dd)-0.5d0*dd
!         endif

          dxc=dxc*ax
          dyc=dyc*ay
          dzc=dzc*az
          drr=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
c           ---> Distance | R_1 - R_2 | in a.u.
          if (drr .ge. rcc) goto 1100
          if (drr .le. 1.0d-10) then
             write(6,*)'ERROR!(SETDHIJ):JS1,JS2,R=',js1,js2,drr
             stop
          endif
c
          dxc=dxc/drr
          dyc=dyc/drr
          dzc=dzc/drr
c           ---> Direction vector : R_1 - R_2

          rnn=drr*0.529177d0
c            ---> Distance in [A]
c
          do isym=1,4
            dha=dhal(isym)
            dna=dnal(isym)
            rca=rcal(isym)
c
            rat1=rnn0/rnn
            rat2=rnn/rca
            rat3=rnn0/rca
            rat4=dnal0/rnn
c
            fac1=rat1**dnal0
            fac2=rat2**dna
            fac3=rat3**dna
            fac4=1.d0+dna*fac2
c
            dargexp=dnal0*(-fac2+fac3)
            dexpon=dexp(dargexp)
            dkwondd(isym,1)=dha*fac1*dexpon
            dkwondd(isym,2)=-dha*fac1*dexpon*rat4*fac4
            if (rnn .gt. r_cut_tail) then
              dddx=rnn-r_cut_tail
              potij=qc0+qc1*dddx+qc2*dddx*dddx
     +            +qc3*dddx*dddx*dddx
              dphidr=qc1+2.0d0*qc2*dddx+3.0d0*qc3*dddx*dddx
              dkwondd(isym,1)=dha*potij
              dkwondd(isym,2)=dha*dphidr
            endif   
          enddo
c
          dvss0=dkwondd(1,1)/ev4au
          dvsp0=dkwondd(2,1)/ev4au
          dvpp0=dkwondd(3,1)/ev4au
          dvpp1=dkwondd(4,1)/ev4au
          dvspi=dvsp0
          dvspj=dvsp0
c            ---> Slator-Koster parameters in au
c
          fss0=dkwondd(1,2)/ev4au*angst
          fsp0=dkwondd(2,2)/ev4au*angst
          fpp0=dkwondd(3,2)/ev4au*angst
          fpp1=dkwondd(4,2)/ev4au*angst
          fspi=fsp0
          fspj=fsp0
c            ---> Derivative of Slator-Koster parameters in au
c
c
          do ja2=1,nval2
            jg=js2j(ja2,js2)
          do ja1=1,nval1
            ig=js2j(ja1,js1)
c
c           dbx1=dbx(ig)
c           dby1=dby(ig)
c           dbz1=dbz(ig)
c           dbx2=dbx(jg)
c           dby2=dby(jg)
c           dbz2=dbz(jg)
c
            dbs1=r_base(1,ig)
            dbx1=r_base(2,ig)
            dby1=r_base(3,ig)
            dbz1=r_base(4,ig)
!
            dbs2=r_base(1,jg)
            dbx2=r_base(2,jg)
            dby2=r_base(3,jg)
            dbz2=r_base(4,jg)

c
ccc Inner products
c
            ad2=dbx2*dxc+dby2*dyc+dbz2*dzc
c              inner product ( a_2 | d )
            ad1=dbx1*dxc+dby1*dyc+dbz1*dzc
c              inner product  ( a_1 | d )
c
ccc Vector :  a' = a - (ad) d
c
            dbx1b=dbx1-ad1*dxc
            dby1b=dby1-ad1*dyc
            dbz1b=dbz1-ad1*dzc
c
            dbx2b=dbx2-ad2*dxc
            dby2b=dby2-ad2*dyc
            dbz2b=dbz2-ad2*dzc
c
ccc < p_1 | p_2 > parts 
c
            app0=ad1*ad2
c             double inner product : ( a_1 | d ) ( a_2 | d ) 
            app1=dbx1b*dbx2b+dby1b*dby2b+dbz1b*dbz2b
c             inner product : ( a'_1 | a'_2 )
c
            aaa= dbs1*dbs2*dvss0
     +          +dbs1*ad2*dvspi
     +          -dbs2*ad1*dvspj
     +          +app0*dvpp0+app1*dvpp1
c
c           aaa=dvss0+dsqrt(3.0d0)*(ad2-ad1)*dvsp0
c    +            +3.0d0*app0*dvpp0+3.0d0*app1*dvpp1
c           aaa=0.25d0*aaa
c
c
c           aaa0=fss0+dsqrt(3.0d0)*(ad2-ad1)*fsp0
c    +            +3.0d0*app0*fpp0+3.0d0*app1*fpp1
c           aaa0=0.25d0*aaa0
c
            aaa0= dbs1*dbs2*fss0  
     +           +dbs1*ad2*fspi
     +           -dbs2*ad1*fspj
     +           +app0*fpp0+app1*fpp1
c
            aaax=aaa0*(-dxc)*2.0d0
            aaay=aaa0*(-dyc)*2.0d0
            aaaz=aaa0*(-dzc)*2.0d0
c              the factor 2 is that for (ij) and (ji)
c
            if (idngl(ig)+idngl(jg) .ne. 0) then
               aaa=0.0d0
               aaax=0.0d0
               aaay=0.0d0
               aaaz=0.0d0
            endif  
c
            bbbj=0.0d0
            bbbi=0.0d0
!
            bbbj=bbbj+dbs1*(-1.0d0/drr)*dvspi
            bbbi=bbbi-dbs2*(-1.0d0/drr)*dvspj
!            ----> second derivative for (sp sigma) interaction
!
!
            bbbj=bbbj+ad1*(-1.0d0/drr)*(dvpp0-dvpp1)
            bbbi=bbbi+ad2*(-1.0d0/drr)*(dvpp0-dvpp1)
!             ----> second deriv. for (pp sigma) (pp pi) int.
!
            bbbx=(bbbi*dbx1b+bbbj*dbx2b)*2.0d0
            bbby=(bbbi*dby1b+bbbj*dby2b)*2.0d0
            bbbz=(bbbi*dbz1b+bbbj*dbz2b)*2.0d0
c
            ahij4d=aaa
            dhij(ja1,ja2,jsd1,jsv2)=ahij4d
            dfij(1,ja1,ja2,jsd1,jsv2)=aaax+bbbx
            dfij(2,ja1,ja2,jsd1,jsv2)=aaay+bbby
            dfij(3,ja1,ja2,jsd1,jsv2)=aaaz+bbbz
c!!!            if (jsv2 >= 1) then
c!!!               write(6,*)'dhij=',ja1,ja2,jsv1,jsv2,ahij4d
c!!!               write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaax+bbbx
c!!!               write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaay+bbby
c!!!               write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaaz+bbbz
c!!!             endif  
c
          enddo
         enddo
 1100   continue
 1000 continue
c
c
      write(6,*)' ... suscessfully ended (LSES_SET_DHIJ)'
c     stop
c
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
