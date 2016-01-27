!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!zzz  @@@@ elses-hami-gaas-molteni.f @@@@@
!zzz  @@@@@ 2008/04/01 @@@@@
!ccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0725: Prepared
!0730: Continued; 'elses_set_hami_gaas_molteni'
!0809: Continued (NT07C-120p05)
!0813: Continued (NT07C-120p11)
!0815: Continued (NT07C-120p11)
!0907: A serious bug-fix for the rest part (NT97C-120,p75)
!0910: A non-essential bugfix in 'elses_force_molteni'
!ccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0401T.Hoshi: function dvec is replaced by dvec2
!            in which iperiodic is not used. 
!                   (NT08A-123p23, 0.00.14a-wrk05) 
!           : i_pbc_x, i_pbc_y, i_pbc_z are introduced.
!              in subroutine 'elses_set_hami_gaas_molteni'
!                            'elses_set_rest_gaas_Molteni'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Rest part in order-N cost
!     Note that the length unit is Angstrom, not a.u., in this routine.
! 
!
subroutine elses_set_hami_gaas_molteni
  use elses_mod_sel_sys,   only : c_system
  use elses_mod_md_dat,    only : itemd
  use elses_mod_ctrl,      only : i_verbose
  use elses_mod_phys_const,only : ev4au,angst
! use elses_mod_sim_cell,  only : iperiodic,ax,ay,az
  use elses_mod_sim_cell,  only : ax,ay,az,  i_pbc_x, i_pbc_y, i_pbc_z
  use elses_mod_tx,        only : tx,ty,tz,jsei
  use elses_mod_js4jsv,    only : js4jsv,jsv4js
  use elses_mod_jsv4jsd,   only : jsv4jsd,njsd
  use elses_mod_noav,      only : noav
  use elses_arr_dhij,      only : dhij
  use elses_arr_dfij,      only : dfij
  use elses_mod_orb1,      only : nval
  use elses_mod_orb2,      only : js2j,dbx, dby, dbz, idngl
  use elses_mod_multi,     only : ict4h
  use elses_mod_r_base,    only : r_base
  
  implicit none
  integer imode, ipe, npe, ierr
  integer jsv2,js2,njsd2,nss2,nval2
  integer jsd1,jsv1,js1,nss1,nval1,jg,ig,ja1,ja2,ioas1,ioas2,isym, isym2 
  integer ipair
  real(8) eunphys, ddd
  real(8) dxc,dyc,dzc,dd,drr,rnn
  real(8) dha,rat
  real(8) fac1,fac2,dfacdr1,dfacdr2,dargexp,dexpon
  real(8) dvss0,dvspi, dvspj, dvpp0,dvpp1,dvapi, dvapj
  real(8) fss0, fspi, fspj, fpp0, fpp1, fapi, fapj
  real(8) dbx1,dby1,dbz1,dbx2,dby2,dbz2,ad1,ad2
  real(8) dbx1b,dby1b,dbz1b,dbx2b,dby2b,dbz2b
  real(8) dbs1,dbs2,dba1,dba2
  real(8) app0,app1,aaa,aaa0,aaax,aaay,aaaz
  real(8) bbbx, bbby, bbbz, bbbi, bbbj
  real(8) dddx,potij,dphidr
  real(8) ahij4d 
  real(8) delta, ddtanh
  real(8) rcc, es0(2), eps(2), ep0(2), esp3a(2), esp3b(2)
  real(8) r0, rcut, rcut2
  real(8) d_onsite(3,2)
!       ----> d_onesite(i,j) :On-site energy, 
!                 E_s, E_p, E_{s*} for i=1,2,3
!                  Ga, As          for j=1,2
  real(8) d_hop_value(5,4)
!       ----> d_hop_func(i,j) :hopping value (as reference)
!           (ss \sigma), (sp \sigma), 
!               (pp, \sigma), (pp, \pi), (s*p, sigma) (i=1,5)
!            Ga-Ga, Ga-As, As-Ga, As-As for j=1,4
  real*8 d_hop_func(6,2)

  write(6,*)'@@ elses_set_hami_gaas_Molteni'
  write(6,*)'c_system=',c_system
  select case(c_system)
  case('GaAs.Molteni.1993')
  case default
     write(6,*)'ERROR:c_system=',c_system
     stop
  end select
  ierr=1

  ! Paremeter for compatibility
  ict4h=1
  if (ict4h .ne. 1) then
     write(6,*)'ERROR!(SET_DHIJ):ict4h=',ict4h
     stop
  endif

  eunphys=5.0d0

  write(6,*)'   ict4h=',ict4h
  write(6,*)' EUNPHYS=',eunphys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input energy parameters in eV, angstrom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! r0 = 2.45d0                   ! in angstrom
  r0 = 2.44782080379671d0       ! in angstrom
  delta = 0.1d0                 ! in angstrom 
  rcut = 1.3d0*r0               ! in angstrom
  rcut2 = rcut + 3.d0*delta     ! in angstrom

  d_onsite(1,1:2) = (/-2.656d0,-8.343d0/)
  d_onsite(2,1:2) = (/ 3.669d0, 1.041d0/)
  d_onsite(3,1:2) = (/ 6.739d0, 8.591d0/)
!
! d_hop_value(1,1:4) = (/ 0.000d0,  0.000d0,  0.000d0,  0.000d0/)
! d_hop_value(2,1:4) = (/ 0.000d0,  0.000d0,  0.000d0,  0.000d0/)
! d_hop_value(3,1:4) = (/ 0.000d0,  0.000d0,  0.000d0,  0.000d0/)
! d_hop_value(4,1:4) = (/ 0.000d0,  0.000d0,  0.000d0,  0.000d0/)
! d_hop_value(5,1:4) = (/ 0.000d0,  0.000d0,  0.000d0,  0.000d0/)
!
  d_hop_value(1,1:4) = (/-2.000d0, -1.613d0, -1.613d0, -1.420d0/)
  d_hop_value(2,1:4) = (/ 2.100d0,  2.504d0,  1.940d0,  2.100d0/)
  d_hop_value(3,1:4) = (/ 2.200d0,  3.028d0,  3.028d0,  3.100d0/)
  d_hop_value(4,1:4) = (/-0.670d0, -0.781d0, -0.781d0, -0.790d0/)
  d_hop_value(5,1:4) = (/ 2.000d0,  2.082d0,  2.097d0,  1.800d0/)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters in a.u.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  rcc = rcut2 / angst
  es0(1:2) = d_onsite(1,1:2) / ev4au
  eps(1:2) = ( d_onsite(2,1:2) - d_onsite(1,1:2) )/ ev4au
!
  ep0(1:2) = es0(1:2) + eps(1:2)
  esp3a(1:2) = 0.25d0*es0(1:2) + 0.75d0*ep0(1:2)
  esp3b(1:2) = 0.25d0*(ep0(1:2)-es0(1:2))
!
 if (i_verbose >= 1) then
    write(6,*)'  es0(1:2) [eV] =', es0(1:2)*ev4au
    write(6,*)'  ep0(1:2) [eV] =', ep0(1:2)*ev4au
    write(6,*)'esp3a(1:2) [eV] =',esp3a(1:2)*ev4au
    write(6,*)'esp3b(1:2) [eV] =',esp3b(1:2)*ev4au
 endif 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Error checking is added for too small system ( ax < 2 r_cut)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (ax <  2.0d0*rcc) then
     write(6,*)'ERROR!(LSES_SET_HAMI_01):ax,rcc=',ax,rcc
     write(6,*)'This version dose not support too small system'
     stop
  endif

  if (ay < 2.0d0*rcc) then
     write(6,*)'ERROR!(LSES_SET_HAMI_01):ay,rcc=',ay,rcc
     write(6,*)'This version dose not support too small system'
     stop
  endif

  if (az < 2.0d0*rcc) then
     write(6,*)'ERROR!(LSES_SET_HAMI_01):az,rcc=',az,rcc
     write(6,*)'This version dose not support too small system'
     stop
  endif

  do  jsv2=1,noav
     js2=js4jsv(jsv2)
     njsd2=njsd(jsv2,ict4h)
     nss2=jsei(js2)
     nval2=nval(nss2)

     do jsd1=1,njsd2
        dhij(:,:,jsd1,jsv2)=0.0d0
        dfij(:,:,:,jsd1,jsv2)=0.0d0

        jsv1=jsv4jsd(jsd1,jsv2)
        js1=js4jsv(jsv1)
        nss1=jsei(js1)
        nval1=nval(nss1)
!
!       write(6,*)'js1,js2,nss1,nss2=',js1,js2,nss1,nss2
!
        if (js1 .eq. js2) then
           ioas1 = jsei(js1)
!
           if (i_verbose >= 1) then
             ddd=0.0d0
             do ja2=1,nval2
                jg=js2j(ja2,js2)
                do ja1=1,nval1
                  ig=js2j(ja1,js1)
                  dbs1=r_base(1,ig)
                  dbx1=r_base(2,ig)
                  dby1=r_base(3,ig)
                  dbz1=r_base(4,ig)
                  dba1=r_base(5,ig)
!
                  dbs2=r_base(1,jg)
                  dbx2=r_base(2,jg)
                  dby2=r_base(3,jg)
                  dbz2=r_base(4,jg)
                  dba2=r_base(5,jg)
!                
                  ddd=dbs1*dbs2+dbx1*dbx2+dby1*dby2+dbz1*dbz2+dba1*dba2
                  if (ig == jg) ddd=ddd-1.0d0
!                 if (js1 <= 2) then
!                    write(6,*)'js, onsite-non-ortho=',js1,ddd
!                 endif 
                  if (dabs(ddd) .gt. 1.0d-14) then
                    write(6,*)'js, onsite-non-ortho=',js1,ddd
                    stop
                  endif  
                enddo
             enddo   
           endif  
!
!
           do ja2=1,nval2
              jg=js2j(ja2,js2)
              do ja1=1,nval1
                 ig=js2j(ja1,js1)
                 if (ja1 == ja2) then
                    ahij4d=esp3a(ioas1)
                    if (ja1 == 5 ) ahij4d=d_onsite(3,ioas1)/ev4au
                    if (idngl(ig) .eq. 1) ahij4d=eunphys
                 else  
                    ahij4d=esp3b(ioas1)
                    if ((ja1 == 5) .or. (ja2 == 5)) then
                       ahij4d=0.0d0
                    endif   
                    if (idngl(ig) .eq. 1) ahij4d=0.0d0
                    if (idngl(jg) .eq. 1) ahij4d=0.0d0
                 endif
                 dhij(ja1,ja2,jsd1,jsv2)=ahij4d
                 if (i_verbose >= 1) then
!                  if (js1 <= 2) then
!                    write(6,*)'js1,ja1,ja2,H[eV] =',js1,ja1,ja2,ahij4d*ev4au
!                  endif 
                 endif   
              enddo
           enddo
           cycle
        endif
!       -----> on-site terms

        dxc=tx(js2)-tx(js1)
        dyc=ty(js2)-ty(js1)
        dzc=tz(js2)-tz(js1)

        if (i_pbc_x == 1) dxc = modulo(dxc + 0.5d0, 1.0d0) - 0.5d0
        if (i_pbc_y == 1) dyc = modulo(dyc + 0.5d0, 1.0d0) - 0.5d0
        if (i_pbc_z == 1) dzc = modulo(dzc + 0.5d0, 1.0d0) - 0.5d0

!       if (iperiodic == 1) then
!          dd=1.0d0
!          dxc=dxc+0.5d0*dd
!          dxc=dmod(dxc+dd,dd)-0.5d0*dd
!          dyc=dyc+0.5d0*dd
!          dyc=dmod(dyc+dd,dd)-0.5d0*dd
!          dzc=dzc+0.5d0*dd
!          dzc=dmod(dzc+dd,dd)-0.5d0*dd
!       endif

        dxc=dxc*ax
        dyc=dyc*ay
        dzc=dzc*az
        drr=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
!        ---> Distance | R_1 - R_2 | in a.u.

        if (drr >= rcc) cycle
        if (drr <= 1.00d0) then
           write(6,*)'Too near atoms?:',itemd,js1,js2,drr
        endif
        if (drr <= 1.0d-3) then
           stop
        endif
!
        dxc=dxc/drr
        dyc=dyc/drr
        dzc=dzc/drr
!       ---> Direction (unit) vector : R_2 - R_1

        rnn=drr*angst
!       ---> Distance in [A]

        ioas1 = jsei(js1)
        ioas2 = jsei(js2)
        if( ioas1 == 1 .and. ioas2 == 1 ) then     ! Ga-Ga
           ipair = 1
        elseif( ioas1 == 1 .and. ioas2 == 2 ) then ! Ga-As
           ipair = 2
        elseif( ioas1 == 2 .and. ioas2 == 1 ) then ! As-Ga
           ipair = 3
        elseif( ioas1 == 2 .and. ioas2 == 2 ) then ! As-As
           ipair = 4
        else
           write(*,*) 'ERROR!,wrong value of jsei',jsei(js1),jsei(js2)
           stop
        end if

!
!!!!!!!!!Parameter for Ga-Ga case
!  
        if (ipair .eq. 1) then
          dvss0=d_hop_value(1,1)
          dvspi=d_hop_value(2,1)
          dvspj=d_hop_value(2,1)
          dvpp0=d_hop_value(3,1)
          dvpp1=d_hop_value(4,1)
          dvapi=d_hop_value(5,1)
          dvapj=d_hop_value(5,1)
        endif  
!
!!!!!!!!!Parameter for Ga-As case
!
        if (ipair .eq. 2) then
          dvss0=d_hop_value(1,2)
          dvspi=d_hop_value(2,2)
          dvspj=d_hop_value(2,3)
          dvpp0=d_hop_value(3,2)
          dvpp1=d_hop_value(4,2)
          dvapi=d_hop_value(5,2)
          dvapj=d_hop_value(5,3)
        endif  
!
!!!!!!!!!Parameter for As-Ga case
!
        if (ipair .eq. 3) then
          dvss0=d_hop_value(1,3)
          dvspi=d_hop_value(2,3)
          dvspj=d_hop_value(2,2)
          dvpp0=d_hop_value(3,3)
          dvpp1=d_hop_value(4,3)
          dvapi=d_hop_value(5,3)
          dvapj=d_hop_value(5,2)
        endif  
!
!!!!!!!!!Parameter for As-As case
!
        if (ipair .eq. 4) then
          dvss0=d_hop_value(1,4)
          dvspi=d_hop_value(2,4)
          dvspj=d_hop_value(2,4)
          dvpp0=d_hop_value(3,4)
          dvpp1=d_hop_value(4,4)
          dvapi=d_hop_value(5,4)
          dvapj=d_hop_value(5,4)
        endif  
!
!        ---> Slator-Koster parameters
!
        rat = r0/rnn
        fac1 = rat ** 2.d0
        ddtanh=dtanh((rnn-rcut)/delta)
        fac2 = 0.5d0*(1.d0-ddtanh)
!          --> f(x) = 0.5 ( 1 - tanh(x) )
        dfacdr1  = - 2.0d0* fac1 / rnn
        dfacdr2 = -0.5d0/delta*(1.0d0-ddtanh*ddtanh)
!          --> f'(x) = -0.5 ( 1 - tanh^2(x) )
!
        fss0=dvss0/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
        fspi=dvspi/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
        fspj=dvspj/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
        fpp0=dvpp0/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
        fpp1=dvpp1/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
        fapi=dvapi/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
        fapj=dvapj/ev4au*angst*(fac1*dfacdr2+dfacdr1*fac2)
!         ---> Derivative of Slator-Koster parameters in au
!
        dvss0=dvss0/ev4au*fac1*fac2
        dvspi=dvspi/ev4au*fac1*fac2
        dvspj=dvspj/ev4au*fac1*fac2
        dvpp0=dvpp0/ev4au*fac1*fac2
        dvpp1=dvpp1/ev4au*fac1*fac2
        dvapi=dvapi/ev4au*fac1*fac2
        dvapj=dvapj/ev4au*fac1*fac2
!         ---> Derivative of Slator-Koster parameters in au
!
!
!         if (i_verbose >= 1) then
!           if (js1 .eq. 1) then
!             write(6,*)'js1,js2   =',js1,js2
!             write(6,*)' rnn [A]  =',rnn
!             write(6,*)'  fac1    =',fac1
!             write(6,*)'  fac2    =',fac2
!             write(6,*)' dfacdr1  =',dfacdr1
!             write(6,*)' dfacdr2  =',dfacdr2
!             write(6,*)' dvss0[eV]=',dvss0*ev4au
!             write(6,*)' dvspi[ev]=',dvspi*ev4au
!             write(6,*)' dvspj[eV]=',dvspj*ev4au
!             write(6,*)' dvpp0[eV]=',dvpp0*ev4au
!             write(6,*)' dvpp1[eV]=',dvpp1*ev4au
!             write(6,*)' dvapi[eV]=',dvapi*ev4au
!             write(6,*)' dvapj[eV]=',dvapj*ev4au
!           endif   
!         endif   
!
          do ja2=1,nval2
             jg=js2j(ja2,js2)
             do ja1=1,nval1
                ig=js2j(ja1,js1)
!
                dbs1=r_base(1,ig)
                dbx1=r_base(2,ig)
                dby1=r_base(3,ig)
                dbz1=r_base(4,ig)
                dba1=r_base(5,ig)
!
                dbs2=r_base(1,jg)
                dbx2=r_base(2,jg)
                dby2=r_base(3,jg)
                dbz2=r_base(4,jg)
                dba2=r_base(5,jg)
!
                ! Inner products
                ad2=dbx2*dxc+dby2*dyc+dbz2*dzc
!                 ----->inner product ( a_2 | d )
                ad1=dbx1*dxc+dby1*dyc+dbz1*dzc
!                 ----->inner product  ( a_1 | d )

                !  Vector :  b_1 = a_1 - (a_1 | d) d
                dbx1b=dbx1-ad1*dxc
                dby1b=dby1-ad1*dyc
                dbz1b=dbz1-ad1*dzc
!
                !  Vector :  b_2 = a_2 - (a_2 | d) d
                dbx2b=dbx2-ad2*dxc
                dby2b=dby2-ad2*dyc
                dbz2b=dbz2-ad2*dzc
!
                ! < p_1 | p_2 > parts 
                app0=ad1*ad2
!                 ----->double inner product : ( a_1 | d ) ( a_2 | d ) 
                app1=dbx1b*dbx2b+dby1b*dby2b+dbz1b*dbz2b
!                 ----->inner product : ( b_1 | b_2 )

                aaa= dbs1*dbs2*dvss0   & 
                    +dbs1*ad2*dvspi    &
                    -dbs2*ad1*dvspj    &
                    +app0*dvpp0+app1*dvpp1 &
                    +dba1*ad2*dvapi    &
                    -dba2*ad1*dvapj    
!
                aaa0= dbs1*dbs2*fss0   & 
                     +dbs1*ad2*fspi    &
                     -dbs2*ad1*fspj    &
                     +app0*fpp0+app1*fpp1 &
                     +dba1*ad2*fapi    &
                     -dba2*ad1*fapj    
!
!               
                aaax=aaa0*(-dxc)*2.0d0
                aaay=aaa0*(-dyc)*2.0d0
                aaaz=aaa0*(-dzc)*2.0d0
!              the factor 2 is that for (ij) and (ji)
!
                if (idngl(ig)+idngl(jg) .ne. 0) then
                   aaa=0.0d0
                   aaax=0.0d0
                   aaay=0.0d0
                   aaaz=0.0d0
                endif
!
                ahij4d=aaa
                dhij(ja1,ja2,jsd1,jsv2)=ahij4d
!               dfij(1,ja1,ja2,jsd1,jsv2)=aaax
!               dfij(2,ja1,ja2,jsd1,jsv2)=aaay
!               dfij(3,ja1,ja2,jsd1,jsv2)=aaaz
!                   ----> first derivative (w.r.t. position of the js)
!
                bbbj=0.0d0
                bbbi=0.0d0
!
                bbbj=bbbj+dbs1*(-1.0d0/drr)*dvspi
                bbbi=bbbi-dbs2*(-1.0d0/drr)*dvspj
!                   ----> second derivative for (sp sigma) interaction
!
                bbbj=bbbj+dba1*(-1.0d0/drr)*dvapi
                bbbi=bbbi-dba2*(-1.0d0/drr)*dvapj
!                   ----> second derivative for (s*p sigma) interaction
!                    
                bbbj=bbbj+ad1*(-1.0d0/drr)*(dvpp0-dvpp1)
                bbbi=bbbi+ad2*(-1.0d0/drr)*(dvpp0-dvpp1)
!                   ----> second derivative for (pp sigma) (pp pi) interactions
!                    
                bbbx=(bbbi*dbx1b+bbbj*dbx2b)*2.0d0
                bbby=(bbbi*dby1b+bbbj*dby2b)*2.0d0
                bbbz=(bbbi*dbz1b+bbbj*dbz2b)*2.0d0
!
!               write(6,*)'js1,js2,ja1,ja2,drr=',js1,js2,ja1,ja2,drr
!               write(6,*)'aaa =',aaa
!               write(6,*)'aaax=',aaax
!               write(6,*)'aaay=',aaay
!               write(6,*)'aaaz=',aaaz
!               write(6,*)'bbbx=',bbbx
!               write(6,*)'bbby=',bbby
!               write(6,*)'bbbz=',bbbz
                dfij(1,ja1,ja2,jsd1,jsv2)=aaax+bbbx
                dfij(2,ja1,ja2,jsd1,jsv2)=aaay+bbby
                dfij(3,ja1,ja2,jsd1,jsv2)=aaaz+bbbz
!               if (jsv2 <= 1) then
!                 write(6,*)'dhij=',ja1,ja2,jsv1,jsv2,ahij4d
!                 write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaax+bbbx
!                 write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaay+bbby
!                 write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,aaaz+bbbz
!               endif  
!
             enddo
          enddo
          
     end do
  end do

  return
end subroutine elses_set_hami_gaas_molteni
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elses_set_rest_gaas_Molteni
  use elses_mod_sel_sys, only : c_system
  use elses_mod_phys_const, only : angst,ev4au
  use elses_mod_sim_cell,   only : noa, ax, ay, az
  use elses_mod_sim_cell,   only : i_pbc_x, i_pbc_y, i_pbc_z
  use elses_mod_tx,         only : tx, ty, tz
  use elses_mod_foi,        only : foi
  use elses_mod_jsv4jsd,    only : jsv4jsd,njsd
  use elses_mod_noav,       only : noav,noao,nncut
  use elses_mod_ene,        only : ecc, enea

  implicit none
  integer ierr, ishow
  integer js, jsv, ioa, joa, joad, noa1, njsd2, nbkedd
  real(8) alpha, delta
  real(8) phi1, phi2
  real(8) r0, rcut, rcut2
  real(8) erep,rxi,ryi,rzi
  real(8) axang, ayang, azang
  real(8) rd, rd2, poti, potij, dphidr
  real(8) term1, term2, tail, ddtanh
  real(8) dtermdr1, dtermdr2, dtaildr
  real(8) rxx,ryy,rzz,dvec2

  ! working arrays
  integer, allocatable, dimension(:,:) :: l4oij
  integer, allocatable, dimension(:,:) :: nbke
  real(8),  allocatable, dimension(:,:) :: force(:,:)
  real(8),  allocatable, dimension(:) :: f1d

  write(6,*)'@@ elses_set_rest_gaas_Molteni'
  write(6,*)'c_system=',c_system
  select case(c_system)
  case('GaAs.Molteni.1993')
  case default
     write(6,*)'ERROR:c_system=',c_system
     stop
  end select

  ierr = 1

  ! parameters for rest part energy in eV, angstrom
  alpha = 0.3555d0              ! in angstrom
! r0 = 2.45d0                   ! in angstrom
  r0 = 2.44782080379671d0       ! in angstrom
  delta = 0.1d0                 ! in angstrom 
  rcut = 1.3d0*r0               ! in angstrom
  rcut2 = 1.3d0*r0 + 3.d0*delta ! in angstrom
  phi1 = 2.3906d0               ! in eV
  phi2 = 1.2347d0               ! in eV

  ishow=5

  !  a.u. -> Angstrom

  axang=ax*angst
  ayang=ay*angst
  azang=az*angst

  write(6,*) '@@ lses_set_rest_01:ax=',axang, ayang, azang

  ierr=0
  if (axang .le. 2.0d0*rcut) ierr=1
  if (ayang .le. 2.0d0*rcut) ierr=1
  if (azang .le. 2.0d0*rcut) ierr=1
  if (ierr == 1) then
     write(6,*) 'error!(repluc):rcut,axang=',rcut,axang
     write(6,*) 'Too small cell size for cutoff radius?!'
     stop
  endif

  ! checking the consistency
  if  (nncut /= 2) then
     write(6,*)'alloc. error!(REPULC)NNCUT=',nncut
     stop
  endif

  if  (noav /= noa) then
     write(6,*)'alloc. error!(REPULC)noa,noav=',noa,noav
     stop
  endif

  if (.not. allocated(njsd) ) then
     write(6,*)'ERROR!:Not yet allocated: NJSD'
     stop
  endif

  if (.not. allocated(jsv4jsd) ) then
     write(6,*)'ERROR!:Not yet allocated: JSV4JSD'
     stop
  endif

  if (.not. allocated(enea) ) then
     write(6,*)'ERROR!:Not yet allocated: ENEA'
     stop
  endif

  ! allocation and setting working arrays

  allocate (l4oij(noao,noa),stat=ierr)
  if( ierr /= 0) then
     write(6,*)'alloc. error!(L4OIJ):ierr=',ierr
     stop
  endif
  l4oij(:,:)=0

  allocate (force(noa,3),stat=ierr)
  if( ierr /= 0) then
     write(6,*)'alloc. error!(force):ierr=',ierr
     stop
  endif
  force(:,:)=0.0d0

  allocate (f1d(noa),stat=ierr)
  if( ierr /= 0) then
     write(6,*)'alloc. error!(force):ierr=',ierr
     stop
  endif
  f1d(:)=0.0d0

  allocate (nbke(noav,0:nncut),stat=ierr)
  if( ierr /= 0) then
     write(6,*)'alloc. error!(NBKE):ierr=',ierr
     stop
  endif
  nbke(:,:)=0

  do js=1,noa
     jsv=js
     njsd2=njsd(jsv,1)
     nbke(jsv,1)=njsd2
     if ((njsd2 <= 0) .or. (njsd2 > noao)) then
        write(6,*)'Error!(REPULC)js,njsd2=',js,njsd2
        stop
     endif
     l4oij(1:njsd2,jsv)=jsv4jsd(1:njsd2,jsv)
  enddo

  ! initialization of the repulsive energy E_{rep}

  erep=0.d0
  do ioa = 1, noa
     noa1 = nbke(ioa,1)

     rxi = tx(ioa)*axang
     ryi = ty(ioa)*ayang
     rzi = tz(ioa)*azang

     if( ioa < ishow ) then
        nbkedd = nbke(ioa,1)
        write(6,*) 'ioa,nbke=',ioa, nbkedd
     end if

     ! initialization of the repulsive pair potential sum_{j}phi(r_{ij})

     poti=0.d0
     force(1:noa1,1)=0.D0
     force(1:noa1,2)=0.D0
     force(1:noa1,3)=0.D0

     !  calculation of distance between joa and ioa ion
     do joad = 1, nbke(ioa,1)
        joa = l4oij(joad,ioa)
        
        if( ioa == joa ) cycle

        if( ( joa <= 0 ) .or. ( joa > noa ) ) then 
           write(6,*) 'error!(zrepul2)'
           write(6,*) 'joa,joad,ioa,nbke=',joa,joad,ioa,nbke(ioa,1)
           stop
        end if

        rxx=dvec2(tx(ioa),tx(joa),1.0d0,ierr,i_pbc_x)*axang 
        ryy=dvec2(ty(ioa),ty(joa),1.0d0,ierr,i_pbc_y)*ayang 
        rzz=dvec2(tz(ioa),tz(joa),1.0d0,ierr,i_pbc_z)*azang 

        rd2=rxx*rxx+ryy*ryy+rzz*rzz
        if(rd2 > rcut2*rcut2) cycle
        rd =dsqrt(rd2)

        if(rd < 1.0d-10) then
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
        
        !    Calculation of the repulsive pair potential and force
        !    (from joa to ioa)
!
        term1 = 0.5d0*phi1*dexp(-(rd-r0)/alpha)
        term2 = 0.5d0*phi2*(r0/rd)
        ddtanh=dtanh((rd-rcut)/delta)
        tail = 0.5d0*(1.d0-ddtanh)

        dtermdr1 = -term1/alpha
        dtermdr2 = -term2/rd
        dtaildr = -0.5d0/delta*( 1.d0 - ddtanh*ddtanh )

        potij =  (term1+term2) * tail
        dphidr = (dtermdr1+dtermdr2)*tail + (term1+term2)*dtaildr 

        poti=poti+potij
        
        force(joad,1)=-dphidr*(rxx/rd)*2.0d0
        force(joad,2)=-dphidr*(ryy/rd)*2.0d0
        force(joad,3)=-dphidr*(rzz/rd)*2.0d0
!             ---> minus sign : - (d E / d R)
!             ---> factor 2   : counting (ij) and (ji)
!
     end do

     erep = erep + poti
     do  joad=1, noa1
        joa=l4oij(joad,ioa)
        if(joa.eq.ioa) cycle
        foi(ioa,1)=foi(ioa,1)+force(joad,1)*(angst/ev4au)
        foi(ioa,2)=foi(ioa,2)+force(joad,2)*(angst/ev4au)
        foi(ioa,3)=foi(ioa,3)+force(joad,3)*(angst/ev4au)
     end do
     write(*,*) ioa
  end do

  write(*,*) 'Erep [eV,a.u.]= ', erep, erep/ev4au
  ecc = erep/ev4au

  return
end subroutine elses_set_rest_gaas_Molteni
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine elses_force_molteni
!   
      use elses_mod_phys_const,only : ev4au,angst
      use elses_mod_ctrl,      only : i_verbose
!     use elses_mod_sim_cell,  only : noa, iperiodic,ax,ay,az
      use elses_mod_sim_cell,  only : noa, ax,ay,az
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
!
      implicit none
      integer imode, ipe, npe
      real*8  detb, dens2
      integer js,nss,ja
!
      real*8  dnal0,rnn0,rcc,es0,eps,ep0,esp3a,esp3b
      real*8  qc0,qc1,qc2,qc3,r_cut_tail
!
      integer jsd2,jsv2,js2,njsd2,nss2,nval2,js2b
      integer jsd1,jsd1b,jsv1,js1,nss1,njsd1, nval1,jsv1b,js1b
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
!
      real*8  dddx1, dddy1, dddz1
      real*8  dbx0, dby0 , dbz0, pinner
      real*8  fdx1, fdy1 , fdz1, fdx2, fdy2 , fdz2
      real*8  acc0, acc3, acc4
      real*8  dddx2, dddy2, dddz2
      real*8  dddx3, dddy3, dddz3
      real*8  dddx4, dddy4, dddz4
      real*8  ddsum1, ddsum2, ddsum3, ddsum
!
      real*8 dhal(4),dnal(4),rcal(4)
      real*8 dkwondd(4,2)
!       ----> work array for TB parameters
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      imode=1
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,*)'@@ elses_force_molteni'
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     if (nvl .ne. 5) then
!       write(6,*)'ERROR!(SETDHIJ):UNSUPPORTED!:NVL=',nvl
!       stop
!     endif
!
      if (imode .ne. 1) then
        WRITE(6,*)'ERROR!(SETDHIJMP2):not supported imode=',imode
        stop
      endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Plot the result (optional)
!
      ddsum1=0.0d0
      ddsum2=0.0d0
      ddsum3=0.0d0
      do js=1,noa
         if ( i_verbose >= 1 ) then
           if (js .le. 1) then
             write(6,*)'js,foi(tot)=',js,foi(js,1),foi(js,2),foi(js,3)   
           endif
         endif  
         ddsum1=ddsum1+foi(js,1)
         ddsum2=ddsum2+foi(js,2)
         ddsum3=ddsum3+foi(js,3)
      enddo
!
      ddsum=(dabs(ddsum1)+dabs(ddsum2)+dabs(ddsum3))/dble(noa)
!
      if ( i_verbose >= 1 ) then
        write(6,*)'Summed force (rep,x)/NOA=',ddsum1/dble(noa)
        write(6,*)'Summed force (rep,y)/NOA=',ddsum2/dble(noa)
        write(6,*)'Summed force (rep,z)/NOA=',ddsum3/dble(noa)
        write(6,*)'Summed force (rep,+)/NOA=',ddsum/dble(noa)
      endif  
!
      if (ddsum .gt. 1.0d-10) then
         write(6,*)'ERROR?:Summed force (rep,+)/NOA=',ddsum/dble(noa)
         stop
      endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      DETB=0.0D0
      DENS2=0.0D0
!
      ddd0=0.0d0
      ipe=0
      npe=0
!     ipe=omp_get_thread_num()+1
!     npe=omp_get_num_threads()
!     write(6,*)'ipe,npe=',ipe,npe
      do jsv1=1,noav
         js=js4jsv(jsv1)
!        write(6,*)'jsv1,js=',jsv1,js
         if ((js .le. 0) .or. (js .gt. noa)) then
           write(6,*)'error!(fortstv):jsv1,js=',jsv1,js
           stop
         endif
         js1=js
         dddx=0.0d0
         dddy=0.0d0
         dddz=0.0d0
         nss=jsei(js)
         nval1=nval(nss)
!
       njsd1=njsd(jsv1,ict4h)
       do jsd2=1,njsd1
          jsv2=jsv4jsd(jsd2,jsv1)
!         write(6,*)'jsv1,jsv2=',jsv1,jsv2
          if (jsv1 .eq. jsv2) cycle
          if ((jsv2 .le. 0) .or. (jsv2 .gt. noav)) then
            write(6,*)'error!(fortstv)'
            write(6,*)'jsv2,jsd2,jsv1=',jsv2,jsd2,jsv1
            stop
          endif
          js2=js4jsv(jsv2)
          if ((js2 .le. 0) .or. (js2 .gt. noa)) then
            write(6,*)'error!(fortstv):js2,jsv2=',js2,jsv2
            stop
          endif
          nss2=jsei(js2)
          nval2=nval(nss2)
!
          jsd1=0
          njsd2=njsd(jsv2,ict4h)
          do jsd1b=1,njsd2
            jsv1b=jsv4jsd(jsd1b,jsv2)
            if (jsv1b .eq. jsv1) then
               jsd1=jsd1b
               exit
            endif   
          enddo
          if (jsd1 == 0) then 
            write(6,*)'ERROR!(FORSTVMP2):JSV1,JSV2=',JSV1,JSV2
            stop
          endif
!            ---> jsd1 is generated.
!
          do ja=1,nval1
            ja1=ja 
            i=js2j(ja,js)
!           write(6,*)' loop 120:js,js2=',js,js2
            js1b=j2js(i)
            if (js1 .ne. js1b) then
             write(6,*)'error!(fortstv):js1,js1b=',js1,js1b
             stop
            endif
!
          do ja2=1,nval2
!            write(6,*)' loop 130:js,js2=',js,js2
             j=js2j(ja2,js2)
             js2b=j2js(j)
             if (js2 .ne. js2b) then
              write(6,*)'error!(fortstv):js2,js2b=',js2,js2b
              stop
             endif
!
             ddalij=dbij(ja2,ja1,jsd2,jsv1)
             dddx1=2.0d0*ddalij*dfij(1,ja1,ja2,jsd1,jsv2)
             dddy1=2.0d0*ddalij*dfij(2,ja1,ja2,jsd1,jsv2)
             dddz1=2.0d0*ddalij*dfij(3,ja1,ja2,jsd1,jsv2)
!                 ---> 2.0d0 : para-spin factor
!
!           if (jsv2 <= 1) then
!              write(6,*)'dbij=',ja1,ja2,jsv1,jsv2,ddalij
!              write(6,*)'dhij=',ja1,ja2,jsv1,jsv2,dhij(ja1,ja2,jsd1,jsv2)
!              write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,dfij(1,ja1,ja2,jsd1,jsv2)
!              write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,dfij(2,ja1,ja2,jsd1,jsv2)
!              write(6,*)'dfij=',ja1,ja2,jsv1,jsv2,dfij(3,ja1,ja2,jsd1,jsv2)
!            endif  
!
             dddx=dddx+dddx1
             dddy=dddy+dddy1
             dddz=dddz+dddz1
!
          enddo   
          enddo   
!
       enddo
!
       foi(js,1)=foi(js,1)-dddx
       foi(js,2)=foi(js,2)-dddy
       foi(js,3)=foi(js,3)-dddz
!             ---> minus sign : - (d E / d R)
!
      enddo
!
  199 CONTINUE
!
!     do js=1,noa
!        if (js .le. 1) then
!           write(6,*)'js,foi(tot)=',js,foi(js,1),foi(js,2),foi(js,3)   
!        endif
!     enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Plot the result (optional)
!
      ddsum1=0.0d0
      ddsum2=0.0d0
      ddsum3=0.0d0
      do js=1,noa
         if (js .le. 1) then
           write(6,*)'js,foi(tot)=',js,foi(js,1),foi(js,2),foi(js,3)   
         endif
         ddsum1=ddsum1+foi(js,1)
         ddsum2=ddsum2+foi(js,2)
         ddsum3=ddsum3+foi(js,3)
      enddo
      write(6,*)'Summed force (all,x)/NOA=',ddsum1/dble(noa)
      write(6,*)'Summed force (all,y)/NOA=',ddsum2/dble(noa)
      write(6,*)'Summed force (all,z)/NOA=',ddsum3/dble(noa)
      ddsum=(dabs(ddsum1)+dabs(ddsum2)+dabs(ddsum3))/dble(noa)
      write(6,*)'Summed force (all,+)/NOA=',ddsum/dble(noa)
      if (ddsum .gt. 1.0d-10) then
         write(6,*)'ERROR?:Summed force (all,+)/NOA=',ddsum/dble(noa)
         stop
      endif   
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

end subroutine elses_force_molteni

