!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
c  @@@@@ elses-eig-leg2 @@@@@
c  @@@@@ 2008/08/22 @@@@@
c2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0517:Prepared (Nt07B-119,p15) 
c       --> Eigen-state solver for standard eigen-value problem
c0823:A. modi.
c2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0105T.Hoshi; 'i_eig_main_init' is introduced for 'elses-eig_mai'
c                 (NT07E-122p31, v000a-wrk05)
c           ; 'elses_set_eig_leg_atmp2' is changed 
c                 into the standard eig-value problem
c0107T.Hoshi; electronic temperature is set;
c                temp_for_electron --> fb (NT07E-122p39)
c           ; 'call elses_alloc_dsij_dpij' 
c                 for compatibility to routines in NRL calc.
c0108T.Hoshi; 'call elses_enetb(etbtmp)' is added 
c               (NT07E-122p39, v000a-wrk05)
c0822T.Hoshi; Unused module variable 'fb' is deleted (NT08E127p13)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Main routine for eigen-state solver (standard eigen value eq.)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_eig_main(imode,i_eig_main_init)
c
      use elses_mod_ctrl,       only : i_verbose
      use elses_mod_phys_const, only : ev4au
      use elses_mod_sim_cell, only : noa
      use elses_mod_val_elec, only : val_elec_tot
      use elses_mod_noav,     only : noas
c     use elses_mod_eig_leg,  only : fb, tot_elec_eig_leg
      use elses_mod_eig_leg,  only : tot_elec_eig_leg
      use elses_mod_orb2,     only : n_tot_base
c     use elses_mod_elec_cond, only : temp_for_electron
      use elses_mod_ene,      only : etb
      use elses_param_ctl_kr, only : noav_kr_def
c
      implicit none
      integer :: nkeig, imode
      integer :: i_chk_eig_states
      real(8) :: rNelec
      real(8) :: time8, tb8, tb8tot
      integer :: i_eig_main_init
      real(8) :: etbtmp
c
      call tclock(time8)
      tb8=time8
      tb8tot=time8
c
      call elses_set_val_elec_tot
      rNelec=val_elec_tot
      tot_elec_eig_leg=val_elec_tot
      if (i_verbose >= 1) then
        write(6,*)'@@ LSES_EIG_MAIN:rNelec=',rNelec
      endif  
c
      if (rNelec .le. 0.01d0) then
        write(6,*)'ERROR!(LSES_EIG_MAIN):rNelec=',rNelec
        stop
      endif
c   
      if (noa .gt. 2000) then
        write(6,*)'Too large for diagonalization ?'
        write(6,*)'NOA=',noa
        stop
      endif
c   
      if (noas .gt. noa) then
        write(6,*)'Inconsistent pararameter!!'
        write(6,*)'NOAS,NOA=',noas,noa
        stop
      endif   
c
      if (i_eig_main_init == 1) then
        call elses_alloc_eig_leg(n_tot_base)
c          ---> allocate arraies : atmp, atmp2 etc.
ccc$        call elses_alloc_eig_leg_wrk
c          ---> allocate work arraies used in lset_mat_eig
c       call elses_init_alloc_nrl_leg
c          ---> allocate arrays  : dsij, dpij
c             ( unnecesary, used for compatibility to 
c                several routines in NRL calc. )
c       call elses_eig_set_dsij_to_unity
c          ---> set S = I ( for compatibility )
      endif   
c
      call elses_set_eig_leg_atmp2
c         : copy DHIJ(Hamiltonian) to ATMP
c
      nkeig=anint(rNelec/2.0d0)
c     fbd=0.025D0
c     fbd=1.0d0/200.0d0*ev4au
c     fbd=temp_for_electron*ev4au
c     fb=fbd
c     fbd=-1.0d0
c      ---> dummy variable (NOT USED ACTUALLY)
      if (i_verbose >= 1) then
        write(6,*)'rNelec =',rNelec
        write(6,*)'nkeig  =',nkeig
c       write(6,*)'fbd[eV,au] =',fbd,fbd/ev4au
      endif  
c
      call tclock(time8)
      if (i_verbose >= 1) then
        WRITE(6,*)'EIGTIME (befor SETEIG)=',time8-tb8
      endif  
      tb8=time8
c
      imode=1
c     imode=2
c        =1 (standard), =2 (generalized) =3 (overlap-mat)
      call elses_seteig3(nkeig,rNelec,imode)
      call tclock(time8)
      if (i_verbose >= 1) then
        WRITE(6,*)'EIGTIME (SETEIG)     =',time8-tb8
      endif
      tb8=time8
c
      call elses_eig_chem_pot
c          ---> calculate chemical potential and occ. num.
      call tclock(time8)
      if (i_verbose >= 1) then
        WRITE(6,*)'EIGTIME (CHEM_POT)   =',time8-tb8
      endif  
      tb8=time8
c
      call elses_eig_set_dens_mat
c        ----> density matrix calc :
c           f(k,1) : occupation number (k:=1,neig)
c                   should be given before
c           Note:dpij is generated (useless)
c
      call tclock(time8)
      if (i_verbose >= 1) then
        WRITE(6,*)'EIGTIME (CALLIJ4_EIG)=',time8-tb8
      endif  
      tb8=time8
c
c
      noav_kr_def=noa
      etbtmp=0.0d0
      call elses_enetb(etbtmp)
      if (i_verbose >= 1) then
        write(6,*)'ETB=,',etbtmp
      endif  
      etb=etbtmp
c
      call tclock(time8)
      if (i_verbose >= 1) then
        WRITE(6,*)'EIGTIME (ENETB   _EIG)=',time8-tb8
      endif  
      tb8=time8
c
      call tclock(time8)
      if (i_verbose >= 1) then
        WRITE(6,*)'EIGTIME (total)      =',time8-tb8tot
      endif  
c
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Set the matrix ATMP, ATMP2
c          for standard eigen value problem
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_eig_leg_atmp2
      use elses_mod_eig_leg, only : n_base_eig_leg
      use elses_arr_eig_leg, only : atmp, atmp2, idngl2
      use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base,idngl
c
c     use elses_mod_sim_cell, only : noa,nos,iperiodic,ax,ay,az
      use elses_mod_tx,       only : jsei
      use elses_mod_noav,     only : noav
      use elses_mod_js4jsv,   only : js4jsv, jsv4js
      use elses_arr_dhij,     only : dhij
c     use elses_arr_dsij,     only : dsij
c     use elses_arr_dpij,     only : dpij
      use elses_mod_jsv4jsd,  only : jsv4jsd,njsd
      use elses_mod_orb1,     only : nvl, nval
      use elses_mod_multi,    only : ict4h
c 
      implicit none
      integer :: neig0
      integer :: nn, i, j, isum, jsv, js, ja, jg, iii, jsv2
      integer :: njsd2, nss2, nval2, jsd1, jsv1, js1, nss1
      integer :: js2, nval1, ja2, ja1, ig, jj 
c     integer :: IDNGL2(NEIG0)
      real(8) :: time2, tb2, dbigd, dsijd, atmpd
      real(8) :: eta
c
      NEIG0=n_base_eig_leg
      write(6,*)'@@LSES_SET_EIG_LEG_ATMP'
c
c
      if (n_tot_base .ne. n_base_eig_leg) then
        write(6,*)'ERROR!(LSES_SET_EIG_LEG_ATMP)'
        write(6,*)'  n_tot_base=',n_tot_base
        write(6,*)'  n_base_eig_leg=',n_base_eig_leg
        stop
      endif
c
      eta=5.0
      WRITE(6,*)'    NEIG0, ETA=',neig0, eta
c
c     if (.not. allocated(dsij)) then
c       write(6,*)'ERROR(SETMAT2):DSIJ is not allocated!!'
c       stop
c     else
c       write(6,*)'...DSIJ is already allocated. OK!'
c     endif   
c
c     if (.not. allocated(dpij)) then
c       write(6,*)'ERROR(SETMAT2):DPIJ is not allocated!!'
c       stop
c     else
c       write(6,*)'...DPIJ is already allocated. OK!'
c     endif   
c
      CALL TCLOCK(TIME2)
      TB2=TIME2
c
      ATMP(:,:)=0.0D0
      ATMP2(:,:)=0.0D0
c
c
        DO 20 I=1,NEIG0
           ATMP(I,I)=10.0D0*ETA
   20   CONTINUE
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Mark the terminated dangling bonds
c            IDNGL2 = 1 : dangling bond
c                   = 0 : not dangling bond
c
      DO 31 J=1,NEIG0
        IDNGL2(J)=0
   31 CONTINUE
c
      ISUM=0
      DO 32 JSV=1,NOAV
        JS=JS4JSV(JSV)
      DO 32 JA=1,NVL
        JG=JS2J(JA,JS)
        J=JA+(JSV-1)*NVL
        IF (J .GT. NEIG0) THEN
          WRITE(6,*)'ERROR!(SETEIG:32):J,JSV,JA=',J,JSV,JA
          STOP
        ENDIF
        III=IDNGL(JG)
        IDNGL2(J)=III
        IF (III .NE. 0) THEN
          ISUM=ISUM+1
c         WRITE(6,*)'DANGLING=',J,III
        ENDIF
   32 CONTINUE
      WRITE(6,*)'TOTAL DANGLING=',ISUM
c
c     STOP
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Setting the matrix to be diagonalized.
c
      DO 100 JSV2=1,NOAV
         JS2=JS4JSV(JSV2)
c        NJSD2=NJSD(JSV2,ICT4L)
         NJSD2=NJSD(JSV2,ICT4H)
         NSS2=JSEI(JS2)
         NVAL2=NVAL(NSS2)
       DO 100 JSD1=1,NJSD2
         JSV1=JSV4JSD(JSD1,JSV2)
         JS1=JS4JSV(JSV1)
         NSS1=JSEI(JS1)
         NVAL1=NVAL(NSS1)
        DO 100 JA2=1,NVAL2
          JG=JS2J(JA2,JS2)
          J=JA2+(JSV2-1)*NVL
          IF ((J .LE. 0) .OR. (J .GT. NEIG0)) THEN
            WRITE(6,*)'ERROR!(SETEIG1:100):J,NEIG0=',J,NEIG0
            STOP
          ENDIF
        DO 100 JA1=1,NVAL1
          IG=JS2J(JA1,JS1)
          I=JA1+(JSV1-1)*NVL
          IF ((J .LE. 0) .OR. (J .GT. NEIG0)) THEN
            WRITE(6,*)'ERROR!(SETEIG1:100):J,NEIG0=',J,NEIG0
            STOP
          ENDIF
          DBIGD=DHIJ(JA1,JA2,JSD1,JSV2)
c            ---> H + 2 eta rho_{PT}
c         DSIJD=DSIJ(JA1,JA2,JSD1,JSV2)
c            ---> Overlap
          ATMP(I,J)=DBIGD
c         ATMP2(I,J)=DSIJD
  100 CONTINUE
c
      CALL TCLOCK(TIME2)
      WRITE(6,*)'@@@ SETEIG TIME1=',TIME2-TB2
      TB2=TIME2
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      DO 33 JSV=1,NOAV
        JS=JS4JSV(JSV)
      DO 33 JA=1,NVL
        JG=JS2J(JA,JS)
        J=JA+(JSV-1)*NVL
        JJ=IDNGL2(J)
        ATMPD=ATMP(J,J)
        IF (JJ .EQ. 1) ATMPD=10.0D0*ETA
        ATMP(J,J)=ATMPD
c       WRITE(6,*)'J,ATMP',J,JS,ATMP(J,J)
   33 CONTINUE
c        ---> correction for dangling bond
c
 9999 CONTINUE
      write(6,*)'.. ended w/o error :LSES_SET_EIG_LEG_ATMP'
c
      END
