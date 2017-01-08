!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Eigen-state solver 
!       NKEIG  : Highest occpued level
!       rNelec : Electron number
!         (Non-essential, just for plotting HOMO-LUMO gap)
!       NEIG   : # of basis 
!
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_seteig3(NKEIG,rNelec,imode)
      use MNRLAdaptor
      use elses_mod_phys_const, only : ev4au
      use elses_mod_sim_cell, only : noa
      use elses_mod_md_dat,    only : itemd,itemdorg
!     use elses_mod_eig_leg, only: n_base_eig_leg, fb,tot_elec_eig_leg
      use elses_mod_eig_leg, only: n_base_eig_leg, tot_elec_eig_leg
!     use elses_arr_eig_leg, only:atmp,atmp2,atmp3,atmpo,eig2
      use elses_arr_eig_leg, only:atmp,atmp2,eig2
      use M_la_eigen_solver, only: la_eigen_solver_wrapper !(routine)
!     INCLUDE 'zconst.f'
!     INCLUDE 'zconst-VR.f'
!     INCLUDE 'zconst-EIG.f'
      implicit none
!     real(8), intent(in) :: fbd
      integer, intent(in) :: nkeig
      real(8), intent(in) :: rNelec
      integer, intent(in) :: imode
      integer :: nn, j, i, nmesh, iii, k, ifddist
      integer :: kk
!     real(8) :: ACWRK1(NEIG0),ACWRK2(NEIG0)
!     real(8) :: EIG2(NEIG0),EIG3B(NEIG0)
!      integer :: IDNGL2(NEIG0)
      real(8) :: tb2, eig3, edev, etop, ebtm
      real(8) :: esum, esum2, eave1, eave2, egapd, ddd
      real(8) :: fmu, eneig, efermi, dsum, dsum1, dsum2
      real(8) :: ddd1, etmp, eigk, eigk2, fdcoe, xxx, fdcoe1
      real(8) :: derr
      real(8) :: fdfunc
      real(8) :: erfcc
      real(8), parameter :: LDOS_EMIN = -0.5
      real(8), parameter :: LDOS_EMAX = +1.0
      real(8), parameter :: LDOS_SIGMA = 1d-3
      integer, parameter :: LDOS_DIVIDE = 1501
      integer, parameter :: LDOS_ATOM1 = 70
      integer, parameter :: LDOS_ATOM2 = 61
      real(8), allocatable :: weights1(:)
      real(8), allocatable :: weights2(:)
      integer :: err 
      integer :: neig0, neig
      real(8) :: time2
!
      neig0=n_base_eig_leg
      NEIG=NEIG0
!
!
      WRITE(6,*)'@@SETEIG:NKEIG    =',NKEIG
      WRITE(6,*)'        NEIG,NEIG0=',NEIG,NEIG0
!
      write(6,*)' imode=',imode
      write(6,*)'   = 1 (standard), 2 (generazed)'
      write(6,*)'   = 3 (overlap-mat)'
!
      CALL TCLOCK(TIME2)
      TB2=TIME2
!
      IF (NKEIG .GT. NEIG-1) THEN
        WRITE(6,*)'ERROR!(SETEIG):NEIG,NKEIG=',NEIG,NKEIG
        STOP
      ENDIF
!
!     DO 12 J=1,NEIG0
!     DO 12 I=1,NEIG0
!       ATMP3(I,J)=ATMP(I,J)
!       ATMPO(I,J)=ATMP2(I,J)
!  12 CONTINUE
   13 FORMAT('ATMP=',2I6,2F20.10) 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      call la_eigen_solver_wrapper(imode)
!     call elses_eig_mateig(imode)
!       --> ATMP is replaced by eigen state
!           ATMP2 is also broken by Cholesky decomposition.
! IGUCHI begin 06/11/07
!   LDOS calculation and outputs.
!       atmp:  eigenwavefunctions
!       atmp2: overlap matrix
!       atmp3: hamiltonian matrix
!       eig2:  eigenenergies 
!       neig:  dimension of Hamiltonian matrix 
!
      if ((itemd == 1) .or. (mod(itemd + itemdorg, 25) == 1)) then
        allocate(weights1(neig), &
             weights2(neig), stat=err)
        if(err /= 0) then
          write(0, *) "Memory Allocation Error in SETEIG3."
          stop
        end if
        weights1(:) = 0d0
        weights2(:) = 0d0
        ! All
        weights1(:) = 1.0d0
        ! One Atom
        !weights1(nvl*(LDOS_ATOM1-1)+1:nvl*(LDOS_ATOM1-1)+9) = 1.0d0
        !weights2(nvl*(LDOS_ATOM2-1)+1:nvl*(LDOS_ATOM2-1)+9) = 1.0d0

        if (itemd == 1) then
          open(91, FILE="LDOS", status="UNKNOWN")
        else
          open(91, FILE="LDOS", status="UNKNOWN", position="APPEND")
        end if

        write(91, "('# ITEMD=', I0)") itemd + itemdorg
        write(91, "('# ATOM1=', I0)") LDOS_ATOM1 
!        call NrlWritePartialDos(91, atmp, eig2, atmpo, neig, LDOS_EMIN,
!     &      LDOS_EMAX, LDOS_DIVIDE, LDOS_SIGMA, weights1, err)
!       call NrlWriteLocalDos(91, atmp, eig2, atmpo, neig, LDOS_EMIN, &
!            LDOS_EMAX, LDOS_DIVIDE, LDOS_SIGMA, weights1, err)
!        write(91, "('# ATOM2=', I0)") LDOS_ATOM2
!        call NrlWriteCohp(91, atmp, eig2, atmp3, atmpo, neig, 
!     &      LDOS_EMIN, LDOS_EMAX, LDOS_DIVIDE, LDOS_SIGMA, 
!     &      weights1, weights2, err)
        close(91)
        if (err /= 0) then
          write(0, *) "Fails in calculation of LDOS in SETEIG3."
          stop
        end if
        deallocate(weights1)
        deallocate(weights2)
      end if 
! IGUCHI end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      CALL TCLOCK(TIME2)
      WRITE(6,*)'@@@ STEIG TIME1=',TIME2-TB2
      TB2=TIME2
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      if (imode .eq. 3) then
        write(6,*)'Analysis of overlap mat'
        do j=1,neig
          eig3=eig2(j)
          write(6,*)'eig S=',j,eig3
          eig2(j)=eig3/ev4au
        enddo
        edev=0.02d0
        nmesh=6000
        etop=5.0d0
        ebtm=-1.0D0
!       call dosgau2(eig2,neig,etop,ebtm,edev,nmesh)
        stop
      endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     WRITE(76,*)N
      WRITE(6,*)NEIG
!
!
      III=MOD(ITEMD,10)
!     III=1
      DO 85 J=1,NEIG
        EIG3=EIG2(J)*EV4AU
        IF (III .EQ. 1) WRITE(6,86)ITEMD,J,EIG3,EIG2(J)*2.0d0
   85 CONTINUE
   86 FORMAT('LEV[eV,Ry]=',2I6,2F24.16) 
!
      ESUM=0.0D0
      ESUM2=0.0D0
      DO 80 K=1,NEIG
       IF (K .LE. NKEIG) THEN
          ESUM=ESUM+EIG2(K)
       ELSE
          ESUM2=ESUM2+EIG2(K)
       ENDIF
   80 CONTINUE
      EAVE1=ESUM/DBLE(NKEIG)*EV4AU
      EAVE2=ESUM2/DBLE(NEIG-NKEIG)*EV4AU
!
      EDEV=0.1D0
      NMESH=6000
!     EDEV=0.01D0
!     NMESH=50000
!     EDEV=0.01D0*2.0D0*1.36058D0
!     NMESH=25000
!
      ETOP=30.0D0
      EBTM=-30.0D0
!     IF (ITEMD .EQ. 1) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 26) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 51) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 76) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 101) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 126) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 151) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 176) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     IF (ITEMD .EQ. 201) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     ETOP=80.0D0
!     EBTM=30.0D0
!     IF (ITEMD .EQ. 1) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     STOP
!     ETOP=290.0D0
!     EBTM=240.0D0
!     IF (ITEMD .EQ. 1) CALL DOSGAU2(EIG2,NEIG,ETOP,EBTM,EDEV,NMESH)
!     
!     CALL DOSGAU2(EIG2,NEIG,000.0D0,-100.0D0,EDEV,NMESH)
!     CALL DOSGAU2(EIG2,NEIG,100.0D0,000.0D0,EDEV,NMESH)
!     CALL DOSGAU2(EIG2,NEIG,200.0D0,100.0D0,EDEV,NMESH)
!     CALL DOSGAU2(EIG2,NEIG,300.0D0,200.0D0,EDEV,NMESH)
!     CALL DOSGAU2(EIG2,NEIG,400.0D0,300.0D0,EDEV,NMESH)
!     CALL DOSGAU2(EIG2,NEIG,500.0D0,400.0D0,EDEV,NMESH)
!
      EGAPD=EIG2(NKEIG+1)-EIG2(NKEIG)
      WRITE(6,*)'# of tot. states  = ',NEIG
      WRITE(6,*)' N (probed state) = ',NKEIG
      WRITE(6,*)'   ETB (tot )[au] = ',ESUM*2.0D0
      WRITE(6,*)'   ETB (atom)[eV] = ',ESUM*2.0D0/DBLE(NOA)*EV4AU
      WRITE(6,*)'   (bandgap1)[au] = ',ITEMD,EGAPD
      WRITE(6,*)'    (lowest )[au] = ',EIG2(1)
      WRITE(6,*)'         E(N)[au] = ',EIG2(NKEIG)
      WRITE(6,*)'       E(N+1)[au] = ',EIG2(NKEIG+1)
      WRITE(6,*)'    (highest)[au] = ',EIG2(NEIG)
      WRITE(6,*)'   Eave(occ.)[eV] = ',EAVE1
      WRITE(6,*)'   Eave(unoc)[eV] = ',EAVE2
      WRITE(6,*)'   (bandgap1)[eV] = ',ITEMD,EGAPD*EV4AU
      WRITE(6,*)'    (lowest )[eV] = ',EIG2(1)*EV4AU
      WRITE(6,*)'         E(N)[eV] = ',EIG2(NKEIG)*EV4AU
      WRITE(6,*)'       E(N+1)[eV] = ',EIG2(NKEIG+1)*EV4AU
      WRITE(6,*)'    (highest)[eV] = ',EIG2(NEIG)*EV4AU
      WRITE(6,*)'(homolumo) = ',&
           ITEMD,EIG2(NKEIG)*EV4AU,EIG2(NKEIG+1)*EV4AU
!
!     STOP
!
      CALL TCLOCK(TIME2)
      WRITE(6,*)'@@@ SETEIG TIME3=',TIME2-TB2
      TB2=TIME2
!
!
 9999 CONTINUE
!
      end subroutine elses_seteig3
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ (Generalized) eigen-value problem for real-symmetrix matrix
!           using LAPACK 
!
!           imode=1: Generalized eigen-value problem
!           imode=2: Standard eigen-value problem
!
!           A : Hamiltonian matrix 
!           B : Overlap matrix 
!
!
!    * Be careful for the fact that 
!          the matrice A and B is not preserved on exit.
!   
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_eig_mateig(imode)
!
      use M_io_dst_write_log,  only : log_unit
      use elses_mod_phys_const, only : ev4au
      use elses_mod_eig_leg, only : n_base_eig_leg
      use elses_arr_eig_leg, only : A=>atmp, B=>atmp2, EIG_WRK=>eig2
      !          renaming for minimal change of the code
      implicit none
      real(8), allocatable, save :: WORK(:) !! not thread-safe routine
      integer imode, info_ham
      integer n
!
      REAL*8 DDDA,DDDB,DDDC,EPS
      REAL*8 CCC1,CCC2,CCC3
      real*8 tb2, time2
!
      INTEGER INFO,LDA,LDB,ITYPE,LWORK,I,J,II
!
!     REAL*8 EIG2(N)
!     REAL*8 A(N,N), B(N,N),WORK(3*N-1)
!     REAL*8 A2(N*(N+1)/2),B2(N*(N+1)/2)
!     REAL*8 EPSZ,EPST
!
!
      n=n_base_eig_leg
      write(6,*)'@@ elses_eig_mateig:n=',n
!
      if( .not. allocated(WORK) )then
         allocate (WORK(20*n),stat=info) ! large size of work is preferable
         !allocate (WORK(3*n-1),stat=info) ! this is the minimum-size
         if( info .ne. 0) then
            write(6,*)'allocation error!(MAT_ALLOC):info=',info
            stop
         end if
      end if
!
      if ( .not. allocated(a) ) then
        write(6,*)'ERROR(elses_eig_mateig)'
        write(6,*)'  An array is not allocated:atmp'
        stop
      endif   
!
      if ((n .le. 0) .or. (n .gt. 100000)) then
        write(6,*)'ERROR?!(elses_eig_mateig)'
        write(6,*)'  matrix size : n=',n
        stop
      endif   
!
      EPS=1.0D-6
!
!     !EV4AU=2.0D0*13.6058D0
!      1 a.u = (EV4AU) eV
!
!     IF (N2 .NE. N) THEN
!        WRITE(6,*)'ERROR!(MATEIG3)'
!        WRITE(6,*)'  N2,N =',N2,N
!        STOP
!     ENDIF
!
      INFO=0
      ITYPE=1
      LDA=N
      LDB=N
!
!     EPSZ=0.0D0
!     EPST=0.0D0
!
      DDDA=0.0D0
      DDDB=0.0D0
      DO 3 J=1,N
      DO 3 I=1,N
         DDDA=DDDA+DABS(A(I,J))
         DDDB=DDDB+DABS(B(I,J))
    3 CONTINUE
      WRITE(6,*)' Summation |A_{ij}|, |B_{ij}| =',DDDA,DDDB
      WRITE(6,*)' If summation is infinite,'
      WRITE(6,*)'   please check the array size'
!
      DDDA=0.0D0
      DDDB=0.0D0
      DDDC=0.0D0
      DO 5 J=1,N
      DO 5 I=1,N
         CCC1=A(I,J)-A(J,I)
         CCC2=B(I,J)-B(J,I)
         CCC3=B(I,J)
         DDDA=DDDA+DABS(CCC1)
         DDDB=DDDB+DABS(CCC2)
         IF (I .EQ. J) CCC3=CCC3-1.0D0
         DDDC=DDDC+DABS(CCC3)
!        write(6,*)'i,j,A,B=',i,j,a(i,j),b(i,j)
    5 CONTINUE
       WRITE(6,*)' Non-symmetric A,B =',DDDA,DDDB
       WRITE(6,*)' Non-orthogonalty B =',DDDC
!      stop
!
      IF (DDDB+DDDA .GT. EPS) THEN
        DO 15 J=1,N
        DO 15 I=1,N
          CCC1=A(I,J)-A(J,I) 
          IF (DABS(CCC1) .GT. 1.0D-10) THEN
           WRITE(6,16)I,J,A(I,J)*ev4au,A(J,I)*ev4au
          ENDIF
   15   CONTINUE
        WRITE(6,*)'ERROR!(MATEIG)'
        STOP
      ENDIF
   16 FORMAT('I,J,=',2I10,2E20.8)
!         
!
!     II=0
!     DO 10 J=1,N
!     DO 10 I=1,J
!        II=II+1
!        A2(II)=(A(I,J)+A(J,I))/2.0D0
!        B2(II)=(B(I,J)+B(J,I))/2.0D0
!  10 CONTINUE
!
      CALL TCLOCK(TIME2)
      TB2=TIME2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LWORK=size(WORK,1)

      if( LWORK < 3*N-1 ) stop '*** ERROR : too small work size'

      select case (imode)
         case(1)
            !  Standard Eigen-value equation (A-e)x=0
            !    ---> for LAPACK
            WRITE(6,*)'Go to LAPACK DSYEV(standard eigen)'
            CALL DSYEV("V","U",N,A,LDA,EIG_WRK,WORK,LWORK,INFO)
         case(2)
            !  Generalized Eigen-value equation (A-eB)x=0
            !    ---> for LAPACK
            WRITE(6,*)'Go to LAPACK DSYGV(genelized eigen)'
            CALL DSYGV(ITYPE,"V","U",N,A,LDA,B, &
                 LDB,EIG_WRK,WORK,LWORK,INFO)
            if(INFO /= 0) then
               info_ham=info
               call CalculateEivenvaluesOfOverlap
               info=info_ham
               write(*,*)'INFO:Error in eigen-value solver'
               write(*,*)'INFO:The cutoff may be too small'
            end if
         case(3)
            !  calculate the eigen values of B
            !    ---> for LAPACK
               call CalculateEivenvaluesOfOverlap
         case default
            write(6,'("elses_eig_mateig: error, imode is out of range.")')
         end select

         ! working array size reconfiguration
         write   (*,'("* WORK size LAPACK suggests : ",I8)') nint(WORK(1))
         write   (*,'("* current WORK size         : ",I8)') LWORK
         if(nint(WORK(1)) > LWORK) then
            write(*,'("* performance is degraded due to small WORK size")')
            LWORK=nint(WORK(1))*2
            ! the facor 2 is a work-around for eroneous mkl v10
            deallocate(work)
            allocate(work(LWORK))
            write(*,'("* WORK size is enlarged   to  ",I8)') LWORK
         end if
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call tclock(time2)
      select case (imode)
         case(1,3)
           write(*,'(a,f25.10)') 'TIME for LAPACK (DSYEV) =',TIME2-TB2
           if (log_unit > 0) then
             write(log_unit,'(a,f25.10)') 'TIME for LAPACK (DSYEV) =',TIME2-TB2
           endif   
         case(2)
           write(*,'(a,f25.10)') '  TIME for LAPACK (DSYGV) =',TIME2-TB2
           if (log_unit > 0) then
             write(log_unit,'(a,f25.10)') 'TIME for LAPACK (DSYGV) =',TIME2-TB2
           endif   
         case default
           write(*,'("elses_eig_mateig: error, imode is out of range.")')
           stop
      end select      
      tb2=time2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Settle the sign freedom for eigen vector 
!      Convert v --> -v, if necessary.
!        so as to satisfy sum_i v(i) > 0 
!       or, so as to set eigen vectors uniquely. 
!
      write(*,'(a)') 'Settle the sign freedom for eigen vectors'
      if (log_unit > 0) then
         write(log_unit,'(a)') 'Settle the sign freedom for eigen vectors'
      endif   
!
!$omp  parallel &
!$omp& default(shared)&
!$omp& private(j)
!$omp do schedule(static)
      do j=1,size(a,2)
        if (sum(a(:,j)) < 0.0d0) a(:,j) = -a(:,j)
      enddo   
!$omp end do
!$omp end parallel 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
 9999 CONTINUE
      IF (INFO .NE. 0) THEN
       WRITE(6,*)' MATEIG3:INFO =',INFO
       STOP
      ENDIF
!
      contains
        subroutine CalculateEivenvaluesOfOverlap
          WRITE(6,*)'Go to LAPACK DSYEV(standard eigen)'
          CALL DSYEV("V","U",N,B,LDB,EIG_WRK,WORK,LWORK,INFO)
          WRITE(6,*)'return from LAPACK:info =',INFO
          write(6,'("Eigenvalues of overlap matrix:")')
          write(6,'(8ES10.3)') EIG_WRK
        end subroutine CalculateEivenvaluesOfOverlap

      end subroutine elses_eig_mateig
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set chemical potential and occupation number
!       FBD    : temperature in eV
!
!       IFDDIST=-1 : FermiDirac
!              = 0 : Error fn. 
!
!       NOTE: eig2(1:n) is assumed to be sorted.
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_eig_chem_pot
      use M_qm_domain,         only : i_verbose    !(unchanged)  
      use M_io_dst_write_log,  only : log_unit !(unchanged)
      use elses_mod_phys_const, only : ev4au !(parameter)
      use  elses_mod_md_dat,    only : itemd !(unchanged)
      use elses_mod_eig_leg, only:tot_elec_eig_leg !(unchanged)
      use elses_arr_eig_leg, only:eig2  !(unchanged)
      use elses_arr_eig_leg, only:f_occ !(CHANGED)
      use elses_mod_elec_cond,  only : temp_for_electron !(unchanged)
      use M_qm_domain,  only : chemical_potential !(CHANGED)
      implicit none
!
      integer :: ifddist
      integer n, jj
      integer ierr
!
      integer ipe,npe
      integer omp_get_thread_num
      integer omp_get_num_threads
      real*8  rNelec
      real*8  tb2, time2
!
      integer iloop, nloopmax, iconv
!     integer jsv, js, nss, nval2, ja, nrecl
      real*8  xtemp, xbeta, ddd1, ddd2, ddemin, ddemax
      real*8  xmu0, xmumn, xmumx, xidos0, xidosE, tdossum, tensum
      real*8  rEd, yyy, xexp, xqinit, err
      real*8  err_cri
      integer level_homo, level_lumo
!
      nloopmax=1000
!     n=n_base_eig_leg
      n=size(eig2,1)
      IFDDIST=-1
      rNelec=tot_elec_eig_leg
      err_cri=1.0d-12 ! original err_cri=1.0d-15
      level_homo=max(nint(rNelec)/2, 1)
      level_lumo=min(level_homo+1,   n)
!
      if (i_verbose >= 1) then
        if (log_unit >0) then 
          write(log_unit,*)'@@ELSES_EIG_CHEM_POT:N=',n
          write(log_unit,*)'   rNelec =',rNelec
          write(log_unit,*)'   err_cri=',err_cri
        endif  
      endif  
!
      CALL TCLOCK(TIME2)
      TB2=TIME2
!
      if ( n < level_homo ) then
        write(*,*) 'ERROR(elses_eig_chem_pot):n=',n
        stop
      endif
!
      if (.not. allocated(f_occ)) then
        allocate(f_occ(n), stat=ierr) 
        if (ierr /= 0) then
          write(*,*) 'Alloc. ERROR(elses_eig_chem_pot):f_occ'
          stop
        endif
      else
        if (size(f_occ,1) /= n) then
          write(*,*) 'ERROR(elses_eig_chem_pot):size(f_occ,1)=', size(f_occ,1)
          stop
        endif
      endif
!
!     if (fb .le. 1.0d-10) then
!       write(6,*) 'zero temperature !'
!       stop
!     endif
!
!     xtemp=fb/ev4au
      xtemp=temp_for_electron
      xbeta=1.0d0/xtemp
      if (i_verbose >= 1) then
        if (log_unit >0) then 
          write(log_unit,*) ' Temp [au,eV] =',xtemp,xtemp*ev4au
          write(log_unit,*) '   type(-1:FermiDi, 0:ErrorFn) =',IFDDIST
          write(log_unit,*) ' eig2(1)[au,eV]=',eig2(1),eig2(1)*ev4au
          write(log_unit,*) ' eig2(N)[au,eV]=',eig2(n),eig2(n)*ev4au
          write(log_unit,*) ' homo and lumo levels =',level_homo, level_lumo
          write(log_unit,'(a,i10,3f20.10)') ' homo-lumo-gap[eV]=', & 
&                            itemd, eig2(level_homo)*ev4au, eig2(level_lumo)*ev4au, &
&                                  (eig2(level_lumo)-eig2(level_homo))*ev4au
          write(log_unit,'(a,i10,3f20.10)') ' homo-lumo-gap[au]=', &
&                            itemd, eig2(level_homo),       eig2(level_lumo), & 
&                                   eig2(level_lumo)-eig2(level_homo)
        endif   
      endif   
      ddemin=eig2(1)
      ddemax=eig2(n)
!
!
      if (ddemax-ddemin .le. 1.0d-10) then
        write(*,*) 'ERROR(CHEM_POT)'
        write(*,*) '  ddemin=',ddemin
        write(*,*) '  ddemax=',ddemax
        stop
      endif   
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ The bisection method
!
      xmu0=ddemin
      xmumn=ddemin
      xmumx=ddemax
!
      iconv=0
      do 1000 iloop=1,nloopmax
        xidos0= 0.0d0
        xidosE= 0.0d0
        tdossum=0.0d0
        tensum =0.0d0
!
!      OMP-NOTE : Shared scalors : xmu0,xbeta
!$omp  parallel &
!$omp& default(shared) &
!$omp& private(ipe,npe,jj)&
!$omp& private(rEd,yyy,xexp)&
!$omp& reduction(+ : xidos0)&
!$omp& reduction(+ : xidosE)&
!$omp& reduction(+ : tdossum)&
!$omp& reduction(+ : tensum)
       ipe=0
       npe=0
!      ipe=omp_get_thread_num()+1
!      npe=omp_get_num_threads()
!      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
       do jj=1,n
         rEd=eig2(jj)
         yyy=xbeta*(rEd-xmu0)
         if(yyy .gt. 100.d0) then
           xexp=0.0d0
         elseif(yyy.lt.-100.d0) then
           xexp=1.0d0
         else
           xexp=1.0d0/(dexp(yyy)+1.0d0)
         endif
         tdossum=tdossum +2.0d0
         tensum =tensum  +2.0d0*rEd
         xidos0=xidos0   +2.0d0*xexp
         xidosE=xidosE   +2.0d0*xexp*rEd
         if (iconv .eq. 1) f_occ(jj)=xexp
       enddo
!$omp end do
!$omp end parallel 
!
!
        xqinit=rNelec
        if (dabs(xqinit) .le. 1.0d-10) then
          write(*,*)'ERROR!(CHEM_POT):xqinit=',xqinit
          stop
        endif   
        err=(xqinit-xidos0)/xqinit
        if (i_verbose >= 1) then
          if (log_unit > 0) then 
            if (iloop .le. 3) write(log_unit,*)'Bisec.',iloop,xmu0,err
          endif
        endif  
        if(abs(err) .le. err_cri) then
          if(iconv .eq. 1) then 
            if (i_verbose >= 1) then
              if (log_unit > 0)  write(log_unit,1092)iloop,xmu0*ev4au,xidos0
            endif 
            goto 1001
          else
            iconv=1
          endif
        endif   
!
        if(err .gt.  0.0d0) then
          xmumn=xmu0 
          xmumx=xmumx
          xmu0=(xmumn+xmumx)*0.5d0
        else
          xmumn=xmumn
          xmumx=xmu0 
          xmu0=(xmumn+xmumx)*0.5d0
        endif
!       
 1000  continue
       write(*,*)'ERROR(elses_eig_chem_pot)'
       write(*,*)' Chemical potential is NOT'
       write(*,*)'   converged in the bisection loop'
       write(*,*)' Relative error (d N / N ) =',err
       write(*,*)' Criteria : err_cri=',err_cri
       write(*,*)'--> The criteria may be too severe'
       stop
 1001  continue
       if (log_unit > 0) then
         write(log_unit, 1091) itemd,iloop,xmu0,xmu0*ev4au
       endif  
 1091  format('chem_pot [au,eV]=',2I10,2F30.15)
       chemical_potential=xmu0
!
 1092  format('Bisec. Conv:',I10,2F30.15)
!      stop
!ccccccccccccccccccccccccccccccccccccccccccccccccc
!    Confirm the elec. num. 
!    
       xidos0=0.0d0
       do jj=1,n
         xidos0=xidos0+2.0d0*f_occ(jj)
       enddo   
       if (i_verbose >= 1) then
         if (log_unit > 0) write(log_unit,*)'Elec. num=',xidos0, rNelec
       endif  
       if (dabs(xidos0-rNelec)/rNelec .ge. 1.0d-8) then
         write(6,*)'ERROR(CHEM_POT):xidos0=',xidos0
         stop
       endif   
!
!
      end subroutine elses_eig_chem_pot
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   @@ Generate L-matrix with fractional occupation
!        with diagonalization method
!
!      BIJ(i,j) = \SUM_{k} C_{i,j} C_{j,k} f_k
!      PIJ(i,j) = \SUM_{k} C_{i,j} C_{j,k} f_k e_k
!
!            ATMP(i,k) is used as eigen states
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_eig_set_dens_mat
!
      use elses_mod_ctrl,     only : i_verbose
      use elses_mod_eig_leg, only: n_base_eig_leg, tot_elec_eig_leg 
      use elses_mod_noav,    only: noav,noas
      use elses_mod_orb1,    only: nvl
      use elses_mod_js4jsv,  only: js4jsv
      use elses_mod_jsv4jsd, only: jsv4jsd,njsd
      use elses_mod_multi,   only: ict4h
      use elses_mod_orb2,    only :js2j
      use elses_arr_dbij,    only: dbij
      use elses_arr_dpij,    only: dpij
      use elses_arr_dsij,    only: dsij
      use elses_arr_eig_leg, only: atmp, eig2, f_occ
!
      implicit none
      integer :: ipe,npe, neig_k, iii, j, k, ierr
      integer :: jsv2, js2, ja2, jsd1, jsv1, js1, ja1
      integer :: i
      integer :: neig0, neig
!      integer  omp_get_thread_num
!      integer  omp_get_num_threads
      real*8, allocatable :: atmp4(:,:)
      real*8, allocatable :: atmp5(:,:)
      real*8, allocatable :: atmp6(:,:)
      real*8, allocatable :: atmp7(:,:)
      real*8, allocatable :: atmp8(:,:)
      real(8) :: time8, tb8, ddcri, ddd, ddsum, rNelec
      real(8) :: fdcoe, ddd1, ddd2, dpijd, dsijd
!
!
      if (i_verbose >= 1) then
        write(6,*)'@@@ LSES_EIG_SET_DENS_MAT'
        write(6,*)'  DBIJ and DPIJ are calculated !!'
      endif  
!
!
      call tclock(time8)
      tb8=time8
      neig0=n_base_eig_leg
      neig=neig0
      neig_k=neig
      rNelec=tot_elec_eig_leg
!
      if (i_verbose >= 1) then
        write(6,*)'noav=',noav
        write(6,*)'noas=',noas
        write(6,*)' nvl=',nvl
        write(6,*)'neig=',neig
      endif  
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   @@ Trivial Checking
!
      if (rNelec .le. 1.0d-8) then
        write(6,*)'ERROR(LSES_EIG_SET_DENS_MAT)'
        write(6,*)' rNelec=',rNelec
        stop
      endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   @@ Set the upper limit of the summation (non-essential)
!
      ddcri=100.0d0
!          ---> dummy setting
!
      if (ddcri .lt. 1.0d0) then
        if (i_verbose >= 1) then
          write(6,*)'set neig_k:criteria=',ddcri
        endif  
        iii=neig
        do k=1,neig
          ddd=f_occ(k)
          if (ddd .lt. ddcri) iii=k
        enddo  
        neig_k=iii
      endif  
!
      if (i_verbose >= 1) then
        write(6,*)' ---> neig_k=',neig_k
      endif  
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   @@ Allocation 
!
      allocate (atmp4(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(atmp4):ierr=',ierr
        stop
      endif
      atmp4(:,:)=0.0d0
!
      allocate (atmp5(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(atmp5):ierr=',ierr
        stop
      endif
      atmp5(:,:)=0.0d0
!
      allocate (atmp6(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(atmp6):ierr=',ierr
        stop
      endif
      atmp6(:,:)=0.0d0
!
      allocate (atmp7(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(atmp7):ierr=',ierr
        stop
      endif
      atmp7(:,:)=0.0d0
!
      allocate (atmp8(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(atmp7):ierr=',ierr
        stop
      endif
      atmp8(:,:)=0.0d0
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Checking 
!
      IF (NEIG .NE. NEIG0) THEN
        WRITE(6,*)'ERROR!(CALLIJ6):NEIG,NEIG0=',NEIG,NEIG0
        STOP
      ENDIF
!
      if (.not. allocated(dpij)) then
        write(6,*)'ERROR(CALLIJ6):DPIJ is not allocated!!'
        stop
      else
        if (i_verbose >= 1) then
          write(6,*)'...DPIJ is already allocated. OK!'
        endif  
      endif   
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Initial clear
!
      dbij(:,:,:,:)=0.0d0
      dpij(:,:,:,:)=0.0d0
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ6: time(10) =',time8-tb8
      endif  
      tb8=time8
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Check the elec. number (optional)
!
      ddsum=0.0d0
      do k=1,neig
        fdcoe=f_occ(k)
        ddsum=ddsum+fdcoe*2.0d0
      enddo  
      if (i_verbose >= 1) then
        WRITE(6,*)'Elec num. (from occ. num)=',ddsum
      endif  
      if (dabs(ddsum-rNelec)/rNelec .gt. 1.0d-8) then
        write(6,*)'ERROR(LSES_EIG_SET_DENS_MAT'
        write(6,*)'   ddsum=',ddsum
        write(6,*)'  rNelec=',rNelec
        stop
      endif    
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ6: time(20) =',time8-tb8
      endif  
      tb8=time8
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   @@ Calculation ATMP7, ATMP8
!
      if (i_verbose >= 1) then
        write(6,*)'Calculation of ATMP7, ATMP8'
      endif  
!
      atmp6=transpose(atmp)
!         ---> traspose of C(j,k) 
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ6: time(22) =',time8-tb8
      endif  
      tb8=time8
!
      do j=1,neig
         atmp7(:,j)=atmp6(:,j)*f_occ(:)
         atmp8(:,j)=atmp6(:,j)*eig2(:)*f_occ(:)
      enddo   
!       ATMP7: D(k,j) = f(k) C(j,k)
!       ATMP8: Q(k,j) = f(k) E(k) C(j,k)
!
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ6: time(30) =',time8-tb8
      endif  
      tb8=time8
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calc. of L-matrix and P-matrix as dbij and dpij
!
!      L(i,j) = sum_k C(i,k) f(k) C(j,k)
!             = sum_k C(i,k) D(k,j) 
!
!      P(i,j) = sum_k C(i,k) f(k) E(k) C(j,k)
!             = sum_k C(i,k) Q(k,j) 
!
      if (i_verbose >= 1) then
        write(6,*)'calc.of L-matrix and P-matrix'
      endif  
!
!
!$omp  parallel &
!$omp& default(shared)&
!$omp& private(ipe,npe,ierr)&
!$omp& private(jsv2,js2,ja2,j)&
!$omp& private(jsd1,jsv1,js1,ja1,i)&
!$omp& private(ddd1,ddd2)
      ipe=0
      npe=0
!     ipe=omp_get_thread_num()+1
!     npe=omp_get_num_threads()
!     write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
      do jsv2=1,noav
         js2=js4jsv(jsv2)
        do ja2=1,nvl
           j=js2j(ja2,js2)
           do jsd1=1,njsd(jsv2,ict4h)
              jsv1=jsv4jsd(jsd1,jsv2)
              js1=js4jsv(jsv1)
              do ja1=1,nvl
                i=js2j(ja1,js1)
                ddd1=dot_product(atmp6(1:neig_k,i),atmp7(1:neig_k,j))
                ddd2=dot_product(atmp6(1:neig_k,i),atmp8(1:neig_k,j))
                dbij(ja1,ja2,jsd1,js2)=ddd1
                dpij(ja1,ja2,jsd1,js2)=ddd2
              enddo   
           enddo  
        enddo  
      enddo  
!$omp end do
!$omp end parallel 
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ5: time(40) =',time8-tb8
      endif  
      tb8=time8
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Deallocation
!
      if (i_verbose >= 1) then
        write(6,*)'Deallc of ATMP4, ATMP5 and ATMP6'
      endif  
!
      deallocate(atmp4,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'deallc error!(atmp4):ierr=',ierr
        stop
      endif
!
      deallocate(atmp5,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'deallc error!(atmp5):ierr=',ierr
        stop
      endif
!
      deallocate(atmp6,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'deallc error!(atmp6):ierr=',ierr
        stop
      endif
!
      deallocate(atmp7,stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'deallc error!(atmp7):ierr=',ierr
        stop
      endif
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ5: time(50) =',time8-tb8
      endif  
      tb8=time8
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Check that Tr[PS] = E (optional)
!
!     goto 9999
!  
      if (i_verbose >= 1) then
        write(6,*)' Check that Tr[PS] = E'
      endif  
!
      ddsum=0.0d0  
      do jsv2=1,noav
           js2=js4jsv(jsv2)
        do jsd1=1,njsd(jsv2,ict4h)
           jsv1=jsv4jsd(jsd1,jsv2)
           js1=js4jsv(jsv1)
           do ja2=1,nvl
             do ja1=1,nvl
               dpijd=dpij(ja1,ja2,jsd1,js2)
               dsijd=dsij(ja1,ja2,jsd1,js2)
               ddsum=ddsum+dpijd*dsijd
             enddo  
           enddo  
        enddo  
      enddo  
      ddsum=ddsum*2.0d0
!
      if (i_verbose >= 1) then
        write(6,*)'ETB as Tr[PS] =',ddsum
      endif  
!
      call tclock(time8)
      if (i_verbose >= 1) then
        write(6,*)'CALLIJ5: time(60) =',time8-tb8
      endif  
      tb8=time8
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 9999 continue
!
      end subroutine elses_eig_set_dens_mat
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

