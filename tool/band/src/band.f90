!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
program band_calc
  use diag
  use m_plot_data, only : plot_and_mask_matrix, plot_eigen_value
  implicit none

  ! ===== Variables for command line options
  character(len=32) :: argv
  integer, parameter :: ID_MODE_LINE = 1
  integer, parameter :: ID_MODE_GRID = 2
  integer :: id_mode = ID_MODE_LINE
  logical :: BinaryIO=.false.

  ! ===== Variables for Ecc calculations
  real(DOUBLE_PRECISION), allocatable :: EccPerAtom(:)
  real(DOUBLE_PRECISION):: Ecc
  character(len=*),parameter :: FileNameEcc="EccPerAtom.txt"
  logical :: EccEnabled

  real(DOUBLE_PRECISION), allocatable :: S(:,:,:,:), H(:,:,:,:), densS(:,:), densH(:,:)
  integer, allocatable  :: jsv4jsd(:,:), njsd(:,:), nval(:)
  integer, parameter    :: IUNIT=10, JUNIT=11, JUNIT2=12
  integer, parameter    :: ict=1
  integer, parameter    :: IWT=0
  real(DOUBLE_PRECISION):: pi,atominUnit,aLat
  integer :: jsv1, jsv2, njsd2, nval1, nval2
  integer :: jsd1, j1, j2, m, i
  integer :: nMaxOrb, nMaxNeigh, noav, nt, int1,int2,int3,int4
  integer :: nOrbs
  integer :: k,ICHECK
  integer, allocatable :: sumNval(:),sumNvalInUnit(:)
  real(DOUBLE_PRECISION) :: double1, double2
  real(DOUBLE_PRECISION), allocatable :: eigval(:)
  character(len=18) :: string
  character(len=*),parameter :: fileNameHamOvl = "OverlapAndHamiltonian.txt"
  ! ===== Variables read from "ForBandCalc.txt"
  character(len=*),parameter :: fileNameParms  = "ForBandCalc.txt"
  integer                :: numPosInUnit, numUnitcells, whereIsOrigin
  real(DOUBLE_PRECISION),allocatable :: posInUnit(:,:), posUnitOrg(:,:)
  ! ===== Inverse k-Space
  integer                :: itotpt
  real(DOUBLE_PRECISION) :: transVec(3,3), volume
  real(DOUBLE_PRECISION) :: recipVec(3,3)
  real(DOUBLE_PRECISION),allocatable :: SLine(:,:),ELine(:,:)
  real(DOUBLE_PRECISION),allocatable :: StartLine(:,:),EndLine(:,:)
  ! ===== Symmetry Lines
  character(len=*),parameter :: fileNameSymLine = "SymLine.txt"
  integer,allocatable :: NumPts(:)
  integer             :: NumLine,iline,nline,ipoint
  ! ===== Grid points
  character(len=*),parameter :: fileNameGridPoint = "Grid.txt"
  integer             :: ika, ikb, ikc
  integer             :: Nka, Nkb, Nkc
  ! ===== Solving Hk-ESk=0 on SymmetricLines 
  real(DOUBLE_PRECISION)    :: ak(3), ak0(3), kdis, kdisSum
  real(DOUBLE_PRECISION),allocatable::eigenVal(:)
  integer                   :: iDimHam 
  complex(DOUBLE_PRECISION), allocatable :: Sk(:,:),S2K(:,:)
  complex(DOUBLE_PRECISION), allocatable, target :: Hk(:,:)
  complex(DOUBLE_PRECISION), pointer     :: Ck(:,:)
  real(DOUBLE_PRECISION), allocatable :: SCk(:,:)
  complex(DOUBLE_PRECISION) :: aim
  ! =====
  integer                   :: atom2,atom22,atom1,atom11,Unit1
  ! ===== output file 
  character(len=*),parameter :: fileNameOutput= "EigenEnergy.txt"
  character(len=*),parameter :: fileNameOutput2= "EigenWave.txt"
  !
  ! =========================================================================== !
  ! Names of variables should be refered in the Manual of ELSES-package         !
  ! nMaxOrb         : #orbitals                                                 !
  ! nt              : internal parameter for checking the input data            !
  ! nMaxNeib        : #neighboring atoms                                        !
  ! noav            : #atoms (Number Of Atoms as the Variables in QM part)      ! 
  ! Correspondence of arguments : S(l1,l2,n1,n2) and H(l1,l2,n1,n2)             !
  ! njsd(I,L)       : #booking atoms for I-th atom with L-th booking radius     !
  ! ict             : =1                                                        !
  ! nval            : #valence orbitals                                         !
  ! nOrbs=sum(nval) : #total valence orbitals = size of dense S,H matrices      !
  ! =========================================================================== !
  !

  
  if( command_argument_count() > 0 ) then
     do j1=1, command_argument_count()
        call get_command_argument(j1,argv)
        !
        if( argv(1:1) == "-" ) then
           if(index(argv(2:),"Bin")>0) then
              BinaryIO=.true.
           else if( index(argv(2:),"line")>0 ) then
              id_mode = ID_MODE_LINE
           else if( index(argv(2:),"grid")>0 ) then
              id_mode = ID_MODE_GRID
           else
              write(*,*) "unknown argument:", trim(argv)
              write(*,*) "possible argument : -line , -grid or -Bin"
              stop
           end if
        end if
     end do
  end if

  pi = 4.d0*atan(1.d0)
  aim=cmplx(0.0d0,1.0d0,DOUBLE_PRECISION)
  atominUnit=27.211d0

  !
  !---- Read Hamiltonian and Overlap matrices (sparse form)----------------
  if(.not. BinaryIO)then
     !! ASCII(ORIGINAL) mode
     open(IUNIT, file=fileNameHamOvl, status='old')
     read(IUNIT,*)
     read(IUNIT,*) nMaxOrb, nt, nMaxNeigh, noav
     if(nt /= nMaxOrb) stop 'format error: int1 /= int2'
     allocate(S(nMaxOrb,nMaxOrb,nMaxNeigh,noav))
     allocate(H(nMaxOrb,nMaxOrb,nMaxNeigh,noav))
     allocate(jsv4jsd(nMaxNeigh,noav))
     allocate(njsd(noav,ict))
     allocate(nval(noav))
     ! ----
     read(IUNIT,*)
     do jsv2=1, noav
        read(IUNIT,'(I8,A18,I8)') nt, string, njsd2
        if(nt /= jsv2) stop 'format error: jsv2'
        njsd(jsv2,ict)=njsd2
        read(IUNIT,*)  (jsv4jsd(jsd1,jsv2) ,jsd1=1, njsd2)!# of neighboring atoms;jsv2=central atom 
     end do
     read(IUNIT,*)
     read(IUNIT,*)  (nval(jsv2),jsv2=1,noav)             !total number of orbitals of each atom 
     read(IUNIT,*)
     do jsv2=1, noav
        njsd2=njsd(jsv2,ict)
        nval2=nval(jsv2)
        do jsd1=1, njsd2
           jsv1=jsv4jsd(jsd1,jsv2)
           nval1=nval(jsv1)
           do j2=1, nval2
              do j1=1, nval1
                 read(IUNIT,*) int1,int2,int3,int4, double1, double2
                 if(int1 /= j1)   stop 'format error: int1'
                 if(int2 /= j2)   stop 'format error: int2'
                 if(int3 /= jsd1) stop 'format error: int3'
                 if(int4 /= jsv2) stop 'format error: int4'
                 S(j1,j2,jsd1,jsv2)=double1
                 H(j1,j2,jsd1,jsv2)=double2
              end do
           end do
        end do
     end do
     close(IUNIT)
     inquire(FILE=fileNameEcc,EXIST=EccEnabled)
     if(EccEnabled) then
        open(IUNIT, file=fileNameEcc, status='old')
        read(IUNIT,*)
        read(IUNIT,*) m
        if(m /= noav) stop 'size mismatch # of atom /= noav'
        allocate(EccPerAtom(noav))
        read(IUNIT,*)
        do m=1, noav
           read(IUNIT,*) EccPerAtom(m)
        end do
        close(IUNIT)
     else
        write(*,*) 'Skip reading Ecc due to backward compatibility'
     end if
  else
     !! Binary(fast) mode
     open(IUNIT, file=fileNameHamOvl, status='old',FORM='UNFORMATTED')
     read(IUNIT) nMaxOrb, nt, nMaxNeigh, noav
     if(nt /= nMaxOrb) stop 'format error: int1 /= int2'
     allocate(njsd(noav,ict))
     allocate(jsv4jsd(nMaxNeigh,noav))
     allocate(nval(noav))
     allocate(S(nMaxOrb,nMaxOrb,nMaxNeigh,noav))
     allocate(H(nMaxOrb,nMaxOrb,nMaxNeigh,noav))
     read(IUNIT) njsd(:,ict) !!! The range of 2nd argument of njsd : 0:nncut !!!
     read(IUNIT) jsv4jsd
     read(IUNIT) (nval(j1),j1=1,noav)
     read(IUNIT) S
     read(IUNIT) H
     close(IUNIT)
     inquire(FILE=fileNameEcc,EXIST=EccEnabled)
     if(EccEnabled) then
        open(IUNIT, file=fileNameEcc, status='old',FORM='UNFORMATTED')
        read(IUNIT) m
        if(m /= noav) stop 'size mismatch # of atom /= noav'
        allocate(EccPerAtom(noav))
        read(IUNIT) EccPerAtom(:)
        close(IUNIT)
     end if
  end if

  !

  !----- read some parameters from "ForBandCalc.txt" ----------------------
  open(IUNIT,file=fileNameParms)
  read(IUNIT,*)                         ! skip comment
  read(IUNIT,*) transVec
  read(IUNIT,*)                         ! skip comment
  read(IUNIT,*) numPosInUnit
  read(IUNIT,*)                         ! skip comment
  allocate(posInUnit(3,numPosInUnit))
  read(IUNIT,*) posInUnit
  read(IUNIT,*)                         ! skip comment
  read(IUNIT,*) numUnitcells
  read(IUNIT,*)                         ! skip comment
  read(IUNIT,*)                         ! skip comment
  allocate(posUnitOrg(3,numUnitcells))
  read(IUNIT,*) posUnitOrg
  read(IUNIT,*)                         ! skip comment
  read(IUNIT,*)                         ! skip comment
  read(IUNIT,*) whereIsOrigin
  close(IUNIT) 
  !
  aLat=sqrt(transVec(1,1)**2+transVec(2,1)**2+transVec(3,1)**2)

  call VecPro(transVec,recipVec)
  volume = dot_product(transVec(:,1),recipVec(:,1))
  recipVec(:,1) = recipVec(:,1) * (2.0d0*pi/volume)
  recipVec(:,2) = recipVec(:,2) * (2.0d0*pi/volume)
  recipVec(:,3) = recipVec(:,3) * (2.0d0*pi/volume)

  write(6,*) 'primitive lattive vectors : '
  write(6,*) '  a1=',transVec(:,1)
  write(6,*) '  a2=',transVec(:,2)
  write(6,*) '  a3=',transVec(:,3)
  write(6,*) 'aLat=|a_1|=', aLat
  write(6,*) 'whereIsOrigin=',whereIsOrigin
  write(6,*) 'posUnitOrg(:,whereIsOrigin)='
  write(6,*)  posUnitOrg(:,whereIsOrigin)
  !
  !------ transform Hamiltonian and Overlap matrix int dens form -------
  !------    from "internal proc. sparse2dens"
  allocate(sumNval(noav))
  k=0
  do jsv1=1, noav
     sumNval(jsv1)=k                 ! sequencial number 
     k=k+nval(jsv1)
  end do
  !------ Hamiltonian and Overlap matrix ; densH and densS (symmetric)
  nOrbs=sum(nval)
  allocate(densS(nOrbs,nOrbs))
  allocate(densH(nOrbs,nOrbs))
  call sparse2dens(S,densS,ICHECK)
  if(ICHECK.ne.0) write(6,*) 'densS is asymmetric.'
  call sparse2dens(H,densH,ICHECK)
  if(ICHECK.ne.0) write(6,*) 'densM is asymmetric.'
  ! diagonalization: eigen values and vectors
  allocate(eigval(nOrbs))
  !   write(*,'(18ES8.1)') ((densS(j1,j2),j2=1,18),j1=1,18) ! Si-Unit cell
  !   write(*,'(10ES8.1)') ((densS(j1,j2),j2=1,10),j1=1,10)
  !   call diag_sym(nOrbs,densS,eigval)
  !   write(*,'(8ES10.3)') eigval

  select case( id_mode )
  case( ID_MODE_LINE )
     !------ Definition of Symmetric Lines --------------------------------
     open(IUNIT,file=fileNameSymLine)
     read(IUNIT,*) NumLine
     allocate(NumPts(NumLine))
     allocate(SLine(3,NumLine))
     allocate(ELine(3,NumLine))
     allocate(StartLine(3,NumLine))
     allocate(EndLine(3,NumLine))
     write(6,*) 'SLine and Eline should be given in a unit of 2*pi/aLat' 
     do nline=1,NumLine
        read(IUNIT,*) NumPts(nline),(SLine(iline,nline),iline=1,3),&
             (ELine(iline,nline),iline=1,3)
        StartLine(:,nline)=SLine(:,nline)*2.d0*pi/aLat
        EndLine(:,nline)  =ELine(:,nline)*2.d0*pi/aLat
     end do
     close(IUNIT)

  case( ID_MODE_GRID )
     !------ Definition of Symmetric Lines --------------------------------
     open(IUNIT,file=fileNameGridPoint)
     read(IUNIT,*) Nka, Nkb, Nkc
     close(IUNIT)

  end select
  !
  !------ Solve Hk-ESk=0 on SymmetricLines ----------------------------- 
!!! ---  Origin-Cell and # of orbitals
  write(6,*) 'Origin-Cell and # of orbitals'
  allocate(sumNvalInUnit(numPosInUnit))
  k=0
  Unit1=whereIsOrigin
  write(6,*) '  numPosInUnit',numPosInUnit
  ! S.Y Sep27, 2013, added in order to claculate Ecc.
  Ecc=0d0
  do atom1=1,numPosInUnit
     sumNvalInUnit(atom1)=k
     atom11=atom1+numPosInUnit*(Unit1-1)
     k=k+nval(atom11)
     write(6,*) '  sumNvalInUnit(atom1)',sumNvalInUnit(atom1)
     write(6,*) '  atom11', atom11
     write(6,*) '  k+nval(atom11)',k
     Ecc = Ecc + EccPerAtom(atom11)
  end do
  if(EccEnabled) write(6,*) '  Ecc/UnitCell = ', Ecc
  iDimHam=k
!!! ---
  ! ---- Set Hamiltonian and Overlap Matrices in k-representation ----------
  ! write(6,*) '#### start- eigenvalue calculation ####'  
  allocate(Hk(iDimHam,iDimHam))
  allocate(Sk(iDimHam,iDimHam))
  allocate(S2k(iDimHam,iDimHam))
  allocate(SCk(iDimHam,iDimHam))
  allocate(eigenVal(iDimHam))


  select case( id_mode )
  case( ID_MODE_LINE )
     open(JUNIT,file=fileNameOutput)

     kdisSum=0.0d0
     write(6,*) 'eigen-value in eV'
     write(6,*)'    ak-x      ak-y      ak-z       Delta-ak          eigen-values ',iDimHam
     itotpt=0
     do nline=1,NumLine 
        ! write(6,*) 'nline=', nline, '      StartLine, EndLine = :'
        ! write(6,'(2X,3E12.5,2X,3E12.5)')(StartLine(k,nline),k=1,3),(EndLine(k,nline),k=1,3)
        do ipoint=1, NumPts(nline)
           itotpt=itotpt+1
           Hk(:,:)=0.0d0
           Sk(:,:)=0.0d0
           ak(:)=StartLine(:,nline) &
                +(EndLine(:,nline)-StartLine(:,nline)) &
                *(ipoint-1)/(NumPts(nline)-1.d0)
           if(ipoint.eq.1) ak0(:)=ak(:)
           kdis=dsqrt((ak(1)-ak0(1))**2+(ak(2)-ak0(2))**2+(ak(3)-ak0(3))**2)
           ak0(:)=ak(:)
           kdisSum=kdisSum+kdis
           ! write(6,*) '  ak=', ak(:)
           ! write(6,*) 'numPosInUnit=',numPosInUnit
           call HkAndSkSetting(ak)
           if(IWT.eq.1) then
              write(6,*) 'Hk='
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(1,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(2,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(3,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(4,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(5,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(6,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(7,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Hk(8,:)
              write(6,*) 'Sk='
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(1,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(2,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(3,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(4,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(5,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(6,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(7,:)
              write(6,'(2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3,1X,2E10.3)') Sk(8,:)
           end if
           ! ---- Solve the eigenvalues ----------------------------------------------
           ! write(6,*) 'iDimHam=',iDimHam
!          call plot_and_mask_matrix(iDimHam, Hk, Sk)
           call  diag_gene_herm(iDimHam,Hk,Sk,eigenVal) ! doesn't preserve Hk & S2k
!          call plot_eigen_value(iDimHam, eigenVal)
           eigenVal(:)=eigenVal(:)*atominUnit
           write(*,1000)itotpt,ak(:),kdisSum,(eigenVal(j1),j1=1,min(iDimHam,20))
           write(JUNIT,1000)itotpt,ak(:),kdisSum,(eigenVal(j1),j1=1,min(iDimHam,20))
        end do
     end do

  case( ID_MODE_GRID )

     open(JUNIT,file=fileNameOutput)
     open(JUNIT2,file=fileNameOutput2)

     write(JUNIT,'(4I4)') Nka, Nkb, Nkc, min(iDimHam,20)
     write(JUNIT2,'(5I4)') Nka, Nkb, Nkc, min(iDimHam,20), iDimHam

     write(6,*) 'eigen-value in eV'
     write(6,*)' grid indeces,   ak-x      ak-y      ak-z     eigen-values ',iDimHam
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              Hk(:,:)=0.0d0
              Sk(:,:)=0.0d0

              ak(:) = recipVec(:,1)*(dble(ika-1)/(Nka-1)) &
              +       recipVec(:,2)*(dble(ikb-1)/(Nkb-1)) &
              +       recipVec(:,3)*(dble(ikc-1)/(Nkc-1))
!!$              ak(:) = recipVec(:,1)*(dble(ika-1)/Nka+0.5d0/Nka) &
!!$              +       recipVec(:,2)*(dble(ikb-1)/Nkb+0.5d0/Nkb) &
!!$              +       recipVec(:,3)*(dble(ikc-1)/Nkc+0.5d0/Nkc)

              call HkAndSkSetting(ak)
              ! ---- Solve the eigenvalues ----------------------------------------------
              S2k=Sk
              call  diag_gene_herm(iDimHam,Hk,S2k,eigenVal) ! doesn't preserve Hk & S2k
              Ck => Hk
              eigenVal(:)=eigenVal(:)*atominUnit
              write(*,    1001) ika, ikb, ikc, ak(:), (eigenVal(j1),j1=1,min(iDimHam,20))
              write(JUNIT,1001) ika, ikb, ikc, ak(:), (eigenVal(j1),j1=1,min(iDimHam,20))

              do m=1, min(iDimHam,20)
                 write(JUNIT2,'(3I4,1E12.3$)') ika,ikb,ikc,eigenVal(m)

                 do i=1, iDimHam
                    SCk(i,m) = cdabs(sum(Sk(i,:)*Ck(:,m)))**2
                 end do
                 SCk(:,m) = SCk(:,m) / sum(SCk(:,m))

                 do i=1, iDimHam
                    write(JUNIT2,'(1F12.6$)') SCk(i,m)
                 end do
                 write(JUNIT2,*)
              end do

           end do ! ikc
        end do ! ikb
     end do ! ika

  end select

  close(JUNIT)
  close(JUNIT2)

1000 format(I4,3E10.3,3X,E12.5,5X,20ES14.5)
!1000 format(1I4,3E12.3,3X,E14.5,5X,20ES14.5) !! reverted Jan 12, 2013
1001 format(3I4,3E12.3,5X,20ES14.5)
  !---------------------------------------------------------------------! 
contains
  subroutine HkAndSkSetting(ak)
    implicit none
    real(DOUBLE_PRECISION)    :: r2(3),r1(3),r12(3),ak(3),r10(3),r20(3),distance,z1,z2
    complex(DOUBLE_PRECISION) :: phFac
    integer                   :: jO2,jO1,curOrb2,curOrb22,curOrb1,curOrb11,atom111,i1,i2
    ! write(6,*) 'ak=', ak(:)

    do atom2=1,numPosInUnit                       !central atom (reduced into Unit)
       atom22=atom2+numPosInUnit*(whereIsOrigin-1)!central atom (absolute #)
       curOrb22=sumNval(atom22)                   !#orbit (absolute atom #)
       curOrb2 =sumNvalInUnit(atom2)              !#orbit (reduced atom #)
       r20(:)=PosInUnit(:,atom2)+posUnitOrg(:,whereIsOrigin)
       r2(:) =transVec(:,1)*r20(1)+transVec(:,2)*r20(2)+transVec(:,3)*r20(3)
       ! write(6,*) '  atom22 and atom2     =',atom22,atom2
       ! write(6,*) '  curOrb22 and curOrb2 =',curOrb22,curOrb2
       ! write(6,*) '  r2=',r2(:)
    ! S.Y Sep27, 2013, added in order to claculate Ecc.
       Ecc=Ecc+EccPerAtom(atom22)
   do atom1=1,njsd(whereIsOrigin,ict)         !neighbor atom (sequential)
          atom11=jsv4jsd(atom1,atom22)             !neighbor atom (absolute #)
          Unit1=(atom11-1)/numPosInUnit+1
          atom111=atom11-numPosInUnit*(Unit1-1)    !neighbor atom (equivalent site in Uni)
          curOrb11  =sumNval(atom11)               !#orbit (absolute atom #)
          curOrb1   =sumNvalInUnit(atom111)        !#orbit (reduced atom #)
          r10(:)=PosInUnit(:,atom111)+posUnitOrg(:,Unit1)
          r1(:) =transVec(:,1)*r10(1)+transVec(:,2)*r10(2)+transVec(:,3)*r10(3)
          r12(:)=r1(:)-r2(:)  
          distance=dsqrt(r12(1)**2+r12(2)**2+r12(3)**2)
          phFac=cdexp(aim*dot_product(ak,r12))
          ! write(6,100)     atom11,atom111,distance
100       format('      atom11,atom111,distance=', 2I5,D13.5)
          ! write(6,*) '    atom11,atom111 =',atom11,atom111
          ! write(6,*) '    curOrb11,curOrb1=',curOrb11,curOrb1
          ! write(6,*) '    r1=',r1(:)
          ! write(6,*) '    |r12|=',distance
          ! write(6,*) '    phFax=',phFac 
          do jO2=1,nval(atom22)
             do jO1=1,nval(atom11)
                Hk(curOrb1+jO1,curOrb2+jO2)=Hk(curOrb1+jO1,curOrb2+jO2) &
                     + phFac*densH(curOrb11+jO1,curOrb22+jO2)
                Sk(curOrb1+jO1,curOrb2+jO2)=Sk(curOrb1+jO1,curOrb2+jO2) &
                     + phFac*densS(curOrb11+jO1,curOrb22+jO2)
             end do
          end do
       end do
    end do
    !
    do i2=1,iDimHam
       do i1=1,iDimHam
          z1=real(Hk(i1,i2),DOUBLE_PRECISION)
          z2=imag(Hk(i1,i2))
          !! if(dabs(z1).le.1.0d-7) z1=0.0d0
          !! if(dabs(z2).le.1.0d-7) z2=0.0d0
          Hk(i1,i2)=cmplx(z1,z2,DOUBLE_PRECISION)
          z1=real(Sk(i1,i2),DOUBLE_PRECISION)
          z2=imag(Sk(i1,i2))
          !! if(abs(z1).le.1.0d-7) z1=0.0d0
          !! if(abs(z2).le.1.0d-7) z2=0.0d0
          Sk(i1,i2)=cmplx(z1,z2,DOUBLE_PRECISION)
       end do
    end do
  end subroutine HkAndSkSetting
  !
  !---------------------------------------------------------------------

  subroutine sparse2dens(mat,densMat,ICHECK)
    real(DOUBLE_PRECISION) :: mat(:,:,:,:), densMat(:,:), acheck
    intent(IN) :: mat
    intent(OUT):: densMat
    integer    :: curOrb1, curOrb2, i1,i2,ICHECK
    !   integer    :: k
    !   integer, allocatable :: sumNval(:)
    !   allocate(sumNval(noav))
    !   k=0
    !   do jsv1=1, noav
    !      sumNval(jsv1)=k                 !sequencial number 
    !      k=k+nval(jsv1)
    !   end do
    do jsv2=1, noav                    !cetral atom
       njsd2=njsd(jsv2,ict)
       nval2=nval(jsv2)
       curOrb2=sumNval(jsv2)
       do jsd1=1, njsd2                !nisd2=total number of neighboring atoms 
          jsv1=jsv4jsd(jsd1,jsv2)      !neighboring atom (absokute #) 
          curOrb1=sumNval(jsv1)
          nval1=nval(jsv1)
          do j2=1, nval2
             do j1=1, nval1
                densMat(curOrb1+j1,curOrb2+j2)=mat(j1,j2,jsd1,jsv2)
             end do
          end do
       end do
    end do
    ! check routine ====== start
    ICHECK=0
    do i1=1,nOrbs
       do i2=1,nOrbs
          acheck=densMat(i1,i2)-densMat(i2,i1)
          if(dabs(acheck).gt.1.0d-7) ICHECK=ICHECK+1
       end do
    end do
    ! check routine ====== end
    !deallocate(sumNval)
  end subroutine sparse2dens
  !---------------------------------------------------------------------
  subroutine VecPro(aVec,bVec)
    implicit none
    real(DOUBLE_PRECISION) :: aVec(3,3),bVec(3,3)
    bVec(:,1)=outer_prod(aVec(:,2),aVec(:,3))
    bVec(:,2)=outer_prod(aVec(:,3),aVec(:,1))
    bVec(:,3)=outer_prod(aVec(:,1),aVec(:,2))
  end subroutine VecPro
  !---------------------------------------------------------------------
  function outer_prod(u,v)
    implicit none
    real(8):: outer_prod(3)
    real(8),intent(in)::u(3),v(3)
    outer_prod(1)=u(2)*v(3)-u(3)*v(2)
    outer_prod(2)=u(3)*v(1)-u(1)*v(3)
    outer_prod(3)=u(1)*v(2)-u(2)*v(1)
  end function outer_prod
  !---------------------------------------------------------------------
!!$      function dot_prod(u,v)
!!$      implicit none
!!$      real(8):: dot_prod
!!$      real(8),intent(in)::u(3),v(3)
!!$      dot_prod =u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
!!$      end function dot_prod
  !---------------------------------------------------------------------
end program band_calc


