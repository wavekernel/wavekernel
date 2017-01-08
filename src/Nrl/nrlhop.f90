!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! A module for calculating off-site elements in NRL scheme
!
module MNrlHopping 
  use MSlaterKosterTable
  use MSlaterKosterDerivative
  use MNrlSystem
  use MNrlPairMatrix

  implicit none

  private

  !!
  ! public types
  !
  public TNrlHoppingData
  public TNrlHoppingBond
  public TNrlHoppingMatrix
  public TNrlHopBuilder


  !!
  ! public procedures related to TNrlHopBuilder
  !
  public NrlInitialize
  public NrlRelease
  public NrlGetHopBond
  public NrlGetMatrix
  public NrlGetDerivative
  public NrlPrepareForRead
  public NrlReadFromUnit
  public NrlAllocateHoppingMatrixXML
  public NrlHoppingMatrixLoadXML

  !!
  ! public constants related to TNrlHopBuilder
  !
  public NRLIO_XB_HOPPING
  public NRLIO_NRL_HOPPING
  public NRLIO_MULTI_POLE

  !!
  ! the dependency of the hopping by the distance
  ! (h1 + h2 (r/r0) + h3 (r/r0)^2) exp(-g2((r0/rc) - 1))
  ! 
  ! r0, rc, delta is a parameter for cutoff function
  !
  type TNrlHoppingData
    real(8) :: h(3)
    real(8) :: g2
    real(8) :: r0
    real(8) :: rc
    real(8) :: delta
  end type

  !!
  ! the hopping of each orbital combination
  !
  type TNrlHoppingBond
    type(TNrlHoppingData) :: m(4)
    integer :: m_max
  end type

  !!
  ! the hopping of each species combination.
  !
  type TNrlHoppingMatrix
    integer :: species1
    integer :: species2
    integer :: orbital_count1
    integer :: orbital_count2             
    type(TNrlHoppingBond), pointer :: bonds(:, :)
  end type

  !!
  ! Stores the data for hopping parameter.
  !
  type TNrlHopBuilder
    private
    integer :: species_count
    real(8) :: cutoff

    type(TNrlHoppingMatrix), pointer :: matrices(:, :)
  end type
 
  type TFileIO
    integer :: port
    integer :: mode
  end type

  integer, parameter :: NRLIO_XB_HOPPING  = 1, &
                        NRLIO_NRL_HOPPING = 2, &
                        NRLIO_MULTI_POLE  = 3

  ! public aliases
  interface NrlInitialize
    module procedure InitHopBuilder
  end interface

  interface NrlRelease
    module procedure ReleaseHopBuilder
  end interface

  interface NrlGetMatrix
    module procedure GetMatrixOfAtomPair
  end interface

  interface NrlGetHopBond
    module procedure GetHopBond
  end interface

  interface NrlPrepareForRead
    module procedure PrepareHopParameters
  end interface

  interface NrlReadFromUnit
    module procedure ReadHopParameters
  end interface

  interface NrlGetDerivative
    module procedure GetDerivativeOfAtomPair
  end interface

  interface NrlAllocateHoppingMatrixXML
    module procedure AllocateHoppingMatrixXML
  end interface

  interface NrlHoppingMatrixLoadXML
    module procedure HoppingMatrixLoadXML
  end interface

contains

  !!
  ! set cut-off length of bond 
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetBondCutOff(iodata, ir0, irc, idelta)
    type(TNrlHoppingData), intent(inout) :: iodata
    real(8),     intent(in)    :: ir0
    real(8),     intent(in)    :: irc
    real(8),     intent(in)    :: idelta

    iodata%r0    = ir0
    iodata%rc    = irc
    iodata%delta = idelta
  end subroutine

  !!
  ! return hopping without cut-off
  !
  real(8) function GetBareHopping(idata, ir) result (oresult)
    type(TNrlHoppingData), intent(in) :: idata
    real(8),     intent(in) :: ir

    real(8) :: the_r

    the_r = ir / idata%r0

    oresult = (idata%h(1) + the_r * (idata%h(2) + the_r * idata%h(3))) &
              * exp(-idata%g2 * (the_r - 1d0))
  end function

  !!
  ! return the derivative of hoppng function without cut-off
  !
  real(8) function GetBareHoppingDerivative(idata, ir) result (oresult)
    type(TNrlHoppingData), intent(in) :: idata
    real(8),     intent(in) :: ir

    real(8) :: the_r

    the_r = ir / idata%r0

    oresult =   (-(idata%h(1) + the_r*(idata%h(2) + the_r*idata%h(3)))*idata%g2 &
                 + idata%h(2) + 2 * the_r * idata%h(3)) / idata%r0            &
               * exp(-idata%g2 * (the_r - 1d0))
  end function

  !!
  ! return hopping e.g. Vspp
  !
  real(8) function GetHopping(idata, ir) result (oresult)
    use MNrlCutOff
    type(TNrlHoppingData), intent(in) :: idata
    real(8),     intent(in) :: ir

    oresult =  GetBareHopping(idata, ir) &
             * NrlCutOffFunction(ir, idata%rc, idata%delta)
  end function

  !!
  ! return the derivative of the Hopping function
  !
  real(8) function GetHoppingDerivative(idata, ir) result (oresult)
    use MNrlCutOff
    type(TNrlHoppingData), intent(in) :: idata
    real(8),     intent(in) :: ir

    oresult = GetBareHopping(idata, ir)                            &
              * NrlCutOffDerivative(ir, idata%rc, idata%delta)     &
            + GetBareHoppingDerivative(idata, ir)                  &
              * NrlCutOffFunction(ir, idata%rc, idata%delta)
   end function

  !!
  ! Initializes TNrlHoppingMatrix object
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitHoppingMatrix(iomatrix, il_count1, il_count2)
    type(TNrlHoppingMatrix), intent(inout)         :: iomatrix
    integer,      intent(in) :: il_count1, il_count2
    integer :: the_err

    allocate(iomatrix%bonds(il_count1, il_count2), STAT=the_err)
    if(the_err .ne. 0) stop 'MEMORY ALLOCATION ERROR in Hop'

    iomatrix%orbital_count1 = il_count1
    iomatrix%orbital_count2 = il_count2
  end subroutine

  !!
  ! Return hopping interaction
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetHopBond(ibond, ir, obond)
    type(TNrlHoppingBond), intent(in) :: ibond
    real(8), intent(in) :: ir
    type(TSKBond) :: obond

    integer :: m

    do m = 1, ibond%m_max + 1
      obond%v(m) = GetHopping(ibond%m(m), ir)
    end do
  end subroutine

  !!
  ! interface routine to calculate position derivative
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine PositionDerivative(ibond, ilm1, ilm2, ipair, oresult)
    type(TNrlHoppingBond), intent(in) :: ibond
    type(TNrlAtomPair), intent(in) :: ipair
    integer, intent(in) :: ilm1, ilm2
    real(8), intent(out) :: oresult(:)

    type(TSKBond) :: the_bond

    call GetHopBond(ibond, ipair%distance, the_bond)
    call SlaterKosterDerivative(the_bond, ilm1, ilm2, &
                                ipair%l, ipair%m, ipair%n, oresult)

    oresult(1:3) = oresult(1:3) / ipair%distance
  end subroutine

  !!
  ! return the derivative hopping interaction
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetHoppingDerivativeList(ibond, ilm1, ilm2, ir, obond)
    type(TNrlHoppingBond), intent(in) :: ibond
    integer, intent(in) :: ilm1, ilm2
    real(8), intent(in) :: ir
    type(TSKBond), intent(out) :: obond
    integer :: the_lmin, m

    the_lmin = LesserL(ilm1, ilm2)

    do m = 1, the_lmin + 1
      obond%v(m) = GetHoppingDerivative(ibond%m(m), ir)
    end do
  end subroutine

  !!
  ! Return the potential part of orbital force.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine PotentialDerivative(ibond, ilm1, ilm2, ipair, oresult)
    type(TNrlHoppingBond), intent(in) :: ibond
    type(TNrlAtomPair), intent(in) :: ipair
    integer, intent(in) :: ilm1, ilm2
    real(8), intent(inout) :: oresult(:)

    type(TSKBond) :: the_bond
    real(8) :: the_hop

    call GetHoppingDerivativeList(ibond, ilm1, ilm2, ipair%distance, the_bond)

    the_hop = SlaterKosterCoefficient(the_bond, ilm1, ilm2,&
                                      ipair%l, ipair%m, ipair%n)

    oresult(1) = ipair%l * the_hop
    oresult(2) = ipair%m * the_hop
    oresult(3) = ipair%n * the_hop
  end subroutine

  !!
  ! Finalizes TNrlHoppingMatrix object 
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseHoppingMatrix(iomatrix)
    type(TNrlHoppingMatrix), intent(inout) :: iomatrix

    !-- it is always allocated if it has been initialized 
    deallocate(iomatrix%bonds)
  end subroutine

  !!
  ! return the minumum angular momentum between given two orbital
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetL(imatrix, iorbital1, iorbital2, ol1, ol2) 
    type(TNrlHoppingMatrix), intent(in)  :: imatrix
    integer,      intent(in)  :: iorbital1, iorbital2
    integer,      intent(out) :: ol1, ol2

    ol1 = iorbital1 - 1
    ol2 = iorbital2 - 1
  end subroutine

  !!
  ! Calculates l-index from lm-index.
  !
  integer function LFromLM(ilm) result(oresult)
    integer, intent(in) :: ilm
    integer :: inc, lmmax
    oresult = 0
    lmmax = 0
    inc = 1
    do while(lmmax < ilm)
      oresult = oresult + 1
      inc = inc + 2
      lmmax = lmmax + inc
    end do    
  end function

  !!
  ! Calculates lesser l-index from two lm-index
  !
  integer function LesserL(ilm1, ilm2) result(oresult)
    integer, intent(in) :: ilm1, ilm2
    integer :: min_lm

    min_lm = min(ilm1, ilm2)
    oresult = LFromLM(min_lm)
  end function

  !!
  ! Sets cut-off parameters for all the combination of orbials.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetCutOffFunction(iomatrix, ir0, irc, idelta)
    type(TNrlHoppingMatrix), intent(inout) :: iomatrix
    real(8),       intent(in)    :: ir0, irc, idelta

    integer :: i, j
    integer :: the_l1, the_l2, the_lmin, m

    do i = 1, iomatrix%orbital_count1
      do j = 1, iomatrix%orbital_count2

        call GetL(iomatrix, i, j, the_l1, the_l2)
        the_lmin = min(the_l1, the_l2)

        do m = 1, the_lmin + 1
          call SetBondCutOff(iomatrix%bonds(i, j)%m(m), ir0, irc, idelta)
        end do
      end do
    end do
  end subroutine

  !!
  ! Initialize routine for TNrlHopBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitHopBuilder(iohop)
    type(TNrlHopBuilder), intent(inout) :: iohop

    nullify(iohop%matrices)
  end subroutine

  !!
  ! Finalize routine for TNrlHopBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseHopBuilder(iohop)
    type(TNrlHopBuilder), intent(inout) :: iohop

    call DeallocateHopBuilder(iohop)
  end subroutine

  !!
  ! Allocation routine for TNrlHoppingMatrix
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AllocateHopBuilder(iohop, ispecies_count, iorbital_count)
    type(TNrlHopBuilder), intent(inout) :: iohop
    integer, intent(in) :: ispecies_count
    integer, intent(in) :: iorbital_count(:)

    integer :: i, j, the_err
    
    call DeallocateHopBuilder(iohop)

    allocate(iohop%matrices(ispecies_count, ispecies_count), STAT=the_err)
    if(the_err .ne. 0) stop 'MEMORY ALLOCATION ERROR in AllocateHopBuilder'

    iohop%species_count = ispecies_count

    do i = 1, ispecies_count
      do j = 1, ispecies_count
        call InitHoppingMatrix(iohop%matrices(j, i), iorbital_count(j), iorbital_count(i))
      end do
    end do
  end subroutine

  !!
  ! Deallocation of main part of TNrlHopBuilder
  ! The consistensy in the object is remained.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine DeallocateHopBuilder(iohop)
    type(TNrlHopBuilder), intent(inout) :: iohop

    integer :: i, j, the_count

    if (associated(iohop%matrices)) then

      the_count = iohop%species_count
      do i = 1, the_count
        do j = 1, the_count
          call ReleaseHoppingMatrix(iohop%matrices(j, i))
        end do
      end do

      deallocate(iohop%matrices)
    end if
  end subroutine

  !!
  ! Returns a hopping matrix of given atom combination.
  ! 
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetMatrixOfAtomPair(iohop, ipair, iomatrix)
    type(TNrlHopBuilder), intent(inout) :: iohop
    type(TNrlAtomPair), intent(in) :: ipair
    type(TNrlPairMatrix), intent(inout) :: iomatrix

    type(TSKBond), allocatable :: the_bonds(:,:)
    type(TNrlHoppingMatrix), pointer :: the_hop_matrix
    integer :: i, j, il1, il2, lm1, lm2, il_count1, il_count2
    integer :: the_species1, the_species2, the_err

    call NrlSetSize(iomatrix, ipair%orbit_count2, ipair%orbit_count1)

    the_species1 = ipair%species1
    the_species2 = ipair%species2 

    the_hop_matrix => iohop%matrices(the_species1, the_species2)

    ! This does not care about orbital index and angular momentum
    ! correspondence, that is, quick hack.
    il_count1 = ipair%il_size1
    il_count2 = ipair%il_size2

    allocate(the_bonds(il_count1, il_count2), stat=the_err)
    if (the_err /= 0) stop "Memory allocation Error in GetMatrixOfAtomPair"

    ! Gets hopping parameter of each bonds such as Vss, Vsp, ..
    do il1 = 1, il_count1
      do il2 = 1, il_count2
        call GetHopBond(the_hop_matrix%bonds(il1, il2), ipair%distance, &
                           the_bonds(il1, il2))
      end do 
    end do

!print "(3F12.4)", the_hop_matrix%bonds(3, 3)%m(2)%h(1:3)
    ! Gets hopping matrix such as Vss, Vsx, ..    
    do i = 1, ipair%orbit_count1
      lm1 = ipair%row_to_lm1(i)
      il1 = ipair%row_to_il1(i)
      do j = 1, ipair%orbit_count2
        lm2 = ipair%row_to_lm2(j)
        il2 = ipair%row_to_il2(j)
        iomatrix%e(j, i) = &
           SlaterKosterCoefficient(the_bonds(il1, il2), lm1, lm2, &
                                   ipair%l, ipair%m, ipair%n)
      end do
    end do
    deallocate(the_bonds)
  end subroutine

  !!
  ! Calculates the derivative of hopping between two atom.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetDerivativeOfAtomPair(iobuilder, ipair, iomatrix)
    type(TNrlHopBuilder), intent(inout) :: iobuilder
    type(TNrlAtomPair), intent(in) :: ipair
    type(TNrlPairForceMatrix), intent(inout) :: iomatrix

    type(TSKBond), allocatable :: the_bonds(:,:)
    type(TNrlHoppingMatrix), pointer :: the_matrix
    integer :: i, j, lm1, lm2, il1, il2
    real(8) :: d1(1:3), d2(1:3)

    the_matrix => iobuilder%matrices(ipair%species1, ipair%species2)

    do i = 1, ipair%orbit_count1
      lm1 = ipair%row_to_lm1(i)
      il1 = ipair%row_to_il1(i)
      do j = 1, ipair%orbit_count2
        lm2 = ipair%row_to_lm2(j)
        il2 = ipair%row_to_il2(j)
        call PositionDerivative(the_matrix%bonds(il1, il2), &
                                   lm1, lm2, ipair, d1(:))
        call PotentialDerivative(the_matrix%bonds(il1, il2), &
                                    lm1, lm2, ipair, d2(:))       

        ! dr_ij/dR_i = - r_ij/r
        iomatrix%e(1:3, j, i) = - (d1(1:3) + d2(1:3))
      end do
    end do
  end subroutine

  !!
  ! Reads cut-off function
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadCutOffFunction(iohop, ifile)
    type(TNrlHoppingMatrix), intent(inout) :: iohop
    type(TFileIO), intent(in) :: ifile

    real(8) :: r0, rc, delta
    integer :: the_err

    select case(ifile%mode)
    case(NRLIO_NRL_HOPPING) 
      read(ifile%port, *, iostat=the_err) rc, delta
      if(the_err /= 0) goto 800
      r0 = 1d0
    case default
      read(ifile%port, *, iostat=the_err) r0, rc, delta
      if(the_err /= 0) goto 800
    end select

    call SetCutOffFunction(iohop, r0, rc, delta)    
    return
  800 continue
    write(0,*) "Read error in ReadCutOffFunction in Hop"
    stop
  end subroutine

  !!
  ! Inverts hopping data in one line stream. e. g. V_pds
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InvertHoppingData(iodata)
    type(TNrlHoppingData), intent(inout) :: iodata

    iodata%h(1) = -iodata%h(1)
    iodata%h(2) = -iodata%h(2)
    iodata%h(3) = -iodata%h(3)
  end subroutine

  !!
  ! Reads header part of hopping energy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHopHeader(iunit, omode)
    integer, intent(in) :: iunit
    integer, intent(out) :: omode

    integer :: the_err
    character(len=32) :: the_header
    read(iunit, '(A32)', iostat=the_err) the_header
    if (the_err /= 0) goto 800

    if(the_header /= 'HOP NRL') goto 810
 
    omode = NRLIO_NRL_HOPPING   
    return

  800 continue
    write(0, *) "There is no-header for hop part."
    stop
  810 continue
    write(0, *) "Hopping header "//trim(the_header)//" is not supported."
    stop
  end subroutine

  !!
  ! Reads Hopping Data in one line stream. e. g. V_pds
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHoppingData(iodata, ifile)
    type(TNrlHoppingData), intent(inout) :: iodata
    type(TFileIO),      intent(in)    :: ifile

    real(8) :: h1, h2, h3, g2, e, f, fbar, g

    select case(ifile%mode)
    case(NRLIO_XB_HOPPING, NRLIO_MULTI_POLE)
      read(ifile%port, *) h1, h2, g2
      iodata%h(1) = h1
      iodata%h(2) = h2
      iodata%h(3) = 0d0
      iodata%g2 = g2
    case(NRLIO_NRL_HOPPING) 
      read(ifile%port, *) g, e, f, fbar
      ! temporal modification
      ! e = e * 0.5d0
      ! f = f * 0.5d0
      ! fbar = fbar * 0.5d0
      iodata%h(1) = e    * exp(-g * g)
      iodata%h(2) = f    * exp(-g * g)
      iodata%h(3) = fbar * exp(-g * g)
      iodata%g2 = g * g
    case default
      write(0, *) "LOGIC ERROR in ReadHoppingData: There is no mode ", &
                  ifile%mode
      stop
    end select
  end subroutine

  !!
  ! Inverse sign of the bond
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InvertBond(iobond)
    type(TNrlHoppingBond), intent(inout) :: iobond
    integer :: i

    do i = 1, iobond%m_max + 1
      call InvertHoppingData(iobond%m(i))
    end do
  end subroutine
    
  !!
  ! Read Hopping Bond e.g. from sigma to pi of V_pd.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHoppingBond(iobond, ifile, ilmin)
    type(TNrlHoppingBond), intent(inout) :: iobond
    type(TFileIO), intent(in) :: ifile
    integer,    intent(in)    :: ilmin
    integer :: i
    iobond%m_max = ilmin

    do i = 1, ilmin + 1
      call ReadHoppingData(iobond%m(i), ifile)
    end do
  end subroutine

  !!
  ! Read matrix data from stream of given file port
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHoppingMatrix(iomatrix, ifile)
    type(TNrlHoppingMatrix), intent(inout) :: iomatrix
    type(TFileIO),        intent(in)    :: ifile

    integer             :: the_orbital1, the_orbital2
    integer             :: i, j, the_l1, the_l2

    the_orbital1 = iomatrix%orbital_count1
    the_orbital2 = iomatrix%orbital_count2

    do i = 1, the_orbital1
      do j = i, the_orbital2
        call GetL(iomatrix, i, j, the_l1, the_l2)
        call ReadHoppingBond(iomatrix%bonds(i, j), ifile, min(the_l1, the_l2))
        iomatrix%bonds(j, i) = iomatrix%bonds(i, j)

        ! sign is negative if the parity of the sum of angular momentum is odd.
        select case(ifile%mode)
        case(NRLIO_XB_HOPPING, NRLIO_NRL_HOPPING)
          if (mod(the_l1 + the_l2, 2) == 1) then
            call InvertBond(iomatrix%bonds(j, i))
          end if
        end select
      end do
    end do
!print "(A,2I4,3F12.4)", "Read Bond", 1, 1, iomatrix%bonds(1, 1)%m(1)%h(1:3)
  end subroutine

  !!
  ! Prepares before reading file.
  ! Before calling ReadHopParameters, you must call this routine
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine PrepareHopParameters(iohop, ispecies_count, iorbital_count)
    type(TNrlHopBuilder), intent(inout) :: iohop
    integer, intent(in) :: ispecies_count
    integer, intent(in) :: iorbital_count(:)

    call AllocateHopBuilder(iohop, ispecies_count, iorbital_count)
  end subroutine
   
  !!
  ! interface routine to read hopping parameters in NRL scheme.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHopParameters(iohop, iunit)
    type(TNrlHopBuilder), intent(inout) :: iohop
    integer,       intent(in)    :: iunit

    type(TFileIO)   :: the_file
    integer :: the_mode
    integer :: i, j

    call ReadHopHeader(iunit, the_mode)

    the_file%port = iunit
    the_file%mode = the_mode

    do i = 1, iohop%species_count
      do j = i, iohop%species_count
        call ReadCutOffFunction(iohop%matrices(j, i), the_file) 
        call ReadHoppingMatrix(iohop%matrices(j, i), the_file)
      end do
    end do

    !one-line skip
    read(iunit, *)
  end subroutine

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AllocateHoppingMatrixXML(iohop, the_species_count, the_orbital_count)
    implicit none
    integer, intent(in) :: the_species_count
    integer, intent(in) :: the_orbital_count(the_species_count)
    integer i, j, ierr 
    type(TNrlHopBuilder), intent(inout) :: iohop

    ! deallocate
    if( associated(iohop%matrices) ) then
       do i = 1, the_species_count
          do j = 1, the_species_count
             deallocate( iohop%matrices(j,i)%bonds )
          end do
       end do
       deallocate(iohop%matrices)
    end if
    
    ! allocate
    allocate( iohop%matrices(the_species_count,the_species_count), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocation matrices in sbrt allocate_hamiltonian'
    end if
    
    do i = 1, the_species_count
       do j = 1, the_species_count
          allocate( iohop%matrices(j,i)%bonds(the_orbital_count(j),the_orbital_count(i) ), &
               stat=ierr )
          if( ierr /= 0 ) then
             write(*,*) 'Error : allocation bond in sbrt allocate_hamiltonian'
          end if
          iohop%matrices(j,i)%orbital_count1 = the_orbital_count(j)
          iohop%matrices(j,i)%orbital_count2 = the_orbital_count(i)
       end do
    end do

  end subroutine AllocateHoppingMatrixXML

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine HoppingMatrixLoadXML(iohop, vnode,delta,cutoff,is,js)
    use flib_dom
    use elses_xml_misc
    implicit none
    type(fnode), target, intent(in) :: vnode
    type(fnode), pointer :: node
    integer, intent(in) :: is, js
    real(8), intent(in) :: delta, cutoff
    
    real(8) rvalue(4)
    character(len=256) value, unit, bond
    type(TNrlHopBuilder), intent(inout) :: iohop
    
    node => vnode

    unit  = getAttribute(node,"unit")
    if( unit == "" ) then
       unit = "a.u."
    endif

    value = getChildValue(node)
    read(unit=value,fmt=*) rvalue(2:4), rvalue(1)
    rvalue(1:4) = 2.d0*rvalue(1:4) * XML_TO_AU(unit) ! Ry -> a.u.->Ry
    
    bond  = getAttribute(node,"bond")
    select case(bond)
    case("sss")
       iohop%matrices(is,js)%bonds(1,1)%m_max = 0
       iohop%matrices(is,js)%bonds(1,1)%m(1)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(1,1)%m(1)%rc = cutoff - 5.d0*delta 
       iohop%matrices(is,js)%bonds(1,1)%m(1)%delta = delta
       iohop%matrices(is,js)%bonds(1,1)%m(1)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,1)%m(1)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,1)%m(1)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,1)%m(1)%g2 = rvalue(1)*rvalue(1)
       
    case("sps")
       iohop%matrices(is,js)%bonds(1,2)%m_max = 0
       iohop%matrices(is,js)%bonds(1,2)%m(1)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(1,2)%m(1)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(1,2)%m(1)%delta = delta
       iohop%matrices(is,js)%bonds(1,2)%m(1)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,2)%m(1)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,2)%m(1)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,2)%m(1)%g2 = rvalue(1)*rvalue(1)
       
       iohop%matrices(is,js)%bonds(2,1)%m(1) = iohop%matrices(is,js)%bonds(1,2)%m(1)
       iohop%matrices(is,js)%bonds(2,1)%m_max = iohop%matrices(is,js)%bonds(1,2)%m_max
       iohop%matrices(is,js)%bonds(2,1)%m(1)%h(1) = -rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,1)%m(1)%h(2) = -rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,1)%m(1)%h(3) = -rvalue(4)*exp(-rvalue(1)*rvalue(1))
       
    case("sds")
       iohop%matrices(is,js)%bonds(1,3)%m_max = 0
       iohop%matrices(is,js)%bonds(1,3)%m(1)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(1,3)%m(1)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(1,3)%m(1)%delta = delta
       iohop%matrices(is,js)%bonds(1,3)%m(1)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,3)%m(1)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,3)%m(1)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(1,3)%m(1)%g2 = rvalue(1)*rvalue(1)
       
       iohop%matrices(is,js)%bonds(3,1)%m(1) = iohop%matrices(is,js)%bonds(1,3)%m(1)
       iohop%matrices(is,js)%bonds(3,1)%m_max = iohop%matrices(is,js)%bonds(1,3)%m_max
         
    case("pps")
       iohop%matrices(is,js)%bonds(2,2)%m_max = 1
       iohop%matrices(is,js)%bonds(2,2)%m(1)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(2,2)%m(1)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(2,2)%m(1)%delta = delta
       iohop%matrices(is,js)%bonds(2,2)%m(1)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,2)%m(1)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,2)%m(1)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,2)%m(1)%g2 = rvalue(1)*rvalue(1)
       
    case("ppp")
       iohop%matrices(is,js)%bonds(2,2)%m_max = 1
       iohop%matrices(is,js)%bonds(2,2)%m(2)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(2,2)%m(2)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(2,2)%m(2)%delta = delta
       iohop%matrices(is,js)%bonds(2,2)%m(2)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,2)%m(2)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,2)%m(2)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,2)%m(2)%g2 = rvalue(1)*rvalue(1)
       
    case("pds")
       iohop%matrices(is,js)%bonds(2,3)%m_max = 1
       iohop%matrices(is,js)%bonds(2,3)%m(1)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(2,3)%m(1)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(2,3)%m(1)%delta = delta
       iohop%matrices(is,js)%bonds(2,3)%m(1)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,3)%m(1)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,3)%m(1)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,3)%m(1)%g2 = rvalue(1)*rvalue(1)
       
       iohop%matrices(is,js)%bonds(3,2)%m(1) = iohop%matrices(is,js)%bonds(2,3)%m(1)
       iohop%matrices(is,js)%bonds(3,2)%m_max = iohop%matrices(is,js)%bonds(2,3)%m_max
       iohop%matrices(is,js)%bonds(3,2)%m(1)%h(1) = -rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,2)%m(1)%h(2) = -rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,2)%m(1)%h(3) = -rvalue(4)*exp(-rvalue(1)*rvalue(1))
       
    case("pdp")
       iohop%matrices(is,js)%bonds(2,3)%m_max = 1
       iohop%matrices(is,js)%bonds(2,3)%m(2)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(2,3)%m(2)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(2,3)%m(2)%delta = delta
       iohop%matrices(is,js)%bonds(2,3)%m(2)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,3)%m(2)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,3)%m(2)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(2,3)%m(2)%g2 = rvalue(1)*rvalue(1)
       
       iohop%matrices(is,js)%bonds(3,2)%m(2) = iohop%matrices(is,js)%bonds(2,3)%m(2)
       iohop%matrices(is,js)%bonds(3,2)%m_max = iohop%matrices(is,js)%bonds(2,3)%m_max
       iohop%matrices(is,js)%bonds(3,2)%m(2)%h(1) = -rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,2)%m(2)%h(2) = -rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,2)%m(2)%h(3) = -rvalue(4)*exp(-rvalue(1)*rvalue(1))
    case("dds")
       iohop%matrices(is,js)%bonds(3,3)%m_max = 2
       iohop%matrices(is,js)%bonds(3,3)%m(1)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(3,3)%m(1)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(3,3)%m(1)%delta = delta
       iohop%matrices(is,js)%bonds(3,3)%m(1)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(1)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(1)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(1)%g2 = rvalue(1)*rvalue(1)
       
    case("ddp")
       iohop%matrices(is,js)%bonds(3,3)%m_max = 2
       iohop%matrices(is,js)%bonds(3,3)%m(2)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(3,3)%m(2)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(3,3)%m(2)%delta = delta
       iohop%matrices(is,js)%bonds(3,3)%m(2)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(2)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(2)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(2)%g2 = rvalue(1)*rvalue(1)
       
    case("ddd")
       iohop%matrices(is,js)%bonds(3,3)%m_max = 2
       iohop%matrices(is,js)%bonds(3,3)%m(3)%r0 = 1.d0
       iohop%matrices(is,js)%bonds(3,3)%m(3)%rc = cutoff-5.d0*delta 
       iohop%matrices(is,js)%bonds(3,3)%m(3)%delta = delta
       iohop%matrices(is,js)%bonds(3,3)%m(3)%h(1) = rvalue(2)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(3)%h(2) = rvalue(3)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(3)%h(3) = rvalue(4)*exp(-rvalue(1)*rvalue(1))
       iohop%matrices(is,js)%bonds(3,3)%m(3)%g2 = rvalue(1)*rvalue(1)
       
    end select

    return
  end subroutine HoppingMatrixLoadXML
  
end module
