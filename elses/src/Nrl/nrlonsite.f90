!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!
! This module treats the onsite element in NRL scheme.
!
! Notes for the notation 
!   Deriv : an abbreviation of "derivative"
!   Prime : an alias of "1st derivative"
!
module MNrlOnSite
  use MNrlSystem
  use MNrlPairMatrix

  implicit none

  private

  !!
  ! public types
  !
  public TNrlOnSiteBuilder

  !!
  ! public procedures for TNrlOnSiteBuilder
  !
  public NrlInitialize
  public NrlRelease
  public NrlGetMatrix
  public NrlGetDerivative
  public NrlPrepareForRead
  public NrlReadFromUnit
  public NrlSetIdentity
  public NrlAllocateOnsiteXML
  public NrlOnsiteSetXML
  public NrlOnsiteLoadXML
  public NrlOnsiteSetIdentityXML

  !!
  ! public constants
  !
  public NRLIO_XB_ONSITE
  public NRLIO_NRL_ONSITE
  public NRLIO_NRLR2_ONSITE

  type TNrlOnSiteDensity
    private
    real(8) :: a
    real(8) :: gamma
    real(8) :: delta
    real(8) :: rc
    real(8) :: r0
  end type

  type TNrlOnSiteEnergyParams
    real(8) :: a
    real(8) :: b
    real(8) :: c
    real(8) :: d
    ! for on-site correction
    real(8) :: e
    real(8) :: gamma
    real(8) :: rho_c
  end type

  type TNrlOnSiteEnergy
    private
    integer :: orbital_count
    type(TNrlOnSiteEnergyParams), pointer :: params(:)
  end type

  type TFileParam
    integer :: port
    integer :: mode
  end type

  type TNrlOnSiteBuilder
    private
    integer :: mode
    integer :: species_count
    type(TNrlOnSiteEnergy),  pointer :: energies(:)
    type(TNrlOnSiteDensity), pointer :: densities(:,:)
    real(8), pointer :: density_cache(:)
  end type

  integer, parameter :: NRLIO_XB_ONSITE    = 1, &
                        NRLIO_NRL_ONSITE   = 2, &
                        NRLIO_NRLR2_ONSITE = 3

  integer, parameter :: NRL_ONSITE_STANDARD = 1, &
                        NRL_ONSITE_CONSTANT = 2

  ! public aliases
  interface NrlInitialize
    module procedure InitOnSiteBuilder
  end interface

  interface NrlRelease
    module procedure ReleaseOnSiteBuilder
  end interface

  interface NrlGetMatrix
    module procedure GetOnSiteOfAtom
  end interface

  interface NrlGetDerivative
    module procedure GetOnSiteDerivative
  end interface

  interface NrlPrepareForRead
    module procedure PrepareForRead
  end interface

  interface NrlReadFromUnit
    module procedure ReadFromUnit
  end interface

  interface NrlSetIdentity
    module procedure SetIdentity
  end interface

  interface NrlAllocateOnsiteXML
    module procedure AllocateOnsiteXML
  end interface

  interface NrlOnsiteSetXML
    module procedure OnsiteSetXML
  end interface

  interface NrlOnsiteLoadXML
    module procedure OnsiteLoadXML
  end interface

  interface NrlOnsiteSetIdentityXML
    module procedure OnsiteSetIdentityXML
  end interface

contains

  !!
  ! Initialize routine for TOnSiteDeisity
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitOnSiteDensity(iodensity)
    type(TNrlOnSiteDensity), intent(inout) :: iodensity

    iodensity%a     = 0d0
    iodensity%gamma = 0d0
    iodensity%delta = 1.0d0
    iodensity%rc    = 5.9d0
    iodensity%r0    = 1d0
  end subroutine
      
  !!
  ! Finalize routine for TNrlOnSiteDensity
  ! Currently, this is only a dummy routine.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseOnSiteDensity(iodensity)
    type(TNrlOnSiteDensity), intent(inout) :: iodensity
  end subroutine

  !!
  ! the density on an atom without cut-off
  !
  real(8) function GetBareOnSiteDensity(idensity, ir) result (oresult) 
    type(TNrlOnSiteDensity), intent(in) :: idensity
    real(8), intent(in) :: ir

    oresult = idensity%a * exp(-idensity%gamma * (ir / idensity%r0 - 1d0))
  end function

  !!
  ! the derivative of the density on an atom without cut-off 
  !
  real(8) function GetBareOnSiteDensityDerivative(idensity, ir) &
       result (oresult)
    type(TNrlOnSiteDensity), intent(in) :: idensity
    real(8),       intent(in) :: ir

    oresult = - idensity%a * idensity%gamma / idensity%r0   &
                 * exp(-idensity%gamma * (ir / idensity%r0 - 1d0))
  end function
  
  !!
  ! routine for the calculation of the density on the atom
  !
  real(8) function GetOnSiteDensity(idensity, ir) result (oresult)
    use MNrlCutOff
    type(TNrlOnSiteDensity), intent(in) :: idensity
    real(8),       intent(in) :: ir

    oresult = GetBareOnSiteDensity(idensity, ir)  &
               * NrlCutOffFunction(ir, idensity%rc, idensity%delta)
  end function

  !!
  ! the derivative of the onsite density
  !
  real(8) function GetDensityDistanceDeriv(idensity, ir) result (oresult)
    use MNrlCutOff
    type(TNrlOnSiteDensity), intent(in) :: idensity
    real(8), intent(in) :: ir

    oresult =   GetBareOnSiteDensityDerivative(idensity, ir)            &
                * NrlCutOffFunction(ir, idensity%rc, idensity%delta)    &
              + GetBareOnSiteDensity(idensity, ir)                      &
                * NrlCutOffDerivative(ir, idensity%rc, idensity%delta) 
  end function

  !!
  ! Read parameters for on-site energy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadOnSiteEnergyParams(ioparams, ifile)
    type(TNrlOnSiteEnergyParams), intent(inout) :: ioparams
    type(TFileParam),     intent(in)    :: ifile
    real(8) :: a, b, c, d, e, rho_c, gamma

    select case(ifile%mode)
    case(NRLIO_XB_ONSITE, NRLIO_NRL_ONSITE)
      read(ifile%port, *) a, b, c, d

      ioparams%a = a
      ioparams%b = b
      ioparams%c = c
      ioparams%d = d

      ioparams%e     = 0d0
      ioparams%rho_c = 100.0d0
      ioparams%gamma = 1.0d0

    case(NRLIO_NRLR2_ONSITE)
      read(ifile%port, *) a, b, c, d, e, rho_c, gamma

      ioparams%a = a
      ioparams%b = b
      ioparams%c = c
      ioparams%d = d

      ioparams%e     = e
      ioparams%rho_c = rho_c
      ioparams%gamma = gamma
    case default
      write(0, *) "Error!: The mode ", ifile%mode, "is not supported."
      write(0, *) "(in ReadOnSiteEnergyParams)"
    end select
  end subroutine


  !!
  ! Set Identity to on-site energy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetOnSiteEnergyParamsIdentity(ioparams)
    type(TNrlOnSiteEnergyParams), intent(inout) :: ioparams

    ioparams%a = 1d0
    ioparams%b = 0d0
    ioparams%c = 0d0
    ioparams%d = 0d0
  end subroutine

  !!
  ! Initialize routine for TNrlOnSiteEnergy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitOnSiteEnergy(ionsite, iorbital_count)
    type(TNrlOnSiteEnergy), intent(inout) :: ionsite
    integer,     intent(in)    :: iorbital_count

    integer :: oerr

    allocate(ionsite%params(iorbital_count), STAT = oerr)
    if (oerr /= 0) stop 'MEMORY ALLOCATION ERROR in OnSite'  

    ionsite%orbital_count = iorbital_count
  end subroutine

  !!
  ! Finalize routine for TNrlOnSiteEnergy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseOnSiteEnergy(ionsite)
    type(TNrlOnSiteEnergy) :: ionsite

    deallocate(ionsite%params)
  end subroutine

  !!
  ! Calculates onsite Hamiltonian element from onsite density
  !
  real(8) function GetOnSiteEnergy(ionsite, iorbital, idensity) &
       result (oresult)
    type(TNrlOnSiteEnergy), intent(in) :: ionsite
    integer,     intent(in) :: iorbital
    real(8),      intent(in) :: idensity

    type(TNrlOnSiteEnergyParams), pointer :: the_c
    real(8)                     :: the_p

    the_c => ionsite%params(iorbital)
    the_p = idensity**(2d0/3d0)

    oresult = ((the_c%d  * the_p + the_c%c) * the_p &
              + the_c%b) * the_p + the_c%a 
    if (idensity > the_c%rho_c) then
      print *, "Density is OverValue", idensity, the_c%rho_c
      oresult = oresult + the_c%e * (idensity - the_c%rho_c) ** the_c%gamma
    end if
    
  end function

  !!
  ! the derivative of the on-site energy polynomial of the given orbital. 
  !
  real(8) function GetOnSitePolynomialDeriv(ionsite, il, idensity) &
       result (oresult)
    type(TNrlOnSiteEnergy), intent(in) :: ionsite
    integer, intent(in) :: il
    real(8), intent(in) :: idensity

    type(TNrlOnSiteEnergyParams), pointer :: the_c
    real(8) :: the_p

    the_c => ionsite%params(il)
    the_p = idensity**(2d0/3d0)

    oresult = ((2d0     * the_c%d  * the_p + 4d0/3d0 * the_c%c) * the_p &
              + 2d0/3d0 * the_c%b) * the_p / idensity

    if (idensity > the_c%rho_c) then
      oresult = oresult + &
        the_c%e * the_c%gamma * (idensity - the_c%rho_c)**(the_c%gamma - 1)
    end if
  end function

  !!
  ! Initialize routine for TNrlOnSiteBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine InitOnSiteBuilder(ioonsite)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite

    ioonsite%mode = NRL_ONSITE_STANDARD
    nullify(ioonsite%energies)
    nullify(ioonsite%densities)
    nullify(ioonsite%density_cache)
  end subroutine

  !!
  ! Finalize routine for TNrlOnSiteBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReleaseOnSiteBuilder(ioonsite)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite

    call DeallocateOnSiteBuilder(ioonsite)
  end subroutine

  !!
  ! Deallocates memory for parameters in TNrlOnSiteBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine DeallocateOnSiteBuilder(ioonsite)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer :: i, j

    do i = 1, ioonsite%species_count
      ! if ioonsite%energies are allocated, %densities is also allocated.
      if (associated(ioonsite%energies)) then
         call ReleaseOnSiteEnergy(ioonsite%energies(i))
         do j = 1, ioonsite%species_count
           call ReleaseOnSiteDensity(ioonsite%densities(j, i))
         end do
      end if
    end do
    ioonsite%species_count = 0
    if(associated(ioonsite%energies)) deallocate(ioonsite%energies)
    if(associated(ioonsite%densities)) deallocate(ioonsite%densities)
    if(associated(ioonsite%density_cache)) deallocate(ioonsite%density_cache)
  end subroutine

  !!
  ! Memory allocation routine for TNrlOnSiteBuilder
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AllocateOnSiteBuilder(ioonsite, ispecies_count, iorbital_count)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer, intent(in) :: ispecies_count
    integer, intent(in) :: iorbital_count(:)

    integer :: i, j, the_err, the_count, the_orbital_count

    call DeallocateOnSiteBuilder(ioonsite)

    the_count = ispecies_count
    ioonsite%species_count = the_count

    allocate(ioonsite%energies(the_count),  &
             ioonsite%densities(the_count, the_count), &
             STAT=the_err)
    if (the_err /= 0) stop 'MEMORY ALLOCATION ERROR in AllocateNrlOnSite'

    do i = 1, the_count
      the_orbital_count = iorbital_count(i)
      call InitOnSiteEnergy(ioonsite%energies(i), the_orbital_count)

      do j = 1, the_count
        call InitOnSiteDensity(ioonsite%densities(j, i))
      end do
    end do
  end subroutine

  !!
  ! Calculates density on the given atom.
  !
  real(8) function GetDensityOnAtom(ionsite, iatom_index) &
       result (oresult)
    type(TNrlOnSitebuilder),    intent(in) :: ionsite
    integer,             intent(in) :: iatom_index

    type(TNrlAtomIterator) :: the_iterator
    type(TNrlAtomPair) :: the_pair
    integer :: the_sp1, the_sp2, jsv2
    real(8) :: the_distance

    call Initialize(the_iterator)

    call Reset(the_iterator, iatom_index)
    oresult = 0d0
    do while(GetNextAtomPair(the_iterator, the_pair))
      the_sp1      = the_pair%species1
      the_sp2      = the_pair%species2
      the_distance = the_pair%distance
      if(the_distance < 1d-12) cycle
      oresult = oresult + GetOnSiteDensity(ionsite%densities(the_sp1, the_sp2),  the_distance) 
    end do
  end function

  !!
  ! Get on-site term of a given atom to Hamiltonian matrix. 
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetOnSiteOfAtom(ioonsite, iatom_index, iomatrix, oerr)
    use MNrlVariables, only: jsei, noav 
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer, intent(in) :: iatom_index
    type(TNrlPairMatrix), intent(inout) :: iomatrix
    integer, intent(out) :: oerr

    real(8) :: the_density
    real(8) :: the_energy
    integer :: l, the_species, row, the_row_size, the_l_size
    character(len=80) :: the_log

    if(.not.associated(ioonsite%density_cache)) then
      allocate(ioonsite%density_cache(noav), stat=oerr)
      if(oerr /= 0) return
    end if

    the_species = jsei(iatom_index)
    the_l_size = ioonsite%energies(the_species)%orbital_count
    the_row_size = the_l_size * the_l_size
    call NrlSetSize(iomatrix, the_row_size, the_row_size)
    iomatrix%e(1:the_row_size, 1:the_row_size) = 0d0

    select case(ioonsite%mode)
    case(NRL_ONSITE_STANDARD)
      the_density = GetDensityOnAtom(ioonsite, iatom_index)
      ioonsite%density_cache(iatom_index) = the_density

      do l = 1, the_l_size
        the_energy = GetOnSiteEnergy(ioonsite%energies(the_species), l, &
                                     the_density)
! print *, l, the_energy
        do row = (l - 1) * (l - 1) + 1, l * l
          iomatrix%e(row, row) =  the_energy
        end do
      end do
    case(NRL_ONSITE_CONSTANT)
      do l = 1, the_l_size
        the_energy = GetOnSiteEnergy(ioonsite%energies(the_species), l, 0.0d0)
        do row = (l - 1) * (l - 1) + 1, l * l
          iomatrix%e(row, row) =  the_energy
        end do
      end do
    end select  
    oerr = 0
  end subroutine

  !!
  ! Calculate the derivative of the matrix.
  ! Return d H_ii / d R_i
  !   i: atom1
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetOnSiteDerivativeOnSite(ioonsite, ipair, iomatrix)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    type(TNrlAtomPair), intent(in) :: ipair
    type(TNrlPairForceMatrix), intent(inout) :: iomatrix

    type(TNrlOnSiteDensity), pointer :: the_density
    type(TNrlOnSiteEnergy), pointer :: the_energy
    type(TNrlAtomIterator) :: the_iterator
    type(TNrlAtomPair) :: the_pair
    integer :: i, il, lm, the_il_size, the_stat
    real(8) :: the_onsite_prime, the_dens_r, the_onsite_dens
    real(8), allocatable :: the_deriv(:, :)

    call NrlSetSize(iomatrix, ipair%orbit_count1, ipair%orbit_count1)
    iomatrix%e(:, :, :) = 0d0

    allocate(the_deriv(3, ipair%il_size1), stat=the_stat)
    if(the_stat /= 0) goto 800

    call Reset(the_iterator, ipair%atom1)
    the_energy  => ioonsite%energies(ipair%species1)
    do while(GetNextAtomPair(the_iterator, the_pair))
      if(the_pair%atom1 == the_pair%atom2) cycle
      the_density => ioonsite%densities(the_pair%species1, the_pair%species2)

      the_dens_r = GetDensityDistanceDeriv(the_density, the_pair%distance) 
      do il = 1, the_pair%il_size1
        the_onsite_dens = &
          GetOnSitePolynomialDeriv(the_energy, il, &
                                   ioonsite%density_cache(the_pair%atom1))

        ! r_ij / R_i = - r_ij / r
        the_deriv(1, il) = - the_onsite_dens * the_dens_r * the_pair%l
        the_deriv(2, il) = - the_onsite_dens * the_dens_r * the_pair%m
        the_deriv(3, il) = - the_onsite_dens * the_dens_r * the_pair%n
      end do
     
      do i = 1, ipair%orbit_count1
        il = ipair%row_to_il1(i)
        iomatrix%e(1:3, i, i) = iomatrix%e(1:3, i, i) + the_deriv(1:3, il)  
      end do
    end do

    deallocate(the_deriv)
    return
  800 continue
    write(0,*) "Memory allocation error in OnSiteDerivative"
    stop
  end subroutine

  !!
  ! Calculate the derivative of the matrix.
  ! Return H_jj / r_i
  !   i: atom1
  !   j: atom2
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine GetOnSiteDerivative(ioonsite, ipair, iomatrix)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    type(TNrlAtomPair), intent(in) :: ipair
    type(TNrlPairForceMatrix), intent(inout) :: iomatrix

    type(TNrlOnSiteDensity), pointer :: the_density
    type(TNrlOnSiteEnergy), pointer :: the_energy
    integer :: i, il, lm, the_il_size, the_stat
    real(8) :: the_onsite_prime, the_dens_r, the_onsite_dens
    real(8), allocatable :: the_deriv(:, :)

    call NrlSetSize(iomatrix, ipair%orbit_count2, ipair%orbit_count2)
    iomatrix%e(:,:,:) = 0d0
    if(ioonsite%mode == NRL_ONSITE_CONSTANT) return

    if (ipair%atom1 == ipair%atom2) then
      call GetOnSiteDerivativeOnSite(ioonsite, ipair, iomatrix)
      return
    end if

    allocate(the_deriv(3, ipair%il_size2), stat=the_stat)
    if(the_stat /= 0) goto 800

    the_energy  => ioonsite%energies(ipair%species2)
    the_density => ioonsite%densities(ipair%species2, ipair%species1)

    the_dens_r = GetDensityDistanceDeriv(the_density, ipair%distance) 
    do il = 1, ipair%il_size1
      the_onsite_dens = &
        GetOnSitePolynomialDeriv(the_energy, il, &
                            ioonsite%density_cache(ipair%atom2))

      ! r_ij / R_i = - r_ij / r
      the_deriv(1, il) = - the_onsite_dens * the_dens_r * ipair%l
      the_deriv(2, il) = - the_onsite_dens * the_dens_r * ipair%m
      the_deriv(3, il) = - the_onsite_dens * the_dens_r * ipair%n
    end do

    iomatrix%e(:, :, :) = 0d0
    do i = 1, ipair%orbit_count2
      il = ipair%row_to_il2(i)
      iomatrix%e(1:3, i, i) = the_deriv(1:3, il)  
    end do

    deallocate(the_deriv)
    return
  800 continue
    write(0,*) "Memory allocation error in OnSiteDerivative"
    stop
  end subroutine

  !!
  ! Reads header part of on-site energy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadOnSiteHeader(iunit, omode)
    integer, intent(in) :: iunit
    integer, intent(out) :: omode

    integer :: the_err
    character(len=32) :: the_header
    read(iunit, '(A32)', iostat=the_err) the_header
    if (the_err /= 0) goto 800

    select case(the_header)
    case('ONSITE NRL')
      omode = NRLIO_NRL_ONSITE   
    case('ONSITE NRLR2')
      omode = NRLIO_NRLR2_ONSITE 
    case default
      goto 810
    end select 
 
    return

  800 continue
    print *, "There is no-header for on-site part."
    stop
  810 continue
    print *, "Onsite header "//trim(the_header)//" is not supported."
    stop
  end subroutine

  !!
  ! Read density for on-site energy
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadOnSiteDensity(ioonsite, ifile)
    type(TNrlOnSiteDensity), intent(inout) :: ioonsite
    type(TFileParam),     intent(in)    :: ifile

    real(8) :: r0, gamma, rc, delta
    real(8) :: lambda

    select case(ifile%mode)
    case(NRLIO_XB_ONSITE) 
      read(ifile%port, *) r0, gamma, rc, delta
      ioonsite%a     = 1.0d0
      ioonsite%r0    = r0
      ioonsite%gamma = gamma
      ioonsite%rc    = rc
      ioonsite%delta = delta
    case(NRLIO_NRL_ONSITE, NRLIO_NRLR2_ONSITE)
      read(ifile%port, *) lambda, rc, delta
      ioonsite%a     = exp(- lambda * lambda)
      ioonsite%r0    = 1.0d0
      ioonsite%gamma = lambda * lambda
      ioonsite%rc    = rc
      ioonsite%delta = delta
    case default
      write(0, *) "Not Supported Density Mode: ", ifile%mode
      stop
    end select
  end subroutine


  !!
  ! Read onsite energy shift parameters.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadOnSiteEnergy(ioonsite, ifile)
    type(TNrlOnSiteEnergy), intent(inout) :: ioonsite
    type(TFileParam),     intent(in)    :: ifile

    integer :: i

    do i = 1, ioonsite%orbital_count
      call ReadOnSiteEnergyParams(ioonsite%params(i), ifile)
    end do
  end subroutine

  !!
  ! Make onsite energy shift parameters 1.0.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetOnSiteEnergyIdentity(ioonsite)
    type(TNrlOnSiteEnergy), intent(inout) :: ioonsite

    integer :: i

    do i = 1, ioonsite%orbital_count
      call SetOnSiteEnergyParamsIdentity(ioonsite%params(i))
    end do
  end subroutine

  !!
  ! Prepares TNrlOnSiteBuilder object for reading file.
  ! Before calling ReadOnSiteParameters, you must call this routine.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine PrepareForRead(ioonsite, ispecies_count, iorbital_count)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer, intent(in) :: ispecies_count
    integer, intent(in) :: iorbital_count(:)

    call AllocateOnSiteBuilder(ioonsite, ispecies_count, iorbital_count)
  end subroutine

  !!
  ! Interface routine to read onsite parameters from the stream of given 
  ! file port.
  ! Before calling this routine, you must call PrepareForRead
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadFromUnit(ioonsite, iunit)  
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer,                 intent(in)    :: iunit
    integer :: the_mode

    type(TFileParam) :: file
    integer :: i, j

    ioonsite%mode = NRL_ONSITE_STANDARD

    call ReadOnSiteHeader(iunit, the_mode)

    file%port = iunit
    file%mode = the_mode

    do i = 1, ioonsite%species_count
      call ReadOnSiteEnergy(ioonsite%energies(i), file)
    end do

    do i = 1, ioonsite%species_count
      do j = 1, ioonsite%species_count
        call ReadOnSiteDensity(ioonsite%densities(i, j), file)
      end do
    end do

    !--skip one line--!
    read(iunit, *)
  end subroutine

  !!
  ! Interface routine to set identity parameter to onsite.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SetIdentity(ioonsite)
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer :: i

    ioonsite%mode = NRL_ONSITE_CONSTANT

    do i = 1, ioonsite%species_count
      call SetOnSiteEnergyIdentity(ioonsite%energies(i))
    end do
  end subroutine

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AllocateOnsiteXML(ioonsite, the_orbital_count, the_species_count)
    implicit none
    integer, intent(in) :: the_species_count
    integer, intent(in) :: the_orbital_count(the_species_count)
    integer i, j, ierr

    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite

    ioonsite%species_count = the_species_count
    
    ! --onsite
    if( associated(ioonsite%energies ) )then
       do i = 1, the_species_count
          deallocate( ioonsite%energies(i)%params )
       end do
       deallocate( ioonsite%energies )
    end if
    
    if( associated(ioonsite%densities ) ) &
         deallocate( ioonsite%densities )
    
    if( associated(ioonsite%density_cache ) ) &
         deallocate( ioonsite%density_cache )
    
    ! allocate
    ! --onsite
    allocate( ioonsite%energies(the_species_count), &
         ioonsite%densities(the_species_count, the_species_count), &
         stat=ierr)
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocation of onsite energies in sbrt allocate_onsiteXML'
    end if
    
    do i = 1, the_species_count
       
       ioonsite%energies(i)%orbital_count = 3  ! s,p,d (defalut) !!
   
       allocate( ioonsite%energies(i)%params(the_orbital_count(i)+1), stat=ierr) ! d=t2g,eg
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation of onsite params in sbrt allocate_onsiteXML'
       end if
       do j = 1, the_species_count
          ioonsite%densities(j,i)%a = 0.d0
          ioonsite%densities(j,i)%gamma = 0.d0
          ioonsite%densities(j,i)%delta = 1.d0
          ioonsite%densities(j,i)%rc = 5.9d0
          ioonsite%densities(j,i)%r0 = 1.d0
       end do
    end do
    
    return
  end subroutine AllocateOnsiteXML

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine OnsiteSetXML(ioonsite, lambda, cutoff, delta, is,js)
    implicit none
    integer, intent(in) :: is, js
    real(8), intent(in) :: lambda, cutoff, delta
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite

    ioonsite%densities(is,js)%rc = cutoff
    ioonsite%densities(is,js)%r0 = 1.d0
    ioonsite%densities(is,js)%delta = delta
    ioonsite%densities(is,js)%a = exp(-lambda*lambda)
    ioonsite%densities(is,js)%gamma = lambda*lambda

    return
  end subroutine OnsiteSetXML

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine OnsiteLoadXML(ioonsite,vnode,the_orbital_count,the_species_count,is)
    use flib_dom
    use elses_xml_misc
    implicit none
    type(fnode), target, intent(in) :: vnode
    type(fnode), pointer :: node
    integer, intent(in) :: is, the_species_count
    integer, intent(in) :: the_orbital_count(the_species_count)
    real(8) rvalue(4)
    character(len=256) value, unit, orbital
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    
    node => vnode

    unit  = getAttribute(node,"unit")
    if( unit == "" ) then
       unit = "a.u."
    endif
    
    value = getChildValue(node)
    read(unit=value,fmt=*) rvalue(1:4)
    rvalue(1:4) = 2.d0*rvalue(1:4)*XML_TO_AU(unit) ! Ry->a.u.->Ry
    orbital  = getAttribute(node,"orbital")
    select case(orbital)
    case("s")
       ioonsite%energies(is)%params(1)%a = rvalue(1)
       ioonsite%energies(is)%params(1)%b = rvalue(2)
       ioonsite%energies(is)%params(1)%c = rvalue(3)
       ioonsite%energies(is)%params(1)%d = rvalue(4)
       ioonsite%energies(is)%params(1)%e = 0.d0
       ioonsite%energies(is)%params(1)%rho_c = 100.d0
       ioonsite%energies(is)%params(1)%gamma = 1.d0
    case("p")
       ioonsite%energies(is)%params(2)%a = rvalue(1)
       ioonsite%energies(is)%params(2)%b = rvalue(2)
       ioonsite%energies(is)%params(2)%c = rvalue(3)
       ioonsite%energies(is)%params(2)%d = rvalue(4)
       ioonsite%energies(is)%params(2)%e = 0.d0
       ioonsite%energies(is)%params(2)%rho_c = 100.d0
       ioonsite%energies(is)%params(2)%gamma = 1.d0
    case("t2g")
       ioonsite%energies(is)%params(3)%a = rvalue(1)
       ioonsite%energies(is)%params(3)%b = rvalue(2)
       ioonsite%energies(is)%params(3)%c = rvalue(3)
       ioonsite%energies(is)%params(3)%d = rvalue(4)
       ioonsite%energies(is)%params(3)%e = 0.d0
       ioonsite%energies(is)%params(3)%rho_c = 100.d0
       ioonsite%energies(is)%params(3)%gamma = 1.d0
    case("eg")
       ioonsite%energies(is)%params(4)%a = rvalue(1)
       ioonsite%energies(is)%params(4)%b = rvalue(2)
       ioonsite%energies(is)%params(4)%c = rvalue(3)
       ioonsite%energies(is)%params(4)%d = rvalue(4)
       ioonsite%energies(is)%params(4)%e = 0.d0
       ioonsite%energies(is)%params(4)%rho_c = 100.d0
       ioonsite%energies(is)%params(4)%gamma = 1.d0
    end select
    return
  end subroutine OnsiteLoadXML
  
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine OnsiteSetIdentityXML(ioonsite, the_orbital_count, the_species_count, is)
    implicit none
    integer, intent(in) :: is 
    type(TNrlOnSiteBuilder), intent(inout) :: ioonsite
    integer, intent(in) :: the_species_count
    integer, intent(in) :: the_orbital_count(the_species_count)
    integer i

    do i = 1, the_orbital_count(is)
       ioonsite%energies(is)%params(i)%a = 1.d0
       ioonsite%energies(is)%params(i)%b = 0.d0
       ioonsite%energies(is)%params(i)%c = 0.d0
       ioonsite%energies(is)%params(i)%d = 0.d0
       ioonsite%energies(is)%params(i)%e = 0.d0
       ioonsite%energies(is)%params(i)%rho_c = 100.d0
       ioonsite%energies(is)%params(i)%gamma = 1.d0
    end do
    
    return
  end subroutine OnsiteSetIdentityXML
  
end module
