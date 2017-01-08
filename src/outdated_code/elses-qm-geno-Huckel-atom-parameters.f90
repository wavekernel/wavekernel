module M_qm_geno_Huckel_atom_params
  use M_qm_domain, only: i_verbose
  implicit none
  integer,parameter :: MAX_ANGULAR_M=2       ! Maximum number of atomic orbitals. 1+3+5, s,p,d
  integer,parameter :: NUMBER_OF_ORBITALS = (MAX_ANGULAR_M+1)*(MAX_ANGULAR_M+1)
  integer,parameter :: NUMBER_OF_TERM_VOIP=MAX_ANGULAR_M+1 ! d0 + d1 Q + d2 Q^2
  integer,parameter :: NUMBER_OF_CONF_VOIP=3 ! For d-orbital
  integer,parameter :: DOUBLE_PRECISION=kind(1d0)
  real(DOUBLE_PRECISION), parameter :: SMALL_N=1d-2

  ! Definition of a set of atomic parameters
  ! mainly consists of the generalized Huckel calculation.
  ! The mass and the chemical hardness are included in the set.
  type THuckelAtomParameters
     character(len=4):: name
     character(len=1):: orbital  (NUMBER_OF_ORBITALS) ! size=(# of orbitals)
     integer         :: principal(NUMBER_OF_ORBITALS) ! principal quantum number
     integer         :: angular  (NUMBER_OF_ORBITALS) ! angular momentum (l)
     integer         :: highest_occupied_orbital
     integer         :: num_val_orb    ! # of valence orbitals
     real(DOUBLE_PRECISION)    :: num_val_elec   ! # of valence electrons
     real(DOUBLE_PRECISION)    :: zeta     (NUMBER_OF_ORBITALS) ! STO screening parameter
     real(DOUBLE_PRECISION)    :: VOIP (NUMBER_OF_TERM_VOIP,NUMBER_OF_ORBITALS)

     ! If and only if orbital(k)=='d', then
     ! k-th component of the following variables are valid.
     ! For example, assuming that
     ! orbital(1)='s', orbital(2)='p', orbital(5)='d',
     ! then zeta2(5:9) is valid but zeta2(1:4) are invalid.
     real(DOUBLE_PRECISION)    :: zeta2    (NUMBER_OF_ORBITALS)
     real(DOUBLE_PRECISION)    :: c1       (NUMBER_OF_ORBITALS)
     real(DOUBLE_PRECISION)    :: c2       (NUMBER_OF_ORBITALS)
     real(DOUBLE_PRECISION)    :: &
          VOIP_D(NUMBER_OF_TERM_VOIP,NUMBER_OF_CONF_VOIP,NUMBER_OF_ORBITALS)
     ! Rescaling factor & rescaled quantities used in repulsive part calc.
     real(DOUBLE_PRECISION),dimension(NUMBER_OF_ORBITALS) &
          :: repulsive_rescaling, &
          rescaled_zeta, rescaled_zeta2, rescaled_c1, rescaled_c2
     ! end ('d' specific parameter)

     ! extra parameters
     real(DOUBLE_PRECISION)    :: mass
     real(DOUBLE_PRECISION)    :: chemical_hardness
                          ! nealy = (Hubbard U of valence electron)
     ! Mass and chemical hardness is not a parameter used in the generalized Huckel formalism.
     ! But it is convenient to treat them here, because most other atom dependent parameters
     ! are set here.
     real(DOUBLE_PRECISION)    :: initial_charge     !(default 0)
     real(DOUBLE_PRECISION)    :: initial_occupation(NUMBER_OF_ORBITALS)
     !  = b^0, determined by "Aufbau Principle" (old version)  -> 
     !                       initial_diagonal_elements (new version)
     real(DOUBLE_PRECISION)    :: initial_diagonal_elements(NUMBER_OF_ORBITALS)
     ! (averaged Coulomb interaction in the free atom $U$ in table (2.2))
     ! Note for ICONC compatibility
     !      The default atom in ICONC is neutral.
     !      In the charge self-consistency mode,  atomic charge (Q)
     !     (s, p, d separately, when d-orbital exists) is calculated as
     !      Mulliken charge, then the $H_{ii}$ = -VOIP(Q).
     !      In the default (charge non-self-consisteency) mode,
     !      Coulomb parameter $U$ in the table 2.2 is used as $H_{ii}$, instead.
     ! 
     !      We choose $H_{ii}$ as follows. If non-zero initial charge
     !     (s, p, d separately, when d-orbital exists) is given,
     !      $H_{ii}$ = -VOIP(Q). If initial charge == 0, $H_{ii}$=$U$, instead.

  end type THuckelAtomParameters

  type(THuckelAtomParameters),allocatable:: AtomParameter(:)
  ! allocate(AtomParameter(NOS))
  !How to access respective parameter?
  ! For example, if you want to write the chemical_hardness
  ! of the d-orbital on the second atom specie,
  ! write(*,*) AtomParameter(2)%chemical_hardness

  logical :: UseVOIP, EmulateICON

contains

  subroutine GetChemicalHardness(chemical_hardness)
    use M_qm_domain
    real(DOUBLE_PRECISION), intent(out) :: chemical_hardness(:)
    integer :: jsv, nss
    if(size(chemical_hardness,1) /= noav) &
       stop 'GetChemicalHardness: Inconsistent size (chemical_hardness)'
    do jsv=1, noav
       nss=atm_element(jsv)
       chemical_hardness(jsv)=Atomparameter(nss)%chemical_hardness
    end do
  end subroutine GetChemicalHardness

  subroutine GetInitialENum(initial_e_num)
    use M_qm_domain
    real(DOUBLE_PRECISION), intent(out) :: initial_e_num(:)
    integer :: jsv, nss
    !write(*,*)  noav
    !write(*,*)  size(initial_e_num,1)
    if(size(initial_e_num,1) /= noav) &
       stop 'GetInitialENum: Inconsistent size (initial_e_num)'
    do jsv=1, noav
       nss=atm_element(jsv)
       initial_e_num(jsv)=Atomparameter(nss)%num_val_elec-Atomparameter(nss)%initial_charge
    end do
  end subroutine GetInitialENum

  subroutine GetInitialOcc(initial_occ)
    use M_qm_domain
    real(DOUBLE_PRECISION), intent(out) :: initial_occ(:,:)
    integer :: jsv, nss, iorb, max_nval
    if(size(initial_occ,2) /= noav) &
       stop 'GetInitialOcc: Inconsistent size (initial_occ(,:))'
    max_nval=maxval(nval)
    if(size(initial_occ,1) /= max_nval) &
       stop 'GetInitialOcc: Inconsistent size (initial_occ(:,))'
    do jsv=1, noav
       nss=atm_element(jsv)
       do iorb=1, nval(nss)
          initial_occ(iorb,jsv)=Atomparameter(nss)%initial_occupation(iorb)
       end do
       do iorb=iorb, max_nval
          initial_occ(iorb,jsv)=0d0
       end do
    end do
  end subroutine GetInitialOcc

  subroutine SetAtomParameters(NOS,elem_name,MASS,num_val_orb,num_val_elec)
    use elses_mod_phys_const
    use flib_dom
    use elses_xml_misc
    use elses_xml_config, only: config

    integer         :: NOS,     num_val_orb(:)
    real(DOUBLE_PRECISION)    :: MASS(:), num_val_elec(:)
    character(len=*):: elem_name(:)
    intent(in)      :: NOS, elem_name
    intent(out)     :: num_val_orb, num_val_elec, MASS
    logical :: Load_XML
    integer :: i, j
    real(DOUBLE_PRECISION):: dtmp

    !checking the consistency of the arguments
    if(NOS /= size(MASS)) &
         stop 'SetAtomParameters: inconsistent size of MASS'
    if(NOS /= size(elem_name)) &
         stop 'SetAtomParameters: inconsistent size of elem_name'

    allocate(AtomParameter(NOS))

    do i=1, NOS
       AtomParameter(i)%principal    =0  ! (uninitialized value)
       AtomParameter(i)%num_val_orb  =0  ! (uninitialized value)
       AtomParameter(i)%num_val_elec =0d0
       AtomParameter(i)%VOIP         =0d0
       AtomParameter(i)%VOIP_D       =0d0
       AtomParameter(i)%zeta         =0d0
       AtomParameter(i)%zeta2        =0d0
       AtomParameter(i)%c1           =0d0
       AtomParameter(i)%c2           =0d0
       AtomParameter(i)%repulsive_rescaling      = 1d0
       AtomParameter(i)%initial_charge           = 0d0
       AtomParameter(i)%initial_occupation       =-1d0       ! (uninitialized value)
       AtomParameter(i)%initial_diagonal_elements= Huge(0d0) ! (uninitialized value)
    end do

    EmulateICON=config%calc%genoOption%Emulate_ICON
    UseVOIP    =config%calc%genoOption%Use_VOIP
    if( i_verbose >= 1 ) then
       write(*,'("SetAtomParametes: EmulateICON = ",L1,", UseVOIP = ",L1)') EmulateICON, UseVOIP
    end if

    ! set atom parameters
    do i=1,NOS
       Load_XML=.false.
       if( config%system%structure%velement(i)%filename /= "" ) then
          Load_XML=Load_XML_AtomParameter( AtomParameter(i),&
               config%system%structure%velement(i)%filename )
       end if
       ! *** tentative implementation ***
       if( .not. Load_XML ) call load_default(AtomParameter(i),elem_name(i))
       call setDefaultElecConf(AtomParameter(i),elem_name(i))

       if( .not. EmulateICON ) call renormalizeDoubleZeta(AtomParameter(i))

      ! setting chemical hardness (*** tentative value ***)
      dtmp=-Huge(0d0) ! candidate of HOMO level
      do j=1, NUMBER_OF_ORBITALS
         if(AtomParameter(i)%initial_diagonal_elements(j) /= Huge(0d0) &
              .and.  AtomParameter(i)%initial_occupation(j) > 0d0 &
              .and.  AtomParameter(i)%initial_diagonal_elements(j) > dtmp ) then
            dtmp  =  AtomParameter(i)%initial_diagonal_elements(j)
            AtomParameter(i)%highest_occupied_orbital=j
         end if
      end do
      AtomParameter(i)%chemical_hardness = - dtmp
      ! Above expression of chemical hardness is nonsensical, but we don't know other estimation.
      ! Definition of chemical hardness is as follows.
      !    d^2 / (dN)^2 E is nearly equal to 1/2{[E(N+1) - E(N)]-[E(N)-E(N-1)]}= 1/2 {[-A] - [-I]}
      ! Here, A and I stands for electron affinity and ionization energy, respectively.
      ! Unfortunately, we have few hints for estimating A and I -- single-electron atomic levels.
      ! Therefore, we use naive approximation A = (LUMO level), I = (HOMO level).
      ! But these approximation is totally nonsense, because electron-electron interaction is omitted
      ! though we would like to know electron-electron interaction.
      ! Consequently, if above naive approximation is applied to open shell system,
      ! the value of chemical hardness is equal to zero!
      ! 

      ! rescaling of repulsive part
      if(NUMBER_OF_ORBITALS /=9 ) stop 'SetAtomParameters: cannot rescale'
      if(AtomParameter(i)%num_val_orb > 4) then
         AtomParameter(i)%rescaled_zeta(1:4)= &
              AtomParameter(i)%repulsive_rescaling(1:4)*&
              AtomParameter(i)%zeta(1:4)
         AtomParameter(i)%rescaled_zeta(5:9)= &
              AtomParameter(i)%zeta(5:9)
         AtomParameter(i)%rescaled_zeta2(5:9)= &
              AtomParameter(i)%repulsive_rescaling(5:9)*&
              AtomParameter(i)%zeta2(5:9)
         AtomParameter(i)%rescaled_zeta2(1:4)=0d0
         AtomParameter(i)%rescaled_c1(1:4)   =0d0
         AtomParameter(i)%rescaled_c2(1:4)   =0d0
         do j=5, 9
            dtmp=normDoubleZeta(&
                 AtomParameter(i)%principal(j),&
                 AtomParameter(i)%rescaled_zeta(j),&
                 AtomParameter(i)%rescaled_zeta2(j),&
                 AtomParameter(i)%c1(j),&
                 AtomParameter(i)%c2(j))
            AtomParameter(i)%rescaled_c1(j)=&
                 AtomParameter(i)%c1(j)/dtmp
            AtomParameter(i)%rescaled_c2(j)=&
                 AtomParameter(i)%c2(j)/dtmp
         end do
      else
         AtomParameter(i)%rescaled_zeta= &
              AtomParameter(i)%repulsive_rescaling*&
              AtomParameter(i)%zeta
      end if

      if( i_verbose >=1 ) call write_information(AtomParameter(i))
    end do

    MASS(:)=AtomParameter(:)%mass
    num_val_orb (:)=AtomParameter(:)%num_val_orb
    num_val_elec(:)=AtomParameter(:)%num_val_elec

    if( i_verbose >= 200 )then
       write(*,'("SetAtomParametes: ")')
       write(*,'(" size of MASS     =",I4)') size(MASS)
       write(*,'(" size of elem_name=",I4)') size(elem_name)
       write(*,'(" MASS        (:)=",ES22.14)') MASS(:)
       write(*,'(" num_val_orb (:)=",I3)')      num_val_orb(:)
       write(*,'(" num_val_elec(:)=",ES22.14)') num_val_elec(:)
    end if

  contains

    function Load_XML_AtomParameter( AtomParameter, filename ) result(succeeded)
      implicit none
      type(THuckelAtomParameters), intent(out) :: AtomParameter
      character(len=*), intent(in) :: filename ! element file name
      logical :: succeeded

      type(fnode), pointer     :: document_node
      type(fnode), pointer     :: element_node
      type(fnode), pointer     :: node, sub_node
      character(len=256)       :: value
      real(DOUBLE_PRECISION),allocatable :: linearArray(:)
      integer :: nsize(3)

      inquire( file=filename, exist=succeeded )

      if( .not. succeeded ) then
         write(*,'("Load_XML_AtomParameter: *** Warning!!! *** Failed to open data file ",A)') trim(filename)
         return
      end if

      document_node => parsefile(filename)
      call normalize(document_node)

      ! get <element> node
      element_node => getFirstElementByTagName(document_node,"element")
      if( .not. associated(element_node) ) then
         call XML_error("<element> not found")
      endif

      ! get type attribute
      value = getAttribute(element_node,"type")
      if( value /= "Huckel" ) then
         call XML_error("<element type> is not Huckel")
      endif

      ! get name attribute
      value = getAttribute(element_node,"name")
      if( value == "" ) then
         call XML_error("<element name> not found")
      endif
      read(unit=value,fmt=*) AtomParameter%name

      ! get <mass> node
      node => getFirstElementByTagName(element_node,"mass")
      if( .not. associated(node) ) then
         AtomParameter%mass = 0
      else
         value = getChildValue(node)
         read(unit=value,fmt=*) AtomParameter%mass
      endif

      ! get <chemical_hardness> node
      node => getFirstElementByTagName(element_node,"chemical_hardness")
      if( .not. associated(node) ) then
         AtomParameter%chemical_hardness = 0.d0
      else
         value = getChildValue(node)
         read(unit=value,fmt=*) AtomParameter%chemical_hardness
      endif

      ! get <initial_charge> node
      node => getFirstElementByTagName(element_node,"initial_charge")
      if( .not. associated(node) ) then
         AtomParameter%initial_charge = 0
      else
         value = getChildValue(node)
         read(unit=value,fmt=*) AtomParameter%initial_charge
      endif

      ! get <principal_quantum_number> node
      node => getFirstElementByTagName(element_node,"principal_quantum_number")
      if( associated(node) ) then
         call load_iArray( AtomParameter%principal, node )
      end if

      ! get <initial_diagonal_elements> node
      node => getFirstElementByTagName(element_node,"initial_diagonal_elements")
      if( associated(node) ) then
         call load_dArray( AtomParameter%initial_diagonal_elements, node )
         where(AtomParameter%initial_diagonal_elements /= Huge(0d0))
            AtomParameter%initial_diagonal_elements = &
                 AtomParameter%initial_diagonal_elements/eV4au
         end where
      end if

      ! get <initial_occupation> node
      node => getFirstElementByTagName(element_node,"initial_occupation")
      if( associated(node) ) then ! If no occupation is supplied in the XML file, use default values.
         if( getChildValue(node) == "default" )then
            stop
         end if
         call load_dArray( AtomParameter%initial_occupation, node )
      end if

      ! get <zeta> node
      node => getFirstElementByTagName(element_node,"zeta")
      if( .not. associated(node) ) then
         AtomParameter%zeta = 0.d0
      else
         call load_dArray( AtomParameter%zeta, node )
      end if

      ! get <zeta2> node
      node => getFirstElementByTagName(element_node,"zeta2")
      if( associated(node) ) then
         call load_dArray( AtomParameter%zeta2, node )
      end if

      ! get <c1> node
      node => getFirstElementByTagName(element_node,"c1")
      if( associated(node) ) then
         call load_dArray( AtomParameter%c1, node )
      end if

      ! get <c2> node
      node => getFirstElementByTagName(element_node,"c2")
      if( associated(node) ) then
         call load_dArray( AtomParameter%c2, node )
      end if

      ! get <VOIP> node
      node => getFirstElementByTagName(element_node,"VOIP")
      if( associated(node) ) then
         nsize(1:2)=(/size(AtomParameter%VOIP,1),size(AtomParameter%VOIP,2)/)
         allocate(linearArray(nsize(1)*nsize(2)))
         call load_dArray( linearArray, node )
         AtomParameter%VOIP(:,:) = reshape(linearArray,(/nsize(1),nsize(2)/))/eV4au
         deallocate(linearArray)
      end if

      ! get <VOIP_D> node
      node => getFirstElementByTagName(element_node,"VOIP_D")
      if( associated(node) ) then
         nsize(1:3)=(/&
              size(AtomParameter%VOIP_D,1),&
              size(AtomParameter%VOIP_D,2),&
              size(AtomParameter%VOIP_D,3)/)
         allocate(linearArray(nsize(1)*nsize(2)*nsize(3)))
         call load_dArray( linearArray, node )
         AtomParameter%VOIP_D(:,:,:) = &
              reshape(linearArray,(/nsize(1),nsize(2),nsize(3)/))/eV4au
         deallocate(linearArray)
      end if

      ! get <repulsive_rescaling> node
      node => getFirstElementByTagName(element_node,"repulsive_rescaling")
      if( associated(node) ) then
         call load_dArray( AtomParameter%repulsive_rescaling, node )
      end if

      call destroyNode(document_node)

    end function Load_XML_AtomParameter

    subroutine replaceAll(ch,ioString)
      character(len=*),intent(inout) :: ioString
      character(len=1),intent(in) :: ch
      integer :: k
      k=index(ioString,ch)
      do while( k > 0  )
         ioString( k:k ) = " "
         k=index(ioString,ch)
      end do
    end subroutine replaceAll

    subroutine normalizeString( string )
      character(len=*),intent(inout) :: string
      if(index('<',string)>0) stop &
           "normalizeString: Extra tag exists. &
           &Possibly obsoletely-formatted atom file."
      call replaceAll(achar(13),string) ! Carriage Retern
      call replaceAll(achar(10),string) ! Line Feed
    end subroutine normalizeString

    subroutine load_cArray( array, node )
      implicit none
      character(len=*), dimension(:), intent(out) :: array
      type(fnode), pointer :: node
      character(len=256)   :: value
      if( associated(node) ) then
         value = getChildValue(node)
         call normalizeString(value)
         read(unit=value, fmt=*, END=999) array
999      continue
      end if
    end subroutine load_cArray

    subroutine load_iArray( array, node )
      implicit none
      integer, dimension(:), intent(out) :: array
      type(fnode), pointer :: node
      character(len=256)   :: value
      if( associated(node) ) then
         value = getChildValue(node)
         call normalizeString(value)
         read(unit=value, fmt=*, END=999) array
999      continue
      end if
    end subroutine load_iArray

    subroutine load_dArray( array, node )
      implicit none
      real(DOUBLE_PRECISION), dimension(:), intent(out) :: array
      type(fnode), pointer :: node
      character(len=256)   :: value
      if( associated(node) ) then
         value = getChildValue(node)
         call normalizeString(value)
         read(unit=value, fmt=*, END=999) array
999      continue
      end if
    end subroutine load_dArray

    subroutine setDefaultElecConf(ioAtomParameters,elemName)
      type(THuckelAtomParameters), intent(inout) :: ioAtomParameters
      character(len=*), intent(in) :: elemName
      integer, parameter :: NUM_SPEC=103
      integer, parameter :: AUFBAU_MAX_L=6
      character(len=2) :: Table(NUM_SPEC)=(/ &
      'H ',                              'He',&
      'Li','Be','B ','C ','N ','O ','F ','Ne',&
      'Na','Mg','Al','Si','P ','S ','Cl','Ar',&
      'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',&
      'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',&
      'Cs','Ba',  'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',&
      'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',&
      'Fr','Ra',  'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' /)
      character(len=1) :: character(0:MAX_ANGULAR_M)=(/'s','p','d'/)
      integer :: i, j, k, n, atomic_number, valence(0:AUFBAU_MAX_L) ! higher shell must be closed.
      integer :: maxL, nValElec
      logical :: found
      found = .false.
      do i=1, NUM_SPEC     ! first try for longest match
         if(elemName(1:2)==Table(i)) then
            found = .true.
            atomic_number=i
            exit
         end if
      end do
      if(.not. found)then ! second try
         do i=1, NUM_SPEC
            if(elemName(1:1)==Table(i)) then
               found = .true.
               atomic_number=i
               exit
            end if
         end do
      end if
      if(.not. found) stop 'setDefaultElecConf: not found'

      ! Aufbau principle
      valence=0
      maxL=0
      k=0
      outer: do i=1, AUFBAU_MAX_L+1   ! i=n+l
         inner: do j=(i+1)/2-1, 0, -1 ! j=l, n=i-j
            n=i-j
            valence(j)=n
            k=k+(2*j+1)
            maxL=max(maxL,j)
            if(2*k >= atomic_number) exit outer
         end do inner
      end do outer

      if( valence(4) > 0 )then
         stop "setDefaultElecConf: Not supported (l > 3)"
      else if( valence(3) > 0 ) then
         if( j == maxL ) then
            if( 2*k > atomic_number ) stop "setDefaultElecConf: open f-shell"
            if( 2*k < atomic_number ) stop "setDefaultElecConf: logic error 1"
         else
            if( maxL /= 3 ) stop "setDefaultElecConf: logic error 2"
            write(*,'("setDefaultElecConf: closed f-shell exist")') 
         end if
         maxL = 2
      end if

      nValElec = atomic_number-2*k ! -(# of vacant orbit)
      if( valence(1) < valence(0) .and. valence(1) > 0 ) then
                                         ! Def. p.q.n : principal quantum number
         valence(1) = valence(0)         ! If  p.q.n of s-orbit > p.q.n of p-orbit,
         nValElec = nValElec - 2*(2*1+1) ! p-orbit is treated as core-orbit.
         if(  ioAtomParameters%initial_diagonal_elements(1) > &
              ioAtomParameters%initial_diagonal_elements(2) ) then
            stop "setDefaultElecConf: unexpected condition. " ! Parameter check
         end if
      end if
      do i=0, maxL
         nValElec = nValElec + 2*(2*i+1)
      end do

      ! If num_val_orb is not initialized, count an availabe number of diagonal elements
      ! and set this number to num_val_orb.
      n = ioAtomParameters%num_val_orb
      if( n == 0 ) then
         do i=1, NUMBER_OF_ORBITALS
            ! If initial_diagonal_elements are initialized, that orbital is available.
            ! Huge(0d0) : uninitialized value
            if(ioAtomParameters%initial_diagonal_elements(i) /= Huge(0d0) )&
               n = n+1
         end do
         ioAtomParameters%num_val_orb=n
      end if

      j=int(sqrt(n+SMALL_N))-1
      if( j > maxL )then
         write(*,'("setDefaultElecConf: (input max L) > (max L after aufbau principle)")')
         maxL=j
      else if ( j < maxL ) then
         write(*,'("setDefaultElecConf: (input max L) < (max L after aufbau principle)")')
         write(*,'(" There may be missing values in your configuration files.")')
         if( maxL == 2 .and. j==1 .and. nValElec > 10) then
            write(*,'(" Try to continue calculation assuming closed d-shell exists.")')
            nValElec = nValElec - 2*(2*maxL+1)
            maxL=1
         else
            stop "setDefaultElecConf: Can't continue."
         end if
      end if

      ! p.q.n, angular mom., orbital character label
      do i=0, maxL
         ! initialize p.q.n, if it is not initialized before.
         if(  ioAtomParameters%principal(i*i+1) == 0 ) then
            if( valence(i) == 0 ) stop "setDefaultElecConf: <principal_quantum_number> may be absent."
            ioAtomParameters%principal(i*i+1:(i+1)*(i+1)) = valence(i)
         end if
         do j=i*i+1, (i+1)*(i+1)
            ioAtomParameters%angular(j)=i
            ioAtomParameters%orbital(j)=character(i)
         end do
      end do

      if( sum(ioAtomParameters%initial_occupation(1:n)) < 0d0 ) &
           call setDefaultOccupation(ioAtomParameters, nValElec)
      ! set numbers of valence electrons
      ioAtomParameters%num_val_elec=sum(ioAtomParameters%initial_occupation(1:n))

      if( i_verbose >= 200 )then
         write(*,'("setDefaultElecConf:")')
         write(*,'(" Atomic Number = ",I3)') atomic_number
         write(*,'(" Atomic Symbol = ",A2)') Table(atomic_number)
      end if
    end subroutine setDefaultElecConf

    subroutine setDefaultOccupation(ioAtomParameters, nValElec)
      type(THuckelAtomParameters), intent(inout) :: ioAtomParameters
      integer, intent(in) :: nValElec
      integer :: i, j, k, m, n
      integer ::                order(MAX_ANGULAR_M+1)
      real(DOUBLE_PRECISION) :: level(MAX_ANGULAR_M+1)

      n=ioAtomParameters%num_val_orb
      k=int(sqrt(n+SMALL_N))
      if( n /= k*k ) stop 'setDefaultOccupation: Unsupported num_val_orb'
      ! k=maxL+1

      ioAtomParameters%initial_occupation(1:n)=0d0
      ! initialize initial_occupation of relvant orbit

      do i=1, k
         level(i)=ioAtomParameters%initial_diagonal_elements((i-1)*(i-1)+1)
      end do

      call dhsort(k,level,order)

      n=0
      do i=1, k
         m=order(i)-1 ! L
         if(n+2*(2*m+1) <= nValElec-ioAtomParameters%initial_charge) then
            do j=m*m+1, (m+1)*(m+1)
               ioAtomParameters%initial_occupation(j) = 2d0
            end do
            n=n+2*(2*m+1)
         else
            do j=m*m+1, (m+1)*(m+1)
               ioAtomParameters%initial_occupation(j) = &
                    (nValElec-ioAtomParameters%initial_charge-n)/(2*m+1)
            end do
            exit
         end if
      end do
        
    end subroutine SetDefaultOccupation

    function normDoubleZeta(n,zeta,zeta2,c1,c2)
      integer :: n
      real(DOUBLE_PRECISION) :: normDoubleZeta, zeta, zeta2, c1, c2
      intent(in) :: n, zeta, zeta2, c1, c2
      normDoubleZeta=(4*zeta*zeta2/((zeta+zeta2)*(zeta+zeta2)))**(n+0.5d0)
      normDoubleZeta=sqrt(c1*c1+c2*c2+2*c1*c2*normDoubleZeta)
    end function normDoubleZeta

    subroutine renormalizeDoubleZeta(ioAtomParameters)
      type(THuckelAtomParameters), intent(inout) :: ioAtomParameters
      real(DOUBLE_PRECISION) :: f, zeta, zeta2, c1, c2
      integer :: i, j, n
!!! Aug8,2009 SY. Renormalization of c1, c2 : To avoid minus repulsive energy.
!!! See the expression of RepulsiveEnergy carefully and find ill-normalized
!!! c1 & c2 cause erroneous repulsive energy.
      i=2 ! L of d-orbital
      if( ioAtomParameters%num_val_orb < (i+1)*(i+1) ) return
      ! Hereafter, we assume d-orbital exists.
      do j=i*i+1, (i+1)*(i+1)
         n    =ioAtomParameters%principal(j)
         zeta =ioAtomParameters%zeta(j)
         zeta2=ioAtomParameters%zeta2(j)
         c1   =ioAtomParameters%c1(j)
         c2   =ioAtomParameters%c2(j)
         f = normDoubleZeta(n,zeta,zeta2,c1,c2)
         if(abs(f - 1d0) > 1d-2) then
            write(*,'("n    =",I1)') n
            write(*,'("c1,c2=",2ES22.14)') c1,c2
            write(*,'("z1,z2=",2ES22.14)') zeta, zeta2
            write(*,'("sqrt(z1.z2)     =", ES22.14)')f
            write(*,'("Suggested c1,c2 =",2ES22.14)')c1/f,c2/f
            stop &
           'renormalizeDoubleZeta in SetAtomParameters:&
           & erroneous normalization of double zeta'
!         else if(abs(f - 1d0) < 1d-12)
!            return
         end if
         ioAtomParameters%c1(j)=c1/f
         ioAtomParameters%c2(j)=c2/f
      end do
    end subroutine renormalizeDoubleZeta
    
    subroutine setStandard(ioAtomParameters,&
      mass,zeta,initial_diagonal_elements,&
         VOIP,zeta2,c1,c2,VOIP_D)
      type(THuckelAtomParameters), intent(inout) :: ioAtomParameters
      real(DOUBLE_PRECISION) :: mass, zeta(:), initial_diagonal_elements(:), zeta2, c1, c2
      real(DOUBLE_PRECISION) :: VOIP(:,:), VOIP_D(:,:,:)
      intent(in) :: mass,zeta,initial_diagonal_elements,VOIP,zeta2,c1,c2,VOIP_D
      optional   :: VOIP,zeta2,c1,c2,VOIP_D
      integer    :: maxL, i, j, k

      ! set mass & maxL
      ioAtomParameters%mass=mass
      maxL=size(initial_diagonal_elements,1)-1

      ! parameter check
      if( UseVOIP .and. .not. (present(VOIP) .or. present(VOIP_D)) ) stop &
           'setStandard in setAtomParameters : UserVOIP but no VOIP parameters are given.'
      if( present(VOIP) .and. present(VOIP_D) ) stop &
           'setStandard in setAtomParameters : VOIP & VOIP_D are not compatible.'
      if( maxL==2 .and. present(VOIP) ) stop &
           'setStandard in setAtomParameters : &
           &d-orbital presents but VOIP parameters for sp-system are given.'
      if( maxL==2 .and. .not. present(zeta2) ) stop &
           'setStandard  in setAtomParameters: d-orbital presents but no zeta2.'
      if( maxL+1 /= size(zeta,1) ) stop 'setStandard: wrong size of zeta.'
      if( maxL+1 /= size(initial_diagonal_elements,1) ) stop &
           'setStandard  in setAtomParameters: wrong size of initial_diagonal_elements.'

      ! set the other parameters
      do i=0, maxL
         k=i+1
         do j=i*i+1, (i+1)*(i+1)
            ioAtomParameters%zeta(j)=zeta(k)
            ioAtomParameters%initial_diagonal_elements(j)=initial_diagonal_elements(k)/eV4au
            if(present(VOIP))ioAtomParameters%VOIP(:,j)=VOIP(:,k)/eV4au
         end do
      end do
      if(present(zeta2))then
         if(.not. present(c1)) stop 'setStandard in setAtomParameters : zeta2 without c1'
         if(.not. present(c2)) stop 'setStandard in setAtomParameters : zeta2 without c2'
         i=2
         do j=i*i+1, (i+1)*(i+1)
            ioAtomParameters%zeta2(j)=zeta2
            ioAtomParameters%c1(j)=c1
            ioAtomParameters%c2(j)=c2
         end do
         if(present(VOIP_D))then
            do j=i*i+1, (i+1)*(i+1)
               ioAtomParameters%VOIP_D(:,:,j)=VOIP_D(:,:,i)
            end do
         end if
      end if
      if(present(VOIP_D) .and. .not. present(zeta2)) &
           stop'setStandard in setAtomParameters : VOIP_D without zeta2'

    end subroutine setStandard

    subroutine load_default(ioAtomParameters,elem_name)
      character(len=*),intent(in)                :: elem_name
      type(THuckelAtomParameters), intent(inout) :: ioAtomParameters
      select case(elem_name)
      case('H')
         call setStandard(ioAtomParameters, 1.0079d0,&
              (/1.3d0/),&
              (/-13.6d0/),&
              RESHAPE((/ (/ 13.618d0, 27.18d0, 13.6d0 /) /),(/3,1/)) )
      case('Li')
         call setStandard(ioAtomParameters, 6.941d0, &
              (/0.650d0,0.650d0/),&
              (/-5.4d0,-3.5d0/),&
              RESHAPE((/ (/ 3.42440d0,9.42287d0,5.38095d0 /), &
                         (/ 3.43934d0,7.32752d0,3.52117d0 /) /),(/3,2/) ) )
      case('Be')
         call setStandard(ioAtomParameters, 9.012182d0, &
              (/0.650d0,0.650d0/),&
              (/-10.0d0,-6.0d0/),&
              RESHAPE((/ (/ 3.4269510d0,12.4357138d0,9.3112879d0 /), &
                         (/ 3.4418290d0,10.0799960d0,5.9388910d0 /) /),(/3,2/)) )
      case('C')
         call setStandard(ioAtomParameters, 12.0107d0, &
              (/1.710d0,1.625d0/),&
              (/-21.4d0,-11.4d0/),&
              RESHAPE((/ (/ 3.4653860d0,17.5563030d0,19.4160805d0 /), &
                         (/ 3.4653860d0,14.6550490d0,10.6379290d0 /) /),(/3,2/)) )
      case('N')
         call setStandard(ioAtomParameters, 14.0067d0, &
              (/2.140d0,1.950d0/),&
              (/-26.0d0,-13.4d0/),&
              RESHAPE((/ (/ 3.4914230d0,20.1103970d0,25.5657463d0 /), &
                         (/ 3.4914230d0,16.5148280d0,13.1920233d0 /) /),(/3,2/)) )
      case('O')
         call setStandard(ioAtomParameters, 15.9994d0, &
              (/2.575d0,2.275d0/),&
              (/-28.2d0,-12.4d0/),&
              RESHAPE((/ (/ 3.4653860d0,22.8876650d0,32.3353348d0 /), &
                         (/ 3.4641460d0,18.5667820d0,15.7957125d0 /) /),(/3,2/)) )
      case('F')
         call setStandard(ioAtomParameters, 18.9984032d0, &
              (/3.010d0,2.425d0/),&
              (/-40.0d0,-18.1d0/),&
              RESHAPE((/ (/ 3.4802640d0,25.5037530d0,40.1216049d0 /), &
                         (/ 3.4629060d0,20.5195480d0,18.6473713d0 /) /),(/3,2/)) )
      case('Na')
         call setStandard(ioAtomParameters, 22.98977d0, &
              (/0.733d0,0.733d0/),&
              (/-5.1d0,-3.0d0/),&
              RESHAPE((/ (/ 1.6341248d0,8.4309931d0,5.0833931d0 /), &
                         (/ 1.6527227d0,6.1248687d0,2.9632462d0 /) /),(/3,2/)) )
      case('Mg')
         call setStandard(ioAtomParameters, 24.305d0, &
              (/1.100d0,1.100d0/),&
              (/-26.0d0,-13.4d0/),&
              RESHAPE((/ (/ 1.6279255d0,9.6956415d0,7.5940928d0 /), &
                         (/ 1.6291650d0,9.7452360d0,8.8893860d0 /) /),(/3,2/)) )
      case('Al')
         call setStandard(ioAtomParameters, 26.981538d0, &
              (/1.167d0,1.167d0/),&
              (/-12.3d0,-6.5d0/),&
              RESHAPE((/ (/ 1.6428040d0,11.0346820d0,11.2578545d0 /), &
                         (/ 1.6477633d0,8.8153467d0,5.9326911d0 /) /),(/3,2/)) )
      case('Al3')
         call setStandard(ioAtomParameters, 26.981538d0, &
              (/1.550d0,1.550d0/),&
              (/-16.31d0,-10.000d0/))
      case('Si')
         call setStandard(ioAtomParameters, 28.0855d0, &
              (/1.383d0,1.383d0/),&
              (/-17.3d0,-9.2d0/),&
              RESHAPE((/ (/ 1.6217260d0,12.3861200d0,14.8286285d0 /), &
                         (/ 1.6142870d0,10.1295900d0,7.7490745d0 /) /),(/3,2/)) )
      case('P')
         call setStandard(ioAtomParameters, 30.973761d0, &
              (/1.600d0,1.300d0/),&
              (/-18.6d0,-14.0d0/),&
              RESHAPE((/ (/ 1.7692690d0,13.2292190d0,18.7713566d0 /), &
                         (/ 1.8907740d0,10.4023570d0,10.1171913d0 /) /),(/3,2/)) )
      case('S')
         call setStandard(ioAtomParameters, 32.065d0, &
              (/2.283d0,1.817d0/),&
              (/-20.5d0,-11.4d0/),&
              RESHAPE((/ (/ 1.5163390d0,15.3741640d0,20.6683311d0 /), &
                         (/ 1.6328850d0,12.2125420d0,11.5802164d0 /) /),(/3,2/)) )
      case('Cl')
         call setStandard(ioAtomParameters, 35.453d0, &
              (/2.183d0,1.733d0/),&
              (/-36.3d0,-14.2d0/),&
              RESHAPE((/ (/ 1.6985970d0,15.7089230d0,25.2681808d0 /), &
                         (/ 1.6725600d0,13.1796260d0,13.6879644d0 /) /),(/3,2/)) )
      case('Fe')
         call setStandard(ioAtomParameters, 55.845d0, &
              (/1.900d0,1.900d0,5.350d0/),&
              (/-9.10d0,-5.32d0,-12.60d0/),&
              zeta2=2.000d0, c1=.550d0, c2=.626d0,&
              VOIP_D = reshape(&
              (/ (/ 0.9110000d0,7.9160000d0,7.1040000d0, &
                    0.9110000d0,9.0570000d0,8.4680000d0, &
                    0.9110000d0,8.3500000d0,10.092000d0 /), &
                 (/ 0.9050000d0,6.2980000d0, 3.707000d0, &
                    0.9050000d0,7.1660000d0, 4.922000d0, &
                    0.9050000d0,7.1660000d0, 4.996000d0 /), &
                 (/ 1.711000d0,10.687000d0, 5.195000d0, &
                    1.711000d0,12.584000d0, 8.678000d0, &
                    1.711000d0,12.634000d0,10.067000d0 /) /), (/3,3,3/)) )
      case('FeM')
         call setStandard(ioAtomParameters, 55.845d0, &
              (/1.900d0,1.900d0,5.350d0/),&
              (/-7.600d0,-3.800d0,-9.200d0/),&
              zeta2=1.800d0, c1=.537d0, c2=.668d0 )
      case('Ni')
         call setStandard(ioAtomParameters, 58.6934d0, &
              (/1.825d0,1.125d0,5.750d0/),&
              (/-9.17d0,-5.15d0,-13.49d0/),&
              zeta2=2.000d0, c1=.568d0, c2=.630d0,&
              VOIP_D = reshape(&
              (/ (/ 0.9110000d0, 8.5600000d0, 7.5380000d0, &
                    0.9110000d0, 9.5500000d0, 8.9630000d0, &
                    0.9110000d0, 9.3790000d0,10.6620000d0 /), &
                 (/ 0.9856000d0, 6.5520000d0, 3.8930000d0, &
                    0.9856000d0, 7.9030000d0, 5.1330000d0, &
                    0.9856000d0, 7.9030000d0, 5.0700000d0/), &
                 (/ 1.7600000d0,11.8400000d0, 5.9000000d0, &
                    1.7600000d0,13.7240000d0,10.0300000d0, &
                    1.7600000d0,13.4140000d0,11.8890000d0  /) /), (/3,3,3/)) )
      case('NiM')
         call setStandard(ioAtomParameters, 58.6934d0, &
              (/2.100d0,2.100d0,5.750d0/),&
              (/-7.80d0,-3.70d0,-9.90d0/),&
              zeta2=2.000d0, c1=.568d0, c2=.629d0 )
      case('Cu')
         call setStandard(ioAtomParameters, 63.546d0, &
              (/2.200d0,2.200d0,5.95d0/),&
              (/-11.4d0,-6.06d0,-14.00d0/),&
              zeta2=2.300d0, c1=.593d0, c2=.574d0,&
              VOIP_D = reshape(&
              (/ (/ 0.9420000d00,9.4590000d0,7.7240000d0, &
                    0.9420000d00,9.6830000d0,9.1740000d0, &
                    0.9420000d00,9.9550000d0,10.860000d0 /), &
                 (/ 0.4340000d0,7.2530000d0,3.987000d0, &
                    0.4340000d0,9.9930000d0,5.281000d0, &
                    0.4340000d0,9.9930000d0,5.033000d0/), &
                 (/ 1.773000d0,16.117000d0, 2.480000d0, &
                    1.773000d0,14.307000d0,10.662000d0, &
                    1.773000d0,13.699000d0,12.856000d0  /) /), (/3,3,3/)) )
      case('Ag')
         call setStandard(ioAtomParameters, 107.8682d0, &
              (/1.697d0,1.201d0,3.912d0/),&
              (/-7.57d0,-3.80d0,-12.70d0/),&
              zeta2=1.545d0, c1=.825d0, c2=.329d0,&
              VOIP_D = reshape(&
              (/ (/ 0.5500000d00,8.3900000d0,7.5800000d0, &
                    0.3700000d00,8.8800000d0,8.8000000d0, &
                    0.3100000d00,9.7100000d0,10.230000d0 /), &
                 (/ 0.7700000d0,6.4600000d0,3.830000d0, &
                    1.1800000d0,6.8600000d0,8.120000d0, &
                    1.1800000d0,6.8600000d0,4.760000d0/), &
                 (/ 3.900000d0,25.600000d0, 0.000000d0, &
                    0.460000d0,12.660000d0,12.770000d0, &
                    0.810000d0,11.670000d0,14.490000d0  /) /), (/3,3,3/)) )
      case('Ag1')
         call setStandard(ioAtomParameters, 107.8682d0, &
              (/1.842d0,1.492d0,3.921d0/),&
              (/-9.24d0,-7.12d0,-13.68d0/),&
              zeta2=1.565d0, c1=.823d0, c2=.328d0)
      case('Pt')
         call setStandard(ioAtomParameters, 195.084d0, &
              (/1.972d0,1.333d0,4.084d0/),&
              (/-9.08d0,-5.48d0,-12.59d0/),&
              zeta2=1.840d0, c1=.798d0, c2=.352d0)
      case('Au')
         call setStandard(ioAtomParameters, 196.966569d0, &
              (/2.602d0,2.584d0,6.163d0/),&
              (/-10.92d0,-5.55d0,-15.07d0/),&
              zeta2=2.794d0, c1=.644d0, c2=.536d0)
      case default
         write(*,'("Atomic symbol:",A4)') elem_name
         stop 'load_default: Thre is no default parameters for that atom.'
      end select
    end subroutine load_default

    subroutine write_information(iAtomParameters)
      type(THuckelAtomParameters), intent(in) :: iAtomParameters
      integer :: j,k,kh
      j=iAtomParameters%num_val_orb
      kh=iAtomParameters%highest_occupied_orbital
      write(*,'("SetAtomParameters:")')
      write(*,'(" Total Number of atomic orbitals =",I8)')  j
      write(*,'(" Orbital            =",10A6)')  (iAtomParameters%orbital(k),k=1, j)
      write(*,'(" Principal          =",10I6)')   iAtomParameters%principal(1:j)
      write(*,'(" Angular            =",10I6)')   iAtomParameters%angular(1:j)
      write(*,'(" zeta      =",4ES10.3)') AtomParameter(i)%zeta(1:j)
      write(*,'(" zeta (Rep)=",4ES10.3)') AtomParameter(i)%rescaled_zeta(1:j)
      if(j >= 5)then
         write(*,'("   <--- For d-orbitals")')
         write(*,'(" zeta      =",5ES10.3)') AtomParameter(i)%zeta(5:j)
         write(*,'(" zeta2     =",5ES10.3)') AtomParameter(i)%zeta2(5:j)
         write(*,'(" c1        =",5ES10.3)') AtomParameter(i)%c1(5:j)
         write(*,'(" c2        =",5ES10.3)') AtomParameter(i)%c2(5:j)
         write(*,'(" zeta (Rep)=",5ES10.3)') AtomParameter(i)%rescaled_zeta(5:j)
         write(*,'(" zeta2(Rep)=",5ES10.3)') AtomParameter(i)%rescaled_zeta2(5:j)
         write(*,'(" c1   (Rep)=",5ES10.3)') AtomParameter(i)%rescaled_c1(5:j)
         write(*,'(" c2   (Rep)=",5ES10.3)') AtomParameter(i)%rescaled_c2(5:j)
         write(*,'("        For d-orbitals --->")')
      end if
      write(*,'(" initial_occupation =",10F6.3)') iAtomParameters%initial_occupation(1:j)
      write(*,'(" Highest Occupied orbital = ",I1,A1)')&
           iAtomParameters%principal(kh), iAtomParameters%orbital(kh)
      write(*,'(" chemical_hardness = ",ES21.14)')&
           iAtomParameters%chemical_hardness 
      write(*,'(" diag. val. = ", ES21.14)') &
           (iAtomParameters%initial_diagonal_elements(k),k=1, j)
      if(j<=4) then
         write(*,'(" VOIP       = ",3ES21.14)') &
              (iAtomParameters%VOIP(:,k),k=1, j)
      else
         write(*,'(" VOIP_D     = ",3ES21.14)') &
              (iAtomParameters%VOIP_D(:,:,k),k=1, j)
      end if
    end subroutine write_information

!***********************************************************************
!     Heap Sort Program
!         written by S. Yamamoto, Dec 12, 1998
!         last revised Jun 23, 2003
! When there are some elements which have the same value,
! the original order isn't conserved.
! This program doesn't need an array for the temporary work area.
!     Inputs
!     N:   number of list elements
!     value: given list of values
!     Outputs (On exit, values are not altered.)
!     index: (value(index(i)),i=1,N) is a non-decreasing series.
!***********************************************************************
!**   Heap sort is not suitable for a recursive algorithm, becase    ***
!** a hierarchy of heap spreads over whole range of storage.         ***
!***********************************************************************
    subroutine dhsort(N,value,index)
      implicit none
      integer, parameter :: knd=kind(0d0)
      integer N,i
      integer index(N),p,q,r,s,t
      real(kind=knd)::  value(N),u
      !     Making Heap ( Ordered Binary Tree )
      !     (i)Every node has greater value than those of its descendents.
      !     (ii)When index(p) is parent node, its child nodes are 
      !     index(2*p) and index(2*p+1).
      index(1)=1
      do i=2, N
         !     A heap is stored in index(1:i-1).
         !     appending the i-th node to the heap
         p=i
         !     If the value of the parent node is less than that of the current
         !     node, exchange their positions with each other.
         do while(p .gt. 1)
            q=p/2
            r=index(q)
            !     If all the ancestors of node p is greater than node p,
            !     quit inner loop.
            if(value(r).ge.value(i)) exit
            !     If not, exchange the positions with each other.
            index(p)=r
            p=q
         end do
         index(p)=i
      end do
      do i=N, 2, -1
         !     Moving the last node to the root.
         r=index(i)
         u=value(r)
         !     Deleting the root node.
         !     (Now, the size of heap decreases by one.)
         index(i)=index(1)
         p=1
         q=2
         !     To keep the rule of heap, the root node (previously,
         !     the last node) will be stored in index(p).
         !     index(q) and index(q+1) are the children of index(p).
         do while(q+1 .lt. i)
            s=index(q)
            t=index(q+1)
            if(value(t).ge.value(s)) then
               !     Node q+1 is greater than node q.
               !     Node q should have the smallest value among two children,
               !     because we exchange the node q with parent node p later.
               q=q+1
               s=t
            end if
            !     If all the descendents of the node p is smaller than the node p,
            !     quit inner loop.
            if(u.ge.value(s)) exit
            !     If not, exchange the positions.
            index(p)=s
            p=q
            q=2*q
         end do
         !     for the case that the last node have only one child
         if(q+1.eq.i) then
            s=index(q)
            if(u.lt.value(s)) then
               index(p)=s
               p=q
            end if
         end if
         index(p)=r
      end do
    end subroutine dhsort
  end subroutine SetAtomParameters

end module M_qm_geno_Huckel_atom_params
