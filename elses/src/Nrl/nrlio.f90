!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Note from 200705
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!0527: A bugfix (information from Dr. Iguchi)
!      (old)   oresult = .false.
!      (fixed) oresult = NRL_IO_OVERLAP_NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module MNrlIO
  use MNrlMatrix
  use MNrlOnsite
  use MNrlHopping
  use flib_dom
  use elses_xml_misc
  use M_config

  implicit none

  private

  public NrlReadFromUnit

  ! module variable
  integer :: m_mode

  ! constants
  integer, parameter :: IO_NRL = 1

  interface NrlReadFromUnit
    module procedure ReadTBParametersFromUnit
  end interface

  interface NrlReadFromUnit
    module procedure ReadTBParametersFromXML
  end interface

  integer, parameter :: NRL_IO_OVERLAP_NONE = 0
  integer, parameter :: NRL_IO_OVERLAP      = 1
  integer, parameter :: NRL_IO_OVERLAP_EX   = 2
  
contains

  !!
  ! Read tight-binding parameter file header and set file reading mode.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHeaderLine(iunit)
    integer, intent(in) :: iunit
    integer :: the_stat
    character(len=32) :: the_header

    read(iunit, "(A32)", iostat=the_stat) the_header
    if(the_stat /= 0) goto 800

    select case(the_header)
    case('TYPE NRL', 'TYPE XB')
      m_mode = IO_NRL
    case default
      goto 810
    end select
    return

  ! On exceptional cases 
  800 continue
    print *, "Can not read a header line in tight-binding parameter file."
    stop
  810 continue
    print *, "Not supported header line in tight-binding parameter file."
    stop
  end subroutine


  !!
  ! Reads the number of atomic species. Usually, this is one.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadSpeciesCount(iunit, ocount)
    integer, intent(in) :: iunit
    integer, intent(out) :: ocount

    integer :: the_stat

    read(iunit, *, iostat=the_stat) ocount
    if(the_stat /= 0) goto 800

    read(iunit, *, iostat = the_stat)
    if(the_stat /= 0) goto 810
    return

  ! On exceptional case
  800 continue
    print *, "Can not read the number of atomic species"
    stop
  810 continue
    print *, "Unexpected end of file in tight-binding parameter file."
    stop
  end subroutine

  !!
  ! Checks the orbital information is assumed one.
  ! Now only supported mode is "spd 11"
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadBasicInfo(iunit, ol_count, ovalence, omass)
    integer, intent(in) :: iunit
    integer, intent(out) :: ol_count
    real(8), intent(out) :: ovalence
    real(8), intent(out) :: omass

    integer :: the_stat
    character(len=32) :: the_line, the_orbital_info
    real(8) :: the_valence
 
    read(iunit, "(A32)", iostat=the_stat) the_line
    if(the_stat /= 0) goto 810
    read(the_line, *, iostat=the_stat) the_orbital_info, ovalence, omass 
    if(the_stat /= 0) goto 820

    select case(the_orbital_info)
    case("spd")
      ol_count = 3
    case default
      goto 800
    end select

    ! skip one line
    read(iunit, *, iostat=the_stat)
    if(the_stat /= 0) goto 810
    return

  !On exceptional case
  800 continue
    write(0,*) "Not supported orbital information."
    write(0,*) "Please check the tight-binding paramter file."
    stop

  810 continue
    write(0,*) "Unexpected end of file in tight-binding parameter file."
    stop

  820 continue
    write(0,*) "Illegal format of tight-binding parameter file."
    stop
  end subroutine

  !!
  ! Reads and checks the cut-off length.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadGlobalCutOff(io_builder, iunit)
    type(TNrlMatrixBuilder), intent(inout) :: io_builder
    integer, intent(in) :: iunit
    real(8) :: the_cutoff
    integer :: the_err

    read(iunit, *, iostat=the_err) the_cutoff
    if(the_err /= 0) goto 810
    print *, "GLOBAL CUTOFF LENGTH", the_cutoff
    io_builder%cut_off = the_cutoff

    read(iunit, *, iostat=the_err)
    if(the_err /= 0) goto 810

    return
  ! On Exception
  810 continue
    print *, "Unexpected end of file in ReadGlobalCutOff"
    stop
  end subroutine
  
  !!
  ! Reads header part of tight-binding potential file.
  ! @param ovalence the number of valence electron in each atomic species
  ! @param omass    the mass of each atomic species   
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadHeader(iunit, ospecies_count, &
                        oorbital_count, ovalence, omass) 
    integer, intent(in) :: iunit
    integer, intent(out) :: ospecies_count
    integer, pointer :: oorbital_count(:)
    real(8), pointer :: ovalence(:)
    real(8), intent(out) :: omass(:)
    integer :: i, the_err
    real(8) :: the_valence, the_mass

    call ReadHeaderLine(iunit)
    call ReadSpeciesCount(iunit, ospecies_count)
    allocate(oorbital_count(ospecies_count), &
             ovalence(ospecies_count), &
             stat=the_err) 
    if(the_err /= 0) goto 900

    do i = 1, ospecies_count
      call ReadBasicInfo(iunit, oorbital_count(i), ovalence(i), omass(i))
    end do

    return
  900 continue
    stop "Memory allocation error in ReadHeader"
  end subroutine

  !!
  ! Read one line and returns true if it is overlap header.
  ! Otherwise, the line is push back.
  !
  function GetOverlapSection(iunit) result(oresult)
    integer :: oresult
    integer, intent(in) :: iunit

    integer :: the_stat
    character(len=80) :: the_header

    read(iunit, "(A80)", iostat=the_stat) the_header
    if(the_stat /= 0) then
!     oresult = .false.
      oresult = NRL_IO_OVERLAP_NONE
      return
    end if

    select case(the_header) 
    case("OVERLAP")
      oresult = NRL_IO_OVERLAP
      read(iunit, *, iostat=the_stat)
    case("OVERLAP EX")
      oresult = NRL_IO_OVERLAP_EX
      read(iunit, *, iostat=the_stat)
    case default
      oresult = NRL_IO_OVERLAP_NONE
      backspace(iunit)
    end select
  end function

  !!
  ! Reads on-site, hopping, intra-atom correction parameters.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadTBParametersFromUnit(iobuilder, iunit, oelec_count, awt)
    type(TNrlMatrixBuilder), intent(inout) :: iobuilder
    integer, intent(in) :: iunit
    integer :: the_species_count
    integer, pointer :: the_orbital_count(:)
    real(8) delta, lambda
    real(8), pointer :: oelec_count(:)
    real(8), intent(out) :: awt(:)
    
    nullify(the_orbital_count)
    call ReadHeader(iunit, the_species_count, &
                    the_orbital_count, oelec_count, awt)

    ! Read cut-off length for global construction.
    call ReadGlobalCutoff(iobuilder, iunit)

    ! Hamiltonian matrix -- on-site
    call NrlPrepareForRead(iobuilder%onsite, the_species_count, &
                            the_orbital_count)
    call NrlReadFromUnit(iobuilder%onsite, iunit)

    ! Hamiltonian matrix -- off-site
    call NrlPrepareForRead(iobuilder%hop, the_species_count, &
                           the_orbital_count)
    call NrlReadFromUnit(iobuilder%hop, iunit)
 
    select case(GetOverlapSection(iunit))
    case(NRL_IO_OVERLAP)
      ! Overlap matrix -- on-site
      call NrlPrepareForRead(iobuilder%overlap_onsite, the_species_count, &
                             the_orbital_count)
      call NrlSetIdentity(iobuilder%overlap_onsite)

      ! Overlap matrix -- off-site
      call NrlPrepareForRead(iobuilder%overlap, the_species_count, &
                             the_orbital_count)
      call NrlReadFromUnit(iobuilder%overlap, iunit)
    case(NRL_IO_OVERLAP_EX)
      ! Overlap matrix -- on-site
      call NrlPrepareForRead(iobuilder%overlap_onsite, the_species_count, &
                             the_orbital_count)
      call NrlReadFromUnit(iobuilder%overlap_onsite, iunit)

      ! Overlap matrix -- off-site
      call NrlPrepareForRead(iobuilder%overlap, the_species_count, &
                             the_orbital_count)
      call NrlReadFromUnit(iobuilder%overlap, iunit)
    end select

    if(associated(the_orbital_count)) deallocate(the_orbital_count)
  end subroutine

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ReadTBParametersFromXML(dnode, iobuilder, filename, oelec_count, awt)
    implicit none
    type(fnode), target      :: dnode
    type(fnode), pointer     :: document_node
    type(fnode), pointer     :: element_node
    type(fnode), pointer     :: classical_node
    type(fnode), pointer     :: mass_node
    type(fnode), pointer     :: ioniccharge_node
    type(fnode), pointer     :: pair_node
    type(fnode), pointer     :: rcut_node
    type(fnode), pointer     :: rscreen_node
    type(fnode), pointer     :: numorb_node
    type(fnode), pointer     :: occupancy_node
    type(fnode), pointer     :: lambda_node
    type(fnodeList), pointer     :: onsite_node
    type(fnodeList), pointer     :: hamiltonian_node
    type(fnodeList), pointer     :: overlap_node
    type(TNrlMatrixBuilder), intent(inout) :: iobuilder
    character(len=*), intent(in) :: filename
    real(8), intent(out) :: awt(:)

    integer :: i, j, ispecies, jspecies, nlength, ierr
    integer :: the_species_count
    integer, allocatable :: the_orbital_count(:)
    real(8), pointer :: oelec_count(:)
    real(8) lambda, cutoff, delta
    integer, allocatable, dimension(:) ::  ivalue
    real(8), allocatable, dimension(:) ::  rvalue
    character(len=256)       :: value
    logical        :: ex

    document_node => dnode
    the_species_count = config%system%structure%nelement

    call initialize
    call NrlAllocateOnsiteXML(iobuilder%onsite, the_orbital_count, the_species_count)
    call NrlAllocateOnsiteXML(iobuilder%overlap_onsite, the_orbital_count, the_species_count)
    call NrlAllocateHoppingMatrixXML(iobuilder%hop, the_species_count, the_orbital_count)
    call NrlAllocateHoppingMatrixXML(iobuilder%overlap, the_species_count, the_orbital_count)

    ! get <element> node
    element_node => getFirstElementByTagName(document_node,"element")

    if( .not. associated(element_node) ) then
       call XML_error("<element> not found")
    endif

    ! get name attribute
    value = getAttribute(element_node,"name")
    if( value == "" ) then
       call XML_error("<element name> not found")
    endif
    do i = 1, config%system%structure%nelement
       if( config%system%structure%velement(i)%name == value ) then
          ispecies = i
       end if
    end do

    ! get <classical_machanics> node
    classical_node => getFirstElementByTagName(element_node,"classical_mechanics")
    if( .not. associated(classical_node) ) then
       call XML_error("<classical> not found")
    endif

    ! get <mass> node
    mass_node => getFirstElementByTagName(classical_node,"mass")
    if( .not. associated(mass_node) ) then
       call XML_error("<mass> not found")
    else
       allocate( rvalue(1), stat=ierr )
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation rvalue<mass_node> in sbrt Read TBParametersFromXML'
       end if
       call real_load( mass_node, rvalue, 1 )
       config%system%structure%velement(ispecies)%classic%mass = rvalue(1)
       awt(ispecies) = rvalue(1)
       write(*,*) awt(ispecies)
       deallocate( rvalue )
    endif

    ioniccharge_node => getFirstElementByTagName(classical_node,"ionic_charge")
    if( .not. associated(ioniccharge_node) ) then
       ! set ionit_charge = 0.d0
    else
       allocate( rvalue(1), stat=ierr)
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation rvalue<ioniccharge_node> in sbrt Read TBParametersFromXML'
       end if
       call real_load( ioniccharge_node, rvalue, 1 )
       ! set rvalue(1) => ioniccharge
       deallocate( rvalue )
    endif

    ! get <pair> node
    pair_node => getFirstElementByTagName(element_node,"tight_binding_pair")
    if( .not. associated(pair_node) ) then
       call XML_error("<pair> not found")
    endif

    ! get name attribute
    value = getAttribute(pair_node,"type")
    if( value == "" ) then
       call XML_error("<type> not found")
    endif

    select case(value)
    case("self","SELF","Self")
       jspecies = ispecies
    case default
       do i = 1, config%system%structure%nelement
          if( config%system%structure%velement(i)%name == value ) then
             jspecies = i
          end if
       end do
    end select

    ! get <rcut> node
    rcut_node => getFirstElementByTagName(pair_node,"rcut")
    if( .not. associated(rcut_node) ) then
       call XML_error("<rcut> not found")
    else
       allocate( rvalue(1), stat=ierr )
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation rvalue<rcut_node> in sbrt Read TBParametersFromXML'
       end if
       call real_load( rcut_node, rvalue, 1 )
       iobuilder%cut_off = rvalue(1)
       cutoff = rvalue(1)
       deallocate( rvalue )
    endif

    ! get <rscreen> node
    rscreen_node => getFirstElementByTagName(pair_node,"rscreen")
    if( .not. associated(rscreen_node) ) then
       call XML_error("<rscreen> not found")
    else
       allocate( rvalue(1), stat=ierr )
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation rvalue<rscreen_node> in sbrt Read TBParametersFromXML'
       end if
       call real_load( rscreen_node, rvalue, 1 )
       delta = rvalue(1)
       deallocate( rvalue )
    endif

    ! get <numorb> node
    numorb_node => getFirstElementByTagName(pair_node,"num_orb")
    if( .not. associated(numorb_node) ) then
       call XML_error("<numorb> not found")
    else
       allocate( ivalue(1), stat=ierr )
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation ivalue<numorb_node> in sbrt Read TBParametersFromXML'
       end if
       call integer_load( numorb_node, ivalue, 1 )
       config%system%structure%velement(ispecies)%quantum%orbital = ivalue(1)
       deallocate( ivalue )
    endif

    ! get <occupancy> node
    occupancy_node => getFirstElementByTagName(pair_node,"occupancy")
    if( .not. associated(occupancy_node) ) then
       call XML_error("<occupancy> not found")
    else
       allocate( rvalue(3), stat=ierr )
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation rvalue<occupancy_node> in sbrt Read TBParametersFromXML'
       end if
       call real_load( occupancy_node, rvalue, 3 )
       oelec_count(ispecies) = sum(rvalue(1:3))
       config%system%structure%velement(ispecies)%classic%charge = oelec_count(ispecies)
       deallocate( rvalue )
    endif

    ! get <lambda> node
    lambda_node => getFirstElementByTagName(pair_node,"lambda")
    if( .not. associated(lambda_node) ) then
       call XML_error("<lambda> not found")
    else
       allocate( rvalue(1), stat=ierr )
       if( ierr /= 0 ) then
          write(*,*) 'Error : allocation rvalue<lambda_node> in sbrt Read TBParametersFromXML'
       end if
       call real_load( lambda_node, rvalue, 1 )
       lambda = rvalue(1)
       deallocate( rvalue )
    end if

    call NrlOnsiteSetXML(iobuilder%onsite,lambda,cutoff,delta,ispecies,jspecies)

    ! get <onsite> node
    onsite_node => getElementsByTagName(pair_node,"onsite")
    if( .not. associated(onsite_node) ) then
       call XML_error("<onsite> not found")
    else
       nlength = getLength(onsite_node)
       do i = 1, nlength
          call NrlOnsiteLoadXML(iobuilder%onsite,item(onsite_node,i-1),the_orbital_count,the_species_count,ispecies)
       end do
       call NrlOnsiteSetIdentityXML(iobuilder%overlap_onsite, the_orbital_count, the_species_count, ispecies)
    endif
    
    ! get <hamiltonian> node
    hamiltonian_node => getElementsByTagName(pair_node,"hamiltonian")
    if( .not. associated(hamiltonian_node) ) then
       call XML_error("<hamiltonian> not found")
    else
       nlength = getLength(hamiltonian_node)
       do i = 1, nlength
          call NrlHoppingMatrixLoadXML(iobuilder%hop,item(hamiltonian_node,i-1),delta,cutoff,ispecies,jspecies)
       end do
    endif

    ! get <overlap> node
    overlap_node => getElementsByTagName(pair_node,"overlap")
    nlength = getLength(overlap_node)
    if( .not. associated(overlap_node) ) then
       call XML_error("<overlap> not found")
    else
       nlength = getLength(overlap_node)
       do i = 1, nlength
          call NrlHoppingMatrixLoadXML(iobuilder%overlap,item(overlap_node,i-1),delta,cutoff,ispecies,jspecies)
       end do
    endif

  contains

    subroutine initialize
      implicit none

      the_species_count = config%system%structure%nelement
      
      allocate( the_orbital_count(the_species_count), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'Error : allocation of oelec_count in sbrt ReadTBParametersFromXML'
         stop
      end if

      the_orbital_count(:) = 3  ! s,p,d (defalut) !!
      allocate( oelec_count(the_species_count), stat=ierr )
      if( ierr /= 0 ) then
         write(*,*) 'Error : allocation of oelec_count in sbrt initialize'
         stop
      end if
      
      return
    end subroutine Initialize

    subroutine integer_load(vnode, ivalue, size)
      implicit none
      type(fnode), target, intent(in) :: vnode
      type(fnode), pointer :: node
      integer, intent(in) :: size 
      integer, intent(out) :: ivalue(size)
      character(len=256) value, unit

      node => vnode

      unit  = getAttribute(node,"unit")
      if( unit == "" ) then
         unit = "a.u."
      endif

      value = getChildValue(node)
      read(unit=value,fmt=*) ivalue(1:size)
      return
    end subroutine integer_load

    subroutine real_load(vnode,rvalue,size)
      implicit none
      type(fnode), target, intent(in) :: vnode
      type(fnode), pointer :: node
      integer, intent(in) :: size 
      real(8), intent(out) :: rvalue(size)
      character(len=256) value, unit

      node => vnode

      unit  = getAttribute(node,"unit")
      if( unit == "" ) then
         unit = "a.u."
      endif

      value = getChildValue(node)
      read(unit=value,fmt=*) rvalue(1:size)
      rvalue(1:size) = rvalue(1:size) * XML_TO_AU(unit)

      return
    end subroutine real_load

  end subroutine ReadTBParametersFromXML
end module MNrlIO
