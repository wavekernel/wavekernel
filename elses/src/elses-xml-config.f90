!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module elses_xml_element
  use M_config,  only : config
  use M_element
  use flib_dom
  use elses_xml_misc
  implicit none

  type :: DB_type
     integer :: size = 0
     type(element_type),pointer,dimension(:) :: data
  end type DB_type
  type(DB_type), save :: DB

  private

  public :: getElements, element_load

contains
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine getElements( size, data )
    integer :: size
    type(element_type),pointer,dimension(:) :: data

!    allocate( data(DB%size) )
    size = DB%size
    data => DB%data

    return 
  end subroutine getElements

  function DB_find( filename ) result(datum)
    character(len=*), intent(in)    :: filename
    type(element_type), pointer :: datum

    integer :: i

    datum => null()
    do i=1, DB%size
       if( DB%data(i)%filename == filename ) then
          datum => DB%data(i)
          exit
       end if
    end do

    return
  end function DB_find

  function DB_pushback() result(datum)
    type(element_type), pointer :: datum
    type(element_type),pointer,dimension(:) :: TMP
    integer :: alloc_size
    integer, parameter :: init_alloc_size=8

    if( .not. associated(DB%data) ) then
       allocate( DB%data(init_alloc_size) )
       DB%size = 0
    end if

    alloc_size = size(DB%data,1)
    if( DB%size == alloc_size ) then
       allocate( TMP(alloc_size) )
       TMP(1:alloc_size) = DB%data(1:alloc_size)
       deallocate(DB%data)
       allocate(DB%data(alloc_size*2))
       DB%data(1:alloc_size)=TMP(1:alloc_size)
       deallocate(TMP)
       write(*,'("DB_pushback: Too small size(DB%data,1). Data sized is enlarged by a factor of 2. ")')
    end if

    DB%size = DB%size+1

    datum => DB%data(DB%size)

    return
  end function DB_pushback


  ! parse <element> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine element_load( element, filename, name, type )
    implicit none
    type(element_type),pointer   :: element
    character(len=*), intent(in) :: filename

!   type(fnode), pointer     :: document_node
!   type(fnode), pointer     :: element_node
!   type(fnode), pointer     :: classic_node
!   type(fnode), pointer     :: quantum_node
!   character(len=256)       :: value
    character(len=*), optional :: name
    character(len=*), optional:: type
!   logical        :: ex

    element => DB_find( filename )

    if( associated(element) ) then
       return
    end if

    element => DB_pushback()

    element%quantum%type = ''
    if( present( name ) )  element%name = name
    if( present( type ) ) element%quantum%type = type
    element%filename = filename

    select case( element%quantum%type )
    case('NRL','geno', '')

    case('GaAs.Molteni.1993', 'Xu.1992', 'Kwon.1994')
       select case(filename)

!---------- BEGIN added by T. Hoshi, 30. Jul 2007 -----
!       case("Au.xml")
!          element%filename = filename
!          element%name = "Au"
!          element%quantum%type = "NRL"
!---------- END: added by T. Hoshi, 30. Jul 2007 -----

!---------- BEGIN added by T. Hoshi, 17. Jul 2007 -----
       case("Ga.xml")
          element%name = "Ga"
          element%classic%mass = 69.723 * XML_TO_AU("amu")
          element%classic%charge = 3.00
          element%quantum%type = "GaAs.Molteni.1993"
          element%quantum%orbital = 5
          element%quantum%interaction_radius =  3.5 * XML_TO_AU("angstrom")

       case("As.xml")
          element%name = "As"
          element%classic%mass = 74.92159 * XML_TO_AU("amu")
          element%classic%charge = 5.00
          element%quantum%type = "GaAs.Molteni.1993"
          element%quantum%orbital = 5
          element%quantum%interaction_radius =  3.5 * XML_TO_AU("angstrom")

!---------- END: added by T. Hoshi, 17. Jul 2007 -----

       case("C.xml")
          element%name = "C"
          
          element%classic%mass = 12.01 * XML_TO_AU("amu")
          element%classic%charge = 4.00
          element%quantum%type = "Xu.1992"
          element%quantum%orbital = 4
          element%quantum%interaction_radius =  2.60 * XML_TO_AU("angstrom")
          element%quantum%DNAL0 = 2.00
          element%quantum%RNN0  = 1.53632900
          
          element%quantum%DHAL(1) = -5.000
          element%quantum%DHAL(2) =  4.700
          element%quantum%DHAL(3) =  5.500
          element%quantum%DHAL(4) = -1.550
          
          element%quantum%DNAL(1) =  6.500
          element%quantum%DNAL(2) =  6.500
          element%quantum%DNAL(3) =  6.500
          element%quantum%DNAL(4) =  6.500
          
          element%quantum%RCAL(1) =  2.180
          element%quantum%RCAL(2) =  2.180
          element%quantum%RCAL(3) =  2.180
          element%quantum%RCAL(4) =  2.180
          
          element%quantum%restpart = 0.00
          
       case("Si.xml")
          element%name = "Si"
          
          element%classic%mass = 28.0855 * XML_TO_AU("amu")
          element%classic%charge = 4.00
          element%quantum%type = "Kwon.1994"
          element%quantum%orbital = 4
          element%quantum%interaction_radius =  4.16 * XML_TO_AU("angstrom")
          element%quantum%DNAL0 = 2.00
          element%quantum%RNN0  = 2.360352
          
          element%quantum%DHAL(1) = -2.038
          element%quantum%DHAL(2) =  1.745
          element%quantum%DHAL(3) =  2.750
          element%quantum%DHAL(4) = -1.075
          
          element%quantum%DNAL(1) =  9.500
          element%quantum%DNAL(2) =  8.500
          element%quantum%DNAL(3) =  7.500
          element%quantum%DNAL(4) =  7.500
          
          element%quantum%RCAL(1) =  3.400
          element%quantum%RCAL(2) =  3.550
          element%quantum%RCAL(3) =  3.700
          element%quantum%RCAL(4) =  3.700
          
          element%quantum%restpart = 0.00
          
       case default
          write(*,*) "Sorry! not supported element ", filename
          stop
       end select

    case default

       write(*,*) "Sorry! not supported type", element%quantum%type
       stop

!!$ We disabled this function until we decide to read data from files
!!$ 
!!$    inquire( file=filename, exist=ex )
!!$    if( .not. ex ) then
!!$       write(*,'(a,a)') '# Error!: element_load : can not open file ', trim(filename)
!!$       stop
!!$    end if
!!$
!!$    document_node => parsefile(filename)
!!$    call normalize(document_node)
!!$
!!$    element%filename = filename
!!$
!!$    ! get <element> node
!!$    element_node => getFirstElementByTagName(document_node,"element")
!!$    if( .not. associated(element_node) ) then
!!$       call XML_error("<element> not found")
!!$    endif
!!$    
!!$    ! get name attribute
!!$    value = getAttribute(element_node,"name")
!!$    if( value == "" ) then
!!$       call XML_error("<element name> not found")
!!$    endif
!!$    read(unit=value,fmt=*) element%name
!!$
!!$    ! get <classic> node
!!$    classic_node => getFirstElementByTagName(element_node,"classic")
!!$    if( .not. associated(classic_node) ) then
!!$       call XML_error("<classic> not found")
!!$    endif
!!$
!!$    call classic_load( element%classic, classic_node )
!!$
!!$    ! get <quantum> node
!!$    quantum_node => getFirstElementByTagName(element_node,"quantum")
!!$    if( .not. associated(quantum_node) ) then
!!$       call XML_error("<quantum> not found")
!!$    endif
!!$
!!$    call quantum_load( element%quantum, quantum_node )
!!$
!!$    call destroyNode(document_node)

    end select
    return
  end subroutine element_load

  ! parse <classic> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine classic_load( classic, classic_node )
    implicit none
    type(classic_type), intent(out) :: classic
    type(fnode), pointer     :: classic_node

    type(fnode), pointer     :: node
    character(len=256)       :: value, unit

    ! get <mass> node
    node => getFirstElementByTagName(classic_node,"mass")
    if( .not. associated(node) ) then
       call XML_error("<mass> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) classic%mass

    unit  = getAttribute(node,"unit")
    if( unit == "" ) then
       unit = "amu"
    endif
    classic%mass = classic%mass * XML_TO_AU(unit)


    ! get <charge> node
    node => getFirstElementByTagName(classic_node,"charge")
    if( .not. associated(node) ) then
       call XML_error("<charge> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) classic%charge

    return
  end subroutine classic_load

  ! parse <quantum> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine quantum_load( quantum, quantum_node )
    implicit none
    type(quantum_type), intent(out) :: quantum
    type(fnode), pointer     :: quantum_node

    type(fnode), pointer     :: node
    character(len=256)       :: value, unit

    ! get name attribute
    value = getAttribute(quantum_node,"type")
    if( value == "" ) then
       call XML_error("<quantum type> not found")
    endif
    read(unit=value,fmt=*) quantum%type


    ! get <orbital> node
    node => getFirstElementByTagName(quantum_node,"orbital")
    if( .not. associated(node) ) then
       call XML_error("<orbital> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%orbital

    ! get <interaction_radius> node
    node => getFirstElementByTagName(quantum_node,"interaction_radius")
    if( .not. associated(node) ) then
       call XML_error("<interaction_radius> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%interaction_radius

    unit  = getAttribute(node,"unit")
    if( unit == "" ) then
       unit = "a.u."
    endif
    quantum%interaction_radius = quantum%interaction_radius * XML_TO_AU(unit)


    ! get <DNAL0> node
    node => getFirstElementByTagName(quantum_node,"DNAL0")
    if( .not. associated(node) ) then
       call XML_error("<DNAL0> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%DNAL0

    ! get <RNN0> node
    node => getFirstElementByTagName(quantum_node,"RNN0")
    if( .not. associated(node) ) then
       call XML_error("<RNN0> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%RNN0

    ! get <DHAL> node
    node => getFirstElementByTagName(quantum_node,"DHAL")
    if( .not. associated(node) ) then
       call XML_error("<DHAL> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%DHAL(1:4)

    ! get <DNAL> node
    node => getFirstElementByTagName(quantum_node,"DNAL")
    if( .not. associated(node) ) then
       call XML_error("<DNAL> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%DNAL(1:4)

    ! get <RCAL> node
    node => getFirstElementByTagName(quantum_node,"RCAL")
    if( .not. associated(node) ) then
       call XML_error("<RCAL> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%RCAL(1:4)

    ! get <restpart> node
    node => getFirstElementByTagName(quantum_node,"restpart")
    if( .not. associated(node) ) then
       call XML_error("<restpart> not found")
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) quantum%restpart

    return
  end subroutine quantum_load


end module elses_xml_element

module elses_xml_structure
  use M_config,  only : config
  use M_structure
  use flib_dom
  use elses_xml_misc
  use elses_xml_element
  implicit none

  private
  public :: structure_loadXYZ, structure_load, structure_element_load, &
            structure_saveXYZ, structure_save, &
            unitcell_load, heatbath_load, &
            unitcell_save, heatbath_save

contains

  ! load structure data from XYZ
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine structure_loadXYZ( structure, filename, cell_info )
  !! NOTE(T.Hoshi, 2010May) 
  !!    : This routine is called ONLY from tool/src/elses-xml-generate.f90
    implicit none
    type(structure_type), intent(out) :: structure
    character(len=*),      intent(in) :: filename
    real(8),                optional  :: cell_info(3) ! cell info

    integer, parameter :: fd = 1
    logical        :: ex

    integer :: j
    character(len=256) :: name
    type(atom_type), pointer :: atom
    integer :: ierr
    character(len=1024) :: chara_wrk
    real(8)             :: cell_info_wrk(3)
!
    cell_info_wrk(3) = -1.0d0 ! dummy value
!
    inquire( file=filename, exist=ex )
    if( .not. ex ) then
       write(*,'(a,a)') '# Error!: structure_load : can not open file ', trim(filename) 
       stop
    end if
    open(fd,file=filename)

    read(fd,'(a)',end=200,iostat=ierr) chara_wrk
    if (len_trim(chara_wrk) > 1024) then
      write(*,*)'ABORT: File read error; (structure_loadXYZ)'
      write(*,*)'first line=', trim(chara_wrk)
      stop
    endif   
!
!   write(*,*)'first line =',trim(chara_wrk)
!   write(*,*)'len_trim    =',len_trim(chara_wrk)
!
    read(chara_wrk,*,iostat=ierr) structure%natom, cell_info_wrk(1:3)
    if (ierr == 0) then
       if (present(cell_info)) then
!        write(*,*) 'INFO:cell info appears in the XYZ file and will be set in the XML file'
!        write(*,*) 'INFO:cell info [A]  =', cell_info_wrk(1:3) 
         cell_info(1:3)=cell_info_wrk(1:3)*XML_TO_AU("angstrom") ! cell info in au
!        write(*,*) 'INFO:cell info [au] =', cell_info(1:3) 
       else
         write(*,*) 'INFO:cell info appears in the XYZ file but is ignored.'
       endif   
    else   
!     write(*,*) 'INFO:no cell info appears in XYZ file'
      read(chara_wrk,*,iostat=ierr) structure%natom
      if (ierr /= 0) then
        write(*,'(a,a)') 'File name : ', trim(filename) 
        write(*,*)'ABORT: File read error; (structure_loadXYZ)'
        write(*,*)'  for structure%natom'
        stop
      endif   
    endif   
!
    allocate( structure%vatom(structure%natom) )

    read(fd,*,iostat=ierr) structure%name
    if (ierr /= 0) then
      write(*,'(a,a)') 'File name : ', trim(filename) 
      write(*,*)'ABORT: File read error; (structure_loadXYZ)'
      write(*,*)'  for structure%name'
      stop
    endif   

    structure%mdstep = 0
    structure%unitcell%vectorA(1) = XML_TO_AU("angstrom") ! dummy value 
    structure%unitcell%vectorA(2) = 0.d0                  ! dummy value
    structure%unitcell%vectorA(3) = 0.d0                  ! dummy value
    structure%unitcell%vectorB(1) = 0.d0                  ! dummy value
    structure%unitcell%vectorB(2) = XML_TO_AU("angstrom") ! dummy value
    structure%unitcell%vectorB(3) = 0.d0                  ! dummy value
    structure%unitcell%vectorC(1) = 0.d0                  ! dummy value 
    structure%unitcell%vectorC(2) = 0.d0                  ! dummy value
    structure%unitcell%vectorC(3) = XML_TO_AU("angstrom") ! dummy value
    structure%unitcell%set = .true.

    do j=1, structure%natom
       atom => structure%vatom(j)
       read(fd,*,end=200) name, atom%position(1:3)

       call element_load( atom%element, trim(name) // ".xml", name )

       atom%position = &
            + atom%position(1) * structure%unitcell%vectorA &
            + atom%position(2) * structure%unitcell%vectorB &
            + atom%position(3) * structure%unitcell%vectorC 

       atom%class = ""
       atom%motion = "free"
       atom%velocity_set = .false.
       atom%velocity = 0.d0

       atom%force_set = .false.
       atom%force = 0.d0

       atom%ncustumize = 0
    end do

!100 continue
   
    close(fd)
    
    return
200 continue ! error_block

    call XML_error("something wrong in the XYZ-formated file")
    
    close(fd)
    stop

  end subroutine structure_loadXYZ


  ! save structure data in XYZ
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine structure_saveXYZ( structure, filename, append )
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: append

    integer, parameter :: fd = 71
    logical, save :: first = .true.

    integer :: j

    close(fd)
    if( (.not. present(append)) .or. (.not. append) ) then
       open(fd,file=filename)
    elseif( first ) then
       open(fd,file=filename)
       first = .false.
    else
       open(fd,file=filename,position='append')
    end if

    write(fd,*) structure%natom
    write(fd,*) trim(structure%name), " mdstep=", structure%mdstep

    do j=1, structure%natom
       write(fd,'(a4,3e23.15)') &
            structure%vatom(j)%name, &
            structure%vatom(j)%position(1:3)*XML_AU_TO("angstrom")
    end do

    close(fd)

    return
  end subroutine structure_saveXYZ

  ! parse <structure> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine structure_load( structure, filename )
    implicit none
    type(structure_type), intent(out) :: structure
    character(len=*), intent(in) :: filename

    type(fnode), pointer     :: document_node
!   type(fnode), pointer     :: base_node
    type(fnode), pointer     :: structure_node
    type(fnode), pointer     :: unitcell_node
    type(fnode), pointer     :: heatbath_node
    type(fnodeList), pointer :: vatom_node

    integer                  :: j
    character(len=256)       :: value
    logical        :: ex

    integer              :: i_verbose, log_unit

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    inquire( file=filename, exist=ex )
    if( .not. ex ) then
       write(*,'(a,a)') '# Error!: structure_load : can not open file ', trim(filename) 
       stop
    end if

    document_node => parsefile(filename)
    call normalize(document_node)

    ! get <structure> node
    structure_node => getFirstElementByTagName(document_node,"structure")

    if( .not. associated(structure_node) ) then
       call XML_error("<structure> not found")
    endif

    ! get name attribute
    value = getAttribute(structure_node,"name")
    if( value == "" ) then
       call XML_error("<structure name> not found")
    endif
    structure%name = value
    if (log_unit > 0) write(log_unit,'(a,a)') 'INFO-XML:structure name = ', trim(structure%name)
!
    ! get mdstep attribute
    value = getAttribute(structure_node,"mdstep")
    if( value == "" ) then
       structure%mdstep = 0
    else
       read(unit=value,fmt=*) structure%mdstep
    end if

    ! get <unitcell> node
    unitcell_node => getFirstElementByTagName(structure_node,"unitcell")
    if( .not. associated(unitcell_node) ) then
       structure%unitcell%set = .false.
    else
       call unitcell_load( structure%unitcell, unitcell_node )
       structure%unitcell%set = .true.
    endif

    ! get <heatbath> node
    heatbath_node => getFirstElementByTagName(structure_node,"heatbath")
    if( .not. associated(heatbath_node) ) then
       call heatbath_default( structure%heatbath )
    else
       call heatbath_default( structure%heatbath )
       call heatbath_load( structure%heatbath, heatbath_node )
    endif

    ! get <atom> nodes
    vatom_node => getElementsByTagName(structure_node,"atom")
    structure%natom = getLength(vatom_node)
    if( structure%natom == 0 ) then
      if (log_unit > 0) then 
        write(log_unit,'(a)') 'INFO-XML-WARN:<atom> tag is not found in the structure XML file'
        write(log_unit,'(a)') 'INFO-XML-WARN:structure%vatom(structure%natom) is not allcoated'
      endif  
    else
      if (log_unit > 0) write(log_unit,'(a)')'INFO-XML:allocate structure%vatom(structure%natom) at structure_load'
      allocate( structure%vatom(structure%natom) )
      do j=1, structure%natom
         call atom_load( structure%vatom(j), item(vatom_node,j-1), structure%unitcell )
      end do
    endif
!
    call destroyNode(document_node)

    return
  end subroutine structure_load

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine structure_element_load( structure, element_node )
    implicit none
    type(structure_type), intent(out) :: structure
    type(fnode), pointer     :: element_node

    integer :: j
    character(len=8)       :: name
    character(len=64)      :: model, filename
    integer              :: i_verbose, log_unit

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    name = getAttribute(element_node,"name")
    filename = getAttribute(element_node,"filename")
    model = getAttribute(element_node,"model")

    do j=1, structure%natom
       if( name == structure%vatom(j)%name .or. name == '') then
          if( filename == '' ) then
             ! Following comments are added by SY Jul 16, 2009.
            if (log_unit > 0) then 
              write(log_unit,'(a)') '"filename" is required in "element_load",&
                 & even if it is not used actually.'
              write(log_unit,'(a)') 'Therefore, now we assume that "atom name"&
                 & gives "filename".'
            endif  
            call element_load( structure%vatom(j)%element, &
                  trim(structure%vatom(j)%name) // ".xml" , name, model )
             ! name, model added by SY Nov28, 2008
          else
             call element_load( structure%vatom(j)%element, &
                  filename, name, model )
          end if
       end if
    end do

    call getElements( structure%nelement, structure%velement )

    return
  end subroutine structure_element_load


  ! parse <heatbath> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine heatbath_load( heatbath, heatbath_node )
    implicit none
    type(heatbath_type), intent(out) :: heatbath
    type(fnode), pointer     :: heatbath_node
    type(fnode), pointer     :: node
    character(len=256)       :: value, unit

    integer              :: i_verbose, log_unit

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    ! get <massperatom> node
    node => getFirstElementByTagName(heatbath_node,"massperatom")
    if( .not. associated(node) ) then
       call XML_error("<massperatom> not found")
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       endif
       value = getChildValue(node)
       read(unit=value,fmt=*) heatbath%massperatom
       heatbath%massperatom = heatbath%massperatom * XML_TO_AU(unit)
    endif

    ! get <position> node
    node => getFirstElementByTagName(heatbath_node,"position")
    if( .not. associated(node) ) then
       heatbath%position = 0.d0
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       endif
       value = getChildValue(node)
       read(unit=value,fmt=*) heatbath%position
       heatbath%position = heatbath%position * XML_TO_AU(unit)
    endif

    ! get <velocity> node
    node => getFirstElementByTagName(heatbath_node,"velocity")
    if( .not. associated(node) ) then
       heatbath%velocity = 0.d0
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       endif
       value = getChildValue(node)
       read(unit=value,fmt=*) heatbath%velocity
       heatbath%velocity = heatbath%velocity * XML_TO_AU(unit)
    endif

    return
  end subroutine heatbath_load

  ! parse <unitcell> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine unitcell_load( unitcell, unitcell_node )
    implicit none
    type(unitcell_type), intent(out) :: unitcell
    type(fnode), pointer     :: unitcell_node

    type(fnodeList), pointer :: vvector_node
    type(fnode), pointer     :: node
    character(len=256)       :: value, unit
    integer                  :: ierr
    integer              :: i_verbose, log_unit

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    ! get length scale
    node => &
    getFirstElementByTagName(unitcell_node,"myLength")
    if( associated(node) ) then
       value = getChildValue(node)
       read(value,*,iostat=ierr) unitcell%myLength
       if (ierr /= 0) then
         write(*,*)'ERROR in the input XML file:<unitcell><myLength>'
         stop
       endif    
       unit = getAttribute(node, "unit")
       if ( unit /= "" ) then
          unitcell%myLength = unitcell%myLength &
               * XML_TO_AU(unit)
          if (log_unit > 0)  write(log_unit,'(a,f10.5)') & 
&                   'INFO-XML:Optional tag detected : mylength [au] =',unitcell%myLength
       end if
    else
       unitcell%myLength = 1d0
    endif

    ! get <vector> nodes
    vvector_node => getElementsByTagName(unitcell_node,"vector")
    if( getLength(vvector_node) /= 3 ) then
       call XML_error("not sufficient <vector>s")
    endif

    ! A-axis
    node => item(vvector_node,0)
    unit  = getAttribute(node,"unit")
    value = getChildValue(node)
    read(unit=value,fmt=*,iostat=ierr) unitcell%vectorA(1:3)
    if (ierr /= 0) then
      write(*,*)'ERROR in the input XML file:<unitcell><vector>'
      stop
    endif    
    select case(unit)
    case("")
       unit = "a.u."
    case("myLength")
       unitcell%vectorA = unitcell%vectorA * unitcell%myLength
       unit = "a.u."
    end select
    unitcell%vectorA = unitcell%vectorA * XML_TO_AU(unit)

    ! B-axis
    node => item(vvector_node,1)
    unit  = getAttribute(node,"unit")
    value = getChildValue(node)
    read(unit=value,fmt=*,iostat=ierr) unitcell%vectorB(1:3)
    if (ierr /= 0) then
      write(*,*)'ERROR in the input XML file:<unitcell><vector>'
      stop
    endif    
    select case(unit)
    case("")
       unit = "a.u."
    case("myLength")
       unitcell%vectorB = unitcell%vectorB * unitcell%myLength
       unit = "a.u."
    end select
    unitcell%vectorB = unitcell%vectorB * XML_TO_AU(unit)

    ! C-axis
    node => item(vvector_node,2)
    unit  = getAttribute(node,"unit")
    value = getChildValue(node)
    read(unit=value,fmt=*,iostat=ierr) unitcell%vectorC(1:3)
    if (ierr /= 0) then
      write(*,*)'ERROR in the input XML file:<unitcell><vector>'
      stop
    endif    
    select case(unit)
    case("")
       unit = "a.u."
    case("myLength")
       unitcell%vectorC = unitcell%vectorC * unitcell%myLength
       unit = "a.u."
    end select
    unitcell%vectorC = unitcell%vectorC * XML_TO_AU(unit)

    return
  end subroutine unitcell_load

  ! parse <atom> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine atom_load( atom, atom_node, unitcell )
    implicit none
    type(atom_type), intent(out) :: atom
    type(fnode), pointer     :: atom_node
    type(unitcell_type), intent(in) :: unitcell

    type(fnodeList), pointer :: vcustumize_node
    type(fnode), pointer     :: node
    character(len=256)       :: value, unit
    integer :: i
    real(8) :: la, lb, lc ! Lengths of unitcell vectors in a.u.
!
    ! get element attribute
    value = getAttribute(atom_node,"element")
    if( value == "" ) then
       call XML_error("<atom element> not found")
    else
       atom%name = value
       atom%element => null()
    endif

    ! get class attribute
    value = getAttribute(atom_node,"class")
    if( value == "" ) then
       atom%class = ""
    else
       atom%class = value
    end if

    ! get md attribute
    value = getAttribute(atom_node,"motion")
    if( value == "" ) then
       atom%motion = "free"
    else
       atom%motion = value
    end if

    ! get group_id attribute
    value = getAttribute(atom_node,"group_id")
    if( value == "" ) then
       atom%group_id_set = .false.
       atom%group_id     = -1
    else
       atom%group_id_set = .true.
       read(unit=value,fmt=*) atom%group_id
    end if

    ! get <position> node
    node => getFirstElementByTagName(atom_node,"position")
    if( .not. associated(node) ) then
       call XML_error("<position> not found")
    endif
    unit  = getAttribute(node,"unit")
    if( unit == "" ) then
       unit = "a.u."
    endif
    value = getChildValue(node)
    read(unit=value,fmt=*) atom%position(1:3)

    la = dsqrt(dot_product(unitcell%vectorA,unitcell%vectorA))
    lb = dsqrt(dot_product(unitcell%vectorB,unitcell%vectorB))
    lc = dsqrt(dot_product(unitcell%vectorC,unitcell%vectorC))

    if( unit == "internal" ) then
       if( .not. unitcell%set ) then
          call XML_error("<unitcell> not found")
       end if
       atom%position(1) = atom%position(1)*la
       atom%position(2) = atom%position(2)*lb
       atom%position(3) = atom%position(3)*lc
!      atom%position = &
!           + atom%position(1) * unitcell%vectorA &
!           + atom%position(2) * unitcell%vectorB &
!           + atom%position(3) * unitcell%vectorC 
    else
       atom%position = atom%position * XML_TO_AU(unit)
    end if

    ! get <velocity> node
    node => getFirstElementByTagName(atom_node,"velocity")
    if( .not. associated(node) ) then
       atom%velocity_set = .false.
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       endif
       value = getChildValue(node)
       read(unit=value,fmt=*) atom%velocity(1:3)
       atom%velocity = atom%velocity * XML_TO_AU(unit)
       atom%velocity_set = .true.
    endif

    ! get <force> node
    node => getFirstElementByTagName(atom_node,"force")
    if( .not. associated(node) ) then
       atom%force_set = .false.
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       endif
       value = getChildValue(node)
       read(unit=value,fmt=*) atom%force(1:3)
       atom%force = atom%force * XML_TO_AU(unit)

       atom%force_set = .true.
    endif
!
    ! get <population> node
    node => getFirstElementByTagName(atom_node,"population")
    if( .not. associated(node) ) then
       atom%population_set = .false.
       atom%population=0.0d0
    else
       value = getChildValue(node)
       atom%population_set = .true.
       read(unit=value,fmt=*) atom%population
    endif
!
    ! get <population_guess> node
    node => getFirstElementByTagName(atom_node,"population_guess")
    if( .not. associated(node) ) then
       atom%population_guess_set = .false.
       atom%population_guess=0.0d0
    else
       value = getChildValue(node)
       atom%population_guess_set = .true.
       read(unit=value,fmt=*) atom%population_guess
    endif
!
    ! get <custumize> nodes
    vcustumize_node => getElementsByTagName(atom_node,"custumize")
    atom%ncustumize = getLength(vcustumize_node)
    if( atom%ncustumize > 0 ) then
       allocate( atom%vcustumize(atom%ncustumize) )

       do i=1, atom%ncustumize
          node => item(vcustumize_node,i-1)

          value = getAttribute(node,"base")
          read(value,*) atom%vcustumize(i)%base(1:4)

          value = getAttribute(node,"status")
          atom%vcustumize(i)%status = value
       end do
    end if

    return
  end subroutine atom_load

  ! save structure data in a file
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine structure_save( structure, filename, append )
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: append
    integer, parameter :: fd = 1
    integer :: j
    logical, save :: first = .true.

    if( (.not. present(append)) .or. (.not. append) ) then
       open(fd,file=filename)
    elseif( first ) then
       open(fd,file=filename)
       first = .false.
    else
       open(fd,file=filename,position='append')
    end if

    write(fd,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(fd,'(a,a,a,I20,a)') '<structure name="', trim(structure%name), &
         '" mdstep="', structure%mdstep, '">'

    call unitcell_save( fd, structure%unitcell )
    call heatbath_save( fd, structure%heatbath )

    do j=1, structure%natom
       call atom_save( fd, structure%vatom(j), structure%unitcell )
    end do

    write(fd,*) '</structure>'

    close(fd)

    return
  end subroutine structure_save

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine unitcell_save( fd, unitcell )
    implicit none
    integer, intent(in) :: fd
    type(unitcell_type), intent(in) :: unitcell

    write(fd,*) ""
    write(fd,*) "<unitcell>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', unitcell%vectorA(1:3), "</vector>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', unitcell%vectorB(1:3), "</vector>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', unitcell%vectorC(1:3), "</vector>"
    write(fd,*) "</unitcell>"
    write(fd,*) ""

    return
  end subroutine unitcell_save

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine heatbath_save( fd, heatbath )
    implicit none
    integer, intent(in) :: fd
    type(heatbath_type), intent(in) :: heatbath

    write(fd,*) ""
    write(fd,*) "<heatbath>"
    write(fd,'(a27,1e23.15,a14)') '  <massperatom unit="a.u.">', heatbath%massperatom, "</massperatom>"
    write(fd,'(a24,1e23.15,a14)') '  <position unit="a.u.">', heatbath%position, "</position>"
    write(fd,'(a24,1e23.15,a14)') '  <velocity unit="a.u.">', heatbath%velocity, "</velocity>"
    write(fd,*) "</heatbath>"
    write(fd,*) ""

    return
  end subroutine heatbath_save

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine atom_save( fd, atom, unitcell )
    implicit none
    integer, intent(in) :: fd
    type(atom_type), intent(in) :: atom
    type(unitcell_type), intent(in) :: unitcell
    integer :: k
    real(8) :: a, b, c
    real(8) :: la, lb, lc
!
    la = dsqrt(dot_product(unitcell%vectorA,unitcell%vectorA))
    lb = dsqrt(dot_product(unitcell%vectorB,unitcell%vectorB))
    lc = dsqrt(dot_product(unitcell%vectorC,unitcell%vectorC))
!
    a = atom%position(1) / la
    b = atom%position(2) / lb
    c = atom%position(3) / lc
!
!   a = atom%position(1) / unitcell%vectorA(1)
!   b = atom%position(2) / unitcell%vectorB(2)
!   c = atom%position(3) / unitcell%vectorC(3)

    write(fd,*) '<atom element="', trim(atom%name), '" ', &
         'class="', trim(atom%class), '" ', &
         'motion="', trim(atom%motion), '">'
    write(fd,'(a28,3e23.15,a11)') '  <position unit="internal">', a, b, c, '</position>'
    write(fd,'(a24,3e23.15,a11)') '  <velocity unit="a.u.">', atom%velocity(1:3), '</velocity>'
    write(fd,'(a24,3e23.15,a11)') '  <force    unit="a.u.">', atom%force(1:3), '</force>'

    if( atom%ncustumize > 0 ) then
       do k=1, atom%ncustumize
          write(fd,'(a19,4f16.10,a10,a16,a4)') &
               '  <custumize base="', atom%vcustumize(k)%base(1:4), &
               '" status="',  trim(atom%vcustumize(k)%status), &
               '" />'
       end do
    end if

    write(fd,*) '</atom>'
    write(fd,*) ""

    return
  end subroutine atom_save


end module elses_xml_structure

module elses_xml_config
  use M_config
  use flib_dom
  use elses_xml_misc
  use elses_xml_structure
  implicit none

  private
  public :: config_load

contains

  ! parse <config> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine config_load( config, filename )
    implicit none
    type(config_type), intent(out) :: config
    character(len=*),  intent(in)  :: filename
    character(len=256)   :: cwd ! current working dir name for backup
    character(len=256)   :: nwd ! current working dir name of config.xml

    type(fnode), pointer :: document_node
    type(fnode), pointer :: config_node
    type(fnode), pointer :: system_node
    type(fnode), pointer :: calc_node
    type(fnode), pointer :: output_node

    integer              :: ierr
    character(len=256)   :: value
    logical              :: ex
    integer              :: i_verbose, log_unit
    real(8), parameter   :: elses_xml_version_current=5.01d0
!
    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit
!
    inquire( file=trim(filename), exist=ex )
    if( .not. ex ) then
       write(*,'(a,a)') '# Error!: config_load : can not open file ', &
            trim(filename)
       stop
    end if

    document_node => parsefile(filename)
    call normalize(document_node)

    call getcwd( cwd )
    nwd = filename(1:scan(filename,"/\\",back=.true.))
    call chdir(nwd)

    ! get <config> node
    config_node => getFirstElementByTagName(document_node,"config")
    if( .not. associated(config_node) ) then
       call XML_error("<config> not found")
    endif

    ! get name attribute
    value = getAttribute(config_node,"name")
    if( value == "" ) then
       call XML_error("<config name> not found")
    endif
    config%name = value
    if (log_unit > 0) write(log_unit,'(a,a)')'INFO-XML:config name=',trim(config%name)

    ! get elses_xml_version attribute
    value = getAttribute(config_node,"elses_xml_version")
    if ( value == "" ) then
      config%elses_xml_version=1.0d0
    else
      if ( value == "current" ) then
        config%elses_xml_version=elses_xml_version_current
        if (log_unit > 0) write(log_unit,'(a)')'INFO-XML:ELSES XML version= (current)'
      else
        read(unit=value,fmt=*,iostat=ierr) config%elses_xml_version
        if ( ierr /= 0 ) then
          if (log_unit > 0) write(log_unit,'(a,a)')'ERROR-XML(elses_xml_version):',trim(value)
          write(*,'(a,a)')'ERROR-XML(elses_xml_version):',trim(value)
          stop
        endif
      endif
    endif
    if (log_unit > 0) write(log_unit,'(a,f10.4)')'INFO-XML:ELSES XML version=',config%elses_xml_version

    ! get <system> node
    system_node => getFirstElementByTagName(config_node,"system")
    if( .not. associated(system_node) ) then
       call XML_error("<system> not found")
    else
       call system_load( config%system, system_node, config%elses_xml_version )
    end if

    ! get <calc> node
    calc_node => getFirstElementByTagName(config_node,"calc")
    if( .not. associated(calc_node) ) then
       call calc_default( config%calc )
    else
       call calc_default( config%calc )
       call calc_load( config%calc, calc_node )
    end if

    ! get <output> node
    output_node => getFirstElementByTagName(config_node,"output")
    if( .not. associated(output_node) ) then
       call output_default( config%output )
    else
       call output_default( config%output )
       call output_load( config%output, output_node )
    endif

    call destroyNode(document_node)

    call chdir(cwd)

    return
  end subroutine config_load



  ! parse <system> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine system_load( system, system_node, elses_xml_version )
    use M_sax_parser,      only : struc_load_sax
    implicit none
    type(system_type), intent(out) :: system
    type(fnode), pointer :: system_node
    real(8),           intent(in)  :: elses_xml_version


    type(fnodeList), pointer :: vtarget_node
    type(fnodeList), pointer :: velement_node

    integer              :: i
    type(fnode), pointer :: node
    type(fnode), pointer :: node_wrk
    character(len=256)   :: value, unit
    integer              :: ierr
!
    character(len=8)       :: name
    character(len=64)      :: model, filename
!
    integer              :: i_verbose, log_unit
    real(8)              :: time_wrk, time_wrk_prev, time_period
!
    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    ! get <split_input_file> node
    node => getFirstElementByTagName(system_node,"split_input_file")
    if( .not. associated(node) ) then
      node => getFirstElementByTagName(system_node,"splitted_input_file")
    endif  
    if( .not. associated(node) ) then
      system%structure%split%set                = .false.
      system%structure%split%file_index         = -1 ! dummy value
      system%structure%split%number_of_files    = -1 ! dummy value
      system%structure%split%atoms_in_this_file = -1 ! dummy value
      system%structure%split%atom_initial       = -1 ! dummy value
      system%structure%split%atom_final         = -1 ! dummy value
    else
      system%structure%split%set = .true.
      value = getAttribute(node,"number_of_files")
      read(unit=value, fmt=*) system%structure%split%number_of_files
      if (system%structure%split%number_of_files == 0) system%structure%split%set = .false.
    endif
    if (system%structure%split%set) then
      if (log_unit >  0) write(log_unit,'(a,i10)')'INFO-XML-SPLIT: Optional attribute detected : number_of_files =', & 
&                              system%structure%split%number_of_files
    endif
!   
    ! get <temperature> node
    node => getFirstElementByTagName(system_node,"temperature")
    if( .not. associated(node) ) then
!      call XML_error("<temperature> not found")
       if (log_unit >  0) write(log_unit,'(a)')'INFO-XML-WARNING : NO tag detected : temperature'
       system%temperature = -1.0d0 ! dummy value
    else   
       unit  = getAttribute(node,"unit")
       if( unit=="" ) then
          unit="kelvin"
       endif
       value = getChildValue(node)
       read(unit=value,fmt=*) system%temperature
       system%temperature = system%temperature * XML_TO_AU(unit)
    endif

    ! get <element> node
    velement_node => getElementsByTagName(system_node,"element")
    system%structure%nelement = getLength(velement_node)
!
    ! get <use_matom> node
    node => getFirstElementByTagName(system_node,"use_matom")
    if( associated(node) ) then
      system%structure%use_matom = .true.
      system%structure%use_vatom = .false.
      if (log_unit >  0) write(log_unit,'(a)')'INFO-XML: Optional tag : use_matom'
    else
      system%structure%use_matom = .false.
      system%structure%use_vatom = .true.
    endif
!   
    ! get <cluster> node
    node => getFirstElementByTagName(system_node,"cluster")
    if( .not. associated(node) ) then
       call XML_error("<cluster> not found")
    else
       ! get number_of_atoms attribute
       value = getAttribute(node,"number_of_atoms")
       if( value == "" ) then
         system%structure%natom = 0  ! dummy value 
       else
         read(unit=value,fmt=*) system%structure%natom
         if (log_unit >  0) write(log_unit, '(a,i10)') 'INFO-XML: Optional attribute detected : number_of_atom =', & 
&                              system%structure%natom
       endif   
       ! get parser attribute
       value = getAttribute(node,"parser")
       if (( value == "SAX" ) .or. ( value == "sax" )) then
         system%structure%parser='sax'
       else
         system%structure%parser='dom'
       endif   
       if (log_unit >  0) write(log_unit, '(a,a)')'INFO-XML:XML parser= ', trim(system%structure%parser)
       if (system%structure%parser == 'dom') then
         if (log_unit >  0) then 
           write(log_unit, '(a,a)')'WARNING:DOM PARSER is used for reading the input XML file.'
           write(log_unit, '(a,a)')'WARNING:DOM PARSER may be impractical for large systems, such as 100,000-atom systems.'
         endif
       endif  
!
       ! get tag_dump attribute
       value = getAttribute(node,"tag_dump")
       if ( value == "on" ) then
         system%structure%tag_dump= .true.
         if (log_unit >  0) write(log_unit, '(a)')'INFO-XML: Optional attribute detected : tag_dump'
       else
         system%structure%tag_dump= .false.
       endif   
       ! get read_mode attribute
       value = getAttribute(node,"read_mode")
       if ( value == "" ) then
         system%structure%read_mode= "default"
       else
         value=trim(value) 
         read(unit=value,fmt=*) system%structure%read_mode
         if (log_unit >  0) write(log_unit, '(a,a)')'INFO-XML: Optional attribute detected : read_mode=', &
&                           system%structure%read_mode
       endif   
!
       ! get structure attribute
       call get_system_clock_time_loc(time_wrk)
       time_wrk_prev=time_wrk
       value = getAttribute(node,"structure")
       if( value == "" ) then
          call XML_error("<cluster structure> not found")
       else
          if (system%structure%parser == "sax") then
            call struc_load_sax( value )
          else
            call structure_load( system%structure, value )
            do i = 1, system%structure%nelement
              call structure_element_load( system%structure, item(velement_node,i-1) )
            end do
          endif  
       end if
       call get_system_clock_time_loc(time_wrk)
       call measure_clock_time_period_loc(time_wrk, time_wrk_prev, time_period)
       if (log_unit >  0) write(log_unit, '(a,f10.5)')'TIME:Read_structure_XML_file = ', time_wrk-time_wrk_prev
       time_wrk_prev=time_wrk
!
       ! get cutoff_radius attribute
       value = getAttribute(node,"cutoff_radius")
       if( value /= "" ) then
          read(unit=value,fmt=*) system%cutoff_radius
          if (log_unit >  0) write(log_unit, '(a,f10.5)')'INFO-XML:Optional tag detected : cutoff_radius [au] =', & 
&                              system%cutoff_radius
       else
          system%cutoff_radius=huge(0d0) ! default value of the cutoff radius is huge(0d0)
       end if
    end if
!
    ! get <element> node, additional procedures only in SAX mode
!
    if (system%structure%parser == "sax") then
      if (log_unit >  0) write(log_unit, '(a,i5)') 'nelement is set to be:nelement=',system%structure%nelement
      allocate (config%system%structure%velement(system%structure%nelement), stat=ierr )
      if (ierr /= 0) then
        write(*,*) 'Alloc error. velement'
        stop 
      endif   
      do i = 1, system%structure%nelement
        name = getAttribute(item(velement_node,i-1), "name")
        filename = getAttribute(item(velement_node,i-1), "filename")
        model = getAttribute(item(velement_node,i-1), "model")
        if (name == '') then
          write(*,*) 'ERROR(minnsing element name)'
          stop
        else
          system%structure%velement(i)%name=name
        endif   
        if (filename == '') then
          write(*,*) 'ERROR(minnsing element filename)'
          stop
        else
          system%structure%velement(i)%filename=filename
        endif   
        if (model == '') then
          write(*,*) 'ERROR(minnsing element model)'
          stop
        else
          system%structure%velement(i)%quantum%type=model
        endif   
      enddo   
    endif  
!
    ! get <boundary> node
    node => getFirstElementByTagName(system_node,"boundary")
    if( .not. associated(node) ) then
       if (elses_xml_version+1.0d-6 > 5.01d0) then 
         if (log_unit >  0) write(log_unit, '(a)') 'ERROR in config XML file: boundary tag missing'
         write(*,*) 'ERROR in config XML file: boundary tag missing'
         stop
       else
         system%boundary%periodic_x = .true.
         system%boundary%periodic_y = .true.
         system%boundary%periodic_z = .true.
         if (log_unit >  0) write(log_unit, '(a)') 'WARNING-XML: boundary tag missing. Periodic boundary condition is assumed'
       endif
    else
       value = getAttribute(node,"x")
       if( value == "nonperiodic" ) then
          system%boundary%periodic_x = .false.
       else
          system%boundary%periodic_x = .true.
       end if

       value = getAttribute(node,"y")
       if( value == "nonperiodic" ) then
          system%boundary%periodic_y = .false.
       else
          system%boundary%periodic_y = .true.
       end if

       value = getAttribute(node,"z")
       if( value == "nonperiodic" ) then
          system%boundary%periodic_z = .false.
       else
          system%boundary%periodic_z = .true.
       end if
    end if

    ! get <heatbath_mass> node
    node => getFirstElementByTagName(system_node,"heatbath_mass")
    if ( associated(node) ) then
      value = getAttribute(node,"mode")
      if (log_unit >  0) then 
        write(log_unit, '(a,a)') 'INFO-XML:Optional tag detected : heatbath mass mode=',trim(value)
      endif       
      if (trim(value) /= "off") then
        node_wrk => getFirstElementByTagName(node,"mass_per_atom")
        if ( associated(node_wrk) ) then
          unit  = getAttribute(node_wrk,"unit")
          if( unit /= "" ) unit="a.u." 
          endif
          value = getChildValue(node_wrk)
          read(unit=value,fmt=*) system%structure%heatbath%massperatom
          system%structure%heatbath%massperatom = system%structure%heatbath%massperatom * XML_TO_AU(unit)
          if (log_unit >  0) then 
           write(log_unit, '(a,f10.5)') 'INFO-XML:Optional tag detected : heatbath mass per atom [au]=', &
&             system%structure%heatbath%massperatom
           write(log_unit, '(a)') 'INFO-XML:The tag of heatbath mass per atom in the structure XML file IS IGNORED'
        endif
      endif  
    endif

    ! get <target> nodes
    vtarget_node => getElementsByTagName(system_node,"target")
    system%ntarget = getLength(vtarget_node)
    allocate( system%vtarget(system%ntarget) )

    do i=1, system%ntarget
       call target_load( system%vtarget(i), item(vtarget_node,i-1) )
    end do

    return
  end subroutine system_load

  ! parse <target> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine target_load( target, target_node )
    implicit none
    type(target_type), intent(out) :: target
    type(fnode), pointer :: target_node

    type(fnode), pointer :: node
    character(len=256)   :: value, unit

    ! get class attribute
    value = getAttribute(target_node,"class")
    if( value == "" ) then
       target%class = ""
       target%class_set = .false.
    else
       read(unit=value,fmt=*) target%class
       target%class_set = .true.
    end if

    ! get xrange attribute
    value = getAttribute(target_node,"xrange")
    if( value == "" ) then
       target%xrange_set = .false.
    else
       unit  = getAttribute(target_node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       read(unit=value,fmt=*) target%xrange(1:2)
       target%xrange = target%xrange * XML_TO_AU(unit)
       target%xrange_set = .true.
    end if

    ! get yrange attribute
    value = getAttribute(target_node,"yrange")
    if( value == "" ) then
       target%yrange_set = .false.
    else
       unit  = getAttribute(target_node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       read(unit=value,fmt=*) target%yrange(1:2)
       target%yrange = target%yrange * XML_TO_AU(unit)
       target%yrange_set = .true.
    end if

    ! get zrange attribute
    value = getAttribute(target_node,"zrange")
    if( value == "" ) then
       target%zrange_set = .false.
    else
       unit  = getAttribute(target_node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       read(unit=value,fmt=*) target%zrange(1:2)
       target%zrange = target%zrange * XML_TO_AU(unit)
       target%zrange_set = .true.
    end if

    ! get trange attribute
    value = getAttribute(target_node,"trange")
    if( value == "" ) then
       target%trange_set = .false.
    else
       unit  = getAttribute(target_node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       read(unit=value,fmt=*) target%trange(1:2)
       target%trange = target%trange * XML_TO_AU(unit)
       target%trange_set = .true.
    end if


    ! get <position> node
    node => getFirstElementByTagName(target_node,"position")
    if( .not. associated(node) ) then
       target%position_set = .false.
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) target%position(1:3)
       target%position = target%position * XML_TO_AU(unit)
       target%position_set = .true.
    end if

    ! get <velocity> node
    node => getFirstElementByTagName(target_node,"velocity")
    if( .not. associated(node) ) then
       target%velocity_set = .false.
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) target%velocity(1:3)
       target%velocity = target%velocity * XML_TO_AU(unit)
       target%velocity_set = .true.
    end if

    ! get <temperature> node
    node => getFirstElementByTagName(target_node,"temperature")
    if( .not. associated(node) ) then
       target%temperature_set = .false.
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) target%temperature
       target%temperature = target%temperature * XML_TO_AU(unit)
       target%temperature_set = .true.
    end if

    ! get <mass> node
    node => getFirstElementByTagName(target_node,"mass")
    if( .not. associated(node) ) then
       target%mass_set = .false.
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) target%mass
       if( unit /= "amu" ) then
          target%mass = target%mass * XML_TO_AU(unit)
       end if
       target%mass_set = .true.
    end if

    ! get <motion> node
    node => getFirstElementByTagName(target_node,"motion")
    if( .not. associated(node) ) then
       target%motion_set = .false.
    else
       value = getChildValue(node)
       read(unit=value,fmt=*) target%motion
       target%motion_set = .true.
    end if

    return
  end subroutine target_load


  ! parse <calc> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_load( calc, calc_node )
    implicit none
    type(calc_type), intent(out) :: calc
    type(fnode), pointer :: calc_node

    type(fnode), pointer :: limit_node
    type(fnode), pointer :: dynamics_node
    type(fnode), pointer :: solver_node
    type(fnode), pointer :: distributed_node
    type(fnode), pointer :: genoOption_node
    type(fnode), pointer :: snapshot_node
    character(len=256)   :: value

    ! added by S.Y Nov 19, 2008
    ! get type attribute
    type(fnode), pointer :: optimization_node
    type(fnode), pointer :: cell_change_node
    type(fnode), pointer :: calc_check_node
    type(fnode), pointer :: calc_force_node
    type(fnode), pointer :: calc_integer_elec_num_node
    type(fnode), pointer :: calc_virial_pressure_node
    type(fnode), pointer :: interaction_range_node
    type(fnode), pointer :: work_node1
    type(fnode), pointer :: work_node2
    type(fnode), pointer :: group_id_node
!
    integer              :: log_unit
!
    log_unit=config%calc%distributed%log_unit
!
    value = getAttribute(calc_node,"mode")
    calc%mode = trim(adjustl(value))
    select case (calc%mode)
    case ("optimization")
       calc%optimization%scheme       = "steepest_descent"
       calc%optimization%max_num_iter = 100
       calc%optimization%sd_ratio     = 1d-1
       calc%optimization%convergence_mode   = "off"
       calc%optimization%energy_convergence = -1.0d0
       calc%optimization%force_convergence  = -1.0d0
       optimization_node => getFirstElementByTagName(calc_node,"optimization")
       if( associated(optimization_node) ) &
            call optimization_load( calc%optimization, optimization_node )
       ! return  ! In optimization mode, there are no more tags in this layer.
    case ("dynamics")
       calc%mode = "dynamics" 
    case ("conversion", "file conversion", "file_conversion")
       calc%mode = "conversion"
    case ("snapshot", "snapshots", "given_snapshot", "given_snapshots")
       calc%mode = "snapshot"
    case ("cell_change_only")
       calc%mode = "cell_change_only"
    case ("matrix_generation")
       calc%mode = "matrix_generation"
! S.Y Sep 21, 2013 two cases "calc%mode = 'matrix_generation'" and "calc%mode = 'write_Hamiltonian'" 
!      should be merged into single case. But before merging, compatibility check should be done.
!      related file : elses-qm-engine.f90
    case ("write_Hamiltonian")
       calc%mode = "optimization"
       calc%optimization%scheme       = "none"
       calc%optimization%max_num_iter = 0
       calc%optimization%convergence_mode   = "off"
    case default
       write(*,*) 'mode=',trim(calc%mode)
       call XML_error("No mode or unknown mode in <calc> tag")
    end select

    ! get <genoOption> node
    genoOption_node => getFirstElementByTagName(calc_node,"genoOption")
    calc%genoOption%Use_VOIP               = .false.
    calc%genoOption%genoMethod             = "STANDARD" !! "ICON" "NOLINEAR"
    calc%genoOption%CSC_method             = "ELSTNER"  !! "TB0"
    calc%genoOption%CSC_max_loop_count     = 0
    calc%genoOption%CSC_charge_convergence = 1d-6
    calc%genoOption%CSC_charge_mixing_ratio= 1d-1
    calc%genoOption%CSC_mode_for_tuning_mixing_ratio= 2
    calc%genoOption%HML_K                  = 1.75d0
    calc%genoOption%HML_small_delta        = 0.35d0 / XML_TO_AU("angstrom")
    calc%genoOption%HML_kappa              = 1.00d0
    calc%genoOption%HML_constant_K         = .false.
    calc%genoOption%set_S_as_I             = .false.
    if( associated(genoOption_node) ) &
       call genoOption_load( calc%genoOption, genoOption_node )

    ! get <limit> node
    limit_node => getFirstElementByTagName(calc_node,"limit")
    if( .not. associated(limit_node) ) then
       call limit_default( calc%limit )
    else
       call limit_default( calc%limit )
       call limit_load( calc%limit, limit_node )
    end if

    ! get <dynamics>
    dynamics_node => getFirstElementByTagName(calc_node,"dynamics")
    if( .not. associated(dynamics_node) ) then
       call dynamics_default( calc%dynamics )
    else
       call dynamics_default( calc%dynamics )
       call dynamics_load( calc%dynamics, dynamics_node )
    end if

    ! get <solver> node
    solver_node => getFirstElementByTagName(calc_node,"solver")
    if( .not. associated(solver_node) ) then
       call solver_default( calc%solver )
    else
       call solver_default( calc%solver )
       call solver_load( calc%solver, solver_node )
    end if

    ! get <distributed> node
    distributed_node => getFirstElementByTagName(calc_node,"distributed")
    call distributed_default( calc%distributed )
    if( associated(distributed_node) ) then
       call distributed_load( calc%distributed, distributed_node )
    end if

    ! get <snapshot> node
    snapshot_node => getFirstElementByTagName(calc_node,"snapshot")
    call snapshot_default( calc%snapshot )
    if( associated(snapshot_node) ) then
      call snapshot_load( calc%snapshot, snapshot_node )
    end if

    ! get <cell_change> node
    cell_change_node => getFirstElementByTagName(calc_node,"cell_change")
    call cell_change_default( calc%cell_change )
    if( associated(cell_change_node) ) then
     call cell_change_load( calc%cell_change, cell_change_node )
    end if

    ! get <calc_check> node
    calc_check_node => getFirstElementByTagName(calc_node,"calc_check")
    call load_calc_check( calc%calc_check, calc_check_node )

    ! get <calc_force_node> node
    calc_force_node => getFirstElementByTagName(calc_node,"calc_force")
    calc%calc_force_mode="on" ! default setting
    if ( associated(calc_force_node) ) then
      value = getAttribute(calc_force_node,"mode")
      if ( value == "off" ) then
        calc%calc_force_mode="off" 
      end if
      if (log_unit > 0) write(log_unit,'(a,a)') & 
        & 'INFO-XML:Optional tag: calc_force : mode=', trim(calc%calc_force_mode)
    end if

    ! get <use_integer_elec_num> node
    calc_integer_elec_num_node => getFirstElementByTagName(calc_node,"use_integer_elec_num")
    if ( associated(calc_integer_elec_num_node) ) then
      calc%use_integer_elec_num =.true. 
      if (getAttribute(calc_integer_elec_num_node,"mode") == "off") calc%use_integer_elec_num =.false. 
    else
      calc%use_integer_elec_num =.false. 
    endif    
    if (calc%use_integer_elec_num) then
      if (log_unit > 0) write(log_unit,'(a,a)') & 
        & 'INFO-XML:Optional tag: calc%use_integer_elec_num'
    endif   

    ! get <calc_virial_pressure> node
    calc_virial_pressure_node => getFirstElementByTagName(calc_node,"calc_virial_pressure")
    if ( associated(calc_virial_pressure_node) ) then
      calc%calc_virial_pressure =.true. 
      if (getAttribute(calc_virial_pressure_node,"mode") == "off") calc%calc_virial_pressure =.false. 
    else
      calc%calc_virial_pressure =.false. 
    endif    
    if (calc%calc_virial_pressure) then
      if (log_unit > 0) write(log_unit,'(a,a)') & 
        & 'INFO-XML:Optional tag: calc%calc_virial_pressure'
    endif   

    ! get <interaction_range> node
    calc%interaction_range%cutoff_rest_cellmax = .false. ! default value
    calc%interaction_range%cutoff_rest         =  -1.0d0 ! dummy value
    calc%interaction_range%cutoff_rest_non_vdW =  -1.0d0 ! dummy value
    interaction_range_node => getFirstElementByTagName(calc_node,"interaction_range")
    if ( associated(interaction_range_node) ) then
      work_node1 => getFirstElementByTagName(interaction_range_node,"cutoff_rest")
      if ( associated(work_node1) ) then
        if (getAttribute(work_node1,"mode") == "cell_max") then
          calc%interaction_range%cutoff_rest_cellmax = .true.
          if (log_unit > 0) write(log_unit,'(a,a)') & 
            & 'INFO-XML:Optional tag: calc%interaction_range%cutoff_rest = cell_max mode' 
        endif  
      endif   
    endif   

    ! get <use_group_id> node
    group_id_node => getFirstElementByTagName(calc_node,"use_group_id")
    calc%use_group_id = .false.
    if ( associated(group_id_node) ) then
      calc%use_group_id =.true.
      if (getAttribute(group_id_node,"mode") == "off") calc%use_group_id =.false.
    endif
    if ( calc%use_group_id ) then
       if (log_unit > 0) write(log_unit,'(a)') & 
         & 'INFO-XML:Optional tag: use_group_id'
    endif

    ! get <constraint_w_group> node
    work_node1 => getFirstElementByTagName(calc_node,"constraint_w_group")
    calc%constraint_w_group = .false.
    if ( associated(work_node1) ) then
      calc%constraint_w_group =.true.
      if (getAttribute(work_node1,"mode") == "off") calc%constraint_w_group =.false.
    endif
    if ( calc%constraint_w_group ) then
       if (log_unit > 0) write(log_unit,'(a)') & 
         & 'INFO-XML:Optional tag: constraint_w_group'
    endif

    return
  end subroutine calc_load


  ! parse <optimization> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine optimization_load( optimization, optimization_node )
    implicit none
    type(optimization_type), intent(out) :: optimization
    type(fnode), pointer                 :: optimization_node
    !
    type(fnode), pointer :: node        ! temporal 
    character(len=256)   :: value, unit
    integer              :: i_verbose, log_unit
    real(8)              :: unit_conv
    character(len=256)   :: value_chara
    real(8)              :: value_real

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    ! get scheme attribute 
    value = getAttribute(optimization_node,"scheme")
    if( value == "" ) then
       optimization%scheme = "steepest_descent"
    end if
    node  => getFirstElementByTagName(optimization_node,"max_num_iter")
    if( associated(node) ) then
       value =  getChildValue(node)
       read(unit=value,fmt=*) optimization%max_num_iter
       if (log_unit > 0) write(log_unit,'(a,i10)') 'INFO-XML:optimization%max_num_iter=', optimization%max_num_iter
    end if
    node  => getFirstElementByTagName(optimization_node,"sd_ratio")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) optimization%sd_ratio
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'INFO-XML:optimization%sd_ratio    =', optimization%sd_ratio
    end if
    node  => getFirstElementByTagName(optimization_node,"convergence_criteria")
    if ( associated(node) ) then
      value = getAttribute(node,"mode")
      if ( value /= "" ) then
        read(unit=value,fmt=*) optimization%convergence_mode
        if (log_unit > 0) write(log_unit,'(a,a)') & 
          & 'INFO-XML:Optional tag: optimization%convergence_mode=', optimization%convergence_mode
        unit  = getAttribute(node,"unit")
        if( unit == "" ) then
          unit_conv=1.0d0 
        else
          unit_conv=XML_TO_AU(unit)
        endif
        value_chara = getChildValue(node)
        read(unit=value_chara,fmt=*) value_real
        select case(trim(optimization%convergence_mode))
          case ("energy_per_atom")
            optimization%energy_convergence =  value_real*unit_conv
            if (log_unit > 0) write(log_unit,'(a,f30.20)') & 
                & 'INFO-XML:Optional tag: optimization%convergence_energy [au] =', optimization%energy_convergence
          case ("force_per_atom")
            optimization%force_convergence  =  value_real*unit_conv
            if (log_unit > 0) write(log_unit,'(a,f30.20)') & 
                & 'INFO-XML:Optional tag: optimization%convergence_force [au]  =', optimization%force_convergence
          case ("off", "none")
          case default
            write(*,*)'ERROR:optimization%convergence_mode=',trim(optimization%convergence_mode)
            stop
        end select
      endif   
    end if
  end subroutine optimization_load

  ! parse <genoOption> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine genoOption_load( genoOption, genoOption_node )
    implicit none
    type(genoOption_type), intent(out) :: genoOption
    type(fnode), pointer :: genoOption_node 
    !
    type(fnode), pointer :: node
    type(fnode), pointer :: node_wrk
    character(len=256)   :: value
    integer              :: i_verbose, log_unit
    character(len=256)   :: chara_mode
    integer              :: ierr
    integer              :: chara_length

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit

    node => getFirstElementByTagName(genoOption_node,"Use_VOIP")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%Use_VOIP
    end if

    node => getFirstElementByTagName(genoOption_node,"geno_method")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%genoMethod
    end if

    node => getFirstElementByTagName(genoOption_node,"CSC_method")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%CSC_METHOD
    end if

    node => getFirstElementByTagName(genoOption_node,"CSC_max_loop_count")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%CSC_MAX_LOOP_COUNT
    end if

    node => getFirstElementByTagName(genoOption_node,"CSC_charge_convergence")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%CSC_charge_convergence
    end if

    node => getFirstElementByTagName(genoOption_node,"HML_constant_K")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%HML_constant_K
    end if

    if(genoOption%HML_constant_K) then
       write(*,'(" INFO-XML:genoOption_load: The value of K is fixed.")')
       node => getFirstElementByTagName(genoOption_node,"HML_K")
       if( associated(node) ) then
          value = getChildValue(node)
          read(unit=value,fmt=*) genoOption%HML_K
       end if
       if (log_unit > 0) write(log_unit,'(" INFO-XML:genoOption_load: genoOption%HML_K=",ES14.7)') &
&           genoOption%HML_K
    else
       node => getFirstElementByTagName(genoOption_node,"HML_small_delta")
       if( associated(node) ) then
          value = getChildValue(node)
          read(unit=value,fmt=*) genoOption%HML_small_delta
          genoOption%HML_small_delta = genoOption%HML_small_delta &
               / XML_TO_AU("angstrom")
       end if

       node => getFirstElementByTagName(genoOption_node,"HML_kappa")
       if( associated(node) ) then
          value = getChildValue(node)
          read(unit=value,fmt=*) genoOption%HML_kappa
       end if
       if (log_unit > 0) then 
         write(log_unit,'("INFO-XML:genoOption_load: genoOption%HML_small_delta=",ES14.7)') &
           & genoOption%HML_small_delta
         write(log_unit,'("INFO-XML:genoOption_load: genoOption%HML_kappa      =",ES14.7)') &
           & genoOption%HML_kappa
       endif  
    end if

    node => getFirstElementByTagName(genoOption_node,"CSC_charge_mixing_ratio")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%CSC_charge_mixing_ratio
    end if

    node => getFirstElementByTagName(genoOption_node,"CSC_mode_for_tuning_mixing_ratio")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) genoOption%CSC_mode_for_tuning_mixing_ratio
       if (log_unit > 0) write(log_unit,'(a,i10)') & 
&               'INFO-XML: Optional tag detected : CSC_mode_for_tuning_mixing_ratio =', &
&                              genoOption%CSC_mode_for_tuning_mixing_ratio
    end if
!
!   ! get <set_S_as_I> node
    node => getFirstElementByTagName(genoOption_node,"set_S_as_I")
    if( associated(node) ) then
      genoOption%set_S_as_I = .true.
      if (getAttribute(node,"mode") == "off") genoOption%set_S_as_I = .false.
      if ( genoOption%set_S_as_I ) then
        if (log_unit > 0) write(log_unit,*)'INFO-XML: optional tag: set_S_as_I =', genoOption%set_S_as_I
      endif   
    endif
!   
    ! get <van_der_Waals> node
    genoOption%vanDerWaals        =  .false.
    genoOption%vanDerWaals_lambda = - 1.0d0 ! dummy value
    node => getFirstElementByTagName(genoOption_node,"van_der_Waals")
    if( associated(node) ) then
       value = getChildValue(node)
       chara_mode  = getAttribute(node,"mode")
       genoOption%vanDerWaals = .true.
       if( chara_mode == "off" ) then
         genoOption%vanDerWaals = .false.
       endif
       if ( genoOption%vanDerWaals ) then
          if (log_unit > 0) write(log_unit,*)'INFO-XML: optional tag: van_der_Waals (vdW) =', genoOption%vanDerWaals
          node_wrk => getFirstElementByTagName(node,"vdW_lambda")
          chara_mode  = getAttribute(node_wrk,"mode")
          value = getChildValue(node_wrk)
          if (trim(value) /= "") then
            if (trim(chara_mode) /= "default") read(unit=value,fmt=*) genoOption%vanDerWaals_lambda
            if (log_unit > 0) write(log_unit,*)'INFO-XML: optional tag: vdW_lambda =', genoOption%vanDerWaals_lambda
          endif   
       endif
    endif

  end subroutine genoOption_load

  ! parse <limit> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine limit_load( limit, limit_node )
    implicit none
    type(limit_type), intent(out) :: limit
    type(fnode), pointer :: limit_node
    type(fnode), pointer :: node
    character(len=256)   :: value, unit
    integer :: lu
    
    lu=config%calc%distributed%log_unit

    ! get <time> node
    node => getFirstElementByTagName(limit_node,"time")
    if( associated(node) ) then
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "second"
       end if

       value = getChildValue(node)
       read(unit=value,fmt=*) limit%time

       select case(unit)
       case("second","seconds","sec","s")
       case("minute","minutes","min","m")
          limit%time = limit%time * 60.0
       case("hour","hours","h")
          limit%time = limit%time * 60*60
       case("day","days","d")
          limit%time = limit%time * 60*60*24
       end select
       if (lu > 0) write(lu,*) 'INFO-XML:Optional tag detected : (time limit)[sec]  =', &
&                                         limit%time
       if (lu > 0) write(lu,*) 'INFO-XML:Optional tag detected : (time limit)[min]  =', &
&                                         limit%time/60.0d0
       if (lu > 0) write(lu,*) 'INFO-XML:Optional tag detected : (time limit)[h  ]  =', &
&                                         limit%time/60.0d0/60.0d0
    end if

    ! get <memory> node
    node => getFirstElementByTagName(limit_node,"memory")
    write(*,*)'read memory'
    if (associated(node) ) then
       unit  = getAttribute(node,"unit")
       if ( unit == "" ) then
          unit = "GB"
       end if
       if ( unit /= "GB" ) then
         write(*,*)'ERROR!:Unknown unit in memory limit=',trim(unit)
         stop
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) limit%memory
       if (lu > 0) write(lu,'(a,f15.10)') 'INFO-XML:Optional tag detected : (memory limit)[GB] =', &
&                                         limit%memory
!      select case(unit)
!      case("byte","bytes","B")
!      case("kbyte","kbytes","kB")
!         limit%memory = limit%memory * 1024
!      case("Mbyte","Mbytes","MB")
!         limit%memory = limit%memory * 1024*1024
!      case("Gbyte","Gbytes","GB")
!         limit%memory = limit%memory * 1024*1024*1024
!      end select
    end if

    return
  end subroutine limit_load

  ! parse <dynamics> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine dynamics_load( dynamics, dynamics_node )
    implicit none
    type(dynamics_type), intent(out) :: dynamics
    type(fnode), pointer :: dynamics_node

    type(fnode), pointer :: node
    character(len=256)   :: value, unit

    ! get scheme attribute
    value = getAttribute(dynamics_node,"scheme")
    if( value == "" ) then
       dynamics%scheme = ""
    else
       dynamics%scheme = value
    end if

    ! get <delta> node
    node => getFirstElementByTagName(dynamics_node,"delta")
    if( .not. associated(node) ) then
       call XML_error("<delta> not found")
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) dynamics%delta
       dynamics%delta = dynamics%delta * XML_TO_AU(unit)
    end if

    ! get <total> node
    node => getFirstElementByTagName(dynamics_node,"total")
    if( .not. associated(node) ) then
       call XML_error("<total> not found")
    else
       unit  = getAttribute(node,"unit")
       if( unit == "" ) then
          unit = "a.u."
       end if
       value = getChildValue(node)
       read(unit=value,fmt=*) dynamics%total
       dynamics%total = dynamics%total * XML_TO_AU(unit)
    end if

    return
  end subroutine dynamics_load


  ! parse <solver> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine solver_load( solver, solver_node )
    implicit none
    type(solver_type), intent(out) :: solver
    type(fnode), pointer :: solver_node

    type(fnode), pointer :: node
    type(fnode), pointer :: node_tmp, node_tmp2
    character(len=256)   :: value
    integer              :: i_verbose, log_unit

    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit
!
    ! get scheme attribute
    value = getAttribute(solver_node,"scheme")
    if( value == "" ) then
       solver%scheme = ""
    else
       solver%scheme = value
       if (log_unit > 0) write(log_unit,'(2a)')'INFO-XML:Optional tag detected : solver scheme = ',solver%scheme
    end if

    ! get <eigen_mpi> node
    solver%eigen_mpi%SEP_solver = "default"        ! dummy value (default)
    solver%eigen_mpi%GS_transformation = "default" ! dummy value (default)
    solver%eigen_mpi%blocksize     = -1                ! dummy value (default)
    solver%eigen_mpi%level_lowest  = -1                ! dummy value (default)
    solver%eigen_mpi%level_highest = -1                ! dummy value (default)
    node => getFirstElementByTagName(solver_node,"eigen_mpi")
    if ( associated(node) ) then
      node_tmp => getFirstElementByTagName(node,"SEP_solver")
      if ( associated(node_tmp) ) then
        value = getChildValue(node_tmp)
        if ( value /= "" ) then
          solver%eigen_mpi%SEP_solver = trim(adjustl(value))
          if (log_unit > 0) write(log_unit,'(2a)') &
           'INFO-XML:Optional tag : eigen_mpi%SEP_solver =', trim(solver%eigen_mpi%SEP_solver)
        endif
      end if
!
      node_tmp => getFirstElementByTagName(node,"GS_transformation")
      if ( associated(node_tmp) ) then
        value = getChildValue(node_tmp)
        if ( value /= "" ) then
          solver%eigen_mpi%GS_transformation = trim(adjustl(value))
          if (log_unit > 0) write(log_unit,'(2a)') &
             'INFO-XML:Optional tag : eigen_mpi%GS_transformation =', trim(solver%eigen_mpi%GS_transformation)
        endif
      end if
!
      node_tmp => getFirstElementByTagName(node,"blocksize")
      if ( associated(node_tmp) ) then
        value = getChildValue(node_tmp)
        if ( value /= "" ) then
          read(unit=value,fmt=*) solver%eigen_mpi%blocksize
          if (log_unit > 0) write(log_unit,'(a,i10)') &
           'INFO-XML:Optional tag : eigen_mpi%blocksize =', solver%eigen_mpi%blocksize
        endif
      end if
!
      node_tmp => getFirstElementByTagName(node,"level_lowest")
      if ( associated(node_tmp) ) then
        value = getChildValue(node_tmp)
        if ( value /= "" ) then
          read(unit=value,fmt=*) solver%eigen_mpi%level_lowest
          if (log_unit > 0) write(log_unit,'(a,i10)') &
           'INFO-XML:Optional tag : eigen_mpi%level_lowest  =', solver%eigen_mpi%level_lowest
        endif
      end if
!
      node_tmp => getFirstElementByTagName(node,"level_highest")
      if ( associated(node_tmp) ) then
      value = getChildValue(node_tmp)
        if ( value /= "" ) then
          read(unit=value,fmt=*) solver%eigen_mpi%level_highest
          if (log_unit > 0) write(log_unit,'(a,i10)') &
           'INFO-XML:Optional tag : eigen_mpi%level_highest =', solver%eigen_mpi%level_highest
        endif
      end if
    endif  
!
    ! get <projection> node
    node => getFirstElementByTagName(solver_node,"projection")
    if( .not. associated(node) ) then
!      call XML_error("<projection> not found")
    else
       value = getChildValue(node)
       read(unit=value,fmt=*) solver%projection
       if (log_unit > 0) write(log_unit,'(a,i10)')'INFO-XML:Optional tag detected : solver projection = ', & 
&                                                     solver%projection
    end if

    ! get <dimension> node
    node => getFirstElementByTagName(solver_node,"dimension")
    if( .not. associated(node) ) then
!      call XML_error("<dimension> not found")
    else
       value = getChildValue(node)
       read(unit=value,fmt=*) solver%dimension
       if (log_unit > 0) write(log_unit,'(a,i10)')'INFO-XML:Optional tag detected : solver dimension = ', & 
&                                                     solver%dimension
    end if

    ! get <mode_for_large_memory> node
    node => getFirstElementByTagName(solver_node,"mode_for_large_memory")
    if( .not. associated(node) ) then
!      call XML_error("<mode_for_large_memory> not found")
    else
       value = getChildValue(node)
       read(unit=value,fmt=*) solver%mode_for_large_memory
       if (log_unit > 0) write(log_unit,'(a,i10)')'INFO-XML:Optional tag detected : solver mode_for_large_memory = ', &
&                             solver%mode_for_large_memory
    end if

    ! get <inner_cg_loop> node
    node => getFirstElementByTagName(solver_node,"inner_cg_loop")
    if (associated(node) ) then
       node_tmp => getFirstElementByTagName(node,"max_iteration")
      if( associated(node_tmp) ) then
        value = getChildValue(node_tmp)
       read(unit=value,fmt=*) solver%inner_cg_loop%max_iteration
       if (log_unit > 0) write(log_unit,'(a,i10)') & 
&               'INFO-XML:Optional tag detected : solver inner_cg_loop max_iteration   = ', &
&                             solver%inner_cg_loop%max_iteration
      endif
       node_tmp => getFirstElementByTagName(node,"convergence_eps")
      if( associated(node_tmp) ) then
        value = getChildValue(node_tmp)
       read(unit=value,fmt=*) solver%inner_cg_loop%convergence_eps
       if (log_unit > 0) write(log_unit,'(a,f10.5)') & 
&                  'INFO-XML:Optional tag detected : solver inner_cg_loop convergence_eps = ', &
&                             solver%inner_cg_loop%convergence_eps
      endif
    end if

    ! get <mArnoldi_q> node
    node => getFirstElementByTagName(solver_node,"mArnoldi_q")
    if( .not. associated(node) ) then
!      call XML_error("<mode_for_large_memory> not found")
    else
       value = getChildValue(node)
       read(unit=value,fmt=*) solver%mArnoldi_q
       if (log_unit > 0) write(log_unit,'(a,i10)') & 
&                   'INFO-XML:Optional tag detected : solver mArnoldi_q = ', solver%mArnoldi_q
    end if

    ! get <mode_for_suggest_projection> node
    node => getFirstElementByTagName(solver_node,"mode_for_suggest_projection")
    if( associated(node) ) then
       value = getChildValue(node)
       read(unit=value,fmt=*) solver%mode_for_suggest_projection
       if (log_unit > 0) write(log_unit,'(a,a)') & 
&                   'INFO-XML:Optional tag detected : solver mode_for_suggest_projection = ', & 
&                                                          trim(solver%mode_for_suggest_projection)
    end if

    ! get <flexible_cutoff> node
    node => getFirstElementByTagName(solver_node,"flexible_cutoff")
    if (associated(node) ) then
      value = getAttribute(node,"mode")
      if( value /= "on" ) then
        solver%flexible_cutoff%set = .false. 
      else   
        solver%flexible_cutoff%set = .true. 
        if (log_unit > 0) write(log_unit,'(a,f10.5)') & 
&                 'INFO-XML:Optional tag detected : solver%flexible_cutoff'
        node_tmp => getFirstElementByTagName(node,"flexible_cutoff_01")
        if( associated(node_tmp) ) then
          value = getChildValue(node_tmp)
          read(unit=value,fmt=*) solver%flexible_cutoff%flexible_cutoff_01
          if (log_unit > 0) write(log_unit,'(a,f10.5)') & 
&                   'INFO-XML:Optional tag detected : solver%flexible_cutoff%flexible_cutoff_01  = ', &
&                             solver%flexible_cutoff%flexible_cutoff_01
          node_tmp2 => getFirstElementByTagName(node,"flexible_cutoff_02")
          if( associated(node_tmp2) ) then
            value = getChildValue(node_tmp2)
            read(unit=value,fmt=*) solver%flexible_cutoff%flexible_cutoff_02
            if (log_unit > 0) write(log_unit,'(a,f10.5)') & 
&                       'INFO-XML:Optional tag detected : solver%flexible_cutoff%flexible_cutoff_02  = ', &
&                               solver%flexible_cutoff%flexible_cutoff_02
          endif   
        endif
      endif
    end if

    return
  end subroutine solver_load
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parse <distributed> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine distributed_load( distributed, distributed_node )
    implicit none
    type(distributed_type), intent(out) :: distributed
    type(fnode), pointer :: distributed_node

    type(fnode), pointer :: node
    type(fnode), pointer :: node_wrk1
    character(len=256)   :: value
    integer              :: log_unit

    log_unit=config%calc%distributed%log_unit

    ! get scheme attribute
    value = getAttribute(distributed_node,"mode")
    if( value == "off" ) then
       distributed%set = .false.
       distributed%global_dens_mat = .false.
       distributed%global_ham_mat = .false.
       return
    else
       distributed%set = .true.
       if (log_unit > 0) write(log_unit,'(a)')'INFO-XML:Optional tag detected : distributed mode is ON'
    end if

    ! get <global_density_matrix> node
    node => getFirstElementByTagName(distributed_node,"global_density_matrix")
    if( .not. associated(node) ) then
       distributed%global_dens_mat = .false.
    else
       distributed%global_dens_mat = .true.
       if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional tag detected : global_density_matrix'
    end if

    ! get <global_hamiltonian_matrix> node
    node => getFirstElementByTagName(distributed_node,"global_hamiltonian_matrix")
    if( .not. associated(node) ) then
       distributed%global_ham_mat = .false.
    else
       distributed%global_ham_mat = .true.
       if (log_unit > 0) write(log_unit,'(a)')'INFO-XML:Optional tag detected : global_hamiltonian_matrix'
    end if

    ! get <use_mpi_barrier> node
    node => getFirstElementByTagName(distributed_node,"use_mpi_barrier")
    distributed%use_mpi_barrier = .true.
    if ( associated(node) ) then
       value = getAttribute(node,"mode")
       if ( value == "on" ) then
         distributed%use_mpi_barrier = .true.
         if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional tag detected : use_mpi_barrier is on'
       endif  
       if ( value == "off" ) then
         distributed%use_mpi_barrier = .false.
         if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional tag detected : use_mpi_barrier is off'
       endif  
    end if
!
    ! get <micro_cell_booking> node  ! default values are set in 'distributed_default'
    node => getFirstElementByTagName(distributed_node,"micro_cell_booking")
    if ( associated(node) ) then
      distributed%micro_cell_booking%set = .true.
      if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional XML tag detected : micro_cell_booking'
!
      value = getAttribute(node,"mode")
      if( value /= "" ) then
        distributed%micro_cell_booking%mode = value 
        if (log_unit > 0) write(log_unit,'(a,a16)')'INFO-XML:Optional tag:micro_cell_booking%mode=', & 
&                                                  distributed%micro_cell_booking%mode
      endif  
!
      node_wrk1 => getFirstElementByTagName(node, "cell_number_x")
      if (associated(node_wrk1)) then
        value = getChildValue(node_wrk1)
        read(unit=value,fmt=*) distributed%micro_cell_booking%cell_number(1)
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional :micro_cell_booking%cell_number_x=', & 
&                                                  distributed%micro_cell_booking%cell_number(1)
      endif  
!
      node_wrk1 => getFirstElementByTagName(node, "cell_number_y")
      if (associated(node_wrk1)) then
        value = getChildValue(node_wrk1)
        read(unit=value,fmt=*) distributed%micro_cell_booking%cell_number(2)
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:micro_cell_booking%cell_number_y=', & 
&                                                  distributed%micro_cell_booking%cell_number(2)
      endif  
!
      node_wrk1 => getFirstElementByTagName(node, "cell_number_z")
      if (associated(node_wrk1)) then
        value = getChildValue(node_wrk1)
        read(unit=value,fmt=*) distributed%micro_cell_booking%cell_number(3)
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optionaltag:micro_cell_booking%cell_number_z=', & 
&                                                  distributed%micro_cell_booking%cell_number(3)
      endif  
!
      node_wrk1 => getFirstElementByTagName(node, "search_range_x")
      if (associated(node_wrk1)) then
        value = getChildValue(node_wrk1)
        read(unit=value,fmt=*) distributed%micro_cell_booking%search_range(1)
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:micro_cell_booking%search_range_x=', & 
&                                                  distributed%micro_cell_booking%search_range(1)
      endif  
!
      node_wrk1 => getFirstElementByTagName(node, "search_range_y")
      if (associated(node_wrk1)) then
        value = getChildValue(node_wrk1)
        read(unit=value,fmt=*) distributed%micro_cell_booking%search_range(2)
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:micro_cell_booking%search_range_y=', & 
&                                                  distributed%micro_cell_booking%search_range(2)
      endif  
!
      node_wrk1 => getFirstElementByTagName(node, "search_range_z")
      if (associated(node_wrk1)) then
        value = getChildValue(node_wrk1)
        read(unit=value,fmt=*) distributed%micro_cell_booking%search_range(3)
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:micro_cell_booking%search_range_z=', & 
&                                                  distributed%micro_cell_booking%search_range(3)
      endif  
!
    endif
!
    return
  end subroutine distributed_load
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parse <snapshot> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine snapshot_load( snapshot, snapshot_node )
    implicit none
    type(snapshot_type), intent(out) :: snapshot
    type(fnode), pointer :: snapshot_node
    type(fnode), pointer :: node
    character(len=256)   :: value
    integer              :: log_unit

    log_unit=config%calc%distributed%log_unit
!
    value = getAttribute(snapshot_node,"initial")
    if( value /= "" ) then 
      read(unit=value,fmt=*) snapshot%initial  ! default value is set in snapshot_default 
      if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:snapshot%initial =', snapshot%initial
    endif
!  
    value = getAttribute(snapshot_node,"number")
    if( value /= "" ) then 
      read(unit=value,fmt=*) snapshot%number   ! default value is set in snapshot_default 
      if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:snapshot%final   =', snapshot%number
    endif
!  
    value = getAttribute(snapshot_node,"interval")
    if( value /= "" ) then
      read(unit=value,fmt=*) snapshot%interval ! default value is set in snapshot_default 
      if (log_unit > 0) write(log_unit,'(a,i5)')'INFO-XML:Optional tag:snapshot%interval=', snapshot%interval
    endif
!  
    return
  end subroutine snapshot_load
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parse <cell_change> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine cell_change_load( cell_change, cell_change_node )
    implicit none
    type(cell_change_type), intent(out) :: cell_change
    type(fnode), pointer :: cell_change_node
    type(fnode), pointer :: node
    character(len=256)   :: value
    integer              :: log_unit

    log_unit=config%calc%distributed%log_unit
!
!   if (log_unit > 0) write(log_unit,*) 'cell_change_load start'
!
    value = getAttribute(cell_change_node,"mode")
    if( value == "off" ) then 
      cell_change%set = .false.
      return
    else
      cell_change%set = .true.
      if (log_unit > 0) write(log_unit,*)'INFO-XML:Optional tag:cell_change%set =', cell_change%set
    endif
!
    value = getAttribute(cell_change_node,"scheme")
    if( value /= "" ) then 
      read(unit=value,fmt=*) cell_change%scheme  ! default value is set in cell_change_default 
      if (log_unit > 0) write(log_unit,*) 'INFO-XML:Optional tag:cell_change%scheme =', cell_change%scheme
    endif
!  
    value = getAttribute(cell_change_node,"filename")
    if( value /= "" ) then 
      read(unit=value,fmt=*) cell_change%filename  ! default value is set in cell_change_default 
      if (log_unit > 0) write(log_unit,'(a,a)')'INFO-XML:Optional tag:cell_change%filename =', cell_change%filename
    endif
!
    node  => getFirstElementByTagName(cell_change_node,"max_num_iter")
    if( associated(node) ) then
       value =  getChildValue(node)
       read(unit=value,fmt=*) cell_change%max_num_iter
      if (log_unit > 0) write(log_unit,'(a,i5)') & 
&               'INFO-XML:Optional tag:cell_change%max_num_iter =', cell_change%max_num_iter
    else
      if (log_unit > 0) write(log_unit,'(a)') & 
&               'INFO-XML:No tag for cell_change%max_num_iter. '
       cell_change%max_num_iter = -1 ! dummy value
    end if
!
!   write(*,*) 'cell_change_load end'
!
    return
  end subroutine cell_change_load
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parse <calc_check> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine load_calc_check( calc_check, calc_check_node )
    implicit none
    type(calc_check_type), intent(out) :: calc_check
    type(fnode), pointer :: calc_check_node
    type(fnode), pointer :: node_tmp
    type(fnode), pointer :: node_wrk1
    character(len=256)   :: value, unit
    integer              :: log_unit
    real(8)              :: unit_conv, value_real
    real(8)              :: warning_level_default, abort_level_default
!
    log_unit=config%calc%distributed%log_unit
!
    warning_level_default = 0.5d0 * XML_TO_AU("angstrom") 
    abort_level_default   = 0.1d0 * XML_TO_AU("angstrom")
!
!   calc_check%set = .false.
!   calc_check%short_atom_pair_distance%set = .false.
!   calc_check%short_atom_pair_distance%warning_level = 1.0d10 ! dummmy value
!   calc_check%short_atom_pair_distance%abort_level   = 1.0d10 ! dummmy value
!
    calc_check%set = .true.
    calc_check%short_atom_pair_distance%set = .true.
    calc_check%short_atom_pair_distance%warning_level = warning_level_default
    calc_check%short_atom_pair_distance%abort_level   = abort_level_default ! dummmy value
!
    if (.not. associated(calc_check_node)) then 
      if (log_unit > 0) write(log_unit,'(a,2f20.10)')  & 
&           'INFO:Default setting:short_atom_pair_distance:warning_level [a.u.], [A] =', &
&             calc_check%short_atom_pair_distance%warning_level,  &
&             calc_check%short_atom_pair_distance%warning_level/XML_TO_AU("angstrom")
      if (log_unit > 0) write(log_unit,'(a,2f20.10)')  & 
&           'INFO:Default setting:short_atom_pair_distance:abort_level [a.u.], [A] =', &
&             calc_check%short_atom_pair_distance%abort_level,  &
&             calc_check%short_atom_pair_distance%abort_level/XML_TO_AU("angstrom")

      return
    endif   
!
    value = getAttribute(calc_check_node,"mode")
    if ( value == "off" ) then 
      calc_check%set = .false.
      if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional tag:calc_check : MODE=OFF'
      return
    else
      calc_check%set = .true.
      if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional tag:calc_check'
    endif
!
    node_wrk1 => getFirstElementByTagName(calc_check_node, "short_atom_pair_distance")
    if (associated(node_wrk1)) then
      value = getAttribute(node_wrk1,"mode")
      if ( value == "off" ) then 
        calc_check%short_atom_pair_distance%set = .false.
        return
      else
        calc_check%short_atom_pair_distance%set = .true.
        if (log_unit > 0) write(log_unit,'(a)') 'INFO-XML:Optional tag:calc_check:short_atom_pair_distance'
      endif
!
      node_tmp  => getFirstElementByTagName(node_wrk1,"warning_level")
      if( associated(node_tmp) ) then
         value =  getChildValue(node_tmp)
         read(unit=value,fmt=*) value_real
         unit  = getAttribute(node_tmp,"unit")
         if( unit == "" ) unit = "a.u."
         calc_check%short_atom_pair_distance%warning_level = value_real * XML_TO_AU(unit)
         if (log_unit > 0) write(log_unit,'(a,2f20.10)')  & 
&           'INFO-XML:Optional tag:calc_check:short_atom_pair_distance:warning_level [a.u.], [A] =', &
&             calc_check%short_atom_pair_distance%warning_level,  &
&             calc_check%short_atom_pair_distance%warning_level/XML_TO_AU("angstrom")
      endif   
      node_tmp  => getFirstElementByTagName(node_wrk1,"abort_level")
      if( associated(node_tmp) ) then
         value =  getChildValue(node_tmp)
         read(unit=value,fmt=*) value_real
         unit  = getAttribute(node_tmp,"unit")
         if( unit == "" ) unit = "a.u."
         calc_check%short_atom_pair_distance%abort_level = value_real * XML_TO_AU(unit)
         if (log_unit > 0) write(log_unit,'(a,2f20.10)')  & 
&           'INFO-XML:Optional tag:calc_check:short_atom_pair_distance:abort_level [a.u.],[A]    =', &
&             calc_check%short_atom_pair_distance%abort_level,  &
&             calc_check%short_atom_pair_distance%abort_level/XML_TO_AU("angstrom")
      endif   
!
    endif   
!
!   stop 'stop manually'
!
    return
  end subroutine load_calc_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parse <output> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine output_load( output, output_node )
    implicit none
    type(output_type), intent(out) :: output
    type(fnode), pointer :: output_node

    type(fnode), pointer :: node
    character(len=256)   :: value, value_mode
    integer              :: log_unit

    log_unit=config%calc%distributed%log_unit

    ! get <main> node
    node => getFirstElementByTagName(output_node,"main")
    output%main%set      = .false.
    output%main%filename = "Output.txt"
    if ( associated(node) ) then
      value       = getAttribute(node,"filename") 
      value_mode  = getAttribute(node,"mode") 
      if ( trim(value_mode) /= "default_name" ) then
        if ( trim(value) /= "" ) then
          output%main%filename = trim(value)
          output%main%set      = .true.
          if (log_unit>0) write(log_unit,'(a,a)') & 
            &     'INFO-XML:Optional tag detected : output%main; filename=', &
            &     trim(output%main%filename)
        endif   
      endif
    end if

    ! get <restart> node
    node => getFirstElementByTagName(output_node,"restart")
    output%restart%first_write = .true. 
    output%restart%append_mode = "off"
    if( .not. associated(node) ) then
       output%restart%set = .true.
       output%restart%filename = "restart.xml"
       output%restart%interval = 10
       output%restart%split    = .false.
       output%restart%atom_id_is_added  = .false.
       output%restart%number_of_split_files  = -1
    else
       call file_load( output%restart, node )
    end if

    ! get <wavefunction> node
    node => getFirstElementByTagName(output_node,"wavefunction")
    if( .not. associated(node) ) then
       output%wavefunction%set = .false.
       output%wavefunction%filename = "output_wavefunction.txt"
       output%wavefunction%interval = -1
    else
       call file_load( output%wavefunction, node )
    end if

    ! get <position> node
    node => getFirstElementByTagName(output_node,"position")
    if( .not. associated(node) ) then
       output%position%set = .false.
    else
       call file_load( output%position, node )
    end if

    ! get <velocity> node
    node => getFirstElementByTagName(output_node,"velocity")
    if( .not. associated(node) ) then
       output%velocity%set = .false.
    else
      if (log_unit>0) write(log_unit,'(a)')'INFO-XML:An optional tag is found : config%output%velocity'
      call file_load( output%velocity, node )
    end if

    ! get <atom_charge> node
    node => getFirstElementByTagName(output_node,"atom_charge")
    if( .not. associated(node) ) then
       output%atom_charge%set = .false.
       output%atom_charge%filename = "output_atom_charge.txt"
       output%atom_charge%interval = -1
    else
      call file_load( output%atom_charge, node )
      if (log_unit>0) write(log_unit,'(a)')'INFO-XML:An optional tag is found : config%output%atom_charge'
    end if

    ! get <atom_energy> node
    node => getFirstElementByTagName(output_node,"atom_energy")
    if( .not. associated(node) ) then
      output%atom_energy%set = .false.
      output%atom_energy%filename = "output_atom_energy.txt"
      output%atom_energy%interval = -1
    else
      call file_load( output%atom_energy, node )
      if (log_unit>0) write(log_unit,'(a)') 'INFO-XML:An optional tag is found : config%output%atom_energy'
    end if

    ! get <density_of_states> node
    node => getFirstElementByTagName(output_node,"density_of_states")
    if( .not. associated(node) ) then
       output%density_of_states%set = .false.
    else
       call file_load( output%density_of_states, node )
    end if

    ! get <bond_list> node
    node => getFirstElementByTagName(output_node,"bond_list")
    if( .not. associated(node) ) then
       output%bond_list%set = .false.
       output%bond_list%filename = "output_bond_list.txt"
       output%bond_list%interval = -1
    else
       call file_load( output%bond_list, node )
       if (log_unit>0) write(log_unit,'(a)') 'INFO-XML:A optional tag is found : config%output%bond_list'
    end if

    ! get <matrix_hamiltonian> node
    node => getFirstElementByTagName(output_node,"matrix_hamiltonian")
    if ( .not. associated(node) ) then
      call file_load_matrices( "hamiltonian" ) ! default setting
    else
      call file_load_matrices( "hamiltonian", node )
    endif
   
    ! get <matrix_overlap> node
    node => getFirstElementByTagName(output_node,"matrix_overlap")
    if ( .not. associated(node) ) then
      call file_load_matrices( "overlap" ) ! default setting
    else
      call file_load_matrices( "overlap", node )
    endif

    ! get <eigen_level> node
    node => getFirstElementByTagName(output_node,"eigen_level")
    if ( .not. associated(node) ) then
      call file_load_eigen_level ! default setting
    else
      call file_load_eigen_level( node )
    endif

    return
  end subroutine output_load

  ! parse <file> node
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine file_load( file, file_node )
    implicit none
    type(file_type), intent(out) :: file
    type(fnode), pointer :: file_node

!   type(fnode), pointer :: node
    character(len=256)   :: value
    integer              :: log_unit

    log_unit=config%calc%distributed%log_unit

    ! get directory name attribute
    value=getAttribute(file_node,"dirname")
    file%dirname=adjustL(value)

    ! get filename attribute
    value=getAttribute(file_node,"filename")
    file%filename=adjustL(value)
    
    ! get history attribute
    value=getAttribute(file_node,"history")
    file%history = adjustL(value)
    if( file%history == "" ) file%history = "last"
    
    ! get format attribute
    value=getAttribute(file_node,"format")
    file%format = adjustL(value)
    
    ! get interval attribute
    value = getAttribute(file_node,"interval")
    if( value /= "" ) then
       read(unit=value,fmt=*) file%interval
    end if

    ! get format attribute
    value=getAttribute(file_node,"method")
    file%method = adjustL(value)
    if( file%method == "" ) file%method = "all"

    ! get append_mode attribute
    value=getAttribute(file_node,"append_mode")
    file%append_mode = adjustL(value)
    if( file%append_mode == "" ) file%append_mode = "off"

    ! get cell_info attribute  (For xyz-formatted position file)
    value=getAttribute(file_node,"cell_info")
    file%cell_info = adjustL(value)
    if( file%cell_info == "" ) then
      file%cell_info = "off"
    else
      if (log_unit>0) write(log_unit,*)'INFO-XML:optional attirbute:cell_info=',trim(file%cell_info)
    endif

    ! get split attribute
    value = getAttribute(file_node,"split")
    if ( value == "on" ) then
      file%split= .true.
      if (log_unit>0) write(log_unit,*)'INFO-XML:optional attirbute:split=',file%split
    else
      file%split= .false.
    endif

    ! get number_of_split_files attribute
    value = getAttribute(file_node,"number_of_split_files")
    if( value /= "" ) then
       read(unit=value,fmt=*) file%number_of_split_files
    end if

    ! get atom_id_is_added attribute
    value=getAttribute(file_node,"atom_id_is_added")
    if ( value == "on" ) then
      file%atom_id_is_added = .true.
      if (log_unit>0) write(log_unit,*)'INFO-XML:optional attirbute:atom_id_is_added=',file%atom_id_is_added
    else
      file%atom_id_is_added = .false.
    endif
   
    file%set = .true.

    return
  end subroutine file_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine file_load_matrices(mode, file_node)
    use M_config,   only  : config
    implicit none
    character(len=*), intent(in)   :: mode
    type(fnode), pointer, optional :: file_node
!
    character(len=256)   :: value
    logical :: mode_for_H
    integer              :: log_unit
    integer              :: value_int ! integer value
!
    log_unit=config%calc%distributed%log_unit
!
!   write(*,*)'@@ file_load_matrices'
!
    select case (mode)
      case("hamiltonian") 
        mode_for_H = .true.
      case("overlap") 
        mode_for_H = .false.
      case default
        write(*,*) 'ERROR:file_load_matrices' 
        stop
    end select    
!
    if (mode_for_H) then  
      config%output%matrix_hamiltonian%set      = .false.
      config%output%matrix_hamiltonian%filename = ""
      config%output%matrix_hamiltonian%mode     = "last"
      config%output%matrix_hamiltonian%format   = "MatrixMarket_sym"
      config%output%matrix_hamiltonian%unit     = "a.u."
      config%output%matrix_hamiltonian%interval = 0
    else
      config%output%matrix_overlap%set          = .false.
      config%output%matrix_overlap%filename     = ""
      config%output%matrix_overlap%mode         = "last"
      config%output%matrix_overlap%format       = "MatrixMarket_sym"
      config%output%matrix_overlap%unit         = "a.u."
      config%output%matrix_overlap%interval     = 0
    endif   
!
    if ( .not. present(file_node) ) then
      if (log_unit>0) write(log_unit,*)'.... skipped : file_load_matrices '
      return 
    endif  
!
    if (mode_for_H) then  
      config%output%matrix_hamiltonian%set = .true.
    else
      config%output%matrix_overlap%set     = .true.
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"filename")
    value=trim(adjustL(value))
    if (value /= "") then
      if (mode_for_H) then  
        config%output%matrix_hamiltonian%filename = trim(value)
        if (log_unit>0) write(log_unit,'(a,a)') &
           &   'INFO-XML:An optional tag is found:config%output%matrix_hamiltonian%filename=',trim(value)
      else
        config%output%matrix_overlap%filename     = trim(value)
        if (log_unit>0) write(log_unit,'(a,a)') &
           &   'INFO-XML:An optional tag is found:config%output%matrix_overlap%filename=',trim(value)
      endif   
    else
      write(*,*)'Stop:file_load_matrices:no filename'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"mode")
    value=trim(adjustL(value))
    select case(value)
      case ("off") 
        return 
      case ("", "on","last","default") 
        value="last"
      case ("periodic")  
        value="periodic"
      case default
        write(*,*)'Stop:file_load_matrices:invalid mode'
        stop
    end select   
!
    if (mode_for_H) then  
      config%output%matrix_hamiltonian%mode = trim(value)
      if (log_unit>0) write(log_unit,'(a,a)') &
         & 'INFO-XML:An optional tag is found:config%output%matrix_hamiltonian%mode=',trim(value)
    else
      config%output%matrix_overlap%mode     = trim(value)
      if (log_unit>0) write(log_unit,'(a,a)') &
         & 'INFO-XML:An optional tag is found:config%output%overlap%mode=', trim(value)
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"format")
    value=trim(adjustL(value))
    select case(value)
      case ("default", "MatrixMarket", "MatrixMarket_symmetric", "MatrixMarket_sym") 
        value="MatrixMarket_sym"
      case ("MatrixMarket_general", "MatrixMarket_gen") 
        value="MatrixMarket_gen"
      case default
        write(*,*)'Stop:file_load_matrices:invalid format'
        stop
    end select   
!
    if (mode_for_H) then  
      config%output%matrix_hamiltonian%format = trim(value)
      if (log_unit>0) write(log_unit,'(a,a)') &
        &  'INFO-XML:An optional tag is found:config%output%matrix_hamiltonian%format=', trim(value)
    else
      config%output%matrix_overlap%format     = trim(value)
      if (log_unit>0) write(log_unit,'(a,a)') &
        &  'INFO-XML:An optional tag is found:config%output%matrix_overlap%format=', trim(value)
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"interval")
    value=trim(adjustL(value))
    if (value /= "") then
      read(unit=value,fmt=*) value_int
      if (mode_for_H) then  
        config%output%matrix_hamiltonian%interval = value_int
        if (log_unit>0) write(log_unit,'(a,i10)') &
          &  'INFO-XML:An optional tag is found:config%output%matrix_hamiltonian%interval =', value_int
      else
        config%output%matrix_overlap%interval     = value_int
        if (log_unit>0) write(log_unit,'(a,i10)') &
          &  'INFO-XML:An optional tag is found:config%output%matrix_overlap%interval     =', value_int
      endif   
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"unit")
    value=trim(adjustL(value))
    if (value /= "") then
      select case(value)
        case ("au", "a.u.") 
          value="a.u."
        case default
          write(*,*)'Stop:file_load_matrices:invalid unit:', trim(adjustL(value))
          stop
      end select   
!
      if (mode_for_H) then  
        config%output%matrix_hamiltonian%unit = trim(adjustL(value))
        if (log_unit>0) write(log_unit,'(a,a10)') &
          &  'INFO-XML:An optional tag is found:config%output%matrix_hamiltonian%unit =', trim(adjustL(value))
      else
        config%output%matrix_overlap%unit     = trim(adjustL(value))
        if (log_unit>0) write(log_unit,'(a,a10)') &
          &  'INFO-XML:An optional tag is found:config%output%matrix_overlap%unit     =', trim(adjustL(value))
      endif   

    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  end subroutine file_load_matrices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine file_load_eigen_level(file_node)
    use M_config,   only  : config
    implicit none
    type(fnode), pointer, optional :: file_node
!
    character(len=256)   :: value
    logical :: mode_for_H
    integer              :: log_unit

    log_unit=config%calc%distributed%log_unit
!
!   write(*,*)'@@ file_load_matrices'
!
    config%output%eigen_level%set      = .false.
    config%output%eigen_level%filename = ""
    config%output%eigen_level%mode     = "last"
    config%output%eigen_level%format   = "full"
    config%output%eigen_level%unit     = "a.u."
    config%output%eigen_level%interval = 0
!
    if ( .not. present(file_node) ) then
      if (log_unit > 0) write(log_unit,*)'.... skipped : file_load_eigen_level '
      return 
    endif  
!
    config%output%eigen_level%set = .true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"filename")
    value=trim(adjustL(value))
    if (value /= "") then
      config%output%eigen_level%filename = trim(value)
      if (log_unit > 0) write(log_unit,'(a,a)') & 
        &    'INFO-XML:An optional tag is found:config%output%eigen_level%filename=',trim(value)
    else
      write(*,*)'Stop:file_load_eigen_level:no filename'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"mode")
    value=trim(adjustL(value))
    select case(value)
      case ("off") 
        return 
      case ("", "on","last","default") 
        value="last"
      case ("periodic")  
        value="periodic"
      case default
        write(*,*)'Stop:file_load_eigen_level:invalid mode'
        stop
    end select   
!
    config%output%eigen_level%mode = trim(value)
    if (log_unit > 0) write(log_unit,'(a,a)') & 
      &       'INFO-XML:An optional tag is found:config%output%eigen_level%mode=',trim(value)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"unit")
    value=trim(adjustL(value))
    select case(value)
      case ("a.u.")
      case ("eV") 
      case default
        write(*,*)'Stop:file_load_eigen_level:invalid unit'
        stop
    end select   
!
    config%output%eigen_level%unit = trim(value)
    if (log_unit > 0) write(log_unit,'(a,a)') & 
      &      'INFO-XML:An optional tag is found:config%output%eigen_level%unit=', trim(value)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value=getAttribute(file_node,"format")
    value=trim(adjustL(value))
    select case(value)
      case ("default", "full") 
        value="full"
      case ("level_only") 
      case default
        write(*,*)'Stop:file_eigen_level:invalid format'
        stop
    end select   
!
    config%output%eigen_level%format = trim(value)
    if (log_unit > 0) write(log_unit,'(a,a)') & 
      &      'INFO-XML:An optional tag is found:config%output%eigen_level%format=', trim(value)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    value = getAttribute(file_node,"interval")
    if( value /= "" ) then
       read(unit=value,fmt=*) config%output%eigen_level%interval
       if (log_unit > 0) write(log_unit,'(a,i10)') & 
      &      'INFO-XML:An optional tag is found:config%output%eigen_level%interval=', config%output%eigen_level%interval
    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  end subroutine file_load_eigen_level
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine measure_clock_time_period_loc(time_present, time_previous, time_period)
!      ----> Measure the clock-time period 
!             with the correction for the possible reset of the system clock
!
     implicit none
     real(kind=8), intent(in)  :: time_present, time_previous
     real(kind=8), intent(out) :: time_period
     integer :: count, rate, max
!
     time_period = time_present - time_previous
!
     if (time_period < 0) then
       call system_clock(count, rate, max)
       time_period = time_period + dble(max)/dble(rate)
     endif   
!
   end subroutine measure_clock_time_period_loc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine get_system_clock_time_loc(elapse_time)
!      ----> Get the system_clock time in sec.
!
     implicit none
     real(kind=8), intent(out) :: elapse_time
     integer :: count, rate
!
     call system_clock(count, rate)
     elapse_time=dble(count)/dble(rate)
!
   end subroutine get_system_clock_time_loc

end module elses_xml_config
