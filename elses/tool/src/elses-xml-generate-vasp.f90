!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

program elses_xml_generate_vasp
  use flib_dom
  use M_structure
  use elses_xml_structure
  use elses_xml_misc
  implicit none

  integer total_num_atom

  type :: generate_type
     integer              :: num_species
     integer, pointer, dimension(:) :: num_atom
     real(8)              :: scale
     character(len=64)    :: comment
     character(len=64)    :: name
     character(len=64)    :: mode
     type(unitcell_type)  :: unitcell
     type(heatbath_type)  :: heatbath
  end type generate_type

  type(generate_type) :: generate
  type(atom_type), pointer, dimension(:) :: atom
  type(element_type), pointer, dimension(:) :: element

  call main

contains

  subroutine main
    implicit none
    integer       :: iargc  
    character(len=50) flnm_poscar_in
    character(len=50) flnm_condition_in
    character(len=50) flnm_xml_out

    if( iargc() /= 3 ) then
       write(*,'(a)') "# Usage of this utility program"
       write(*,'(a)') "#   elses-xml-generate input_generate.xml POSCAR output_structures.xml"
       stop
    end if
    
    call getarg(1,flnm_condition_in)
    call getarg(2,flnm_poscar_in)
    call getarg(3,flnm_xml_out)

    call condition_loadVASP(flnm_condition_in)
    call structure_loadVASP(flnm_poscar_in)
    call generate_output(flnm_xml_out)
    
    call deallocate_all
    return
  end subroutine main

  subroutine condition_loadVASP(filename)
    implicit none
    integer, parameter ::  fd = 1
    integer nspecies
    character(len=*) filename
    character(len=256) buffer
    character(len=256) value

    type(fnode), pointer     :: document_node
    type(fnode), pointer     :: element_node
    type(fnode), pointer     :: name_node
    type(fnode), pointer     :: number_node
    type(fnode), pointer     :: heatbath_node

    document_node => parsefile(filename)
    call normalize(document_node)
    
    ! get <element> node
    element_node => getFirstElementByTagName(document_node,"element")
    if( .not. associated(element_node) ) then
       call XML_error("<element> not found")
    endif

    number_node => getFirstElementByTagName(element_node,"number")
    if( .not. associated(number_node) ) then
       call XML_error("<number> not found")
    else
       value = getChildValue(number_node)
       read(unit=value,fmt=*) generate%num_species
       nspecies = generate%num_species
       allocate( generate%num_atom(generate%num_species) )
       allocate( element(generate%num_species) )
    endif

    name_node => getFirstElementByTagName(element_node,"name")
    if( .not. associated(name_node) ) then
       call XML_error("<name> not found")
    else
       value = getChildValue(name_node)
       read(unit=value,fmt=*) element(1:nspecies)%name
    endif
    
    ! get <heatbath> node
    heatbath_node => getFirstElementByTagName(document_node,"heatbath")
    if( .not. associated(heatbath_node) ) then
       call XML_error("<heatbath> not found")
       call heatbath_default( generate%heatbath )
    else
       call heatbath_default( generate%heatbath )
       call heatbath_load( generate%heatbath, heatbath_node )
    endif

    return
  end subroutine condition_loadVASP

  subroutine structure_loadVASP(filename)
    implicit none
    integer i, line
    integer, parameter :: fd = 1 
    character(1)  ctmp(3)
    character(265) buffer
    character(len=*) filename

    open(fd,file=filename)
    write(*,*) filename
    line = 1
    do
       read(fd,'(a)',end=1000) buffer
       if(buffer == "" ) cycle
       if( line == 1 ) then
          read(buffer,*) generate%comment
       elseif( line == 2 ) then
          read(buffer,*) generate%scale
       elseif( line == 3 ) then
          read(buffer,*) generate%unitcell%vectorA(1:3)
          generate%unitcell%vectorA(1:3) = generate%unitcell%vectorA(1:3) * generate%scale * XML_TO_AU("angstrom") 
       elseif( line == 4 ) then
          read(buffer,*) generate%unitcell%vectorB(1:3)
          generate%unitcell%vectorB(1:3) = generate%unitcell%vectorB(1:3) * generate%scale * XML_TO_AU("angstrom")
       elseif( line == 5 ) then
          read(buffer,*) generate%unitcell%vectorC(1:3)
          generate%unitcell%vectorC(1:3) = generate%unitcell%vectorC(1:3) * generate%scale * XML_TO_AU("angstrom")
       elseif( line == 6 ) then
          read(buffer,*) generate%num_atom(1:generate%num_species)
          total_num_atom = sum(generate%num_atom(1:generate%num_species))

          allocate( atom( total_num_atom ) )

       elseif( line == 8 ) then
          read(buffer,*) generate%mode

       elseif( line > 8 ) then
          read(buffer,*) atom(line-8)%position(1:3), ctmp(1:3)
          ! temporaty
          if( line-8 <= generate%num_atom(1) ) then
             atom(line-8)%element => element(1)
          else
             do i = 1, generate%num_species-1
                if( line-8 > generate%num_atom(i) .and. line-8 <= sum(generate%num_atom(1:i+1)) ) then 
                   atom(line-8)%element => element(i+1)
                end if
             end do
          end if

          do i = 1, 3
             select case(ctmp(i))
             case('T','t')
                atom(line-8)%motion = 'free'
             case('F','f')
                atom(line-8)%motion = 'fix'
             end select
          end do
       end if
       line = line + 1
    end do
    
1000 continue

    return
  end subroutine structure_loadVASP

  subroutine generate_output(filename)
    implicit none
    integer i
    integer, parameter :: fd = 1
    real(8) a, b, c
    character(len=*), intent(in) :: filename

    open(fd,file=filename)

    write(fd,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(fd,'(a,a,a)') '<structure name="', trim(generate%comment), '" mdstep="0">'

    write(fd,*) ""

    call unitcell_save( fd, generate%unitcell )
    call heatbath_save( fd, generate%heatbath )

    do i=1, total_num_atom

       write(fd,*) '<atom element="', trim(atom(i)%element%name), '" ', &
!            'class="', trim(class), '" ', &
            'motion="', trim(atom(i)%motion), '">'
       select case(generate%mode)
       case('Direct','direct','DIRECT')
          a = atom(i)%position(1)
          b = atom(i)%position(2)
          c = atom(i)%position(3)
          write(fd,'(a28,3e23.15,a11)') '  <position unit="internal">', a, b, c, '</position>'
       case('cart','Cartesian','cartesian','CARTESIAN')
          a = atom(i)%position(1) * XML_TO_AU("angstrom")
          b = atom(i)%position(2) * XML_TO_AU("angstrom")
          c = atom(i)%position(3) * XML_TO_AU("angstrom")
          write(fd,'(a28,3e23.15,a11)') '  <position unit="a.u.">', a, b, c, '</position>'
       end select

       if( atom(i)%velocity_set ) then
          write(fd,'(a24,3e23.15,a11)') '  <velocity unit="a.u.">', atom(i)%velocity(1:3), '</velocity>'
       end if
       if( atom(i)%force_set ) then
          write(fd,'(a24,3e23.15,a11)') '  <force    unit="a.u.">', atom(i)%force(1:3), '</force>'
       end if
       
       write(fd,*) '</atom>'
       write(fd,*) ""
       
    end do
    
    write(fd,*) '</structure>'
    
    close(fd)
    
    return
  end subroutine generate_output

  subroutine deallocate_all
    implicit none

    deallocate( generate%num_atom )
    deallocate( element )
    deallocate( atom )

    return
  end subroutine deallocate_all

end program elses_xml_generate_vasp

