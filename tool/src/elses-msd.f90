!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

program elses_msd
  use flib_dom
  use M_structure
  use elses_xml_structure
  use elses_xml_misc
  implicit none

  type step_type
     type(structure_type) :: structure
  end type step_type

  type state_type
     integer step_start
     integer step_end
     integer interval
  end type state_type

  type condition_type 
     type( state_type ) initial
     type( state_type ) final
     integer nsample
     real(8) time_step  ! default in fsec
     character(256) time_unit
     character(256) length_unit
  end type condition_type

  real(8), allocatable, dimension(:) :: msd

  type(step_type), pointer, dimension(:) :: step
  type(condition_type) condition

  call main

  stop

contains

  subroutine main
    implicit none
    integer i, j, ierr
    integer, parameter :: num = 6
    integer npi, npo
    integer iargc  
    character(len=50) format
    character(len=50) flnm1_in, flnm2_in
    character(len=50) flnm_out
    character(len=50) flnmi1, flnmi2, flnmo
    character(len=50) tmpflnmi, tmpflnmo
    character(len=6) stepname

    if( iargc() /= 3 ) then
       write(*,'(a)') "# Usage of this utility program"
       write(*,'(a)') "#   ./elses-msd condition.xml coordinate.xml output_msd.dat"
       stop
    end if

    call getarg(1,flnm1_in)
    call getarg(2,flnm2_in)
    call getarg(3,flnm_out)

    allocate( step(2), stat=ierr )
    if( ierr /= 0 ) stop 'Error : allocate step'

    npi = len_trim(flnm2_in)
    tmpflnmi=''
    do i = 1, len(flnm2_in)-len('.')
       if(flnm2_in(i:i)=='.') then
          npi=i-1
          tmpflnmi = flnm2_in(i:len_trim(flnm2_in))
       end if
    end do
    flnmi1 =''
    flnmi1(1:npi) = flnm2_in(1:npi)
    if( tmpflnmi=='') tmpflnmi='.xml'
    flnmi1(npi+num+1:npi+num+len_trim(tmpflnmi)) = trim(tmpflnmi)
    flnmi2 = flnmi1

    npo = len_trim(flnm_out)
    tmpflnmo =''
    do i = 1, len(flnm_out)-len('.')
       if(flnm_out(i:i)=='.') then
          npo=i-1
          tmpflnmo = flnm_out(i:len_trim(flnm_out))
       end if
    end do
    flnmo =''
    if( tmpflnmo=='') tmpflnmo='.dat'
    flnmo(1:npo) = flnm_out(1:npo)
    flnmo(npo+num+1:npo+num+len_trim(tmpflnmo)) = trim(tmpflnmo)

    write(format,'(a2,i1,a1,i1,a1)')  '(i',num,'.',num,')'

    call condition_loadMSD( flnm1_in )
    condition%nsample = (condition%final%step_end-condition%final%step_start)/condition%final%interval + 1
    do i = condition%initial%step_start, condition%initial%step_end, condition%initial%interval
       write(stepname,format) i
       flnmi1(npi+1:npi+num) = trim(stepname)
       flnmo(npo+1:npo+num) = trim(stepname)
       call structure_load( step(1)%structure, flnmi1 )

       allocate( msd(condition%nsample), stat=ierr )
       if( ierr /= 0 ) then
          stop 'Error : allocate msd in sbrt main'
       end if

       do j = condition%final%step_start, condition%final%step_end, condition%final%interval
          write(stepname,format) j
          flnmi2(npi+1:npi+num) = trim(stepname)

          call structure_load( step(2)%structure, flnmi2 )

          call calc_msd(step(1)%structure, step(2)%structure, msd(j) )

       end do

       call output_msd(trim(flnmo))

       deallocate( msd )
    end do

    deallocate(step)

    return
  end subroutine main

  subroutine condition_loadMSD( filename )
    implicit none
    character(len=*) filename
    character(len=256) value, unit
    type(fnode), pointer     :: document_node
    type(fnode), pointer     :: condition_node
    type(fnode), pointer     :: step_node, time_step_node
    type(fnode), pointer     :: initial_node, final_node
    type(fnode), pointer     :: istart_node, iend_node, iinterval_node
    type(fnode), pointer     :: fstart_node, fend_node, finterval_node
    type(fnode), pointer     :: length_unit_node

    document_node => parsefile(filename)
    call normalize(document_node)

    condition_node => getFirstElementByTagName(document_node,"condition")
    if( .not. associated(condition_node) ) then
       call XML_error("<condition> not found")
    endif

    ! get <step> node
    condition%initial%step_start = 1
    condition%initial%step_end = 1
    condition%initial%interval = 1
    condition%final%step_start = 1
    condition%final%step_end = 1
    condition%final%interval = 1

    step_node => getFirstElementByTagName(condition_node,"step")
    if( associated(step_node) ) then
       time_step_node => getFirstElementByTagName(step_node,"timestep")
       if( associated(time_step_node) ) then
          unit  = getAttribute(time_step_node,"unit")
          if( unit == "" ) then
             unit = "fsec"
          endif
          condition%time_unit = unit
          value = getChildValue(time_step_node)
          read(unit=value,fmt=*) condition%time_step
       end if

       initial_node => getFirstElementByTagName(step_node,"initial")
       if( associated(initial_node) ) then
          istart_node => getFirstElementByTagName(initial_node,"start")
          if( associated( istart_node ) )then
             value = getChildValue(istart_node)
             read(unit=value,fmt=*) condition%initial%step_start
          end if
          iend_node => getFirstElementByTagName(initial_node,"end")
          if( associated( iend_node ) )then
             value = getChildValue(iend_node)
             read(unit=value,fmt=*) condition%initial%step_end
          end if
          iinterval_node => getFirstElementByTagName(initial_node,"interval")
          if( associated( iinterval_node ) )then
             value = getChildValue(iinterval_node)
             read(unit=value,fmt=*) condition%initial%interval
          end if
       end if

       final_node => getFirstElementByTagName(step_node,"final")
       if( associated(final_node) ) then
          fstart_node => getFirstElementByTagName(final_node,"start")
          if( associated( fstart_node ) )then
             value = getChildValue(fstart_node)
             read(unit=value,fmt=*) condition%final%step_start
          end if
          fend_node => getFirstElementByTagName(final_node,"end")
          if( associated( fend_node ) )then
             value = getChildValue(fend_node)
             read(unit=value,fmt=*) condition%final%step_end
          end if
          finterval_node => getFirstElementByTagName(final_node,"interval")
          if( associated( finterval_node ) )then
             value = getChildValue(finterval_node)
             read(unit=value,fmt=*) condition%final%interval
          end if
       end if
    end if

    length_unit_node => getFirstElementByTagName(condition_node,"length_unit")
    if( associated(length_unit_node) ) then
       value = getChildValue(length_unit_node)
       read(unit=value,fmt=*) condition%length_unit
    end if

    call destroyNode(document_node)
    return
  end subroutine condition_loadMSD

  subroutine calc_msd( str1, str2, distsum2 )
    implicit none
    integer i, index, ierr
    integer ii, jj, kk
    integer num_atom
    real(8), allocatable, dimension(:,:) ::  r1, r2
    real(8) unitcell(3,3)
    real(8) dist2, distsum2
    real(8) r12(3), rshift(3)
    character(len=8), allocatable, dimension(:) ::  name
    type(structure_type) :: str1 
    type(structure_type) :: str2 

    if( str1%natom /= str2%natom ) then
       write(*,*) 'Error : numbers of atoms are different between str1/str2'
       write(*,*) 'str1%natom=', str1%natom
       write(*,*) 'str2%natom=', str2%natom
       stop
    end if
    num_atom = str1%natom

    allocate( r1(num_atom,3), r2(num_atom,3), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocate r1/r2 in sbrt calc_rdf'
       stop
    end if
    
    allocate( name(num_atom), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocate name in sbrt calc_rdf'
       stop
    end if

    unitcell(1,1:3) = str2%unitcell%vectorA*XML_AU_TO(condition%length_unit)
    unitcell(2,1:3) = str2%unitcell%vectorB*XML_AU_TO(condition%length_unit)
    unitcell(3,1:3) = str2%unitcell%vectorC*XML_AU_TO(condition%length_unit)

    do i = 1, num_atom
       r1(i,1:3) = str1%vatom(i)%position(1:3)*XML_AU_TO(condition%length_unit)
       r2(i,1:3) = str2%vatom(i)%position(1:3)*XML_AU_TO(condition%length_unit)
       name(i) = str1%vatom(i)%name
    end do

    distsum2 = 0.d0
    do i = 1, num_atom
       rshift(1:3) = dble(ii)*unitcell(1,1:3) &
            +dble(jj)*unitcell(2,1:3)+dble(kk)*unitcell(3,1:3)
       r12(1:3) = r1(i,:) - r2(i,:)
       dist2 = dot_product(r12(:),r12(:))
       distsum2 = distsum2 + dist2
    end do

    return
  end subroutine calc_msd

  subroutine output_msd(filename)
    implicit none
    integer i
    character(len=*) filename

    open(1,file=filename)

    write(1,*) ' time(', trim(condition%time_unit), ')    ',  'msd' 
    do i = 1, condition%nsample
       write(1,'(f10.6,f20.12)') condition%time_step*dble(i), msd(i) 
    end do

    close(1)
    return
  end subroutine output_msd

end program elses_msd
