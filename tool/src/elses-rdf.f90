!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

program elses_rdf
  use flib_dom
  use M_structure
  use elses_xml_structure
  use elses_xml_misc
  implicit none

  type step_type
     type(structure_type) :: structure
  end type step_type

  type condition_type 
     integer start_step, end_step, interval
     integer ncell(3)
     integer ndiv
     real(8) rmin, rmax, rdiv
     character(8) pair_element1, pair_element2
     character(256) unit
  end type condition_type

  real(8), allocatable, dimension(:) :: rdf
  real(8), allocatable, dimension(:) :: rdf_r

  type(step_type), pointer, dimension(:) :: step
  type( condition_type) condition

  call main

  stop
contains

  subroutine main
    implicit none
    integer i, ierr
    integer npi, npo
    integer  iargc  
    integer, parameter :: num = 6
    character(len=50) format
    character(len=50) flnm1_in, flnm2_in
    character(len=50) flnm_out
    character(len=50) flnmi, flnmo
    character(len=50) tmpflnmi, tmpflnmo
    character(len=6) stepname

    if( iargc() /= 3 ) then
       write(*,'(a)') "# Usage of this utility program"
       write(*,'(a)') "#   ./elses-rdf condition.xml coordinate.xml output-msr.dat"
       stop
    end if

    call getarg(1,flnm1_in)
    call getarg(2,flnm2_in)
    call getarg(3,flnm_out)

    allocate( step(1), stat=ierr )
    if( ierr /= 0 ) stop 'Error : allocate step'

    npi = len_trim(flnm2_in)
    tmpflnmi=''
    do i = 1, len(flnm2_in)-len('.')
       if(flnm2_in(i:i)=='.') then
          npi=i-1
          tmpflnmi = flnm2_in(i:len_trim(flnm2_in))
       end if
    end do
    flnmi =''
    if( tmpflnmi=='') tmpflnmi='.xml'
    flnmi(1:npi) = flnm2_in(1:npi)
    flnmi(npi+num+1:npi+num+len_trim(tmpflnmi)) = trim(tmpflnmi)
    
    npo = len_trim(flnm_out)
    tmpflnmo=''
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
    
    call condition_loadRDF( flnm1_in )

    do i = condition%start_step, condition%end_step, condition%interval
       write(stepname,format) i
       flnmi(npi+1:npi+num) = trim(stepname)
       flnmo(npo+1:npo+num) = trim(stepname)
       call structure_load( step(1)%structure, flnmi )

       write(*,*) 'calculating rdf of step:', i
       call calc_rdf( step(1)%structure )
       call output_rdf(trim(flnmo))

       deallocate(rdf, rdf_r)
    end do

    deallocate(step)
    return
  end subroutine main

  subroutine condition_loadRDF(filename)
    implicit none
    integer i, j, k
    integer, parameter ::  fd = 1
    integer nspecies
    character(len=*) filename
    character(len=256) buffer
    character(len=256) value

    type(fnode), pointer     :: document_node
    type(fnode), pointer     :: condition_node
    type(fnode), pointer     :: dummycell_node
    type(fnode), pointer     :: ncellx_node, ncelly_node, ncellz_node
    type(fnode), pointer     :: rdfrange_node
    type(fnode), pointer     :: step_node, start_node, end_node, interval_node
    type(fnode), pointer     :: rmax_node, rmin_node
    type(fnode), pointer     :: rdiv_node, ndiv_node
    type(fnode), pointer     :: pair_node
    type(fnode), pointer     :: output_node, unit_node

    document_node => parsefile(filename)
    call normalize(document_node)
    
    ! get <condition> node
    condition_node => getFirstElementByTagName(document_node,"condition")
    if( .not. associated(condition_node) ) then
       call XML_error("<condition> not found")
    endif

    ! get <dummycell> node
    dummycell_node => getFirstElementByTagName(condition_node,"dummycell")
    condition%ncell(1:3) = 1
    if( associated(dummycell_node) ) then
       ncellx_node => getFirstElementByTagName(dummycell_node,"ncellx")
       if( associated( ncellx_node ) )then
          value = getChildValue(ncellx_node)
          read(unit=value,fmt=*) condition%ncell(1) 
       end if
       ncelly_node => getFirstElementByTagName(dummycell_node,"ncelly")
       if( associated( ncelly_node ) )then
          value = getChildValue(ncelly_node)
          read(unit=value,fmt=*) condition%ncell(2) 
       end if
       ncellz_node => getFirstElementByTagName(dummycell_node,"ncellz")
       if( associated( ncellz_node ) )then
          value = getChildValue(ncellz_node)
          read(unit=value,fmt=*) condition%ncell(3) 
       end if
    endif

    ! get <step> node
    condition%start_step = 1
    condition%end_step = 1
    condition%interval = 1
    step_node => getFirstElementByTagName(condition_node,"step")
    if( associated(step_node) ) then
       start_node => getFirstElementByTagName(step_node,"start")
       if( associated( start_node ) )then
          value = getChildValue(start_node)
          read(unit=value,fmt=*) condition%start_step
       end if
       end_node => getFirstElementByTagName(step_node,"end")
       if( associated( end_node ) )then
          value = getChildValue(end_node)
          read(unit=value,fmt=*) condition%end_step
       end if
       interval_node => getFirstElementByTagName(step_node,"interval")
       if( associated( interval_node ) )then
          value = getChildValue(interval_node)
          read(unit=value,fmt=*) condition%interval
       end if
    end if

    ! get <output> node
    condition%rmin = 0.d0
    condition%rmax = -1.d0
    condition%rdiv = 0.1d0
    condition%ndiv = int(condition%rmax/condition%rdiv)+1
    condition%unit = "a.u."

    output_node => getFirstElementByTagName(condition_node,"output")
    if( associated(output_node) ) then

       rmax_node => getFirstElementByTagName(output_node,"rmax")
       if( associated( rmax_node ) )then
          value = getChildValue(rmax_node)
          read(unit=value,fmt=*) condition%rmax
       end if

       rmin_node => getFirstElementByTagName(output_node,"rmin")
       if( associated( rmin_node ) )then
          value = getChildValue(rmin_node)
          read(unit=value,fmt=*) condition%rmin
       end if

       rdiv_node => getFirstElementByTagName(output_node,"rdiv")
       ndiv_node => getFirstElementByTagName(output_node,"ndiv")
       if( associated(rdiv_node) .and. associated(ndiv_node) ) then
          call XML_error("<rdiv> and <ndiv> cannot be set simultaneously")
       elseif( associated( rdiv_node ) )then
          value = getChildValue(rdiv_node)
          read(unit=value,fmt=*) condition%rdiv
          condition%ndiv = int(condition%rmax/condition%rdiv)+1
       elseif( associated( ndiv_node ) )then
          value = getChildValue(ndiv_node)
          read(unit=value,fmt=*) condition%ndiv
       end if

       condition%pair_element1 = 'all'
       condition%pair_element2 = 'all'
       pair_node => getFirstElementByTagName(output_node,"pair")
       if( associated( pair_node ) )then
          value = getChildValue(pair_node)
          j=1; k=1
          condition%pair_element1 = ''
          do i = 1, len_trim(value)
             if(value(i:i)=='-') then
                condition%pair_element2 = ''
                j = 1
                k = k+1
                cycle
             end if
             if( value(i:i) /= ' ' ) then
                select case(k)
                case(1)
                   condition%pair_element1(j:j) = value(i:i)
                case(2)
                   condition%pair_element2(j:j) = value(i:i)
                case default
                   stop 'setting of <pair> is wrong!!'
                end select
                j = j + 1
             end if
          end do
       end if
       unit_node => getFirstElementByTagName(output_node,"unit")
       if( associated( unit_node ) )then
          value = getChildValue(unit_node)
          if( value == "" ) then
             value = "a.u."
          endif
          read(unit=value,fmt=*) condition%unit
       end if
       
    endif

    call destroyNode(document_node)

    return
  end subroutine condition_loadRDF

  subroutine calc_rdf( str )
    implicit none

    integer i, j, index, ierr
    integer ii, jj, kk
    integer num_atom
    real(8) dist, dist_min
    real(8) r12(3)
    real(8) rshift(3)
    real(8) unitcell(3,3)
    real(8), allocatable, dimension(:,:) ::  r
    character(len=8), allocatable, dimension(:) ::  name
    type(structure_type) :: str 

    num_atom = str%natom

    allocate( r(num_atom,3), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocate r in sbrt calc_rdf'
       stop
    end if

    allocate( name(num_atom), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocate name in sbrt calc_rdf'
       stop
    end if

    unitcell(1,1:3) = str%unitcell%vectorA*XML_AU_TO(condition%unit)
    unitcell(2,1:3) = str%unitcell%vectorB*XML_AU_TO(condition%unit)
    unitcell(3,1:3) = str%unitcell%vectorC*XML_AU_TO(condition%unit)

    do i = 1, num_atom
       r(i,1:3) = str%vatom(i)%position(1:3)*XML_AU_TO(condition%unit)
       name(i) = str%vatom(i)%name
    end do
    
    if( condition%rmax < 0.d0 ) then
       condition%rmax = 0.5d0 * sqrt((sum(unitcell(1:3,1))**2.d0) &
            + (sum(unitcell(1:3,2))**2.d0) + (sum(unitcell(1:3,3))**2.d0) )
       condition%ndiv = int(condition%rmax/condition%rdiv)+1
    end if

    allocate( rdf(condition%ndiv), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocate rdf in sbrt calc_rdf'
       stop
    end if

    allocate( rdf_r(condition%ndiv), stat=ierr )
    if( ierr /= 0 ) then
       write(*,*) 'Error : allocate rdf_r in sbrt calc_rdf'
       stop
    end if

    do i = 1, condition%ndiv
       rdf_r(i) = condition%rdiv * ( dble(i-1) + 0.5d0 )
    end do

    rdf(:) = 0.d0
    do i = 1, num_atom
       do j = 1, num_atom
          if( ((condition%pair_element1 /= name(i)) .or. (condition%pair_element2 /= name(j))) &
               .and.((condition%pair_element1 /= name(j)) .or. (condition%pair_element2 /= name(i))) &
               .and. ((condition%pair_element1 /= 'all') .or. (condition%pair_element2 /= 'all'))) cycle

          dist_min = -1.d0
          do ii = -condition%ncell(1), condition%ncell(1)
             do jj = -condition%ncell(2), condition%ncell(2)
                do kk = -condition%ncell(3), condition%ncell(3)
                   rshift(1:3) = dble(ii)*unitcell(1,1:3) &
                        +dble(jj)*unitcell(2,1:3)+dble(kk)*unitcell(3,1:3)
                   r12(1:3) = r(i,:) - ( r(j,:) + rshift )
                   dist = sqrt(dot_product(r12(:),r12(:)))
                   if( dist_min < 0.d0 .or. dist < dist_min ) dist_min = dist
                end do
             end do
          end do

          dist = dist_min
          if( dist > condition%rmax ) cycle
          index = int(dist/condition%rdiv)

          rdf(index) = rdf(index) + 1.d0

       end do
    end do

    deallocate( r )
    deallocate( name )

    return

  end subroutine calc_rdf

  subroutine output_rdf(filename)
    implicit none
    integer i
    character(len=*) filename
    open(1,file=filename)

    do i = 1, condition%ndiv
       write(1,'(2f20.12)') rdf_r(i), rdf(i) 
    end do

    close(1)
    return
  end subroutine output_rdf

end program elses_rdf
