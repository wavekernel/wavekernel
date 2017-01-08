!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

program elses_xml_generate
! use M_lib_phys_const, only : angst
  use flib_dom
  use M_structure
! use elses_xml_misc
  use elses_xml_structure
  implicit none
  real(8), parameter :: angst=0.529177d0
  type :: cluster_type
     type(structure_type)  :: structure
     character(len=16)     :: class
!    logical               :: center_set
!    real(8), dimension(3) :: center
!    logical               :: motion_set
!    real(8), dimension(3) :: motion
!    logical               :: translation_set
!    real(8), dimension(3) :: translation
!    logical               :: boltzmann_set
!    real(8)               :: boltzmann
  end type cluster_type

  type :: generate_type
     character(len=64)    :: name
     type(unitcell_type)  :: unitcell
     type(heatbath_type)  :: heatbath
     integer              :: ncluster
!    type(cluster_type),pointer,dimension(:) :: vcluster
     logical              :: atom_id_is_added_flag ! set in generate_load
     integer, allocatable :: group_id_booking(:)   ! set in read_group_file
     integer              :: n_group               ! set in read_group_file
     integer              :: max_num_group_atoms   ! set in read_group_file
     logical              :: boundary_condition_set ! set in generate_load
     character(len=64)    :: boundary_condition(3)  ! set in generate_load 
!
  end type generate_type
!
  type :: atom_info_type
     real(8)              :: position_angst(3)
     character(len=16)    :: element_name
     character(len=16)    :: motion
  end type atom_info_type
!
  type :: structure_info_type
     type(atom_info_type), allocatable :: atom_info(:)
     integer              :: total_number_of_atoms
     real(8)              :: cell_size_angst(3)
  end type structure_info_type
!
  type(structure_info_type) :: structure_info
!
  call main

contains

  ! the main function of this program
  subroutine main
    implicit none
    type(generate_type) :: generate
    character(64) :: filename_in
    character(64) :: filename_out
    character(64) :: filename_xyz
    character(64) :: filename_fixed_atom_list
    integer       :: iargc  
    character(64) :: argv3
    integer       :: ierr
    integer       :: n_split
!    
    call get_info_from_command_argument(filename_xyz, filename_in, filename_out, & 
&            filename_fixed_atom_list, n_split, ierr)
!
    if ( ierr /= 0 ) then
      write(*,'(a)') "# Usage of this utility program (1) for single output file"
      write(*,'(a)') "#   elses-xml-generate input_generate.xml output_structures.xml"
      write(*,'(a)') "# Usage of this utility program (2) for splitted output file"
      write(*,'(a)') "#   elses-xml-generate input_generate.xml output_structures.xml -split=XX"
      write(*,'(a)') "#   where XX is the number of the splitted files" 
      write(*,'(a)') "# Usage of this utility program (3) "
      write(*,'(a)') "#   elses-xml-generate -xyz_file=input_xyz.xyz -setting_file=input_generate.xml "
      write(*,'(a)') "#                        -output_file=output_structures.xml -split=XX"
      write(*,'(a)') "#   where XX is the number of the splitted files" 
      stop
    end if
!    
    write(*,'(a)')'@@@ elses-xml-generate'
!
    call generate_load( generate, filename_in, filename_xyz )
!
    call set_fixed_atom ( filename_fixed_atom_list )
!
    call generate_output_split( generate, filename_out, n_split )
!
    write(*,'(a)')'.... elses-xml-generate ends without error'
!
  end subroutine main


  ! parse <generate> node
  subroutine generate_load( generate, filename, filename_xyz )
    implicit none
    type(generate_type), intent(out) :: generate
    character(len=*),intent(in)   :: filename
    character(len=*),intent(in)   :: filename_xyz
    character(len=256) :: cwd ! current working dir name for backup
    character(len=256) :: nwd ! current working dir name of config.xml

    type(fnode), pointer     :: document_node
    type(fnode), pointer     :: generate_node
    type(fnode), pointer     :: unitcell_node
    type(fnode), pointer     :: heatbath_node
    type(fnodeList), pointer :: vcluster_node
    type(fnode), pointer     :: cluster_node
    type(fnode), pointer     :: group_node
    type(fnode), pointer     :: atom_id_node
    type(fnode), pointer     :: boundary_node
!   type(fnode), pointer     :: node

    integer                  :: i, j, ierr
!   character(len=256)       :: value, unit
    character(len=256)       :: value
    character(len=256)       :: filename_wrk
    character(len=256)       :: name_of_xyz_file
    character(len=256)       :: chara_wrk
    logical        :: ex
    logical        :: cell_info_in_xyz_file
    real(8)        :: cell_info(3)
!
    cell_info(1:3)   = -1.0d0 ! dummy value
    name_of_xyz_file = ''   ! dummy value
!
!   write(*,*)'INFO:subroutine generate_load'
!
    inquire( file=filename, exist=ex )
    if( .not. ex ) then
       write(*,'(a,a)') '# Error!: generate_load : can not open file ', trim(filename)
       stop
    end if

    document_node => parsefile(filename)
    call normalize(document_node)

    call getcwd( cwd )
    nwd = filename(1:scan(filename,"/\\",back=.true.))
    call chdir(nwd)

    ! get <generate> node
    generate_node => getFirstElementByTagName(document_node,"generate")
    if( .not. associated(generate_node) ) then
      write(*,*)"ERROR:<generate> not found"
      stop 
    endif

    ! get name attribute
    value = getAttribute(generate_node,"name")
    if( value == "" ) then
      write(*,*)"ERROR:<generate name> not found"
      stop
    endif
    generate%name = value
!
    if ( len_trim(filename_xyz) /= 0) then 
       name_of_xyz_file = trim(filename_xyz)
       write(*,'(a,a)')'INFO:XYZ file name=',trim(name_of_xyz_file)
    endif
!
    ! get <cluster> nodes for multiple tags
    !    (only for the compatibility to the older code
    vcluster_node => getElementsByTagName(generate_node,"cluster")
    generate%ncluster = getLength(vcluster_node)
    if( generate%ncluster == 0 ) then
      write(*,*)"INFO:<cluster> not found in the setting XML file"
    endif
    if ( generate%ncluster > 1  ) then
      write(*,*)'ERROR:ncluster =', generate%ncluster 
    endif   
!
    ! get <cluster> nodes for single tag
    cluster_node => getFirstElementByTagName(generate_node,"cluster")
    if ( len_trim(filename_xyz) == 0) then 
      if ( associated(cluster_node) ) then
       name_of_xyz_file = getAttribute(cluster_node,"structure")
       write(*,'(a,a)')'INFO:XYZ file name=',trim(name_of_xyz_file)
      endif
    endif
    cell_info_in_xyz_file = .false. ! dummy value
    call get_info_from_xyz_file(name_of_xyz_file, cell_info_in_xyz_file)

    generate%unitcell%vectorA(:)=0.0d0   ! dummy value
    generate%unitcell%vectorB(:)=0.0d0   ! dummy value
    generate%unitcell%vectorC(:)=0.0d0   ! dummy value

    generate%unitcell%vectorA(1)=-1.0d0  ! dummy value
    generate%unitcell%vectorB(2)=-1.0d0  ! dummy value
    generate%unitcell%vectorC(3)=-1.0d0  ! dummy value

    ! get <unitcell> node
    unitcell_node => getFirstElementByTagName(generate_node,"unitcell")
    if ( associated(unitcell_node) ) then
      call unitcell_load( generate%unitcell, unitcell_node )
      write(*,'(a)')'INFO:<unitcell> tag was found in the XML file and will be used for the unit cell info.' 
      structure_info%cell_size_angst(1)=generate%unitcell%vectorA(1)*angst
      structure_info%cell_size_angst(2)=generate%unitcell%vectorB(2)*angst
      structure_info%cell_size_angst(3)=generate%unitcell%vectorC(3)*angst
    else  
      write(*,'(a)')'INFO:No <unitcell> tag was found in the XML file' 
      write(*,'(a)')'INFO:The data in the XYZ file will be used for the unit cell info, if exist.' 
    endif
!
    ! get <boundary> node
    boundary_node => getFirstElementByTagName(generate_node,"boundary")
    generate%boundary_condition_set = .false. ! dummy value
    generate%boundary_condition(1:3)=''       ! dummy value
    if(associated(boundary_node) ) then
      generate%boundary_condition_set = .true.
      value=getAttribute(boundary_node,"mode")
      if (trim(value) == 'off') then
        generate%boundary_condition_set = .false. 
      else
        do j=1,3
          if (j==1) chara_wrk=getAttribute(boundary_node,"x")
          if (j==2) chara_wrk=getAttribute(boundary_node,"y")
          if (j==3) chara_wrk=getAttribute(boundary_node,"z")
!         write(*,*)'INFO:j, chara_wrk=',j, trim(chara_wrk)
          select case(trim(chara_wrk))
            case('nonperiodic')
              generate%boundary_condition(j)='nonperiodic' ! non-periodic BC
            case('periodic')
              generate%boundary_condition(j)='periodic'    ! periodic BC
            case default
              write(*,*)'ERROR:<boundary> tag'
              stop
          end select   
        enddo
      endif  
    endif   
!
    ! get <heatbath> node
    heatbath_node => getFirstElementByTagName(generate_node,"heatbath")
    if( .not. associated(heatbath_node) ) then
       generate%heatbath%massperatom = 25.0
       generate%heatbath%position = 0.0
       generate%heatbath%velocity = 0.0
       write(*,*)'INFO:No value is specified in the XML file  : heatbath%massperatom '
       write(*,*)'INFO:The default value is set               : heatbath%massperatom =', generate%heatbath%massperatom
    else
       call heatbath_load( generate%heatbath, heatbath_node )
       write(*,*)'INFO:The value is specified in the XML file : heatbath%massperatom =', generate%heatbath%massperatom
    endif
!
!   allocate( generate%vcluster(generate%ncluster) )
!   do i=1, generate%ncluster
!     if (cell_info_in_xyz_file) then 
!       call cluster_load( generate%vcluster(i), item(vcluster_node,i-1), cell_info )
!       write(*,*)'INFO:cell info in XYZ file =',cell_info(1:3)
!       ierr=0
!       if (cell_info(1) < -1.0d-10) ierr=1
!       if (cell_info(2) < -1.0d-10) ierr=1
!       if (cell_info(3) < -1.0d-10) ierr=1
!       if (ierr == 1) then
!         write(*,*)'ERROR:No unitcell info in the XML and XYZ files'
!         stop
!       else
!         write(*,'(a,f20.10)')'INFO:unitcell info is obtained from the XYZ file: AX [au] =',cell_info(1)
!         write(*,'(a,f20.10)')'INFO:unitcell info is obtained from the XYZ file: AY [au] =',cell_info(2)
!         write(*,'(a,f20.10)')'INFO:unitcell info is obtained from the XYZ file: AZ [au] =',cell_info(3)
!       endif   
!       generate%unitcell%vectorA(1)=cell_info(1)
!       generate%unitcell%vectorB(2)=cell_info(2)
!       generate%unitcell%vectorC(3)=cell_info(3)
!     else
!       call cluster_load( generate%vcluster(i), item(vcluster_node,i-1) )
!     endif   
!   end do
!
    ! get <atom_id_is_added> node
    generate%atom_id_is_added_flag = .false.
    atom_id_node => getFirstElementByTagName(generate_node,"atom_id_is_added")
    if( associated(atom_id_node) ) then
      value        =''
      value        = getAttribute(atom_id_node,"mode")
!     write(*,*)'info:',trim(value)
      if (trim(value) /= "off") then 
        generate%atom_id_is_added_flag = .true.
        write(*,*)'INFO: <atom_id_is_added> tag is detected'
      endif  
    endif   

    ! get <group_id> node
    generate%n_group=-1  ! dummy setting
    group_node => getFirstElementByTagName(generate_node,"group_id_is_added")
    if( associated(group_node) ) then
      value        =''
      value        = getAttribute(group_node,"mode")
      if (trim(value) /= "off") then
        filename_wrk = ''
        filename_wrk = getAttribute(group_node,"file")
        write(*,*)'INFO: <group_id_is_added> tag is detected: file=',trim(filename_wrk)
        call read_group_file(trim(filename_wrk),generate)
      endif  
    endif   


    call destroyNode(document_node)

    call chdir(cwd)

  end subroutine generate_load
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  ! parse <cluster> node
  subroutine cluster_load( cluster, cluster_node, cell_info )
    implicit none
    type(cluster_type), intent(out) :: cluster
    type(fnode), pointer     :: cluster_node

    type(fnode), pointer     :: node
    character(len=256)       :: value, unit
    integer :: j
    real(8),       optional  :: cell_info(3)
!
    ! get class attribute
    value = getAttribute(cluster_node,"class")
    if( value == "" ) then
       cluster%class = ""
    else
       read(unit=value,fmt=*) cluster%class
    end if

    ! get cluster attribute
    value = getAttribute(cluster_node,"structure")

    if( value == "" ) then
       write(*,*) '<cluster structure> not found'
       stop
    else if( index( value, ".xml" )+3 == len(trim(value)) ) then
       call structure_load( cluster%structure, value )
    else if( index( value, ".xyz" )+3 == len(trim(value)) ) then
       if ( present(cell_info) ) then 
         call structure_loadXYZ( cluster%structure, value, cell_info )
!        write(*,*)'after structure_loadXYZ : cell_info=',cell_info(1:3)
       else
         call structure_loadXYZ( cluster%structure, value )
       endif   
    end if

    return
  end subroutine cluster_load
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  ! output <generate> node in the splitted files
  subroutine generate_output_split( generate, filename, n_split )
    implicit none
    type(generate_type), intent(in) :: generate
    character(len=*),    intent(in) :: filename
    integer,             intent(in) :: n_split
    character(len=64)               :: filename_header
    character(len=64)               :: filename_wrk
    integer                         :: lenf
    integer                         :: noa
    integer, parameter              :: i_cluster = 1
    logical, parameter              :: atom_id_is_added     = .true.
    logical, parameter              :: atom_id_is_not_added = .false.
    integer                         :: j_split
    integer                         :: atom_ini, atom_fin, atom_num 
    character(len=64)               :: chara_wrk
!
    lenf=len_trim(filename)
    filename_header=filename(1:lenf-4)
!   write(*,*)'filename_header =',trim(filename_header)
!
    filename_wrk=filename_header(1:lenf-4)//'_basic.xml'
!   write(*,*)'filename_basic  =',trim(filename_wrk)
    call generate_output(generate, filename_wrk, 0, 0, atom_id_is_added, 0,0 )
!
    noa=structure_info%total_number_of_atoms
!   noa=generate%vcluster(i_cluster)%structure%natom
    j_split=0 ! dummy value
!
    if (n_split == 0) then 
      call generate_output(generate, filename, 1, noa, atom_id_is_not_added, j_split, n_split)
      return
    endif  
!
    if (n_split == 1) then 
      call generate_output(generate, filename, 1, noa, atom_id_is_added, j_split, n_split)
      return
    endif  
!
    do j_split=0, n_split-1
      call a_x_b_divided_by_c_i8(j_split, noa, n_split, atom_ini)
      atom_ini = atom_ini +1
!       -->  atom_ini = j_split * noa / n_split + 1
      call a_x_b_divided_by_c_i8(j_split+1, noa, n_split, atom_fin)
!       -->  atom_fin = (j_split+1) * noa / n_split
      atom_num = atom_fin - atom_ini + 1
      write(chara_wrk, '(i6.6)') j_split
      filename_wrk=filename_header(1:lenf-4)//'_'//trim(chara_wrk)//'.xml'
      call generate_output(generate, trim(filename_wrk), atom_ini, atom_num, atom_id_is_added, j_split, n_split)
    enddo
!   
  end subroutine generate_output_split
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  ! output <generate> node
  subroutine generate_output( generate, filename, atom_ini, atom_num, atom_id_is_added, j_split, n_split )
    implicit none
    type(generate_type), intent(in) :: generate
    character(len=*),    intent(in) :: filename
    integer,             intent(in) :: atom_ini, atom_num
    logical,             intent(in) :: atom_id_is_added
    integer,             intent(in) :: j_split, n_split
    logical, parameter :: total_number_of_atoms_is_added = .true.
    integer, parameter :: fd = 1
    real(8), parameter :: max_cell_length = 1.0d10 ! used only for error checking
    integer :: i, j, k
    type(cluster_type), pointer :: cluster
    type(atom_type), pointer :: atom
    real(8) :: a, b, c
    real(8) :: ax, ay, az
    real(8) :: pos(3)
    character(len=256) :: class
    integer            :: atom_fin
    integer, parameter :: i_cluster = 1
    integer            :: noa
    character(len=64)  :: chara_wrk1, chara_wrk2, chara_wrk3, chara_wrk4, chara_wrk5
    character(len=1024) :: chara_wrk_long
    logical             :: atom_id_is_added_wrk
    logical             :: group_id_is_added_wrk
    integer             :: nog, noga
    integer             :: group_id
    integer             :: ierr
!
!   write(*,*)'generate_output:filename=', trim(filename)
!   write(*,*)' generate%unitcell%vectorA(1)=', generate%unitcell%vectorA(1)
!   write(*,*)' generate%unitcell%vectorB(2)=', generate%unitcell%vectorB(2)
!   write(*,*)' generate%unitcell%vectorC(3)=', generate%unitcell%vectorC(3)
!
    atom_id_is_added_wrk = atom_id_is_added
    if (generate%atom_id_is_added_flag) atom_id_is_added_wrk = .true.
!
    if (allocated(generate%group_id_booking)) then
      group_id_is_added_wrk = .true. 
    else
      group_id_is_added_wrk = .false. 
    endif
!
    noa=structure_info%total_number_of_atoms
!   noa=generate%vcluster(i_cluster)%structure%natom
    nog =generate%n_group
    noga=generate%max_num_group_atoms
!
    if (group_id_is_added_wrk) then
      ierr=0 
      if (nog < 1)   ierr=1
      if (nog > noa) ierr=1
      if (noga < 1)   ierr=1
      if (noga > noa) ierr=1
      if (ierr /=0) then
        write(*,*)'ERROR:n_group=',nog
        write(*,*)'ERROR:noga   =',noga
        stop
      endif   
      if (.not. allocated(generate%group_id_booking)) then
        write(*,*)'ERROR:group_id_booking is not allocated'
        stop
      endif
    endif
!
    atom_fin=atom_ini+atom_num-1
    if (atom_fin > noa) then
      write(*,*)'ERROR:atom_fin,noa=',atom_fin,noa
      stop
    endif
!
    if (n_split >= 2 ) then
      write(*,*)          'generate_output:filename=',trim(filename)
      write(*,'(a,3i10)')'generate_output :atom_ini, atom_fin, atom_num =',atom_ini, atom_fin, atom_num
    endif  
!
    open(fd,file=filename)
!
    write(fd,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(fd,'(a,a,a)') '<structure name="', trim(generate%name), '" mdstep="0">'
!
    write(fd,*) ""
!
!   call unitcell_save( fd, generate%unitcell )
!
    if (total_number_of_atoms_is_added) then
      write(chara_wrk1, *) noa
      chara_wrk1=adjustl(chara_wrk1)
      chara_wrk_long = ' <total_number_of_atoms> '//trim(chara_wrk1)//' </total_number_of_atoms>'
      write(fd,'(a)') trim(chara_wrk_long)
      write(fd,*) ""
    endif
!
    ax = structure_info%cell_size_angst(1)/angst
    ay = structure_info%cell_size_angst(2)/angst
    az = structure_info%cell_size_angst(3)/angst
!
    ierr=0
    if ((ax > max_cell_length) .or. (ax < 0.0d0)) ierr=1
    if ((ay > max_cell_length) .or. (ay < 0.0d0)) ierr=1
    if ((az > max_cell_length) .or. (az < 0.0d0)) ierr=1
!
    if (ierr /= 0) then
      write(*,*)'ERROR(elses-xml-generate):ax=',ax
      write(*,*)'ERROR(elses-xml-generate):ay=',ay
      write(*,*)'ERROR(elses-xml-generate):az=',az
      stop
    endif
!

!
    write(fd,*) ""
    write(fd,*) "<unitcell>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', ax, 0.0d0, 0.0d0, "</vector>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', 0.0d0, ay, 0.0d0, "</vector>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', 0.0d0, 0.0d0, az, "</vector>"
    write(fd,*) "</unitcell>"
    write(fd,*) ""
!
    if (generate%boundary_condition_set) then
      chara_wrk1=trim(generate%boundary_condition(1))
      chara_wrk2=trim(generate%boundary_condition(2))
      chara_wrk3=trim(generate%boundary_condition(3))
      write(fd,*) ""
      chara_wrk_long = ' <boundary x="'//trim(chara_wrk1)//'" y="'//trim(chara_wrk2) & 
&                         //'" z="'//trim(chara_wrk3)//'" />'
      write(fd,'(a)') trim(chara_wrk_long)
      write(fd,*) ""
    endif
!
    call heatbath_save( fd, generate%heatbath )
!
    if (group_id_is_added_wrk) then
      write(fd,'(a)') ' <group_id_is_added> '
      write(chara_wrk1, *) nog
      chara_wrk1=adjustl(chara_wrk1)
      chara_wrk_long = '   <number_of_groups> '//trim(chara_wrk1)//' </number_of_groups>'
      write(fd,'(a)') trim(chara_wrk_long)
      write(chara_wrk1, *) noga
      chara_wrk1=adjustl(chara_wrk1)
      chara_wrk_long = '   <max_number_of_group_atoms> '//trim(chara_wrk1)//' </max_number_of_group_atoms>'
      write(fd,'(a)') trim(chara_wrk_long)
      write(fd,'(a)') ' </group_id_is_added> '
      write(fd,*) ""
    endif
!
    if (n_split >= 2 ) then
      write(chara_wrk1, *) j_split
      chara_wrk1=adjustl(chara_wrk1)
      write(chara_wrk2, *) n_split
      chara_wrk2=adjustl(chara_wrk2)
      write(chara_wrk3, *) atom_ini
      chara_wrk3=adjustl(chara_wrk3)
      write(chara_wrk4, *) atom_fin
      chara_wrk4=adjustl(chara_wrk4)
      write(chara_wrk5, *) atom_num
      chara_wrk5=adjustl(chara_wrk5)
      chara_wrk_long = ' <split file_index="'//trim(chara_wrk1)//'" number_of_files="'//trim(chara_wrk2) & 
&                         //'" atoms_in_this_file="'//trim(chara_wrk5) &
&                         //'" atom_initial="'//trim(chara_wrk3)//'" atom_final="'//trim(chara_wrk4)//'" />'
      write(fd,'(a)') trim(chara_wrk_long)
      write(fd,*) ""
    endif   
!
    do j=atom_ini, atom_fin
      pos(1:3) = structure_info%atom_info(j)%position_angst(1:3) / structure_info%cell_size_angst(1:3)
      group_id = -1
      if (allocated(generate%group_id_booking)) then
        group_id =  generate%group_id_booking(j) 
      endif
      class = ''

      if (atom_id_is_added_wrk) then
        write(chara_wrk1, *) j
        chara_wrk1=adjustl(chara_wrk1)
        if (group_id_is_added_wrk) then
          write(chara_wrk2, *) group_id
          chara_wrk1=adjustl(chara_wrk1)
          chara_wrk2=adjustl(chara_wrk2)
          write(fd,*) '<atom element="', trim(structure_info%atom_info(j)%element_name), '" ', &
           'id="', trim(chara_wrk1), '" ', 'group_id="', trim(chara_wrk2), '" ',  &
           'motion="', trim(structure_info%atom_info(j)%motion), '">'
        else
          write(fd,*) '<atom element="', trim(structure_info%atom_info(j)%element_name), '" ', &
           'id="', trim(chara_wrk1), '" ', 'class="', trim(class), '" ', &
           'motion="', trim(structure_info%atom_info(j)%motion), '">'
        endif   
      else
        write(fd,*) '<atom element="', trim(structure_info%atom_info(j)%element_name), '" ', &
           'motion="', trim(structure_info%atom_info(j)%motion), '">'
      endif   
!
      write(fd,'(a28,3e23.15,a11)') '  <position unit="internal">', pos(1:3), '</position>'

      write(fd,*) '</atom>'
      write(fd,*) ""

    end do

    write(fd,*) '</structure>'

    close(fd)

    return
  end subroutine generate_output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set number of splitted files
  subroutine set_split( argv3, n )
    implicit none
    character(len=*),intent(in)  :: argv3
    integer,         intent(out) :: n
    integer                      :: ierr
!
    n=0
!
    if( argv3(1:7) /= "-split=" ) then
      write(*,*) 'ERROR:argv3=',trim(argv3)
      stop 
    endif   
!
    read(unit=argv3(8:), fmt=*, iostat=ierr) n
    if (ierr /= 0) then
      write(*,*) 'ERROR:argv3(8:)=',trim(argv3(8:))
      stop
    endif   
!
  end subroutine set_split
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine a_x_b_divided_by_c_i8(a, b, c, d)
!      integer(8) calculation of d = a*b/c 
!
   implicit none
   integer,           intent(in)  :: a, b, c
   integer,           intent(out) :: d
   integer(kind=8) :: a_i8, b_i8, c_i8, d_i8
!
   a_i8=a
   b_i8=b
   c_i8=c
   d_i8=a_i8*b_i8/c_i8
!
   d=d_i8
!
  end subroutine a_x_b_divided_by_c_i8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_group_file(filename,generate)
!
   implicit none
   character(len=*),  intent(in)     :: filename
   type(generate_type),intent(inout) :: generate
   character(len=10240) :: chara_wrk
   integer              :: i, j, n_max_comment_lines
   integer, parameter   :: fd = 2
   integer              :: ierr
   integer              :: noa      ! number of atoms
   integer              :: nog      ! number of groups
   integer              :: noga     ! max. number of group atoms
   integer              :: atom_id  ! atom_id
   integer              :: group_id ! group_id
!  logical, parameter   :: verbose_mode = .true.
   logical, parameter   :: verbose_mode = .false.
!
   open(fd,file=filename,status='old')
!
   n_max_comment_lines = 1000   ! possible maximum number of comment lines
!
   do j=1, n_max_comment_lines
     read(fd, '(a)', iostat=ierr) chara_wrk 
     if (ierr /=0) exit
     if (index(chara_wrk,"#") /= 1) exit  ! ignore the comment line
     if (verbose_mode) write(*,*)'INFO:comment line:', trim(chara_wrk)
   enddo   
!
   read(chara_wrk, * , iostat=ierr) noa, nog, noga
!
   if (ierr /= 0) then
     write(*,*)'ERROR:read_group_file' 
     stop
   endif   
!
   ierr=0
   if (noa  < 1) ierr =1
   if (nog  < 1) ierr =1
   if (noga < 1) ierr =1
   if (noga > noa) ierr =1
   if (ierr /= 0) then
     write(*,*)'ERROR:read_group_file:noa, nog, noga=', noa, nog, noga
     stop
   endif
!
   write(*,*)'INFO:number of atoms   = ', noa
   write(*,*)'INFO:number of groups  = ', nog
   write(*,*)'INFO:maximum number of group atoms  = ', noga
   generate%n_group=nog
   generate%max_num_group_atoms=noga
!
   allocate(generate%group_id_booking(noa),stat=ierr)
   if (ierr /= 0) then
     write(*,*)'Alloc. error: group_id_booking'
     stop
   endif   
!
   do j=1, noa+n_max_comment_lines
     read(fd, '(a)', iostat=ierr) chara_wrk 
     if (ierr /= 0) exit
     if (index(chara_wrk,"#") == 1) cycle  ! ignore the comment line
!    write(*,*)'j,line   = ',j, trim(chara_wrk)
     read(chara_wrk, * , iostat=ierr) atom_id, group_id
     if (ierr /= 0) exit
     ierr=0
     if (atom_id  < 1  ) ierr=1
     if (atom_id  > noa) ierr=1
     if (group_id < 1  ) ierr=1
     if (group_id > nog) ierr=1
     if (ierr /=0) then
       write(*,*)'ERROR:read_group_file:atom_id, group_id=',atom_id, group_id
       stop 
     endif   
     generate%group_id_booking(atom_id)=group_id
     if (verbose_mode) write(*,*)'INFO:atom ID, group ID = ',atom_id, group_id
   enddo
   if (verbose_mode) write(*,*)'INFO:The file end is found.'
!
   close(fd)
!
  end subroutine read_group_file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_info_from_xyz_file( filename, cell_info_in_xyz_file)
   implicit none
   character(len=*),    intent(in)    :: filename
   logical         ,    intent(inout) :: cell_info_in_xyz_file
   integer, parameter :: fd = 1
   integer            :: noa
   real(8)            :: ax, ay, az  ! cell sizes in Agnstrom
   integer            :: ierr, j
   character(len=1024)  :: chara_wrk
   character(len=8)     :: elem_name_wrk ! work array
   real(8)              :: data_wrk(3)   ! work array
!
   write(*,*)'INFO:get_basic_info: filename=',trim(filename)
!
   if (len_trim(filename) == 0) then
     write(*,*) 'ERROR(get_info_from_xyz_file): No XYZ filename '
     stop
   endif
!
   open(fd,file=filename,status='old')
!
   structure_info%total_number_of_atoms = -1     ! dummy value
   structure_info%cell_size_angst(3)    = -1.0d0 ! dummy value
   read(fd, '(a)') chara_wrk  ! read the first line
!  write(*,*)'INFO:chara_wrk        =', trim(chara_wrk)
   read(chara_wrk, *, iostat=ierr) noa, ax, ay, az
   if (ierr == 0) then
     cell_info_in_xyz_file = .true. 
!    write(*,*)'INFO:number of atoms        =', noa
!    write(*,*)'INFO:cell sizes [Angstrom]  =', ax, ay, az
     structure_info%total_number_of_atoms = noa
     structure_info%cell_size_angst(1)    = ax
     structure_info%cell_size_angst(2)    = ay
     structure_info%cell_size_angst(3)    = az
   else   
     read(chara_wrk, *, iostat=ierr) noa
     cell_info_in_xyz_file = .false. 
     structure_info%total_number_of_atoms = noa
!    write(*,*)'INFO:number of atoms        =', noa
   endif   
!
   allocate(structure_info%atom_info(noa), stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR:alloc. error for structure_info%atom_info:noa=',noa
     stop
   endif
!
   read(fd, '(a)') chara_wrk   ! read the second (comment) line 
!
   do j=1,noa
     read(fd, *) elem_name_wrk, data_wrk(1:3)
!    write(*,*)'INFO:j, name=', j, trim(elem_name_wrk)
!    write(*,*)'INFO:j, data(1:3)=', j, data_wrk(1:3)
     structure_info%atom_info(j)%element_name = trim(elem_name_wrk)
     structure_info%atom_info(j)%position_angst(1:3)=data_wrk(1:3)
   enddo
!
   structure_info%atom_info(:)%motion = 'free'    ! default setting
!
   close(fd)
  end subroutine get_info_from_xyz_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_fixed_atom (filename)
   implicit none
   character(len=*),    intent(in)    :: filename
   integer, parameter :: fd = 1
   integer            :: j_atom
   integer            :: ierr, j, n_list, n_list_max
!
    n_list_max=structure_info%total_number_of_atoms+1
!
   if (len_trim(filename) == 0) then 
     write(*,*)'INFO:set_fixed_atom ... is ignored'
     return
   else
     write(*,*)'INFO:set_fixed_atom: filename =',trim(filename)
   endif
!
   open(fd,file=filename,status='old')
!
   do j=1,n_list_max
     if (j == n_list_max) then
       write(*,*)'ERROR(set_fixed_atom);Too many list ?:j=',j
     endif
     read(fd, *, iostat=ierr) j_atom
     if (ierr /= 0) then 
       n_list=j-1
       exit
     endif
     write(*,*)'j=',j
     structure_info%atom_info(j_atom)%motion = 'fixed'  
   enddo
!
   write(*,*)'INFO:set_fixed_atom: number of fixed atoms =',n_list
!
   close(fd)
  end subroutine set_fixed_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine get_info_from_command_argument(filename_xyz, filename_setting, & 
&               filename_output, filename_fixed_atom_list, n_split, ierr_result)
   implicit none
    character(len=*), intent(out) :: filename_xyz, filename_setting, filename_output, filename_fixed_atom_list
    integer,          intent(out) :: n_split
    integer,          intent(out) :: ierr_result
    integer             :: iargc  
    character(len=1024) :: arg, chara_header, chara_wrk
    integer             :: i, len_arg, len_header, ierr
    logical             :: debug=.false.

!    
    n_split=0             ! default (non-split) setting
    filename_xyz     =''  ! dummy value
    filename_setting =''  ! dummy value
    filename_output  =''  ! dummy value
    filename_fixed_atom_list  =''  ! dummy value
    ierr_result      = 0  ! dummy value
!
    write(*,'(a)')'@@@ get_info_from_command_argument'
!
    ierr=1
    if ( iargc() == 2 ) ierr=0
    if ( iargc() == 3 ) ierr=0
    if ( iargc() == 4 ) ierr=0
!
    if (ierr /=0) then 
     ierr_result=ierr
     return
    endif
!
    do i=1, iargc()
      call getarg(i, arg)
      len_arg=len_trim(arg)
!
      chara_header='-xyz_file='
      if (index(arg, trim(chara_header)) > 0) then
        len_header=len_trim(chara_header)
        chara_wrk=trim(arg(len_header+1:len_arg))
        if (debug) write(*,*)' INFO:command argument : ', i, trim(chara_wrk)
        read(chara_wrk, *, iostat=ierr)  filename_xyz
        if (ierr /= 0) stop 'ERROR (get_command_argument) : filename_xyz'
        cycle
      endif
!
      chara_header='-setting_file='
      if (index(arg, trim(chara_header)) > 0) then
        len_header=len_trim(chara_header)
        chara_wrk=trim(arg(len_header+1:len_arg))
        if (debug) write(*,*)' INFO:command argument : ', i, trim(chara_wrk)
        read(chara_wrk, *, iostat=ierr)  filename_setting
        if (ierr /= 0) stop 'ERROR (get_command_argument) : filename_xyz'
        cycle
      endif
!
      chara_header='-output_file='
      if (index(arg, trim(chara_header)) > 0) then
        len_header=len_trim(chara_header)
        chara_wrk=trim(arg(len_header+1:len_arg))
        if (debug) write(*,*)' INFO:command argument : ', i, trim(chara_wrk)
        read(chara_wrk, *, iostat=ierr)  filename_output
        if (ierr /= 0) stop 'ERROR (get_command_argument) : filename_xyz'
        cycle
      endif
!
      chara_header='-split='
      if (index(arg, trim(chara_header)) > 0) then
        len_header=len_trim(chara_header)
        chara_wrk=trim(arg(len_header+1:len_arg))
        if (debug) write(*,*)' INFO:command argument : ', i, trim(chara_wrk)
        read(chara_wrk, *, iostat=ierr)  n_split
        if (ierr /= 0) stop 'ERROR (get_command_argument) : filename_xyz'
        cycle
      endif
!
      chara_header='-fixed_atom_list='
      if (index(arg, trim(chara_header)) > 0) then
        len_header=len_trim(chara_header)
        chara_wrk=trim(arg(len_header+1:len_arg))
        if (debug) write(*,*)' INFO:command argument : ', i, trim(chara_wrk)
        read(chara_wrk, *, iostat=ierr)  filename_fixed_atom_list
        if (ierr /= 0) stop 'ERROR (get_command_argument) : filename_fixed_atom_list'
        cycle
      endif
!
      write(*,*)'INFO: non-header argument : ', trim(arg)
!
      if (len_trim(filename_setting) == 0) then 
        filename_setting=trim(arg)
        cycle
      endif
!
      if (len_trim(filename_output)  == 0) then 
        filename_output=trim(arg)
        cycle
      endif
!
    enddo
!
    write(*,*)'INFO:filename of the xyz file     =', trim(filename_xyz)
    write(*,*)'INFO:filename of the setting file =', trim(filename_setting)
    write(*,*)'INFO:filename of the output  file =', trim(filename_output)
    write(*,*)'INFO:filename of the fixed atom list file =', trim(filename_fixed_atom_list)
    if (n_split > 0) then
      write(*,*)'INFO: the number of split files   =', n_split
    else
      write(*,*)'INFO: the single (non-split) file will be generated'
    endif
!    
  end subroutine get_info_from_command_argument
!

!
end program elses_xml_generate
