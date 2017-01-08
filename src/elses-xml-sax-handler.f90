!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_sax_handler
!
  use M_config, only : config  !(CHANGED)
  use elses_xml_misc, only : XML_TO_AU !(function)
  implicit none
! logical, parameter :: tag_dump = .true. ! true only for debugging
! logical, parameter :: tag_dump = .false. ! true only for debugging
!
  logical :: in_structure_first
  logical :: in_structure_history
  logical :: in_structure, in_unitcell, in_vector, in_heatbath, in_atom 
  logical :: in_myLength
!
  logical :: in_heatbath_massperatom
  logical :: in_heatbath_position
  logical :: in_heatbath_velocity
  logical :: in_atom_position
  logical :: in_atom_velocity
  logical :: in_atom_population
  logical :: in_atom_population_guess
  integer :: cell_vector_counter
  integer :: counter_atom
  real(8) :: unit_conv
  logical :: in_atom_position_internal
  logical :: in_atom_force
  logical :: in_split
  integer :: atom_id
!
  private
  public :: in_structure_first
  public :: begin_element
  public :: pcdata_chunk
  public :: end_element
  public :: cell_vector_counter
  public :: counter_atom
  public :: xml_flag_init
  public :: xml_flag_check

contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xml_flag_init
    implicit none
    in_structure_history    = .false.
    in_structure            = .false.
    in_unitcell             = .false.
    in_vector               = .false.
    in_heatbath             = .false. 
    in_atom                 = .false.
    in_myLength             = .false.
    in_heatbath_massperatom = .false.
    in_heatbath_position    = .false.
    in_heatbath_velocity    = .false.
    in_atom_position        = .false.
    in_atom_velocity        = .false.
    in_atom_population      = .false.
    in_atom_population_guess  = .false.
    in_atom_position_internal = .false.
    in_atom_force             = .false.
    in_split                  = .false.
    atom_id      = -1
  end subroutine xml_flag_init
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xml_flag_check(error_count)
    implicit none
    integer, intent(out) :: error_count 
    error_count=0
    if (in_structure) error_count=error_count+1
    if (in_unitcell) error_count=error_count+1
    if (in_vector) error_count=error_count+1
    if (in_heatbath) error_count=error_count+1 
    if (in_atom) error_count=error_count+1
    if (in_myLength) error_count=error_count+1
    if (in_heatbath_massperatom) error_count=error_count+1
    if (in_heatbath_position) error_count=error_count+1
    if (in_heatbath_velocity) error_count=error_count+1
    if (in_atom_position) error_count=error_count+1
    if (in_atom_velocity) error_count=error_count+1
    if (in_atom_population) error_count=error_count+1
    if (in_atom_population_guess) error_count=error_count+1
    if (in_atom_position_internal) error_count=error_count+1
    if (in_atom_force) error_count=error_count+1
    if (in_split) error_count=error_count+1
    if (error_count > 0) then
      write(*,*)'ERROR(xml_flag_check):error_count=' ,error_count
      stop
    endif   
  end subroutine xml_flag_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine begin_element(name, attributes)
!
    use flib_sax, only : dictionary_t, get_value
    implicit none
    character(len=*), intent(in) :: name
    type(dictionary_t), intent(in) :: attributes
    integer :: status, ierr
    character(len=1024) :: chara_wrk
    real(8) :: value_wrk
    integer :: value_int_wrk
    logical :: tag_dump
!
    tag_dump=config%system%structure%tag_dump
!
    select case(name)
      case('structure_history')
       in_structure_history = .true.
       write(*,*)'ERROR(XML-SAX):<structure_history> was found in the XML file and is not supported'
       stop
!
      case('structure')
       in_structure = .true.
       if (.not. in_structure_first) then
         write(*,*)'ERROR(XML-SAX):not supported now:in_structure_first=',in_structure_first
         stop
       endif   
       call get_value(attributes,"name", chara_wrk, status)
       if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
         write(*,*)'ERROR(XML-SAX):No structure name:status=',status
         stop
       endif   
       config%system%structure%name=trim(chara_wrk)
       if (tag_dump) write(*,'(a,a)')'INFO-XML-SAX:structure%name=',config%system%structure%name
       call get_value(attributes,"mdstep", chara_wrk, status)
       if (trim(chara_wrk) =='') then
         write(*,*)'ERROR(XML-SAX):No mdstep'
         stop
       endif  
       read(chara_wrk,*,iostat=ierr) config%system%structure%mdstep
       if (ierr /= 0) then
         write(*,*)'ERROR(XML-SAX):mdstep can not be obtained'
         stop
       endif   
       if (tag_dump) write(*,'(a,i10)')'INFO-XML-SAX:structure%mdstep=',config%system%structure%mdstep
!
      case('unitcell')
        in_unitcell = .true. 
!
      case('split')
        in_split = .true. 
        call get_value(attributes,"number_of_files", chara_wrk, status)
       if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
          write(*,*)'ERROR(XML-SAX):No number_of_files attirbute:status=',status
          stop
        endif   
        read(chara_wrk,*,iostat=ierr) value_int_wrk
        if (config%system%structure%split%number_of_files /= value_int_wrk) then
          write(*,*)'ERROR:INFO-XML-SAX:unmatched: number_of_files=', value_int_wrk
          stop
        else  
          if (tag_dump) write(*,*)'number_of_files=',value_int_wrk
        endif   
!
        call get_value(attributes,"file_index", chara_wrk, status)
        if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
          write(*,*)'ERROR(XML-SAX):No file_index attirbute:status=',status
          stop
        endif   
        read(chara_wrk,*,iostat=ierr) config%system%structure%split%file_index
        if (tag_dump) write(*,*)'file_index =', config%system%structure%split%file_index
!
        call get_value(attributes,"atoms_in_this_file", chara_wrk, status)
        if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
          write(*,*)'ERROR(XML-SAX):No atoms_in_this_file attirbute:status=',status
          stop
        endif   
        read(chara_wrk,*,iostat=ierr) config%system%structure%split%atoms_in_this_file
        if (tag_dump) write(*,*)'atoms_in_this_file =', config%system%structure%split%atoms_in_this_file
!
        call get_value(attributes,"atom_initial", chara_wrk, status)
        if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
          write(*,*)'ERROR(XML-SAX):No atom_initial attirbute:status=',status
          stop
        endif   
        read(chara_wrk,*,iostat=ierr) config%system%structure%split%atom_initial
        if (tag_dump) write(*,*)'atom_initial =', config%system%structure%split%atom_initial
!
        call get_value(attributes,"atom_final", chara_wrk, status)
        if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
          write(*,*)'ERROR(XML-SAX):No atom_final attirbute:status=',status
          stop
        endif   
        read(chara_wrk,*,iostat=ierr) config%system%structure%split%atom_final
        if (tag_dump) write(*,*)'atom_final =', config%system%structure%split%atom_final
!
      case('atom')
        in_atom = .true. 
        counter_atom=counter_atom+1
        if (tag_dump) write(*,*)'counter_atom (p) =',counter_atom
        if (counter_atom > config%system%structure%natom) then
          write(*,*)'ERROR(XML-SAX):counter_atom=',counter_atom
          stop
        endif   
!
        call get_value(attributes,"id", chara_wrk, status)
        if (status /= 0) then
          atom_id=counter_atom
        else  
          if (trim(chara_wrk) =='') then
            write(*,*)' ERROR:reading atom_id;chara_wrk=',trim(chara_wrk)
          endif   
          read(chara_wrk,*,iostat=ierr) atom_id
          if (ierr /= 0) then
            write(*,*)' ERROR:reading atom_id:ierr=',ierr
            write(*,*)' ERROR:reading atom_id:chara_wrk=',trim(chara_wrk)
            stop
          endif   
          if (config%system%structure%split%set) then
            if (atom_id < config%system%structure%split%atom_initial) then
              write(*,*)' ERROR:myrank, reading atom_id, initial=', config%calc%distributed%myrank,  &
&                        atom_id, config%system%structure%split%atom_initial
              stop
            endif   
            if (atom_id > config%system%structure%split%atom_final) then
              write(*,*)' ERROR:reading atom_id, final  =',atom_id, config%system%structure%split%atom_final
              write(*,*)' ERROR:myrank=',config%calc%distributed%myrank
              stop
            endif
          endif
        endif   
!
        call get_value(attributes,"element", chara_wrk, status)
        if ( (trim(chara_wrk) =='') .or. (status /=  0) ) then
          write(*,*)'ERROR(XML-SAX):No atom element tag:status=',status
          stop
        endif   
        if (config%system%structure%use_matom) then
          if ((atom_id < 1) .or. (atom_id > size(config%system%structure%matom))) then
            write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
            stop
          endif
          config%system%structure%matom(atom_id)%name=trim(chara_wrk)
          if (tag_dump) write(*,*)'  element name =',config%system%structure%matom(counter_atom)%name
        endif  
        if (config%system%structure%use_vatom) then
          if ((atom_id < 1) .or. (atom_id > size(config%system%structure%vatom))) then
            write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
            stop
          endif
          config%system%structure%vatom(atom_id)%name=trim(chara_wrk)
          if (tag_dump) write(*,*)'  element name =',config%system%structure%vatom(counter_atom)%name
        endif  

!
        call get_value(attributes,"motion", chara_wrk, status)
        if (config%system%structure%use_vatom) then
          if (trim(chara_wrk) =='') then
            config%system%structure%vatom(atom_id)%motion='free' 
          else   
            config%system%structure%vatom(atom_id)%motion=trim(chara_wrk)
          endif   
          if (tag_dump) write(*,*)'  motion name =',config%system%structure%vatom(counter_atom)%motion
        endif  
!
!       chara_wrk=''
        call get_value(attributes,"group_id", chara_wrk, status)
        if (config%system%structure%use_vatom) then
          if (status /= 0) then
            config%system%structure%vatom(atom_id)%group_id=-1
          else   
            read(chara_wrk,*,iostat=ierr) config%system%structure%vatom(atom_id)%group_id
            if (ierr /= 0) then
              write(*,*)'ERROR(sax XML handler for group_id): ',trim(chara_wrk)
              stop
            endif
            if (tag_dump) write(*,*)'  group_id    =',config%system%structure%vatom(atom_id)%group_id
          endif   
        endif  
!
      case('myLength')
        in_myLength = .true. 
        if (in_unitcell) then
          call get_value(attributes,"unit", chara_wrk, status)
          if (trim(chara_wrk) == '') then
            unit_conv = 1d0 
          else
            if (tag_dump) write(*,*) 'INFO-XML:optional tag:mylenth:unit=',trim(chara_wrk)
            unit_conv = XML_TO_AU(trim(chara_wrk))
          endif   
        else
          write(*,*)'ERROR(XML-SAX):Wrong XML structure:myLength'
          stop
        endif  
!
      case('vector')
        in_vector = .true.
        if (in_unitcell) then
          call get_value(attributes,"unit", chara_wrk, status)
          if (trim(chara_wrk) == '') then
            write(*,*)'ERROR(XML-SAX):No unit in the unitcell tag'
            stop
          else
           if (tag_dump) write(*,*) 'call XML_TO_AU for unitcell vector:',trim(chara_wrk)
           if (trim(chara_wrk) == 'myLength') then
             unit_conv=config%system%structure%unitcell%myLength
           else
             unit_conv=XML_TO_AU(trim(chara_wrk))
           endif 
          endif  
        endif  
!
      case('heatbath')
        in_heatbath = .true. 
!
      case('massperatom')
        if (in_heatbath) then 
          in_heatbath_massperatom = .true. 
          if (tag_dump) write(*,*)'massperatom tag'
        else
          write(*,*)'ERROR(XML-SAX):begin_element:tag error:massperatom'
          stop
        endif   
        call get_value(attributes,"unit", chara_wrk, status)
        if (trim(chara_wrk) == '') chara_wrk="a.u."
        if (tag_dump) write(*,*) 'call XML_TO_AU for heatbath massperatom:',trim(chara_wrk)
        unit_conv=XML_TO_AU(trim(chara_wrk))
!
      case('position')
        if (tag_dump) write(*,*)'position tag'
        if (in_heatbath) then 
          in_heatbath_position = .true.
          call get_value(attributes,"unit", chara_wrk, status)
          if (trim(chara_wrk) == '') chara_wrk="a.u."
          if (tag_dump) write(*,*) 'call XML_TO_AU for heatbath position:',trim(chara_wrk)
          unit_conv=XML_TO_AU(trim(chara_wrk))
        else
          if (in_atom) then 
            in_atom_position = .true.
            in_atom_position_internal = .false.
            call get_value(attributes,"unit", chara_wrk, status)
            unit_conv=1.0d0
            if (trim(chara_wrk) == '') then
              write(*,*)'ERROR(XML-SAX):No unit in the unitcell tag'
              stop
            else
              if (trim(chara_wrk) == 'internal') then 
                if (.not. config%system%structure%unitcell%set) then
                  write(*,*)'ERROR(XML-SAX):unitcell it not set, when atom position is loaded' 
                  stop 
                endif   
                unit_conv=1.0d0
                in_atom_position_internal = .true.
              else   
                if (tag_dump) write(*,*) 'call XML_TO_AU for atom position:',trim(chara_wrk)
                unit_conv=XML_TO_AU(trim(chara_wrk))
              endif  
            endif  
          else
            write(*,*)'ERROR(XML-SAX):tag error:position'
            stop
          endif   
        endif
!
      case('velocity')
        if (tag_dump) write(*,*)'velocity tag'
        if (in_heatbath) then 
          in_heatbath_velocity = .true.
          call get_value(attributes,"unit", chara_wrk, status)
          if (trim(chara_wrk) == '') chara_wrk="a.u."
          if (tag_dump) write(*,*) 'call XML_TO_AU for heatbath velocity:',trim(chara_wrk)
          unit_conv=XML_TO_AU(trim(chara_wrk))
        else
          if (in_atom) then 
            in_atom_velocity = .true.
          else
            write(*,*)'ERROR(XML-SAX):tag error:velocity'
            stop
          endif   
        endif
!
      case('population')
       if (in_atom) then 
         in_atom_population = .true.
         if (tag_dump) write(*,*)'INFO-XML:tag_dump: atom%population'
       else
         write(*,*)'ERROR(XML-SAX):tag error:population'
         stop
       endif   
!
      case('population_guess')
       if (in_atom) then 
         in_atom_population_guess = .true.
         if (tag_dump) write(*,*)'INFO-XML:tag_dump: atom%population_guess'
       else
         write(*,*)'ERROR(XML-SAX):tag error:population_guess'
         stop
       endif   
!
      case('force')
       if (in_atom) then 
         in_atom_force = .true.
         if (tag_dump) write(*,*)'INFO-XML:tag_dump: atom%force'
       else
         write(*,*)'ERROR(XML-SAX):tag error:force'
         stop
       endif   
!
      case('total_number_of_atoms','boundary','group_id_is_added','number_of_groups', 'max_number_of_group_atoms')
      write(*,'(a,a)')'INFO-XML-SAX:the beginning of the following tag was found and ignored;', trim(name)
!
      case default
      write(*,*)'ERROR:INFO-XML-SAX:unknown tag is found:', name
      stop
!
    end select
!
  end subroutine begin_element
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pcdata_chunk(chunk)
!
    implicit none
    character(len=*), intent(in) :: chunk
    integer :: status, ierr
    real(8) :: value_wrk
    real(8) :: value_wrk1, value_wrk2, value_wrk3
    real(8) :: atom_pos_wrk(3)
    logical :: tag_dump
!
    tag_dump=config%system%structure%tag_dump
!
    if (.not. in_structure) return
! 
    if ((in_unitcell) .and. (in_myLength)) then
      read(chunk,*,iostat=ierr) value_wrk
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX): myLength can not be obtained'
        stop
      endif   
      config%system%structure%unitcell%myLength=value_wrk*unit_conv
      if (tag_dump) write(*,*)'myLength (a.u.)=',config%system%structure%unitcell%myLength
    endif
!  
    if ((in_unitcell) .and. (in_vector)) then
      read(chunk,*,iostat=ierr) value_wrk1, value_wrk2, value_wrk3
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):unit_cell length can not be obtained'
!       write(31,*) chunk
        stop
      endif   
      cell_vector_counter=cell_vector_counter+1
      select case (cell_vector_counter)
        case (1)
          config%system%structure%unitcell%vectorA(1) = value_wrk1*unit_conv
          config%system%structure%unitcell%vectorA(2) = 0.0d0
          config%system%structure%unitcell%vectorA(3) = 0.0d0
          if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:cell length Lx [au] =',config%system%structure%unitcell%vectorA(1)
          if (tag_dump) write(*,'(a,2f20.10)')'INFO-XML-SAX:cell length Lx [au]2=',value_wrk1, unit_conv
        case (2)
          config%system%structure%unitcell%vectorB(1) = 0.0d0
          config%system%structure%unitcell%vectorB(2) = value_wrk2*unit_conv               
          config%system%structure%unitcell%vectorB(3) = 0.0d0
          if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:cell length Ly [au] =',config%system%structure%unitcell%vectorB(2)
          if (tag_dump) write(*,'(a,2f20.10)')'INFO-XML-SAX:cell length Ly [au]2=',value_wrk2, unit_conv
        case (3)
          config%system%structure%unitcell%vectorC(1) = 0.0d0
          config%system%structure%unitcell%vectorC(2) = 0.0d0
          config%system%structure%unitcell%vectorC(3) = value_wrk3*unit_conv               
          config%system%structure%unitcell%set        = .true.
          if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:cell length Lc [au] =',config%system%structure%unitcell%vectorC(3)
          if (tag_dump) write(*,'(a,2f20.10)')'INFO-XML-SAX:cell length Lc [au]2=',value_wrk3, unit_conv
        case default 
          write(*,*)'Error(XML-SAX):Too many unit cell vectors' 
          stop
      end select
    endif  
!
    if ((in_heatbath) .and. (in_heatbath_massperatom)) then
      read(chunk,*,iostat=ierr) value_wrk1
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):massperatom is not set'
!       write(31,*) chunk
        stop
      endif   
      config%system%structure%heatbath%massperatom=value_wrk1*unit_conv
      if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:massperatom=',config%system%structure%heatbath%massperatom
    endif   
!
    if ((in_heatbath) .and. (in_heatbath_position)) then
      read(chunk,*,iostat=ierr) value_wrk1
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):heatbath position is not set'
!       write(31,*) chunk
        stop
      endif   
      config%system%structure%heatbath%position=value_wrk1*unit_conv
      if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:heatbath position=',config%system%structure%heatbath%position
    endif   
!
    if ((in_heatbath) .and. (in_heatbath_velocity)) then
      read(chunk,*,iostat=ierr) value_wrk1
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):heatbath velocity is not set'
        write(31,*) chunk
        stop
      endif   
      config%system%structure%heatbath%velocity=value_wrk1*unit_conv
      if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:heatbath velocity=',config%system%structure%heatbath%velocity
    endif   
!
    if ((in_atom) .and. (in_atom_position)) then
      read(chunk,*,iostat=ierr) value_wrk1, value_wrk2, value_wrk3
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):heatbath velocity is not set'
        write(31,*) chunk
        stop
      endif   
      if (config%system%structure%use_matom) then
        if ((atom_id < 1) .or. (atom_id > size(config%system%structure%matom))) then
          write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
          stop
        endif
        config%system%structure%matom(atom_id)%position(1)=value_wrk1*unit_conv
        config%system%structure%matom(atom_id)%position(2)=value_wrk2*unit_conv
        config%system%structure%matom(atom_id)%position(3)=value_wrk3*unit_conv
        if (tag_dump) write(*,'(a,3f20.10)')'INFO-XML-SAX:position(1:3)=', &
&                   config%system%structure%matom(atom_id)%position(1:3)
      endif
      if (config%system%structure%use_matom) then
        if ((atom_id < 1) .or. (atom_id > size(config%system%structure%matom))) then
          write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
          stop
        endif
        config%system%structure%matom(atom_id)%position(1)=value_wrk1*unit_conv
        config%system%structure%matom(atom_id)%position(2)=value_wrk2*unit_conv
        config%system%structure%matom(atom_id)%position(3)=value_wrk3*unit_conv
        if (tag_dump) write(*,'(a,3f20.10)')'INFO-XML-SAX:position(1:3)=', &
&                   config%system%structure%matom(atom_id)%position(1:3)
      endif
      if (config%system%structure%use_vatom) then
        if ((atom_id < 1) .or. (atom_id > size(config%system%structure%vatom))) then
          write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
          stop
        endif
        config%system%structure%vatom(atom_id)%position(1)=value_wrk1*unit_conv
        config%system%structure%vatom(atom_id)%position(2)=value_wrk2*unit_conv
        config%system%structure%vatom(atom_id)%position(3)=value_wrk3*unit_conv
        if (tag_dump) write(*,'(a,3f20.10)')'INFO-XML-SAX:position(1:3)=', &
&                   config%system%structure%vatom(atom_id)%position(1:3)
      endif
      if (in_atom_position_internal) then
        if (.not. config%system%structure%unitcell%set) then
          write(*,*)'ERROR(XML-SAX):unitcell it not set, when atom position is set'
          stop 
        endif   
        if (config%system%structure%use_matom) then
          atom_pos_wrk(1:3)=config%system%structure%matom(atom_id)%position(1:3)
          atom_pos_wrk(1)=atom_pos_wrk(1)*config%system%structure%unitcell%vectorA(1)
          atom_pos_wrk(2)=atom_pos_wrk(2)*config%system%structure%unitcell%vectorB(2)
          atom_pos_wrk(3)=atom_pos_wrk(3)*config%system%structure%unitcell%vectorC(3)
          config%system%structure%matom(atom_id)%position(1:3)=atom_pos_wrk(1:3)
        endif  
        if (config%system%structure%use_vatom) then
          atom_pos_wrk(1:3)=config%system%structure%vatom(atom_id)%position(1:3)
          atom_pos_wrk(1)=atom_pos_wrk(1)*config%system%structure%unitcell%vectorA(1)
          atom_pos_wrk(2)=atom_pos_wrk(2)*config%system%structure%unitcell%vectorB(2)
          atom_pos_wrk(3)=atom_pos_wrk(3)*config%system%structure%unitcell%vectorC(3)
          config%system%structure%vatom(atom_id)%position(1:3)=atom_pos_wrk(1:3)
        endif  
        in_atom_position_internal=.false.
      endif   
!
    endif   
!
    if ((in_atom) .and. (in_atom_velocity)) then
      read(chunk,*,iostat=ierr) value_wrk1, value_wrk2, value_wrk3
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):atom velocity is not set'
        write(31,*) chunk
        stop
      endif   
      if (config%system%structure%use_vatom) then
        if ((atom_id < 1) .or. (atom_id > size(config%system%structure%vatom))) then
          write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
          stop
        endif
        config%system%structure%vatom(atom_id)%velocity_set = .true. 
        config%system%structure%vatom(atom_id)%velocity(1)=value_wrk1*unit_conv
        config%system%structure%vatom(atom_id)%velocity(2)=value_wrk2*unit_conv
        config%system%structure%vatom(atom_id)%velocity(3)=value_wrk3*unit_conv
        if (tag_dump) write(*,'(a,3f20.10)')'INFO-XML-SAX:velocity(1:3)=', &
&                   config%system%structure%vatom(atom_id)%velocity(1:3)
      endif  
!
    endif   
!
    if ((in_atom) .and. (in_atom_force)) then
      read(chunk,*,iostat=ierr) value_wrk1, value_wrk2, value_wrk3
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):atom force is not set'
        write(31,*) chunk
        stop
      endif   
      if (config%system%structure%use_vatom) then
        if ((atom_id < 1) .or. (atom_id > size(config%system%structure%vatom))) then
          write(*,*)'ERROR(XML-SAX):atom_id=',atom_id
          stop
        endif
        config%system%structure%vatom(atom_id)%force_set = .true. 
        config%system%structure%vatom(atom_id)%force(1)=value_wrk1*unit_conv
        config%system%structure%vatom(atom_id)%force(2)=value_wrk2*unit_conv
        config%system%structure%vatom(atom_id)%force(3)=value_wrk3*unit_conv
        if (tag_dump) write(*,'(a,3f20.10)')'INFO-XML-SAX:velocity(1:3)=', &
&                   config%system%structure%vatom(atom_id)%velocity(1:3)
      endif  
!
    endif   

    if ((in_atom) .and. (in_atom_population)) then
      read(chunk,*,iostat=ierr) value_wrk1
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):atom_population is not set'
        write(*,*) chunk
        stop
      endif   
      if (config%system%structure%use_vatom) then
        config%system%structure%vatom(counter_atom)%population_set = .true. 
        config%system%structure%vatom(counter_atom)%population     = value_wrk1
        if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:population=', &
&                   config%system%structure%vatom(counter_atom)%population
      endif  
    endif   
!
    if ((in_atom) .and. (in_atom_population_guess)) then
      read(chunk,*,iostat=ierr) value_wrk1
      if (ierr /= 0) then
        write(*,*)'ERROR(XML-SAX):atom_population_guess is not set'
        write(*,*) chunk
        stop
      endif   
      if (config%system%structure%use_vatom) then
        config%system%structure%vatom(counter_atom)%population_guess_set = .true. 
        config%system%structure%vatom(counter_atom)%population_guess     = value_wrk1
        if (tag_dump) write(*,'(a,f20.10)')'INFO-XML-SAX:population_guess=', &
&                   config%system%structure%vatom(counter_atom)%population_guess
      endif  
    endif   
!
   end subroutine pcdata_chunk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine end_element(name)
     implicit none
     character(len=*), intent(in) :: name
     logical :: tag_dump
     tag_dump=config%system%structure%tag_dump
     select case(name)
       case('structure_history') 
         in_structure_history= .false. 
       case('structure') 
         in_structure= .false. 
       case('unitcell') 
         in_unitcell= .false. 
       case('myLength') 
         in_myLength = .false. 
       case('vector') 
         in_vector= .false. 
      case('heatbath')
        in_heatbath = .false. 
      case('split')
        in_split = .false. 
      case('atom')
        in_atom = .false. 
      case('massperatom')
        if (in_heatbath) then 
          in_heatbath_massperatom = .false.
        else
          write(*,*)'ERROR(XML-SAX);end_element:tag error:massperatom'
          stop
        endif   
      case('position')
        if (in_heatbath) then 
          in_heatbath_position = .false.
        else
          if (in_atom) then 
            in_atom_position = .false.
          else
            write(*,*)'ERROR(XML-SAX):end_element:tag error:position'
            stop
          endif   
        endif
      case('velocity')
        if (in_heatbath) then 
          in_heatbath_velocity = .false.
        else
          if (in_atom) then 
            in_atom_velocity = .false.
          else
            write(*,*)'ERROR(XML-SAX):end_element:tag error:velocity'
            stop
          endif   
        endif
      case('population')
        if (in_atom) then 
            in_atom_population = .false.
        else
          write(*,*)'ERROR(XML-SAX):end_element:tag error:population'
          stop
        endif   
!
      case('population_guess')
        if (in_atom) then 
            in_atom_population_guess = .false.
        else
          write(*,*)'ERROR(XML-SAX):end_element:tag error:population_guess'
          stop
        endif   

      case('force')
        if (in_atom) then 
            in_atom_force = .false.
        else
          write(*,*)'ERROR(XML-SAX):end_element:tag error:force'
          stop
        endif   

      case('total_number_of_atoms','boundary','group_id_is_added','number_of_groups', 'max_number_of_group_atoms')
      write(*,'(a,a)')'INFO-XML-SAX:the end of the following tag was found and ignored;', trim(name)

      case default
        write(*,*)'ERROR:INFO-XML-SAX:unknown tag end is found:', name
        stop

     end select   
!
   end subroutine end_element
!
end module M_sax_handler


