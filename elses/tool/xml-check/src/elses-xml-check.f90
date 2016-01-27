!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_xml_checker_hander

  implicit none
  logical :: unknown_tag 
  logical :: incorrect_hierarchy
!
  logical :: in_heatbath
  logical :: in_atom
  logical :: in_output
!
  private
  public :: begin_element
  public :: pcdata_chunk
  public :: end_element
!
  public :: unknown_tag
  public :: incorrect_hierarchy

  public :: in_heatbath
  public :: in_atom
  public :: in_output
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine begin_element(name, attributes)
!
    use flib_sax, only : dictionary_t, get_value
    implicit none
    character(len=*), intent(in) :: name
    type(dictionary_t), intent(in) :: attributes
    integer :: status, ierr
    character(len=1024) :: chara_wrk
    logical :: check_hierarchy
!
    select case(name)
!!--------------------------------------------------
!!------------tags in the configuration file--------
!!--------------------------------------------------
!!------------config tag------------
      case('config')
       write(*,*)'tag starts:config'  
!!------------system tag------------
      case('system')
       write(*,*)'tag starts:system'  
      case('cluster')
       write(*,*)'tag starts:cluster'  
      case('splitted_input_file')
       write(*,*)'tag starts:splitted_input_file'
      case('boundary')
       write(*,*)'tag starts:boundary'  
      case('temperature')
       write(*,*)'tag starts:temperature'  
      case('element')
       write(*,*)'tag starts:element'  
!!------------calc tag------------
      case('calc')
       write(*,*)'tag starts:calc'  
!!------------calc/genoOption tag------------
      case('genoOption')
       write(*,*)'tag starts:genoOption'  
      case('CSC_method')
       write(*,*)'tag starts:CSC_method'  
      case('CSC_max_loop_count')
       write(*,*)'tag starts:CSC_max_loop_count'
      case('CSC_charge_convergence')
       write(*,*)'tag starts:CSC_charge_convergence'
      case('CSC_charge_mixing_ratio')
       write(*,*)'tag starts:CSC_charge_mixing_ratio'
      case('HML_kappa')
       write(*,*)'tag starts:HML_kappa'
      case('HML_small_delta')
       write(*,*)'tag starts:HML_small_delta'
      case('HML_constant_K')
       write(*,*)'tag starts:HML_constant_K'
      case('HML_K')
       write(*,*)'tag starts:HML_K'
!!------------calc/optimization tag------------
      case('optimization')
       write(*,*)'tag starts:optimization'
      case('sd_ratio')
       write(*,*)'tag starts:sd_ratio'
      case('max_num_iter')
       write(*,*)'tag starts:max_num_iter'
      case('convergence_criteria')
       write(*,*)'tag starts:convergence_criteria'
!!------------calc/dynamics tag------------
      case('dynamics')
       write(*,*)'tag starts:dynamics'
      case('delta')
       write(*,*)'tag starts:delta'
      case('total')
       write(*,*)'tag starts:total'
!!------------calc/solver tag------------
      case('solver')
       write(*,*)'tag starts:solver'
!!------------calc/limit tag------------
      case('limit')
       write(*,*)'tag starts:limit'
      case('time')
       write(*,*)'tag starts:time'
!!------------output tag------------
      case('output')
       in_output = .true.  
       write(*,*)'tag starts:output'
      case('restart')
       write(*,*)'tag starts:restart'
      case('wavefunction')
       write(*,*)'tag starts:wavefunction'
!!--------------------------------------------------
!!------------tags in the structure file------------
!!--------------------------------------------------
      case('structure')
       write(*,*)'tag starts:structure'
      case('split')
       write(*,*)'tag starts:split'
!!------------unitcell tag------------
      case('lengthScale')
       write(*,*)'tag starts:lengthScale'
      case('unitcell')
       write(*,*)'tag starts:unitcell'
      case('vector')
       write(*,*)'tag starts:vector'
!!------------heatbath tag------------
      case('heatbath')
       in_heatbath = .true.
       write(*,*)'tag starts:heatbath'
      case('massperatom')
       write(*,*)'tag starts:massperatom'
!!------------atom tag------------
      case('atom')
       in_atom = .true.
       write(*,*)'tag starts:atom'
!!------------multiple-defined tag:position ------
!!------------  (defined in output, atom and heatbath tags) ------
      case('position')
       check_hierarchy = .false.
       if (in_output) then 
         check_hierarchy = .true.
         write(*,*)'tag starts:position in output'
       endif  
       if (in_atom) then 
         check_hierarchy = .true.
         write(*,*)'tag starts:position in atom'
       endif  
       if (in_heatbath) then 
         check_hierarchy = .true.
         write(*,*)'tag starts:position in heatbath'
       endif  
       if (.not. check_hierarchy) then 
         write(*,*)'incorrect_hierarchy for position'
         incorrect_hierarchy = .true.
       endif  
!!------------multiple-defined tag:velocity ------
!!------------  (defined in atom and heatbath tags) ------
      case('velocity')
       check_hierarchy = .false.
       if (in_atom) then 
         check_hierarchy = .true.
         write(*,*)'tag starts:velocity in atom'
       endif  
       if (in_heatbath) then 
         check_hierarchy = .true.
         write(*,*)'tag starts:velocity in heatbath'
       endif  
       if (.not. check_hierarchy) then 
         write(*,*)'incorrect_hierarchy for velocity'
         incorrect_hierarchy = .true.
       endif  
!!--------------------------------------------------
!!------------tags in the 'generate' file-----------
!!--------------------------------------------------
      case('generate')
       write(*,*)'tag starts:generate'
!!--------------------------------------------------
!!-------------------unknown tag--------------------
!!--------------------------------------------------
      case default
       write(*,*)'UNKNOWN TAG!!! :',trim(name)
       unknown_tag = .true.
    end select
!
  end subroutine begin_element
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pcdata_chunk(chunk)
!
    implicit none
    character(len=*), intent(in) :: chunk
    integer :: status, ierr
!
    if (trim(chunk) /= '') then
      write(*,*)'tag content=',trim(chunk)
    endif  
!
  end subroutine pcdata_chunk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine end_element(name)
     character(len=*), intent(in) :: name
!
     select case(name)
       case('heatbath')
        in_heatbath = .false.
        write(*,*)'tag ends:heatbath'
       case('atom')
        in_atom = .false.
        write(*,*)'tag ends:atom'
       case('output')
        in_atom = .false.
        write(*,*)'tag ends:outputm'
       case default
        write(*,*)'tag ends:',trim(name)
     end select   
!
  end subroutine end_element
!

end module M_xml_checker_hander

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program xml_checker
    use flib_sax,      only : open_xmlfile, xml_parse !(routine)
    use flib_sax,      only : xml_t                   !(type) 
    use M_xml_checker_hander, only : begin_element, end_element, pcdata_chunk !(routine)
    use M_xml_checker_hander, only : unknown_tag         !(CHANGED)
    use M_xml_checker_hander, only : incorrect_hierarchy !(CHANGED)
    use M_xml_checker_hander, only : in_heatbath         !(CHANGED)
    use M_xml_checker_hander, only : in_atom             !(CHANGED)
    use M_xml_checker_hander, only : in_output           !(CHANGED)
!
    implicit none
    character(len=64) :: filename
    integer       :: iargc  
    type(xml_t) :: fxml ! XML file object (opaque)
    integer :: ierr
    logical :: tag_is_ok
!    
    write(*,'(a)') 'XML tag checker'
!
    if( iargc() /= 1 ) then
       write(*,'(a)') "# Usage of this utility program"
       write(*,'(a)') "#   xml-check test.xml"
       stop
    end if
!    
    call getarg(1,filename)
!    
    call open_xmlfile(trim(filename),fxml,iostat=ierr)
    if (ierr /=0) then
      write(*,*)'ERROR(structure_load_sax):file not found;',trim(filename)
    else
      write(*,*)'scanned file name=',trim(filename)
    endif
!
    tag_is_ok           = .true.
    unknown_tag         = .false.
    incorrect_hierarchy = .false.
!
    in_heatbath = .false.
    in_atom     = .false.
    in_output   = .false.
!
    write(*,*)'---------------Tag check starts--------------------'
!
    call xml_parse(fxml, begin_element_handler=begin_element, &
&                      end_element_handler=end_element, &
&                      pcdata_chunk_handler=pcdata_chunk )
!
    write(*,*)'---------------Tag check ends--------------------'
!
    write(*,*)'scanned file name=',trim(filename)
!
    if (unknown_tag) then
      write(*,*)'RESULT:ERROR:unknown tag is found!'
      tag_is_ok = .false.
    endif   
!
    if (incorrect_hierarchy) then
      write(*,*)'RESULT:ERROR: incorrect hierarchy'
      tag_is_ok = .false.
    endif   
!
    if (tag_is_ok) then
      write(*,*)'RESULT: OK !'
    endif   
!
end program xml_checker

