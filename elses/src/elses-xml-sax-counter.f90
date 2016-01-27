!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_sax_counter
!
  use M_config, only : config  !(CHANGED)
  use elses_xml_misc, only : XML_TO_AU !(function)
  implicit none
  logical, parameter :: tag_dump = .false. ! true only for debugging 
!
  logical :: in_structure_first
  logical :: in_structure, in_atom 
!
  integer :: counter_atom
  integer :: counter_element
!
  integer :: j
  logical :: booked
!
  integer, parameter :: max_element_number = 100
  character(len=16)  :: element_store(max_element_number)
!
  private
  public :: in_structure_first
  public :: begin_element
  public :: pcdata_chunk
  public :: end_element
!
  public :: counter_atom
  public :: counter_element

contains
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
!
    select case(name)
      case('structure')
       in_structure = .true.
       if (.not. in_structure_first) then
         write(*,*)'ERROR(XML-SAX):not supported now:in_structure_first=',in_structure_first
         stop
       endif   
!
      case('atom')
        in_atom = .true. 
        counter_atom=counter_atom+1
        if (tag_dump) write(*,*)'atom tag (c)'
        call get_value(attributes,"element", chara_wrk, status)
        if (trim(chara_wrk) == '') then
          write(*,*)'No element attribute in atom tag' 
          stop
        endif   
        if (counter_element > max_element_number) then
          write(*,*)'ERROR(begin_element): counter_element=', counter_element
          stop
        endif
        booked = .false. 
        if (counter_element > 0) then 
          do j=1, counter_element
            if (trim(chara_wrk) == element_store(j)) then
              booked = .true. 
              exit
            endif   
          enddo   
        endif  
        if (.not. booked) then
          counter_element=counter_element+1
          element_store(counter_element)=trim(chara_wrk)
          if (tag_dump) write(*,*)'element =',counter_element, element_store(counter_element)
        endif   
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
!
!      Do nothing
!
  end subroutine pcdata_chunk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine end_element(name)
    character(len=*), intent(in) :: name
!
    select case(name)
      case('structure') 
        in_structure= .false. 
     case('atom')
       in_atom = .false. 
    end select   
!
  end subroutine end_element
!
end module M_sax_counter



