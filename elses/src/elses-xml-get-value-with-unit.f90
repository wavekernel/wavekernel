module M_xml_get_value_with_unit
  use flib_dom
  use M_lib_phys_const, only : convert_unit ! routine
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
  private
  public xml_get_value_w_unit
  public xml_get_value_w_unit_detail
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xml_get_value_w_unit(node, tag_name, value_au,comment)
    implicit none
    type(fnode), pointer                   :: node
    character(len=*),       intent(in)     :: tag_name
    real(DOUBLE_PRECISION), intent(out)    :: value_au
    character(len=*),       intent(inout)  :: comment
!
    type(fnode), pointer                :: node_wrk
    character(len=256)                  :: value_chara, unit_name
    logical, parameter                  :: debug_mode = .true.
!
    real(DOUBLE_PRECISION)      :: value_org
    character(len=256)          :: unit_org
!
    value_org = 0.0d0 ! dummy 
    unit_org  = ''
    call xml_get_value_w_unit_detail(node, tag_name, value_au, value_org, unit_org, comment)
!
  end subroutine xml_get_value_w_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xml_get_value_w_unit_detail(node, tag_name, value_au, value_org, unit_org, comment)
    implicit none
    type(fnode), pointer                   :: node
    character(len=*),       intent(in)     :: tag_name
    real(DOUBLE_PRECISION), intent(out)    :: value_au
    real(DOUBLE_PRECISION), intent(out)    :: value_org
    character(len=*),       intent(inout)  :: unit_org
    character(len=*),       intent(inout)  :: comment
!
    type(fnode), pointer                :: node_wrk
    character(len=256)                  :: value_chara, unit_name
    logical, parameter                  :: debug_mode = .false.
!
    comment=''
!
    node_wrk => getFirstElementByTagName(node,tag_name)
    if ( .not. associated(node_wrk) ) then
      comment='error_no_tag'
      return
    endif
!
    value_chara = getChildValue(node_wrk) 
    unit_name   = getAttribute(node_wrk,"unit")
    if (debug_mode) write(*,*) 'value_chara, unit=', trim(value_chara), ' ', trim(unit_name)
    if ( unit_name == "" ) then
      write(*,*)'ERROR(xml_get_value_with_unit):no unit is specified' 
      write(*,*) 'value_chara = ', trim(value_chara)
      stop
    endif   
!
    unit_org=unit_name                                 ! the original unit  
    read (unit=value_chara, fmt=*) value_org           ! value in the original unit
    value_au=value_org*convert_unit('from',unit_name)  ! value in au
!
  end subroutine xml_get_value_w_unit_detail
!
end module M_xml_get_value_with_unit
