module m_dom_utils

  use m_dom_types
  use m_dom_element
  use m_dom_document
  use m_dom_node
  use m_dom_namednodemap
  use m_dom_debug
  use m_strings

!  use flib_wxml ! Mizuho-IR

  public :: dumpTree
  public :: xmlize

  private

CONTAINS

  subroutine dumpTree(startNode)

    type(fnode), pointer :: startNode   

    character(len=50) :: indent = " "
    integer           :: indent_level
!    type(string)      :: s ! Mizuho-IR

    indent_level = 0

    call dump2(startNode)

  contains

    recursive subroutine dump2(input)
      type(fnode), pointer :: input
      type(fnode), pointer :: temp     
      temp => input
      do while(associated(temp))
!         s = getNodeName(temp) ! Mizuho-IR
         write(*,'(3a,i3)') indent(1:indent_level), &
!                        char(s), " of type ", &
                        char(temp%nodeName), " of type ", & ! Mizuho-IR
                        getNodeType(temp)
         if (hasChildNodes(temp)) then
            indent_level = indent_level + 3
!            call dump2(getFirstChild(temp)) ! Mizuho-IR
            call dump2(temp%firstChild) ! Mizuho-IR
            indent_level = indent_level - 3
         endif
!         temp => getNextSibling(temp) ! Mizuho-IR
         temp => temp%nextSibling ! Mizuho-IR
      enddo

    end subroutine dump2

  end subroutine dumpTree


!----------------------------------------------------------------

  subroutine xmlize(startNode,fname)
    type(fnode), pointer :: startNode   
    character(len=*), intent(in) :: fname
    !SY-from here
    print *, "Warning xmlize is disabled"
    print *, associated(startNode), fname
    !SY-to here
! Mizuho-IR
!!$
!!$    type(xmlf_t)  :: xf
!!$    type(string)  :: s, sv, sn       ! to avoid memory leaks
!!$
!!$    call xml_OpenFile(fname,xf)
!!$    call dump_xml(startNode)
!!$    call xml_Close(xf)
!!$
!!$  contains
!!$
!!$    recursive subroutine dump_xml(input)
!!$      type(fnode), pointer         :: input
!!$!
!!$!     Just this node and its descendants, no siblings. 
!!$!     Of course, the document node has only children...
!!$!
!!$      type(fnode), pointer         :: node, attr
!!$      type(fnamedNodeMap), pointer :: attr_map
!!$      integer  ::  i
!!$
!!$      node => input
!!$      do
!!$         if (.not. associated(node)) exit
!!$         select case (getNodeType(node))
!!$
!!$          case (DOCUMENT_NODE)
!!$
!!$             call xml_AddXMLDeclaration(xf)
!!$!             if (hasChildNodes(node)) call dump_xml(getFirstChild(node)) ! Mizuho-IR
!!$             if (hasChildNodes(node)) call dump_xml(node%firstChild) ! Mizuho-IR
!!$
!!$          case (ELEMENT_NODE)
!!$
!!$             s = getNodeName(node)
!!$             call xml_NewElement(xf,char(s))
!!$             attr_map => getAttributes(node)
!!$             do i = 0, getLength(attr_map) - 1
!!$                attr => item(attr_map,i)
!!$                sn = getNodeName(attr)
!!$                sv = getNodeValue(attr)
!!$                call xml_AddAttribute(xf, char(sn),char(sv))
!!$             enddo
!!$!             if (hasChildNodes(node)) call dump_xml(getFirstChild(node)) ! Mizuho-IR
!!$             if (hasChildNodes(node)) call dump_xml(node%firstChild) ! Mizuho-IR
!!$             s = getNodeName(node)
!!$             call xml_EndElement(xf,char(s))
!!$
!!$          case (TEXT_NODE)
!!$             
!!$             s = getNodeValue(node)
!!$             call xml_AddPcdata(xf,char(s))
!!$
!!$          case (CDATA_SECTION_NODE)
!!$             
!!$             s = getNodeValue(node)
!!$             call xml_AddCdataSection(xf,char(s))
!!$
!!$          case (COMMENT_NODE)
!!$             
!!$             s = getNodeValue(node)
!!$             call xml_AddComment(xf,char(s))
!!$
!!$        end select
!!$        if (associated(node,StartNode)) exit  ! In case we request the 
!!$                                              ! dumping of a single element,
!!$                                              ! do not do siblings
!!$!        node => getNextSibling(node) ! Mizuho-IR
!!$        node => node%nextSibling ! Mizuho-IR
!!$     enddo
!!$
!!$    end subroutine dump_xml
!!$
  end subroutine xmlize

end module m_dom_utils
