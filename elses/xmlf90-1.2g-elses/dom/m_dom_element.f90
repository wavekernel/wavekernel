module m_dom_element

use m_dom_types
use m_dom_namednodemap
use m_dom_nodelist
use m_dom_attribute
use m_dom_document
use m_dom_debug
use m_dom_node
use m_strings

private

  !-------------------------------------------------------   
  ! METHODS FOR ELEMENT NODES
  !-------------------------------------------------------   
  public :: getTagName
  public :: getElementsByTagName
  public :: getFirstElementByTagName ! Mizuho-IR
  public :: getAttribute
  public :: getAttributeNode
  public :: setAttribute
  public :: setAttributeNode
  public :: removeAttribute
  public :: normalize       !--- combines adjacent text nodes ---!

CONTAINS

  !-----------------------------------------------------------
  !  METHODS FOR ELEMENT NODES
  !-----------------------------------------------------------
  function getTagName(element)

    type(fnode), intent(in) :: element   
    type(string)            :: getTagName

    if (element % nodeType == ELEMENT_NODE) then
       getTagName = element % nodeName 
    else
       getTagName = ''
    endif

  end function getTagName

  !-----------------------------------------------------------
  function getElementsByTagName(element, tag) result(nodelist)
    type(fnode), pointer         :: element
    character(len=*), intent(in) :: tag
    type(fnodeList), pointer     :: nodelist 

!    type(fnode), pointer        :: np ! Mizuho-IR

    nodelist => null()

!    np => element ! Mizuho-IR
    if (dom_debug) print *, "Going into search for tag: ", trim(tag)

!    call search(np) ! Mizuho-IR
    call search(element,0) ! Mizuho-IR

    CONTAINS

!    recursive subroutine search(np) ! Mizuho-IR
    recursive subroutine search(mizuho,level) ! Mizuho-IR
    type(fnode), pointer        :: mizuho ! Mizuho-IR
    type(fnode), pointer        :: np
    integer                     :: level ! Mizuho-IR
!    type(string)                :: name ! Mizuho-IR

    np => mizuho ! Mizuho-IR
    !
    ! Could replace the calls to helper methods by direct lookups of node 
    ! components to make it faster.
    ! 
    do
       if (.not. associated(np)) exit
       select case(np%nodeType)

          case(DOCUMENT_NODE) 
             ! special case ... search its children 
!             if (hasChildNodes(np)) call search(getFirstChild(np)) ! Mizuho-IR
             if (hasChildNodes(np)) call search(np%firstChild,level+1) ! Mizuho-IR
             ! will exit for lack of siblings
          case(ELEMENT_NODE)

!             name = getNodeName(np) ! Mizuho-IR

             if (dom_debug) print *, "exploring node: ", char(np%nodeName)
             if ((tag == "*") .or. (tag == np%nodeName)) then
                call append(nodelist,np)
                if (dom_debug) print *, "found match ", nodelist%length
             endif
!             if (hasChildNodes(np)) call search(getFirstChild(np)) ! Mizuho-IR
             if ( associated(np%firstChild) ) call search(np%firstChild,level+1) ! Mizuho-IR

          case default
             
             ! do nothing

        end select

        if (associated(np,element)) exit  ! no siblings of element...
!        np => getNextSibling(np) ! Mizuho-IR
        if( level > 0 ) then ! Mizuho-IR
           np => np%nextSibling ! Mizuho-IR
        endif ! Mizuho-IR

     enddo
    
    end subroutine search

  end function getElementsByTagName

  !-----------------------------------------------------------

!----------------------------------------------------------------------
! Mizuho-IR

  function getFirstElementByTagName(element, tag) result(node)
    type(fnode), pointer         :: element
    character(len=*), intent(in) :: tag
    type(fnodeList), pointer     :: nodelist 
    type(fnode), pointer         :: node

    nodelist => getElementsByTagName( element, tag )
    if( .not. associated(nodelist) ) then
       node => null()
    else
       node => item(nodelist,0)
    end if
  end function getFirstElementByTagName

!----------------------------------------------------------------------

  function getAttribute(element, name)
    
    type(fnode), intent(in) :: element
    character(len=*), intent(in) :: name
    type(string)                 :: getAttribute

    type(fnode), pointer :: nn

    getAttribute = ""  ! as per specs, if not found
    if (element % nodeType /= ELEMENT_NODE) RETURN
    nn => getNamedItem(element%attributes,name)
    if (.not. associated(nn)) RETURN
    
    getAttribute = nn%nodeValue

        
  end function getAttribute

  !-----------------------------------------------------------

  function getAttributeNode(element, name)
    
    type(fnode), intent(in) :: element
    type(fnode), pointer    :: getAttributeNode
    character(len=*), intent(in) :: name

    getAttributeNode => null()     ! as per specs, if not found
    if (element % nodeType /= ELEMENT_NODE) RETURN
    getAttributeNode => getNamedItem(element%attributes,name)

  end function getAttributeNode
  
  !-----------------------------------------------------------

  subroutine setAttributeNode(element, newattr)
    type(fnode), pointer :: element
    type(fnode), pointer :: newattr

    type(fnode), pointer :: dummy

    if (element % nodeType /= ELEMENT_NODE) then
       if (dom_debug) print *, "not an element node in setAttributeNode..."
       RETURN
    endif

    dummy => setNamedItem(element%attributes,newattr)
     
  end subroutine setAttributeNode

!-------------------------------------------------------------------
  subroutine setAttribute(element, name, value)
    type(fnode), pointer :: element
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    type(fnode), pointer      :: newattr

    newattr => createAttribute(name)
    call setValue(newattr,value)
    call setAttributeNode(element,newattr)

  end subroutine setAttribute

  !-----------------------------------------------------------

  subroutine removeAttribute(element, name)
    type(fnode), pointer :: element
    character(len=*), intent(in) :: name

    type(fnode), pointer :: dummy

    if (element % nodeType /= ELEMENT_NODE) RETURN
    if (.not. associated(element%attributes)) RETURN

    dummy => removeNamedItem(element%attributes,name)
     
  end subroutine removeAttribute

  !-----------------------------------------------------------
  recursive subroutine normalize(element)
    type(fnode), pointer         :: element

    type(fnode), pointer        :: np, ghost
    logical                     :: first

    type(fnode), pointer        :: head

    first = .true.  ! next Text node will be first

    if (dom_debug) print *, "Normalizing: ", trim(element%nodeName)
    np => element%firstChild
    ! 
    do
       if (.not. associated(np)) exit
       select case(np%nodeType)

          case(TEXT_NODE) 
             if (first) then
                if (dom_debug) print *, "normalize: found first in chain"
                head => np
                first = .false.
!                np => getNextSibling(np) ! Mizuho-IR
                np => np%nextSibling ! Mizuho-IR
             else                    ! a contiguous text node
                if (dom_debug) print *, "normalize: found second in chain"
!                head%nodeValue = head%nodeValue // np%nodeValue ! Mizuho-IR
                call append_to_string( head%nodeValue, np%nodeValue ) ! Mizuho-IR

                head%nextSibling => np%nextSibling
                if (associated(np,np%parentNode%lastChild)) then
                   np%parentNode%lastChild => head
                   head%nextSibling => null()
                else
                   np%nextSibling%previousSibling => head
                endif
                ghost => np
!                np => getNextSibling(np) ! Mizuho-IR
                np => np%nextSibling ! Mizuho-IR
                call destroyNode(ghost)
             endif

          case(ELEMENT_NODE)

             first = .true.
             if (dom_debug) print *, "element sibling: ", trim(np%nodeName)
             if (hasChildNodes(np)) call normalize(np)
!             np => getNextSibling(np) ! Mizuho-IR
             np => np%nextSibling ! Mizuho-IR

          case default
             
             ! do nothing, just mark that we break the chain of text nodes
             if (dom_debug) print *, "other sibling: ", trim(np%nodeName)
             first = .true.
!             np => getNextSibling(np) ! Mizuho-IR
             np => np%nextSibling ! Mizuho-IR

        end select

     enddo

    end subroutine normalize


end module m_dom_element
