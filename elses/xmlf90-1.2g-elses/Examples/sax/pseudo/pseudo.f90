program pseudo
!
! Example driver for the construction of data structures from an xml document
!
use flib_sax
use m_pseudo      ! Defines begin_element, end_element, pcdata_chunk

integer      :: iostat
type(xml_t)  :: fxml

call open_xmlfile("pseudo.xml",fxml,iostat)
if (iostat /=0) stop "Cannot open file"

call xml_parse(fxml, &
                begin_element,end_element,pcdata_chunk,verbose=.false.)

end program pseudo














