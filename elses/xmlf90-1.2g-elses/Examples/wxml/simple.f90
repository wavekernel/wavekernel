program simple

use flib_wxml

type(xmlf_t) :: xf

integer :: age = 34
real, dimension(20)  :: x
real, dimension(20,20)  :: y

call xml_OpenFile("simple.xml",xf, indent=.true.)

call xml_AddXMLDeclaration(xf,"UTF-8")
call xml_NewElement(xf,"john")
call xml_AddAttribute(xf,"age",str(age))
call xml_NewElement(xf,"peter")
call xml_NewElement(xf,"tim")
call xml_AddAttribute(xf,"age","37")
call xml_AddAttribute(xf,"weight",str(123.45,"(f7.3)"))
call xml_AddAttribute(xf,"cholesterol",str(167.0,format="(f8.0)"))
call xml_EndElement(xf,"tim")
call xml_AddPcdata(xf,"Ping-pong")
call xml_AddPcdata(xf,"champion", line_feed=.false.)
call xml_AddPcdata(xf," in 2004", space=.false., line_feed=.false.)
call xml_NewElement(xf,"data")
call xml_AddAttribute(xf,"units","eV")
call random_number(x)
call random_number(y)
call xml_AddArray(xf,x)
call xml_AddArray(xf,reshape(y,(/ 400 /)))
call xml_EndElement(xf,"data")
call xml_EndElement(xf,"peter")
call xml_EndElement(xf,"john")

call xml_Close(xf)

end program simple
