!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_get_atom_param

  private 
  public :: get_orb_num
!
contains

 function get_orb_num(elem_name_in) result(orb_num)
   use M_config
   use M_qm_domain,          only : nos, nval
   implicit none
   character(len=*), intent(in)  :: elem_name_in
   character(len=64)             :: elem_name_wrk
   integer                       :: orb_num
   integer                       :: nss
!
   orb_num=-1 ! dummy value
   do nss=1,nos
     elem_name_wrk=trim(config%system%structure%velement(nss)%name)
     if (trim(elem_name_in) == trim(elem_name_wrk)) then
       orb_num=nval(nss)
     endif
   enddo
!
   if (orb_num == -1) then
     write(*,*)'ERROR(get_orb_num): elem_name=', trim(elem_name_in)
     stop
   endif
!
 end function get_orb_num

end module M_get_atom_param

