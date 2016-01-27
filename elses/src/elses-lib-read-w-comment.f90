!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_read_w_comment
!
   implicit none
!
   private
!
   public :: read_w_comment
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read a file with comment
!
  subroutine read_w_comment(fd, comment_chara, result_chara, end_of_file)
   implicit none
   integer,             intent(in)  :: fd
   character(len=*),    intent(in)  :: comment_chara
   character(len=*),  intent(inout) :: result_chara
   integer, parameter               :: max_comment_line=1024
   logical,             intent(out) :: end_of_file
   character(len=10240)             :: chara_wrk
   integer  :: j, ierr
!
   result_chara = ''
   end_of_file  = .false.
!
   do j = 1, max_comment_line
     if (j == max_comment_line) then
       write(*,*)'ERROR(read_w_comment) j=', j
       stop
     endif
     read(fd, '(a)', iostat=ierr) chara_wrk 
     if (ierr /= 0) then 
       end_of_file = .true.
       chara_wrk   = ''
       exit
     endif
     if (index(chara_wrk,comment_chara) == 1) then 
       cycle
     else
       exit
     endif
   enddo
!
   if ( len_trim(chara_wrk) > len(result_chara) ) then
     write(*,*)'ERROR(read_w_comment) chara_wrk=', trim(chara_wrk)
     stop
   endif
!
   result_chara=trim(chara_wrk)
!
  end subroutine read_w_comment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_read_w_comment
