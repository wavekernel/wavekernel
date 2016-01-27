!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_pair_list_for_cohp
!
  use M_qm_domain,          only : i_verbose !(unchanged)
  use elses_mod_phys_const, only : ev4au     !(parameter)
!
  character(len=*), parameter :: file_name='input_pair_list_for_cohp.txt'
!
  private
!
! Public routines
   public set_pair_num_for_cohp
   public set_pair_list_for_cohp
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_pair_num_for_cohp(num_pair_list)
!
    use elses_mod_file_io, only : vacant_unit
    use M_qm_domain,       only : njsd, noav
    implicit none
    integer, intent(out) :: num_pair_list
!
    integer :: iunit, ierr
    logical :: file_exist
    integer :: ipair, atm_index1, atm_index2
    integer :: line_count, add_line_count
    integer, parameter :: line_count_max = 100000000
    integer, parameter :: ict4h = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) write(*,*)'@@ set_pair_num_for_cohp'
!
    inquire (file=trim(file_name), exist=file_exist)
!
    if (file_exist .eqv. .false.) then
      if (i_verbose >= 1) write(*,*)'The file does not exist:',file_name
      num_pair_list=0  ! dummy data
      return
    endif
!
    iunit=vacant_unit()
!
    open(iunit, file=file_name, status='unknown')
!
!
    add_line_count=0
    do line_count=1,line_count_max
      if (line_count >= line_count_max) then
        write(*,*)'ERROR(set_pair_num_for_cohp):line_count=',line_count 
        stop 
      endif   
      read(iunit,*,iostat=ierr) atm_index1, atm_index2
      if (ierr /= 0) then
        num_pair_list=line_count-1
        exit
      endif  
      write(*,*)'pair list=',atm_index1, atm_index2
      if ((atm_index1 < 0) .or. (atm_index1 > noav)) then
        write(*,*)'ERROR(set_pair_num_for_cohp):atm_index1=',atm_index1
        stop
      endif   
      if ((atm_index2 < 0) .or. (atm_index2 > noav)) then
        write(*,*)'ERROR(set_pair_num_for_cohp):atm_index2=',atm_index2
        stop
      endif   
      if (atm_index2 == 0) then
        add_line_count=add_line_count+njsd(atm_index1,ict4h)-2 
!            ! Note : '-2' exclcudes the atom itself.
      endif   
    enddo   
!
    num_pair_list=num_pair_list+add_line_count
!
    if (i_verbose >= 1) then 
      write(*,*)' num_pair_list=',num_pair_list
    endif  
!
!   stop
!
    close(iunit)
!
   end subroutine set_pair_num_for_cohp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_pair_list_for_cohp(num_pair_list, pair_list)
!
    use elses_mod_file_io, only : vacant_unit
    use M_qm_domain,       only : njsd, noav, jsv4jsd
!
    implicit none
    integer, intent(in)    :: num_pair_list
    integer, intent(inout) :: pair_list(:,:)
!
    integer :: iunit, ierr
    logical :: file_exist
    integer :: ipair, atm_index1, atm_index2, atm_index2b
    integer :: num_pair_list_dummy
    integer :: line_count, jsd
    integer, parameter :: ict4h = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    inquire (file=trim(file_name), exist=file_exist)
!
    if (file_exist .eqv. .false.) then
      if (i_verbose >= 1) then 
        write(*,*)'ERROR(set_pair_list_for_cohp)The file does not exist:',file_name
        stop
      endif  
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    iunit=vacant_unit()
!
    open(iunit, file=file_name, status='unknown')
!
    if (i_verbose >= 1) then 
      write(*,*)' open file for pair list:',file_name
    endif  
!
!   read(iunit,*) num_pair_list_dummy
!
!   if (num_pair_list /= num_pair_list_dummy) then
!     write(*,*)'ERROR:num_pair_list=', num_pair_list, num_pair_list_dummy 
!     stop
!   endif
!
    if (i_verbose >= 1) then 
      write(*,*)' num_pair_list=',num_pair_list
    endif  
!
!   allocate(pair_list(2,num_pair_list), stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in alloc:pair_list'
!
    line_count=0
    do ipair=1,num_pair_list
      read(iunit,*,iostat=ierr) atm_index1, atm_index2
      if (ierr /= 0) then
        exit
      endif   
      if (i_verbose >= 1) then 
         write(*,*)' pair : =',ipair, atm_index1,atm_index2
      endif  
      if (atm_index1 == atm_index2) then
        write(*,*)'ERROR(set_pair_list_for_cohp):atm_index1,2=',atm_index1,atm_index2 
        stop
      endif   
      if (atm_index2 /= 0) then
        line_count=line_count+1
        if (line_count > size(pair_list,2)) then
          write(*,*)'ERROR(set_pair_list_for_cohp):line_count=',line_count
          stop
        endif   
        pair_list(1,line_count)=atm_index1
        pair_list(2,line_count)=atm_index2
      else
        do jsd=1,njsd(atm_index1,ict4h) 
          atm_index2b=jsv4jsd(jsd,atm_index1) 
          if (atm_index1 == atm_index2b) cycle
          line_count=line_count+1
          if (line_count > size(pair_list,2)) then
            write(*,*)'ERROR(set_pair_list_for_cohp):line_count=',line_count
            stop
          endif   
          pair_list(1,line_count)=atm_index1
          pair_list(2,line_count)=atm_index2b
        enddo
      endif   
    enddo 
!
!
    close(iunit)
!
   end subroutine set_pair_list_for_cohp
!
end module M_pair_list_for_cohp




