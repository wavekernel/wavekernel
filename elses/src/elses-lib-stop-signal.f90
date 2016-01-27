!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_stop_signal
!
  implicit none
!
  private
  public :: detect_stop_signal
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine detect_stop_signal(stop_signal, elapse_time)
!
    use M_config,          only : config !(unchanged)
    use elses_mod_file_io,    only : vacant_unit
    implicit none
    real(8), intent(in)  :: elapse_time
    integer, intent(out) :: stop_signal
    integer              :: unit_num, ierr
    character(len=*), parameter :: filename="00_stop_signal.txt"
    integer              :: tmp_variable
    logical              :: file_exists
    real(8)              :: time_limit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Detect the stop signal from the time limit  
!
    time_limit=config%calc%limit%time
    if (time_limit > -1.0d-10) then
      if (elapse_time > time_limit) then 
        stop_signal=1
        write(*,*)'Warning:The calculation will stop, because of the time limit' 
        write(*,*)' elapse time, limit (sec) =', elapse_time, time_limit
        return
      endif  
    endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Detect the stop signal from the file
!
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
      stop_signal=0
      return
    endif
!
    unit_num=vacant_unit()
!
    open(unit_num, file=filename)
!
    read(unit_num,*,iostat=ierr) tmp_variable
!
    if (ierr /= 0) then 
      stop_signal=0
    else
      stop_signal=tmp_variable
      write(*,*)'Warning:The calculation will stop, because of the presentce of the file'
      write(*,*)'  the file name=', filename
    endif   
!
    close(unit_num)
!
  end subroutine detect_stop_signal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_stop_signal



