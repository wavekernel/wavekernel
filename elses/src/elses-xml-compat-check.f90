!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_xml_compat_chk
!
  use M_config,           only: config    !(unchanged)
!
  implicit none
  integer  :: i_verbose, log_unit
!
  private
  public check_xml_compat
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine check_xml_compat
!
    implicit none
    integer i_v, lu
    integer ierr
!
    i_v = config%option%verbose
    lu  = config%calc%distributed%log_unit
!
    if (i_v >= 1) then
      if (lu > 0) write(lu,'(a)') '@@ check_xml_compat:check XML file for compatibility'
    endif
!
    if (config%system%structure%use_matom) then
      ierr=1
      if (trim(config%calc%mode) == "cell_change_only" ) ierr=0
      if (trim(config%calc%mode) == "conversion" ) ierr=0
      if (ierr /= 0) then
        write(*,*)'ERROR:(check_xml_compat)incompatible setting'
        write(*,*)'  use_matom  = ', config%system%structure%use_matom
        write(*,*)'  calc%mode  = ', trim(config%calc%mode)
        stop
      endif
    endif
!
    if (config%calc%distributed%set) then
      call check_xml_compat_dst
    else
      call check_xml_compat_non_dst
    endif  
!
!   stop 'stop manually'
!
  end subroutine check_xml_compat
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine check_xml_compat_dst
    implicit none
    integer :: ierr
!
    if (log_unit > 0) write(log_unit,'(a)') '@@ check_xml_compat in DST workflow'
!
    if (trim(config%system%structure%velement(1)%quantum%type) /= 'geno') then
      write(*,*) 'ERROR(XML compatibility) : system =', trim(config%system%structure%velement(1)%quantum%type) 
      stop
    endif
!
    ierr=1
    if (trim(config%calc%solver%scheme) == 'gKrylov_A') ierr = 0
    if (trim(config%calc%solver%scheme) == 'eigen_mpi') ierr = 0
!
    if (ierr == 1) then
      write(*,*) 'ERROR(XML compatibility) : solver =', trim(config%calc%solver%scheme)
      stop
    endif
!
  end subroutine check_xml_compat_dst
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine check_xml_compat_non_dst
    implicit none
!
    if (log_unit > 0) write(log_unit,'(a)') '@@ check_xml_compat in non-DST workflow'
!
    if (trim(config%calc%calc_force_mode) == "off") then
      write(*,*) 'ERROR(XML compatibility) : calc_force_mode =', trim(config%calc%calc_force_mode)
      stop
    endif   
!
  end subroutine check_xml_compat_non_dst
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module M_xml_compat_chk


