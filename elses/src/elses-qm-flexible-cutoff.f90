!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_flexible_cutoff
!
  private
  public flexible_cutoff_ini
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine flexible_cutoff_ini
    use M_config,    only : config !(unchanged) 
    implicit none
    integer              :: i_verbose, log_unit
!
    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit
!
    if (.not. config%calc%solver%flexible_cutoff%set) then
      if (I_verbose >=1) then
        if (log_unit > 0) then
          write(log_unit, '(a)') '@@ flexible_cutoff_ini .... is skipped'
        endif   
      endif   
      return
    endif   
!
    if (I_verbose >=1) then
      if (log_unit > 0) then
        write(log_unit, '(a)') '@@ flexible_cutoff_ini'
      endif   
    endif   
!
  end subroutine flexible_cutoff_ini
!
end module M_qm_flexible_cutoff



