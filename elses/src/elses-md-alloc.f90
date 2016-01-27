!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_alloc
!
  use M_config,           only: config    !(unchanged)
!
  implicit none
  integer  :: i_verbose, log_unit
!
  private
  public array_alloc_for_md
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine array_alloc_for_md(n)
!
    use elses_mod_tx,       only : tx,  ty,  tz,  jsei
    use elses_mod_txp,      only : txp, typ, tzp
    use elses_mod_foi,      only : foi
    use elses_mod_foiold,   only : foiold
    use elses_mod_iflag,    only : iflag
    implicit none
    integer, intent(in) :: n
    integer             :: ierr
    logical             :: flag_for_alloc
!
    i_verbose=config%option%verbose
    log_unit  = config%calc%distributed%log_unit
!
    allocate (jsei(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:jsei'
    if (ierr /= 0) then
      write(*,*)'Alloc error:jsei'
      stop
    endif
!
    allocate (tx(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:tx'
    if (ierr /= 0) then
      write(*,*)'Alloc error:tx'
      stop
    endif
!
    allocate (ty(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:ty'
    if (ierr /= 0) then
      write(*,*)'Alloc error:ty'
      stop
    endif
!
    allocate (tz(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:tz'
    if (ierr /= 0) then
      write(*,*)'Alloc error:tz'
      stop
    endif
!
    allocate (txp(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:txp'
    if (ierr /= 0) then
      write(*,*)'Alloc error:txp'
      stop
    endif
!
    allocate (typ(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:typ'
    if (ierr /= 0) then
      write(*,*)'Alloc error:typ'
      stop
    endif
!
    allocate (tzp(n),stat=ierr)
    if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:tzp'
    if (ierr /= 0) then
      write(*,*)'Alloc error:tzp'
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    flag_for_alloc = .true.
    if (config%calc%distributed%set) then
      if (trim(config%calc%calc_force_mode) == 'off') flag_for_alloc = .false.
    endif
    if ( flag_for_alloc ) then
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:foi'
      allocate (foi(n,3),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc error:foi'
        stop
      endif
      foi(:,:)=0.0d0
    else  
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:foi:IGNORED'
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    flag_for_alloc = .true.
    if (config%calc%distributed%set) then
      if ( trim(config%calc%mode) /= "dynamics" ) flag_for_alloc = .false.
    endif
    if ( flag_for_alloc ) then
      allocate (foiold(n,3),stat=ierr)
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:foiold'
      if (ierr /= 0) then
        write(*,*)'Alloc error:foiold'
        stop
      endif
      foiold(:,:)=0.0d0
    else  
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:foiold:IGNORED'
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    flag_for_alloc = .true.
    if (config%calc%distributed%set) then
      flag_for_alloc = .false.
      if ( trim(config%calc%mode) /= "dynamics" ) flag_for_alloc = .true.
      if ( trim(config%calc%mode) /= "optimization" ) flag_for_alloc = .true.
    endif
    if ( flag_for_alloc ) then
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:iflag'
      allocate (iflag(n),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc error:iflag'
        stop
      endif
      iflag(:)=1
    else
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:allocation:iflag:IGNORED'
    endif  
!
  end subroutine array_alloc_for_md
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_md_alloc



