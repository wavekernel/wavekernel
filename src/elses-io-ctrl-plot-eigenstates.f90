!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_io_ctrl_output_eigen
!
  use M_qm_domain,          only : i_verbose !(unchanged)
!
   private
!
! Public routines
   public set_filename_for_save_wfn
   public init_for_plot_wavefunction 
!        ! A wrapper routine for OLD interface (used for compatibility)
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_filename_for_save_wfn(filename_wrk)
!
    use M_config,  only : config
    use elses_arr_eig_leg, only : atmp
    use elses_mod_md_dat,  only : final_iteration !(unchanged)
    use elses_mod_md_dat,  only : itemd           !(unchanged)
    use M_qm_domain,       only : c_system        !(unchanged)
    implicit none
    character(len=100),  intent(out) :: filename_wrk
    character(len=32)                :: chara_wrk
    integer                          :: lenf
    integer                          :: log_unit
    logical                          :: save_now
!
    log_unit=config%calc%distributed%log_unit
!
    filename_wrk=""   ! dummy setting
!
    if (log_unit > 0) write(log_unit,*)' @@ init_for_plot_wavefunction:GENO compatible:c_system=',trim(c_system)
    if (log_unit > 0) write(log_unit,*)' @@ init_for_plot_wavefunction:interval=',config%output%wavefunction%interval
!
    if ( .not. config%output%wavefunction%set ) return
!

!
!   write(*,*)' wavefunction%interval=', config%output%wavefunction%interval
!
    if (config%output%wavefunction%interval < 0 ) then 
      write(*,*)' ERROR!:wavefunction%interval=', config%output%wavefunction%interval
      stop
    endif
!   
!
    save_now = .false.
    if (final_iteration) then 
      save_now = .true.
    else   
      if (config%output%wavefunction%interval > 0 ) then 
       if ( mod(itemd-1,config%output%wavefunction%interval) ==  0 ) save_now = .true.
      endif
    endif     
!
    if (.not. save_now) return
!
    if (.not. allocated(atmp)) return
!   
    filename_wrk=config%output%wavefunction%filename
    if (trim(filename_wrk) == "") then
      write(*,*)'ERROR(init_for_plot_wavefunction)'
      stop
    endif   
!
    if (config%output%wavefunction%interval > 0 ) then 
      lenf=len_trim(filename_wrk) 
      if (filename_wrk(lenf-3:lenf)=='.txt') then
        filename_wrk=filename_wrk(1:lenf-4) 
      endif   
      write(chara_wrk, '(i8.8)') config%system%structure%mdstep
      filename_wrk=trim(filename_wrk)//'_'//trim(chara_wrk)//'.txt'
    endif   
!
    if (i_verbose >= 1) then 
      if (log_unit > 0) write(log_unit,*)' save wfn: itemd, filename=',itemd, filename_wrk
    endif  
!
!   stop
!
   end subroutine set_filename_for_save_wfn
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine init_for_plot_wavefunction(imode,filename_wrk)
!
    use M_config,  only : config
    use elses_arr_eig_leg, only : atmp
    implicit none
    integer,             intent(out) :: imode
    character(len=100),  intent(out) :: filename_wrk
!
     call set_filename_for_save_wfn(filename_wrk)
!
     if (trim(filename_wrk) /= '') then 
       imode=1
     else
       imode=0 
     endif  
!
   end subroutine init_for_plot_wavefunction
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module M_io_ctrl_output_eigen
