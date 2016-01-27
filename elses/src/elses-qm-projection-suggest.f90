!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_proj_suggest_radius
!
  use M_io_dst_write_log,  only : log_unit !(unchanged)
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
!
   implicit none
!
   private
!
! Public routine
   public :: suggest_proj_radius
   public :: suggest_proj_radius_wrap
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine suggest_proj_radius(noa_proj, mode, proj_radius_suggested)
!
    use M_config !(unchanged)
    use M_qm_domain,        only : noav, ax, ay, az    !(unchanged)
    use elses_mod_sel_sys,  only : r_cut_book          !(unchanged)
    use M_lib_phys_const,   only : pi                  !(unchanged)
!
    implicit none
    integer,                intent(in) :: noa_proj
    character(len=32),      intent(out):: mode
    real(DOUBLE_PRECISION), intent(out):: proj_radius_suggested
!
    real(DOUBLE_PRECISION) :: volume
    real(DOUBLE_PRECISION), parameter :: alpha=0.7d0 ! a factor
    real(DOUBLE_PRECISION), parameter :: proj_radius_min=2.0d0 ! minimum value (for error checking)
!
    mode=config%calc%solver%mode_for_suggest_projection
!
!   if (i_verbose >=0) then
!     if (i_verbose >=1) write(*,*)'@@ suggest_proj_radius: mode=',trim(mode) 
!   endif  
!
    if (log_unit > 0) then
      if (i_verbose >=1) write(log_unit,*)'@@ suggest_proj_radius: mode=',trim(mode) 
    endif
!   
    if ( trim(mode) == "default") mode="from_density"
!
    if ( trim(mode) == "upto_v0.03.08" ) then
      proj_radius_suggested=r_cut_book*2.0d0 
      return
    endif
!   
    if ( trim(mode) == "from_density" ) then
      volume=ax*ay*az
      proj_radius_suggested= alpha*      &
&                (volume*dble(noa_proj)/dble(noav)*3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)
      if ( proj_radius_suggested < proj_radius_min ) then
        write(*,*)'ERROR(suggest_proj_radius):proj_radius_suggested=',proj_radius_suggested 
        write(*,*)'ERROR(suggest_proj_radius):proj_radius_min      =',proj_radius_min
        write(*,*)'ERROR(suggest_proj_radius):noa_proj             =',noa_proj
        stop
      endif
      return
    endif
!
    write(*,*)'ERROR(suggest_proj_radius):mode=',trim(mode)
    stop
!
  end subroutine suggest_proj_radius
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine suggest_proj_radius_wrap(proj_radius_suggested)
    use M_config                                !(unchanged)
    use elses_param_ctl_kr, only : noak_min_def !(unchanged)
    implicit none
    real(DOUBLE_PRECISION), intent(out):: proj_radius_suggested
!
    integer           :: noa_proj
    character(len=32) :: suggestion_mode
!
    suggestion_mode=trim(config%calc%solver%mode_for_suggest_projection)
    noa_proj=noak_min_def
!
    call suggest_proj_radius(noa_proj, suggestion_mode, proj_radius_suggested)
!
  end subroutine suggest_proj_radius_wrap
!
end module M_qm_proj_suggest_radius

