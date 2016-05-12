!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_output_levels
!
! use M_io_dst_write_log,  only : log_unit
  implicit none
! integer, allocatable :: called_first(:)
!
   private
!
! Public routines
   public qm_output_levels
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Output for eigen levels
!
  subroutine qm_output_levels(filename_wrk, unit_conv, append_mode, format_mode)
!
    use M_qm_domain ,  only : i_verbose    !(unchanged)
    use M_config,            only : config !(unchanged)
    use elses_arr_eig_leg,    only: eig2, f_occ, atmp !(unchanged)
    use elses_mod_file_io, only : vacant_unit !(function)
    use M_qm_domain ,  only : total_electron_number  !(unchanged)
    use M_output_participation, only : gene_participation !(routine)
!
    implicit none
    character(len=*), intent(in) :: filename_wrk
    logical,          intent(in) :: append_mode
    real(8),          intent(in) :: unit_conv
    character(len=*), intent(in) :: format_mode
!
    integer      :: jj, ierr, iunit, step_count, lu
    character(len=32) :: format_mode_wrk
    logical           :: part_ratio_is_on
    real(kind=8), allocatable    :: part_ratio(:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    lu = config%calc%distributed%log_unit
    part_ratio_is_on = .false.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Show the title (non-essential)
!
    if (i_verbose >= 1) then
      if (lu > 0) then 
        write(lu,*)'@@ qm_output_levels'
        write(lu,*)'  append_mode = ',append_mode
        write(lu,*)'  total_electron_number = ',total_electron_number
      endif  
    endif 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    select case(format_mode)
      case ("level_only") 
        format_mode_wrk="level_only" 
      case ("full") 
        format_mode_wrk="full" 
        part_ratio_is_on = .true.
      case default
        write(*,*)'ERROR(qm_output_levels):format_mode=',format_mode
        stop
    end select   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checking the allocation status
!
    if (.not. allocated(eig2)) then
      write(*,*)'INFO:(qm_output_levels):eig2 is not allocated --> skipped'
      return
    endif
!
    if (.not. allocated(f_occ)) then
      write(*,*)'INFO:(qm_output_levels):f_occ is not allocated --> skipped'
      return
    endif
!
    iunit=vacant_unit()
!
    if ( append_mode ) then
      if (i_verbose >= 1) then
        if (lu > 0) write(lu,*)' The output file is initialized : ',trim(filename_wrk)
      endif
      open(iunit, file=filename_wrk, status='unknown', position='append')
    else
      open(iunit, file=filename_wrk, status='unknown')
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculation of participation ratio (optional)KChecking the allocation status
!
    if (part_ratio_is_on) then
      allocate(part_ratio(size(eig2)), stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc error in part_ratio'
        stop
      endif
      if (trim(config%calc%solver%scheme) == 'eigen_mpi') then
        ! Participation ratio is already calculated in eig_solver_center() by set_pratio_mpi()
        part_ratio(:) = atmp(:, 1)
      else
        call gene_participation(part_ratio,lu)
      end if
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plot the eigen levels
!
    step_count=config%system%structure%mdstep
    do jj=1,size(eig2)
      select case(format_mode_wrk)
        case ("level_only") 
          write(iunit,'(i10,f20.12)') jj, eig2(jj)*unit_conv
        case ("full") 
          if (allocated(part_ratio)) then
            write(iunit,'(2i10,3f30.20)') step_count, jj, eig2(jj)*unit_conv, f_occ(jj), part_ratio(jj)
          else
            write(iunit,'(2i10,2f30.20)') step_count, jj, eig2(jj)*unit_conv, f_occ(jj)
          endif
        case default
          write(*,*)'ERROR(qm_output_levels,level):format_mode=',format_mode
          stop
      end select
    enddo
!
    close(iunit)
!
  end subroutine qm_output_levels
!
end module M_qm_output_levels

