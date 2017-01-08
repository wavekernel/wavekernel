!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_io_dst_write_log
!
 use M_lib_dst_info,   only : log_unit    !(CHANGED in set_dst_write_log_file)
 use M_qm_domain,      only : i_verbose   !(uncanged)
 implicit none  
 logical ::  log_file_is_set  ! set in set_dst_write_log_file
!
  private
  public :: log_file_is_set, log_unit
  public :: set_dst_write_log_file
!
  contains
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set the variable log_unit
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine set_dst_write_log_file
!
    use M_config,          only : config      !(unchanged) (only config%option%log_node_number, config%calc%distributed%log_unit)
    use M_lib_dst_info,    only : output_unit !(CHANGED)
    use elses_mod_file_io, only : vacant_unit !(function)
    use M_lib_dst_info,    only : mpi_is_active, myrank, root_node !(unchanged)
!   use M_lib_dst_info,    only : output_filename  !(unchanged)
    use M_lib_dst_info,    only : true_mpi_wrapper !(unchanged)
    use M_00_v_info,       only : version_info     !(unchanged)
    implicit none
    character(len=32)    :: log_node_filename
    character(len=32)    :: myrank_chara
    integer              :: max_node_for_log
    character(len=8)  :: chara_date
    character(len=10) :: chara_time
!
    log_file_is_set = .true.
    max_node_for_log  = config%option%log_node_number - 1
!
    write(myrank_chara, '(i6.6)') myrank
    log_node_filename=trim(config%option%output_dir)//'log-node'//trim(myrank_chara)//'.txt'
!
    if (myrank > max_node_for_log) then
      log_file_is_set = .false.
      log_node_filename=trim(config%option%output_dir)//'log-nodeXXXXXX.txt'  ! dummy 
    endif
!   
    output_unit = -1
!   if (root_node) then
!     output_unit=vacant_unit()
!     open (output_unit, file=output_filename, status='unknown')
!   endif
!  
    log_unit=-1
    if (log_file_is_set) then 
      log_unit=vacant_unit()
      open (log_unit, file=log_node_filename, status='unknown')
    endif
!  
    config%calc%distributed%log_unit=log_unit 
!
!   if (i_verbose >= 0) then 
!     write(*,'(a)') '@@ ELSES Project'
!     write(*,'(a)') trim(version_info)
!   endif
!  
    if (log_unit > 0) then 
      write(log_unit,'(a,a)') '@@ ', trim(version_info)
      if (true_mpi_wrapper) then
        write(log_unit,'(a)') 'INFO-MPI: True  MPI wrapper' 
      else
        write(log_unit,'(a)') 'INFO-MPI: Dummy MPI wrapper' 
      endif   
    endif  
!
    if (log_unit > 0) then 
      call date_and_time(chara_date, chara_time)
      write(log_unit,'(13a)')'Date: ', chara_date(1:4), ' ', chara_date(5:6), ' ', chara_date(7:8), '; ', &
&                             'Time: ', chara_time(1:2), ' ', chara_time(3:4), ' ', chara_time(5:6)
    endif   
!
!   if (output_unit > 0) then 
!     write(output_unit,'(a)') '@@ Output of ELSES'
!     write(output_unit,'(a)') trim(version_info)
!   endif  
!
  end subroutine set_dst_write_log_file
!    
end module M_io_dst_write_log

