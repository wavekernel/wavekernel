!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_io_update_snapshot
!
   use M_qm_domain,      only : i_verbose, DOUBLE_PRECISION !(unchanged)
   use M_io_dst_write_log,  only : log_unit
!
   implicit none
!
   private
   public :: update_snapshot
!
 contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine update_snapshot
!
    use M_qm_domain,          only : ax, ay, az           !(CHNAGED)
    use elses_mod_txp,        only : txp, typ, tzp        !(CHANGED)
!
    use M_config  !(unchanged)(only : )
    use M_io_dst_system_load, only : file_name_structure  !(unchanged)
    use elses_mod_md_dat,     only : itemd                !(unchanged)
!
!   use elses_mod_file_io,  only : vacant_unit  !(function)
!   use M_lib_phys_const,   only : angst        !(parameter)
    use M_io_read_xyz_file, only : get_basic_info_from_xyz !(routine)
    use M_io_read_xyz_file, only : get_position_from_xyz   !(routine)
!   use M_md_motion,        only : md_update_tx            !(routine)
!
    implicit none
!   character(len=*), parameter :: file_name_structure='input_structure.txt'
!   character(len=32) :: chara1
    logical :: file_exist
    integer ierr
    integer atm_index, elm_index
!
    character(len=4) :: name_wrk
!   real(DOUBLE_PRECISION) :: ddx, ddy, ddz
!
    integer num_atoms, num_species
    character(len=4)   :: species(100)
    character(len=100) :: header
    real(DOUBLE_PRECISION), allocatable   :: atm_data_position(:,:)
    character(len=4),       allocatable   :: atm_data_elem(:)
    real(DOUBLE_PRECISION)                :: cell_data(3)
!
    integer :: j
!
    integer :: snap_index, snap_initial, snap_interval
!
    snap_initial  = config%calc%snapshot%initial
    snap_interval = config%calc%snapshot%interval
!
    snap_index=snap_initial+itemd*snap_interval
!
    if (i_verbose >= 0) write(*,*)'@@ update-snapshot:snap_index=',snap_index
    if (log_unit > 0) write(log_unit,*)'@@ update-snapshot:snap_index=',snap_index
!
    inquire (file=file_name_structure, exist=file_exist)
    if (.not. file_exist) then 
      write(*,*)'ERROR(update-snapshot):File does not exist:',trim(file_name_structure)
      stop
    endif
!
    call get_basic_info_from_xyz(snap_index, file_name_structure, num_atoms, num_species, header, species, cell_data)
!                          --------> cell_data(1:3) : unit in a.u.
!
    if (cell_data(1) > 0.0d0) then
      ax=cell_data(1)
      ay=cell_data(2)
      az=cell_data(3)
      if (i_verbose >= 0) write(*,*)'INFO:update-snapshot:ax(au)=',ax
      if (i_verbose >= 0) write(*,*)'INFO:update-snapshot:ay(au)=',ay
      if (i_verbose >= 0) write(*,*)'INFO:update-snapshot:az(au)=',az
      if (log_unit > 0) write(log_unit,*) 'INFO:update-snapshot:ax(au)=',ax
      if (log_unit > 0) write(log_unit,*) 'INFO:update-snapshot:ay(au)=',ay
      if (log_unit > 0) write(log_unit,*) 'INFO:update-snapshot:az(au)=',az
    endif  
!
    allocate(atm_data_elem(num_atoms), stat=ierr )
    if (ierr /= 0) stop 'Alloc ERROR(atm_data_elem)'
!
    allocate(atm_data_position(3, num_atoms), stat=ierr )
    if (ierr /= 0) stop 'Alloc ERROR(atm_data_position)'
!
    call get_position_from_xyz(snap_index, file_name_structure, species, atm_data_position, atm_data_elem)
!                         --------> atm_data_position(1:3,num_atoms) : unit in a.u.
!
    do atm_index=1, num_atoms
      txp(atm_index)=atm_data_position(1,atm_index)/ax
      typ(atm_index)=atm_data_position(2,atm_index)/ay
      tzp(atm_index)=atm_data_position(3,atm_index)/az
    enddo
!
!   call md_update_tx
!
    deallocate(atm_data_elem, stat=ierr )
    if (ierr /= 0) stop 'Dealloc ERROR(atm_data_elem)'
!
    deallocate(atm_data_position, stat=ierr )
    if (ierr /= 0) stop 'Dealloc ERROR(atm_data_position)'
!
    write(*,*)'.... end. : dst_structure_load_from_txt'
    if (log_unit > 0) write(log_unit,*)'.... end. : dst_structure_load_from_txt'
!
!   stop
!
  end subroutine update_snapshot
!
end module M_io_update_snapshot

