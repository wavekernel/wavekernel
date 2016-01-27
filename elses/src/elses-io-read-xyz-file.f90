!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_io_read_xyz_file
!
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
!
   private
   public get_basic_info_from_xyz
   public get_position_from_xyz
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine get_basic_info_from_xyz(snap_index, file_name, num_atoms, num_species, header, species, cell_data_wrk)
!
    use elses_mod_file_io,    only : vacant_unit !(function)
    use M_lib_phys_const,     only : angst       !(unchanged)
    implicit none
    integer,          intent(in)    :: snap_index
    integer,          intent(out)   :: num_atoms
    integer,          intent(out)   :: num_species
    character(len=*), intent(in)    :: file_name
    character(len=4), intent(inout) :: species(:)
    character(len=*), intent(inout) :: header
    real(DOUBLE_PRECISION), intent(out) :: cell_data_wrk(3)
    character(len=4)                :: chara2      ! Element name should be less than 4 characters 
    integer                         :: unit_num
    integer                         :: atm_index, elm_index, spe_index
    integer                         :: snp_index_wrk
    integer                         :: ierr, ierr2
    integer                         :: cell_info_value
!
    if (i_verbose >=1) then
      if (log_unit > 0) then
        write(log_unit,*)'@@ get_num_atoms_from_xyz:snap_index, file_name=',snap_index, file_name
      endif  
    endif  
!
    unit_num = vacant_unit()
    cell_data_wrk(1:3)=-1.0d0   ! Dummy data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check whether the cell info is included or not.
!    
    open(unit_num, file=file_name)

    cell_info_value=0
    call read_first_line_of_xyz(cell_info_value, unit_num, num_atoms, cell_data_wrk)
    rewind(unit=unit_num)    
!
    if (i_verbose >=1) then
      write(*,*)'cell_info_value=',cell_info_value
      if (log_unit > 0) then
        write(log_unit,*)'cell_info_value=',cell_info_value
      endif
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
    do snp_index_wrk=0, snap_index
!
      if (snp_index_wrk /= snap_index) then
        call read_first_line_of_xyz(cell_info_value, unit_num, num_atoms, cell_data_wrk)
      else
        call read_first_line_of_xyz(cell_info_value, unit_num, num_atoms, cell_data_wrk)
        if (i_verbose >=1) then
          write(*,*)'cell size [A]=',cell_data_wrk(1:3)
          if (log_unit > 0) then
            write(log_unit,*)'cell size [A]=',cell_data_wrk(1:3)
          endif
        endif  
        cell_data_wrk(1:3)=cell_data_wrk(1:3)/angst  ! cell data in a.u.
      endif
!
      read(unit_num,*) header
!
      num_species=0
      do atm_index=1,num_atoms
        read(unit_num, *) chara2
        if (snp_index_wrk /= snap_index) cycle
        elm_index=0
        do spe_index=1,num_species
          if (trim(chara2) == trim(species(spe_index))) then 
            elm_index=spe_index
            exit
          endif
        enddo
        if (elm_index == 0) then
          num_species=num_species+1
          species(num_species)=trim(chara2)
!         write(*,*)' element =', num_species, trim(chara2)
        endif
      enddo
!
    enddo
!
    close(unit_num)
!
  end subroutine get_basic_info_from_xyz
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine get_position_from_xyz(snap_index, file_name, species, atm_data_position, atm_data_elem)
!
!
    use elses_mod_file_io,    only : vacant_unit !(function)
    use M_lib_phys_const,     only : angst       !(unchanged)
    implicit none
    integer,                intent(in)    :: snap_index
    character(len=*),       intent(in)    :: file_name
    character(len=4),       intent(in)    :: species(:)
    real(DOUBLE_PRECISION), intent(inout) :: atm_data_position(:,:)
    character(len=4),       intent(inout) :: atm_data_elem(:)
    character(len=4)                      :: chara2      ! Element name should be less than 4 characters 
    character(len=100)                    :: header      ! Header of the file
    integer :: num_atoms, num_species
    integer :: unit_num
    integer :: atm_index
    integer :: snp_index_wrk
    integer :: ierr, ierr2
    real(DOUBLE_PRECISION) :: dddx, dddy, dddz
!   real(DOUBLE_PRECISION) :: cell_data_wrk(3)
!
    atm_data_position(:,:)=-1.0d0 ! dummy data
    atm_data_elem(:)=''           ! dummy data
!
    num_atoms=size(atm_data_position,2)
    write(*,*)' num_atoms   = ',num_atoms
    num_species=size(species,1)
    write(*,*)' num_species = ',num_species
!
!   if (snap_index /= 0) then
!     write(*,*)'ERROR(get_basic_info_from_xyz):snap_index=',snap_index
!     stop
!   endif
!
    unit_num = vacant_unit()
!
    open(unit_num, file=file_name)
!
    do snp_index_wrk=0, snap_index
!
      read(unit_num,*) num_atoms
      read(unit_num,*) header
!
!     if (snp_index_wrk == snap_index) then
!       write(*,*)'INFO:cell_data [A] =', cell_data_wrk(1:3)
!       cell_data(1:3)=cell_data_wrk(1:3)/angst
!       write(*,*)'  header = ', trim(header)
!     endif
!
      do atm_index=1,num_atoms
        read(unit_num, *) chara2, dddx, dddy, dddz
        if (snp_index_wrk /= snap_index) cycle
        atm_data_elem(atm_index)=trim(chara2)
        atm_data_position(1,atm_index)=dddx/angst
        atm_data_position(2,atm_index)=dddy/angst
        atm_data_position(3,atm_index)=dddz/angst
      enddo
!
    enddo
!
    close(unit_num)
!
  end subroutine get_position_from_xyz
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine read_first_line_of_xyz(cell_info_value, unit_num, num_atoms, cell_data_wrk)
!
!       cell_info_value =  0 : check whether the cell info is included or not
!       cell_info_value =  1 : read the number of atoms and the cell info
!       cell_info_value = -1 : read the number of atoms
!
    implicit none
    integer,          intent(inout)    :: cell_info_value
    integer,          intent(in)       :: unit_num
    integer,          intent(out)      :: num_atoms
    real(8),          intent(out)      :: cell_data_wrk(3)
    integer                            :: ierr
!
    num_atoms=0                ! dummy data
    cell_data_wrk(1:3)=-10.d0  ! dummy data
!
    if (cell_info_value == 0) then 
      read(unit_num, *, iostat=ierr) num_atoms, cell_data_wrk(1:3)
      if (ierr == 0) then
        cell_info_value = 1
      else  
        cell_info_value = -1
      endif  
      return
    endif  
!
    if (cell_info_value == 1) then 
      read(unit_num, *, iostat=ierr) num_atoms, cell_data_wrk(1:3)
      if (ierr /= 0) then
        write(*,*)'ERROR(read_first_line_of_xyz):cell_info_value, ierr=' ,cell_info_value, ierr
        stop
      endif  
    endif  
!
    if (cell_info_value == -1) then 
      read(unit_num, *, iostat=ierr) num_atoms
      if (ierr /= 0) then
        write(*,*)'ERROR(read_first_line_of_xyz):cell_info_value, ierr=' ,cell_info_value, ierr
        stop
      endif  
    endif  
!
  end subroutine read_first_line_of_xyz
!
end module M_io_read_xyz_file







