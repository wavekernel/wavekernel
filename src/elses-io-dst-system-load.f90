!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_io_dst_system_load
!
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
!
   implicit none
   character(len=*), parameter :: file_name_structure='input_structure.txt'
!
   private
   public file_name_structure
   public dst_structure_load_from_txt
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine dst_structure_load_from_txt
!
    use M_config  !(CHANGED)
    use elses_mod_file_io,  only : vacant_unit  !(function)
    use M_lib_phys_const,   only : angst        !(parameter)
    use M_io_read_xyz_file, only : get_basic_info_from_xyz !(routine)
    use M_io_read_xyz_file, only : get_position_from_xyz   !(routine)
    implicit none
!   character(len=32) :: chara1
    logical :: file_exist
    integer ierr
    integer atm_index, elm_index
!
    character(len=4) :: name_wrk
!   real(DOUBLE_PRECISION) :: ddx, ddy, ddz
!
    integer snap_index, num_atoms, num_species
    character(len=4)   :: species(100)
    character(len=100) :: header
    real(DOUBLE_PRECISION), allocatable   :: atm_data_position(:,:)
    character(len=4),       allocatable   :: atm_data_elem(:)
    real(DOUBLE_PRECISION)                :: cell_data(3)
!
    integer j
!
    snap_index=config%calc%snapshot%initial
!
    write(*,*)'@@ dst_structure_load_from_txt:snap_index=',snap_index
    if (log_unit > 0) write(log_unit,*)'@@ dst_structure_load_from_txt:snap_index=',snap_index
!
!
    inquire (file=file_name_structure, exist=file_exist)
    if (.not. file_exist) then 
      write(*,*)' INFO:File does not exist:',trim(file_name_structure)
      if (log_unit > 0) write(log_unit,*)' INFO:File does not exist:',trim(file_name_structure)
      return
    else
      write(*,*)' INFO:EXYZ file exists:',trim(file_name_structure)
    endif
!
!   snap_index=0
    call get_basic_info_from_xyz(snap_index, file_name_structure, num_atoms, num_species, header, species, cell_data)
!                          --------> cell_data(1:3) : unit in a.u.
!
    write(*,*)' exyz file info:  header       = ', trim(header)
    write(*,*)' exyz file info:  num_atoms    = ', num_atoms
    write(*,*)' exyz file info:  num_species  = ', num_species
    write(*,*)' exya file info:  cell_data(au) =', cell_data(1:3)
    do j=1, num_species
      write(*,*)' exyz file info:  element list =',j, species(j)
    enddo
!
    config%system%structure%natom    = num_atoms
    config%system%structure%nelement = num_species
!
    if (cell_data(1) < 0.0d0) then
      write(*,*)'INFO: No meaningful cell data is found in the exyz file'
      if (config%system%structure%unitcell%set .eqv. .false.) then
        write(*,*)'ERROR:config%system%structure%unitcell%set =',config%system%structure%unitcell%set 
        stop
      endif
    else
      write(*,*)'INFO: A meaningful cell data is found in xyz file'
      config%system%structure%unitcell%set        = .true.
      config%system%structure%unitcell%vectorA(:) = 0.0d0
      config%system%structure%unitcell%vectorB(:) = 0.0d0
      config%system%structure%unitcell%vectorC(:) = 0.0d0
      config%system%structure%unitcell%vectorA(1) = cell_data(1)
      config%system%structure%unitcell%vectorB(2) = cell_data(2)
      config%system%structure%unitcell%vectorC(3) = cell_data(3)
    endif
!
    if (config%system%structure%unitcell%set .eqv. .false.) then
      write(*,*)'ERROR(structure_load):config%system%structure%unitcell%set =',config%system%structure%unitcell%set 
      stop
    endif
!
    write(*,*)'before setting : natom   =',config%system%structure%natom
    write(*,*)'before setting : nelement=',config%system%structure%nelement
    if (log_unit > 0) write(log_unit,*)'before setting : natom   =',config%system%structure%natom
    if (log_unit > 0) write(log_unit,*)'before setting : nelement=',config%system%structure%nelement
!
    write(*,*)'INFO:allocate structure%vatom(structure%natom) at dst_structure_load_from_txt'
    if (log_unit > 0) then
      write(log_unit,*)'INFO:allocate structure%vatom(structure%natom) at dst_structure_load_from_txt'
    endif
!
    allocate( config%system%structure%vatom(config%system%structure%natom), stat=ierr )
    if (ierr /= 0) stop 'Alloc ERROR(vatom)'
!
    write(*,*)'INFO:allocate structure%nelement at dst_structure_load_from_txt'
    if (log_unit > 0) then
      write(log_unit,*) 'INFO:allocate structure%nelement at dst_structure_load_from_txt'
    endif
!
    allocate( config%system%structure%velement(config%system%structure%nelement), stat=ierr )
    if (ierr /= 0) stop 'Alloc ERROR(velement)'
!
    do elm_index=1,config%system%structure%nelement
      name_wrk=trim(species(elm_index))
      config%system%structure%velement(elm_index)%name=trim(name_wrk)
      config%system%structure%velement(elm_index)%filename=trim(name_wrk)//'.xml'
      config%system%structure%velement(elm_index)%quantum%type = 'geno'
    enddo
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
    do atm_index=1,config%system%structure%natom
      config%system%structure%vatom(atm_index)%name       =trim(atm_data_elem(atm_index))
      config%system%structure%vatom(atm_index)%position(1:3)=atm_data_position(1:3,atm_index) ! value in au
      config%system%structure%vatom(atm_index)%class      =""
      config%system%structure%vatom(atm_index)%motion     ="free"
      config%system%structure%vatom(atm_index)%velocity_set             = .false.
      config%system%structure%vatom(atm_index)%force_set                = .false.
      config%system%structure%vatom(atm_index)%population_set           = .false.
      config%system%structure%vatom(atm_index)%population_guess_set     = .false.
    enddo
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
  end subroutine dst_structure_load_from_txt
!
end module M_io_dst_system_load





