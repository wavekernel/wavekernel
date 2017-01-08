!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_save_struct
!
   use M_io_dst_write_log, only : log_unit !(unchanged)
   implicit none
   logical, allocatable :: called_first_for_velocity_save(:)
!
   public :: elses_md_save_struct
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for saving structures  
!
  subroutine elses_md_save_struct(final_iteration)
!
    use M_config,           only : config
    use elses_mod_ctrl,     only : i_verbose
    use M_save_restart_xml, only : save_restart_xml !(routine)
    use M_lib_dst_info,     only : mpi_is_active, myrank, nprocs
    implicit none
    logical, intent(in) :: final_iteration
    logical             :: save_at_this_node
    logical             :: split_save_mode
!
    if (config%system%structure%use_vatom) then
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,*)'@@ elses_md_save_struct'
      endif 
    else
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,*)'@@ elses_md_save_struct is SKIPPED:use_vatom=',config%system%structure%use_vatom
      endif 
      return
    endif   
!
    if (.not. config%output%restart%set) then
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,*)'@@ elses_md_save_struct is SKIPPED:restart set=',config%output%restart%set
      endif
      return
    endif   
!
    save_at_this_node = .false.
    split_save_mode   = .false.
!
    if (mpi_is_active) then 
      if (myrank == 0) save_at_this_node = .true.
      if (config%output%restart%split) split_save_mode = .true.
    else
      save_at_this_node = .true.
    endif   
!
    if ( (save_at_this_node) .or. (split_save_mode) ) then
      call elses_xml_save_txp
!       ---> Copy structure data in config%system%structure%.....
!
      call save_restart_xml(final_iteration)
!       ---> Save restart XML file
!
    endif  
!
    if (save_at_this_node) then
!
      if ( trim(config%calc%mode) == "dynamics" ) then
        call save_velocity_data
      endif   
!
      call save_position_data
!       ---> Save position data for molecular visualization tools
    endif  
!
  end subroutine elses_md_save_struct
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine for saving velocity data for analysis
!
  subroutine save_velocity_data
!
    use M_config,           only : config     !(unchanged)
    use M_qm_domain,        only : i_verbose  !(unchanged)
    use elses_mod_md_dat,   only : itemd      !(unchanged)
    use elses_mod_file_io,    only : vacant_unit  !(function)
    use elses_mod_vel,      only : velx, vely, velz !(unchanged)
    implicit none
    integer :: ierr, fd, j
    logical :: flag_for_init
    character(len=100) filename_wrk
!
    if (i_verbose >= 1) then
      write(*,*)'@@ save_velocity_data'
    endif 
!
    if (config%output%velocity%interval == 0 ) then
      if (i_verbose >= 1) then
        write(*,*)'....is skipped, since (interval)=0'
      endif 
      return
    endif   
!
    if (.not. allocated(velx)) return
    if (.not. allocated(vely)) return
    if (.not. allocated(velz)) return
!
    if ( config%output%velocity%filename == "") return
    if (mod(itemd-1,config%output%velocity%interval) /= 0 ) return
!
    fd=vacant_unit()
    filename_wrk=trim(config%output%velocity%filename)
!
    if (.not. allocated(called_first_for_velocity_save)) then
      flag_for_init = .true. 
      allocate(called_first_for_velocity_save(1), stat=ierr)
      if (ierr /= 0) stop 'ERROR in alloc. (init_flag_for_velocity_save)'
      called_first_for_velocity_save(1) = .false.
    else
      flag_for_init = .false.
    endif
!
    if( flag_for_init ) then
       open(fd,file=filename_wrk)
    else
       open(fd,file=filename_wrk, position='append')
    end if
!
    write(fd,'(i10)') config%system%structure%natom 
    write(fd,*) trim(config%system%structure%name), " mdstep=", & 
&                                    config%system%structure%mdstep
!
    do j=1, config%system%structure%natom
       write(fd,'(a4,3e23.15)') &
&           config%system%structure%vatom(j)%name, &
&           config%system%structure%vatom(j)%velocity(1:3)
    end do
!
    close(fd)
!
  end subroutine save_velocity_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine for saving position data for visualization tools
!
  subroutine save_position_data

    use M_config,           only : config
    use M_qm_domain,        only : i_verbose, i_pbc_x, i_pbc_y, i_pbc_z
    use elses_mod_md_dat,   only : itemd
    implicit none
    integer :: lenf
    integer :: i_pbc
    logical :: flag_for_init
    character(len=100) ctemd
!   character(len=100) filename_org
    character(len=100) filename_xyz
    character(len=100) filename_axsf
    character(len=100) filename_xsf
    character(len=100) filename_pdb
    logical :: plot_for_xyz
    logical :: plot_for_axsf
    logical :: plot_for_xsf
    logical :: plot_for_pdb
    logical :: cell_info_mode
!
    if (i_verbose >= 1) then
      write(*,*)'@@ save_position_data'
    endif 
!
    if (config%output%position%interval == 0 ) then
      if (i_verbose >= 1) then
        write(*,*)'....is skipped, since (interval)=0'
      endif 
      return
    endif   
!
    if ( config%output%position%filename == "") return
    if (mod(itemd-1,config%output%position%interval) /= 0 ) return
!
!
    if (itemd == 1) then 
      flag_for_init = .true. 
    else
      flag_for_init = .false.
    endif   
!
    lenf = len_trim(config%output%position%filename)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Determine the file type to be generated
!
    plot_for_xyz  = .false.
    plot_for_axsf = .false.
    plot_for_xsf  = .false.
    plot_for_pdb  = .false.
    i_pbc = i_pbc_x*i_pbc_y*i_pbc_z

    filename_xyz  = ''
    filename_axsf = ''

    if((config%output%position%filename(lenf-3:lenf)=='.xyz') &
&         .or. (config%output%position%format =='xyz' ) ) then
       plot_for_xyz = .true.
       filename_xyz(1:lenf-4)=config%output%position%filename(1:lenf-4)
       filename_xyz(lenf-3:lenf)='.xyz'
    endif
!
    if((config%output%position%filename(lenf-4:lenf)=='.axsf') &
&         .or. (config%output%position%format =='axsf' ) ) then
       plot_for_axsf = .true.
       filename_axsf(1:lenf-5)=config%output%position%filename(1:lenf-5)
       filename_axsf(lenf-4:lenf+1)='.axsf'
    end if
!
    if((config%output%position%filename(lenf-3:lenf)=='.xsf') &
&         .or. (config%output%position%format =='xsf' ) ) then
        write(ctemd, '(i10.10)') itemd-1
        plot_for_xsf = .true.
    end if
!
    if((config%output%position%filename(lenf-3:lenf)=='.pdb') &
&         .or. (config%output%position%format =='pdb' ) ) then
        plot_for_pdb = .true.
    end if
!
    if(config%output%position%filename(lenf-3:lenf)=='.all') then
       plot_for_xyz  = .true.
       plot_for_xsf  = .true.
       plot_for_axsf = .true.
       plot_for_pdb  = .true.
       filename_xyz(1:lenf-4)=config%output%position%filename(1:lenf-4)
       filename_xyz(lenf-3:lenf)='.xyz'
       filename_axsf(1:lenf-4)=config%output%position%filename(1:lenf-4)
       filename_axsf(lenf-3:lenf+1)='.axsf'
      if (i_verbose >= 1) then
         write(*,*)' INFO:Experimental mode: Files in all formats are generated'
         write(*,*)' INFO:  ----> .xyz, .pdb, .xsf, .axsf'
      endif   
    end if
!
    if (lenf > 11) then
      if (config%output%position%filename(lenf-8:lenf)=='.xyz_axsf') then
         plot_for_xyz  = .true.
         plot_for_axsf = .true.
         filename_xyz(1:lenf-9)      = config%output%position%filename(1:lenf-9)
         filename_xyz(lenf-8:lenf-5)  = '.xyz'
         filename_axsf(1:lenf-9)     = config%output%position%filename(1:lenf-9)
         filename_axsf(lenf-8:lenf-4) = '.axsf'
         if (i_verbose >= 1) then
           write(*,*)' INFO:Experimental mode: two formats are generated; xyz and axsf'
         endif   
      endif   
    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Genaration of the output files
!
     if (plot_for_xyz) then
       write(*,*)'filename_xyz=',trim(filename_xyz)
       if (trim(config%output%position%cell_info) == "on") then
         cell_info_mode = .true. 
       else
         cell_info_mode = .false. 
       endif
       if (i_verbose >= 1) then
         write(*,*)' INFO:Whether cell_sizes are written in xyz=', cell_info_mode
       endif   
       call save_struct_xyz( config%system%structure, &
&                                  trim(filename_xyz), flag_for_init, cell_info_mode )
     endif  
!
     if (plot_for_xsf) then
      filename_xsf=''
      filename_xsf(1:lenf-4)=config%output%position%filename(1:lenf-4)
      write(ctemd, '(i10.10)') itemd-1
      if (i_pbc == 0) then
        call save_struct_xsf_npbc( config%system%structure, &
&              filename_xsf(1:lenf-4)//trim(ctemd)//'.xsf', flag_for_init )
      else  
        call save_struct_xsf( config%system%structure,      &
&              filename_xsf(1:lenf-4)//trim(ctemd)//'.xsf', flag_for_init )
      endif
     endif 
!
      if (plot_for_axsf) then
       write(*,*)'filename_wrk=',trim(filename_axsf)
       if (i_pbc == 0) then
         call save_struct_axsf_npbc( config%system%structure, &
&              trim(filename_axsf), flag_for_init )
       else
         call save_struct_axsf( config%system%structure,      &
&              trim(filename_axsf), flag_for_init )
       endif
      endif
!
      if (plot_for_pdb) then  
        if (flag_for_init) call elses_xml_vis
        filename_pdb=''
        filename_pdb(1:lenf-4)=config%output%position%filename(1:lenf-4)
        write(*,*)'filename_pdb=',trim(filename_pdb)
        call elses_vis_rasmol(itemd-1, trim(filename_pdb))
      endif
!
  end subroutine save_position_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save structure data in XYZ format
!
  subroutine save_struct_xyz( structure, filename, flag_for_init, cell_info_mode)
    use M_structure,          only : structure_type
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_phys_const, only : angst
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: flag_for_init
    logical, intent(in), optional:: cell_info_mode
    real(8)                      :: cell_size_angst(3)
    logical :: cell_info_mode_wrk
    integer ::fd
    integer :: j
!
    fd=vacant_unit()
!    
    if (present(cell_info_mode)) then 
       cell_info_mode_wrk = cell_info_mode
    else
       cell_info_mode_wrk = .false.
    endif   
!
    if( flag_for_init ) then
       open(fd,file=filename)
    else
       open(fd,file=filename,position='append')
    end if
!
    if (cell_info_mode_wrk) then
      cell_size_angst(1)=structure%unitcell%vectorA(1)*angst
      cell_size_angst(2)=structure%unitcell%vectorB(2)*angst
      cell_size_angst(3)=structure%unitcell%vectorC(3)*angst
      write(fd,'(i10,3f30.20)') structure%natom, cell_size_angst(1:3)
    else
      write(fd,'(i10)') structure%natom 
    endif
!
    write(fd,*) trim(structure%name), " mdstep=", structure%mdstep
!
    do j=1, structure%natom
       write(fd,'(a4,3e23.15)') &
            structure%vatom(j)%name, &
            structure%vatom(j)%position(1:3)*angst
    end do
!
    close(fd)
    return
!
  end subroutine save_struct_xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save structure data in AXSF format
! 20081222 Iitaka
!
  subroutine save_struct_axsf( structure, filename, flag_for_init )
    use M_config
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_phys_const, only : angst
    use elses_mod_md_dat,     only : itemd, itemdmx
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: flag_for_init
!  In future, these flags should be provided via xml
    logical, parameter :: flag_for_variable_cell=.true.
    logical, parameter :: flag_for_saving_force=.true.
    logical, parameter :: flag_for_system_name=.true.
!   logical, parameter :: flag_for_atomic_num=.true.
    logical, parameter :: flag_for_atomic_num=.false.
!
    integer ::fd
    integer :: j
    integer :: anum 
!       ----> atomic number 
    integer :: nstep ! total animation step
    integer :: istep ! current animation step
!
    if (config%output%position%interval == 0 ) then
       write(*,*)'ERROR:(save_struct_axsf):interval=',config%output%position%interval
       stop
    endif   
!
    nstep=(itemdmx-1)/config%output%position%interval+1
    istep=(itemd-1)/config%output%position%interval+1
! variable cell
    if(flag_for_variable_cell) then
!
    fd=vacant_unit()
!    
    if( flag_for_init ) then
       open(fd,file=filename)
    else
       open(fd,file=filename,position='append')
    end if

    if( flag_for_init ) then
      write(fd,*) '# TOTAL MD STEP ', itemdmx
      write(fd,*) 'ANIMSTEPS', nstep
      write(fd,*) 'CRYSTAL'
    endif
    if (flag_for_system_name) write(fd,*) '# SYSTEM :',trim(structure%name)
    write(fd,*) '# MD STEP ', structure%mdstep
    write(fd,*) 'PRIMVEC',istep
    write(fd,'(3e23.15)') structure%unitcell%vectorA(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorB(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorC(1:3)*angst
    write(fd,*) 'CONVVEC',istep
    write(fd,'(3e23.15)') structure%unitcell%vectorA(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorB(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorC(1:3)*angst
    write(fd,*) 'PRIMCOORD',istep
    write(fd,*) structure%natom, 1
    if(flag_for_saving_force)then
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,6e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst &
&               , structure%vatom(j)%force(1:3)/angst
        else
          write(fd,'(a4,6e23.15)') &
&           structure%vatom(j)%name, &
&           structure%vatom(j)%position(1:3)*angst, &
&           structure%vatom(j)%force(1:3)/angst
        endif   
      end do
    else
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,3e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst
        else
          write(fd,'(a4,3e23.15)') &
&           structure%vatom(j)%name, &
&           structure%vatom(j)%position(1:3)*angst
        endif  
      end do
    endif
!
    close(fd)
! end of variable cell
    else
! fixed cell
!
    fd=vacant_unit()
!    
    if( flag_for_init ) then
       open(fd,file=filename)
    else
       open(fd,file=filename,position='append')
    end if

    if( flag_for_init ) then
      write(fd,*) '# TOTAL MD STEP ', itemdmx
      write(fd,*) 'ANIMSTEPS', nstep
      write(fd,*) 'CRYSTAL'
      if (flag_for_system_name) write(fd,*) '# SYSTEM :',trim(structure%name)
      write(fd,*) 'PRIMVEC'
      write(fd,'(3e23.15)') structure%unitcell%vectorA(1:3)*angst
      write(fd,'(3e23.15)') structure%unitcell%vectorB(1:3)*angst
      write(fd,'(3e23.15)') structure%unitcell%vectorC(1:3)*angst
      write(fd,*) 'CONVVEC'
      write(fd,'(3e23.15)') structure%unitcell%vectorA(1:3)*angst
      write(fd,'(3e23.15)') structure%unitcell%vectorB(1:3)*angst
      write(fd,'(3e23.15)') structure%unitcell%vectorC(1:3)*angst
    endif
    write(fd,*) '# MD STEP ', structure%mdstep
    write(fd,*) 'PRIMCOORD',istep
    write(fd,*) structure%natom, 1
    if(flag_for_saving_force)then
      do j=1, structure%natom
        call get_atomic_num(trim(structure%vatom(j)%name),anum)
        write(fd,'(I4,6e23.15)') &
&         anum, structure%vatom(j)%position(1:3)*angst &
&             , structure%vatom(j)%force(1:3)/angst 
      end do
    else
      do j=1, structure%natom
        call get_atomic_num(trim(structure%vatom(j)%name),anum)
        write(fd,'(I4,3e23.15)') &
&         anum, structure%vatom(j)%position(1:3)*angst
!        write(fd,'(a4,3e23.15)') &
!&         structure%vatom(j)%name, &
!&         structure%vatom(j)%position(1:3)*angst
      end do
    endif
!
    close(fd)
! end of fixed cell
    endif
    return
!
  end subroutine save_struct_axsf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save structure data in AXSF format in the non-periodic boundary conditions
!
  subroutine save_struct_axsf_npbc( structure, filename, flag_for_init )
    use M_config
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_phys_const, only : angst
    use elses_mod_md_dat,     only : itemd, itemdmx
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: flag_for_init
!
    logical, parameter :: flag_for_saving_force=.true.
    logical, parameter :: flag_for_system_name=.true.
!   logical, parameter :: flag_for_atomic_num=.true.
    logical, parameter :: flag_for_atomic_num=.false.
!
    integer ::fd
    integer :: j
    integer :: anum 
!       ----> atomic number 
    integer :: nstep ! total animation step
    integer :: istep ! current animation step
!
    if (config%output%position%interval == 0 ) then
       write(*,*)'ERROR:(save_struct_axsf):interval=',config%output%position%interval
       stop
    endif   
!
    nstep=(itemdmx-1)/config%output%position%interval+1
    istep=(itemd-1)/config%output%position%interval+1
!
    fd=vacant_unit()
!    
    if( flag_for_init ) then
       open(fd,file=filename)
    else
       open(fd,file=filename,position='append')
    end if

    if( flag_for_init ) then
      write(fd,*) '# TOTAL MD STEP ', itemdmx
      write(fd,*) 'ANIMSTEPS', nstep
    endif
    if (flag_for_system_name) write(fd,*) '# SYSTEM :',trim(structure%name)
    write(fd,*) '# MD STEP ', structure%mdstep
    write(fd,*) 'ATOMS',istep
    if(flag_for_saving_force)then
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,6e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst &
&               , structure%vatom(j)%force(1:3)/angst
        else
          write(fd,'(a4,6e23.15)') &
&           structure%vatom(j)%name, &
&           structure%vatom(j)%position(1:3)*angst, &
&           structure%vatom(j)%force(1:3)/angst
        endif   
      end do
    else
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,3e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst
        else
          write(fd,'(a4,3e23.15)') &
&           structure%vatom(j)%name, &
&           structure%vatom(j)%position(1:3)*angst
        endif  
      end do
    endif
!
    close(fd)
!
  end subroutine save_struct_axsf_npbc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save structure data in XSF format
!
  subroutine save_struct_xsf( structure, filename, flag_for_init )
    use M_structure,          only : structure_type
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_phys_const, only : angst, au_mass
    use elses_mod_md_dat,     only : itemdmx
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: flag_for_init
!  In future, this flag should be provided via xml
    logical, parameter :: flag_for_saving_force=.true.
    logical, parameter :: flag_for_system_name=.true.
!   logical, parameter :: flag_for_atomic_num=.true.
    logical, parameter :: flag_for_atomic_num=.false.
    integer ::fd
    integer :: j
    integer :: anum 
!       ----> atomic number 
!
    fd=vacant_unit()
!    
    if( flag_for_init ) then
       open(fd,file=filename)
    else
       open(fd,file=filename,position='append')
    end if

    write(fd,*) 'CRYSTAL'
    if (flag_for_system_name) write(fd,*) '# SYSTEM :',trim(structure%name)
    write(fd,*) 'PRIMVEC'
    write(fd,'(3e23.15)') structure%unitcell%vectorA(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorB(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorC(1:3)*angst
    write(fd,*) 'CONVVEC'
    write(fd,'(3e23.15)') structure%unitcell%vectorA(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorB(1:3)*angst
    write(fd,'(3e23.15)') structure%unitcell%vectorC(1:3)*angst
    write(fd,*) 'PRIMCOORD'
    write(fd,*) structure%natom, 1
    if(flag_for_saving_force)then
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,6e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst &
&               , structure%vatom(j)%force(1:3)/angst
        else
          write(fd,'(a4,6e23.15)') &
&           structure%vatom(j)%name, &
&           structure%vatom(j)%position(1:3)*angst, &
&           structure%vatom(j)%force(1:3)/angst
        endif  
      end do
    else
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,3e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst
        else  
          write(fd,'(a4,3e23.15)') &
&         structure%vatom(j)%name, &
&         structure%vatom(j)%position(1:3)*angst
        endif  
      end do
    endif
!
    close(fd)
    return
!
  end subroutine save_struct_xsf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save structure data in XSF format in non-periodic boundary conditions
!
  subroutine save_struct_xsf_npbc( structure, filename, flag_for_init )
    use M_structure,          only : structure_type
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_phys_const, only : angst, au_mass
    use elses_mod_md_dat,     only : itemdmx
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: flag_for_init
!  In future, this flag should be provided via xml
    logical, parameter :: flag_for_saving_force=.true.
    logical, parameter :: flag_for_system_name=.true.
!   logical, parameter :: flag_for_atomic_num=.true.
    logical, parameter :: flag_for_atomic_num=.false.
    integer ::fd
    integer :: j
    integer :: anum 
!       ----> atomic number 
!
    fd=vacant_unit()
!    
    if( flag_for_init ) then
       open(fd,file=filename)
    else
       open(fd,file=filename,position='append')
    end if

    write(fd,*) 'ATOMS'
    if(flag_for_saving_force)then
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,6e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst &
&               , structure%vatom(j)%force(1:3)/angst
        else
          write(fd,'(a4,6e23.15)') &
&           structure%vatom(j)%name, &
&           structure%vatom(j)%position(1:3)*angst, &
&           structure%vatom(j)%force(1:3)/angst
        endif  
      end do
    else
      do j=1, structure%natom
        if (flag_for_atomic_num) then
          call get_atomic_num(trim(structure%vatom(j)%name),anum)
          write(fd,'(I4,3e23.15)') &
&           anum, structure%vatom(j)%position(1:3)*angst
        else  
          write(fd,'(a4,3e23.15)') &
&         structure%vatom(j)%name, &
&         structure%vatom(j)%position(1:3)*angst
        endif  
      end do
    endif
!
    close(fd)
    return
!
  end subroutine save_struct_xsf_npbc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get atomic number from element name
!
  subroutine get_atomic_num(element_name, atom_num)
    use M_lib_element_database, only : get_atomic_num_new => get_atomic_num
    character(len=*), intent(in) :: element_name
    integer,          intent(out) :: atom_num
!
    call get_atomic_num_new(element_name, atom_num)
    return
!
  end subroutine get_atomic_num
!
end module M_md_save_struct
