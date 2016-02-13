!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_cell_change
!
   use elses_mod_ctrl,       only : i_verbose
   use elses_mod_file_io,    only : vacant_unit
   use M_io_dst_write_log,   only : log_unit !(unchanged)
!
   integer, parameter   :: DOUBLE_PRECISION=8
   real(DOUBLE_PRECISION), allocatable :: cell_parameter_store(:)

   public :: elses_md_cell_change
!
 contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cell change 
!
  subroutine elses_md_cell_change(for_next_step)
!
    use M_config,             only : config      
!       !(CHANGED, only config%system%structure%unitcell, config%system%cutoff_radius)
!       !(CHANGED, only config%calc%cell_change%max_num_iter
    use elses_mod_sim_cell,   only : ax, ay, az  !(CHANGED)
    use elses_mod_sim_cell,   only : i_pbc_x, i_pbc_y, i_pbc_z  !(CHANGED)
    use elses_mod_sel_sys,    only : r_cut_book          !(CHANGED)
    use M_qm_domain,          only : i_verbose, c_system !(unchanged)
    use elses_mod_md_dat,     only : itemdmx             !(CHANGED)
    use elses_mod_vel,        only : velx, vely, velz    !(CHANGED)
    use M_md_velocity_routines, only : calc_kinetic_energy !(routine)
    use elses_mod_phys_const, only : angst
!
    implicit none
    logical, intent(in) :: for_next_step
    logical :: file_exist
    integer :: line_count_max
    character(len=64) :: file_name_wrk
    character(len=32)  :: chara1
    integer :: iunit,ierr
    integer :: j, n
    integer :: imode
    integer :: step_count
    real(DOUBLE_PRECISION) :: cell_size_ref(3)
    real(DOUBLE_PRECISION) :: cell_scale(3)
    real(DOUBLE_PRECISION) :: dddx, dddy, dddz
    real(DOUBLE_PRECISION) :: r_cut_book_old
    logical                :: rescale_cutoff
!
    real(DOUBLE_PRECISION) :: ax_old, ay_old, az_old
    logical                :: flag_for_scale
    logical                :: velocity_scaling
!
    logical                :: check_kinetic_energy
    real(DOUBLE_PRECISION) :: kinetic_energy, kinetic_energy_old
!
    character(len=1024)    :: chara_unit_wrk, chara_unit
    character(len=10240)   :: chara_read_line
    real(DOUBLE_PRECISION) :: unit_value
!
    logical                :: called_first
!

    kinetic_energy     = 0.0d0
    kinetic_energy_old = 0.0d0
!
    called_first = .false.
    if (config%calc%cell_change%cell_sizes_in_xml_file(1) < 0.0d0) called_first = .true.
    if (config%calc%cell_change%cell_sizes_in_xml_file(2) < 0.0d0) called_first = .true.
    if (config%calc%cell_change%cell_sizes_in_xml_file(3) < 0.0d0) called_first = .true.
    if (called_first) then
      config%calc%cell_change%cell_sizes_in_xml_file(1)=config%system%structure%unitcell%vectorA(1)
      config%calc%cell_change%cell_sizes_in_xml_file(2)=config%system%structure%unitcell%vectorB(2)
      config%calc%cell_change%cell_sizes_in_xml_file(3)=config%system%structure%unitcell%vectorC(3)
    endif
!
    rescale_cutoff=.false.
    if (c_system == 'geno')    rescale_cutoff= .true.
    if (c_system == 'NRL_leg') rescale_cutoff= .true.
!
    velocity_scaling     = .false.
    check_kinetic_energy = .false.
    if ( (config%calc%mode == "dynamics") .and. (for_next_step) ) then
      velocity_scaling     = .true.
      check_kinetic_energy = .true.
    endif   
!
    if (config%calc%mode == "cell_change_only") then
      if (.not. config%calc%cell_change%set) then
        write(*,*)'ERROR(elses_md_cell_change):cell_change_only mode'
        write(*,*)'ERROR(elses_md_cell_change):config%calc%cell_change%set=',config%calc%cell_change%set 
        write(*,*)'# The cell_change tag may be missing in the cell-change-only mode'
        stop
      endif   
    endif   
!
    if (.not. config%calc%cell_change%set) return
!
    if (config%calc%cell_change%scheme /= "from_file")  then
      return 
    endif   
!
    file_name_wrk=trim(config%calc%cell_change%filename)
!
!   file_name='input_cell_change.txt'
!
    if (for_next_step) then
      step_count=config%system%structure%mdstep+1
    else
      step_count=config%system%structure%mdstep
    endif  
!
    if (i_verbose >= 1) then
       if (log_unit > 0) write(log_unit,*)'@@ elses_md_cell_change:step_count=', step_count
    endif   
!
    inquire (file=trim(file_name_wrk), exist=file_exist)
    if (i_verbose >= 1) then 
      if (log_unit > 0) write(log_unit,*)'file_exist=',file_exist
    endif   
!
    if (file_exist .eqv. .false.) then
      if (config%calc%mode == "cell_change_only ") then
        write(*,*) 'ERROR:cell_change_only : file is missing:',trim(file_name_wrk)
        stop 
      endif   
    endif   
!
    if (file_exist .eqv. .false.) then
      if (i_verbose >= 1) then 
        if (log_unit > 0) write(log_unit,*)'The file of input_cell_change.txt does not exist'
      endif   
      return
    else
      if (i_verbose >= 1) then 
        if (log_unit > 0) write(log_unit,*)'The file of input_cell_change.txt EXITS'
      endif   
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the old kinetic energy, if specified
!
    if (check_kinetic_energy) then
      call calc_kinetic_energy(kinetic_energy_old)
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the reference cell sizes
!
    iunit=vacant_unit()
    open(iunit, file=file_name_wrk,status='old')
!
    read(iunit,'(a)') chara_read_line
!
    if (index(chara_read_line, 'default') > 0) then
      cell_size_ref(1:3)=config%calc%cell_change%cell_sizes_in_xml_file(1:3)
      unit_value=1.0d0
      if (i_verbose >= 1) then
        if (log_unit > 0) then 
          write(log_unit,*) 'INFO:cell_size_ref(default): set by those in the structure XML file'
        endif
      endif
    else
      read(chara_read_line,*,iostat=ierr) cell_size_ref(1:3), chara_unit_wrk
!
      if (ierr == 0) then
        select case(trim(adjustl(chara_unit_wrk))) 
          case ('angst', 'Angst', 'angstrom', 'Angstrom', 'A')
            chara_unit='angst'
          case ('au', 'a.u.')
            chara_unit='au'
          case default
            write(*,*)'ERROR(cell size reading):unit name = ', trim(chara_unit_wrk)
            stop
        end select
      else
        chara_unit='au'
        read(chara_read_line,*,iostat=ierr) cell_size_ref(1:3)
        if (ierr /=0) then
          write(*,*)'ERROR(cell size reading)'
          write(*,*)'line = ', trim(chara_read_line)
          stop
        endif
      endif
!
      if (trim(chara_unit) == 'angst') then 
        unit_value=angst
      else
        unit_value=1.0d0
      endif
      cell_size_ref(:)=cell_size_ref(:)/unit_value
!
      if (i_verbose >= 1) then
        if (log_unit > 0) then 
          write(log_unit,*) '  cell_size_ref: unit in the file = ', trim(chara_unit)
        endif
       endif
!
    endif
!
    if (i_verbose >= 1) then
      if (log_unit > 0) then 
        write(log_unit,*) '  cell_size_ref x (au)=',cell_size_ref(1)
        write(log_unit,*) '  cell_size_ref y (au)=',cell_size_ref(2)
        write(log_unit,*) '  cell_size_ref z (au)=',cell_size_ref(3)
      endif  
    endif  
!
!
    do j=1,3
      if (dabs(cell_size_ref(j)) < 1.0d1) then
       write(*,'(a)') 'ERROR(elses_md_cell_change):cell_size_ref'
       write(*,'(a,i10,f30.20)') 'ERROR:Too small size of cell lengths: j, cell_length (au) =', & 
&                                j, cell_size_ref(j)
       stop
      endif   
    enddo
!   
    cell_scale(:)=1.0d0
    line_count_max=1000000
    n=0
    do j=0,line_count_max
      if (j == line_count_max) then
        write(*,*)'ERROR(elses_md_cell_change)'
        write(*,*)' .. reaches line_count_max'
        stop
      endif   
      read(iunit,fmt=*, iostat=ierr) n, dddx, dddy, dddz
      if (ierr /= 0) exit
      if (n <= step_count) then 
        cell_scale(1)=dddx
        cell_scale(2)=dddy
        cell_scale(3)=dddz
      endif  
    enddo   
!    
    if (config%calc%mode == "cell_change_only") then
      if (config%calc%cell_change%max_num_iter < 0) then
        config%calc%cell_change%max_num_iter = n+1 
        itemdmx=config%calc%cell_change%max_num_iter
        if (i_verbose >=1) then
          if (log_unit > 0) then 
            write(log_unit,'(a,i10)')' INFO:Set the maximum number of iteration in cell-change-only mode:', & 
&                                     config%calc%cell_change%max_num_iter
            write(log_unit,'(a,i10)')' INFO:Set itemdmx (=config%calc%cell_change%max_num_iter):', itemdmx 
          endif  
        endif 
      endif   
    endif   
!
    close(iunit)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the cell sizes
!
    ax_old=ax
    ay_old=ay
    az_old=az

    ax=cell_size_ref(1)*cell_scale(1)
    ay=cell_size_ref(2)*cell_scale(2)
    az=cell_size_ref(3)*cell_scale(3)
!
    config%system%structure%unitcell%vectorA(1)=ax
    config%system%structure%unitcell%vectorB(2)=ay
    config%system%structure%unitcell%vectorC(3)=az
!
    if (i_verbose >=1) then
      if (log_unit > 0) then 
        write(log_unit,'(a,i10,3f15.8)')' INFO:cell size [scaled]=', step_count, cell_scale(1:3)
        write(log_unit,'(a,i10,3f15.8)')' INFO:cell size [au    ]=', step_count, ax, ay, az
      endif  
    endif   
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Detect the scale procedure
!
    flag_for_scale = .false.
    if ( abs(ax/ax_old-1.0d0) > 1.0d-6 ) flag_for_scale = .true.
    if ( abs(ay/ay_old-1.0d0) > 1.0d-6 ) flag_for_scale = .true.
    if ( abs(az/az_old-1.0d0) > 1.0d-6 ) flag_for_scale = .true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scale the velocity so as to conserve the kinetic energy
!
    if ( (velocity_scaling) .and. (flag_for_scale) ) then
      if (allocated(velx) .eqv. .false.) then
        write(*,*) 'Alloc. Error (elses_md_cell_change) : velx'
        stop
      endif   
      if (allocated(vely) .eqv. .false.) then
        write(*,*) 'Alloc. Error (elses_md_cell_change) : vely'
        stop
      endif   
      if (allocated(velz) .eqv. .false.) then
        write(*,*) 'Alloc. Error (elses_md_cell_change) : velz'
        stop
      endif   
      velx(:)=velx(:)/(ax/ax_old)
      vely(:)=vely(:)/(ay/ay_old)
      velz(:)=velz(:)/(az/az_old)
      if (log_unit > 0) then 
        write(log_unit,'(a,i10,6f15.8)')' INFO:Change cell size [au](new/old)=', & 
&                       step_count, ax, ay, az, ax_old, ay_old, az_old
      endif  
    endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the new kinetic energy, if specified
!
    if (check_kinetic_energy) then
      call calc_kinetic_energy(kinetic_energy)
      if (abs(kinetic_energy) < 1.0d-10) then
        write(*,*)' ERROR(elses_md_cell_change):kinetic energy=',kinetic_energy
        stop
      endif   
      if (i_verbose >= 1) then
        if (log_unit > 0) then 
           write(log_unit,'(a,i10,2f15.8)')' INFO:Change cell size:KE(new, old)=', & 
&                       step_count, kinetic_energy, kinetic_energy_old
        endif  
      endif  
      if (abs(kinetic_energy/kinetic_energy_old)-1.0d0 > 1.0d-10) then
        write(*,*)' ERROR(elses_md_cell_change):kinetic energy    =',kinetic_energy
        write(*,*)' ERROR(elses_md_cell_change):kinetic energy_old=',kinetic_energy_old
        stop
      endif   
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return, if the maximum cutoff is not modified
!
    if (.not. rescale_cutoff) then
      if (i_verbose >=1) then
        if (log_unit > 0) then 
          write(log_unit,'(a)')' Modification of cutoff is not carried out.'
        endif  
      endif  
      return
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the maximum cutoff 
!
    if (config%calc%mode == "cell_change_only") then
      if ( i_pbc_x*i_pbc_y*i_pbc_z == 1 ) then
        if (trim(config%calc%solver%scheme) == 'eigen') then
          r_cut_book_old=r_cut_book
          r_cut_book=min(ax,ay,az)/2.01d0
         if (log_unit > 0) write(log_unit,'(a,2f15.8)') & 
&                ' INFO:Modify the cutoff(cell_change_only, eigen-solver ):new,old [au] =', r_cut_book, r_cut_book_old
        endif   
      endif  
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modify the cutoff with the PBC 
!
    if ( i_pbc_x == 1 ) then 
      if (ax/2.01d0 < r_cut_book) then ! 0.01 is a tolerance factor
        r_cut_book_old=r_cut_book
        r_cut_book=min(r_cut_book,ax/2.01d0)
        if (log_unit > 0) write(log_unit,'(a,2f15.8)') & 
&                ' INFO:Modify the cutoff:new,old [au] =', r_cut_book, r_cut_book_old
      endif   
    endif   
!
    if ( i_pbc_y == 1 ) then 
      if (ay/2.01d0 < r_cut_book) then ! 0.01 is a tolerance factor
        r_cut_book_old=r_cut_book
        r_cut_book=min(r_cut_book,ay/2.01d0)
        if (log_unit > 0) write(log_unit,'(a,2f15.8)') & 
&                ' INFO:Modify the cutoff:new,old [au] =', r_cut_book, r_cut_book_old
      endif   
    endif   
!
    if ( i_pbc_z == 1 ) then 
      if (az/2.01d0 < r_cut_book) then ! 0.01 is a tolerance factor
        r_cut_book_old=r_cut_book
        r_cut_book=min(r_cut_book,az/2.01d0)
        if (log_unit > 0) write(log_unit,'(a,f15.8)')' Modify the cutoff [au] =', r_cut_book
        if (log_unit > 0) write(log_unit,'(a,2f15.8)') & 
&                ' INFO:Modify the cutoff:new,old [au] =', r_cut_book, r_cut_book_old
      endif   
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    config%system%cutoff_radius=r_cut_book
!
  end subroutine elses_md_cell_change
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cell change (OLD routine, 2011.12.29, T. Hoshi)
!
  subroutine elses_md_cell_change_old
!
    use elses_mod_sim_cell,   only : ax, ay, az
    use elses_mod_md_dat,     only : itemd
!
    implicit none
    logical :: file_exist
    integer :: line_count_max
    character(len=16) :: file_name
    character(len=32)  :: chara1
    integer :: iunit,ierr
    integer :: j, n
    integer :: imode
    real(DOUBLE_PRECISION) :: scale_x, scale_y, scale_z
    real(DOUBLE_PRECISION) :: ddd
!
!
    file_name='Cell.txt'
!
    if (i_verbose >= 1) then
       write(*,*)'@@ elses_md_cell_change:itemd=',itemd
    endif   
!
    inquire (file=trim(file_name), exist=file_exist)
    if (i_verbose >= 1) write(*,*)'file_exist=',file_exist
!
    if (file_exist .eqv. .false.) then
      if (i_verbose >= 1) write(*,*)'The file of Cell.txt does not exist'
      return
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Store the initial cell sizes at the first time
!
    if ( .not. allocated(cell_parameter_store)) then
       if (i_verbose >= 1) then
         write(*,*)'  cell sizes are stored as original'
         write(*,*)'   original ax [au]=',ax
         write(*,*)'   original ay [au]=',ay
         write(*,*)'   original az [au]=',az
       endif  
       allocate (cell_parameter_store(3),stat=ierr)
       cell_parameter_store(1)=ax
       cell_parameter_store(2)=ay
       cell_parameter_store(3)=az
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the scales from the file 
!
    iunit=vacant_unit()
    open(iunit, file=file_name,status='old')
!
    read(iunit,'(a)') chara1
    if (i_verbose >= 1) then
      write(*,*)'  info from file:',trim(chara1)
    endif  
    read(iunit,'(a)') chara1
!
    if (trim(chara1) == 'aaa') then
      write(*,*)'  info from file:mode=',trim(chara1)
      imode=1
    else
      write(*,*)'  info from file:mode=',trim(chara1)
      write(*,*)'  invalid mode'
      stop
    endif
!    
    if (imode == 1) then
!
      ddd=-1.0d0
      scale_x=ddd
      scale_y=ddd
      scale_z=ddd
      line_count_max=1000000
      do j=1,line_count_max
!        write(*,*)'j=',j
         if (j == line_count_max) then
           write(*,*)'ERROR(elses_md_cell_change)'
           write(*,*)' .. reaches line_count_max'
           stop
         endif   
         read(iunit,fmt=*, iostat=ierr) n,ddd
!        if (ierr > 0) then
!           write(*,*)'ERROR(elses_md_cell_change)'
!           write(*,*)'Error in reading the file'
!           stop
!        endif   
         if (ierr /= 0) exit
!        write(*,*)'n,ddd=',n,ddd
         if (n <= itemd) then 
           scale_x=ddd
           scale_y=ddd
           scale_z=ddd
         endif  
      enddo   
      if ((ddd < 0.0d0) .or. (scale_x < 0.0d0)) then
        write(*,*)'ERROR(elses_md_cell_change)'
        write(*,*)'  --> The scale is not set'
        write(*,*)'  --> Check the input file '
        stop
      endif   
!
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the cell sizes
!
    ax=cell_parameter_store(1)*scale_x
    ay=cell_parameter_store(2)*scale_y
    az=cell_parameter_store(3)*scale_z
!
    if (i_verbose >=1) then
      write(*,'(a,i10,3f15.8)')' Cell size [scaled]=', itemd,scale_x, scale_y, scale_z
      write(*,'(a,i10,3f15.8)')' Cell size [au    ]=', itemd,ax,ay,az
    endif   
!
    close(iunit)
!
  end subroutine elses_md_cell_change_old
!
end module M_md_cell_change
