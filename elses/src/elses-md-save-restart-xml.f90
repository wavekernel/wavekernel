!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module M_save_restart_xml
!
  use elses_mod_ctrl,     only : i_verbose !(unchanged)
  use M_io_dst_write_log, only : log_unit  !(unchanged)
  use M_md_dst,          only : mpi_is_active, myrank, nprocs !(unchanged)
  private
  public :: save_restart_xml

contains

  subroutine save_restart_xml(final_iteration)
    use M_config,           only : config
    use elses_mod_md_dat,   only : itemd, itemdmx
    implicit none
    logical, intent(in) :: final_iteration
    integer :: lenf
    logical :: flag_for_save
!   logical :: append_mode
    logical :: final_write
    character(len=1024) :: filename_wrk
    character(len=1024) :: filename_basic
    character(len=1024) :: chara_wrk
!
    if (i_verbose >= 1) then
!     write(*,*)'@@ save_restart_xml'
      if (log_unit > 0) write(log_unit,*)'@@ save_restart_xml'
      if (mpi_is_active) then
!       write(*,*)'myrank=',myrank
        if (log_unit > 0) write(log_unit,*)'myrank=',myrank
      endif
    endif 
!
    if (config%output%restart%interval == 0 ) then
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,*)'....is skipped, since (interval)=0'
      endif 
      return
    endif   
!
    lenf = len_trim(config%output%restart%filename)
!
    flag_for_save=.false.
    final_write  =.false.
    if( config%output%restart%filename /= "") then
      if ( mod(itemd-1,config%output%restart%interval) == 0 ) then
         flag_for_save=.true.
      endif   
      if ( final_iteration ) then
         flag_for_save=.true.
         final_write  =.true.
      endif   
    endif   
!
    if (i_verbose >= 1) then
      if (log_unit > 0) write(log_unit,*)'flag_for_save=',flag_for_save
    endif  
!
    filename_basic=trim(adjustL(config%output%restart%filename))
    if (config%output%restart%mode == "sequential" ) then
      write(chara_wrk, '(i8.8)') config%system%structure%mdstep
      filename_basic=filename_basic(1:lenf-4)//'_step'//trim(adjustl(chara_wrk))//'.xml'
      filename_basic=trim(adjustL(filename_basic))
    endif
!
    if (log_unit > 0) write(log_unit,*)'filename_restart=', trim(adjustL(filename_basic))
!
    if( flag_for_save ) then
       call save_restart_xml_main( config%system%structure, &
&                 filename_basic, final_write)
    end if
!
!
  end subroutine save_restart_xml
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine save_restart_xml_main( structure, filename, final_write )
    use elses_mod_file_io, only : vacant_unit  !(function)
    use M_config,           only : config
!        !! MAY BE CHANGED : config%output%restart%first_write
    use M_structure,   only : structure_type !(type)
    use M_structure,   only : atom_type !(type)
    use M_qm_domain,   only :e_num_on_atom !(unchanged)
    use M_md_dst_get_atom_range, only : a_x_b_divided_by_c_i8 !(function)
    implicit none
    type(structure_type), intent(in) :: structure
    character(len=*),     intent(in) :: filename
    logical,              intent(in) :: final_write
    logical, parameter :: total_number_of_atoms_is_added = .true.
    type(atom_type),         pointer :: atom
    integer :: fd
    integer :: ierr
    integer :: j
    integer :: stg_index, max_stg_index ! stage index ( = 0, 1, 2.. )
    real(8) :: population_value, population_guess_value
!
    character(len=64)   :: filename_to_be_saved
    character(len=64)   :: filename_header
    integer :: lenf
    character(len=64)   :: chara_wrk1, chara_wrk2, chara_wrk3, chara_wrk4, chara_wrk5
    character(len=1024) :: chara_wrk_long
    integer             :: j_ini, j_fin, atom_id, group_id_wrk
    integer             :: j_tot_ini, j_tot_fin
    integer             :: split_index ! ( = 0, 1, 2..., number_of_split_files -1)
!
    fd=vacant_unit()
!
    filename_to_be_saved = filename
    j_ini     = 1
    j_fin     = structure%natom
    j_tot_ini = 1
    j_tot_fin = structure%natom
    max_stg_index = 0
!
    if (i_verbose >= 1) then
      if (log_unit > 0) write(log_unit,*)'flag_for_split_save=', & 
&         config%output%restart%split    
    endif  
!
    if (config%output%restart%split) then
      config%output%restart%atom_id_is_added = .true.
      if (config%output%restart%number_of_split_files <= 0) then
        config%output%restart%number_of_split_files = nprocs 
      endif
      if (config%output%restart%number_of_split_files > structure%natom) then
        write(*,*)'ERROR(save_restart_xml_main):number_of_split_files, natom=', & 
&                      config%output%restart%number_of_split_files, structure%natom
        stop
      endif
      if (nprocs <= 0) then
        write(*,*)'ERROR(save_restart_xml_main):nprocs=',nprocs
        stop
      endif   
      max_stg_index = ( config%output%restart%number_of_split_files - 1 ) / nprocs
      if (log_unit > 0) write(log_unit,*) 'split_output:max_stg_index=',max_stg_index
    endif
!
    if (i_verbose >= 1) then
      if (log_unit > 0) write(log_unit,*)'max_stg_index=', max_stg_index
    endif  
!
    do stg_index=0, max_stg_index
!
      split_index= myrank + nprocs * stg_index
!
      if ( split_index < 0 ) then
        write(*,*)'ERROR(save_restart_xml_main):split_index=',split_index 
        stop
      endif   
!
      if (config%output%restart%split) then
        if ( split_index > config%output%restart%number_of_split_files-1 ) then 
          if (i_verbose >= 1) then
            if (log_unit > 0) write(log_unit,*)'NOTE:split save is not carried out at this node'
          endif  
          cycle
        endif  
      endif
!
      if (config%output%restart%split) then
        call a_x_b_divided_by_c_i8(split_index, structure%natom , config%output%restart%number_of_split_files, j_ini)
        j_ini = j_ini +1
        call a_x_b_divided_by_c_i8(split_index+1, structure%natom , config%output%restart%number_of_split_files, j_fin)
        lenf=len_trim(filename) 
        if (filename(lenf-3:lenf) /= '.xml') then
          write(*,*) 'ERROR(save_restart_xml_main): filename=',trim(filename)
          stop
        endif
        filename_header=filename(1:lenf-4)
        write(chara_wrk1, '(i6.6)') split_index
!       write(*,*)'INFO-XML-SPLIT:chara_wrk1=',trim(chara_wrk1)
        chara_wrk2=trim(filename_header)//'_'//trim(chara_wrk1)//'.xml'
        filename_to_be_saved = trim(chara_wrk2)
        if ((i_verbose >= 1) .and. (log_unit > 0)) then 
          write(log_unit,*) 'INFO-XML-SPLIT:restat file name=',trim(filename_to_be_saved), j_ini, j_fin
        endif  
      endif   
!
      if (config%output%restart%first_write) then
        open(fd,file=filename_to_be_saved, iostat=ierr)
        if (ierr /= 0) then
          write(*,*)'ERROR(save_restart_xml_main):open:ierr=',ierr
          stop
        endif   
        write(fd,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
        if (config%output%restart%append_mode=='on') then
          write(fd,'(a)') '<structure_history>'
        endif   
        close(fd)
      endif
!   
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'   save_restart_xml:append_mode=',config%output%restart%append_mode
      endif 
!
!
      if (config%output%restart%append_mode=='on') then
        open(fd,file=filename_to_be_saved,position='append', iostat=ierr)
        if (ierr /= 0) then
          write(*,*)'ERROR(save_restart_xml_main):open:ierr=',ierr
          stop
        endif   
      else
        open(fd,file=filename_to_be_saved, iostat=ierr)
        if (ierr /= 0) then
          write(*,*)'ERROR(save_restart_xml_main):open:ierr=',ierr
          stop
        endif   
        write(fd,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
      endif
!
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'Write name tag in file'
      endif 
!
      write(fd,'(a,a,a,I20,a)') '<structure name="', trim(structure%name), &
&        '" mdstep="', structure%mdstep, '">'
!
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'Write unitcell tag in file'
      endif 
!
      if (total_number_of_atoms_is_added) then
         write(fd,*) ""
         write(chara_wrk1, *) structure%natom
         chara_wrk1=adjustl(chara_wrk1)
         chara_wrk_long = ' <total_number_of_atoms> '//trim(chara_wrk1)//' </total_number_of_atoms>'
         write(fd,'(a)') trim(chara_wrk_long)
         write(fd,*) ""
      endif
!
      call unitcell_save( fd, structure%unitcell )
!
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'Write heatbath tag in file'
      endif  
!
      call heatbath_save( fd, structure%heatbath )
!
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'Write atom tag in file'
      endif  
!
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'Write atom tag in file'
      endif  
!
      if (config%output%restart%split) then
        write(chara_wrk1, *) split_index 
        write(chara_wrk2, *) config%output%restart%number_of_split_files
        write(chara_wrk3, *) j_fin - j_ini + 1
        write(chara_wrk4, *) j_ini
        write(chara_wrk5, *) j_fin
        write(fd,'(15a)') ' <split file_index="', trim(adjustl(chara_wrk1)), '" ',  &
&                         'number_of_files="', trim(adjustl(chara_wrk2)), '" ',  &
&                      'atoms_in_this_file="', trim(adjustl(chara_wrk3)), '" ',  &
&                           'atom_initial="', trim(adjustl(chara_wrk4)), '" ',  &
&                             'atom_final="', trim(adjustl(chara_wrk5)), '" />'
        write(fd,'(a)') ''
      endif
!   
      do j=j_ini, j_fin
        if (i_verbose >= 1) then 
          if (j < 5) then
            if ((i_verbose >= 1) .and. (log_unit > 0)) write(log_unit,*)'atom id=',j
          endif  
        endif
        atom => config%system%structure%vatom(j) 
        if (atom%population_guess_set) then
          population_guess_value=atom%population_guess
          if ((population_guess_value < -100d0) .or. (population_guess_value > 100.0d0)) then
            write(*,*)'Error: Unphysical population on atom ?:j,population_guess=',j,population_guess_value
            stop
          endif   
        endif   
        if (atom%population_set) then
          population_value=atom%population
          if ((population_value < -100d0) .or. (population_value > 100.0d0)) then
            write(*,*)'Error: Unphysical population on atom ?:j,population_guess=',j,population_value
            stop
          endif   
        endif   
        if (config%output%restart%atom_id_is_added) then
          atom_id=j 
          group_id_wrk=atom%group_id
        else  
          atom_id=-1 
          group_id_wrk=-1
        endif   
        call atom_save( fd, structure%vatom(j), structure%unitcell, &
&                       atom%population_set, atom%population_guess_set, &
&                      atom%population, atom%population_guess, atom_id, group_id_wrk)
      end do
!
      if ((i_verbose >= 1) .and. (log_unit > 0)) then 
        write(log_unit,*)'Write atom tag in file...ends'
      endif  
!
      write(fd,*) '</structure>'
!
      if (config%output%restart%append_mode=='on') then
        if (final_write) write(fd,'(a)') '</structure_history>'
      endif  
!
      close(fd)
!
    enddo  ! loop end for stg_index
!
    config%output%restart%first_write = .false.
!  
    return
  end subroutine save_restart_xml_main

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine unitcell_save( fd, unitcell )
    use M_structure,   only : unitcell_type !(type)
    implicit none
    integer, intent(in) :: fd
    type(unitcell_type), intent(in) :: unitcell

    write(fd,*) ""
    write(fd,*) "<unitcell>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', unitcell%vectorA(1:3), "</vector>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', unitcell%vectorB(1:3), "</vector>"
    write(fd,'(a23,3e23.15,a9)') '  <vector unit="a.u.">', unitcell%vectorC(1:3), "</vector>"
    write(fd,*) "</unitcell>"
    write(fd,*) ""

    return
  end subroutine unitcell_save

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine heatbath_save( fd, heatbath )
    use M_structure,   only : heatbath_type !(type)
    implicit none
    integer, intent(in) :: fd
    type(heatbath_type), intent(in) :: heatbath

    write(fd,*) ""
    write(fd,*) "<heatbath>"
    write(fd,'(a27,1e23.15,a14)') '  <massperatom unit="a.u.">', heatbath%massperatom, "</massperatom>"
    write(fd,'(a24,1e23.15,a14)') '  <position unit="a.u.">', heatbath%position, "</position>"
    write(fd,'(a24,1e23.15,a14)') '  <velocity unit="a.u.">', heatbath%velocity, "</velocity>"
    write(fd,*) "</heatbath>"
    write(fd,*) ""

    return
  end subroutine heatbath_save

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine atom_save( fd, atom, unitcell, population_set, population_guess_set, & 
&                       population, population_guess, atom_id, group_id)
    use M_structure,   only : atom_type, unitcell_type !(type)
    implicit none
    integer, intent(in) :: fd
    type(atom_type), intent(in) :: atom
    type(unitcell_type), intent(in) :: unitcell
    logical,  intent(in) :: population_set
    logical,  intent(in) :: population_guess_set
    real(8),  intent(in) :: population
    real(8),  intent(in) :: population_guess
    integer,  intent(in) :: atom_id    ! NOTE: 'atom_id = -1' means the atom_id is not added.
    integer,  intent(in) :: group_id   ! NOTE: 'atom_id = -1' means the atom_id is not added.
!   integer :: k
    real(8) :: a, b, c
    real(8) :: la, lb, lc
    character(len=64)    :: chara_wrk
    character(len=64)    :: chara_wrk_class
    character(len=64)    :: chara_wrk_atom_id
    character(len=64)    :: chara_wrk_group_id
    character(len=64)    :: chara_sep
!
    la = dsqrt(dot_product(unitcell%vectorA,unitcell%vectorA))
    lb = dsqrt(dot_product(unitcell%vectorB,unitcell%vectorB))
    lc = dsqrt(dot_product(unitcell%vectorC,unitcell%vectorC))
!
    a = atom%position(1) / la
    b = atom%position(2) / lb
    c = atom%position(3) / lc
!
!   a = atom%position(1) / unitcell%vectorA(1)
!   b = atom%position(2) / unitcell%vectorB(2)
!   c = atom%position(3) / unitcell%vectorC(3)
!
    chara_sep='"'
!
    if (trim(atom%class) == '') then
      chara_wrk_class='' 
    else  
      chara_wrk=trim(adjustl(atom%class)) 
      chara_wrk_class='class="'//trim(chara_wrk)//chara_sep
    endif  
!
    if (atom_id == -1) then
      chara_wrk_atom_id='' 
    else
      write(chara_wrk,*) atom_id
      chara_wrk_atom_id='id="'//trim(adjustl(chara_wrk))//chara_sep
    endif  
!
    if (group_id == -1) then
      chara_wrk_group_id='' 
    else
      write(chara_wrk,*) group_id
      chara_wrk_group_id='group_id="'//trim(adjustl(chara_wrk))//chara_sep
    endif  
!
    write(fd,'(a,a,a,a,a,a,a,a,a,a,a,a)') ' <atom element="', trim(atom%name), '" ', &
&        trim(adjustl(chara_wrk_class)), ' ', trim(adjustl(chara_wrk_atom_id)), ' ', trim(adjustl(chara_wrk_group_id)), ' ', &
&        'motion="', trim(atom%motion), '">'
!      
!   if (atom_id == -1) then
!     write(fd,*) '<atom element="', trim(atom%name), '" ', &
!&        'class="', trim(atom%class), '" ', &
!&        'motion="', trim(atom%motion), '">'
!    else
!     write(chara_wrk,*) atom_id
!     write(fd,'(a,a,a,a,a,a,a,a,a,a,a,a)') '<atom element="', trim(atom%name), '" ', &
!&        'id="', trim(adjustl(chara_wrk)), '" ', 'class="', trim(atom%class), '" ', &
!&        'motion="', trim(atom%motion), '">'
!    endif   
!
    write(fd,'(a28,3e23.15,a11)') '  <position unit="internal">', a, b, c, '</position>'
    write(fd,'(a24,3e23.15,a11)') '  <velocity unit="a.u.">', atom%velocity(1:3), '</velocity>'
    write(fd,'(a24,3e23.15,a11)') '  <force    unit="a.u.">', atom%force(1:3), '</force>'

    if( population_set ) then
       write(fd,'(a,f30.20,a)') '  <population>  ', population, '  </population>'
    endif
   
    if( population_guess_set ) then
       write(fd,'(a,f30.20,a)') '  <population_guess>  ', population_guess, '  </population_guess>'
    endif

!   if( atom%ncustumize > 0 ) then
!      do k=1, atom%ncustumize
!         write(fd,'(a19,4f16.10,a10,a16,a4)') &
!              '  <custumize base="', atom%vcustumize(k)%base(1:4), &
!              '" status="',  trim(atom%vcustumize(k)%status), &
!              '" />'
!      end do
!   end if

    write(fd,'(a)') ' </atom>'
    write(fd,*) ""

    return
  end subroutine atom_save

end module M_save_restart_xml



