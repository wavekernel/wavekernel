!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_sax_parser
!
  private
  public :: struc_load_sax

contains

  subroutine struc_load_sax(filename)
    use M_config,      only : config ! for config%option%verbose
    use flib_sax,      only : open_xmlfile, xml_parse, close_xmlfile !(routine)
    use flib_sax,      only : xml_t                   !(type) 
!
    use M_sax_counter, only : in_structure_first_c   => in_structure_first !(CHANGED)
    use M_sax_counter, only : counter_atom, counter_element                !(CHANGED)
    use M_sax_counter, only : begin_element_c => begin_element             !(routine)
    use M_sax_counter, only : end_element_c   => end_element               !(routine)
    use M_sax_counter, only : pcdata_chunk_c  => pcdata_chunk              !(routine)
!
    use M_sax_handler, only : in_structure_first  !(CHANGED)
    use M_sax_handler, only : counter_atom_p  => counter_atom  !(CHANGED)
    use M_sax_handler, only : cell_vector_counter !(CHANGED)
    use M_sax_handler, only : begin_element !(routine)
    use M_sax_handler, only : end_element   !(routine)
    use M_sax_handler, only : pcdata_chunk  !(routine)
    use M_sax_handler, only : xml_flag_init, xml_flag_check  !(routine)
!   use M_sax_data_sync, only : struc_data_sync              !(routine)
!
    use M_config, only : config             !(CHANGED in slave routine) 
!
!   use M_lib_dst_info,    only : log_unit  !(unchanged)
!   use M_wall_clock_time, only : get_system_clock_time, measure_clock_time_period !(routine)
!
    implicit none
    character(len=*), intent(in) :: filename
    type(xml_t) :: fxml ! XML file object (opaque)
    integer :: ierr
    integer :: i_verbose
    real(8) :: time_wrk, time_wrk_previous, time_period
    integer :: error_count
!
    character(len=128)           :: filename_header
    character(len=128)           :: filename_to_be_read
!
    logical :: read_at_this_node
    character(len=128)           :: filename_basic
    integer :: lenf
    logical :: file_exist
    character(len=256)           :: chara_wrk, chara_wrk1, chara_wrk2
    logical                      :: generate_basic_file
    integer                      :: myrank, nprocs, log_unit
!
    integer                      :: stg_index, max_stg_index ! stage index ( = 0, 1, 2.. )
    integer                      :: j_tot_ini, j_tot_fin, j_ini, j_fin, atom_id
    integer                      :: split_index
    integer                      :: j
!
    i_verbose = config%option%verbose 
    myrank    = config%calc%distributed%myrank
    nprocs    = config%calc%distributed%nprocs
    log_unit  = config%calc%distributed%log_unit
!
    filename_to_be_read=trim(filename)
!   write(*,*)'filename_to_be_read (a) =',trim(filename_to_be_read)
!
    j_ini     = 1
    j_fin     = config%system%structure%natom
    j_tot_ini = 1
    j_tot_fin = config%system%structure%natom
    max_stg_index = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call get_system_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
    if (trim(config%system%structure%read_mode) == 'default') then
      config%system%structure%read_mode='redundant'
    endif   
!
    if (i_verbose >= 0) then
      if (log_unit > 0) write(log_unit,*) '@@ structure_load_sax:read_mode=', config%system%structure%read_mode
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the header of filename 
!
    lenf = len_trim(filename) 
    if (filename(lenf-3:lenf) /= '.xml') then
      write(*,*) 'ERROR:structure XML filename=',trim(filename)
      stop
    endif
    filename_header=filename(1:lenf-4)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the 'basic' filename (filename_basic) for 'root_only' and 'split' mode
!
    generate_basic_file = .false. 
    if (trim(config%system%structure%read_mode) == 'root_only') generate_basic_file = .true.
    if (trim(config%system%structure%read_mode) == 'split')     generate_basic_file = .true.
!
    if ( generate_basic_file) then
      lenf = len_trim(filename) 
      if (filename(lenf-3:lenf) /= '.xml') then
        write(*,*) 'ERROR:structure XML filename=',trim(filename)
        stop
      endif
      filename_basic=''
      filename_basic(1:lenf-4)=filename(1:lenf-4)
      filename_basic(lenf-3:lenf+6)='_basic.xml'
      if (i_verbose >= 0) then
        if (log_unit > 0) write(log_unit,*) 'filename_basic=',trim(filename_basic)
      endif  
      inquire (file=trim(filename_basic), exist=file_exist)
      if (.not. file_exist) then
        write(*,*)'ERROR(struc_load_sax):file is missing:',trim(filename_basic)
        stop
      endif   
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the filename in the 'root_only' mode
!
    if (trim(config%system%structure%read_mode) == 'root_only') then
      if (config%calc%distributed%root_node) then
        filename_to_be_read=trim(filename)
      else   
        filename_to_be_read=trim(filename_basic)
      endif   
    endif
!   
    if (i_verbose >= 0) then
      if (log_unit > 0) write(log_unit,*) 'filename_to_be_read (b) =',trim(filename_to_be_read)
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the filename in the 'split' mode
!
    if (trim(config%system%structure%read_mode) == 'split') then
      if (.not. config%system%structure%split%set) then
        write(*,*)'ERROR: split%set =',config%system%structure%split%set
        stop
      endif   
!
      if (config%system%structure%split%number_of_files <= 1) then
        write(*,*)'ERROR: split%number_of_files =',config%system%structure%split%number_of_files
        stop
      endif   
!
      max_stg_index = ( config%system%structure%split%number_of_files - 1 ) / nprocs
!
    endif   
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call get_system_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the dummy variables before parsing
!
    config%system%structure%unitcell%set        = .false.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the number of elements and the number of atoms for allocation
!
!
    if ( (config%system%structure%natom > 0) .and. (config%system%structure%nelement > 0) ) then
      counter_atom    = config%system%structure%natom 
      counter_element = config%system%structure%nelement
      if (i_verbose >= 1) then
        if (log_unit > 0) then 
          write(log_unit,*) 'INFO: number of atoms found in the config file:   =',counter_atom
          write(log_unit,*) 'INO: number of elements found in the config file:=',counter_element
        endif  
      endif  
    else   
      if (trim(config%system%structure%read_mode) /= 'redundant') then
        write(*,*)'ERROR(struc_load_sax):read_mode=',trim(config%system%structure%read_mode) 
        write(*,*)'# If read_mode=roor_only, '
        write(*,*)'#  the number of atoms should be given in the config. XML file' 
        stop
      endif   
      counter_atom    = 0
      counter_element = 0
!
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,*) 'filename_to_be_read (c) =',trim(filename_to_be_read)
      endif  
!
      call open_xmlfile(trim(filename_to_be_read),fxml,iostat=ierr)
      if (ierr /=0) then
        write(*,*)'ERROR(structure_load_sax):myrank, open_xmlfile(1);',myrank, trim(filename_to_be_read)
        write(*,*)'ERROR(structure_load_sax): ierr=',ierr
        write(*,*)'ERROR:The structure XML file may be missing or have a trouble'
        stop
      endif
!
!     write(*,*)'filename_to_be_read (d) =',trim(filename_to_be_read)
!
      call get_system_clock_time(time_wrk)
      call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
      if (i_verbose >= 0) then 
        write(*,'(a,f15.10)')'TIME:open_xml_file= ', time_period
      endif  
      time_wrk_previous=time_wrk
!
      in_structure_first_c = .true.
      call xml_flag_init
      call xml_parse(fxml, begin_element_handler=begin_element_c, &
&                      end_element_handler=end_element_c, &
&                      pcdata_chunk_handler=pcdata_chunk_c )
!
      call close_xmlfile(fxml)
      error_count=0
      call xml_flag_check(error_count)
      if (i_verbose >= 1) then
        if (log_unit > 0) write(log_unit,'(a,i5)')'INFO:result for xml_flag_check:count=',error_count
      endif    
!
      if (counter_atom <= 0) then
        write(*,*) 'ERROR:struc_load_sax:after xml_parse_c:counter_atom=',counter_atom
        write(*,*) 'ERROR:The structure XML file may be wrong.'
        stop
      endif
!
      if (counter_element <= 0) then
        write(*,*) 'ERROR*struc_load_sax:after xml_parse_c:counter_element=',counter_element
        write(*,*) 'ERROR:The structure XML file may be wrong.'
        stop
      endif
!
      if (i_verbose >= 1) then
        if (log_unit > 0) then
          write(log_unit,*)'number of counted atoms   =', counter_atom
          write(log_unit,*)'number of counted element =', counter_element
        endif  
      endif  
!
      config%system%structure%natom    = counter_atom
!
      if (config%system%structure%nelement /= counter_element) then
        write(*,*) 'ERROR*struc_load_sax:after xml_parse_c:nelement, counter_element=', & 
&                     config%system%structure%nelement, counter_element
        write(*,*) 'ERROR:The structure XML file may be wrong.'
        stop
      endif   
!
      call get_system_clock_time(time_wrk)
      call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
      if (i_verbose > 0) then 
        if (log_unit > 0) write(log_unit,'(a,f15.10)')'TIME:XML parse for atom count= ', time_period
      endif
      time_wrk_previous=time_wrk
!
    endif  
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (counter_atom <= 0) then
      write(*,*) 'ERROR(struc_load_sax:befor alloc. vatom):counter_atom=',counter_atom
      write(*,*) 'ERROR:The structure XML file may be wrong.'
      stop
    endif
!   
    if (config%system%structure%use_matom) then
      allocate (config%system%structure%matom(counter_atom), stat=ierr )
      if (ierr /= 0) then
        stop 'Alloc error. matom (structure_load_sax)' 
      endif   
      config%system%structure%matom(:)%name                 = ""      ! Default setting
      do j=1,counter_atom
        config%system%structure%matom(j)%position(1:3)        = 0.0d0   ! Default setting
      enddo  
    endif  
!
    if (config%system%structure%use_vatom) then
      allocate (config%system%structure%vatom(counter_atom), stat=ierr )
      if (ierr /= 0) then
        stop 'Alloc error. vatom (structure_load_sax)' 
      endif   
      do j=1,counter_atom
        config%system%structure%vatom(j)%position(1:3)      = 0.0d0   ! Default setting
        config%system%structure%vatom(j)%velocity(1:3)      = 0.0d0   ! Default setting
        config%system%structure%vatom(j)%force(1:3)         = 0.0d0   ! Default setting
      enddo  
      config%system%structure%vatom(:)%name                 = ""      ! Default setting
      config%system%structure%vatom(:)%class                = ""      ! Default setting
      config%system%structure%vatom(:)%motion               = "free"  ! Default setting
      config%system%structure%vatom(:)%velocity_set         = .false. ! Default setting
      config%system%structure%vatom(:)%force_set            = .false. ! Default setting
      config%system%structure%vatom(:)%population_set       = .false. ! Default setting
      config%system%structure%vatom(:)%population_guess_set = .false. ! Default setting
      config%system%structure%vatom(:)%group_id             =  -1     ! Dummy   setting
      config%system%structure%vatom(:)%group_id_set         = .false. ! Default setting
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the structure data
!
    call get_system_clock_time(time_wrk)
    call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
    time_wrk_previous=time_wrk
!
    if (i_verbose > 0) then 
      if (log_unit > 0) write(log_unit,*)' root_node=', config%calc%distributed%root_node
    endif
!
    read_at_this_node = .true.
!
!   if (trim(config%system%structure%read_mode) == 'root_only') then
!     if (.not. config%calc%distributed%root_node) then
!        read_at_this_node = .false.
!     endif      
!   endif   
!
    if (read_at_this_node) then
!
      do stg_index=0, max_stg_index 
! 
        if (trim(config%system%structure%read_mode) == 'split') then
!
          split_index= myrank + nprocs * stg_index
! 
          if ( split_index < 0 ) then
            write(*,*)'ERROR(save_restart_xml_main):split_index=',split_index
            stop
          endif  
!
          if (stg_index /= 0) then
            if ( split_index > config%system%structure%split%number_of_files - 1 ) then 
              cycle 
            endif  
          endif   
!
          if ( split_index > config%system%structure%split%number_of_files - 1 ) then
            filename_to_be_read=trim(filename_basic)
!           write(*,*)'filename_to_be_read (e) =',trim(filename_to_be_read)
          else
            call a_x_b_divided_by_c_i8(split_index, config%system%structure%natom, & 
&                    config%system%structure%split%number_of_files, j_ini)
            j_ini = j_ini +1
            call a_x_b_divided_by_c_i8(split_index+1, config%system%structure%natom, &
&                    config%system%structure%split%number_of_files, j_fin)
            lenf=len_trim(filename) 
            if (filename(lenf-3:lenf) /= '.xml') then
              write(*,*) 'ERROR(sax-parser): filename=',trim(filename)
              stop
            endif
            filename_header=filename(1:lenf-4)
            write(chara_wrk1, '(i6.6)') split_index
            chara_wrk2=trim(filename_header)//'_'//trim(chara_wrk1)//'.xml'
            filename_to_be_read = trim(chara_wrk2)
!           write(*,*)'INFO-XML-SPLIT:restat file name=',trim(filename_to_be_read), j_ini, j_fin
          endif  
!
          if (trim(filename_to_be_read) == '') then
            write(*,*)'ERROR(struc_load_sax):empty file name for filename_to_be_read'
            stop
          endif   
!
        endif ! end for if (trim(config%system%structure%read_mode) == 'split') then
!
!       write(*,*)'filename_to_be_read (f) =',trim(filename_to_be_read)
!
        if (i_verbose > 0) then 
          if (log_unit > 0) write(log_unit,*)'INFO-XML-SPLIT:file to be read = :', &
&                 myrank, trim(filename_to_be_read)
        endif
!
        inquire (file=trim(filename_to_be_read), exist=file_exist)
        if (.not. file_exist) then
          write(*,*)'ERROR(struc_load_sax):file is missing:',myrank, trim(filename_to_be_read)
          stop
        endif   
!
        call open_xmlfile(trim(filename_to_be_read),fxml,iostat=ierr)
        if (ierr /=0) then
          write(*,*)'ERROR(structure_load_sax):myrank, open_xmlfile(2);',myrank, trim(filename_to_be_read)
          write(*,*)'ERROR(structure_load_sax): ierr=',ierr
          write(*,*)'ERROR:The structure XML file may be missing or have a trouble'
          stop
        endif
!
        call get_system_clock_time(time_wrk)
        call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
        if (i_verbose >= 0) then 
          write(*,'(a,f15.10)')'TIME:open_xml_file for str. data = ', time_period
        endif  
        time_wrk_previous=time_wrk
!
        in_structure_first  = .true.
        cell_vector_counter = 0
        counter_atom_p      = 0
!
        call xml_flag_init
        call xml_parse(fxml, begin_element_handler=begin_element, &
&                          end_element_handler=end_element, &
&                          pcdata_chunk_handler=pcdata_chunk )
!
        call close_xmlfile(fxml)
        error_count=0
        call xml_flag_check(error_count)
        if (i_verbose >= 1) then
          if (log_unit > 0) write(log_unit,*) 'INFO:result for xml_flag_check=',error_count
        endif    
!
      enddo ! end of do_loop with stg_index
!  
    endif  ! for if (read_at_this_node)
!  
    call get_system_clock_time(time_wrk)
    call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
    if (i_verbose >= 0) then
      write(*,'(a,i5)')'INFO:result for xml_flag_check:count=',error_count
    endif  
    time_wrk_previous=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the data
!
    if (read_at_this_node) then
!
      if (.not. config%system%structure%unitcell%set) then
        write(*,*) 'ERROR(struc_load_sax):unitcell vectors are not set'
        stop
      endif   
!
      if (trim(filename_to_be_read) == trim(filename)) then
        if (counter_atom_p /= config%system%structure%natom) then
          write(*,*) 'ERROR(struc_load_sax):wrong number of atoms:counter_atom_p=',counter_atom_p
          write(*,*) '# The number of atoms in the config. XML file may be wrong'
          stop
        endif  
      endif   
!
      if (trim(filename_to_be_read) == trim(filename_basic)) then
        if (counter_atom_p /= 0) then
          write(*,*) 'ERROR(struc_load_sax):wrong number of atoms:counter_atom_p=',counter_atom_p
          write(*,*) '# The "basic" structure XML file may be wrong'
          stop
        endif  
      endif   
!
    endif  
!
!   stop
!
  end subroutine struc_load_sax
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine get_system_clock_time(elapse_time)
!      ----> Get the system_clock time in sec.
!
     implicit none
     real(kind=8), intent(out) :: elapse_time
     integer :: count, rate
!
     call system_clock(count, rate)
     elapse_time=dble(count)/dble(rate)
!
   end subroutine get_system_clock_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine measure_clock_time_period(time_present, time_previous, time_period)
!      ----> Measure the clock-time period 
!             with the correction for the possible reset of the system clock
!
     implicit none
     real(kind=8), intent(in)  :: time_present, time_previous
     real(kind=8), intent(out) :: time_period
     integer :: count, rate, max
!
     time_period = time_present - time_previous
!
     if (time_period < 0) then
       call system_clock(count, rate, max)
       time_period = time_period + dble(max)/dble(rate)
     endif   
!
   end subroutine measure_clock_time_period
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine a_x_b_divided_by_c_i8(a, b, c, d)
!      integer(8) calculation of d = a*b/c 
!
   implicit none
   integer,           intent(in)  :: a, b, c
   integer,           intent(out) :: d
   integer(kind=8) :: a_i8, b_i8, c_i8, d_i8
!
   a_i8=a
   b_i8=b
   c_i8=c
   d_i8=a_i8*b_i8/c_i8
!
   d=d_i8
!
  end subroutine a_x_b_divided_by_c_i8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_sax_parser

