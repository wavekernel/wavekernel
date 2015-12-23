!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_cohp_dstm_plot
!
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_unit   !(unchanged)
  implicit none
  integer, allocatable :: file_unit_cohp(:)           ! set in 'prep_cohp_dst'
  real(DOUBLE_PRECISION), allocatable :: cohp_cri(:)  ! set in 'prep_cohp_dst'
!
  character(len=*), parameter :: filename_cohp_head = 'output-cohp-head.txt'
  character(len=*), parameter :: filename_cohp_body = 'output-cohp-body.txt'
  character(len=*), parameter :: filename_cohp_data = 'output-cohp-data-dump.txt'
!  
  private
!
! Public routines
  public prep_cohp_dst
  public calc_cohp_dstm_plot
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Preparation of COHP calculation 
!          in the DST work flow only at the root node
!
  subroutine prep_cohp_dst
!
    use M_md_dst,           only : myrank       !(unchanged)
    use M_qm_domain,        only : noav         !(unchanged)
    use elses_mod_file_io,  only : vacant_unit  !(function)
    use M_config,           only : config !(unchanged) 
!                        (only config%output%bond_list%set, config%output%bond_list%interval)
    use elses_mod_md_dat,    only : itemd  !(unchanged)
    implicit none
    integer :: ierr, iunit1, iunit2
    real(DOUBLE_PRECISION) :: cohp_cri_wrk
!
    if (allocated(file_unit_cohp)) return
    if ( .not. config%output%bond_list%set ) return
    if (mod(itemd-1,config%output%bond_list%interval) /= 0 ) return
!
    cohp_cri_wrk = -0.001d0 ! criteria for COHP in eV
!
    if (i_verbose >= 1) then
      write(*,*)'@@ prep_cohp_dst:cohp_criteria (eV)=', cohp_cri_wrk 
    endif   
!
!   allocate(file_unit_cohp(2), stat=ierr)
!
    allocate(file_unit_cohp(1), stat=ierr)
    if (ierr /= 0) stop 'Alloc Error (prep_cohp_dst)'
    file_unit_cohp(:)=-1
!
    allocate(cohp_cri(1), stat=ierr)
    if (ierr /= 0) stop 'Alloc Error (cohp_cri)'
    cohp_cri(1)=cohp_cri_wrk
!
    iunit1=vacant_unit()
    file_unit_cohp(1)=iunit1
    open(iunit1, file=filename_cohp_data, status='unknown')
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Write the file header only at the root node
!
    if (myrank == 0) then
      write(iunit1,'(a)') trim(filename_cohp_data)
      write(iunit1,'(i10,f20.10)') noav, cohp_cri_wrk
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   iunit2=vacant_unit()
!   file_unit_cohp(2)=iunit2
!   open(iunit2, file=filename_cohp_body, status='unknown')
!
!   write(*,*) 'INFO:COHP file unit1=',iunit1, iunit2
!
!   write(iunit1,'(a)') trim(filename_cohp_head)
!   write(iunit1,*) noav
!
!   write(iunit2,'(a)') trim(filename_cohp_body)
!
  end subroutine prep_cohp_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the partial trace
!       such as local energy for given atom
!
  subroutine calc_cohp_dstm_plot & 
&      (atm_index,orb_index,cohp_loc,jsv4jsk,booking_list_dstm,booking_list_dstm_len)
!    
    use M_config,           only : config !(unchanged) 
!                        (only config%output%bond_list%set, config%output%bond_list%interval)
!                        (only config%system%structure%mdstep)
    use elses_mod_md_dat,     only : itemd  !(unchanged)
!
!   use M_lib_get_core_index, only : get_core_index !(routine) 
!   use elses_mod_file_io,    only : vacant_unit    !(routine) 
    use M_lib_phys_const,     only : ev4au          !(unchanged)
    implicit none
    integer,                   intent(in)  :: atm_index, orb_index
    real(DOUBLE_PRECISION),    intent(in)  :: cohp_loc(:)
    integer,                   intent(in)  :: jsv4jsk(:)
    integer,                   intent(in)  :: booking_list_dstm(:,:)
    integer,                   intent(in)  :: booking_list_dstm_len(:)
!
    integer :: jsv2, jsd1, jsv1, jsk1, jsk2, ja2
    integer :: iunit1, iunit2
    integer :: step_count
    real(DOUBLE_PRECISION) :: cohp_in_ev, cohp_cri_wrk
    logical                :: cohp_cri_is_set
    logical                :: cohp_plot
!
!   character(len=32) :: file_header_name
!   character(len=32) :: file_number
!   character(len=70) :: file_name
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (.not. allocated(file_unit_cohp)) return
!
    if ( .not. config%output%bond_list%set ) return
    if (mod(itemd-1,config%output%bond_list%interval) /= 0 ) return
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    iunit1=file_unit_cohp(1)
!   iunit2=file_unit_cohp(2)
!
    if ((iunit1 <= 0) .or. (iunit1 >= 100)) then
      write(*,*)'ERROR(calc_cohp_dstm_plot):iunit1=',iunit1 
      stop
    endif
!   
!   if ((iunit2 <= 0) .or. (iunit2 >= 100)) then
!     write(*,*)'ERROR(calc_cohp_dstm_plot):iunit2=',iunit2
!     stop
!   endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    step_count=config%system%structure%mdstep
    jsv2=atm_index
    ja2 =orb_index
    jsk2=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    cohp_cri_is_set = .false.
    cohp_cri_wrk=1.0d10       ! dummy value
    if (allocated(cohp_cri)) then
      cohp_cri_is_set = .true.
      cohp_cri_wrk=cohp_cri(1)
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   if (ja2 == 1) then
!     write(iunit1,'(a,3i10)')'cohp-head ', step_count, jsv2, booking_list_dstm_len(jsk2)
!   endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   call get_core_index(core_index)
!   file_header_name='output-bond-split-'
!   write(file_number, '(i7.7)') core_index
!   file_name = trim(file_header_name)//trim(file_number)//'.txt'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Plot the COHP
!
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      jsv1=jsv4jsk(jsk1)
      cohp_in_ev=cohp_loc(jsd1)*ev4au
      cohp_plot = .true.
      if (cohp_cri_is_set) then
        if (cohp_in_ev > cohp_cri_wrk) cohp_plot= .false.
      endif   
      if (cohp_plot) then
         write(iunit1,'(a,4i10,f20.10)')'cohp-body ', step_count, jsv2, ja2, jsv1, cohp_in_ev
      endif  
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  end subroutine calc_cohp_dstm_plot
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_cohp_dstm_plot
