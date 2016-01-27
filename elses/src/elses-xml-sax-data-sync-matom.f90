!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_sax_data_sync_matom
!
  private
  public :: struc_data_sync_matom

contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Synchronize the structure data among nodes
!
  subroutine struc_data_sync_matom
!
    use M_config,           only : config    !(CHANGED)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r2 !(routine)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_i1 !(routine)
    use M_lib_element_database, only : get_atomic_num, get_element_name !(routine)
!
    implicit none
    integer :: error_count
    integer :: i_verbose, ierr
    integer :: myrank_wrk, noa_wrk, nos_wrk
    integer :: n_item_in_array
    real(8), allocatable :: wrk_array_real(:,:)
    integer, allocatable :: wrk_array_int(:)
    integer :: j, i_show
!   logical :: data_included(4)  ! for velocity, force, population, population_guess
    integer :: log_unit
    integer :: atomic_num
!
!   logical, parameter :: data_dump_for_debug = .true.
    logical, parameter :: data_dump_for_debug = .false.
    character(len=8)   :: element_name_wrk
!
    i_verbose = config%option%verbose 
!
    i_show=5
!
    noa_wrk     = config%system%structure%natom
    nos_wrk     = config%system%structure%nelement
    myrank_wrk  = config%calc%distributed%myrank
    log_unit    = config%calc%distributed%log_unit
!
    if (trim(config%system%structure%read_mode) == 'default') then
      config%system%structure%read_mode='redundant' 
    endif   
!
    select case(trim(config%system%structure%read_mode))
      case ('redundant') 
        return 
      case ('root_only') 
      case ('split') 
      case default
        write(*,*)'ERROR(struc_data_sync):read_mode=',trim(config%system%structure%read_mode)
        stop 
    end select   
!
    if ( log_unit > 0 ) then
      write(log_unit,'(a,3i10)')'@@ struc_data_sync:myrank_wrk, noa_wrk,nos_wrk=', myrank_wrk, noa_wrk, nos_wrk
    endif  
!
    allocate(wrk_array_real(3,noa_wrk), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in alloc:wrk_array_real'
      stop
    endif   
!
    allocate(wrk_array_int(noa_wrk), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in alloc:wrk_array_int'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the element data among the nodes
!
    wrk_array_int(:)=0
!
!$omp parallel default(shared) &
!$omp& private (j, atomic_num, element_name_wrk) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      if (trim(config%system%structure%matom(j)%name) /= '') then
        element_name_wrk=trim(config%system%structure%matom(j)%name) 
        call get_atomic_num(element_name_wrk, atomic_num)  
      else
        atomic_num=0
      endif   
      wrk_array_int(j)=atomic_num
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then 
          write(log_unit,*)' j, atomic_num, name (befor) =', j, atomic_num, config%system%structure%matom(j)%name
        endif  
      endif  
    enddo  
!$omp end do
!$omp end parallel
!
    call mpi_wrapper_allreduce_i1(wrk_array_int)
!
!$omp parallel default(shared) &
!$omp& private (j, atomic_num, element_name_wrk) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      atomic_num=wrk_array_int(j)
      call get_element_name(atomic_num, element_name_wrk)
      config%system%structure%matom(j)%name=element_name_wrk
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then 
          write(log_unit,*)' j, atomic_num, name (after) =', j, atomic_num, config%system%structure%matom(j)%name
        endif  
      endif  
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then 
        do j=1,min(i_show,noa_wrk)
          write(log_unit,*)' j, name (after sync) =', j, config%system%structure%matom(j)%name
        enddo  
      endif  
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the position data among the nodes
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then
        do j=1,min(i_show,noa_wrk)
          write(log_unit,'(a,2i10,3f20.10)')' j,pos(befor sync)=',myrank_wrk, j, config%system%structure%matom(j)%position(1:3)
        enddo   
      endif  
    endif
!
    wrk_array_real(:,:)=0.0d0
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      wrk_array_real(1:3,j)=config%system%structure%matom(j)%position(1:3)
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then
        write(log_unit,*)' mpi_command (position) start :myrank=',myrank_wrk
      endif  
    endif  
!
    call mpi_wrapper_allreduce_r2(wrk_array_real)
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then
        write(log_unit,*)' mpi_command (position) ends  :myrank=',myrank_wrk
      endif  
    endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      config%system%structure%matom(j)%position(1:3)=wrk_array_real(1:3,j)
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then
        do j=1,min(i_show,noa_wrk)
          write(log_unit,'(a,2i10,3f20.10)')' j,pos(after sync)=',myrank_wrk, j, config%system%structure%matom(j)%position(1:3)
        enddo   
      endif  
    endif
!
!   return
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    deallocate(wrk_array_real, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in dealloc:wrk_array_real'
      stop
    endif   
!
    deallocate(wrk_array_int, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in dealloc:wrk_array_int'
      stop
    endif   
!
  end subroutine struc_data_sync_matom
!
!
end module M_sax_data_sync_matom
