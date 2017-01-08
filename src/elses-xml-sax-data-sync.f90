!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_sax_data_sync
!
  private
  public :: struc_data_sync

contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Synchronize the structure data among nodes
!
  subroutine struc_data_sync
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
    logical :: data_included(4)  ! for velocity, force, population, population_guess
    integer :: log_unit
    integer :: atomic_num
!
!   logical, parameter :: data_dump_for_debug = .true.
    logical, parameter :: data_dump_for_debug = .false.
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
!$omp& private (j, atomic_num) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      if (trim(config%system%structure%vatom(j)%name) /= '') then
        call get_atomic_num(config%system%structure%vatom(j)%name, atomic_num)  
      else
        atomic_num=0
      endif   
      wrk_array_int(j)=atomic_num
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then 
          write(log_unit,*)' j, atomic_num, name (befor) =', j, atomic_num, config%system%structure%vatom(j)%name
        endif  
      endif  
    enddo  
!$omp end do
!$omp end parallel
!
    call mpi_wrapper_allreduce_i1(wrk_array_int)
!
!$omp parallel default(shared) &
!$omp& private (j, atomic_num) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      atomic_num=wrk_array_int(j)
      call get_element_name(atomic_num, config%system%structure%vatom(j)%name)
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then 
          write(log_unit,*)' j, atomic_num, name (after) =', j, atomic_num, config%system%structure%vatom(j)%name
        endif  
      endif  
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then 
        do j=1,min(i_show,noa_wrk)
          write(log_unit,*)' j, name (after sync) =', j, config%system%structure%vatom(j)%name
        enddo  
      endif  
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the motion-status data among the nodes
!
    wrk_array_int(:)=0
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then 
        do j=1,min(i_show,noa_wrk)
          write(log_unit,*)' j, motion (befor sync) =', j, config%system%structure%vatom(j)%motion
        enddo  
      endif  
    endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      if (trim(config%system%structure%vatom(j)%motion) == 'free') then
        wrk_array_int(j)=0
      else
        wrk_array_int(j)=1
      endif   
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then 
        do j=1,min(i_show,noa_wrk)
          write(log_unit,*)' j, motion-array (befor mpi-allreduce) =', j, wrk_array_int(j)
        enddo  
      endif  
    endif  
!
    call mpi_wrapper_allreduce_i1(wrk_array_int)
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then 
        do j=1,min(i_show,noa_wrk)
          write(log_unit,*)' j, motion-array (after mpi-allreduce) =', j, wrk_array_int(j)
        enddo  
      endif  
    endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
    do j=1,noa_wrk
      if (wrk_array_int(j) == 0) then
        config%system%structure%vatom(j)%motion='free'
      else
        config%system%structure%vatom(j)%motion='fixed'
      endif  
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then 
        do j=1,min(i_show,noa_wrk)
          write(log_unit,*)' j, motion (after sync) =', j, config%system%structure%vatom(j)%motion
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
          write(log_unit,'(a,2i10,3f20.10)')' j,pos(befor sync)=',myrank_wrk, j, config%system%structure%vatom(j)%position(1:3)
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
      wrk_array_real(1:3,j)=config%system%structure%vatom(j)%position(1:3)
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
      config%system%structure%vatom(j)%position(1:3)=wrk_array_real(1:3,j)
    enddo  
!$omp end do
!$omp end parallel
!
    if (data_dump_for_debug) then
      if ( log_unit > 0 ) then
        do j=1,min(i_show,noa_wrk)
          write(log_unit,'(a,2i10,3f20.10)')' j,pos(after sync)=',myrank_wrk, j, config%system%structure%vatom(j)%position(1:3)
        enddo   
      endif  
    endif
!
!   return
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call inquire_data_kind(data_included)
!
    if ( log_unit > 0 ) then
      write(log_unit,*)'data_included:myrank, vel, force, popu, popu_guess)',myrank_wrk, data_included(1:4)
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the velocity data among the nodes, if exist
!
    if (data_included(1)) then
!
      if ( log_unit > 0 ) then
        write(log_unit,*)'INFO-XML-SAX:data sync: velocity '
      endif  
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,3f20.10)')' j,vel(befor sync)=',myrank_wrk, j, config%system%structure%vatom(j)%velocity(1:3)
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
        wrk_array_real(1:3,j)=config%system%structure%vatom(j)%velocity(1:3)
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (velocity) start :myrank=',myrank_wrk
        endif  
      endif  
!
      call mpi_wrapper_allreduce_r2(wrk_array_real)
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (velocity) ends  :myrank=',myrank_wrk
        endif  
      endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
      do j=1,noa_wrk
        config%system%structure%vatom(j)%velocity_set = .true.
        config%system%structure%vatom(j)%velocity(1:3)=wrk_array_real(1:3,j)
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,3f20.10)')' j,vel(after sync)=',myrank_wrk, j, config%system%structure%vatom(j)%velocity(1:3)
          enddo   
        endif  
      endif
!
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the force data among the nodes, if exist
!
    if (data_included(2)) then
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,3f20.10)')' j,force(befor sync)=',myrank_wrk, j, config%system%structure%vatom(j)%force(1:3)
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
        wrk_array_real(1:3,j)=config%system%structure%vatom(j)%force(1:3)
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (force) start :myrank=',myrank_wrk
        endif  
      endif  
!
      call mpi_wrapper_allreduce_r2(wrk_array_real)
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (force) ends  :myrank=',myrank_wrk
        endif  
      endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
      do j=1,noa_wrk
        config%system%structure%vatom(j)%force_set = .true.
        config%system%structure%vatom(j)%force(1:3)=wrk_array_real(1:3,j)
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,3f20.10)')' j,force(after sync)=',myrank_wrk, j, config%system%structure%vatom(j)%force(1:3)
          enddo   
        endif  
      endif
!
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    deallocate(wrk_array_real, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in dealloc:wrk_array_real'
      stop
    endif   
!
    allocate(wrk_array_real(1,noa_wrk), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in alloc:wrk_array_real'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the population data among the nodes, if exist
!
    if (data_included(3)) then
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,f20.10)')' j,population(befor sync)=',myrank_wrk, j, config%system%structure%vatom(j)%population
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
        wrk_array_real(1,j)=config%system%structure%vatom(j)%population
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (population) start :myrank=',myrank_wrk
        endif  
      endif  
!
      call mpi_wrapper_allreduce_r2(wrk_array_real)
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (populaton) ends  :myrank=',myrank_wrk
        endif  
      endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
      do j=1,noa_wrk
        config%system%structure%vatom(j)%population_set = .true.
        config%system%structure%vatom(j)%population     = wrk_array_real(1,j)
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,f20.10)')' j,force(after sync)=',myrank_wrk, j, config%system%structure%vatom(j)%population
          enddo   
        endif  
      endif
!
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sync the population_guess data among the nodes, if exist
!
    if (data_included(4)) then
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,f20.10)')' j,population(befor sync)=',& 
&               myrank_wrk, j, config%system%structure%vatom(j)%population_guess
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
        wrk_array_real(1,j)=config%system%structure%vatom(j)%population_guess
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (population) start :myrank=',myrank_wrk
        endif  
      endif  
!
      call mpi_wrapper_allreduce_r2(wrk_array_real)
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          write(log_unit,*)' mpi_command (populaton_guess) ends  :myrank=',myrank_wrk
        endif  
      endif  
!
!$omp parallel default(shared) &
!$omp& private (j) 
!$omp  do schedule(static)
      do j=1,noa_wrk
        config%system%structure%vatom(j)%population_guess_set = .true.
        config%system%structure%vatom(j)%population_guess     = wrk_array_real(1,j)
      enddo  
!$omp end do
!$omp end parallel
!
      if (data_dump_for_debug) then
        if ( log_unit > 0 ) then
          do j=1,min(i_show,noa_wrk)
            write(log_unit,'(a,2i10,f20.10)')' j,force(after sync)=', & 
&                myrank_wrk, j, config%system%structure%vatom(j)%population_guess
          enddo   
        endif  
      endif
!
    endif   
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
  end subroutine struc_data_sync
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine inquire_data_kind(data_included)
    use M_config,           only : config    !(unchanged)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_i1 !(routine)
    implicit none
    logical, intent(out) :: data_included(4) ! for velocity, force, population, population_guess
    integer              :: work_array(4)
    logical              :: wrong_mode
!
    wrong_mode = .true.
    if (trim(config%system%structure%read_mode) == 'root_only') wrong_mode = .false.
    if (trim(config%system%structure%read_mode) == 'split') wrong_mode = .false.

    if (wrong_mode) then
      write(*,*) 'ERROR(inquire_data_kind):mode=',trim(config%system%structure%read_mode)
      stop
    endif
!   
    work_array(:)=0
!
    if (config%calc%distributed%root_node) then
      if (config%system%structure%vatom(1)%velocity_set)         work_array(1)=1
      if (config%system%structure%vatom(1)%force_set)            work_array(2)=1
      if (config%system%structure%vatom(1)%population_set)       work_array(3)=1
      if (config%system%structure%vatom(1)%population_guess_set) work_array(4)=1
    endif
!
    call mpi_wrapper_allreduce_i1(work_array)
!
    data_included(:)= .false.
!
    if (work_array(1) /= 0) data_included(1) = .true.
    if (work_array(2) /= 0) data_included(2) = .true.
    if (work_array(3) /= 0) data_included(3) = .true.
    if (work_array(4) /= 0) data_included(4) = .true.
!
  end subroutine inquire_data_kind
!
end module M_sax_data_sync
