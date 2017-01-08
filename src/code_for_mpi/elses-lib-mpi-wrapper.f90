!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_mpi_wrapper !(TRUE ROUTINES)
!
  use M_lib_dst_info, only  : log_unit  
  use mpi
  implicit none
  logical, parameter :: measure_wtime     = .true.
  logical, parameter :: sum_up_wtime      = .true.
  logical, parameter :: barrier_everytime = .true.
!
  real(kind(1d0)), allocatable :: total_allreduce_time(:)
  real(kind(1d0)), allocatable :: total_barrier_time(:)
!
  real(kind(1d0)), allocatable :: time_origin(:)  ! used in get_wall_clock_time
!
  private
  public :: mpi_wrapper_initial
  public :: mpi_wrapper_final
!
  public :: mpi_wrapper_allreduce_r1
  public :: mpi_wrapper_allreduce_r2
  public :: mpi_wrapper_allreduce_r3
  public :: mpi_wrapper_allreduce_r4
  public :: mpi_wrapper_allreduce_r5
!
  public :: mpi_wrapper_allreduce_i2
  public :: mpi_wrapper_allreduce_i1
  public :: mpi_wrapper_allreduce_i1c
  public :: mpi_wrapper_allreduce_r0
!
  public :: mpi_wrapper_allreduce_minmax_r1
!
  public :: mpi_wrapper_barrier_time
  public :: mpi_wrapper_allreduce_chk
!
  public :: sum_up_wtime
  public :: total_allreduce_time
  public :: total_barrier_time
!
  public :: mpi_wrapper_wtime
  public :: get_wall_clock_time
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_wtime(time_data)
    implicit none
    real(8), intent(out) :: time_data
!
    time_data=mpi_wtime()
!
  end subroutine mpi_wrapper_wtime
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wall_clock_time(elapse_time)
!      ----> Get the elapse wall-clock time in sec.
!             from the 'time origin'
!
    implicit none
    real(8), intent(out) :: elapse_time
    real(8)              :: time_data
    integer              :: ierr
!
    time_data=mpi_wtime()
!
    if ( .not. allocated(time_origin)) then
      allocate(time_origin(1), stat=ierr)
      time_origin(1)=time_data
    endif
!
    elapse_time=time_data-time_origin(1)
!
  end subroutine get_wall_clock_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_initial(myrank, nprocs, mpi_is_active, true_routine)
!
    implicit none
    integer, intent(out) :: myrank
    integer, intent(out) :: nprocs
    logical, intent(out) :: mpi_is_active
    logical, intent(out) :: true_routine
!
    integer :: myrank_for_plot_info
    integer :: ierr
!
!
    myrank_for_plot_info=8
!
    true_routine = .true.
!
    call mpi_init(ierr)
    if (ierr .ne. 0) then
       stop 'ERROR:mpi_init'
    endif
!
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    if (ierr .ne. 0) then
       stop 'ERROR:mpi_size'
    endif
!
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
    if (ierr .ne. 0) then
       stop 'ERROR:mpi_rank'
    endif
!
    if (myrank < myrank_for_plot_info) then
      write(*,*)'INFO-MPI:mpi_wrapper_initial(true routine)'
      write(*,*)'INFO-MPI:myrank, nprocs=',myrank,nprocs
    endif  
!
    mpi_is_active = .true.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preparation for summing up the time
!
    if (sum_up_wtime) then
!
      allocate (total_allreduce_time(0:nprocs-1), stat=ierr)
      if( ierr .ne. 0 ) then
         write(6,*)'alloc. error!(total_allreduce_time):ierr=',ierr
         stop
      endif
      total_allreduce_time(:)=0.0d0
!
      allocate (total_barrier_time(0:nprocs-1), stat=ierr)
      if( ierr .ne. 0 ) then
         write(6,*)'alloc. error!(total_barrier_time):ierr=',ierr
         stop
      endif
      total_barrier_time(:)=0.0d0
!
    endif   
!    
  end subroutine mpi_wrapper_initial
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_final
!
    implicit none
    integer :: ierr
!
    call mpi_finalize(ierr)
    if (ierr .ne. 0) then
       stop 'ERROR:mpi_final'
    endif
!
    if (log_unit > 0 )then
      write(log_unit,'(a)') '..... ELSES ended successfully (with MPI)'
      write(*,'(a)')        '..... ELSES ended successfully (with MPI)'
    endif   
!
  end subroutine mpi_wrapper_final
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r5(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:,:,:,:)
    integer :: mat_size1, mat_size2, mat_size3, mat_size4, mat_size5, mat_size
    integer :: ierr
    real(8), allocatable :: loc_array(:,:,:,:,:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size1=size(wrk_array,1)
    mat_size2=size(wrk_array,2)
    mat_size3=size(wrk_array,3)
    mat_size4=size(wrk_array,4)
    mat_size5=size(wrk_array,5)
    mat_size=mat_size1*mat_size2*mat_size3*mat_size4*mat_size5
!
    allocate (loc_array(mat_size1, mat_size2, mat_size3, mat_size4, mat_size5),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (mpi_wrapper_allreduce_r5)'
    loc_array=0.0d0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r5'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_r5 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_r5)'
!
  end subroutine mpi_wrapper_allreduce_r5
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r4(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:,:,:)
    integer :: mat_size1, mat_size2, mat_size3, mat_size4, mat_size
    integer :: ierr
    real(8), allocatable :: loc_array(:,:,:,:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size1=size(wrk_array,1)
    mat_size2=size(wrk_array,2)
    mat_size3=size(wrk_array,3)
    mat_size4=size(wrk_array,4)
    mat_size=mat_size1*mat_size2*mat_size3*mat_size4
!
    allocate (loc_array(mat_size1, mat_size2, mat_size3, mat_size4),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (mpi_wrapper_allreduce_r4)'
    loc_array=0.0d0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r4'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_r4 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_r4)'
!
  end subroutine mpi_wrapper_allreduce_r4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r3(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:,:)
    integer :: mat_size1, mat_size2, mat_size3, mat_size
    integer :: ierr
    real(8), allocatable :: loc_array(:,:,:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size1=size(wrk_array,1)
    mat_size2=size(wrk_array,2)
    mat_size3=size(wrk_array,3)
    mat_size=mat_size1*mat_size2*mat_size3
!
    allocate (loc_array(mat_size1, mat_size2, mat_size3),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (mpi_wrapper_allreduce_r3)'
    loc_array=0.0d0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r3'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_r3 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_r3)'
!
  end subroutine mpi_wrapper_allreduce_r3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r2(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:,:)
    integer :: mat_size1, mat_size2, mat_size
    integer :: ierr
    real(8), allocatable :: loc_array(:,:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size1=size(wrk_array,1)
    mat_size2=size(wrk_array,2)
    mat_size=mat_size1*mat_size2
!
    allocate (loc_array(mat_size1, mat_size2),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (mpi_wrapper_allreduce_r2)'
    loc_array=0.0d0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r2'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_r2 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_r2)'
!
  end subroutine mpi_wrapper_allreduce_r2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r1(wrk_array)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:)
    integer :: mat_size
    integer :: ierr
    real(8), allocatable :: loc_array(:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size=size(wrk_array,1)
!
    allocate (loc_array(mat_size),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (mpi_wrapper_allreduce_r1)'
    loc_array=0.0d0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r1'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_r1 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_r1)'
!
  end subroutine mpi_wrapper_allreduce_r1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_i2(wrk_array)
!
    implicit none
    integer, intent(inout) :: wrk_array(:,:)
    integer :: mat_size1, mat_size2, mat_size
    integer :: ierr
    integer,   allocatable :: loc_array(:,:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size1=size(wrk_array,1)
    mat_size2=size(wrk_array,2)
    mat_size=mat_size1*mat_size2
!
    allocate (loc_array(mat_size1, mat_size2),stat=ierr)
    if (ierr /= 0) then 
      write(*,*) 'Abort:ERROR in alloc (mpi_wrapper_allreduce_i2):size=', mat_size1, mat_size2
      stop
    endif  
    loc_array=0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_i2'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_i2 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_i2)'
!
  end subroutine mpi_wrapper_allreduce_i2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_i1(wrk_array)
!
    implicit none
    integer, intent(inout) :: wrk_array(:)
    integer :: mat_size1
    integer :: ierr
    integer,   allocatable :: loc_array(:)
    real(8) :: time_wrk, time_wrk_prev
!
    mat_size1=size(wrk_array,1)
!
    allocate (loc_array(mat_size1),stat=ierr)
    if (ierr /= 0) then 
      write(*,*) 'Abort:ERROR in alloc (mpi_wrapper_allreduce_i1):size=',mat_size1
      stop
    endif  
    loc_array=0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(wrk_array,loc_array,mat_size1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_i1'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_i1 =',mat_size1, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_i1)'
!
  end subroutine mpi_wrapper_allreduce_i1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_i1c(array_in,array_out)  ! Conversion mode
!
    implicit none
    integer, intent(in)  :: array_in(:)
    integer, intent(out) :: array_out(:)
    integer :: mat_size1
    integer :: ierr
    integer :: j
    real(8) :: time_wrk, time_wrk_prev, time_check
!
    mat_size1=size(array_in,1)
!
    array_out(:)=0
!
!   if (myrank == 0) then
!     do j=1,mat_size1
!       write(99,'(a,2i15)') 'befor_mpi',j, array_in(j)
!     enddo
!   endif
!
!
    if (measure_wtime) then
      call mpi_wrapper_allreduce_chk(time_check)  ! Routine for checking 
      time_wrk_prev=mpi_wtime()
    endif
!
    call mpi_allreduce(array_in,array_out,mat_size1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_i1c'
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,f20.10)')    'TIME:mpi_allreduce_chk =',time_check
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_i1c =',mat_size1, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
!   if (myrank == 0) then
!     do j=1,mat_size1
!       write(99,'(a,2i15)') 'after_mpi',j, array_out(j)
!     enddo
!   endif
!
  end subroutine mpi_wrapper_allreduce_i1c
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_r0(wrk_variable)
!
    implicit none
    real(8), intent(inout) :: wrk_variable
    real(8) :: loc_array(1)
!
    loc_array(1)=wrk_variable
!
    call mpi_wrapper_allreduce_r1(loc_array)
!
    wrk_variable=loc_array(1)
!
  end subroutine mpi_wrapper_allreduce_r0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_minmax_r1(wrk_array,mode_chara)
!
    implicit none
    real(8), intent(inout) :: wrk_array(:)
    integer :: mat_size
    integer :: ierr, imode
    real(8), allocatable :: loc_array(:)
    character(len=*), intent(in) :: mode_chara
    real(8) :: time_wrk, time_wrk_prev
!
    imode=0
    if (mode_chara == 'max') imode=1
    if (mode_chara == 'min') imode=2
    if (imode == 0) then
      write(*,*)'ERROR(mpi_wrapper_allreduce_minmax_r1)'
      stop
    endif   
!
    mat_size=size(wrk_array,1)
!
    allocate (loc_array(mat_size),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (mpi_wrapper_allreduce_r1)'
    loc_array=0.0d0
!
    if (measure_wtime) then
      if (barrier_everytime) call mpi_wrapper_barrier_time_no_out
      time_wrk_prev=mpi_wtime()
    endif
!
    if (imode  == 1) then
      call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_max,mpi_comm_world,ierr)
      if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r1'
    endif  
!      
    if (imode  == 2) then
      call mpi_allreduce(wrk_array,loc_array,mat_size,mpi_real8,mpi_min,mpi_comm_world,ierr)
      if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_allreduce_r1'
    endif  
!
    if (measure_wtime) then
      time_wrk=mpi_wtime()
      if (log_unit > 0 )then
        write(log_unit,'(a,i15,f20.10)')'TIME:mpi_allreduce_minmax_r1 =',mat_size, time_wrk-time_wrk_prev
      endif
      if (sum_up_wtime) call add_allreduce_time(time_wrk-time_wrk_prev)
    endif
!
    wrk_array=loc_array
!
    deallocate (loc_array,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (mpi_wrapper_allreduce_r1)'
!
  end subroutine mpi_wrapper_allreduce_minmax_r1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_barrier_time(time_check)  ! Routine for checking 
!
    use M_config,    only : config !(unchanged)
    implicit none
    real(8), intent(out) :: time_check
    integer :: ierr
    real(8) :: time_wrk, time_wrk_prev
!
    time_wrk_prev=mpi_wtime()
!
    if (config%calc%distributed%use_mpi_barrier) then
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_barrier_time'
    endif  
!
    time_wrk=mpi_wtime()
    time_check=time_wrk-time_wrk_prev
!
    if (sum_up_wtime) call add_barrier_time(time_check)
!
  end subroutine mpi_wrapper_barrier_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_barrier_time_no_out  ! Routine for checking 
!
    use M_config,    only : config !(unchanged)
    implicit none
    integer :: ierr
    real(8) :: time_wrk, time_wrk_prev
!
    time_wrk_prev=mpi_wtime()
!
    if (config%calc%distributed%use_mpi_barrier) then
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if (ierr /= 0) stop 'Abort:ERROR in mpi_wrapper_barrier_time'
    endif  
!
    time_wrk=mpi_wtime()
!
    if (sum_up_wtime) call add_barrier_time(time_wrk-time_wrk_prev)
!
  end subroutine mpi_wrapper_barrier_time_no_out 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_wrapper_allreduce_chk(time_check)  ! Routine for checking 
!
    implicit none
    real(8), intent(out) :: time_check
    call mpi_wrapper_barrier_time(time_check)
!
  end subroutine mpi_wrapper_allreduce_chk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_allreduce_time(wtime)  ! local routine
!
    use M_lib_dst_info, only  : myrank, nprocs
    implicit none
    real(8), intent(in) :: wtime
    logical, parameter  :: debug_mode = .true.
    integer :: array_size
!
    if (debug_mode) then
!
      if ( (myrank < 0 ) .or. (myrank > nprocs-1) ) then
        write(*,*)'ERROR(add_allreduce_time):myrank, nprocs=',myrank, nprocs
        stop  
      endif   
!
      if (.not. allocated(total_allreduce_time)) then
        write(*,*)'ERROR(add_allreduce_time):not allocated: total_allreduce_time' 
        stop  
      endif   
!
      array_size=size(total_allreduce_time,1)
      if (myrank > array_size-1) then
        write(*,*)'ERROR(add_allreduce_time):myrank, array_size=',myrank, array_size
        stop  
      endif   
!
    endif   
!
!
    total_allreduce_time(myrank)=total_allreduce_time(myrank)+wtime
!
  end subroutine add_allreduce_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_barrier_time(wtime)  ! local routine
!
    use M_lib_dst_info, only  : myrank, nprocs
    implicit none
    real(8), intent(in) :: wtime
    logical, parameter  :: debug_mode = .true.
    integer :: array_size
!
    if (debug_mode) then
!
      if ( (myrank < 0 ) .or. (myrank > nprocs-1) ) then
        write(*,*)'ERROR(barrier_time):myrank, nprocs=',myrank, nprocs
        stop  
      endif   
!
      if (.not. allocated(total_barrier_time)) then
        write(*,*)'ERROR(add_barrier_time):not allocated: total_allreduce_time' 
        stop  
      endif   
!
      array_size=size(total_barrier_time,1)
      if (myrank > array_size-1) then
        write(*,*)'ERROR(add_barrier_time):myrank, array_size=',myrank, array_size
        stop  
      endif   
!
    endif   
!
!
    total_barrier_time(myrank)=total_barrier_time(myrank)+wtime
!
  end subroutine add_barrier_time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module M_lib_mpi_wrapper
