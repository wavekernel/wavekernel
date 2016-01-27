!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_dst
!
  use M_qm_domain,    only: i_verbose !(CHANGED!)
  use M_lib_dst_info, only: mpi_is_active, myrank, nprocs, root_node
  use M_00_v_info,    only: version_info
  implicit none
  logical, parameter :: print_on_all_nodes = .false.
! logical, parameter :: print_on_all_nodes = .true.
!           ! If true, data will be printed on the all nodes (FOR DEBUGGING)
!
! logical :: mpi_is_active !(CHANGED)
! integer :: myrank   ! node index (CHANGED!!) (DIFFERENT VALUE AMONG THE NODES)
! integer :: nprocs   ! number of nodes (CHANGED!!)
!
  logical :: print_on_this_node ! (CHANGED!!) (DIFFERENT VALUE AMONG THE NODES)
!             ! If true, data will be printed on this node.
!  
! integer, allocatable :: dst_atm_list(:)     !  (CHANGED!) (DIFFERENT VALUE AMONG THE NODES)
!                ! distributed atom list 
! integer, allocatable :: len_dst_atm_list(:) !  (CHANGED!) (DIFFERENT VALUE AMONG THE NODES)
!                ! len_dst_atm_list(1)        : length of distributed atom list
!                ! len_dst_atm_list(n) (n>1)  : not used now
! integer, allocatable :: dst_atm_list_rev(:)     !  (CHANGED!) (DIFFERENT VALUE AMONG THE NODES)
!                ! reverse list for the distributed atoms
!
! integer, allocatable :: jsv4jsk_dst(:,:)     !  (CHANGED!) (DIFFERENT VALUE AMONG THE NODES)
! integer, allocatable :: noak_dst(:)          !  (CHANGED!) (DIFFERENT VALUE AMONG THE NODES)
! real(8), allocatable :: rcut_kry_dst(:)      !  (CHANGED!) (DIFFERENT VALUE AMONG THE NODES)
!
  private
  public :: mpi_is_active, myrank, nprocs
  public :: print_on_this_node
!
! public :: dst_atm_list, len_dst_atm_list, dst_atm_list_rev
!
  public :: set_dst_initial, set_dst_final
  public :: show_mpi_info
!
! public :: jsv4jsk_dst, noak_dst, rcut_kry_dst
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_dst_initial
!
   use M_config !(unchanged)
   use M_lib_dst_info,    only : true_mpi_wrapper    ! (CHANGED!)
   use M_lib_mpi_wrapper, only : mpi_wrapper_initial ! (routine) 
   use M_qm_domain,       only : i_verbose
   implicit none
   logical :: true_routine
!
!  if (i_verbose >= 0) then 
!    write(*,*)'@@ set_dst_initial'
!    write(*,*)'   config%option%verbose =', config%option%verbose
!  endif   
!
   call mpi_wrapper_initial(myrank, nprocs, mpi_is_active, true_mpi_wrapper)
!
   if (myrank == 0) then
     root_node = .true. 
   else
     root_node = .false. 
   endif
!   
   config%calc%distributed%mpi_is_active = mpi_is_active
   config%calc%distributed%myrank        = myrank
   config%calc%distributed%nprocs        = nprocs
   config%calc%distributed%root_node     = root_node
!
   if (.not. mpi_is_active) then
     print_on_this_node = .true.
     i_verbose = config%option%verbose
     return
   endif
!
   if (print_on_all_nodes) then
     print_on_this_node = .true.
     i_verbose = config%option%verbose
   else
     if (myrank == 0) then 
       print_on_this_node = .true.
       i_verbose = config%option%verbose
     else
       print_on_this_node = .false.
       config%option%verbose = -1
       i_verbose = config%option%verbose
     endif   
   endif
!
!
!  write(*,*)'INITIAL:config%calc%distributed%myrank   =',config%calc%distributed%myrank
!  write(*,*)'INITIAL:config%calc%distributed%root_node=',config%calc%distributed%root_node
!
!  if (i_verbose >= 0) then 
!    write(*,*)  '@@ ELESE Project'
!    write(*,*)  trim(version_info)
!    write(*,*)  'INFO-MPI:   mpi_is_active or not : ', mpi_is_active
!    write(*,*)  'INFO-MPI:   nprocs               =' , nprocs
!    write(*,*)  'INFO-MPI:   print_on_all_nodes   =', print_on_all_nodes
!    if (true_routine) then
!      write(*,*)'INFO-MPI:mpi_wrapper:true routine'
!    else  
!      write(*,*)'INFO-MPI:mpi_wrapper:dummy routine'
!    endif   
!  endif   
!
  end subroutine set_dst_initial
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_dst_final
!
   use M_lib_mpi_wrapper, only : mpi_wrapper_final
   implicit none
!
!  if (i_verbose >= 0) then 
!    write(*,*)'@@ set_dst_final:mpi_is_active=',mpi_is_active
!  endif  
!
   call mpi_wrapper_final
!     
!  if (i_verbose >= 0) then 
!    write(*,*)' .... ended (set_mpi_final)'
!  endif  
!
  end subroutine set_dst_final
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine show_mpi_info
!
   use M_lib_dst_info, only: log_unit !(unchanged)
   implicit none
!
!  if (i_verbose >= 0) then
!    write(*,*)'INFO-MPI:   mpi_is_active or not : ',mpi_is_active
!    write(*,*)'INFO-MPI:   myrank, nprocs       =', myrank, nprocs
!  endif  
!
   if (log_unit > 0) then
     write(log_unit,*)'INFO-MPI:   mpi_is_active or not : ',mpi_is_active
     write(log_unit,*)'INFO-MPI:   myrank, nprocs       =', myrank, nprocs
   endif
!
  end subroutine show_mpi_info
!   
end module M_md_dst
