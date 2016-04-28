!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_main
!
  private
  public :: elses_md_main_01
!
  contains
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Main routine for MD
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine elses_md_main_01
!
    use elses_mod_md_dat, only : first_iteration, final_iteration   !(CHANGED)
!
    use M_config,          only : config !(unchanged)
    use elses_mod_sel_sys, only : c_system !(unchanged)
    use elses_mod_md_dat, only : itemdmx   !(unchanged)
    use M_md_motion,      only : elses_md_motion, elses_ini_set_velocity !(routine)
!   use M_md_verlet,      only : md_motion_verlet_velocity               !(routine)
    use M_md_velocity_dst,  only : md_motion_verlet_velocity_dst         !(routine)
    use M_md_save_struct, only : elses_md_save_struct                    !(routine)
    use M_md_output,      only : elses_md_output_main                    !(routine)
    use M_wall_clock_time,only : get_elapse_wall_clock_time              !(routine)
    use M_qm_engine,      only : elses_qm_engine                         !(routine)
    use M_ini_load_geno,   only : ini_load_geno                          !(routine)
    use M_qm_engine_dst,   only : qm_engine_dst                          !(routine)
    use M_lib_dst_info,     only : mpi_is_active, myrank, nprocs !(unchanged)
    use M_lib_dst_info,     only : log_unit         !(unchanged)
    use M_lib_dst_info,     only : root_node        !(unchanged)
!   use M_io_dst_write_log, only : log_file_is_set  !(unchanged)
!   use M_io_dst_write_log, only : set_dst_write_log_file    !(routine)
    use M_md_dst,           only : show_mpi_info !(routine)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_barrier_time !(routine)
    use M_lib_stop_signal,  only : detect_stop_signal !(routine)
    use M_lib_rep_mpi_time, only : get_lap_of_allreduce_time     !(routine) 
    use M_lib_rep_mpi_time, only : get_lap_of_barrier_time       !(routine) 
    use M_lib_rep_mpi_time, only : record_allreduce_barrier_time !(routine)
    use M_wall_clock_time,  only : measure_clock_time_period     !(routine)
    use M_wall_clock_time,  only : show_reset_period_of_clock    !(routine)
    use M_md_motion,        only : optimization_check            !(routine)
    use M_qm_output_matrices, only : output_levels_matrices      !(routine)
    use M_qm_geno_output,     only : output_for_eigen_solver     !(routine)
    use M_output_wfn_charge,  only : output_for_wfn_charge       !(routine)
!
    use M_lib_timer_in_thread,    only : init_timer_in_thread       !(routine)
    use M_lib_timer_in_thread,    only : matvec_timer_in_thread     !(CHANGED)
!
    implicit none
!   real*8 tb0, tb, ptime, entime, entime_bak
    integer istop, itemd2
    real(8) :: elapse_time, elapse_time_previous, time_period
    integer :: stop_signal
!   logical :: first_iteration
!   logical :: final_iteration
    logical :: dst_calculation
    character(len=32) :: scheme_mode
    real(8) :: time_wrk, time_wrk_previous, time_check ! work variable for measuring the time
    integer :: imode
    logical, parameter :: mpi_time_check = .true.
    integer itemdmx_result
    real(8) :: elapse_time_mpi_allreduce, elapse_time_mpi_barrier
    logical, parameter :: record_total_mpi_time_is_active = .false.
    logical :: converged_in_optimization
    integer :: ierr
    real(8) :: matvec_timer_wrk
    integer :: ipe
!
!   log_file_is_set = .false.    ! dummy setting
!   log_unit        = 0          ! dummy setting
!   call set_dst_write_log_file
!    ----> set log_file_is_set, log_unit
!
    if (log_unit > 0) call show_mpi_info
!    ----> write the mpi info into the specified logfile
!
    if (root_node) then
      write(*,*) '-----------------------'
      write(*,*) '|ELSES standard output|'
      write(*,*) '-----------------------'
      call show_reset_period_of_clock(time_period) 
      write(*,*) 'INFO:The reset period of the system clock (sec) =', time_period
      write(*,*) 'INFO:The largest INTEGER-type number: huge(0)   =', huge(0)
    endif   
!
    call get_elapse_wall_clock_time(elapse_time)
    elapse_time_previous=elapse_time
!
     call get_elapse_wall_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
    call get_lap_of_allreduce_time(elapse_time_mpi_allreduce)
    call get_lap_of_barrier_time(elapse_time_mpi_barrier)
!
!   entime_bak=0.0d0
!   call tclock(ptime)
!   tb0=ptime
!   entime_bak=0.0d0
!
!   call elses_md_init
!
    if (mpi_time_check) then
      call mpi_wrapper_barrier_time(time_check)
      if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(bef ini_load_geno)= ',time_check
    endif
!
    call get_elapse_wall_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
    call ini_load_geno
!         ---> Initial setting 
!
    if (mpi_time_check) then
      call mpi_wrapper_barrier_time(time_check)
      if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(aft ini_load_geno)= ',time_check
    endif
!
    call get_elapse_wall_clock_time(time_wrk)
    call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
    if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:ini_load_geno = ',time_period
    time_wrk_previous=time_wrk
!
    call elses_ini_set_velocity
!         ---> Set the initial velocity, if necessary
!
    if (mpi_time_check) then
      call mpi_wrapper_barrier_time(time_check)
      if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(aft ini_set_velocity)= ',time_check
    endif
!
    call get_elapse_wall_clock_time(time_wrk)
    call measure_clock_time_period(time_wrk, time_wrk_previous, time_period)
    if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME: ini_set_velocity= ',time_period
    time_wrk_previous=time_wrk
!
    istop=0
    itemd2=0
!   call elses_mdsave_leg(istop,entime,itemd2)
!         ---> load the data, if needed.
!
    if (mpi_time_check) then
      call mpi_wrapper_barrier_time(time_check)
!     write(*,'(a,f20.10)')'TIME:mpi_check(bef)= ',time_check
      if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(bef)= ',time_check
    endif
!
    call get_elapse_wall_clock_time(elapse_time)
    call get_lap_of_allreduce_time(elapse_time_mpi_allreduce)
    call get_lap_of_barrier_time(elapse_time_mpi_barrier)
!
    matvec_timer_wrk=0.0d0
    if (root_node) then
      print '(A,I15,F24.10)', 'elaps-time(befor  MDloop)=', & 
&            itemd2,elapse_time-elapse_time_previous
    endif  
    if (log_unit > 0) write(log_unit,'(A,I15,4F15.6)')  'elaps-time(befor  MDloop)=', &
&            itemd2,elapse_time-elapse_time_previous, & 
&                 elapse_time_mpi_allreduce, elapse_time_mpi_barrier, matvec_timer_wrk

    elapse_time_previous=elapse_time
!
!   write(*,*)'config%option%test_mode=', trim(adjustl(config%option%test_mode))
!
    if (trim(adjustl(config%option%test_mode)) == 'initial') then
      if (log_unit > 0) write(log_unit, '(a)') 'INFO:TEST MODE ONLY FOR INITIAL PROCEDURE !'
      if (root_node)    write(*,'(a)')         'INFO:TEST MODE ONLY FOR INITIAL PROCEDURE !'
      itemdmx=0
    else
      if (log_unit > 0) write(log_unit, *) '--------MD loop starts:itemdmx=',itemdmx
      if (root_node)    write(*,*)         '--------MD loop starts:itemdmx=',itemdmx
    endif  
!
    itemdmx_result=itemdmx
    do itemd2=1,itemdmx
!
      final_iteration = .false.
      first_iteration = .false.
      if (itemd2 == 1) first_iteration = .true.
      if (itemd2 == itemdmx) final_iteration = .true.
!
!     call tclock(ptime)
!     tb=ptime
!         ---> Measure the time
!
!     call elses_md_save_xml
!
!     if (itemd2 .ge. 2) then
!       istop=0
!       call tclock2(entime,entime_bak,tb0)
!       call elses_mdsave_leg(istop,entime,itemd2)
!     endif  
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
!      write(*,'(a,f20.10)')'TIME:mpi_check(ini)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(ini)= ',time_check
     endif
!
     scheme_mode=trim(config%calc%solver%scheme)
     if (scheme_mode == 'ekrgl') scheme_mode='gKrylov_A'
     config%calc%solver%scheme=scheme_mode
     if (root_node) then
       write(*,*)'scheme_mode=',scheme_mode
       write(*,*)'c_system   =',c_system 
     endif  
!
     dst_calculation = config%calc%distributed%set
!    if (c_system /= 'geno') then
!      write(*,*)'Warning:NO DST calculation, because system=',c_system
!      dst_calculation = .false.
!    endif   
!
     if (root_node) then
       write(*,*)'INFO:dst_calculation=',dst_calculation
       write(*,*)'INFO:mpi_is_active  =',mpi_is_active
     endif  
!
     if (mpi_is_active) then
       if ((nprocs > 1) .and. (.not. dst_calculation)) then
         write(*,*) 'WARNING:Non-DST calcuation with nprocs > 1 : nprocs=', nprocs
       endif   
     endif   
!
     call get_elapse_wall_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
     if ( dst_calculation ) then
       if (root_node) write(*,*)'INFO:Distributed calculation mode'
       call init_timer_in_thread
       call qm_engine_dst
     else
       if (root_node) write(*,*)'INFO:NON-distributed calculation mode'
       call elses_qm_engine  
     endif  
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
!      write(*,'(a,f20.10)')'TIME:mpi_check(qme)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(qme)= ',time_check
     endif
!
      call get_elapse_wall_clock_time(time_wrk)
!     write(*,*)'TIME:qm_engine = ',time_wrk-time_wrk_previous
      if (log_unit > 0) write(log_unit,*) 'TIME:qm_engine = ',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk 
!
      call detect_stop_signal(stop_signal, time_wrk)
!         ---> Detect the stop signal, if any
!
      if (stop_signal /= 0)  final_iteration = .true.
!
      call optimization_check(converged_in_optimization)
!         ---> Check the convergence (only in the optimization scheme)
      if (converged_in_optimization) final_iteration = .true.
!
      call output_levels_matrices
!         ---> Plot the levels and matrices
!
      call output_for_eigen_solver
!         ---> Plot the wavefunction and so on (only for the eigen solver)
!
      call output_for_wfn_charge
!         ---> Plot the wavefunction charge
!
      imode = 0
      if (first_iteration) imode=1
!     call md_motion_verlet_velocity(imode)
      call md_motion_verlet_velocity_dst(imode)
!         ---> Update the velocity (only in the dynamics mode)
!
!     stop 'STOP MANUALLY'
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
!      write(*,'(a,f20.10)')'TIME:mpi_check(vel)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(vel)= ',time_check
     endif
!
      call elses_md_save_struct(final_iteration)
!
!     if (mpi_is_active) then
!       if (myrank == 0) then
!         write(*,*)'INFO-MPI:elses_md_save_struct is called only at the root node (myrank=0)'
!         if (log_unit > 0) write(log_unit,*) 'INFO-MPI:elses_md_save_struct is called only at the root node (myrank=0)'
!         call elses_md_save_struct(final_iteration)
!        endif
!      else
!        call elses_md_save_struct(final_iteration)
!      endif
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
!      write(*,'(a,f20.10)')'TIME:mpi_check(svs)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(svs)= ',time_check
     endif
!
!     if (stop_signal /= 0) then
!       write(*,*)'.....stop by user signal (in 00_stop_signal.txt)'
!       stop
!     endif   
!
      call elses_md_output_main
!         ---> Output (energy and so on)
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
!      write(*,'(a,f20.10)')'TIME:mpi_check(out)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out)= ',time_check
     endif
!
      if (.not. final_iteration) call elses_md_motion
!         ---> Update the atom positions
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
!      write(*,'(a,f20.10)')'TIME:mpi_check(mov)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov)= ',time_check
     endif
!
!
     call get_lap_of_allreduce_time(elapse_time_mpi_allreduce)
     call get_lap_of_barrier_time(elapse_time_mpi_barrier)
!
     call get_elapse_wall_clock_time(elapse_time)
     if (root_node) then
       print '(A,I15,F24.10)', 'elaps-time(lap    MDloop)=', &
&            itemd2,elapse_time-elapse_time_previous
     endif 
!
     matvec_timer_wrk=0.0d0
     if (allocated(matvec_timer_in_thread)) then 
       matvec_timer_wrk=sum(matvec_timer_in_thread(:))/size(matvec_timer_in_thread,1)
       do ipe=1, size(matvec_timer_in_thread,1)
        if (log_unit > 0) then 
          write(log_unit,*) 'INFO:matvec time per thread =', itemd2, ipe, matvec_timer_in_thread(ipe)
        endif
       enddo
     else
       matvec_timer_wrk=0.0d0
     endif
!
     if (log_unit > 0) write(log_unit,'(A,I15, 4F15.6)') 'elaps-time(lap    MDloop)=', &
&            itemd2,elapse_time-elapse_time_previous, & 
&                 elapse_time_mpi_allreduce, elapse_time_mpi_barrier, matvec_timer_wrk
      elapse_time_previous=elapse_time
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(fin)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(fin)= ',time_check
     endif
!
     if (final_iteration) then 
       itemdmx_result=itemd2 
       exit
     endif 
!
    enddo
!
    call get_elapse_wall_clock_time(elapse_time)
!
    if (root_node) then
      print '(A,I30)',  '.....Total MD loops ends; # MD step =', &
&                                 itemdmx_result
      print '(A,F30.5)','.....Total Simulation time (sec ) =',    &
&                                 elapse_time
      print '(A,F30.5)','.....Total Simulation time (hour) =',    &
&                                 elapse_time/3600.0d0
      print '(A,F30.5)','.....Total Simulation time (day ) =',    &
&                                 elapse_time/3600.0d0/24.0d0
      if (itemdmx_result /= itemdmx) then
        print '(A)', 'INFO: The simulation ends, before the maximum iteration,'
        print '(A)', 'INFO:   ... for the convergence or the stop signal from user'
      endif   
    endif  
!
   if (log_unit > 0) then 
     write(log_unit,'(A,I30)')   '.....Total MD loops ended; # MD step =', itemdmx_result
     write(log_unit,'(A,F30.5)') '.....Total Simulation time (sec )    =', elapse_time
     write(log_unit,'(A,F30.5)') '.....Total Simulation time (hour)    =', elapse_time/3600.0d0
     write(log_unit,'(A,F30.5)') '.....Total Simulation time (day )    =', elapse_time/3600.0d0/24.0d0
     if (itemdmx_result /= itemdmx) then
       write(log_unit,'(A)') 'INFO: The simulation ends, before the maximum iteration,'
       write(log_unit,'(A)') 'INFO:   ... for the convergence or the stop signal from user'
     endif  
   endif
!   
   if (record_total_mpi_time_is_active)  call record_allreduce_barrier_time('output-total-consumed-mpi-time.txt')
!
  end subroutine elses_md_main_01
!
end module M_md_main
