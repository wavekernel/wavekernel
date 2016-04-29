!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_output
!
   use elses_mod_ctrl,      only : i_verbose     !(unchanged)
   use M_00_v_info,         only : version_info  !(unchanged)
   use M_lib_dst_info,      only : log_unit      !(unchanged)
   use M_lib_dst_info,      only : output_unit   !(CHANGED in setting_for_main_output)
!
!  NOTE: The following module variables are set in elses_md_output_set
!
   implicit none
   integer :: step_count
!
   integer, allocatable :: unit_num_for_output(:)
!        ---> Unit number for the output file 
!
   character(len=*), parameter :: filename_for_output="Output.txt"
!        ---> File name of the output file

   private
   public :: elses_md_output_main
   public :: setting_for_main_output
!
 contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for output
!
   subroutine elses_md_output_main
!
     use M_config,             only : config !(unchanged)
     use M_output_atom_charge, only : output_atom_charge     !(routine)
     use M_output_atom_energy, only : output_atom_energy     !(routine)
     use M_qm_domain,          only : dhij, dsij, dbij, dpij !(unchanged)
     use elses_mod_file_io,    only : vacant_unit  !(function)
     use M_config,             only : config !(unchanged)
!    use M_qm_output_matrices, only : qm_output_matrices     !(routine)
!    use M_qm_output_levels,   only : qm_output_levels       !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_barrier_time  !(routine)
     use M_md_virial_pressure, only : plot_virial_pressure   !(routine)
!
     implicit none
     integer :: i_global_mat
     logical, parameter :: plot_hamiltonian = .false.
     logical, parameter :: plot_levels      = .false.
!    character(len=8) :: mode_for_plot_matrices
     logical          :: mpi_time_check
     logical          :: small_output
     logical          :: flag_for_init
     integer :: unit_num, ierr
     real(8) :: time_wrk, time_wrk_previous, time_check
!
     if (config%calc%distributed%set .eqv. .true.) then
       mpi_time_check = .true.
       small_output   = .true.
     else
       mpi_time_check = .false.
       small_output   = .false.
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Initial setting for main output file
!
     if (config%calc%distributed%root_node) then
       if (output_unit <= 0) then
         write(*,*)'ERROR(elses_md_output_main):output_unit=',output_unit
         stop
       endif   
     endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    if (.not. allocated(unit_num_for_output)) then
!      allocate (unit_num_for_output(1),stat=ierr)
!      if ( ierr .ne. 0) then
!        write(*,*)'Alloc. error(elses_md_output_main):ierr=',ierr
!        stop
!      endif  
!      unit_num_for_output(1)=vacant_unit()
!      flag_for_init = .true.
!    else
!      flag_for_init = .false.
!    endif  
!
!    unit_num=output_unit
!
     step_count=config%system%structure%mdstep
!     
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out1)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out1)= ',time_check
     endif
!
!    if (log_file_is_set) then
!      unit_num=log_unit
!    else
!      unit_num=vacant_unit()       
!    endif
!
     if (i_verbose >= 1) then
       write(*,'("@@ elses_md_output_main: step_count = ",I10)') step_count
       if (log_unit > 0) write(*,'("@@ elses_md_output_main: step_count = ",I10)') step_count
     endif  
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out2)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out2)= ',time_check
     endif
!
!
     if (output_unit > 0) then 
       write(output_unit,'("------------------------------------------------------------------")')
       write(output_unit,'("Output; step_count=",I10)') step_count
     endif
!  
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out3)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out3)= ',time_check
     endif
!
     call output_energy
!      --> output energy and so on
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out4)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out4)= ',time_check
     endif
!
     call md_output_compat
!      --> call the old routine for compatibility
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out5)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out5)= ',time_check
     endif
!
     call output_csc_convergence
!
     if (small_output) return
!
     call output_force_ave_max
!      --> output force amplitude
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out6)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out6)= ',time_check
     endif
!
     call output_atom_charge
!        --> output for atom charge (optional)
!
     call plot_virial_pressure
!        --> plot the virial pressure (optional)
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(out7)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(out7)= ',time_check
     endif
!
!    if (plot_hamiltonian) then 
!      mode_for_plot_matrices='H'
!      call qm_output_matrices(mode_for_plot_matrices)
!      call qm_output_matrices
!    endif
!        --> output for Hamiltonian and overlap matrices, if defined
!
!    if (plot_levels) then 
!      call qm_output_levels
!    endif
!        --> output for eigen levels
!
     i_global_mat=1
     if (.not. allocated(dhij)) i_global_mat=0
     if (.not. allocated(dsij)) i_global_mat=0
     if (.not. allocated(dbij)) i_global_mat=0
     if (.not. allocated(dpij)) i_global_mat=0
!
     if (i_global_mat == 1) then
       call output_atom_energy
!        --> output for atom energy (optional)
!    else
!      write(*,*)'INFO-WARN:No global matrices are defined'
!      write(*,*)'INFO-WARN:SKIPPED  output_atom_energy'
     endif  
!
     if (log_unit > 0) then
       write(log_unit,'("------------------------------------------------------------------")')
     endif  
!
     if (i_verbose >= 1) then
       call plot_detailed_kin_energy
     endif
!
!    close(unit_num)
!
   end subroutine elses_md_output_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting for main output file
!   --> Set step_count, flag_for_unit, unit_num, local_switches
!
   subroutine setting_for_main_output
!
     use M_config,            only : config  !(unchanged)
     use M_lib_dst_info,      only: myrank, nprocs !(unchanged)
     use elses_mod_file_io,   only : vacant_unit   !(function)
     implicit none
     integer :: ierr
     character(len=32) :: output_filename
     integer :: omp_get_thread_num
     integer :: omp_get_num_threads
     integer :: ipe, npe
     character(len=8)  :: chara_date
     character(len=10) :: chara_time
!
     if (.not. config%calc%distributed%root_node) return 
!   
     output_unit=vacant_unit() 
     output_filename=trim(config%output%main%filename)
     if (trim(output_filename) == "") then
       write(*,*)'ERROR(elses_md_output_main):empty output_filename'
       stop
     endif   
     open (output_unit, file=output_filename, status='unknown')
     write(output_unit,'(a,a)') '@@ Main output : ', trim(version_info)
     write(output_unit,'(a,a)') 'INFO: config_name = ', trim(config%name)
!
!$omp parallel default(shared) &
!$omp& private(ipe,npe)
    ipe=-1
    npe=-1
!$  ipe=omp_get_thread_num()+1
!$  npe=omp_get_num_threads()
!
    if (npe == -1) then 
      write(output_unit,'(a,2i10,a)')'INFO-MPI-OMP: P_MPI, P_OMP=', nprocs, npe+2, ' (NO-OMP or OMP-STUB)'
    else  
      if (ipe == 1) write(output_unit,'(a,2i10)')  'INFO-MPI-OMP: P_MPI, P_OMP=', nprocs, npe
    endif  
!
!$omp end parallel
!
    call date_and_time(chara_date, chara_time)
    write(output_unit,'(13a)')'Date: ', chara_date(1:4), ' ', chara_date(5:6), ' ', chara_date(7:8), '; ', &
&                             'Time: ', chara_time(1:2), ' ', chara_time(3:4), ' ', chara_time(5:6) 
!
   end subroutine setting_for_main_output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine for plotting energy
!
   subroutine output_energy
!
     use elses_mod_md_dat,     only : e_kin
     use elses_mod_ene,        only : etb, ecc,ecsc
     use elses_mod_sim_cell,   only : noa
     use elses_mod_phys_const, only : ev4au
     use M_md_motion,          only : sd_energy_diff !(unchanged)
     use M_qm_domain,          only : chemical_potential !(unchanged) 
     implicit none
     integer :: unit_num
     real(8) :: con_unit
!      ---> constant for unit conversion
!
     con_unit=ev4au/dble(noa)
!
     if (output_unit > 0) then
       unit_num=output_unit
     else   
       return 
     endif
!   
     write(unit_num,*)'Output energy'
!
     write(unit_num,'(" Band         Energy : EBD  [au]: ",F30.15)') etb
     write(unit_num,'(" ECSC         Energy : ECSC [au]: ",F30.15)') ecsc
     write(unit_num,'(" Core-core    Energy : ECC  [au]: ",F30.15)') ecc
     write(unit_num,'(" Kinetic      Energy : EKE  [au]: ",F30.15)') e_kin
     write(unit_num,'(" EBD+ECSC+ECC Energy : ETOT [au]: ",F30.15)') etb+ecsc+ecc
     write(unit_num, &
& '(" Energy summary (explan.):  step_count        EBD+ECSC+ECC      EBD+ECSC+ECC+EKE  EBD &
&              ECSC              ECC                EKE")')

     if (abs(etb)+abs(ecsc)+abs(ecc)+abs(e_kin) > 1.0d7) then
       write(unit_num,'(" Energy summary (au     ): ",I12,6F23.8)') step_count, & 
&         etb+ecsc+ecc, etb+ecsc+ecc+e_kin, etb, ecsc, ecc, e_kin
     else
       write(unit_num,'(" Energy summary (au     ): ",I12,6F18.8)') step_count, & 
&         etb+ecsc+ecc, etb+ecsc+ecc+e_kin, etb, ecsc, ecc, e_kin
     endif   
     write(unit_num,'(" Energy summary (eV/atom): ",I12,6F18.8)') step_count, & 
&         (etb+ecsc+ecc)*con_unit, (etb+ecsc+ecc+e_kin)*con_unit,  &
&         etb*con_unit, ecsc*con_unit, ecc*con_unit, e_kin*con_unit
     if (log_unit > 0) then
       write(log_unit,'(" Energy summary (eV/atom): ",I12,6F18.8)') step_count, & 
&         (etb+ecsc+ecc)*con_unit, (etb+ecsc+ecc+e_kin)*con_unit,  &
&         etb*con_unit, ecsc*con_unit, ecc*con_unit, e_kin*con_unit
     endif
!
     write(unit_num,'(a,i10,2f30.20)') ' Chemical potential [au, eV]=', & 
&                                        step_count, chemical_potential, chemical_potential*ev4au
     if (allocated(sd_energy_diff)) then
        write(unit_num,'(a,i10,f30.20)') ' SD energy difference per atom [eV]=', step_count, sd_energy_diff(1)*con_unit
     endif   
!
   end subroutine output_energy
!
   subroutine output_force_ave_max
     use elses_mod_phys_const, only : ev4au, angst
     use elses_mod_sim_cell,   only : noa
     use elses_mod_foi,    only : foi
     implicit none
     real(8) :: ddave, ddmax, ddd, ddd1, ddd2, ddd3, ddave_au, ddmax_au
     integer :: js, js_max
     integer :: unit_num
!
     if (output_unit > 0) then
       unit_num=output_unit
     else   
       return 
     endif
!
     write(unit_num,"(' Output Force Amplitude ')")
     ddave=0.0d0
     ddmax=0.0d0
     js_max=0
     do js=1,noa
      ddd=0.0d0
      ddd1=foi(js,1)
      ddd2=foi(js,2)
      ddd3=foi(js,3)
      ddd=ddd1*ddd1+ddd2*ddd2+ddd3*ddd3
      ddave=ddave+ddd
      if (ddd > ddmax) then
        ddmax=ddd
        js_max=js
      endif   
     enddo
     ddave=dabs(ddave/dble(noa))
     ddmax=dabs(ddmax)
     ddave_au=ddave
     ddmax_au=ddmax
     ddave=ddave*ev4au/angst
     ddmax=ddmax*ev4au/angst
     write(unit_num,'(" Force Amp. Average [au] [eV/A]= ",2F30.15)') ddave_au, ddave
     write(unit_num,'(" Force Amp. Max     [au] [eV/A]= ",2F30.15)') ddmax_au, ddmax
     write(unit_num,'(" The atom that gives the max. force amp. = ",I10)') js_max
     write(unit_num,"(' Force_summary(ave[eV/A],max[eV/A],atom for max)=',I10,2F30.20,I10)") step_count,ddave,ddmax,js_max
   end subroutine output_force_ave_max
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ OLD output routine for compatibility
!
   subroutine md_output_compat
!
     use M_config,             only : config  !(unchanged)(only config%system%temperature)
     use elses_mod_phys_const, only : ev4au
     use elses_mod_sim_cell,   only : noa
!    use elses_mod_thermo, only : tempk0,amq,thbold,vhb
     use elses_mod_thermo, only : amq,thbold,vhb
     use elses_mod_md_dat, only : itemd, itemdorg,dtmd, e_kin
     use elses_mod_ene,     only : etb, ecc
!
     implicit none
     integer itemd4, inose
     real(8) dtmd2, eki, etb2, eadd
     real(8) ehb, ehb1, ehb2, pkin2
!
     real(8) eptadd1, eptadd2, fb, entr, tentr
!        ----> dummy only for combatibility to the old code
     real(8) tempk0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     tempk0=config%system%temperature
     dtmd2=dtmd*dtmd
     inose=1
     eki=e_kin
     etb2=etb
!
     eadd=0.0d0
     eptadd1=0.0d0
     eptadd2=0.0d0
     fb=0.0d0
     entr=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if (i_verbose >=0) then
       write(*,*)'@@ md_output_compat'
     endif 
!
     itemd4=itemd
     if (itemdorg .ne. 0) itemd4=itemd+itemdorg-2
     entr=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     EHB=0.0D0
     EHB1=0.0D0
     EHB2=0.0D0
     PKIN2=0.0D0
!
     if (inose /= 0) then
!
       EHB1=0.5d0*VHB*VHB/DTMD2/AMQ
       EHB2=3.0D0*DBLE(NOA)*TEMPK0*THBOLD
       EHB=EHB1+EHB2
       if (i_verbose >=0) then
         WRITE(*,*)'THB,VHB at t =',THBOLD,VHB
         WRITE(*,'(a,3f25.10)') 'EHB,1,2=',EHB,EHB1,EHB2
       endif
       PKIN2=3.0D0/2.0D0*DBLE(NOA)*TEMPK0
       if (i_verbose >=0) then
         write(*,'(a,i8,3f25.10)') 'PKIN,d_PKIN,VHB=',ITEMD,PKIN2,EKI-PKIN2,VHB
       endif  
!
!         PKIN2: given kinetic-energy (temperature)
!         EKI  : kinetic-energy of particles
!         VHB  : (damping coeficient) x dt
!         When ( EKI > PKIN2 ), VHB must be positive
!         When ( EKI < PKIN2 ), VHB must be negative
!        
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      TENTR=FB/EV4AU*ENTR
!
      if (i_verbose >=0) then
        WRITE(6,*)'    ETB2=',ETB2
        WRITE(6,*)'    ECC =',ECC
        WRITE(6,*)'EPTADD1 =',EPTADD1
        WRITE(6,*)'EPTADD2 =',EPTADD2
        WRITE(6,*)'    EKI =',EKI
        WRITE(6,*)'    EHB =',EHB
        WRITE(6,*)'    EADD=',EADD
        WRITE(6,*)' FB [eV]=',FB
        WRITE(6,*)' FB [au]=',FB/EV4AU
        WRITE(6,*)'   ENTR =',ENTR
        WRITE(*,'(a,2i8,7f18.6)')'IONMOV3',ITEMD,step_count,ETB2+ECC+EPTADD1+EPTADD2, &
&          ETB2,ECC,EKI,EPTADD1+EPTADD2,EHB,EADD
      endif  
!
!
   end subroutine md_output_compat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine for plotting the CSC convergence info.
!
   subroutine output_csc_convergence
     use M_config,             only : config !(unchanged)
     use M_qm_domain,   only : csc_ite_converged, csc_dq_converged !(unchanged)
     implicit none
     integer unit_num
!
     if (output_unit > 0) then
       unit_num=output_unit
     else   
       return 
     endif
!
     if (config%calc%genoOption%CSC_max_loop_count == 0) return
!
     write(output_unit,'(a, 2i10,f20.10)') ' CSC convergence info (step_count, CSC ite., dq) :',  &
&                            config%system%structure%mdstep, csc_ite_converged, csc_dq_converged

!
   end subroutine output_csc_convergence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine plot_detailed_kin_energy
     use M_config,               only : config !(unchanged)
     use M_md_velocity_routines, only : calc_kinetic_energy !(routine)
     implicit none
     integer :: lu
     real(8) :: kinetic_energy
     real(8) :: k_e_component(3)
!
     lu=log_unit
!
     if (i_verbose >=1) then
       if (lu >0) write(lu,*)'@@ plot_detailed_kin_energy'
       call calc_kinetic_energy(kinetic_energy, k_e_component)
       if (lu >0) write(lu,'(a,i10,4f20.10)')'E_kin_t,x,y,z=', &
&          config%system%structure%mdstep, kinetic_energy,k_e_component(1:3)
     endif
!
   end subroutine plot_detailed_kin_energy
!
end module M_md_output

