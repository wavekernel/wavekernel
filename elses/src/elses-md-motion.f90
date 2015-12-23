!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_motion
!
   use elses_mod_ctrl,   only: i_verbose
   use M_io_dst_write_log,  only : log_unit
!
   implicit none
   real(8), allocatable :: sd_energy_diff(:)  ! SD energy difference, given in optimization_check
!
   private
   public :: sd_energy_diff
   public :: elses_md_motion
   public :: elses_md_motion_ini_set
   public :: elses_ini_set_velocity
   public :: md_update_tx
   public :: get_atom_dist
   public :: optimization_check
!
 contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial parameter setting for atom motion
!
   subroutine elses_md_motion_ini_set
     use M_config, only: config
     use elses_mod_md_dat, only: dtmd, itemdmx
     implicit none
     real(8) :: dtmd_tmp, tot_sim_time
     real(8) :: sd_ratio_local
     character(len=20) :: motion_type
!
     motion_type=trim(config%calc%mode) ! now this value is set at elses-xml-config.f90
     dtmd = -1.0d0 ! dummy value
!
!
     if (i_verbose >= 1) then
       write(*,'("@@ elses_md_motion_ini_set: mode = ",A)') config%calc%mode
     endif  
     !
     select case(trim(motion_type))
     case ("dynamics")
       tot_sim_time    = config%calc%dynamics%total ! a.u.
       dtmd_tmp        = config%calc%dynamics%delta ! a.u.
       if (dtmd_tmp < 1.0d-10 ) then
        write(*,*)'ERROR(elses_md_motion_ini_set): dtmd=',dtmd_tmp
        write(*,*)' This value should be positive'
        stop
       endif
       itemdmx=nint(tot_sim_time/dtmd_tmp)+1
       dtmd=dtmd_tmp
       if (itemdmx <= 0 ) then
        write(*,*)'ERROR(elses_md_motion_ini_set): itemdmx=',itemdmx
        write(*,*)' Total MD step:This value should be larger than zero'
        stop
       endif
       if (i_verbose >= 1 ) then
         write(*,*)'Total MD step =',itemdmx
         write(*,*)'Time slice dt =',dtmd
       endif
     case ("optimization")
       itemdmx=config%calc%optimization%max_num_iter+1
       sd_ratio_local=config%calc%optimization%sd_ratio
       if (itemdmx <= 0 ) then
         write(*,*)'ERROR(elses_md_motion_ini_set): itemdmx =',itemdmx
         write(*,*)' Maximum Iteration number : This value should be larger than zero'
         stop
       endif
       if (sd_ratio_local <= 1.0d-10 ) then
         write(*,*)'Warning(elses_md_motion_ini_set): sd_ratio =',sd_ratio_local
!        write(*,*)' Invalid value: this value should be positive'
!        stop
       endif
       if (i_verbose >= 1 ) then
         write(*,*)'Maximum iteration number in optimization mode=',itemdmx
       endif
     case ("conversion")
       itemdmx=1  ! dummy value (so as to be compatible to other routines; T.Hoshi)
       dtmd=1.0d0 ! dummy value (so as to be compatible to other routines; T.Hoshi)
       config%calc%dynamics%delta = dtmd
     case ("matrix_generation")
       write(*,*)'matrix generation'
       itemdmx=1  ! dummy value (so as to be compatible to other routines; T.Hoshi)
       dtmd=1.0d0 ! dummy value (so as to be compatible to other routines; T.Hoshi)
       config%calc%dynamics%delta = dtmd
     case ("snapshot")
       itemdmx=config%calc%snapshot%number
     case ("cell_change_only")
       itemdmx=config%calc%cell_change%max_num_iter
     case default
       write(*,'("     config%calc%mode = ", A)') motion_type
       stop "elses_md_motion: unknown motion_type"
     end select
     !
   end subroutine elses_md_motion_ini_set
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial setting of velocity
!
   subroutine elses_ini_set_velocity
     use M_config, only: config
!    use elses_mod_vel,      only : velx, vely, velz
     use M_md_velocity_routines, only : calc_initial_velocity
     character(len=32) :: motion_type
!
     motion_type=trim(config%calc%mode)
!
     if (motion_type == "dynamics") then 
       if (i_verbose >= 1) write(*,*)'@@ elses_ini_set_velocity'
       call calc_initial_velocity
     else  
       if (i_verbose >= 1) write(*,*)'@@ elses_ini_set_velocity....is skipped'
     endif  
!   
!
!    select case(motion_type)
!    case ("dynamics")
!      call calc_initial_velocity
!      call elses_vel_ini_chk
!    case ("optimization", "snapshot")
!      if (i_verbose >= 1) then
!        write(*,*)' INFO : set velocity to be zero'
!      endif
!      velx(:)=0.0d0
!      vely(:)=0.0d0
!      velz(:)=0.0d0
!    case default
!      stop "elses_md_motion: unknown motion_type"
!    end select
!       
   end subroutine elses_ini_set_velocity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for atom motion
!
   subroutine elses_md_motion
     use M_config, only: config
     use M_md_cell_change, only : elses_md_cell_change      !(routine)
     use M_md_constraint,  only : elses_md_constraint       !(routine)
     use M_md_verlet,      only : md_motion_verlet_position !(routine)
     use M_qm_population,  only : update_population_guess   !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_barrier_time !(routine)
     use M_io_update_snapshot, only : update_snapshot       !(routine)
     use M_gid_constraint,   only : add_constraint_w_groups !(routine)
!
     implicit none
     character(len=32) :: motion_type
     logical, parameter :: mpi_time_check = .true.
     real(8) :: time_check
     logical, parameter :: cell_change_for_next_step = .true.
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(mov1)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov1)= ',time_check
     endif
!
     motion_type=trim(config%calc%mode) !now this value is set at elses-xml-config.f90
     write(*,'("@@ ELSES_MD_MOTION, config%calc%mode = ",A)')config%calc%mode
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(mov2)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov2)= ',time_check
     endif
!
     select case(motion_type)
     case ("dynamics")
        call md_motion_verlet_position
        !       --> Atom motion by Newton eq.
     case ("optimization")
        if (i_verbose >= 1) then
          write(*,'("     config%calc%optimization%scheme = ",A)')      config%calc%optimization%scheme
          write(*,'("                            sd_ratio =",ES21.14)') config%calc%optimization%sd_ratio
          write(*,'("                        max_num_iter =",I8)')      config%calc%optimization%max_num_iter
        endif  
        call steepest_descent(config%calc%optimization%sd_ratio)
          !       --> Atom motion by steepest descent
     case ("snapshot")
        call update_snapshot
        call md_update_tx
          !       --> Update the snapshot data
     case ("cell_change_only")
        write(*,*)'cell_change_only mode'
     case default
        write(*,'("     config%calc%mode = ", A)') config%calc%mode
        stop "elses_md_motion: unknown motion_type"
     end select
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(mov3)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov3)= ',time_check
     endif
!
     if (config%system%structure%use_vatom) then
       call update_population_guess
!       ---> update config%system%structure%vatom(:)%population_guess for the next CSC loop
     endif  
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(mov4)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov4)= ',time_check
     endif
!
     call elses_md_cell_change(cell_change_for_next_step)
!      ---> Change the simulation (peridic) cell sizes (experimental)
!            only if the file 'Cell.txt' exist.
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(mov5)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov5)= ',time_check
     endif
!
     call elses_md_constraint
!      ---> Constraint motion
!
     call add_constraint_w_groups
!      ---> Add the constraint on atom positions with group id (Optional)
!
     if (mpi_time_check) then
       call mpi_wrapper_barrier_time(time_check)
       write(*,'(a,f20.10)')'TIME:mpi_check(mov6)= ',time_check
       if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:mpi_check(mov6)= ',time_check
     endif
!
   end subroutine elses_md_motion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the convergence in the optimization mode 
!
  subroutine optimization_check(converged)
!
    use M_config, only: config
    use elses_mod_sim_cell,   only : noa
    use elses_mod_iflag,      only : iflag
    use elses_mod_md_dat, only : itemd
    use elses_mod_foi,    only : foi
    use M_lib_phys_const, only : ev4au
!
    implicit none
    logical, intent(out) :: converged
    integer  js, js_max, ierr
    real(8)  ddave, ddmax, ddd, ddd1, ddd2, ddd3
    real(8)  da1, da2, da3, ddrr1, ddrr2, ddrr3
    real(8)  ddsum, ddsumx, ddsumy, ddsumz 
    real(8)  sd_ratio
!
!
    converged = .false.  ! dummy value
!
    if (i_verbose >= 1) then
      if (log_unit>0) write(log_unit,*)'@@ optimization_check'
    endif  
!
    if (trim(config%calc%mode) /= 'optimization') then
      if (log_unit>0) write(log_unit,*)'    ..is skipped'
      return
    endif   
!
    if (.not. allocated(iflag)) then
      write(*,*) 'ERROR(optimization_check):not allocated:iflag' 
      stop
    endif   
!
    sd_ratio = config%calc%optimization%sd_ratio
!
    if (sd_ratio < -1.0d-10) then
      write(*,*)'ERROR(optimization_check):sd_ratio=',sd_ratio
      stop
    endif
!
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
    if (i_verbose >= 1) then
      if (log_unit>0) write(log_unit,"('force_ave,max=',I10,2F30.20,I10)")itemd,ddave,ddmax,js_max
    endif  
!
!   
    da1=sd_ratio
    da2=sd_ratio
    da3=sd_ratio
!
    ddsum=0.0d0
    ddsumx=0.0d0
    ddsumy=0.0d0
    ddsumz=0.0d0
    do js=1,noa
     ddrr1=foi(js,1)*da1
     ddrr2=foi(js,2)*da2
     ddrr3=foi(js,3)*da3
     if (iflag(js) == 1) then
       ddsumx=ddsumx+foi(js,1)*ddrr1
       ddsumy=ddsumy+foi(js,2)*ddrr2
       ddsumz=ddsumz+foi(js,3)*ddrr3
     endif  
    enddo   
!
    ddsum=ddsumx+ddsumy+ddsumz
!
    if (.not. allocated(sd_energy_diff)) then
      allocate(sd_energy_diff(4), stat=ierr) 
      if (ierr /= 0) then
        write(*,*)'Alloc. ERROR(optimization_check)' 
        stop
      endif   
    endif   
!
    sd_energy_diff(1)=ddsum
    sd_energy_diff(2)=ddsumx
    sd_energy_diff(3)=ddsumy
    sd_energy_diff(4)=ddsumz
!
    if (i_verbose >= 1) then
      if (log_unit>0) then 
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference (tot) [au] =',  itemd, ddsum
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference ( x ) [au] =',  itemd, ddsumx
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference ( y ) [au] =',  itemd, ddsumy
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference ( z ) [au] =',  itemd, ddsumz
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference per atom (tot) [eV] =', & 
                        & itemd, ddsum/dble(noa)*ev4au
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference per atom (x)   [eV] =', & 
                        & itemd, ddsumx/dble(noa)*ev4au
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference per atom (y)   [eV] =', & 
                        & itemd, ddsumy/dble(noa)*ev4au
        write(log_unit,'(a,i10,f30.20)') '  SD energy difference per atom (z)   [eV] =', & 
                        & itemd, ddsumz/dble(noa)*ev4au
      endif  
    endif  
!
    if (trim(config%calc%optimization%convergence_mode) == 'energy_per_atom') then
      if (config%calc%optimization%energy_convergence < 0.0d0) then
        if (log_unit > 0) write(log_unit,*) & 
             &  'ERROR(optimization_check) energy_convergence=',config%calc%optimization%energy_convergence
        stop
      endif   
      if (config%calc%optimization%energy_convergence > dabs(ddsum)/dble(noa)) then
        converged = .true.
        if (log_unit > 0) write(log_unit,*) & 
             &  'INFO:optimization_check: converged:criteria, energy difference [eV]=', &
             &   config%calc%optimization%energy_convergence*ev4au, dabs(ddsum)/dble(noa)*ev4au
      endif   
    endif   
!
  end subroutine optimization_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine for simple steepest descent (SD) method
!
  subroutine steepest_descent(sd_ratio)
!
    use elses_mod_ctrl,       only : i_verbose
    use elses_mod_mass,       only : amm
    use elses_mod_sim_cell,   only : noa, ax, ay, az
    use elses_mod_iflag,      only : iflag
    use elses_mod_txp,        only : txp, typ, tzp
    use elses_mod_md_dat, only : itemd
    use elses_mod_foi,    only : foi
!   use elses_mod_foiold, only : foiold
!
    implicit none
    real(8), intent(in) :: sd_ratio
    integer  js, js_max
    real(8)  ddave, ddmax, ddd, ddd1, ddd2, ddd3
    real(8)  da0, da1, da2, da3, ddrr1, ddrr2, ddrr3
    real(8)  ddsum
    logical :: update_atom_position
!
!
    da0=sd_ratio
    if (i_verbose >= 1) then
      write(*,*)'@@ Steepest Descent Method:Ratio=',da0
    endif  
!
    if (.not. allocated(iflag)) then
      write(*,*) 'ERROR(steepest_descent):not allocated:iflag' 
      stop
    endif   
!
    if (sd_ratio < -1.0d-10) then
      update_atom_position = .false.
      if (i_verbose >= 1) then
        write(*,*)'INFO:The atom positions are not updated (sd_ratio < 0)'
      endif  
    else
      update_atom_position = .true.
    endif
!
!   foiold(:,:)=foi(:,:)
!
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
    write(*,"('force_ave,max=',I10,2F30.20,I10)")itemd,ddave,ddmax,js_max
!
    da1=da0
    da2=da0
    da3=da0
!
    ddsum=0.0d0
    do js=1,noa
     ddrr1=foi(js,1)*da1
     ddrr2=foi(js,2)*da2
     ddrr3=foi(js,3)*da3
     if (iflag(js) == 1) then
       if (update_atom_position) then
         txp(js)=txp(js)+ddrr1/ax
         typ(js)=typ(js)+ddrr2/ay
         tzp(js)=tzp(js)+ddrr3/az
       endif
       ddsum=ddsum+foi(js,1)*ddrr1 &
&                 +foi(js,2)*ddrr2 &
&                 +foi(js,3)*ddrr3 
     endif  
!    if (i_verbose >= 1) then
!      write(44,"('SD: dx,dy,dz=',2I10,3F20.10)")itemd,js,ddrr1,ddrr2,ddrr3
!    endif  
    enddo   
!
    if (i_verbose >= 1) then
      write(*,*)'  SD: Energy Difference [au] =', ddsum
    endif  
!
    call elses_gene_tx
!
  end subroutine steepest_descent
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Generate TX from TXP
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine md_update_tx
    use elses_mod_ctrl,     only : i_verbose
    use elses_mod_sim_cell, only : noa, i_pbc_x, i_pbc_y, i_pbc_z
    use elses_mod_tx,       only : tx, ty, tz, jsei
    use elses_mod_txp,      only : txp, typ, tzp
    implicit none
    integer :: ioa
!
    if (i_verbose >= 1) then
      WRITE(6,*) '@@ GENETX:NOA,I_PBC_X=',noa,i_pbc_x, i_pbc_y, i_pbc_z
    endif  
!
    do ioa=1,noa
!
      if (i_pbc_x == 1) then
        tx(ioa)=txp(ioa)-dble(int(txp(ioa)))
        if(tx(ioa).lt.0.d0) then
          tx(ioa)=tx(ioa)+1.d0
        endif    
      else
        tx(ioa)=txp(ioa)
      endif   
!
      if (i_pbc_y == 1) then
        ty(ioa)=typ(ioa)-dble(int(typ(ioa)))
        if(ty(ioa).lt.0.d0) then
          ty(ioa)=ty(ioa)+1.d0
        endif
      else
        ty(ioa)=typ(ioa)
      endif   
!
      if (i_pbc_z == 1) then
        tz(ioa)=tzp(ioa)-dble(int(tzp(ioa)))
        if(tz(ioa).lt.0.d0) then
          tz(ioa)=tz(ioa)+1.d0
        endif
      else
       tz(ioa)=tzp(ioa)
      endif   
!
    enddo
!
  end subroutine md_update_tx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Funtion to get the distance between two atoms 
!         in a.u. 
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  function get_atom_dist(js1,js2) result(dd)
    use elses_mod_sim_cell,   only : i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az           !(unchanged)
    use elses_mod_tx,         only : tx, ty, tz                !(unchanged)
!
    implicit none
    real(8) ::  dxc, dyc, dzc, dd
    integer ::  ierr, js1, js2
!
    ierr=0
    if ((js1 <= 0) .or. (js1 > noa)) ierr=1
    if ((js2 <= 0) .or. (js2 > noa)) ierr=1
    if (ierr == 1) then
       write(*,*)'ERROR( get_atom_dist)'
       write(*,*)'js1, js2,noa=',js1,js2,noa
       stop
    endif
!   
    dxc=tx(js2)-tx(js1)
    dyc=ty(js2)-ty(js1)
    dzc=tz(js2)-tz(js1)
    if (i_pbc_x == 1) dxc = modulo(dxc + 0.5d0, 1.0d0) - 0.5d0
    if (i_pbc_y == 1) dyc = modulo(dyc + 0.5d0, 1.0d0) - 0.5d0
    if (i_pbc_z == 1) dzc = modulo(dzc + 0.5d0, 1.0d0) - 0.5d0
!
    dxc=dxc*ax
    dyc=dyc*ay
    dzc=dzc*az
    dd=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)
!
  end function get_atom_dist

end module M_md_motion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elses_gene_tx
  use M_md_motion, only : md_update_tx
  call md_update_tx
end subroutine elses_gene_tx
