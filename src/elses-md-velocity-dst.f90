!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_velocity_dst
!
  use M_config                                    ! (unchanged)
  use M_qm_domain,        only : i_verbose
  use M_io_dst_write_log, only : log_unit
  implicit none
  real(8), allocatable :: vel_dst(:,:)    ! DIFFERENT AMONG NODES !!
  real(8), allocatable :: pos_dst(:,:)    ! DIFFERENT AMONG NODES !!
  integer, allocatable :: iflag_dst(:)    ! DIFFERENT AMONG NODES !!
  integer :: atm_index_ini, atm_index_fin, noa_dst ! DIFFERENT AMONG NODES !!
!                ----> determined in copy_velocity_to_vel_dst
  integer :: noa_mobile   ! number of mobile atoms. determined in copy_velocity_to_vel_dst 
!
  private
  public :: convert_velocity
  public :: calc_kinetic_energy_dst
  public :: md_motion_verlet_velocity_dst
  public :: md_motion_verlet_position_dst
!
! public :: calc_initial_velocity
! public :: calc_total_momentum
! public :: adjust_velocity
! public :: allocate_velocity
! public :: calc_kin_ene_atom
!
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Convert velocity
!       imode=1 : (velx, vely, velz) --> vel_dst 
!       imode=2 :       vel_dst      --> (velx, vely, velz)
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine convert_velocity(imode, mpi_elapse_time)
!    
    use M_config                                                ! (unchanged)
    use elses_mod_sim_cell, only : noa                          !(unchanged)
    use elses_mod_vel,      only : velx, vely, velz             !(CHANGED)
    use elses_mod_iflag,    only : iflag                        !(unchanged)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r1     !(routine)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_wtime            !(routine)
    use M_md_dst_get_atom_range, only : dst_get_atm_index_range !(routine)
    implicit none
    integer, intent(in)  :: imode
    real(8), optional    :: mpi_elapse_time
    integer :: i_v, lu
    integer :: ierr, dst_atm_index, atm_index
    logical, parameter :: debug_mode = .false.
!   logical, parameter :: debug_mode = .true.
    real(8)            :: timer_wrk1, timer_wrk2
!
    i_v = config%option%verbose_level
    lu  = config%calc%distributed%log_unit
!
    if (lu > 0) write(lu,*) '@@@ convert_velocity:imode=', imode
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Initial procedures
!
    if (.not. allocated(vel_dst)) then
      call dst_get_atm_index_range(atm_index_ini,atm_index_fin)
      noa_dst = atm_index_fin - atm_index_ini + 1
      noa_mobile = sum(iflag(1:noa))
      allocate(vel_dst(3,noa_dst), stat=ierr) 
      if (ierr /= 0) then
        write(*,*)'ALLOC ERROR:copy_velocity_to_vel_dst'
        stop
      endif
      allocate(iflag_dst(noa_dst), stat=ierr) 
      if (ierr /= 0) then
        write(*,*)'ALLOC ERROR:copy_velocity_to_vel_dst'
        stop
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

    if (debug_mode) then
      ierr=0
      if (atm_index_ini < 1) ierr=1
      if (atm_index_ini > noa) ierr=1
      if (atm_index_fin < 1) ierr=1
      if (atm_index_fin > noa) ierr=1
      if (atm_index_fin < atm_index_ini+1) ierr=1
      if (ierr == 1) then
        write(*,*)'ERROR(copy_velocity_to_velx):atm_index_ini, atm_index_fin=',atm_index_ini, atm_index_fin
        stop
      endif
      if (.not. allocated(velx)) ierr=1
      if (.not. allocated(vely)) ierr=1
      if (.not. allocated(velz)) ierr=1
      if (ierr == 1) then
        write(*,*)'ERROR(copy_velocity_to_velx):velx is not allocated'
        stop
      endif
    endif
!
    timer_wrk1=0.0d0
    timer_wrk2=0.0d0
!
    if (imode == 1) then
      vel_dst(:,:)=0.0d0
      do dst_atm_index=1,noa_dst
        atm_index = atm_index_ini - 1 + dst_atm_index
        vel_dst(1,dst_atm_index)=velx(atm_index)
        vel_dst(2,dst_atm_index)=vely(atm_index)
        vel_dst(3,dst_atm_index)=velz(atm_index)
        iflag_dst(dst_atm_index)=iflag(atm_index)
      enddo
    else
      velx(:)=0.0d0
      vely(:)=0.0d0
      velz(:)=0.0d0
      do dst_atm_index=1,noa_dst
        atm_index = atm_index_ini - 1 + dst_atm_index
        velx(atm_index)=vel_dst(1,dst_atm_index)
        vely(atm_index)=vel_dst(2,dst_atm_index)
        velz(atm_index)=vel_dst(3,dst_atm_index)
      enddo
      call mpi_wrapper_wtime(timer_wrk1)
      call mpi_wrapper_allreduce_r1(velx)
      call mpi_wrapper_allreduce_r1(vely)
      call mpi_wrapper_allreduce_r1(velz)
      call mpi_wrapper_wtime(timer_wrk2)
    endif
!
    if (present(mpi_elapse_time)) mpi_elapse_time=timer_wrk2-timer_wrk1

  end subroutine convert_velocity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of kinetic energy of atom
!      Output : kinetic_energy [au]
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_kinetic_energy_dst(kinetic_energy)
!
    use elses_mod_md_dat, only : itemd, dtmd          !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r0  !(routine)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_wtime         !(routine)
!
    implicit none
    logical, parameter :: debug_mode = .false.
!   logical, parameter :: debug_mode = .true.
    real(8), intent(out)           :: kinetic_energy    ! total kinetic energy
    real(8)  ddsum,ddd
    real(8)  mass
    integer lu, i_v
    integer :: ierr, dst_atm_index, atm_index
!
    i_v = config%option%verbose_level
    lu  = config%calc%distributed%log_unit
!
    if (.not. allocated(vel_dst)) then
      write(*,*)'ERROR(calc_kinetic_energy_dst):vel_dst is not allocated'
      stop
    endif
!
    kinetic_energy=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      if (lu > 0) then 
        write(lu,*)'@@ calc_kinetic_energy_dst'
      endif  
    endif  
!
    if (allocated(AMM) .eqv. .false.) then
      write(6,*)'ERROR!:Not yet allocated: AMM'
      stop
    endif 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    ddsum=0.0d0
!$omp  parallel default(shared) &
!$omp& private (dst_atm_index, atm_index, mass) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) &
!$omp& reduction (+ : ddsum) 
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
      atm_index = atm_index_ini - 1 + dst_atm_index
      mass=1.0d0/amm(atm_index)
      ddsum=ddsum+mass*(vel_dst(1,dst_atm_index)**2*ax*ax &
&                          +vel_dst(2,dst_atm_index)**2*ay*ay &
&                          +vel_dst(3,dst_atm_index)**2*az*az) 
    enddo    
!$omp end do
!$omp end parallel
    ddsum=ddsum*0.5d0/(dtmd*dtmd)
    call mpi_wrapper_allreduce_r0(ddsum)
!
    kinetic_energy=ddsum
!
  end subroutine calc_kinetic_energy_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Molecular Dynamics Part with heat bath
!       Velocity Verlet Method
!           written by T. Hoshi and H. OHTANI
!
!        when inose =1 --> canonical ensemble 
!                    0 --> micro-canonical
!
!          ivelini = 1 : The velocities are not updated
!                         (for the first MD iteration)
!          ivelini = 2 : The posisions are not updated
!                         (for the last MD iteration)
!           otherwise  : The posisions and velocites are updated
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine md_motion_verlet_velocity_dst(i_switch)
!
    use M_config,             only : config           !(unchanged)(only config%system%temperature)
    use M_qm_domain,          only : noav, ax,ay,az   !(unchanged)
    use elses_mod_phys_const, only : ev4au,ev2kel     !(parameter)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(CHANGED)
    use elses_mod_iflag,      only : iflag            !(unchanged)
!   use elses_mod_thermo,     only : tempk0,amq,thb,thbold,vhb,vhbold !(changed in the slaved routine)
    use elses_mod_thermo,     only : amq,thb,thbold,vhb,vhbold !(changed in the slaved routine)
    use elses_mod_md_dat,     only : itemd, itemdorg,itemdmx, dtmd !(unchanged)
    use elses_mod_md_dat,     only : e_kin                         !(CHANGED!)
    use elses_mod_foi,        only : foi
    use elses_mod_foiold,     only : foiold
    use M_lib_mpi_wrapper,    only : mpi_wrapper_wtime         !(routine)
!   use M_md_velocity_routines, only : calc_total_momentum !(routine)
!   use M_md_velocity_routines, only : adjust_velocity     !(routine)
!
    implicit none
    integer, intent(in) :: i_switch
    integer  ivelini, inose, ierr, myrank, iadjmom2
    integer  noa_def
    real(8)  dtmd2, err, tempk
!   integer  ioa
    real(8)  eki3, eki, ekion, sc1, sc2, ddvel, acc0
    real(8)  pkin, pkin2, dthb, x
    real(8)  ddrr1, ddrr2, ddrr3
    real(8)   :: total_momentum(3)
    real(8)   :: tempk0
!
    integer :: imode
    integer :: dst_atm_index, atm_index
    real(8) :: kinetic_energy_wrk
    real(8) :: mpi_elapse_time
    real(8) :: timer_now, timer_prev
    real(8) :: iflag_d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call mpi_wrapper_wtime(timer_now)
    timer_prev=timer_now
!
    tempk0=config%system%temperature
!
    if ( trim(config%calc%mode) /= "dynamics" ) then
      if( i_verbose >= 1 )then
        if (log_unit > 0) write(log_unit,*) '@@ md_motion_verlet ... is skipped'
      endif
      return
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (.not. allocated(iflag)) then
      write(*,*)'ERROR(md_motion_verlet_velocity):not allocated:iflag'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    myrank=1
    inose=1
    if (i_verbose >= 1) myrank=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note : in the micro-canonical ensemble,
!  all the heat-bath parameters are set to be zero.
!        ( vhb,vhbold,thb,thbold )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Caution
!      velocity:  velx := (velocity)*(time step)
!      foi: atomic unit
!      velx: reduced length unit
!      tx: reduced length unit
!      ekion : atomic unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) '@@ md_motion_verlet_velocity_dst:itemd,noa_dst=',itemd,noa_dst
    endif
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:ini1:', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Select the mode (ivelini)
!
    ivelini=0
!
    if (i_switch == 1) ivelini=1
!!
    if(i_verbose >= 1)then
      if (log_unit > 0) write(log_unit,*) 'ivelini,inose =',ivelini,inose
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set the heatbath variables 
!
    if (inose == 1) then
      if (ivelini /=1) call md_heatbath_dst
    else
      vhb=0.0d0
      vhbold=0.0d0
      thb=0.0d0
      thbold=0.0d0
    endif
!
    dtmd2=dtmd*dtmd
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:bath:', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  vleocity(t-dt) -> velocity(t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if(ivelini /= 1) then
!
      eki3=e_kin
      ekion=0.d0
      sc1=1.0d0/(1.0d0+0.5d0*vhb)
      sc2=1.0d0-0.5d0*vhbold

!$omp  parallel default(none) &
!$omp& shared  (iflag_dst, amm, foi, foiold, vel_dst) &
!$omp& private (dst_atm_index, atm_index, ddvel, iflag_d) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) &
!$omp& firstprivate (sc1, sc2, ax, ay, az, dtmd2) 
!$omp  do schedule(static)
      do dst_atm_index=1,noa_dst
        atm_index = atm_index_ini - 1 + dst_atm_index
        iflag_d = dble(iflag_dst(dst_atm_index))
        ddvel=0.5d0*dtmd2*amm(atm_index)*(foi(atm_index,1)+foiold(atm_index,1))/ax
        vel_dst(1,dst_atm_index)=sc1*(sc2*vel_dst(1,dst_atm_index)+ddvel)*iflag_d
        ddvel=0.5d0*dtmd2*amm(atm_index)*(foi(atm_index,2)+foiold(atm_index,2))/ay
        vel_dst(2,dst_atm_index)=sc1*(sc2*vel_dst(2,dst_atm_index)+ddvel)*iflag_d
        ddvel=0.5d0*dtmd2*amm(atm_index)*(foi(atm_index,3)+foiold(atm_index,3))/az
        vel_dst(3,dst_atm_index)=sc1*(sc2*vel_dst(3,dst_atm_index)+ddvel)*iflag_d
      enddo
!$omp end do
!$omp end parallel
!
    endif   
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:vel :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Check or adjust the velocity 
!         so as to keep the total momentum to be zero
!            the first variable =0 : check
!                               =1 : adjust
!
    if(ivelini /= 1) then
!
      iadjmom2=1
!     iadjmom2=0
!
      select case (iadjmom2)
        case (0)
          call calc_total_momentum_dst(total_momentum)
          if (i_verbose >= 1) then
            if (log_unit > 0) write(log_unit,*)'(check) total P_x =',total_momentum(1)
            if (log_unit > 0) write(log_unit,*)'(check) total P_y =',total_momentum(2)
            if (log_unit > 0) write(log_unit,*)'(check) total P_z =',total_momentum(3)
          endif
        case (1)
          call adjust_velocity_dst
      end select
!
    endif  
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:adj :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call calc_kinetic_energy_dst(kinetic_energy_wrk)
    e_kin=kinetic_energy_wrk
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:k_ene:', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
    ekion=e_kin/dble(noa_mobile)
    tempk=2.d0/3.d0*ekion
    if(i_verbose >= 1)then
      if (log_unit > 0) write(log_unit,*) 'ekion per noa (a.u.)= ', ekion
      if (log_unit > 0) write(log_unit,*) 'ekion per noa (eV  )= ', ekion*ev4au
      if (log_unit > 0) write(log_unit,*) 'ekion per noa (Kelv)= ', ekion*ev4au*ev2kel
      if (log_unit > 0) write(log_unit,*) 'tempk,tempk0,ratio  = ',  tempk,tempk0,dble(tempk/tempk0)
    endif 
    eki=ekion*dble(noa_mobile)
    e_kin=eki
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:k_ene2:', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
    if (.not. config%calc%distributed%dst_bench_mode) then
      imode=2
      call convert_velocity(imode, mpi_elapse_time)
      if (log_unit > 0) write(log_unit,*) 'TIME:mpi_time for convert velocity =', mpi_elapse_time
!
      call mpi_wrapper_wtime(timer_now)
      if( i_verbose >= 1 )then
        if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:conv:', timer_now-timer_prev
      endif
      timer_prev=timer_now
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Check the consistency of the velocity of heatbath 
!
    if ((ivelini /= 1) .and. (inose /= 0)) then
!
       pkin2=3.0d0*dble(noa_mobile)*tempk0
       x=vhbold/dtmd+0.5d0*dtmd*amq*(2.0d0*eki+2.0d0*eki3-2.0d0*pkin2)
       err=vhb-x*dtmd
       if(i_verbose >= 1)then
         if (log_unit > 0) write(log_unit,*)' comparison for vhb = ',vhb,x*dtmd
       endif
!
    endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  In the last md-iteration, only velocity V(t) is generated
!    thbold --> thb(t)
!
      thbold=thb
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:chk :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
  end subroutine md_motion_verlet_velocity_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Update the position of atom and the heat bath
!     using  Velocity Verlet Method
!
!        when inose =1 --> canonical ensemble 
!                    0 --> micro-canonical
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine md_motion_verlet_position_dst
!
    use M_config     !(unchanged)(only config%system%temperature)
    use M_qm_domain,          only : noav, ax,ay,az   !(unchanged)
    use elses_mod_phys_const, only : ev4au,ev2kel     !(parameter)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(unchanged)
    use elses_mod_iflag,      only : iflag            !(unchanged)
    use elses_mod_txp,        only : txp, typ, tzp    !(CHANGED)
!   use elses_mod_thermo,     only : tempk0,amq,thb,thbold,vhb,vhbold !(unchanged)
    use elses_mod_thermo,     only : amq,thb,thbold,vhb,vhbold !(unchanged)
    use elses_mod_md_dat,     only : itemd, itemdorg,itemdmx, dtmd, e_kin !(unchanged) 
    use elses_mod_foi,        only : foi        !(unchanged)
    use elses_mod_foiold,     only : foiold     !(CHANGED)
!   use M_md_velocity_routines, only : calc_total_momentum, adjust_velocity !(routine)
    use M_lib_mpi_wrapper,    only : mpi_wrapper_wtime         !(routine)
    use M_lib_mpi_wrapper,    only : mpi_wrapper_allreduce_r1  !(routine)
!
    implicit none
    integer  inose, ierr, myrank, iadjmom2
!   integer  ivelini, inose, ierr, myrank, iadjmom2
!   integer  noa_def
!   real(8)  dtmd2, err, tempk
    real(8)  dtmd2
!   real(8), allocatable :: foi2x(:),foi2y(:),foi2z(:)
    integer  ioa
    real(8)  eki3, eki, ekion, sc1, sc2, ddvel, acc0
    real(8)  pkin, pkin2, dthb
    real(8)  ddrr1, ddrr2, ddrr3
!
    real(8)  tempk0
!
    real(8)  :: foi2_wrk(3)
    real(8)  :: cell_length(3)
    real(8)  :: dr(3)
!
    integer :: i_v, lu
    integer :: dst_atm_index, atm_index
    real(8) :: mpi_elapse_time, timer_wrk1, timer_wrk2
    real(8) :: timer_now, timer_prev
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call mpi_wrapper_wtime(timer_now)
    timer_prev=timer_now
!
    i_v = config%option%verbose_level
    lu  = config%calc%distributed%log_unit
!
    inose=1
!
    cell_length(1)=ax
    cell_length(2)=ay
    cell_length(3)=az
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
    if (.not. allocated(iflag)) then
      write(*,*)'ERROR(md_motion_verlet_position):not allocated:iflag'
      stop
    endif   
!
    if (i_v > 0)then
      if (lu > 0) write(lu,*) '@@@ md_motion_verlet_position_dst:noa_dst=',noa_dst
    endif
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:ini :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Preparation of local variable and arrays
!
    tempk0=config%system%temperature
!
    dtmd2=dtmd*dtmd
!
    if (.not. allocated(pos_dst)) then 
      allocate (pos_dst(3,noa_dst),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(pos_dst):ierr=',ierr
        stop
      endif
    endif
    pos_dst(:,:)=0.0d0
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:alloc pos_dst :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!$omp  parallel default(none) &
!$omp& shared  (pos_dst, txp, typ, tzp) &
!$omp& private (dst_atm_index, atm_index) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) 
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
      atm_index = atm_index_ini - 1 + dst_atm_index
      pos_dst(1,dst_atm_index)=txp(atm_index)
      pos_dst(2,dst_atm_index)=typ(atm_index)
      pos_dst(3,dst_atm_index)=tzp(atm_index)
    enddo   
!$omp end do
!$omp end parallel
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:copy to pos_dst :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Force including the interaction with the heat bath
!    (note that VHB have already updated at HEATBATH)
!        force = foi(t) - m v_eta(t) * v(t)
!        position(t) -> position(t+dt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!$omp  parallel default(none) &
!$omp& shared  (amm, foi, vel_dst, iflag, pos_dst) &
!$omp& private (dst_atm_index, atm_index, foi2_wrk, dr, acc0 ) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst, cell_length, vhb, dtmd2)
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
      atm_index = atm_index_ini - 1 + dst_atm_index
      acc0=vhb/dtmd2/amm(atm_index)
      foi2_wrk(1:3)=foi(atm_index,1:3)-vel_dst(1:3,dst_atm_index)*acc0*cell_length(1:3)
      dr(1:3)=dble(iflag(atm_index))*(vel_dst(1:3,dst_atm_index)  &
&            +0.5d0*dtmd2*amm(atm_index)*foi2_wrk(1:3)/cell_length(1:3))
      pos_dst(1:3,dst_atm_index)=pos_dst(1:3,dst_atm_index)+dr(1:3)
    enddo   
!$omp end do
!$omp end parallel
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:upd pos_dst :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Copy the variable : pos_dst --> (txp, typ, tzp) and (tx, ty, tz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    txp(:)=0.0d0
    typ(:)=0.0d0
    tzp(:)=0.0d0
!
!$omp  parallel default(none) &
!$omp& shared  (pos_dst, txp, typ, tzp) &
!$omp& private (dst_atm_index, atm_index) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) 
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
      atm_index = atm_index_ini - 1 + dst_atm_index
      txp(atm_index)=pos_dst(1,dst_atm_index)
      typ(atm_index)=pos_dst(2,dst_atm_index)
      tzp(atm_index)=pos_dst(3,dst_atm_index)
    enddo   
!$omp end do
!$omp end parallel
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:copy from pos_dst :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
    call mpi_wrapper_allreduce_r1(txp)
    call mpi_wrapper_allreduce_r1(typ)
    call mpi_wrapper_allreduce_r1(tzp)
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:conv pos_dst :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
    foiold(:,:)=foi(:,:)
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:copy foi :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (.not. config%calc%distributed%dst_bench_mode) then
      call elses_gene_tx
!       ---> Periodic boundary condition, if you like
!
      call mpi_wrapper_wtime(timer_now)
      if( i_verbose >= 1 )then
        if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:gene tx :', timer_now-timer_prev
      endif
      timer_prev=timer_now
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  heatbath(t) -> heatbath(t+dt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (inose /= 0) then 
!
!     ekion=0.0d0
!     do ioa=1,noav
!       if (amm(ioa) >= 1.0d-10) then
!         ekion=ekion+velx(ioa)*velx(ioa)/amm(ioa)*ax*ax &
!&                    +vely(ioa)*vely(ioa)/amm(ioa)*ay*ay &
!&                    +velz(ioa)*velz(ioa)/amm(ioa)*az*az
!        endif
!      enddo  
!     ekion=0.5d0*ekion/dtmd2
!
      ekion=e_kin
      pkin=2.0d0*ekion
      pkin2=3.0d0*dble(noa_mobile)*tempk0
!
      pkin2=3.0d0*dble(noa_mobile)*tempk0
      dthb=vhb+0.5d0*dtmd2*amq*(pkin-pkin2)
      thb=thb+dthb
!
      if(i_verbose > 0)then
        if (log_unit > 0) write(log_unit,*) 'pkin2=',pkin2
        if (log_unit > 0) write(log_unit,*) 'dpkin,vhb,dthb=',pkin-pkin2,vhb,dthb
        if (log_unit > 0) write(log_unit,*) 'position of heatbath=',thb,thbold
      endif
!
    endif  
!
    call mpi_wrapper_wtime(timer_now)
    if( i_verbose >= 1 )then
      if (log_unit > 0) write(log_unit,*) 'TIME:md_motion_verlet_velocity_dst:bath :', timer_now-timer_prev
    endif
    timer_prev=timer_now
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  end subroutine md_motion_verlet_position_dst
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Determine the velocity of the heat bath
!          using Newton-Raphson method
!           partial compatible to Miyata's code
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine md_heatbath_dst
!
    use M_config,           only : config           !(unchanged)(only config%system%temperature)
!   use elses_mod_ctrl,     only : i_verbose
    use elses_mod_sim_cell, only : noa, ax, ay, az
    use elses_mod_mass,     only : amm
    use elses_mod_vel, only : velx, vely, velz
!   use elses_mod_thermo,   only : tempk0,amq,thb,thbold,vhb,vhbold
    use elses_mod_thermo,   only : amq,thb,thbold,vhb,vhbold
    use elses_mod_md_dat,   only : itemd, itemdorg,itemdmx, dtmd, e_kin
    use elses_mod_foi,    only : foi
    use elses_mod_foiold, only : foiold
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r0  !(routine)
!
    implicit none
    integer ivelini
    real(8) dtmd2, ams, x0, x, ekion, pkin, pkin2
    real(8) dsum, tkin, ac0, accex, accey, accez
    real(8) velx2, vely2, velz2
    real(8) velx3, vely3, velz3
    real(8) ac1, ddd, ddd1, etafnc, etafnd, detafnc
    real(8) dconv_cri
    integer iii,iii_max
    real(8) tempk0
!
    integer :: ierr, dst_atm_index, atm_index
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    tempk0=config%system%temperature
!
    dtmd2=dtmd*dtmd
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Select the mode
!
    ivelini=0
!
!   if ((itemdorg .eq. 0) .and. (itemd .eq. 1)) then
!     ivelini=1
!   endif   
!
!   if (itemd .eq. itemdmx) then
!     ivelini=2
!     if (itemdmx == 1) then
!        foiold(:,:)=foi(:,:)
!        if (i_verbose >= 1) then
!          write(*,*)'NOTE:foiold=foi, (itemdmx=1)'
!        endif
!     endif   
!   endif   
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
    if (i_verbose >= 1) then 
      if (log_unit > 0) write(log_unit,*) '@@ md_heatbath_dst:itemd,ivelini,amq=', itemd,ivelini,amq
    endif
!
    if (amq .lt. 1.0d-10) then 
      write(6,*) 'ERROR!(LSES_MD_HEATBATH)'
      write(6,*) ' amq=',amq
      stop
    endif
!
    ams=1.0d0/amq
!        Mass of heat bath
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   if (ivelini .eq. 1) then
!      vhb=0.0d0
!      vhbold=0.0d0
!      thb=0.0d0
!      thbold=0.0d0
!      return
!   endif
!
    x0=vhb/dtmd
!        OLD velocity of heat bath 
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   ekion=0.0d0
!   do ioa=1,noa
!      if (amm(ioa) .lt. 1.0d-10) cycle
!          ekion=ekion+velx(ioa)*velx(ioa)/amm(ioa)*ax*ax &
!&                    +vely(ioa)*vely(ioa)/amm(ioa)*ay*ay &
!&                    +velz(ioa)*velz(ioa)/amm(ioa)*az*az
!   enddo
!  ekion=0.5d0*ekion/dtmd2
!   eki=ekion
!
    ekion=e_kin
    pkin=2.0d0*ekion
    pkin2=3.0d0*dble(noa_mobile)*tempk0
!     
!      ekion : kinetic energy in a.u. : T = \sum_I 0.5 M_I (V_I)^2
!
    dsum=0.0d0
    ac0=1.0d0-0.5d0*dtmd*x0
!$omp  parallel default(none) &
!$omp& shared  (iflag_dst, foi, foiold, amm, vel_dst) &
!$omp& private (dst_atm_index, atm_index) &
!$omp& private (accex, accey, accez) &
!$omp& private (velx2, vely2, velz2) &
!$omp& private (velx3, vely3, velz3) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) &
!$omp& firstprivate (ac0, dtmd, ax, ay, az) &
!$omp& reduction (+ : dsum) 
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
       if ( iflag_dst(dst_atm_index) == 0 ) cycle
       atm_index = atm_index_ini - 1 + dst_atm_index
       accex=(foi(atm_index,1)+foiold(atm_index,1))*amm(atm_index)
       accey=(foi(atm_index,2)+foiold(atm_index,2))*amm(atm_index)
       accez=(foi(atm_index,3)+foiold(atm_index,3))*amm(atm_index)
       velx2=vel_dst(1,dst_atm_index)*ax/dtmd
       vely2=vel_dst(2,dst_atm_index)*ay/dtmd
       velz2=vel_dst(3,dst_atm_index)*az/dtmd
       velx3=ac0*velx2+0.5d0*dtmd*accex
       vely3=ac0*vely2+0.5d0*dtmd*accey
       velz3=ac0*velz2+0.5d0*dtmd*accez
       dsum=dsum+(velx3)**2/amm(atm_index)
       dsum=dsum+(vely3)**2/amm(atm_index)
       dsum=dsum+(velz3)**2/amm(atm_index)
     enddo
!$omp end do
!$omp end parallel
     call mpi_wrapper_allreduce_r0(dsum)
     tkin=dsum
!    
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Initial guess for new variable
!       assuming the ionic kinetic-energies are same
!                between the old and the present configuration
!
     x=x0+0.5d0*dtmd*amq*2.0d0*(pkin-pkin2)
     if (i_verbose >= 1) then 
       if (log_unit > 0) write(log_unit,*) 'guess of vhb,prev,dpkin=',x,x0,pkin-pkin2
     endif
!
     dconv_cri=1.0d-12
     iii_max=5000
     do iii=1,iii_max
       ac1=1.d0+0.5d0*dtmd*x
       ddd=pkin+tkin/(ac1**2)-2.d0*3.0d0*dble(noa)*tempk0
       ddd1=tkin/(ac1**3)*dtmd
       etafnc=x-x0-dtmd/(2.d0*ams)*ddd
       etafnd=1.d0+dtmd/(2.d0*ams)*ddd1
       detafnc= dabs(etafnc)
       if (i_verbose >= 1) then 
         if (log_unit > 0) write(log_unit,'(a,i10,4f30.20)')' iii,x,diff=',iii,x,etafnc,etafnd,ac1
       endif
       x=x-etafnc/etafnd
       if (detafnc .lt. dconv_cri) exit
       if (iii == iii_max) then
         write(*,*) 'ABORT(heat bath):do not converge'
         write(*,*) 'Criteria for convergence maybe too tight'
         write(*,*) '   (particularly in MP calc)'
         write(*,*) '     dconv_cri=',dconv_cri
         write(*,*) 'Detailed INFO:x0     =',x0
         write(*,*) 'Detailed INFO:dtmd   =',dtmd
         write(*,*) 'Detailed INFO:pkin   =',pkin
         write(*,*) 'Detailed INFO:tkin   =',tkin
         write(*,*) 'Detailed INFO:noa    =',noa
         write(*,*) 'Detailed INFO:tempk0 =',tempk0
         write(*,*) 'Detailed INFO:ac1    =',ac1
         write(*,*) 'NOTE:The parameter ac1 is a stalibity parameter'
         stop
        endif
      enddo
      if (i_verbose >= 1) then 
        if (log_unit > 0) write(log_unit,'(a,i10,3f30.20)') 'velocity of heatbath (new,old) =',iii,x,x0,ac1
      endif
      vhbold=vhb
!
      vhb=x*dtmd
!
!     stop 'STOP MANUALLY (md_heatbath_dst)'
!
  end subroutine md_heatbath_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of total momemtum
!      Output : total_momemtum(3)
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_total_momentum_dst(total_momentum)
!
!   use elses_mod_phys_const, only : ev4au,ev2kel     !(unchanged)
    use elses_mod_md_dat,     only : dtmd             !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
!   use elses_mod_vel,        only : velx, vely, velz !(unchanged)
    use M_lib_mpi_wrapper,    only : mpi_wrapper_allreduce_r1  !(routine)
!
    implicit none
    real(8),          intent(out) :: total_momentum(3)
    integer                       :: js
    real(8)                       :: mass
!
    integer :: ierr, dst_atm_index, atm_index
!
    total_momentum(1:3)=0.0d0
!
!$omp  parallel default(none) &
!$omp& shared  (amm, vel_dst) &
!$omp& private (dst_atm_index, atm_index, mass) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) &
!$omp& firstprivate (ax, ay, az, dtmd) &
!$omp& reduction (+ : total_momentum) 
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
      atm_index = atm_index_ini - 1 + dst_atm_index
      mass=1.0d0/amm(atm_index)
      total_momentum(1)=total_momentum(1)+mass*vel_dst(1,dst_atm_index)*ax
      total_momentum(2)=total_momentum(2)+mass*vel_dst(2,dst_atm_index)*ay
      total_momentum(3)=total_momentum(3)+mass*vel_dst(3,dst_atm_index)*az
    enddo
!$omp end do 
!$omp end parallel
!
    total_momentum(1:3)=total_momentum(1:3)/dtmd
!
    call mpi_wrapper_allreduce_r1(total_momentum)
!
  end subroutine calc_total_momentum_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Adjustment of the velocity so that (total momemtum) = 0
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine adjust_velocity_dst
!
    use elses_mod_sim_cell,   only : noa, ax, ay, az   !(unchanged)
    use elses_mod_md_dat,     only : dtmd              !(unchanged)
    use elses_mod_mass,       only : amm               !(unchanged)
    use elses_mod_iflag,      only : iflag             !(unchanged)
    implicit none
    real(8)                       :: total_momentum(3)
    real(8)                       :: d_px, d_py, d_pz
    real(8)                       :: mass_inv
    integer                       :: i_v, lu, imode
!
    integer :: ierr, dst_atm_index, atm_index
!
    i_v = config%option%verbose_level
    lu  = config%calc%distributed%log_unit
!
    if (i_v >= 1) then
      if (lu > 0) write(lu,*)' @@ Adjust velocity so that (total momemtum) = 0'
    endif
!   
    call calc_total_momentum_dst(total_momentum)
!
    if (i_v >= 1) then
      if (lu > 0) then
        write(lu,*)'(befor) total P_x =',total_momentum(1)
        write(lu,*)'(befor) total P_y =',total_momentum(2)
        write(lu,*)'(befor) total P_z =',total_momentum(3)
      endif
    endif
!
!
    if (i_v >= 1) then
      if (lu > 0) then
        write(lu,*)' (the number of the movable atoms) =',noa_mobile
      endif
    endif
!
    if ((noa_mobile < 1) .or. (noa_mobile > noa)) then
      write(*,*)'ERROR:noa_mobile=',noa_mobile
      stop
    endif
!   
    d_px = total_momentum(1)/dble(noa_mobile)*dtmd/ax
    d_py = total_momentum(2)/dble(noa_mobile)*dtmd/ay
    d_pz = total_momentum(3)/dble(noa_mobile)*dtmd/az
!
!$omp  parallel default(none) &
!$omp& shared  (amm, vel_dst, iflag) &
!$omp& private (dst_atm_index, atm_index) &
!$omp& firstprivate (atm_index_ini, atm_index_fin, noa_dst) &
!$omp& firstprivate (mass_inv, d_px, d_py, d_pz) 
!$omp  do schedule(static)
    do dst_atm_index=1,noa_dst
      atm_index = atm_index_ini - 1 + dst_atm_index
      mass_inv=amm(atm_index)
      vel_dst(1,dst_atm_index)=(vel_dst(1,dst_atm_index)-d_px*mass_inv)*dble(iflag(atm_index))
      vel_dst(2,dst_atm_index)=(vel_dst(2,dst_atm_index)-d_py*mass_inv)*dble(iflag(atm_index))
      vel_dst(3,dst_atm_index)=(vel_dst(3,dst_atm_index)-d_pz*mass_inv)*dble(iflag(atm_index))
    enddo      
!$omp end do
!$omp end parallel
!
    call calc_total_momentum_dst(total_momentum)
!
    if (i_v >= 1) then
      if (lu > 0) then
        write(lu,*)'(after) total P_x =',total_momentum(1)
        write(lu,*)'(after) total P_y =',total_momentum(2)
        write(lu,*)'(after) total P_z =',total_momentum(3)
      endif
    endif
!
  end subroutine adjust_velocity_dst
!
!
end module M_md_velocity_dst
