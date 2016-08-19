!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_velocity_routines
!
  use M_config           ! (unchanged)
  use M_qm_domain,      only : i_verbose
  integer :: noa_mobile ! set in calc_kinetic_energy

  private
! public :: set_zero_velocity_for_fixed_atoms
  public :: calc_kinetic_energy
  public :: calc_initial_velocity
  public :: calc_total_momentum
  public :: adjust_velocity
  public :: allocate_velocity
  public :: calc_kin_ene_atom
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of kinetic energy of atom
!      Output : kinetic_energy [au]
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_kinetic_energy(kinetic_energy, k_e_component)
!
    use elses_mod_md_dat, only : itemd, dtmd          !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(unchanged)
!
    implicit none
    integer js
    real(8), intent(out)           :: kinetic_energy    ! total kinetic energy
    real(8), intent(out), optional :: k_e_component(3)  ! x,y,z components 
    real(8)  ddsum,ddd
    real(8)  ddsum_compo(3)
    logical :: velocity_is_allocated
    integer lu
!
    lu = config%calc%distributed%log_unit
!
    velocity_is_allocated = .true.
!
    if (present(k_e_component)) then
      k_e_component(:)=0.0d0
    endif  
!
    if (allocated(velx) .eqv. .false.) velocity_is_allocated = .false.
    if (allocated(vely) .eqv. .false.) velocity_is_allocated = .false.
    if (allocated(velz) .eqv. .false.) velocity_is_allocated = .false.
!
    if (.not. velocity_is_allocated) then
      kinetic_energy=0.0d0
      if (i_verbose >= 1) then
        if (lu > 0) then 
          write(lu,*)'@@ calc_kinetic_energy...is skipped'
        endif  
      endif  
      return
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      if (lu > 0) then 
        write(lu,*)'@@ calc_kinetic_energy'
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
    ddsum_compo(:)=0.d00
!$omp  parallel default(shared) &
!$omp& private (js) &
!$omp& reduction (+ : ddsum_compo) 
!$omp  do schedule(static)
    do js=1,noa
      if (dabs(amm(js)) .le. 1.0D-10) then
        write(*,*)'ERROR!(CALC_KIN_ENE):js,amm=',js,amm(js)
        stop
      endif   
      ddsum_compo(1)=ddsum_compo(1)+velx(js)*velx(js)/amm(js)*ax*ax 
      ddsum_compo(2)=ddsum_compo(2)+vely(js)*vely(js)/amm(js)*ay*ay
      ddsum_compo(3)=ddsum_compo(3)+velz(js)*velz(js)/amm(js)*az*az
    enddo    
!$omp end do
!$omp end parallel
    ddsum_compo(:)=ddsum_compo(:)*0.5d0/(dtmd*dtmd)
    ddsum=sum(ddsum_compo(:))
!
    kinetic_energy=ddsum
    if (present(k_e_component)) then
      k_e_component(:)=ddsum_compo(:)
    endif  
!
  end subroutine calc_kinetic_energy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of kinetic energy for each atom
!      Output : kinetic_energy [au]
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_kin_ene_atom(jsv, kinetic_energy, k_e_component)
!
    use elses_mod_md_dat,     only : dtmd      !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(unchanged)
!
    implicit none
    integer, intent(in)            :: jsv
    real(8), intent(out)           :: kinetic_energy    ! total kinetic energy
    real(8), intent(out), optional :: k_e_component(3)  ! x,y,z components 
    real(8)  ddsum,ddd
    real(8)  ddsum_compo(3)
    logical :: velocity_is_allocated
    integer js
!
    ddsum_compo(:)=0.d00
    js=jsv
!
    kinetic_energy=0.0d0
    if (present(k_e_component)) then
      k_e_component(:)=0.0d0
    endif  
!
    velocity_is_allocated = .true.
    if (allocated(velx) .eqv. .false.) velocity_is_allocated = .false.
    if (allocated(vely) .eqv. .false.) velocity_is_allocated = .false.
    if (allocated(velz) .eqv. .false.) velocity_is_allocated = .false.
    if (.not. velocity_is_allocated) return
!
    if ( (js < 1) .or. (js > noa) ) then
      write(*,*)'ERROR(calc_kinetic_energy_atom):jsv=',js
      stop
    endif
!
    if (dabs(amm(js)) .le. 1.0D-10) then
      write(*,*)'ERROR!(CALC_KIN_ENE):js,amm=',js,amm(js)
      stop
    endif   
!
    ddsum_compo(1)=ddsum_compo(1)+velx(js)*velx(js)/amm(js)*ax*ax 
    ddsum_compo(2)=ddsum_compo(2)+vely(js)*vely(js)/amm(js)*ay*ay
    ddsum_compo(3)=ddsum_compo(3)+velz(js)*velz(js)/amm(js)*az*az

    ddsum_compo(:)=ddsum_compo(:)*0.5d0/(dtmd*dtmd)
    ddsum=sum(ddsum_compo(:))
!
    kinetic_energy=ddsum
    if (present(k_e_component)) then
      k_e_component(:)=ddsum_compo(:)
    endif  
!
!
  end subroutine calc_kin_ene_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of total momemtum
!      Output : total_momemtum(3)
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_total_momentum(total_momentum)
!
!   use elses_mod_phys_const, only : ev4au,ev2kel     !(unchanged)
    use elses_mod_md_dat,     only : dtmd             !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(unchanged)
    use elses_mod_iflag,      only : iflag            !(unchanged)
    implicit none
    real(8),          intent(out) :: total_momentum(3)
    integer                       :: js
    real(8)                       :: mass
!
    total_momentum(1:3)=0.0d0
!
    if (.not. allocated(iflag)) then
      write(*,*)'ERROR(calc_total_momentum):not allocated : iflag' 
      stop
    endif   
!
!$omp  parallel default(shared) &
!$omp& private (js,mass) &
!$omp& reduction (+ : total_momentum)
!$omp  do schedule(static)
    do js=1,noa
      mass=1.0d0/amm(js)
      if (iflag(js) == 1) then
        total_momentum(1)=total_momentum(1)+mass*velx(js)*ax/dtmd
        total_momentum(2)=total_momentum(2)+mass*vely(js)*ay/dtmd
        total_momentum(3)=total_momentum(3)+mass*velz(js)*az/dtmd
      endif    
    enddo      
!$omp end do 
!$omp end parallel
!
  end subroutine calc_total_momentum
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of total momemtum as absolute value
!      Output : total_momemtum(3)
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_total_momentum_abs(total_momentum)
!
!   use elses_mod_phys_const, only : ev4au,ev2kel     !(unchanged)
    use elses_mod_md_dat,     only : dtmd             !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(unchanged)
    use elses_mod_iflag,      only : iflag            !(unchanged)
    implicit none
    real(8),          intent(out) :: total_momentum(3)
    integer                       :: js
    real(8)                       :: mass
!
    if (.not. allocated(iflag)) then
      write(*,*)'ERROR(calc_total_momentum_abs):not allocated : iflag' 
      stop
    endif   
!
    total_momentum(1:3)=0.0d0
!
!$omp  parallel default(shared) &
!$omp& private (js,mass) &
!$omp& reduction (+ : total_momentum)
!$omp  do schedule(static)
    do js=1,noa
      mass=1.0d0/amm(js)
      if (iflag(js) == 1) then
        total_momentum(1)=total_momentum(1)+mass*dabs(velx(js))*ax/dtmd
        total_momentum(2)=total_momentum(2)+mass*dabs(vely(js))*ay/dtmd
        total_momentum(3)=total_momentum(3)+mass*dabs(velz(js))*az/dtmd
      endif    
    enddo      
!$omp end do 
!$omp end parallel
!
  end subroutine calc_total_momentum_abs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Adjustment of the velocity so that (total momemtum) = 0
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine adjust_velocity
!
    use elses_mod_sim_cell,   only : noa, ax, ay, az   !(unchanged)
    use elses_mod_md_dat,     only : dtmd              !(unchanged)
    use elses_mod_mass,       only : amm               !(unchanged)
    use elses_mod_iflag,      only : iflag             !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz  !(unchanged)
    implicit none
    integer                       :: noac, js
    real(8)                       :: total_momentum(3)
    real(8)                       :: total_momentum_abs(3)
    real(8)                       :: d_px, d_py, d_pz
    real(8)                       :: dddx, dddy, dddz
    real(8)                       :: mass
    integer lu
!
    lu = config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*)' @@ Adjust velocity so that (total momemtum) = 0 : adjust_velocity'
    endif
!   
    call calc_total_momentum(total_momentum)
!
    call calc_total_momentum_abs(total_momentum_abs)
!
    dddx=total_momentum(1)/total_momentum_abs(1)
    dddy=total_momentum(2)/total_momentum_abs(2)
    dddz=total_momentum(3)/total_momentum_abs(3)
!
    if (i_verbose >= 1) then
      if (lu > 0) then
        write(lu,*)'(befor) total P_x =',total_momentum(1), dabs(dddx)
        write(lu,*)'(befor) total P_y =',total_momentum(2), dabs(dddy)
        write(lu,*)'(befor) total P_z =',total_momentum(3), dabs(dddz)
      endif
    endif
!
    noac=noa_mobile
!
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*)' (the number of the movable atoms) =',noac
    endif
!
    if (noac <= 0) then
      write(*,*)'ERROR:NOAC=',noac
      stop
    endif
!   
    d_px = total_momentum(1)/dble(noac)
    d_py = total_momentum(2)/dble(noac)
    d_pz = total_momentum(3)/dble(noac)

!$omp  parallel default(shared) &
!$omp& private (js,mass)
!$omp  do schedule(static)      
    do js=1,noa
      if (iflag(js) == 1) then
        mass=1.0d0/amm(js)
        velx(js)=(velx(js)-d_px*dtmd/ax/mass)*dble(iflag(js))
        vely(js)=(vely(js)-d_py*dtmd/ay/mass)*dble(iflag(js))
        velz(js)=(velz(js)-d_pz*dtmd/az/mass)*dble(iflag(js))
      endif    
    enddo      
!$omp end do
!$omp end parallel
!
    call calc_total_momentum(total_momentum)
!
    dddx=total_momentum(1)/total_momentum_abs(1)
    dddy=total_momentum(2)/total_momentum_abs(2)
    dddz=total_momentum(3)/total_momentum_abs(3)
!
    if (i_verbose >= 1) then
      if (lu > 0) then
        write(lu,*)'(after) total P_x =',total_momentum(1), dabs(dddx)
        write(lu,*)'(after) total P_y =',total_momentum(2), dabs(dddy)
        write(lu,*)'(after) total P_z =',total_momentum(3), dabs(dddz)
      endif
    endif
!
  end subroutine adjust_velocity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Check anc calculate the initial velocity
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_initial_velocity
!
!      NOTE:Input parameters: 
!          tempk0     : Temperature
!          amm(1:noa) : Mass of atoms
!
    use M_config,             only : config            !(unchanged)(only config%system%temperature)
    use elses_mod_phys_const, only : ev4au, ev2kel     !(parameter)
    use elses_mod_sim_cell,   only : noa, ax, ay, az   !(unchanged)
    use elses_mod_md_dat,     only : dtmd              !(unchanged)
    use elses_mod_mass,       only : amm               !(unchanged)
    use elses_mod_iflag,      only : iflag             !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz  !(CHANGED)
!   use elses_mod_thermo, &
!               only : tempk0,thb,thbold,vhb,vhbold
    use elses_mod_thermo,     only : thb,thbold,vhb,vhbold
!
    implicit none
    integer iseed, myrank, ioa, js, mm
    integer noac, iflagd
    real(8) sig, sum, rnd
    real(8) vxsum, vysum, vzsum, ekion, tempk, tempr
    real(8) :: kinetic_energy
    real(8) :: total_momentum(3)
    real(8) :: tempk0
    integer lu
    logical :: forced_initialization 
    logical :: parameter_checking
!
    noa_mobile=sum(iflag(1:noa))
!
    if (config%calc%distributed%dst_bench_mode) then
       forced_initialization = .true.
    else
       forced_initialization = .false.
    endif
!
    parameter_checking    = .false.
!
    lu = config%calc%distributed%log_unit
!
    tempk0 = config%system%temperature
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      if (lu > 0) then
        write(lu,*)'@@ calc_initial_velocity:checking initial velocity'
        write(lu,*)' temperature[au]:tempk0=',tempk0
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Parameter checking
!
    if (parameter_checking) then 
      if (tempk0 .le. -1.0d-10) then
        write(*,*)'ERROR!(calc_initial_velocity)'
        stop
      endif   
!
!$omp  parallel default(shared) &
!$omp& private (ioa)
!$omp  do schedule(static) 
      do ioa=1,noa
        if (dabs(amm(ioa)) .le. 1.0d-10) then
          write(*,*)'ERROR!(calc_initial_velocity)'
          write(*,*)'ioa, amm=',ioa,amm(ioa)
          stop
        endif   
      enddo   
!$omp end do
!$omp end parallel
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Detect the meaningful velocity data in the XML file
!
    if (.not. forced_initialization) then 
!
      call calc_kinetic_energy(kinetic_energy)
!
      call calc_total_momentum(total_momentum)
!
      if (i_verbose >= 1) then
        if (lu > 0) then
          write(lu,*)'Kinetic eneregy = ',kinetic_energy
          write(lu,*)' total momemtum in x direction =',total_momentum(1)
          write(lu,*)' total momemtum in y direction =',total_momentum(2)
          write(lu,*)' total momemtum in z direction =',total_momentum(3)
        endif
      endif
!
      myrank=0
!
      tempk=2.d0/3.d0*kinetic_energy/dble(noa_mobile)
!
      if(i_verbose >= 1) then
        if (lu > 0) write(lu,*) 'tempk (a.u., eV)= ', tempk, tempk*ev4au
      endif
!
      if (dabs(tempk0) .lt. 1.0d-10) then
        write(*,*)'ERROR:zero temperature:tempk0=',tempk0
        stop
      endif   
!
      tempr=tempk/tempk0
!
      if (dabs(tempr) .gt. 1.0d-10) then
        if (lu > 0) then
          write(lu,*) 'INFO:Optional XML tag detected : Meaningful velocity data (will be used)'
          write(lu,*) ' The temperature of the input velocity data [eV, kel]=', & 
  &                    tempk*ev4au, tempk*ev4au*ev2kel
          write(lu,*) ' thb, thbold=',thb, thbold
          write(lu,*) ' vhb, vhbold=',vhb, vhbold
          return
        endif
      endif   
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      if (lu > 0) then
        write(lu,*) '  INFO: Velocity is not set in the input file'
        write(lu,*) '  INFO: Heatbath parameters are set to be zero'
      endif
    endif
!
    vhb=0.0d0
    vhbold=0.0d0
    thb=0.0d0
    thbold=0.0d0
!
    if (i_verbose >= 1) then
      if (lu > 0) write(lu,*) '  ---> Velocity will be generated from the maxwell distribution'
    endif
    call set_velocity_from_maxwell
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call calc_kinetic_energy(kinetic_energy)
    tempk=2.d0/3.d0*kinetic_energy/dble(noa_mobile)
    tempr=tempk/tempk0
!
    if(i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'tempk (a.u., eV)= ', tempk, tempk*ev4au
      if (lu > 0) write(lu,*) 'tempk/tempk0 = ',tempr
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Adjust the velocity so that (Total Momemtum) = 0
!
    call adjust_velocity    
!
    call calc_kinetic_energy(kinetic_energy)
    tempk=2.d0/3.d0*kinetic_energy/dble(noa_mobile)
    tempr=tempk/tempk0
!
    if(i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'tempk (a.u., eV)= ', tempk, tempk*ev4au
      if (lu > 0) write(lu,*) 'tempk/tempk0 = ',tempr
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Scale the velocities to gurantee that
!       the total kinetic energy 
!             is exactly equal to the temper.
!
    if(i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'AFTER SCALING'
    endif
!      
    tempr=dsqrt(tempr)
!
    velx(:)=velx(:)/tempr
    vely(:)=vely(:)/tempr
    velz(:)=velz(:)/tempr
!
!   do ioa=1, noa
!     if (lu > 0) write(lu,*) 'velx,y,z=', ioa, velx(ioa), vely(ioa), velz(ioa)
!   enddo   
!
    call calc_kinetic_energy(kinetic_energy)
    tempk=2.d0/3.d0*kinetic_energy/dble(noa_mobile)
    tempr=tempk/tempk0
!
    if(i_verbose >= 1) then
      if (lu > 0) write(lu,*) 'tempk (a.u., eV)= ', tempk, tempk*ev4au
      if (lu > 0) write(lu,*) 'tempk/tempk0 = ',tempr
    endif
!
  end subroutine calc_initial_velocity
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Initialization of velocity
!          written by Hiroaki OHTANI; modified by T.Hoshi
!          ONLY FOR FREE ATOMS
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine set_velocity_from_maxwell
!
!     NOTE:Input parameters: 
!           tempk0     : Temperature
!           amm(1:noa) : Mass of atoms
!
    use M_config,             only : config           !(unchanged)(only config%system%temperature)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use elses_mod_md_dat,     only : dtmd             !(unchanged)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_iflag,      only : iflag            !(unchagend)
    use elses_mod_vel,        only : velx, vely, velz !(CHANGED)
!   use elses_mod_thermo, &
!&             only : tempk0,thb,thbold,vhb,vhbold
!
    use M_lib_random_num, only : rndini, rndu
    implicit none
    integer iseed, jinit30, myrank, ioa, js, mm
!   integer noac, iflagd
    real(8) dtmd2, sig, sum, rnd
    real(8) vxsum, vysum, vzsum, ekion, tempk, tempr
    real(8) tempk0
    integer lu
!
    lu = config%calc%distributed%log_unit
!
    tempk0=config%system%temperature
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      if (lu > 0) then
        write(lu,*)'@@ set_velocity_from_maxwell'
        write(lu,*)' temperature[au]:tempk0=',tempk0
        write(lu,*)'  noa_mobile =',noa_mobile
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Parameter checking
!
    if (tempk0 .le. -1.0d-10) then
      write(*,*)'ERROR!(set_velocity_from_maxwell)'
      stop
    endif   
!
!$omp  parallel default(shared) &
!$omp& private (ioa)
!$omp  do schedule(static)
    do ioa=1,noa
      if (dabs(amm(ioa)) .le. 1.0d-10) then
        write(*,*)'ERROR!(set_velocity_from_maxwell)'
        write(*,*)'ioa, amm=',ioa,amm(ioa)
        stop
      endif   
    enddo   
!$omp end do
!$omp end parallel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Preperation for rundom number generator
!
    iseed=13
!      : Seed for rundom number generator (should be odd)
    call rndini(iseed)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Generation of Velocity under Maxwell distribution
!
    velx(:)=0.0d0
    vely(:)=0.0d0
    velz(:)=0.0d0
!!
    do ioa=1,noa
!
       if (iflag(ioa) /= 1) cycle
!
       sig=tempk0*amm(ioa)
       sum = 0.d0
       do mm = 1, 12
         call rndu(rnd) 
         sum=sum+rnd
       enddo
       velx(ioa)=(sum-6.d0)*dsqrt(sig)
       velx(ioa)=velx(ioa)*dtmd/ax
!
       sum = 0.d0
       do mm = 1, 12
         call rndu(rnd) 
         sum=sum+rnd
       enddo  
       vely(ioa)=(sum-6.d0)*dsqrt(sig)
       vely(ioa)=vely(ioa)*dtmd/ay
!
       sum = 0.d0
       do mm = 1, 12
         call rndu(rnd) 
         sum=sum+rnd
       enddo  
       velz(ioa)=(sum-6.d0)*dsqrt(sig)
       velz(ioa)=velz(ioa)*dtmd/az
!
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  end subroutine set_velocity_from_maxwell
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Allocate the velocity
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine allocate_velocity
!
    use M_qm_domain,          only : i_verbose !(unchanged)
    use elses_mod_sim_cell,   only : noa  !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(CHANGED)
!
    implicit none
    integer :: ierr
!
    if (i_verbose >= 1) then
      write(*,*)'@@ allocate_velocity'
    endif  
!
    allocate (velx(noa),stat=ierr)
    if( ierr /= 0 ) stop 'Alloc. Error(velx)'
    velx(:)=0.0d0
!
    allocate (vely(noa),stat=ierr)
    if( ierr /= 0 ) stop 'Alloc. Error(vely)'
    vely(:)=0.0d0
!
    allocate (velz(noa),stat=ierr)
    if( ierr /= 0 ) stop 'Alloc. Error(velz)'
    velz(:)=0.0d0
!
  end subroutine allocate_velocity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_md_velocity_routines
