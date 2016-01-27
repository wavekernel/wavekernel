!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_verlet
!
  use M_qm_domain,      only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_unit                  !(unchanged)
! 
  private
  public :: md_motion_verlet_velocity
  public :: md_motion_verlet_position
!
  contains
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
  subroutine md_motion_verlet_velocity(i_switch)
!
    use M_config,             only : config           !(unchanged)(only config%system%temperature)
    use M_qm_domain,          only : noav, ax,ay,az   !(unchanged)
    use elses_mod_phys_const, only : ev4au,ev2kel     !(parameter)
    use elses_mod_mass,       only : amm              !(unchanged)
    use elses_mod_vel,        only : velx, vely, velz !(CHANGED)
    use elses_mod_iflag,      only : iflag            !(unchanged)
!   use elses_mod_thermo,     only : tempk0,amq,thb,thbold,vhb,vhbold !(changed in the slaved routine)
    use elses_mod_thermo,     only : amq,thb,thbold,vhb,vhbold !(changed in the slaved routine)
    use elses_mod_md_dat,     only : itemd, itemdorg,itemdmx, dtmd, e_kin
    use elses_mod_foi,        only : foi
    use elses_mod_foiold,     only : foiold
    use M_md_velocity_routines, only : calc_total_momentum, adjust_velocity !(routine)
!
    implicit none
    integer, intent(in) :: i_switch
    integer  ivelini, inose, ierr, myrank, iadjmom2
    integer  noa_def
    real(8)  dtmd2, err, tempk
    integer  ioa
    real(8)  eki3, eki, ekion, sc1, sc2, ddvel, acc0
    real(8)  pkin, pkin2, dthb, x
    real(8)  ddrr1, ddrr2, ddrr3
    real(8)   :: total_momentum(3)
    real(8)   :: tempk0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      if (log_unit > 0) write(log_unit,*) '@@ md_motion_verlet:itemd=',itemd
    endif
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
      if (ivelini /=1) call md_heatbath
    else
      vhb=0.0d0
      vhbold=0.0d0
      thb=0.0d0
      thbold=0.0d0
    endif
!
    dtmd2=dtmd*dtmd
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
      do ioa=1,noav
! 
       if (iflag(ioa) .eq. 1) then 
        ddvel=dble(iflag(ioa))*(0.5d0*dtmd2*amm(ioa)*(foi(ioa,1)+foiold(ioa,1))/ax)
        velx(ioa)=sc1*(sc2*velx(ioa)+ddvel)
!
        ddvel=dble(iflag(ioa))*(0.5d0*dtmd2*amm(ioa)*(foi(ioa,2)+foiold(ioa,2))/ay)
        vely(ioa)=sc1*(sc2*vely(ioa)+ddvel)
!
        ddvel=dble(iflag(ioa))*(0.5d0*dtmd2*amm(ioa)*(foi(ioa,3)+foiold(ioa,3))/az)
        velz(ioa)=sc1*(sc2*velz(ioa)+ddvel)
!
        if (dabs(amm(ioa)) .gt. 1.0D-10) then
          ekion=ekion+velx(ioa)*velx(ioa)/amm(ioa)*ax*ax &
&                    +vely(ioa)*vely(ioa)/amm(ioa)*ay*ay &
&                    +velz(ioa)*velz(ioa)/amm(ioa)*az*az
        endif  
!
       endif
! 
      enddo
!
      ekion=0.5d0*ekion/dtmd2/dble(noav)
      tempk=2.d0/3.d0*ekion
      if(myrank.eq.0)then
        if (log_unit > 0) write(log_unit,*) 'ekion per noa (a.u.)= ', ekion
        if (log_unit > 0) write(log_unit,*) 'ekion per noa (eV  )= ', ekion*ev4au
        if (log_unit > 0) write(log_unit,*) 'ekion per noa (Kelv)= ', ekion*ev4au*ev2kel
        if (log_unit > 0) write(log_unit,*) 'tempk,tempk0,ratio  = ',  tempk,tempk0,dble(tempk/tempk0)
      endif 
      eki=ekion*dble(noav)
      e_kin=eki
!
    endif   
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
          call calc_total_momentum(total_momentum)
          if (i_verbose >= 1) then
            if (log_unit > 0) write(log_unit,*)'(check) total P_x =',total_momentum(1)
            if (log_unit > 0) write(log_unit,*)'(check) total P_y =',total_momentum(2)
            if (log_unit > 0) write(log_unit,*)'(check) total P_z =',total_momentum(3)
          endif
        case (1)
          call adjust_velocity
      end select
!
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Check the consistency of the velocity of heatbath 
!
    if ((ivelini /= 1) .and. (inose /= 0)) then
!
       pkin2=3.0d0*dble(noav)*tempk0
       x=vhbold/dtmd+0.5d0*dtmd*amq*(2.0d0*eki+2.0d0*eki3-2.0d0*pkin2)
       err=vhb-x*dtmd
       if(myrank.eq.0)then
         if (log_unit > 0) write(log_unit,*)' comparison for vhb = ',vhb,x*dtmd
       endif
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  In the last md-iteration, only velocity V(t) is generated
!    thbold --> thb(t)
!
      thbold=thb
!
  end subroutine md_motion_verlet_velocity
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
  subroutine md_motion_verlet_position
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
!
    implicit none
    integer  inose, ierr, myrank, iadjmom2
!   integer  ivelini, inose, ierr, myrank, iadjmom2
!   integer  noa_def
!   real(8)  dtmd2, err, tempk
    real(8)  dtmd2
    real(8), allocatable :: foi2x(:),foi2y(:),foi2z(:)
    integer  ioa
    real(8)  eki3, eki, ekion, sc1, sc2, ddvel, acc0
    real(8)  pkin, pkin2, dthb
    real(8)  ddrr1, ddrr2, ddrr3
!   real(8)   :: total_momentum(3)
    real(8)  tempk0
!
!
    inose=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (.not. allocated(iflag)) then
      write(*,*)'ERROR(md_motion_verlet_position):not allocated:iflag'
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Preparation of local variable and arrays
!
    tempk0=config%system%temperature
!
    dtmd2=dtmd*dtmd
!
    allocate (foi2x(noav),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'alloc. error!(foi2x:ierr=',ierr
      stop
    endif
!
    allocate (foi2y(noav),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'alloc. error!(foi2y:ierr=',ierr
      stop
    endif
!
    allocate (foi2z(noav),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'alloc. error!(foi2z):ierr=',ierr
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Force including the interaction with the heat bath
!    (note that VHB have already updated at HEATBATH)
!        force = foi(t) - m v_eta(t) * v(t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    do ioa=1,noav
      foi2x(ioa)=foi(ioa,1)
      foi2y(ioa)=foi(ioa,2)
      foi2z(ioa)=foi(ioa,3)
      if (amm(ioa) >= 1.0d-10) then
         acc0=vhb/dtmd2/amm(ioa)
         foi2x(ioa)=foi2x(ioa)-velx(ioa)*acc0*ax
         foi2y(ioa)=foi2y(ioa)-vely(ioa)*acc0*ay
         foi2z(ioa)=foi2z(ioa)-velz(ioa)*acc0*az
      endif   
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  heatbath(t) -> heatbath(t+dt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (inose /= 0) then 
!
      ekion=0.0d0
      do ioa=1,noav
        if (amm(ioa) >= 1.0d-10) then
          ekion=ekion+velx(ioa)*velx(ioa)/amm(ioa)*ax*ax &
&                    +vely(ioa)*vely(ioa)/amm(ioa)*ay*ay &
&                    +velz(ioa)*velz(ioa)/amm(ioa)*az*az
        endif
      enddo  
!
      ekion=0.5d0*ekion/dtmd2
      pkin=2.0d0*ekion
      pkin2=3.0d0*dble(noav)*tempk0
!
      pkin2=3.0d0*dble(noav)*tempk0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  position(t) -> position(t+dt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    do ioa=1,noav
!
      ddrr1=dble(iflag(ioa))*(velx(ioa)+0.5d0*dtmd2*amm(ioa)*foi2x(ioa)/AX)
      ddrr2=dble(iflag(ioa))*(vely(ioa)+0.5d0*dtmd2*amm(ioa)*foi2y(ioa)/AY)
      ddrr3=dble(iflag(ioa))*(velz(ioa)+0.5d0*dtmd2*amm(ioa)*foi2z(ioa)/AZ)
!
      txp(ioa)=txp(ioa)+ddrr1
      foiold(ioa,1)=foi(ioa,1)
!
      typ(ioa)=typ(ioa)+ddrr2
      foiold(ioa,2)=foi(ioa,2)
!
      tzp(ioa)=tzp(ioa)+ddrr3
      foiold(ioa,3)=foi(ioa,3)
!
    enddo  
!
    deallocate (foi2x,stat=ierr)
    if( ierr .ne. 0) stop 'dealloc. error!:foi2x'
!
    deallocate (foi2y,stat=ierr)
    if( ierr .ne. 0) stop 'dealloc. error!:foi2y'
!
    deallocate (foi2z,stat=ierr)
    if( ierr .ne. 0) stop 'dealloc. error!:foi2z'
!
    call elses_gene_tx
!       ---> Periodic boundary condition, if you like
!
  end subroutine md_motion_verlet_position
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Determine the velocity of the heat bath
!          using Newton-Raphson method
!           partial compatible to Miyata's code
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine md_heatbath
!
    use M_config,           only : config           !(unchanged)(only config%system%temperature)
    use elses_mod_ctrl,     only : i_verbose
    use elses_mod_sim_cell, only : noa, ax, ay, az
    use elses_mod_mass,     only : amm
    use elses_mod_vel, only : velx, vely, velz
!   use elses_mod_thermo,   only : tempk0,amq,thb,thbold,vhb,vhbold
    use elses_mod_thermo,   only : amq,thb,thbold,vhb,vhbold
    use elses_mod_md_dat,   only : itemd, itemdorg,itemdmx, dtmd, e_kin
    use elses_mod_foi,    only : foi
    use elses_mod_foiold, only : foiold
!
    implicit none
    integer ivelini
    real(8) dtmd2, ams, x0, x, ekion, pkin, pkin2
    real(8) dsum, tkin, ac0, accex, accey, accez
    real(8) velx2, vely2, velz2
    real(8) velx3, vely3, velz3
    real(8) ac1, ddd, ddd1, etafnc, etafnd, detafnc
    real(8) dconv_cri
    integer ioa,iii,iii_max
    real(8) tempk0
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
      if (log_unit > 0) write(log_unit,*) '@@ md_heatbath:itemd,ivelini,amq=', itemd,ivelini,amq
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
    ekion=0.0d0
    do ioa=1,noa
       if (amm(ioa) .lt. 1.0d-10) cycle
          ekion=ekion+velx(ioa)*velx(ioa)/amm(ioa)*ax*ax &
&                    +vely(ioa)*vely(ioa)/amm(ioa)*ay*ay &
&                    +velz(ioa)*velz(ioa)/amm(ioa)*az*az
   enddo
    ekion=0.5d0*ekion/dtmd2
!   eki=ekion
    pkin=2.0d0*ekion
    pkin2=3.0d0*dble(noa)*tempk0
!     
!      ekion : kinetic energy in a.u. : T = \sum_I 0.5 M_I (V_I)^2
!
    dsum=0.0d0
    ac0=1.0d0-0.5d0*dtmd*x0
    do ioa=1,noa
       if (amm(ioa) .lt. 1.0d-10) cycle
       accex=(foi(ioa,1)+foiold(ioa,1))*amm(ioa)
       accey=(foi(ioa,2)+foiold(ioa,2))*amm(ioa)
       accez=(foi(ioa,3)+foiold(ioa,3))*amm(ioa)
       velx2=velx(ioa)*ax/dtmd
       vely2=vely(ioa)*ay/dtmd
       velz2=velz(ioa)*az/dtmd
       velx3=ac0*velx2+0.5d0*dtmd*accex
       vely3=ac0*vely2+0.5d0*dtmd*accey
       velz3=ac0*velz2+0.5d0*dtmd*accez
       dsum=dsum+(velx3)**2/amm(ioa)
       dsum=dsum+(vely3)**2/amm(ioa)
       dsum=dsum+(velz3)**2/amm(ioa)
     enddo
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
  end subroutine md_heatbath


end module M_md_verlet
