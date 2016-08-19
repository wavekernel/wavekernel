!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_ini_load_vatom
!
  private
  public set_structure_data_from_vatom
  public set_atm_position_from_tx
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set atm_position from (tx,ty,tz)
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine set_atm_position_from_tx
!
    use M_config !(unchanged)
    use M_qm_domain,        only : noav  !(CHANGED)
    use M_qm_domain,        only : atm_position, atm_element !(CHANGED)
    use elses_mod_sim_cell, only : noa   ! (unchanged)
    use elses_mod_tx,       only : tx, ty, tz, jsei
    use elses_mod_txp,      only : txp, typ, tzp
    use M_qm_domain,        only : make_booking_list !(routine)
    use elses_mod_sim_cell,   only : i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
!
    implicit none
    integer :: jsv,js, nss, ierr
    logical :: calledFirst
    real(8) :: time_wrk, time_wrk_previous
    integer :: i_verbose
    integer :: i_pbc_all
!
    i_pbc_all= i_pbc_x*i_pbc_y*i_pbc_z 
!       = 1  if the PBC is imposed on the three directions
!       = 0  otherwise
    i_verbose=config%option%verbose
    noav=config%system%structure%natom
!
    if (i_verbose >= 1) then
      write(*,*)'@@ set_atm_position_from_tx'
      write(*,*)'  Set noav=',noav
    endif   
!
    if (noav <= 0) then
      write(*,*)'ERROR(set_atm_position_from_tx):noav=', noav
      stop
    endif   
!
    if( allocated(atm_position) ) then
       if( size(atm_position,1) .ne. 3 ) &
            stop 'set_atm_position_from_tx: allocation size mismatch atm_position 1'
       if( size(atm_position,2) .ne. noav ) &
            stop 'set_atm_position_from_tx: allocation size mismatch atm_position 2'
       if( .not. allocated(atm_element) ) &
            stop 'set_atm_position_from_tx: inconsistent state of allocation 1'
       if( size(atm_element,1) .ne. noav ) &
            stop 'set_atm_position_from_tx: allocation size mismatch atm_element'
       calledFirst=.false.
    else
       allocate (atm_position(3,noav),stat=ierr)
       if( ierr .ne. 0 ) then
          write(6,*)'alloc. error!(atom_position):ierr=',ierr
          stop
       endif
       if( allocated(atm_element) ) &
            stop 'set_atm_position_from_tx: inconsistent state of allocation 2'
       allocate (atm_element(noav),stat=ierr)
       if( ierr .ne. 0) then
          write(6,*)'alloc. error!(atom_element):ierr=',ierr
          stop
       endif
       calledFirst=.true.
    end if
!
    if (config%calc%distributed%dst_bench_mode) then
      if (i_pbc_all == 1) then 
!$omp  parallel default(none) &
!$omp& shared (jsei, atm_position, txp, typ, tzp, atm_element, noav) &
!$omp& private (jsv) 
!$omp  do schedule(static)
        do jsv=1,noav
          atm_position(1,jsv)=modulo(txp(jsv), 1.0d0)
          atm_position(2,jsv)=modulo(typ(jsv), 1.0d0)
          atm_position(3,jsv)=modulo(tzp(jsv), 1.0d0)
          atm_element(jsv)=jsei(jsv)
        enddo
!$omp  end do
!$omp  end parallel
      else
        write(*,*)'ERROR(set_atm_position_from_tx):unmatched in the dst_bench mode'
        stop
      endif
      return
    endif
!
    do jsv=1,noav
       js=jsv  !!! ASSUMED (Memo by T. Hoshi, 2010.Jan.02)
       nss=jsei(js)
       atm_position(1,jsv)=tx(js)
       atm_position(2,jsv)=ty(js)
       atm_position(3,jsv)=tz(js)
       atm_element(jsv)=nss
    enddo
!
  end subroutine set_atm_position_from_tx
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set strcuture data
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine set_structure_data_from_vatom
!
    use M_config !(unchanged)
!   use elses_mod_tx,       only : tx, ty, tz, jsei    !(CHANGED)
    use elses_mod_tx,       only : jsei                !(CHANGED)
    use elses_mod_txp,      only : txp, typ, tzp       !(CHANGED)
    use elses_mod_md_dat,   only : dtmd                !(unchagend)   
!   use elses_mod_sim_cell, only : noa,ax,ay,az, i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
    use elses_mod_sim_cell, only : ax,ay,az            !(unchanged)
    use elses_mod_vel,      only : velx, vely, velz    !(CHANGED) 
    use elses_mod_foi,      only : foi                 !(CHANGED)
!   use elses_mod_md_dat,   only : itemdorg            !(CHANGED)
!   use elses_mod_thermo,   only : amq, thb,thbold,vhb,vhbold !(CHANGED)
    use elses_mod_iflag,    only : iflag                      !(CHANGED)
!   use elses_param_ctl_kr,  only : nreclc_def, noak_min_def  !(CHANGED) 
!   use M_md_motion, only : elses_md_motion_ini_set           !(routine)
!
    implicit none
    integer :: j, k
    type(atom_type), pointer :: atom
!   real(8) :: r_mem_limit
!   real(8) :: ntlimit
!   integer :: noak_min_tmp, nrecl_tmp
!   real(8) :: r_hb_mass_per_atom
    real(8) :: dtmd_wrk
!      
    integer :: log_unit, i_verbose
!    
    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit
    dtmd_wrk=dtmd
!
    if (i_verbose >=1 ) then
      if (log_unit > 0 ) then
        write(log_unit,'(a)') '@@ set_structure_data_from_vatom'
      endif   
    endif  
!
    do j=1, config%system%structure%natom
      atom => config%system%structure%vatom(j)
!
!
      txp(j) = atom%position(1) /ax
      typ(j) = atom%position(2) /ay
      tzp(j) = atom%position(3) /az
!
      jsei(j)=0   ! dummy setting
!
      do k = 1, config%system%structure%nelement
        if( atom%name == config%system%structure%velement(k)%name ) then
           jsei(j)= k               
        end if
      end do
!
      if (jsei(j) == 0) then
         write(*,*)'ERROR(set_structure_data_from_vatom):j,jsei(j)=',j,jsei(j)
         write(*,*) '#ELSES: Stop by error:'
         write(*,*) '#This may be due to mismatch'
         write(*,*) '#  between the configuration XML file'
         write(*,*) '#  and the structure XML file'
         write(*,*) '# The element tag is required '
         write(*,*) '#  for each atom species'
         write(*,*) '#  in the configuration file.'
         stop
       endif   
!
!      if (allocated(iflag)) then
!        if( atom%motion == "free" ) then
!           iflag(j)=1
!        else
!           iflag(j)=0
!           if (log_unit > 0 ) then
!             write(log_unit,'(a,i15)') 'INFO:fixed atom: j=',j
!           endif   
!        end if
!      endif  
!         
       if (allocated(iflag)) then
         select case (atom%motion)
           case ("free") 
             iflag(j)=1
!            if (log_unit > 0 ) then
!              write(log_unit,'(a,i15)') 'INFO:free  atom: j=',j
!            endif   
           case ("fixed") 
             iflag(j)=0
             if (log_unit > 0 ) then
               write(log_unit,'(a,i15)') 'INFO:fixed atom: j=',j
             endif   
           case default
             write(*,*) 'ERROR:: atom id, motion=',j, trim(atom%motion)
             stop
         end select
       endif  
!         
       if( atom%velocity_set ) then
          if ( allocated(velx) ) velx(j) = atom%velocity(1)*dtmd_wrk/ax
          if ( allocated(vely) ) vely(j) = atom%velocity(2)*dtmd_wrk/ay
          if ( allocated(velz) ) velz(j) = atom%velocity(3)*dtmd_wrk/az
       end if
!
       if (allocated(foi)) then
         if( atom%force_set ) then
            foi(j,:) = atom%force(:)
         end if
       endif  
!
    end do
!
    call elses_gene_tx
!
  end subroutine set_structure_data_from_vatom
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module M_ini_load_vatom


