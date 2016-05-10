!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_ini_load_matom
!
  private
  public set_structure_data_from_matom
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Set strcuture data
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine set_structure_data_from_matom
!
    use M_config !(unchanged)
    use M_qm_domain,        only : atm_position, atm_element, noav !(CHANGED)
!   use elses_mod_tx,       only : jsei                !(CHANGED)
!   use elses_mod_txp,      only : txp, typ, tzp       !(CHANGED)
    use elses_mod_md_dat,   only : dtmd                !(unchagend)   
    use elses_mod_sim_cell, only : ax,ay,az            !(unchanged)
!
    implicit none
    integer :: j, k
    type(mini_atom_type), pointer :: atom
    real(8) :: dtmd_wrk
    integer :: ierr
    real(8) :: txpd, typd, tzpd
    real(8) :: tx,   ty,   tz
    integer :: n_element_wrk
!      
    integer :: log_unit, i_verbose
!    
    i_verbose=config%option%verbose
    log_unit=config%calc%distributed%log_unit
    dtmd_wrk=dtmd
    noav=config%system%structure%natom
    n_element_wrk=config%system%structure%nelement
!
!   write(log_unit,*) 'modulo( 2.1d0, 1.0d0) =', modulo( 1.1d0, 1.0d0)
!   write(log_unit,*) 'modulo(-1.1d0, 1.0d0) =', modulo(-0.1d0, 1.0d0)
!   write(log_unit,*) 'modulo( 1.1d0, 1.0d0) =', modulo( 1.1d0, 1.0d0)
!   write(log_unit,*) 'modulo(-0.1d0, 1.0d0) =', modulo(-0.1d0, 1.0d0)
!
    if (i_verbose >=1 ) then
      write(*,'(a)') '@@ set_structure_data_from_matom'
      write(*,*) 'INFO:nelement=', n_element_wrk
      write(*,*) 'INFO:Set noav = ', noav
      if (log_unit > 0 ) then
        write(log_unit,'(a)') '@@ set_structure_data_from_matom'
        write(log_unit,*) 'INFO:nelement=', n_element_wrk
        write(log_unit,*) 'INFO:Set noav = ', noav
      endif   
    endif  
!
    if ((n_element_wrk <= 0) .or. (n_element_wrk > 100)) then
      write(*,*)'ERROR(set_structure_data_from_matom):n_element=',n_element_wrk
      stop
    endif
!
    if (noav <= 0) then
      write(*,*)'ERROR(set_structure_data_from_matom):noav=',noav
      stop
    endif
!
    if (allocated(atm_position)) then
      write(*,*)'ERROR(set_structure_data_from_matom):atom_position is already allocated'
      stop
    else
      allocate(atm_position(3,noav),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'ERROR(set_structure_data_from_matom):Alloc. error: atom_position'
        stop
      endif   
      if (i_verbose >=1 ) then
        if (log_unit > 0 ) then
          write(log_unit,'(a,f20.10)')'INFO:Alloc. of atm_position : size [GB]  =', 8.0d0*3.0d0*dble(noav)/1.0d9
        endif
      endif  
    endif   
!
    if (allocated(atm_element)) then
      write(*,*)'ERROR(set_structure_data_from_matom):atom_position is already allocated'
      stop
    else
      allocate(atm_element(noav),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'ERROR(set_structure_data_from_matom):Alloc. error: atom_element'
        stop
      endif   
      atm_element(:)=0   ! dummy setting
      if (i_verbose >=1 ) then
        if (log_unit > 0 ) then
          write(log_unit,'(a,f20.10)')'INFO:Alloc. of atm_element  : size [GB]  =', 4.0d0*dble(noav)/1.0d9
        endif
      endif  
    endif   
!
    if (dabs(ax) < 1.0d-10) then
      write(*,*)'ERROR(set_structure_data_from_matom):ax=',ax
      stop
    endif   
!
    if (dabs(ay) < 1.0d-10) then
      write(*,*)'ERROR(set_structure_data_from_matom):ay=',ay
      stop
    endif   
!
    if (dabs(az) < 1.0d-10) then
      write(*,*)'ERROR(set_structure_data_from_matom):az=',az
      stop
    endif   
!
!$omp  parallel &
!$omp& default(none) &
!$omp& shared(atm_position, atm_element, config) &
!$omp& private(j, txpd, typd, tzpd, tx, ty, tz) &
!$omp& private(k, atom) &
!$omp& firstprivate(ax, ay, az)
!$omp  do schedule(static)
    do j=1, config%system%structure%natom
      atom => config%system%structure%matom(j)
!
      if (trim(atom%name) == '') then
        write(*,*) 'ERROR(set_structure_data_from_matom):empty name:j=',j
        stop
      endif   
!
      txpd = atom%position(1) /ax
      typd = atom%position(2) /ay
      tzpd = atom%position(3) /az
!
      tx = modulo(txpd, 1.0d0)
      ty = modulo(typd, 1.0d0)
      tz = modulo(tzpd, 1.0d0)
!
      atm_position(1,j)=tx
      atm_position(2,j)=ty
      atm_position(3,j)=tz
!
      do k = 1, config%system%structure%nelement+1
        if (k == config%system%structure%nelement+1) then
           write(*,*)'ERROR(set_structure_data_from_matom):j,atm_element(j)=',j, atm_element(j)
           write(*,*)'ERROR(set_structure_data_from_matom): element name  )=',trim(atom%name)
           write(*,*) '#ELSES: Stop by error:'
           write(*,*) '#This may be due to mismatch'
           write(*,*) '#  between the configuration XML file'
           write(*,*) '#  and the structure XML file'
           write(*,*) '# The element tag is required '
           write(*,*) '#  for each atom species'
           write(*,*) '#  in the configuration file.'
           stop
        endif
        if( trim(atom%name) == trim(config%system%structure%velement(k)%name) ) then
           atm_element(j)= k               
           exit
        end if
      end do
!
    end do
!$omp end do
!$omp end parallel
!
    if (log_unit > 0 ) then
       write(log_unit,'(a)') 'INFO:dealloc. of config%system%structure%matom'
    endif   
!
    deallocate(config%system%structure%matom, stat=ierr) 
    if (ierr /= 0) then
      write(*,*) 'ERROR in dealloc: config%system%structure%matom'
      stop
    endif   

!   stop 'Stop manually'
!
  end subroutine set_structure_data_from_matom
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module M_ini_load_matom


