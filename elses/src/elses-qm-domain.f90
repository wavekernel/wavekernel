!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_domain
!  geno w.g history
!  Sep.16, 2008, S. Yamamoto, Added CSC arrays.
!  Oct.11, 2008, T. Hoshi, Added energy density matrix (dpij)
!
   use elses_mod_ctrl,     only : i_verbose
   use elses_mod_sel_sys,  only : c_system
   use elses_mod_jsv4jsd,  only : jsv4jsd, njsd
   use elses_mod_orb1,     only : nval
   use elses_mod_noav,     only : noav
   use elses_mod_sim_cell, only : NOS, ax, ay, az, i_pbc_x, i_pbc_y, i_pbc_z
   use elses_mod_elem_name,only : element_name => elem_name 
!
   use elses_mod_elec_cond,  only : temp_for_electron
   use elses_mod_val_elec,   only : total_electron_number => val_elec_tot
!
   use elses_arr_dhij,     only : dhij ! hamiltonian matrix
   use elses_arr_dsij,     only : dsij ! overlap matrix
   use elses_arr_dbij,     only : dbij ! density matrix
   use elses_arr_dpij,     only : dpij ! energy density matrix
!
   implicit none
   integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
   real(DOUBLE_PRECISION), allocatable :: atm_position(:,:)
!
   real(DOUBLE_PRECISION), allocatable :: atm_force(:,:), atm_force_csc(:,:)
   real(DOUBLE_PRECISION), allocatable :: atm_force_tb0(:,:)
!
   integer, allocatable :: atm_element(:)
!    
   real(DOUBLE_PRECISION) :: chemical_potential
!
   real(DOUBLE_PRECISION), allocatable :: ham_tb0(:,:,:,:)
!
   real(DOUBLE_PRECISION), allocatable :: ham_csc(:,:,:,:)
!
   real(DOUBLE_PRECISION), allocatable ::  gamma_csc(:,:)       ! \gamma(atom1, atom2)
   real(DOUBLE_PRECISION), allocatable :: dgamma_csc(:,:,:)
   real(DOUBLE_PRECISION), allocatable :: fgamma_csc(:,:)
!
   real(DOUBLE_PRECISION), allocatable :: e_num_on_basis(:,:) ! (Mulliken Charge) / (par atomic orbital), (orbital,atom) 
   real(DOUBLE_PRECISION), allocatable :: previous_e_num_on_basis(:,:) ! For charge self-consistent (CSC) loop 
   real(DOUBLE_PRECISION), allocatable :: delta_e_num(:), e_num_on_atom(:), tau_csc(:), previous_e_num_on_atom(:)
!
   real(DOUBLE_PRECISION), allocatable :: ddhij(:,:,:,:,:), ddsij(:,:,:,:,:), dham_tb0(:,:,:,:,:)
!
   integer                :: csc_ite_converged
   real(DOUBLE_PRECISION) :: csc_dq_converged
!
   logical, allocatable   :: CalledFirst(:)
!
   public
!  private :: make_booking_list
!
   logical  :: S_is_I
!
   real(DOUBLE_PRECISION) :: E_vdW_only   ! Energy only for vdW interaction
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_domain_setting
!
!   module variables unchanged
!      --> 
!   module variables changed
!      --> 
!
     use elses_mod_sim_cell, only : noa   ! (unchanged)
     use elses_mod_tx,       only : tx, ty, tz, jsei
!
     implicit none
     integer :: neig0, jsv,js, nss, ierr
!    logical :: calledFirst

     noav=noa
     if (i_verbose >= 1) then
       write(*,*)'  Set noav=',noav
     endif   

     if( allocated(atm_position) ) then
        if( size(atm_position,1) .ne. 3 ) &
             stop 'qm_domain_setting: allocation size mismatch atm_position 1'
        if( size(atm_position,2) .ne. noav ) &
             stop 'qm_domain_setting: allocation size mismatch atm_position 2'
        if( .not. allocated(atm_element) ) &
             stop 'qm_domain_setting: inconsistent state of allocation 1'
        if( size(atm_element,1) .ne. noav ) &
             stop 'qm_domain_setting: allocation size mismatch atm_element'
     else
        allocate (atm_position(3,noav),stat=ierr)
        if( ierr .ne. 0 ) then
           write(6,*)'alloc. error!(atom_position):ierr=',ierr
           stop
        endif
        if( allocated(atm_element) ) &
             stop 'qm_domain_setting: inconsistent state of allocation 2'
        allocate (atm_element(noav),stat=ierr)
        if( ierr .ne. 0) then
           write(6,*)'alloc. error!(atom_element):ierr=',ierr
           stop
        endif
     end if
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
     if (.not. allocated(CalledFirst)) then
       call qm_domain_ini_setting
       allocate(CalledFirst(1))
     endif  
!
     atm_force(1:3,1:noav)=0.0d0
!
     call  make_booking_list
!
   end subroutine qm_domain_setting
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_domain_ini_setting
!       called only at the first MD step
!
!   module variables unchanged
!      --> noa, 
!   module variables changed here
!      --> noav, noao, noas, noab, noab0, nncut
!
     use elses_mod_sel_sys,  only : r_cut_book !(unchanged)
     use elses_mod_sim_cell, only : noa ! (unchanged)
     use elses_mod_noav,     only : noao,noab,noas,noab0,nncut ! (CHANGED)
     use M_config,           only : config !(unchanged)
!
     implicit none
     integer ierr
     integer nval_max
     real(DOUBLE_PRECISION) :: cutoff_radius
     integer length_of_list
     logical :: non_orderN_memory
!
     if (i_verbose >= 1) then
       write(*,*)'@@ qm_domain_ini_setting'
       write(*,*)'  solver scheme=',trim(config%calc%solver%scheme)
     endif
!
     non_orderN_memory = .false. 
     if (trim(config%calc%solver%scheme) == 'eigen') then
        non_orderN_memory = .true. 
     endif   
!
     call qm_domain_ini_compat
!
     nval_max=maxval(nval)
     nncut=2
!
     cutoff_radius=r_cut_book

     call get_length_of_list(cutoff_radius,length_of_list)
!
     if (i_verbose >= 1) then
       write(*,*)'INFO:length of list =',length_of_list
       if (non_orderN_memory) then
         write(*,*)'INFO: NOAO=NOAV (Non-order-N memory for eigen solver)'
       else  
         write(*,*)'INFO: NOAO is set to be min (noav, 2*length_of_list) at the first MD step'
       endif  
     endif
!
     if (non_orderN_memory) then
       noao  = noav
     else  
       noao  = min(length_of_list*2,noav)
     endif  
!
     noas  = noao
     noab  = noao
     noab0 = noas
!
!    noao=min(2000,noav)
!    noas=min(2000,noav)
!    noab=min(2000,noav)
!    noab0=noas
!
     call elses_init_alloc2b
!      --> allocation : 
!           js4jsv(noav), jsv4js(noa), jsv4js(noa), jsv4js(noa),
!           dbij(nvl,nvl,noas,noav), dhij(nvl,nvl,noas,noav),
!           dhij(nvl,nvl,noas,noav), 
!           js4jsv(noav), jsv4js(noa), 
!
     
     call elses_init_alloc_nrl_leg
!      --> allocation : dsij, dpij
!          set : js4jsv(noav), jsv4js(noa)
!
!
     allocate (atm_force(3,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(atom_force):ierr=',ierr
        stop
     endif
!
     if (c_system == "geno") then
!
       allocate (ddhij(3,nval_max,nval_max,noas,noav),stat=ierr)
        if( ierr .ne. 0) then
          write(6,*)'alloc. error!(ddhij):ierr=',ierr
          stop
       endif
!
       allocate (ddsij(3,nval_max,nval_max,noas,noav),stat=ierr)
        if( ierr .ne. 0) then
          write(6,*)'alloc. error!(ddsij):ierr=',ierr
          stop
       endif
!
       call qm_domain_allocate_csc (nval_max, noas, noav)
!
     endif
!  
   end subroutine qm_domain_ini_setting
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Initial allocation for CSC method
!       called only at the first MD step
!
   subroutine qm_domain_allocate_csc (nval_max, noas, noav)
!
     use M_config, only : config !(unchanged)
     implicit none
     integer, intent(in) :: nval_max, noas, noav
     integer :: ierr
     integer :: n_csc_loop
!
!    real(DOUBLE_PRECISION) :: cutoff_radius
!    integer length_of_list
!
     n_csc_loop=config%calc%genoOption%CSC_max_loop_count
!
     allocate (ham_tb0(nval_max, nval_max, noas, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(ham_tb0):ierr=',ierr
        stop
     endif
!
     allocate (ham_csc(nval_max, nval_max, noas, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(ham_csc):ierr=',ierr
        stop
     endif
!
     if (n_csc_loop >0 ) then
       write(6,*)'INFO: gamma_csc is allocated (and may be a large array)'
       allocate (gamma_csc(noav, noav),stat=ierr)
       if( ierr .ne. 0) then
         write(6,*)'alloc. error!(gamma_csc):ierr=',ierr
         stop
       endif
       allocate (dgamma_csc(3,noav,noav),stat=ierr)
       if( ierr .ne. 0) then
         write(6,*)'alloc. error!(dgamma_csc):ierr=',ierr
         stop
       endif
     else
       write(6,*)'INFO: gamma_csc is NOT allocated'
     endif
!
     allocate (e_num_on_basis(nval_max, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(e_num_on_basis):ierr=',ierr
        stop
     endif
!
     allocate (previous_e_num_on_basis(nval_max, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(previous_e_num_on_basis):ierr=',ierr
        stop
     endif
!
     allocate (e_num_on_atom(noav),stat=ierr)
      if( ierr .ne. 0) then
          write(6,*)'alloc. error!(e_num_on_atom):ierr=',ierr
          stop
      end if
!
      allocate (previous_e_num_on_atom(noav),stat=ierr)
      if( ierr .ne. 0) then
         write(6,*)'alloc. error!(previous_e_num_on_atom):ierr=',ierr
         stop
      end if
!
     allocate (delta_e_num(noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(delta_e_num):ierr=',ierr
        stop
     endif
!
     allocate (tau_csc(noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(tau_csc):ierr=',ierr
        stop
     endif
!
     allocate (dham_tb0(3,nval_max,nval_max,noas,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(dham_tb0):ierr=',ierr
        stop
     endif
!
     allocate (atm_force_csc(3,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(atm_force_csc):ierr=',ierr
        stop
     endif

     allocate (fgamma_csc(noas,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(fgamma_csc):ierr=',ierr
        stop
     endif

   end subroutine qm_domain_allocate_csc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Booking the neighboring list 
!
!         Output : jsv4jsd(noao,noav)
!                : njsd(noav,0:nncut)
!
!         non-order-N calculation (temporaly) !!!
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
   subroutine make_booking_list
!
    use M_config,           only : config !(unchanged)
!     only:  config%calc%calc_check%set
!            config%calc%calc_check%short_atom_pair_distance%set
!            config%calc%calc_check%short_atom_pair_distance%warning_level
!            config%calc%calc_check%short_atom_pair_distance%abort_level
!
     use elses_mod_sel_sys,    only : r_cut_book
     use elses_mod_noav,       only : noao, nncut
!
     implicit none
     real(DOUBLE_PRECISION) :: rcut_bk
!
     integer                :: ishow,ict
     real(DOUBLE_PRECISION) :: rcut1
!
!    integer :: ipe,npe ! Variable for debugging. Currently not available
     integer                :: jsv1,jsv2,js1d,js1b,jsd
     real(DOUBLE_PRECISION) :: dvecx,dvecy,dvecz,w
!
     real(DOUBLE_PRECISION) :: min_atom_pair_distance ! minimum atom-pair distance
     integer                :: jsv1_min,jsv2_min      ! atom pair for the minimum distance 
!
     ishow=5
!
     ict=1
!
     if (i_verbose >= 1) then
       write(6,*)'@@@ make_booking_list'
       write(6,*)'   noav=',noav
       write(6,*)'r_cut_book =',r_cut_book
     endif  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Check the parapmeters
!
     if (nncut /= 2) then 
       write(6,*)'ERROR!(BKMATV)nncut=',nncut
       stop
     endif
!
     rcut_bk=r_cut_book
     if (rcut_bk <= 1.0d0) then 
       write(6,*)'ERROR!(BKMATV)rcut_bk=',rcut_bk
       stop
     endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     if( rcut_bk .lt. huge(1d0) ) then
        rcut1=rcut_bk*1.001d0
     else
        rcut1=rcut_bk
     end if
!       -->  radius  for cutoff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!$omp  parallel & 
!$omp& default(shared) &
!$omp& private(dvecx,dvecy,dvecz,w) & 
!$omp& private(jsv1,jsv2,js1d,js1b,jsd)
!-- Not available -- $omp& private(ipe,npe) & 
!    ipe=omp_get_thread_num()+1
!    npe=omp_get_num_threads()
!    write(6,*)'ipe,npe=',ipe,npe
!$omp  do schedule(static)
     do jsv2=1,noav
!
       jsv4jsd(:,jsv2)=0
!         
       jsv4jsd(1,jsv2)=jsv2
       njsd(jsv2,0)=1
!          ----> zero-th shell 
!
       js1d=1
       js1b=1
       ict=1
!   
       jsv1_min=0                      ! dummu value
       min_atom_pair_distance=1.0d100  ! dummy value
!
       do jsv1=1,noav
!        dd=1.0d0
         dvecx=atm_position(1,jsv2)-atm_position(1,jsv1)
         dvecy=atm_position(2,jsv2)-atm_position(2,jsv1)
         dvecz=atm_position(3,jsv2)-atm_position(3,jsv1)
!
         if (i_pbc_x == 1) dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
         if (i_pbc_y == 1) dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
         if (i_pbc_z == 1) dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0
!
         dvecx=dvecx*ax
         dvecy=dvecy*ay
         dvecz=dvecz*az
         w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
         if (jsv1 /= jsv2) then
           if (w < min_atom_pair_distance) then
             jsv1_min=jsv1
             min_atom_pair_distance=w 
           endif   
         endif
!
         if ((jsv1 /= jsv2) .and. (w .lt. rcut1)) then
           js1d=js1d+1
           if (js1d .gt. noao) then
             write(6,*)'error!(bkmatv2):js1d,noao=',js1d,noao
             stop
           endif   
           jsv4jsd(js1d,jsv2)=jsv1
         endif
       enddo  
       njsd(jsv2,ict)=js1d
!
       if ((config%calc%calc_check%set) .and. (config%calc%calc_check%short_atom_pair_distance%set)) then
         if (min_atom_pair_distance < config%calc%calc_check%short_atom_pair_distance%abort_level  ) then
           write(*,'(a,2i10, f20.10)') 'ABORT   :Too short atom-pair distance:', &
&                                        jsv2, jsv1_min, min_atom_pair_distance
           stop
         endif
         if (min_atom_pair_distance < config%calc%calc_check%short_atom_pair_distance%warning_level) then
           write(*,'(a,2i10, f20.10)') 'WARNING :Too short atom-pair distance:', &
&                                        jsv2, jsv1_min, min_atom_pair_distance
         endif
       endif   
!
     enddo
!$omp end do
!$omp end parallel 
!     
      write(*,*)' ...booking is done'
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Plot the result (optional)
!
     if (i_verbose .ge. 1) then
       do jsv2=1,noav
         if (jsv2 .le. ishow) then
           ict=1
           write(6,*)'jsv,njsd=',jsv2,njsd(jsv2,ict)
         endif   
       enddo
     endif  
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   end subroutine make_booking_list
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine generate_foi
!
!   module variables unchanged
!      --> js4jsv
!   module variables changed
!      --> foi
!
     use elses_mod_js4jsv,   only : js4jsv 
     use elses_mod_foi,      only : foi
!
     implicit none
     integer jsv,js
!
     if (i_verbose >= 1) then
       write(*,*)'@@ generate_foi'
     endif   
!
     foi(:,:)=0.0d0
!
     do jsv=1,noav
       js=js4jsv(jsv)
       foi(js,1:3)=atm_force(1:3,jsv)
     enddo
!   
   end subroutine generate_foi
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calculate Mulliken charge :  q_i = sum_j S_ij Rho_ji   
!         (S_ij, Rho_ij are symmetric matrices)
!
subroutine mulliken
!
! Module variables unchanged
!   ---> noav, njsd, atm_element, ksv4jsd, dbij, dsij
! Module variables changed
!   ---> None
!
! use elses_arr_dsij,     only : dsij
!
  use elses_mod_phys_const, only : para_spin_factor
!
  implicit none
  integer, parameter :: ict=1
  integer :: jsv1, jsv2, njsd2, nss1, nss2, nval1, nval2
  integer :: jsd1, ja1, ja2
  real(DOUBLE_PRECISION) :: d_mulliken, dbijd, dsijd, e_num_temp
  integer :: ierr
!
  ierr=0
  if (.not. allocated(dsij)) ierr=1
  if (.not. allocated(dbij)) ierr=1
  if (.not. allocated(e_num_on_basis)) ierr=1
  if (.not. allocated(e_num_on_atom)) ierr=1
  if (ierr /= 0) then
    write(*,*)'ERROR:mulliken'
    stop
  endif
!
  do jsv2=1,noav
    njsd2=njsd(jsv2,ict)
    nss2=atm_element(jsv2)
    nval2=nval(nss2)
    e_num_temp=0d0
    do ja2=1,nval2
      d_mulliken=0.0d0
      do jsd1=1,njsd2
        jsv1=jsv4jsd(jsd1,jsv2)
        nss1=atm_element(jsv1)
        nval1=nval(nss1)
        do ja1=1,nval1
          dbijd=dbij(ja1,ja2,jsd1,jsv2)
          dsijd=dsij(ja1,ja2,jsd1,jsv2)
          d_mulliken=d_mulliken+para_spin_factor*dbijd*dsijd
!            : para_spin_factor(=2.0d0) is multiplied
        enddo
      enddo
      e_num_on_basis(ja2,jsv2) = d_mulliken ! Nov 30, 2008
      e_num_temp=e_num_temp+d_mulliken
    enddo
    e_num_on_atom(jsv2)=e_num_temp
  enddo
!
  if (i_verbose >= 10) then
     do jsv2=1,noav
        njsd2=njsd(jsv2,ict)
        nss2=atm_element(jsv2)
        nval2=nval(nss2)
        do ja2=1,nval2
           write(*,'(" mulliken: jsv2,ja2,mulliken charge (Occ. for each orb.)=",&
                &2I8,ES21.14)')  jsv2,ja2, e_num_on_basis(ja2,jsv2)
        end do
           write(*,'(" mulliken: jsv2,ja2,mulliken charge (Occ. for each atom)=",&
                &I8,"        ",ES21.14)')  jsv2, e_num_on_atom(jsv2)
     end do
     write(*,'(" mulliken: Total electron number=",ES21.14)') sum(e_num_on_atom)
  endif  
!
end subroutine mulliken

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Subroutine for compatibility to legacy ('non-geno') code
!
subroutine qm_domain_ini_compat
   use elses_mod_sim_cell, only : noa
   use elses_mod_orb1,     only : nvl
   use elses_mod_eig_leg,  only : n_base_eig_leg
   implicit none
   integer   :: neig0
!
   if (i_verbose >= 1) then
     write(*,*)'@@ qm_domain_ini_compat'
     write(*,*)'c_system=',c_system
   endif  
!
   if (c_system == "NRL_leg") then
     nval(:)=9
     nvl=9
     n_base_eig_leg=nvl*noa
     if (i_verbose >= 1) then
       write(*,*)'  nvl =',nvl
       write(*,*)'  n_base_eig_leg =',n_base_eig_leg
     endif  
     neig0=n_base_eig_leg
!    call elses_alloc_eig_leg(neig0)
!    call elses_alloc_eig_leg_wrk
   endif
!
end subroutine qm_domain_ini_compat
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Get the length of the booking list for a given radius
!       
!
!         non-order-N calculation (temporaly) !!!
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
   subroutine get_length_of_list(cutoff_radius,length_of_list)
!
     use M_config,           only : config !(unchanged)
     implicit none
     integer, intent(out) :: length_of_list
     real(DOUBLE_PRECISION), intent(in)  :: cutoff_radius
!
     integer, allocatable :: number_of_members(:)
     integer i_show, ierr
     integer jsv1, jsv2, jsd
     real(DOUBLE_PRECISION) ::  rcut1
     real(DOUBLE_PRECISION) ::  dvecx, dvecy, dvecz, w
     integer                :: log_unit
!
     log_unit=config%calc%distributed%log_unit
!
     if (i_verbose >= 1) then
       if (log_unit > 0) write(log_unit, *) '@@@ get_length_of_list (non-order-N): cutoff= ', cutoff_radius
     endif  
!
     i_show=5
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Check the parapmeters
!
     if (cutoff_radius < 0.0d0) then 
       write(6,*)'ERROR:get_length_of_list'
       stop
     endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     allocate (number_of_members(noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error :ierr=',ierr
        stop
     endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     if( cutoff_radius .lt. huge(1d0) ) then
        rcut1=cutoff_radius*1.001d0
     else
        rcut1=cutoff_radius
     end if
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!$omp  parallel & 
!$omp& default(shared) &
!$omp& private(dvecx,dvecy,dvecz,w) & 
!$omp& private(jsv1,jsv2,jsd)
!$omp  do schedule(static)
     do jsv2=1,noav
!
       jsd=0
!   
       do jsv1=1,noav
         dvecx=atm_position(1,jsv2)-atm_position(1,jsv1)
         dvecy=atm_position(2,jsv2)-atm_position(2,jsv1)
         dvecz=atm_position(3,jsv2)-atm_position(3,jsv1)
!
         if (i_pbc_x == 1) dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
         if (i_pbc_y == 1) dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
         if (i_pbc_z == 1) dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0
!
         dvecx=dvecx*ax
         dvecy=dvecy*ay
         dvecz=dvecz*az
         w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
         if (w .lt. rcut1) then
           jsd=jsd+1
         endif
       enddo  
       number_of_members(jsv2)=jsd
!
     enddo
!$omp end do
!$omp end parallel 
!     
     length_of_list=maxval(number_of_members)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Plot the result (optional)
!
     if (i_verbose >= 1) then
       if (log_unit > 0) then 
         do jsv2=1,noav
           if (jsv2 .le. i_show) then
             write(log_unit, '(a,2i10)') 'jsv,njsd=',jsv2,number_of_members(jsv2)
           endif  
         enddo
       endif  
     endif  
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Final deallocation
!
     deallocate (number_of_members,stat=ierr)
     if( ierr .ne. 0) then
        write(*,*)'alloc. error :ierr=',ierr
        stop
     endif
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   end subroutine get_length_of_list
!
end module M_qm_domain
