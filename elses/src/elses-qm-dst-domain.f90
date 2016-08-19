!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_domain_dst
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_wall_clock_time,  only : get_system_clock_time !(routine)
   implicit none
   logical :: global_dens_mat
   logical :: global_ham_mat 
   logical :: overlap_is_on  
!
   private
!
! Public module variable
   public  :: global_dens_mat
   public  :: global_ham_mat 
   public  :: overlap_is_on  
!
! Public routine
   public :: qm_domain_setting_dst
   public :: proj_init_end_dst
   public :: set_projection_dst_cell
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_domain_setting_dst
!
!   module variables unchanged
!      --> 
!   module variables changed
!      --> 
!
     use M_config,           only : config !(unchanged) 
     use M_qm_domain,        only : noav  !(CHANGED)
     use M_qm_domain,        only : atm_position, atm_element !(CHANGED)
     use elses_mod_sim_cell, only : noa   ! (unchanged)
     use elses_mod_tx,       only : tx, ty, tz, jsei
     use M_qm_domain,        only : make_booking_list !(routine)
     use M_md_dst_cell,      only : initial_mcell_booking !(routine)
     use M_md_dst_update_cell_bk, only : update_cell_booking !(routine)
     use M_md_dst_update_cell_bk, only : check_cell_booking !(routine)
     use M_qm_dstm_ini_setting,   only : qm_domain_ini_setting_dstm !(routine)
     use M_ini_load_vatom,        only : set_atm_position_from_tx   !(routine)
!
     implicit none
     integer :: jsv,js, nss, ierr
     logical :: calledFirst
     real(8) :: time_wrk, time_wrk_previous
!
     noav=noa
     if (i_verbose >= 1) then
       write(*,*)'@@ qm_domain_setting_dst'
       write(*,*)'   global_dens_mat =',global_dens_mat
       write(*,*)'   global_ham_mat  =',global_ham_mat
       write(*,*)'   overlap_is_on   =',overlap_is_on
     endif   
!
     if (i_verbose >= 1) then
       write(*,*)'  Set noav=',noav
     endif   
!
     if (config%system%structure%use_vatom) then
       call set_atm_position_from_tx 
!        ------> Set atm_position, atm_element
!                 from tx, ty, tz, jsei
     endif   
!
     call get_system_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
     call initial_mcell_booking
!        ------> Set n_cell, n_cell_x, n_cell_y, n_cell_z
!                   cell_for_atom(:), num_of_atoms_in_cell(:), atom_list_in_cell(:,:)
!
     call check_cell_booking
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a',time_wrk-time_wrk_previous
     if (log_unit > 0) then 
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a',time_wrk-time_wrk_previous
     endif
     time_wrk_previous=time_wrk
!
     call qm_domain_ini_setting_dstm
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:b',time_wrk-time_wrk_previous
     if (log_unit > 0) then 
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:b',time_wrk-time_wrk_previous
     endif
     time_wrk_previous=time_wrk
!
!    stop 'stop manually'
!
   end subroutine qm_domain_setting_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Get the distance for an atom pair
!
   subroutine get_distance(jsv1, jsv2, dist)
!
     use M_qm_domain, only : noav, ax, ay, az             !(unchanged)
     use M_qm_domain, only : i_pbc_x, i_pbc_y, i_pbc_z    !(unchanged)
     use M_qm_domain, only : atm_position                 !(unchanged)
!
     implicit none
     integer, intent(in)                  :: jsv1, jsv2
     real(DOUBLE_PRECISION), intent(out)  :: dist
!
     real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz
!
     if (jsv1 == jsv2) then
      dist=0.0d0
      return
     endif
!   
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
     dist=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
   end subroutine get_distance
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_init_end_dst(imode)
!
    use elses_arr_kry_glo,  only : rcut_kry,noak_str,jsv4jsk_str !(CHANGED)
    use M_qm_domain,        only : noav, ax, ay, az              !(unchanged)
    use elses_param_ctl_kr, only : noak_min_def                  !(unchanged)
    use elses_mod_sel_sys,  only : r_cut_book                    !(unchanged)
    use M_wall_clock_time,  only : get_system_clock_time         !(routine)
    implicit none
!   integer :: imode, ierr
    integer :: imode
    real(8) :: time_wrk, time_wrk_previous ! work variable for measuring the time
!   real(8) :: proj_radius
!   integer :: atm_index
!   integer :: size1
!   integer :: num_atom_list
!
    if (i_verbose >=1) write(*,*)'@@ qm_solver_projection'
!
    call get_system_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
    call elses_set_param_ctl_kr
!
    call get_system_clock_time(time_wrk)
    if (log_unit > 0) write(log_unit,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(1)  =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
    if (imode == 1) then
!TEST call elses_alloc_kry_glo(imode)
!
      call get_system_clock_time(time_wrk)
      if (log_unit > 0) write(log_unit,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(2)  =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
!     if (.not. allocated(jsv4jsk_str)) then
!       size1=noak_min_def*2
!       allocate(jsv4jsk_str(size1, noav),stat=ierr)
!       if (ierr /= 0) stop 'ERROR(proj_init_end_dst):jsv4jsk_str'
!     endif  
!
!     do atm_index=1,noav
!      proj_radius=r_cut_book
!      num_atom_list=noak_min_def
!      call set_projection_dst_cell(atm_index, info, num_atom_list, proj_radius, jsv4jsk_str(:, atm_index))
!      if (info /= 0) then
!        write(*,*)'ERROR(proj_init_end_dst):atm_index, info=',atm_index, info
!        stop
!      endif
!      rcut_kry(atm_index)=proj_radius
!      noak_str(atm_index)=num_atom_list
!     enddo
!
!     call get_system_clock_time(time_wrk)
!     write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(3)  =',time_wrk-time_wrk_previous
!     time_wrk_previous=time_wrk
!
!     do atm_index=1,noav
!       write(*,*)'atm_index, noak_str, rcut_kry=',atm_index, noak_str(atm_index), rcut_kry(atm_index)
!     enddo
!
!     write(*,'(a,3i10,f10.5)')'noav_kr,noak_max,min,ave=', &
!&                             noav, maxval(noak_str), minval(noak_str), dble(sum(noak_str))/dble(noav)
!
!
!     call elses_alloc_jsv4jsk_str(imode)
!
      call get_system_clock_time(time_wrk)
      if (log_unit > 0)  write(log_unit,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(4)  =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
!     call get_system_clock_time(time_wrk)
!     write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(5)  =',time_wrk-time_wrk_previous
!     time_wrk_previous=time_wrk
!
    endif  
!
    if (imode == 2) then
      call elses_alloc_kry_glo(imode)
!     call elses_alloc_jsv4jsk_str(imode)
    endif
!
!   stop 'stop manually'
!  
  end subroutine proj_init_end_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine set_projection_dst_cell(atm_index, info, ntry, num_atom_list, proj_radius, atom_list)
!
!   use elses_arr_kry_glo,  only : rcut_kry,noak_min,noak_str  !(CHANGED!)
!
    use M_config,           only : config !(unchanged) 
!     only:  config%calc%calc_check%set 
!            config%calc%calc_check%short_atom_pair_distance%set 
!            config%calc%calc_check%short_atom_pair_distance%warning_level                                                      
!            config%calc%calc_check%short_atom_pair_distance%abort_level  
!
    use M_qm_domain,   only : noav, ax, ay, az                 !(unchanged)
    use M_qm_domain,   only : i_pbc_x, i_pbc_y, i_pbc_z        !(unchanged)
    use M_md_dst_cell, only : one_cell_setting                 !(unchanged)
    use M_md_dst_cell, only :  n_cell_x,  n_cell_y,  n_cell_z  !(unchanged)
    use M_md_dst_cell, only : dn_cell_x, dn_cell_y, dn_cell_z  !(unchanged)
!   use elses_mod_sel_sys,    only : r_cut_book                !(unchanged)
!   use elses_param_ctl_kr,   only : noak_min_def              !(unchanged)
    use M_qm_proj_list_dstm, only : set_projection_list_dstm   !(routine)
!
    implicit none
    logical, parameter :: debug_mode = .false.
!
    integer, intent(in)     :: atm_index
    integer, intent(out)    :: info
    integer, intent(out)    :: ntry
    integer, intent(inout)  :: num_atom_list 
!                               ! Threashhold value  (input)
!                               ! Resultant value (output)
    real(8), intent(inout)  :: proj_radius
    integer, intent(out)    :: atom_list(:)
!
    integer :: jdx, jdy, jdz
    integer :: imode, ierr
    integer, allocatable :: atm_list_wrk(:)
!   integer, allocatable :: atm_list_wrk2(:)
!   integer :: length_of_list, size1, ntrymax
    integer :: length_of_list, ntrymax
    real(8) :: max_distance
!   real(8) :: rcut_book_wrk, rcut1, dist
    real(8) :: rcut_book_wrk, dist
!   real(8) :: rcut_book_tmp
    integer :: atm_index2
!   integer :: itry, n_count, j_atom
    integer :: n_count, j_atom
    integer :: num_atom_list_min
!   integer :: iexit
    logical, allocatable :: booked(:)
    real(8), allocatable :: atm_distance_list_wrk(:) ! List for atom-pair distance for the given atom seed
!
    real(8) :: nna_distance        ! nearest neighbor atomic distance          (for checking)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameter setting
!
    num_atom_list_min=num_atom_list
    rcut_book_wrk=proj_radius
    ntrymax=1000
    info=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting for the cell-booking method
!
    if (.not. one_cell_setting) then
      jdx=dn_cell_x                      !! TEMPORAL SETTING!
      jdy=dn_cell_y                      !! TEMPORAL SETTING!
      jdz=dn_cell_z                      !! TEMPORAL SETTING!
!     size1=noav                         !! TEMPORAL SETTING! 
!
      if ((jdx < 1) .or. (jdx*2+1 > n_cell_x)) then
        write(*,*)'jdx*2 , n_cell_x=',jdx*2 , n_cell_x
        stop
      endif
!   
      if ((jdy < 1) .or. (jdy*2+1 > n_cell_y)) then
        write(*,*)'jdy*2 , n_cell_y=',jdy*2 , n_cell_y
        stop
      endif
!
      if ((jdz < 1) .or. (jdz*2+1 > n_cell_z)) then
        write(*,*)'jdz*2 , n_cell_z=',jdz*2 , n_cell_z
        stop
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write(*,*)'@@ set_projection_dst_cell'
!   write(*,*)'  rcut1=',rcut1
!   write(*,*)'  jdx, jdy, jdz =',jdx, jdy, jdz
!   write(*,*)'  num_atom_list_min=',num_atom_list_min
!   write(*,*)'  rcut_book_wrk    =',rcut_book_wrk
!
    if (debug_mode) then 
      allocate (booked(noav),stat=ierr)
      if (ierr /= 0) stop 'Alloc. error(set_rcut_kry_dst_cell):booked'
    endif
!
    if (rcut_book_wrk > 0.5d0*min(ax,ay,az)) then
      if ( i_pbc_x + i_pbc_y + i_pbc_z /= 0 ) then
        write(*,*)'ERROR:Too large r_cut_book: =',rcut_book_wrk
        write(*,*)'ERROR: (rcut_book_wrk > 0.5d0*min(ax,ay,az))'
        write(*,*)'  i_pbcx,y,z=', i_pbc_x, i_pbc_y, i_pbc_z
        stop
      else
        if (log_unit > 0) then
          write(log_unit,*)'WARNING:Large r_cut_book ? : =',rcut_book_wrk
          write(log_unit,*)'WARANIG: (rcut_book_wrk > 0.5d0*min(ax,ay,az))'
          write(log_unit,*)'  i_pbcx,y,z=', i_pbc_x, i_pbc_y, i_pbc_z
        endif
      endif   
    endif
!
    if ((num_atom_list_min < 1) .or. (num_atom_list_min > noav)) then
      write(*,*)'ERROR; num_atom_list_min =', num_atom_list_min
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set length_of_list, max_distance
!
    if (one_cell_setting) then
      max_distance=0.5d0*min(ax, ay, az)
      if ( i_pbc_x + i_pbc_y + i_pbc_z == 0 ) max_distance=huge(1.0d0)
      length_of_list=noav
    else
      imode=1
      call set_atom_list_from_cell(imode, atm_index, jdx, jdy, jdz, & 
&                                  length_of_list, max_distance)
!      ----> Set length_of_list : number of atoms in the neighboring M cells
!                                     with M = (jdx+1)*(jdy+1)*(jdz+1) 
!                 max_distance  : projection radius should be less than max_distance
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      if (atm_index < 5) then
       if (log_unit > 0) then
         write(log_unit,*)' length of prebooking list=', atm_index, length_of_list, num_atom_list_min
       endif
      endif
    endif
!
    if (length_of_list < num_atom_list_min) then
      write(*,*)'ERROR(set_projection_dst_cell): atm_index= ', atm_index
      write(*,*)'  length_of_list    =', length_of_list 
      write(*,*)'  num_atom_list_min =', num_atom_list_min
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set atm_list_wrk
!
    allocate (atm_list_wrk(length_of_list),stat=ierr)
    if (ierr /= 0) stop 'Alloc. error(set_projection_dst_cell):atm_list_wrk'
!
    if (one_cell_setting) then
      if (length_of_list /= noav) then
        write(*,*)'ERROR(set_projection_dst_cell):length_of_list, noav=',length_of_list, noav
        stop
      endif
      atm_list_wrk(1)=atm_index
      n_count=1
      do j_atom=1,noav
        if (j_atom /= atm_index) then
          n_count=n_count+1
          atm_list_wrk(n_count)=j_atom
        endif
      enddo
    else
      call set_atom_list_from_cell(imode, atm_index, jdx, jdy, jdz, & 
&                                  length_of_list, max_distance, atm_list_wrk)
    endif
!      ----> Set atm_list_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set atm_distance_list_wrk
!
    allocate (atm_distance_list_wrk(length_of_list),stat=ierr)
    if (ierr /= 0) stop 'Alloc. error(set_projection_dst_cell):atm_distance_list_wrk'
!
    do j_atom=1,length_of_list
      atm_index2=atm_list_wrk(j_atom)
      if ((atm_index2 < 1) .or. (atm_index2 > noav)) then
        write(*,*)'ERROR(set_rcut_kry_dst_cell):j_atom,atm_index2=',j_atom,atm_index2
        stop
      endif   
      call get_distance(atm_index, atm_index2, dist)
      atm_distance_list_wrk(j_atom)=dist
    enddo
!      ----> Set atm_distance_list_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the nearest nighbor atom (optional)
!
    if ((config%calc%calc_check%set) .and. (config%calc%calc_check%short_atom_pair_distance%set)) then
      nna_distance=minval(atm_distance_list_wrk(2:length_of_list))
      j_atom      =minloc(atm_distance_list_wrk(2:length_of_list),1)
!     write(*,'(a,2i10, f20.10)') 'INFO :Short atom-pair distance:', & 
!&                                       atm_index, atm_list_wrk(j_atom), nna_distance
      if (nna_distance < config%calc%calc_check%short_atom_pair_distance%abort_level  ) then
        write(*,'(a,2i10, f20.10)') 'ABORT   :Too short atom-pair distance [au]:', & 
&                                        atm_index, atm_list_wrk(j_atom), nna_distance
        stop
      endif  
      if (nna_distance < config%calc%calc_check%short_atom_pair_distance%warning_level) then
        write(*,'(a,2i10, f20.10)') 'WARNING :Too short atom-pair distance [au]:', & 
&                                        atm_index, atm_list_wrk(j_atom), nna_distance
      endif  
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the projection list
!
    ntry=ntrymax                    !(input value, will be changed in the routine)
    proj_radius=rcut_book_wrk       !(input value, will be changed in the routine)
    num_atom_list=num_atom_list_min !(input value, will be changed in the routine)
    call set_projection_list_dstm(atm_index, num_atom_list, info, ntry, proj_radius, &
&                              max_distance, atom_list, atm_list_wrk, atm_distance_list_wrk)
!                           ---> Set proj_radius, num_of_atm_list, atom_list  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    deallocate (atm_list_wrk, stat=ierr)
    if (ierr /= 0) stop 'Dealloc. error(set_projection_dst_cell):atm_list_wrk'
!
    deallocate (atm_distance_list_wrk, stat=ierr)
    if (ierr /= 0) stop 'Dealloc. error(set_projection_dst_cell):atm_distance_list_wrk'

  end subroutine set_projection_dst_cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine set_atom_list_from_cell(imode, atm_index, jdx, jdy, jdz, length_of_list, max_distance, atom_list) 
!
    use M_qm_domain,   only : noav, ax, ay, az          !(unchanged)
    use M_qm_domain,   only : i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
    use M_qm_domain,   only : atm_position              !(unchanged)
    use M_md_dst_cell, only : n_cell, n_cell_x, n_cell_y, n_cell_z !(unchanged)
    use M_md_dst_cell, only : cell_for_atom        !(unchanged)
    use M_md_dst_cell, only : num_of_atoms_in_cell !(unchanged)
    use M_md_dst_cell, only : atom_list_in_cell    !(unchanged)
!
    implicit none
    integer, intent(in)    :: imode, atm_index, jdx, jdy, jdz
    integer, intent(out)   :: length_of_list
    real(8), intent(out)   :: max_distance
    integer, optional, intent(out)   :: atom_list(:)
!
    integer :: j_cell, j_cell_x, j_cell_y, j_cell_z
    integer :: n_count, jx, jy, jz
    integer :: k, kx, ky, kz, k_atom
    integer :: atm_index2, size1
    real(8) :: rx, ry, rz
    real(8) :: dx1, dx2, dy1, dy2, dz1, dz2, dx, dy, dz
!    
    if (imode == -1) then
      write(*,*)'STOP;set_atom_list_from_cell:imode=',imode
      stop
    endif
!
    if ((atm_index < 1) .or. (atm_index > noav)) then
      write(*,*)'ERROR(set_atom_list_from_cell):atm_index=',atm_index
      stop
    endif
!
    j_cell=cell_for_atom(atm_index)
    if ((j_cell < 0) .or. (j_cell > n_cell-1)) then
      write(*,*)'ERROR(set_atom_list_from_cell):atm_index,j_cell=',atm_index,j_cell
      stop
    endif
    j_cell_z = j_cell/(n_cell_x*n_cell_y)
    j_cell_y = (j_cell-j_cell_z*n_cell_x*n_cell_y)/n_cell_x
    j_cell_x = j_cell - j_cell_z*n_cell_x*n_cell_y - j_cell_y*n_cell_x
!
    rx = atm_position(1,atm_index)*ax
    ry = atm_position(2,atm_index)*ay
    rz = atm_position(3,atm_index)*az
!
    dx1=dble(j_cell_x+1+jdx)*(ax/dble(n_cell_x))
    dx2=dble(j_cell_x  -jdx)*(ax/dble(n_cell_x))
    dy1=dble(j_cell_y+1+jdy)*(ay/dble(n_cell_y))
    dy2=dble(j_cell_y  -jdy)*(ay/dble(n_cell_y))
    dz1=dble(j_cell_z+1+jdz)*(az/dble(n_cell_z))
    dz2=dble(j_cell_z  -jdz)*(az/dble(n_cell_z))
!
    dx=min(dabs(dx1-rx),dabs(dx2-rx))
    dy=min(dabs(dy1-ry),dabs(dy2-ry))
    dz=min(dabs(dz1-rz),dabs(dz2-rz))
!
    max_distance=min(dx,dy,dz)
!
!   write(*,*)'atm_index   =', atm_index
!   write(*,*)'max_distance=', max_distance
!   write(*,*)'     j_cell =',j_cell
!   write(*,*)' j_cell_x, j_cell_y, j_cell_z =', j_cell_x, j_cell_y, j_cell_z
!   write(*,*)' jdx, jdy, jdz =', jdx, jdy, jdz
!   write(*,*)' ax/n_cell_x =', ax/dble(n_cell_x)
!   write(*,*)' ay/n_cell_y =', ay/dble(n_cell_y)
!   write(*,*)' az/n_cell_z =', az/dble(n_cell_z)
!   write(*,*)'        tx  =', atm_position(1,atm_index)
!   write(*,*)'        ty  =', atm_position(2,atm_index)
!   write(*,*)'        tz  =', atm_position(3,atm_index)
!   write(*,*)'        rx  =', rx
!   write(*,*)'        ry  =', ry
!   write(*,*)'        rz  =', rz
!   write(*,*)'       dx1  =', dx1
!   write(*,*)'       dx2  =', dx2
!   write(*,*)'       dy1  =', dy1
!   write(*,*)'       dy2  =', dy2
!   write(*,*)'       dz1  =', dz1
!   write(*,*)'       dz2  =', dz2
!   write(*,*)'        dx  =', dx
!   write(*,*)'        dy  =', dy
!   write(*,*)'        dz  =', dz
!
    if (present(atom_list)) then
      atom_list(:)=0
      size1=size(atom_list,1)
      atom_list(1)=atm_index
    endif
!
    n_count=1
!    
    do jz=-jdz,jdz
      kz=jz+j_cell_z
      if (kz < 0)          kz = kz + n_cell_z
      if (kz > n_cell_z-1) kz = kz - n_cell_z
      do jy=-jdy,jdy
        ky=jy+j_cell_y
        if (ky < 0)          ky = ky + n_cell_y
        if (ky > n_cell_y-1) ky = ky - n_cell_y
        do jx=-jdx,jdx
          kx=jx+j_cell_x
          if (kx < 0)          kx = kx + n_cell_x
          if (kx > n_cell_x-1) kx = kx - n_cell_x
          k=kx+ky*n_cell_x+kz*n_cell_x*n_cell_y
          do k_atom=1,num_of_atoms_in_cell(k)
            atm_index2=atom_list_in_cell(k_atom,k)
            if (atm_index /= atm_index2) then
              n_count=n_count+1
              if (present(atom_list)) then
                if (n_count > size1) then
                  write(*,*)'ERROR(set_atom_list_from_cell)'
                  write(*,*)'atm_index, size1, n_count=',atm_index, size1, n_count
                  write(*,*)' jx, jy, jz=',jx, jy, jz
                  write(*,*)' kx, ky, kz=',kx, ky, kz
                  stop
                endif
                atom_list(n_count)=atm_index2
              endif
            endif
          enddo
        enddo
      enddo
    enddo  
!
    length_of_list=n_count
!
  end subroutine set_atom_list_from_cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_domain_dst

