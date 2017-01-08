!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_proj_cell
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION      !(unchanged)
  use M_io_dst_write_log, only : lu => log_unit                  !(unchanged)
  use M_wall_clock_time, only : get_system_clock_time      !(routine)
  use M_md_dst,          only : set_dst_final              !(routines)
  implicit none
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  DST arrays for the member list  (DIFFERENT VALUES AMONG THE NODES)
!
  integer, allocatable :: dst_atm_list(:)     !  set in 'set_dst_atm_list'
!                ! the member list of the distributed atoms on each node
  integer, allocatable :: len_dst_atm_list(:) !  set in 'set_dst_atm_list'
!                ! len_dst_atm_list(1)        : length of the 'core' member list
!                ! len_dst_atm_list(n) (n>1)  : not used now
  integer, allocatable :: dst_atm_list_rev(:) !  set in 'set_dst_atm_list'
!                ! reverse list for the distributed atoms
!
! |
! |  ex. do dst_atm_index=1,len_dst_atm_list(n)   (n is one for two)
! |        atm_index=dst_atm_list(dst_atm_index)
! |      enddo
! |
! |  ex. do atm_index=1,noav
! |        dst_atm_index=dst_atm_list_rev(atm_index)
! |      enddo
! |
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  DST arrays for the booking list  (DIFFERENT VALUES AMONG THE NODES)
!
  integer, allocatable :: jsv4jsk_dst(:,:)     ! set in 'set_dst_atm_list'
  integer, allocatable :: noak_dst(:)          ! set in 'set_dst_atm_list'
  real(8), allocatable :: rcut_kry_dst(:)      ! set in 'set_dst_atm_list'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  private
!
! PUBLIC ARRAYS
  public :: dst_atm_list, len_dst_atm_list, dst_atm_list_rev
  public :: jsv4jsk_dst, noak_dst, rcut_kry_dst
!
! PUBLIC ROUTINES
  public :: set_dst_atm_list
  public :: proj_get_mat_size_dst
  public :: proj_get_list_dst
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_dst_atm_list
!
   use M_config,           only : config  !(unchanged)
!                            only used for config%calc%solver%projection_list_length
   use M_md_dst,           only : myrank, nprocs     !(unchanged)
   use M_qm_domain,        only : noav                !(unchanged)
   use elses_param_ctl_kr, only : noak_min_def        !(unchanged)
!  use elses_mod_sel_sys,  only : r_cut_book          !(unchanged)
!
   use M_qm_domain_dst,    only : set_projection_dst_cell !(routine)
   use M_md_dst_get_atom_range, only : dst_get_atm_index_range !(routine)
   use M_qm_proj_suggest_radius, only : suggest_proj_radius    !(routine)
   use M_dstm_reordering,        only : reorder_booked_atoms   !(routine)
!
   implicit none
   logical, parameter :: flag_for_reordering = .false.
   integer :: ierr
   integer :: atm_index_ini, atm_index_fin, atm_index, jj
   integer :: jsk, jsv, dst_atm_index
   logical, allocatable :: booked_for_core(:)
   integer, allocatable :: booked_for_projection(:)
!
   integer :: size1, size2
   integer :: info, num_atom_list
   real(8) :: proj_radius
!  integer :: dst_atm_index_max
   logical, parameter :: second_booking = .false.
   integer :: ntry
   real(DOUBLE_PRECISION) :: time_wrk, time_wrk_previous
   logical :: flag_for_allocate
!
   real(DOUBLE_PRECISION) :: proj_radius_suggested
   character(len=32) :: suggestion_mode
!
   if (i_verbose >= 0) then 
     if (lu>0) write(lu,*)'@@ set_dst_atm_list'
   endif  
!
   if ((nprocs < 1) .or. (nprocs > 1000000)) then
     write(*,*)'ERROR:nprocs=',nprocs
     stop
   endif
!   
   if ((myrank < 0) .or. (myrank > nprocs-1)) then
     write(*,*)'ERROR:myrank=',nprocs
     stop
   endif
!   
   if (.not. allocated(len_dst_atm_list)) then
     allocate(len_dst_atm_list(2),stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. of len_dst_atm_list'
   endif
   len_dst_atm_list(:)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set len_dst_atm_list(:)
!
!  atm_index_ini=myrank*noav/nprocs+1
!  atm_index_fin=(myrank+1)*noav/nprocs
   call dst_get_atm_index_range(atm_index_ini,atm_index_fin)
   if (lu >0) write(lu,'(a,3i10)')'INFO-MPI: procs, ini, fin=',myrank, atm_index_ini, atm_index_fin
!
   jj=0
   do atm_index=atm_index_ini, atm_index_fin
     jj=jj+1
   enddo
   len_dst_atm_list(1)=jj
!
   if (jj /= len_dst_atm_list(1)) then
     write(*,*)'ERROR(set_dst_atm_list):jj=',jj, len_dst_atm_list(1)
     stop
   endif
!
   len_dst_atm_list(2)=len_dst_atm_list(1) ! dummy setting for compatibility
   if (lu >0) write(lu,'(a,i10)')'INFO-MPI: len_dst_atm_list(1)=',len_dst_atm_list(1)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate dst_atm_list, rcut_kry_dst, noak_dst
!
   flag_for_allocate = .false.
   if (.not. allocated(dst_atm_list)) then 
     flag_for_allocate = .true. 
   else
     if (size(dst_atm_list,1) < len_dst_atm_list(2)) flag_for_allocate = .true.
   endif
!
   if (flag_for_allocate) then
     size1=int(dble(len_dst_atm_list(2))*1.2d0)
!
     if (allocated(dst_atm_list)) then 
       deallocate(dst_atm_list,stat=ierr)
       if (ierr /= 0) stop 'ERROR in dealloc. of dst_atm_list'
     endif
     allocate(dst_atm_list(size1),stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. of dst_atm_list'
!
     if (allocated(rcut_kry_dst)) then
       deallocate(rcut_kry_dst,stat=ierr)
       if (ierr /= 0) stop 'ERROR in dealloc. of rcut_kry_dst'
     endif
     allocate(rcut_kry_dst(size1),stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. of rcut_kry_dst'
!
     if (allocated(noak_dst)) then
       deallocate(noak_dst,stat=ierr)
       if (ierr /= 0) stop 'ERROR in dealloc. of noak_dst'
     endif
     allocate(noak_dst(size1),stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. of noak_dst'
!
   endif
!
   dst_atm_list(:)=0
   rcut_kry_dst(:)=0
   noak_dst(:)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set dst_atm_list(:)
!
   jj=0
   do atm_index=atm_index_ini, atm_index_fin
     jj=jj+1
     dst_atm_list(jj)=atm_index
   enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate jsv4jsk_dst(:,:)
!
   flag_for_allocate = .false.
   if (.not. allocated(jsv4jsk_dst)) then 
     flag_for_allocate = .true. 
   else
     if (size(jsv4jsk_dst,2) < len_dst_atm_list(2)) flag_for_allocate = .true.
   endif
!
   size1=min(noav, noak_min_def*2)
   size2=int(dble(len_dst_atm_list(2))*1.2d0)
!
   select case(config%calc%solver%projection_list_length)
     case (-1,-2)
       size1=min(noav, noak_min_def*2) ! default setting
       if (lu >0) write(lu,'(a,i10)')'INFO-DST: projection_list_length (default)   =', size1
     case (-3)
       size1=noav                      ! 'full'  setting
       if (lu >0) write(lu,'(a,i10)')'INFO-DST: projection_list_length (full)      =', size1
     case default
       size1=config%calc%solver%projection_list_length
       if ((size1 < noak_min_def) .or. (size1 > noav)) then
         write(*,*)'ERROR:projection_list_length=',size1
         stop
       endif
       if (lu >0) write(lu,'(a,i10)')'INFO-DST: projection_list_length (specified) =', size1
   end select
!

!
   if (flag_for_allocate) then
     if (allocated(jsv4jsk_dst)) then 
       deallocate(jsv4jsk_dst,stat=ierr)
       if (ierr /= 0) stop 'ERROR in dealloc. of jsv4jsk_dst'
     endif
     allocate(jsv4jsk_dst(size1, size2),stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. of noak_dst'
     jsv4jsk_dst(:,:)=0
   endif
!
!
   call suggest_proj_radius(noak_min_def, suggestion_mode, proj_radius_suggested)
!
!  if (i_verbose >= 0) then 
!    write(*,*) 'INFO:suggested projecion radius:      mode =',trim(suggestion_mode)
!    write(*,*) 'INFO:suggested projecion radius: value(au) =', proj_radius_suggested 
!  endif  
!
   if (lu >0) then
     write(lu,*) 'INFO:suggested projecion radius:      mode =',trim(suggestion_mode)
     write(lu,*) 'INFO:suggested projecion radius: value(au) =', proj_radius_suggested 
   endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting jsv4jsk_dst(:,:), rcut_kry_dst(:), noak_dst(:)
!
   call get_system_clock_time(time_wrk)
   time_wrk_previous=time_wrk
!
!$omp  parallel default(shared) &
!$omp& private (dst_atm_index, atm_index) &
!$omp& private (proj_radius, num_atom_list, info, ntry)
!$omp  do schedule(static)
   do dst_atm_index=1,len_dst_atm_list(1)
     atm_index=dst_atm_list(dst_atm_index)
     proj_radius=proj_radius_suggested
!    proj_radius=r_cut_book
     num_atom_list=noak_min_def
     call set_projection_dst_cell(atm_index, info, ntry, num_atom_list, proj_radius, jsv4jsk_dst(:, dst_atm_index))
     if (info /= 0) then
       write(*,*)'ERROR(proj_init_end_dst):atm_index, info=',atm_index, info
       stop
     endif
     if (i_verbose >= 1) then 
       if (atm_index < 10) then 
         if (lu>0) write(lu,'(a,3i10)')'  atm_index, ntry, num_atom_list =',atm_index, ntry, num_atom_list
       endif
     endif  
     rcut_kry_dst(dst_atm_index)=proj_radius
     noak_dst(dst_atm_index)=num_atom_list
   enddo
!$omp end do
!$omp end parallel
!
   call get_system_clock_time(time_wrk)
!  write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:dst_atm_list(main) =',time_wrk-time_wrk_previous
   if (lu > 0) write(lu,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:dst_atm_list(main) =', & 
&                         time_wrk-time_wrk_previous
   time_wrk_previous=time_wrk
!
   if (flag_for_reordering) then
!
!$omp  parallel default(shared) &
!$omp& private (dst_atm_index, atm_index) &
!$omp& private (num_atom_list)
!$omp  do schedule(static)
   do dst_atm_index=1,len_dst_atm_list(1)
     num_atom_list=noak_dst(dst_atm_index)
     atm_index=dst_atm_list(dst_atm_index)
     call reorder_booked_atoms( atm_index, jsv4jsk_dst(1:num_atom_list,dst_atm_index) )
   enddo
!$omp end do
!$omp end parallel
!
     if (lu > 0) write(lu,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:reordering       =', & 
&                         time_wrk-time_wrk_previous
     time_wrk_previous=time_wrk
!
   endif

   if (.not. second_booking) then 
     if (lu >0)  write(lu,*)'.... scond booking is skipped'
!    stop 'Stop manually'
     return
   endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preparation for the second shell booking
!
   if (second_booking) then
     if (.not. allocated(booked_for_core)) then
       allocate(booked_for_core(noav),stat=ierr)
       if (ierr /= 0) stop 'ERROR in alloc. of booked_for_core'
     endif
     booked_for_core(:)= .false.
!   
     jj=0
     do atm_index=atm_index_ini, atm_index_fin
       jj=jj+1
       dst_atm_list(jj)=atm_index
       booked_for_core(atm_index)= .true.
     enddo
!
     if (.not. allocated(booked_for_projection)) then
       allocate(booked_for_projection(noav),stat=ierr)
       if (ierr /= 0) stop 'ERROR in alloc. of booked_for_projection'
     endif
     booked_for_projection(:)= 0
!
     if (.not. allocated(dst_atm_list_rev)) then
       allocate(dst_atm_list_rev(noav),stat=ierr)
       if (ierr /= 0) stop 'ERROR in alloc. of dst_atm_list_rev'
     endif
     dst_atm_list_rev(:)=0
   endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second shell booking
!
   do dst_atm_index=1,len_dst_atm_list(1)
     atm_index=dst_atm_list(dst_atm_index)
     do jsk=1,noak_dst(dst_atm_index)
       jsv=jsv4jsk_dst(jsk,dst_atm_index)
       if (.not. booked_for_core(jsv)) then
         booked_for_projection(jsv)=booked_for_projection(jsv)+1
       endif
     enddo   
   enddo  
!
   jj=len_dst_atm_list(1)
   do atm_index=1,noav
     if (.not. booked_for_core(atm_index)) then
       if ( booked_for_projection(atm_index) >0   ) then
         jj=jj+1
!        write(*,'(a,3i10)')'myrank, jj, atm_index=',myrank, jj, atm_index
         if (jj > size(dst_atm_list,1)) then
           write(*,*)'error in booking dst_atm_list : jj=',jj
           stop
         endif   
         dst_atm_list(jj)=atm_index
       endif   
     endif
   enddo   
   len_dst_atm_list(2)=jj
   if (lu >0) write(lu,'(a,i10)')'INFO-MPI: len_dst_atm_list(2)=',len_dst_atm_list(2)
!
   dst_atm_list_rev(:)=0
   do dst_atm_index=1,len_dst_atm_list(2)
     atm_index=dst_atm_list(dst_atm_index)
     if ((atm_index > noav) .or. (atm_index < 1)) then
       write(*,*)'ERROR(set_dst_atm_list):dst_atm_list=',dst_atm_index, atm_index
       stop
     endif
     if (dst_atm_list_rev(atm_index) /= 0) then
       write(*,*)'ERROR(set_dst_atm_list):dst_atm_list=',dst_atm_index, atm_index
       write(*,*)'dst_atm_list_rev(atm_index) =',dst_atm_list_rev(atm_index)
       stop
     endif
     dst_atm_list_rev(atm_index)=dst_atm_index
   enddo
!
   deallocate(booked_for_core,stat=ierr)
   if (ierr /= 0) stop 'ERROR in dealloc. of booked_for_core'
!
   deallocate(booked_for_projection,stat=ierr)
   if (ierr /= 0) stop 'ERROR in dealloc. of booked_for_projection'
!
   if (i_verbose >= 0) then 
     if (lu>0) write(lu,*)' ... ended (set_dst_atm_lst)'
   endif  
!
  end subroutine set_dst_atm_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_get_mat_size_dst(dst_atm_index,mat_size,num_atom_proj)
!
    use elses_mod_js4jsv,     only : js4jsv                ! (unchanged)
!   use elses_mod_tx,         only : jsei                  ! (unchanged)
    use M_qm_domain,          only : noav, nval            ! (unchanged)
    use M_qm_domain,          only : atm_element           ! (unchanged)
!
    implicit none
    integer, intent(in)  :: dst_atm_index
    integer, intent(out) :: mat_size
    integer, intent(out) :: num_atom_proj
!
    integer :: atm_index
    integer :: jj, jsk, jsv, js, nss, nval2, jsk_self
!
!    
    if ((dst_atm_index <1) .or. (dst_atm_index > len_dst_atm_list(1))) then
      write(*,*)'ERROR(proj_get_mat_size_dst):dst_atm_index=', dst_atm_index
      write(*,*)'ERROR(proj_get_mat_size_dst):len_dst_atm_list(1)=', len_dst_atm_list(1)
      stop
    endif
!
    atm_index=dst_atm_list(dst_atm_index)
!
    if ((atm_index <1) .or. (atm_index > noav)) then
      write(*,*)'ERROR(proj_get_mat_size_dst):atm_index=',atm_index
      stop
    endif
!
    num_atom_proj=noak_dst(dst_atm_index)
!       ---> number of atom for projection
!
    jj=0
    jsk_self=0
    do jsk=1,num_atom_proj
      jsv=jsv4jsk_dst(jsk,dst_atm_index)
      if ((jsv .le. 0) .or. (jsv .gt. noav)) then
        write(6,*)'ERROR!(qm_proj_get_mat_size):jsv,noav=',jsv,noav
        stop
      endif   
      if (jsv == atm_index) then
        if (jsk_self == 0) then
          jsk_self=jsk
        else
          write(6,*)'ERROR!(qm_proj_get_mat_size):jsv,jsk_self=',jsv,jsk_self
          stop
        endif   
      endif   
      js=js4jsv(jsv)
!     nss=jsei(js)
      nss=atm_element(js)
      if ((nss <= 0) .or. (nss > 100)) then
        write(*,*)'ERROR(proj_get_mat_size_dst):nss=',nss
        stop
      endif   
      nval2=nval(nss)
      jj=jj+nval2
    enddo  
    mat_size=jj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Check that The 'self' atom is set to be the first element
!              ( jsv4jsk_str(1,atm_index) = atm_index )
!
    if (jsk_self /= 1) then
      write(*,*)'ERROR!(qm_proj_get_mat_size):jsk_self=',jsk_self
      stop
    endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Trivial checking
!
    if ((mat_size <= 0) .or. (mat_size > num_atom_proj*100)) then
      write(*,*)'ERROR!(qm_proj_get_mat_size)'
      write(*,*)'atm_index, mat_size=',atm_index,mat_size
      stop
    endif   
!
  end subroutine proj_get_mat_size_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_get_list_dst(dst_atm_index,jsv4jsk,jjkset)
!
!   use elses_arr_kry_glo,    only : jsv4jsk_str, noak_str ! (unchanged)
    use elses_mod_js4jsv,     only : js4jsv                ! (unchanged)
!   use elses_mod_tx,         only : jsei                  ! (unchanged)
    use M_qm_domain,          only : noav, nval            ! (unchanged)
    use M_qm_domain,          only : atm_element           ! (unchanged)
!
    implicit none
    integer, intent(in)      :: dst_atm_index
    integer, intent(out)     :: jsv4jsk(:)
    integer, intent(out)     :: jjkset(:)
!   integer, intent(out)     :: jsk4jsv(:)
!
    integer :: atm_index
    integer :: num_atom_proj
    integer :: jsk, jsv, jj, js, nss, nval2
    integer :: nval_upper_limit
!
    nval_upper_limit=100
!
    if ((dst_atm_index <1) .or. (dst_atm_index > len_dst_atm_list(1))) then
      write(*,*)'ERROR(proj_get_mat_size_dst):dst_atm_index=', dst_atm_index
      write(*,*)'ERROR(proj_get_mat_size_dst):len_dst_atm_list(1)=', len_dst_atm_list(1)
      stop
    endif
!
    atm_index=dst_atm_list(dst_atm_index)
!
    if ((atm_index <1) .or. (atm_index > noav)) then
      write(*,*)'ERROR(proj_get_mat_size_dst):atm_index=',atm_index
      stop
    endif
!
    num_atom_proj=size(jsv4jsk,1)
    if (num_atom_proj /= noak_dst(dst_atm_index)) then
      write(*,*)'ERROR(proj_get_list):num_atom_proj=',num_atom_proj
      stop
    endif
!
    jsv4jsk(1:num_atom_proj)=jsv4jsk_dst(1:num_atom_proj,dst_atm_index)
!   jsk4jsv(:)=0
!
!   do jsk=1, num_atom_proj
!     jsv=jsv4jsk(jsk)
!     if ((jsv .le. 0) .or. (jsv .gt. noav)) then
!        write(6,*)'ERROR!(SETKRY):jsk,jsv=',jsk,jsv
!        stop
!     endif   
!     jsk4jsv(jsv)=jsk
!   enddo
!
    jj=0
    do jsk=1, num_atom_proj
      jjkset(jsk)=jj
      jsv=jsv4jsk(jsk)
      if ((jsv .le. 0) .or. (jsv .gt. noav)) then
        write(*,*)'ERROR!(SETKRY:1100):JSK,JSV=',jsk,jsv
        stop
      endif   
      js=js4jsv(jsv)
!     nss=jsei(js)
      nss=atm_element(js)
      if ((nss <= 0) .or. (nss > 100)) then
        write(*,*)'ERROR(proj_get_mat_size_dst):nss=',nss
        stop
      endif   
      nval2=nval(nss)
      if ((nval2 .le. 0) .or. (nval2 .gt. nval_upper_limit)) then
        write(*,*)'ERROR!(SETKRY:1100):JS,NVAL2=',js,nval2
        stop
      endif   
      jj=jj+nval2
     enddo   
!
  end subroutine proj_get_list_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_dst_proj_cell

