!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_proj_list_dstm
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_wall_clock_time,  only : get_system_clock_time !(routine)
   use elses_mod_file_io,  only : vacant_unit  !(function)
   implicit none
   logical :: global_dens_mat
   logical :: global_ham_mat 
   logical :: overlap_is_on  
!
   private
!
! Public routine
   public :: set_projection_list_dstm
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine set_projection_list_dstm(atm_index, length_of_atm_list, info, ntry, proj_radius, & 
&                       max_distance, atm_list, atm_list_wrk, atm_distance_list_wrk)
!
    use M_qm_domain,        only : noav, ax, ay, az !(unchanged)
    use M_qm_domain,        only : i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
!
    implicit none
    logical, parameter :: debug_mode = .false.
!   logical, parameter :: debug_mode = .true.
!
    integer, intent(in)     :: atm_index
    integer, intent(inout)  :: length_of_atm_list
!                            ! input : Minimum number of atom (used as num_atm_list_min)
!                            ! output: Resultant value        
    integer, intent(out)    :: info
    integer, intent(inout)  :: ntry
!                               ! input : Maximum trial (used as ntrymax)
!                               ! output: Resultant value        
    real(8), intent(in)     :: max_distance
!
!   integer, intent(out)    :: length_of_atm_list 
!                               ! output: Resultant value
    real(8), intent(inout)  :: proj_radius
!                               ! input : Unit for radius (used as rcut_book_wrk)
!                               ! output: Resultant value 
    integer, intent(out)    :: atm_list(:)
!                               ! output: Resultant atom list
    integer, intent(in)     :: atm_list_wrk(:)
!                               ! input : Atom list for candidate
    real(8), intent(in)     :: atm_distance_list_wrk(:) ! List for atom-pair distance for the given atom seed

!
    integer :: size1
    integer :: ierr
!   integer, allocatable :: atm_list_wrk(:)
    integer :: length_of_atm_list_wrk
    integer :: ntrymax
    real(8) :: rcut_book_wrk, dist
    real(8) :: rcut_book_tmp
    integer :: atm_index2
    integer :: itry, n_count, j_atom
    integer :: iexit
    integer :: iunit
    integer :: num_atm_list_min
    logical, allocatable :: booked(:)
    logical :: reach_max_length
    integer :: i_pbc_sum
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the local variables as the inputs
!
    i_pbc_sum = i_pbc_x + i_pbc_y + i_pbc_z
!      ----> zero, if ( i_pbc_x, i_pbc_y, i_pbc_z) = ( 0, 0, 0 )
!
    ntrymax=ntry
    length_of_atm_list_wrk=size(atm_list_wrk,1)
    rcut_book_wrk=proj_radius
    num_atm_list_min=length_of_atm_list
!
    if (size(atm_distance_list_wrk,1) /= length_of_atm_list_wrk) then
      write(*,*)'ERROR(set_projection_list_dstm)'
      write(*,*)'size =',length_of_atm_list_wrk, atm_distance_list_wrk
      stop
    endif
!
    if (debug_mode) then
      write(*,*)'set_projection_list_dstm:debug mode'
      size1=noav
      allocate (booked(size1), stat=ierr) 
      if (ierr /= 0) then
        write(*,*)'Alloc error(set_projection_dstm:booked)'
        stop
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set num_atm_list, proj_radius
!
    iexit=0
    reach_max_length=.false.
    do itry=1,ntrymax
      ntry=itry
      if (reach_max_length) then
        write(*,*)'ERROR(set_rcut_kry_dst_cell):reach_max_length=',reach_max_length
        write(*,*)'ERROR(set_rcut_kry_dst_cell):itry=',itry
        write(*,*)'ERROR(set_rcut_kry_dst_cell):rcut_book_tmp=',rcut_book_tmp
        write(*,*)'ERROR(set_rcut_kry_dst_cell):length_of_atm_list_wrk=',length_of_atm_list_wrk
        write(*,*)'ERROR(set_rcut_kry_dst_cell):num_atm_list_min =',num_atm_list_min
        write(*,*)'ERROR(set_rcut_kry_dst_cell):n_count =',n_count
        stop
      endif   
      rcut_book_tmp=rcut_book_wrk*(1.0d0+dble(itry-1)/20.0d0)
      if (i_pbc_sum /= 0) then
        if (rcut_book_tmp > 0.47d0*min(ax,ay,az)) then
          rcut_book_tmp = 0.47d0*min(ax,ay,az)
          reach_max_length=.true.
        endif  
      endif   
      if (debug_mode) then 
        booked(:)=.false.
        booked(atm_index)=.true.
      endif
      n_count=1
      atm_list(n_count)=atm_index
      do j_atom=1,length_of_atm_list_wrk
        atm_index2=atm_list_wrk(j_atom)
        dist=atm_distance_list_wrk(j_atom)
        if (dist < rcut_book_tmp) then 
          if (atm_index2 /= atm_index) then
            if (debug_mode) then 
              if (booked(atm_index2) .eqv. .true.) then
                write(*,*)'ERROR(set_rcut_kry_dst_cell):already booked'
                stop
              endif
            endif
            n_count=n_count+1
            if (n_count > size(atm_list,1)) then
              iexit=1 !!! Error detection !!
              write(*,*)'ERROR(set_rcut_kry_dst_cell):n_count=',n_count, size(atm_list,1)
              write(*,*)'ERROR(set_rcut_kry_dst_cell):itry=',itry
              write(*,*)'ERROR(set_rcut_kry_dst_cell):rcut_book_tmp=',rcut_book_tmp
              write(*,*)'ERROR(set_rcut_kry_dst_cell):length_of_atm_list_wrk=',length_of_atm_list_wrk
              write(*,*)'ERROR(set_rcut_kry_dst_cell):num_atm_list_min =',num_atm_list_min
              write(*,*)'ERROR: Too small projection_list_length ? or Too large cut_off ?'
              stop
            else  
              atm_list(n_count)=atm_index2
              if (debug_mode) booked(atm_index2) = .true.
            endif
          endif   
        endif  
      enddo ! for j_atom loop
!
      if (rcut_book_tmp > max_distance) then
        write(*,*)'ERROR:rcut_book_tmp, max_distance=',rcut_book_tmp, max_distance
        write(*,*)'         atm_index=',atm_index
        write(*,*)'length_of_atm_list    =',length_of_atm_list
        write(*,*)'length_of_atm_list_wrk=',length_of_atm_list_wrk
        write(*,*)'              itry=',itry
        write(*,*)'  rcut_book_wrk=',rcut_book_wrk
        write(*,*)'        n_count=',n_count
        iunit=vacant_unit()
        open (iunit,file='debug-info.txt',status='unknown')
          do j_atom=1,length_of_atm_list
            write(iunit,*)j_atom, atm_list_wrk(j_atom), atm_distance_list_wrk(j_atom)
          enddo
        close (iunit)
        stop
      endif
!
      if (n_count > noav) then
        write(*,*)'ERROR(set_rcut_kry_dst_cell):n_count,noav=',n_count,noav
        stop
      endif   
!
      if (n_count > num_atm_list_min) then
        if (iexit == 1) info=n_count   !!!! Error detection
        length_of_atm_list=n_count
        proj_radius=rcut_book_tmp
        exit
      endif   
!
    enddo ! for itry loop
!
  end subroutine set_projection_list_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_proj_list_dstm
