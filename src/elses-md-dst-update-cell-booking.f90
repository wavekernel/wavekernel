!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_dst_update_cell_bk
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_unit !(unchanged)
  use M_wall_clock_time,  only : get_system_clock_time !(routine)
  implicit none
!
  private
!
  public :: update_cell_booking
  public :: check_cell_booking
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Update cell booking
!
   subroutine update_cell_booking
!
     use M_md_dst_cell,     only : num_of_atoms_in_cell !(CHANGED)
     use M_md_dst_cell,     only : atom_list_in_cell    !(CHANGED)
     use M_md_dst_cell,     only : n_cell, n_cell_x, n_cell_y, n_cell_z  !(unchanged)
     use M_md_dst_cell,     only : one_cell_setting     !(unchanged)
     use M_md_dst,          only : myrank, nprocs       !(unchanged)
     use M_md_dst_cell,     only : get_neighbor_cell_list !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i1  !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i1c !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i2  !(routine)
     use M_md_dst_cell,     only : get_cell_index_for_atom   !(routine)
     use M_md_dst_get_atom_range, only : dst_get_index_range !(routine)
!
     implicit none
!
!    integer :: j_cell, jx, jy, jz
!    integer :: cell_index, sorted_index
!    real(8) :: rx, ry, rz
!    real(8) :: rcut1
!    integer :: max_num_of_atom_in_cell, size1
!
!    integer, allocatable :: booked_atoms_in_cell(:)
!    logical :: CalledFirst
!
!    integer, allocatable :: atom_list_wrk_omp(:,:)
!    integer, allocatable :: n_count_omp(:)
!    integer :: number_of_omp_threads, id_of_my_omp_thread, i_thread
!    integer :: omp_get_num_threads, omp_get_thread_num
!
     integer :: ierr
     integer :: j_cell_ini, j_cell_fin
     integer :: j_cell, jx, jy, jz
     integer :: jdx, jdy, jdz
     integer, allocatable :: cell_list(:)
     integer :: size_for_cell_list, k
     integer :: j_cell_old, j_cell_new
     integer :: atm_index_cell, atm_index
     logical :: result_value
     integer :: n_count
     integer, allocatable :: atom_list_in_cell_wrk(:,:)
     integer, allocatable :: num_of_atoms_in_cell_wrk(:)
     real(DOUBLE_PRECISION) :: time_wrk, time_wrk_prev
!
     if (one_cell_setting) then
       if (i_verbose >=0) write(*,*)'@@ update_cell_booking:one_cell_setting... skipped'
       return
     endif
!
     if (i_verbose >=0) write(*,*)'@@ update_cell_booking'
!
     jdx=1
     jdy=1
     jdz=1
     size_for_cell_list=(2*jdx+1)*(2*jdy+1)*(2*jdz+1)
!
     call dst_get_index_range(n_cell, j_cell_ini, j_cell_fin)
     j_cell_ini=j_cell_ini-1
     j_cell_fin=j_cell_fin-1
!    j_cell_ini=myrank*n_cell/nprocs
!    j_cell_fin=(myrank+1)*n_cell/nprocs-1
!       ----> decomposition of the loop for j_cell = 0, n_cell-1
!
     call get_system_clock_time(time_wrk)
     time_wrk_prev=time_wrk
!
     allocate( atom_list_in_cell_wrk &
&              (size(atom_list_in_cell,1), 0:size(atom_list_in_cell,2)-1 ), stat=ierr)
     if (ierr /= 0) then
       stop 'Alloc error in atom_list_in_cell_wrk'
     endif
!
     allocate(num_of_atoms_in_cell_wrk( 0:size(num_of_atoms_in_cell,1)-1 ), stat=ierr)
     if (ierr /= 0) then
       stop 'Alloc error in num_of_atoms__in_cell_wrk'
     endif
!
     num_of_atoms_in_cell_wrk(:)=0
     atom_list_in_cell_wrk(:,:)=0
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a1',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a1',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!$omp   parallel default(shared) &
!$omp&  private (j_cell, cell_list, ierr) &
!$omp&  private (n_count, k, j_cell_old, atm_index_cell, atm_index, j_cell_new)
!$omp   do schedule(static)
     do j_cell = j_cell_ini, j_cell_fin
!
       allocate(cell_list(size_for_cell_list), stat=ierr)
       if (ierr /= 0) then
         write(*,*)'Error in alloc. of cell_list'
         stop
       endif
!
       call get_neighbor_cell_list(j_cell, jdx, jdy, jdz, cell_list)
!        ----> get the neighbor cell list
!
       n_count=0
       do k = 1, size_for_cell_list
         j_cell_old = cell_list(k)
         if ((j_cell_old < 0) .or. (j_cell_old > n_cell-1)) then
           write(*,*)'ERROR(update_cell_booking):j_cell_old=',j_cell_old
           stop
         endif
!
         do atm_index_cell = 1,num_of_atoms_in_cell(j_cell_old)
           atm_index = atom_list_in_cell(atm_index_cell, j_cell_old)
           call get_cell_index_for_atom(atm_index, j_cell_new)
!            -----> Get the present cell index for the atom
           if (j_cell == j_cell_new) then
             n_count=n_count+1
             if (n_count > size(atom_list_in_cell,1)) then
               write(*,*)'ERROR(update_cell_booking):j_cell, n_count=',j_cell,n_count
               stop
             endif
             if (j_cell < 10 ) then
!              write(*,'(a,4i10)')'j_cell, j_cell_old, atm_index, n_count=', &
!&                                     j_cell, j_cell_old, atm_index, n_count
             endif
             atom_list_in_cell_wrk(n_count, j_cell)=atm_index
           endif
         enddo
!
       enddo
       num_of_atoms_in_cell_wrk(j_cell)=n_count
!       write(*,'(a,2i10)')'j_cell, num_of_atoms_in_cell_wrk=', & 
!&                          j_cell, num_of_atoms_in_cell_wrk(j_cell)
!
       deallocate(cell_list, stat=ierr)
       if (ierr /= 0) then
         write(*,*)'Error in dealloc. of cell_list'
         stop
       endif
!
     enddo  ! end of j_cell loop
!$omp end do
!$omp end parallel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a2',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a2',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
     atom_list_in_cell(:,:)=0
     atom_list_in_cell(:,:)=atom_list_in_cell_wrk(:,:)
!
!    num_of_atoms_in_cell(:)=0
!    num_of_atoms_in_cell(:)=num_of_atoms_in_cell_wrk(:)
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3a',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3a',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
!    call mpi_wrapper_allreduce_i1(num_of_atoms_in_cell)
     call mpi_wrapper_allreduce_i1c(num_of_atoms_in_cell_wrk, num_of_atoms_in_cell)
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3bn',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,i10)')'size for allreduce =',size(num_of_atoms_in_cell_wrk,1)
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3bn',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
     call mpi_wrapper_allreduce_i2(atom_list_in_cell)
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3c',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3c',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
     deallocate(atom_list_in_cell_wrk, stat=ierr)
     if (ierr /= 0) then
       stop 'Dealloc error in atom_list_in_cell_wrk'
     endif
!
     deallocate(num_of_atoms_in_cell_wrk, stat=ierr)
     if (ierr /= 0) then
       stop 'Dealloc error in num_of_atoms_in_cell_wrk'
     endif
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3d',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a3d',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     call set_cell_for_atom
!       ------> Set cell_for_atom(1:noav)
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
!    do j_cell = 0, n_cell-1
!      write(*,*)'num_of_atoms_in_cell=',j_cell, num_of_atoms_in_cell(j_cell)
!    enddo
!
   end subroutine update_cell_booking
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_cell_for_atom
!
     use M_md_dst_cell,     only : cell_for_atom        !(CHANGED)
     use M_md_dst_cell,     only : num_of_atoms_in_cell !(unchanged)
     use M_md_dst_cell,     only : atom_list_in_cell    !(unchanged)
     use M_md_dst_cell,     only : n_cell, n_cell_x, n_cell_y, n_cell_z  !(unchanged)
     use M_md_dst_cell,     only : one_cell_setting     !(unchanged)
     use M_md_dst,          only : myrank, nprocs       !(unchanged)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i1  !(routine)
     use M_md_dst_get_atom_range, only : dst_get_index_range !(routine)

!
     implicit none
     integer :: j_cell_ini, j_cell_fin, ierr
     integer :: atm_index_cell, atm_index, j_cell
     logical, allocatable :: flag_for_set_cell(:)
     real(DOUBLE_PRECISION) :: time_wrk, time_wrk_prev
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if (one_cell_setting) then
       write(*,*)'ERROR(set_cell_for_atom):one_cell_setting=',one_cell_setting
       stop
     endif
!
     if (i_verbose >=1) write(*,*)'@@ set_cell_for_atom'
!
     call get_system_clock_time(time_wrk)
     time_wrk_prev=time_wrk
!
     call dst_get_index_range(n_cell, j_cell_ini, j_cell_fin)
     j_cell_ini=j_cell_ini-1
     j_cell_fin=j_cell_fin-1
!    j_cell_ini=myrank*n_cell/nprocs
!    j_cell_fin=(myrank+1)*n_cell/nprocs-1
!       ----> decomposition of the loop for j_cell = 0, n_cell-1
     cell_for_atom(:)=0
!
     allocate(flag_for_set_cell(size(cell_for_atom)), stat=ierr)
     if (ierr /= 0) stop 'ERROR!(set_cell_for_atom):flag_for_cell'
     flag_for_set_cell(:)=.false.
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4a',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4a',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!$omp   parallel default(shared) &
!$omp&  private (j_cell, atm_index_cell, atm_index)
!$omp   do schedule(static)
     do j_cell = j_cell_ini, j_cell_fin
!
       do atm_index_cell = 1,num_of_atoms_in_cell(j_cell)
         atm_index = atom_list_in_cell(atm_index_cell, j_cell)
         if (flag_for_set_cell(atm_index) .eqv. .true.) then
           write(*,*)'ERROR(update_cell_booking:flag):atm_index=',atm_index
           stop
         endif
         flag_for_set_cell(atm_index)=.true.
         cell_for_atom(atm_index)=j_cell
       enddo
!
     enddo  ! end of j_cell loop
!$omp end do
!$omp end parallel
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4b',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4b',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
     deallocate(flag_for_set_cell, stat=ierr)
     if (ierr /= 0) stop 'Dealloc error!(set_cell_for_atom):flag_for_cell'
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4c',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4c',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
     call mpi_wrapper_allreduce_i1(cell_for_atom)
!
     call get_system_clock_time(time_wrk)
!    write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4d',time_wrk-time_wrk_prev
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:a4d',time_wrk-time_wrk_prev
     endif
     time_wrk_prev=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine set_cell_for_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine check_cell_booking  !Routine for debugging
!
     use M_md_dst_cell,     only : cell_for_atom        !(unchanged)
     use M_md_dst_cell,     only : num_of_atoms_in_cell !(unchanged)
     use M_md_dst_cell,     only : atom_list_in_cell    !(unchanged)
     use M_md_dst_cell,     only : n_cell               !(unchanged)
!    use M_md_dst_cell,     only : one_cell_setting     !(unchanged)
     use M_qm_domain,       only : noav  !(unchanged)
     implicit none
     integer :: atm_index, j_cell, ierr
!
!    if (one_cell_setting) then
!      write(*,*)'ERROR(check_cell_booking):one_cell_setting=',one_cell_setting
!      stop
!    endif
!
     if (i_verbose >=2) then 
       write(*,'(a)')'@@ Check_cell_booking (for debugging)(Non-parallel routine)'
     else
       if (i_verbose >=1) write(*,'(a)') &
&         '@@ Check_cell_booking (for debugging) ...is skipped. Set (verbose level) > =2 for checking.'
       return
     endif  
!
     ierr=0
!
     do atm_index=1,noav
       j_cell=cell_for_atom(atm_index) 
       if ((j_cell < 0) .or. (j_cell > n_cell-1)) then
         write(*,'(a,2i10)')'ERROR!(check_cell_booking):atm_index,j_cell=',atm_index, j_cell
         ierr=1
       endif
     enddo
!
     if (sum(num_of_atoms_in_cell) /= noav) then
       write(*,'(a,i10)')'ERROR!(check_cell_booking):sum=',sum(num_of_atoms_in_cell)
       ierr=1
     endif
!
     if (ierr == 1) then
       stop 'Stop in check_cell_booking'
     endif
!
   end subroutine check_cell_booking
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module M_md_dst_update_cell_bk
