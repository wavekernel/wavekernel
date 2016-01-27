!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_set_bk_list_w_cell
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_wall_clock_time,  only : get_system_clock_time !(routine)
   implicit none
!
   private
!
! Public module variable
!
! Public routine
   public :: make_bk_list_w_cell
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Set the booking list with the cell technique
!
   subroutine make_bk_list_w_cell(imode, cutoff_radius, length_of_list, booking_list)
!
     use M_qm_domain, only : noav, ax, ay, az             !(unchanged)
     use M_qm_domain, only : i_pbc_x, i_pbc_y, i_pbc_z    !(unchanged)
     use M_qm_domain, only : atm_position                 !(unchanged)
     use M_md_dst,    only : myrank, nprocs               !(unchanged)
!
     use M_md_dst_cell, only : n_cell, n_cell_x, n_cell_y, n_cell_z !(unchanged)
     use M_md_dst_cell, only : cell_for_atom, num_of_atoms_in_cell  !(unchanged)
     use M_md_dst_cell, only : one_cell_setting                     !(unchanged)
!
     use M_md_dst_cell, only : num_of_atoms_in_cell, atom_list_in_cell !(routine)
     use M_md_get_distance, only : get_distance                        !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i1  !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i2  !(routine)
     use M_md_dst_get_atom_range, only : dst_get_atm_index_range !(routine)
!
     implicit none
     integer, intent(in)    :: imode
     integer, intent(inout) :: length_of_list(:)
     integer, optional      :: booking_list(:,:)
!
     real(DOUBLE_PRECISION), intent(in)  :: cutoff_radius
     real(DOUBLE_PRECISION) ::  dx,  dy,  dz 
!
     real(DOUBLE_PRECISION) :: rcut1, dist
     integer                :: jdx, jdy, jdz, ierr
     integer                :: j, jx, jy, jz, ix, iy, iz
     integer                :: k, kx, ky, kz
     integer                :: j_atom, k_atom, atm_index1, atm_index2, atm_index
     integer, allocatable   :: num_of_booked_atoms(:)
     real(DOUBLE_PRECISION) :: time_wrk, time_wrk_previous
     integer                :: atm_index_ini, atm_index_fin
!
!
     if (i_verbose >=1) then
       write(*,*)'@@ make_bk_list_w_cell : imode=',imode
       write(*,*)' Warning : Only PBC is supported now'
       if (one_cell_setting) write(*,*)'one_cell_setting = ',one_cell_setting
       write(*,*)'Stop manually(make_bk_list_w_cell)'
       stop
     endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Error checking 
!
     if (size(length_of_list,1) /= noav) &
&       stop 'ERROR(set_booking_list_w_cell):wrong array size'
!
     select case(imode)
      case(1) 
        if (i_verbose >=1) write(*,*)' Mode 1: Set the array size for booking list'
      case(2)
        if (i_verbose >=1) write(*,*)' Mode 2: Construct the booking list'
        if (.not. present(booking_list)) then 
          write(*,'(a,i5)')'ERROR(set_booking_list_w_cell):booking list is missing:imode=',imode
          stop
        endif
        if (size(booking_list,2) /= noav) then
          write(*,'(a,i5)') 'ERROR(set_booking_list_w_cell):wrong array size of booking_list:imode=',imode
          stop
        endif
      case default
        stop 'ERROR in set_booking_list_w_cell (mode setting)'
     end select  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Parameter setting
!
     if( cutoff_radius .lt. huge(1d0) ) then
       rcut1=cutoff_radius*1.001d0
     else
       rcut1=cutoff_radius
     endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocate working arrays
!
     allocate (num_of_booked_atoms(noav),stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. : num_of_booked_atoms'
     num_of_booked_atoms(:)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if (i_verbose >=1) then
      write(*,*)'@@ initial_cell_booking:myrank,nprocs=',myrank,nprocs
      write(*,*)'noav=',noav
      write(*,*)'n_cell_x=',n_cell_x
      write(*,*)'n_cell_y=',n_cell_y
      write(*,*)'n_cell_z=',n_cell_z
      write(*,*)' rcut1  =',rcut1
      write(*,*)' ax     =',ax
      write(*,*)' ay     =',ay
      write(*,*)' az     =',az
     endif
!
     if (i_pbc_x*i_pbc_y*i_pbc_z /= 1) then
      write(*,*)'ERROR:non periodic system is not suppored in Cell decomposition'
      stop
     endif
!
     jdx=int(aint(rcut1/ax*dble(n_cell_x)))+1
     jdy=int(aint(rcut1/ay*dble(n_cell_y)))+1
     jdz=int(aint(rcut1/az*dble(n_cell_z)))+1
!
     write(*,'(a,i10)')'jdx =',jdx
     write(*,'(a,i10)')'jdy =',jdy
     write(*,'(a,i10)')'jdz =',jdz
!
     ierr=0
     if (jdx > n_cell_x/2) ierr=1
     if (jdy > n_cell_y/2) ierr=1
     if (jdz > n_cell_z/2) ierr=1
     if (ierr == 1) then
      write(*,*)'ERROR in get_length_of_list_cell:Too small cell'
      stop
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     call get_system_clock_time(time_wrk)
     time_wrk_previous=time_wrk
!
!    if (imode == 1) then
!      num_of_booked_atoms(:)=0
!    endif
!  
!    if (imode == 2) then
!      booking_list(:,:)=0
!      num_of_booked_atoms(:)=1
!      do atm_index1=1,noav
!        booking_list(1,atm_index1)=atm_index1
!      enddo
!    endif
!
     num_of_booked_atoms(:)=0
     if (imode == 2) booking_list(:,:)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    atm_index_ini=myrank*noav/nprocs+1
!    atm_index_fin=(myrank+1)*noav/nprocs
!
     call dst_get_atm_index_range(atm_index_ini,atm_index_fin)
!
!$omp  parallel default(shared) &     
!$omp& private (atm_index1, j) &
!$omp& private (jz, jy, jx, iz, iy, ix, kz, ky, kx, k, k_atom, atm_index2) &
!$omp& private (dist)
!$omp  do schedule(static)
     do atm_index1=atm_index_ini, atm_index_fin
!
       num_of_booked_atoms(atm_index1)=1
       if (imode == 2) then 
         booking_list(:,atm_index1)=0
         booking_list(1,atm_index1)=atm_index1
       endif
!
       j=cell_for_atom(atm_index1)
       if ((j < 0) .or. (j > n_cell-1)) then
         write(*,*)'ERROR(make_bk_list_w_cell):atm_index1,j_cell=',atm_index1,j
       endif
!
      jz=j/(n_cell_x*n_cell_y)
      jy=(j-jz*n_cell_x*n_cell_y)/n_cell_x
      jx= j - jz*n_cell_x*n_cell_y - jy*n_cell_x
!     write(*,'(a,4i7)')'j,jx,jy,jz =',j,jx,jy,jz
!
      do iz=-jdz,jdz
       kz=jz+iz
       if (kz < 0)          kz = kz + n_cell_z
       if (kz > n_cell_z-1) kz = kz - n_cell_z
       do iy=-jdy,jdz
        ky=jy+iy
        if (ky < 0)          ky = ky + n_cell_y
        if (ky > n_cell_y-1) ky = ky - n_cell_y
        do ix=-jdx,jdx
         kx=jx+ix
         if (kx < 0)          kx = kx + n_cell_x
         if (kx > n_cell_x-1) kx = kx - n_cell_x
         k=kx+ky*n_cell_x+kz*n_cell_x*n_cell_y
!        write(*,'(a,4i7)')'k,kx,ky,kz=',k,kx,ky,kz
!
         do k_atom=1,num_of_atoms_in_cell(k)
           atm_index2=atom_list_in_cell(k_atom,k)
           if ((atm_index2 <= 0) .or. (atm_index2 > noav)) then
            write(*,*)'ERROR:get_length_of_list_cell'
            write(*,*)' k_atom, atm_index2 =', k_atom, atm_index2
            stop
           endif   
!
           call get_distance(atm_index1, atm_index2, dist)
           if (dist < rcut1) then 
             if (atm_index1 /= atm_index2) then
               num_of_booked_atoms(atm_index1)=num_of_booked_atoms(atm_index1)+1
               if (imode == 2) then
                 booking_list(num_of_booked_atoms(atm_index1),atm_index1)=atm_index2
               endif
             endif  
           endif
!
         enddo ! end of k_atom loop
!
        enddo ! end of ix loop
       enddo  ! end of iy loop
      enddo   ! end of iz loop
     enddo   ! end of atm_index1 loop
!$omp end do
!$omp end parallel
!
     call get_system_clock_time(time_wrk)
     write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:bw:OMP',time_wrk-time_wrk_previous
     if (log_unit > 0) then
       write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:bw:OMP',time_wrk-time_wrk_previous
     endif
     time_wrk_previous=time_wrk
!
     call mpi_wrapper_allreduce_i1(num_of_booked_atoms)
     length_of_list(:)=num_of_booked_atoms(:)
!
     if (imode == 2) then
       call mpi_wrapper_allreduce_i2(booking_list)
     endif
!
     call get_system_clock_time(time_wrk)
     if (imode == 1) then
       write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:bw1:MPI',time_wrk-time_wrk_previous
       if (log_unit > 0) then
         write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:bw1:MPI',time_wrk-time_wrk_previous
       endif
     endif
     if (imode == 2) then
       write(*,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:bw2:MPI',time_wrk-time_wrk_previous
       if (log_unit > 0) then
         write(log_unit,'(a,f20.10)')'TIME:qm_engine_dst:qm_domain_setting_dst:bw2:MPI',time_wrk-time_wrk_previous
       endif
     endif
     time_wrk_previous=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Checking (optional)
!
     if (imode == 1) then
!      do atm_index=1,noav
!       write(*,*)'atom, length_of_list=',atm_index, length_of_list(atm_index)
!      enddo
       write(*,*)'max and min length of the list=',maxval(length_of_list), minval(length_of_list)
     endif
!
     if (imode == 2) then
       do atm_index=1,noav
        if (length_of_list(atm_index) /= num_of_booked_atoms(atm_index)) then
          write(*,*)'ERROR in length_of_list=', atm_index, length_of_list(atm_index), & 
&                                               num_of_booked_atoms(atm_index)
          stop
        endif
       enddo
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
   end subroutine make_bk_list_w_cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_md_set_bk_list_w_cell


