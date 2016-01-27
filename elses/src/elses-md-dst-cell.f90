!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_dst_cell
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
  integer :: n_cell, n_cell_x, n_cell_y, n_cell_z   ! set in 'initial_cell_booking'
  integer :: dn_cell_x, dn_cell_y, dn_cell_z        ! set in 'initial_cell_booking' 
  integer, allocatable :: cell_for_atom(:)
  integer, allocatable :: num_of_atoms_in_cell(:)
  integer, allocatable :: atom_list_in_cell(:,:)
  logical :: one_cell_setting
  logical :: mcell_booking_initial                  ! true for initial creation of booking list
!
  private
!
  public :: one_cell_setting
  public :: mcell_booking_initial 
!
  public :: n_cell, n_cell_x, n_cell_y, n_cell_z
  public :: dn_cell_x, dn_cell_y, dn_cell_z
  public :: cell_for_atom, num_of_atoms_in_cell
  public :: atom_list_in_cell
!
  public :: initial_mcell_booking
! public :: get_j_cell
! public :: get_j_cell_xyz
  public :: get_neighbor_cell_list
  public :: get_cell_index_for_atom
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Initial cell booking
!
   subroutine initial_mcell_booking
!
     use M_config,          only : config   !(unchanged)(for micro_cell_booking)
     use M_lib_phys_const,  only : angst                     !(unchanged)
     use M_qm_domain,       only : noav, ax, ay, az          !(unchanged)
     use M_qm_domain,       only : i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)
!    use M_qm_domain,       only : atm_position              !(unchanged)
     use M_md_dst,          only : myrank, nprocs            !(unchanged)
     use elses_mod_sel_sys, only : r_cut_book                !(unchanged)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i1  !(routine)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_i2  !(routine)
     use M_md_dst_get_atom_range, only : dst_get_atm_index_range !(routine)
     use M_md_dst_get_atom_range, only : dst_get_index_range     !(routine)
     use M_lib_dst_info,          only : log_unit            !(unchanged)
!
     implicit none
!
     integer :: atm_index, ierr
     integer :: j_cell, jx, jy, jz
     integer :: atm_index_ini, atm_index_fin
     integer :: cell_index, sorted_index
!    real(8) :: rx, ry, rz
!    real(8) :: rcut1
     integer :: max_num_of_atom_in_cell, size1
!
     integer, allocatable :: booked_atoms_in_cell(:)
!    logical :: CalledFirst
!    real(8) :: qx, qy, qz
!
     integer, allocatable :: atom_list_wrk_omp(:,:)
     integer, allocatable :: n_count_omp(:)
     integer :: j_cell_ini, j_cell_fin, n_count, j
     integer :: number_of_omp_threads, id_of_my_omp_thread, i_thread
     integer :: omp_get_num_threads, omp_get_thread_num
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the initial flag
!
     if(allocated(cell_for_atom)) mcell_booking_initial = .true.
!
     if (i_verbose >=0) then
       write(*,*)'@@ initial_mcell_booking:myrank,nprocs=',myrank,nprocs
       write(*,*)'noav=',noav
!      write(*,*)'r_cut_book=',r_cut_book
     endif
!
     call set_paramters_for_booking
!       ----> one_cell_setting, n_cell_x, n_cell_y, n_cell_z, n_cell, 
!             dn_cell_x, dn_cell_y, dn_cell_z
!
     if (one_cell_setting) then
       call initial_cell_booking_one_cell
!      write(*,*)'Stop manually'
       return
     endif
!
     if (i_verbose >= 0) then
       write(*,'(a,i10)')       '  INFO-CELL:n_cell=', n_cell
       write(*,'(a,i10,2f20.10)')'  INFO-CELL:n_cell_x, ax/(n_cell_x)[au,A]=', n_cell_x, ax/dble(n_cell_x), ax/dble(n_cell_x)*angst
       write(*,'(a,i10,2f20.10)')'  INFO-CELL:n_cell_y, ay/(n_cell_y)[au,A]=', n_cell_y, ay/dble(n_cell_y), ax/dble(n_cell_y)*angst
       write(*,'(a,i10,2f20.10)')'  INFO-CELL:n_cell_z, az/(n_cell_z)[au,A]=', n_cell_z, az/dble(n_cell_z), ax/dble(n_cell_z)*angst
       write(*,'(a,i10)')       '  INFO-CELL:  dn_cell_x =', dn_cell_x
       write(*,'(a,i10)')       '  INFO-CELL:  dn_cell_y =', dn_cell_y
       write(*,'(a,i10)')       '  INFO-CELL:  dn_cell_z =', dn_cell_z
     endif  
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate the arrays
!
     if (.not. allocated(cell_for_atom)) then
       allocate(cell_for_atom(noav), stat=ierr)
       if ( ierr /=  0 ) then
         write(*,*)'ERROR in alloc. : cell_for_atom' 
         stop
       endif   
     endif  
     cell_for_atom(:)=0 ! dummy setting
!
     if(mcell_booking_initial) then
       deallocate(num_of_atoms_in_cell, stat=ierr)
       if ( ierr /=  0 ) then
        write(*,*)'ERROR in dealloc. : num_of_atoms_in_cell'
        stop
       endif   
     endif   
!
     if (.not. allocated(num_of_atoms_in_cell)) then
       allocate(num_of_atoms_in_cell(0:n_cell-1), stat=ierr)
       if ( ierr /=  0 ) then
        write(*,*)'ERROR in alloc. : num_of_atoms_in_cell'
        stop
       endif   
     endif  
     num_of_atoms_in_cell(:)=0 ! dummy setting0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     call dst_get_atm_index_range(atm_index_ini, atm_index_fin)
     if (log_unit > 0) then
       write(log_unit,'(a,3i10)')'INFO-MPI(initial_mcell_booking): myrank, ini, fin=',myrank, atm_index_ini, atm_index_fin
     endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set cell_for_atom(1:noav)
!
     cell_for_atom(:)=0 ! Initial clear (essential for the MPI allreduce command)
!
!$omp  parallel default(shared) &
!$omp& private (atm_index, j_cell) 
!$omp  do schedule(static)
     do atm_index=atm_index_ini, atm_index_fin
!
       call get_cell_index_for_atom(atm_index, j_cell)
       cell_for_atom(atm_index)=j_cell
!
     enddo   
!$omp end do
!$omp end parallel
!
     call mpi_wrapper_allreduce_i1(cell_for_atom)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check the MPI parameters 
!
     if (nprocs > n_cell) then
       write(*,*)'ERROR(initial_cell_booking):nprocs,n_cell=',nprocs,n_cell
       stop
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set num_of_atoms_in_cell(0:n_cell-1)
!
     num_of_atoms_in_cell(:)=0 ! Initial clear (essential for the MPI allreduce command)
!
     call dst_get_index_range(n_cell, j_cell_ini, j_cell_fin)
     j_cell_ini=j_cell_ini-1 
     j_cell_fin=j_cell_fin-1 
!    j_cell_ini=myrank*n_cell/nprocs
!    j_cell_fin=(myrank+1)*n_cell/nprocs-1
!       ----> decomposition of the loop for j_cell = 0, n_cell-1
!
     if (log_unit >0) then 
       write(log_unit,'(a,3i10)')'INFO-MPI:(initial_mcell_booking) myrank, j_cell_ini, fin=',myrank, j_cell_ini, j_cell_fin
     endif
!
     do j_cell = j_cell_ini, j_cell_fin
!
       n_count=0
!$omp  parallel default(shared) &
!$omp& private (atm_index) &
!$omp& reduction (+ : n_count)
!$omp  do schedule(static)
       do atm_index=1,noav
         if ( cell_for_atom(atm_index) == j_cell) n_count=n_count+1
       enddo
!$omp end do
!$omp end parallel
      num_of_atoms_in_cell(j_cell)=n_count
!
     enddo  
!
     call mpi_wrapper_allreduce_i1(num_of_atoms_in_cell)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check num_of_atoms_in_cell(0:n_cell-1)
!
     if (sum(num_of_atoms_in_cell) /= noav) then
       write(*,*)'ERROR(initial_cell_booking)'
       write(*,*)' sum of atoms in cell (should be noav)=',sum(num_of_atoms_in_cell)
       write(*,*)' Debug info: debug-info.txt, debug-info2.txt'
       open (99, file='debug-info.txt',status='unknown')
       open (98, file='debug-info2.txt',status='unknown')
       do atm_index = 1, noav
         write(99,*) atm_index, cell_for_atom(atm_index)
       enddo
       do j_cell = 0 , n_cell-1
         write(98,*) j_cell, num_of_atoms_in_cell(j_cell)
       enddo
       stop
     else
       if (i_verbose >=0) write(*,*)'Check : sum of atoms in cell (should be noav)...OK!'
     endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocation of atom_list_in_cell(1:size1, 0:n_cell
!
     max_num_of_atom_in_cell=maxval(num_of_atoms_in_cell)
     size1=int(max_num_of_atom_in_cell*1.5) ! NOTE: The factor of 1.5 is tolerance factor
     if (i_verbose >= 1) write(*,*)'max. number of the atoms in a cell, allocated size =',max_num_of_atom_in_cell, size1
!
     if ( allocated(atom_list_in_cell) ) then
       if ( size1 < size(atom_list_in_cell,1) ) then
         deallocate(atom_list_in_cell, stat=ierr)
         if ( ierr /=  0 ) then
           write(*,*)'ERROR in dealloc. : atom_list_in_cell'
           stop
         endif   
       endif
     endif
!
     if ( .not. allocated(atom_list_in_cell) ) then
       allocate(atom_list_in_cell(size1, 0:n_cell-1), stat=ierr)
       if ( ierr /=  0 ) then
         write(*,*)'ERROR in alloc. : atom_list_in_cell'
         stop
       endif   
     endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allcoate the work array for the OMP loop
!
     number_of_omp_threads=1
!$omp parallel default(shared)
!$     number_of_omp_threads=omp_get_num_threads()
!$omp end parallel
!    write(*,*)'number_of_omp_threads=',number_of_omp_threads
!
     allocate(atom_list_wrk_omp(size(atom_list_in_cell,1),0:number_of_omp_threads-1), stat=ierr)
     if ( ierr /=  0 ) then
       write(*,*)'ERROR in alloc. : atom_list_wrk_omp'
       stop
     endif   
     atom_list_wrk_omp(:,:)=0
!
     allocate(n_count_omp(0:number_of_omp_threads-1), stat=ierr)
     if ( ierr /=  0 ) then
       write(*,*)'ERROR in alloc. : atom_list_wrk_omp'
       stop
     endif   
     n_count_omp(:)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set atom_list_in_cell(size1, 0:n_cell-1)
!
      atom_list_in_cell(:,:)=0
!
      do j_cell = j_cell_ini, j_cell_fin
!
        atom_list_wrk_omp(:,:)=0
!$omp   parallel default(shared) &
!$omp&  private (atm_index, id_of_my_omp_thread) 
        id_of_my_omp_thread=0
!$      id_of_my_omp_thread=omp_get_thread_num()
!$      if (number_of_omp_threads /= omp_get_num_threads()) then
!$        write(*,*)'ERROR(OMP loop in initial cell booking)=', number_of_omp_threads, omp_get_num_threads()
!$        stop
!$      endif
        n_count_omp(id_of_my_omp_thread)=0
!$omp   do schedule(static)
        do atm_index=1,noav
          if ( cell_for_atom(atm_index) == j_cell) then 
            n_count_omp(id_of_my_omp_thread)=n_count_omp(id_of_my_omp_thread)+1
            atom_list_wrk_omp(n_count_omp(id_of_my_omp_thread),id_of_my_omp_thread)=atm_index
          endif
        enddo
!$omp end do
!$omp end parallel
!
        do i_thread=0,number_of_omp_threads-1
          if (i_thread == 0) then 
            n_count=0
          else
            n_count=sum(n_count_omp(0:i_thread-1))
          endif
          do j=1,n_count_omp(i_thread)
            atom_list_in_cell(n_count+j,j_cell)=atom_list_wrk_omp(j,i_thread)
          enddo
        enddo
!
     enddo  
!
     call mpi_wrapper_allreduce_i2(atom_list_in_cell)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallcoate the work array for the OMP loop
!
     deallocate(atom_list_wrk_omp,stat=ierr)
     if ( ierr /=  0 ) then
       write(*,*)'ERROR in dealloc. : atom_list_wrk'
       stop
     endif   
!
     deallocate(n_count_omp,stat=ierr)
     if ( ierr /=  0 ) then
       write(*,*)'ERROR in dealloc. : n_count_omp'
       stop
     endif   
!
!    allocate(booked_atoms_in_cell(0:n_cell-1), stat=ierr)
!    if ( ierr /=  0 ) then
!      write(*,*)'ERROR in alloc. : booked_atom_in_cell'
!    endif   
!
!    atom_list_in_cell(:,:)=0
!    booked_atoms_in_cell(:)=0
!    do atm_index=1,noav ! Not parallelizable
!      j_cell=cell_for_atom(atm_index)
!      booked_atoms_in_cell(j_cell)=booked_atoms_in_cell(j_cell)+1
!      atom_list_in_cell(booked_atoms_in_cell(j_cell), j_cell)=atm_index
!    enddo   
!
!    deallocate(booked_atoms_in_cell, stat=ierr)
!    if ( ierr /=  0 ) then
!      write(*,*)'ERROR in dealloc. : booked_atom_in_cell'
!      stop
!    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    stop 'stop manually'
!
   end subroutine initial_mcell_booking
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Set the parameter for booking; 
!    ----> one_cell_setting, n_cell_x, n_cell_y, n_cell_z, n_cell, dn_cell_x, dn_cell_y, dn_cell_z
!
   subroutine set_paramters_for_booking
     use M_config,          only : config      !(unchanged)(for micro_cell_booking)
     use M_qm_domain,       only : ax, ay, az  !(unchanged)
     use M_qm_domain,       only : i_pbc_x, i_pbc_y, i_pbc_z  !(unchanged)
     use elses_mod_sel_sys, only : r_cut_book  !(unchanged)
     implicit none
     real(8) :: rcut1
     real(8) :: qx, qy, qz
!
     one_cell_setting = .false.
!
     if ( i_pbc_x + i_pbc_y + i_pbc_z /= 3 ) then
       one_cell_setting = .true. 
       if (i_verbose >= 1) then
         write(*,*)'Warning(DST):one_cell_setting is used for non-periodic system' 
       endif   
       return
     endif   
!
     rcut1=r_cut_book
     if (rcut1 > max(ax,ay,az)*10d0) then
       write(*,*)'ERROR(initial_cell_booking):Not suppoted:rcut1=',rcut1
       write(*,*)'Too large cutoff distance or non-periodic case'
       stop
     endif
!   
     qx = log(ax/(2.0d0*rcut1)) / log(2.0d0)
     qy = log(ay/(2.0d0*rcut1)) / log(2.0d0)
     qz = log(az/(2.0d0*rcut1)) / log(2.0d0)
!
     n_cell_x=2**(nint(aint(qx)))         ! TEMPORAL SETTING
     n_cell_y=2**(nint(aint(qy)))         ! TEMPORAL SETTING
     n_cell_z=2**(nint(aint(qz)))         ! TEMPORAL SETTING
     n_cell=n_cell_x*n_cell_y*n_cell_z
!
     rcut1=min(ax,ay,az)/4.0d0                 ! TEMPORAL SETTING
     rcut1=min(rcut1, r_cut_book*4.0d0)        ! TEMPORAL SETTING
!
     dn_cell_x=int(aint(rcut1/ax*dble(n_cell_x)))+1 ! TEMPORAL SETTING
     dn_cell_y=int(aint(rcut1/ay*dble(n_cell_y)))+1 ! TEMPORAL SETTING
     dn_cell_z=int(aint(rcut1/az*dble(n_cell_z)))+1 ! TEMPORAL SETTING
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the parameters from the XML tag, if available
!
     if (config%calc%distributed%micro_cell_booking%set .eqv. .true.) then
       if (config%calc%distributed%micro_cell_booking%mode == "inactive") then
         if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML) micro_cell_booking :inactive tag'
       else   
!
         if (config%calc%distributed%micro_cell_booking%mode == "one_cell") then
           one_cell_setting = .true. 
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML)one_cell_setting=',one_cell_setting
         endif   
!
         if (config%calc%distributed%micro_cell_booking%cell_number(1) > 0) then
           n_cell_x=config%calc%distributed%micro_cell_booking%cell_number(1) 
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML)  n_cell_x=',n_cell_x
         endif   
!
         if (config%calc%distributed%micro_cell_booking%cell_number(2) > 0) then
           n_cell_y=config%calc%distributed%micro_cell_booking%cell_number(2)
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML)  n_cell_y=',n_cell_y
         endif   
!
         if (config%calc%distributed%micro_cell_booking%cell_number(3) > 0) then
           n_cell_z=config%calc%distributed%micro_cell_booking%cell_number(3)
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML)  n_cell_z=',n_cell_z
         endif   
!
         n_cell=n_cell_x*n_cell_y*n_cell_z
!
         if (config%calc%distributed%micro_cell_booking%search_range(1) > 0) then
           dn_cell_x=config%calc%distributed%micro_cell_booking%search_range(1) 
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML) dn_cell_x=',dn_cell_x
         endif   
!
         if (config%calc%distributed%micro_cell_booking%search_range(2) > 0) then
           dn_cell_y=config%calc%distributed%micro_cell_booking%search_range(2) 
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML) dn_cell_y=',dn_cell_y
         endif   
!
         if (config%calc%distributed%micro_cell_booking%search_range(3) > 0) then
           dn_cell_z=config%calc%distributed%micro_cell_booking%search_range(3) 
           if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking, from XML) dn_cell_z=',dn_cell_z
         endif   
!
       endif  
     endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if (dn_cell_x*2+1 > n_cell_x) then
       if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking):dn_cell_x, n_cell_x=',dn_cell_x, n_cell_x
       one_cell_setting = .true.
!      stop
     endif
!
     if (dn_cell_y*2+1 > n_cell_y) then
       if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking):dn_cell_y, n_cell_y=',dn_cell_y, n_cell_y
       one_cell_setting = .true.
!      stop
     endif
!
     if (dn_cell_z*2+1 > n_cell_z) then
       if (i_verbose >=0) write(*,*)'INFO(initial_cell_booking):dn_cell_z, n_cell_z=',dn_cell_z, n_cell_z
       one_cell_setting = .true.
!      stop
     endif
!
   end subroutine set_paramters_for_booking
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_j_cell(j_cell, j_cell_x, j_cell_y, j_cell_z)
     implicit none
     integer, intent(out)   :: j_cell
     integer, intent(in)    :: j_cell_x, j_cell_y, j_cell_z
!
     if (one_cell_setting) then
       write(*,*)'ERROR(get_j_cell):one_cell_setting=',one_cell_setting
       stop
     endif
!
     j_cell = j_cell_x + j_cell_y * n_cell_x + j_cell_z * n_cell_x * n_cell_y
!
   end subroutine get_j_cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_j_cell_xyz(j_cell, j_cell_x, j_cell_y, j_cell_z)
     implicit none
     integer, intent(in)    :: j_cell
     integer, intent(out)   :: j_cell_x, j_cell_y, j_cell_z
!
     if (one_cell_setting) then
       write(*,*)'ERROR(get_j_cell):one_cell_setting=',one_cell_setting
       stop
     endif
!
     j_cell_z = j_cell/(n_cell_x*n_cell_y)
     j_cell_y = (j_cell-j_cell_z*n_cell_x*n_cell_y)/n_cell_x
     j_cell_x = j_cell - j_cell_z*n_cell_x*n_cell_y - j_cell_y*n_cell_x
!
   end subroutine get_j_cell_xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_neighbor_cell_list(j_cell, jdx, jdy, jdz, cell_list)
    implicit none
    integer, intent(in)    :: j_cell, jdx, jdy, jdz
    integer, intent(out)   :: cell_list(:)
    integer :: jj, jj_max
    integer :: j_cell_x, j_cell_y, j_cell_z
    integer :: jx, jy, jz, kx, ky, kz, k
!
     if (one_cell_setting) then
       write(*,*)'ERROR(get_j_cell):one_cell_setting=',one_cell_setting
       stop
     endif
!
    cell_list(1)=j_cell
!
    call get_j_cell_xyz(j_cell, j_cell_x, j_cell_y, j_cell_z)
!      ------------> Get the cell adress of (j_cell_x, j_cell_y, j_cell_z)
!
    jj_max=size(cell_list,1)
!
    jj=1
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
          if (k == j_cell) cycle
          jj=jj+1
          if (jj > jj_max) then
            write(*,'(a,2i10)')'ERROR(get_neighbor_cell_list):j_cell,jj=',j_cell, jj
            stop
          endif
          cell_list(jj)=k
        enddo
      enddo
    enddo  

!
   end subroutine get_neighbor_cell_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_cell_index_for_atom(atm_index, j_cell)
!
    use M_qm_domain,       only : atm_position              !(unchanged)
    implicit none
    integer, intent(in)    :: atm_index
    integer, intent(out)   :: j_cell
    real(DOUBLE_PRECISION) :: rx, ry, rz
    integer                :: jx, jy, jz
    logical, parameter     :: check_mode = .true.
!
    if (one_cell_setting) then
      write(*,*)'ERROR(get_j_cell):one_cell_setting=',one_cell_setting
      stop
    endif
!
    rx=atm_position(1,atm_index)
    ry=atm_position(2,atm_index)
    rz=atm_position(3,atm_index)
!
    if (check_mode) then
      if ((rx < 0.0d0) .or. (rx > 1.0)) then
        write(*,*)' ERROR:j,rx=',atm_index,rx
        stop
      endif
      if ((ry < 0.0d0) .or. (ry > 1.0)) then
        write(*,*)' ERROR:j,ry=',atm_index,ry
        stop
      endif
      if ((rz < 0.0d0) .or. (rz > 1.0)) then
        write(*,*)' ERROR:j,ry=',atm_index,rz
        stop
      endif
    endif
!
    jx=int(aint(rx*dble(n_cell_x)))
    jy=int(aint(ry*dble(n_cell_y)))
    jz=int(aint(rz*dble(n_cell_z)))
!
    j_cell=jx+jy*n_cell_x+jz*n_cell_x*n_cell_y
!
  end subroutine get_cell_index_for_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Initial cell booking for the one cell situation
!
   subroutine initial_cell_booking_one_cell
!
     use M_qm_domain,   only : noav !(unchanged)
!     
     implicit none
     integer :: ierr
     logical :: CalledFirst
!
     n_cell=1
     n_cell_x=1
     n_cell_y=1
     n_cell_z=1
     dn_cell_x=0
     dn_cell_y=0
     dn_cell_z=0
!
     if( allocated(cell_for_atom)) then
       CalledFirst = .false.
     else
       CalledFirst = .true.
     endif
!   
     if (CalledFirst) then
       allocate(cell_for_atom(noav), stat=ierr)
       if ( ierr /=  0 ) then
         write(*,*)'ERROR in alloc. : cell_for_atom' 
         stop
       endif   
     endif  
     cell_for_atom(:)=0
!
     if (i_verbose >=0) write(*,'(a,i10)')       '  INFO-CELL(ONE CELL):n_cell=', n_cell
!
     if (CalledFirst) then
       allocate(num_of_atoms_in_cell(0:n_cell-1), stat=ierr)
       if ( ierr /=  0 ) then
        write(*,*)'ERROR in alloc. : num_of_atoms_in_cell'
        stop
       endif   
     endif  
     num_of_atoms_in_cell(0)=noav
!
   end subroutine initial_cell_booking_one_cell
!
end module M_md_dst_cell
