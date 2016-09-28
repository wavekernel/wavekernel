!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_solver_gkrylov_dst
!
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_unit                    !(unchanged)
! use M_io_dst_write_log, only : log_file_is_set, log_unit   !(unchanged)
  use M_md_dst,           only : set_dst_final               !(routines)
  use M_qm_domain,        only : c_system, S_is_I            !(unchanged)
!
  implicit none
  logical                     :: dst_micro_mat_is_active
!  
  private
!
! Public routines
  public qm_solver_gkrylov_dst
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine qm_solver_gkrylov_dst(scheme_mode)
!
    use M_config,           only : config !(unchaged
!                                         !  except config%calc%solver%mode_for_large_memory, config%calc%calc_force_mode)
    use elses_mod_orb2,     only : js2j ! (unchanged)
    use M_qm_domain,        only : i_verbose, noav, atm_element, nval, &
&                                    total_electron_number, chemical_potential ! (unchanged)
    use M_qm_projection,    only : proj_get_mat_size, &
&                                  proj_get_list, convert_dm_and_edm, get_interac_list_num_proj_atom ! (routines)
    use M_la_gkrylov_hist,  only : check_proj_list, alloc_for_hist, set_s_inv_e_j_for_hist ! (routines)
    use elses_mod_md_dat,    only : itemd, itemdmx
    use M_la_gkrylov_output, only : calc_eigen_in_gkrylov
!
    use M_la_gkrylov_main,  only : gkrylov_main !(routine)
!
    use M_qm_dst_proj_cell,           only : set_dst_atm_list !(routine)
    use M_qm_dst_proj_cell,           only : dst_atm_list, len_dst_atm_list !(CHANGED in set_dst_atm_list)
    use M_md_dst,           only : myrank !(unchanged)
!
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r1, mpi_wrapper_allreduce_r2, & 
&                                  mpi_wrapper_allreduce_r3, mpi_wrapper_allreduce_r4, &
&                                  mpi_wrapper_allreduce_i2, mpi_wrapper_allreduce_r0, mpi_wrapper_final !(routine)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_barrier_time !(routine)
!
!    use M_qm_dst_global_weight, only : set_global_weight_mat, set_chem_pot_by_global_weight_mat, &
!&                                      plot_ldos_by_global_weight_mat  !(routine)
!
!   use M_qm_dst_global_dens_mat, only : set_global_density_mat !(routine)
!
    use M_qm_domain,        only : atm_position                 !(unchanged)
    use M_qm_domain,        only : atm_force_tb0, atm_force_csc !(CHANGED)
    use elses_mod_sel_sys,  only : r_cut_book
!
    use elses_mod_ene,             only : etb, ecc   !(CHANGED)
!
    use M_qm_domain,               only : e_num_on_atom, e_num_on_basis !(CHANGED)
    use M_qm_domain_dst,           only : global_dens_mat !(unchanged)
!
    use M_qm_geno_dst,   only : set_ham_tb0_and_overlap_dst !(routine)
    use M_qm_dst_global_ham_mat, only : set_ham_tb0_and_overlap_global !(routine)
!
    use M_md_dst,        only : set_dst_final !(routine)
!
    use M_qm_domain_dst, only : global_ham_mat !(unchanged)
    use M_qm_geno_CSC_dst, only : set_hamiltonian_csc_dst !(routine)
    use M_qm_geno_CSC_dst, only : set_atm_force_csc_dst   !(routine)
    use M_qm_geno_CSC_dst, only : set_CSC_parameters_geno_dst !(routine)
!
    use M_wall_clock_time, only : get_system_clock_time   !(routine)
!
    use M_qm_dst_proj_cell,      only : proj_get_mat_size_dst        !(routine)
    use M_qm_dst_proj_cell,      only : proj_get_list_dst            !(routine)
    use M_md_get_distance,       only : get_distance                 !(routine)
    use M_qm_dst_micro_mat,      only : set_size_for_dstm_micro_mat  !(routine)
    use M_qm_dst_micro_mat,      only : set_ham_tb0_and_overlap_dstm !(routine)
    use M_la_gkrylov_main_dstm,  only : gkrylov_main_dstm            !(routine)
    use M_qm_charge_from_wt_dst, only : calc_charge_from_wt_dst      !(routine)
    use M_qm_partial_trace_dstm, only : calc_partial_trace_dstm      !(routine)
    use M_qm_partial_trace_dstm, only : calc_partial_tb0_force_dstm  !(routine)
    use M_qm_partial_trace_dstm, only : set_booking_list_rev1_dstm   !(routine)
    use M_qm_partial_trace_dstm, only : set_interac_list_dstm        !(routine)
    use M_qm_geno_rest_dstm,     only : calc_geno_rest_dstm          !(routine)
!
    use M_cohp_dstm_plot,     only : file_unit_cohp       !(unchanged)
    use M_qm_cohp_dstm,          only : calc_cohp_dstm       !(routine)
!
    use M_qm_chk_allocated,      only : chk_allocated_dst    !(routine)
!
    implicit none
    character(len=32), intent(in) :: scheme_mode
!
    integer :: atm_index, orb_index, mat_size, num_energy_points
!   integer :: atm_index, orb_index, mat_size, num_energy_points, j_src
    integer :: num_atom_proj
    integer :: nss2, nval2, m, n, ierr, i, k
    integer :: prc_index
    integer :: imode
    integer :: kr_dim_max_input, i_kr_hst_str
    integer :: i_show
!
    real(8) :: elec_num, Nelec
!
    integer, allocatable :: jsv4jsk(:)
!   integer, allocatable :: jsk4jsv(:)
    integer, allocatable :: jjkset(:)
    real(8), allocatable :: b(:)
!
    integer :: m_int
    real(8), allocatable :: s_inv_e_j_wrk(:)
    real(8), allocatable :: dm_wrk(:,:)        
!               dm_wrk(:,1) : Density matrix in the compressed format 
!               dm_wrk(:,2) : Energy density matrix in the compress format 
!
    real(8), allocatable :: u_hst_wrk(:,:)
    real(8), allocatable :: v_mat_kr(:,:)
    real(8), allocatable :: eig_wrk(:)
    real(8), allocatable :: u_b_hst_wrk(:)
    real(8), allocatable :: wt_kr_wrk(:)
!
    integer              :: kr_dim_max
!
    integer              :: dst_atm_index
!
    integer, allocatable :: kr_dim_dst(:,:)
!   real(8), allocatable :: s_inv_e_j_dst(:,:,:)
!
    real(8), allocatable :: wt_kr_dst(:,:,:)
    real(8), allocatable :: eig_kr_dst(:,:,:)
!
    real(8) :: cutoff_radius
!
    real(8) :: xmu
!
    logical :: global_weight_mat_is_active
!   logical :: global_density_mat_is_active
!
    integer              :: alloc_size_dst
!
    real(8) :: ene_tmp, ene_tmp_rep
    real(8),  allocatable ::  dst_work_array(:,:)
!
    integer              :: size1
    integer              :: jsv, ja, jsv2, ja2 
    integer              :: n_csc_loop
    integer              :: nval_max
!
    integer              :: size_for_int_atoms         ! Used for DSTM matrix size                  ! DIFFERENT AMONG THEADS!
    integer, allocatable :: booking_list_dstm(:,:)     ! Booking list for DSTM matrix               ! DIFFERENT AMONG THEADS!
    integer, allocatable :: booking_list_rev1_dstm(:)  ! Reverse booking list for DSTM matrix       ! DIFFERENT AMONG THEADS!
    integer, allocatable :: booking_list_dstm_len(:)   ! Length of the booking list for DSTM matrix ! DIFFERENT AMONG THEADS!
    real(8), allocatable :: ham_tot_dstm(:,:,:,:)      ! DSTM matrix for total Hamiltoninan         ! DIFFERENT AMONG THEADS!
    real(8), allocatable :: overlap_dstm(:,:,:,:)      ! DSTM matrix for overlap                    ! DIFFERENT AMONG THEADS!
    real(8), allocatable :: ham_tb0_dstm(:,:,:,:)      ! DSTM matrix for tb0 Hamiltoninan           ! DIFFERENT AMONG THEADS!
    real(8), allocatable :: d_overlap_dstm(:,:,:,:,:)  ! DSTM matrix for derivative of overlap      ! DIFFERENT AMONG THEADS!
    real(8), allocatable :: d_ham_tb0_dstm(:,:,:,:,:)  ! DSTM matrix for derivative of tb0 Ham mat  ! DIFFERENT AMONG THEADS!
!
    logical, parameter :: mpi_time_check = .true.
    real(8) :: time_wrk, time_wrk_previous, time_check ! work variable for measuring the time
!
    logical :: calc_charge
    logical, parameter :: dst_matrices = .false.
!
    real(DOUBLE_PRECISION), allocatable :: atm_energy_dst(:,:)
    real(DOUBLE_PRECISION), allocatable :: local_energy_on_basis(:,:)
    real(DOUBLE_PRECISION), allocatable :: atm_tb0_force_dstm(:,:)
    integer                             :: jsk
    real(8)                             :: memory_size
    real(DOUBLE_PRECISION), allocatable :: atm_wrk_force_omp(:,:,:)
    integer                             :: number_of_omp_threads, id_of_my_omp_thread
    integer                             :: omp_get_num_threads, omp_get_thread_num
    logical                             :: calc_rest_part
    logical, parameter                  :: calc_atm_energy_by_2nd_def = .false.
!   logical, parameter                  :: calc_cohp                  = .false.
    logical                             :: calc_cohp
!
    real(8), allocatable :: e_num_on_atom_dst(:)    ! DIFFERENT VALUES AMONG NODES
    real(8), allocatable :: e_num_on_basis_dst(:,:) ! DIFFERENT VALUES AMONG NODES
!
    integer              :: v_level  ! verbose level
    integer              :: lu       ! log unit
!
    v_level = config%option%verbose_level
    lu      = config%calc%distributed%log_unit
    n_csc_loop = config%calc%genoOption%CSC_max_loop_count
    nval_max=maxval(nval)
    if (allocated(atm_force_tb0))  atm_force_tb0(:,:)=0.0d0
!
    if (allocated(file_unit_cohp)) then
      calc_cohp = .true.
    else
      calc_cohp = .false.
    endif
!
    if (v_level >= 1) then
      if (lu > 0) then
        write(lu,*)'@@ qm_solver_gkrylov_dst:scheme mode,n_csc_loop=',trim(scheme_mode),n_csc_loop
        write(lu,*)'  nval_max=',nval_max
      endif
    endif  
!
    if (config%calc%solver%mode_for_large_memory == 1) then
      config%calc%solver%mode_for_large_memory=0
      if (i_verbose >= 0) then 
       if (lu>0) write(lu,*)'Warning:Now mode_for_large_memory is forced to be zero in DST calculation'
      endif
    endif
!
    if (n_csc_loop == 0) then
      calc_rest_part = .true.
      calc_charge = .false.
      if (config%output%atom_charge%set) calc_charge = .true.
      if (i_verbose >= 0) then
        if (lu>0) write(lu,*)'INFO:calc_charge =', calc_charge
      endif
    else
      calc_rest_part = .false.
      calc_charge = .true.
      if (i_verbose >= 0) then
        if (lu>0) write(lu,*)'INFO:calc_rest_part=',calc_rest_part
      endif
    endif
!
    if (v_level  >= 1) then 
      if (lu >0) then
        write(lu,*)'INFO:calc_atm_energy_by_2nd_def= ', calc_atm_energy_by_2nd_def
        write(lu,*)'INFO:calc_cohp                 = ', calc_cohp
      endif
    endif   
!
    dst_micro_mat_is_active = .true.
    if (v_level >= 1) then
      if (lu >0) then
        write(lu,*)' dst_micro_mat_is_active= ',dst_micro_mat_is_active
      endif
    endif  
!
!
!   global_weight_mat_is_active = .true. 
    global_weight_mat_is_active = .false. 
!
    call get_system_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
!   imode=1
!   call proj_init_end_dst(imode)
    call proj_init_compat ! Routine for compatibility to the old code
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:init_end    =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:init_end    =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
    call set_dst_atm_list
!      ------> setting of dst_atm_list, len_dst_atm_list
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:dst_atm_list =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:dst_atm_list =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
!   if (global_ham_mat) then
!     call set_ham_tb0_and_overlap_global
!   endif  
!
    if (n_csc_loop >=1 ) then
      call set_CSC_parameters_geno_dst
      call set_hamiltonian_csc_dst
    endif  
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc1         =',time_wrk-time_wrk_previous
    if (log_unit > 0 ) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc1         =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
!
!   call check_proj_list !! COMMENTED OUT Jun. 2011, T.Hoshi
!   call alloc_for_hist  !! COMMENTED OUT Jun. 2011, T.Hoshi
!
!   imode=1
!   call set_s_inv_e_j_for_hist(imode) !! COMMENTED OUT Jun. 2011, T.Hoshi
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc2         =',time_wrk-time_wrk_previous
    if (log_unit > 0 ) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc2         =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
    kr_dim_max_input=config%calc%solver%dimension
!
    if (v_level >=1) then
      if (lu > 0) then
        write(lu,'(a,i10)')'INFO:kr_dim_max_input    =', kr_dim_max_input
        write(lu,'(a,i10)')'INFO:maxval(nval)        =', maxval(nval)
        write(lu,'(a,i10)')'INFO:len_dst_atm_list(1) =', len_dst_atm_list(1)
      endif
    endif       
!
    if (.not. allocated(kr_dim_dst)) then
      allocate (kr_dim_dst(maxval(nval),len_dst_atm_list(1)),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (kr_dim_dst)'
      kr_dim_dst(:,:)=kr_dim_max_input
      if (i_verbose >=1) then
        if (lu>0) write(lu,'(a,f20.10)')'INFO:Alloc. of kr_dim_dst           : size [GB]  =', & 
&                                4.0d0*dble(maxval(nval))*dble(len_dst_atm_list(1))/1.0d9
      endif  
    endif
!   
    allocate (wt_kr_dst(kr_dim_max_input, maxval(nval), len_dst_atm_list(1)), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (wt_kr_dst)'
    wt_kr_dst(:,:,:)=0.0d0
    if (i_verbose >=1) then
      if (lu>0) write(lu,'(a,f20.10)')'INFO:Alloc. of wt_kr_dst              : size [GB]  =', & 
&                                8.0d0*dble(kr_dim_max_input)*dble(maxval(nval))*dble(len_dst_atm_list(1))/1.0d9
    endif  
!
    allocate (eig_kr_dst(kr_dim_max_input, maxval(nval), len_dst_atm_list(1)), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (wt_kr_dst)'
    eig_kr_dst(:,:,:)=0.0d0
    if (i_verbose >=1) then
      if (lu>0) write(lu,'(a,f20.10)')'INFO:Alloc. of eig_kr_dst             : size [GB]  =', & 
&                                8.0d0*dble(kr_dim_max_input)*dble(maxval(nval))*dble(len_dst_atm_list(1))/1.0d9
    endif  
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc3         =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc3         =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
    allocate (atm_energy_dst(len_dst_atm_list(1),0:2), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (atm_energy_dst)'
    atm_energy_dst(:,:)=0.0d0
    if (i_verbose >=1) then
      if (lu>0) write(lu,'(a,f20.10)')'INFO:Alloc. of atm_energy_dst         : size [GB]  =', & 
&                                8.0d0*3.0d0*dble(len_dst_atm_list(1))/1.0d9
    endif  

!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc4         =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc4         =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc6         =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc6         =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
!   call set_dst_final
!   stop 'Stop manually:before prc_index loop'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    number_of_omp_threads=1
!$omp parallel default(shared) 
!$  number_of_omp_threads=omp_get_num_threads()
    if (i_verbose >= 1) then 
      if (lu>0) write(lu,*)'number_of_omp_threads=',number_of_omp_threads
    endif
!$omp end parallel
!
    if (trim(config%calc%calc_force_mode) /= "off") then 
      allocate (atm_wrk_force_omp(3,noav,0:number_of_omp_threads-1),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (atm_wrk_force_omp)'
    else   
      if (log_unit > 0) write(log_unit,'(a)')'INFO:array is not allocated: atm_wrk_force_omp  '
    endif  
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc7         =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:msc7         =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
!   call set_dst_final
!   stop 'Stop manually:prc_index'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do prc_index=1,2
!
      if (lu > 0)  write(lu,*)' prc_index = ',prc_index
!
!     call set_dst_final
!     stop 'Stop manually:prc_index'
!
      Nelec=0.0d0 
!
      call get_system_clock_time(time_wrk)
      time_wrk_previous=time_wrk
!
      if (lu > 0)  write(lu,*)'len_dst_atm_list(1)=',len_dst_atm_list(1)
!
      if (v_level >= 1) then 
        if (lu > 0)  write(lu,*)'call chk_allocated_dst:prc_index=', prc_index
        call chk_allocated_dst
      endif  
!
!     call set_dst_final
!     stop 'Stop manually:prc_index:2'
!
!$omp  parallel default(shared) &
!$omp& private (atm_index, orb_index, i_show) &
!$omp& private (elec_num, mat_size, num_atom_proj, nss2, nval2, ierr) &
!$omp& private (jsv4jsk, jjkset, b) &
!$omp& private (s_inv_e_j_wrk, dm_wrk, m_int, u_hst_wrk) &
!$omp& private (kr_dim_max_input, i_kr_hst_str) &
!$omp& private (v_mat_kr, eig_wrk, u_b_hst_wrk, wt_kr_wrk, kr_dim_max) &
!$omp& private (size_for_int_atoms) &
!$omp& private (booking_list_dstm, booking_list_rev1_dstm, booking_list_dstm_len) &
!$omp& private (overlap_dstm, ham_tot_dstm, ham_tb0_dstm) &
!$omp& private (d_overlap_dstm, d_ham_tb0_dstm) &
!$omp& private (local_energy_on_basis) &
!$omp& private (atm_tb0_force_dstm, jsk) &
!$omp& private (memory_size) &
!$omp& private (id_of_my_omp_thread) &
!$omp& reduction (+ : Nelec) 
       id_of_my_omp_thread=0
!$     id_of_my_omp_thread=omp_get_thread_num()
!$     if (number_of_omp_threads /= omp_get_num_threads()) then
!$       write(*,*)'ERROR(Main OMP loop)=', number_of_omp_threads, omp_get_num_threads()
!$       stop
!$     endif
!     write(*,*)'OMP loop start:id,num_theads=',id_of_my_omp_thread, number_of_omp_threads
      if (v_level >= 1) then 
       if (lu > 0) then 
         write(lu,'(a,2i10)')'INFO:OMP_id, P_OMP in Krylov loop =', & 
&                              omp_get_thread_num(), omp_get_num_threads()
       endif
      endif
      if (allocated(atm_wrk_force_omp)) atm_wrk_force_omp(:,:,id_of_my_omp_thread)=0.0d0
!$omp  do schedule(static)
      do dst_atm_index=1,len_dst_atm_list(1)
!
!       write(*,*)'OMP loop : dst_atm_index =',dst_atm_index
!       cycle
!
        atm_index=dst_atm_list(dst_atm_index)
        if ((atm_index < 1) .or. (atm_index > noav)) then
          write(*,*)'ERROR:atm_index=',atm_index
          stop
        endif
        memory_size=0.0d0
!
        kr_dim_max_input=config%calc%solver%dimension
        i_kr_hst_str=config%calc%solver%mode_for_large_memory
!
        i_show=0
        if (i_verbose >= 1) then
          if (atm_index <= 2) then
            i_show=1
          endif  
        endif
!   
        call get_elec_num_for_atom(atm_index, elec_num)
        Nelec=Nelec+elec_num ! summing up of electron numbers
!
        call proj_get_mat_size_dst(dst_atm_index, mat_size, num_atom_proj)
        if (i_show >= 1) then
          if (lu>0) then
            write(lu,*)'atm_index     =',atm_index
            write(lu,*)'mat_size      =',mat_size
            write(lu,*)'num_atom_proj =',num_atom_proj
          endif
        endif  
        nss2=atm_element(atm_index)
        nval2=nval(nss2)
!
        if (i_show >= 1) then
          if (lu>0) then
            write(lu,*)'nss2          =',nss2
            write(lu,*)'nval2         =',nval2
          endif
        endif  
!
        allocate (jsv4jsk(num_atom_proj),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (jsv4jsk)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of jsv4jsk           : size [GB]  =', &
&                                          4.0d0*dble(num_atom_proj)/1.0d9
        endif   
!
!       allocate (jsk4jsv(noav),stat=ierr)
!       if (ierr /= 0) stop 'Abort:ERROR in alloc (jsk4jsv)'
!
        allocate (jjkset(num_atom_proj),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (jjkset)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of jjkset            : size [GB]  =', &
&                                          4.0d0*dble(num_atom_proj)/1.0d9
        endif   
!
        allocate (b(mat_size),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (b)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of b                 : size [GB]  =', & 
&                                          8.0d0*dble(mat_size)/1.0d9
        endif   
!
        allocate (s_inv_e_j_wrk(mat_size),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (s_inv_e_j_wrk)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of s_inv_ej          : size [GB]  =', & 
&                                          8.0d0*dble(mat_size)/1.0d9
        endif   
!
        call proj_get_list_dst(dst_atm_index, jsv4jsk, jjkset)
!           ---> Set jsv4jsk, jjkset
!
        call set_size_for_dstm_micro_mat(atm_index, jsv4jsk, num_atom_proj, size_for_int_atoms)
!           ---> Set size_for_int_atoms
!
        allocate (booking_list_dstm(size_for_int_atoms, num_atom_proj),stat=ierr)
        if( ierr /= 0) stop 'ERROR in alloc (booking_list_dstm)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of booking_list_dstm : size [GB]  =', & 
&                              4.0d0*dble(size_for_int_atoms)*dble(num_atom_proj)/1.0d9
        endif   
!
        allocate (booking_list_rev1_dstm(num_atom_proj),stat=ierr)
        if( ierr /= 0) stop 'ERROR in alloc (booking_list_rev1_dstm)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of booking_list_dstm : size [GB]  =', & 
&                              4.0d0*dble(num_atom_proj)/1.0d9
        endif   
!
        allocate (booking_list_dstm_len(num_atom_proj),stat=ierr)
        if( ierr /= 0) stop 'ERROR in alloc (booking_list_dstm_len)'
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of booking_list_dstm_len : size [GB]  =', & 
&                              4.0d0*dble(num_atom_proj)/1.0d9
        endif   
!
        allocate (ham_tb0_dstm(nval_max, nval_max, size_for_int_atoms, num_atom_proj), stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (ham_tb0_dstm)'
        memory_size=memory_size+8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of ham_tb0_dstm          : size [GB]  =', & 
&                            8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)/1.0d9
        endif   
!
        allocate (ham_tot_dstm(nval_max, nval_max, size_for_int_atoms, num_atom_proj), stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (ham_tot_dstm)'
        memory_size=memory_size+8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of ham_tot_dstm          : size [GB]  =', & 
&                            8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)/1.0d9
        endif   
!
        allocate (d_ham_tb0_dstm(3, nval_max, nval_max, size_for_int_atoms, num_atom_proj), stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (d_ham_tot_dstm)'
        memory_size=memory_size+8.0d0*3.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of ham_tb0_dstm          : size [GB]  =', & 
&                            8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)/1.0d9
        endif   
!
        if (S_is_I .eqv. .false.) then
          allocate (overlap_dstm(nval_max, nval_max, size_for_int_atoms, num_atom_proj), stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in alloc (overlap_dstm)'
          memory_size=memory_size+8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)
          if (i_show >= 1) then
            if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of overlap_dstm          : size [GB]  =', & 
&                            8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)/1.0d9
          endif   
!
          allocate (d_overlap_dstm(3, nval_max, nval_max, size_for_int_atoms, num_atom_proj), stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in alloc (d_overlap_dstm)'
          memory_size=memory_size+8.0d0*3.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)
          if (i_show >= 1) then
            if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of d_overlap_dstm        : size [GB]  =', & 
&                      3.0d0*8.0d0*dble(nval_max)*dble(nval_max)*dble(size_for_int_atoms)*dble(num_atom_proj)/1.0d9
          endif   
        endif  
!
        if (i_show >= 1) then
          if (lu>0) write(lu,'(a,f20.10)')'memory_size [GB]   =', memory_size/1.0d9
        endif  
!
        call set_ham_tb0_and_overlap_dstm(atm_index, jsv4jsk, num_atom_proj, &
&                booking_list_dstm, booking_list_dstm_len,  overlap_dstm, ham_tb0_dstm, d_overlap_dstm, d_ham_tb0_dstm)
!         ------> Set booking_list_dstm, booking_list_dstm, 
!                     overlap_dstm, overlap_dstm, overlap_dstm, overlap_dstm
!
        call set_booking_list_rev1_dstm(booking_list_rev1_dstm, booking_list_dstm, booking_list_dstm_len, num_atom_proj)
!         ------> Set booking_list_rev1_dstm(num_atom_proj)
!
        m_int=0 ! Dummy value
        call set_interac_list_dstm(jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset, m_int)
!         ------> Set m_int
!
        ham_tot_dstm(:,:,:,:)=ham_tb0_dstm(:,:,:,:)
!        
        if ((m_int <= 0) .or. (m_int > 1000000)) then
          write(*,*)'ERROR(qm_solver_gkrylov_dst):m_int=',m_int
          stop
        endif
!
        if (calc_rest_part) then
          if (allocated(atm_wrk_force_omp)) then 
            call calc_geno_rest_dstm(atm_index, jsv4jsk, num_atom_proj,  &
&                 booking_list_dstm, booking_list_dstm_len,            & 
&                 atm_energy_dst(dst_atm_index,0), atm_wrk_force_omp(:,:,id_of_my_omp_thread))
          else  
            call calc_geno_rest_dstm(atm_index, jsv4jsk, num_atom_proj,  &
&                 booking_list_dstm, booking_list_dstm_len,            & 
&                 atm_energy_dst(dst_atm_index,0))
          endif  
        endif
!
        allocate (dm_wrk(m_int,2),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (dm_wrk)'
!
        allocate (u_hst_wrk(m_int,kr_dim_max_input),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (u_hst_wrk)'
!
        allocate (v_mat_kr(kr_dim_max_input,kr_dim_max_input),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (v_mat_kr)'
!
        allocate (eig_wrk(kr_dim_max_input),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (eig_wrk)'
!
        allocate (u_b_hst_wrk(kr_dim_max_input),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (eig_wrk)'
!
        allocate (wt_kr_wrk(kr_dim_max_input),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (wt_kr_wrk)'
!
        allocate (local_energy_on_basis(nval2,2),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (local_energy_on_basis)'
        local_energy_on_basis(:,:)=0.0d0
!
        if (trim(config%calc%calc_force_mode) /= "off") then 
          allocate (atm_tb0_force_dstm(3,num_atom_proj), stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in alloc (atm_tb0_force_dstm)'
          if (i_show >= 1) then
            if (lu>0) write(lu,'(a,f20.10)')'INFO(in OMP loop):Alloc. of atm_tb0_force_dstm    : size [GB]  =', & 
&                        3.0d0*8.0d0*dble(num_atom_proj)/1.0d9
          endif   
        endif  
!
        do orb_index=1,nval2
          if (i_show >= 1) then
            if (lu>0) write(lu,*)"atom=",atm_index,"orbit=",orb_index 
          endif  
!         j_src=js2j(orb_index,atm_index)
          b(:)=0.0d0
          b(orb_index)=1.0d0
          s_inv_e_j_wrk(:)=b(:)
          dm_wrk(:,:)=0.0d0
          u_hst_wrk(:,:)=0.0d0
          v_mat_kr(:,:)=0.0d0
          eig_wrk(:)=0.0d0
          u_b_hst_wrk(:)=0.0d0
          wt_kr_wrk(:)=0.0d0
          if (i_kr_hst_str == 1) then
            if (prc_index == 1) then
              kr_dim_dst(orb_index, dst_atm_index)=kr_dim_max_input
              wt_kr_dst(:,orb_index,dst_atm_index)=0.0d0
              eig_kr_dst(:,orb_index,dst_atm_index)=0.0d0
            else
              wt_kr_wrk(1:kr_dim_max_input)=wt_kr_dst(1:kr_dim_max_input,orb_index,dst_atm_index)
              kr_dim_max=kr_dim_dst(orb_index, dst_atm_index)
            endif  
          endif  
          call gkrylov_main_dstm(dst_atm_index, atm_index, orb_index, b, jsv4jsk, jjkset, prc_index, scheme_mode, & 
&                    s_inv_e_j_wrk,dm_wrk,u_hst_wrk, v_mat_kr, eig_wrk, u_b_hst_wrk, wt_kr_wrk, kr_dim_max, &
&                    booking_list_dstm, booking_list_dstm_len, overlap_dstm, ham_tot_dstm, id_of_my_omp_thread )
!
          if (prc_index == 1) then
            wt_kr_dst(1:kr_dim_max_input,orb_index,dst_atm_index) =wt_kr_wrk(1:kr_dim_max_input)
            eig_kr_dst(1:kr_dim_max_input,orb_index,dst_atm_index)=eig_wrk(1:kr_dim_max_input)
            kr_dim_dst(orb_index, dst_atm_index)=kr_dim_max
          endif  
!         if (allocated(s_inv_e_j_dst)) then
!           s_inv_e_j_dst(1:mat_size,orb_index,dst_atm_index)=s_inv_e_j_wrk(1:mat_size)
!         endif
!
          if (prc_index == 2) then
             call calc_partial_trace_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&                   dm_wrk(:,1), ham_tb0_dstm, local_energy_on_basis(orb_index,1), & 
&                        jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset)
             if (i_show >=1) then 
               if (lu>0) write(lu,*)'local energy on basis 1=',local_energy_on_basis(orb_index,1)
             endif
!
             if (calc_atm_energy_by_2nd_def) then
               call calc_partial_trace_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&                   dm_wrk(:,2), overlap_dstm, local_energy_on_basis(orb_index,2), & 
&                        jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset)
               if (i_show >=1) then 
                 if (lu>0) write(lu,*)'local energy on basis 2=',local_energy_on_basis(orb_index,2)
               endif
             endif  
!
             if (calc_cohp) then
               call calc_cohp_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&                   dm_wrk(:,1), ham_tb0_dstm, local_energy_on_basis(orb_index,1), & 
&                        jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset, id_of_my_omp_thread)
               if (i_show >=1) then 
                  if (lu>0) write(lu,*)'local energy on basis 3=',local_energy_on_basis(orb_index,1)
               endif
             endif  
!
             if (allocated(atm_tb0_force_dstm)) then
               atm_tb0_force_dstm(:,:)=0.0d0
               call calc_partial_tb0_force_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&                  dm_wrk(:,1), dm_wrk(:,2), d_ham_tb0_dstm, d_overlap_dstm,                 &
&                  atm_tb0_force_dstm, jsv4jsk, booking_list_dstm, booking_list_rev1_dstm, booking_list_dstm_len, jjkset)
               if (allocated(atm_wrk_force_omp)) then
                 do jsk=1,num_atom_proj
                   atm_wrk_force_omp(:,jsv4jsk(jsk),id_of_my_omp_thread)=atm_wrk_force_omp(:,jsv4jsk(jsk),id_of_my_omp_thread) &
&                                                                     +atm_tb0_force_dstm(:,jsk)
                 enddo
               endif  
             endif  
!
           endif  
        enddo ! loop end for orb_index
!
        if (prc_index == 2) then
!         write(*,*)'atm_energy1=',dst_atm_index,sum(local_energy_on_basis(:,1))
!         write(*,*)'atm_energy2=',dst_atm_index,sum(local_energy_on_basis(:,2))
          atm_energy_dst(dst_atm_index, 1) = sum(local_energy_on_basis(1:nval2,1))
          if (calc_atm_energy_by_2nd_def) atm_energy_dst(dst_atm_index, 2) = sum(local_energy_on_basis(1:nval2,2))
        endif  
!
        deallocate (dm_wrk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (dm_wrk)'
!
        deallocate (u_hst_wrk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (u_hst_wrk)'
!
        deallocate (v_mat_kr,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (v_mat_kr)'
!
        deallocate (eig_wrk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (eig_wrk)'
!
        deallocate (u_b_hst_wrk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (u_b_hst_wrk)'
!
        deallocate (wt_kr_wrk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (wt_kr_wrk)'
!
        deallocate (local_energy_on_basis, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (local_energy_on_basis)'
!
        if (allocated(atm_tb0_force_dstm)) then
          deallocate (atm_tb0_force_dstm, stat=ierr)
          if (ierr /= 0) stop 'Abort:ERROR in dealloc (atm_tb0_force_dstm)'
        endif  
!
        deallocate(jsv4jsk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate(jjkset,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate(b, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate(s_inv_e_j_wrk,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
        deallocate (booking_list_dstm, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (booking_list_dstm)'
!
        deallocate (booking_list_rev1_dstm, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (booking_list_rev1_dstm)'
!
        deallocate (booking_list_dstm_len, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (booking_list_dstm_len)'
!
        deallocate (overlap_dstm,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (overlap_dstm)'
!
        deallocate (ham_tb0_dstm,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (ham_tb0_dstm)'
!
        deallocate (ham_tot_dstm,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (ham_tot_dstm)'
!
        deallocate (d_overlap_dstm,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (d_overlap_dstm)'
!
        deallocate (d_ham_tb0_dstm,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (d_ham_tb0_dstm)'
!
      enddo ! loop end for atm_index 
!
!$omp end do
!$omp end parallel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (mpi_time_check) then
      call mpi_wrapper_barrier_time(time_check)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:barrier   =',time_check
      if (log_unit > 0) write(log_unit,'(a,f20.10)') 'TIME:qm_solver_gkrylov_dst:barrier   =',time_check
    endif
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:OMP loop    =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5,i4)') & 
&        'TIME:qm_solver_gkrylov_dst:MAIN loop    =',time_wrk-time_wrk_previous, prc_index
    time_wrk_previous=time_wrk
!
    call mpi_wrapper_allreduce_r0(Nelec)
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:MPI comm    =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:MPI comm    =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (config%calc%use_integer_elec_num) then
      if (log_unit > 0) write(log_unit,'(a,f40.20)') &
&           'INFO:integer_elec_num corrction (befor) : Nelec=',Nelec
       Nelec=anint(Nelec)
      if (log_unit > 0) write(log_unit,'(a,f40.20)') &
&           'INFO:integer_elec_num corrction (after) : Nelec=',Nelec
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (prc_index == 2) then ! (calculation of energy)
!
      call get_system_clock_time(time_wrk)
      time_wrk_previous=time_wrk
!
      ene_tmp=0.0d0
      ene_tmp_rep=0.0d0
!$omp  parallel default(shared) &
!$omp& private (dst_atm_index)  &
!$omp& reduction (+ : ene_tmp, ene_tmp_rep)
!$omp  do schedule(static)
      do dst_atm_index=1,len_dst_atm_list(1)
        ene_tmp     = ene_tmp     + atm_energy_dst(dst_atm_index,1)
        ene_tmp_rep = ene_tmp_rep + atm_energy_dst(dst_atm_index,0)
      enddo
!$omp end do
!$omp end parallel
!
      call get_system_clock_time(time_wrk)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:ENE_OMP     =',time_wrk-time_wrk_previous
      if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:ENE_OMP     =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      call mpi_wrapper_allreduce_r0(ene_tmp)
      call mpi_wrapper_allreduce_r0(ene_tmp_rep)
!
      call get_system_clock_time(time_wrk)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:ENE_MPI     =',time_wrk-time_wrk_previous
      if (log_unit > 0 ) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:ENE_MPI     =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      etb=ene_tmp
      ecc=ene_tmp_rep
      if (i_verbose >=0) then
        if (lu>0) then
          write(lu,*)'ETB0 (from atm_energy_dstm)=',etb
          write(lu,*)'Erest(from atm_energy_dstm)=',ecc
        endif
      endif  
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (prc_index == 2) then ! (calculation of force)
!
      call get_system_clock_time(time_wrk)
      time_wrk_previous=time_wrk
!
      if (allocated(atm_force_tb0)) atm_force_tb0(:,:)=0.0d0
      if (allocated(atm_wrk_force_omp) .and. allocated(atm_force_tb0)) then
!$omp  parallel default(shared) &
!$omp& private (atm_index, k)  
!$omp  do schedule(static)
       do atm_index=1,noav
        do k=0, number_of_omp_threads-1
          atm_force_tb0(:,atm_index) = atm_force_tb0(:,atm_index) + atm_wrk_force_omp(:,atm_index,k)
        enddo
       enddo
!$omp end do
!$omp end parallel
      endif
!
      call get_system_clock_time(time_wrk)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:FORCE_a     =',time_wrk-time_wrk_previous
      if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:FORCE_a     =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      if (allocated(atm_force_tb0)) call mpi_wrapper_allreduce_r2(atm_force_tb0)
!
      call get_system_clock_time(time_wrk)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:FORCE_OMP    =',time_wrk-time_wrk_previous
      if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:FORCE_OMP   =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      if (i_verbose >= 1) then
        if (allocated(atm_force_tb0)) then
          do atm_index=1,min(noav,11)
            if (lu>0) write(lu,'(a,i10,3f20.10)')'tb0 force=',atm_index, atm_force_tb0(1:3,atm_index)
          enddo
        endif  
      endif
!
!     stop 'Stop manually'
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (prc_index == 1) then ! (calculation of chemical potential)
      if (i_verbose > 1) then 
        if (lu>0) write(lu,*)'Nelec=',Nelec
      endif
      if ( dabs( total_electron_number - Nelec)/dabs(Nelec) >= 1.0d-10 ) then
         write(*,*)'ERROR in Nelec : Nelec=',Nelec, total_electron_number
         stop
      endif   
!
      call get_system_clock_time(time_wrk)
      time_wrk_previous=time_wrk
!
      call set_chemical_potential_dst(Nelec, dst_atm_list, len_dst_atm_list, & 
&                                     kr_dim_dst, wt_kr_dst, eig_kr_dst, xmu)
      chemical_potential=xmu
!
      call get_system_clock_time(time_wrk)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:chem pot    =',time_wrk-time_wrk_previous
      if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:chem pot    =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
!     call set_dst_final
!     stop 'Stop manually'
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   calc_charge = .true.
    if (calc_charge) then ! (calculation of charge)
!
      call get_system_clock_time(time_wrk)
      time_wrk_previous=time_wrk
!
      size1=size(e_num_on_basis,1)
!
      allocate (e_num_on_basis_dst(size1,len_dst_atm_list(1)), stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (e_num_on_basis_dst)'
      e_num_on_basis_dst(:,:)=0.0d0
!
      call calc_charge_from_wt_dst(dst_atm_list, len_dst_atm_list, & 
&                                   kr_dim_dst, wt_kr_dst, eig_kr_dst, xmu, e_num_on_basis_dst)
!       -----> e_num_on_basis_dst
!
!     call set_dst_final
!     stop 'Stop manually:after call calc_charge_from_wt_dst'
!
      allocate (e_num_on_atom_dst(len_dst_atm_list(1)), stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (e_num_on_atom_dst)'
!
      do dst_atm_index=1,len_dst_atm_list(1)
        e_num_on_atom_dst(dst_atm_index)=sum(e_num_on_basis_dst(:,dst_atm_index))
      enddo
!
      e_num_on_atom(:)=0.0d0
      e_num_on_basis(:,:)=0.0d0
      do dst_atm_index=1,len_dst_atm_list(1)
        atm_index=dst_atm_list(dst_atm_index)
        e_num_on_basis(:,atm_index)=e_num_on_basis_dst(:,dst_atm_index)
        e_num_on_atom(atm_index)=e_num_on_atom_dst(dst_atm_index)
      enddo  
!
      call mpi_wrapper_allreduce_r1(e_num_on_atom)
      call mpi_wrapper_allreduce_r2(e_num_on_basis)
!       -----> e_num_on_atom, e_num_on_basis
!
      deallocate (e_num_on_atom_dst, stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc (e_num_on_atom_dst)'
!
      deallocate (e_num_on_basis_dst, stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc (e_num_on_basis_dst)'
!
      call get_system_clock_time(time_wrk)
!     write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:charge calc    =',time_wrk-time_wrk_previous
      if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:charge calc    =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    enddo ! loop end for prc_index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call get_system_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
    deallocate (wt_kr_dst, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (wt_kr_dst)'
!
    deallocate (eig_kr_dst, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (eig_kr_dst)'
!
    deallocate (kr_dim_dst,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (kr_dim_dst)'
!
!   if (allocated(s_inv_e_j_dst)) then
!     deallocate (s_inv_e_j_dst, stat=ierr )
!     if (ierr /= 0) stop 'Abort:ERROR in dealloc (s_inv_e_j_dst)'
!   endif  
!
    deallocate (atm_energy_dst, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (atm_energy_dst)'
!
    if (allocated(atm_wrk_force_omp)) then
      deallocate (atm_wrk_force_omp, stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc (atm_wrk_force_omp)'
    endif  
!
!
    call get_system_clock_time(time_wrk)
!   write(*,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:deallc          =',time_wrk-time_wrk_previous
    if (log_unit > 0) write(log_unit,'(a,f20.5)')'TIME:qm_solver_gkrylov_dst:deallc          =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
  end subroutine qm_solver_gkrylov_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Setting chemical potential in the distirbuted system
!
  subroutine set_chemical_potential_dst(Nelec, dst_atm_list, len_dst_atm_list, & 
&                                       kr_dim_dst, wt_kr_dst, eig_kr_dst, xmu)
!
    use M_config,            only : config !(unchaged)
    use M_lib_phys_const,    only : ev4au !(unchaged)
    use M_md_dst,            only : myrank, nprocs !(unchanged)
    use M_qm_domain,         only : i_verbose, noav, atm_element, nval, & 
&                                   temp_for_electron !(unchanged)
    use M_lib_math_func,     only : Fermi_Dirac_Func !(function)
    use M_lib_mpi_wrapper,   only : mpi_wrapper_allreduce_minmax_r1, &
&                                   mpi_wrapper_allreduce_r0 !(routine)
!
    implicit none
    real(8),          intent(in)  :: Nelec
    integer,          intent(in)  :: dst_atm_list(:)
    integer,          intent(in)  :: len_dst_atm_list(:)
    integer,          intent(in)  :: kr_dim_dst(:,:)
    real(8),          intent(in)  :: wt_kr_dst(:,:,:)
    real(8),          intent(in)  :: eig_kr_dst(:,:,:)
    real(8),          intent(out) :: xmu  ! Chemical potential (OUTPUT)
    integer :: dst_atm_index, orb_index, atm_index, kr_dim_max
!
    real(8) :: max_energy(1)
    real(8) :: min_energy(1)
!
    real(8) :: ddemin, ddemax
    real(8) :: xtemp, xbeta
!
    integer :: ipe, npe
    integer :: iloop, nloopmax, irec
    real(8) :: xmu0, N_elec_tmp
    real(8) :: xmumn, xmumx, rRd, rEd, x, xexp, err
    integer :: lu
!
    lu= config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      if (lu>0) write(lu,*)'@@ set_chemical_potential'
    endif
!  
    nloopmax=100
    xtemp=temp_for_electron
!        ----> Temperature in Fermi Distribution
    xbeta=1.0d0/xtemp
!
    if (i_verbose >= 1) then
      if (lu>0) write(lu,*)'xtemp [au,eV]=',xtemp,xtemp*ev4au
    endif  
!
    max_energy(1)=maxval(eig_kr_dst)
    min_energy(1)=minval(eig_kr_dst)
!
    if (i_verbose >=1) then 
      if (lu>0) write(lu,'(a,i10,2f30.20)')'myrank, max,min (loc)=',myrank, max_energy(1), min_energy(1)
    endif
!
    call mpi_wrapper_allreduce_minmax_r1(max_energy,'max')
    call mpi_wrapper_allreduce_minmax_r1(min_energy,'min')
!
    if (i_verbose >=1) then 
      if (lu>0) write(lu,'(a,i10,2f30.20)')'myrank, max,min (glo)=',myrank, max_energy(1), min_energy(1)
    endif
!
    ddemin=min_energy(1)
    ddemax=max_energy(1)
!
    if (i_verbose >= 1) then
      if (lu>0) write(lu,*)' ddemin [au,eV]=',ddemin,ddemin*ev4au
      if (lu>0) write(lu,*)' ddemax [au,eV]=',ddemax,ddemax*ev4au
      if (dabs(ddemin-ddemax) .lt. 1.0d-10) then
        write(*,*)'ERROR!:SET_KRY_CHEM_POT'
        stop
      endif
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ The bisection method
!
      xmu0=ddemin
      xmumn=ddemin
      xmumx=ddemax
!
      do iloop=1,nloopmax
!
       if (iloop == nloopmax) then
         write(*,*)'Bisection was not converged ! '
         stop
       endif
!   
       N_elec_tmp =0.0d0
!
!      OMP-NOTE : Shared scalors : xmu0,xbeta
!$omp  parallel &
!$omp& default(shared) & 
!$omp& private(ipe,npe) &
!$omp& private(dst_atm_index, atm_index, orb_index) &
!$omp& private(irec, rRd, rEd, x, xexp) &
!$omp& reduction(+ : N_elec_tmp)
       ipe=0
       npe=0
!      ipe=omp_get_thread_num()+1
!      npe=omp_get_num_threads()
!      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
        do dst_atm_index=1,len_dst_atm_list(1)
          atm_index=dst_atm_list(dst_atm_index)
!         if ((atm_index < 1) .or. (atm_index > noav)) then
!           write(*,*)'ERROR:atm_index=',atm_index
!           stop
!         endif
          do orb_index=1, nval(atm_element(atm_index))
            do irec=kr_dim_dst(orb_index, dst_atm_index),1,-1
              rRd= wt_kr_dst(irec, orb_index, dst_atm_index)
              rEd=eig_kr_dst(irec, orb_index, dst_atm_index)
              x=xbeta*(rEd-xmu0)
              xexp=Fermi_Dirac_Func(x)
              N_elec_tmp =  N_elec_tmp + 2.0d0*xexp*rRd
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel 
!
        call mpi_wrapper_allreduce_r0(N_elec_tmp)
!
        err=(Nelec-N_elec_tmp)/Nelec
        if (i_verbose >=1) then 
          if (lu>0) write(lu,*)'Bisec.',iloop,xmu0,err
        endif
!
        if (abs(err) < 1.0d-12) exit
!
        if (err >  0.0d0) then
          xmumn=xmu0 
          xmumx=xmumx
          xmu0=(xmumn+xmumx)*0.5d0
        else
          xmumn=xmumn
          xmumx=xmu0 
          xmu0=(xmumn+xmumx)*0.5d0
        endif
!
      enddo ! loop by iloop
!
      if (i_verbose >=0) then
        if (lu>0) write(lu,'(a,i20,2f30.20)')'chem_pot [au,eV]=',iloop,xmu0,xmu0*ev4au
      endif  
      xmu=xmu0
!
!    write(*,*)' Stop Manually'
!    stop
!
  end subroutine set_chemical_potential_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Get the allocation size for density matrix and energy density matrix
!
  subroutine get_max_booking_list_length(max_len_book_list)
!
    use M_qm_domain,     only : noav !(unchanged)
    use M_qm_dst_proj_cell,        only : dst_atm_list, len_dst_atm_list !(unchanged)
    use M_qm_projection, only : get_interac_list_num_proj_atom !(routine)
    implicit none
    integer, intent(out) :: max_len_book_list
    integer :: ierr
    integer :: atm_index, dst_atm_index, m_int
    integer, allocatable :: int_list(:)
!
    allocate ( int_list(len_dst_atm_list(1)) , stat=ierr )
    if (ierr /= 0) stop 'Abort:error in get_alloc_size'
!
    do dst_atm_index=1,len_dst_atm_list(1)
      atm_index=dst_atm_list(dst_atm_index)
!     cutoff_radius=r_cut_book
!     call get_length_of_list_dst(cutoff_radius,length_of_list)
!           ---> get number of atom in the interaction list : length_of_list
      int_list(dst_atm_index)=0
    enddo
!
    max_len_book_list=min( maxval(int_list)*2, noav )
!     ----> two is the tolerance factor
!
  end subroutine get_max_booking_list_length
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ subroutine for obtaining the valence elec. num. for the given atom
!
  subroutine get_elec_num_for_atom(atm_index,elec_num)
!
    use elses_mod_val_elec,     only : val_elec_atm ! (unchanged)
    use M_qm_domain,            only : atm_element  ! (unchanged)
!
    implicit none
    integer,          intent(in)  :: atm_index
    real(8),          intent(out) :: elec_num
    integer :: nss, size1, size2
!
    size1=size(atm_element,1)
    size2=size(val_elec_atm,1)
!
    if ((atm_index <= 0) .or. (atm_index > size1)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'atm_index=',atm_index
      stop
    endif   
!
    nss=atm_element(atm_index)
!
    if ((nss <= 0) .or. (nss > size2)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'nss=',nss
      stop
    endif   
!
    elec_num=val_elec_atm(nss)
!
    if ((elec_num <= 1.0d-10) .or. (elec_num >= 1000.0d0)) then
      write(*,*)'ERROR:get_elec_num_for_atom'
      write(*,*)'  elec_num=',elec_num
      stop
    endif   
!
  end subroutine get_elec_num_for_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Convert the density matrix and energy density matrix 
!        into dm_dst, edm_dst
!
  subroutine convert_dm_and_edm_dst(atm_index, dst_atm_index, orb_index, m_int, & 
&                               dm_loc, edm_loc, dm_dst, edm_dst, jsv4jsk, jsk4jsv, jjkset)
!    
    use M_qm_domain,   only : jsv4jsd, noav, njsd, nval, atm_element !(unchanged)
!
    implicit none
    integer, intent(in)      :: atm_index, dst_atm_index, orb_index, m_int
    real(8), intent(in)      ::  dm_loc(:) ! density matrix for a given basis
    real(8), intent(in)      :: edm_loc(:) ! energy density matrix for a given basis
    integer, intent(in)      :: jsv4jsk(:), jsk4jsv(:), jjkset(:) ! info. of the projection 
    real(8), intent(out)     ::  dm_dst(:,:,:) ! density matrix for a given basis
    real(8), intent(out)     :: edm_dst(:,:,:) ! energy density matrix for a given basis
!
    integer :: jsv2, ja2, ict4h, jj, jsd1, jsv1, jsk1, jjkset1, ja1, jjk1
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set the local variables
!
    jsv2=atm_index
    ja2=orb_index
    ict4h=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Trivial checking for the input quantities (non-essential)
!
    if (size(dm_loc,1) /= m_int) then
       write(*,*)'ERROR:convert_dm_and_edm:size(dm_loc,1)= ',size(dm_loc,1)
       stop
    endif
!
    if (size(edm_loc,1) /= m_int) then
       write(*,*)'ERROR:convert_dm_and_edm:size(edm_loc,1)= ',size(edm_loc,1)
       stop
    endif
!
    if ((jsv2 <=0 ).or. (jsv2 > noav)) then
       write(*,*)'ERROR:convert_dm_and_edm'
       write(*,*)' atm_index, orb_index=',jsv2, ja2
       stop
    endif   
!
    if ((ja2 <=0 ).or. (ja2 > nval(atm_element(jsv2)))) then
       write(*,*)'ERROR:convert_dm_and_edm'
       write(*,*)' atm_index, orb_index=',jsv2, ja2
       stop
    endif   
!
    jj=0
    do jsd1=1,njsd(jsv2,ict4h)
      jsv1=jsv4jsd(jsd1,jsv2)
      jsk1=jsk4jsv(jsv1)
      if (jsk1 == 0) then
        write(*,*)'ERROR(convert_dm_and_edm):jsk1=',jsk1
        stop
      else
        jjkset1=jjkset(jsk1)
        do ja1=1,nval(atm_element(jsv1))
          jjk1=jjkset1+ja1
          jj=jj+1
          if (jj > m_int) then
            write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
            write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
            stop
          endif
           dm_dst(ja1,ja2,jsd1) =  dm_loc(jj)
          edm_dst(ja1,ja2,jsd1) = edm_loc(jj)
        enddo
       endif  
    enddo   
!
    if (m_int /= jj) then
       write(*,*)'ERROR(convert_dm_and_edm):jj,m_int=',jj,m_int
       stop
    endif   
!
  end subroutine convert_dm_and_edm_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_init_compat ! Routine for compatibility to the old code
!
   call elses_set_param_ctl_kr
!
  end subroutine proj_init_compat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_solver_gkrylov_dst
