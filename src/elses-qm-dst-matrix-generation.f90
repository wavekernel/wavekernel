!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_matrix_gene
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
  public qm_dst_mat_gene
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine qm_dst_mat_gene(scheme_mode)
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
        write(lu,*)'@@ qm_dst_mat_gene:scheme mode,n_csc_loop=',trim(scheme_mode),n_csc_loop
      endif
    endif  
!
!   stop 'STOP MANUALLY'
!
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
      if (prc_index == 2) cycle
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
        if (ierr /= 0) then 
          write(*,*)'num_atom_proj=', num_atom_proj
          stop 'Abort:ERROR in alloc (jsv4jsk)'
        endif
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
!       cycle
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
      stop 'STOP MANUALLY'
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
  end subroutine qm_dst_mat_gene
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
end module M_qm_dst_matrix_gene

