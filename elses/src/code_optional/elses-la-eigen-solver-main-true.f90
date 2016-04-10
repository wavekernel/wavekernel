!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_eig_solver_center
  !
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use eigen_test
  use global_variables
  implicit none
  !include 'mpif.h'
  !
  private
  public :: eig_solver_center, set_density_matrix_mpi, set_pratio_mpi
  !
contains
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! @@ eigen solver center
  !     imode = 1 : standard eigen-value problem with real-symmetrix matrix    (A y = e y )
  !     imode = 2 : generalized eigen-value problem with real-symmetrix matrix (A y = e B y)
  !
  !         mat_a      : input  : input matrix A
  !                      output : eigen vectors ( A(:,k) is the k-th eigen vector)
  !         eig_levels : input  : empty
  !                      output : eigen values e ( e(1) =< e(2) =< ...)
  !         mat_b      : input  : input matrix B
  !                      output : not preserved
  !
  subroutine eig_solver_center(imode, log_unit, SEP_solver_in, GS_transformation_in, &
       blocksize_in, level_low_high, eig_levels, desc_eigenvectors, eigenvectors)
    !
    use mpi
    use elses_mod_md_dat, only : final_iteration
    use M_config, only : config
    use M_ext_matrix_data
    use M_lib_mpi_wrapper
    use wp_setting_m
    use wp_main_aux_m
    use M_wavepacket  ! For testing wavepacket_main_ext().
    implicit none

    integer,          intent(in) :: imode
    character(len=*),       intent(in)               :: SEP_solver_in
    character(len=*),       intent(in)               :: GS_transformation_in
    integer,                intent(in)               :: blocksize_in
!                                                        ! will be stored as 'g_block_size',
!                                                        ! a module variable in 'eigen_test'
    integer,                intent(inout)            :: level_low_high(2)
    integer,                intent(in)               :: log_unit
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    integer,                intent(out)              :: desc_eigenvectors(9)
    real(DOUBLE_PRECISION), intent(out), allocatable :: eigenvectors(:, :)
!
    character(len=1024)                              :: eigen_mpi_scheme_wrk
    character(len=1024)                              :: SEP_solver
    character(len=1024)                              :: GS_transformation
    integer                                          :: mat_size
    integer                                          :: context
    integer                                          :: n_procs_row, n_procs_col, my_proc_row, my_proc_col, ierr
    type(process) :: proc
    type(eigenkernel_sparse_matrix) :: matrix_A, matrix_B
    type(eigenpairs_types_union) :: eigenpairs
    logical :: is_first_call_of_wavepacket = .true.  ! Switch initialization and main loop of wavepacket calculation.
    logical :: is_wavepacket_end = .false.  ! Avoid calling finalization twice.
    type(wp_setting_t), save :: setting
    type(wp_state_t), save :: state
    integer :: i, j
    real(DOUBLE_PRECISION) :: wtime_start, wtime_end
    type(fson_value), pointer :: output
    character(len=128) :: timer_output_filename
    integer, parameter :: timer_output_iunit = 31

    call mpi_barrier(mpi_comm_world, ierr)
    g_mpi_wtime_init = mpi_wtime()
    wtime_start = g_mpi_wtime_init
    output => fson_value_create()
    output%value_type = TYPE_OBJECT

    do i = 1, 2
      if (allocated(matrix_data(i)%element_index)) then
        deallocate(matrix_data(i)%element_index)
      end if
      if (allocated(matrix_data(i)%element_data)) then
        deallocate(matrix_data(i)%element_data)
      end if
    end do
    call set_matrix_data

    if (config%calc%wave_packet%mode == 'on' .and. .not. is_wavepacket_end) then
      if (is_first_call_of_wavepacket) then
        call wavepacket_init(setting, state)
        is_first_call_of_wavepacket = .false.  ! wavepacket_init is called only once.
      else
        call wavepacket_replace_matrix(setting, state)  ! Update result of MD step.
      end if
      call wavepacket_main(setting, state)  ! Compute wavepacket dynamics while atoms are fixed.
      if (final_iteration .or. setting%delta_t * (state%i + 1) >= setting%limit_t) then
        call output_fson_and_destroy(setting, state%output, state%split_files_metadata, &
             state%states, state%structures, state%wtime_total)
        is_wavepacket_end = .true.  ! output_fson_and_destroy is called only once.
      end if
    end if

    wtime_end = mpi_wtime()
    call add_event('main:wavepacket', wtime_end - wtime_start)
    wtime_start = wtime_end
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the work variables
!
    SEP_solver=trim(SEP_solver_in)
    GS_transformation=trim(GS_transformation_in)
!
    mat_size=matrix_data(1)%matrix_size
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Default setting
!
    if (trim(SEP_solver) == 'default') SEP_solver='scalapack'
    if (trim(GS_transformation) == 'default') GS_transformation='scalapack'
    if (blocksize_in /= -1) g_block_size=blocksize_in
    if (level_low_high(1) == -1) level_low_high(1)=1
    if (level_low_high(2) == -1) level_low_high(2)=mat_size
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (log_unit >0) then
      write(log_unit,*)'@@ eig_solver_center'
      write(log_unit,*)'  SEP_solver        =', trim(SEP_solver)
      write(log_unit,*)'  GS_transformation =', trim(GS_transformation)
      write(log_unit,*)'  blocksize         =', g_block_size
      write(log_unit,*)'  lowest level      =', level_low_high(1) ! NOT USED
      write(log_unit,*)'  highest level     =', level_low_high(2) ! NOT USEDkkk
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Rename the solver routine names
!
    if (trim(SEP_solver) == 'scalapack')    SEP_solver='ScaLAPACK'
    if (trim(SEP_solver) == 'EigenExa_sx')  SEP_solver='EigenExa'
    if (trim(SEP_solver) == 'EigenExa_s')   SEP_solver='EigenK'
!
    if (trim(GS_transformation) == 'scalapack') GS_transformation='ScaLAPACK'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Select the solver routine
!
    eigen_mpi_scheme_wrk=''
    if ((trim(SEP_solver) == 'ScaLAPACK') .and. (trim(GS_transformation) == 'ScaLAPACK')) then
      eigen_mpi_scheme_wrk='scalapack'
    else if ((trim(SEP_solver) == 'EigenExa') .and. (trim(GS_transformation) == 'ScaLAPACK')) then
      eigen_mpi_scheme_wrk='eigenexa_scalapack'
    else if ((trim(SEP_solver) == 'ScaLAPACK') .and. (trim(GS_transformation) == 'ELPA')) then
      eigen_mpi_scheme_wrk='scalapack_elpa'
    else if ((trim(SEP_solver) == 'ELPA1') .and. (trim(GS_transformation) == 'ELPA')) then
      eigen_mpi_scheme_wrk='elpa1_elpa'
    else if ((trim(SEP_solver) == 'ELPA2') .and. (trim(GS_transformation) == 'ELPA')) then
      eigen_mpi_scheme_wrk='elpa2_elpa'
    else if ((trim(SEP_solver) == 'EigenExa') .and. (trim(GS_transformation) == 'ELPA')) then
      eigen_mpi_scheme_wrk='eigenexa_elpa'
    else if ((trim(SEP_solver) == 'EigenK') .and. (trim(GS_transformation) == 'ELPA')) then
      eigen_mpi_scheme_wrk='eigenk_elpa'
    endif
!
    if (trim(eigen_mpi_scheme_wrk) == '') then
      write(*,*)'ERROR(eig_solver_center):No solver is chosen '
      write(*,*)' SEP solver = "',trim(SEP_solver), '"'
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    !
    if (log_unit >0) then
      write(log_unit,*)'matrix size of A =', matrix_data(1)%matrix_size
      write(log_unit,*)'matrix size of B =', matrix_data(1)%matrix_size
    end if

    if (imode /= 2) stop ! Standard problems not supported
!
    if (log_unit >0) write(log_unit,*) 'solver: ', trim(eigen_mpi_scheme_wrk)

    call convert_sparse_matrix_data_to_eigenkernel_sparse(matrix_data(1), matrix_A)
    call convert_sparse_matrix_data_to_eigenkernel_sparse(matrix_data(2), matrix_B)
!
    select case (trim(eigen_mpi_scheme_wrk))
    case ('scalapack')
      call eig_solver_center_scalapack(matrix_A, matrix_B, level_low_high, eig_levels, &
           desc_eigenvectors, eigenvectors)
    !case ('eigenexa_scalapack')
    !  call eig_solver_center_eigenexa_scalapack(matrix_A, matrix_B, level_low_high, eig_levels, &
    !       desc_eigenvectors, eigenvectors)
    !case ('scalapack_elpa')
    !  call eig_solver_center_scalapack_elpa(matrix_A, matrix_B, level_low_high, eig_levels, &
    !       desc_eigenvectors, eigenvectors)
    !case ('elpa1_elpa')
    !  call eig_solver_center_elpa1_elpa(matrix_A, matrix_B, level_low_high, eig_levels, &
    !  desc_eigenvectors, eigenvectors)
    !case ('elpa2_elpa')
    !  call eig_solver_center_elpa2_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
    !case ('eigenexa_elpa')
    !  call eig_solver_center_eigenexa_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
    !case ('eigenk_elpa')
    !  call eig_solver_center_eigenk_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
    case default
      stop 'ERROR(eig_solver_center): Unknown solver'
    end select

    wtime_end = mpi_wtime()
    call add_event('main:total', wtime_end - wtime_start)

    if (check_master()) then
      call fson_events_add(output)
      if (trim(config%calc%solver%eigen_mpi%timer_output_filename) == '') then
        timer_output_filename = 'timer_output.json'
      else
        timer_output_filename = trim(config%calc%solver%eigen_mpi%timer_output_filename)
      end if
      open(timer_output_iunit, file=trim(timer_output_filename), status='unknown', position='append')
      call fson_print(timer_output_iunit, output)
      close(timer_output_iunit)
    end if

    if (log_unit >0) then
      write(log_unit,*)'eig_solver_center finished'
    end if
  end subroutine eig_solver_center

  subroutine convert_sparse_matrix_data_to_eigenkernel_sparse(matrix_data_, ek_sparse)
    use M_ext_matrix_data
    type(sparce_matrix_data_real_type) :: matrix_data_
    type(eigenkernel_sparse_matrix), intent(out) :: ek_sparse

    ek_sparse%size = matrix_data_%matrix_size
    ek_sparse%num_non_zeros = matrix_data_%num_of_non_zero_elements
    allocate(ek_sparse%suffix(2, ek_sparse%num_non_zeros), ek_sparse%value(ek_sparse%num_non_zeros))
    ek_sparse%suffix(:, :) = matrix_data_%element_index(:, :)
    ek_sparse%value(:) = matrix_data_%element_data(:)
  end subroutine convert_sparse_matrix_data_to_eigenkernel_sparse

  subroutine eig_solver_center_scalapack(matrix_A, matrix_B, level_low_high, eig_levels, &
       desc_eigenvectors, eigenvectors)
    type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
    type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
    integer,                intent(in)               :: level_low_high(2)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    integer,                intent(out)              :: desc_eigenvectors(9)
    real(DOUBLE_PRECISION), allocatable, intent(out) :: eigenvectors(:, :)

    integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), ierr
    type(process) :: proc
    double precision, allocatable :: A_dist(:, :), B_dist(:, :)
    type(eigenpairs_types_union) :: eigenpairs

    dim = matrix_A%size
    num_output_vectors = level_low_high(2) - level_low_high(1) + 1
    call setup_distribution(proc)
    call setup_distributed_matrix('A', proc, dim, dim, desc_A, A_dist)
    call setup_distributed_matrix('B', proc, dim, dim, desc_B, B_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A, A_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, B_dist)
    call reduce_generalized(dim, A_dist, desc_A, B_dist, desc_B)
    call eigen_solver_scalapack_all(proc, desc_A, A_dist, eigenpairs)
    call recovery_generalized(dim, dim, B_dist, desc_B, eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    desc_eigenvectors(:) = eigenpairs%blacs%desc
    if (allocated(eigenvectors)) then
      deallocate(eigenvectors)
    end if
    allocate(eigenvectors(size(eigenpairs%blacs%Vectors, 1), size(eigenpairs%blacs%Vectors, 2)))
    eigenvectors(:, :) = eigenpairs%blacs%Vectors
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_scalapack

  !subroutine eig_solver_center_eigenexa_scalapack(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
  !  integer,                intent(in)               :: level_low_high(2)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_vectors(:, :)
  !
  !  integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), desc_A_re(9), ierr
  !  type(process) :: proc
  !  double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
  !  type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp
  !
  !  dim = matrix_A%size
  !  num_output_vectors = level_low_high(2) - level_low_high(1) + 1
  !  call setup_distribution(proc)
  !  call setup_distributed_matrix('A', proc, dim, dim, desc_A, matrix_A_dist)
  !  call setup_distributed_matrix('B', proc, dim, dim, desc_B, matrix_B_dist)
  !  call setup_distributed_matrix_for_eigenexa(dim, desc_A_re, matrix_A_redist, eigenpairs_tmp)
  !  call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
  !  call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
  !  call reduce_generalized(dim, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
  !  call pdgemr2d(dim, dim, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
  !  deallocate(matrix_A_dist)
  !  call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, dim, eigenpairs_tmp, 'L')
  !  deallocate(matrix_A_redist)
  !  eigenpairs%type_number = 2
  !  allocate(eigenpairs%blacs%values(dim), stat = ierr)
  !  if (ierr /= 0) then
  !    call terminate('eigen_solver, general_scalapack_eigenexa: allocation failed', ierr)
  !  end if
  !  eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
  !  eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
  !  call setup_distributed_matrix('Eigenvectors', proc, dim, dim, &
  !       eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
  !  call pdgemr2d(dim, dim, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
  !       eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
  !       eigenpairs_tmp%blacs%desc(context_))
  !  call recovery_generalized(dim, dim, matrix_B_dist, desc_B, &
  !       eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
  !  call gather_matrix_part(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
  !       1, level_low_high(1), dim, num_output_vectors, 0, 0, eig_vectors)
  !  call mpi_bcast(eig_vectors, dim * num_output_vectors, mpi_double_precision, 0, mpi_comm_world, ierr)
  !  eig_levels(:) = eigenpairs%blacs%values(:)
  !end subroutine eig_solver_center_eigenexa_scalapack
  !
  !subroutine eig_solver_center_scalapack_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
  !  integer,                intent(in)               :: level_low_high(2)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_vectors(:, :)
  !
  !  integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), desc_A_re(9), ierr
  !  type(process) :: proc
  !  double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
  !  type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp
  !
  !  dim = matrix_A%size
  !  num_output_vectors = level_low_high(2) - level_low_high(1) + 1
  !  call setup_distribution(proc)
  !  call solve_with_general_elpa_scalapack(dim, proc, matrix_A, eigenpairs, matrix_B)
  !  call gather_matrix_part(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
  !       1, level_low_high(1), dim, num_output_vectors, 0, 0, eig_vectors)
  !  call mpi_bcast(eig_vectors, dim * num_output_vectors, mpi_double_precision, 0, mpi_comm_world, ierr)
  !  eig_levels(:) = eigenpairs%blacs%values(:)
  !end subroutine eig_solver_center_scalapack_elpa
  !
  !subroutine eig_solver_center_elpa1_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
  !  integer,                intent(in)               :: level_low_high(2)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_vectors(:, :)
  !
  !  integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), desc_A_re(9), ierr
  !  type(process) :: proc
  !  double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
  !  type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp
  !
  !  dim = matrix_A%size
  !  num_output_vectors = level_low_high(2) - level_low_high(1) + 1
  !  call setup_distribution(proc)
  !  call solve_with_general_elpa1(dim, proc, matrix_A, eigenpairs, matrix_B)
  !  call gather_matrix_part(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
  !       1, level_low_high(1), dim, num_output_vectors, 0, 0, eig_vectors)
  !  call mpi_bcast(eig_vectors, dim * num_output_vectors, mpi_double_precision, 0, mpi_comm_world, ierr)
  !  eig_levels(:) = eigenpairs%blacs%values(:)
  !end subroutine eig_solver_center_elpa1_elpa
  !
  !
  !subroutine eig_solver_center_elpa2_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
  !  integer,                intent(in)               :: level_low_high(2)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_vectors(:, :)
  !
  !  integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), desc_A_re(9), ierr
  !  type(process) :: proc
  !  double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
  !  type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp
  !
  !  dim = matrix_A%size
  !  num_output_vectors = level_low_high(2) - level_low_high(1) + 1
  !  call setup_distribution(proc)
  !  call solve_with_general_elpa2(dim, proc, matrix_A, eigenpairs, matrix_B)
  !  call gather_matrix_part(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
  !       1, level_low_high(1), dim, num_output_vectors, 0, 0, eig_vectors)
  !  call mpi_bcast(eig_vectors, dim * num_output_vectors, mpi_double_precision, 0, mpi_comm_world, ierr)
  !  eig_levels(:) = eigenpairs%blacs%values(:)
  !end subroutine eig_solver_center_elpa2_elpa
  !
  !subroutine eig_solver_center_eigenexa_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
  !  integer,                intent(in)               :: level_low_high(2)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_vectors(:, :)
  !
  !  integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), desc_A_re(9), ierr
  !  type(process) :: proc
  !  double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
  !  type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp
  !
  !  dim = matrix_A%size
  !  num_output_vectors = level_low_high(2) - level_low_high(1) + 1
  !  call setup_distribution(proc)
  !  call solve_with_general_elpa_eigenexa(dim, proc, matrix_A, eigenpairs, matrix_B)
  !  call gather_matrix_part(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
  !       1, level_low_high(1), dim, num_output_vectors, 0, 0, eig_vectors)
  !  call mpi_bcast(eig_vectors, dim * num_output_vectors, mpi_double_precision, 0, mpi_comm_world, ierr)
  !  eig_levels(:) = eigenpairs%blacs%values(:)
  !end subroutine eig_solver_center_eigenexa_elpa
  !
  !subroutine eig_solver_center_eigenk_elpa(matrix_A, matrix_B, level_low_high, eig_levels, eig_vectors)
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_A
  !  type(eigenkernel_sparse_matrix), intent(in)      :: matrix_B
  !  integer,                intent(in)               :: level_low_high(2)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
  !  real(DOUBLE_PRECISION), intent(out)              :: eig_vectors(:, :)
  !
  !  integer :: max_block_size, dim, num_output_vectors, desc_A(9), desc_B(9), desc_A_re(9), ierr
  !  type(process) :: proc
  !  double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
  !  type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp
  !
  !  dim = matrix_A%size
  !  num_output_vectors = level_low_high(2) - level_low_high(1) + 1
  !  call setup_distribution(proc)
  !  call solve_with_general_elpa_eigenk(dim, proc, matrix_A, eigenpairs, matrix_B)
  !  call gather_matrix_part(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, &
  !       1, level_low_high(1), dim, num_output_vectors, 0, 0, eig_vectors)
  !  call mpi_bcast(eig_vectors, dim * num_output_vectors, mpi_double_precision, 0, mpi_comm_world, ierr)
  !  eig_levels(:) = eigenpairs%blacs%values(:)
  !end subroutine eig_solver_center_eigenk_elpa

  subroutine set_pratio_mpi()  ! set pratios to atmp(k, 1)
    use elses_arr_eig_leg, only : atmp, desc_eigenvectors, eigenvectors

    integer :: n_tot_base, i, j, n_procs_row, n_procs_col, my_proc_row, my_proc_col
    real(8) :: sum_power4(desc_eigenvectors(cols_)), sum_power2(desc_eigenvectors(cols_)), pratio
    complex(kind(0d0)) :: elem
    integer :: indxg2p

    n_tot_base = desc_eigenvectors(cols_)

    if (allocated(atmp)) then
      deallocate(atmp)
    endif    
    allocate(atmp(n_tot_base, 1))
    
    call blacs_gridinfo(desc_eigenvectors(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    sum_power4(:) = 0d0
    sum_power2(:) = 0d0
    do j = 1, n_tot_base
      if (indxg2p(j, desc_eigenvectors(block_col_), 0, 0, n_procs_col) == my_proc_col) then
        do i = 1, n_tot_base
          if (indxg2p(i, desc_eigenvectors(block_row_), 0, 0, n_procs_row) == my_proc_row) then
            call pdelget('Self', ' ', elem, eigenvectors, i, j, desc_eigenvectors)
            sum_power4(j) = sum_power4(j) + elem ** 4d0
            sum_power2(j) = sum_power2(j) + elem ** 2d0
          end if
        end do
      end if
    end do
    call dgsum2d(desc_eigenvectors(context_), 'All', ' ', 1, n_tot_base, sum_power4, 1, -1, -1)
    call dgsum2d(desc_eigenvectors(context_), 'All', ' ', 1, n_tot_base, sum_power2, 1, -1, -1)

    do j = 1, n_tot_base
      pratio = (sum_power2(j) ** 2d0) / sum_power4(j)
      atmp(j, 1) = pratio
    end do
  end subroutine set_pratio_mpi

  subroutine set_density_matrix_mpi()
    use M_qm_domain, only : i_verbose, &
         &      dhij, dsij, dbij, dpij, njsd, noav, atm_element, nval, jsv4jsd, &
         &      temp_for_electron, chemical_potential, DOUBLE_PRECISION, ddsij
    use elses_arr_eig_leg, only : atmp, eig2, f_occ, desc_eigenvectors, eigenvectors
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base
    use elses_mod_js4jsv,   only : js4jsv
    use elses_mod_multi,    only : ict4h
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r4
!
    implicit none
    integer :: neig_k
    integer :: jsv2, js2, ja2, jsd1
    integer :: jsv1, js1, ja1, j1, j2, k
    integer :: maxNumOrb, numOrb1, numOrb2
    integer :: ierr
    integer :: desc_atmp(desc_size), local_rows, local_cols, nprow, npcol, myprow, mypcol, prow, pcol
    integer :: indxg2p  ! Function.
!
    real(DOUBLE_PRECISION), allocatable, dimension(:,:)&
            :: f_eigenvectors, f_e_eigenvectors, l_matrix, p_matrix
    real(DOUBLE_PRECISION) :: w0, db_elem, dp_elem
    real(DOUBLE_PRECISION) :: EPSILON=1d-14
!
    if (desc_eigenvectors(rows_) /= n_tot_base .or. desc_eigenvectors(cols_) /= n_tot_base) then
      write(*,*) 'ERROR(set_density_matrix_mpi): Inconsistent matrix dimension'
      stop
    end if
    call blacs_gridinfo(desc_eigenvectors(context_), nprow, npcol, myprow, mypcol)

    dbij(:,:,:,:) = 0.0d0
    dpij(:,:,:,:) = 0.0d0

    do neig_k = 1, n_tot_base
       if(f_occ(neig_k) < EPSILON) exit
    end do
    print *, 'neig_k, n_tot_base', neig_k, n_tot_base
!     ----> upper limit of sum_k
!            ( k = 1, 2, ..., neig_k )

    local_rows = size(eigenvectors, 1)
    local_cols = size(eigenvectors, 2)
    allocate(f_eigenvectors(local_rows, local_cols))
    allocate(f_e_eigenvectors(local_rows, local_cols))
    allocate(l_matrix(local_rows, local_cols))
    allocate(p_matrix(local_rows, local_cols))

    f_eigenvectors(:, :) = eigenvectors(:, :)

    do ja1 = 1, neig_k
      w0 = sqrt(f_occ(ja1))
      call pdscal(n_tot_base, w0, f_eigenvectors, 1, ja1, desc_eigenvectors, 1)
    end do
    do ja1 = neig_k + 1, n_tot_base
      call pdscal(n_tot_base, 0d0, f_eigenvectors, 1, ja1, desc_eigenvectors, 1)
    end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calc. of L-matrix and P-matrix as dbij and dpij
!
!      L(i,j) = sum_k C(i,k) f(k) C(j,k)
!             = sum_k C(i,k) D(k,j) 
!
!      P(i,j) = sum_k C(i,k) f(k) E(k) C(j,k)
!             = sum_k C(i,k) Q(k,j) 
!
    if (i_verbose >= 1) then
      write(6,*)'calc.of L-matrix and P-matrix'
    endif

    f_e_eigenvectors(:, :) = f_eigenvectors(:, :)  ! f_eigenvectors is used as temporary copy of f_eigenvectors.
    call pdgemm('N', 'T', n_tot_base, n_tot_base, n_tot_base, 1d0, &
         f_eigenvectors, 1, 1, desc_eigenvectors, &
         f_e_eigenvectors, 1, 1, desc_eigenvectors, 0d0, &
         l_matrix, 1, 1, desc_eigenvectors)    
    do k = 1, neig_k
      call pdscal(n_tot_base, eig2(k), f_e_eigenvectors, 1, k, desc_eigenvectors, 1)
    end do
    call pdgemm('N', 'T', n_tot_base, n_tot_base, n_tot_base, 1d0, &
         f_eigenvectors, 1, 1, desc_eigenvectors, &
         f_e_eigenvectors, 1, 1, desc_eigenvectors, 0d0, &
         p_matrix, 1, 1, desc_eigenvectors)    
    
    do jsv2 = 1, noav
      js2 = js4jsv(jsv2)
      numOrb2 = nval(atm_element(jsv2))
      do jsd1 = 1, njsd(jsv2, ict4h)
        jsv1 = jsv4jsd(jsd1, jsv2)
        js1 = js4jsv(jsv1)
        numOrb1 = nval(atm_element(jsv1))
        do ja2=1, numOrb2
          do ja1=1, numOrb1
            j1 = js2j(ja1, js1)
            j2 = js2j(ja2, js2)
            prow = indxg2p(j1, desc_eigenvectors(block_row_), 0, desc_eigenvectors(rsrc_), nprow)
            pcol = indxg2p(j2, desc_eigenvectors(block_col_), 0, desc_eigenvectors(csrc_), npcol)
            call pdelget('Self', ' ', db_elem, l_matrix, j1, j2, desc_eigenvectors)
            call pdelget('Self', ' ', dp_elem, p_matrix, j1, j2, desc_eigenvectors)
            if (myprow == prow .and. mypcol == pcol) then
              dbij(ja1, ja2, jsd1, js2) = db_elem
              dpij(ja1, ja2, jsd1, js2) = dp_elem
            end if
          enddo
        enddo
      end do
    enddo
    call mpi_wrapper_allreduce_r4(dbij)
    call mpi_wrapper_allreduce_r4(dpij)
    deallocate(f_eigenvectors, f_e_eigenvectors, l_matrix, p_matrix)
  end subroutine set_density_matrix_mpi
end module M_eig_solver_center
