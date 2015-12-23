!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_eig_solver_center
  !
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use eigen_test
  implicit none
  include 'mpif.h'
  !
  private
  public :: eig_solver_center
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
&                blocksize_in, level_low_high_in, mat_a, eig_levels, mat_b)
    !
    use mpi
    use M_config, only : config
    implicit none

    integer,          intent(in) :: imode
    character(len=*),       intent(in)               :: SEP_solver_in
    character(len=*),       intent(in)               :: GS_transformation_in
    integer,                intent(in)               :: blocksize_in
!                                                        ! will be stored as 'g_block_size',
!                                                        ! a module variable in 'eigen_test'
    integer,                intent(in)               :: level_low_high_in(2)
    integer,                intent(in)               :: log_unit
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout), optional  :: mat_b(:,:)
!
    character(len=1024)                              :: eigen_mpi_scheme_wrk
    character(len=1024)                              :: SEP_solver
    character(len=1024)                              :: GS_transformation
    integer                                          :: level_low_high(2)
    integer                                          :: mat_size
    integer                                          :: context
    integer                                          :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    type(process) :: proc
    type(sparse_mat) :: matrix_A, matrix_B
    type(eigenpairs_types_union) :: eigenpairs
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the work variables
!
    SEP_solver=trim(SEP_solver_in)
    GS_transformation=trim(GS_transformation_in)
    level_low_high(1:2)=level_low_high_in(1:2)
!
    mat_size=size(mat_a,1)
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
!     write(log_unit,*)'  lowest level      =', level_low_high(1) ! NOT USED
!     write(log_unit,*)'  highest level     =', level_low_high(2) ! NOT USED
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
    if (imode == 2) then
      if (.not. present(mat_b)) then
        write(*,*)'ERROR(eig_solver_center):No matrix B is missing imode=',imode
        stop
      endif
    endif
    !
    if (log_unit >0) write(log_unit,*)'matrix size of A =', size(mat_a, 1), size(mat_a,2)
    !
    if (present(mat_b)) then
      if (log_unit >0) write(log_unit,*)'matrix size of B =', size(mat_b, 1), size(mat_b,2)
    endif

    if (imode /= 2) stop ! Standard problems not supported
!
    if (log_unit >0) write(log_unit,*) 'solver: ', trim(eigen_mpi_scheme_wrk)
!
    select case (trim(eigen_mpi_scheme_wrk))
    case ('scalapack')
      call eig_solver_center_scalapack(mat_a, eig_levels, mat_b)
    case ('eigenexa_scalapack')
      call eig_solver_center_eigenexa_scalapack(mat_a, eig_levels, mat_b)
    case ('scalapack_elpa')
      call eig_solver_center_scalapack_elpa(mat_a, eig_levels, mat_b)
    case ('elpa1_elpa')
      call eig_solver_center_elpa1_elpa(mat_a, eig_levels, mat_b)
    case ('elpa2_elpa')
      call eig_solver_center_elpa2_elpa(mat_a, eig_levels, mat_b)
    case ('eigenexa_elpa')
      call eig_solver_center_eigenexa_elpa(mat_a, eig_levels, mat_b)
    case ('eigenk_elpa')
      call eig_solver_center_eigenk_elpa(mat_a, eig_levels, mat_b)
    case default
      stop 'ERROR(eig_solver_center): Unknown solver'
    end select
    !
    !  stop 'stop manually'
    !
  end subroutine eig_solver_center

  ! Convert a local dense matrix into a sparse one.
  ! Used to match format of input matrices with EigenTest routines
  subroutine dense_to_sparse(dim, dense, sparse)
    integer, intent(in) :: dim
    double precision, intent(in) :: dense(dim, dim)
    type(sparse_mat), intent(out) :: sparse

    integer :: i, j, num_non_zeros, value_index, ierr

    sparse%size = dim
    num_non_zeros = 0
    do j = 1, dim
      do i = 1, dim
        if (dense(i, j) /= 0.0d0) then
          num_non_zeros = num_non_zeros + 1
        end if
      end do
    end do
    sparse%num_non_zeros = num_non_zeros
    allocate(sparse%suffix(2, num_non_zeros), sparse%value(num_non_zeros), stat = ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR(eig_solver_center):allocation failed, ierr=', ierr
    end if

    value_index = 1
    do j = 1, dim
      do i = 1, dim
        if (dense(i, j) /= 0.0d0) then
          sparse%value(value_index) = dense(i, j)
          sparse%suffix(1, value_index) = i
          sparse%suffix(2, value_index) = j
          value_index = value_index + 1
        end if
      end do
    end do
  end subroutine dense_to_sparse

  subroutine eig_solver_center_scalapack(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size
    integer :: dim, desc_A(9), desc_B(9), ierr
    type(process) :: proc
    double precision, allocatable :: A_dist(:, :), B_dist(:, :)
    type(eigenpairs_types_union) :: eigenpairs

    dim = size(mat_a, 1)
    call setup_distribution(proc)
    call setup_distributed_matrix('A', proc, dim, dim, desc_A, A_dist)
    call setup_distributed_matrix('B', proc, dim, dim, desc_B, B_dist)
    call distribute_global_dense_matrix(mat_a, desc_A, A_dist)
    call distribute_global_dense_matrix(mat_b, desc_B, B_dist)
    call reduce_generalized(dim, A_dist, desc_A, B_dist, desc_B)
    call eigen_solver_scalapack_all(proc, desc_A, A_dist, eigenpairs)
    call recovery_generalized(dim, dim, B_dist, desc_B, eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, dim * dim, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_scalapack

  subroutine eig_solver_center_eigenexa_scalapack(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size, n, desc_A(9), desc_B(9), desc_A_re(9), ierr
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp

    n = size(mat_a, 1)
    call setup_distribution(proc)
    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
    call distribute_global_dense_matrix(mat_a, desc_A, matrix_A_dist)
    call distribute_global_dense_matrix(mat_b, desc_B, matrix_B_dist)
    call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
    call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
    deallocate(matrix_A_dist)
    call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')
    deallocate(matrix_A_redist)
    eigenpairs%type_number = 2
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_scalapack_eigenexa: allocation failed', ierr)
    end if
    eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
    eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         eigenpairs_tmp%blacs%desc(context_))
    call recovery_generalized(n, n, matrix_B_dist, desc_B, &
         eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, n * n, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_eigenexa_scalapack

  subroutine eig_solver_center_scalapack_elpa(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size, n, desc_A(9), desc_B(9), desc_A_re(9), ierr
    type(process) :: proc
    type(sparse_mat) :: sparse_mat_a, sparse_mat_b
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp

    n = size(mat_a, 1)
    call dense_to_sparse(n, mat_a, sparse_mat_a)
    call dense_to_sparse(n, mat_b, sparse_mat_b)
    call setup_distribution(proc)
    call solve_with_general_elpa_scalapack(n, proc, sparse_mat_a, eigenpairs, sparse_mat_b)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, n * n, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_scalapack_elpa

  subroutine eig_solver_center_elpa1_elpa(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size, n, desc_A(9), desc_B(9), desc_A_re(9), ierr
    type(process) :: proc
    type(sparse_mat) :: sparse_mat_a, sparse_mat_b
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp

    n = size(mat_a, 1)
    call dense_to_sparse(n, mat_a, sparse_mat_a)
    call dense_to_sparse(n, mat_b, sparse_mat_b)
    call setup_distribution(proc)
    call solve_with_general_elpa1(n, proc, sparse_mat_a, eigenpairs, sparse_mat_b)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, n * n, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_elpa1_elpa

  subroutine eig_solver_center_elpa2_elpa(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size, n, desc_A(9), desc_B(9), desc_A_re(9), ierr
    type(process) :: proc
    type(sparse_mat) :: sparse_mat_a, sparse_mat_b
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp

    n = size(mat_a, 1)
    call dense_to_sparse(n, mat_a, sparse_mat_a)
    call dense_to_sparse(n, mat_b, sparse_mat_b)
    call setup_distribution(proc)
    call solve_with_general_elpa2(n, proc, sparse_mat_a, eigenpairs, sparse_mat_b)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, n * n, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_elpa2_elpa

  subroutine eig_solver_center_eigenexa_elpa(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size, n, desc_A(9), desc_B(9), desc_A_re(9), ierr
    type(process) :: proc
    type(sparse_mat) :: sparse_mat_a, sparse_mat_b
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp

    n = size(mat_a, 1)
    call dense_to_sparse(n, mat_a, sparse_mat_a)
    call dense_to_sparse(n, mat_b, sparse_mat_b)
    call setup_distribution(proc)
    call solve_with_general_elpa_eigenexa(n, proc, sparse_mat_a, eigenpairs, sparse_mat_b)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, n * n, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_eigenexa_elpa

  subroutine eig_solver_center_eigenk_elpa(mat_a, eig_levels, mat_b)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_a(:,:)
    real(DOUBLE_PRECISION), intent(out)              :: eig_levels(:)
    real(DOUBLE_PRECISION), intent(inout)            :: mat_b(:,:)

    integer :: max_block_size, n, desc_A(9), desc_B(9), desc_A_re(9), ierr
    type(process) :: proc
    type(sparse_mat) :: sparse_mat_a, sparse_mat_b
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs, eigenpairs_tmp

    n = size(mat_a, 1)
    call dense_to_sparse(n, mat_a, sparse_mat_a)
    call dense_to_sparse(n, mat_b, sparse_mat_b)
    call setup_distribution(proc)
    call solve_with_general_elpa_eigenk(n, proc, sparse_mat_a, eigenpairs, sparse_mat_b)
    call gather_matrix(eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, 0, 0, mat_a)
    call mpi_bcast(mat_a, n * n, mpi_double_precision, 0, mpi_comm_world, ierr)
    eig_levels(:) = eigenpairs%blacs%values(:)
  end subroutine eig_solver_center_eigenk_elpa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
end module M_eig_solver_center
