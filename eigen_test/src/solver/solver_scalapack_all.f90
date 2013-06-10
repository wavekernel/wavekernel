!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_scalapack_all
  !
  use time, only : get_wclock_time !(routine)
  use processes, only : layout_procs
  use distribute_matrix, only : gather_matrix, copy_global_dense_matrix_to_local, allgather_row_wise
  !
  implicit none

  type conf_distribution
    integer :: dim, my_rank, n_procs, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: block_size, n_local_row, n_local_col
  end type conf_distribution
  !
  private
  public :: conf_distribution, setup_distribution, distribute_dense, eigen_solver_scalapack_all
  !
contains
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setup_distribution(dim, conf)
    integer, intent(in) :: dim
    type(conf_distribution), intent(out) :: conf
    integer :: info

    ! Functions
    integer :: numroc, iceil

    conf%dim = dim

    call blacs_pinfo(conf%my_rank, conf%n_procs)
    call layout_procs(conf%n_procs, conf%n_procs_row, conf%n_procs_col)
    call blacs_get(-1, 0, conf%context)
    call blacs_gridinit(conf%context, 'R', conf%n_procs_row, conf%n_procs_col)
    call blacs_gridinfo(conf%context, conf%n_procs_row, conf%n_procs_col, &
         conf%my_proc_row, conf%my_proc_col)

    conf%block_size = min(32, dim / max(conf%n_procs_row, conf%n_procs_col))

    if (conf%my_rank == 0) then
       print '( "block size:", I7 )', conf%block_size
       print '( "procs: ", I5, " x ", I5, " (", I5, ")" )', &
            conf%n_procs_row, conf%n_procs_col, conf%n_procs
    end if

    if (conf%my_proc_row >= conf%n_procs_row .or. conf%my_proc_col >= conf%n_procs_col) then
       call blacs_exit(0)
       stop 'out of process grid'
    end if

    conf%n_local_row = max(1, numroc(dim, conf%block_size, &
         conf%my_proc_row, 0, conf%n_procs_row))
    conf%n_local_col = max(1, numroc(dim, conf%block_size, &
         conf%my_proc_col, 0, conf%n_procs_col))
  end subroutine setup_distribution

  subroutine distribute_dense(conf, mat_in, desc, mat)
    type(conf_distribution), intent(in) :: conf
    real(kind(1.d0)), intent(in) :: mat_in(:, :)
    integer, intent(out) :: desc(9)
    real(kind(1.d0)), intent(out), allocatable :: mat(:, :)

    integer :: info

    call descinit(desc, conf%dim, conf%dim, conf%block_size, conf%block_size, &
         0, 0, conf%context, conf%n_local_row, info)
    allocate(mat(1 : conf%n_local_row, 1 : conf%n_local_col))

    call copy_global_dense_matrix_to_local(mat_in, desc, mat)
  end subroutine distribute_dense

  subroutine eigen_solver_scalapack_all(conf, desc_A, A, eigen_level, eigenvectors_global)
    implicit none
    include 'mpif.h'

    type(conf_distribution) :: conf
    integer, intent(in) :: desc_A(9)
    real(kind(1.d0)), intent(in) :: A(:, :)
    real(kind(1.d0)), intent(out) :: eigen_level(:), eigenvectors_global(:, :)

    integer :: ierr, info
    integer :: work_size, iwork_size, diag_size, subdiag_size
    integer :: desc_Eigenvectors(9)

    character(len = 1) :: uplo, compz, side, trans

    real(kind(1.d0)), allocatable :: Eigenvectors(:, :)
    real(kind(1.d0)), allocatable :: diag_local(:), subdiag_local(:)
    real(kind(1.d0)), allocatable :: diag_global(:), subdiag_global(:)
    real(kind(1.d0)), allocatable :: tau(:), work(:), work_print(:)
    integer, allocatable :: iwork(:)

    ! Time
    integer, parameter :: n_intervals = 7
    integer :: i
    real(kind(1.d0)) :: t_intervals(n_intervals)
    real(kind(1.d0)) :: t_init, t_pdsytrd, t_pdsytrd_end, t_pdstedc, &
         t_pdstedc_end, t_pdormtr_end, t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'init   ', &
         'pdsytrd', 'gather1', 'pdstedc', 'pdormtr', 'gather2', 'total  '/)

    ! Functions
    integer :: numroc, iceil

    call get_wclock_time(t_init)

    diag_size = numroc(conf%dim, conf%block_size, conf%my_proc_col, 0, conf%n_procs_col)
    subdiag_size = numroc(conf%dim - 1, conf%block_size, conf%my_proc_col, 0, conf%n_procs_col)
    work_size = max(conf%block_size * (conf%n_local_row + 1), 3 * conf%block_size)

    allocate(diag_local(diag_size))
    allocate(subdiag_local(subdiag_size))
    allocate(tau(diag_size))
    allocate(work(work_size))
    allocate(work_print(conf%block_size))

    call get_wclock_time(t_pdsytrd)

    uplo = 'l'
    call pdsytrd(uplo, conf%dim, A, 1, 1, desc_A, diag_local, subdiag_local, tau, work, work_size, info)
    deallocate(work)
    if (conf%my_rank == 0) then
       print *, 'info(pdsytrd): ', info
    end if

    call get_wclock_time(t_pdsytrd_end)

    allocate(diag_global(conf%dim))
    allocate(subdiag_global(conf%dim - 1))

    call allgather_row_wise(diag_local, conf%context, conf%block_size, diag_global)
    call allgather_row_wise(subdiag_local, conf%context, conf%block_size, subdiag_global)

    call descinit(desc_Eigenvectors, conf%dim, conf%dim, conf%block_size, conf%block_size, &
         0, 0, conf%context, conf%n_local_row, info)

    allocate(Eigenvectors(1 : conf%n_local_row, 1 : conf%n_local_col))
    Eigenvectors(:, :) = 0.0

    work_size = 6 * conf%dim + 2 * conf%n_local_row * conf%n_local_col
    iwork_size = 2 + 7 * conf%dim + 8 * conf%n_procs_col
    allocate(work(work_size))
    allocate(iwork(iwork_size))

    call get_wclock_time(t_pdstedc)

    call pdstedc('i', conf%dim, diag_global, subdiag_global, Eigenvectors, 1, 1, &
         desc_Eigenvectors, work, work_size, iwork, iwork_size, info)
    if (conf%my_rank == 0) then
       print *, 'info(pdstedc): ', info
    end if

    call get_wclock_time(t_pdstedc_end)

    eigen_level(:) = diag_global(:)

    side = 'l'
    work_size = work_size_for_pdormtr(side, uplo, conf%dim, conf%dim, 1, 1, 1, 1, desc_A, desc_Eigenvectors)
    deallocate(work)
    allocate(work(work_size))

    call pdormtr(side, uplo, 'n', conf%dim, conf%dim, A, 1, 1, desc_A, tau, &
         Eigenvectors, 1, 1, desc_Eigenvectors, work, work_size, info)
    if (conf%my_rank == 0) then
       print *, 'info(pdormtr): ', info
    end if

    call get_wclock_time(t_pdormtr_end)

    ! call pdlaprnt(dim, dim, Eigenvectors, 1, 1, desc_Eigenvectors, 0, 0, 'Eigenvectors', 6, work_print)


    call gather_matrix(Eigenvectors, desc_Eigenvectors, 0, 0, eigenvectors_global)

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_pdsytrd - t_init
    t_intervals(2) = t_pdsytrd_end - t_pdsytrd
    t_intervals(3) = t_pdstedc - t_pdsytrd_end
    t_intervals(4) = t_pdstedc_end - t_pdstedc
    t_intervals(5) = t_pdormtr_end - t_pdstedc_end
    t_intervals(6) = t_all_end - t_pdormtr_end
    t_intervals(7) = t_all_end - t_init

    if (conf%my_rank == 0) then
       call MPI_Reduce(MPI_IN_PLACE, t_intervals, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
       print *, 'Elapse time (sec)'
       do i = 1, n_intervals
          print *, ' ', interval_names(i), ':', t_intervals(i)
       end do
    else
       call MPI_Reduce(t_intervals, 0, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    end if

    call blacs_exit(0)

  end subroutine eigen_solver_scalapack_all

  integer function work_size_for_pdormtr(side, uplo, m, n, ia, ja, ic, jc, desc_A, desc_C) result(size)
    character(len = 1), intent(in) :: side, uplo
    integer, intent(in) :: m, n, ia, ja, ic, jc
    integer, intent(in) :: desc_A(9), desc_C(9)
    integer :: context, n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: mb_a, nb_a, mb_c, nb_c
    integer :: iaa, jaa, icc, jcc
    integer :: lcmq, npa0 , mpc0, nqc0
    integer :: iroffa, icoffa, iarow, iroffc, icoffc, icrow, iccol
    integer :: mi, ni
    logical :: is_upper, is_left

    integer :: numroc, indxg2p, ilcm
    logical :: lsame

    is_upper = lsame(uplo, 'U')
    is_left = lsame(side, 'L')

    mb_a = desc_A(5)
    nb_a = desc_A(6)
    mb_c = desc_C(5)
    nb_c = desc_C(6)

    context = desc_A(2)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (is_upper) then
       iaa = ia
       jaa = ja + 1
       icc = ic
       jcc = jc
    else
       iaa = ia + 1
       jaa = ja
       if (is_left) then
          icc = ic + 1
          jcc = jc
       else
          icc = ic
          jcc = jc + 1
       end if
    end if

    if (is_left) then
       mi = m - 1
       ni = n
    else
       mi = m
       ni = n - 1
    end if

    iroffc = mod(icc - 1, mb_c)
    icoffc = mod(jcc - 1, nb_c)
    icrow = indxg2p(icc, mb_c, my_proc_row, desc_C(7), n_procs_row)
    iccol = indxg2p(jcc, nb_c, my_proc_col, desc_C(8), n_procs_col)

    mpc0 = numroc(mi + iroffc, mb_c, my_proc_row, icrow, n_procs_row)
    nqc0 = numroc(ni + icoffc, nb_c, my_proc_col, iccol, n_procs_col)

    if (is_left) then
       size = max((nb_a * (nb_a - 1)) / 2, (nqc0 + mpc0) * nb_a) + nb_a * nb_a
    else
       iroffa = mod(iaa - 1, mb_a)
       icoffa = mod(jaa - 1, nb_a)
       iarow = indxg2p(iaa, mb_a, my_proc_row, desc_A(7), n_procs_row)
       npa0 = numroc(ni + iroffa, mb_a, my_proc_row, iarow, n_procs_row)
       lcmq = ilcm(n_procs_row, n_procs_col) / n_procs_col
       size = max((nb_a * (nb_a - 1)) / 2, &
            (nqc0 + max(npa0 + numroc(numroc(ni + icoffc, nb_a, 0, 0, n_procs_col), &
            nb_a, 0, 0, lcmq), mpc0)) * nb_a) + &
            nb_a * nb_a
    end if
  end function work_size_for_pdormtr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module solver_scalapack_all
