!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_scalapack_all
  use time, only : get_wclock_time
  use distribute_matrix, only : &
       process, get_local_cols, gather_matrix, allgather_row_wise, &
       desc_size, desc_type_, context_, rows_, cols_, block_row_, block_col_, &
       rsrc_, csrc_, local_rows_
  use eigenpairs_types, only : eigenpairs_types_union
  implicit none

  private
  public :: eigen_solver_scalapack_all

contains

  subroutine eigen_solver_scalapack_all(proc, desc_A, A, eigenpairs)
    include 'mpif.h'

    type(process) :: proc
    integer, intent(in) :: desc_A(9)
    double precision, intent(in) :: A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: ierr, info
    integer :: dim, work_size, iwork_size, diag_size, subdiag_size
    integer :: eigenvectors_local_cols

    character(len = 1) :: uplo, side

    double precision, allocatable :: diag_local(:), subdiag_local(:)
    double precision, allocatable :: subdiag_global(:)
    double precision, allocatable :: tau(:), work(:), work_print(:)
    integer, allocatable :: iwork(:)

    ! Time
    integer, parameter :: n_intervals = 7
    integer :: i
    double precision :: t_intervals(n_intervals)
    double precision :: t_init, t_pdsytrd, t_pdsytrd_end, t_pdstedc, &
         t_pdstedc_end, t_pdormtr_end, t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'init   ', &
         'pdsytrd', 'gather1', 'pdstedc', 'pdormtr', 'finish ', 'total  '/)

    ! Functions
    integer :: numroc

    call get_wclock_time(t_init)

    eigenpairs%type_number = 2

    dim = desc_A(rows_)
    diag_size = numroc(dim, desc_A(block_col_), proc%my_proc_col, 0, proc%n_procs_col)
    subdiag_size = numroc(dim - 1, desc_A(block_col_), proc%my_proc_col, 0, proc%n_procs_col)
    work_size = max(desc_A(block_row_) * (desc_A(local_rows_) + 1), 3 * desc_A(block_row_))

    allocate(diag_local(diag_size))
    allocate(subdiag_local(subdiag_size))
    allocate(tau(diag_size))
    allocate(work(work_size))
    allocate(work_print(desc_A(block_row_)))

    call get_wclock_time(t_pdsytrd)

    uplo = 'l'
    call pdsytrd(uplo, dim, A, 1, 1, desc_A, diag_local, subdiag_local, tau, work, work_size, info)
    deallocate(work)
    if (proc%my_rank == 0) then
       print *, 'info(pdsytrd): ', info
    end if

    call get_wclock_time(t_pdsytrd_end)

    ! Diagonal elements of the tridiagonal matrix Initially
    allocate(eigenpairs%blacs%values(dim))

    allocate(subdiag_global(dim - 1))

    call allgather_row_wise(diag_local, proc%context, desc_A(block_col_), &
         eigenpairs%blacs%values)
    call allgather_row_wise(subdiag_local, proc%context, desc_A(block_col_), &
         subdiag_global)

    call descinit(eigenpairs%blacs%desc, dim, dim, desc_A(block_row_), &
         desc_A(block_col_), 0, 0, proc%context, desc_A(local_rows_), info)
    eigenvectors_local_cols = get_local_cols(proc, eigenpairs%blacs%desc)
    allocate(eigenpairs%blacs%Vectors( &
         1 : eigenpairs%blacs%desc(local_rows_), &
         1 : eigenvectors_local_cols))
    eigenpairs%blacs%Vectors(:, :) = 0.0

    work_size = 6 * dim + 2 * eigenpairs%blacs%desc(local_rows_) * &
         eigenvectors_local_cols
    iwork_size = 2 + 7 * dim + 8 * proc%n_procs_col
    allocate(work(work_size))
    allocate(iwork(iwork_size))

    call get_wclock_time(t_pdstedc)

    call pdstedc('i', dim, eigenpairs%blacs%values, subdiag_global, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         work, work_size, iwork, iwork_size, info)
    if (proc%my_rank == 0) then
       print *, 'info(pdstedc): ', info
    end if

    call get_wclock_time(t_pdstedc_end)

    side = 'l'
    work_size = work_size_for_pdormtr(side, uplo, dim, dim, 1, 1, 1, 1, desc_A, eigenpairs%blacs%desc)
    deallocate(work)
    allocate(work(work_size))

    call pdormtr(side, uplo, 'n', dim, dim, A, 1, 1, desc_A, tau, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, work, work_size, info)
    if (proc%my_rank == 0) then
       print *, 'info(pdormtr): ', info
    end if

    call get_wclock_time(t_pdormtr_end)

    ! call eigentest_pdlaprnt(dim, dim, eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, 0, 0, 'Eigenvectors', 6, work_print)

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_pdsytrd - t_init
    t_intervals(2) = t_pdsytrd_end - t_pdsytrd
    t_intervals(3) = t_pdstedc - t_pdsytrd_end
    t_intervals(4) = t_pdstedc_end - t_pdstedc
    t_intervals(5) = t_pdormtr_end - t_pdstedc_end
    t_intervals(6) = t_all_end - t_pdormtr_end
    t_intervals(7) = t_all_end - t_init

    if (proc%my_rank == 0) then
       call MPI_Reduce(MPI_IN_PLACE, t_intervals, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
       print *, 'Elapsed time (sec)'
       do i = 1, n_intervals
          print *, ' ', interval_names(i), ':', t_intervals(i)
       end do
    else
       call MPI_Reduce(t_intervals, 0, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    end if
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
end module solver_scalapack_all
