!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_scalapack_all
  !
  implicit none
  !
  private
  public :: eigen_solver_scalapack_all
  !
contains
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eigen_solver_scalapack_all(mat,eigen_level)
    !
    use time, only : get_wclock_time !(routine)
    use processes, only : layout_procs
    use distribute_matrix, only : gather_matrix, copy_global_dense_matrix_to_local, allgather_row_wise
    implicit none
    include 'mpif.h'

    real(kind(1.d0)), intent(inout) :: mat(:,:)   ! ( n x n ) matrix
    real(kind(1.d0)), intent(out)   :: eigen_level(:)

    integer :: my_rank, n_procs, block_size, n_procs_row, n_procs_col, context
    integer :: my_proc_row, my_proc_col

    integer :: ierr, info
    integer :: dim, n_local_row, n_local_col, work_size, iwork_size
    integer :: diag_size, subdiag_size
    integer :: desc_A(9), desc_Eigenvectors(9)

    character(len = 1) :: uplo, compz, side, trans

    real(kind(1.d0)), allocatable :: A(:,:)
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

    dim = size(mat, 1)

    call blacs_pinfo(my_rank, n_procs)
    call layout_procs(n_procs, n_procs_row, n_procs_col)
    call blacs_get(-1, 0, context)
    call blacs_gridinit(context, 'R', n_procs_row, n_procs_col)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    block_size = min(32, dim / max(n_procs_row, n_procs_col))

    if (my_rank == 0) then
       print '( "block size:", I7 )', block_size
       print '( "procs: ", I5, " x ", I5, " (", I5, ")" )', n_procs_row, n_procs_col, n_procs
    end if

    if (my_proc_row >= n_procs_row .or. my_proc_col >= n_procs_col) then
       call blacs_exit(0)
       stop 'out of process grid'
    end if

    n_local_row = max(1, numroc(dim, block_size, my_proc_row, 0, n_procs_row))
    n_local_col = max(1, numroc(dim, block_size, my_proc_col, 0, n_procs_col))
    call descinit(desc_A, dim, dim, block_size, block_size, 0, 0, context, n_local_row, info)

    allocate(A(1 : n_local_row, 1 : n_local_col))

    call copy_global_dense_matrix_to_local(mat, desc_A, A)

    diag_size = numroc(dim, block_size, my_proc_col, 0, n_procs_col)
    subdiag_size = numroc(dim - 1, block_size, my_proc_col, 0, n_procs_col)
    work_size = max(block_size * (n_local_row + 1), 3 * block_size)

    allocate(diag_local(diag_size))
    allocate(subdiag_local(subdiag_size))
    allocate(tau(diag_size))
    allocate(work(work_size))
    allocate(work_print(block_size))

    call get_wclock_time(t_pdsytrd)

    uplo = 'l'
    call pdsytrd(uplo, dim, A, 1, 1, desc_A, diag_local, subdiag_local, tau, work, work_size, info)
    deallocate(work)
    if (my_rank == 0) then
       print *, 'info(pdsytrd): ', info
    end if

    call get_wclock_time(t_pdsytrd_end)

    allocate(diag_global(dim))
    allocate(subdiag_global(dim - 1))

    call allgather_row_wise(diag_local, context, block_size, diag_global)
    call allgather_row_wise(subdiag_local, context, block_size, subdiag_global)

    call descinit(desc_Eigenvectors, dim, dim, block_size, block_size, 0, 0, context, n_local_row, info)

    allocate(Eigenvectors(1 : n_local_row, 1 : n_local_col))
    Eigenvectors(:, :) = 0.0

    work_size = 6 * dim + 2 * n_local_row * n_local_col
    iwork_size = 2 + 7 * dim + 8 * n_procs_col
    allocate(work(work_size))
    allocate(iwork(iwork_size))

    call get_wclock_time(t_pdstedc)

    call pdstedc('i', dim, diag_global, subdiag_global, Eigenvectors, 1, 1, &
         desc_Eigenvectors, work, work_size, iwork, iwork_size, info)
    if (my_rank == 0) then
       print *, 'info(pdstedc): ', info
    end if

    call get_wclock_time(t_pdstedc_end)

    eigen_level(:) = diag_global(:)

    side = 'l'
    work_size = work_size_for_pdormtr(side, uplo, dim, dim, 1, 1, 1, 1, desc_A, desc_Eigenvectors)
    deallocate(work)
    allocate(work(work_size))

    call pdormtr(side, uplo, 'n', dim, dim, A, 1, 1, desc_A, tau, &
         Eigenvectors, 1, 1, desc_Eigenvectors, work, work_size, info)
    if (my_rank == 0) then
       print *, 'info(pdormtr): ', info
    end if

    call get_wclock_time(t_pdormtr_end)

    ! call pdlaprnt(dim, dim, Eigenvectors, 1, 1, desc_Eigenvectors, 0, 0, 'Eigenvectors', 6, work_print)


    call gather_matrix(Eigenvectors, desc_Eigenvectors, 0, 0, mat)

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_pdsytrd - t_init
    t_intervals(2) = t_pdsytrd_end - t_pdsytrd
    t_intervals(3) = t_pdstedc - t_pdsytrd_end
    t_intervals(4) = t_pdstedc_end - t_pdstedc
    t_intervals(5) = t_pdormtr_end - t_pdstedc_end
    t_intervals(6) = t_all_end - t_pdormtr_end
    t_intervals(7) = t_all_end - t_init

    if (my_rank == 0) then
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
