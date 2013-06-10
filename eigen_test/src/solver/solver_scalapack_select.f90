!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_scalapack_select

  use time, only : get_wclock_time
  use distribute_matrix, only : conf_distribution, gather_matrix, allgather_row_wise
  !
  implicit none
  !
  private
  public :: eigen_solver_scalapack_select
  !
contains
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eigen_solver_scalapack_select(conf, desc_A, A, eigen_level, eigenvectors_global)
    !
    implicit none
    include 'mpif.h'

    type(conf_distribution) :: conf
    integer, intent(in) :: desc_A(9)
    real(kind(1.d0)), intent(in) :: A(:, :)
    real(kind(1.d0)), intent(out)   :: eigen_level(:), eigenvectors_global(:, :)

    integer :: ierr, info
    integer :: work_size, iwork_size
    integer :: desc_Eigenvectors(9)

    real(kind(1.d0)), allocatable :: Eigenvectors(:, :)
    real(kind(1.d0)), allocatable :: work(:), work_print(:)
    integer, allocatable :: iwork(:)

    ! For pdsyevx
    character :: jobz, range
    integer :: n_eigenvalues, n_eigenvectors
    integer, allocatable :: ifail(:), iclustr(:)
    real(kind(1.d0)), allocatable :: eigenvalues(:), gap(:)
    real(kind(1.d0)) :: abstol, orfac

    ! Time
    integer, parameter :: n_intervals = 1
    integer :: i
    real(kind(1.d0)) :: t_intervals(n_intervals)
    real(kind(1.d0)) :: t_init, t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'total'/)

    ! Functions
    integer :: numroc
    real(kind(1.d0)) :: pdlamch

    call get_wclock_time(t_init)

    allocate(Eigenvectors(1 : conf%n_local_row, 1 : conf%n_local_col))
    Eigenvectors(:, :) = 0.0

    work_size = max(3, work_size_for_pdsyevx('V', conf%dim, desc_A, conf%dim))
    iwork_size = 6 * max(conf%dim, conf%n_procs_row * conf%n_procs_col + 1, 4)
    allocate(eigenvalues(conf%dim))
    allocate(work(work_size))
    allocate(iwork(iwork_size))
    allocate(ifail(conf%dim))
    allocate(iclustr(2 * conf%n_procs_row * conf%n_procs_col))
    allocate(gap(conf%n_procs_row * conf%n_procs_col))

    jobz = 'V'
    range = 'A'
    abstol = 2.0 * pdlamch(desc_A(2), 'S')
    orfac = 1.e-3_8
    call pdsyevx(jobz, range, 'L', conf%dim, A, 1, 1, Desc_A, &
         0, 0, 0, 0, abstol, n_eigenvalues, n_eigenvectors, eigenvalues, &
         orfac, Eigenvectors, 1, 1, desc_Eigenvectors, &
         work, work_size, iwork, iwork_size, &
         ifail, iclustr, gap, info)
    if (conf%my_rank == 0) then
      call pdsyevx_report(conf%context, jobz, abstol, orfac, info, &
           n_eigenvalues, n_eigenvectors, ifail, iclustr)
    end if

    eigen_level(:) = eigenvalues(:)

    call gather_matrix(Eigenvectors, desc_Eigenvectors, 0, 0, eigenvectors_global)

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_all_end - t_init

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
  end subroutine eigen_solver_scalapack_select

  integer function work_size_for_pdsyevx(jobz, n, desc, neig) result (size)
    character(len = 1), intent(in) :: jobz
    integer, intent(in) :: n, desc(9), neig

    integer :: context, n_procs_row, n_procs_col, my_proc_row, my_proc_col, block_size
    integer :: nn, np0, mq0, anb, sqnpc, nps, nsytrd_lwopt
    integer :: numroc, iceil, pjlaenv

    context = desc(2)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    block_size = desc(5)
    nn = MAX(n, block_size, 2)
    np0 = numroc(nn, block_size, 0, 0, n_procs_row)
    if (jobz == 'N') then
      size = 5 * n + max( 5 * nn, block_size * (np0 + 1))
    else if (jobz == 'V') then
      mq0 = numroc(max(neig, block_size, 2), block_size, 0, 0, n_procs_col)
      size = 5 * n + MAX( 5 * nn, np0 * mq0 + 2 * block_size * block_size ) + &
           iceil(neig, n_procs_row * n_procs_col) * nn
    else
      stop 'unknown jobz value'
    end if

    anb = pjlaenv(desc(2), 3, 'pdsyttrd', 'l', 0, 0, 0, 0)
    sqnpc = int(sqrt(dble(n_procs_row * n_procs_col)))
    nps = max(numroc(n, 1, 0, 0, sqnpc), 2 * anb)
    nsytrd_lwopt = n + 2 * (anb + 1) * (4 * nps + 2) + (nps + 3) * nps
    size = max(size, 5 * n + nsytrd_lwopt)
  end function work_size_for_pdsyevx

  subroutine pdsyevx_report(context, jobz, abstol, orfac, info, &
       n_eigenvalues, n_eigenvectors, ifail, iclustr)
    integer, intent(in) :: context, info, n_eigenvalues, n_eigenvectors, ifail(:), iclustr(:)
    character, intent(in) :: jobz
    real(kind(1.d0)) :: abstol, orfac

    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: i

    print *, 'PDSYEVX report'
    print *, ' found eigenvalues: ', n_eigenvalues
    print *, ' computed eigenvectors: ', n_eigenvectors
    print *, ' abstol: ', abstol
    print *, ' orfac: ', orfac
    print *, ' info: ', info
    if (info /= 0 .and. jobz == 'V') then
      if (mod(info, 2) /= 0) then
        print *, ' ifail: ', ifail(n_eigenvalues + 1 :)
      end if
      write (*, '(A)', advance = 'no') '  iclustr: '
      call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)
      do i = 1, n_procs_row * n_procs_col
        write (*, '(I6, A, I6)', advance = 'no') iclustr(2 * i - 1), ' - ', iclustr(2 * i)
        if (iclustr(2 * i) /= 0 .and. iclustr(2 * i + 1) == 0) then
          print *
          exit
        else
          write (*, '(A)', advance = 'no') ', '
        end if
      end do
    end if
  end subroutine pdsyevx_report

end module solver_scalapack_select
