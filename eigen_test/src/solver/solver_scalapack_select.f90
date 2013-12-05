!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_scalapack_select
  use descriptor_parameters
  use distribute_matrix, only : &
       get_local_cols, gather_matrix, allgather_row_wise, setup_distributed_matrix
  use eigenpairs_types, only : eigenpairs_types_union
  use processes, only : process
  use time, only : get_wclock_time
  implicit none

  private
  public :: eigen_solver_scalapack_select

contains
  subroutine eigen_solver_scalapack_select(proc, desc_A, A, n_vec, eigenpairs)
    include 'mpif.h'

    type(process) :: proc
    integer, intent(in) :: desc_A(9), n_vec
    double precision, intent(in) :: A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: ierr, info
    integer :: dim, work_size, iwork_size

    double precision, allocatable :: work(:)
    ! double precision, allocatable :: work_print(:)
    integer, allocatable :: iwork(:)

    ! For pdsyevx
    character :: jobz, range
    integer :: n_eigenvalues, n_eigenvectors
    integer, allocatable :: ifail(:), iclustr(:)
    double precision, allocatable :: gap(:)
    double precision :: abstol, orfac

    ! Time
    integer, parameter :: n_intervals = 1
    integer :: i
    double precision :: t_intervals(n_intervals)
    double precision :: t_init, t_all_end
    character(*), parameter :: interval_names(n_intervals) = (/'total'/)

    ! Functions
    double precision :: pdlamch

    call get_wclock_time(t_init)

    eigenpairs%type_number = 2

    dim = desc_A(rows_)

    call setup_distributed_matrix(proc, dim, dim, eigenpairs%blacs%desc, eigenpairs%blacs%Vectors)

    work_size = max(3, work_size_for_pdsyevx('V', dim, desc_A, dim))
    iwork_size = 6 * max(dim, proc%n_procs_row * proc%n_procs_col + 1, 4)
    allocate(eigenpairs%blacs%values(dim))
    allocate(work(work_size))
    allocate(iwork(iwork_size))
    allocate(ifail(dim))
    allocate(iclustr(2 * proc%n_procs_row * proc%n_procs_col))
    allocate(gap(proc%n_procs_row * proc%n_procs_col))

    jobz = 'V'
    range = 'I'
    abstol = 2.0 * pdlamch(desc_A(2), 'S')
    orfac = 0.0d0
    call pdsyevx(jobz, range, 'L', dim, A, 1, 1, Desc_A, &
         0, 0, 1, n_vec, abstol, n_eigenvalues, n_eigenvectors, eigenpairs%blacs%values, &
         orfac, eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         work, work_size, iwork, iwork_size, &
         ifail, iclustr, gap, info)
    if (proc%my_rank == 0) then
      call pdsyevx_report(proc%context, jobz, abstol, orfac, info, &
           n_eigenvalues, n_eigenvectors, ifail, iclustr)
      if (n_eigenvalues < n_vec) then
        print '("[Warning] eigen_solver_scalapack_select: PDSYEVX computed only ", I0, " of ", I0, " requested eigenvalues")', &
             n_eigenvalues, n_vec
      end if
    end if

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_all_end - t_init

    if (proc%my_rank == 0) then
       call MPI_Reduce(MPI_IN_PLACE, t_intervals, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
       print '("elapsed time (sec)")'
       do i = 1, n_intervals
         print '("  ", a, " ", f12.2)', interval_names(i), t_intervals(i)
       end do
    else
       call MPI_Reduce(t_intervals, 0, n_intervals, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    end if
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
      stop '[Error] work_size_for_pdsyevx: Unknown jobz value'
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
    double precision :: abstol, orfac

    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: i

    print '("pdsyevx report")'
    print '("  found eigenvalues: ", i0)', n_eigenvalues
    print '("  computed eigenvectors: ", i0)', n_eigenvectors
    print '("  abstol: ", e23.16)', abstol
    print '("  orfac: ", e23.16)', orfac
    print '("  info: ", i0)', info
    if (info /= 0 .and. jobz == 'V') then
      if (mod(info, 2) /= 0) then
        print *, ' ifail: ', ifail(n_eigenvalues + 1 :)
      end if
      write (*, '("  iclustr: ")', advance = 'no')
      call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)
      do i = 1, n_procs_row * n_procs_col
        write (*, '(i0, " - ", i0)', advance = 'no') iclustr(2 * i - 1), iclustr(2 * i)
        if (iclustr(2 * i) /= 0 .and. iclustr(2 * i + 1) == 0) then
          print *
          exit
        else
          write (*, '(", ")', advance = 'no')
        end if
      end do
    end if
  end subroutine pdsyevx_report
end module solver_scalapack_select
