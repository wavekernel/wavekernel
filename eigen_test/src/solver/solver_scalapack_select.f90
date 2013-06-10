!================================================================
! ELSES version 0.03
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module solver_scalapack_select
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
  subroutine eigen_solver_scalapack_select(mat, eigen_level)
    !
    use time, only : get_wclock_time !(routine)
    use processes, only : layout_procs
    implicit none
    include 'mpif.h'

    real(kind(1.d0)), intent(inout) :: mat(:,:)   ! ( n x n ) matrix
    real(kind(1.d0)), intent(out)   :: eigen_level(:)

    integer :: my_rank, n_procs, block_size, n_procs_row, n_procs_col, context
    integer :: my_proc_row, my_proc_col

    integer :: ierr, info
    integer :: dim, n_local_row, n_local_col, work_size, iwork_size
    integer :: desc_A(9), desc_Eigenvectors(9)

    real(kind(1.d0)), allocatable :: A(:,:)
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
    call descinit(desc_Eigenvectors, dim, dim, block_size, block_size, 0, 0, context, n_local_row, info)

    allocate(A(1 : n_local_row, 1 : n_local_col))
    allocate(Eigenvectors(1 : n_local_row, 1 : n_local_col))

    call copy_global_matrix_to_local(mat, desc_A, A)
    Eigenvectors(:, :) = 0.0

    work_size = max(3, work_size_for_pdsyevx('V', dim, desc_A, dim))
    iwork_size = 6 * max(dim, n_procs_row * n_procs_col + 1, 4)
    allocate(eigenvalues(dim))
    allocate(work(work_size))
    allocate(iwork(iwork_size))
    allocate(ifail(dim))
    allocate(iclustr(2 * n_procs_row * n_procs_col))
    allocate(gap(n_procs_row * n_procs_col))

    jobz = 'V'
    range = 'A'
    abstol = 2.0 * pdlamch(desc_A(2), 'S')
    orfac = 1.e-3_8
    call pdsyevx(jobz, range, 'L', dim, A, 1, 1, Desc_A, &
         0, 0, 0, 0, abstol, n_eigenvalues, n_eigenvectors, eigenvalues, &
         orfac, Eigenvectors, 1, 1, desc_Eigenvectors, &
         work, work_size, iwork, iwork_size, &
         ifail, iclustr, gap, info)
    if (my_rank == 0) then
      call pdsyevx_report(context, jobz, abstol, orfac, info, &
           n_eigenvalues, n_eigenvectors, ifail, iclustr)
    end if

    eigen_level(:) = eigenvalues(:)

    call gather_matrix(Eigenvectors, desc_Eigenvectors, 0, 0, mat)

    call get_wclock_time(t_all_end)

    t_intervals(1) = t_all_end - t_init

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
  end subroutine eigen_solver_scalapack_select

  subroutine gather_matrix(mat, desc, dest_row, dest_col, global_mat)
    real(kind(1.d0)), intent(in) :: mat(:, :)
    integer, intent(in) :: desc(9), dest_row, dest_col
    real(kind(1.d0)), intent(out) :: global_mat(:, :)

    real(kind(1.d0)), allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: m, n, mb, nb, m_local, n_local, context, m_recv, n_recv
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global

    integer :: numroc, iceil

    m = desc(3)
    n = desc(4)
    if (m /= size(global_mat, 1) .or. n /= size(global_mat, 2)) then
       stop 'gather_matrix: illegal matrix size'
    end if
    mb = desc(5)
    nb = desc(6)
    m_local = size(mat, 1)
    n_local = size(mat, 2)
    context = desc(2)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
       buf_size_row = numroc(m, mb, 0, 0, n_procs_row)
       buf_size_col = numroc(n, nb, 0, 0, n_procs_col)
       allocate(recv_buf(buf_size_row, buf_size_col))
    end if

    do pr = 0, n_procs_row - 1
       do pc = 0, n_procs_col - 1
          if (my_proc_row == pr .and. my_proc_col == pc) then
             call dgesd2d(context, m_local, n_local, mat, m_local, dest_row, dest_col)
          end if
          if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
             m_recv = numroc(m, mb, pr, 0, n_procs_row)
             n_recv = numroc(n, nb, pc, 0, n_procs_col)
             call dgerv2d(context, m_recv, n_recv, recv_buf, buf_size_row, pr, pc)
             do br = 1, iceil(iceil(m, mb) - my_proc_row, n_procs_row)
                m1_global = 1 + mb * (n_procs_row * (br - 1) + pr)
                m2_global = min(m1_global + mb - 1, m)
                m1_local = 1 + mb * (br - 1)
                m2_local = m1_local + m2_global - m1_global
                do bc = 1, iceil(iceil(n, nb) - my_proc_col, n_procs_col)
                   n1_global = 1 + nb * (n_procs_col * (bc - 1) + pc)
                   n2_global = min(n1_global + nb - 1, n)
                   n1_local = 1 + nb * (bc - 1)
                   n2_local = n1_local + n2_global - n1_global
                   global_mat(m1_global : m2_global, n1_global : n2_global) = &
                        recv_buf(m1_local : m2_local, n1_local : n2_local)
                end do
             end do
          end if
       end do
    end do

    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
       deallocate(recv_buf)
    end if
  end subroutine gather_matrix

  subroutine copy_global_matrix_to_local(global_mat, desc, local_mat)
    real(kind(1.d0)), intent(in) :: global_mat(:, :) ! assumed to be same in all the procs
    integer, intent(in) :: desc(9)
    real(kind(1.d0)), intent(out) :: local_mat(:, :)

    integer :: m, n, mb, nb, m_local, n_local, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: m_start, n_start, m_end, n_end, m_start_global, n_start_global
    integer :: i, j

    integer :: iceil

    m = desc(3)
    n = desc(4)
    mb = desc(5)
    nb = desc(6)
    m_local = size(local_mat, 1)
    n_local = size(local_mat, 2)
    context = desc(2)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    do j = 1, iceil(iceil(n, nb), n_procs_col)
       n_start = 1 + (j - 1) * nb
       if (n_start > n_local) then
          cycle
       end if
       n_end = min(j * nb, n_local)
       n_start_global = 1 + nb * (n_procs_col * (j - 1) + my_proc_col)
       do i = 1, iceil(iceil(m, mb), n_procs_row)
          m_start = 1 + mb * (i - 1)
          if (m_start > m_local) then
             cycle
          end if
          m_end = min(i * mb, m_local)
          m_start_global = 1 + mb * (n_procs_row * (i - 1) + my_proc_row)
          local_mat(m_start : m_end, n_start : n_end) = &
               global_mat(m_start_global : m_start_global + m_end - m_start, &
               n_start_global : n_start_global + n_end - n_start)
       end do
    end do
  end subroutine copy_global_matrix_to_local

  subroutine allgather_row_wise(array_local, context, block_size, array_global)
    real(kind(1.d0)), intent(in) :: array_local(:)
    integer, intent(in) :: context, block_size
    real(kind(1.d0)), intent(out) :: array_global(:)

    real(kind(1.d0)), allocatable :: send_buf(:)
    integer :: n_global, n_local, max_buf_size
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col, sender_proc_col
    integer :: b, s, s2, e

    integer :: numroc, iceil

    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    n_global = size(array_global)
    max_buf_size = numroc(n_global, block_size, 0, 0, n_procs_col) !iceil(n_global, n_procs_col)
    allocate(send_buf(max_buf_size))

    do sender_proc_col = 0, n_procs_col - 1
       n_local = numroc(n_global, block_size, sender_proc_col, 0, n_procs_col)
       if (my_proc_col == sender_proc_col) then
          if (n_local /= size(array_local)) stop 'wrong local array size'
          send_buf(1 : n_local) = array_local(1 : n_local)
          call dgebs2d(context, 'Row', 'I', max_buf_size, 1, send_buf, max_buf_size)
       else
          call dgebr2d(context, 'Row', 'I', max_buf_size, 1, send_buf, max_buf_size, &
               my_proc_row, sender_proc_col)
       end if
       do b = 1, iceil(n_local, n_procs_col)
          s = 1 + block_size * (n_procs_col * (b - 1) + sender_proc_col)
          e = min(s + block_size - 1, n_global)
          s2 = 1 + block_size * (b - 1)
          array_global(s : e) = send_buf(s2 : s2 + e - s)
       end do
    end do

    deallocate(send_buf)
  end subroutine allgather_row_wise

  subroutine print_mat(name, mat, m, n)
    character(*), intent(in) :: name
    real(kind(1.d0)), intent(in) :: mat(:, :)
    integer :: i, j, m, n
    if (m < 0) then
       m = size(mat, 1)
    end if
    if (n < 0) then
       n = size(mat, 2)
    end if
    do j = 1, n
       do i = 1, m
          print ' (A, "(", I6, ",", I6, ")=", D30.18) ', name, i, j, mat(i, j)
       end do
    end do
  end subroutine print_mat

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


