module distribute_matrix
  use descriptor_parameters
  use processes, only : layout_procs, terminate
  use matrix_io, only : sparse_mat
  implicit none

  type process
    integer :: my_rank, n_procs, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
  end type process

  public :: process, setup_distribution, &
       get_local_cols, setup_distributed_matrix, &
       copy_global_dense_matrix_to_local, copy_global_sparse_matrix_to_local, &
       create_dense_matrix, gather_matrix, allgather_row_wise, &
       desc_size, desc_type_, context_, rows_, cols_, block_row_, block_col_, &
       rsrc_, csrc_, local_rows_

contains

  subroutine setup_distribution(proc)
    type(process), intent(out) :: proc

    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call layout_procs(proc%n_procs, proc%n_procs_row, proc%n_procs_col)
    call blacs_get(-1, 0, proc%context)
    call blacs_gridinit(proc%context, 'R', proc%n_procs_row, proc%n_procs_col)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    if (proc%my_rank == 0) then
       print '( "procs: ", I5, " x ", I5, " (", I5, ")" )', &
            proc%n_procs_row, proc%n_procs_col, proc%n_procs
    end if

    if (proc%my_proc_row >= proc%n_procs_row .or. proc%my_proc_col >= proc%n_procs_col) then
       call blacs_exit(0)
       stop '[Warning] setup_distribution: Out of process grid, process exit'
    end if
  end subroutine setup_distribution


  integer function get_local_cols(proc, desc)
    type(process), intent(in) :: proc
    integer, intent(in) :: desc(desc_size)

    integer :: numroc

    get_local_cols = max(1, numroc(desc(cols_), desc(block_col_), &
         proc%my_proc_col, 0, proc%n_procs_col))
  end function get_local_cols


  subroutine setup_distributed_matrix(proc, rows, cols, desc, mat)
    type(process), intent(in) :: proc
    integer, intent(in) :: rows, cols
    integer, intent(out) :: desc(desc_size)
    double precision, intent(out), allocatable :: mat(:, :)

    integer :: numroc
    integer :: block_size, local_rows, info

    block_size = min(32, rows / max(proc%n_procs_row, proc%n_procs_col))
    if (proc%my_rank == 0) then
      print '( "Creating distributed matrix  M, N, MB, NB: ", &
           &I0, ", ", I0, ", ", I0, ", ", I0 )', rows, cols, block_size, block_size
    end if

    local_rows = max(1, numroc(rows, block_size, &
         proc%my_proc_row, 0, proc%n_procs_row))

    call descinit(desc, rows, cols, block_size, block_size, &
         0, 0, proc%context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('[Error] setup_distributed_matrix: descinit failed')
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols(proc, desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('[Error] setup_distributed_matrix: allocation failed')
    end if

    mat(:, :) = 0.0d0
  end subroutine setup_distributed_matrix


  subroutine create_dense_matrix(mat_in, mat)
    use matrix_io, only : sparse_mat

    type(sparse_mat), intent(in) :: mat_in
    double precision, intent(out), allocatable :: mat(:,:)

    integer :: k, i, j, n

    n = mat_in%size

    allocate (mat(n, n))
    mat(:, :) = 0.0d0

    do k = 1, mat_in%num_non_zeros
      i = mat_in%suffix(1, k)
      j = mat_in%suffix(2, k)
      mat(i, j) = mat_in%value(k)
      if (i /= j) then
        mat(j, i) = mat_in%value(k)
      endif
    enddo
  end subroutine create_dense_matrix


  subroutine gather_matrix(mat, desc, dest_row, dest_col, global_mat)
    use time, only : get_wclock_time

    double precision, intent(in) :: mat(:, :)
    integer, intent(in) :: desc(9), dest_row, dest_col
    double precision, intent(out) :: global_mat(:, :)

    double precision, allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: m, n, mb, nb, m_local, n_local, context, m_recv, n_recv
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global

    integer :: numroc, iceil

    double precision :: t_start, t_end

    call get_wclock_time(t_start)

    m = desc(3)
    n = desc(4)
    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      if (m /= size(global_mat, 1) .or. n /= size(global_mat, 2)) then
        stop '[Error] gather_matrix: Illegal matrix size'
      end if
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
     else
       allocate(recv_buf(0, 0)) ! Only destination process uses receive buffer
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

    call get_wclock_time(t_end)
    ! print *, 'gather time', my_proc_row, my_proc_col, (t_end - t_start)
  end subroutine gather_matrix


  subroutine copy_global_dense_matrix_to_local(global_mat, desc, local_mat)
    double precision, intent(in) :: global_mat(:, :) ! assumed to be same in all the procs
    integer, intent(in) :: desc(9)
    double precision, intent(out) :: local_mat(:, :)

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
  end subroutine copy_global_dense_matrix_to_local


  subroutine copy_global_sparse_matrix_to_local(mat_in, desc, mat)
    type(sparse_mat), intent(in) :: mat_in
    integer, intent(in) :: desc(9)
    double precision, intent(out) :: mat(:,:)

    integer :: k, i, j

    do k = 1, mat_in%num_non_zeros
      i = mat_in%suffix(1, k)
      j = mat_in%suffix(2, k)
      call pdelset(mat, i, j, desc, mat_in%value(k))
      if (i /= j) then
      call pdelset(mat, j, i, desc, mat_in%value(k))
      endif
    enddo
  end subroutine copy_global_sparse_matrix_to_local


  subroutine allgather_row_wise(array_local, context, block_size, array_global)
    double precision, intent(in) :: array_local(:)
    integer, intent(in) :: context, block_size
    double precision, intent(out) :: array_global(:)

    double precision, allocatable :: send_buf(:)
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
          if (n_local /= size(array_local)) then
            stop '[Error] distribute_matrix: Wrong local array size'
          end if
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
end module distribute_matrix
