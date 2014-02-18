module distribute_matrix
  use descriptor_parameters
  use matrix_io, only : sparse_mat
  use processes, only : process, layout_procs, terminate
  implicit none

  public :: get_local_cols, setup_distributed_matrix, &
       distribute_global_dense_matrix, distribute_global_sparse_matrix, &
       convert_sparse_matrix_to_dense, gather_matrix, allgather_row_wise

contains

  integer function get_local_cols(proc, desc)
    type(process), intent(in) :: proc
    integer, intent(in) :: desc(desc_size)

    integer :: numroc

    get_local_cols = max(1, numroc(desc(cols_), desc(block_col_), &
         proc%my_proc_col, 0, proc%n_procs_col))
  end function get_local_cols


  subroutine setup_distributed_matrix(name, proc, rows, cols, desc, mat, &
       block_size)
    character(*), intent(in) :: name
    type(process), intent(in) :: proc
    integer, intent(in) :: rows, cols
    integer, intent(out) :: desc(desc_size)
    double precision, intent(out), allocatable :: mat(:, :)
    integer, intent(in), optional :: block_size

    integer :: numroc
    integer :: actual_block_size, local_rows, info

    if (present(block_size)) then
      actual_block_size = block_size
    else
      actual_block_size = 32
    end if

    actual_block_size = min(actual_block_size, rows / max(proc%n_procs_row, proc%n_procs_col))
    if (proc%my_rank == 0) then
      print '( "Creating distributed matrix ", A, " with M, N, MB, NB: ", &
           &I0, ", ", I0, ", ", I0, ", ", I0 )', name, rows, cols, &
           actual_block_size, actual_block_size
    end if

    local_rows = max(1, numroc(rows, actual_block_size, &
         proc%my_proc_row, 0, proc%n_procs_row))

    call descinit(desc, rows, cols, actual_block_size, actual_block_size, &
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


  subroutine convert_sparse_matrix_to_dense(mat_in, mat)
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
  end subroutine convert_sparse_matrix_to_dense


  subroutine gather_matrix(mat, desc, dest_row, dest_col, global_mat)
    double precision, intent(in) :: mat(:, :)
    integer, intent(in) :: desc(desc_size), dest_row, dest_col
    double precision, intent(out) :: global_mat(:, :)

    double precision, allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: m, n, mb, nb, m_local, n_local, context, m_recv, n_recv
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global

    integer :: numroc, iceil

    m = desc(rows_)
    n = desc(cols_)
    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      if (m /= size(global_mat, 1) .or. n /= size(global_mat, 2)) then
        stop '[Error] gather_matrix: Illegal matrix size'
      end if
    end if
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(mat, 1)
    n_local = size(mat, 2)
    context = desc(context_)
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
  end subroutine gather_matrix


  subroutine distribute_global_dense_matrix(global_mat, desc, local_mat)
    double precision, intent(in) :: global_mat(:, :) ! assumed to be same in all the procs
    integer, intent(in) :: desc(desc_size)
    double precision, intent(out) :: local_mat(:, :)

    integer :: m, n, mb, nb, m_local, n_local, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: m_start, n_start, m_end, n_end, m_start_global, n_start_global
    integer :: i, j

    integer :: iceil

    m = desc(rows_)
    n = desc(cols_)
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(local_mat, 1)
    n_local = size(local_mat, 2)
    context = desc(context_)
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
  end subroutine distribute_global_dense_matrix


  subroutine distribute_global_sparse_matrix(mat_in, desc, mat)
    type(sparse_mat), intent(in) :: mat_in
    integer, intent(in) :: desc(desc_size)
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
  end subroutine distribute_global_sparse_matrix


  ! Supposing arrays of the same size are distributed by the block cyclic manner
  ! in each processor row, broadcast whole pieces of the array then share the
  ! array globally in the processor row.
  !
  ! array_local: local input, dimension LOCc(N) where N is the dimension of the
  ! global array.
  subroutine allgather_row_wise(array_local, context, block_size, array_global)
    double precision, intent(in) :: array_local(:)
    integer, intent(in) :: context, block_size
    double precision, intent(out) :: array_global(:)

    double precision, allocatable :: send_buf(:)
    integer :: n_global, n_local, max_buf_size, ierr
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col, sender_proc_col
    integer :: block_, head_global, tail_global, head_local

    integer :: numroc, iceil

    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    n_global = size(array_global)
    max_buf_size = numroc(n_global, block_size, 0, 0, n_procs_col)
    allocate(send_buf(max_buf_size))

    do sender_proc_col = 0, n_procs_col - 1
      n_local = numroc(n_global, block_size, sender_proc_col, 0, n_procs_col)
      if (my_proc_col == sender_proc_col) then
        if (n_local /= size(array_local)) then
          stop '[Error] distribute_matrix: Wrong local array size'
        end if
        send_buf(1 : n_local) = array_local(1 : n_local)
        call dgebs2d(context, 'Row', ' ', max_buf_size, 1, send_buf, max_buf_size)
      else
        call dgebr2d(context, 'Row', ' ', max_buf_size, 1, send_buf, max_buf_size, &
             my_proc_row, sender_proc_col)
      end if
      ! Slice entries in a process into blocks
      do block_ = 1, iceil(n_local, block_size)
        head_global = 1 + block_size * (n_procs_col * (block_ - 1) + sender_proc_col)
        tail_global = min(head_global + block_size - 1, n_global)
        head_local = 1 + block_size * (block_ - 1)
        array_global(head_global : tail_global) = &
             send_buf(head_local : head_local + tail_global - head_global)
      end do
    end do
  end subroutine allgather_row_wise
end module distribute_matrix
