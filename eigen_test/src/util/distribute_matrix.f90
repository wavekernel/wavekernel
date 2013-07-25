module distribute_matrix
  use processes, only : layout_procs
  use matrix_io, only : sparse_mat
  implicit none

  type processor
    integer :: my_rank, n_procs, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
  end type processor

  type distribution
    integer :: dim, block_size, n_local_row, n_local_col
    integer :: desc(9) ! duplicated information for convenience
  end type distribution

  private
  public :: processor, distribution, &
       setup_distribution, setup_distributed_matrix, &
       copy_global_dense_matrix_to_local, copy_global_sparse_matrix_to_local, &
       create_dense_matrix, gather_matrix, allgather_row_wise

contains
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


  subroutine setup_distributed_matrix(conf, desc, mat)
    type(conf_distribution), intent(in) :: conf
    integer, intent(out) :: desc(9)
    real(kind(1.d0)), intent(out), allocatable :: mat(:, :)

    integer :: info

    call descinit(desc, conf%dim, conf%dim, conf%block_size, conf%block_size, &
         0, 0, conf%context, conf%n_local_row, info)
    allocate(mat(1 : conf%n_local_row, 1 : conf%n_local_col))
    mat(:, :) = 0.0d0
  end subroutine setup_distributed_matrix


  subroutine create_dense_matrix(verbose_level, mat_in, mat)
    use matrix_io, only : sparse_mat
    implicit none
    integer, intent(in) :: verbose_level
    type(sparse_mat), intent(in) :: mat_in
    real(kind(1.d0)), intent(out), allocatable :: mat(:,:)

    integer, allocatable :: mat_chk(:,:)
    integer :: k, i, j, n
    integer :: ierr

    logical :: debug_mode

    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif

    n = mat_in%size

    if (debug_mode) write(*,*)'@@ create_dense_matrix'

    allocate (mat(n, n), stat = ierr)
    if (ierr /= 0) stop 'Stop:error in alloc. for mat'
    mat(:, :) = 0.0d0

    allocate (mat_chk(n, n), stat = ierr)
    if (ierr /= 0) stop 'Stop:error in alloc. for mat'
    mat_chk(:, :) = 0

    do k = 1, mat_in%num_non_zeros
      i = mat_in%suffix(1, k)
      j = mat_in%suffix(2, k)
      mat(i, j) = mat_in%value(k)
      mat_chk(i, j) = mat_chk(i, j) + 1
      if (i /= j) then
        mat(j, i) = mat_in%value(k)
        mat_chk(j, i) = mat_chk(j, i) + 1
      endif
    enddo

    ! Check the multiple record
    do j = 1, n
      do i = 1, n
        if (mat_chk(i, j) > 1) then
          write(*,*)'ERROR(create_dense_matrix): multiple record : i,j =', i, j
          stop
        endif
      enddo
    enddo

    if (verbose_level >= 1)  then
      write(*,*)' No multiple record in mat_value ... OK!'
    endif
  end subroutine create_dense_matrix


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


  subroutine copy_global_dense_matrix_to_local(global_mat, desc, local_mat)
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
  end subroutine copy_global_dense_matrix_to_local


  subroutine copy_global_sparse_matrix_to_local(mat_in, desc, mat)
    type(sparse_mat), intent(in) :: mat_in
    integer, intent(in) :: desc(9)
    real(kind(1.d0)), intent(out) :: mat(:,:)

    integer :: k, i, j
    integer :: ierr

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
end module distribute_matrix



