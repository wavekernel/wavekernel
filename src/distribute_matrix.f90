module wk_distribute_matrix_m
  use mpi
  use wk_event_logger_m
  use wk_descriptor_parameters_m
  use wk_global_variables_m
  use wk_matrix_io_m
  use wk_processes_m, only : check_master, wk_process_t, layout_procs, terminate
  implicit none

  type wk_local_matrix_t
    real(8), allocatable :: val(:, :)
  end type wk_local_matrix_t

  type wk_distributed_block_matrices_t
    ! Global values.
    integer :: m, n, num_blocks, num_procs
    integer, allocatable :: block_to_row(:), block_to_col(:)
    ! Local values.
    integer :: my_rank
    type(wk_local_matrix_t), allocatable :: local_matrices(:)
  end type wk_distributed_block_matrices_t

  public :: get_local_cols_wk, setup_distributed_matrix_complex, setup_distributed_matrix_real, &
       distribute_matrix_real_part, distribute_matrix_complex, &
       gather_matrix_complex, gather_matrix_real_part, gather_matrix_complex_part, &
       read_distribute_eigenvalues, read_distribute_eigenvectors, &
       bcast_sparse_matrix

contains

  integer function get_local_cols_wk(proc, desc)
    type(wk_process_t), intent(in) :: proc
    integer, intent(in) :: desc(desc_size)

    integer :: numroc

    get_local_cols_wk = max(1, numroc(desc(cols_), desc(block_col_), &
         proc%my_proc_col, 0, proc%n_procs_col))
  end function get_local_cols_wk


  ! If there is a process which owns no entries in given block size
  ! configuration, diminish the block size and warn about in.
  subroutine correct_block_size(name, dim, num_procs, block_size)
    character(*), intent(in) :: name
    integer, intent(in) :: dim, num_procs
    integer, intent(inout) :: block_size

    integer :: max_block_size

    max_block_size = max(dim / num_procs, 1)
    if (block_size > max_block_size) then
      if (check_master() .and. trim(name) /= '') then
        write(0, *) '[Warning] correct_block_size: dimension of matrix ' // &
             trim(name) // ' is very small relative to the number of processes'
      end if
      block_size = max_block_size
    end if
  end subroutine correct_block_size


  subroutine setup_distributed_matrix_real(name, proc, rows, cols, &
       desc, mat, is_square_block, block_size_row, block_size_col)
    character(*), intent(in) :: name
    type(wk_process_t), intent(in) :: proc
    integer, intent(in) :: rows, cols

    integer, intent(out) :: desc(desc_size)
    real(8), intent(out), allocatable :: mat(:, :)
    logical, intent(in), optional :: is_square_block
    integer, intent(in), optional :: block_size_row, block_size_col

    integer :: numroc
    integer :: actual_block_size_row, actual_block_size_col, &
         local_rows, info

    if (present(block_size_row)) then
      actual_block_size_row = block_size_row
    else
      actual_block_size_row = g_wk_block_size
      call correct_block_size(name, rows, proc%n_procs_row, actual_block_size_row)
    end if
    if (present(block_size_col)) then
      actual_block_size_col = block_size_col
    else
      actual_block_size_col = g_wk_block_size
      call correct_block_size(name, cols, proc%n_procs_col, actual_block_size_col)
    end if

    if (present(is_square_block)) then
      if (is_square_block) then
        actual_block_size_row = min(actual_block_size_row, actual_block_size_col)
        actual_block_size_col = actual_block_size_row
      end if
    end if

    if (check_master() .and. trim(name) /= '') then
      write (0, '(A, F16.6, 3A, I0, ", ", I0, ", ", I0, ", ", I0)') &
           ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] setup_distributed_matrix_real() ', trim(name), ' with M, N, MB, NB: ', &
           rows, cols, actual_block_size_row, actual_block_size_col
    end if

    local_rows = max(1, numroc(rows, actual_block_size_row, &
         proc%my_proc_row, 0, proc%n_procs_row))

    call descinit(desc, rows, cols, actual_block_size_row, actual_block_size_col, &
         0, 0, proc%context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix_real: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols_wk(proc, desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('setup_distributed_matrix_real: allocation failed', info)
    end if

    mat(:, :) = 0d0
  end subroutine setup_distributed_matrix_real


  subroutine setup_distributed_matrix_complex(name, proc, rows, cols, &
       desc, mat, is_square_block, block_size_row, block_size_col)
    character(*), intent(in) :: name
    type(wk_process_t), intent(in) :: proc
    integer, intent(in) :: rows, cols

    integer, intent(out) :: desc(desc_size)
    complex(kind(0d0)), intent(out), allocatable :: mat(:, :)
    logical, intent(in), optional :: is_square_block
    integer, intent(in), optional :: block_size_row, block_size_col

    integer :: numroc
    integer :: actual_block_size_row, actual_block_size_col, &
         local_rows, info

    if (present(block_size_row)) then
      actual_block_size_row = block_size_row
    else
      actual_block_size_row = g_wk_block_size
      call correct_block_size(name, rows, proc%n_procs_row, actual_block_size_row)
    end if
    if (present(block_size_col)) then
      actual_block_size_col = block_size_col
    else
      actual_block_size_col = g_wk_block_size
      call correct_block_size(name, cols, proc%n_procs_col, actual_block_size_col)
    end if

    if (present(is_square_block)) then
      if (is_square_block) then
        actual_block_size_row = min(actual_block_size_row, actual_block_size_col)
        actual_block_size_col = actual_block_size_row
      end if
    end if

    if (check_master() .and. trim(name) /= '') then
      write (0, '(A, F16.6, 3A, I0, ", ", I0, ", ", I0, ", ", I0)') &
           ' [Event', mpi_wtime() - g_wk_mpi_wtime_init, &
           '] setup_distributed_matrix_complex() ', trim(name), ' with M, N, MB, NB: ', &
           rows, cols, actual_block_size_row, actual_block_size_col
    end if

    local_rows = max(1, numroc(rows, actual_block_size_row, &
         proc%my_proc_row, 0, proc%n_procs_row))

    call descinit(desc, rows, cols, actual_block_size_row, actual_block_size_col, &
         0, 0, proc%context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix_complex: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols_wk(proc, desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('setup_distributed_matrix_complex: allocation failed', info)
    end if

    mat(:, :) = (0d0, 0d0)
  end subroutine setup_distributed_matrix_complex


  subroutine gather_matrix_complex(context, mat_src_desc, mat_src, dest_row, dest_col, mat_dest_lld, mat_dest)
    integer, intent(in) :: context, mat_src_desc(desc_size), dest_row, dest_col, mat_dest_lld
    complex(kind(0d0)), intent(in) :: mat_src(:, :)
    complex(kind(0d0)), intent(out) :: mat_dest(*)

    integer :: m, n, mb, nb, mb_, nb_, bi, bj, i, j, il, jl, k
    integer :: nprow, npcol, myprow, mypcol
    integer :: indxg2l

    call blacs_gridinfo(context, nprow, npcol, myprow, mypcol)
    m = mat_src_desc(rows_)
    n = mat_src_desc(cols_)
    mb = mat_src_desc(block_row_)
    nb = mat_src_desc(block_col_)
    mat_dest(1 : m * n) = (0d0, 0d0)

    do bi = 0, (m - mb * myprow - 1) / (mb * nprow)
      do bj = 0, (n - nb * mypcol - 1) / (nb * npcol)
        i = mb * nprow * bi + mb * myprow + 1  ! Global start index of block to copy.
        j = nb * npcol * bj + nb * mypcol + 1
        mb_ = min(mb, m - i + 1)  ! The last block may be smaller than (mb, nb).
        nb_ = min(nb, n - j + 1)
        il = indxg2l(i, mb, 0, 0, nprow)
        jl = indxg2l(j, nb, 0, 0, npcol)
        k = mat_dest_lld * (j - 1) + i
        call zlacpy('A', mb_, nb_, mat_src(il, jl), &
             mat_src_desc(local_rows_), mat_dest(k), mat_dest_lld)
      end do
    end do
    call zgsum2d(context, 'A', ' ', m, n, mat_dest, mat_dest_lld, dest_row, dest_col)
  end subroutine gather_matrix_complex


  subroutine gather_matrix_real_part(mat, desc, i, j, m, n, dest_row, dest_col, global_mat)
    real(8), intent(in) :: mat(:, :)
    integer, intent(in) :: i, j, m, n, desc(desc_size), dest_row, dest_col
    real(8), intent(out) :: global_mat(:, :)

    real(8), allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: rows, cols, mb, nb, m_local, n_local, context, m_send, n_send
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global
    integer :: ierr
    real(8) :: time_start, time_end

    integer :: numroc, iceil, indxg2l

    time_start = mpi_wtime()

    rows = desc(rows_)
    cols = desc(cols_)
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(mat, 1)
    n_local = size(mat, 2)
    context = desc(context_)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      buf_size_row = numroc(rows, mb, 0, 0, n_procs_row)
      buf_size_col = numroc(cols, nb, 0, 0, n_procs_col)
      allocate(recv_buf(buf_size_row, buf_size_col), stat = ierr)
    else
      allocate(recv_buf(0, 0), stat = ierr) ! Only destination process uses receive buffer
    end if
    if (ierr /= 0) then
      call terminate('gather_matrix_real_part: allocation failed', ierr)
    end if

    do pr = 0, n_procs_row - 1
      do pc = 0, n_procs_col - 1
        m_send = numroc(rows, mb, pr, 0, n_procs_row)
        n_send = numroc(cols, nb, pc, 0, n_procs_col)
        if (my_proc_row == pr .and. my_proc_col == pc) then
          call dgesd2d(context, m_send, n_send, mat, desc(local_rows_), dest_row, dest_col)
        end if
        if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
          call dgerv2d(context, m_send, n_send, recv_buf, buf_size_row, pr, pc)
          do br = 1, iceil(iceil(rows, mb) - pr, n_procs_row)
            m1_global = 1 + mb * (n_procs_row * (br - 1) + pr)
            m2_global = m1_global + mb - 1
            m1_global = max(m1_global, i)
            m2_global = min(m2_global, i + m - 1)
            if (m1_global > m2_global) then
              cycle
            end if
            m1_local = indxg2l(m1_global, mb, 0, 0, n_procs_row)
            m2_local = indxg2l(m2_global, mb, 0, 0, n_procs_row)
            do bc = 1, iceil(iceil(cols, nb) - pc, n_procs_col)
              n1_global = 1 + nb * (n_procs_col * (bc - 1) + pc)
              n2_global = n1_global + nb - 1
              !print *, 'D1', n1_global, n2_global, j, j + n - 1
              n1_global = max(n1_global, j)
              n2_global = min(n2_global, j + n - 1)
              !print *, 'D2', n1_global, n2_global
              if (n1_global > n2_global) then
                cycle
              end if
              n1_local = indxg2l(n1_global, nb, 0, 0, n_procs_col)
              n2_local = indxg2l(n2_global, nb, 0, 0, n_procs_col)
              !print *, 'VC', m1_global - i + 1, m2_global - i + 1, ',', n1_global - j + 1, n2_global - j + 1, ',', &
              !     m1_local, m2_local, ',', n1_local, n2_local
              global_mat(m1_global - i + 1 : m2_global - i + 1, n1_global - j + 1 : n2_global - j + 1) = &
                   recv_buf(m1_local : m2_local, n1_local : n2_local)
            end do
          end do
        end if
      end do
    end do

    time_end = mpi_wtime()
    call add_event('gather_matrix_real_part', time_end - time_start)
  end subroutine gather_matrix_real_part


  subroutine gather_matrix_complex_part(mat, desc, i, j, m, n, dest_row, dest_col, global_mat)
    complex(kind(0d0)), intent(in) :: mat(:, :)
    integer, intent(in) :: i, j, m, n, desc(desc_size), dest_row, dest_col
    complex(kind(0d0)), intent(out) :: global_mat(:, :)

    complex(kind(0d0)), allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: rows, cols, mb, nb, m_local, n_local, context, m_recv, n_recv
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global
    integer :: ierr
    real(8) :: time_start, time_end

    integer :: numroc, iceil, indxg2l

    time_start = mpi_wtime()

    !if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
    !  if (m /= size(global_mat, 1) .or. n /= size(global_mat, 2)) then
    !    call terminate('gather_matrix: illegal matrix size', 1)
    !  end if
    !end if

    rows = desc(rows_)
    cols = desc(cols_)
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(mat, 1)
    n_local = size(mat, 2)
    context = desc(context_)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      buf_size_row = numroc(rows, mb, 0, 0, n_procs_row)
      buf_size_col = numroc(cols, nb, 0, 0, n_procs_col)
      allocate(recv_buf(buf_size_row, buf_size_col), stat = ierr)
    else
      allocate(recv_buf(0, 0), stat = ierr) ! Only destination process uses receive buffer
    end if
    if (ierr /= 0) then
      call terminate('gather_matrix_complex_part: allocation failed', ierr)
    end if

    do pr = 0, n_procs_row - 1
       do pc = 0, n_procs_col - 1
          if (my_proc_row == pr .and. my_proc_col == pc) then
             call zgesd2d(context, m_local, n_local, mat, m_local, dest_row, dest_col)
          end if
          if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
             m_recv = numroc(rows, mb, pr, 0, n_procs_row)
             n_recv = numroc(cols, nb, pc, 0, n_procs_col)
             call zgerv2d(context, m_recv, n_recv, recv_buf, buf_size_row, pr, pc)
             do br = 1, iceil(iceil(rows, mb) - pr, n_procs_row)
                m1_global = max(1 + mb * (n_procs_row * (br - 1) + pr), i)
                m2_global = min(m1_global + mb - 1, i + m - 1)
                if (m1_global > m2_global) then
                  cycle
                end if
                m1_local = indxg2l(m1_global, mb, 0, 0, n_procs_row)
                m2_local = indxg2l(m2_global, mb, 0, 0, n_procs_row)
                do bc = 1, iceil(iceil(cols, nb) - pc, n_procs_col)
                   n1_global = max(1 + nb * (n_procs_col * (bc - 1) + pc), j)
                   n2_global = min(n1_global + nb - 1, j + n - 1)
                   if (n1_global > n2_global) then
                     cycle
                   end if
                   n1_local = indxg2l(n1_global, nb, 0, 0, n_procs_col)
                   n2_local = indxg2l(n2_global, nb, 0, 0, n_procs_col)
                   global_mat(m1_global - i + 1 : m2_global - i + 1, n1_global - j + 1 : n2_global - j + 1) = &
                        recv_buf(m1_local : m2_local, n1_local : n2_local)
                end do
             end do
          end if
       end do
    end do

    time_end = mpi_wtime()
    call add_event('gather_matrix_complex_part', time_end - time_start)
  end subroutine gather_matrix_complex_part


  subroutine gather_matrix_complex_with_pzgemr2d(context, mat_src_desc, mat_src, &
       dest_row, dest_col, mat_dest_lld, mat_dest)
    use mpi
    integer, intent(in) :: context, mat_src_desc(desc_size), dest_row, dest_col, mat_dest_lld
    complex(kind(0d0)), intent(in) :: mat_src(:, :)
    ! Needed to be allocated only on the destination process.
    complex(kind(0d0)), intent(out) :: mat_dest(*)

    integer :: usermap(1, 1), new_context, mat_dest_desc(desc_size), ierr, m, n
    integer :: blacs_pnum
    integer :: nprow, npcol, myrow, mycol

    call blacs_gridinfo(context, nprow, npcol, myrow, mycol )
    m = mat_src_desc(rows_)
    n = mat_src_desc(cols_)
    new_context = context
    usermap(1, 1) = blacs_pnum(context, dest_row, dest_col)
    call blacs_gridmap(new_context, usermap, 1, 1, 1)

    if (dest_row == myrow .and. dest_col == mycol) then
      call descinit(mat_dest_desc, m, n, m, n, 0, 0, new_context, mat_dest_lld, ierr)
    end if
    call mpi_bcast(mat_dest_desc, desc_size, mpi_integer, usermap(1, 1), mpi_comm_world, ierr)

    call pzgemr2d(m, n, mat_src, 1, 1, mat_src_desc, mat_dest, 1, 1, mat_dest_desc, context)
    if (dest_row == myrow .and. dest_col == mycol) then
      call blacs_gridexit(new_context)
    end if
  end subroutine gather_matrix_complex_with_pzgemr2d


  subroutine distribute_matrix_real_part(context, m, n, i, j, src_row, src_col, mat_src_lld, mat_src, &
       mat_dest_desc, mat_dest)
    integer, intent(in) :: context, m, n, i, j, src_row, src_col, mat_src_lld, mat_dest_desc(desc_size)
    real(8), intent(in) :: mat_src(*)
    real(8), intent(out) :: mat_dest(:, :)

    integer :: usermap(1, 1), new_context, mat_src_desc(desc_size), ierr
    integer :: rows, cols, nprow, npcol, myrow, mycol
    integer :: blacs_pnum

    rows = mat_dest_desc(rows_)
    cols = mat_dest_desc(cols_)

    call blacs_gridinfo(context, nprow, npcol, myrow, mycol)
    new_context = context
    usermap(1, 1) = blacs_pnum(context, src_row, src_col)
    call blacs_gridmap(new_context, usermap, 1, 1, 1)

    if (src_row == myrow .and. src_col == mycol) then
      call descinit(mat_src_desc, m, n, m, n, 0, 0, new_context, mat_src_lld, ierr)
    end if
    call mpi_bcast(mat_src_desc, desc_size, mpi_integer, usermap(1, 1), mpi_comm_world, ierr)
    if (src_row /= myrow .or. src_col /= mycol) then
      mat_src_desc(context_) = -1
    end if

    call pdgemr2d(m, n, mat_src, 1, 1, mat_src_desc, &
         mat_dest, i, j, mat_dest_desc, context)
    if (src_row == myrow .and. src_col == mycol) then
      call blacs_gridexit(new_context)
    end if
  end subroutine distribute_matrix_real_part


  subroutine distribute_matrix_complex(context, src_row, src_col, mat_src_lld, mat_src, &
       mat_dest_desc, mat_dest)
    use mpi
    integer, intent(in) :: context, src_row, src_col, mat_src_lld, mat_dest_desc(desc_size)
    complex(kind(0d0)), intent(in) :: mat_src(*)
    complex(kind(0d0)), intent(out) :: mat_dest(:, :)

    integer :: usermap(1, 1), new_context, mat_src_desc(desc_size), ierr
    integer :: m, n, nprow, npcol, myrow, mycol
    integer :: blacs_pnum

    m = mat_dest_desc(rows_)
    n = mat_dest_desc(cols_)

    call blacs_gridinfo(context, nprow, npcol, myrow, mycol)
    new_context = context
    usermap(1, 1) = blacs_pnum(context, src_row, src_col)
    call blacs_gridmap(new_context, usermap, 1, 1, 1)

    if (src_row == myrow .and. src_col == mycol) then
      call descinit(mat_src_desc, m, n, m, n, 0, 0, new_context, mat_src_lld, ierr)
    end if
    call mpi_bcast(mat_src_desc, desc_size, mpi_integer, usermap(1, 1), mpi_comm_world, ierr)

    call pzgemr2d(m, n, mat_src, 1, 1, mat_src_desc, &
         mat_dest, 1, 1, mat_dest_desc, context)
    if (src_row == myrow .and. src_col == mycol) then
      call blacs_gridexit(new_context)
    end if
  end subroutine distribute_matrix_complex


  subroutine distribute_global_sparse_matrix_wk(mat_in, desc, mat)
    type(sparse_mat), intent(in) :: mat_in
    integer, intent(in) :: desc(desc_size)
    real(8), intent(out) :: mat(:, :)

    integer :: k, i, j
    real(8) :: time_start, time_end

    time_start = mpi_wtime()

    do k = 1, mat_in%num_non_zeros
      i = mat_in%suffix(1, k)
      j = mat_in%suffix(2, k)
      call pdelset(mat, i, j, desc, mat_in%value(k))
      if (i /= j) then
      call pdelset(mat, j, i, desc, mat_in%value(k))
      endif
    enddo

    time_end = mpi_wtime()
    call add_event('distribute_global_sparse_matrix', time_end - time_start)
  end subroutine distribute_global_sparse_matrix_wk


  subroutine read_distribute_eigenvalues(filename, dim, full_vecs, col_eigenvalues, full_vecs_desc)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: dim, col_eigenvalues, full_vecs_desc(desc_size)
    complex(kind(0d0)), intent(out) :: full_vecs(:, :)

    integer, parameter :: iunit = 14
    integer :: line, i, ierr
    real(8) :: re
    complex(kind(0d0)) :: eigenvalues(dim)

    if (check_master()) then
      open(iunit, file=filename)
      do line = 1, dim
        read(iunit, *) i, re
        eigenvalues(i) = cmplx(re, 0d0, kind(0d0))
      end do
      close(iunit)
    end if
    call mpi_bcast(eigenvalues, dim, mpi_double_complex, g_wk_master_pnum, mpi_comm_world, ierr)
    do i = 1, dim
      call pzelset(full_vecs, i, col_eigenvalues, full_vecs_desc, eigenvalues(i))
    end do
  end subroutine read_distribute_eigenvalues


  subroutine read_distribute_eigenvectors(dirname, dim, fst_filter, num_filter, &
         Y_filtered, Y_filtered_desc)
    character(len=*), intent(in) :: dirname
    integer, intent(in) :: dim, fst_filter, num_filter, Y_filtered_desc(desc_size)
    complex(kind(0d0)), intent(out) :: Y_filtered(:, :)

    integer, parameter :: iunit = 15, max_num_digits = 8
    integer :: digit, j, i, l, k, len, ierr
    character(1024) :: num_str, filename, str_dummy1, str_dummy2
    real(8) :: re
    complex(kind(0d0)) :: eigenvector(dim)

    do j = fst_filter, fst_filter + num_filter - 1
      if (check_master()) then
        write(num_str, '(I0)') j
        len = len_trim(num_str)
        num_str(max_num_digits - len + 1 : max_num_digits) = num_str(1 : len)
        do digit = 1, max_num_digits - len
          num_str(digit : digit) = '0'
        end do

        filename = trim(dirname) // '/' // trim(num_str) // '.dat'
        open(iunit, file=trim(filename), status='old')
        do i = 1, dim
          ! i == k .and. l == j must be satisfied.
          read(iunit, '(I8, A1, I8, A1, E30.18e3)') k, str_dummy1, l, str_dummy2, re
          eigenvector(i) = cmplx(re, 0d0, kind(0d0))
        end do
        close(iunit)
      end if
      call mpi_bcast(eigenvector, dim, mpi_double_complex, g_wk_master_pnum, mpi_comm_world, ierr)
      do i = 1, dim
        call pzelset(Y_filtered, i, j - fst_filter + 1, Y_filtered_desc, eigenvector(i))
      end do
    end do
  end subroutine read_distribute_eigenvectors


  subroutine bcast_sparse_matrix(root, mat)
    integer, intent(in) :: root
    type(sparse_mat), intent(inout) :: mat

    integer :: my_rank, ierr
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    call mpi_bcast(mat%size, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(mat%num_non_zeros, 1, mpi_integer, root, mpi_comm_world, ierr)
    if (my_rank /= root) then
      allocate(mat%suffix(2, mat%num_non_zeros), mat%value(mat%num_non_zeros), stat=ierr)
      if (ierr /= 0) then
        call terminate('bcast_sparse_matrix: allocation failed', ierr)
      end if
    end if

    time_end = mpi_wtime()
    call add_event('bcast_sparse_matrix:bcast_aux', time_end - time_start_part)
    time_start_part = time_end

    call mpi_bcast(mat%suffix, 2 * mat%num_non_zeros, mpi_integer, root, mpi_comm_world, ierr)

    time_end = mpi_wtime()
    call add_event('bcast_sparse_matrix:bcast_suffix', time_end - time_start_part)
    time_start_part = time_end

    call mpi_bcast(mat%value, mat%num_non_zeros, mpi_double_precision, root, mpi_comm_world, ierr)

    time_end = mpi_wtime()
    call add_event('bcast_sparse_matrix:bcast_value', time_end - time_start_part)
    call add_event('bcast_sparse_matrix', time_end - time_start)
  end subroutine bcast_sparse_matrix
end module wk_distribute_matrix_m
