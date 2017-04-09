module wk_distribute_matrix_m
  use mpi
  use wk_event_logger_m
  use wk_descriptor_parameters_m
  use wk_global_variables_m
  use wk_matrix_io_m
  use wk_processes_m
  use wk_util_m
  implicit none

  type wk_local_matrix_t
    real(8), allocatable :: val(:, :)
  end type wk_local_matrix_t


  ! Index range of g-th group block is
  ! [block_to_row(g) : block_to_row(g + 1) - 1, block_to_col(g) : block_to_col(g + 1) - 1]
  ! (entire matrix is [1 : m, 1 : n]).

  ! OLD INFO
  ! Blocks are assinged to processes in cyclic manner.
  ! Therefore, index conversion can be performed by the ScaLAPACK utility routines as follows.
  ! The g-th group block is stored in local_matrices(l) in p-th process,
  ! where l = indxg2l(g, 1, my_rank, 0, num_procs),
  ! and p = indxg2p(g, 1, my_rank, 0, num_procs).
  ! Eack process has numroc(num_blocks, 1, my_rank, 0, num_procs) blocks.
  ! The l-th local block corresponds to g-th global block, where g = indxl2g(l, 1, my_rank, 0, num_procs),

  ! NEW INFO
  ! Processes are divided into colors.
  ! The size of color (i.e. the number of processes in the color) should be identical among all colors.
  ! The g-th group block is stored in local_matrices(l) in c-th colored processes in column cyclic manner,
  ! where l = indxg2l(g, 1, dummy, dummy, num_colors),
  ! and c = indxg2p(g, 1, dummy, 0, num_colors).
  ! Therefore, size of local_matrices is numroc(num_blocks, 1, my_color_index, 0, num_colors),
  ! and size for rank 1 of local_matrices(l) is rows(g),
  ! and size for rank 2 of local_matrices(l) is numroc(cols(g), 1, my_rank_in_color, 0, g_wk_num_procs_per_color),
  ! where rows(g) = block_to_row(g + 1) - block_to_row(g),
  !       cols(g) = block_to_col(g + 1) - block_to_col(g),
  !   and g = indxl2g(l, 1, my_color_index, 0, num_colors).
  ! To manipulate j-th column, first, g = find_range_index(block_to_col, j) is determined.
  ! Then l and c above are determined and
  ! j_column (column index in the block) = j - block_to_col(g) + 1.
  ! if my_color_index == c and my_rank_in_color == indxg2p(j_local, 1, dummy, 0, g_wk_num_procs_per_color),
  ! local_matrices(l)%val(:, indxg2l(j_local, 1, dummy, dummy, num_procs / num_colors)) is desired physical column.
  type wk_distributed_block_diagonal_matrices_t
    ! Global values.
    ! num_procs in comm should be g_wk_num_procs_per_color * num_colors.
    integer :: m, n, num_blocks
    integer :: comm_master, color_master, num_colors
    integer, allocatable :: block_to_row(:), block_to_col(:), colors(:)
    real(8), allocatable :: dv_block_eigenvalues_full(:), dv_block_eigenvalues_filtered(:)
    ! Local values.
    integer :: my_comm, my_blacs_context
    integer :: my_color_index  ! 0-origin.
    type(wk_local_matrix_t), allocatable :: local_matrices(:)
    integer, allocatable :: descs(:, :)
  end type wk_distributed_block_diagonal_matrices_t

  integer, parameter :: index_kinds = 4

  type wk_sparse_matrix_block_with_group_t
    integer :: block_size = -1  ! block_size == -1 means unallocated
    ! indices_with_group: size is (index_kinds, block_size).
    ! For the first index, (1, 2, 3, 4) is (group col, group row, matrix row, matrix col).
    integer, allocatable :: indices_with_group(:, :)
    real(8), allocatable :: values(:)  ! Size is (block_size).
  end type wk_sparse_matrix_block_with_group_t

  private
  public :: wk_distributed_block_diagonal_matrices_t, &
       get_local_cols_wk, setup_distributed_matrix_complex, setup_distributed_matrix_real, &
       setup_distributed_matrix_with_context_real, &
       distribute_matrix_real_part, distribute_matrix_complex, &
       gather_matrix_complex, gather_matrix_real_part, gather_matrix_complex_part, &
       gather_matrix_complex_with_pzgemr2d, &
       distribute_global_sparse_matrix_wk, distribute_global_sparse_matrix_wk_part, &
       read_distribute_eigenvalues, read_distribute_eigenvectors, &
       bcast_sparse_matrix, diag_to_sparse_matrix, &
       sort_sparse_matrix_by_group, &
       initialize_distributed_block_diagonal_matrices, &
       destroy_distributed_block_diagonal_matrices, &
       check_eq_distributed_block_diagonal_matrices, &
       copy_distributed_block_diagonal_matrices, &
       copy_allocation_distributed_block_diagonal_matrices

contains

  integer function get_local_cols_wk(desc)
    integer, intent(in) :: desc(desc_size)

    integer :: numroc

    get_local_cols_wk = max(1, numroc(desc(cols_), desc(block_col_), &
         g_my_proc_col, 0, g_n_procs_col))
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


  subroutine setup_distributed_matrix_real(name, rows, cols, &
       desc, mat, is_square_block, block_size_row, block_size_col)
    character(*), intent(in) :: name
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
      call correct_block_size(name, rows, g_n_procs_row, actual_block_size_row)
    end if
    if (present(block_size_col)) then
      actual_block_size_col = block_size_col
    else
      actual_block_size_col = g_wk_block_size
      call correct_block_size(name, cols, g_n_procs_col, actual_block_size_col)
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
         g_my_proc_row, 0, g_n_procs_row))

    call descinit(desc, rows, cols, actual_block_size_row, actual_block_size_col, &
         0, 0, g_context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix_real: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols_wk(desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('setup_distributed_matrix_real: allocation failed', info)
    end if

    mat(:, :) = 0d0
  end subroutine setup_distributed_matrix_real


  subroutine setup_distributed_matrix_with_context_real(name, context, rows, cols, &
       desc, mat, is_square_block, block_size_row, block_size_col)
    character(*), intent(in) :: name
    integer, intent(in) :: context, rows, cols

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
      call correct_block_size(name, rows, g_n_procs_row, actual_block_size_row)
    end if
    if (present(block_size_col)) then
      actual_block_size_col = block_size_col
    else
      actual_block_size_col = g_wk_block_size
      call correct_block_size(name, cols, g_n_procs_col, actual_block_size_col)
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
         g_my_proc_row, 0, g_n_procs_row))

    call descinit(desc, rows, cols, actual_block_size_row, actual_block_size_col, &
         0, 0, context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix_real: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols_wk(desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('setup_distributed_matrix_real: allocation failed', info)
    end if

    mat(:, :) = 0d0
  end subroutine setup_distributed_matrix_with_context_real


  subroutine setup_distributed_matrix_complex(name, rows, cols, &
       desc, mat, is_square_block, block_size_row, block_size_col)
    character(*), intent(in) :: name
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
      call correct_block_size(name, rows, g_n_procs_row, actual_block_size_row)
    end if
    if (present(block_size_col)) then
      actual_block_size_col = block_size_col
    else
      actual_block_size_col = g_wk_block_size
      call correct_block_size(name, cols, g_n_procs_col, actual_block_size_col)
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
         g_my_proc_row, 0, g_n_procs_row))

    call descinit(desc, rows, cols, actual_block_size_row, actual_block_size_col, &
         0, 0, g_context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix_complex: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols_wk(desc)), stat = info)
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


  subroutine distribute_global_sparse_matrix_wk_part(mat_in, i, j, m, n, desc, mat)
    type(sparse_mat), intent(in) :: mat_in
    integer, intent(in) :: i, j, m, n, desc(desc_size)
    real(8), intent(out) :: mat(:, :)

    integer :: k, r, c
    real(8) :: time_start, time_end

    time_start = mpi_wtime()

    do k = 1, mat_in%num_non_zeros
      r = mat_in%suffix(1, k)
      c = mat_in%suffix(2, k)
      if (i <= r .and. r < i + m .and. j <= c .and. c < j + n) then
        call pdelset(mat, r, c, desc, mat_in%value(k))
        if (c /= r) then
          call pdelset(mat, c, r, desc, mat_in%value(k))
        endif
      end if
    enddo

    time_end = mpi_wtime()
    call add_event('distribute_global_sparse_matrix_wk_part', time_end - time_start)
  end subroutine distribute_global_sparse_matrix_wk_part


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
      if (allocated(mat%suffix)) then
        deallocate(mat%suffix)
      end if
      if (allocated(mat%value)) then
        deallocate(mat%value)
      end if
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


  ! Generate a sparse diagonal matrix X with X(i, i) = ds(i).
  subroutine diag_to_sparse_matrix(dim, ds, X)
    integer, intent(in) :: dim
    real(8), intent(in) :: ds(dim)
    type(sparse_mat), intent(out) :: X

    integer :: i

    call destroy_sparse_mat(X)
    X%size = dim
    X%num_non_zeros = dim
    allocate(X%value(dim), X%suffix(2, dim))
    X%value(:) = ds(:)
    do i = 1, dim
      X%suffix(1, i) = i
      X%suffix(2, i) = i
    end do
  end subroutine diag_to_sparse_matrix


  ! Returns whether i1 < i2.
  ! (1, 2, 3, 4) is (group col, group row, matrix row, matrix col).
  ! Supposes indices are no identical (i1 /= i2).
  logical function lt_index_with_group(i1, i2)
    integer, intent(in) :: i1(index_kinds), i2(index_kinds)

    integer :: k

    do k = 1, index_kinds
      if (i1(k) < i2(k)) then
        lt_index_with_group = .true.
        return
      else if (i1(k) > i2(k)) then
        lt_index_with_group = .false.
        return
      end if  ! If i1(1) == i2(1), go to next index.
    end do
    call terminate('lt_index_with_group: identical indices found', 1)
  end function lt_index_with_group


  subroutine get_my_initial_block(num_groups, group_indices, A, my_rank, num_procs, my_block)
    integer, intent(in) :: num_groups, group_indices(num_groups + 1), my_rank, num_procs
    type(sparse_mat), intent(in) :: A
    type(wk_sparse_matrix_block_with_group_t), intent(out) :: my_block

    integer :: rank, block_size, block_size_acc, k, r, c, gr, gc
    ! Functions.
    integer :: numroc

    my_block%block_size = numroc(A%num_non_zeros, 1, my_rank, 0, num_procs)
    allocate(my_block%indices_with_group(index_kinds, my_block%block_size))
    allocate(my_block%values(my_block%block_size))

    block_size_acc = 0
    do rank = 0, num_procs -1
      block_size = numroc(A%num_non_zeros, 1, rank, 0, num_procs)
      if (rank == my_rank) then
        do k = 1, block_size
          r = A%suffix(1, block_size_acc + k)  ! Matrix row.
          c = A%suffix(2, block_size_acc + k)  ! Matrix col.
          gr = find_range_index(group_indices, r)  ! Group row.
          gc = find_range_index(group_indices, c)  ! Group col.
          my_block%indices_with_group(1, k) = gc
          my_block%indices_with_group(2, k) = gr
          my_block%indices_with_group(3, k) = r
          my_block%indices_with_group(4, k) = c
          my_block%values(k) = A%value(block_size_acc + k)
        end do
      end if
      block_size_acc = block_size_acc + block_size
    end do
    if (block_size_acc /= A%num_non_zeros) then
      call terminate('get_my_initial_block: non-zero element count mismatch', 1)
    end if
  end subroutine get_my_initial_block


  ! Supposes block_recv is unallocated.
  subroutine send_block(num_blocks, depth, my_block, my_merge_rank, block_recv)
    integer, intent(in) :: num_blocks, depth
    type(wk_sparse_matrix_block_with_group_t), intent(in) :: my_block
    ! my_merge_rank: At the end, on ranks to send a block, becomes -1.
    !                Otherwise (i.e. on receive ranks) unchanged.
    integer, intent(inout) :: my_merge_rank
    type(wk_sparse_matrix_block_with_group_t), intent(out) :: block_recv

    integer :: block_size, rank_recv, rank_send, status(mpi_status_size), ierr, ierr2

    if (my_merge_rank == -1) then  ! Already sent my block.
      return
    end if

    if (mod(my_merge_rank, 2) == 0) then  ! This rank is receiver.
      if (my_merge_rank == num_blocks - 1) then  ! Edge rank (receives dummy empty block).
        block_recv%block_size = 0
        allocate(block_recv%indices_with_group(index_kinds, 0), block_recv%values(0))
      else
        rank_recv = g_my_rank
        rank_send = rank_recv + 2 ** (depth - 1)
        ! Share block size first.
        call mpi_recv(block_size, 1, mpi_integer, rank_send, 0, mpi_comm_world, status, ierr)
        if (ierr /= 0) then
          call terminate ('send_block: block_size receives failed', 1)
        end if
        ! Receive block body.
        block_recv%block_size = block_size
        allocate(block_recv%indices_with_group(index_kinds, block_size))
        allocate(block_recv%values(block_size))
        call mpi_recv(block_recv%indices_with_group, index_kinds * block_size, mpi_integer, &
             rank_send, 2 * depth - 1, mpi_comm_world, status, ierr)
        call mpi_recv(block_recv%values, block_size, mpi_double_precision, &
             rank_send, 2 * depth, mpi_comm_world, status, ierr)
        if (ierr /= 0 .or. ierr2 /= 0) then
          call terminate ('send_block: block body receive failed', 1)
        end if
      end if
    else  ! This rank is sender.
      rank_send = g_my_rank
      rank_recv = rank_send / (2 ** depth) * (2 ** depth)
      ! Share block size first.
      call mpi_send(my_block%block_size, 1, mpi_integer, rank_recv, 0, mpi_comm_world, ierr)
      if (ierr /= 0) then
        call terminate ('send_block: block_size send failed', 1)
      end if
      ! Send block body.
      call mpi_send(my_block%indices_with_group, index_kinds * my_block%block_size, mpi_integer, &
           rank_recv, 2 * depth - 1, mpi_comm_world, ierr)
      call mpi_send(my_block%values, my_block%block_size, mpi_double_precision, &
           rank_recv, 2 * depth, mpi_comm_world, ierr2)
      if (ierr /= 0 .or. ierr2 /= 0) then
        call terminate ('send_block: block body send failed', 1)
      end if
      my_merge_rank = -1  ! Already sent my block.
    end if
  end subroutine send_block


  subroutine destroy_block(b)
    type(wk_sparse_matrix_block_with_group_t), intent(out) :: b

    b%block_size = -1
    if (allocated(b%indices_with_group)) then
      deallocate(b%indices_with_group)
    end if
    if (allocated(b%values)) then
      deallocate(b%values)
    end if
  end subroutine destroy_block


  subroutine copy_block(orig, copy)
    type(wk_sparse_matrix_block_with_group_t), intent(in) :: orig
    type(wk_sparse_matrix_block_with_group_t), intent(out) :: copy

    call destroy_block(copy)
    copy%block_size = orig%block_size
    allocate(copy%indices_with_group(index_kinds, copy%block_size))
    allocate(copy%values(copy%block_size))
    copy%indices_with_group(:, :) = orig%indices_with_group(:, :)
    copy%values(:) = orig%values(:)
  end subroutine copy_block


  ! Merge bl into b.
  subroutine merge_block(bl, b)
    type(wk_sparse_matrix_block_with_group_t), intent(in) :: bl
    type(wk_sparse_matrix_block_with_group_t), intent(inout) :: b

    integer :: n, nl, nr, isl(index_kinds), isr(index_kinds)
    logical :: to_copy_from_left
    type(wk_sparse_matrix_block_with_group_t) :: br

    call copy_block(b, br)
    call destroy_block(b)
    b%block_size = bl%block_size + br%block_size
    allocate(b%indices_with_group(index_kinds, b%block_size))
    allocate(b%values(b%block_size))
    nl = 1
    nr = 1
    do n = 1, b%block_size
      ! Determine which block should be merged in this step.
      if (nl > bl%block_size) then
        to_copy_from_left = .false.
      else if (nr > br%block_size) then
        to_copy_from_left = .true.
      else
        isl = bl%indices_with_group(:, nl)
        isr = br%indices_with_group(:, nr)
        to_copy_from_left = lt_index_with_group(isl, isr)
      end if
      ! One step of merge.
      if (to_copy_from_left) then
        b%indices_with_group(:, n) = bl%indices_with_group(:, nl)
        b%values(n) = bl%values(nl)
        nl = nl + 1
      else
        b%indices_with_group(:, n) = br%indices_with_group(:, nr)
        b%values(n) = br%values(nr)
        nr = nr + 1
      end if
    end do
  end subroutine merge_block


  subroutine final_block_to_sparse_matrix(dim, b, A)
    integer, intent(in) :: dim
    type(wk_sparse_matrix_block_with_group_t), intent(in) :: b
    type(sparse_mat), intent(out) :: A

    integer :: n

    call destroy_sparse_mat(A)
    A%size = dim
    n = b%block_size
    A%num_non_zeros = n
    allocate(A%suffix(2, n), A%value(n))
    A%suffix(1 : 2, 1 : n) = b%indices_with_group(3 : 4, 1 : n)
    A%value(1 : n) = b%values(1 : n)
  end subroutine final_block_to_sparse_matrix


  subroutine final_block_to_block_indices(num_groups, fb, group_nonzero_indices)
    integer, intent(in) :: num_groups
    type(wk_sparse_matrix_block_with_group_t), intent(in) :: fb
    integer, intent(out) :: group_nonzero_indices(num_groups * num_groups + 1)

    integer :: i, gc, gr, gc_next, gr_next

    gc = 0
    gr = 0
    do i = 1, fb%block_size
      gc_next = fb%indices_with_group(1, i)
      gr_next = fb%indices_with_group(2, i)
      if (gc /= gc_next .or. gr /= gr_next) then
        if (gc > gc_next .or. (gc == gc_next .and. gr > gc_next)) then
          call terminate('final_block_to_block_indices: sort is incorrect', 1)
        else
          group_nonzero_indices(num_groups * gc_next + gr_next) = i
        end if
      end if
    end do
    group_nonzero_indices(num_groups * num_groups + 1) = fb%block_size + 1
  end subroutine final_block_to_block_indices


  subroutine sort_sparse_matrix_by_group(num_groups, group_indices, A, group_nonzero_indices)
    integer, intent(in) :: num_groups, group_indices(num_groups + 1)
    type(sparse_mat), intent(inout) :: A
    integer, intent(out) :: group_nonzero_indices(num_groups * num_groups + 1)

    integer :: my_merge_rank, num_blocks, depth, dim, ierr
    type(wk_sparse_matrix_block_with_group_t) :: my_block, block_recv

    my_merge_rank = g_my_rank
    num_blocks = g_n_procs

    call get_my_initial_block(num_groups, group_indices, A, g_my_rank, g_n_procs, my_block)
    depth = 1
    do while (num_blocks > 1)
      ! At the start of the do-loop, block_recv must be unallocated.
      call send_block(num_blocks, depth, my_block, my_merge_rank, block_recv)
      ! If the process received a block (including the edge block receiving dummy empty block), merge it.
      if (my_merge_rank /= -1) then
        call merge_block(block_recv, my_block)
        call destroy_block(block_recv)
        my_merge_rank = my_merge_rank / 2
      end if
      num_blocks = (num_blocks + 1) / 2
      depth = depth + 1
    end do
    if (my_merge_rank == 0) then
      dim = A%size
      call final_block_to_sparse_matrix(dim, my_block, A)
      call final_block_to_block_indices(num_groups, my_block, group_nonzero_indices)
    end if
    call bcast_sparse_matrix(0, A)
    call mpi_bcast(group_nonzero_indices, num_groups ** 2 + 1, mpi_integer, 0, mpi_comm_world, ierr)
    if (ierr /= 0) then
      call terminate('sort_sparse_matrix_by_group: mpi_bcast failed', 1)
    end if
  end subroutine sort_sparse_matrix_by_group


  subroutine initialize_distributed_block_diagonal_matrices(comm_master, m, n, &
       num_blocks, block_to_row, block_to_col, A)
    integer, intent(in) :: comm_master, m, n
    integer, intent(in) :: num_blocks, block_to_row(num_blocks + 1), block_to_col(num_blocks + 1)
    type(wk_distributed_block_diagonal_matrices_t), intent(out) :: A

    integer :: color_index, rank_range_for_color(3, 1), color_tmp, l, g, ierr
    integer :: num_procs_row_per_color, num_procs_col_per_color, rows, cols, num_local_blocks
    integer :: desc_tmp(desc_size)
    ! Functions.
    integer :: numroc, indxl2g

    A%m = m
    A%n = n
    A%num_blocks = num_blocks
    A%comm_master = comm_master
    if (mod(g_n_procs, g_wk_num_procs_per_color) /= 0) then
      call terminate('setup_distributed_matrices: the number of MPI processes must be a multiple of ' // &
           'the number of colors', 1)
    end if
    A%num_colors = g_n_procs / g_wk_num_procs_per_color
    allocate(A%block_to_row(num_blocks + 1), A%block_to_col(num_blocks + 1), A%colors(A%num_colors))
    allocate(A%dv_block_eigenvalues_full(m), A%dv_block_eigenvalues_filtered(n))
    ! Set colors and BLACS process grid.
    call mpi_comm_group(mpi_comm_world, A%color_master, ierr)
    do color_index = 0, A%num_colors - 1
      rank_range_for_color(1, 1) = g_wk_num_procs_per_color * color_index
      rank_range_for_color(2, 1) = g_wk_num_procs_per_color * (color_index + 1) - 1
      rank_range_for_color(3, 1) = 1
      call mpi_group_range_incl(A%color_master, 1, rank_range_for_color, color_tmp, ierr)
      A%colors(color_index + 1) = color_tmp
      if (rank_range_for_color(1, 1) <= g_my_rank .and. g_my_rank <= rank_range_for_color(2, 1)) then
        A%my_color_index = color_index  ! Save the color which my rank belongs to.
        ! my_comm is related to my_color.
        call mpi_comm_create(mpi_comm_world, A%colors(color_index + 1), A%my_comm, ierr)
        if (ierr /= 0) then
          call terminate('initialize_distributed_block_diagonal_matrices: mpi_comm_create failed', ierr)
        end if
        ! Get BLACS context that is related to the communicator my_comm.
        call layout_procs(g_wk_num_procs_per_color, num_procs_row_per_color, num_procs_col_per_color)
        A%my_blacs_context = A%my_comm  ! A%my_blacs_context is rewritten in blacs_gridinit below.
        call blacs_gridinit(A%my_blacs_context, 'R', num_procs_row_per_color, num_procs_col_per_color, ierr)
        if (ierr /= 0) then
          call terminate('initialize_distributed_block_diagonal_matrices: blacs_gridinit failed', ierr)
        end if
      end if
    end do
    ! Set local matrices.
    num_local_blocks = numroc(A%num_blocks, 1, A%my_color_index, 0, A%num_colors)
    allocate(A%local_matrices(num_local_blocks), A%descs(desc_size, num_local_blocks))
    do l = 1, num_local_blocks
      g = indxl2g(l, 1, A%my_color_index, 0, A%num_colors)
      rows = A%block_to_row(g + 1) - A%block_to_row(g)
      cols = A%block_to_col(g + 1) - A%block_to_col(g)
      call setup_distributed_matrix_real('Y_block', rows, cols, &
           desc_tmp, A%local_matrices(l)%val, .false., rows, 1)
      A%descs(1 : desc_size, l) = desc_tmp(1 : desc_size)
    end do
  end subroutine initialize_distributed_block_diagonal_matrices


  subroutine destroy_distributed_block_diagonal_matrices(A)
    type(wk_distributed_block_diagonal_matrices_t), intent(out) :: A

    integer :: l, color_index, ierr

    call mpi_group_free(A%color_master, ierr)
    do color_index = 0, A%num_colors - 1
      call mpi_group_free(A%colors(color_index + 1), ierr)
    end do
    deallocate(A%block_to_row, A%block_to_col, A%colors)
    deallocate(A%dv_block_eigenvalues_full, A%dv_block_eigenvalues_filtered)
    do l = 1, size(A%local_matrices, 1)
      deallocate(A%local_matrices(l)%val)
    end do
    deallocate(A%local_matrices, A%descs)
  end subroutine destroy_distributed_block_diagonal_matrices


  subroutine check_eq_distributed_block_diagonal_matrices(A, B)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: A, B

    if (B%m /= A%m .or. B%n /= A%n .or. B%num_blocks /= A%num_blocks) then
      call terminate('check_eq_distributed_block_diagonal_matrices: inconsistent matrix shape', 1)
    end if
    if (B%comm_master /= A%comm_master .or.  B%color_master /= A%color_master .or. &
         B%num_colors /= A%num_colors .or. any(B%colors(:) /= A%colors(:))) then
      call terminate('check_eq_distributed_block_diagonal_matrices: inconsistent communication setting', 1)
    end if
    if (any(B%block_to_row(:) /= A%block_to_row(:)) .or. any(B%block_to_col(:) /= A%block_to_col(:))) then
      call terminate('check_eq_distributed_block_diagonal_matrices: inconsistent block structure', 1)
    end if
  end subroutine check_eq_distributed_block_diagonal_matrices


  subroutine copy_distributed_block_diagonal_matrices(A, B)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: A
    type(wk_distributed_block_diagonal_matrices_t), intent(out) :: B

    integer :: k, l, r, c
    ! Functions.
    integer :: numroc

    call check_eq_distributed_block_diagonal_matrices(A, B)
    B%dv_block_eigenvalues_full(:) = A%dv_block_eigenvalues_full(:)
    B%dv_block_eigenvalues_filtered(:) = B%dv_block_eigenvalues_filtered(:)
    do l = 1, numroc(A%num_blocks, 1, A%my_color_index, 0, A%num_colors)
      B%local_matrices(l)%val(:, :) = A%local_matrices(l)%val(:, :)
    end do
  end subroutine copy_distributed_block_diagonal_matrices


  subroutine copy_allocation_distributed_block_diagonal_matrices(A, B)
    type(wk_distributed_block_diagonal_matrices_t), intent(in) :: A
    type(wk_distributed_block_diagonal_matrices_t), intent(out) :: B

    integer :: l, r, c, num_local_blocks

    B%m = A%m
    B%n = A%n
    B%num_blocks = A%num_blocks
    B%comm_master = A%comm_master
    B%color_master = A%color_master
    B%num_colors = A%num_colors
    allocate(B%block_to_row(A%num_blocks + 1), B%block_to_col(A%num_blocks + 1), B%colors(A%num_colors))
    allocate(B%dv_block_eigenvalues_full(A%m), B%dv_block_eigenvalues_filtered(A%n))
    num_local_blocks = size(A%local_matrices, 1)
    allocate(B%local_matrices(num_local_blocks))
    do l = 1, num_local_blocks
      r = size(A%local_matrices(l)%val, 1)
      c = size(A%local_matrices(l)%val, 2)
      allocate(B%local_matrices(l)%val(r, c))
    end do
  end subroutine copy_allocation_distributed_block_diagonal_matrices
end module wk_distribute_matrix_m
