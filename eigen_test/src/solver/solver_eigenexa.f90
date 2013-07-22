module solver_eigenexa
  !use time, only : get_wclock_time
  use distribute_matrix, only : conf_distribution, gather_matrix, allgather_row_wise, copy_global_sparse_matrix_to_local
  use matrix_io, only : sparse_mat

  implicit none

  private
  public :: eigen_solver_eigenexa
contains
  subroutine eigen_solver_eigenexa(mat, n_vec, eigenvalues, eigenvectors_global)
    use MPI
    !!use eigen_libs
    implicit none

    type(sparse_mat), intent(in) :: mat
    integer, intent(in) :: n_vec
    double precision, intent(out) :: eigenvalues(:), eigenvectors_global(:, :)

    integer :: num_nodes, num_nodes_row, num_nodes_col
    integer :: my_rank, my_node_row, my_node_col
    integer :: dim, local_size_row, local_size_col, ierr, info, context
    integer :: i, j, i_local, j_local, k, owner_node_row, owner_node_col
    double precision :: d1, d2
    double precision, allocatable :: A(:, :), Eigenvectors(:, :)
    integer :: desc(9)

    ! For descriptor
    !integer, parameter :: dlen_ = 9, dtype_ = 1, context_  = 2, m_ = 3, n_ = 4
    !integer, parameter :: mb_ = 5, nb_ = 6, rsrc_ = 7, csrc_ = 8, lld_ = 9

    !integer :: thread_level

    !call MPI_Init_thread(MPI_THREAD_MULTIPLE, thread_level, ierr)
    !call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    !call MPI_Comm_size(MPI_COMM_WORLD, num_nodes, ierr)

    !!call eigen_init(order='R')

    dim = mat%size

    !!call eigen_get_matdims(dim, local_size_row, local_size_col)

    print *, 'dim, local_size_row, local_size_col : ', dim, local_size_row, local_size_col

    allocate(A(local_size_row, local_size_col), &
         Eigenvectors(local_size_row, local_size_col), &
         stat = info)
    if (info /= 0) then
      stop "Memory exhausted"
    end if

    !!context = eigen_get_blacs_context()
    call descinit(desc, dim, dim, 64, 64, 0, 0, context, local_size_row, info)
    !call copy_global_sparse_matrix_to_local(mat, desc, A)

    !!call eigen_get_procs(num_nodes, num_nodes_col, num_nodes_row)
    !!call eigen_get_id(my_rank, my_node_col, my_node_row)

    !j_2 = eigen_loop_start( 1, num_nodes_col, my_proc_col)
    !j_3 = eigen_loop_end  ( n, num_nodes_col, my_proc_col)
    !i_2 = eigen_loop_start( 1, num_nodes_row, my_proc_row)
    !i_3 = eigen_loop_end  ( n, num_nodes_row, my_proc_row)

    !do i_1 = i_2, i_3
    !  i = eigen_translate_l2g(i_1, num_nodes_row, my_proc_row)
    !  do j_1 = j_2, j_3
    !    j = eigen_translate_l2g(j_1, num_nodes_col, my_proc_col)
    !     a(j_1, i_1) = (n+1-Max(n+1-i,n+1-j))*1.0D+00
    !  end do
    !end do

    do k = 1, mat%num_non_zeros
      i = mat%suffix(1, k)
      j = mat%suffix(2, k)
     !! owner_node_row = eigen_owner_node(i, num_nodes_row, my_node_row)
     !! owner_node_col = eigen_owner_node(j, num_nodes_col, my_node_col)
      if (my_node_row == owner_node_row .and. my_node_col == owner_node_col) then
        !!i_local = eigen_translate_g2l(i, num_nodes_row, my_node_row)
        !!j_local = eigen_translate_g2l(j, num_nodes_col, my_node_col)
        A(i_local, j_local) = mat%value(k)
        print *, 'A(', i_local, ', ', j_local, ') on (', owner_node_row, ', ', owner_node_col, '): ', mat%value(k)
      end if
      if (i /= j) then
        !!owner_node_row = eigen_owner_node(j, num_nodes_row, my_node_row)
        !!owner_node_col = eigen_owner_node(i, num_nodes_col, my_node_col)
        if (my_node_row == owner_node_row .and. my_node_col == owner_node_col) then
          !!i_local = eigen_translate_g2l(j, num_nodes_row, my_node_row)
          !!j_local = eigen_translate_g2l(i, num_nodes_col, my_node_col)
          A(i_local, j_local) = mat%value(k)
        end if
      end if
    enddo

    Eigenvectors(:, :) = 0.0

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    d1 = MPI_WTIME()
    print *, d1

    !!call eigen_sx(dim, n_vec, A, local_size_row, eigenvalues, &
    !!     Eigenvectors, local_size_row, m_forward=8, m_backward=128)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    d2 = MPI_WTIME()
    print *, d2

    !!call gather_matrix(Eigenvectors, desc, 0, 0, eigenvectors_global)

    !!call eigen_free()

    call MPI_Finalize(ierr)

    stop 'ouch'
  end subroutine eigen_solver_eigenexa
end module solver_eigenexa
