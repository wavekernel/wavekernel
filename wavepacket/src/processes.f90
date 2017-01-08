module wp_processes_m
  use mpi
  use wp_global_variables_m
  implicit none

  type wp_process_t
    integer :: my_rank, n_procs, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: n_omp_threads
  end type wp_process_t

  private
  public :: wp_process_t, setup_distribution, print_proc, get_num_procs, layout_procs, &
       print_map_of_grid_to_processes, check_master, terminate, &
       check_nan_scalar, check_nan_vector, check_nan_matrix

contains

  subroutine setup_distribution(proc)
    type(wp_process_t), intent(out) :: proc

    call get_num_procs(proc%n_procs, proc%n_omp_threads)
    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call layout_procs(proc%n_procs, proc%n_procs_row, proc%n_procs_col)
    call blacs_get(-1, 0, proc%context)
    call blacs_gridinit(proc%context, 'R', proc%n_procs_row, proc%n_procs_col)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    if (check_master()) then
      write (0, '(A, F16.6, A, I0, " x ", I0, " (", I0, ")" )') &
           ' [Event', mpi_wtime() - g_wp_mpi_wtime_init, "] BLACS process grid: ", &
           proc%n_procs_row, proc%n_procs_col, proc%n_procs
    end if

    if (proc%my_proc_row >= proc%n_procs_row .or. proc%my_proc_col >= proc%n_procs_col) then
       call blacs_exit(0)
       stop '[Warning] setup_distribution: Out of process grid, process exit'
    end if
  end subroutine setup_distribution


  subroutine print_proc(proc)
    type(wp_process_t), intent(in) :: proc

    print *, 'num_mpi_processes: ', proc%n_procs
    print *, 'num_mpi_processes_row: ', proc%n_procs_row
    print *, 'num_mpi_processes_col: ', proc%n_procs_col
    print *, 'num_omp_threads: ', proc%n_omp_threads
  end subroutine print_proc


  subroutine get_num_procs(num_mpi_procs, num_omp_procs)
    !$ use omp_lib
    integer, intent(out) :: num_mpi_procs, num_omp_procs

    integer :: ierr

    call mpi_comm_size(mpi_comm_world, num_mpi_procs, ierr)
    if (ierr /= 0) then
      call terminate('get_num_procs: mpi_comm_size failed', ierr)
    end if

    num_omp_procs = 1
    !$ num_omp_procs = omp_get_max_threads()
  end subroutine get_num_procs


  subroutine layout_procs(n_procs, n_procs_row, n_procs_col)
    integer, intent(in) :: n_procs
    integer, intent(out) :: n_procs_row, n_procs_col

    n_procs_row = int(sqrt(dble(n_procs + 1)))
    do while (mod(n_procs, n_procs_row) /= 0)
      n_procs_row = n_procs_row - 1
    end do
    n_procs_col = n_procs / n_procs_row
  end subroutine layout_procs


  subroutine map_grid_to_processes(context, num_procs_row, num_procs_col, map)
    integer, intent(in) :: context, num_procs_row, num_procs_col
    integer, intent(out) :: map(num_procs_row, num_procs_col)

    integer :: i, j
    integer :: blacs_pnum ! function

    do j = 0, num_procs_col - 1
      do i = 0, num_procs_row - 1
        map(i + 1, j + 1) = blacs_pnum(context, i, j)
      end do
    end do
  end subroutine map_grid_to_processes


  subroutine print_map_of_grid_to_processes()
    integer :: context, num_procs_row, num_procs_col, proc_row, proc_col
    integer :: my_rank, ierr
    integer, allocatable :: map(:, :)

    call blacs_get(-1, 0, context) ! Get default system context

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (my_rank /= g_wp_master_pnum) return

    call blacs_gridinfo(context, num_procs_row, num_procs_col, proc_row, proc_col)
    allocate(map(num_procs_row, num_procs_col))
    call map_grid_to_processes(context, num_procs_row, num_procs_col, map)

    print '("process numbers in BLACS grid is")'
    do proc_row = 1, num_procs_row
      do proc_col = 1, num_procs_col - 1
        write (*, '(i6, " ")', advance = 'no') map(proc_row, proc_col)
      end do
      print '(i6)', map(proc_row, num_procs_col)
    end do
  end subroutine print_map_of_grid_to_processes


  logical function check_master()
    integer :: my_rank, ierr

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (ierr /= 0) then
      call terminate('check_master: mpi_comm_rank failed', ierr)
    end if

    check_master = (my_rank == g_wp_master_pnum)
  end function check_master


  subroutine terminate(error_message, error_code)
    character(*), intent(in) :: error_message
    integer, intent(in) :: error_code

    integer :: ierr

    write (0, '("[Error] ", a)') error_message
    call mpi_abort(mpi_comm_world, error_code, ierr)
  end subroutine terminate


  subroutine check_nan_scalar(name, x)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: x
    if (x /= x) then
      call terminate('check_nan_scalar: detect NaN in ' // name, 1)
    end if
  end subroutine check_nan_scalar


  subroutine check_nan_vector(name, xs)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: xs(:)
    integer :: i
    do i = 1, size(xs, 1)
      if (xs(i) /= xs(i)) then
        print *, '[Error] NaN at ', i
        call terminate('check_nan_vector: detect NaN in ' // name, 1)
      end if
    end do
  end subroutine check_nan_vector


  subroutine check_nan_matrix(name, xss)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: xss(:, :)
    integer :: i, j
    do j = 1, size(xss, 2)
      do i = 1, size(xss, 1)
        if (xss(i, j) /= xss(i, j)) then
          print *, '[Error] NaN at ', i, j
          call terminate('check_nan_matrix: detect NaN in ' // name, 1)
        end if
      end do
    end do
  end subroutine check_nan_matrix
end module wp_processes_m
