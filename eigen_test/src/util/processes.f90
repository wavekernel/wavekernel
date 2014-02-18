module processes
  implicit none

  type process
    integer :: my_rank, n_procs, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
  end type process

  private
  public :: process, setup_distribution, get_num_procs, layout_procs, &
       print_map_of_grid_to_processes, check_master, terminate

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
      print '("BLACS process grid: ", I0, " x ", I0, " (", I0, ")")', &
           proc%n_procs_row, proc%n_procs_col, proc%n_procs
    end if

    if (proc%my_proc_row >= proc%n_procs_row .or. proc%my_proc_col >= proc%n_procs_col) then
       call blacs_exit(0)
       stop '[Warning] setup_distribution: Out of process grid, process exit'
    end if
  end subroutine setup_distribution


  subroutine get_num_procs(num_mpi_procs, num_omp_procs)
    !$ use omp_lib
    include 'mpif.h'

    integer, intent(out) :: num_mpi_procs, num_omp_procs

    integer :: ierr

    call mpi_comm_size(mpi_comm_world, num_mpi_procs, ierr)
    if (ierr /= 0) then
      write (0, *) '[Error] get_num_procs: mpi_comm_size failed, error code is: ', ierr
      stop
    end if

    num_omp_procs = 1
    !$ num_omp_procs = omp_get_max_threads()
  end subroutine get_num_procs


  subroutine layout_procs(n_procs, n_procs_row, n_procs_col)
    integer, intent(in) :: n_procs
    integer, intent(out) :: n_procs_row, n_procs_col

    integer :: n_procs_tmp, switch = 0, denom
    n_procs_tmp = n_procs

    n_procs_row = 1
    n_procs_col = 1
    do while (n_procs_tmp > 1)
       denom = 2
       do while (mod(n_procs_tmp, denom) /= 0)
          denom = denom + 1
       end do
       if (switch == 0) then
          n_procs_row = n_procs_row * denom
       else
          n_procs_col = n_procs_col * denom
       end if
       n_procs_tmp = n_procs_tmp / denom
       switch = 1 - switch
    end do
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
    include 'mpif.h'

    integer :: context, num_procs_row, num_procs_col, proc_row, proc_col
    integer :: my_rank
    integer, allocatable :: map(:, :)

    call blacs_get(-1, 0, context) ! Get default system context

    call mpi_comm_rank(mpi_comm_world, my_rank)
    if (my_rank /= 0) return

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
    include 'mpif.h'

    integer :: my_rank, ierr

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (ierr /= 0) then
      write (0, *) '[Error] check_master: mpi_comm_rank failed, error code is ', ierr
      stop
    end if

    check_master = (my_rank == 0)
  end function check_master


  subroutine terminate(err_msg)
    include 'mpif.h'

    character(*), intent(in) :: err_msg

    write (0, *) err_msg
    call mpi_finalize()
    stop
  end subroutine terminate
end module processes
