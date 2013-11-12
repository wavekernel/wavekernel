module processes
  implicit none

  private
  public :: layout_procs, check_master, terminate

contains
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

    if (check_master()) then
      write (0, *) err_msg
    end if
    call mpi_finalize()
    stop
  end subroutine terminate
end module processes
