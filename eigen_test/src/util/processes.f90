module processes
  implicit none

  private
  public :: check_master, layout_procs

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


  subroutine check_master(solver_type, is_master)
    character(len=256), intent(in) :: solver_type
    logical, intent(out) :: is_master
    integer :: my_rank, n_procs

    select case (trim(solver_type))
    case ('lapack')
      is_master = .true.
    case default
      call blacs_pinfo(my_rank, n_procs)
      is_master = (my_rank == 0)
    end select
  end subroutine check_master
end module processes
