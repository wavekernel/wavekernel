module M_check_master
  implicit none
  private
  public :: check_master
!
contains
!
  subroutine check_master(solver_type, is_master)
    implicit none
    character(len=256), intent(in) :: solver_type
    logical, intent(out) :: is_master
    integer :: my_rank, n_procs
!
    select case (trim(solver_type))
    case ('lapack')
       is_master = .true.
    case ('scalapack_all')
       call blacs_pinfo(my_rank, n_procs)
       is_master = (my_rank == 0)
    case ('scalapack_select')
       call blacs_pinfo(my_rank, n_procs)
       is_master = (my_rank == 0)
    case default
       write(*,*) 'Error(lib_check_master): unknown solver type ', trim(solver_type)
       stop
    end select
!
  end subroutine check_master
!
end module M_check_master



