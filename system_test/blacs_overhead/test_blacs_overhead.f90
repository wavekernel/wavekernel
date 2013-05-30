program test_mpi_overhead
  implicit none
  include 'mpif.h'
  integer :: my_rank, n_procs, n_procs_row, n_procs_col, context

  call blacs_pinfo(my_rank, n_procs)
  call layout_procs(n_procs, n_procs_row, n_procs_col)
  call blacs_get(-1, 0, context)
  call blacs_gridinit(context, 'R', n_procs_row, n_procs_col)

  print *, 'This is rank ', my_rank
  call sleep(1)
end program test_mpi_overhead

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
