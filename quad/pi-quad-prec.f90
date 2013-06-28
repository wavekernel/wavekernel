program main
  use mpi
  implicit none
  
  integer,parameter::quad_precision=selected_real_kind(30)
  real(kind=quad_precision)::pi,pi_sum
  
  integer::nprocs,myrank,ierr

  call mpi_init(ierr)
  if (ierr /= 0) then
     stop "@@ Init Error"
  end if

  call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  if (ierr /= 0) then
     stop "@@ Size Error"
  end if

  call mpi_comm_rank(mpi_comm_world,myrank,ierr)
  if (ierr /= 0) then
     stop "@@ Rank Error"
  end if

  pi = 4.0_quad_precision*atan(1.0_quad_precision)/nprocs
  write(*,*) nprocs ,myrank, pi

  call mpi_reduce(pi,pi_sum,1,mpi_real16,&
                              mpi_sum,0,mpi_comm_world,ierr)
  if (ierr /= 0) then
     stop "@@ Reduce Error"
  end if

  if (myrank .eq. 0) then
     write(*,*) "** result **",pi_sum
  end if

  call mpi_finalize(ierr)
  if (ierr /= 0) then
     stop "@@ Finalize Error"
  end if
end program main
