!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_omp
!
  private
  public :: show_mpi_omp_info
  public :: show_omp_info
  public :: zero_clear_omp_r4
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Show the MPI and OMP infromation
!
  subroutine show_mpi_omp_info
    use M_qm_domain,      only : i_verbose
    use M_lib_dst_info,   only : log_unit
    use M_lib_dst_info,    only: myrank, nprocs
    implicit none
    integer ipe, npe
    integer omp_get_thread_num
    integer omp_get_num_threads
!
!
    if (log_unit <= 0) return
!
!$omp parallel default(shared) &
!$omp& private(ipe,npe)
    ipe=-1
    npe=-1
!$  ipe=omp_get_thread_num()+1
!$  npe=omp_get_num_threads()
!
    if (npe == -1) then 
      write(log_unit,'(a,2i10,a)')'INFO-MPI-OMP: P_MPI, P_OMP=', nprocs, npe+2, '(NO-OMP or OMP-STUB)'
    else  
      if (ipe == 1) write(log_unit,'(a,2i10)')  'INFO-MPI-OMP: P_MPI, P_OMP=', nprocs, npe
    endif  
!
!$omp end parallel 
!
  end subroutine show_mpi_omp_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Show the OMP infromation
!      i.e. the information of the thread number
!
  subroutine show_omp_info
    implicit none
    integer ipe, npe
    integer omp_get_thread_num
    integer omp_get_num_threads
!
    write(*,*)'@@ show_omp_info'
!
!$omp parallel default(shared) &
!$omp& private(ipe,npe)
    ipe=-1
    npe=-1
!$  ipe=omp_get_thread_num()+1
!$  npe=omp_get_num_threads()
!$  write(*,*)'--->OMP:Thread num, Total num. threads=',ipe,npe
!
    if (npe == -1) then 
      write(*,*)'----> NO-OMP or OMP-STUB'
    endif   
!
!$omp end parallel 
!
  end subroutine show_omp_info
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Zero clear an array with the OMP direvtive
!
!    NOTE : the array should be in the form of a(:, :, :, :)
!       
!
  subroutine zero_clear_omp_r4(a)
    implicit none
    real(8), intent(out) :: a(:,:,:,:)
    integer :: k
!
!$omp parallel default(shared) &
!$omp& private(k)
!$omp do schedule(static)
    do k=lbound(a,4), ubound(a,4)
      a(:,:,:,k)=0.0d0
    enddo
!$omp end do
!$omp end parallel
!
  end subroutine zero_clear_omp_r4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module M_lib_omp
