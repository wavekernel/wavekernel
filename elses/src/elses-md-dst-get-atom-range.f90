!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_dst_get_atom_range
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
  public :: dst_get_atm_index_range
  public :: dst_get_index_range
  public :: a_x_b_divided_by_c_i8

!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dst_get_index_range(n, index_ini, index_fin)
!
   use M_lib_dst_info,    only : myrank, nprocs  !(unchanged)
   implicit none
   integer,           intent(in)  :: n
   integer,           intent(out) :: index_ini, index_fin
   integer :: j, ierr
!
!  Note : This routine is written for 
!     index_ini=myrank*n/nprocs+1
!     index_fin=(myrank+1)*n/nprocs
!     The integer(kind=8) expression is used so as to avoide 
!          the possible overflow of the integer(kind=4) limit (= 2**31 - 1 = 2,147,483,647)
!          for the value of (myrank+1)*n 
!
   call a_x_b_divided_by_c_i8(myrank,   n, nprocs, j)
   index_ini = j +1
!
   call a_x_b_divided_by_c_i8(myrank+1, n, nprocs, j)
   index_fin = j 
!
   ierr=0
   if (index_ini < 1) ierr=1
   if (index_fin < 1) ierr=1
   if (index_ini > n ) ierr=1
   if (index_fin > n ) ierr=1
   if (index_fin < index_ini ) ierr=1
!
   if (myrank == 0) then
     if (index_ini /= 1) ierr=1
   endif
!
   if (myrank == nprocs-1) then
     if (index_fin /= n) ierr=1
   endif
!
   if (ierr == 1) then
     write(*,'(a,4i15)')'ERROR(dst_get_index_range)=',myrank, nprocs, index_ini, index_fin
     stop
   endif
!
  end subroutine dst_get_index_range
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dst_get_atm_index_range(atm_index_ini, atm_index_fin)
!
   use M_md_dst,    only : myrank, nprocs  !(unchanged)
   use M_qm_domain, only : noav            !(unchanged)
   implicit none
   integer,           intent(out) :: atm_index_ini, atm_index_fin
   integer(kind=8) :: myrank_i8, nprocs_i8, noav_i8
   integer(kind=8) :: atm_index_ini_i8, atm_index_fin_i8
!
!  Note : This routine is written for 
!     atm_index_ini=myrank*noav/nprocs+1
!     atm_index_fin=(myrank+1)*noav/nprocs
!     The integer(kind=8) expression is used so as to avoide 
!          the possible overflow of the integer(kind=4) limit (= 2**31 - 1 = 2,147,483,647)
!          for the value of myrank*noav in a huge system.
!
   integer :: ierr
!
   myrank_i8=myrank
   nprocs_i8=nprocs
   noav_i8=noav
   atm_index_ini_i8=myrank_i8*noav_i8/nprocs_i8+1
   atm_index_fin_i8=(myrank_i8+1)*noav_i8/nprocs_i8
!   
   atm_index_ini=atm_index_ini_i8
   atm_index_fin=atm_index_fin_i8
!
   ierr=0
   if (atm_index_ini < 1) ierr=1
   if (atm_index_fin < 1) ierr=1
   if (atm_index_ini > noav ) ierr=1
   if (atm_index_fin > noav ) ierr=1
   if (atm_index_fin < atm_index_ini ) ierr=1
!
   if (myrank == 0) then
     if (atm_index_ini /= 1) ierr=1
   endif
!
   if (myrank == nprocs-1) then
     if (atm_index_fin /= noav) ierr=1
   endif
!
   if (ierr == 1) then
     write(*,'(a,4i15)')'ERROR(dst_get_atm_index_range)=',myrank, nprocs, atm_index_ini, atm_index_fin
     stop
   endif
!
  end subroutine dst_get_atm_index_range
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine a_x_b_divided_by_c_i8(a, b, c, d)
!      integer(8) calculation of d = a*b/c 
!
   implicit none
   integer,           intent(in)  :: a, b, c
   integer,           intent(out) :: d
   integer(kind=8) :: a_i8, b_i8, c_i8, d_i8
!
   a_i8=a
   b_i8=b
   c_i8=c
   d_i8=a_i8*b_i8/c_i8
!
   d=d_i8
!
  end subroutine a_x_b_divided_by_c_i8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_md_dst_get_atom_range

