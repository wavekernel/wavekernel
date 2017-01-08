!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_alloc_charge
!
   use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION
   use M_io_dst_write_log, only : log_unit !(unchanged)
   use M_wall_clock_time,  only : get_system_clock_time !(routine)
   implicit none
!
   private
!
! Public routine
!
   public :: allocate_charge_arrays
   public :: allocate_csc_arrays
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocate_charge_arrays(nval_max, noav)
!       called only at the first MD step
!
     use M_qm_domain,  only : e_num_on_basis, previous_e_num_on_basis !(CHANGED)
     use M_qm_domain,  only : e_num_on_atom, previous_e_num_on_atom   !(CHANGED)
     use M_qm_domain,  only : delta_e_num                             !(CHANGED)
!
     implicit none
     integer, intent(in) :: nval_max, noav
     integer             :: ierr
     real(8)             :: memory_size
!
     if (i_verbose >= 1) then
       write(*,*)'@@ allocate_charge_arrays'
     endif   
!
     memory_size=0.0d0
!
     allocate (e_num_on_basis(nval_max, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(*,*)'alloc. error!(e_num_on_basis):ierr=',ierr
        stop
     endif
     memory_size=memory_size+8.0d0*dble(nval_max*noav)/1.0d9
!
     allocate (previous_e_num_on_basis(nval_max, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(*,*)'alloc. error!(previous_e_num_on_basis):ierr=',ierr
        stop
     endif
     memory_size=memory_size+8.0d0*dble(nval_max*noav)/1.0d9
!
     allocate (e_num_on_atom(noav),stat=ierr)
      if( ierr .ne. 0) then
          write(*,*)'alloc. error!(e_num_on_atom):ierr=',ierr
          stop
      end if
     memory_size=memory_size+8.0d0*dble(noav)/1.0d9
!
      allocate (previous_e_num_on_atom(noav),stat=ierr)
      if( ierr .ne. 0) then
         write(*,*)'alloc. error!(previous_e_num_on_atom):ierr=',ierr
         stop
      end if
     memory_size=memory_size+8.0d0*dble(noav)/1.0d9
!
     allocate (delta_e_num(noav),stat=ierr)
      if( ierr .ne. 0) then
        write(*,*)'alloc. error!(delta_e_num):ierr=',ierr
        stop
     endif
     memory_size=memory_size+8.0d0*dble(noav)/1.0d9
!
     if (i_verbose >= 1) then
       write(*,'(a,f10.5)')' allocated memory size (allocate_chage_arrays) [GB]=',memory_size
     endif   
!
   end subroutine allocate_charge_arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine allocate_csc_arrays(nval_max, noas, noav)
!       called only at the first MD step
!
     use M_qm_domain,  only : ham_tb0, ham_csc, tau_csc !(CHANGED)
     use M_qm_domain,  only : dham_tb0, fgamma_csc, gamma_csc, dgamma_csc !(CHANGED)
!
     implicit none
     integer, intent(in) :: nval_max, noas, noav
     integer             :: ierr
     real(8)             :: memory_size1, memory_size2
!
     if (i_verbose >= 1) then
       write(*,*)'@@ allocate_csc_arrays'
     endif   
!
     memory_size1=0.0d0
     memory_size2=0.0d0
!
     allocate (ham_tb0(nval_max, nval_max, noas, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(ham_tb0):ierr=',ierr
        stop
     endif
     memory_size1=memory_size1+8.0d0*dble(nval_max*nval_max*noas*noav)/1.0d9
!
     allocate (ham_csc(nval_max, nval_max, noas, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(ham_csc):ierr=',ierr
        stop
     endif
     memory_size1=memory_size1+8.0d0*dble(nval_max*nval_max*noas*noav)/1.0d9
!
     allocate (tau_csc(noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(tau_csc):ierr=',ierr
        stop
     endif
     memory_size1=memory_size1+8.0d0*dble(noav)/1.0d9
!
     allocate (dham_tb0(3,nval_max,nval_max,noas,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(dham_tb0):ierr=',ierr
        stop
     endif
     memory_size1=memory_size1+3.0d0*8.0d0*dble(nval_max*nval_max*noas*noav)/1.0d9
!
     allocate (fgamma_csc(noas,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(fgamma_csc):ierr=',ierr
        stop
     endif
     memory_size1=memory_size1+8.0d0*dble(noas*noav)/1.0d9
!
     allocate (gamma_csc(noav, noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(gamma_csc):ierr=',ierr
        stop
     endif
     memory_size2=memory_size2+8.0d0*dble(noav*noav)/1.0d9
!
     allocate (dgamma_csc(3,noav,noav),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(dgamma_csc):ierr=',ierr
        stop
     endif
     memory_size2=memory_size2+3.0d0*8.0d0*dble(noav*noav)/1.0d9
!
     if (i_verbose >= 1) then
       write(*,'(a,f10.5)')' allocated memory size 1 (allocate_csc_arrays) [GB]=',memory_size1
       write(*,'(a,f10.5)')' allocated memory size 2 (allocate_csc_arrays) [GB]=',memory_size2
     endif   
!
   end subroutine allocate_csc_arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_alloc_charge


