!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_charge_from_wt_dst
!
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_file_is_set, log_unit   !(unchanged)
  implicit none
  logical                     :: dst_micro_mat_is_active
!  
  private
!
! Public routines
  public calc_charge_from_wt_dst
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Setting chemical potential in the distirbuted system
!
  subroutine calc_charge_from_wt_dst(dst_atm_list, len_dst_atm_list, & 
&                               kr_dim_dst, wt_kr_dst, eig_kr_dst, xmu, e_num_on_basis_dst)
!
    use M_lib_phys_const,    only : ev4au !(unchaged)
    use M_md_dst,            only : myrank, nprocs !(unchanged)
    use M_qm_domain,         only : i_verbose, noav, atm_element, nval, & 
&                                   temp_for_electron !(unchanged)
    use M_lib_math_func,     only : Fermi_Dirac_Func !(function)
    use M_lib_mpi_wrapper,   only : mpi_wrapper_allreduce_minmax_r1, &
&                                   mpi_wrapper_allreduce_r0 !(routine)
!
    implicit none
    integer,                intent(in)  :: dst_atm_list(:)
    integer,                intent(in)  :: len_dst_atm_list(:)
    integer,                intent(in)  :: kr_dim_dst(:,:)
    real(DOUBLE_PRECISION), intent(in)  :: wt_kr_dst(:,:,:)
    real(DOUBLE_PRECISION), intent(in)  :: eig_kr_dst(:,:,:)
    real(DOUBLE_PRECISION), intent(in)  :: xmu  
    real(DOUBLE_PRECISION), intent(out) :: e_num_on_basis_dst(:,:)
!
    integer :: dst_atm_index, orb_index, atm_index, kr_dim_max
!
    real(8) :: max_energy(1)
    real(8) :: min_energy(1)
!
    real(8) :: ddemin, ddemax
    real(8) :: xtemp, xbeta
!
    integer :: ipe, npe
    integer :: iloop, nloopmax, irec
    real(8) :: xmu0, N_elec_tmp
    real(8) :: xmumn, xmumx, rRd, rEd, x, xexp, err
!
    if (i_verbose >= 1) then
      write(*,*)'@@ calc_charge_from_wt_dst'
    endif
!  
    xtemp=temp_for_electron
!        ----> Temperature in Fermi Distribution
    xbeta=1.0d0/xtemp
!
    if (i_verbose >= 1) then
      write(*,*)'xtemp [au,eV]=',xtemp,xtemp*ev4au
    endif  
!
    e_num_on_basis_dst(:,:)=0.0d0
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      OMP-NOTE : Shared scalors : xmu0,xbeta
!$omp  parallel &
!$omp& default(shared) & 
!$omp& private(ipe,npe) &
!$omp& private(dst_atm_index, atm_index, orb_index) &
!$omp& private(irec, rRd, rEd, x, xexp) &
!$omp& private(N_elec_tmp)
       ipe=0
       npe=0
!      ipe=omp_get_thread_num()+1
!      npe=omp_get_num_threads()
!      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
        do dst_atm_index=1,len_dst_atm_list(1)
          atm_index=dst_atm_list(dst_atm_index)
          if ((atm_index < 1) .or. (atm_index > noav)) then
            write(*,*)'ERROR:atm_index=',atm_index
            stop
          endif
          do orb_index=1, nval(atm_element(atm_index))
            N_elec_tmp =0.0d0
            do irec=kr_dim_dst(orb_index, dst_atm_index),1,-1
              rRd= wt_kr_dst(irec, orb_index, dst_atm_index)
              rEd=eig_kr_dst(irec, orb_index, dst_atm_index)
              x=xbeta*(rEd-xmu)
              xexp=Fermi_Dirac_Func(x)
              N_elec_tmp =  N_elec_tmp + 2.0d0*xexp*rRd
            enddo
            e_num_on_basis_dst(orb_index,dst_atm_index)=N_elec_tmp
        enddo
        enddo
!$omp end do
!$omp end parallel 
!
  end subroutine calc_charge_from_wt_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_charge_from_wt_dst


