!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_CSC_dst
    use M_qm_domain
!
  real(DOUBLE_PRECISION), allocatable :: ham_csc_dst(:,:,:,:)
  real(DOUBLE_PRECISION), allocatable :: d_ham_csc_dst(:,:,:,:,:)
!
  private
  public :: set_CSC_parameters_geno_dst
  public :: set_hamiltonian_csc_dst
! public :: elses_get_gamma 
  public :: set_atm_force_csc_dst
! public :: set_atm_force_csc_prep
! public :: set_atm_force_csc_overlap
! public :: set_atm_force_csc_dgamma 

contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_CSC_parameters_geno_dst 
!
    use M_qm_geno_Huckel_atom_params, only : GetChemicalHardness, GetInitialENum
    use M_qm_geno_CSC               , only : elses_get_gamma
    implicit none
    integer      :: jsv
    real(kind=8) :: chemical_hardness(noav), initial_e_num(noav)
    call GetChemicalHardness(chemical_hardness)
    call GetInitialENum(initial_e_num)
    do jsv=1, noav
       tau_csc(jsv)=chemical_hardness(jsv) * 16d0/5d0
       delta_e_num(jsv)=e_num_on_atom(jsv)-initial_e_num(jsv)
    end do
    call elses_get_gamma
!
  end subroutine set_CSC_parameters_geno_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate correction term for charge selfconsistency
! Results are stored in ham_csc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_hamiltonian_csc_dst
!
    use M_qm_dst_proj_cell,  only : dst_atm_list, len_dst_atm_list !(unchanged)
    use elses_mod_noav,   only : noas !(unchanged)
    use M_qm_geno_dst,    only : ham_tot_dst, overlap_dst !(unchanged)
    implicit none
!
    integer, parameter :: ict=1
    integer :: jsv0, jsv1, jsv2
    integer :: jsd1, k
    real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz
    real(DOUBLE_PRECISION),allocatable:: GammaDotDeltaN(:)
    integer :: nval_max, ierr, dst_atm_index
    integer :: len_dst_atm_list_alloc
!
    if (i_verbose >= 1) then
       write(*,*)'@@ set_hamiltonian_csc_dst'
    endif
!
     nval_max=maxval(nval)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of ham_csc_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     len_dst_atm_list_alloc=size(ham_tot_dst,4)
     if ( allocated(ham_csc_dst)) then
       if (size(ham_csc_dst,4) /= len_dst_atm_list_alloc) then
         deallocate ( ham_csc_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (ham_csc_dst)'         
       endif
     endif
!
     if (.not. allocated(ham_csc_dst)) then
       allocate ( ham_csc_dst( nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (ham_csc_dst)'
     endif  
     ham_csc_dst(:,:,:,:)=0.0d0
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of d_ham_csc_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     len_dst_atm_list_alloc=size(ham_tot_dst,4)
     if ( allocated(d_ham_csc_dst)) then
       if (size(d_ham_csc_dst,5) /= len_dst_atm_list_alloc) then
         deallocate ( d_ham_csc_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (d_ham_csc_dst)'
       endif
     endif
!
     if (.not. allocated(d_ham_csc_dst)) then
       allocate ( d_ham_csc_dst(3, nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (d_ham_csc_dst)'
     endif
     d_ham_csc_dst(:,:,:,:,:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate(GammaDotDeltaN(noav))
!$omp parallel do
    do jsv1=1, noav
       GammaDotDeltaN(jsv1)=dot_product(gamma_csc(:,jsv1),delta_e_num(:))
    end do
!
!   do jsv1=1, noav
!     write(*,*)'GammaDotDeltaN=',jsv1,GammaDotDeltaN(jsv1)
!   enddo
!
    ham_csc=0.0d0
!
!$omp parallel do private(jsv1,k)
    do dst_atm_index=1,len_dst_atm_list(2)
       jsv2=dst_atm_list(dst_atm_index)
       do jsd1=1,njsd(jsv2,ict)
          jsv1=jsv4jsd(jsd1,jsv2)
!
!     off-diagonal elements
!
          if (jsv1 /= jsv2) then
             ! possibly jsd1==1 
             ham_csc_dst(:,:,jsd1,dst_atm_index) = overlap_dst(:,:,jsd1,dst_atm_index) * &
                  (GammaDotDeltaN(jsv1)+GammaDotDeltaN(jsv2))/2
!
!     diagonal elements (jsv1 == jsv2)
!
          else
             do k=1, size(ham_csc_dst,1) ! CAUTION, tentative code SY, Nov21, 2008
                                         ! Orbital dependence is not considered.
                ham_csc_dst(k,k,jsd1,dst_atm_index) = GammaDotDeltaN(jsv1)
             end do
          end if
       enddo ! end jsv1 loop
    enddo ! end dst_atm_index loop
    deallocate(GammaDotDeltaN)
!
    ham_tot_dst(:,:,:,:)=ham_tot_dst(:,:,:,:)+ham_csc_dst(:,:,:,:)
!
  end subroutine set_hamiltonian_csc_dst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_atm_force_csc_dst
     
    use M_qm_geno_CSC, only : set_atm_force_csc_dgamma     !(routine)
    use M_lib_mpi_wrapper,  only : mpi_wrapper_allreduce_r2 !(routine)
    implicit none
    integer :: jsv
!
    atm_force_csc(:,:)=0.0d0           !! initial clearance
    call set_atm_force_csc_prep_dst    !! calculation of fgamma_csc(noas,noav)
!
    call set_atm_force_csc_overlap_dst !! force with the derivative of overlap
    call mpi_wrapper_allreduce_r2(atm_force_csc)
    do jsv=1,noav
      if (jsv < 3) then
        write(*,*)'force_csc (after overlap)',jsv,atm_force_csc(1:3,jsv)
      endif
    enddo
!
    call set_atm_force_csc_dgamma      !! force with the derivative of gamma
    do jsv=1,noav
      if (jsv < 3) then
        write(*,*)'force_csc (after dgamma )',jsv,atm_force_csc(1:3,jsv)
      endif
    enddo
!
  end subroutine set_atm_force_csc_dst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_atm_force_csc_prep_dst ! DEAD COPY of set_atm_force_csc_prep
     
    use elses_mod_phys_const, only : para_spin_factor
      
    implicit none
       
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, jsv3, njsd1, njsd2, nval1, nval2
    integer :: jsd1, jsd2, nss1, nss2
    real(8) :: fgamma
!
    if (i_verbose >= 1) then
      write(*,*)'@@ set_hamiltonian_csc_prep_dst'
      write(*,*)'min,max of gamma =',minval(gamma_csc),maxval(gamma_csc)
      write(*,*)'min,max of delta_e_num =',minval(delta_e_num),maxval(delta_e_num)
    endif
!
    do jsv2=1,noav
       
       njsd2=njsd(jsv2,ict)
       nss2=atm_element(jsv2)
       nval2=nval(nss2)
       
       do jsd1=1,njsd2
          
          jsv1=jsv4jsd(jsd1,jsv2)
          
          njsd1=njsd(jsv1,ict)
          
          
          do jsd2=1, njsd1
             if(jsv2 == jsv4jsd(jsd2,jsv1)) exit
          end do
          
          if( jsd2 == njsd1+1 ) stop 'set_qm_force_geno: no pair'
          
             nss1=atm_element(jsv1)
             nval1=nval(nss1)
             
             if (jsv1 /=  jsv2) then

                fgamma=0.0d0
              
                do jsv3=1,noav
                   fgamma = fgamma+ &
                   (gamma_csc(jsv2,jsv3) + gamma_csc(jsv1,jsv3))*delta_e_num(jsv3)
                end do
                
                fgamma_csc(jsd1,jsv2)=fgamma

             endif
   
       enddo
    enddo   
!
    if (i_verbose >= 1) then
      write(*,*)'min,max of fgamma =',minval(fgamma_csc),maxval(fgamma_csc)
    endif
!
  end subroutine set_atm_force_csc_prep_dst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_atm_force_csc_overlap_dst
     
    use elses_mod_phys_const, only : para_spin_factor      !(parameter)
!   use M_qm_dst_dens_mat, only : dens_mat_dst             !(unchanged)
    use M_qm_dst_proj_cell,    only : dst_atm_list, len_dst_atm_list !(unchanged)
    use M_qm_geno_dst, only : d_overlap_dst                !(unchnaged)
    implicit none
       
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, jsv3, njsd1, njsd2, nval1, nval2
    integer :: jsd1, jsd2, nss1, nss2, ja1, ja2,i
    real(8) :: fgamma

    real(8) :: dvecx,dvecy,dvecz,dist1,dist2
    integer :: dst_atm_index
!
    write(*,*)'@@ set_atm_force_csc_overlap_dst(dummy routine)'
    write(*,*)'ERROR(set_atm_force_csc_overlap_dst SHOULD NOT BE CALLED) '
    stop
!
!   if (i_verbose >= 1) then
!     write(*,*)'@@ set_atm_force_csc_overlap_dst'
!     write(*,*)'min,max of fgamma =',minval(fgamma_csc),maxval(fgamma_csc)
!     write(*,*)'min,max of DM_dst =',minval(dens_mat_dst),maxval(dens_mat_dst)
!     write(*,*)'min,max of dS_dst =',minval(d_overlap_dst),maxval(d_overlap_dst)
!   endif
!
    do dst_atm_index=1,len_dst_atm_list(1)
       jsv2=dst_atm_list(dst_atm_index)
       
       njsd2=njsd(jsv2,ict)
       nss2=atm_element(jsv2)
       nval2=nval(nss2)
       
       do jsd1=1,njsd2
          
          jsv1=jsv4jsd(jsd1,jsv2)
          
          njsd1=njsd(jsv1,ict)
          
          
          do jsd2=1, njsd1
             if(jsv2 == jsv4jsd(jsd2,jsv1)) exit
          end do
          
          if( jsd2 == njsd1+1 ) stop 'set_qm_force_geno: no pair'
          
             nss1=atm_element(jsv1)
             nval1=nval(nss1)
             
             if (jsv1 /=  jsv2) then

                fgamma=fgamma_csc(jsd1,jsv2)
                
                do ja2=1, nval2
                   do ja1=1, nval1
!                     atm_force_csc(:,jsv2)=atm_force_csc(:,jsv2)-&
!&                          dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)*d_overlap_dst(:,ja1,ja2,jsd1,dst_atm_index)*fgamma
                   end do
                end do

             endif
          enddo
       enddo
       
       atm_force_csc = para_spin_factor * atm_force_csc

  end subroutine set_atm_force_csc_overlap_dst
!
end module M_qm_geno_CSC_dst
