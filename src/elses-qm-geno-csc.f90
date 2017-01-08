!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_CSC
    use M_qm_domain
!
  private
  public :: set_CSC_parameters_geno
  public :: set_hamiltonian_csc
  public :: elses_get_gamma 
  public :: set_atm_force_csc
  public :: set_atm_force_csc_prep
  public :: set_atm_force_csc_overlap
  public :: set_atm_force_csc_dgamma 

contains
  subroutine set_CSC_parameters_geno
!  geno w.g history
!  Sep.16, 2008, S. Yamamoto, Added initialization for CSC loop
    use M_qm_geno_Huckel_atom_params, only : GetChemicalHardness, GetInitialENum
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

  end subroutine set_CSC_parameters_geno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function s11 Eq.(1.36) in ELSES Code Overview (CSC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function s11(tau1,tau2,dist)
    !
    implicit none
    real(kind=8) :: s11
    real(kind=8), intent(in) :: tau1, tau2, dist
!
    real(kind=8) :: td, tau122, tau123
!
    if(dist==0.0d0) then 
       s11=0.0d0
    else    
       if (tau1 == tau2) then
          td=tau1*dist
          s11=exp(-td)/dist*  &
               (1.d0+11.d0/16.d0*td +3.d0/16.d0*td**2 +1.d0/48.d0*td**3)
       else ! (tau1 /= tau2)
          tau122=(tau1**2-tau2**2)**2
          tau123=(tau1**2-tau2**2)**3
          S11=exp(-tau1*dist)*( tau2**4*tau1/2.0d0/tau122                  &
               - (tau2**6-3.0d0*tau2**4*tau1**2)/tau123/dist )               &
               + exp(-tau2*dist)*( tau1**4*tau2/2.0d0/tau122                 &
               + (tau1**6-3.0d0*tau1**4*tau2**2)/tau123/dist )
       end if
    end if
  end function s11
!
!
  function ds11(tau1,tau2,dist)

    implicit none
    real(8) :: ds11
    real(8), intent(in) :: tau1, tau2, dist
    real(8) :: tau122, tau123

    real(8) :: epsilon=0.000001

    if(abs(tau1-tau2)<epsilon) then
       ds11 = exp(-tau1*dist)/dist*(-1.0d0/dist - 13.0d0/8.0d0*tau1 - 0.50d0*tau1**2*dist &
            -7.0d0/48.0d0*tau1**3*dist**2 - 1.0d0/48.0d0*tau1**4*dist**3)
    else

       tau122 = (tau1**2-tau2**2)**2
       tau123 = (tau1**2-tau2**2)**3

       ds11 = exp(-tau1*dist)*(-tau2**4*tau1**2/2.0d0/tau122 + (tau1*tau2**6-3.0d0*tau2**4*tau1**3)/tau123/dist &
            + (tau2**6-3.0d0*tau2**4*tau1**2)/tau123/dist**2) &
            + exp(-tau2*dist)*(-tau1**4*tau2**2/2.0d0/tau122 - (tau2*tau1**6-3.0d0*tau1**4*tau2**3)/tau123/dist &
            - (tau1**6-3.0d0*tau1**4*tau2**2)/tau123/dist**2)
    end if


  end function ds11

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gamma in Charge SelfConsistency (CSC) loop
! Results are stored in gamma_csc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine elses_get_gamma
    implicit none
!
!!$    interface
!!$       function s11(tau1,tau2,dist)
!!$         implicit none
!!$         real(kind=8) :: s11
!!$         real(kind=8), intent(in) :: tau1, tau2, dist
!!$       end function s11
!!$    end interface
!
    integer :: jsv1, jsv2, i
    real(kind=8) :: dvecx, dvecy, dvecz
    real(kind=8) :: tau1, tau2, dist
!
    if (i_verbose >= 1) then
       write(*,*) '@@ elses_get_gamma'
    endif
!
    do jsv2=1,noav
!
       tau2 = tau_csc(jsv2)
!
       do jsv1=1,jsv2
          tau1 = tau_csc(jsv1)
!
          if (jsv1 /=  jsv2) then
             dvecx=atm_position(1,jsv2)-atm_position(1,jsv1)
             dvecy=atm_position(2,jsv2)-atm_position(2,jsv1)
             dvecz=atm_position(3,jsv2)-atm_position(3,jsv1)
!
             if (i_pbc_x == 1) dvecx = modulo(dvecx+0.5d0,1.0d0)-0.5d0
             if (i_pbc_y == 1) dvecy = modulo(dvecy+0.5d0,1.0d0)-0.5d0
             if (i_pbc_z == 1) dvecz = modulo(dvecz+0.5d0,1.0d0)-0.5d0
!
             dvecx=dvecx*ax
             dvecy=dvecy*ay
             dvecz=dvecz*az
!
!       Vector R2 - R1 : (dvecx, dvecy, dvecz)
!
             dist=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
!       Distance | R2 - R1 |
!
             gamma_csc(jsv1,jsv2) = 1.0d0/dist-s11(tau1,tau2,dist)
             gamma_csc(jsv2,jsv1) = gamma_csc(jsv1,jsv2)

             dgamma_csc(1,jsv1,jsv2) = dvecx/dist**3 + (dvecx/dist)*ds11(tau1,tau2,dist)
             dgamma_csc(2,jsv1,jsv2) = dvecy/dist**3 + (dvecy/dist)*ds11(tau1,tau2,dist)
             dgamma_csc(3,jsv1,jsv2) = dvecz/dist**3 + (dvecz/dist)*ds11(tau1,tau2,dist)
             dgamma_csc(:,jsv2,jsv1) = -dgamma_csc(:,jsv1,jsv2)      

          else   !  jsv1 == jsv2
             gamma_csc(jsv1,jsv2) = 5.0d0/16.0d0*tau1

             do i=1,3
                dgamma_csc(i,jsv1,jsv2)=0.0d0
             end do
             
          endif
!
       enddo ! end jsv1 loop
    enddo ! end jsv2 loop
!
  end subroutine elses_get_gamma
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate correction term for charge selfconsistency
! Results are stored in ham_csc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_hamiltonian_csc
!
    implicit none
!
    integer, parameter :: ict=1
    integer :: jsv0, jsv1, jsv2
    integer :: jsd1, k
    real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz
    real(DOUBLE_PRECISION),allocatable:: GammaDotDeltaN(:)
!
    if (i_verbose >= 1) then
       write(*,*)'@@ elses_get_ham_csc'
    endif
!
    allocate(GammaDotDeltaN(noav))
!$omp parallel do
    do jsv1=1, noav
       GammaDotDeltaN(jsv1)=dot_product(gamma_csc(:,jsv1),delta_e_num(:))
    end do

    ham_csc=0.0d0
!
!$omp parallel do private(jsv1,k)
    do jsv2=1,noav
       do jsd1=1,njsd(jsv2,ict)
          jsv1=jsv4jsd(jsd1,jsv2)
!
!     off-diagonal elements
!
          if (jsv1 /= jsv2) then
             ! possibly jsd1==1 
             ham_csc(:,:,jsd1,jsv2) = dsij(:,:,jsd1,jsv2) * &
                  (GammaDotDeltaN(jsv1)+GammaDotDeltaN(jsv2))/2
!
!     diagonal elements (jsv1 == jsv2)
!
          else
             do k=1, size(ham_csc,1) ! CAUTION, tentative code SY, Nov21, 2008
                                     ! Orbital dependence is not considered.
                ham_csc(k,k,jsd1,jsv2) = GammaDotDeltaN(jsv1)
             end do
          end if
       enddo ! end jsv1 loop
    enddo ! end jsv2 loop
    deallocate(GammaDotDeltaN)
!
  end subroutine set_hamiltonian_csc

  subroutine set_atm_force_csc
     
    implicit none
       
    atm_force_csc(:,:)=0.0d0       !! initial clearance
    call set_atm_force_csc_prep    !! calculation of fgamma_csc(noas,noav)
    call set_atm_force_csc_overlap !! force with the derivative of overlap
    call set_atm_force_csc_dgamma  !! force with the derivative of gamma

  end subroutine set_atm_force_csc

  subroutine set_atm_force_csc_prep
     
    use elses_mod_phys_const, only : para_spin_factor
      
    implicit none
       
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, jsv3, njsd1, njsd2, nval1, nval2
    integer :: jsd1, jsd2, nss1, nss2
    real(8) :: fgamma


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

  end subroutine set_atm_force_csc_prep

  subroutine set_atm_force_csc_overlap
     
    use elses_mod_phys_const, only : para_spin_factor
      
    implicit none
       
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, jsv3, njsd1, njsd2, nval1, nval2
    integer :: jsd1, jsd2, nss1, nss2, ja1, ja2,i
    real(8) :: fgamma

    real(8) :: dvecx,dvecy,dvecz,dist1,dist2

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

                fgamma=fgamma_csc(jsd1,jsv2)
                
                do ja2=1, nval2
                   do ja1=1, nval1
                                        
                      atm_force_csc(:,jsv2)=atm_force_csc(:,jsv2)-&
                           dbij(ja1,ja2,jsd1,jsv2)*ddsij(:,ja1,ja2,jsd1,jsv2)*fgamma
                      
                   end do
                end do

             endif
          enddo
       enddo
       
       atm_force_csc = para_spin_factor * atm_force_csc

  end subroutine set_atm_force_csc_overlap

  subroutine set_atm_force_csc_dgamma
     
    implicit none
       
    integer, parameter :: ict=1
    integer :: jsv1, jsv2
    real(kind=8),allocatable :: fdgamma(:,:)
       
    allocate(fdgamma(3,noav))
    
    !!!!! force term of dgamma !!!!!!!
    do jsv2=1,noav
       
       fdgamma=0.0d0
         
       do jsv1=1,noav
             
          fdgamma(:,jsv2) = fdgamma(:,jsv2) + dgamma_csc(:,jsv2,jsv1)*delta_e_num(jsv1)
          
       enddo
          
       atm_force_csc(:,jsv2) = atm_force_csc(:,jsv2) - delta_e_num(jsv2)*fdgamma(:,jsv2)
          
    enddo
       !
    deallocate(fdgamma)

  end subroutine set_atm_force_csc_dgamma

end module M_qm_geno_CSC
