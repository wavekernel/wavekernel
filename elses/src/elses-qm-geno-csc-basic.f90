!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_CSC_basic
!
  private
  public gamma_csc_func_plot
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function gamma_csc for off site term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function gamma_csc_func_off_site(tau1,tau2,dist)

    implicit none
    real(kind(1d0)), intent(in) :: tau1, tau2, dist
    real(kind(1d0)) :: gamma_csc_func_off_site
!
    if (dist < 1.0d-10) then
      write(*,*)'ERROR(gamma_csc_func_off_sit):dist=',dist
      stop
    endif
!
    gamma_csc_func_off_site = 1.0d0/dist-s11(tau1,tau2,dist)
!
  end function gamma_csc_func_off_site
!
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
! Plot gamma CSC function for demonstration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gamma_csc_func_plot
!
    use M_config,             only : config       !(unchanged)
    use elses_mod_phys_const, only : ev4au, angst !(parameter)
    use M_qm_geno_Huckel_atom_params, only : AtomParameter !(unchanged)
    use elses_mod_file_io,    only : vacant_unit  !(function)
    implicit none
    integer :: lu, nos, n_elem_pair, iunit
    integer :: j, jmax, i, elm_index1, elm_index2, ierr
    real(kind(1d0)) :: r_min, r_max, dist, tau1, tau2, data_wrk
    real(kind(1d0)), allocatable :: value_data(:,:)
    real(kind(1d0)), allocatable :: dist_data(:)
    character(len=*),  parameter :: filename_wrk='output_csc_function_plot.txt'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input paramter for plot
!
    r_min =  1d0/angst
    r_max = 10.0d0/angst
    jmax  = 900
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (.not. config%calc%distributed%root_node) return
!
    lu=config%calc%distributed%log_unit
    if (lu > 0) then
      write(lu,'(a)') '@@@ gamma_csc_func_plot'
    else  
      return
    endif
!
    iunit=vacant_unit()
    open(iunit,file=filename_wrk,status='unknown')
    nos=config%system%structure%nelement
    n_elem_pair=nos*(nos+1)/2
    write(lu,*)'INFO:nos, n_elem_pair=', nos, n_elem_pair
!
    allocate(dist_data(0:jmax),stat=ierr)
    if (ierr /=0) stop 'Alloc error (gamma_csc_func_plot)'
!
    allocate(value_data(0:n_elem_pair, 0:jmax),stat=ierr)
    if (ierr /=0) stop 'Alloc error (gamma_csc_func_plot)'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    write(lu,*)'INFO:plot the gamma CSC function : filename=',trim(filename_wrk)
    write(iunit,'(a)') '% Plotting the gamma CSC function for element pairs'
!
!   write(*,*) 'Point 1'
!
    do elm_index1 = 1,nos
      data_wrk=Atomparameter(elm_index1)%chemical_hardness
      write(iunit,'(a,i10,a10,2f20.10)') '%   element index, element name, Chem. Hardness [au,eV]=', &
&           elm_index1, config%system%structure%velement(elm_index1)%name, &
&           Atomparameter(elm_index1)%chemical_hardness,   &
&           Atomparameter(elm_index1)%chemical_hardness*ev4au
    enddo
!
!   write(*,*) 'Point 2'
!
    do j=0, jmax
      dist_data(j)=r_min+(r_max-r_min)*dble(j)/dble(jmax)
      value_data(0,j)=1.0d0/dist_data(j)
    enddo
!
!   write(*,*) 'Point 3'
!   stop 'Stop manually'
!
    i=0
    do elm_index1 = 1,nos
      tau1=Atomparameter(elm_index1)%chemical_hardness*16.0d0/5.0d0
      do elm_index2 = elm_index1 ,nos
        tau2=Atomparameter(elm_index2)%chemical_hardness*16.0d0/5.0d0
        i=i+1
         write(iunit,'(a,3i10,2f20.10)')    '%   k: element index i, element index j, tau_i, tau_j=', & 
&                         i, elm_index1, elm_index2, tau1, tau2
        do j=0, jmax
          dist=dist_data(j)
          value_data(i,j)=gamma_csc_func_off_site(tau1,tau2,dist)
        enddo
      enddo
    enddo
!
!   write(*,*) 'Point 4'
!   stop 'Stop manually'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    dist_data(:)=dist_data(:)*angst
    value_data(:,:)=value_data(:,:)*ev4au
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   write(*,*) 'Point 5'
!   stop 'Stop manually'
!   
    write(iunit,'(a)') '%   index, distance[A], bare, gamma(i,j),.. [eV]'
!
    do j=0, jmax
      select case(n_elem_pair)
        case(3) 
          write(iunit,'(i10,5f20.10)') j, dist_data(j), value_data(0:3,j)
        case default
          write(iunit,*)               j, dist_data(j), value_data(:,j)
      end select
    enddo
!
!   write(*,*) 'Point 6'
!   stop 'Stop manually'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    close(iunit)
!
    deallocate(value_data,stat=ierr)
    if (ierr /=0) stop 'Dealloc error (gamma_csc_func_plot)'
!
!   stop 'Stop manually'
!
  end subroutine gamma_csc_func_plot
!
!
end module M_qm_geno_CSC_basic
