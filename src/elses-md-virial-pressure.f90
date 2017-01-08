!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_virial_pressure
!
  private
!
  public :: plot_virial_pressure
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Plot pressure by the virial theorem
!      into the main output file
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine plot_virial_pressure
!
    use M_config,         only : config       ! (unchanged)
    use M_lib_phys_const, only : convert_unit !(function)
    implicit none
    real(8)    :: pressure_value           ! total kinetic energy
    real(8)    :: pressure_component(3,2)  ! x,y,z components 
    integer    :: i_verbose 
    integer    :: lu
    integer    :: step_count
    real(8)    :: GPa_for_au               ! Constant X (1 au = X GPa)
!
    lu          = config%calc%distributed%log_unit
    i_verbose   = config%option%verbose
    step_count  = config%system%structure%mdstep
    GPa_for_au  = convert_unit('into', 'Pa')*1.0d-9
!
    if (.not. config%calc%calc_virial_pressure) then
      if (i_verbose >=1 )then
        if (lu > 0) then
          write(lu,'(a)') '@@ plot_virial_pressure (TEST)..is skipped '
        endif
      endif  
      return
    endif   
!
    if (i_verbose >=1 )then
      if (lu > 0) then
        write(lu,'(a)') '@@ plot_virial_pressure (TEST)'
        write(lu,*) ' note : Pressure unit : 1 au =', GPa_for_au, 'GPa'
      endif  
    endif  
!

!
    call calc_virial_pressure(pressure_value, pressure_component)
!
    if (lu > 0) then
      write(lu,'(a,i10,f20.10)')  'Pressure:P(total)   (GPa) = ', step_count, pressure_value*GPa_for_au
      write(lu,'(a,i10,f20.10)')  'Pressure:P_1(total) (GPa) = ', step_count, sum(pressure_component(:,1))*GPa_for_au
      write(lu,'(a,i10,f20.10)')  'Pressure:P_2(total) (GPa) = ', step_count, sum(pressure_component(:,2))*GPa_for_au
      write(lu,'(a,i10,3f20.10)') 'Pressure:P_1(x,y,z) (GPa) = ', step_count, pressure_component(:,1)*GPa_for_au
      write(lu,'(a,i10,3f20.10)') 'Pressure:P_2(x,y,z) (GPa) = ', step_count, pressure_component(:,2)*GPa_for_au
    endif  
!
  end subroutine plot_virial_pressure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculation of pressure by the virial theorem
!      Output : pressure_value [au] 
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_virial_pressure(pressure_value, pressure_component)
!
!
    use M_config,             only : config           ! (unchanged)
    use elses_mod_sim_cell,   only : noa, ax, ay, az  !(unchanged)
    use elses_mod_foi,        only : foi              !(unchanged)
    use elses_mod_tx,         only : tx, ty, tz       !(unchanged)
!
    use M_md_velocity_routines, only : calc_kinetic_energy !(routine)
!
    implicit none
!
    real(8), intent(out)           :: pressure_value           ! total kinetic energy
    real(8), intent(out), optional :: pressure_component(3,2)  ! x,y,z components 
!
    integer :: j
!
    real(8)  ddsum,ddd
    real(8)  unit_value
    real(8)  k_e_total, k_e_compo(3)
!
    real(8)  ddsum_compo(3,2)
    real(8)  dsum, dsum_check(3)
    real(8)  product_r_dot_f(3)
    real(8)  cell_volume
!
    cell_volume = ax*ay*az 
!
!   write(*,*)'@@ calc_virial_pressure'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate kinetic energy (only in the dynamics mode)
!
    k_e_total=0.0d0
    k_e_compo(:)=0.0d0
    if (trim(config%calc%mode) == 'dynamics') then
      call calc_kinetic_energy(k_e_total, k_e_compo)
!       ---> kinetic energy and its x,y,z components
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Check that the sum of the force vectors is zero
!
    do j=1,3
      dsum_check(j)=sum(foi(:,j))
    enddo 
    dsum=abs(sum(dsum_check(:)))
!
    if (dsum/dble(noa) > 1.0d-12) then
      write(*,*) 'ERROR(calc_virial_pressure):force error =',dsum/dble(noa) 
      stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@  Calculate the inner producct : \sum_i ( r_i, f_i )
!
    product_r_dot_f(1)=ax*dot_product(tx(:),foi(:,1))
    product_r_dot_f(2)=ay*dot_product(ty(:),foi(:,2))
    product_r_dot_f(3)=az*dot_product(tz(:),foi(:,3))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@  Calculate pressure components for the first term 
!           : P_1x = (2/3V) K_x 
!           : P_1y = (2/3V) K_y 
!           : P_1z = (2/3V) K_z 
!           : P_1  = P_1x + P_1y + P_1z
!
    pressure_component(:,1) = (2.0d0/3.0d0/cell_volume) * k_e_compo(:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@  Calculate pressure components for the second term 
!           : P_2x = (1/3V) \sum_i (r_i, f_i)_x
!           : P_2y = (1/3V) \sum_i (r_i, f_i)_y
!           : P_2z = (1/3V) \sum_i (r_i, f_i)_z
!           : P_2  = P_2x + P_2y + P_2z
!
    pressure_component(:,2) = (1.0d0/3.0d0/cell_volume) * product_r_dot_f(:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@  Calculate the total pressue : P = P_1 + P_2 
!
    pressure_value = sum( pressure_component(:,:) )
!
  end subroutine calc_virial_pressure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_md_virial_pressure

