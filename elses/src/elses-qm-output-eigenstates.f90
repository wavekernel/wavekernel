!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_output_eigenstates

  use M_qm_domain,        only : i_verbose
  implicit none
  integer :: unit_num

  public :: output_eigenstates

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  subroutine output_eigenstates(filename, no_wavefunction)

    use M_config, only : config
    use elses_mod_file_io,      only : vacant_unit
!   use M_io_ctrl_output_eigen, only : init_prep_for_plot_eigenstates

    implicit none
    integer :: i, imode
    character(len=*), intent(in) :: filename
    logical,          intent(in), optional :: no_wavefunction
    logical                                :: no_wavefunction_wrk
!
    no_wavefunction_wrk = .false.
    if (present(no_wavefunction)) no_wavefunction_wrk=no_wavefunction
!   
    if (i_verbose >= 1) then
      write(*,*)'@@ output_eigenstates'
      write(*,*)'    filename=',trim(filename)
    endif   

!   call init_prep_for_plot_eigenstates(imode)
!
!   if (imode == 0) then 
!      if (i_verbose >= 1) then
!        write(*,*)' skipped: output_eigenstates'
!      endif   
!      return
!   endif   
!
    if (trim(filename) == '') then
      write(*,*)' ERROR!:output_eigenstates: empty filename'
      stop
    endif   
!
    unit_num=vacant_unit()
    open(unit_num,file=trim(filename))
!
    write(unit_num,'(a)')'# file_format= v0.04.05'
!
    call output_number_of_atom

    call output_cell_information

    do i=1,config%system%structure%natom

       call output_atom_species(i)

       call output_atom_parameter(i)

       call output_atom_position(i)

    end do
    
    if (no_wavefunction_wrk) then
       write(unit_num,'(a)') '# LCAO coefficients .. are not included.'
    else
      call output_coefficient
    endif  

    close(unit_num)    

  end subroutine output_eigenstates
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output number of atom
!
  subroutine output_number_of_atom

    use M_config, only : config
    use elses_mod_orb2, only : n_tot_base

    implicit none

!    write(*,'("number of atom =",2I10)') config%system%structure%natom, n_tot_base

    write(unit_num,'(2I10)') config%system%structure%natom, n_tot_base
    
  end subroutine output_number_of_atom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output unit cell vector
!
  subroutine output_cell_information

    use M_config, only : config
    use M_qm_domain, only : i_pbc_x, i_pbc_y, i_pbc_z !(unchanged)

    implicit none

!    write(*,'("vectorA:",3F16.12)') config%system%structure%unitcell%vectorA(1:3)
!    write(*,'("vectorB:",3F16.12)') config%system%structure%unitcell%vectorB(1:3)
!    write(*,'("vectorC:",3F16.12)') config%system%structure%unitcell%vectorC(1:3)

    write(unit_num,'(3F16.8,i5)') config%system%structure%unitcell%vectorA(1:3), i_pbc_x
    write(unit_num,'(3F16.8,i5)') config%system%structure%unitcell%vectorB(1:3), i_pbc_y
    write(unit_num,'(3F16.8,i5)') config%system%structure%unitcell%vectorC(1:3), i_pbc_z


  end subroutine output_cell_information
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  subroutine output_atom_species(atom_index)

    use M_config, only : config
    use M_md_save_struct, only : get_atomic_num

    implicit none

    integer, intent(in) :: atom_index

    integer :: i, atom_num

    atom_num=-1

    call get_atomic_num(trim(config%system%structure%vatom(atom_index)%name),atom_num)

!    write(*,'(" atom element: ",A2, "  atomic number: ",I2)') &
!         config%system%structure%vatom(atom_index)%name, atom_num

    write(unit_num,'(I2)') atom_num
        
  end subroutine output_atom_species
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  subroutine output_atom_parameter(atom_index)

    use M_config, only : config
    use M_qm_domain, only :atm_element
    use M_qm_geno_Huckel_atom_params, only : AtomParameter

    implicit none

    integer, intent(in) :: atom_index
    integer :: i,j,k
    character(len=50) :: fmt, fmt_t

    if (.not. allocated(AtomParameter) ) then
      write(*,*)'ERROR(output_atom_parameter):AtomParameter is not allocated'
      stop
    endif

    j=AtomParameter(atm_element(atom_index))%num_val_orb

    if (i_verbose >= 10) then    
       write(*,'(" Total Number of atomic orbitals =",I8)')  j
       
       write(fmt,'(I1,"A6")') j
       fmt_t='" Orbital         =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt) (AtomParameter(atm_element(atom_index))%orbital(k),k=1, j)
       write(fmt,'(I1,"I6")') j
       
       fmt_t='" Principal       =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt)   AtomParameter(atm_element(atom_index))%principal(1:j)
       write(fmt,'(I1,"I6")') j
       
       fmt_t='" Angular         =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt)   AtomParameter(atm_element(atom_index))%angular(1:j)
       
       write(fmt,'(I1,"ES10.3")') j
       fmt_t='" zeta            =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt)   AtomParameter(atm_element(atom_index))%zeta(1:j)
       
       write(fmt,'(I1,"ES10.3")') j
       fmt_t='" zeta2           =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt)   AtomParameter(atm_element(atom_index))%zeta2(1:j)
       
       write(fmt,'(I1,"ES10.3")') j
       fmt_t='" c1              =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt)   AtomParameter(atm_element(atom_index))%c1(1:j)
       
       write(fmt,'(I1,"ES10.3")') j
       fmt_t='" c2              =",'
       fmt="(" // trim(fmt_t) // trim(fmt) // ")"
       write(*,fmt)   AtomParameter(atm_element(atom_index))%c2(1:j)
    endif


    write(unit_num,'(I8)')  j
 
    write(fmt,'(I1,"I6")') j
    fmt="(" // trim(fmt) // ")"
    write(unit_num,fmt)   AtomParameter(atm_element(atom_index))%principal(1:j)

    write(fmt,'(I1,"I6")') j
    fmt="(" // trim(fmt) // ")"
    write(unit_num,fmt)   AtomParameter(atm_element(atom_index))%angular(1:j)

    write(fmt,'(I1,"ES10.3")') j
    fmt="(" // trim(fmt) // ")"
    write(unit_num,fmt)   AtomParameter(atm_element(atom_index))%zeta(1:j)

    write(fmt,'(I1,"ES10.3")') j
    fmt="(" // trim(fmt) // ")"
    write(unit_num,fmt)   AtomParameter(atm_element(atom_index))%zeta2(1:j)

    write(fmt,'(I1,"ES10.3")') j
    fmt="(" // trim(fmt) // ")"
    write(unit_num,fmt)   AtomParameter(atm_element(atom_index))%c1(1:j)

    write(fmt,'(I1,"ES10.3")') j
    fmt="(" // trim(fmt) // ")"
    write(unit_num,fmt)   AtomParameter(atm_element(atom_index))%c2(1:j)
 
  end subroutine output_atom_parameter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  subroutine output_atom_position(atom_index)

    use M_config, only : config

    implicit none
    integer :: i
    integer, intent(in) :: atom_index

    write(unit_num,'(3F23.15)') config%system%structure%vatom(atom_index)%position(1:3)
    
  end subroutine output_atom_position
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine output_coefficient

    use elses_mod_orb2, only : n_tot_base
    use elses_arr_eig_leg, only : atmp

    implicit none
    integer :: n_eigen_vectors
    integer :: i,j
!
    n_eigen_vectors = n_tot_base
!
    write(unit_num,'(a)') '# LCAO coefficients'
    write(unit_num,'(I10)') n_eigen_vectors
!
    if (allocated(atmp)  .eqv. .false.)  then
      write(*,*) 'ERROR(output_coefficient): ATMP is not allocated'
      stop
    endif   

    do i=1, n_tot_base
       do j=1, n_tot_base
          write(unit_num,'(I10,F25.15,a,I10)') j, atmp(j,i), '   k=', i
       end do
    end do

  end subroutine output_coefficient
!

end module M_output_eigenstates

