!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module M_qm_population
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_lib_dst_info,     only : log_unit !(unchanged)
!
  implicit none
!
  private
!
! Public routines
!
  public :: initial_guess_for_charge
  public :: save_population_after_convergence
  public :: update_population_guess
!
contains
!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine initial_guess_for_charge
!
    use M_config  !(unchanged : config%system%structure%vatom(:)%population_guess)
    use M_qm_domain, only : noav !(unchanged)
    use M_qm_domain, only : e_num_on_atom, previous_e_num_on_atom !(CHANGED)
    use M_structure, only : atom_type !(type)
    use M_qm_geno,   only : initialize_charge !(routine)
    implicit none
    type(atom_type), pointer :: atom
    integer :: atm_index
    integer :: n_count, i_show
!
    if (i_verbose >= 1) then
      write(*,*) '@@ initial_guess_for_charge'
    endif  
!
    if (.not. allocated(e_num_on_atom)) then
      write(*,*)'Error(initial_guess_for_charge):e_num_on_atom'
      stop
    endif
!
    if (.not. allocated(previous_e_num_on_atom)) then
      write(*,*)'Error(initial_guess_for_charge):previous_e_num_on_atom'
      stop
    endif
!
    call initialize_charge
!    ----> Set e_num_on_atom, previous_e_num_on_atom, as default
!
    i_show=10
    n_count=0
    do atm_index=1, noav
      atom => config%system%structure%vatom(atm_index)
      if (.not. atom%population_guess_set) then
        atom%population_guess     = e_num_on_atom(atm_index)  ! set as the default value
      else
        if ((atom%population_guess < -100.d00) .or. (atom%population_guess > 100.d0)) then
          write(*,*)'ERROR:population_guess=', atm_index, atom%population_guess
          stop
        endif
        n_count=n_count+1
        e_num_on_atom(atm_index)=atom%population_guess
        previous_e_num_on_atom(atm_index)=e_num_on_atom(atm_index)
      endif
      atom%population_guess_set = .true.
      if (i_verbose >= 1) then
        if (atm_index <= i_show) then
          write(*,*)' initial guess for charge=',atm_index, e_num_on_atom(atm_index), atom%population_guess
        endif
      endif
    enddo
!
    if (i_verbose >= 1) then
      write(*,*) 'INFO:initial guess charge is set to be a NON default value:n=',n_count
    endif  
!
!   stop
!
  end subroutine initial_guess_for_charge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine save_population_after_convergence
!
    use M_config  !(CHANGED : config%system%structure%vatom(:)%population)
    use M_qm_domain, only : noav !(unchanged)
    use M_qm_domain, only : e_num_on_atom, previous_e_num_on_atom !(CHANGED)
    use M_structure, only : atom_type !(type)
    use M_qm_geno,   only : initialize_charge !(routine)
    implicit none
    type(atom_type), pointer :: atom
    integer :: atm_index
    integer :: i_show
!   integer :: n_count
!
    i_show=3
!
    if (i_verbose >= 1) then
      if (log_unit > 0) write(log_unit,*) '@@ save_population_after_convergence'
    endif  
!
    if (.not. allocated(e_num_on_atom)) then
      write(*,*)'Error(initial_guess_for_charge):e_num_on_atom'
      stop
    endif
!
    if (.not. allocated(previous_e_num_on_atom)) then
      write(*,*)'Error(initial_guess_for_charge):previous_e_num_on_atom'
      stop
    endif
!
    do atm_index=1, noav
      atom => config%system%structure%vatom(atm_index)
      atom%population_set       = .true.
      atom%population=e_num_on_atom(atm_index)
      if (atm_index <= i_show) then
        if (log_unit > 0) write(log_unit,*) ' save popu=', atm_index, atom%population
      endif
      if ((atom%population < -100.d00) .or. (atom%population > 100.d0)) then
        write(*,*)'Error(population)=', atm_index, atom%population
        stop
      endif
    enddo
!   stop 'stop manually'
!
  end subroutine save_population_after_convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_population_guess
!
    use M_config  !(unchanged : config%system%structure%vatom(:)%population)
    use M_config  !(CHANGED   : config%system%structure%vatom(:)%population_guess)
    use M_qm_domain, only : noav !(unchanged)
    use M_structure, only : atom_type !(type)
    implicit none
    type(atom_type), pointer :: atom
    integer :: atm_index
    integer :: n_count
    integer :: i_show
!
    i_show=3
!
    if (i_verbose >= 1) then
      write(*,*) '@@ update_population_guess for the next CSC loop'
    endif  
!
    n_count=0
    do atm_index=1, noav
      atom => config%system%structure%vatom(atm_index)
      if (atom%population_guess_set) n_count=n_count+1
    enddo
    if (n_count /= noav) then
      write(*,*)'.... is skipped (population_guess is not defined)'
      return
    endif   
!
    if (i_verbose >= 1) then
      write(*,*) '... now the converged population will be used as the initial guess in the next CSC loop'
    endif  
!
    do atm_index=1, noav
      atom => config%system%structure%vatom(atm_index)
      if (.not. atom%population_set) then
        write(*,*)'ERROR:update_population_guess:population_set=',atom%population_set
        stop
      endif
      atom%population_guess_set       = .true.
      atom%population_guess=atom%population
      if (atm_index <= i_show) then
        if (log_unit > 0) write(log_unit,*) ' save popu_guess=', atm_index, atom%population_guess
      endif
      if ((atom%population_guess < -100d0) .or. (atom%population_guess > 100.d0)) then
        write(*,*)'Error(population_guess)=', atm_index, atom%population_guess
        stop
      endif
    enddo
!   stop 'stop manually'
!
  end subroutine update_population_guess
!
end module M_qm_population

