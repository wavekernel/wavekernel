!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note from 2007/07
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!0731: AMM is now set in elses_set_mass2, not here. (NT07B-119,p68)
!    : total elec. num. (eNelec) is not set here.
!    : val_elec_atm(1:nos) is now set here.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
! Adaptor module between NRL and Hoshi Order N program.
!
module MNrlAdaptor
  use param_ham_nrl
  use MNrlMatrix
  use MNrlLDos
  use MNrlRepulsive
  use flib_dom

  implicit none

  private

  public NrlReadTBParameters
  public NrlBuildMatrices
  public NrlBuildForces
  public NrlWriteLocalDos
  public NrlWritePartialDos
  public NrlWriteCohp
  public NrlBuildRepulsiveEnergy
  public NrlAddRepulsiveForce

  integer, parameter :: IO_UNIT = 83
  integer, parameter :: NRL_SUPPORTED_ATOMIC_NUMBER = 79

  interface NrlBuildRepulsiveEnergy
    module procedure BuildRepulsiveEnergy
  end interface

  interface NrlAddRepulsiveForce
    module procedure AddRepulsiveForce
  end interface

  type(TNrlMatrixBuilder) :: g_matrix_builder
  type(TNrlRepulsiveBuilder) :: g_repulsive_builder
  logical :: g_nrl_initialized = .false.

contains

  !!
  ! Reads the tight-binding parameters.
  ! Before calling this routine you must set the number of atoms
  ! and each species i.e., noav and jsei(:) in order to calculate
  ! the number of valence electrons.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlReadTBParameters(ifile_name, rNelec)
    use MNrlVariables, only: awt, noav, amm, aumass, jsei
    use MNrlIO
    use elses_mod_sim_cell, only: nos
    use elses_mod_val_elec, only: val_elec_atm
    character(len=*), intent(in) :: ifile_name
    integer :: i, the_stat
    real(8), pointer :: the_elec_count(:)
    real(8), intent(out) :: rNelec
    character(3) param_type
    logical        :: ex
    type(fnode), pointer     :: document_node

    nullify(the_elec_count)

    if (.not. g_nrl_initialized) then
      call NrlInitialize(g_matrix_builder)
      call NrlInitialize(g_repulsive_builder)
      g_nrl_initialized = .true.
    end if

!     param_type='old'
    param_type='new'
    select case(param_type) 
    case('old')
    !  old formatted parameter file
       open(IO_UNIT, FILE=ifile_name, status='OLD', iostat=the_stat)
       if(the_stat /= 0) goto 900
       call NrlReadFromUnit(g_matrix_builder, IO_UNIT, the_elec_count, awt)
       close(IO_UNIT)
    case('new')
    ! new XML-formatted parameter file 
       inquire( file=ifile_name, exist=ex )
       if( .not. ex ) then
          write(*,'(a,a)') '# Error!: nrladapt : can not open file ', &
               trim(ifile_name) 
          stop
       end if
       document_node => parsefile(ifile_name)
       call normalize(document_node)
       call NrlReadFromUnit(document_node, g_matrix_builder, ifile_name, &
            the_elec_count, awt)
    end select 

    write(*,*)'awt(1)=',awt(1)
!   stop
!
!   ! updates inverse of atomic weights.
!   do i = 1, noav
!     amm(i) = 1.0d0 / (aumass * awt(jsei(i)))
!   end do
!
!   ! updates the total number of valence electron.
!   rNelec = 0.0d0
!   do i = 1, noav
!     rNelec = rNelec + the_elec_count(jsei(i))
!   end do

    rNelec = 0.0d0
    do i=1,nos
      write(*,*) ' i,the_elec_count(i)=',the_elec_count(i)
      val_elec_atm(i)=the_elec_count(i)
      write(*,*) ' i,val_elec_atm(i)=',val_elec_atm(i)
    enddo  
!
    rcut_ham_nrl = g_matrix_builder%cut_off
    ! The below number is atomic number of element.
    i_elem_nrl = NRL_SUPPORTED_ATOMIC_NUMBER

    if(associated(the_elec_count)) deallocate(the_elec_count)
    return
  900 continue
    write(0,*) "Cannot open file: ", ifile_name
    stop
  end subroutine

  !!
  ! Build Hamiltonian and overlap matrices.
  ! Before calling this routine you must set atomic configuration
  ! by some procedures and set tight-binding parameters 
  ! by NrlReadTBParameters.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlBuildMatrices(dhij, dsij)
    real(8), intent(inout) :: dhij(:,:,:,:), dsij(:,:,:,:)
    integer :: the_err
    call NrlBuildHamiltonianMatrix(g_matrix_builder, dhij, the_err)
    if(the_err /= 0) goto 800
    
    ! Converts energy value from Ry to Ha
    dhij(:,:,:,:) = dhij(:,:,:,:) * 0.5d0

    call NrlBuildOverlapMatrix(g_matrix_builder, dsij, the_err)
    if(the_err /= 0) goto 810
    return

  800 continue
    write(0,*) "Fails in building Hamiltonian matrix"
    stop
  810 continue
    write(0,*) "Fails in building overlap matrix"
    stop
  end subroutine

  !!
  ! Calculates the force on each atom.
  ! @param foi force on ion. Note that xyz index is latter one, 
  !            that is, foi(index of atom, x-y-z). 
  ! 
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlBuildForces(foi, dbij, dpij)
    use MNrlVariables, only: noav
    real(8), intent(inout) :: foi(:,:)
    real(8), intent(in) :: dbij(:,:,:,:), dpij(:,:,:,:)
    integer :: i, the_err
    real(8) :: the_force(3)

    integer :: jsv2

    foi(:,:) = 0d0
    do jsv2 = 1, noav   
      call NrlBuildHamiltonianForceOnAtom(g_matrix_builder, dbij, jsv2, the_force, the_err)
      if(the_err /= 0) goto 800
      foi(jsv2, 1:3) = foi(jsv2, 1:3) + the_force(1:3) 
    end do

    do jsv2 = 1, noav   
      call NrlBuildOverlapForceOnAtom(g_matrix_builder, dpij, jsv2, the_force, the_err)
      if(the_err /= 0) goto 810
      ! The "two" is the convert coefficient from Ha to Ry.
      ! Note that the unit of dpij is Ha. 
      foi(jsv2, 1:3) = foi(jsv2, 1:3) + the_force(1:3) * 2d0
    end do
    ! The "two" is a spin degeneracy.
    ! The "half" is the convert coefficient from Ry to Ha.
    ! These two terms are canceled each other.
    ! foi(1:noav, 1:3) = foi(1:noav, 1:3) * 2.0d0 * 0.5d0
    return
  800 continue
    write(0,*) "Fails in building Hamiltonian force"
    stop
  810 continue
    write(0,*) "Fails in building overlap force"
    stop
  end subroutine

  !!
  ! Calculate and output LDOS from given eigenwavefunction, eigenenergy, 
  ! and overlap matrix.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlWriteLocalDos(i_unit, io_eigen, io_energies, io_overlap, &
      i_dimension, i_emin, i_emax, i_count, i_sigma, i_weights, o_err)
    integer, intent(in) :: i_unit
    real(8), intent(inout) :: io_eigen(:, :)
    real(8), intent(inout) :: io_energies(:)
    real(8), intent(inout) :: io_overlap(:, :)
    integer, intent(in) :: i_dimension
    real(8), intent(in) :: i_emin
    real(8), intent(in) :: i_emax
    integer, intent(in) :: i_count
    real(8), intent(in) :: i_sigma
    real(8), intent(in) :: i_weights(:)
    integer, intent(out) :: o_err
    type(TNrlLDosCalculator) :: calculator
    type(TNrlLDos) :: ldos

    integer :: row0

    call Initialize(ldos)
    call Initialize(calculator)
    call AssociateMatrices(calculator, io_eigen, io_energies, io_overlap, &
                            i_dimension, o_err)
    if(o_err /= 0) goto 900

    call SetSigma(calculator, i_sigma)
    call GetLDos(calculator, i_emin, i_emax, i_count, i_weights, ldos, o_err)
    if(o_err /= 0) goto 900

    call WriteTo(ldos, i_unit) 

  ! Finally
  900 continue 
    call Release(calculator)
    call Release(ldos)
    return
  end subroutine

  !!
  ! Calculate and output LDOS from given eigenwavefunction, eigenenergy, 
  ! and overlap matrix.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlWritePartialDos(i_unit, io_eigen, io_energies, io_overlap, &
      i_dimension, i_emin, i_emax, i_count, i_sigma, i_weights, o_err)
    integer, intent(in) :: i_unit
    real(8), intent(inout) :: io_eigen(:, :)
    real(8), intent(inout) :: io_energies(:)
    real(8), intent(inout) :: io_overlap(:, :)
    integer, intent(in) :: i_dimension
    real(8), intent(in) :: i_emin
    real(8), intent(in) :: i_emax
    integer, intent(in) :: i_count
    real(8), intent(in) :: i_sigma
    real(8), intent(in) :: i_weights(:)
    integer, intent(out) :: o_err
    type(TNrlLDosCalculator) :: calculator
    type(TNrlLDos) :: ldos

    real(8), allocatable :: weights(:)
    integer :: row0

    call Initialize(ldos)
    call Initialize(calculator)
    call AssociateMatrices(calculator, io_eigen, io_energies, io_overlap, &
                            i_dimension, o_err)
    if(o_err /= 0) goto 900

    call SetSigma(calculator, i_sigma)
    call GetPartialDos(calculator, i_emin, i_emax, i_count, i_weights, ldos, o_err)
    if(o_err /= 0) goto 900

    call WriteTo(ldos, i_unit) 

  ! Finally
  900 continue 
    call Release(calculator)
    call Release(ldos)
    return
  end subroutine

  !!
  ! Calculate and output LDOS from given eigenwavefunction, eigenenergy, 
  ! and overlap matrix.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine NrlWriteCohp(i_unit, io_eigen, io_energies, io_hamiltonian, &
      io_overlap, i_dimension, i_emin, i_emax, i_count, &
      i_sigma, i_weights1, i_weights2, o_err)
    integer, intent(in) :: i_unit
    real(8), intent(inout) :: io_eigen(:, :)
    real(8), intent(inout) :: io_energies(:)
    real(8), intent(inout) :: io_hamiltonian(:, :)
    real(8), intent(inout) :: io_overlap(:, :)
    integer, intent(in) :: i_dimension
    real(8), intent(in) :: i_emin
    real(8), intent(in) :: i_emax
    integer, intent(in) :: i_count
    real(8), intent(in) :: i_sigma
    real(8), intent(in) :: i_weights1(:)
    real(8), intent(in) :: i_weights2(:)
    integer, intent(out) :: o_err
    type(TNrlLDosCalculator) :: calculator
    type(TNrlLDos) :: ldos

    real(8), allocatable :: weights(:)
    integer :: row0

    call Initialize(ldos)
    call Initialize(calculator)
    call AssociateMatrices(calculator, io_eigen, io_energies, io_overlap, &
                            i_dimension, o_err)
    if(o_err /= 0) goto 900

    call SetSigma(calculator, i_sigma)
    call GetCohp(calculator, i_emin, i_emax, i_count, io_hamiltonian, &
                 i_weights1, i_weights2, ldos, o_err)
    if(o_err /= 0) goto 900

    call WriteTo(ldos, i_unit) 

  ! Finally
  900 continue 
    call Release(calculator)
    call Release(ldos)
    return
  end subroutine

  !!
  ! Adaptor procedure for calculating repulsive energy.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine BuildRepulsiveEnergy(ecc)
    use MNrlVariables, only: noav
    real(8), intent(out) :: ecc
    integer :: jsv2, the_err
    real(8) :: the_energy

    ecc = 0d0
    do jsv2 = 1, noav
      call NrlBuildRepulsiveEnergyOnAtom(g_repulsive_builder, jsv2, the_energy, the_err)
      if(the_err /= 0) goto 800
      ! The "half" is the convert coefficient from Ry to Ha.
      if (dabs(the_energy) > 1.0d-10) print *, "REPULSIVE: ", jsv2, the_energy*0.5d0
      ecc = ecc + the_energy * 0.5d0
    end do
    return
  800 continue
    write(0,*) "Fails in building repulsive energy. ABORT!"
    stop
  end subroutine

  !!
  ! Adaptor procedure for calculating repulsive energy.
  !
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine AddRepulsiveForce(foi)
    use MNrlVariables, only: noav
    real(8), intent(out) :: foi(:,:)
    integer :: jsv2, the_err
    real(8) :: the_force(3)

    do jsv2 = 1, noav
      call NrlBuildRepulsiveForceOnAtom(g_repulsive_builder, jsv2, the_force, the_err) 
      if(the_err /= 0) goto 800
      ! The "half" is the convert coefficient from Ry to Ha.
      foi(jsv2, 1:3) = foi(jsv2, 1:3) + the_force(1:3) * 0.5d0
    end do

    return

  800 continue
    write(0,*) "Fails in building repulsive forces. ABORT!"
    stop
  end subroutine

end module
