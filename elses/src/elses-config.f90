!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module M_element
  type :: classic_type
     real(8) :: mass
     real(8) :: charge
  end type classic_type

  type :: quantum_type
     character(len=64)  :: type
     integer :: orbital
     real(8) :: interaction_radius
     real(8) :: DNAL0
     real(8) :: RNN0
     real(8), dimension(4) :: DHAL
     real(8), dimension(4) :: DNAL
     real(8), dimension(4) :: RCAL
     real(8) :: restpart
  end type quantum_type

  type :: element_type
     character(len=64)  :: filename
     character(len=8)   :: name
     type(classic_type) :: classic
     type(quantum_type) :: quantum
  end type element_type
contains
  ! print <element> data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine element_print( element )
    implicit none
    type(element_type) :: element

    write(*,*) "element%name = ", element%name

    call classic_print( element%classic )
    call quantum_print( element%quantum )

    return
  end subroutine element_print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print <classic> data
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine classic_print( classic )
    implicit none
    type(classic_type) :: classic

    write(*,*) "classic%mass = ",   classic%mass
    write(*,*) "classic%charge = ", classic%charge

    return
  end subroutine classic_print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print <quantum> data
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine quantum_print( quantum )
    implicit none
    type(quantum_type) :: quantum

    write(*,*) "quantum%orbital=",  quantum%orbital
    write(*,*) "quantum%DNAL0=",    quantum%DNAL0
    write(*,*) "quantum%RNN0=",     quantum%RNN0
    write(*,*) "quantum%DHAL=",     quantum%DHAL(1:4)
    write(*,*) "quantum%DNAL=",     quantum%DNAL(1:4)
    write(*,*) "quantum%RCAL=",     quantum%RCAL(1:4)
    write(*,*) "quantum%restpart=", quantum%restpart

    return
  end subroutine quantum_print
end module M_element
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module M_structure
  use M_element

  type :: unitcell_type
     logical               :: set
     real(8), dimension(3) :: vectorA
     real(8), dimension(3) :: vectorB
     real(8), dimension(3) :: vectorC
     real(8)               :: myLength
  end type unitcell_type

  type :: heatbath_type
     real(kind(1d0)) :: massperatom
     real(kind(1d0)) :: position
     real(kind(1d0)) :: velocity
  end type heatbath_type

  type :: custumize_type
     real(8), dimension(4) :: base
     character(len=16)     :: status
  end type custumize_type

  type :: mini_atom_type
     character(len=4)      :: name
     real(8), dimension(3) :: position
  end type mini_atom_type

  type :: atom_type
     character(len=8)      :: name
     type(element_type),pointer :: element
     character(len=16)     :: class
     character(len=8)      :: motion
     real(8), dimension(3) :: position
     real(8), dimension(3) :: velocity
     logical               :: velocity_set
     real(8), dimension(3) :: force
     logical               :: force_set
     real(8)               :: population
     logical               :: population_set
     real(8)               :: population_guess
     logical               :: population_guess_set
     integer               :: group_id
     logical               :: group_id_set
     integer                                     :: ncustumize
     type(custumize_type), pointer, dimension(:) :: vcustumize
  end type atom_type

  type :: split_type
     logical :: set
     integer :: file_index
     integer :: number_of_files
     integer :: atoms_in_this_file
     integer :: atom_initial
     integer :: atom_final
  end type split_type

  type :: structure_type
     character(len=64)    :: name
     integer              :: mdstep
     character(len=8)     :: parser ! "sax", "SAX", "dom" or "DOM"
     logical              :: tag_dump
     character(len=32)    :: read_mode ! "redundant", "root_only", "split" or "default(=redundant)"
     type(unitcell_type)  :: unitcell
     type(heatbath_type)  :: heatbath
     integer                              :: nelement
     type(element_type), pointer, dimension(:) :: velement
     integer                              :: natom
     type(atom_type),pointer,dimension(:) :: vatom
     type(split_type)     :: split
     logical              :: use_matom
     logical              :: use_vatom
     type(mini_atom_type),pointer,dimension(:) :: matom
  end type structure_type


contains

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine heatbath_default( heatbath )
    implicit none
    type(heatbath_type), intent(out) :: heatbath

    heatbath%massperatom = 25.0d0
    heatbath%position = 0.0d0
    heatbath%velocity = 0.0d0

    return
  end subroutine heatbath_default

  ! print <structure> data
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine structure_print( structure )
    implicit none
    type(structure_type) :: structure
    integer :: j

    write(*,*) "name   = ", structure%name
    write(*,*) "mdstep = ", structure%mdstep

    call unitcell_print( structure%unitcell )

    do j=1, structure%natom
       write(*,*) "atom array index = ", j
       call atom_print( structure%vatom(j) )
    end do

    return
  end subroutine structure_print


  ! print <unitcell> data
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine unitcell_print( unitcell )
    implicit none
    type(unitcell_type) :: unitcell

    write(*,*) "unitcell%vectorA = ", unitcell%vectorA(1:3)
    write(*,*) "unitcell%vectorB = ", unitcell%vectorB(1:3)
    write(*,*) "unitcell%vectorC = ", unitcell%vectorC(1:3)
    write(*,*) "unitcell%myLength =", unitcell%myLength

    return
  end subroutine unitcell_print

  ! print <atom> data
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine atom_print( atom )
    implicit none
    type(atom_type) :: atom

    if( associated(atom%element) ) call element_print( atom%element )
    write(*,*) "atom%position = ", atom%position(1:3)
    write(*,*) "atom%velocity = ", atom%velocity(1:3)
    write(*,*) "atom%class    = ", atom%class
    write(*,*) "atom%motion   = ", atom%motion

    return
  end subroutine atom_print

end module M_structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module M_config
  use M_structure

  type :: target_type
     logical             :: class_set
     logical             :: xrange_set
     logical             :: yrange_set
     logical             :: zrange_set
     logical             :: trange_set

     character(len=64)   :: class
     real(8),dimension(2):: xrange ! xmin and xmax
     real(8),dimension(2):: yrange ! ymin and ymax
     real(8),dimension(2):: zrange ! zmin and zmax
     real(8),dimension(2):: trange ! tmin and tmax

     logical             :: temperature_set
     logical             :: mass_set
     logical             :: position_set
     logical             :: velocity_set
     logical             :: motion_set

     real(8)             :: temperature
     real(8)             :: mass
     real(8),dimension(3):: position
     real(8),dimension(3):: velocity
     character(len=8)    :: motion
  end type target_type

  type :: boundary_type
     logical             :: periodic_x
     logical             :: periodic_y
     logical             :: periodic_z
  end type boundary_type

  type :: system_type
     real(8)             :: temperature
     real(8)             :: cutoff_radius ! SY Nov21,2008
     type(structure_type):: structure
     type(boundary_type) :: boundary
     integer             :: ntarget
     type(target_type),pointer,dimension(:) :: vtarget
  end type system_type


  type :: limit_type
     real(8)             :: time
     real(8)             :: memory            ! in GB
     real(8)             :: memory_allocated  ! in GB
  end type limit_type

  type :: dynamics_type
     character(len=16)   :: scheme
     real(8)             :: delta
     real(8)             :: total
  end type dynamics_type

  type :: solver_inner_cg_type
     integer             :: max_iteration
     real(kind(1d0))     :: convergence_eps
  end type solver_inner_cg_type

  type :: flexible_cutoff_type    ! in calc%solver tag
     logical             :: set
     real(kind(1d0))     :: flexible_cutoff_01
     real(kind(1d0))     :: flexible_cutoff_02
  end type flexible_cutoff_type

  type :: eigen_mpi_type
     character(len=128)   :: SEP_solver
     character(len=128)   :: GS_transformation
     integer              :: blocksize
     integer              :: level_lowest
     integer              :: level_highest
  end type eigen_mpi_type

  type :: solver_type
     character(len=16)   :: scheme
     character(len=64)   :: eigen_mpi_scheme
     integer             :: projection
     integer             :: dimension
     integer             :: mode_for_large_memory
     character(len=32)   :: mode_for_suggest_projection
     integer             :: mArnoldi_q
     type(solver_inner_cg_type) :: inner_cg_loop
     type(flexible_cutoff_type) :: flexible_cutoff
     type(eigen_mpi_type) :: eigen_mpi
  end type solver_type

  type :: optimization_type
     character(len=20)   :: scheme      !
     integer             :: max_num_iter! maximum number of iterations
     real(kind(1d0))     :: sd_ratio    ! for "SD" type calculation
     character(len=20)   :: convergence_mode ! "off", "none(=off)",  "energy_per_atom", "force_per_atom"
     real(kind(1d0))     :: energy_convergence
     real(kind(1d0))     :: force_convergence
  end type optimization_type

  type :: snapshot_type
     integer             :: initial  ! initial snapshot (0,1,2...)
     integer             :: number   ! number of snapshots
     integer             :: interval ! interval for loading the snapshot
  end type snapshot_type

  type :: calc_check_item_type    ! used in calc%check
     logical             :: set
     real(kind(1d0))     :: warning_level
     real(kind(1d0))     :: abort_level
  end type calc_check_item_type

  type :: calc_check_type    ! used for calc%check
     logical                    :: set
     type(calc_check_item_type) :: short_atom_pair_distance
  end type calc_check_type

  type :: cell_change_type
     logical             :: set
     character(len=20)   :: scheme      !
     integer             :: max_num_iter! maximum number of iterations
     character(len=64)   :: filename
  end type cell_change_type

  type :: genoOption_type
     logical             :: Use_VOIP
     character(len=8)    :: genoMethod
     character(len=8)    :: CSC_method
     integer             :: CSC_max_loop_count
     real(kind(1d0))     :: CSC_charge_convergence
     real(kind(1d0))     :: CSC_charge_mixing_ratio
     integer             :: CSC_mode_for_tuning_mixing_ratio
     logical             :: HML_constant_K
     real(kind(1d0))     :: HML_small_delta
     real(kind(1d0))     :: HML_kappa
     real(kind(1d0))     :: HML_K
     logical             :: set_S_as_I
     logical             :: vanDerWaals
     real(kind(1d0))     :: vanDerWaals_lambda
  end type genoOption_type

  type :: micro_cell_type
     logical :: set
     character(len=20)   :: mode   ! "one_cell", "default" or "inactive"
     integer :: cell_number(3)
     integer :: search_range(3)
  end type micro_cell_type

  type :: distributed_type
     logical               :: set
     logical               :: global_dens_mat
     logical               :: global_ham_mat
     type(micro_cell_type) :: micro_cell_booking
     logical               :: mpi_is_active
     integer               :: myrank
     integer               :: nprocs
     integer               :: log_unit
     logical               :: root_node
     logical               :: use_mpi_barrier
  end type distributed_type

  type :: interation_range_type
     logical               :: cutoff_rest_cellmax
     real(kind(1d0))       :: cutoff_rest
     real(kind(1d0))       :: cutoff_rest_non_vdW
  end type interation_range_type

  type :: wave_packet_type
     logical               :: set
     character(len=20)     :: mode
     real(kind(1d0))       :: delta_t
     real(kind(1d0))       :: limit_t
     real(kind(1d0))       :: replace_t
     real(kind(1d0))       :: charge_factor_common
     real(kind(1d0))       :: charge_factor_H
     real(kind(1d0))       :: charge_factor_C
     real(kind(1d0))       :: phase_factor_coef
     real(kind(1d0))       :: localize_potential_depth
     real(kind(1d0))       :: temperature
     real(kind(1d0))       :: alpha_delta_min_x
     real(kind(1d0))       :: alpha_delta_max_x
     real(kind(1d0))       :: perturb_interval
     real(kind(1d0))       :: amplitude_print_threshold
     real(kind(1d0))       :: amplitude_print_interval
     integer               :: alpha_delta_index
     integer               :: fst_filter
     integer               :: end_filter
     integer               :: num_steps_per_output_split
     integer               :: num_group_filter_from_homo
     integer               :: localize_start
     integer               :: localize_end
     integer               :: output_interval
     character(len=1024)   :: output_filename
     character(len=1024)   :: h1_type
     character(len=1024)   :: time_evolution_mode
     character(len=1024)   :: init_type
     character(len=1024)   :: filter_mode
     character(len=1024)   :: filter_group_filename
     logical               :: to_multiply_phase_factor
     logical               :: is_output_split
     logical               :: is_overlap_ignored
     logical               :: is_binary_output_mode
     ! Restrictions for setting in ELSES with wavepacket.
     ! - Restarting is not supported.
     ! - Reduction mode is not supported.
     ! - Multistep input and group id are always enabled.
     ! - Reading eigenpairs or initial values are not supported.
  end type wave_packet_type

  type :: calc_type
     character(len=20)   :: mode        ! "dynamics", "optimization" and so on
     type(optimization_type) :: optimization
     type(genoOption_type)   :: genoOption
     type(limit_type)    :: limit
     type(dynamics_type) :: dynamics
     type(solver_type)   :: solver
     type(distributed_type)   :: distributed
     type(snapshot_type)      :: snapshot
     type(cell_change_type)   :: cell_change
     type(calc_check_type)    :: calc_check
     type(wave_packet_type)   :: wave_packet
     character(len=20)        :: calc_force_mode
     logical                  :: use_integer_elec_num
     logical                  :: calc_virial_pressure
     logical                  :: use_group_id
     logical                  :: constraint_w_group
     type(interation_range_type)    :: interaction_range
  end type calc_type

  type :: file_type
     logical             :: set
     character(len=64)   :: dirname
     character(len=64)   :: filename
     character(len=8)    :: history
     character(len=8)    :: format
     character(len=8)    :: method
     integer             :: interval
     character(len=8)    :: append_mode
     logical             :: first_write
     character(len=8)    :: cell_info  ! for xyz-formatted position file
     logical             :: split
     integer             :: number_of_split_files
     logical             :: atom_id_is_added
  end type file_type

  type :: file_matrix_type
     logical             :: set
     character(len=64)   :: filename
     character(len=64)   :: mode  ! "last" "periodic" "default(=last)" "off"
     character(len=64)   :: format
     integer             :: interval
     character(len=8)    :: unit  ! "a.u." or "eV"
     logical             :: first_write
  end type file_matrix_type

  type :: output_type
     type(file_type)     :: main
     type(file_type)     :: restart
     type(file_type)     :: position
     type(file_type)     :: wavefunction
     type(file_type)     :: atom_charge
     type(file_type)     :: atom_energy
     type(file_type)     :: density_of_states
     type(file_type)     :: bond_list
     type(file_type)     :: velocity
     type(file_matrix_type) :: matrix_hamiltonian
     type(file_matrix_type) :: matrix_overlap
     type(file_matrix_type) :: eigen_level
  end type output_type

  type :: option_type
     character(len=256)  :: filename
     integer             :: verbose       ! Old variable
     integer             :: verbose_level ! New variable (8. Feb. 2015)
     integer             :: debug
     character(len=16)   :: functionality
     integer             :: log_node_number
  end type option_type

  type :: config_type
     type(option_type)   :: option
     character(len=64)   :: name
     type(system_type)   :: system
     type(calc_type)     :: calc
     type(output_type)   :: output
     real(8)             :: elses_xml_version
  end type config_type

  type(config_type)      :: config

contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine option_default( option )
    implicit none
    type(option_type)    :: option

    option%filename = "config.xml"
    option%verbose = 0
    option%debug   = 0
    option%functionality = "default"
    option%log_node_number = 8
    return
  end subroutine option_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_default( calc )
    implicit none
    type(calc_type)      :: calc
    call limit_default( calc%limit )
    call dynamics_default( calc%dynamics )
    call solver_default( calc%solver )
    call distributed_default( calc%distributed )
    call snapshot_default( calc%snapshot )
    calc%use_integer_elec_num = .false.
    calc%calc_virial_pressure = .false.
    return
  end subroutine calc_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine limit_default( limit )
    implicit none
    type(limit_type)     :: limit

    limit%time  = -10.0d0 ! dummy value
    limit%memory = -1.0d0 ! dummy value

    limit%memory_allocated = 0.0d0 ! initial value

    return
  end subroutine limit_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine dynamics_default( dynamics )
    implicit none
    type(dynamics_type)  :: dynamics

    dynamics%scheme = "velocity verlet"
    dynamics%delta = 1.00
    dynamics%total = 1000.00

    return
  end subroutine dynamics_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine snapshot_default( snapshot )
    implicit none
    type(snapshot_type)  :: snapshot
    snapshot%initial  = 0
    snapshot%number   = 1  ! dummy value
    snapshot%interval = 1
    return
  end subroutine snapshot_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine solver_default( solver )
    implicit none
    type(solver_type)    :: solver
    solver%scheme = "scheme_default"
    solver%projection = 200
    solver%dimension = 30
    solver%mode_for_large_memory = 1
    solver%mode_for_suggest_projection = "default"
    solver%mArnoldi_q        = 15
    solver%inner_cg_loop%max_iteration = 100
    solver%inner_cg_loop%convergence_eps = -8
    solver%flexible_cutoff%set = .false.
    solver%flexible_cutoff%flexible_cutoff_01 = -1.0d0
    solver%flexible_cutoff%flexible_cutoff_02 = -1.0d0
    return
  end subroutine solver_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine cell_change_default( cell_change )
    implicit none
    type(cell_change_type)    :: cell_change
    cell_change%set          = .false.
    cell_change%max_num_iter = 0
    cell_change%filename     = ""
    return
  end subroutine cell_change_default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine output_default( output )
    implicit none
    type(output_type)    :: output

    output%restart%set = .false.
    output%restart%filename = ""

    output%position%set = .false.
    output%position%filename = ""

    output%density_of_states%set = .false.
    output%density_of_states%filename = ""

    return
  end subroutine output_default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine distributed_default( distributed )
    implicit none
    type(distributed_type)    :: distributed

    distributed%set = .false.
    distributed%global_dens_mat = .false.
    distributed%global_ham_mat = .false.

    distributed%micro_cell_booking%set  = .false.
    distributed%micro_cell_booking%mode = "default"
    distributed%micro_cell_booking%cell_number(:)  = -1
    distributed%micro_cell_booking%search_range(:) = -1

    return
  end subroutine distributed_default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine config_print( config )
    implicit none
    type(config_type)    :: config

    write(*,*) "config%name = ", config%name

    call system_print( config%system )
    call calc_print( config%calc )
    call output_print( config%output )

    return
  end subroutine config_print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine system_print( system )
    implicit none
    type(system_type)    :: system
    integer              :: i

    write(*,*) "system%temperature = ", system%temperature

    call structure_print( system%structure )

    write(*,*) "system%ntarget = ", system%ntarget

    do i=1, system%ntarget
       call target_print( system%vtarget(i) )
    end do

    return
  end subroutine system_print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine target_print( target )
    implicit none
    type(target_type)    :: target

    if( target%class_set ) then
       write(*,*) "target%class = ", target%class
    end if
    if( target%xrange_set ) then
       write(*,*) "target%xrange = ", target%xrange(1:2)
    end if
    if( target%yrange_set ) then
       write(*,*) "target%yrange = ", target%yrange(1:2)
    end if
    if( target%zrange_set ) then
       write(*,*) "target%zrange = ", target%zrange(1:2)
    end if
    if( target%trange_set ) then
       write(*,*) "target%trange = ", target%trange(1:2)
    end if

    if( target%temperature_set ) then
       write(*,*) "target%temperature = ", target%temperature
    end if

    if( target%mass_set ) then
       write(*,*) "target%mass = ", target%mass
    end if

    if( target%position_set ) then
       write(*,*) "target%position = ", target%position(1:3)
    end if

    if( target%velocity_set ) then
       write(*,*) "target%velocity = ", target%velocity(1:3)
    end if

    if( target%motion_set ) then
       write(*,*) "target%motion = ", target%motion
    end if

    return
  end subroutine target_print

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine calc_print( calc )
    implicit none
    type(calc_type)      :: calc

    call limit_print( calc%limit )
    call dynamics_print( calc%dynamics )
    call solver_print( calc%solver )

    return
  end subroutine calc_print

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine limit_print( limit )
    implicit none
    type(limit_type)     :: limit

    write(*,*) "limit%time   = ", limit%time
    write(*,*) "limit%memory = ", limit%memory

    return
  end subroutine limit_print

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine dynamics_print( dynamics )
    implicit none
    type(dynamics_type)  :: dynamics

    write(*,*) "dynamics%scheme = ", trim(dynamics%scheme)
    write(*,*) "dynamics%delta = ", dynamics%delta
    write(*,*) "dynamics%total = ", dynamics%total

    return
  end subroutine dynamics_print

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine solver_print( solver )
    implicit none
    type(solver_type)    :: solver

    write(*,*) "solver%scheme = ", trim(solver%scheme)
    write(*,*) "solver%projection = ", solver%projection
    write(*,*) "solver%dimension  = ", solver%dimension

    return
  end subroutine solver_print

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine output_print( output )
    implicit none
    type(output_type)    :: output

    call file_print( output%restart )
    call file_print( output%position )
    call file_print( output%density_of_states )

    return
  end subroutine output_print

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine file_print( file )
    implicit none
    type(file_type)      :: file

    if( .not. file%set ) return

    write(*,*) "file%filename = ", trim(file%filename)
    write(*,*) "file%history  = ", file%history
    write(*,*) "file%format   = ", file%format
    write(*,*) "file%interval = ", file%interval

    return
  end subroutine file_print

  function target_match( target, atom ) result(match)
    implicit none
    type(target_type), intent(in) :: target
    type(atom_type),   intent(in) :: atom

    logical              :: match

    match = .false.

    if( target%class_set ) then
       if( target%class == atom%class ) then
          match = .true.
       else
          match = .false.
          return
       end if
    end if

    if( target%xrange_set ) then
       if( target%xrange(1) <= atom%position(1) .and. &
            atom%position(1) <= target%xrange(2) ) then
          match = .true.
       else
          match = .false.
          return
       end if
    end if

    if( target%yrange_set ) then
       if( target%yrange(1) <= atom%position(2) .and. &
            atom%position(2) <= target%yrange(2) ) then
          match = .true.
       else
          match = .false.
          return
       end if
    end if

    if( target%zrange_set ) then
       if( target%zrange(1) <= atom%position(3) .and. &
            atom%position(3) <= target%zrange(2) ) then
          match = .true.
       else
          match = .false.
          return
       end if
    end if

    return
  end function target_match

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine target_overwrite( target, atom )
    implicit none
    type(target_type), intent(in)  :: target
    type(atom_type),   intent(out) :: atom

    if( target%temperature_set ) then
!       atom%temperature = target%temperature
    end if

    if( target%mass_set ) then
!       atom%mass = target%mass
    end if

    if( target%position_set ) then
       atom%position = target%position
    end if

    if( target%velocity_set ) then
       atom%velocity = target%velocity
    end if

    if( target%motion_set ) then
       atom%motion = target%motion
    end if

    return
  end subroutine target_overwrite
end module M_config
