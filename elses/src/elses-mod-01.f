!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-mod-01.f @@@@@
czzz  @@@@@ 2008/08/22 @@@@@
cccc2006cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0722: Prepared. (NT116p31): Modules for LSES engine
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0104: (TXP0, TYP0, TZP0) are added (NT07Ap01)
c        --> allocated in elses_alloc_txp0(noa_def)
c    : Module for the Nose thermostat is added (NT07A-p01).
c        --> module 'elses_mod_heat_thermo'
c    : Module for force on atom at the previous timestep
c        --> module 'module elses_mod_foiold'
c    : Module for flag for atoms:  moving(=1) or fixed(=0)
c        --> module 'elses_mod_iflag'
c0324: 'target' is added for the compatibility of NRL-part code.
c          (NT07A-118,p39) ---> jsei
c    : 'awt(:)' is prepared in elses_mod_mass (NT07A-118,p.39)
c0517: Module val_elec_atm is prepared. (NT07B-119,p15)
c        ---> valence electron numbmer per atom, and
c             total number of valence electrons
c0606: Module 'elses_mod_ctrl' is prepared.
!0722: Module 'elses_mod_elem_name' is prepared. (NT07B-119,p60)
!         ----> element name, such as 'Si', 'Ge' 
!    : Initial dummy value is set in elses_alloc_vel
!0805: In 'elses_mod_ctrl', iverbose is prepared. (NT07B-119,p69)
!0823: A modi. in 'elses_alloc_val_elec'
!ccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0108T.Hoshi: FOI, FOIOLD is initizalized (=0.0d0) 
!               after allocation (NT07E-122p39, v000a-wrk05)
!0401T.Hoshi: i_pbc_x, i_pbc_y, i_pbc_z are prepared.
!              in module elses_mod_sim_cell
!               (NT08A-123p23, 0.00.14a-wrk05) 
!0822T.Hoshi: para_spin_factor=2.0d0 is introduced as physical constant
!               (NT08E-127)
!ccc2009cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!01101T.Hoshi: The module for physical constants 
!               are moved into elses-lib-phys-const.f90
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_phys_const
c
        use M_lib_phys_const
c
c       real*8, parameter :: pi=.314159265358979323D1
c
c       real*8, parameter :: ev4au=2.0d0*13.6058d0
c                        --->   1 au = (ev4au) eV
c
c       real*8, parameter :: angst=0.529177d0
c                        --->   1 au = 0.529177 A
c
c       real*8, parameter :: au_fsec=124.06948d0/3.0d0
c                        --->   1 fsec = (AUFSEC) au
c
c       real*8, parameter :: au_mass=1.6605655d-24/9.109534d-28
c
c       real*8, parameter :: ev2kel=1.60217733d0/1.380658d0*10000.0d0
c                        --->   1 eV  = (EV2KEL) kelvin
c
c       real*8, parameter :: ene_j4ev=1.60219d-19
c                        ---> 1 [eV] = 1.60219 x 10^{-19} [J]
c
c       real(8), parameter :: para_spin_factor=2.0d0
c
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_sim_cell
        integer noa
        integer nos
!       integer iperiodic
        integer i_pbc_x, i_pbc_y, i_pbc_z
        real*8  ax,ay,az
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_tx
        real*8,          allocatable :: tx(:),ty(:),tz(:)
        integer, target, allocatable :: jsei(:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_txp
        real*8, allocatable :: txp(:),typ(:),tzp(:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Force on atoms
c
      module elses_mod_foi
        real*8, allocatable :: foi(:,:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Controlling parameters
c
      module elses_mod_ctrl
        integer i_verbose
c         ---> verbose mode for when i_verbose = 0
        real*8  r_mem_limit
c         ----> memory limit in GB
        real*8  r_time_limit
c         ----> time limit in sec
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Force on atoms at the previous time step
c
      module elses_mod_foiold
        real*8, allocatable :: foiold(:,:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Flag for atoms:  moving(=1) or fixed(=0)
c
      module elses_mod_iflag
        integer, allocatable :: iflag(:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Element name
c
      module elses_mod_elem_name
        character(len=8), allocatable :: elem_name(:)
!         --> defined for each element: elem_name(1:nos)
!                Ex. elem_name(1)="Si"
!
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Paramater for mass 
c
      module elses_mod_mass
        real*8, allocatable :: awt(:)
c          --> atomic mass
c               defined for each element: awt(1:nos)
c
        real*8, allocatable :: amm(:)
c          --> (mass)^{-1} in a.u.
c               defined for each atom amm(1:noa)
c
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Paramater for charge
c
      module elses_mod_val_elec
        real(8), allocatable :: val_elec_atm(:)
c          --> valence electron number per atom
c               defined for each element : val_elec_atm(1:nos)
        real(8) :: val_elec_tot
c          --> total number of valence electrons
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_vel
        real*8, allocatable :: velx(:),vely(:),velz(:)
c          --> (velocity)*dt, NOT velocity itself
c                where length is normalized by ax, ay or az.
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Paramater for Heat-bath 
c    (All quantities are written in au.)
c
      module elses_mod_thermo
c       real*8 tempk0,amq,thb,thbold,vhb,vhbold
        real*8 amq,thb,thbold,vhb,vhbold
c             tempk0      : temperature (input parameter)
c             amq         : (mass)^(-1)
c             thb, thbold : cordinate
c             vhb, vhbold : velocity * dt
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Quaitities in MD simulation (as dynamics)
c
      module elses_mod_md_dat
        integer itemd
c          --> iteration number in MD
c
        integer itemdorg
c          --> total iteration number in MD
c              from the original structure
c              (if the calculation is a continued simulation)
c
        integer itemdmx
c          --> maximum iteration number in MD
c
        real*8  dtmd
c          --> time step in MD
c
        real*8  e_kin
c          --> Kinetic energy of atoms in MD
c
        logical first_iteration
c          --> True at the first MD iteration
c
        logical final_iteration
c          --> True at the final MD iteration
c
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of atom mass : AMM(1:NOA_DEF)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_amm(noa_def)
      use elses_mod_mass, only : amm
      implicit none
      integer noa_def,ierr
c
      write(6,*)'@@ ALLOC_AMM:NOA_DEF=',noa_def
c
      allocate (amm(noa_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(AMM):ierr=',ierr
        stop
      endif
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of atom mass : AWT(1:NOS_DEF)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_awt(nos_def)
      use elses_mod_mass, only : awt
      implicit none
      integer nos_def,ierr
c
      write(6,*)'@@ ALLOC_AWT:NOS_DEF=',nos_def
c
      allocate (awt(nos_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(AWT):ierr=',ierr
        stop
      endif
      awt(:)=-1.0d0
c       ---> initial dummy value
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of val_elec_atm(1:NOS_DEF)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_val_elec(nos_def)
      use elses_mod_val_elec, only : val_elec_atm, val_elec_tot
      implicit none
      integer nos_def,ierr
c
      write(6,*)'@@ ALLOC_VAL_ELEC:NOS_DEF=',nos_def
c
      if (allocated(val_elec_atm) .eqv. .true.) then
         write(*,*)'ERROR:.. already allocated'
         stop
      endif   
c
      allocate (val_elec_atm(nos_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(val_elec_atm):ierr=',ierr
        stop
      endif
      val_elec_atm(:)=-1.0d0
      val_elec_tot=-1.0d0
c     ----> dummy values
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of elem_name(1:NOS_DEF)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_elem_name(nos_def)
      use elses_mod_elem_name, only : elem_name
      implicit none
      integer nos_def,ierr
c
      write(6,*)'@@ ALLOC_elem_name:NOS_DEF=',nos_def
c
      allocate (elem_name(nos_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'alloc. error!(val_elec_atm):ierr=',ierr
        stop
      endif
      elem_name(:)="ZZ"
c     ----> dummy values
c
      end
