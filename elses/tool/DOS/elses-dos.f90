!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
! mostly imported from elses-qm-solver-eig-geno.f90

program elses_dos

! use M_qm_domain,          only : i_verbose, noav
  use M_io_ctrl_dos_plot,   only : initial_preparation_for_dos !(routine)
  use M_lib_phys_const,     only : pi !(constant)
! use M_config,             only : config !(unchanged )
  use elses_mod_phys_const, only : ev4au  !(unchanged )
! use elses_arr_eig_leg, only : atmp, atmp2, eig2, f_occ  !(unchanged )
! use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged )
! use elses_mod_js4jsv,   only : js4jsv                  !(unchanged )
  use elses_mod_file_io, only : vacant_unit
  use M_lib_math_func,      only : Fermi_Dirac_Func
  
  implicit none
  integer, parameter :: DOUBLE_PRECISION=kind(1d0)

  integer :: i,j,k, step_count, i_enable, ierr, jj, imode
  integer :: number_energy_mesh
  integer :: atm_index, orb_index
  integer :: iunit
  real(DOUBLE_PRECISION) :: sum_energy, sum_occupation
  real(DOUBLE_PRECISION) :: sum_eigen_levels
  !
  real(DOUBLE_PRECISION) :: width_of_peak
  real(DOUBLE_PRECISION) :: x
  real(DOUBLE_PRECISION) :: step_fn, weight
  real(DOUBLE_PRECISION) :: energy_btm, energy_top, de
  real(DOUBLE_PRECISION), allocatable  :: energy_mesh(:)
  real(DOUBLE_PRECISION), allocatable  :: nos(:)
  real(DOUBLE_PRECISION), allocatable  :: enos(:)
  real(DOUBLE_PRECISION), allocatable  :: pnos(:,:)
! real(DOUBLE_PRECISION), allocatable  :: lnos(:,:)
  real(DOUBLE_PRECISION), allocatable  :: eigen_states(:,:)
  real(DOUBLE_PRECISION) :: nos_prv, nos_tmp, dos_tmp, ene_tmp, ene_tot
  !
  character(len=250) :: string, string2
  !
  real(DOUBLE_PRECISION), allocatable :: atmp(:,:), atmp2(:,:), eig2(:), f_occ(:)
  integer :: i_verbose
  integer :: n_tot_base
!
  i_verbose = 1
!
  write(*,'("Input # of levels")')
  read (*,*)  n_tot_base

  call initial_preparation_for_dos(imode, number_energy_mesh, energy_btm, energy_top, width_of_peak)
  
  if (i_verbose >= 1) then
     write(*,*)'@@ plot_ldos_from_eigen (experimental routine)'
     write(*,*)' i_enable=',i_enable
     write(*,*)' imode  =',imode
     write(*,*)' number of energy mesh  =',number_energy_mesh
     write(*,*)' energy bottom   [eV]   =',energy_btm*ev4au
     write(*,*)' energy top      [eV]   =',energy_top*ev4au
     write(*,*)' width of peak      [eV] =',width_of_peak*ev4au
  endif
  !
  allocate (energy_mesh(0:number_energy_mesh),stat=ierr)
  if (ierr /= 0) stop 'Abort:ERROR in alloc'
  !    
  allocate (nos(0:number_energy_mesh),stat=ierr)
  if (ierr /= 0) stop 'Abort:ERROR in alloc'
  nos(:)=0.0d0
  !
  allocate (enos(0:number_energy_mesh),stat=ierr)
  if (ierr /= 0) stop 'Abort:ERROR in alloc'
  enos(:)=0.0d0
  !
  allocate (pnos(0:number_energy_mesh,n_tot_base),stat=ierr)
  if (ierr /= 0) stop 'Abort:ERROR in alloc:pnos'
  pnos(:,:)=0.0d0
  !
! allocate (lnos(0:number_energy_mesh,noav),stat=ierr)
! if (ierr /= 0) stop 'Abort:ERROR in alloc:lnos'
! lnos(:,:)=0.0d0
  !
  de=(energy_top-energy_btm)/dble(number_energy_mesh)
  !      -----> de : energy interval [au]
  !
  do jj=0,number_energy_mesh
     energy_mesh(jj)=energy_btm+de*dble(jj)
     !      -----> energy mesh in eV
  enddo

  allocate (f_occ(n_tot_base))
  allocate (eig2 (n_tot_base))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! load eigen levels
  iunit=vacant_unit()
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! @@ Eigen levels
  !
  open(iunit, file='levels.txt', status='old',PAD='YES')
  !
  write(*,*)'read eigen levels from levels.txt'
  !
  do while (.not. isFoundInStream(iunit,"Orbital energies and kinetic energies"))
  end do

  read(iunit,*) string ! skip single data line
  do k=1,n_tot_base
     read(iunit,*) string, string2, eig2(k) ! skip , eig2(k)*ev4au
     jj=len_trim(string2)
     select case(string2(jj:jj))
     case("O")
        f_occ(k)=1d0
     case("V")
        f_occ(k)=0d0
     case default
        stop 'Illeagal format (OCCUPATION)'
     end select
  enddo
  !    
  close(iunit)
  !
  !
  write(*,*)' number of total bases      =', n_tot_base
  write(*,*)' Eigen level (lowest)  [eV] = ',eig2(1)*ev4au
  write(*,*)' Eigen level (highest) [eV] = ',eig2(n_tot_base)*ev4au
  !
  if (i_verbose >= 40) then
     write(*,*)'Eigen levels (eV) and occupation number'
  endif
  !
  sum_energy=0.0d0
  sum_occupation=0.0d0
  sum_eigen_levels=0.0d0
  do k=1,n_tot_base
     if (i_verbose >= 40) then
        write(*,'(a,2I10,2F15.5)')'eigen_level:',step_count,k, eig2(k)*ev4au, f_occ(k)
     endif
     sum_energy=sum_energy+2.0d0*eig2(k)*f_occ(k)
     sum_eigen_levels=sum_eigen_levels+2.0d0*eig2(k)
     sum_occupation=sum_occupation+2.0d0*f_occ(k)
  enddo
  write(*,*)'ETB (as Sum_i E_i f_i) [au]=',sum_energy
  write(*,*)'ETB (as Sum_i E_i f_i) [eV]=',sum_energy*ev4au
  write(*,*)'Sum of eigen levels (as 2 Sum_i E_i ) [au]=',sum_eigen_levels
  write(*,*)'Sum of eigen levels (as 2 Sum_i E_i ) [eV]=',sum_eigen_levels*ev4au
  write(*,*)'Total elec. number (check): N_elec=',sum_occupation
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! @@ Total DOS
  !
  if (imode >= 1) then
     !
     open(iunit, file='output_tdos.txt', status='unknown')
     !
     write(*,*)'Save TDOS into output_tdos.txt'
     nos(:)=0.0d0
     enos(:)=0.0d0
     !$omp  parallel default(shared) &
     !$omp& private(k, jj) &
     !$omp& private(x, step_fn)
     !$omp  do schedule(static)
     do jj=0,number_energy_mesh
        do k=1,n_tot_base
           x=(eig2(k)-energy_mesh(jj))/width_of_peak
           step_fn=Fermi_Dirac_Func(x)
           nos(jj)  =  nos(jj) + step_fn*2.0d0
           enos(jj) = enos(jj) + step_fn*2.0d0*eig2(k) ! enos in au
        enddo
     enddo
     !$omp end do
     !$omp end parallel
     
     !
     nos_prv=nos(0)
     !   ene_tot=0.0d0
     do jj=0,number_energy_mesh
        nos_tmp=nos(jj)
        dos_tmp=(nos_tmp-nos_prv)/(de*ev4au)
        !     ene_tmp=energy_mesh(jj)
        !     ene_tot=ene_tot+ene_tmp*dos_tmp*de
        write(iunit,'(a,i7,4f20.10)')'tdos=', jj, energy_mesh(jj)*ev4au, dos_tmp, nos(jj), enos(jj)*ev4au
        nos_prv=nos_tmp
     enddo
     !    
     close(iunit)
     !
  endif
  !   
  
  contains

    logical function isFoundInStream(iunit,searchString)
      implicit none
      integer          :: iunit
      character(len=*) :: searchString
      character(len=132)   :: string
      intent(in)       :: iunit, searchString
      read(iunit,'(A132)')string
      isFoundInStream = isFoundInString(trim(searchString),string)
      if(isFoundInStream)write(*,*)string
    end function isFoundInStream

    logical function isFoundInString(searchString,string)
      implicit none
      character(len=*),intent(in):: searchString,string
      if( index(string,searchString) ==0 )then
         isFoundInString = .false.
      else
         isFoundInString = .true.
      end if
    end function isFoundInString
  
end program elses_dos

