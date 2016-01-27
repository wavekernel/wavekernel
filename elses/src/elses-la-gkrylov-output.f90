!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_gkrylov_output
!
  implicit none
!
!
  private
  public :: calc_eigen_in_gkrylov
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Plot eigen state in the gKrylov method
!
!
  subroutine calc_eigen_in_gkrylov
!
    use M_la_krgl_main                  !(unchagend)
    use M_lib_phys_const,     only : pi !(constant)
    use M_config,     only : config !(unchanged )
    use elses_mod_phys_const, only : ev4au  !(unchanged )
!
!   use elses_arr_eig_leg, only : atmp, atmp2, eig2, f_occ  !(unchanged )
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged )
!
    use elses_mod_file_io, only : vacant_unit
    use M_qm_domain,      only : i_verbose, noav, nval, temp_for_electron, &
&                                atm_element, chemical_potential, total_electron_number !(unchanged)
!
    use M_lib_math_func,      only : Fermi_Dirac_Func !(function)
    use M_io_ctrl_dos_plot,   only : initial_preparation_for_dos !(routine)
    use elses_arr_eig_leg,    only : atmp !(CHANGED)
    use M_output_eigenstates, only : output_eigenstates  !(routine)
    use M_io_ctrl_output_eigen, only : init_for_plot_wavefunction !(routine)
    use M_la_gkrylov_hist, only : u_b_hst_str,  v_mat_kr_str !(unchanged)
!
    implicit none
    integer :: i,j,k, step_count, ierr, jj
    integer :: number_energy_mesh
    integer :: atm_index, orb_index
    integer :: iunit1, iunit2
    real(8) :: sum_energy, sum_occupation
    real(8) :: width_of_peak
    real(8) :: x
    real(8) :: step_fn, weight
    real(8) :: energy_btm, energy_top, de
    real(8), allocatable  :: energy_mesh(:)
    real(8), allocatable  :: nos(:)
!
    integer :: imode
    integer :: kr_dim_max, kr_dim_max_input
    integer :: jsv, nss, nval2, ja, nrecl, irec
    real(8) :: rRd, rEd
    real(8) :: nos_prv, nos_tmp, dos_tmp
    real(8) :: enos_prv, enos_tmp, edos_tmp
!
    real(8) :: energy_inf, nos_inf, enos_inf
!
    real(8), allocatable  :: eigen_vector(:,:)
    real(8), allocatable  :: eigen_level(:)
    real(8), allocatable  :: eigen_vector0(:,:)
    real(8), allocatable  :: eigen_vector1(:,:)
    integer(8), allocatable  :: mesh_index_for_level(:,:)
    integer :: min_eigen_level,  max_eigen_level, num_eigen_level
    integer :: nos_int
    real(8) :: ddd, ddmin
    character(len=100) :: filename_wfn
    integer :: imode2
    integer :: kr_dim
    real(8), allocatable  :: v_b_hst_str(:,:,:)
    integer :: nval_max
    integer :: al
    integer :: kr_dim_max_tmp
    real(8) :: ene1, ene2, ddd1, ddd2, ddmax
    real(8), allocatable :: wt_kr_k(:,:,:,:)
    integer, allocatable :: indices_max_weight(:,:)
    integer, allocatable :: sign_factor(:,:,:,:)  ! (+1) or (-1) for wavefuntion
    integer :: al_m, orb_m, atm_m
    real(8) :: result_value
    real(8) :: ratio_criteria
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    ratio_criteria=0.00000001d0
!
    min_eigen_level=max(1, int(total_electron_number/2.0d0)-15)
    max_eigen_level=min(n_tot_base, int(total_electron_number/2.0d0)+16)
!
!   min_eigen_level=6
!   max_eigen_level=7
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    num_eigen_level=max_eigen_level-min_eigen_level+1
    kr_dim_max=config%calc%solver%dimension
    kr_dim_max_input=config%calc%solver%dimension
    nval_max=maxval(nval)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      write(*,*)'@@ calc_eigen_in_gkrylov'
      write(*,*)' total electron number =',total_electron_number
      write(*,*)' total basis    number =',n_tot_base
      write(*,*)' ratio_criteria        =',ratio_criteria
      write(*,*)' min_eigen_level       =',min_eigen_level
      write(*,*)' max_eigen_level       =',max_eigen_level
      write(*,*)' num_eigen_level       =',num_eigen_level
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call initial_preparation_for_dos(imode, number_energy_mesh, energy_btm, energy_top, width_of_peak)
    if (imode == 0) then
      if (i_verbose >= 1) write(*,*)'.... is skipped'
      return
    endif   
!
    if (i_verbose >= 1) then
       write(*,*)' number of energy mesh  =',number_energy_mesh
       write(*,*)' energy bottom   [eV]   =',energy_btm*ev4au
       write(*,*)' energy top      [eV]   =',energy_top*ev4au
       write(*,*)' width of peak      [eV] =',width_of_peak*ev4au
       write(*,*)' chemical potential [au] =',chemical_potential
       write(*,*)' chemical potential [eV] =',chemical_potential*ev4au
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate (energy_mesh(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!    
    allocate (nos(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    nos(:)=0.0d0
!
    allocate (eigen_vector0(n_tot_base,num_eigen_level),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    eigen_vector0(:,:)=0.0d0
!
    allocate (eigen_vector1(n_tot_base,num_eigen_level),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    eigen_vector1(:,:)=0.0d0
!
    allocate (eigen_vector(n_tot_base,num_eigen_level),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    eigen_vector(:,:)=0.0d0
!
    allocate (eigen_level(n_tot_base),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    nos(:)=0.0d0
!
    allocate (mesh_index_for_level(n_tot_base,2),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    mesh_index_for_level(:,:)=0
!
    de=(energy_top-energy_btm)/dble(number_energy_mesh)
!      -----> de : energy interval [au]
!
    do jj=0,number_energy_mesh
      energy_mesh(jj)=energy_btm+de*dble(jj)
!      -----> energy mesh in au
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ NOS calculation
!
!
    energy_inf=10.0d60
!
    nos_inf=0.0d0
    enos_inf=0.0d0
    do jsv=1,noav
      nss=atm_element(jsv)
      nval2=nval(nss)
      do ja=1,nval2
        nrecl=kr_dim_max
        do irec=nrecl,1,-1
          rRd=wt_kr(irec,ja,jsv)*2.0d0
          rEd=eig_kr(irec,ja,jsv)
          do jj=0,number_energy_mesh
            x=(rEd-energy_mesh(jj))/width_of_peak
            step_fn=Fermi_Dirac_Func(x)
            nos(jj)  =  nos(jj)+step_fn*rRd
          enddo
          x=(rEd-energy_inf)/width_of_peak
          step_fn=Fermi_Dirac_Func(x)
          nos_inf  = nos_inf  + step_fn*rRd
        enddo
      enddo
    enddo
!
    write(*,*) 'nos_inf =',nos_inf, energy_inf
!   write(*,*) 'enos_inf=',enos_inf, energy_inf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Determine the (quasi) eigen levels ; e_k
!      mesh_index_for_level(k,1) : mesh index for the integer
!                 N( e_k ) / 2 = k
!      mesh_index_for_level(k,2) : mesh index for the half integer
!                 N( e_k ) / 2 = k- (1/2)
!
!          the factor two in l.h.s is the para-spin factor
!
!          
!
    do k=1,n_tot_base
      nos_int=0
      ddmin=1.0d10 ! dummy value
      do jj=0, number_energy_mesh-1
        ddd=dabs(nos(jj)/2.0d0-dble(k))
        if (ddd < ddmin) then
          ddmin=ddd
          nos_int=jj
        endif   
      enddo
      mesh_index_for_level(k,1)=nos_int
    enddo
!   
    do k=1,n_tot_base
      nos_int=0
      ddmin=1.0d10 ! dummy value
      do jj=0, number_energy_mesh-1
        ddd=dabs(nos(jj)/2.0d0-dble(k)+0.5d0)
        if (ddd < ddmin) then
          ddmin=ddd
          nos_int=jj
        endif   
      enddo
      mesh_index_for_level(k,2)=nos_int
    enddo
!
!   do k=1, n_tot_base
!     write(*,*) 'mesh for level=',k,mesh_index_for_level(k,1:2)
!   enddo
!
    iunit1=vacant_unit()
    open(iunit1, file='output_nos.txt', status='unknown')
    do jj=0, number_energy_mesh-1
      write(iunit1,'(i10,2f20.10)') jj, energy_mesh(jj)*ev4au, nos(jj)/2.0d0
    enddo
    close(iunit1)
!
    iunit1=vacant_unit()
    open(iunit1, file='output_eigen_levels.txt', status='unknown')
    do k=1, n_tot_base
       jj=mesh_index_for_level(k,2)
       write(iunit1,'(i10,f20.10)') k, energy_mesh(jj)*ev4au
    enddo   
    close(iunit1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Allocation 
!
    if (allocated(u_b_hst_str)  .eqv. .false.)  then
       write(*,*) 'ERROR:not allcoated:u_b_hst_str'
       stop
    endif
!
    if (allocated(v_b_hst_str)  .eqv. .true.)  then
       write(*,*) 'ERROR:already allcoated:v_b_hst_str'
       stop
    endif
    allocate(v_b_hst_str(kr_dim_max_input, nval_max, noav), stat=ierr) 
    if (ierr /= 0) then
      write(*,*) 'Alloc. Error'
      stop
    endif  
!
!
    if (allocated(kr_dim_str)  .eqv. .false.)  then
       write(*,*) 'ERROR:not allcoated:kr_dim_str'
       stop
    endif
!
    if (allocated(atmp)  .eqv. .false.)  then
       allocate(atmp(n_tot_base, n_tot_base), stat=ierr) 
       if (ierr /= 0) then
         write(*,*) 'Alloc. Error'
         stop
       endif  
    endif
    atmp(:,:)=0.0d0
!
    allocate(wt_kr_k(kr_dim_max_input, nval_max, noav, num_eigen_level), stat=ierr) 
    if (ierr /= 0) then
      write(*,*) 'Alloc. Error'
      stop
    endif  
!
    allocate(sign_factor(kr_dim_max_input, nval_max, noav, num_eigen_level), stat=ierr) 
    if (ierr /= 0) then
      write(*,*) 'Alloc. Error'
      stop
    endif  
    sign_factor(:,:,:,:)=0
!
    allocate(indices_max_weight(3, num_eigen_level), stat=ierr) 
    if (ierr /= 0) then
      write(*,*) 'Alloc. Error'
      stop
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot weight (optional)
!
!   do atm_index=1,noav
!     do orb_index=1,nval(atm_element(atm_index))
!       do al=1,kr_dim_str(orb_index,atm_index)
!         write(*,'(a,3i10,f30.20)')' wt_kr=', al, orb_index, atm_index, wt_kr(al, orb_index, atm_index)
!       enddo
!       write(*,*)'sum of wt_k, ja, jsr= ', orb_index, atm_index, sum(wt_kr(:,orb_index,atm_index))
!     enddo
!   enddo

!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Weight calculation for k-th eigen state
!
    wt_kr_k(:,:,:,:)=-1.0d0
    do k=min_eigen_level, max_eigen_level
      ene1=energy_mesh(mesh_index_for_level(k,1))
      if (k == 1) then 
        ene2=-1.0d10 ! dummy value ( -infinity )
      else
        ene2=energy_mesh(mesh_index_for_level(k-1,1))
      endif
      do atm_index=1,noav
        do orb_index=1,nval(atm_element(atm_index))
          do al=1,kr_dim_str(orb_index,atm_index)
            ddd1=Fermi_Dirac_Func((-ene1+eig_kr(al,orb_index,atm_index))/width_of_peak)  ! smoothed step func.
            ddd2=Fermi_Dirac_Func((-ene2+eig_kr(al,orb_index,atm_index))/width_of_peak)  ! smoothed step func.
            ddd=wt_kr(al, orb_index, atm_index)*(ddd1-ddd2)
            wt_kr_k(al, orb_index, atm_index, k-min_eigen_level+1)=ddd
            if (dabs(ddd) > 1.0d-10 ) then
              write(*,'(a,4i10,3f30.20)')'k,al,ja,jsv, wt_kr_k1 =', &
&                k, al, orb_index, atm_index, ddd, ddd1, ddd2
              write(*,'(a,4i10,3f30.20)')'k,al,ja,jsv, wt_kr_k2 =', &
&                k, al, orb_index, atm_index, wt_kr_k(al, orb_index, atm_index, k-min_eigen_level+1)
            endif   
          enddo
        enddo
      enddo
    enddo  
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set indices_max_weight: the indices for the maximum weight
!
    indices_max_weight(:,:)=0
    do k=min_eigen_level, max_eigen_level
      ddmax=-1.0d0
      do atm_index=1,noav
        do orb_index=1,nval(atm_element(atm_index))
          do al=1,kr_dim_str(orb_index,atm_index)
            ddd=wt_kr_k(al, orb_index, atm_index, k-min_eigen_level+1)
!           write(*,'(a,4i10,f20.10)') ' k, al, ja, jsv, ddd=',k,al, orb_index, atm_index, ddd
            if (ddd >= ddmax) then
              ddmax=ddd
              indices_max_weight(1,k-min_eigen_level+1)=al
              indices_max_weight(2,k-min_eigen_level+1)=orb_index
              indices_max_weight(3,k-min_eigen_level+1)=atm_index
            endif   
          enddo
        enddo
      enddo
      write(*,'(a,4i10,f20.10)') ' k, indices of (al,ja,jsv) = ', k, indices_max_weight(1:3,k-min_eigen_level+1), ddmax
    enddo  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the sign for the most important component
!
    sign_factor(:,:,:,:)=0
    do k=min_eigen_level, max_eigen_level
      al        = indices_max_weight(1,k-min_eigen_level+1)
      orb_index = indices_max_weight(2,k-min_eigen_level+1)
      atm_index = indices_max_weight(3,k-min_eigen_level+1)
      sign_factor(al, orb_index, atm_index, k-min_eigen_level+1)=1
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the sign for the important components
!
    do k=min_eigen_level, max_eigen_level
      write(*,*)'k = ', k
      al_m  = indices_max_weight(1,k-min_eigen_level+1)
      orb_m = indices_max_weight(2,k-min_eigen_level+1)
      atm_m = indices_max_weight(3,k-min_eigen_level+1)
      ddmax=wt_kr_k(al_m, orb_m, atm_m, k-min_eigen_level+1)
      do atm_index=1,noav
        do orb_index=1,nval(atm_element(atm_index))
          do al=1,kr_dim_str(orb_index,atm_index)
            if ( sign_factor(al, orb_index, atm_index, k-min_eigen_level+1) /= 0 ) cycle
            if ( wt_kr_k(al, orb_index, atm_index, k-min_eigen_level+1) < ddmax*ratio_criteria ) cycle
            call get_v_v_product(al, orb_index, atm_index, al_m, orb_m, atm_m, result_value)
            write(*,'(a,3i10,f30.20)')'al, orb, atm, v_v_product = ', al, orb_index, atm_index, result_value
            if (result_value > 0) then
              sign_factor(al, orb_index, atm_index, k-min_eigen_level+1)=1
            else
              sign_factor(al, orb_index, atm_index, k-min_eigen_level+1)=-1
            endif   
          enddo
        enddo
      enddo
    enddo
!
!   stop 'Stop manually T.H.'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Determine the (quasi) eigen states
!
    do atm_index=1,noav
      nss=atm_element(atm_index)
      nval2=nval(nss)
      do orb_index=1,nval2
        kr_dim_max_tmp=kr_dim_str(orb_index,atm_index)
        do al=1,kr_dim_max_tmp
         v_b_hst_str(al, orb_index, atm_index)= &
&             dot_product( v_mat_kr_str(1:kr_dim_max_tmp, al, orb_index,atm_index),  &
&                          u_b_hst_str(1:kr_dim_max_tmp,     orb_index, atm_index) )
        enddo
      enddo
    enddo
!  
    do atm_index=1,noav
      nss=atm_element(atm_index)
      nval2=nval(nss)
      do orb_index=1,nval2
        j=js2j(orb_index,atm_index)
        do k=min_eigen_level, max_eigen_level
           ene1=energy_mesh(mesh_index_for_level(k,1))
           if (k == 1) then 
             ene2=-1.0d10 ! dummy value ( -infinity )
           else
             ene2=energy_mesh(mesh_index_for_level(k-1,1))
           endif
           ddd=0.0d0
           do al=1,kr_dim_str(orb_index,atm_index)
             ddd1=Fermi_Dirac_Func((-ene1+eig_kr(al,orb_index,atm_index))/width_of_peak)
             ddd2=Fermi_Dirac_Func((-ene2+eig_kr(al,orb_index,atm_index))/width_of_peak)
             ddd=ddd+v_b_hst_str(al, orb_index, atm_index)*(ddd1-ddd2)* &
&                  sign_factor(al, orb_index, atm_index, k-min_eigen_level+1)
           enddo  
           atmp(j,k)=ddd
        enddo
      enddo
    enddo  
!
    call init_for_plot_wavefunction(imode2, filename_wfn)
    if (imode2 == 0) then
      if (i_verbose >= 1) write(*,*)'.... plotting wfn. is skipped (imode2=0)'
      return
    endif   
    if (imode2 == 1) call output_eigenstates(filename_wfn)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Deallocation
!
    deallocate (nos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
  end subroutine calc_eigen_in_gkrylov
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Inner product, (v,v)
!
  subroutine get_v_v_product(al_0, orb_0, atm_0, al_m, orb_m, atm_m, result_value)
!
    use elses_mod_orb2,    only : n_tot_base        !(unchanged)
    implicit none
    integer,  intent(in)  :: al_0, orb_0, atm_0, al_m, orb_m, atm_m 
    real(8),  intent(out) :: result_value
    real(8), allocatable  :: v1(:), v2(:)
    integer               :: ierr
!
!
    allocate (v1(n_tot_base),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    v1(:)=0.0d0
!
    allocate (v2(n_tot_base),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    v2(:)=0.0d0
!
    call get_kr_eig_vec(v1, al_m, orb_m, atm_m)
    call get_kr_eig_vec(v2, al_0, orb_0, atm_0)
!
    result_value=dot_product(v1,v2)
!
!   write(*,'(a,6i10,f30.20)')'al1, ja1, jsv1, al2, ja2, jsv2, result_value = ', &
!&                         al_0, orb_0, atm_0, al_m, orb_m, atm_m, result_value
!
  end subroutine get_v_v_product
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Get eigen vector of the Krylov space, as a global vector
!
  subroutine get_kr_eig_vec(v, al_m, orb_m, atm_m)

    use M_la_gkrylov_hist, only : u_hst_str, v_mat_kr_str          !(unchanged)
    use elses_mod_orb2,    only : j2js,j2ja,js2j,n_tot_base        !(unchanged)
    use M_qm_domain,       only : jsv4jsd, njsd, atm_element, nval !(unchanged)
    use M_la_krgl_main,    only : kr_dim_str                       !(unchanged)
    implicit none
    integer,  intent(in)    :: al_m, orb_m, atm_m 
    real(8),  intent(inout) :: v(:)
    integer               :: j_src
    integer               :: ict4h
    integer               :: jsd1, jsv1, ja1, jj, j
    integer               :: kr_dim_max
!
    ict4h=1
!
    kr_dim_max=kr_dim_str(orb_m,atm_m)
    jj=0
    j_src=js2j(orb_m, atm_m)
    do jsd1=1,njsd(atm_m, ict4h)
      jsv1=jsv4jsd(jsd1, atm_m)
      do ja1=1,nval(atm_element(jsv1))
        jj=jj+1
        j=js2j(ja1, jsv1)
        v(j)=dot_product(v_mat_kr_str(1:kr_dim_max, al_m, orb_m, atm_m), & 
&                          u_hst_str(jj, 1:kr_dim_max, orb_m, atm_m))
      enddo
    enddo
!    
  end subroutine get_kr_eig_vec

end module M_la_gkrylov_output


