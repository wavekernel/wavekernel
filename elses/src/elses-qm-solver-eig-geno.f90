!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_solver_eig_geno
!
!
  use M_qm_domain, only : i_verbose, &
&      dhij, dsij, dbij, dpij, njsd, noav, atm_element, nval, jsv4jsd, &
&      temp_for_electron, chemical_potential, DOUBLE_PRECISION, ddsij
!
   private
!
! Public routines
   public qm_solver_eig_geno
   public plot_eigen_levels
   public plot_spectrum_from_eigen
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine qm_solver_eig_geno
!
     use M_config,     only : config !(unchanged )
     use elses_mod_orb2,     only : n_tot_base !(unchanged)
     use elses_arr_eig_leg,  only : atmp       !(CHANGED)
     use  M_wall_clock_time, only : get_system_clock_time
     use elses_mod_eig_leg,  only : n_base_eig_leg !(CHANGED)
     use M_eig_solver_center, only : set_density_matrix_mpi
!
     implicit none
     integer i_init
     real(8) :: elapse_time, elapse_time_bak
     integer lu
!    integer neig0, nval_max
     logical :: legacy_workflow ! with full matrix 
!
     lu=config%calc%distributed%log_unit
!
     if (trim(config%calc%solver%scheme) == 'eigen_mpi') then
       legacy_workflow = .false.
     else
       legacy_workflow = .true.
     endif
!
     if (i_verbose >= 1) then
       if (legacy_workflow) then
         if (lu > 0) write(lu,"(a)") '@@ qm_solver_eig_geno:legacy_workflow'
       else
         if (lu > 0) write(lu,"(a)") '@@ qm_solver_eig_geno:Non_legecy_workflow'
       endif
     endif   
!
     if (i_verbose >= 1) then
       call get_system_clock_time(elapse_time)
       elapse_time_bak=elapse_time
     endif   
!
     n_base_eig_leg=n_tot_base
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Set full matrices  (ONLY FOR legacy workflow)
!
     if (legacy_workflow) then
!
       i_init=0
       if (allocated(atmp) .eqv. .false.) i_init=1
!
       if (i_init == 1) then
!      if (i_verbose >= 1) then
!        write(*,*)' n_tot_base=',n_tot_base
!      endif   
         call elses_alloc_eig_leg(n_tot_base)
!         ----> allocations of working matrices 
!               in eigen-state solvers
       endif  
!
       if (i_verbose >= 1) then
         call get_system_clock_time(elapse_time)
         if (lu > 0) write(lu,"(a,f20.10)") ' TIME:qm_solver_eig_geno:initial alloc=',elapse_time-elapse_time_bak
         elapse_time_bak=elapse_time
       endif
!
       call copy_to_full_matrices
!      --> Copy the Hamiltonian and overlap matrices
!              into full-matrix arraies (atmp, atmp2)
!        INPUT : dhij, dsij 
!        OUTPUT: atmp (n,n)  (Hamiltonian, as full matrix)
!        OUTPUT: atmp2(n,n)  (Overlap,     as full matrix)
!
       if (i_verbose >= 1) then
         call get_system_clock_time(elapse_time)
         if (lu > 0) write(lu,"(a,f20.10)") ' TIME:qm_solver_eig_geno:copy_to_full =',elapse_time-elapse_time_bak
         elapse_time_bak=elapse_time
       endif   
!
     endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Generalized eigen value problem 
!     with matrix size of n = n_tot_base
!
     call set_eigen_states
!      --> Set the eigen states 
!        OUTPUT: eig2(n)    (Eigen levels)
!        OUTPUT (normal mode): atmp
!        OUTPUT (eigen_mpi mode): desc_eigenvectors(:), eigenvectors(:, :)  (Distributed eigen vectors A(i,k))
!
     if (i_verbose >= 1) then
       call get_system_clock_time(elapse_time)
       write(*,"(a,f20.10)") '  TIME:qm_solver_eig_geno:set_eigen    =',elapse_time-elapse_time_bak
       if (lu > 0) write(lu,"(a,f20.10)") '  TIME:qm_solver_eig_geno:set_eigen    =',elapse_time-elapse_time_bak
       elapse_time_bak=elapse_time
     endif   
!
     call elses_eig_chem_pot
!      --> Set the chemical potential and occupation number 
!               (legacy code)
!        INPUT : temp_for_electron (elenctron temperature)
!        OUTPUT: f_occ(n)
!
     if (i_verbose >= 1) then
       call get_system_clock_time(elapse_time)
       if (lu > 0) write(lu,"(a,f20.10)") ' TIME:qm_solver_eig_geno:set_chem_pot =',elapse_time-elapse_time_bak
       elapse_time_bak=elapse_time
     endif   
!
     if (trim(config%calc%solver%scheme) == 'eigen_mpi') then
       call set_density_matrix_mpi()
     else
       call set_density_matrix
     end if
!      --> Set the density matrix
!        OUTPUT: dbij : density matrix
!        OUTPUT: dpij : energy density matrix
!
     if (i_verbose >= 1) then
       call get_system_clock_time(elapse_time)
       if (lu > 0) write(lu,"(a,f20.10)") ' TIME:qm_solver_eig_geno:set_dens mat =',elapse_time-elapse_time_bak
       elapse_time_bak=elapse_time
     endif   
!
     if (i_verbose >= 60) then
       call plot_eigen_levels
     endif  
!      --> Plot eigen levels and so on (optional)
!
     if (i_verbose >= 1) then
       call get_system_clock_time(elapse_time)
       if (lu > 0) write(lu,"(a,f20.10)") ' TIME:qm_solver_eig_geno:plot eigen   =',elapse_time-elapse_time_bak
     endif   
!
   end subroutine qm_solver_eig_geno
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Plot eigen levels and so on, if they are defined (optional)
!
  subroutine plot_eigen_levels
    use M_config,     only : config !(unchanged )
    use elses_mod_phys_const, only : ev4au  !(unchanged )
    use elses_arr_eig_leg, only : atmp, eig2, f_occ  !(unchanged )
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged )
    implicit none
    integer :: i,j,k, step_count, i_enable
    real(8) :: sum_energy, sum_occupation
    real(8) :: sum_eigen_levels
!
    i_enable=1
    if (allocated(eig2)  .eqv. .false.) i_enable=0
    if (allocated(f_occ) .eqv. .false.) i_enable=0
    if (allocated(atmp)  .eqv. .false.) i_enable=0
!         
    if (i_enable == 0) then
       if (i_verbose >= 1) then
         write(*,*)'@@ plot_eigen_levels...skipped'
       endif   
       return
    endif   
!   
    step_count=config%system%structure%mdstep
    if (i_verbose >= 1) then
      write(*,*)'@@ plot_eigen_levels'
    endif   
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
    write(*,*)'Total elec. number (check): N_elec=',sum_occupation
    write(*,*)'Sum of eigen levels (as 2 Sum_i E_i ) [au]=',sum_energy
    write(*,*)'Sum of eigen levels (as 2 Sum_i E_i ) [eV]=',sum_energy*ev4au
!
!
    if (i_verbose >= 100) then
      write(*,*)'Plot coefficents: i,  A(i,k)'
      do k=1,n_tot_base
        write(*,*)'Coefficient of ',k,'-th eigen state'      
        do i=1,n_tot_base
          write(*,*)i,atmp(i,k)
        enddo   
      enddo   
    endif
!
  end subroutine plot_eigen_levels
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Copy the Hamiltonian and overlap matrices
!         into full-matrix arraies (atmp, atmp2)
!
  subroutine copy_to_full_matrices
!
    use elses_arr_eig_leg, only : atmp, atmp2
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base
    use elses_mod_js4jsv,   only : js4jsv
    use elses_mod_multi,    only : ict4h
!
    implicit none
    integer :: jsv2, nss2, js2, nval2, ja2, jsd1
    integer :: jsv1, nss1, js1, nval1, ja1
    integer :: ig, jg
!
    if (i_verbose >= 1) then
      write(*,*)'@@ copy_to_full_matrices'
    endif   
!
    if (.not. allocated(dsij)) then
      write(6,*)'ERROR:DSIJ is not allocated'
      stop
    endif
!
    atmp(:,:)=0.0d0
    atmp2(:,:)=0.0d0
!
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      js2=js4jsv(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd(jsv2,ict4h)
        jsv1=jsv4jsd(jsd1,jsv2)
        nss1=atm_element(jsv1)
        js1=js4jsv(jsv1)
        nval1=nval(nss1)
        do ja2=1,nval2
          jg=js2j(ja2,js2)
          if ((jg .le. 0) .or. (jg .gt. n_tot_base)) then
            write(6,*)'ERROR:JS2J:j=',jg
            stop
          endif
          do ja1=1,nval1
            ig=js2j(ja1,js1)
            if ((ig .le. 0) .or. (ig .gt. n_tot_base)) then
              write(6,*)'ERROR:JS2J:j=',ig
              stop
            endif
            atmp(ig,jg)=dhij(ja1,ja2,jsd1,jsv2)
!              ---> Hamiltonian
            atmp2(ig,jg)=dsij(ja1,ja2,jsd1,jsv2)
!              ---> Overlap
          enddo
        enddo
      enddo
    enddo
!
  end subroutine copy_to_full_matrices
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Set eigen states with 'elses_eig_mateig' (legacy code)
!    based on 'elses_seteig3'
!
  subroutine set_eigen_states
!
!   use elses_arr_eig_leg, only:atmp,atmp2,atmp3,atmpo
!
    use M_config,          only : config !(unchanged)
    use elses_arr_eig_leg
    implicit none
    integer :: imode
    logical :: flag_for_dump
!
!   flag_for_dump = .true.  ! 'true' is chosen only for a debug stage
    flag_for_dump = .false.  ! 'true' is chosen only for a debug stage
!
    if (trim(config%calc%solver%scheme) == 'eigen_mpi') then
      imode=2
      call eig_mateig_wrapper(imode)
    else   
      imode=2
      call elses_eig_mateig(imode)
    endif  
!      --> ATMP is replaced by eigen state
!      --> ATMP2 is also broken by Cholesky decomposition.
!      OUTPUT: atmp(:,:) : eigen states
!      OUTPUT: eig2(:)   : eigen levels
!
    if (flag_for_dump) call dump_eigen_pair
!
  end subroutine set_eigen_states
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Set eigen states compatible to the legecy code of 'elses_eig_mateig' 
!
  subroutine eig_mateig_wrapper(imode)
!   use elses_arr_eig_leg,      only : mat_a=>atmp, mat_b=>atmp2, eig_levels=>eig2 ! CHANGED
    use elses_mod_orb2,     only : n_tot_base
    use elses_arr_eig_leg,      only : atmp, eig_levels=>eig2, desc_eigenvectors, eigenvectors ! CHANGED
    use M_eig_solver_center, only : eig_solver_center ! routine
    use M_config,          only : config !(unchanged)                                                                        
    implicit none
    integer :: imode
    character(len=128) :: SEP_solver_wrk
    character(len=128) :: GS_transformation_wrk
    integer            :: blocksize_wrk
    integer            :: level_low_high(2)
    integer            :: log_unit_wrk
    integer            :: vec_size, ierr

    if (allocated(eig_levels)) then
      deallocate(eig_levels)
    end if    
    allocate(eig_levels(n_tot_base))    
!
    vec_size=size(eig_levels,1)
!
    log_unit_wrk=config%calc%distributed%log_unit
    SEP_solver_wrk=trim(config%calc%solver%eigen_mpi%SEP_solver)
    GS_transformation_wrk=trim(config%calc%solver%eigen_mpi%GS_transformation)
    blocksize_wrk=config%calc%solver%eigen_mpi%blocksize
    level_low_high(1)=config%calc%solver%eigen_mpi%level_lowest
    level_low_high(2)=config%calc%solver%eigen_mpi%level_highest
    imode=2
!   solver_scheme_wrk='scalapack'
    call eig_solver_center(imode, log_unit_wrk, trim(SEP_solver_wrk), & 
         trim(GS_transformation_wrk), blocksize_wrk, level_low_high, eig_levels, &
         desc_eigenvectors, eigenvectors)

!    if ((level_low_high(1) == 1) .and. (level_low_high(2) == vec_size )) then
!      if (.not. allocated(atmp)) then 
!        allocate(atmp(vec_size, vec_size),stat=ierr) 
!        if (ierr /= 0) then
!          stop 'Alloc error. atmp2 (eig_mateig_wrapper)' 
!        endif
!      endif
!      ierr=0
!      if (size(atmp,1) /= size(eig_vectors,1)) ierr=1
!      if (size(atmp,2) /= size(eig_vectors,2)) ierr=1
!      if (ierr /= 0) then
!         write(*,*) 'ERROR:imcompatible matrix size (eig_mateig_wrapper)'
!         write(*,*) ' size(atmp2,1)=', size(atmp,1)
!         write(*,*) ' size(atmp2,2)=', size(atmp,2)
!         write(*,*) ' size(eig_vectors,1)=', size(eig_vectors,1)
!         write(*,*) ' size(eig_vectors,2)=', size(eig_vectors,2)
!         stop
!      endif
!      atmp(:,:)=eig_vectors(:,:)
!    endif
    
!
!   write(*,*)'@@ eig_mateig_wrapper'
!    stop 'stop manually'
!
  end subroutine eig_mateig_wrapper
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Dump the eigen value and vectors
!
  subroutine dump_eigen_pair
    use M_lib_dst_info,     only : mpi_is_active, myrank, nprocs               !(unchanged)
    use elses_arr_eig_leg,  only : mat_a=>atmp, mat_b=>atmp2, eig_levels=>eig2 !(unchanged)
    use elses_mod_file_io,  only : vacant_unit                                 !(function)
    implicit none
    integer :: i, j, unit_num_e, unit_num_v
!
    if (myrank /= 0) then 
      write(*,*)'@@ Dump eigen pair (for debugging) .. is ignored (non-root node)'
      return
    endif   
!
    write(*,*)'@@ Dump eigen pair (for debugging) '
    write(*,*)'INFO: mpi_is_active  =', mpi_is_active
    write(*,*)'INFO: myrank, nprocs =', myrank, nprocs
!
    unit_num_e=vacant_unit()
    open(unit_num_e, file='test_data_eigen_levels.txt', status='unknown')
!
    unit_num_v=vacant_unit()
    open(unit_num_v, file='test_data_eigen_vectors.txt', status='unknown')
!
    write(*,*)'unit_num_e, unit_num_v=', unit_num_e, unit_num_v
!
!
    do j=1,size(mat_a,2)
      write(unit_num_e, '(i10, f30.20)') j, eig_levels(j)
      do i=1,size(mat_a,1)
        write(unit_num_v, '(2i10, f30.20)') i, j, mat_a(i,j)
      enddo   
    enddo   
!
    close(unit_num_e)
    close(unit_num_v)
!
    stop 'Stop (test mode) : dump_eigen_pair'
!
!
  end subroutine dump_eigen_pair
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Generate density matrix and energy density matrix
!    base on 'elses_eig_set_dens_mat'
!
!      BIJ(i,j) = sum_{k} A{i,k} A{j,k} f_k
!      PIJ(i,j) = sum_{k} A{i,k} A{j,k} f_k e_k
!
  subroutine set_density_matrix
    use elses_arr_eig_leg, only : atmp, eig2, f_occ
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base
    use elses_mod_js4jsv,   only : js4jsv
    use elses_mod_multi,    only : ict4h
!
    implicit none
    integer :: neig, neig0, neig_k
    integer :: jsv2, js2, ja2, jsd1
    integer :: jsv1, js1, ja1
    integer :: maxNumOrb, numOrb1, numOrb2
    integer :: ierr
!
    real(DOUBLE_PRECISION), allocatable, dimension(:,:)&
            :: atmp6, atmp7, atmp8, atmp9
    real(DOUBLE_PRECISION) :: w0,w1,w2,w3
    real(DOUBLE_PRECISION) :: EPSILON=1d-14
!
    if (.not. allocated(atmp)) then
      write(*,*) 'ERROR(set_density_matrix):atmp is not allocated'
      stop
    endif
!
!ccccccccccccccccccccccccccccccccccccccccccc
! @@ Local parameter setting
!
    neig0=n_tot_base
    neig=neig0
!ccccccccccccccccccccccccccccccccccccccccccc
! @@ Initial clearance
!
    dbij(:,:,:,:)=0.0d0
    dpij(:,:,:,:)=0.0d0
!
!ccccccccccccccccccccccccccccccccccccccccccc
!   @@ Allocation 
!
    do ja1=1, neig0
       if(f_occ(ja1) < EPSILON) exit
    end do
    neig_k=min(ja1,neig0)
!     ----> upper limit of sum_k
!            ( k = 1, 2, ..., neig_k )
!
    maxNumOrb=maxval(nval)
    allocate (atmp6(neig_k,neig0))
    !atmp6=transpose(atmp)
!$omp parallel do private(w0,w1,w2,w3)
    do ja1=1, neig_k-3, 4
       w0=sqrt(f_occ(ja1))
       w1=sqrt(f_occ(ja1+1))
       w2=sqrt(f_occ(ja1+2))
       w3=sqrt(f_occ(ja1+3))
       do ja2=1, neig0
          atmp6(ja1,  ja2)=w0*atmp(ja2,ja1)
          atmp6(ja1+1,ja2)=w1*atmp(ja2,ja1+1)
          atmp6(ja1+2,ja2)=w2*atmp(ja2,ja1+2)
          atmp6(ja1+3,ja2)=w3*atmp(ja2,ja1+3)
       end do
    end do
!$omp end parallel do
!$omp parallel do private(w0)
    do ja1=neig_k-mod(neig_k,4)+1, neig_k
       w0=sqrt(f_occ(ja1))
       do ja2=1, neig0
          atmp6(ja1,  ja2)=w0*atmp(ja2,ja1)
       end do
    end do
!$omp end parallel do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calc. of L-matrix and P-matrix as dbij and dpij
!
!      L(i,j) = sum_k C(i,k) f(k) C(j,k)
!             = sum_k C(i,k) D(k,j) 
!
!      P(i,j) = sum_k C(i,k) f(k) E(k) C(j,k)
!             = sum_k C(i,k) Q(k,j) 
!
      if (i_verbose >= 1) then
        write(6,*)'calc.of L-matrix and P-matrix'
      endif  
!
!$omp  parallel default(shared) &
!$omp& private(     jsv2,js2,ja2,numOrb2) &
!$omp& private(jsd1,jsv1,js1,ja1,numOrb1) &
!$omp& private(atmp7,atmp8,atmp9)
      allocate (atmp7(neig_k,maxNumOrb))
      allocate (atmp8(maxNumOrb,maxNumOrb))
      allocate (atmp9(maxNumOrb,maxNumOrb))
!$omp do
      do jsv2=1,noav
         js2=js4jsv(jsv2)
         numOrb2=nval(atm_element(jsv2))
         do ja2=1, numOrb2
            atmp7(1:neig_k,ja2)=atmp6(1:neig_k,js2j(ja2,js2))*eig2(1:neig_k)
         end do
         do jsd1=1,njsd(jsv2,ict4h)
            jsv1=jsv4jsd(jsd1,jsv2)
            js1=js4jsv(jsv1)
            numOrb1=nval(atm_element(jsv1))
            ! Assumption : js2j(:,js1) is a set of contiguous and ascending integers
            call dgemm('T','N',numOrb1,numOrb2,neig_k,&
                 1d0,atmp6(1,js2j(1,js1)),neig_k,atmp6(1,js2j(1,js2)),&
                 neig_k,0d0,atmp8,maxNumOrb)
            call dgemm('T','N',numOrb1,numOrb2,neig_k,&
                 1d0,atmp6(1,js2j(1,js1)),neig_k,atmp7,neig_k,0d0,atmp9,maxNumOrb)
            do ja2=1, numOrb2
               do ja1=1, numOrb1
                  dbij(ja1,ja2,jsd1,js2)=atmp8(ja1,ja2)
                  dpij(ja1,ja2,jsd1,js2)=atmp9(ja1,ja2)
               enddo
            enddo
               ! original code
!!$               do ja2=1, numOrb2
!!$                  do ja1=1, numOrb1
!!$                     dbij(ja1,ja2,jsd1,js2)=&
!!$                          dot_product(atmp6(1:neig_k,js2j(ja1,js1)),&
!!$                                      atmp6(1:neig_k,js2j(ja2,js2)))
!!$                     dpij(ja1,ja2,jsd1,js2)=&
!!$                          dot_product(atmp6(1:neig_k,js2j(ja1,js1)),&
!!$                                      atmp7(1:neig_k,ja2))
!!$                  enddo
!!$               enddo
         end do
      enddo
!$omp end do
      deallocate(atmp7,atmp8,atmp9)
!$omp end parallel
!
  end subroutine set_density_matrix
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Plot LDOS from eigen states (optional)
!
!
  subroutine plot_spectrum_from_eigen
    use M_qm_domain, only : nos !(unchanged)
    use M_io_ctrl_dos_plot,   only : initial_preparation_for_dos !(routine)
    use M_lib_phys_const,     only : pi !(constant)
    use M_config,             only : config !(unchanged )
    use elses_mod_phys_const, only : ev4au  !(unchanged )
    use elses_arr_eig_leg, only : atmp, atmp2, eig2, f_occ  !(unchanged )
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged )
    use elses_mod_js4jsv,   only : js4jsv                  !(unchanged )
    use elses_mod_file_io, only : vacant_unit
    use M_lib_math_func,      only : Fermi_Dirac_Func
    use M_qm_geno,            only : set_hamiltonian_and_overlap_geno
    use M_qm_domain,          only : ham_tb0, dham_tb0, DOUBLE_PRECISION
    use M_pair_list_for_cohp, only : set_pair_num_for_cohp  !(routine)
    use M_pair_list_for_cohp, only : set_pair_list_for_cohp !(routine)
!
    implicit none
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
    real(DOUBLE_PRECISION), allocatable  :: tnos(:)
    real(DOUBLE_PRECISION), allocatable  :: enos(:)
    real(DOUBLE_PRECISION), allocatable  :: pnos(:,:)
    real(DOUBLE_PRECISION), allocatable  :: lnos(:,:)
    real(DOUBLE_PRECISION), allocatable  :: tpnos(:,:,:)
    real(DOUBLE_PRECISION), allocatable  :: eigen_states(:,:)
    real(DOUBLE_PRECISION), allocatable  :: S_eigen_states(:,:)
    real(DOUBLE_PRECISION) :: nos_prv, nos_tmp, dos_tmp, ene_tmp, ene_tot
    real(DOUBLE_PRECISION) :: nos_prv_l(0:2), nos_tmp_l(0:2), dos_tmp_l(0:2)
!
!   For variables for COHP
    integer :: jsv2, nss2, nval2, jsd1, jsv1, nss1, nval1, ja2, ja1, ict4h
    integer :: pair_kind, num_pair_kind
    integer :: js2, js1, jg, ig, jj_prev
    real(DOUBLE_PRECISION), allocatable  :: itcohp(:,:)
    real(DOUBLE_PRECISION) :: eig_k
    integer :: ite_cohp, ite_cohp_max
    real(DOUBLE_PRECISION), allocatable  :: bondMask(:)
    real(DOUBLE_PRECISION), allocatable  :: dsijBK(:,:,:,:), ddsijBK(:,:,:,:,:) 
!
    integer :: num_pair_list, ipair
    integer, allocatable  :: atm_pair_list(:,:)
    character(len=32)     :: num_pair_list_chara
    character(len=32)     :: format_cohp_output
!
    integer :: nval_max, elm_index
    integer :: l_sym
!
    nval_max=maxval(nval)
!
    i_enable=1
    if (.not. allocated(eig2))  i_enable=0
    if (.not. allocated(f_occ)) i_enable=0
    if (.not. allocated(atmp))  i_enable=0
    if (.not. allocated(atmp2)) i_enable=0
!         
    call initial_preparation_for_dos(imode, number_energy_mesh, energy_btm, energy_top, width_of_peak)
    if (imode == 0) i_enable=0
!
    if (i_enable == 0) then
       if (i_verbose >= 1) then
         write(*,*)'@@ plot_spectrum_from_eigen ... skipped'
       endif   
       return
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate (eigen_states(n_tot_base,n_tot_base),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'allocation error'
      stop
    endif
!
    allocate (S_eigen_states(n_tot_base,n_tot_base),stat=ierr)
    if( ierr .ne. 0) then
      write(6,*)'allocation error'
      stop
    endif
!
    eigen_states(:,:)=atmp(:,:)
    S_eigen_states(:,:)=0.0d0
!
     call copy_to_full_matrices
!      --> Copy the Hamiltonian and overlap matrices
!              into full-matrix arraies (atmp, atmp2)
!        INPUT : dhij, dsij 
!        OUTPUT: atmp (n,n)  (Hamiltonian, as full matrix)
!        OUTPUT: atmp2(n,n)  (Overlap,     as full matrix)
!
    S_eigen_states=matmul(atmp2,eigen_states)
!
!   
!
    if (i_verbose >= 1) then
       write(*,*)'@@ plot_ldos_from_eigen (experimental routine)'
       write(*,*)' i_enable=',i_enable
       write(*,*)' imode  =',imode
       write(*,*)' number of energy mesh  =',number_energy_mesh
       write(*,*)' energy bottom   [eV]   =',energy_btm*ev4au
       write(*,*)' energy top      [eV]   =',energy_top*ev4au
       write(*,*)' width of peak      [eV] =',width_of_peak*ev4au
       write(*,*)' chemical potential [au] =',chemical_potential
       write(*,*)' chemical potential [eV] =',chemical_potential*ev4au
    endif   
!
!   if (temp_for_electron <= 1.0d-10) then
!     write(*,*)' ERROR(plot_ldos_from_eigen)'
!     write(*,*)' temperature for electron (au)=',temp_for_electron
!     stop
!   endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate (energy_mesh(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!    
    allocate (tnos(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    tnos(:)=0.0d0
!
    allocate (enos(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    enos(:)=0.0d0
!
    allocate (pnos(0:number_energy_mesh,n_tot_base),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc:pnos'
    pnos(:,:)=0.0d0
!
    allocate (lnos(0:number_energy_mesh,noav),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc:lnos'
    lnos(:,:)=0.0d0
!
    allocate (tpnos(0:number_energy_mesh, 0:2, nos),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc:lnos'
    tpnos(:,:,:)=0.0d0
!
    de=(energy_top-energy_btm)/dble(number_energy_mesh)
!      -----> de : energy interval [au]
!
    do jj=0,number_energy_mesh
      energy_mesh(jj)=energy_btm+de*dble(jj)
!      -----> energy mesh in eV
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    step_count=config%system%structure%mdstep
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
    iunit=vacant_unit()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Eigen levels
!
    open(iunit, file='output_eigen_levels.txt', status='unknown')
!
    write(*,*)'Save eigen levels into output_eigen_levels.txt'
!
    do k=1,n_tot_base
       write(iunit,'(i10,2f20.10)') k, eig2(k)*ev4au, f_occ(k)
    enddo
!    
    close(iunit)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Total DOS
!
   if (imode >= 1) then
!
    open(iunit, file='output_tdos.txt', status='unknown')
!
    write(*,*)'Save TDOS into output_tdos.txt'
!
    write(*,*) 'TDOS Result : Chem_pot [eV] =',chemical_potential*ev4au
!   write(iunit,'(a,f20.10)') 'TDOS Result : Chem_pot [eV] =',chemical_potential*ev4au
!   write(iunit,'(a)') 'tdos=: energy index, energy, DOS, IDOS'
!
!   nos(:)=0.0d0
!   do k=1,n_tot_base
!     do jj=0,number_energy_mesh
!       x=(energy_mesh(jj)-eig2(k))/width_of_peak
!       step_fn=datan(x)/pi+0.5d0
!       nos(jj)=nos(jj)+step_fn*2.0d0
!      enddo   
!   enddo   
!
    tnos(:)=0.0d0
    enos(:)=0.0d0
!$omp  parallel default(shared) &
!$omp& private(k, jj) &
!$omp& private(x, step_fn)
!$omp  do schedule(static)
    do jj=0,number_energy_mesh
      do k=1,n_tot_base
        x=(eig2(k)-energy_mesh(jj))/width_of_peak
        step_fn=Fermi_Dirac_Func(x)
        tnos(jj) = tnos(jj) + step_fn*2.0d0
        enos(jj) = enos(jj) + step_fn*2.0d0*eig2(k) ! enos in au
      enddo   
    enddo
!$omp end do
!$omp end parallel

!
    nos_prv=tnos(0)
!   ene_tot=0.0d0
    do jj=0,number_energy_mesh
      nos_tmp=tnos(jj)
      dos_tmp=(nos_tmp-nos_prv)/(de*ev4au)
!     ene_tmp=energy_mesh(jj)
!     ene_tot=ene_tot+ene_tmp*dos_tmp*de
      write(iunit,'(a,i7,4f20.10)')'tdos=', jj, energy_mesh(jj)*ev4au, dos_tmp, tnos(jj), enos(jj)*ev4au
      nos_prv=nos_tmp
    enddo   
!    
    close(iunit)
!
  endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ total COHP 
!
  if ((imode == 3 ) .or. (imode == 4 )) then
!
   if (imode == 3) ite_cohp_max=1
   if (imode == 4) ite_cohp_max=2
!  
   do ite_cohp=1,ite_cohp_max
!
     if (ite_cohp == 2) then
       allocate(dsijBK(size(dsij,1),size(dsij,2),size(dsij,3),size(dsij,4)))
       allocate(ddsijBK(3,size(dsij,1),size(dsij,2),size(dsij,3),size(dsij,4)))
       allocate(bondMask(size(dhij,1)))
        dsijBK = dsij
       ddsijBK =ddsij
       bondMask(:)=1d0
       bondMask(4)=0d0 !! pp-\pi is masked
       !! bondMask/ss-\sigma,n.a,pp-\sigma,pp-\pi,n.a,n.a,dd-\sigma,dd-\pi,dd-\delta/
       call set_hamiltonian_and_overlap_geno(bondMask)
       dhij(:,:,:,:)   = ham_tb0(:,:,:,:)
       dsij=dsijBK
       ddsij=ddsijBK
       deallocate(bondMask)
       deallocate(dsijBK)
       deallocate(ddsijBK)
     endif  
!
     call copy_to_full_matrices
!      --> Copy the Hamiltonian and overlap matrices
!              into full-matrix arraies (atmp, atmp2)
!        INPUT : dhij, dsij 
!        OUTPUT: atmp (n,n)  (Hamiltonian, as full matrix)
!        OUTPUT: atmp2(n,n)  (Overlap,     as full matrix)
!
!
    call set_pair_num_for_cohp(num_pair_list)
    if (i_verbose >=1) then
      write(*,*)'num_pair_list=',num_pair_list
    endif   
!
    if (num_pair_list == 0) then
      num_pair_kind=3
      allocate (itcohp(0:number_energy_mesh,num_pair_kind),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc'
      if (ite_cohp == 2) then
        open(iunit, file='output_tcohp_wo_pppi.txt', status='unknown')
        write(*,*)'Save TCOHP into output_tcohp_wo_pppi.txt'
      else
        open(iunit, file='output_tcohp.txt', status='unknown')
        write(*,*)'Save TCOHP into output_tcohp.txt'
      endif   
    else
      allocate (itcohp(0:number_energy_mesh,num_pair_list),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc:itcohp'
      if (allocated(atm_pair_list)) then
        deallocate (atm_pair_list,stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc:atm_pair_list'
      endif
      allocate (atm_pair_list(2,num_pair_list),stat=ierr)
      if (ierr /= 0) then
        write(*,*) 'Abort:ERROR in alloc:atm_pair_list:um_pair_list=',num_pair_list
        stop
      endif
      call set_pair_list_for_cohp(num_pair_list, atm_pair_list)
      if (ite_cohp == 2) then
        open(iunit, file='output_cohp_wo_pppi.txt', status='unknown')
        write(*,*)'Save TCOHP into output_cohp_wo_pppi.txt'
      else
        open(iunit, file='output_cohp.txt', status='unknown')
        write(*,*)'Save TCOHP into output_cohp.txt'
      endif   
    endif   
    itcohp(:,:)=0.0d0
!
!
    ict4h=1
    do jsv2=1,noav
      js2=js4jsv(jsv2)
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      write(*,*)'jsv2, njsd(jsv2,ict4h)=',jsv2,njsd(jsv2,ict4h)
      do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         js1=js4jsv(jsv1)
!        write(*,*)'jsv2, jsd1, jsv1=',jsv2,jsd1,jsv1
         if (jsv1 == jsv2) cycle
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         if (num_pair_list == 0) then
           pair_kind=0
           if ((nss1 == 1) .and. (nss2 == 1)) pair_kind = 1
           if ((nss1 == 1) .and. (nss2 == 2)) pair_kind = 2
           if ((nss1 == 2) .and. (nss2 == 1)) pair_kind = 2
           if ((nss1 == 2) .and. (nss2 == 2)) pair_kind = 3
           if (pair_kind == 0) then
             write(*,*)'STOP in COHP calc. : pair_kind=',pair_kind
             stop
           endif   
         else
           pair_kind=0 
           do ipair=1,num_pair_list
             if ((atm_pair_list(1,ipair)==jsv1) .and. (atm_pair_list(2,ipair)==jsv2)) pair_kind=ipair
             if ((atm_pair_list(1,ipair)==jsv2) .and. (atm_pair_list(2,ipair)==jsv1)) pair_kind=ipair
           enddo   
         endif   
         if (pair_kind == 0) then 
           cycle 
         else
           write(*,*)'pair_kind, jsv1, jsv2=',pair_kind, jsv1, jsv2
         endif   
!        write(*,*)'jsv1, pair_kind=',jsv1,pair_kind
         do ja2=1,nval2
           jg=js2j(ja2,js2)
           do ja1=1,nval1
             ig=js2j(ja1,js1)
             do k=1,n_tot_base
!              if ((k <= 0) .or. (k > n_tot_base)) then
!                write(*,*)'ERROR(a)!:k=',k
!                stop
!              endif   
               weight=eigen_states(ig,k)*atmp(ig,jg)*eigen_states(jg,k)
               eig_k=eig2(k)
!$omp  parallel default(shared) &
!$omp& private(jj) &
!$omp& private(x, step_fn)
!$omp  do schedule(static)
               do jj=0,number_energy_mesh
!                if ((k <= 0) .or. (k > n_tot_base)) then
!                  write(*,*)'ERROR!:k=',k
!                  stop
!                endif   
                 x=(eig_k-energy_mesh(jj))/width_of_peak
!                write(*,*)'jj,k,eig,x=',jj,k,eig2(k),x
                 step_fn=Fermi_Dirac_Func(x)
!                write(*,*)'cohp-data=',jj,x,weight,step_fn
                 itcohp(jj,pair_kind)=itcohp(jj,pair_kind)+weight*step_fn  ! ICOHP in au
               enddo   
!$omp end do
!$omp end parallel
             enddo  
           enddo
         enddo
      enddo   
    enddo
!
    write(num_pair_list_chara,*) 2*num_pair_list+1
!
    num_pair_list_chara=adjustl(num_pair_list_chara)
    format_cohp_output='(a,i7,'//trim(adjustl(num_pair_list_chara))//'f20.10)'
    write(*,*) 'format_cohp_output=',trim(format_cohp_output)
!
    jj_prev=0
    do jj=0,number_energy_mesh
      if (jj == 0) then
        jj_prev = 0
      else
        jj_prev = jj-1
      endif   
      write(iunit,format_cohp_output)'cohp=', jj, energy_mesh(jj)*ev4au,  &
&           (itcohp(jj,1:num_pair_list)-itcohp(jj_prev,1:num_pair_list))/de, itcohp(jj,1:num_pair_list)*ev4au
    enddo   

    deallocate (itcohp)
!
    close(iunit)
!
    open(iunit, file='output_cohp_pair_list.txt', status='unknown')
    do ipair=1,size(atm_pair_list,2)
      write(iunit,'(3i10)')ipair, atm_pair_list(1:2,ipair) 
    enddo   
    close(iunit)
!
   enddo
!
!  write(*,*)'Stop manually (T.Hoshi)'
!  stop
!
  endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Partial DOS
!
  if (imode ==2 ) then
!
    open(iunit, file='output_pdos.txt', status='unknown')
!
!
    write(iunit,'(a)') 'pdos=: atom index, orb. index, energy index, energy, PDOS, IPDOS'
!
    pnos(:,:)=0.0d0
    do j=1,n_tot_base
      do k=1,n_tot_base
        weight=eigen_states(j,k)*S_eigen_states(j,k)
!$omp  parallel default(shared) &
!$omp& private (jj, x, step_fn)
!$omp  do schedule(static)
        do jj=0,number_energy_mesh
          x=(energy_mesh(jj)-eig2(k))/width_of_peak
          step_fn=datan(x)/pi+0.5d0
          pnos(jj,j)=pnos(jj,j)+step_fn*weight*2.0d0
        enddo
!$omp end do
!$omp end parallel
      enddo  
    enddo
!

    do j=1,n_tot_base
      write(iunit,'(a,2f20.10)') 'PDOS Result : Hami_diag, Chem_pot [eV] =',atmp(j,j)*ev4au, chemical_potential*ev4au
      nos_prv=pnos(0,j)
      atm_index=j2js(j)
      orb_index=j2ja(j)
      do jj=0,number_energy_mesh
        nos_tmp=pnos(jj,j)
        dos_tmp=(nos_tmp-nos_prv)/(de*ev4au)
        write(iunit,'(a,3i7,3f20.10)')'pdos=',atm_index,orb_index,jj,energy_mesh(jj)*ev4au,dos_tmp,nos_tmp
        nos_prv=nos_tmp
      enddo   
    enddo
!  
    close(iunit)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Local DOS
!
    open(iunit, file='output_ldos.txt', status='unknown')
!
    write(iunit,'(a)') 'ldos=: atom index, energy index, energy, LDOS, ILDOS'
!
    lnos=0.0d0
    do j=1,n_tot_base
      atm_index=j2js(j)
      lnos(:,atm_index)=lnos(:,atm_index)+pnos(:,j)
    enddo
!  
    do atm_index=1,noav
      nos_prv=lnos(0,atm_index)
      do jj=0,number_energy_mesh
        nos_tmp=lnos(jj,atm_index)
        dos_tmp=(nos_tmp-nos_prv)/(de*ev4au)
        write(iunit,'(a,2i7,3f20.10)')'ldos=',atm_index,jj,energy_mesh(jj)*ev4au,dos_tmp,nos_tmp
        nos_prv=nos_tmp
      enddo  
    enddo
!
    close(iunit)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Total TPDOS
!
    open(iunit, file='output_tpdos.txt', status='unknown')
!
    write(iunit,'(a)') 'tpdos=: element index, energy index, energy, TPDOS(s,p,d), TPNOS(s,p,d)'
!
    tpnos(:,:,:)=0.0d0
    do j=1,n_tot_base
      atm_index=j2js(j)
      orb_index=j2ja(j)
      elm_index=atm_element(atm_index)
      write(*,*)' tpnos-calc=',j, atm_index, elm_index, orb_index
      select case(orb_index)
        case (1)     ! s orbital
          l_sym=0
        case (2:4)   ! p orbitals
          l_sym=1
        case (5:9)   ! d orbitals
          l_sym=2
        case default
          write(*,*)'ERROR in select(plot_spectrum_from_eigen)' 
          stop
      end select   
      tpnos(:, l_sym, elm_index) = tpnos(:, l_sym, elm_index) + pnos(:,j)
    enddo
!
    do elm_index=1,nos
      nos_prv_l(0:2)=tpnos(0, 0:2, elm_index)
      do jj = 0, number_energy_mesh
        nos_tmp_l(0:2)=tpnos(jj, 0:2, elm_index)
        dos_tmp_l(0:2)=(nos_tmp_l(0:2)-nos_prv_l(0:2))/(de*ev4au)
        write(iunit,'(a,2i7,7f20.10)')'tpdos=',elm_index, jj, energy_mesh(jj)*ev4au, &
&                                         dos_tmp_l(0:2), nos_tmp_l(0:2)  
        nos_prv_l(0:2)=nos_tmp_l(0:2)
      enddo  
    enddo
!
    close(iunit)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Deallocation
!
    deallocate (tnos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
    deallocate (pnos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
    deallocate (lnos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
    deallocate (tpnos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
  endif
!
  end subroutine plot_spectrum_from_eigen
!
!
end module M_qm_solver_eig_geno

