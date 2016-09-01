!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_output
!
! Module variables : changed 
!   --> ETB (in qm_geno_output)
!
  use M_qm_domain
!
   private
!
! Public routines
   public qm_geno_output
   public qm_calc_etb_atom
   public qm_calc_ecsc_atom
   public qm_calc_ecsc
   public qm_calc_ecsc_list
   public output_for_eigen_solver
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Output routine
!      Module variables : changed 
!        --> ETB (as 2 Tr [rho H_{TB0}])
!
  subroutine qm_geno_output
!
    use M_config, only : config
    use M_qm_solver_eig_geno, only : plot_eigen_levels  ! (routine)
    use elses_mod_ene,     only : etb,ecsc ! (CHANGED)
!   use elses_mod_md_dat, only : itemd, itemdmx
    use M_io_ctrl_output_eigen, only : set_filename_for_save_wfn ! (routine)
    use M_output_eigenstates, only : output_eigenstates          ! (routine)
!   use M_io_ctrl_output_eigen, only : init_for_plot_wavefunction
    use M_qm_output_icohp, only : qm_output_icohp ! (subroutine)
    use M_qm_domain,       only : dbij !(unchanged)
    use elses_mod_md_dat,  only : final_iteration !(unchanged)
!
!   use M_output_participation, only : calc_participation
!
    implicit none
    integer      :: imode, n_csc_loop
    real(kind=8) :: value_of_etb0, dddd, value_of_ecsc
    integer      :: imode_icohp
    integer      :: log_unit
!
    character(len=8) :: CSC_METHOD
    character(len=100) :: filename_wfn
    CSC_METHOD = config%calc%genoOption%CSC_method
    n_csc_loop = config%calc%genoOption%CSC_max_loop_count
    log_unit   = config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_geno_output'
      write(*,*)'  CSC_METHOD  =',CSC_METHOD
      write(*,*)'  N_CSC_loop  =',n_csc_loop
    endif  
!
    if (.not. allocated(dbij)) then
      write(*,*)'WARN-INFO : global dens mat is not allocated'
      write(*,*)'WARN-INFO : qm_geno_output is skipped'
      return
    endif
!   
    if (.not. allocated(dsij)) then
      write(*,*)'WARN-INFO : global overlap mat is not allocated'
      write(*,*)'WARN-INFO : qm_geno_output is skipped'
      return
    endif
!   
    if (.not. allocated(dhij)) then
      write(*,*)'WARN-INFO : global hamiltonian mat is not allocated'
      write(*,*)'WARN-INFO : qm_geno_output is skipped'
      return
    endif
!   
!   if (i_verbose >= 1) then
!     write(*,'(a,i10,f20.10)')'  sum(dbij)   =',itemd,sum(dbij)
!     write(*,'(a,i10,f20.10)')'  sum(dsij)   =',itemd,sum(dsij)
!     write(*,'(a,i10,f20.10)')'  sum(dhij)   =',itemd,sum(dhij)
!   endif  
!
!   imode=3
!   call qm_calc_etb(value_of_etb0,imode)
!   write(6,*)'ETB0 as 2 Tr[rho H_{TB0}] =',value_of_etb0
!   write(6,*)' --> set as ETB'
!
    ecsc=0.0d0
    if(CSC_METHOD == "ELSTNER" .and. n_csc_loop >0) then
       imode=3
       call qm_calc_etb(value_of_etb0,imode)
       etb=value_of_etb0
       if (i_verbose >= 1) then
         write(6,*)'ETB0 as 2 Tr[rho H_{TB0}] =',value_of_etb0
         write(6,*)' --> set as ETB'
       endif  
       call qm_calc_ecsc(value_of_ecsc)
       ecsc=value_of_ecsc
       if (i_verbose >= 1) then
         write(6,*) 'n_csc_loop=',n_csc_loop
         write(6,*) 'Ecsc as 0.5*sum_{alpha,beta}gamma*dq_{alpha}*dq_{beta}=', value_of_ecsc
         write(6,*) '--> set as ECSC'
       endif  
    else
       imode=1
       call qm_calc_etb(value_of_etb0,imode)
       etb=value_of_etb0
       if (i_verbose >= 1) then
         write(6,*)'ETB0 as 2 Tr[rho H_{TB}] =',value_of_etb0
         write(6,*)' --> set as ETB (non CSC case)'
       endif  
    end if
!
    if (i_verbose >= 40) then
      call plot_eigen_levels
!        ---> plot eigen levels, if they are defined (optional, eigen solver only)
    endif   
!
    if (i_verbose >= 1) then
      dddd=0.0d0
      imode=2
      call qm_calc_etb(dddd,imode)
      write(6,*)'ETB  as 2 Tr[pi S] =',dddd
      imode=1
      call qm_calc_etb(dddd,imode)
      write(6,*)'ETB  as 2 Tr[rho H] =',dddd
    endif
!  
    imode_icohp=1
    call qm_output_icohp(imode_icohp)  
!        ---> plot ICOHP into file, at specific steps.
!
!   if (.not. final_iteration) then
!     call set_filename_for_save_wfn(filename_wfn)
!     if (trim(filename_wfn) /= '') call output_eigenstates(filename_wfn)
!   endif
!
  end subroutine qm_geno_output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine output_for_eigen_solver
!
    use M_config,               only : config  !(unchanged)
    use elses_mod_md_dat,       only : final_iteration !(unchanged)
    use M_io_ctrl_output_eigen, only : set_filename_for_save_wfn ! (routine)
    use M_output_eigenstates,   only : output_eigenstates
    use M_output_participation, only : calc_participation
    use M_qm_solver_eig_geno,   only : plot_spectrum_from_eigen  ! (routine)
    use M_lib_dst_info,         only : root_node        !(unchanged)
    implicit none
    character(len=100) :: filename_wfn
    integer            :: log_unit
!
    log_unit = config%calc%distributed%log_unit
!
    if (i_verbose >= 1) then
      if (log_unit > 0) write(log_unit,*) '@@@ output_for_eigen_solver'
    endif
!   
    call set_filename_for_save_wfn(filename_wfn)
!
    if (i_verbose >= 1) then
      if (log_unit > 0) write(log_unit,*) '  filename_wfn=',filename_wfn
    endif
!
    if (trim(filename_wfn) /= '') call output_eigenstates(filename_wfn)
!
    if (final_iteration) then
      !if (i_verbose >= 1) then  ! Comment out because i_verbose is not global and distributed operation fails.
      if (trim(filename_wfn) /= '') call calc_participation
      !endif
    endif
    if (root_node) then
      if (final_iteration) then
        call plot_spectrum_from_eigen
!          ---> plot LDOS from eigen states,
!               if they are defined (optional, eigen solver only)
      endif
    end if
!
  end subroutine output_for_eigen_solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calculate electronic energy in trace form
!           ETB = 2 Tr [rho H]        ( in imode = 1 )
!           ETB = 2 Tr [pi S]         ( in imode = 2 )
!           ETB = 2 Tr [rho H_{TB0}]  ( in imode = 3 )
!
!   Used module variable ( unchanged )
!               : dbij (as rho)
!               : dpij (as pi)
!               : dhij (as H)
!               : dsij (as S)
!               : ham_tb0 (as H_{TB0})
!
  subroutine qm_calc_etb(value_of_etb,imode)
!
!
    implicit none
    integer      :: imode
    integer, parameter  :: ict4h=1
    real(kind=8) :: value_of_etb
    integer      :: jsv2, nss2, nval2, ja2, jsd1
    integer      :: jsv1, nss1, nval1, ja1
    real(kind=8) :: ddsum, ddd1, ddd2
!
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_calc_etb:imode=',imode
    endif
!  
    ddsum=0.0d0  
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         do ja2=1,nval2
           do ja1=1,nval1
             if (imode == 1) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=dhij(ja1,ja2,jsd1,jsv2)
!              write(*,*)ja1,ja2,jsd1,jsv2,ddd1,ddd2
             endif  
             if (imode == 2) then
               ddd1=dpij(ja1,ja2,jsd1,jsv2)
               ddd2=dsij(ja1,ja2,jsd1,jsv2)
             endif
             if (imode == 3) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=ham_tb0(ja1,ja2,jsd1,jsv2)
             endif  
             ddsum=ddsum+ddd1*ddd2
           enddo  
         enddo  
      enddo  
    enddo  
    ddsum=ddsum*2.0d0
    value_of_etb=ddsum
!
      if (i_verbose >= 1) then
        if (imode == 1) then
          write(6,*)'ETB as 2 Tr[rho H] =',ddsum
        endif   
        if (imode == 2) then
          write(6,*)'ETB as 2 Tr[pi S]  =',ddsum
        endif   
        if (imode == 3) then
          write(6,*)'ETB as 2 Tr[rho H_{TB0}] =',ddsum
        endif   
      endif  

  end subroutine qm_calc_etb
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine qm_calc_ecsc(value_of_ecsc,plot_atom_energy)
    
    use M_config, only : config !(unchanged)
    use M_qm_domain, only : e_num_on_atom !(unchanged)
    use M_lib_phys_const,    only : ev4au !(unchaged)
!
    implicit none
    real(8), intent(out) :: value_of_ecsc
    logical, intent(in), optional :: plot_atom_energy
    integer :: jsv1,jsv2
    real(8) :: value_of_ecsc_atm, value_of_ecsc_atm_onsite
    integer :: lu

    lu = config%calc%distributed%log_unit
    value_of_ecsc = 0.0d0

    do jsv1=1,noav
      value_of_ecsc_atm=0.0d0
      call qm_calc_ecsc_atom(jsv1, value_of_ecsc_atm,value_of_ecsc_atm_onsite)
      value_of_ecsc = value_of_ecsc + value_of_ecsc_atm
      if (present(plot_atom_energy)) then
        if (plot_atom_energy) then 
          if (lu > 0) then
            write(lu,'(a,i10,3f20.10)')'atom index, charge, CSC energy, CSC onsite energy [eV] =',jsv1, &
&              e_num_on_atom(jsv1), value_of_ecsc_atm*ev4au, value_of_ecsc_atm_onsite*ev4au
          endif
        endif
      endif
    end do
    
  end subroutine qm_calc_ecsc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate electronic energy for each atom
!           ETB_i = 2 ( rho H )_{ii} 
!
  subroutine qm_calc_etb_atom(atm_index,value_of_etb)
!
    implicit none
    integer,      intent(in)  :: atm_index
    real(kind=8), intent(out) :: value_of_etb
    integer      :: imode
    integer, parameter  :: ict4h=1
    integer      :: jsv2, nss2, nval2, ja2, jsd1
    integer      :: jsv1, nss1, nval1, ja1
    real(kind=8) :: ddsum, ddd1, ddd2
!
    imode=1
    ddsum=0.0d0  
    jsv2=atm_index
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         do ja2=1,nval2
           do ja1=1,nval1
             if (imode == 1) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=dhij(ja1,ja2,jsd1,jsv2)
!              write(*,*)ja1,ja2,jsd1,jsv2,ddd1,ddd2
             endif  
             if (imode == 2) then
               ddd1=dpij(ja1,ja2,jsd1,jsv2)
               ddd2=dsij(ja1,ja2,jsd1,jsv2)
             endif
             if (imode == 3) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=ham_tb0(ja1,ja2,jsd1,jsv2)
             endif  
             ddsum=ddsum+ddd1*ddd2
           enddo  
         enddo  
      enddo  
    ddsum=ddsum*2.0d0
    value_of_etb=ddsum
!
  end subroutine qm_calc_etb_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine qm_calc_ecsc_atom(atm_index, value_of_ecsc, value_of_ecsc_onsite)
    
    use M_qm_domain, only: gamma_csc, noav, delta_e_num !(unchanged)
    use M_config, only : config !(unchanged)

    implicit none
    integer,  intent(in)   :: atm_index
    real(8),  intent(out)  :: value_of_ecsc, value_of_ecsc_onsite
    integer :: jsv1,jsv2, imode


    value_of_ecsc        = 0.0d0
    value_of_ecsc_onsite = 0.0d0
!
    imode=1
!   if (config%calc%genoOption%CSC_method == 'ELSTNER') imode=0
    if (config%calc%genoOption%CSC_max_loop_count <= 0) imode=0
    if (imode == 0) return
!
    jsv1=atm_index
    do jsv2=1,noav
      value_of_ecsc = value_of_ecsc + 0.5*gamma_csc(jsv1,jsv2)*delta_e_num(jsv1)*delta_e_num(jsv2)
    end do
!
    jsv2=jsv1
    value_of_ecsc_onsite=0.5*gamma_csc(jsv1,jsv2)*delta_e_num(jsv1)*delta_e_num(jsv2)
!
  end subroutine qm_calc_ecsc_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine qm_calc_ecsc_list(ecsc_atom_list)
!    
    use M_qm_domain, only: noav

    implicit none
    real(8), intent(inout)   :: ecsc_atom_list(:,:)
    integer :: jsv1
    real(8) :: value_of_ecsc_atm, value_of_ecsc_atm_onsite
!
    if (size(ecsc_atom_list,1) /= noav) then
      write(*,*)'ERROR(qm_calc_ecsc_list):size1=', size(ecsc_atom_list,1), noav
      stop
    endif
!
    if (size(ecsc_atom_list,2) /= 2) then
      write(*,*)'ERROR(qm_calc_ecsc_list):size1=', size(ecsc_atom_list,2), noav
      stop
    endif
!
    do jsv1=1,noav
      call qm_calc_ecsc_atom(jsv1, value_of_ecsc_atm, value_of_ecsc_atm_onsite)
      ecsc_atom_list(jsv1,1)=value_of_ecsc_atm
      ecsc_atom_list(jsv1,2)=value_of_ecsc_atm_onsite
    enddo
!    
  end subroutine qm_calc_ecsc_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_geno_output

