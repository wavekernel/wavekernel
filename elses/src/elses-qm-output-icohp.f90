!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_output_icohp
!
    implicit none
    integer, allocatable :: called_first(:)
!
   private
!
! Public routines
   public qm_output_icohp
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Output for ICOHP
!       imode = 1 : for full Hamiltonian
!       imode = 2 : for TB0  Hamiltonian
!       imode = 3 : for (pp-pi)-interaction 
!
  subroutine qm_output_icohp(imode)
!
    use elses_mod_file_io, only : vacant_unit !(function)
    use M_lib_phys_const, only : ev4au               !(unchanged) 
    use M_qm_domain ,  only : i_verbose, dhij, dbij, ham_tb0, atm_element, nval, &
&                             jsv4jsd, njsd, noav, atm_element !(unchanged) 
    use elses_mod_js4jsv,    only : js4jsv !(unchanged)
    use M_config,            only : config !(unchanged)
    use elses_mod_md_dat,    only : itemd  !(unchanged)
    use M_qm_geno,           only : set_hamiltonian_and_overlap_geno !(routine)
!
    use M_qm_domain ,          only : c_system              !(unchanged) 
    use M_qm_output_rest_pair, only : plot_rest_energy_pair !(routine)
!
    implicit none
!   character(len=*), parameter :: filename="output_icohp.txt"
    integer      :: imode
    integer, parameter  :: ict4h=1
!   real(kind=8) :: value_of_etb
    integer      :: jsv2, nss2, nval2, ja2, jsd1
    integer      :: jsv1, nss1, nval1, ja1
    real(kind=8) :: ddsum, ddd1, ddd2
!   real(kind=8) :: dx, dy, dz
    integer      :: n_count
    integer      :: js1, js2
    integer      :: iunit, ierr
    logical      :: flag_for_first_write
    character(len=32) :: filename_wrk
    real(kind=8), allocatable :: bondmask(:)
    real(kind=8), allocatable :: ham_tb0_part(:,:,:,:)
    real(kind=8) :: ddsum2
!
!   write(*,*)'@@ qm_output_icohp(TEST):imode,step_count=',imode, config%system%structure%mdstep
!
    if ( .not. config%output%bond_list%set ) return
    if (mod(itemd-1,config%output%bond_list%interval) /= 0 ) return
    filename_wrk=trim(config%output%bond_list%filename)
!
    if (.not. allocated(called_first)) then
      flag_for_first_write = .true.
      allocate(called_first(1), stat=ierr)
      if (ierr /=0 ) stop 'Alloc. error (qm_output_icohp):called_first'
    else
      flag_for_first_write = .false.
    endif   
!
    if (i_verbose >= 1) then
      write(*,*)'@@ qm_output_icohp:imode,step_count=',imode, config%system%structure%mdstep
    endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Partial interaction Hamiltonian
!
    if (imode == 3) then

      if (i_verbose >= 1) then
        write(*,*)'INFO: pi-COHP calculation'
      endif
!
      if (size(dhij,1) < 4) then
        write(*,*)'ERROR(qm_output_icohp):size(dhij)= ', size(dhij,1)
        stop
      endif
!
      allocate(bondMask(size(dhij,1)))
      bondMask(:)=1d0
      bondMask(4)=0d0   ! The pp-pi interaction is masked
!       ! bondMask = (ss-\sigma,n.a,pp-\sigma,pp-\pi,n.a,n.a,dd-\sigma,dd-\pi,dd-\delta)
      call set_hamiltonian_and_overlap_geno(bondMask)
!       !    ----> dsij and ham_tb0 are replaced by the pi-interaction components
      allocate(ham_tb0_part(size(dhij,1), size(dhij,2), size(dhij,3), size(dhij,4)), stat=ierr)
      if (ierr /= 0) stop 'Alloc. error (qm_output_icohp):ham_tb0_pi'
      ham_tb0_part(:,:,:,:)=ham_tb0(:,:,:,:)
      call set_hamiltonian_and_overlap_geno
!       !    ----> dsij and ham_tb0 are restored as the true ones.
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iunit=vacant_unit()
!
    if ( flag_for_first_write ) then
      if (i_verbose >= 1) then
        write(*,*)' The output ICOHP file is initialized : ',trim(filename_wrk)
      endif
      open(iunit, file=filename_wrk, status='unknown')
    else
      open(iunit, file=filename_wrk, status='unknown', position='append')
    endif   
!
    write(*,*)'CHECK-INFO:sum(dhij)=',sum(dhij)
    write(*,*)'CHECK-INFO:sum(dbij)=',sum(dbij)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Count the number of the atom pairs
!  
    n_count=0
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd(jsv2,ict4h)
        jsv1=jsv4jsd(jsd1,jsv2)
        n_count=n_count+1
!       if (jsv1 /= jsv2) n_count=n_count+1
      enddo
    enddo
    if (i_verbose >= 1) then
      write(*,*)'   n_count=',n_count
    endif
!
    write(iunit,'(a)')'%%MatrixMarket matrix coordinate real general'
    write(iunit,'(a,i10)')'% step_count=', config%system%structure%mdstep
    write(iunit,*)noav, noav, n_count
!
    ddsum2=0.0d0
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      js2=js4jsv(jsv2)
      do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         js1=js4jsv(jsv1)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         ddsum=0.0d0
         do ja2=1,nval2
           do ja1=1,nval1
             if (imode == 1) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=dhij(ja1,ja2,jsd1,jsv2)
             endif  
             if (imode == 2) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=ham_tb0(ja1,ja2,jsd1,jsv2)
             endif  
             if (imode == 3) then
               ddd1=dbij(ja1,ja2,jsd1,jsv2)
               ddd2=ham_tb0(ja1,ja2,jsd1,jsv2)-ham_tb0_part(ja1,ja2,jsd1,jsv2)
!              write(*,'(a,4i10,2f20.10)')'COHP-comp=',ja1, ja2, jsv1,jsv2,ddd1,ddd2
             endif  
             ddsum=ddsum+ddd1*ddd2
           enddo  
         enddo  
         ddsum=ddsum*2.0d0*ev4au
         ddsum2=ddsum2+ddsum
         write(iunit,'(2i10,f20.10)')jsv1, jsv2, ddsum
!        if (jsv1 /= jsv2) write(iunit,'(2i10,f20.10)')jsv1, jsv2, ddsum
      enddo  
    enddo  
!
    if (allocated(bondMask))  deallocate(bondMask)
    if (allocated(ham_tb0_part)) deallocate(ham_tb0_part)
!
    close(iunit)
!
    if (i_verbose >= 1) then
      write(*,'(a,2f20.10)')'INFO:ICOHP total sum (should be ETB) [au,eV] =', ddsum2/ev4au, ddsum2
    endif
!
    if (c_system == 'geno') call plot_rest_energy_pair
!
  end subroutine qm_output_icohp
!
end module M_qm_output_icohp

