!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!czzz  @@@@ elses-ana-cohp.f90 @@@@@
!czzz  @@@@@ 2008/01/02 @@@@@
!ccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0102T.Hoshi: Prepared (NT07E-122,p29)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Bond list analysis with ICOHP
!   imode=0 : data output for bond analysis
!              into the extended molden format
!                (d_cri is ignored)
!   imode=1 : bond analysis with ICOHP
!              d_cri : critical energy in eV
!              into the PDB format
!   imode=2 : bond analysis with atom distance
!              d_cri : critical length in Angst
!              into the PDB format
!
!   If n_bond = -1, no bond list is generated.
!              and only n_bond is calculated as # bonds
!   If n_bond /= -1 nor 0, 
!           the equation of (n_bond = n_bond_count) is checked
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_ana_icohp(imode,d_cri,fname,n_bond)
   use elses_mod_ctrl,     only : i_verbose
   use elses_mod_phys_const, only : angst, ev4au
   use elses_mod_md_dat,     only : itemd
   use elses_mod_tx,       only : jsei
   use elses_mod_js4jsv,   only : js4jsv 
   use elses_mod_jsv4jsd,  only : jsv4jsd, njsd
   use elses_mod_orb1,     only : nval
   use elses_arr_dbij,     only : dbij
   use elses_arr_dhij,     only : dhij
   use elses_mod_multi,    only : ict4h, ict4l
   use elses_mod_orb2,     only : js2j
   use elses_mod_file_io,  only : vacant_unit
   use elses_arr_dhij_cohp, only : dhij_cohp
   use M_md_motion,     only : get_atom_dist
!
   use elses_param_ctl_kr, only : noav_kr_def
!
   implicit none
   integer :: imode
   real(8) ::  etbtmp
!  real(8) ::  get_atom_dist
   real(8) ::  d_cri
!
   integer :: noav_kr
   integer :: jsv2,      js2, ja2, njsd2, nss, nval2, jg
   integer :: jsv1,jsd1, js1, ja1,             nval1, ig
   integer :: i_partial
!
   real(8) :: dsum, dlijd, ahijd
   real(8) :: d_icohp, drr
   real(8) :: ahijd2, d_icohp2
!
   integer :: iunit
   integer :: n_bond
   integer :: n_bond_count
   character(len=*), intent(in) :: fname
!
   noav_kr=noav_kr_def
!
   iunit=vacant_unit()
   if (i_verbose >= 1) then 
     write(*,*)'@@ elses_ana_icohp'
     write(*,*)' imode, d_cri=',imode,d_cri
     write(*,*)' fname=',trim(fname)
     write(*,*)' iunit=',iunit
     write(*,*)' n_bond=',n_bond
   endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Set dhij_cohp 
!        as Hamiltonian with partial interaction
!
   call elses_alloc_dhij_cohp(1)
   i_partial=4
   call elses_set_hami_01(i_partial)
   if (i_verbose >= 1) then 
     write(*,*)' i_partial=',i_partial
   endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   open (iunit,file=fname,position='append',status='old')
!  write(iunit,*)'test'
!
   dsum=0.0D0
   n_bond_count=0
   do jsv2=1,noav_kr
      js2=js4jsv(jsv2)
!     if ((js2 .le. 0) .or. (js2 .gt. noa)) then
!       write(6,*)'error!(enetb):jsv2,js2=',jsv2,js2
!       stop
!     endif
      njsd2=njsd(jsv2,ict4h)
      nss=jsei(js2)
      nval2=nval(nss)
      do jsd1=1,njsd2
        jsv1=jsv4jsd(jsd1,jsv2)
        if ((jsv1 .le. 0) .or. (jsv1 .gt. noav_kr)) then
          write(6,*)'error!(enetb):jsv1,noav=',jsv1,noav_kr
          stop
        endif
        js1=js4jsv(jsv1)
!       if ((js1 .le. 0) .or. (js1 .gt. noa)) then
!         write(6,*)'error!(enetb):jsv1,js1=',jsv1,js1
!         stop
!       endif
        nss=jsei(js1)
        nval1=nval(nss)
        drr=get_atom_dist(js1,js2)
        d_icohp=0.0d0
        d_icohp2=0.0d0
        do ja2=1,nval2
          jg=js2j(ja2,js2)
          do ja1=1,nval1
            ig=js2j(ja1,js1)
            dlijd =dbij(ja1,ja2,jsd1,jsv2)
            ahijd =dhij(ja1,ja2,jsd1,jsv2)
            ahijd2=dhij_cohp(ja1,ja2,jsd1,jsv2)
            dsum=dsum+2.0d0*dlijd*ahijd
            d_icohp =d_icohp+2.0d0*dlijd*ahijd
            d_icohp2=d_icohp2+2.0d0*dlijd*ahijd2
          enddo
        enddo   
        drr=drr*angst
        d_icohp =d_icohp *ev4au
        d_icohp2=d_icohp2*ev4au
        if (js1 /= js2) then
           if (imode == 0) then
             if (n_bond /= -1) write(iunit,"(2I5,3F15.8)")js2,js1, drr, d_icohp, d_icohp2
             n_bond_count=n_bond_count+1
           endif  
           if (imode == 1) then
              if (d_icohp <= d_cri) then
                 write(iunit,"('CONECT',2I5)")js2,js1
                 n_bond_count=n_bond_count+1
              endif
           endif  
           if (imode == 2) then
              if (drr <= d_cri) then
                 if (n_bond /= -1) write(iunit,"('CONECT',2I5)")js2,js1
                 n_bond_count=n_bond_count+1
              endif
           endif  
           if (imode == 3) then
              if (d_icohp2 <= d_cri) then
                 if (n_bond /= -1) write(iunit,"('CONECT',2I5)")js2,js1
                 n_bond_count=n_bond_count+1
              endif
             endif  
        endif  
     enddo
   enddo
!
   close(iunit)
!
   if ((n_bond /=0) .and. (n_bond /= -1)) then
     if (n_bond /= n_bond_count) then 
       write(6,*)'ERROR(ELSES_ANA_COHP)'
       write(6,*)'n_bond, n_bond_count=',n_bond,n_bond_count
       stop
     endif   
   endif   
!
   if (i_verbose >= 1) then 
     write(6,*)'n_bond_count=',n_bond_count
     write(6,*)'ETB,N_elec  =',DSUM
!          ---> energy and electron number
   endif  
!
   n_bond=n_bond_count

   call elses_alloc_dhij_cohp(2)
!        ---> deallocation of dhij_cohp
!
end subroutine elses_ana_icohp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c @@ Allocation or deallocation of dhij_cohp
!       imode = 1 : allocation
!       imode = 2 : deallocation
!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_alloc_dhij_cohp(imode)
   use elses_mod_ctrl,      only : i_verbose
   use elses_arr_dhij_cohp, only : dhij_cohp
   use elses_mod_noav,      only : noav,noas
   use elses_mod_orb1,      only : nvl
!
   implicit none
   integer :: imode, ierr, nddd1
   real(8) ::  dsum,ddd1
!
   if (imode == 2) then
     if (i_verbose >= 1) then 
       write(6,*)'@@ DeaLLOC_DHIJ_COHP'
     endif  
     deallocate(dhij_cohp,stat=ierr)
     if( ierr .ne. 0) then
       write(6,*)'allocation error!(MAT_ALLOC_DHIJ_COHP):ierr=',ierr
       stop
     endif
     goto 9999
   endif  
!
   if (i_verbose >= 1) then 
      write(6,*)'@@ ALLOC_DHIJ_COHP'
      write(6,*)'  noav=',noav
      write(6,*)'  noas=',noas
      write(6,*)'   nvl=',nvl
   endif  
!
!
   ierr=0
   if ((noav <= 0) .or. (noav > 1000000)) ierr=1
   if ((noas <= 0) .or. (noas > noas)) ierr=1
   if ((nvl <= 0) .or. (nvl > 9)) ierr=1
!
   if (ierr .ne. 0) then
     write(6,*)'ERROR?(ALLOC_DHIJ_COHP):matrix size?'
     stop
   endif   
!
   dsum=0.0d0
!    ---> total size (GB)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Allocation of dhij_cohp
!
   if (i_verbose >= 1) then 
      write(6,*)'Allocation of dhij_cohp'
   endif
!
   allocate (dhij_cohp(nvl,nvl,noas,noav),stat=ierr)
   if( ierr .ne. 0) then
     write(6,*)'allocation error!(MAT_ALLOC_DHIJ_COHP):ierr=',ierr
     stop
   endif
   dhij_cohp(:,:,:,:)=0.0d0
!
   nddd1=8*nvl*nvl*noas*noav
   ddd1=dble(nddd1)/1.0d9
   dsum=dsum+ddd1
!
   if (i_verbose >= 1) then 
      write(6,*)'... DHIJ_COHP      : size (GB)=',ddd1
   endif
!
9999 continue
!
end subroutine elses_alloc_dhij_cohp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
