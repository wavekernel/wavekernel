!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Sepecify the visualization tool
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_xml_vis
      use M_config
      implicit none
      integer :: n_interval2
      character(len=64) flnm
!
      n_interval2=config%output%position%interval
      flnm=trim(config%output%position%filename)
      call elses_vis_init(n_interval2,flnm)
!
      end subroutine elses_xml_vis
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!@@ Read the temperature
!
!      !! Copyright (C) ELSES. 2007-2016 all rights reserved
!     subroutine elses_xml_load_temperature
!     use M_config
!     use elses_mod_phys_const, only : ev2kel, ev4au
!     use elses_mod_thermo,   only : tempk0
!     implicit none 
!     real(8) :: tempkel
!
!     write(*,*)'@@ LSES_LOAD_TEMPERATURE'
!
!
!      tempkel=600.0d0
!        ---> temperature in kelvin
!
!      write(*,*)' tempkel=',tempkel
!
!      tempk0 = config%system%temperature
!
!      tempk0=tempkel/ev2kel/ev4au
!      ---> temperature in au
!     write(*,*)' tempk0=',tempk0
!     stop
!
!     end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_xml_load_ctrl
      use M_config
!
      implicit none
!     integer :: i, j, js
!     type(atom_type), pointer :: atom
      real(8) :: r_mem_limit
!     real*8  :: r_mem_limit, dtmd_tmp, tot_sim_time
      integer :: ntlimit
!     integer :: ntlimit, noak_min_tmp, nrecl_tmp
!
      ntlimit = config%calc%limit%time ! sec
      r_mem_limit = config%calc%limit%memory  ! byte
!
      write(*,*)'ntlimit     =',ntlimit
      write(*,*)'r_mem_limit =',r_mem_limit
!
      end subroutine elses_xml_load_ctrl
!      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_xml_save_txp
      use M_config
!     use elses_mod_sel_sys, only : c_system, r_cut_book
!     use elses_mod_phys_const, only : angst, ev2kel, ev4au, au_fsec
!      use elses_mod_sim_cell, only : noa
!     use elses_mod_sim_cell, only : noa,nos,iperiodic,ax,ay,az
!     use elses_mod_noav,     only : noav,noao,noab,noas,noab0,nncut
!     use elses_mod_orb1,     only : nvl
      use elses_mod_tx,       only : tx, ty, tz, jsei
      use elses_mod_txp,      only : txp, typ, tzp
!     use elses_mod_txp0,     only : txp0, typ0, tzp0
!     use elses_mod_foi,      only : foi
      use elses_mod_md_dat,   only : dtmd, itemdmx
      use elses_mod_md_dat,   only : itemd,dtmd, itemdmx
!     use elses_mod_mass,     only : awt
!     use elses_mod_thermo, 
!    +          only : tempk0,amq,thb,thbold,vhb,vhbold
!
!     use elses_mod_sim_cell, only : noa,iperiodic,ax,ay,az
      use elses_mod_sim_cell, only : noa,ax,ay,az
      use elses_mod_vel,      only : velx, vely, velz
      use elses_mod_foi,      only : foi
      use elses_mod_md_dat,   only : itemdorg
      use elses_mod_thermo, 
     +          only : amq, thb,thbold,vhb,vhbold
      use elses_mod_iflag, only : iflag

      implicit none
      integer :: j
      type(atom_type), pointer :: atom
!     integer :: ntlimit, noak_min_tmp, nrecl_tmp
      real(8) :: r_hb_mass_per_atom
!    
!     config%system%structure%mdstep = itemdorg+itemd-1

      r_hb_mass_per_atom=1.0d0/amq/dble(noa)

      config%system%structure%heatbath%massperatom 
     +         = r_hb_mass_per_atom
      config%system%structure%heatbath%position = thb
      config%system%structure%heatbath%velocity = vhb
!
      config%system%structure%unitcell%vectorA(1)=ax
      config%system%structure%unitcell%vectorB(2)=ay
      config%system%structure%unitcell%vectorC(3)=az
!
      write(*,*) '@@ XML_SAVE_TXP:NOA=',noa
!
!     if ( .not. allocated(config%system%structure%vatom) ) then
!       write(*,*)'ERROR(elses_xml_save_txp):Not allocated:vatom'
!       stop
!     endif
!
      write(*,*)'natom=',config%system%structure%natom
!
      do j=1, config%system%structure%natom
!        write(*,*)'j=',j
!        write(*,*)'txp=',txp(j),typ(j),tzp(j)
         atom => config%system%structure%vatom(j)

         atom%position(1) = txp(j)*ax
         atom%position(2) = typ(j)*ay
         atom%position(3) = tzp(j)*az

         atom%motion = "free"
         if (allocated(iflag)) then
           if( iflag(j) == 1 ) then
              atom%motion = "free"
           else if( iflag(j) == 0 ) then
              atom%motion = "fix"
           end if
         endif  
!         
         if ( allocated(velx) ) then
           atom%velocity(1) = velx(j)*(ax/dtmd)
         else  
           atom%velocity(1) = 0.0d0
         endif   
!
         if ( allocated(vely) ) then
           atom%velocity(2) = vely(j)*(ay/dtmd)
         else  
           atom%velocity(2) = 0.0d0
         endif   
!
         if ( allocated(velz) ) then
           atom%velocity(3) = velz(j)*(az/dtmd)
         else  
           atom%velocity(3) = 0.0d0
         endif   
!
         if (allocated(foi)) then
           atom%force(:) = foi(j,:)
         endif  
!
      end do
!
      write(*,*)'...routine(elses_xml_save_txp) is ended. '
!
      end subroutine elses_xml_save_txp
!      
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@ Allocate and set various parameters for compatibility
c         input : noa, nos, nvl
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_init_compat2
      use elses_mod_sim_cell, only : noa, nos
      use elses_mod_tx,       only : jsei
      use elses_mod_orb1,     only : nvl,nval
      use elses_mod_orb2,     only : j2js,j2ja,js2j,n_tot_base,
     +                                dbx,dby,dbz,idngl
      use elses_mod_multi,    only : ict4h,ict4l,intf
      use elses_mod_r_base,   only : r_base
c
      implicit none
!     integer nos_def,noa_def
      integer noa_def
      integer js,ja,iii,nss,nval2,j,n_tot_base_def,nvl_def,imode
!     integer jb
      real(8) :: dbx1,dby1,dbz1,ddd,dd, dbs1, dbs2
c
      write(6,*)'@@ elses_init_compat2:NOA,NOS,NVL=',noa,nos,nvl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Error checking
c
      if (noa*nos*nvl .le. 0) then
         write(6,*)'ERROR!(elses_init_compat)'
         stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Allocate : NVAL
c
!     nos_def=nos
!     imode=1
!     call elses_alloc_nval(nos_def,imode)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Setting NVAL,JSEI for one-component system
c
!     nval(1)=nvl
!     do js=1,noa
!       jsei(js)=1
!     enddo  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Setting ICT4H, ICT4L, INTF(NOA) 
c         for compatibility to multi-solver scheme
c
      ict4h=1
      ict4l=1
      noa_def=noa
      call elses_alloc_intf(noa_def,imode)
c         --> alloc: intf(noa)
      intf(:)=0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Count up the total number of bases 
c
      iii=0
      do js=1,noa
        nss=jsei(js)
        nval2=nval(nss)
        iii=iii+nval2
      enddo   
      n_tot_base=iii
      write(6,*)' total number of bases=',n_tot_base
c
      imode=1
      nvl_def=nvl
      n_tot_base_def=n_tot_base
      noa_def=noa
      call elses_alloc_mat2(noa_def,n_tot_base_def,nvl_def,imode)
c         --> alloc: j2js,j2ja,js2j,dbx,dby,dbz,idngl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      j=0
      dbx(:)=0.0d0
      dby(:)=0.0d0
      dbz(:)=0.0d0
      do js=1,noa
!       jsei(js)=1
       do ja=1,nvl
        j=j+1 
        j2js(j)=js
        j2ja(j)=ja
        js2j(ja,js)=j
        idngl(j)=0
        if (ja .eq. 1) then
          dbx1=1.0d0 
          dby1=1.0d0 
          dbz1=1.0d0 
        endif   
        if (ja .eq. 2) then
          dbx1=-1.0d0 
          dby1=-1.0d0 
          dbz1=1.0d0 
        endif   
        if (ja .eq. 3) then
          dbx1=1.0d0 
          dby1=-1.0d0 
          dbz1=-1.0d0 
        endif   
        if (ja .eq. 4) then
          dbx1=-1.0d0 
          dby1=1.0d0 
          dbz1=-1.0d0 
        endif   
c
        ddd=dbx1*dbx1+dby1*dby1+dbz1*dbz1
        dd=dsqrt(ddd)
        dbx(j)=dbx1/dd
        dby(j)=dby1/dd
        dbz(j)=dbz1/dd
c
       enddo
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   @@ R_BASE is prepared.
c
      write(*,*)'R_BASE is prepared.'
      call elses_alloc_r_base
c         ----> allocate r_base(nvl,n_tot_base)
c
      j=0
      do js=1,noa
       do ja=1,nvl
        j=j+1 
c
        if (ja .le. 4) then
          dbs1=0.5d0
          dbx1=0.5d0*dsqrt(3.0d0)*dbx(j)
          dby1=0.5d0*dsqrt(3.0d0)*dby(j)
          dbz1=0.5d0*dsqrt(3.0d0)*dbz(j)
          dbs2=0.0d0
        endif   
c
        if (ja .eq. 5) then
          dbs1=0.0d0
          dbx1=0.0d0
          dby1=0.0d0
          dbz1=0.0d0
          dbs2=1.0d0
        endif   
c
        ddd=dbs1*dbs1+dbx1*dbx1+dby1*dby1+dbz1*dbz1+dbs2*dbs2
        dd=dsqrt(ddd)
c
        r_base(1,j)=dbs1/dd
        r_base(2,j)=dbx1/dd
        r_base(3,j)=dby1/dd
        r_base(4,j)=dbz1/dd
        if (nvl .eq. 5) r_base(5,j)=dbs2/dd
c
       enddo
      enddo 
c
c     do j=1,n_tot_base
c       do jb=1,nvl
c         write(6,*)'jb,j, r_base(jb,j)=',jb,j,r_base(jb,j)
c       enddo
c     enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(*,*)'....ended:LSES_INIT_COMPAT'
c     stop
c
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Mass setting : from awt(nos) to amm(noa)
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_set_mass2
      use elses_mod_phys_const, only : au_mass
      use elses_mod_tx,         only : jsei
      use elses_mod_mass,       only : amm
      use elses_mod_sim_cell,   only : noa, nos
      use elses_mod_mass,       only : awt
      implicit none
      integer noa_def, js, nss
      real(8) :: awt0
c
      noa_def=noa
      write(*,*)'@@ LSES_SET_MASS2'
c
c     call elses_alloc_amm(noa_def)
c
      if (allocated(amm) .eqv. .false.) then
         write(*,*)'ERROR:LSES_SET_MASS2'
         write(*,*)' AMM is not allocated'
         stop
      endif   
c
      do js=1,noa
        nss=jsei(js)
        awt0=awt(nss)
!
        if ((nss .le. 0) .or. (nss .gt. nos)) then
          write(*,*)'ERROR?:LSES_SET_MASS2:js,nss=',js,nss
          stop
        endif
        if (awt0 .le. 1.0d-10) then
          write(*,*)'ERROR?:LSES_SET_MASS2:awt0=',awt0
          write(*,*)'                       nss=',nss
          stop
        endif
!
        amm(js)=1.0d0/(au_mass*awt0)
        if (js .le. 5) then
          write(*,*)'js,amm=',js,amm(js)
        endif  
      enddo
c
      write(*,*)'....ended:LSES_SET_MASS'
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

