!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
czzz  @@@@ elses-mod4.f @@@@@
czzz  @@@@@ 2009/06/27 @@@@@
cccc2006cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0727: Prepared. (NT116p33) based on zmodule.f on pt40/work04
cccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0216T.Hoshi:(NT07E-122p69) Modi. on WRITE command
c        in subroutine elses_alloc_kry_rKnstr
!ccc2009cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0627T.Hoshi; Modi. on  'elses_alloc_ini_txp' (plus-20090627)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_param_ctl_kr
        integer noav_kr_def
        integer noak_min_def
        integer nreclc_def
        integer iKnstr
c
        integer nrec_kr_def
c           ---> used only for allocation command
        integer nval_kr_def
c           ---> used only for allocation command
c
        integer nkene_def
        real*8  ddemin_def
        real*8  ddemax_def
        real*8  rNelec_kr_def
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_arr_kry_glo
c
        integer nrecmx
c
        real*8, allocatable :: rA(:,:,:)
        real*8, allocatable :: rB(:,:,:)
        real*8, allocatable :: rE(:,:,:)
        real*8, allocatable :: rR(:,:,:)
        real*8, allocatable :: rW(:,:,:)
        real*8, allocatable :: rxmu(:,:)
        real*8, allocatable :: rRho(:,:,:)
        real*8, allocatable :: rFocc(:,:,:)
        real*8, allocatable :: rRNave(:,:)
        real*8, allocatable :: rRNdia(:,:)
        integer,allocatable :: nrecg(:,:)
c
        real*8, allocatable :: rcut_kry(:)
        integer,allocatable :: noak_str(:)
        integer,allocatable :: noak_min(:)
c
        real*8, allocatable :: rKnstr(:,:,:,:)
c           vector : < i | K_n^(j) > in KR space 
c
        integer,allocatable :: jsv4jsk_str(:,:)
c           used as jsv1 = jsv4jsk_str(jsk,jsv2)
c
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@@ Initial allocation or final deallocation 
c            for the Krylov subspace
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_kry_glo(imode)
      use elses_param_ctl_kr, only : noav_kr_def, nrec_kr_def,
     +              nval_kr_def
      use elses_arr_kry_glo, only : rA,rB,rE,rR,rW,rxmu,rRho,rFocc,
     +               rRNave,rRNdia,nrecg,rcut_kry,noak_str,noak_min
c
c
      implicit none
      integer noav_def,nrec_def,nvl_def
      integer imode,ierr
c
c
      noav_def=noav_kr_def
      nvl_def=nval_kr_def
      nrec_def=nrec_kr_def
      write(6,*)'@@ ALLOC_KRY2:IMODE=',imode
      write(6,*)'           noav_def=',noav_def
      write(6,*)'            nvl_def=',nvl_def
      write(6,*)'           nrec_def=',nrec_def
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  @@ Parameter checking
c
      ierr=0
      if (noav_def .le. 0) ierr=1
      if (noav_def .gt. 100000000) ierr=1
      if (nrec_def .le. 0) ierr=1
      if (nrec_def .gt. 100000) ierr=1
      if (nvl_def .le. 0) ierr=1
      if (nvl_def .gt. 20) ierr=1
      if (ierr .ne. 0) then
        write(6,*)'ERROR!(Parameter inconsistency)'
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (imode .ne. 1) goto 1999
c
       allocate (rA(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rA:ierr=',ierr
         stop
       endif
       rA(:,:,:)=0
c
       allocate (rB(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rB:ierr=',ierr
         stop
       endif
       rB(:,:,:)=0
c
       allocate (rE(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rE:ierr=',ierr
         stop
       endif
       rE(:,:,:)=0
c
       allocate (rR(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)eR:ierr=',ierr
         stop
       endif
       rR(:,:,:)=0
c
       allocate (rW(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rW:ierr=',ierr
         stop
       endif
       rW(:,:,:)=0
c
       allocate (nrecg(nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)nrecg:ierr=',ierr
         stop
       endif
       nrecg(:,:)=0
c
       allocate (rxmu(nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rxmu:ierr=',ierr
         stop
       endif
       rxmu(:,:)=-10000.0d0
c
       allocate (rRho(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rRho:ierr=',ierr
         stop
       endif
       rRho(:,:,:)=0.0d0
c
       allocate (rRNave(nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rRNave:ierr=',ierr
         stop
       endif
       rRNave(:,:)=0.0d0
c
       allocate (rRNdia(nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rRNdia:ierr=',ierr
         stop
       endif
       rRNdia(:,:)=0.0d0
c
       allocate (rFocc(nrec_def,nvl_def,noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)rFocc:ierr=',ierr
         stop
       endif
       rFocc(:,:,:)=0.0d0
c
       allocate(noak_str(noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(ALLOC_KRY)NOAK:ierr=',ierr
         stop
       endif
       noak_str(:)=0.0d0
c
       allocate(noak_min(noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'alloc. error!(ALLOC_KRY)NOAK:ierr=',ierr
         stop
       endif
       noak_min(:)=0.0d0
c
       allocate (rcut_kry(noav_def),stat=ierr)
       if(ierr .ne. 0) then
         write(6,*)'allocation error!(ALLOC_KRY)r_cut_kry:ierr=',ierr
         stop
       endif
       rcut_kry(:)=0.0d0
c
c
c     stop
c
 1999 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (imode .ne. 2) goto 2999
c
       deallocate(rA)
       deallocate(rB)
       deallocate(rE)
       deallocate(rR)
       deallocate(rW)
       deallocate(nrecg)
       deallocate(rxmu)
       deallocate(rRho)
       deallocate(rRNave)
       deallocate(rRNdia)
       deallocate(rFocc)
       deallocate(noak_str)
       deallocate(noak_min)
       deallocate(rcut_kry)
c
 2999 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 9999 continue
c
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@@ Initial allocation or final deallocation 
c            for the Krylov subspace
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_kry_rKnstr(imode)
      use elses_mod_ctrl,     only : i_verbose
      use elses_param_ctl_kr, only : noav_kr_def, nrec_kr_def,
     +                              nval_kr_def
      use elses_arr_kry_glo, only : rKnstr 
c
      implicit none
      integer imode,ierr
      integer njjt_def,nrec_def,nvl_def,noav_def
      real(8) :: ddd
c
      njjt_def=50*nval_kr_def
      nrec_def=nrec_kr_def
      nvl_def=nval_kr_def
      noav_def=noav_kr_def
c
      ddd=dble(njjt_def)*dble(nrec_def)*dble(nvl_def)*dble(noav_def)
      if (i_verbose >= 1) then
        write(6,*)'@@ elses_alloc_kry_rKnstr'
        write(6,*)'   njjt_def =',njjt_def
        write(6,*)'   nrec_def =',nrec_def
        write(6,*)'   nvl_def  =',nrec_def
        write(6,*)'   noav_def =',noav_def
        write(6,*)'   number of elements  =',ddd
        write(6,*)'... rKnstr     : size (GB)=',8.0d0*ddd/1.0d9
      endif  
c
      allocate(rKnstr(njjt_def,nrec_def,nvl_def,noav_def),
     +                stat=ierr)
      if(ierr .ne. 0) then
        write(6,*)'alloc. error!(SETKRY)rKnstr:ierr=',ierr
        stop
      endif
      rKnstr(:,:,:,:)=0.0d0
c
      end    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c @@@ Allocation or deallocation of jsv4jsk_str
c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_jsv4jsk_str(imode)
      use elses_param_ctl_kr,only : noav_kr_def
      use elses_arr_kry_glo, only : jsv4jsk_str
      use elses_arr_kry_glo, only : noak_str
c
      implicit none
      integer imode, ierr
      integer noak_max
      integer jsv, noak, noav_kr

c
      if (imode .ne. 1) goto 1999
c
      noav_kr=noav_kr_def
      write(6,*)' Alloc of jsv4jsk_str: noav_kr=',noav_kr
c
      noak_max=-1
      do jsv=1,noav_kr
        noak=noak_str(jsv)
        if (noak .gt. noak_max) noak_max=noak
      enddo   
      write(6,*)' noak_max=',noak_max
c
      if (noak_max .gt. noav_kr) then
        write(6,*)' ERROR noak_max_def=',noak_max
        stop
      endif
c
      allocate(jsv4jsk_str(noak_max, noav_kr),stat=ierr)
      if(ierr .ne. 0) then
        write(6,*)'alloc. error!(JSV4JSK):ierr=',ierr
        stop
      endif
      jsv4jsk_str(:,:)=-1
c
      write(6,*)'  ..alloc. is done successfully'
      goto 9999
c
 1999 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
      if (imode .eq. 2) then
        write(6,*)' Dealloc of jsv4jsk_str'
        deallocate(jsv4jsk_str,stat=ierr)
        if(ierr .ne. 0) then
          write(6,*)'dealloc. error!(JSV4JSK):ierr=',ierr
          stop
        endif
        write(6,*)'  ..dealloc. is done successfully'
      endif  
c
 9999 continue

      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
