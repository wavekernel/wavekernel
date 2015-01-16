!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
czzz  @@@@ elses-mod-nrl.f @@@@@
czzz  @@@@@ 2007/03/26 @@@@@
cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0736: Prepared. (NT07A-118,p.39): Modules for NRL-related variables
c    :   'param_ham_nrl' : not essential, 
c                          --> should be deleted near future.
c    : 'elses_arr_dpij', 'elses_arr_dsij' prepared. 
c        ---> allocated in elses_alloc_dsij_dpij
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module param_ham_nrl
        integer i_elem_nrl
        real*8 rcut_ham_nrl
      end module
c
      module elses_arr_dpij
        real*8, allocatable :: dpij(:,:,:,:)
      end module
c
      module elses_arr_dsij
        real*8, allocatable :: dsij(:,:,:,:)
      end module
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of DPIJ, DSIJ
c
      !! Copyright (C) ELSES. 2007-2015 all rights reserved
      subroutine elses_alloc_dsij_dpij
     +           (nvl_def,noas_def,noav_def,imode)
      use elses_arr_dsij,  only : dsij
      use elses_arr_dpij,  only : dpij
c
      implicit none
      integer nvl_def, noas_def, noao_def, noav_def
      integer imode, ierr, nddd1
      real*8 dsum,ddd1
c
      write(6,*)'@@ ALLOC_MAT_DSIJ_DPIJ'
      write(6,*)'   :matrix allocation :imode=',imode
      write(6,*)'  noav_def=',noav_def
      write(6,*)'  noas_def=',noas_def
      write(6,*)'   nvl_def=',nvl_def
c
      ierr=0
c     if ((noav_def .le. 0) .or. (noav_def .gt. 1000000)) ierr=1
      if (noav_def .le. 0) ierr=1
      if ((nvl_def .le. 0) .or. (nvl_def .gt. 9)) ierr=1
      if ((noas_def .le. 0) .or. (noas_def .gt. noav_def)) ierr=1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dsum=0.0d0
c       ---> total size (GB)
c
      nddd1=8*nvl_def*nvl_def*noas_def*noav_def
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
      write(6,*)'... DSIJ      : size (GB)=',ddd1
c
      nddd1=8*nvl_def*nvl_def*noas_def*noav_def
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
      write(6,*)'... DPIJ     : size (GB)=',ddd1
c
      write(6,*)'... total     : size (GB)=',dsum
c
      if (ierr .ne. 0) then
        write(6,*)'ERROR?(ALLOC_MOD2):matrix size?'
        stop
      endif   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of dsij
c
      write(6,*)'Allocation of dsij'
      if (allocated(dsij) .eqv. .true.) then
        write(6,*)'ERROR:dsij is already allcoated'
        stop
      endif   
      allocate (dsij(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC)dsij:ierr=',ierr
        stop
      endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of dpij
c
      write(6,*)'Allocation of dpij'
      if (allocated(dpij) .eqv. .true.) then
        write(6,*)'ERROR:dpij is already allcoated'
        stop
      endif   
      allocate (dpij(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC):dpij:ierr=',ierr
        stop
      endif
c     stop
c
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
