!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!czzz  @@@@ elses-mod-nrl.f @@@@@
!czzz  @@@@@ 2007/03/26 @@@@@
!cccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c0736: Prepared. (NT07A-118,p.39): Modules for NRL-related variables
!c    :   'param_ham_nrl' : not essential, 
!c                          --> should be deleted near future.
!c    : 'elses_arr_dpij', 'elses_arr_dsij' prepared. 
!c        ---> allocated in elses_alloc_dsij_dpij
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
      module param_ham_nrl
        integer i_elem_nrl
        real*8 rcut_ham_nrl
      end module
!c
      module elses_arr_dpij
        real*8, allocatable :: dpij(:,:,:,:)
      end module
!c
      module elses_arr_dsij
        real*8, allocatable :: dsij(:,:,:,:)
      end module
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of DPIJ, DSIJ
!c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
      subroutine elses_alloc_dsij_dpij(nvl_def,noas_def,noav_def,imode)
      use M_config,        only : config ! CHANGED : config%calc%limit%memory_allocated
      use elses_arr_dsij,  only : dsij
      use elses_arr_dpij,  only : dpij
!c
      implicit none
      integer nvl_def, noas_def, noao_def, noav_def
      integer imode, ierr, nddd1
      real*8 dsum,ddd1
      integer :: lu
      logical :: alloc_dpij
!c
      alloc_dpij = .true.
      lu = config%calc%distributed%log_unit
!
      if (config%calc%mode == 'matrix_generation') alloc_dpij = .false.
!
      if (lu > 0) then
        write(lu,*)'@@ ALLOC_MAT_DSIJ_DPIJ'
        write(lu,*)'   :matrix allocation :imode=',imode
        write(lu,*)'  noav_def=',noav_def
        write(lu,*)'  noas_def=',noas_def
        write(lu,*)'   nvl_def=',nvl_def
      endif  
!c
      ierr=0
!c    if ((noav_def .le. 0) .or. (noav_def .gt. 1000000)) ierr=1
      if (noav_def .le. 0) ierr=1
      if ((nvl_def .le. 0) .or. (nvl_def .gt. 9)) ierr=1
      if ((noas_def .le. 0) .or. (noas_def .gt. noav_def)) ierr=1
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
      dsum=0.0d0
!c       ---> total size (GB)
!c
      ddd1=dble(nvl_def)*dble(nvl_def)*dble(noas_def)*dble(noav_def)
      ddd1=8.0d0*ddd1/1.0d9
      dsum=dsum+ddd1
      if (lu > 0) write(lu,*)'... DSIJ      : size (GB)=',ddd1
!c
      if ( alloc_dpij ) then
        ddd1=dble(nvl_def)*dble(nvl_def)*dble(noas_def)*dble(noav_def)
        ddd1=8.0d0*ddd1/1.0d9
        dsum=dsum+ddd1
        if (lu > 0) write(lu,*)'... DPIJ     : size (GB)=',ddd1
      endif
!c
      if (lu > 0) write(lu,*)'... total     : size (GB)=',dsum
!
      config%calc%limit%memory_allocated=config%calc%limit%memory_allocated+dsum
      if (lu > 0) write(lu,'(a,f20.8)')'INFO-MEMORY:estimated size [GB]=', &
&                               config%calc%limit%memory_allocated
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Check the memory size
!
      if ( config%calc%limit%memory > 0.0d0 ) then
        if (config%calc%limit%memory_allocated > config%calc%limit%memory) then
          write(*,'(a,f20.10)')'ERROR!:(allocated memory size) [GB] =', config%calc%limit%memory_allocated
          write(*,'(a,f20.10)')'           (memory size limit) [GB] =', config%calc%limit%memory
          stop
        endif
      endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of dsij
!c
      if (lu > 0) write(lu,*)'Allocation of dsij'
      if (allocated(dsij) .eqv. .true.) then
        write(*,*)'ERROR:dsij is already allcoated'
        stop
      endif   
      allocate (dsij(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
      if( ierr .ne. 0) then
        write(*,*)'allocation error!(MAT_ALLOC)dsij:ierr=',ierr
        stop
      endif
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of dpij
!c
      if ( .not. alloc_dpij ) return
!
      if (lu > 0) write(lu,*)'Allocation of dpij'
      if (allocated(dpij) .eqv. .true.) then
        write(*,*)'ERROR:dpij is already allcoated'
        stop
      endif   
      allocate (dpij(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
      if( ierr .ne. 0) then
        write(*,*)'allocation error!(MAT_ALLOC):dpij:ierr=',ierr
        stop
      endif
!c     stop
!c
!
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
