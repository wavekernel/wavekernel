!================================================================
! ELSES version 0.05
! Copyright (C) ELSES. 2007-2015 all rights reserved
!================================================================
c  @@@@@ elses-eig-leg @@@@@
c  @@@@@ 2008/01/05 @@@@@
c2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0326: Prepared (Nt07A-118,p39) 
c       --> Eigen-state solver for generalized eigen-value problem
c        Based on zeigenl2.f 
c2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c0105: A modi. on comment
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      module elses_mod_eig_leg
        integer n_base_eig_leg
c         : number of basis 
        real*8 :: fb
c         : temperature (level-broadening) parameter in eV
        real*8 :: tot_elec_eig_leg
c         : Total electron number
      end module
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module elses_arr_eig_leg
         real*8,  allocatable :: atmp(:,:)
         real*8,  allocatable :: atmp2(:,:)
!        real*8,  allocatable :: atmp3(:,:)
!        real*8,  allocatable :: atmpo(:,:)
         real*8,  allocatable :: eig2(:)
         real*8,  allocatable :: f_occ(:)
         integer, allocatable :: idngl2(:)
      end module
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Allocation of matrix and set n_base_eig_leg
c
      !! Copyright (C) ELSES. 2007-2015 all rights reserved
      subroutine elses_alloc_eig_leg(neig0)
      use elses_mod_eig_leg, only : n_base_eig_leg
      use elses_arr_eig_leg, only : atmp, atmp2, eig2, idngl2, f_occ
!     use elses_arr_eig_leg, only : atmp, atmp2, atmp3, atmpo, 
!    +                             eig2, idngl2, f_occ
      implicit none
      integer neig0, nddd1, ierr
      real*8 dsum, ddd1
c
      write(6,*)'@@ Allocation for EIG_LEG'
c
      if ((neig0 .le.0) .or. (neig0 .gt. 100000)) then
        write(6,*)'ERROR?(Allocation for EIG_LEG)'
        write(6,*)'neig0=',neig0
        stop
      endif
      n_base_eig_leg=neig0
      write(6,*)'n_base_eig_leg is set to=',n_base_eig_leg
c
c
      dsum=0.0d0
c
      allocate (atmp(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
        stop
      endif
      nddd1=8*neig0*neig0
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
c
      allocate (atmp2(neig0,neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
        stop
      endif
      nddd1=8*neig0*neig0
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
c
c     allocate (atmp3(neig0,neig0),stat=ierr)
c     if( ierr .ne. 0) then
c       write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
c       stop
c     endif
c     nddd1=8*neig0*neig0
c     ddd1=dble(nddd1)/1.0d9
c     dsum=dsum+ddd1
c
c     allocate (atmpo(neig0,neig0),stat=ierr)
c     if( ierr .ne. 0) then
c       write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
c       stop
c     endif
c     nddd1=8*neig0*neig0
c     ddd1=dble(nddd1)/1.0d9
c     dsum=dsum+ddd1
c
      allocate (idngl2(neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
        stop
      endif
      nddd1=8*neig0
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
c
      allocate (eig2(neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
        stop
      endif
      nddd1=8*neig0*neig0
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
c
      allocate (f_occ(neig0),stat=ierr)
      if( ierr .ne. 0) then
        write(6,*)'allocation error!(MAT_ALLOC):ierr=',ierr
        stop
      endif
      nddd1=8*neig0*neig0
      ddd1=dble(nddd1)/1.0d9
      dsum=dsum+ddd1
c
      write(6,*)'  total matrix size in GB =',dsum
c
      end subroutine elses_alloc_eig_leg
