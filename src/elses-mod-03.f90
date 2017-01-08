! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  @@ Select system (for Hamiltonian)
!c
module elses_mod_sel_sys
  character(len=32) c_system
  real(8) r_cut_book
end module
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
module elses_mod_orb1
  integer  nvl
  integer, target, allocatable :: nval(:)
end module
!
module elses_mod_orb2
  integer  n_tot_base
  integer, allocatable :: j2js(:),j2ja(:),js2j(:,:)
  real*8,  allocatable :: dbx(:), dby(:), dbz(:), idngl(:)
end module
!
module elses_mod_r_base
   real*8,  allocatable :: r_base(:,:)
!        ---> ratio of basis components : r_base(isym,jj)
!             isym = 1,2,3,4,5 : |s>, |px>, |py>, |pz>, |s*>
!               jj = 1, n_tot_base
end module
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
module elses_arr_dbij
  real*8, allocatable :: dbij(:,:,:,:)
end module
!
module elses_arr_dhij
  real*8, allocatable :: dhij(:,:,:,:)
end module
!
module elses_arr_dfij
  real*8, allocatable :: dfij(:,:,:,:,:)
end module
!c
module elses_arr_dbij2
  real*8, allocatable :: dbij2(:,:,:,:)
!c          ----> used only as work array in elses_sym_dbij
end module
!c
!cccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
module arr_kry_ldos
  integer ite_ldos_calc
  integer nmesh
  real*8  etop,ebtm,dev,de
  real*8,  allocatable :: dldos(:,:)
!c        -----> LDOS for each atom :D(E,na)
!c                NOTE:  na=0 for the total DOS
!c
end module
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Module for Hamiltonian with partial interaction
!c        (so as to calculate COHP )
!
module elses_arr_dhij_cohp
  real*8, allocatable :: dhij_cohp(:,:,:,:)
end module
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of DBJI, DBIJ2, DHIJ
!c
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_alloc_mat1(nvl_def,noas_def,noav_def,imode)
  use M_config,        only : config ! CHANGED : config%calc%limit%memory_allocated
  use elses_arr_dbij,  only : dbij
  use elses_arr_dbij2, only : dbij2
  use elses_arr_dhij,  only : dhij
!c
  implicit none
  integer nvl_def, noas_def, noao_def, noav_def
  integer imode, ierr, nddd1
  real*8 dsum,ddd1
  logical :: alloc_dbij, alloc_dbij2 
  integer :: lu
!c
  alloc_dbij  = .true.
  alloc_dbij2 = .true.
!
  lu=config%calc%distributed%log_unit
!
  if (config%calc%mode == 'matrix_generation') then
    alloc_dbij  = .false.
    alloc_dbij2 = .false.
  endif
!
  if (lu >0) then
    write(lu,*)'@@ ALLOC_MAT1:matrix allocation :imode=',imode
    write(lu,*)'  noav_def=',noav_def
    write(lu,*)'  noas_def=',noas_def
    write(lu,*)'   nvl_def=',nvl_def
  endif
!c
  ierr=0
  if ((noav_def .le. 0) .or. (noav_def .gt. 3000000)) ierr=1
  if ((nvl_def .le. 0) .or. (nvl_def .gt. 9)) ierr=1
  if ((noas_def .le. 0) .or. (noas_def .gt. noav_def)) ierr=1
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
  dsum=0.0d0
!c       ---> total size (GB)
!c
  if (alloc_dbij) then
    ddd1=8.d0*dble(nvl_def)*dble(nvl_def)*dble(noas_def) &
&          *dble(noav_def)/1.0d9
    dsum=dsum+ddd1
    if (lu >0) write(lu,*)'... DBIJ      : size (GB)=',ddd1
  endif
!
  if (alloc_dbij2) then
    ddd1=8.d0*dble(nvl_def)*dble(nvl_def)*dble(noas_def) &
&          *dble(noav_def)/1.0d9
    dsum=dsum+ddd1
    if (lu >0) write(lu,*)'... DBIJ2     : size (GB)=',ddd1
  endif
!
  ddd1=8.d0*dble(nvl_def)*dble(nvl_def)*dble(noas_def) &
&          *dble(noav_def)/1.0d9
  dsum=dsum+ddd1
  if (lu >0) write(lu,*)'... DHIJ      : size (GB)=',ddd1
!c
  if (lu >0) write(lu,*)'... total     : size (GB)=',dsum
!c
   config%calc%limit%memory_allocated=config%calc%limit%memory_allocated+dsum
   if (lu > 0) write(lu,'(a,f20.8)')'INFO-MEMORY:estimated size [GB]=', &
&                               config%calc%limit%memory_allocated
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Check the memory size (optional)
!
   if ( config%calc%limit%memory > 0.0d0 ) then
     if (config%calc%limit%memory_allocated > config%calc%limit%memory) then
       write(*,'(a,f20.10)')'ERROR!:(allocated memory size) [GB] =', config%calc%limit%memory_allocated
       write(*,'(a,f20.10)')'           (memory size limit) [GB] =', config%calc%limit%memory
       stop
     endif
   endif
!
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of dbij
!c
  if (alloc_dbij) then
    if (lu >0) write(lu,*)'Allocation of dbij'
    allocate (dbij(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
    if( ierr .ne. 0) then
      write(*,*)'allocation error!(MAT_ALLOC;dbij):ierr=',ierr
      stop
    endif
  endif  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of dbij2
!c
  if (alloc_dbij2) then
    if (lu >0) write(lu,*)'Allocation of dbij2 and (only for SYM_DBIJ)'
    allocate (dbij2(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
    if( ierr .ne. 0) then
      write(*,*)'allocation error!(MAT_ALLOCdbij2):ierr=',ierr
      stop
    endif
  endif  
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of dhij
!c
  if (lu >0) write(lu,*)'Allocation of dhij'
  allocate (dhij(nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(MAT_ALLOC;dhij):ierr=',ierr
    stop
  endif
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     stop
!c
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Allocation of nval
!c
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_alloc_nval(nos_def,imode)
  use elses_mod_orb1, only : nval
  implicit none
  integer nos_def, imode, ierr
!c
  write(*,*)'@@ LSES_ALLOC_NVAL:nos,imode=',nos_def,imode
!c
  allocate (nval(nos_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(NVAL):ierr=',ierr
    stop
  endif
  nval(:)=-1
!c      ---> initial dummy value
!c
end
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  @@ Allocation of j2js,j2ja,js2j,dbx,dby,dbz,idngl
!c
     !! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_alloc_mat2(noa_def,n_tot_base_def,nvl_def,imode)
  use elses_mod_orb2, only : j2js,j2ja,js2j,n_tot_base,dbx,dby,dbz,idngl
!c
  implicit none
  integer noa_def, n_tot_base_def,nvl_def,imode,ierr
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Checking parameters
!c
  if ((nvl_def .le. 0) .or. (nvl_def .gt. 1000)) then
    write(*,*)'ERROR(LSES_ALLOC_MAT2):nvl_def=',nvl_def
    stop
  endif   
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
  allocate (j2js(n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(J2JS):ierr=',ierr
    stop
  endif
!c
  allocate (j2ja(n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(J2JA):ierr=',ierr
    stop
  endif
!c
  allocate (js2j(nvl_def,n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(JS2J):ierr=',ierr
    stop
  endif
!c
  allocate (dbx(n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(DBX):ierr=',ierr
    stop
  endif
!c
  allocate (dby(n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(DBY):ierr=',ierr
    stop
  endif
!c
  allocate (dbz(n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(DBZ):ierr=',ierr
    stop
  endif
!c
  allocate (idngl(n_tot_base_def),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'allocation error!(IDNGL):ierr=',ierr
    stop
  endif
!c
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  @@ Allocation of dfij
!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_alloc_dfij(nvl_def, noas_def, noav_def,imode)
  use M_config,        only : config !(unchanged)
  use elses_arr_dfij, only : dfij
  implicit none
  integer nvl_def, noas_def, noav_def, imode, ierr
  integer nddd1
  real*8  ddd1
  logical :: alloc_dfij
  integer :: lu
!
  alloc_dfij = .true. 
!
  lu=config%calc%distributed%log_unit
!
  if (config%calc%mode == 'matrix_generation') then
    alloc_dfij  = .false.
  endif
!
  if ( .not. alloc_dfij ) return
!
  if (imode .eq. 1) then
    if (lu >0) write(lu,*)'Allocation of dfij'
    allocate (dfij(3,nvl_def,nvl_def,noas_def,noav_def),stat=ierr)
    if( ierr .ne. 0) then
      write(*,*)'Alloc. error!(DFIJ):ierr=',ierr
      stop
    endif
    ddd1=8.d0*3.0d0*dble(nvl_def)*dble(nvl_def)*dble(noas_def) &
&                  *dble(noav_def)/1.0d9
    if (lu >0) write(lu,*)'... DFIJ      : size (GB)=',ddd1
  endif
!c  
  if (imode .eq. 2) then
     if (lu >0) write(lu,*)'Deallocation of dfij'
     deallocate (dfij,stat=ierr)
     if( ierr .ne. 0) then
       write(*,*)'Dealloc. error!(DFIJ):ierr=',ierr
       stop
     endif
  endif   
!c
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  @@ Allocation of r_base
!c
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_alloc_r_base
  use M_config,         only : config !(unchanged)
  use elses_mod_orb1,   only : nvl
  use elses_mod_orb2,   only : n_tot_base
  use elses_mod_r_base, only : r_base
!c
  integer ierr, lu
!c
  lu=config%calc%distributed%log_unit
!
  if (lu >0) then
    write(*,*)'@@ Allocation of r_base'
    write(*,*)'  nvl, n_tot_bae=',nvl,n_tot_base
  endif
!c
  allocate (r_base(nvl,n_tot_base),stat=ierr)
  if( ierr .ne. 0) then
    write(*,*)'Alloc. error!(R_BASE):ierr=',ierr
    stop
  endif
  r_base(:,:)=0.0d0
!c
end subroutine elses_alloc_r_base


      
