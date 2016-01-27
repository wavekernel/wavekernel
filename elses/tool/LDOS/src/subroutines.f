!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
************************************************************************
      module common_constants
*     define knd=kind(1d0)
*            SMALL_EPSILON
************************************************************************
      implicit none
      integer,  parameter :: knd=kind(1d0)
      real(knd),parameter :: SMALL_EPSILON=1d-10
      !real(knd),parameter :: Ry2eV=13.6058d0, eV2Ry=1d0/Ry2eV
      end module common_constants

************************************************************************
      module parallel_operation
*     utilities for parallel operation
************************************************************************
      implicit none
      integer :: NPROC=0        ! # of processors
      private
      public :: get_np_id, init_workshare, pdprod
      public :: get_NPROC, check_NPROC

      contains

      logical function NPROC_initialized()
      implicit none
      NPROC_initialized=.not.(NPROC.eq.0)
      end function NPROC_initialized

      subroutine get_NPROC(np)
      implicit none
      integer,intent(out)::np
      if(.not. NPROC_initialized()) call set_NPROC
      np=NPROC
      end subroutine get_NPROC

      subroutine check_NPROC(np)
      implicit none
      integer,intent(in) :: np
      if(.not. NPROC_initialized()) call set_NPROC
      if(np.ne.NPROC) then
         write(*,*) 'check_NPROC: np mismatch', np, NPROC
         stop
      end if
      end subroutine check_NPROC

      subroutine set_NPROC
      implicit none
      integer :: np, id, NPROC_OLD
      if(NPROC_initialized())then
         write(*,*) 'set_NPROC: *** WARNING ***'
         write(*,*) 'redundant call of set_NPROC'
         write(*,*) 'the old value of NPROC=',NPROC
      end if
      NPROC_OLD=NPROC
      np=1 ! if OMP is disabled, NPROC=1
C$OMP PARALLEL
C$OMP SINGLE
c     for parallelism, save # of procs. to NPROC
C$    call get_np_id(np,id)
C$OMP END SINGLE
C$OMP END PARALLEL
      NPROC=np
      write(*,*) 'set_NPROC: NPROC=',NPROC
      if(NPROC_OLD .ne. 0 .and. NPROC.ne.NPROC_OLD)
     $     STOP 'set_NPROC: NPROC changed'
      end subroutine set_NPROC

      subroutine set_workshare_boundary(N,np,id,lb,ub)
      implicit none
c     input   N: size of problem
c            np: # of processors
c            id: process id number (stating from 1, unlike omp_get_thread_num)
c     output lb: lower bound of index
c            ub: upper bound of index
      integer, intent(in)  :: N, np, id
      integer, intent(out) :: lb, ub
      integer              :: share
      share=(N+np-1)/np
      lb=(id-1)*share+1
      ub=min(id*share,N)
      end subroutine set_workshare_boundary

      subroutine get_np_id(np,id)
C$    use omp_lib
      implicit none
      integer, intent(out) :: np, id
      np=1  ! for non-openmp operation
      id=1  ! for non-openmp operation
C$    if(.not. omp_in_parallel())stop'init_workshare: not in parallel'
C$    np=omp_get_num_threads()
C$    id=omp_get_thread_num()+1
      end subroutine get_np_id

      subroutine init_workshare(N,np,id,lb,ub)
      implicit none
c     This routine should be called from omp-parallel region.
c     input  N: size of problem
c     output lb: lower bound of index
c            ub: upper bound of index
c            np: # of processors
c            id: process id number (stating from 1, unlike omp_get_thread_num)
      integer, intent(in)  :: N
      integer, intent(out) :: np, id, lb, ub
      call get_np_id(np,id)
      call set_workshare_boundary(N,np,id,lb,ub)
      end subroutine init_workshare

      function pdprod(v1,v2)
      use common_constants
      implicit none
      real(knd):: pdprod, s
      real(knd),intent(in)::v1(:),v2(:)
      integer  :: i, nd, np, id, lower_bound, upper_bound
      nd=min(size(v1),size(v2))
      s=0d0
C$OMP PARALLEL DEFAULT(PRIVATE), SHARED(nd,v1,v2), REDUCTION(+:s)
      call init_workshare(nd,np,id,lower_bound,upper_bound)
      do i=lower_bound, upper_bound
         s=s+v1(i)*v2(i)
      end do
C$OMP END PARALLEL
      pdprod=s
      end function pdprod
      end module parallel_operation

************************************************************************
      module sort
*     Sorting Algorithms
************************************************************************
      private
      public :: imsort, imsorts, lmsort

      contains

************************************************************************
*     Merge Sort Program (recursive version)
*         written by S. Yamamoto, June 20, 2003
*         modified on Aug 15, 2005
*         (useless conditional branch was removed)
*         modified on Feb 25, 2007
*         (iw(N) -> allocatable)
*         parallelized on March 13, 2007           
* When there are some elements which have the same value,
* the original order isn't conserved.
* This program doesn't need the array for the temporary work area.
*     Inputs
*     N:   number of elements
*     value: given list of values
*     Outputs (On exit, values are not altered.)
*     index: (value(index(i)),i=1,N) is a non-decreasing series.
************************************************************************
      subroutine imsort(N,value,index)
      use parallel_operation, only: init_workshare, get_nproc
      implicit none
      integer, intent(in) :: N
      integer, intent(out):: index(N)
      integer :: i, p, q, np, id, lb, ub, m, mh
                                ! {l,u}b:{lower,upper}bound
      integer, allocatable :: iw(:), lub(:)
                                ! iw: work area, lub: list of ub
      integer, intent(in) ::value(N)


      call get_NPROC(np)

      allocate(iw(N))

      select case (np)
      case(:0)
         stop 'imsort: Illegal # of procs.'
      case(1)

*     serial branch

         do i=1, N
            index(i)=i
         end do
         call imsub(N,value(1),index(1),iw(1))

*     end of serial branch

      case(2:)

*     parallel branch

      allocate(lub(np))

C$OMP PARALLEL PRIVATE(np,id,lb,ub,i,p,q,m,mh)
C$OMP.SHARED(N,index,value,iw)
         call init_workshare(N,np,id,lb,ub)
c     N : # of works
c     np: # of processors
c     id: id number (from 1 to np)
c     lb: lower bound of work
c     ub: upper bound of work
         lub(id)=ub

         do i=lb, ub
            index(i)=i
         end do
         call imsub(ub-lb+1,value(1),index(lb),iw(lb))
         p=1
         q=p
         do while (p .lt. np)
C$OMP BARRIER
            if(iand(id-1,q) .eq. 0)then
               if(id+p .le. np) then
                  m =lub(min(id+2*p-1,np))-lb+1
                  mh=lub(id+p-1)-lb+1
                  call imerge(m,mh,value(1),index(lb),iw(lb))
               end if
            end if
            p=p*2
            q=q+p
         end do
C$OMP END PARALLEL

      deallocate(lub)

*     end of parallel branch

      end select

      deallocate(iw)
      end subroutine imsort


      subroutine imsorts(N,value,index) ! serial version
      implicit none
      integer, intent(in) :: N
      integer, intent(out):: index(N)
      integer :: i
      integer, allocatable:: iw(:)
      integer, intent(in) :: value(N)
      allocate(iw(N))
      do i=1, N
         index(i)=i
      end do
      call imsub(N,value(1),index(1),iw(1))
      deallocate(iw)
      end subroutine imsorts


      recursive subroutine imsub(m,value,index,iw)
      implicit none
      integer,intent(in)   :: m
      integer,intent(inout):: index(*), iw(*)
      integer,intent(in) :: value(*)
      integer :: mh
      if(m.le.1) return
      mh=m/2
*     Divide & Conquer
      call imsub(mh,  value(1),index(1),   iw(1)   )
      call imsub(m-mh,value(1),index(mh+1),iw(mh+1))
      call imerge(m,mh,value(1),index(1),iw(1))
      end subroutine imsub

      subroutine imerge (m,mh,value,index,iw)
      implicit none
      integer,intent(in) :: m,mh
      integer,intent(inout):: index(*)
      integer,intent(out)  :: iw(*)
      integer,intent(in) :: value(*)
      integer :: i,j,k
*     merging two halves
      i=1 ; j=mh+1 ; k=1
      do while( (i .le. mh) .and. (j .le. m) )
         if( value(index(i)) .le. value(index(j)) ) then
            iw(k)=index(i)
            i=i+1
         else
            iw(k)=index(j)
            j=j+1
         end if
         k=k+1
      end do
*      When no elements are left in the upper half after merging, the
*     rest of the lower half should be moved to the end of the list.
      do i=i, mh
         index(m-mh+i)=index(i)
      end do
*     If some elements are left in upper half,
*     there is nothing to be done to sort, any more.

*     copying the merged (sorted) list to the original storage ...
      do i=1, k-1
         index(i)=iw(i)
      end do
      end subroutine imerge


************************************************************************
*     Merge Sort Program (recursive version)
*         written by S. Yamamoto, June 20, 2003
*         modified on Aug 15, 2005
*         (useless conditional branch was removed)
*         modified on Feb 25, 2007
*         (iw(N) -> allocatable)
*         parallelized on March 13, 2007           
* When there are some elements which have the same value,
* the original order isn't conserved.
* This program doesn't need the array for the temporary work area.
*     Inputs
*     N:   number of elements
*     value: given list of values
*     Outputs (On exit, values are not altered.)
*     index: (value(index(i)),i=1,N) is a non-decreasing series.
************************************************************************
      subroutine lmsort(N,value,index)
      use parallel_operation, only: init_workshare, get_nproc
      implicit none
      integer, intent(in) :: N
      integer, intent(out):: index(N)
      integer :: i, p, q, np, id, lb, ub, m, mh
                                ! {l,u}b:{lower,upper}bound
      integer, allocatable :: iw(:), lub(:)
                                ! iw: work area, lub: list of ub
      integer(8), intent(in):: value(N)


      call get_NPROC(np)

      allocate(iw(N))

      select case (np)
      case(:0)
         stop 'lmsort: Illegal # of procs.'
      case(1)

*     serial branch

         do i=1, N
            index(i)=i
         end do
         call lmsub(N,value(1),index(1),iw(1))

*     end of serial branch

      case(2:)

*     parallel branch

      allocate(lub(np))

C$OMP PARALLEL PRIVATE(np,id,lb,ub,i,p,q,m,mh)
C$OMP.SHARED(N,index,value,iw)
         call init_workshare(N,np,id,lb,ub)
c     N : # of works
c     np: # of processors
c     id: id number (from 1 to np)
c     lb: lower bound of work
c     ub: upper bound of work
         lub(id)=ub

         do i=lb, ub
            index(i)=i
         end do
         call lmsub(ub-lb+1,value(1),index(lb),iw(lb))
         p=1
         q=p
         do while (p .lt. np)
C$OMP BARRIER
            if(iand(id-1,q) .eq. 0)then
               if(id+p .le. np) then
                  m =lub(min(id+2*p-1,np))-lb+1
                  mh=lub(id+p-1)-lb+1
                  call lmerge(m,mh,value(1),index(lb),iw(lb))
               end if
            end if
            p=p*2
            q=q+p
         end do
C$OMP END PARALLEL

      deallocate(lub)

*     end of parallel branch

      end select

      deallocate(iw)
      end subroutine lmsort

      recursive subroutine lmsub(m,value,index,iw)
      implicit none
      integer,intent(in)   :: m
      integer,intent(inout):: index(*), iw(*)
      integer(8),intent(in):: value(*)
      integer :: mh
      if(m.le.1) return
      mh=m/2
*     Divide & Conquer
      call lmsub(mh,  value(1),index(1),   iw(1)   )
      call lmsub(m-mh,value(1),index(mh+1),iw(mh+1))
      call lmerge(m,mh,value(1),index(1),iw(1))
      end subroutine lmsub

      subroutine lmerge (m,mh,value,index,iw)
      implicit none
      integer,intent(in) :: m,mh
      integer,intent(inout):: index(*)
      integer,intent(out)  :: iw(*)
      integer(8),intent(in):: value(*)
      integer :: i,j,k
*     merging two halves
      i=1 ; j=mh+1 ; k=1
      do while( (i .le. mh) .and. (j .le. m) )
         if( value(index(i)) .le. value(index(j)) ) then
            iw(k)=index(i)
            i=i+1
         else
            iw(k)=index(j)
            j=j+1
         end if
         k=k+1
      end do
*      When no elements are left in the upper half after merging, the
*     rest of the lower half should be moved to the end of the list.
      do i=i, mh
         index(m-mh+i)=index(i)
      end do
*     If some elements are left in upper half,
*     there is nothing to be done to sort, any more.

*     copying the merged (sorted) list to the original storage ...
      do i=1, k-1
         index(i)=iw(i)
      end do
      end subroutine lmerge

      end module sort


      module file_io

      private
      public :: vacant_unit

      contains

      integer function vacant_unit()
      integer, parameter :: MAX_UNIT=99
      integer            :: IUNIT
      logical            :: not_found
*     debug
c!!!      do IUNIT=1, MAX_UNIT
c!!!         open(IUNIT)
c!!!      end do
      do IUNIT=1, MAX_UNIT
         inquire(IUNIT,OPENED=not_found)
         if(.not. not_found) exit
      end do
      if(not_found) then
         write(*,'("file_io%vacant_unit: all UNIT is busy.")')
         stop
      end if
c!!!      write(*,'("file_io%vacant_unit: IUNIT=",I3)') IUNIT
      vacant_unit=IUNIT
      end function vacant_unit

      end module file_io

      module sparse_matrix
      use common_constants
      use parallel_operation, only: init_workshare, get_np_id,
     $     get_NPROC
      implicit none
      type sparsemat
      integer :: dim, nelem
c     fortran 95 does not allow allocatable array in derived type
      real(knd), pointer :: val(:)
      integer,   pointer :: id1(:),id2(:)
cAL      real(knd), allocatable :: val(:)
cAL      integer,   allocatable :: id1(:),id2(:)
*     following variables are needed for parallel calc.
      integer :: NPROC
      integer,   pointer :: st_p(:)
cAL      integer,   allocatable :: st_p(:)
      end type sparsemat

      contains

      logical function check_matrix(X)
      type(sparsemat),intent(in)::x
      integer :: i
      logical :: flag
      flag=.true.
      if(.not. associated(X%val))then
cAL      if(.not. allocated(X%val))then
         write(*,'("check_matrix: %val is not allocated")')
         flag=.false.
      end if
      if(.not. associated(X%id1))then
cAL      if(.not. allocated(X%id1))then
         write(*,'("check_matrix: %id1 is not allocated")')
         flag=.false.
      end if
      if(.not. associated(X%id2))then
cAL      if(.not. allocated(X%id2))then
         write(*,'("check_matrix: %id2 is not allocated")')
         flag=.false.
      end if
      do i=1, X%nelem
         if(X%id1(i) <= 0 .or. X%id1(i) > X%dim)then !Fixed May 9, 2008
            write(*,'("check_matrix: id1(",I10,")=",I10)')
     $           i, X%id1(i)
            flag=.false.
         end if
         if(X%id2(i).le.0 .or. X%id2(i) > X%dim)then
            write(*,'("check_matrix: id2(",I10,")=",I10)')
     $           i, X%id2(i)
            flag=.false.
         end if
      end do
      check_matrix=flag
      end function check_matrix

c!!!      logical function check_order_1(X)
c!!!      type(sparsemat),intent(in)::x
c!!!      integer :: i, np, id, lb, ub, M
c!!!      logical :: flag
c!!!      M=X%nelem
c!!!      flag=.true.
c!!!C$OMP PARALLEL PRIVATE(np,id,lb,ub) REDUCTION(.and.:flag)
c!!!      call init_workshare(M,np,id,lb,ub)
c!!!      do i=lb, min(ub,M-1)
c!!!         if(X%id1(i) > X%id1(i+1)) then
c!!!            write(*,'("i=",I10,", X%id1(i)=",I10,
c!!!     $           "X%id1(i+1)=",I10)')i, X%id1(i), X%id1(i+1)
c!!!            flag=.false.
c!!!         end if
c!!!      end do
c!!!C$OMP END PARALLEL
c!!!      check_order_1=flag
c!!!      end function check_order_1

      subroutine allocate_matrix(N,M,X)
      type(sparsemat),intent(out)::x
      integer, intent(in) :: N, M
      integer :: ierr
      if(associated(X%val)) then
cAL      if(allocated(X%val)) then
         if(size(X%val).ne.M .or. X%dim.ne.N)then
            deallocate(X%val)
            deallocate(X%id1)
            deallocate(X%id2)
         end if
      end if
      if(.not. associated(X%val))then
cAL      if(.not. allocated(X%val))then
         allocate(X%val(M),STAT=ierr)
         if(ierr.ne.0) call allocation_error("val",ierr)
         allocate(X%id1(M))
         if(ierr.ne.0) call allocation_error("id1",ierr)
         allocate(X%id2(M))
         if(ierr.ne.0) call allocation_error("id2",ierr)
         X%dim  =N
         X%nelem=M
      end if
      contains
      subroutine allocation_error(VNAME,nerr)
      character(len=*):: VNAME
      integer         :: nerr
      intent(in)      :: VNAME, nerr
      write(*,'("allocate_matrix: allocation error (",A10,")")')
     $     VNAME(1:min(len(VNAME),10))
      write(*,'("   error status=",I4)') nerr
      end subroutine allocation_error
      end subroutine allocate_matrix

      subroutine parallel_setup(X)
      use sort, only: imsort, imsorts
      type(sparsemat),intent(inout) :: X
      integer, allocatable :: idx(:), upper_bound(:)
      integer, allocatable :: iw1(:), iw2(:)
      real(knd),allocatable:: wal(:)
      integer :: np, id, lb, ub, N, M, i, jd

      call get_NPROC(np)
      if(associated(X%st_p)) then
cAL      if(allocated(X%st_p)) then
         if(np==X%NPROC) return
         deallocate(X%st_p)
         write(*,'("parallel_setup: relocation occured")')
         write(*,'("   old size=",I2,", news ize=",I2)') X%NPROC, np
      end if
      N=X%dim
      M=X%nelem
      allocate(X%st_p(np+1))
      allocate(upper_bound(0:np))
      allocate(idx(M))
      allocate(iw1(M))
      allocate(iw2(M))
      allocate(wal(M))
      call imsort(M,X%id1,idx)
      X%NPROC=np
      X%st_p(1)=1
      X%st_p(np+1)=M+1
      upper_bound=0
C$OMP PARALLEL PRIVATE(np,id,lb,ub,jd,i)
      call init_workshare(N,np,id,lb,ub)
      upper_bound(id)=ub
      call init_workshare(M,np,id,lb,ub)
      do i=lb, ub
         iw1(i)=X%id1(idx(i))
         iw2(i)=X%id2(idx(i))
         wal(i)=X%val(idx(i))
      end do
      jd=1
C$OMP BARRIER
*     partitioning id1 into np regions.
*     in i-th region,  id1 satisfies
*     upper_bound(i-1)<id1<=upper_bound(i)
      do i=lb, min(ub,M-1)
         do while( iw1(i+1) > upper_bound(jd) )
            jd=jd+1
         end do
c     now iw1(i+1) <= upper_bound(jd)
         if(  iw1(i)  <= upper_bound(jd-1) )then
c            write(*,*) 'id=',id, 'jd=',jd, 'i=',i
            X%st_p(jd)=i+1
         end if
      end do
C$OMP BARRIER
      lb=X%st_p(id)
      ub=X%st_p(id+1)-1
      i=ub-lb+1
      call imsorts(i,iw2(lb:ub),idx(lb:ub))
      do i=lb, ub
         X%id1(i)=iw1(idx(i)+lb-1)
         X%id2(i)=iw2(idx(i)+lb-1)
         X%val(i)=wal(idx(i)+lb-1)
      end do
      X%id1(lb:ub)=iw1(lb:ub)
      X%id2(lb:ub)=iw2(lb:ub)
      X%val(lb:ub)=wal(lb:ub)
C$OMP END PARALLEL
c!!!      write(*,*) upper_bound
c!!!      write(*,*) X%st_p
c!!!      write(*,*) X%id1(X%st_p(1))
c!!!      write(*,*) X%id1(X%st_p(2)),X%id1(X%st_p(2)-1)
c!!!      write(*,*) X%id1(X%st_p(3)),X%id1(X%st_p(3)-1)
c!!!      write(*,*) X%id1(X%st_p(4)),X%id1(X%st_p(4)-1)
c!!!      pause
      deallocate(wal)
      deallocate(iw1)
      deallocate(iw2)
      deallocate(idx)
      deallocate(upper_bound)
      end subroutine parallel_setup

      subroutine apply_matrix(X,vec,nvec)
      use parallel_operation, only: get_np_id, init_workshare
      type(sparsemat), intent(inout) :: X
      real(knd), intent(in)  ::  vec(:)
      real(knd), intent(out) :: nvec(:)
      integer :: np, id, lb, ub, N, i
      N=X%dim
      if(N > size(vec)) then
         write(*,'("apply_matrix: too small vec")')
         write(*,'("N=",I10,", size(vec )=",I10)')N, size(vec)
         stop
      end if
      if(N > size(nvec)) then
         write(*,'("apply_matrix: too small nvec")')
         write(*,'("N=",I10,", size(nvec)=",I10)')N, size(nvec)
         stop
      end if

c     Most simple implementation
c      nvec=0d0
c      do i=1, X%nelem
c         nvec(X%id1(i))=nvec(X%id1(i))+X%val(i)*vec(X%id2(i))
c      end do
c      return
c     end Most simple implementation

      call parallel_setup(X)

C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call init_workshare(N,np,id,lb,ub)
      nvec(lb:ub)=0d0
C$OMP END PARALLEL      

C$OMP PARALLEL PRIVATE(np,id,lb,ub)
      call get_np_id(np,id)
c!!!      call init_workshare(N,np,id,lb,ub) ! private(np,id,lb,ub)
c!!!      write(*,*) X%st_p
      do i=X%st_p(id), X%st_p(id+1)-1
c!!!         if( X%id1(i)< lb .or. X%id1(i) > ub)then
c!!!            write(*,*) X%id1(i),X%id2(i),lb,ub
c!!!            stop
c!!!         end if
         nvec(X%id1(i))=nvec(X%id1(i))+X%val(i)*vec(X%id2(i))
      end do
c!!!      if(id==1)write(*,*) X%id1(X%st_p(id):X%st_p(id+1)-1)
c!!!      if(id==1)write(*,*) X%id2(X%st_p(id):X%st_p(id+1)-1)
C$OMP END PARALLEL
      end subroutine apply_matrix

      end module sparse_matrix


************************************************************************
      module hamiltonian
*     Needed in the module COCG
************************************************************************
      use sparse_matrix
      implicit none
      type(sparsemat)   :: Ham

      type list_of_base
         integer :: number_of_base
         integer, pointer :: index_of_atom_orb_pair(:,:) !(orb, atom)
      end type list_of_base

      type(list_of_base):: base
      
      contains

      subroutine read_hamiltonian(FNAME)
      use file_io
      character(len=*),intent(in) :: FNAME
      integer :: IUNIT, ndim, nelem, mnel_row
      if(index(FNAME,".dat").ne.0) then
         call read_binary_hamiltonian
         goto 99
      end if
      IUNIT=vacant_unit()
      open(IUNIT,FILE=FNAME,STATUS="OLD")
      read(IUNIT,*) ndim
      write(*,'("DIM=",I8)') ndim
      read(IUNIT,*) nelem
      write(*,'("# of non-zero elements=",I10)') nelem
      call allocate_matrix(ndim,nelem,Ham)
      call read_hamiltonian_elements(Ham)
      call read_hamiltonian_base(base)
      close(IUNIT)
 99   continue
      if(.not. check_matrix(Ham)) stop  'read_hamiltonian: error(M)'
      if(.not. check_hamiltonian()) stop'read_hamiltonian: error(H)'
      return

      contains

      subroutine read_binary_hamiltonian
      integer :: SPARSE_CALC,ns,TotOcc,TotSz,norb,nsp, i

      IUNIT=vacant_unit()
      open(IUNIT,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      read(IUNIT) SPARSE_CALC,ns,ndim,TotOcc,TotSz,norb,nsp
                                ! ns=#(upper non-diagonal block)
      nelem=2*ns+ndim
      if(SPARSE_CALC.ne.1) stop 'wrong input'
      call allocate_matrix(ndim,nelem,Ham)
      read(IUNIT) (Ham%id1(i),Ham%id2(i),i=1,ns)
      read(IUNIT) (Ham%val(i),i=1,ns)
      read(IUNIT) (Ham%val(i),i=2*ns+1,2*ns+ndim)
      Ham%id2(ns+1:2*ns)=Ham%id1(1:ns)
      Ham%id1(ns+1:2*ns)=Ham%id2(1:ns)
      Ham%val(ns+1:2*ns)=Ham%val(1:ns)
      write(*,'("Minimum value of diagonal part",ES22.14)')
     $     minval(Ham%val(2*ns+1:2*ns+ndim))
      write(*,'("Maximum value of diagonal part",ES22.14)')
     $     maxval(Ham%val(2*ns+1:2*ns+ndim))
      do i=1, ndim
         Ham%id1(2*ns+i)=i
         Ham%id2(2*ns+i)=i
      end do
      close(IUNIT)
      end subroutine read_binary_hamiltonian

      logical function check_hamiltonian()
      use sort, only: lmsort
      integer(8),allocatable::val(:)
      integer   ,allocatable::idx(:)
      integer :: i, k1, k2, ND, N, count
      ND=Ham%dim
      N =Ham%nelem
      allocate(val(N))
      allocate(idx(N))
      count=0
      do i=1, N
         if(Ham%id1(i) > Ham%id2(i))then
            count=count+1
            val(i)=Ham%id1(i)
            val(i)=(val(i)-1)*ND+Ham%id2(i)
         else
            val(i)=Ham%id2(i)
            val(i)=(val(i)-1)*ND+Ham%id1(i)
         end if
      end do
      write(*,'("check_hamiltonian:")')
      write(*,'("   Dim   =                      ",I10)') ND
      write(*,'("   Nelem =                      ",I10)') N
      write(*,'("   # of Ham%id1(i) > Ham%id2(i)=",I10)')count
      write(*,'("   (Ham%Nelem-Ham%dim)/2       =",I10)')(N-ND)/2
      check_hamiltonian=(count == (N-ND)/2)
      call lmsort(N,val,idx)
      count=0
      do i=1, N
         k1=idx(i)
         if(Ham%id1(k1)==Ham%id2(k1))cycle
         if(count==1) then
            count=0
            cycle
         end if
         k2=idx(i+1)
         if(  Ham%id1(k1)/=Ham%id2(k2).or.
     $        Ham%id2(k1)/=Ham%id1(k2)    ) then
            write(*,*) i, Ham%id1(k1), Ham%id2(k1),
     $                    Ham%id1(k2), Ham%id2(k2)
            write(*,*) val(idx(i-2:i+2))
            check_hamiltonian=.false.
            exit
         end if
         if( abs(Ham%val(k1)-Ham%val(k2)) .gt. 1d-10 )then
            write(*,*) i, Ham%id1(k1), Ham%id2(k1),
     $                    Ham%id1(k2), Ham%id2(k2)
            write(*,*)    Ham%val(k1), Ham%val(k2)
            check_hamiltonian=.false.
            exit
         end if
         count=1
      end do
      write(*,'("   symmetric=",L1)') check_hamiltonian
      end function check_hamiltonian


      subroutine read_hamiltonian_elements(H)
      type(sparsemat),intent(inout)::H
      integer   :: i, j, k, p, irow, inum
      real(knd) :: val_r, val_i
      logical   :: in_order

      in_order=.true.
      irow=0
      do p=1, nelem
         read(IUNIT,*) i, j, k, val_r, val_i
         if(i /= p)then
            write(*,'("read_hamiltonian: missing elements p=",
     $           I10)') p
            stop
         end if
         if(j < irow) then
            if(.not.(inum==1 .and. H%id1(p-1)==H%id2(p-1))) then
c     diagonal elements is stored in leading line of the
c     column data block.
c!!!               write(*,'("read_hamiltonian: *** WARNING ***",
c!!!     $              " row index decreased. p=",I8)') p
               in_order=.false.
            end if
         else if(j > irow)then
            inum=1
            irow=j
         else ! j=irow
            inum=inum+1
            if(H%id2(p-1).gt.k)then
               if(.not.(inum==2 .and. H%id1(p-1)==H%id2(p-1)))then
c     diagonal elements is stored in leading line of the
c     column data block.
c!!!                  write(*,'("read_hamiltonian: WARNING!",
c!!!     $                 " column index decreased. p=",I8)') p
                  in_order=.false.
               end if
            end if
         end if
         if(val_i /= 0d0)then
            write(*,'("read_sparse_hamiltonian: WARNING!",
     $           " non-zero imaginary part. p=",I8)') p
         end if
         H%id1(p)=j
         H%id2(p)=k
         H%val(p)=val_r
      end do
      end subroutine read_hamiltonian_elements

      subroutine read_hamiltonian_base(base)
      implicit none
      type(list_of_base),intent(out)::base
      integer,allocatable:: idx(:), atom(:), orbital(:)
      integer :: j, k, m
      base%number_of_base=Ham%dim
      allocate(idx    (base%number_of_base))
      allocate(atom   (base%number_of_base))
      allocate(orbital(base%number_of_base))
      do j=1, base%number_of_base
         read(IUNIT,*,END=900)idx(j), atom(j), orbital(j)
         if(idx(j)/=j) then
            write(*,'(A)') '*** Warning! ***'
            write(*,'(A)') '   The orbital index jumped.'
c            stop 'read_hamiltonian_base: error in Ser.#'
         end if
      end do
 900  continue
      if(j == 1) then
         write(*,'("   Fallback to compatibility-mode.")')
         goto 910
      else if(j /= base%number_of_base+1) then
         write(*,'("read_hamiltonian_base: *** Warning ***")')
         write(*,'("   Insufficient base data.")')
      end if
      k=maxval(atom)
      m=maxval(orbital)
      if(associated(base%index_of_atom_orb_pair))
     $     call base_error('base%index_of_atom_orb_pair') 
      allocate(base%index_of_atom_orb_pair(m,k))
      base%index_of_atom_orb_pair=0
      do j=1, base%number_of_base
         base%index_of_atom_orb_pair(orbital(j),atom(j))=idx(j)
      end do
      goto 910
 910  continue
      deallocate(idx)
      deallocate(atom)
      deallocate(orbital)
      end subroutine read_hamiltonian_base

      subroutine base_error(array_name)
      character(len=*) :: array_name
      write(*,'(A)') 'read_hamiltonian_base: ', array_name,
     $     'is already allocated'
      stop
      end subroutine base_error
      end subroutine read_hamiltonian

      subroutine apply_hamiltonian(N,vec,nvec)
      integer,  intent(in) :: N
      real(knd),intent(in) :: vec(N)
      real(knd),intent(out):: nvec(N)
      call apply_matrix(Ham,vec,nvec)
      end subroutine apply_hamiltonian
      end module hamiltonian


************************************************************************
      module control
*     manages the variables in input file
************************************************************************
      use common_constants
      implicit none

      type control_parameters
      real(knd) :: eps_ref, eps_shift
      real(knd) :: Esample, enmin, enmax
      real(knd) :: escale, denergy ! derived variables
      real(knd) :: gamma, gammaXmesh
      integer   :: nitemax, iswEscale, nenemax, iswgamma
      integer   :: target_i, target_j, target_k
      logical   :: compatibility_mode
      integer   :: number_of_orbitals_to_be_shown
      integer,pointer:: list_of_atom_orbitals_pairs(:,:)
      integer,pointer:: list_of_orbitals_to_be_shown(:)
      end type control_parameters

      contains

      subroutine normalize_parameters(x)
      type(control_parameters),intent(inout):: x

      if(x%iswEscale /= 0) then
c     when input and output data are in unit of eV
c!!!         x%escale  = Ry2eV*2
c!!!         x%Esample = x%Esample /x%escale
c!!!         x%enmin   = x%enmin   /x%escale
c!!!         x%enmax   = x%enmax   /x%escale
c!!!         x%gamma   = x%gamma   /x%escale
         write(*,'("normalize_parameter: Energy scale conversion",
     $        " is obsolescent. Turn off the ""iswEscale"".")')
         stop
      else
         x%escale=1d0
      endif

      x%denergy=(x%enmax-x%enmin)/(x%nenemax-1)
      if( x%iswgamma == 1 ) then
         x%gamma=x%gammaXmesh*x%denergy
         write(*,*) '#### gamma is replaced by denergy*gammaXmesh ####'
         write(*,*) 'gamma=',x%gamma
      end if


      end subroutine normalize_parameters

      subroutine read_parameters(FNAME,input)
      use file_io
      character(len=*),intent(in)          :: FNAME
      type(control_parameters),intent(out) :: input
      integer :: IUNIT, j,k, atom1, orb1, atom2, orb2

      input%compatibility_mode=.true. ! for comatibility

      IUNIT=vacant_unit()
      open(IUNIT,FILE=FNAME,STATUS="OLD")

      read(IUNIT,*) input%eps_ref, input%eps_shift, input%nitemax
      read(IUNIT,*) input%iswEscale, input%Esample
      read(IUNIT,*) input%enmin, input%enmax, input%nenemax
      read(IUNIT,*) input%iswgamma, input%gamma, input%gammaXmesh
      read(IUNIT,*) input%target_i, input%target_j, input%target_k
      if(input%target_i > 0) then
         ! compatibility mode
         input%number_of_orbitals_to_be_shown=1
         allocate(input%list_of_orbitals_to_be_shown(1))
         input%list_of_orbitals_to_be_shown = input%target_i
      else if(input%target_i < 0) then
         ! advanced mode
         input%compatibility_mode=.false.
         read(IUNIT,*) k
         input%number_of_orbitals_to_be_shown=k
         allocate(input%list_of_atom_orbitals_pairs(2,k))
         allocate(input%list_of_orbitals_to_be_shown(k))
         do j=1, input%number_of_orbitals_to_be_shown
            read(IUNIT,*) atom1, orb1 !, atom2, orb2
            input%list_of_atom_orbitals_pairs(1,j)=atom1
            input%list_of_atom_orbitals_pairs(2,j)=orb1
         end do
      end if
      close(IUNIT)

      input%eps_ref  =exp(log(10d0)*input%eps_ref)
      input%eps_shift=exp(log(10d0)*input%eps_shift)

      call normalize_parameters(input)

      write(6,602) input%eps_ref, input%eps_shift, input%nitemax
      write(6,604) input%iswEscale, input%Esample
      write(6,606) input%enmin, input%enmax, input%nenemax
      write(6,608) input%iswgamma, input%gamma, input%gammaXmesh
      write(6,610) input%target_i, input%target_j, input%target_k
 602  format("eps_ref,eps_shift,nitemax= ",2G14.6,I5)
 604  format("iswEscale,Esample=         ",I5,G14.6)
 606  format("enmin,enmax,nenemax=       ",2G14.6,I5)
 608  format("iswgamma,gamma,gammaXmesh= ",I5,2G14.6)
 610  format("target_i,target_j,target_k=",3I5)

      end subroutine read_parameters

      subroutine set_list_of_orbitals_to_be_shown(input)
      use hamiltonian
      implicit none
      type(control_parameters),intent(out) :: input
      integer :: j
      if(input%compatibility_mode) return
      do j=1, input%number_of_orbitals_to_be_shown
         input%list_of_orbitals_to_be_shown(j)=
     $        base%index_of_atom_orb_pair(
     $        input%list_of_atom_orbitals_pairs(2,j),
     $        input%list_of_atom_orbitals_pairs(1,j))
      end do
      end subroutine set_list_of_orbitals_to_be_shown
      end module control
