!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================

module M_lib_random_num
! Routine for random number generator
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=8
  integer, parameter   :: ip = 521
  integer, parameter   :: iq =  32
  integer              :: ir(ip)
  integer              :: iptr
!
  private
  public :: rndini
  public :: rndu
!
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initial routine for random number generator
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine rndini(iseed)
!
    implicit none
    integer, parameter     :: nq    = ip-iq
    integer, parameter     :: nbit  = 32
    integer, parameter     :: nbit1 = 31
    integer, parameter     :: ia    = 69069
    integer, intent(inout) :: iseed
    integer                :: iw(ip)
    integer                :: i,j,ih,mj,ii,ij
!
    if(iseed.le.0) then
      write(*,*) 'rndini: iseed must be positive.'
      stop
    end if
!
    if(mod(iseed,2).eq.0) then
      write(*,*) 'rndini: iseed must be odd.'
      stop
    end if
!
    do i=1,ip
      iseed=iseed*ia
      iw(i)=isign(1,iseed)
      if (i .le. 10) write(6,*)'i,iw=',i,iw(i)
    enddo  
!
    do j=1,ip
      ih=mod((j-1)*nbit,ip)
      mj=0
      do i=1,nbit1
        ii=mod(ih+i-1, ip)+1
        mj=2*mj+(iw(ii)-1)/(-2)
        ij=mod(ii+nq-1,ip)+1
        iw(ii)=iw(ii)*iw(ij)
      enddo  
      ir(j)=mj
      ii=mod(ih+nbit1,ip)+1
      ij=mod(ii+nq-1, ip)+1
      iw(ii)=iw(ii)*iw(ij)
    enddo  
!
    iptr=0
!
    return
!
  end subroutine rndini
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random number generator
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine rndu(rnd)
!
    implicit none
    integer, parameter   :: ip = 521
    real(DOUBLE_PRECISION), parameter  :: ra    =2.0d0**(-31)
    real(DOUBLE_PRECISION), parameter  :: rb    =2.0d0**(-32)
    integer :: jptr
    real(DOUBLE_PRECISION) :: rnd
!
    iptr=iptr+1
    if(iptr.gt.ip) iptr=1
    jptr=iptr-iq
    if(jptr.le.0) jptr=jptr+ip
!
    ir(iptr)=ieor(ir(iptr),ir(jptr))
    rnd=dble(ir(iptr))*ra+rb
!
    return
  end subroutine rndu
!
!
end module M_lib_random_num
