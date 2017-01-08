       subroutine CSTAB_get_optdim(n_min, n_unroll,
     &                             delta_L1, delta_L2, n_opt)
       implicit none

       integer                :: n_min, n_unroll, delta_L1, delta_L2
       integer                :: n_opt, n_delta
       integer                :: n_opt2, n_delta2

       include 'CSTAB.h'

       integer                :: i,j,k


       n_opt=n_min
       do

          n_opt   = (n_opt-1)/L1_WINDOW+1
          n_opt   = (n_opt/2)*2+1
          n_opt   = n_opt*L1_WINDOW

          n_delta = 0

          do i=1,INT((n_unroll*1.2-1.0)/L1_WAY+1)

             k=mod(i*n_opt+L1_LSIZE/2,L1_LSIZE)-L1_LSIZE/2
             if(abs(k)<=delta_L1/2)then
                n_delta=(delta_L1/2-k-1)/i+1
                goto 10000
             end if

          end do

          do i=1,INT((n_unroll*1.2-1.0)/L2_WAY+1)

             k=mod(i*n_opt+L2_LSIZE/2,L2_LSIZE)-L2_LSIZE/2
             if(abs(k)<=delta_L2/2)then
                n_delta=(delta_L2/2-k-1)/i+1
                goto 10000
             end if

          end do

10000     continue

          if(n_delta==0)exit
          n_opt = n_opt + n_delta

       end do

       return
       end subroutine

       subroutine CSTAB_adjust_base(a,b,offset)
       implicit none

       real(8)                :: a(*), b(*)
       integer                :: offset

       include 'CSTAB.h'


       call get_delta(a(1),b(1),offset)
       offset=(offset/8)

       if(offset>0)then
         offset=MOD(L2_LSIZE-MOD(+offset,L2_LSIZE),L2_LSIZE)
       else
         offset=MOD(L2_LSIZE+MOD(-offset,L2_LSIZE),L2_LSIZE)
       endif

       return
       end subroutine
!
       subroutine CSTAB_adjust_page(a,b,offset)
       implicit none

       real(8)                :: a(*), b(*)
       integer                :: offset

       include 'CSTAB.h'


       call get_delta(a(1),b(1),offset)
       offset=(offset/8)
       if(offset>0)then
         offset=MOD(PAGE_LSIZE-MOD(+offset,PAGE_LSIZE),PAGE_LSIZE)
       else
         offset=MOD(PAGE_LSIZE+MOD(-offset,PAGE_LSIZE),PAGE_LSIZE)
       endif

       return
       end subroutine
!
       subroutine CSTAB_round_offset(offset)
       implicit none

       integer                :: offset

       include 'CSTAB.h'

       offset=MOD(offset,L2_LSIZE)

       return
       end subroutine
!
