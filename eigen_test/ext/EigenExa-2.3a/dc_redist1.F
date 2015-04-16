      subroutine dc_redist1(n, NB, a, b, ldm, wk, lwk)
      use MPI
      use eigen_devel
      implicit NONE

      integer, intent(IN)    :: n, NB, ldm, lwk
      real(8), intent(IN)    :: a(ldm, *)
      real(8), intent(OUT)   :: b(ldm, *)
      real(8), intent(INOUT) :: wk(*)

      real(8), pointer       :: wk1(:, :), wk2(:, :)

      integer :: i,j,j0,k,k0,l,lx,IERR
      integer :: iNQ
      integer :: iblk_, jblk_
      integer :: iblk, jblk
      integer :: idist, ir_size, is_size
      integer :: his_rank, her_rank
      integer, pointer :: ir_sz(:, :, :)
      integer :: temp(2)

      integer :: NQ, NBQ


      iblk_ = (n-1)/y_nnod+1
      iblk  = (iblk_-1)/NB+1
      jblk_ = (n-1)/x_nnod+1
      jblk  = (jblk_-1)/NB+1

!
! Calculate the size of required buffer for re-distribution;
! divide the Block width by NQ and get the buffer shape
! as (1:ldm, lx) where lx is the max number of vectors to be
! transfered such that ldm * lx * 2 <= lwk (lwk is the size of
! working vector passed from the callee procedure.

      NQ = 0
      lx = 0
!
      if ( TRD_inod == 1 ) then
      do iNQ=1,NB

         NQ = iNQ
         NBQ = (NB-1)/NQ+1

! Emulate the communication and evaluate data size
         lx = 0
         do i=1,iblk
         do j0=1,NB,NBQ
         do idist = 1, y_nnod-1

! stride communication with the upper and the lower processess

            his_rank = MOD(y_inod-1-idist+y_nnod,y_nnod)+1
            her_rank = MOD(y_inod-1+idist+y_nnod,y_nnod)+1

! count up the number of vectors which should be sent to her_rank
            l = 0
            do j=j0,MIN(NB,j0+NBQ-1)
              K=j+((y_inod-1)+(i-1)*y_nnod)*NB
              if ( K <= n .AND. MOD(K-1,y_nnod)+1 == her_rank ) then
                 l = l + 1
              end if
            end do
            lx = MAX(l, lx)

! count up the number of vectors which should be recieved from his_rank
            l = 0
            do j=j0,MIN(NB,j0+NBQ-1)
              K=j+((his_rank-1)+(i-1)*y_nnod)*NB
              if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
                 l = l + 1
              end if
            end do
            lx = MAX(l, lx)

         end do
         end do
         end do

         if ( ldm*lx*2 <= lwk ) then
            exit
         end if

      end do
      end if

      temp(1) = NQ
      temp(2) = lx
      call MPI_Bcast( temp(1), 2, MPI_INTEGER, 0,
     $                       TRD_COMM_WORLD, IERR )
      NQ  = temp(1)
      NBQ = (NB-1)/NQ+1
      lx  = temp(2)


      allocate ( ir_sz(y_nnod, iblk, NQ) )

      do i=1,iblk
      do j0=1,NB,NBQ
      do idist = 1, y_nnod-1

! Find up the number of elements which should be recieved from his_rank
         his_rank = MOD(y_inod-1-idist+y_nnod,y_nnod)+1
         l = 0
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              l = l + 1
           end if
         end do
         ir_sz(idist,i,(j0-1)/NBQ+1) = l

      end do
      end do
      end do


      if ( ldm*lx*2 <= lwk ) then
! The most eco case
         call dc_redist1_sub ( n, NB, a, b, ldm,
     &                       wk(1), wk(1+ldm*lx),
     &                       ir_sz, y_nnod, iblk, NBQ )
      else if ( ldm*lx <= lwk ) then
! One buffer is allocated
         allocate ( wk2(1:ldm, 1:lx) )
         call dc_redist1_sub ( n, NB, a, b, ldm,
     &                       wk, wk2(1,1),
     &                       ir_sz, y_nnod, iblk, NBQ )
         deallocate ( wk2 )
      else
! Both buffers are allocated
         allocate ( wk1(1:ldm, 1:lx) )
         allocate ( wk2(1:ldm, 1:lx) )
         call dc_redist1_sub ( n, NB, a, b, ldm,
     &                       wk1(1,1), wk2(1,1),
     &                       ir_sz, y_nnod, iblk, NBQ )
         deallocate ( wk1 )
         deallocate ( wk2 )
      end if


      deallocate ( ir_sz )

      return
      end subroutine dc_redist1


      subroutine dc_redist1_sub(n, NB, a, b, ldm,
     $               wk1, wk2, ir_sz, l1, l2, NBQ)
      use eigen_devel
      implicit NONE

      integer, intent(IN)    :: n, NB, ldm, l1, l2, NBQ
      real(8), intent(IN)    :: a(ldm,*)
      real(8), intent(OUT)   :: b(ldm,*)
      real(8), intent(INOUT) :: wk1(ldm,*), wk2(ldm,*)
      integer, intent(IN)    :: ir_sz(l1,l2,*)

      integer :: i,j,j0,j1,j2,k,k0,k1,k2,l,L0,lx
      integer :: iblk_, jblk_
      integer :: iblk,  jblk
      integer :: idist, ir_size, is_size
      integer :: iq_r, iq_s, his_rank, her_rank

      integer, pointer :: l_(:)


      iblk_ = (n-1)/y_nnod+1
      iblk  = (iblk_-1)/NB+1
      jblk_ = (n-1)/x_nnod+1
      jblk  = (jblk_-1)/NB+1


      allocate ( l_(0:NB) )


      do i=1,iblk
      do j0=1,NB,NBQ

         K1 = MOD((y_inod-1)*NB,y_nnod)
         j2 = MOD(y_inod-K1+y_nnod-1,y_nnod)+1
         j2 = j2 + ((j0-1)/y_nnod)*y_nnod
         if ( j2 < j0 ) j2 = j2 + y_nnod

!$OMP PARALLEL DO PRIVATE(j,K,K0,l)
#if 0
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              K0=(K-1)/y_nnod+1
              b(1:ldm, K0) = a(1:ldm, j+NB*(i-1))
           end if
         end do
#else
         do j=j2,MIN(NB,j0+NBQ-1),y_nnod
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n ) then
              K0=(K-1)/y_nnod+1
              b(1:ldm, K0) = a(1:ldm, j+NB*(i-1))
           end if
         end do
#endif
!$OMP END PARALLEL DO

      do idist = 1, y_nnod-1

         her_rank = MOD(y_inod-1+idist+y_nnod,y_nnod)+1
         his_rank = MOD(y_inod-1-idist+y_nnod,y_nnod)+1

         ir_size=ldm*ir_sz(idist,i,(j0-1)/NBQ+1)
         if ( ir_size > 0 ) then
         call irecv_dbl ( wk2, ir_size, his_rank, iq_r, y_COMM_WORLD )
         end if

!
! l_ contains the offset which represents XXXXXX ....
!
         K1 = MOD((y_inod-1)*NB,y_nnod)
         j1 = MOD(her_rank-K1+y_nnod-1,y_nnod)+1
         j1 = j1 + ((j0-1)/y_nnod)*y_nnod
         if ( j1 < j0 ) j1 = j1 + y_nnod

         K2 = MOD((his_rank-1)*NB,y_nnod)
         j2 = MOD(y_inod-K2+y_nnod-1,y_nnod)+1
         j2 = j2 + ((j0-1)/y_nnod)*y_nnod
         if ( j2 < j0 ) j2 = j2 + y_nnod

         l = 0; l_(0)=0
#if 0
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == her_rank ) then
              l = l + 1
           end if
           l_(j-j0+1) = l
         end do
#else
         l=0
         do j=j1,MIN(NB,j0+NBQ-1),y_nnod
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           l_(j-j0) = l
           if ( K <= n ) then
              l = l + 1
           end if
         end do
#endif
         lx = l
!$OMP PARALLEL DO PRIVATE(j,K,K0,l)
#if 0
         do j=j0,MIN(NB,j0+NBQ-1)
           l = l_(j-j0)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == her_rank ) then
              l = l + 1
              wk1(1:ldm, l) = a(1:ldm, j+NB*(i-1))
           end if
         end do
#else
         do j=j1,MIN(NB,j0+NBQ-1),y_nnod
           l = l_(j-j0)
           K=j+((y_inod-1)+(i-1)*y_nnod)*NB
           if ( K <= n ) then
              l = l + 1
              wk1(1:ldm, l) = a(1:ldm, j+NB*(i-1))
           end if
         end do
#endif
!$OMP END PARALLEL DO

         is_size=ldm*lx
         if ( is_size > 0 ) then
         call isend_dbl ( wk1, is_size, her_rank, iq_s, y_COMM_WORLD )
         end if

         if ( ir_size > 0 ) then
         call wait_dbl ( iq_r )

         l = 0; l_(0)=0
#if 0
         do j=j0,MIN(NB,j0+NBQ-1)
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              l = l + 1
           end if
           l_(j-j0+1) = l
         end do
#else
         do j=j2,MIN(NB,j0+NBQ-1),y_nnod
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           l_(j-j0) = l
           if ( K <= n ) then
              l = l + 1
           end if
         end do
#endif
         lx = l
!$OMP PARALLEL DO PRIVATE(j,K,K0,l)
#if 0
         do j=j0,MIN(NB,j0+NBQ-1)
           l = l_(j-j0)
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           if ( K <= n .AND. MOD(K-1,y_nnod)+1 == y_inod ) then
              l = l + 1
              K0=(K-1)/y_nnod+1
              b(1:ldm, K0) = wk2(1:ldm, l)
           end if
         end do
#else
         do j=j2,MIN(NB,j0+NBQ-1),y_nnod
           l = l_(j-j0)
           K=j+((his_rank-1)+(i-1)*y_nnod)*NB
           if ( K <= n ) then
              l = l + 1
              K0=(K-1)/y_nnod+1
              b(1:ldm, K0) = wk2(1:ldm, l)
           end if
         end do
#endif
!$OMP END PARALLEL DO
         end if

         if ( is_size > 0 ) then
         call wait_dbl ( iq_s )
         end if

      enddo

      end do
      end do


      deallocate ( l_ )


      return
      end subroutine dc_redist1_sub

