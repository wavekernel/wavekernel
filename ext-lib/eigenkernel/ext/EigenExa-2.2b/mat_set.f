       subroutine  mat_set( n, a, nm, mtype )
       use MPI
       use eigen_libs
       implicit NONE

       integer, intent(IN)    :: n, nm, mtype
       real(8)                :: a(1:nm,*)

! Parameters BLACS array descritor(the position of entry tags), etc
       INTEGER, PARAMETER     :: BLOCK_CYCLIC_2D = 1
       INTEGER, PARAMETER     :: DLEN_  = 9
       INTEGER, PARAMETER     :: DTYPE_ = 1
       INTEGER, PARAMETER     :: CTXT_  = 2
       INTEGER, PARAMETER     :: M_     = 3
       INTEGER, PARAMETER     :: N_     = 4
       INTEGER, PARAMETER     :: MB_    = 5
       INTEGER, PARAMETER     :: NB_    = 6
       INTEGER, PARAMETER     :: RSRC_  = 7
       INTEGER, PARAMETER     :: CSRC_  = 8
       INTEGER, PARAMETER     :: LLD_   = 9

       integer                :: DESCA( DLEN_ )
       integer                :: NPROW, NPCOL

       real(8), pointer       :: as(:,:), w(:)
       real(8)                :: t, s, hi_, hj_

       integer                :: COMM, x_COMM, y_COMM
       integer                :: nnod, x_nnod, y_nnod
       integer                :: inod, x_inod, y_inod
       integer                :: i, i_1, i_2, i_3
       integer                :: j, j_1, j_2, j_3
       integer                :: k
       integer                :: nx, info, ICTXT, ierror
       integer, pointer       :: iseed(:)

       external               :: DESCINIT, PDTRAN

       real(8), parameter     :: ZERO= 0.0D0
       real(8), parameter     :: ONE = 1.0D0
       real(8)                :: PAI, EPS, EPS2, EPS4, theta

       integer :: n1, n2, ne


!          if ( mtype < 0 .OR. mtype > 10 ) then
!             print*,"Illeagal Matrix Type is specified."
!             call  MPI_Abort( COMM, MPI_ERR_OTHER, ierror )
!          end if

          call  eigen_get_comm ( COMM, x_COMM, y_COMM )
          call  eigen_get_procs( nnod, x_nnod, y_nnod )
          call  eigen_get_id   ( inod, x_inod, y_inod )

          j_2 = eigen_loop_start( 1, x_nnod, x_inod )
          j_3 = eigen_loop_end  ( n, x_nnod, x_inod )
          i_2 = eigen_loop_start( 1, y_nnod, y_inod )
          i_3 = eigen_loop_end  ( n, y_nnod, y_inod )

          PAI  = 4*ATAN(1D0)
          EPS  = machine_epsilon()
          EPS2 = SQRT(EPS)
          EPS4 = SQRT(EPS2)

          if ( mtype == 0 ) then

             do i_1 = i_2, i_3
                i = eigen_translate_l2g( i_1, y_nnod, y_inod )
                do j_1 = j_2, j_3
                   j = eigen_translate_l2g( j_1, x_nnod, x_inod )
                   a(j_1, i_1) = MIN(i,j)*ONE
                end do
             end do

          end if

          if ( mtype == 1 ) then

             do i_1 = i_2, i_3
                i = eigen_translate_l2g( i_1, y_nnod, y_inod )
                do j_1 = j_2, j_3
                   j = eigen_translate_l2g( j_1, x_nnod, x_inod )
                   if ( i == j ) then
                      a(j_1, i_1) = -7.2D+00
                   else
                      a(j_1, i_1) = -3.0D+00/(i-j)**2
                   end if
                end do
             end do

          end if

          if ( mtype == 2 ) then

             NPROW = x_nnod
             NPCOL = y_nnod

             ICTXT = eigen_get_blacs_context( )

             CALL  DESCINIT( DESCA, n, n, 1, 1, 0, 0, ICTXT, nm, INFO )

             nx = (n-1)/NPCOL+1
             allocate( as(1:nm,nx) )

             CALL  RANDOM_SEED( size = i )
             allocate( iseed(i) )
             iseed(1:i) = inod
             CALL  RANDOM_SEED( put = iseed )
             deallocate( iseed )

             do i_1 = i_2, i_3
                do j_1 = j_2, j_3
                   CALL RANDOM_NUMBER( t )
                   as(j_1, i_1) = t
                   a (j_1, i_1) = t
                end do
             end do

             CALL  PDTRAN( n, n,
     &                    ONE, as, 1, 1, DESCA, ONE, a, 1, 1, DESCA )

             deallocate( as )

          end if

          if ( mtype == 3 ) then

             do i_1 = i_2, i_3
                i = eigen_translate_l2g( i_1, y_nnod, y_inod )
                do j_1 = j_2, j_3
                   j = eigen_translate_l2g( j_1, x_nnod, x_inod )
                   a(j_1, i_1) = (n+1-MAX(i,j))*ONE
                end do
             end do

          end if


          if ( 4 <= mtype .AND. mtype <= 10 ) then
             allocate ( w(1:n) )
             call  w_set( n, w, mtype )
             call  helmert_trans( n, a, nm, w, i_2, i_3, j_2, j_3 )
             deallocate ( w )

          end if

       return
       contains

       subroutine  helmert_trans( n, a, lda, w, i_2, i_3, j_2, j_3 )
       implicit NONE
       integer, intent(IN)    :: n, lda, i_2, i_3, j_2, j_3
       real(8), intent(OUT)   :: a(lda, *)
       real(8), intent(IN)    :: w(*)

       integer                :: i_1, j_1, i, j
       real(8)                :: hi_, hj_, t, s, scale
       real(8), pointer       :: w_(:), hi(:), hj(:), ht(:)
       integer, pointer       :: iseed(:)


          allocate ( w_(1:n), hi(1:n), hj(1:n), ht(1:n) )

          scale = 0D0
          do k=1,n
             scale = MAX(scale, ABS(w(k)))
          end do
          if ( scale < 1D0 ) then
             scale = 1D0
          end if
          w_(1:n) = w(1:n) / scale

          CALL RANDOM_SEED( size = i )
          allocate( iseed(i) )
          iseed(1:i) = 0
          CALL RANDOM_SEED( put = iseed )
          deallocate( iseed )
          do k=1,n
             CALL RANDOM_NUMBER( t )
             i = INT(t*n)+1
             CALL RANDOM_NUMBER( t )
             j = INT(t*n)+1
             t = w_(i); w_(i) = w_(j); w_(j) = t
          enddo

          do i=1,n
             ht(i) = SQRT(DBLE(i))
          end do

!$OMP PARALLEL DO PRIVATE(i,j,k, i_1,j_1, s,t, hi,hj, hi_,hj_)
          do i_1 = i_2, i_3
             i = eigen_translate_l2g( i_1, y_nnod, y_inod )

             if ( i==1 ) then
                hi_ = ht(n) ! SQRT(DBLE(n))
                hi(1:n)   = 1D0 / hi_
             else
                s   = DBLE(i-1)
!                hi_ = s * SQRT(1D0+1D0/s)
                hi_ = ht(i-1) * ht(i) ! SQRT(s) * SQRT(s+1D0)
                hi(1:i-1) =  1D0 / hi_
                hi(i)     = -s / hi_
                hi(i+1:n) =  0D0
             end if

             do j_1 = j_2, j_3
                j = eigen_translate_l2g( j_1, x_nnod, x_inod )

             if ( j==1 ) then
                hj_ = ht(n) ! SQRT(DBLE(n))
                hj(1:n)   = 1D0 / hj_
             else
                s   = DBLE(j-1)
!                hj_ = s * SQRT(1D0+1D0/s)
                hj_ = ht(j-1)*ht(j) ! SQRT(s) * SQRT(s+1D0)
                hj(1:j-1) =  1D0 / hj_
                hj(j)     = -s / hj_
                hj(j+1:n) =  0D0
             end if

                t = 0D0
                do k=1,n
                   t = t + w_(k) * hi(k) * hj(k)
                end do
                a(j_1, i_1) = t ! / (hi_ * hj_)

             end do
          end do
!$OMP END PARALLEL DO

          a(j_2:j_3, i_2:i_3) = a(j_2:j_3, i_2:i_3) * scale

          deallocate( w_, hi, hj, ht )


       end subroutine  helmert_trans

       end subroutine  mat_set

       subroutine  print_mat_name_list( )
       implicit NONE
       integer :: mtype

       do  mtype = 0, 10
          call  print_mat_name( mtype )
       end do

       end subroutine  print_mat_name_list

       subroutine  print_mat_name( mtype )
       implicit NONE
       integer, intent(IN)    :: mtype
       character*(100)        :: message

          message = " (unkonwn or invalid)"
          select case (mtype)
          case (0)
             message = " (Frank matrix)"
          case (1) 
             message = " (Toeplitz matrix)"
          case (2)
             message = " (Random matrix)"
          case (3)
             message = " (Frank matrix 2)"
          case (4)
             message = " (W: 0, 1, ..., n-1)"
          case (5)
             message = " (W: sin(PAI*5*i/(n-1)+EPS^1/4)^3)"
          case (6)
             message = " (W: MOD(i,5)+MOD(i,2))"
          case (7)
             message = " (W: same as Frank matrix)"
          case (8)
             message = " (W: Uniform Distribution, [0,1))"
          case (9)
             message = " (W: Gauss Distribution, m=0,s=1)"
          case (10)
             message = " (W: Read from the data file 'W.dat')"
          end select

          print*, "Matrix type = ", mtype, trim(message)

       end subroutine  print_mat_name

       subroutine  w_set( n, w, mtype )
       use MPI
       use eigen_libs
       implicit NONE

       integer, intent(IN)    :: n, mtype
       real(8)                :: w(1:n)

       real(8), parameter     :: ZERO = 0.0D0
       real(8), parameter     :: ONE  = 1.0D0

       real(8)                :: PAI, EPS, EPS2, EPS4
       real(8)                :: ax, bx, x, y, z, theta, s, t
       real(8), pointer       :: ww(:)
       integer, pointer       :: iseed(:)

       integer                :: i, j, k, kk, l


          if ( mtype == 1 .OR. mtype == 2 ) then
          return
          end if

          PAI  = 4*ATAN(ONE)
          EPS  = machine_epsilon( )
          EPS2 = SQRT(EPS)
          EPS4 = SQRT(EPS2)

          if ( mtype == 0 .OR. mtype == 3 .OR. mtype == 7 ) then
             do i = 1, n
                j = n-i
                theta = PAI*(2*j+1)/(2*n+1)
                w(i) = 5D-1/(1D0-COS(theta))
             end do
          end if

          if ( mtype == 4 ) then
             do i = 1, n
                w(i) = DBLE(i-1)
             end do
          end if

          if ( mtype == 5 ) then
             do i = 1, n
                theta = PAI*5*i/(n-1)+EPS4
                w(i) = SIN(theta)**3
             end do
          end if

          if ( mtype == 6 ) then
             do i = 1, n
                w(i) = DBLE(MOD(i, 5)+MOD(i, 2))
             end do
          end if

          if ( mtype == 8 ) then
             CALL RANDOM_SEED( size = i )
             allocate( iseed(i) )
             iseed(1:i) = 0
             CALL RANDOM_SEED( put = iseed )
             deallocate( iseed )

             do i = 1, n
                CALL RANDOM_NUMBER( t )
                w(i) = t
             end do
          end if

          if ( mtype == 9 ) then
             CALL RANDOM_SEED( size = i )
             allocate( iseed(i) )
             iseed(1:i) = 0
             CALL RANDOM_SEED( put = iseed )
             deallocate( iseed )

             do i = 1, n
                CALL RANDOM_NUMBER( t )
                CALL RANDOM_NUMBER( s )
                w(i) = SQRT(-2*LOG(s)) * SIN(2*PAI*s)
             end do
          end if

          if ( mtype == 10 ) then
             open(11, FILE='W.dat')
             read(11, *) w(1:n)
             close(11)
          end if


       return
       end subroutine  w_set

