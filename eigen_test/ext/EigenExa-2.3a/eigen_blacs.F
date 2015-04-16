       module eigen_blacs

       implicit NONE   
!--
          public :: eigen_blacs_init
          public :: eigen_blacs_exit
          public :: eigen_get_blacs_context
!--
          integer, private :: BLACS_ICONTXT_FOR_EIGEN
!--
       contains
!--
          subroutine eigen_blacs_init( TRD_COMM_WORLD,
     &                                 x_nnod, y_nnod,
     &                                 GRID_major )
          use MPI
          implicit NONE   

          integer, intent(in)       ::  TRD_COMM_WORLD
          integer, intent(in)       ::  x_nnod, y_nnod
          character*(*), intent(in) ::  GRID_major

          integer, pointer ::  tmpgrid(:,:)
          integer, pointer ::  kk0(:), kk1(:)
          integer          ::  group0, group1
          integer          ::  i, j, k, ierr


             call BLACS_GET( 0, 0, BLACS_ICONTXT_FOR_EIGEN )

             allocate( tmpgrid( 1:x_nnod, 1:y_nnod ) )
             allocate( kk0( 1:x_nnod ) )
             allocate( kk1( 1:x_nnod ) )
             call MPI_Comm_group( MPI_COMM_WORLD, group0, ierr )
             call MPI_Comm_group( TRD_COMM_WORLD, group1, ierr )

             if ( Grid_major == 'R' ) then
                do j = 1, y_nnod
                do i = 1, x_nnod
                   kk1(i) = j+(i-1)*y_nnod-1
                enddo
                   k = x_nnod
                   call MPI_Group_translate_ranks( group1, k, kk1,
     $                                             group0, kk0, ierr )
                do i = 1, x_nnod
                   tmpgrid(i, j) = kk0(i)
                end do
                end do
             else
                do j = 1, y_nnod
                do i = 1, x_nnod
                   kk1(i) = i+(j-1)*x_nnod-1
                enddo
                   k = x_nnod
                   call MPI_Group_translate_ranks( group1, k, kk1,
     $                                             group0, kk0, ierr )
                do i = 1, x_nnod
                   tmpgrid(i, j) = kk0(i)
                end do
                end do
             endif

             call BLACS_GRIDMAP( BLACS_ICONTXT_FOR_EIGEN,
     $                           tmpgrid, x_nnod, x_nnod, y_nnod )

             call MPI_Group_free( group0, ierr )
             call MPI_Group_free( group1, ierr )
             deallocate( tmpgrid )
             deallocate( kk0 )
             deallocate( kk1 )


          return
          end subroutine  eigen_blacs_init
*-
          subroutine eigen_blacs_exit( )
          implicit NONE   


             call BLACS_GRIDEXIT( BLACS_ICONTXT_FOR_EIGEN )
! commented out on ver 2.2a
!             call BLACS_EXIT( 1 )


          return
          end subroutine  eigen_blacs_exit
!--
          integer function eigen_get_blacs_context( )
          implicit NONE   


             eigen_get_blacs_context = BLACS_ICONTXT_FOR_EIGEN


          return
          end function  eigen_get_blacs_context
*-

       end module eigen_blacs

