      module  eigen_devel
!$    use OMP_LIB
      use MPI, only : MPI_COMM_WORLD, MPI_THREAD_SINGLE
      implicit NONE


         public ! default is public

!  CONSTANT PARAMETERS
         real(8), parameter     :: ZERO   =  0.0D+00
         real(8), parameter     :: HALF   =  0.5D+00
         real(8), parameter     :: ONE    =  1.0D+00
         real(8), parameter     :: TWO    =  2.0D+00
         real(8), parameter     :: THREE  =  3.0D+00
         real(8), parameter     :: FOUR   =  4.0D+00
         real(8), parameter     :: FIVE   =  5.0D+00
         real(8), parameter     :: SIX    =  6.0D+00
         real(8), parameter     :: SEVEN  =  7.0D+00
         real(8), parameter     :: EIGHT  =  8.0D+00
         real(8), parameter     :: NINE   =  9.0D+00
         real(8), parameter     :: TEN    =  1.0D+01
!  PREFIX MINUS_ MEANS THE NEGATIVE SIGN
         real(8), parameter     :: MHALF  = -0.5D+00
         real(8), parameter     :: MONE   = -1.0D+00
         real(8), parameter     :: MTWO   = -2.0D+00


!$       integer                :: MPI_THREAD_MODE = MPI_THREAD_SINGLE
         character*(8), private :: Process_Grid_Major
         logical, private       :: Eigen_initialized_flag = .FALSE.

         integer                :: TRD_inod       = 0
         integer                :: TRD_nnod       = 0
         integer                :: TRD_COMM_WORLD = MPI_COMM_WORLD

         integer                :: x_inod, x_nnod, x_COMM_WORLD,
     &                             y_inod, y_nnod, y_COMM_WORLD,
     &                             z_inod, z_nnod, z_COMM_WORLD,
     &                             w_inod, w_nnod, w_COMM_WORLD,
     &                             n_common, diag_0, diag_1

         real(8)                :: Barrier_Overhead_x
         real(8)                :: Barrier_Overhead_y
         real(8)                :: Reduce_Overhead_x
         real(8)                :: Reduce_Overhead_y
         real(8)                :: Bcast_Overhead_x
         real(8)                :: Bcast_Overhead_y

         real(8), pointer       :: tmp_save(:)

         integer, pointer       :: p0_(:), q0_(:)

         integer, parameter     :: nsx = 480
         integer, parameter     :: nsm = 256 ! 64 +32
         integer, parameter     :: ns0 = nsm * nsm + 6
         integer, parameter     :: MBAND = 2

!$       integer(8) :: TRBK_MASK_FULL = 0

         real(8), pointer       :: rt_timer(:,:)
         integer                :: rt_ptr(0:16)

         real(8)                :: time_bcast, time_reduce,
     &                             time_redist, time_gather
         integer                :: items_bcast, items_reduce,
     &                             items_redist, items_gather
         real(8)                :: time_bcast_ (1:16),
     &                             time_reduce_(1:16),
     &                             time_redist_(1:16),
     &                             time_gather_(1:16)


      contains


         real(8) function eigen_get_wtime( ) result(ret)
!$       use OMP_LIB
         use MPI
         implicit NONE

#ifdef _OPENMP
!$          ret = omp_get_wtime( )
#else
            ret = MPI_Wtime( )
#endif

         return
         end function  eigen_get_wtime

         subroutine eigen_set_initialized( )

            Eigen_initialized_flag = .TRUE.

         return
         end subroutine eigen_set_initialized

         subroutine eigen_unset_initialized( )

            Eigen_initialized_flag = .FALSE.

         return
         end subroutine eigen_unset_initialized

         subroutine eigen_get_initialized( flag )
         logical :: flag

            flag = Eigen_initialized_flag

         return
         end subroutine eigen_get_initialized

         subroutine eigen_set_grid_major( Major )
         character*(*) Major

            Process_Grid_Major(1:1) = Major(1:1)

         return
         end subroutine eigen_set_grid_major

         subroutine eigen_get_grid_major( Major )
         character*(*) Major

            Major(1:1) = Process_Grid_Major(1:1)

         return
         end subroutine eigen_get_grid_major

         subroutine  eigen_timer_reset( num_bcast, num_reduce,
     &                                  num_redist, num_gather)
         integer, intent(in), optional  ::  num_bcast, num_reduce
         integer, intent(in), optional  ::  num_redist, num_gather


            time_bcast  = ZERO
            time_reduce = ZERO
            time_redist = ZERO
            time_gather = ZERO

            if ( PRESENT(num_bcast) ) then
               items_bcast  = num_bcast
            else
               items_bcast   = 0
            end if

            if ( PRESENT(num_reduce) ) then
               items_reduce = num_reduce
            else
               items_reduce  = 0
            end if

            if ( PRESENT(num_redist) ) then
               items_redist = num_redist
            else
               items_redist  = 0
            end if

            if ( PRESENT(num_gather) ) then
               items_gather = num_gather
            else
               items_gather  = 0
            end if

            time_bcast_ (1:items_bcast ) = ZERO
            time_reduce_(1:items_reduce) = ZERO
            time_redist_(1:items_redist) = ZERO
            time_gather_(1:items_gather) = ZERO

 
         return
         end subroutine eigen_timer_reset

         subroutine eigen_timer_print( messages )
         character*(*) :: messages
         real(8)       :: total_time

#if TIMER_PRINT
            if ( TRD_inod == 1 ) then

               print*, "COMM_STAT / [ ",  messages, " ]"

               print*, "   BCAST  :: ", time_bcast," // ",
     &                                  time_bcast_(1:items_bcast),
     &                            " ;; ", time_bcast_(
     &                   items_bcast+1:items_bcast+items_redist)

               print*, "   REDUCE :: ", time_reduce," // ",
     &                                  time_reduce_(1:items_reduce)

               print*, "   REDIST :: ", time_redist," // ",
     &                                  time_redist_(1:items_redist)

               print*, "   GATHER :: ", time_gather," // ",
     &                                  time_gather_(1:items_gather),
     &                            " ;; ", time_gather_(
     &                   items_gather+1:items_gather+items_redist)

               total_time = time_bcast + time_reduce
     &                                 + time_redist
     &                                 + time_gather
               print*, "   Total :: ", total_time

            end if
#endif
               items_bcast   = 0
               items_reduce  = 0
               items_redist  = 0
               items_gather  = 0


         return
         end subroutine  eigen_timer_print

!$       subroutine reduction_sum( x )
!$       implicit NONE
!$       real(8), intent(INOUT) :: x
!$       integer :: j
!$
!$          if ( omp_get_num_threads() > 1 ) then
!$OMP BARRIER
!$             tmp_save(omp_get_thread_num()) = x
!$OMP BARRIER
!$OMP MASTER
!$             x = ZERO
!$             do j=0,omp_get_num_threads()-1
!$                x = x + tmp_save(j)
!$             end do
!$             tmp_save(0) = x
!$OMP END MASTER
!$OMP BARRIER
!$             x = tmp_save(0)
!$OMP BARRIER
!$          end if
!$
!$       return
!$       end subroutine reduction_sum

!$       subroutine sync_other_than_master_init( TRBK_lock, TRBK_mask )
!$       use OMP_LIB
!$       use MPI
!$       integer(kind=omp_lock_kind) :: TRBK_lock
!$       integer(8)                  :: TRBK_mask(1:2)
!$
!$          integer :: i
!$
!$          TRBK_MASK_FULL = 0
!$          if ( omp_get_num_threads() > 64 ) then
!$             TRBK_MASK_FULL = omp_get_num_threads() - 1
!$          else
!$             do i=1,omp_get_num_threads()-1
!$                TRBK_MASK_FULL = ibset( TRBK_MASK_FULL, i )
!$             end do
!$          end if
!$
!$          call omp_init_lock( TRBK_lock )
!$
!$          call omp_set_lock( TRBK_lock )
!$          TRBK_mask(1) = 0
!$          TRBK_mask(2) = 0
!$          call omp_unset_lock( TRBK_lock )
!$
!$       end subroutine  sync_other_than_master_init

!$       subroutine sync_other_than_master( TRBK_lock, TRBK_mask )
!$       use OMP_LIB
!$       integer(kind=omp_lock_kind) :: TRBK_lock
!$       integer(8)                  :: TRBK_mask(1:2)
!$
!$       integer(8) :: T
!$       integer    :: local_size, local_rank
!$
!$          local_size = omp_get_num_threads()
!$          local_rank = omp_get_thread_num()
!$
!$          if ( local_size == 1 ) return
!$          if ( local_rank == 0 ) return
!$
!$          call omp_set_lock( TRBK_lock )
!$          if ( local_size > 64 ) then
!$             T = TRBK_mask(1) + 1
!$          else
!$             T = ibset( TRBK_mask(1), local_rank )
!$          end if
!$          if ( T == TRBK_MASK_FULL ) TRBK_mask(2) = T
!$          TRBK_mask(1) = T
!$          call omp_unset_lock( TRBK_lock )
!$
!$          do
!$             call omp_set_lock( TRBK_lock )
!$             T = TRBK_mask(1)
!$             call omp_unset_lock( TRBK_lock )
!$             if ( T == TRBK_MASK_FULL ) exit
!$          end do
!$
!$          call omp_set_lock( TRBK_lock )
!$          if ( local_size > 64 ) then
!$             T = TRBK_mask(2) - 1
!$          else
!$             T = ibclr( TRBK_mask(2), local_rank )
!$          end if
!$          if ( T == 0 ) TRBK_mask(1) = T
!$          TRBK_mask(2) = T
!$          call omp_unset_lock( TRBK_lock )
!$
!$          do
!$             call omp_set_lock( TRBK_lock )
!$             T = TRBK_mask(2)
!$             call omp_unset_lock( TRBK_lock )
!$             if ( T == 0 ) exit
!$          end do
!$
!$       end subroutine  sync_other_than_master

      end module  eigen_devel

