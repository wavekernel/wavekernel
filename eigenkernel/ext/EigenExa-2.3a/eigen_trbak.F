      module eigen_TRBAK
!$    use OMP_LIB
      implicit NONE
!=
         real(8) :: TRBK_TIME_BCAST
         real(8) :: TRBK_TIME_REDUC
         real(8) :: TRBK_TIME_FR
         real(8) :: TRBK_TIME_TRBK1
         real(8) :: TRBK_TIME_TRBK1_
         real(8) :: TRBK_TIME_TRBK1x
         real(8) :: TRBK_TIME_TRBK1x_
         real(8) :: TRBK_TIME_TRBK2

         real(8) :: TRBK_TIME_REDUC_OVERHEAD_x

!$       integer(kind=omp_lock_kind) :: TRBK_LOCK
!$       integer(8)                  :: TRBK_MASK(1:2)

         integer :: DO_OVERLAP_BCAST_LEVEL
         integer :: TRBK_TIME_COUNTER
         integer :: TRBK_TIME_INTERVAL
         integer :: TRBK_TIME_NEXT
         integer :: TRBK_SWITCHED
!=
         real(8), pointer :: TRBK_buf(:)
!=
      end module eigen_TRBAK
