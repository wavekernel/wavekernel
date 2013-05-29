module M_create_dense_matrix
  implicit none
  private 
  public :: create_dense_matrix
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine create_dense_matrix(verbose_level, matrix_size, num_non_zeros, mat_value, mat_suffix, mat, mat_bak, eig)
    implicit none
    integer,             intent(in) :: verbose_level
    integer,             intent(in) :: matrix_size
    integer,             intent(in) :: num_non_zeros
    real(kind(1.d0)), allocatable   :: mat_value(:,:)
    integer,          allocatable   :: mat_suffix(:,:)
    real(kind(1.d0)), allocatable   :: mat(:,:)
    real(kind(1.d0)), allocatable   :: mat_bak(:,:)
    real(kind(1.d0)), allocatable   :: eig(:)
    integer,          allocatable   :: mat_chk(:,:)
!
    integer :: k, i, j, n
    integer :: ierr
!
    logical :: debug_mode
!
    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif   
!
    n=matrix_size
!
    if (debug_mode) write(*,*)'@@ create_dense_matrix'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Allocate matrices
!
    allocate(mat(n,n),stat=ierr)
    if (ierr /= 0) stop 'Stop:error in alloc. for mat'
    mat(:,:)=0.0d0
!                                                                                                                                   
    allocate(mat_bak(n,n),stat=ierr)
    if (ierr /= 0) stop 'Stop:error in alloc. for mat'
    mat_bak(:,:)=0.0d0
!                                                                                                                                   
    allocate(eig(n),stat=ierr)
    if (ierr /= 0) stop 'Stop:error in alloc. for eig'
    eig(:)=0.0d0
!                                                                                                                                   
    allocate(mat_chk(n,n),stat=ierr)
    if (ierr /= 0) stop 'Stop:error in alloc. for mat'
    mat_chk(:,:)=0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the element of the dense matrix
!
    do k=1, num_non_zeros
      i=mat_suffix(1,k) 
      j=mat_suffix(2,k) 
      mat(i,j)=mat_value(k,1)
      mat_chk(i,j)=mat_chk(i,j)+1
      if (i /= j) then
        mat(j,i)=mat_value(k,1)
        mat_chk(j,i)=mat_chk(j,i)+1
      endif   
    enddo   
!
    mat_bak(:,:)=mat(:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Check the multiple record 
!
    do j=1,n
      do i=1,n
        if (mat_chk(i,j) > 1) then 
          write(*,*)'ERROR(create_dense_matrix): multiple record : i,j =', i,j
          stop
        endif   
      enddo   
    enddo   
    if (verbose_level >= 1)  then
      write(*,*)' No multiple record in mat_value ... OK!'
    endif  
!
!
  end subroutine create_dense_matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_create_dense_matrix



