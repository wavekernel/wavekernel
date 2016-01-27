!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_eig_standard
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
  public :: eig_standard
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eig_standard(v_mat, eig, a_mat, m_ini, m_fin, msg)
!
!     Standard eigen-value problem : A v_k = e_k v_k
!          for k:= m_ini, m_ini+1, .... m_fin
!
!        NOTE: The parameters m_ini and m_fin are ignored tentatively.
!
!          
!
    implicit none
    real(8),          intent(out)   :: v_mat(:,:)      ! eigen vector
    real(8),          intent(out)   :: eig(:)          ! eigen value
    real(8),          intent(in)    :: a_mat(:,:)      ! input matrix
    integer,          intent(in)    :: m_ini,m_fin     ! range of the target levels
    character(len=*), intent(inout) :: msg             ! message (optional)
!
    integer :: n,ierr
    real(8), allocatable ::  a_wrk(:,:)
    real(8), allocatable ::  wrk_vec(:)
!
    integer :: info, lda, lwork
    real(8), allocatable ::  work(:)
    integer :: i_check_eigen, i,j
    real(8) ::ddd
    character(len=72) :: msg_wrk
!
    if (trim(msg) == 'debug') then
      i_check_eigen=1
    else 
      i_check_eigen=0
    endif 
!
    n=size(a_mat, 1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Checking the input parameters
!
    if (n < 1) then
      write(*,*)'ERROR(eig_standard):matrix size=', n
      stop
    endif   
!
    if (len(msg) < 1) then
      write(*,*)'ERROR(eig_standard):len(msg)=',len(msg)
      stop
    endif
    msg=''
!
    ierr=0
    if (size(v_mat, 1) /= n) ierr=1 
    if (size(v_mat, 2) /= n) ierr=1 
    if (size(eig,   1) /= n) ierr=1 
    if (size(a_mat, 1) /= n) ierr=1 
    if (size(a_mat, 2) /= n) ierr=1 
!
    if (ierr /=0) then
      write(msg_wrk,'(a,i20)') 'Unmatched matrix size:n=', n
      if (len(msg) < len_trim(msg_wrk)) then
        msg='E' 
      else
        msg=trim(msg_wrk)
      endif   
      write(*,*) 'ERROR(eig_standard):n=', n
      write(*,*) 'ERROR(eig_standard) size(v_mat)=', size(v_mat,1), size(v_mat,2)
      write(*,*) 'ERROR(eig_standard) size(eig)  =', size(eig,1)
      write(*,*) 'ERROR(eig_standard) size(a_mat)=', size(a_mat,1), size(a_mat,2)
      return
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    lda=n
    lwork=3*n-1
    info=0
!
    allocate (a_wrk(n,n),stat=ierr)
    if (ierr /= 0) stop 'Alloc error (a_wrk):eig_standard'
    allocate (work(3*n-1),stat=ierr)
    if (ierr /= 0) stop 'Alloc error (work):eig_standard'
!
    a_wrk(:,:)=a_mat(:,:)
    eig(:)=0.0d0
!
    call dsyev("V","U",n,a_wrk,lda,eig,work,lwork,info)
    if (info /= 0) then
      write(*,*)'LAPACK ERROR:info=',info
      msg='LAPACK ERROR'
      return
    endif   
!
    v_mat(:,:)=a_wrk(:,:)
!
    if (i_check_eigen == 1) then
      allocate (wrk_vec(n),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc'
      do j=1,n
        wrk_vec(:)=0.0d0
        wrk_vec(:)=matmul(a_mat(:,:), v_mat(:,j))
        wrk_vec(:)=wrk_vec(:)-eig(j)*v_mat(:,j)
        ddd=dot_product(wrk_vec(:),wrk_vec(:))
        if (dabs(ddd) >= 1.0d-10) then
          write(*,*)'ERROR(calc_eig): j, residual =',j, ddd
          stop
        endif
!       write(*,*)' residual =',j, ddd
      enddo   
      deallocate (wrk_vec,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
    endif
!
    if (i_check_eigen == 1) then
      do j=1,n
       do i=1,n
         ddd=dot_product(v_mat(i,:),v_mat(j,:))
         if (i == j) ddd=ddd-1.0d0
         if (dabs(ddd) >= 1.0d-10) then
           write(*,*)'ERROR(calc_eig): j, residual in completeness =',i,j, ddd
           stop
         endif
       enddo
      enddo   
    endif
!
    deallocate (a_wrk,stat=ierr)
    if (ierr /= 0) stop 'Dealloc error (a_wrk):eig_standard'
!
    deallocate (work,stat=ierr)
    if (ierr /= 0) stop 'Dealloc error (work) :eig_standard'
!
    if (allocated(wrk_vec)) then
      deallocate (wrk_vec,stat=ierr)
      if (ierr /= 0) stop 'Dealloc error (wrk_vec) :eig_standard'
    endif   
!
  end subroutine eig_standard
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_eig_standard
