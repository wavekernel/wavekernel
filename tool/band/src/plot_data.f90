!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module m_plot_data
!
  implicit none
  integer,  parameter :: DOUBLE_PRECISION=kind(1d0)
!
  private
  public  plot_and_mask_matrix
  public  plot_eigen_value
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine plot_and_mask_matrix(n,H,S)
    implicit none
    integer,                   intent(in) :: n
    complex(DOUBLE_PRECISION), intent(inout) :: H(n,n)   
    complex(DOUBLE_PRECISION), intent(inout) :: S(n,n)   
!
    complex(DOUBLE_PRECISION) :: H_bak(n,n)   
    complex(DOUBLE_PRECISION) :: S_bak(n,n)   
    integer,      allocatable :: member_i(:)
    integer,      allocatable :: member_j(:)
!
    real(DOUBLE_PRECISION), parameter     :: value_mask=1.0d-15 ! (input) 
!
    integer :: i,j,k, n_count, i_mat
    integer :: ierr, iunit
!
    complex(DOUBLE_PRECISION) :: mat_value
    real(DOUBLE_PRECISION)    :: abs_mat_value
!
    iunit=30
!
    write(*,*)'n=',n
!
    H_bak(:,:)=H(:,:)
    S_bak(:,:)=S(:,:)
    H(:,:)=0.0d0
    S(:,:)=0.0d0
!
    allocate (member_i(n*n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error:member_i'
      stop
    endif   
!
    allocate (member_j(n*n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'Alloc. error:member_j'
      stop
    endif   
!    
!   do j=1,n
!     do i=1,n   
!       write(*,*)'H,S=',i,j,H(i,j), S(i,j)
!     enddo
!   enddo
!  
    do i_mat=1,2
!
      n_count=0 
      do i=1,n
        do j=i,n   
          if (i_mat == 1) then
            mat_value=H_bak(i,j) 
          else   
            mat_value=S_bak(i,j) 
          endif
          abs_mat_value=abs(mat_value) 
          if (i /= j) then
            if ( abs_mat_value < value_mask ) cycle
          endif   
          n_count=n_count+1
          member_i(n_count)=i
          member_j(n_count)=j
        enddo
      enddo
!
      write(iunit,'(a)') '%%MatrixMarket matrix coordinate complex hermitian'
      if (i_mat == 1) then
        write(iunit,'(a)') &
&        '%Matrix data "BNZ30" A, from ELSES Matrix Library (http://www.elses.jp/matrix), T. Hoshi, 30. Apr. 2013'
      else  
        write(iunit,'(a)') &
&        '%Matrix data "BNZ30" B, from ELSES Matrix Library (http://www.elses.jp/matrix), T. Hoshi, 30. Apr. 2013'
      endif
!  
      write(iunit,'(3i20)') n, n, n_count
!
      write(*,*)'n_count, n*(n+1)/2 =', n_count, n*(n+1)/2
      do k=1,n_count
        i=member_i(k) 
        j=member_j(k) 
        if (i_mat == 1) then
          mat_value=H_bak(i,j)
        else
          mat_value=S_bak(i,j)
        endif
        if (i_mat == 1) then
          if (i == j ) then 
            H(i,i) = real(H_bak(i,i),kind(0d0))
          else
            H(i,j)=H_bak(i,j)
            H(j,i)=conjg(H(i,j))
          endif  
          write(iunit,'(2i10,2f30.18)') i, j, real(H(i,j),kind(0d0)), dimag(H(i,j))
        else
          if (i == j ) then 
            S(i,i) = real(S_bak(i,i),kind(0d0))
          else
            S(i,j)=S_bak(i,j)
            S(j,i)=conjg(S(i,j))
          endif  
          write(iunit,'(2i10,2f30.18)') i, j, real(S(i,j),kind(0d0)), dimag(S(i,j))
        endif
      enddo
!
    enddo
!
!   stop 'Stop manually'
!
  end subroutine plot_and_mask_matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine plot_eigen_value(n,eig)
    implicit none
    integer,                 intent(in) :: n
    real(DOUBLE_PRECISION),  intent(in) :: eig(n)   
    integer :: k, iunit
!
    iunit=31
    do k=1, n
      write(iunit,'(i10,f20.12)') k, eig(k)
    enddo
!
  end subroutine plot_eigen_value
!
end module m_plot_data
