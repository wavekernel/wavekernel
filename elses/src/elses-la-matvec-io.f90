!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_matvec_io
!
  implicit none
  integer m
  integer nz
!
  real(8), allocatable :: val(:)
  integer, allocatable :: col_ind(:)
  integer, allocatable :: row_ptr(:)
  integer, allocatable :: diag_ptr(:)
!
  real(8), allocatable :: val_s(:)
  integer, allocatable :: col_ind_s(:)
  integer, allocatable :: row_ptr_s(:)
  integer, allocatable :: diag_ptr_s(:)
!
  integer, allocatable :: flag_for_init(:)
!
  private
  public :: matvec_init
  public :: matvec_mul
  public :: get_interaction_list_number
  public :: set_interaction_list
  public :: calc_nj
  public :: inner_product_with_s
  public :: inner_product_with_h


!
  contains
!
  subroutine matvec_init(m2,nz2)
!    ---> Set all the module variables
!
    implicit none
    integer, intent(in)  :: m2, nz2
    real(8), allocatable :: dat(:)
    real(8), allocatable :: dat_s(:)
    integer              :: ierr
    integer              :: i
!
    write(*,*)'Initial setting for matrix'
    write(*,*)'  size=',m2
    write(*,*)'  Number of non-zero elements=',nz2
!
    m=m2
    nz=nz2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocation
!
    allocate (flag_for_init(1),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (val(nz),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (col_ind(nz),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (row_ptr(m+1),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (diag_ptr(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (val_s(nz),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (col_ind_s(nz),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (row_ptr_s(m+1),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (diag_ptr_s(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (dat(nz),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (dat_s(nz),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the data of matrix H
!
    open(unit=10,file='./Mat_H/ROW_PTR',status='old') 
    read(10,*) (row_ptr(i),i=1,m+1)

    open(unit=11,file='./Mat_H/COL_IND',status='old') 
    read(11,*) (col_ind(i),i=1,nz)

    open(unit=12,file='./Mat_H/VAL',status='old') 
    read(12,*) (dat(i),i=1,nz)

    open(unit=13,file='./Mat_H/DIAG_PTR',status='old') 
    read(13,*) (diag_ptr(i),i=1,m)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the data of matrix S
!
    open(unit=14,file='./Mat_S/ROW_PTR_S',status='old') 
    read(14,*) (row_ptr_s(i),i=1,m+1)

    open(unit=15,file='./Mat_S/COL_IND_S',status='old') 
    read(15,*) (col_ind_s(i),i=1,nz)

    open(unit=16,file='./Mat_S/VAL_S',status='old') 
    read(16,*) (dat_s(i),i=1,nz)

    open(unit=17,file='./Mat_S/DIAG_PTR_S',status='old') 
    read(17,*) (diag_ptr_s(i),i=1,m)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate val(1:m) and val_s(:)

    val(:) = dcmplx(dat(:),0.0d0)
    val_s(:) = dcmplx(dat_s(:),0.0d0)

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
!    
!
  end subroutine matvec_init
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct the number of the 'interaction list' 
!  : The list of the indices 'i' for the non-zero components 
!              of H(i,j_src) or S(i,j_src) 
!
  subroutine get_interaction_list_number(j_src,list_number)
    implicit none
    integer,       intent(in)   :: j_src
    integer,       intent(out)  :: list_number
    list_number=row_ptr(j_src+1)-row_ptr(j_src)
  end subroutine get_interaction_list_number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct the 'interaction list' 
!  : The list of the indices 'i' for the non-zero components 
!              of H(i,j_src) or S(i,j_src) 
!
  subroutine set_interaction_list(j_src,list_member)
!
    implicit none
    integer,       intent(in)   :: j_src
    integer,       intent(out)  :: list_member(:)
!
    integer   :: i,k
    integer   :: list_number, list_number2
    integer   :: index_in_global_list
!
    list_number=size(list_member,1)
    list_number2=row_ptr(j_src+1)-row_ptr(j_src)
!
    if (list_number /= list_number2) then
      write(*,*)'Abort(set_interaction_list):list_number=',list_number,list_number2
      stop
    endif
!
!   write(*,*)'j_src=',j_src
    do k=1,list_number
       index_in_global_list=row_ptr(j_src)+k-1
       i=col_ind(index_in_global_list)
       list_member(k)=i
!      write(*,*)'k,i=',k,i
    enddo
!
  end subroutine set_interaction_list

  subroutine calc_nj(j_src,dos,nj)
    implicit none
    real(8)      ::dos(:,:),nj(:)
    integer      ::list_number,n,ind,i,j_src
    
    list_number=size(dos,1)
    n=size(dos,2)
    write(*,*)list_number,n

    nj(:)=0.0D0

    do i=1,n
       do ind=1,list_number
          nj(i)=nj(i)+dos(ind,i)*val_s(row_ptr_s(j_src)+ind-1)
       enddo
       write(80+j_src,*)i,nj(i)
    enddo

  end subroutine calc_nj




!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with H
!   mat_kind='H' : (vect_out) = H (vect_in)
!   mat_kind='S' : (vect_out) = S (vect_in)
!
  subroutine matvec_mul(vect_in,vect_out,mat_kind)
!    ---> Set all the module variables
    implicit none
    complex(8),       intent(in)  :: vect_in(:)
    complex(8),       intent(out) :: vect_out(:)
    character(len=*),  intent(in)  :: mat_kind
    integer :: m_in, m_out
    integer :: ierr, i, j
    complex(8) :: tmp
!
    m_in=size(vect_in,1)
    m_out=size(vect_out,1)
!
    if (.not. allocated(flag_for_init) ) stop 'Abort:Not initialized'
    if (m /= m_in) stop 'Abort:Size mismatch;m,m_in'
    if (m /= m_out) stop 'Abort:Size mismatch;m,m_out'
!
    if ((mat_kind /='H') .and. (mat_kind /='S')) stop 'Abort:Wrong mode'
!
    vect_out(:)=0.0d0
!
    if (mat_kind == 'H') then
      do i=1,m
       tmp=0.0d0
       do j=row_ptr(i),row_ptr(i+1)-1
         tmp=tmp+val(j)*vect_in(col_ind(j))
       enddo
       vect_out(i)=tmp
      enddo   
    endif  
!
    if (mat_kind == 'S') then
      do i=1,m
       tmp=0.0d0
       do j=row_ptr_s(i),row_ptr_s(i+1)-1
         tmp=tmp+val_s(j)*vect_in(col_ind(j))
       enddo
       vect_out(i)=tmp
      enddo   
    endif  
!
  end subroutine matvec_mul
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inner product with the j-th column of the S matrix
!     within the 'compressed' format
!
  subroutine inner_product_with_s(vect_in,result_value,j_src)
!    ---> Set all the module variables
    implicit none
    real(8),       intent(in)  :: vect_in(:)
    real(8),       intent(out) :: result_value
    integer,          intent(in)  :: j_src
    integer                       :: i, i_run, m_in, m_s
    real(8) :: tmp
!
    m_in=size(vect_in,1)
    m_s=row_ptr_s(j_src+1)-row_ptr_s(j_src)
!
    if (.not. allocated(flag_for_init) ) stop 'Abort:Not initialized'
    if (m_in /= m_s) stop 'Abort:Size mismatch;m_in,m_s'
!
    result_value=0.0d0
!
    tmp=0.0d0
    do i_run=0,m_s-1
      i=row_ptr(j_src)+i_run
      tmp=tmp+val_s(i)*vect_in(i_run+1)
    enddo
    result_value=tmp
!
  end subroutine inner_product_with_s

  subroutine inner_product_with_h(vect_in,result_value,j_src)
!    ---> Set all the module variables
    implicit none
    real(8),       intent(in)  :: vect_in(:)
    real(8),       intent(out) :: result_value
    integer,          intent(in)  :: j_src
    integer                       :: i, i_run, m_in, m_h
    real(8) :: tmp
!
    m_in=size(vect_in,1)
    m_h=row_ptr(j_src+1)-row_ptr(j_src)
!
    if (.not. allocated(flag_for_init) ) stop 'Abort:Not initialized'
    if (m_in /= m_h) stop 'Abort:Size mismatch;m_in,m_s'
!
    result_value=0.0d0
!
    tmp=0.0d0
    do i_run=0,m_h-1
      i=row_ptr(j_src)+i_run
      tmp=tmp+val(i)*vect_in(i_run+1)
    enddo
    result_value=tmp
!
  end subroutine inner_product_with_h

end module M_la_matvec_io


