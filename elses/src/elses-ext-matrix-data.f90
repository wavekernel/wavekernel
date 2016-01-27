module M_ext_matrix_data
!
  use M_qm_domain ,  only : i_verbose    !(unchanged)
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  type :: sparce_matrix_data_real_type
     integer                      :: matrix_size
     integer                      :: num_of_non_zero_elements
     integer,         allocatable :: element_index(:,:)
     real(DOUBLE_PRECISION), allocatable :: element_data(:)
  end type sparce_matrix_data_real_type
!
  type(sparce_matrix_data_real_type) :: matrix_data(2)
!            matrix_data(1) : H
!            matrix_data(2) : S
!
  private
!
  public sparce_matrix_data_real_type
  public matrix_data
  public set_matrix_data
  public plot_matrix_data
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ test routine for quantum dynamics
!
  subroutine test_for_quantum_dynamics
!
    if (i_verbose >= 1) then
      write(*,*)'@@ test_for_quantum_dynamics'
    endif
!
    call set_matrix_data
    call plot_matrix_data
!
  end subroutine test_for_quantum_dynamics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set matrix data array for external library
!
  subroutine set_matrix_data
!
    use M_config,            only : config !(unchanged)
    use elses_mod_orb2,      only : js2j   !(unchanged)
    use elses_mod_js4jsv,    only : js4jsv !(unchanged)
!
!   use elses_mod_file_io, only : vacant_unit !(function)
    use M_qm_domain ,  only : dhij, dsij, atm_element, nval, &
&                             jsv4jsd, njsd, noav, atm_element !(unchanged)
!
    implicit none
!
    real(DOUBLE_PRECISION) :: min_value
    character(len=32) :: format_mode
!
    integer, parameter  :: ict4h=1
!
    integer      :: jj, prc_index, ierr
    integer      :: js2, jsv2, nss2, nval2, ja2, jj2, jsd1
    integer      :: js1, jsv1, nss1, nval1, ja1, jj1
    integer      :: mat_size, num_of_non_zero
    real(DOUBLE_PRECISION) :: data_values(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control parameters
!
    format_mode='sym'
    min_value = 1.0d-16
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    if (i_verbose >= 1) then
      write(*,*)'@@ set_matrix_data'
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Trivial Checking
!
    if (.not. allocated(dhij)) then
      write(*,*)'ERROR(set_matrix_data):dhij is not allocated'
      stop
    endif
!
    if (.not. allocated(dsij)) then
      write(*,*)'ERROR(set_matrix_data):dhij is not allocated'
      stop
    endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count up the matrix size
!
    jj=0
    do jsv2=1,noav
      jj=jj+nval(atm_element(jsv2))
    enddo
    mat_size=jj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main loop
!   prc_index = 1 : Count up the number of non-zero elements
!             = 2 : Store the matrix into the array
!
   do prc_index=1,2
     jj=0
     do jsv2=1,noav
       nss2=atm_element(jsv2)
       nval2=nval(nss2)
       js2=js4jsv(jsv2)
       do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         js1=js4jsv(jsv1)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         js1=js4jsv(jsv1)
         do ja2=1,nval2
           jj2=js2j(ja2,js2)
           do ja1=1,nval1
             jj1=js2j(ja1,js1)
             data_values(1)=dhij(ja1,ja2,jsd1,jsv2)
             data_values(2)=dsij(ja1,ja2,jsd1,jsv2)
             if ( dabs(data_values(2)) < min_value ) cycle
             if ( trim(format_mode) == "sym" ) then
               if (jj1 < jj2) cycle
             endif
             jj=jj+1
             if (prc_index == 2) then
               call set_matrix_data_element(jj, jj1, jj2, data_values)
             endif
          enddo
        enddo
       enddo
     enddo
     if (prc_index == 1) then
       num_of_non_zero=jj
       call allocate_matrix_data(mat_size, num_of_non_zero)
     endif
   enddo
!
  end subroutine set_matrix_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Allocate matrix data array
!
  subroutine allocate_matrix_data(m,nnz)
    implicit none
    integer, intent(in) :: nnz, m
    integer             :: k, ierr
!
    do k=1,2
      matrix_data(k)%matrix_size=m
      matrix_data(k)%num_of_non_zero_elements=nnz
      allocate(matrix_data(k)%element_index(2,nnz), stat=ierr )
      allocate(matrix_data(k)%element_data(nnz),    stat=ierr )
      if (ierr /= 0) then
        write(*,*) 'Alloc. error (allocate_matrix_data) :matrix_data '
        stop
      endif
    enddo
!
  end subroutine allocate_matrix_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Set matrix data element
!
  subroutine set_matrix_data_element(k, i, j, data_values)
    implicit none
    integer, intent(in)                :: i, j, k
    real(DOUBLE_PRECISION), intent(in) :: data_values(2)
    integer                            :: mat_kind
    integer                            :: m, nnz, ierr
    logical,  parameter                :: debug_mode = .true.
!
    if (debug_mode) then
     do mat_kind=1,2
       m   = matrix_data(mat_kind)%matrix_size
       nnz = matrix_data(mat_kind)%num_of_non_zero_elements
       ierr=0
       if ((i < 1) .or. (i > m)) ierr=1
       if ((j < 1) .or. (j > m)) ierr=1
       if ((k < 1) .or. (k > nnz)) ierr=1
       if (ierr == 1) then
         write(*,*) 'ERROR(set_matrix_data_element):k,i,j=',k,i,j
         stop
       endif
     enddo
    endif
!
    do mat_kind=1,2
      matrix_data(mat_kind)%element_index(1,k)=i
      matrix_data(mat_kind)%element_index(2,k)=j
      matrix_data(mat_kind)%element_data(k)=data_values(mat_kind)
    enddo
!
  end subroutine set_matrix_data_element
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Plot matrix data element
!
  subroutine plot_matrix_data
    use elses_mod_file_io,  only : vacant_unit  !(function)
    implicit none
    character(len=32)      :: filename_wrk
    integer                :: i, j, k
    real(DOUBLE_PRECISION) :: data_value
    integer                :: mat_kind
    integer                :: m, nnz, ierr
    integer                :: iunit
    logical,  parameter    :: debug_mode = .true.
!
    iunit=vacant_unit()
!
    do mat_kind=1,2
      m   = matrix_data(mat_kind)%matrix_size
      nnz = matrix_data(mat_kind)%num_of_non_zero_elements
      if (mat_kind == 1) filename_wrk='matrix_H.txt'
      if (mat_kind == 2) filename_wrk='matrix_S.txt'
      open(iunit, file=filename_wrk, status='unknown')
      write(iunit,'(a)')'%%MatrixMarket matrix coordinate real symmetric'
      write(iunit,*) m, m, nnz
      do k=1,nnz
        i=matrix_data(mat_kind)%element_index(1,k)
        j=matrix_data(mat_kind)%element_index(2,k)
        data_value = matrix_data(mat_kind)%element_data(k)
        write(iunit,'(2i10,f30.20)') i, j, data_value
      enddo
      close(iunit)
    enddo
!
  end subroutine plot_matrix_data
!
end module M_ext_matrix_data
