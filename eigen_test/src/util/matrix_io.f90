module matrix_io
  implicit none
  private
  public :: read_matrix_file, print_matrix
!
contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_matrix_file(verbose_level, mtx_filename, matrix_type, mat_size, num_non_zeros, mat_value, mat_suffix)
    implicit none
    integer,             intent(in) :: verbose_level
    character(len=256),  intent(in) :: mtx_filename
    character(len=256),  intent(in) :: matrix_type
    integer,             intent(out) :: mat_size
    integer,             intent(out) :: num_non_zeros
    real(kind(1.d0)), allocatable    :: mat_value(:,:)
    integer,          allocatable    :: mat_suffix(:,:)
    integer,             parameter :: unit_num=1
    integer :: ierr
    logical :: debug_mode
!
    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif
!
    if (debug_mode) write(*,'(a)')'@@ read_matrix_file'
!
    if (debug_mode) write(*,'(a,a)')'  filename=', trim(mtx_filename)
!
    call check_file_name(verbose_level, mtx_filename)
!
    open(unit_num,file=mtx_filename)
!
    call read_matrix_file_header(verbose_level, unit_num, matrix_type, mat_size, num_non_zeros)
!
    allocate (mat_suffix(2, num_non_zeros),  stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR in allocation : mat_suffix'
      stop
    endif
!
    if (trim(matrix_type) == 'real_symmetric') then
      allocate (mat_value(num_non_zeros,1),  stat=ierr)
      if (ierr /= 0) then
        write(*,*)'ERROR in allocation : mat_value'
        stop
      endif
    else
      write(*,*) 'ERROR:unsuported matrix type = ',trim(matrix_type)
    endif
!
    call read_matrix_file_value(verbose_level, unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
!
    close(unit_num)
!
  end subroutine read_matrix_file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_file_name(verbose_level, mtx_filename)
    implicit none
    integer,            intent(in) :: verbose_level
    character(len=256), intent(in) :: mtx_filename
!
    logical :: debug_mode
!
    logical :: file_exist
    integer :: k
!
    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif
!
    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_header'
!
    if (debug_mode) write(*,'(a,a)')'  filename=', trim(mtx_filename)
    if (debug_mode) write(*,*)'verbose level =', verbose_level
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Checking the filename
!
    inquire (file=trim(mtx_filename), exist=file_exist)
    if (.not. file_exist) then
      write(*,*)' ERROR:No matrix file: filename=',trim(mtx_filename)
      stop
    else
      if (debug_mode) write(*,*)' INFO: Matrix file founded : filename=',trim(mtx_filename)
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Check the extension
!
    k=len_trim(mtx_filename)
    if ( mtx_filename(k-3:k) /= '.mtx' ) then
      write(*,*)'ERROR:Wrong file : filename=', trim(mtx_filename)
      stop
    else
      if (debug_mode) write(*,*)'The extension of the is checked: filename=',trim(mtx_filename)
    endif
!
  end subroutine check_file_name
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_matrix_file_header(verbose_level, unit_num, matrix_type, mat_size, num_non_zeros)
    implicit none
    integer,            intent(in)  :: verbose_level
    integer,            intent(in)  :: unit_num
    character(len=256), intent(in)  :: matrix_type
    integer,            intent(out) :: mat_size
    integer,            intent(out) :: num_non_zeros
!
    integer,            parameter  :: max_line_count = 100
!
    logical :: keyword_real_exist
    logical :: keyword_symmetric_exist
    integer             :: line_count
    integer             :: ierr
    integer             :: mat_size2
!
    logical :: debug_mode
!
    character(len=1024) :: chara_wrk
!
    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif
!
    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_header'
!
    if (debug_mode) write(*,*)'verbose level =', verbose_level
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Read the first line
!
    read(unit_num,'(a)') chara_wrk
    if (debug_mode) write(*,*)'first line=',trim(chara_wrk)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Check the header of the first line
!
    if (chara_wrk(1:14) /= '%%MatrixMarket') then
      write(*,*)'ERROR:Fist line :', trim(chara_wrk)
      stop
    else
      if (debug_mode) write(*,*)'INFO:header of file =', chara_wrk(1:14)
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Seach the keyword of 'real' on the first line
!
    keyword_real_exist = .false.
    if ( index(chara_wrk, 'real') /= 0 ) keyword_real_exist = .true.
    if ( index(chara_wrk, 'Real') /= 0 ) keyword_real_exist = .true.
    if ( index(chara_wrk, 'REAL') /= 0 ) keyword_real_exist = .true.
!
    if (keyword_real_exist) then
      if (debug_mode) write(*,'(a)') '  INFO:keyword found : real'
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Search the keyword of 'symmetic' on the first line
!
    keyword_symmetric_exist = .false.
    if ( index(chara_wrk, 'symmetric') /= 0 ) keyword_symmetric_exist = .true.
    if ( index(chara_wrk, 'Symmetric') /= 0 ) keyword_symmetric_exist = .true.
    if ( index(chara_wrk, 'SYMMETRIC') /= 0 ) keyword_symmetric_exist = .true.
!
    if (keyword_symmetric_exist) then
      if (debug_mode) write(*,'(a)') '  INFO:keyword found : symmetric'
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Check the file type
!
    if (trim(matrix_type) == 'real_symmetric') then
      if (.not. keyword_real_exist) then
        write(*,*)'ERROR:keyword_real_exist =', keyword_real_exist
        stop
      endif
      if (.not. keyword_symmetric_exist) then
        write(*,*)'ERROR:keyword_symmetric_exist =', keyword_symmetric_exist
        stop
      endif
    else
      write(*,*)'ERROR:unsupported matrix type :', trim(matrix_type)
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot the comment lines
!
    do line_count=1, max_line_count
      read(unit_num,'(a)') chara_wrk
      if (index(chara_wrk, '%') /= 1) exit
      if (verbose_level >= 1) write(*,'(a,a)')'comment line:',trim(chara_wrk)
    enddo
!
    read(chara_wrk,*,iostat=ierr) mat_size, mat_size2, num_non_zeros
    if (ierr /=0) then
      write(*,*)'ERROR : matrix size or number of non zero elements :', trim(chara_wrk)
      stop
    endif
!
    if (mat_size /= mat_size2) then
      write(*,*)'ERROR : matrix size info :', mat_size, mat_size2
      stop
    endif
!
    if (num_non_zeros <= 0) then
      write(*,*)'ERROR:num_non_zeros =',num_non_zeros
      stop
    endif
!
    if (mat_size <= 0) then
      write(*,*)'ERROR:mat_size =',mat_size
      stop
    endif
!
    if ( num_non_zeros > mat_size*mat_size ) then
      write(*,*)'ERROR : imcompatible mat_size, number of non zero elements =', mat_size, num_non_zeros
      stop
    endif
!
    if (debug_mode) write(*,*)'mat_size      =', mat_size
    if (debug_mode) write(*,*)'num_non_zeros =', num_non_zeros
!
!
  end subroutine read_matrix_file_header
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_matrix_file_value(verbose_level, unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
!
    implicit none
    integer,            intent(in)  :: verbose_level
    integer,            intent(in)  :: unit_num
    integer,            intent(in)  :: mat_size
    integer,            intent(in)  :: num_non_zeros
    real(kind(1.d0)),  allocatable, intent(inout) :: mat_value(:,:)
    integer, allocatable, intent(inout) :: mat_suffix(:,:)
!
    integer :: line_count
    logical :: debug_mode
    integer :: i, j, ierr
    real(kind(1.d0))    :: value_wrk
!
    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif
!
    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_value'
    if (debug_mode) write(*,'(a)')' ONLY REAL-SYMMETRIC MATRIX is supported'
!
    if (debug_mode) write(*,*)'size(mat_value,1)=',size(mat_value,1)
    if (debug_mode) write(*,*)'size(mat_value,2)=',size(mat_value,2)
!
    if (debug_mode) write(*,*)'size(mat_suffix,1)=',size(mat_suffix,1)
    if (debug_mode) write(*,*)'size(mat_suffix,2)=',size(mat_suffix,2)
!
    if (size(mat_value,1) /= num_non_zeros) then
      write(*,*)'ERROR(read_matrix_file_value): size(mat_value,1)=',size(mat_value,1)
      stop
    endif
!
    if (size(mat_value,2) /= 1) then
      write(*,*)'ERROR(read_matrix_file_value): size(mat_value,2)=',size(mat_value,2)
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    do line_count=1, num_non_zeros
      read(unit_num, *, iostat=ierr) i, j, value_wrk
      if (ierr /= 0) then
        write(*,*)'ERROR: file reading for matrix value'
        stop
      endif
!
!     if (debug_mode) write(*,*)'i,j, matrix_value=',i, j, value_wrk
!
      if ( (i < 1) .or. ( i > mat_size ) ) then
        write(*,*)'ERROR: file reading for matrix suffix'
        stop
      endif
!
      if ( (j < 1) .or. ( j > mat_size ) ) then
        write(*,*)'ERROR: file reading for matrix suffix'
        stop
      endif
!
      mat_value(line_count, 1) = value_wrk
      mat_suffix(1,line_count) = i
      mat_suffix(2,line_count) = j
    enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   if (num_non_zeros /= size(mat_value,1)) then
!   endif
!
  end subroutine read_matrix_file_value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine print_matrix(name, mat, m, n)
    character(*), intent(in) :: name
    real(kind(1.d0)), intent(in) :: mat(:, :)
    integer :: i, j, m, n
    if (m < 0) then
       m = size(mat, 1)
    end if
    if (n < 0) then
       n = size(mat, 2)
    end if
    do j = 1, n
       do i = 1, m
          print ' (A, "(", I6, ",", I6, ")=", D30.18) ', name, i, j, mat(i, j)
       end do
    end do
  end subroutine print_matrix

end module matrix_io
