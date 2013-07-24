module matrix_io
  use command_argument, only : matrix_info
  implicit none

  type sparse_mat
    integer :: size, num_non_zeros
    real(kind(1.d0)), allocatable :: value(:)
    integer, allocatable :: suffix(:, :)
  end type sparse_mat

  private
  public :: read_matrix_file, print_matrix, sparse_mat

contains
  subroutine read_matrix_file(filename, info, matrix)
    implicit none

    character(*), intent(in) :: filename
    type(matrix_info), intent(in) :: info
    type(sparse_mat), intent(out) :: matrix

    integer, parameter :: iunit = 8
    integer :: ierr
    integer dummy1, dummy2

    allocate(matrix%suffix(2, info%entries), matrix%value(info%entries), &
         stat = ierr)

    if (ierr /= 0) then
      stop 'ERROR read_matrix_filematrix_io : allcation failed'
    end if

    open(iunit, file = filename)

    ! read_matrix_file_header is added to skip comment lines
    call read_matrix_file_header(0, iunit, 'real_symmetric', dummy1, dummy2)
    call read_matrix_file_value(0, iunit, info%rows, &
         info%entries, matrix%value, matrix%suffix)

    close(iunit)
  end subroutine read_matrix_file


  subroutine read_matrix_file_header(verbose_level, unit_num, matrix_type, mat_size, num_non_zeros)
    implicit none
    integer, intent(in) :: verbose_level
    integer, intent(in) :: unit_num
    character(len=256), intent(in) :: matrix_type
    integer, intent(out) :: mat_size
    integer, intent(out) :: num_non_zeros

    integer, parameter :: max_line_count = 100

    logical :: keyword_real_exist
    logical :: keyword_symmetric_exist
    integer :: line_count
    integer :: ierr
    integer :: mat_size2

    logical :: debug_mode

    character(len=1024) :: chara_wrk

    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif

    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_header'

    if (debug_mode) write(*,*)'verbose level =', verbose_level

    ! Read the first line
    read(unit_num,'(a)') chara_wrk
    if (debug_mode) write(*,*)'first line=',trim(chara_wrk)

    ! Check the header of the first line
    if (chara_wrk(1:14) /= '%%MatrixMarket') then
      write(*,*)'ERROR:Fist line :', trim(chara_wrk)
      stop
    else
      if (debug_mode) write(*,*)'INFO:header of file =', chara_wrk(1:14)
    endif

    ! Seach the keyword of 'real' on the first line
    keyword_real_exist = .false.
    if (index(chara_wrk, 'real') /= 0) keyword_real_exist = .true.
    if (index(chara_wrk, 'Real') /= 0) keyword_real_exist = .true.
    if (index(chara_wrk, 'REAL') /= 0) keyword_real_exist = .true.

    if (keyword_real_exist) then
      if (debug_mode) write(*,'(a)') '  INFO:keyword found : real'
    endif

    ! Search the keyword of 'symmetic' on the first line
    keyword_symmetric_exist = .false.
    if (index(chara_wrk, 'symmetric') /= 0) keyword_symmetric_exist = .true.
    if (index(chara_wrk, 'Symmetric') /= 0) keyword_symmetric_exist = .true.
    if (index(chara_wrk, 'SYMMETRIC') /= 0) keyword_symmetric_exist = .true.

    if (keyword_symmetric_exist) then
      if (debug_mode) write(*,'(a)') '  INFO:keyword found : symmetric'
    endif

    ! Plot the comment lines
    do line_count = 1, max_line_count
      read(unit_num,'(a)') chara_wrk
      if (index(chara_wrk, '%') /= 1) exit
      if (verbose_level >= 1) write(*,'(a,a)')'comment line:',trim(chara_wrk)
    enddo

    read (chara_wrk, *, iostat=ierr) mat_size, mat_size2, num_non_zeros
    if (ierr /= 0) then
      write(*,*)'ERROR : matrix size or number of non zero elements :', trim(chara_wrk)
      stop
    endif

    if (mat_size /= mat_size2) then
      write(*,*)'ERROR : matrix size info :', mat_size, mat_size2
      stop
    endif

    if (num_non_zeros <= 0) then
      write(*,*)'ERROR:num_non_zeros =',num_non_zeros
      stop
    endif

    if (mat_size <= 0) then
      write(*,*)'ERROR:mat_size =',mat_size
      stop
    endif

    if (num_non_zeros > mat_size * mat_size) then
      write(*,*)'ERROR : imcompatible mat_size, number of non zero elements =', mat_size, num_non_zeros
      stop
    endif

    if (debug_mode) write(*,*)'mat_size      =', mat_size
    if (debug_mode) write(*,*)'num_non_zeros =', num_non_zeros
  end subroutine read_matrix_file_header


  subroutine read_matrix_file_value(verbose_level, unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
    implicit none

    integer, intent(in) :: verbose_level
    integer, intent(in) :: unit_num
    integer, intent(in) :: mat_size
    integer, intent(in) :: num_non_zeros
    real(kind(1.d0)), allocatable, intent(inout) :: mat_value(:)
    integer, allocatable, intent(inout) :: mat_suffix(:,:)

    integer :: line_count
    logical :: debug_mode
    integer :: i, j, ierr
    real(kind(1.d0)) :: value_wrk

    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif

    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_value'
    if (debug_mode) write(*,'(a)')' ONLY REAL-SYMMETRIC MATRIX is supported'

    if (debug_mode) write(*,*)'size(mat_value)=',size(mat_value)

    if (debug_mode) write(*,*)'size(mat_suffix,1)=',size(mat_suffix,1)
    if (debug_mode) write(*,*)'size(mat_suffix,2)=',size(mat_suffix,2)

    if (size(mat_value, 1) /= num_non_zeros) then
      write(*,*)'ERROR(read_matrix_file_value): size(mat_value,1)=', size(mat_value,1)
      stop
    endif

    do line_count = 1, num_non_zeros
      read(unit_num, *, iostat=ierr) i, j, value_wrk
      if (ierr /= 0) then
        write(*,*)'ERROR: file reading for matrix value'
        stop
      endif

      ! if (debug_mode) write(*,*)'i,j, matrix_value=',i, j, value_wrk

      if (i < 1 .or. i > mat_size) then
        write(*,*)'ERROR: file reading for matrix suffix'
        stop
      endif

      if (j < 1 .or. j > mat_size) then
        write(*,*)'ERROR: file reading for matrix suffix'
        stop
      endif

      mat_value(line_count) = value_wrk
      mat_suffix(1,line_count) = i
      mat_suffix(2,line_count) = j
    enddo

    ! if (num_non_zeros /= size(mat_value,1)) then
    ! endif
  end subroutine read_matrix_file_value


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
