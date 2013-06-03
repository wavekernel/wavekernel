module M_command_argument
  implicit none
  private
  public :: read_command_argument
!
contains
!
  subroutine read_command_argument(verbose_level, matrix_filename, eqn_type, matrix_type, solver_type, n_vec, n_check_vec)
    implicit none
    integer,            parameter   :: num_filename = 2
    integer,            intent(out) :: verbose_level
    character(len=256), intent(out) :: eqn_type
    character(len=256), intent(out) :: matrix_type
    character(len=256), intent(out) :: solver_type
    integer, intent(out) :: n_vec, n_check_vec
!
    integer :: j,n
    integer :: k
!
    integer :: count_filename
    character(len=256) :: chara_wrk
    character(len=256) :: matrix_filename(num_filename)
!
    logical :: debug_mode
!
    debug_mode         = .false.
!   debug_mode         = .true.
!
    verbose_level      = 0
    matrix_filename(1) = 'matrix_A.mtx'
    matrix_filename(2) = 'matrix_B.mtx'
    eqn_type           = 'standard'         ! (fixed)
    matrix_type        = 'real_symmetric'   ! (fixed)
    solver_type        = 'lapack'
!
!   if (verbose_level >= 100) then
!     debug_mode = .true.
!   else
!     debug_mode = .false.
!   endif
!
!   if (debug_mode) write(*,'(a)') '@@ read command arguments'
!
    n=command_argument_count()
    if (debug_mode) write(*,*)'number of command arguments = ',n
    if (n == 0) then
      write(*,*)'ERROR:No command argument'
      stop
    endif
!
    count_filename=0
    do j=1, n
      call get_command_argument(j,chara_wrk)
      if (debug_mode) write(*,*)'message= ',trim(chara_wrk)
      if ( chara_wrk(1:8) == '-solver=') then
        read(unit=chara_wrk(9:),fmt=*) solver_type
        cycle
      else if (chara_wrk(1:7) == '-n_vec=') then
        read(unit=chara_wrk(8:), fmt=' (I10) ') n_vec
        cycle
      else if (chara_wrk(1:13) == '-n_check_vec=') then
        read(unit=chara_wrk(14:), fmt=' (I10) ') n_check_vec
        cycle
      endif
      if ( chara_wrk(1:1) == '-') then
        select case (chara_wrk(2:))
          case('verbose')
            verbose_level=1
          case('debug')
            verbose_level=100
            debug_mode = .true.
!         case('solver=':'solver=:')
!           if (debug_mode) write(*,*)' solver is detected'
!           read(unit=chara_wrk(9:),fmt=*) solver_type
          case default
            write(*,*)'ERROR:unknown option: ', trim(chara_wrk)
            stop
        end select
      else
        k=len_trim(chara_wrk)
        if ( chara_wrk(k-3:k) /= '.mtx' ) then
          write(*,*)'ERROR: unknown file extension :', chara_wrk
          stop
        endif
        count_filename=count_filename+1
        if (count_filename > size(matrix_filename,1)) then
          write(*,*)'ERROR:too many filename'
          stop
        else
          matrix_filename(count_filename)=trim(chara_wrk)
          if (debug_mode) write(*,*)'filename =', count_filename, trim(matrix_filename(count_filename))
        endif
      endif
    enddo
!

!
  end subroutine read_command_argument
!
end module M_command_argument


