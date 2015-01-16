module M_gid_constraint
!
  implicit none
  integer, parameter   :: DOUBLE_PRECISION=kind(1d0)
!
  private
  public :: add_constraint_w_groups
!
  contains
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Add constraints on motion with groups
!
  subroutine add_constraint_w_groups
!
    use M_config,             only : config  ! unchanged
    use M_gid_basic_routines, only : fix_group_center !(routine)
    use M_lib_read_w_comment, only : read_w_comment   !(routine)
    use elses_mod_file_io,    only : vacant_unit      !(function)
    implicit none
    character(len=*), parameter :: filename='input_group_id_control.txt'
    integer :: lu, vl
    integer :: fd
    logical :: end_of_file
    logical :: file_exist
    character(len=10240) :: line_chara
    character(len=10240) :: chara_mode
    integer :: j, command_num, ierr, gid
!
    vl = config%option%verbose  ! verbose level
    lu = config%calc%distributed%log_unit
    fd = vacant_unit()
!
    if ( .not. config%calc%constraint_w_group ) then 
      if ( vl > 0 ) then
        if (lu > 0) write(lu,'(a)') '@@@ add_constraint_w_groups ...is skipped'
      endif
      return
    endif
!
    inquire(file=trim(filename), exist=file_exist)
    if (.not. file_exist) then
      write(*,*)'ERROR(add_constraint_w_groups):filename=', trim(filename)
      stop
    endif
!
    open(fd,file=filename,status='old')
!
    if (lu > 0) write(lu,'(a)') '@@@ add_constraint_w_groups'
!
    call read_w_comment(fd, '#', line_chara, end_of_file)
!   write(*,*)'line = ', trim(line_chara)
    if (end_of_file) then
      write(*,*)'ERROR(add_constraint_w_groups):end_of_file'
      stop
    endif
    read(line_chara,*,iostat=ierr) command_num
    if (ierr /= 0) then
      write(*,*)'ERROR(add_constraint_w_groups):ierr=',ierr
      stop
    endif
!
    if (lu > 0) write(lu,'(a,i10)')'command_num = ', command_num
!
    do j=1,command_num
      call read_w_comment(fd, '#', line_chara, end_of_file)
!     write(*,*)'line = ', trim(line_chara)
      if (end_of_file) then
        write(*,*)'ERROR(add_constraint_w_groups):end_of_file'
        stop
      endif
      read(line_chara,*,iostat=ierr) gid, chara_mode
      if (ierr /= 0) then
        write(*,*)'ERROR(add_constraint_w_groups):ierr=',ierr
        stop
      endif
      select case(trim(chara_mode))
        case ('fix_group_center')
          if (lu > 0) write(lu,'(a,i10)') " call fix_group_center(xyz):gid=",gid
          call fix_group_center(gid,'xyz')
        case ('fix_group_center_x')
          if (lu > 0) write(lu,'(a,i10)') " call fix_group_center(x):gid=",gid
          call fix_group_center(gid,'x')
        case ('fix_group_center_y')
          if (lu > 0) write(lu,'(a,i10)') " call fix_group_center(y):gid=",gid
          call fix_group_center(gid,'y')
        case ('fix_group_center_z')
          if (lu > 0) write(lu,'(a,i10)') " call fix_group_center(z):gid=",gid
          call fix_group_center(gid,'z')
        case default
          write(*,*)'ERROR(add_constraint_w_groups):mode=',trim(chara_mode)
          stop
      end select
    enddo
!
!
    close(fd)
!

  end subroutine add_constraint_w_groups
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_gid_constraint


