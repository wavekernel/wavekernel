! mk_supercell_xyz.f90 (build 20160215)
!
program main
  implicit none
  real(8),allocatable :: position(:,:)
  character(len=3),allocatable :: element_name(:)
  real(8) :: ax, ay, az
  integer :: noa, noa_new
  integer :: ix, iy, iz, k
  real(8) :: ddx, ddy, ddz
  integer :: nx, ny, nz
  character (len=8) :: chara1
  character (len=64) :: filename_org, filename_new
  character (len=1024) :: comment_line
  character(len=*), parameter :: filename_cell_info = 'input_supercell_size.txt'
! character(len=*), parameter :: filename_org       = 'structure_org.xyz'
! character(len=*), parameter :: filename_new       = 'structure_new.xyz'
  integer, parameter :: fd_cell_info = 29
  integer, parameter :: fd_org       = 30
  integer, parameter :: fd_new       = 31
!
  write(*,'(a)') '@@@ mk_supecell_xyz.f90 (build 20160215)'
  write(*,'(a)') ' NOTE: This code uses GET_COMMAND_ARGUMENT defined in Fortran2003'
!
  call read_parameters(nx, ny, nz, filename_org, filename_new, comment_line)
! 
  open (fd_cell_info, file=filename_cell_info)
  open (fd_org, file=filename_org)
  open (fd_new, file=filename_new)
!
  read(fd_org,*) noa, ax, ay, az
  write(*,'(a)')          '@@@ Input file : structure_org.xyz'
  write(*,'(a,i10)')      ' -->      N        =', noa
  write(*,'(a,3f20.10)')  ' -->  ax,ay,az [A] =', ax, ay, az
!
  noa_new=noa*nx*ny*nz
  write(*,'(a)')          '@@@ Output file : structure_new.xyz'
  write(*,'(a,i10)')      ' -->      N        =', noa_new
  write(*,'(a,3f20.10)')  ' -->  ax,ay,az [A] =', ax*dble(nx), ay*dble(ny), az*dble(nz)
!
  write(*,'(a)')          '@@@ Generate the output file'
!
  write(fd_new,'(i15, 3f20.10)') noa_new, ax*dble(nx), ay*dble(ny), az*dble(nz)
  write(fd_new,'(a)') trim(comment_line)
!
  allocate (position(3,noa))
  allocate (element_name(noa))
!
  read(fd_org,*) chara1
  do k=1,noa
    read(fd_org,*) element_name(k), position(1:3,k)
  enddo
!
  do iz=0,nz-1
    do iy=0,ny-1
      do ix=0,nx-1
        do k=1,noa
          ddx=position(1,k)+ax*dble(ix)
          ddy=position(2,k)+ay*dble(iy)
          ddz=position(3,k)+az*dble(iz)
          write(fd_new,'(a2,3f30.20)') element_name(k), ddx, ddy, ddz
        enddo
      enddo
     enddo
  enddo
!
  write(*,'(a)') '.... mk_supecell_xyz.f90 (build 20160215) ended without error'
!
end program

subroutine read_parameters(nx, ny, nz, filename_org, filename_new, comment_line)
  implicit none
  integer, intent(out)              :: nx, ny, nz
  character(len=*), intent(out)     :: filename_org, filename_new, comment_line
!
  integer i, len_arg, ierr
  character(len=1024) :: arg, chara_wrk, comment_line_default
!
  logical :: debug=.false.
! logical :: debug=.true.
!
  nx=-1  ! dummy value
  ny=-1  ! dummy value
  nz=-1  ! dummy valu
  filename_org = '' 
  filename_new = '' 
  comment_line = '' 
!
  write(*,*) '@@@ Read the paramters from the argument'
  i=0
!
  do 
    call get_command_argument(i, arg)
    len_arg=len_trim(arg)
    if (len_arg == 0) exit
    if (debug) write(*,*) 'INFO:command line argument:',i, trim(arg)
!
    if (index(arg, '-nx=') > 0) then
      chara_wrk=trim(arg(5:len_arg))
      read(chara_wrk, *, iostat=ierr) nx
      if (ierr /= 0) stop 'ERROR in reading NX from the command argument'
      if (debug) write(*,*) '  nx =',nx
    endif
!
    if (index(arg, '-ny=') > 0) then
      chara_wrk=trim(arg(5:len_arg))
      read(chara_wrk, *) ny
      if (debug) write(*,*) '  ny =',ny
    endif
!
    if (index(arg, '-nz=') > 0) then
      chara_wrk=trim(arg(5:len_arg))
      read(chara_wrk, *) nz
      if (debug) write(*,*) '  nz =',nz
    endif
!
    if (index(arg, '-filename_input=') > 0) then
      chara_wrk=trim(arg(17:len_arg))
      read(chara_wrk, *, iostat=ierr) filename_org
      if (ierr /= 0) stop 'ERROR in reading FILENAME_INPUT from the command argument'
      if (debug) write(*,*) '  filename_input =',trim(filename_org)
    endif
!
    if (index(arg, '-filename_output=') > 0) then
      chara_wrk=trim(arg(18:len_arg))
      read(chara_wrk, *, iostat=ierr) filename_new
      if (ierr /= 0) stop 'ERROR in reading FILENAME_OUTPUT from the command argument'
      if (debug) write(*,*) '  filename_output =',trim(filename_new)
    endif
!
    if (index(arg, '-comment=') > 0) then
      chara_wrk=trim(arg(10:len_arg))
      read(chara_wrk, *, iostat=ierr) comment_line
      if (ierr /= 0) stop 'ERROR in reading COMMENT from the command argument'
      if (debug) write(*,*) '  comment_line =',trim(comment_line)
    endif
!
    i=i+1
  enddo
!
  if (len_trim(comment_line) == 0) then
    write(comment_line,'(a,3i10)') & 
&         '# generated by mk_supercell_xyz.f90(build 20160215) : nx, ny, nz =', nx, ny, nz
  endif
!
  write(*,*)'INFO:supercell sizes  = ', nx, ny, nz
  write(*,*)'INFO:filename input   = ', trim(filename_org)
  write(*,*)'INFO:filename output  = ', trim(filename_new)
  write(*,*)'INFO:comment line     = ', trim(comment_line)
!
  ierr=0
  if (nx == -1) ierr=1
  if (ny == -1) ierr=1
  if (nz == -1) ierr=1
  if (len_trim(filename_org) == 0) ierr=1
  if (len_trim(filename_new) == 0) ierr=1
  if (len_trim(comment_line) == 0) ierr=1
  if (ierr == 1) stop 'ERROR in reading from the command argument'
!
end subroutine read_parameters
