!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
program make_group
  implicit none
!
  type :: atom_info_type
     real(8)              :: position_angst(3)
     character(len=16)    :: element_name
     character(len=16)    :: motion
  end type atom_info_type
!
  type :: structure_info_type
     type(atom_info_type), allocatable :: atom_info(:)
     integer              :: total_number_of_atoms
     real(8)              :: cell_size_angst(3)
  end type structure_info_type
!
  type :: bond_info_type
     integer              :: total_number
     real(8), allocatable :: cutoff(:)
     character(len=8), allocatable :: element(:,:)
  end type bond_info_type
!
  type(structure_info_type) :: structure_info
!
  type(bond_info_type)      :: bond_info
!
  integer, allocatable      :: hop_mat(:,:)  ! hopping matrix
  integer, allocatable      :: group_id(:)
!
! real(8), parameter :: cutoff_C_H = 1.3d0
! real(8), parameter :: cutoff_C_C = 1.6d0
!
  call main
!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the main function of this program
!
  subroutine main
   implicit none
   character(len=80)  :: file_sample
   character(len=80)  :: file_bond_info
   character(len=80)  :: file_result_group
   character(len=80)  :: filename_list(3)
   logical            :: cell_info_in_xyz_file
!
   filename_list(:)=''
   call read_options(filename_list)
!
   file_sample       = trim(adjustl(filename_list(1)))
   file_bond_info    = trim(adjustl(filename_list(2)))
   file_result_group = trim(adjustl(filename_list(3)))
!
   call get_bond_info(file_bond_info)
!
   call get_info_from_xyz_file(file_sample, cell_info_in_xyz_file)
!
   write(*,*)'noa=', structure_info%total_number_of_atoms
!
   call set_hop_mat
!
   call set_group_id
!
   call save_group_id(file_result_group) ! Save the group ID information into file
!
  end subroutine main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine save_group_id(filename)
   implicit none
   integer, parameter :: fd = 2
   character(len=*), intent(in)  :: filename
   character(len=256)  :: chara_wrk
   integer             :: i, ierr
   integer             :: noa
   integer             :: id, id_max
   integer, allocatable :: mem_group(:) ! The number of member atoms for each group
!
   write(*,*)'INFO:save_groupid: filename=',trim(filename)
!
   id_max=maxval(group_id(:))
   write(*,*)'Number of groups=', id_max
   allocate(mem_group(id_max), stat=ierr)
   if (ierr /= 0) stop 'ERROR:alloc. error : mem_group'
   mem_group(:)=0
!
   noa = structure_info%total_number_of_atoms
   do i=1,noa
     id=group_id(i)
     mem_group(id)=mem_group(id)+1
   enddo
!
   do id=1,id_max
     write(*,*)'INFO:mem_group =', id, mem_group(id)
   enddo
!
   open(fd,file=filename,status='unknown')
!
   chara_wrk='# Sample'
   write(fd,'(a)') trim(chara_wrk)
   chara_wrk='# number of atoms, number of groups, maximum number of group atoms'
   write(fd,'(a)') trim(chara_wrk)
   write(fd,'(3i10)') noa, id_max, maxval(mem_group)
   chara_wrk='#'
   write(fd,'(a)') trim(chara_wrk)
   do i=1,noa
     id=group_id(i)
     write(fd,'(2i10)') i, id
   enddo
   close(fd)
!
  end subroutine save_group_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_group_id
    implicit none
    integer :: noa
    integer :: try_max
    integer :: itry, ierr
    integer :: i_seed
    integer, allocatable :: vect_wrk(:)
    integer :: id, id_max
    integer :: i
    integer :: count, count_old
!
    write(*,'(a)')'INFO:set_group_id'
!
    noa = structure_info%total_number_of_atoms
!
   allocate(group_id(noa), stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR:alloc. error : group_id'
     stop
   endif
   group_id(:)=0
!
   allocate(vect_wrk(noa), stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR:alloc. error : vect_wrk'
     stop
   endif
   vect_wrk(:)=0
!
    try_max=noa
    id_max=noa
    do id=1, id_max
      call find_vacant_comp(i_seed)
      if (i_seed == -1) exit
      write(*,*)'id, i_seed=',id, i_seed
      vect_wrk(:)=0
      vect_wrk(i_seed)=1
      count_old=0
      do itry=1,try_max
        vect_wrk = matmul(hop_mat, vect_wrk)
        call normalize_vect(vect_wrk)
        count = sum(vect_wrk)
        write(*,*)'count=',id, count
        if (count == count_old) exit
        count_old=count
      enddo  
      do i=1,noa
        if (vect_wrk(i) /= 0) then
          if (group_id(i) /= 0) then
            write(*,*)'ERROR, group_id'
            stop
          endif
          group_id(i) = id
        endif
      enddo
    enddo
!
  end subroutine set_group_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine normalize_vect(vect_wrk)
    implicit none
    integer, intent(inout) :: vect_wrk(:)
    integer                :: j
!
    do j=1, size(vect_wrk,1)
      if (vect_wrk(j) /= 0) vect_wrk(j)=1
    enddo
!
  end subroutine normalize_vect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine find_vacant_comp(i_vacant_comp)
    implicit none
    integer, intent(out) :: i_vacant_comp ! vacant component
    integer              :: j

    i_vacant_comp = -1 ! dummy value
!
    do j=1,structure_info%total_number_of_atoms
      if (group_id(j) == 0) then
        i_vacant_comp=j
        exit
      endif   
    enddo
!
!   if (i_vacant_comp == -1) then
!     write(*,*)'ERROR(find_vacant_comp)'
!     stop
!   endif
!
  end subroutine find_vacant_comp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_hop_mat
    implicit none
    integer :: noa
    integer :: i,j 
    integer :: counter
    integer :: ierr
    real(8) :: cutoff_wrk, distance
    character(len=16) :: elem_name_i
!
    write(*,'(a)')'INFO:set_hop_mat only for the non-periodic case'
    noa = structure_info%total_number_of_atoms
!
   allocate(hop_mat(noa,noa), stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR:alloc. error : hop_mat'
     stop
   endif
   hop_mat(:,:)=0
!
   do i=1,noa
     hop_mat(i,i)=1
   enddo
!
   do i=1,noa
     do j=1,noa
       if (i == j) cycle
       call set_cutoff(i,j,cutoff_wrk)
       call get_distance(i,j,distance)
!      write(*,'(a,2i10,2f10.5)')'i, j, dist, cutoff=',i,j,distance, cutoff_wrk
       if (distance < cutoff_wrk) hop_mat(i,j)=1
       if (distance < cutoff_wrk) hop_mat(j,i)=1
     enddo
   enddo
!
!  do i=1,noa
!    do j=1,noa
!      write(*,'(a,3i10)')'i, j, hop_mat=',i,j,hop_mat(i,j)
!    enddo
!  enddo
!
   do i=1,noa
     counter=0
     elem_name_i=trim(structure_info%atom_info(i)%element_name)
     do j=1,noa
       if (i == j) cycle
       if (hop_mat(i,j) /=0) counter=counter+1
     enddo
     write(*,'(a,i10,a,a,a,i10)')'i, num. of bonds=',i,' ', trim(elem_name_i), ' ', counter
   enddo
!
  end subroutine set_hop_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_cutoff(i,j,distance)
   implicit none
   integer,    intent(in)   :: i, j
   real(8),    intent(out)  :: distance
   real(8)                  :: cutoff_wrk
   character(len=16) :: elem_name_i
   character(len=16) :: elem_name_j
!
   integer           :: k
   character(len=16) :: elem_name_a
   character(len=16) :: elem_name_b
   logical           :: bond_set
!
   elem_name_i=trim(structure_info%atom_info(i)%element_name)
   elem_name_j=trim(structure_info%atom_info(j)%element_name)
!
   cutoff_wrk = -1.0d0 ! dummy value
!
   do k=1, bond_info%total_number
     bond_set = .false.
     elem_name_a=bond_info%element(1,k)
     elem_name_b=bond_info%element(2,k)
     if ( (trim(elem_name_i) == trim(elem_name_a)) .and. (trim(elem_name_j) == trim(elem_name_b)) ) bond_set = .true.
     if ( (trim(elem_name_i) == trim(elem_name_b)) .and. (trim(elem_name_j) == trim(elem_name_a)) ) bond_set = .true.
     if (bond_set) then
       cutoff_wrk = bond_info%cutoff(k)
     endif
   enddo
!
   distance = cutoff_wrk
!
  end subroutine set_cutoff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_distance(i,j,distance)
   implicit none
   integer,    intent(in)   :: i, j
   real(8),    intent(out)  :: distance
   real(8)                  :: vect(3)
!
   vect(:)= structure_info%atom_info(i)%position_angst(:) &
&         - structure_info%atom_info(j)%position_angst(:)
!
   distance = dsqrt(dot_product(vect,vect))
!
  end subroutine get_distance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_info_from_xyz_file( filename, cell_info_in_xyz_file)
   implicit none
   character(len=*),    intent(in)    :: filename
   logical         ,    intent(inout) :: cell_info_in_xyz_file
   integer, parameter :: fd = 1
   integer            :: noa
   real(8)            :: ax, ay, az  ! cell sizes in Agnstrom
   integer            :: ierr, j
   character(len=1024)  :: chara_wrk
   character(len=8)     :: elem_name_wrk ! work array
   real(8)              :: data_wrk(3)   ! work array
!
   write(*,*)'INFO:get_basic_info: filename=',trim(filename)
   open(fd,file=filename,status='old')
!
   structure_info%total_number_of_atoms = -1     ! dummy value
   structure_info%cell_size_angst(3)    = -1.0d0 ! dummy value
   read(fd, '(a)') chara_wrk  ! read the first line
!  write(*,*)'INFO:chara_wrk        =', trim(chara_wrk)
   read(chara_wrk, *, iostat=ierr) noa, ax, ay, az
   if (ierr == 0) then
     cell_info_in_xyz_file = .true. 
!    write(*,*)'INFO:number of atoms        =', noa
!    write(*,*)'INFO:cell sizes [Angstrom]  =', ax, ay, az
     structure_info%total_number_of_atoms = noa
     structure_info%cell_size_angst(1)    = ax
     structure_info%cell_size_angst(2)    = ay
     structure_info%cell_size_angst(3)    = az
   else   
     read(chara_wrk, *, iostat=ierr) noa
     cell_info_in_xyz_file = .false. 
     structure_info%total_number_of_atoms = noa
!    write(*,*)'INFO:number of atoms        =', noa
   endif   
!
   allocate(structure_info%atom_info(noa), stat=ierr)
   if (ierr /= 0) then
     write(*,*)'ERROR:alloc. error for structure_info%atom_info:noa=',noa
     stop
   endif
!
   read(fd, '(a)') chara_wrk   ! read the second (comment) line 
!
   do j=1,noa
     read(fd, *) elem_name_wrk, data_wrk(1:3)
!    write(*,*)'INFO:j, name=', j, trim(elem_name_wrk)
!    write(*,*)'INFO:j, data(1:3)=', j, data_wrk(1:3)
     structure_info%atom_info(j)%element_name = trim(elem_name_wrk)
     structure_info%atom_info(j)%position_angst(1:3)=data_wrk(1:3)
   enddo
!
   structure_info%atom_info(:)%motion = 'free'    ! default setting
!
   close(fd)
  end subroutine get_info_from_xyz_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_bond_info(filename)
   implicit none
   character(len=*),    intent(in)    :: filename
   integer :: i, n, ierr
!
!
   write(*,'(a)')'@@ get_bond_info'
   open(30, file=filename, status='old')
!
   read(30,*) n
   write(*,'(a,i5)') ' total number of bond info =', n
   bond_info%total_number=n
!
   allocate(bond_info%cutoff(n), stat=ierr)
   if (ierr /= 0) stop 'Alloc. error in cutoff'
!
   allocate(bond_info%element(2,n), stat=ierr)
   if (ierr /= 0) stop 'Alloc. error in member_atom'
!
   do i=1,n
     read(30,*) bond_info%element(1:2,i), bond_info%cutoff(i)
     write(*,*) trim(bond_info%element(1,i)), ' ', trim(bond_info%element(2,i)), ' ',  bond_info%cutoff(i)
   enddo
!
   close(30)
!
  end subroutine get_bond_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_options(filename_list)
   implicit none
   character(len=80), intent(inout) :: filename_list(3)
   integer                          :: i, counter
   character(len=80)                :: argc
!
   counter = iargc()
!
   if (counter /= 3) then
     write(*,*)' ERROR(read_options):counter=', counter
     stop
   endif
!
   do i=1,counter
     call getarg(i,argc)
     filename_list(i) = trim(adjustl(argc))
   enddo
!
  end subroutine
!
end program make_group

