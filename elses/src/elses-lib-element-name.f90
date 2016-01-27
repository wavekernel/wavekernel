!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_element_database
!
!  NOTE: One should call "make_db_element_name" first.
!
   implicit none
   integer, parameter            :: db_length = 200
   character(len=8), allocatable :: db_element_name(:)
!
   private
   public :: get_element_name
   public :: get_atomic_num
   public :: make_db_element_name
   public :: show_db_element_name
!
   public :: db_element_name
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get element name from atomic number
!
  subroutine get_element_name(atomic_num, element_name)
   implicit none
   integer,          intent(in)  :: atomic_num
   character(len=8), intent(out) :: element_name
!
   if ( (atomic_num < 1) .or. (atomic_num > size(db_element_name,1)) ) then
     write(*,*) 'ERROR(make_db_element_name):atomic_num=', atomic_num
     stop
   endif   
!
   if (.not. allocated(db_element_name)) then
     write(*,*) 'ERROR:db_element_name is not allocated'
     stop
   endif   
!
   element_name=trim(db_element_name(atomic_num))
!
   if (trim(element_name) == '') then
     write(*,*)'ERROR(get_element_name):atomic_num=', atomic_num
     stop
   endif   
!
  end subroutine get_element_name
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get atomic number from element name
!
  subroutine get_atomic_num(element_name, atomic_num)
   implicit none
   character(len=*), intent(in)   :: element_name
   integer,          intent(out)  :: atomic_num
   character(len=8)               :: element_name_wrk
   integer                        :: j
!
   if (trim(element_name) == "") then
     stop 'ERROR(get_atomic_num)' 
   endif   
!
   if (.not. allocated(db_element_name)) then
     write(*,*) 'ERROR:db_element_name is not allocated'
     stop
   endif   
!
   atomic_num = -1
!
   do j=1,size(db_element_name,1)
     element_name_wrk=trim(db_element_name(j))
!    write(*,*)'try =', element_name_wrk, element_name
     if (trim(element_name_wrk) == trim(element_name)) then
       atomic_num=j 
       exit
     endif   
   enddo   
!
   if (atomic_num == -1) then
     write(*,*)'ERROR(get_atomic_num):element_name=',trim(element_name)
     stop
   endif   
!
  end subroutine get_atomic_num
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Show the database of element name (optional)
!
  subroutine show_db_element_name
!
   implicit none
   integer  :: j
   character(len=8) :: name_wrk
!
   if (.not. allocated(db_element_name)) then
     write(*,*) 'ERROR:db_element_name is not allocated'
     stop
   endif   
!
   write(*,*) '@@ show_db_element_name'
   do j=1,db_length
     name_wrk=trim(db_element_name(j))
     if (name_wrk /= '') write(*,*) j, name_wrk
   enddo
!   
  end subroutine show_db_element_name
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make the database for element name, if it is not made
!
  subroutine make_db_element_name
!
   implicit none
   integer  :: ierr
!
   if (allocated(db_element_name)) return
!  
   allocate(db_element_name(db_length), stat=ierr) 
   if (ierr /= 0) then
     stop 'Alloc error(db_element_name)'
   endif   
!    
   db_element_name(:)  = ''
!
   db_element_name(1)  = 'H'
   db_element_name(2)  = 'He'
   db_element_name(3)  = 'Li'
   db_element_name(4)  = 'Be'
   db_element_name(5)  = 'B'
   db_element_name(6)  = 'C'
   db_element_name(7)  = 'N'
   db_element_name(8)  = 'O'
   db_element_name(9)  = 'F'
   db_element_name(10) = 'Ne'
   db_element_name(11) = 'Na'
   db_element_name(12) = 'Mg'
   db_element_name(13) = 'Al'
   db_element_name(14) = 'Si'
   db_element_name(15) = 'P'
   db_element_name(16) = 'S'
   db_element_name(17) = 'Cl'
   db_element_name(18) = 'Ar'
   db_element_name(19) = 'K'
   db_element_name(20) = 'Ca'
   db_element_name(21) = 'Sc'
   db_element_name(22) = 'Ti'
   db_element_name(23) = 'V'
   db_element_name(24) = 'Cr'
   db_element_name(25) = 'Mn'
   db_element_name(26) = 'Fe'
   db_element_name(27) = 'Co'
   db_element_name(28) = 'Ni'
   db_element_name(29) = 'Cu'
   db_element_name(30) = 'Zn'
   db_element_name(31) = 'Ga'
   db_element_name(32) = 'Ge'
   db_element_name(33) = 'As'
   db_element_name(34) = 'Se'
   db_element_name(35) = 'Br'
   db_element_name(36) = 'Kr'
   db_element_name(37) = 'Rb'
   db_element_name(38) = 'Sr'
   db_element_name(39) = 'Y'
   db_element_name(40) = 'Zr'
   db_element_name(41) = 'Nb'
   db_element_name(42) = 'Mo'
   db_element_name(43) = 'Tc'
   db_element_name(44) = 'Ru'
   db_element_name(45) = 'Rh'
   db_element_name(46) = 'Pd'
   db_element_name(47) = 'Ag'
   db_element_name(48) = 'Cd'
   db_element_name(49) = 'In'
   db_element_name(50) = 'Sn'
   db_element_name(51) = 'Sb'
   db_element_name(52) = 'Te'
   db_element_name(53) = 'I'
   db_element_name(54) = 'Xe'
   db_element_name(55) = 'Cs'
   db_element_name(56) = 'Ba'
!
   db_element_name(57) = 'La'
   db_element_name(58) = 'Ce'
   db_element_name(59) = 'Pr'
   db_element_name(60) = 'Nd'
   db_element_name(61) = 'Pm'
   db_element_name(62) = 'Sm'
   db_element_name(63) = 'Eu'
   db_element_name(64) = 'Gd'
   db_element_name(65) = 'Tb'
   db_element_name(66) = 'Dy'
   db_element_name(67) = 'Ho'
   db_element_name(68) = 'Er'
   db_element_name(69) = 'Tm'
   db_element_name(70) = 'Yb'
   db_element_name(71) = 'Lu'
!
   db_element_name(72) = 'Hf'
   db_element_name(73) = 'Ta'
   db_element_name(74) = 'W'
   db_element_name(75) = 'Re'
   db_element_name(76) = 'Os'
   db_element_name(77) = 'Ir'
   db_element_name(78) = 'Pt'
   db_element_name(79) = 'Au'
   db_element_name(80) = 'Hg'
   db_element_name(81) = 'Tl'
   db_element_name(82) = 'Pb'
   db_element_name(83) = 'Bi'
   db_element_name(84) = 'Po'
   db_element_name(85) = 'At'
   db_element_name(86) = 'Rn'
!
  end subroutine make_db_element_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_element_database

! program main
!  use M_lib_element_database
!  implicit none
!  integer :: j , atomic_num
!  character(len=8) :: element_name_wrk
!
! do j=1,100
!   call get_element_name(j, element_name_wrk)
!   if (trim(element_name_wrk) /= '') then
!     call get_atomic_num(element_name_wrk, atomic_num)
!     write(*,*)'test', j, trim(element_name_wrk), atomic_num
!   endif   
! enddo   
!
! end program main
