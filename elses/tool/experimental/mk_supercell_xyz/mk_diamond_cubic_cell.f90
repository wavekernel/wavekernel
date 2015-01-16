! Simple tool for generating diamond structure in the eight-atom cubic cell 
!   input file  : (none)
!   output file : diamond_in_cubic_cell.xyz
!                ( structure in the extended XYZ format )
!
program main
  implicit none
  integer :: noa
  real(8) :: ax, ay, az, dnn
  real(8)          :: pos(3,8)
  integer :: j
  integer, parameter          :: fd=31
  character(len=*), parameter :: filename='diamond_in_cubic_cell.xyz'
!
  noa=8
!
  dnn=1.536329d0   
!  ---> n.n atomic distance in Angstrom
!   Ref. Xu et al, J. Phys.: Condens. Matter 4, 6047-6054 (1992) 
!
  ax=dnn*4.0d0/dsqrt(3.0d0)   ! lattice constant in Angstrom
  ay=ax
  az=ax
!
  pos(1,1)=0.0d0
  pos(2,1)=0.0d0
  pos(3,1)=0.0d0
!
  pos(1,2)=0.25d0
  pos(2,2)=0.25d0
  pos(3,2)=0.25d0
!
  pos(1,3)=0.5d0
  pos(2,3)=0.5d0
  pos(3,3)=0.0d0
!
  pos(1,4)=0.75d0
  pos(2,4)=0.75d0
  pos(3,4)=0.25d0
!
  pos(1,5)=0.5d0
  pos(2,5)=0.0d0
  pos(3,5)=0.5d0
!
  pos(1,6)=0.75d0
  pos(2,6)=0.25d0
  pos(3,6)=0.75d0
!
  pos(1,7)=0.0d0
  pos(2,7)=0.5d0
  pos(3,7)=0.5d0
!
  pos(1,8)=0.25d0
  pos(2,8)=0.75d0
  pos(3,8)=0.75d0
!
  open(fd,file=filename)
  write(fd,'(i15, 3f30.20)') noa, ax, ay, az
  write(fd,'(a)')'Diamond structure in the eight-atom cubic cell'
!
  do j=1,8
    write(fd,'(a,3f30.20)') 'C  ', pos(1,j)*ax, pos(2,j)*ay, pos(3, j)*az
  enddo
!   
end program
