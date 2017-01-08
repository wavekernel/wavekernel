!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
program calculate_msd_main

  implicit none
  
  integer :: io

  integer :: n

  integer :: md_step
  integer :: t0       !!! index of t=0

  integer :: natom          !!! number of atom 
  integer :: total_md_step  !!! total number of md step
  integer :: tmax           !!! total number of time step
  integer :: it0            !!! interval to sample
  integer :: t0max          !!! the maximum number of t=0
  real(8) :: dt             !!! interval [fs] between two successive md step

  real(8), allocatable :: x0(:,:,:) !!! atom position for given t=0
  integer, allocatable :: time0(:)  !!! time for given t=0


  integer, allocatable :: ntime(:) !!! number of samples for time i
  real(8), allocatable :: msd(:)   !!! mean square displacement
  real(8), allocatable :: msd_x(:), msd_y(:), msd_z(:)   !!! mean square displacement for each component

  character(len=*),parameter :: filename_00='position.xyz'

  character(len=100) :: structure_name
 
  real(8),allocatable :: atom_position(:,:)



!----------------------------------------------------------------------
  
  call read_parameter_file
 
  open(10,file=filename_00,status='old',IOSTAT=io)
  if(io>0) then
     write(*,'("!!! ERROR !!! The file postion.xyz does not exist")')
  end if

  !!! initialization !!!
  allocate(ntime(tmax))
  allocate(msd(tmax))
  allocate(msd_x(tmax))
  allocate(msd_y(tmax))
  allocate(msd_z(tmax))
  allocate(time0(t0max))
  allocate(atom_position(natom,3))
  allocate(x0(t0max,natom,3))

  atom_position(:,:)=0.0
  x0(:,:,:)=0.0

  ntime(:)=0
  msd(:)=0.0d0
  msd_x(:)=0.0d0
  msd_y(:)=0.0d0
  msd_z(:)=0.0d0
  t0=0
  !!!!!!!!!!!!!!!!!!!!!!

  do md_step=0, total_md_step-1

     !!! write(*,'("step:",I10)') md_step

     read(10,*,iostat=io) n
 
     if(io<0) then
        write(*,'("file is end")')
        exit
     end if

     !!! write(*,'("number of atom: ",I10)') n
 
     read(10,*) structure_name  
     !!! write(*,'(A100)') structure_name 

     call read_atom_position
     
     call calc_msd
        
  end do

  deallocate(atom_position) 
  deallocate(x0)
  deallocate(time0)


  call output_msd

  deallocate(ntime)
  deallocate(msd)
  deallocate(msd_x)
  deallocate(msd_y)
  deallocate(msd_z)

  close(10)

contains
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read calculation conditions from the input file
! 
  subroutine read_parameter_file

    implicit none

    character(len=*),parameter :: filename_01='input_msd.txt'


    open(10,file=filename_01,status='old',iostat=io)
    if(io>0) then
       write(*,'("!!! ERROR !!! The file input_msd_vacf.txt does not exist")')
    end if

    read(10,*) natom
    read(10,*) total_md_step
    read(10,*) tmax
    read(10,*) it0
    read(10,*) t0max
    read(10,*) dt

    write(*,'("natom        :",I10)') natom
    write(*,'("total_md_step:",I10)') total_md_step
    write(*,'("tmax         :",I10)') tmax
    write(*,'("it0          :",I10)') it0
    write(*,'("t0max        :",I10)') t0max
    write(*,'("dt [fs]      :",F10.5)') dt

    close(10)

  end subroutine read_parameter_file
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Read atom position from positio.xyz file
!
  subroutine read_atom_position

    implicit none

    integer :: i
    character(len=4) :: atom_name

    do i=1, natom
       read(10,*) atom_name,atom_position(i,1:3)
       !!!write(*,'(a4,3e23.15)') atom_name,atom_position(i,1:3)
    end do

  end subroutine read_atom_position
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate mean square displacement
!
  subroutine calc_msd

    implicit none

    integer :: i,t,tt0
    integer :: delt  !!! t-t0

    if(mod(md_step,it0)==0) then

       t0=t0+1

       tt0=mod(t0-1,t0max)+1
       time0(tt0)=md_step

       !!! write(*,'("t0         :",I10)') t0
       !!! write(*,'("tt0        :",I10)') tt0
       !!! write(*,'("time0(tt0) :",I10)') time0(tt0)

       do i=1,natom
          x0(tt0,i,1:3)=atom_position(i,1:3)
          !!! write(*,'("x0(",I10,"):",3e23.15)') i,atom_position(i,1:3)
       end do

    end if

    do t=1,min(t0,t0max)

       delt=md_step-time0(t)+1

       !!! write(*,'("md_step   :",I10)') md_step
       !!! write(*,'("t         :",I10)') t
       !!! write(*,'("time0(t)  :",I10)') time0(t)
       !!! write(*,'("delt      :",I10)') delt

       if(delt <= tmax) then

          ntime(delt)=ntime(delt)+1

          do i=1,natom
 
            msd(delt)=msd(delt) &
                  +(atom_position(i,1)-x0(t,i,1))**2 &
                  +(atom_position(i,2)-x0(t,i,2))**2 &
                  +(atom_position(i,3)-x0(t,i,3))**2

            msd_x(delt)=msd_x(delt) + (atom_position(i,1)-x0(t,i,1))**2
            msd_y(delt)=msd_y(delt) + (atom_position(i,2)-x0(t,i,2))**2
            msd_z(delt)=msd_z(delt) + (atom_position(i,3)-x0(t,i,3))**2

            !!! write(*,'("atom index:",I8," delt:",I8," msd:",F23.15)') i,delt,msd(delt)

          end do

       end if

    end do

  end subroutine calc_msd
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output mean square displacement
!
  subroutine output_msd

    implicit none

    integer :: i
    real(8) :: time

    character(len=*),parameter :: filename_02="output_msd.txt"



    open(12,file=filename_02,status='replace')


    !!! write(*,'("***** MSD *****")')

    do i=1,tmax

       time=dt*(i-0.5)

       msd(i)=msd(i)/(natom*ntime(i))
       msd_x(i)=msd_x(i)/(natom*ntime(i))
       msd_y(i)=msd_y(i)/(natom*ntime(i))
       msd_z(i)=msd_z(i)/(natom*ntime(i))

       !!! write(*,'(I8, F23.15)') i, msd(i)
       write(12,'(5F23.15)') time,msd(i),msd_x(i),msd_y(i),msd_z(i)

    end do

    close(12) 

  end subroutine output_msd
!
end program calculate_msd_main
