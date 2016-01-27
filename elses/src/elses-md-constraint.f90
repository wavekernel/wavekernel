!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_constraint
!
   use elses_mod_ctrl,       only : i_verbose
   use elses_mod_file_io,    only : vacant_unit
!
   integer, parameter   :: DOUBLE_PRECISION=8
   real(DOUBLE_PRECISION), allocatable :: cell_parameter_store(:)

   public :: elses_md_constraint
!
 contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MD constraint motion
!
  subroutine elses_md_constraint
!
!   use elses_mod_sim_cell,   only : ax, ay, az
    use elses_mod_md_dat,     only : itemd
    use elses_mod_txp,        only : txp, typ, tzp
    use elses_mod_sim_cell,   only : noa
!
    implicit none
    logical :: file_exist
    integer :: line_count_max
    character(len=64) :: file_name
    character(len=32)  :: component
    integer :: iunit,ierr
    integer :: number_of_groups, number_of_members, jj 
    integer :: grp_index, atm_index
    integer :: noa2
    integer, allocatable :: member_of_group(:)
!
    real(DOUBLE_PRECISION) :: average_position
    real(DOUBLE_PRECISION) :: ddave, ddshift
!
!
    file_name='input-md-constraint.txt'
!
    if (i_verbose >= 1) then
       write(*,*)'@@  elses_md_constraint'
    endif   
!
    inquire (file=trim(file_name), exist=file_exist)
    if (i_verbose >= 1) write(*,*)'file_exist=',file_exist
!
    if (file_exist .eqv. .false.) then
      if (i_verbose >= 1) write(*,*)'The input file does not exist'
      return
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Atom movement for the constraint
!
    iunit=vacant_unit()
    open(iunit, file=file_name,status='old')
!
    read(iunit,*)noa2
    if (noa2 /= noa) then
      write(*,*)'Unsupported : noa2=',noa2
      stop
    endif   
!
    read(iunit,*)number_of_groups
!
    if (i_verbose >= 1) then
      write(*,*)'number of groups =', number_of_groups
    endif  
!
    do grp_index=1,number_of_groups
!
      read(iunit,*)component, number_of_members, average_position
      write(*,*)'component         =', trim(component)
      write(*,*)'number of members =', number_of_members
      write(*,*)'average position  =', average_position
!
      if (trim(component) /= 'z') then
        write(*,*)'unsupported format:',trim(component)
        stop
      endif   
!
      if ((number_of_members <= 0) .or. (number_of_members > noa)) then
        write(*,*)'unsupported format:number_of_members=',number_of_members
        stop
      endif   
!
      allocate(member_of_group(number_of_members),stat=ierr)
      if (ierr  /= 0) then
        write(*,*)'alloc. error'
        stop
      endif   
!
      do jj=1, number_of_members
        read(iunit,*) atm_index
        write(*,*)'   atm_index =',atm_index
        member_of_group(jj)=atm_index
        if ((atm_index <= 0) .or. (atm_index > noa)) then
          write(*,*)'unsupported format:atm_index=',atm_index
          stop
        endif   
      enddo   
!
      ddave=0.0d0
      do jj=1, number_of_members
        ddave=ddave+tzp(member_of_group(jj))/number_of_members
      enddo   
      write(*,'(a,2f16.10)')'befor: ddave, ave_pos =',ddave, average_position
!
      ddshift=ddave-average_position
      if (dabs(ddshift) > 0.1d0) then
        write(*,*)'Error(md constraint):shift=',ddshift
        stop
      endif   
!
      do jj=1, number_of_members
        tzp(member_of_group(jj))=tzp(member_of_group(jj))-ddshift
      enddo   
!
      ddave=0.0d0
      do jj=1, number_of_members
        ddave=ddave+tzp(member_of_group(jj))/number_of_members
      enddo   
      write(*,'(a,2f16.10)')'after: ddave, ave_pos =',ddave, average_position
!
      if (dabs(ddave-average_position) > 1.0d-10) then
        write(*,*)'Error(md constraint):=',ddave, average_position
        stop
      endif   
!
      deallocate(member_of_group)
!
    enddo
!   
    call elses_gene_tx
!
!   write(*,*)'stop manually'
!   stop
!
    close(iunit)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  end subroutine elses_md_constraint
!
end module M_md_constraint
