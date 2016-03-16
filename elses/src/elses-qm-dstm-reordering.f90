!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_dstm_reordering
!
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_file_is_set, log_unit   !(unchanged)
  implicit none
!  
  private
!
  public reorder_booked_atoms
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the partial trace
!       such as local energy for given atom
!
!  CAUTION : Do NOT reorder the first element of the booking list !!
!            The first element is one the 'center' atom
!            The booking list is that for the neighbor atoms of the 'center' atom.
!
  subroutine reorder_booked_atoms(center_atom_id, booking_list)
!    
    use M_md_get_distance, only : get_vector_for_atom_pair !(routine)
    use M_qm_domain,      only : noav   !(unchanged)
    use M_lib_sort,       only : heap_sort => dhsort !(routine)
!
    implicit none
    logical, parameter                        :: use_heap_sort = .true.
    logical, parameter                        :: debug_mode    = .false.
!   logical, parameter                        :: debug_mode    = .true.
    integer,                   intent(in)     :: center_atom_id
    integer,                   intent(inout)  :: booking_list(:)
!
    integer, allocatable  :: booking_list_org(:)
    integer, allocatable  :: sort_list(:)
!
    real(DOUBLE_PRECISION), allocatable :: loc_coord(:,:) 
                  ! the 'local' coodinates in which the 'center' atom is located at (0,0,0)
    real(DOUBLE_PRECISION), allocatable :: distance(:) 
                  ! the distance of the member atoms from the 'center' atom
    real(DOUBLE_PRECISION), parameter :: eps = 1.0d-10
    integer :: n, ierr
!
    integer :: j, atom_id, i
!
    real(DOUBLE_PRECISION) :: vec(3)
    real(DOUBLE_PRECISION) :: diff
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Initial setting and backup the original booking list.
!
    n = size(booking_list,1)
!
    if (debug_mode) then
      if (n <= 0) then
        write(*,*)'ERROR(reorder_booked_atoms):n=',n
        stop
      endif
    endif
!
    allocate (sort_list(n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ALloc. error sort_list'
      stop
    endif
!
    allocate (booking_list_org(n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ALloc. error booking_list_old'
      stop
    endif
    booking_list_org(:)=booking_list
!
    if (debug_mode) then
      do j=1,n
        atom_id=booking_list_org(j)
        if ((atom_id < 1) .or. (atom_id > noav)) then
          write(*,*)'ERROR(reorder_booked_atoms) in booking list:',j,atom_id
          stop
        endif
      enddo
    endif
!
    allocate (loc_coord(3,n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ALloc. error loc_coord'
      stop
    endif
!
    allocate (distance(n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ALloc. error distance'
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Initial setting and backup the original booking list.
!
    write(*,*)'Rerodering:center, size=', center_atom_id, n
!
    do j=1,n
      atom_id=booking_list_org(j)
      call get_vector_for_atom_pair(center_atom_id, atom_id, vec(1), vec(2), vec(3))
      loc_coord(1:3,j)=vec(1:3)
      distance(j)=dot_product(vec, vec)
    enddo
!
    if (debug_mode) then
      do j=1,n
        write(*,*)'old order :', j, booking_list_org(j), distance(j)
      enddo
    endif
!
    if (use_heap_sort) then
      call heap_sort(n,distance(:),sort_list(:))
    else
      call sort_routine_test(distance(:),sort_list(:))
    endif
!
    if (debug_mode) then
      if ( sort_list(1) /= 1 ) then
        write(*,*) 'ERROR(reorder_booked_atoms):sort_list(1) =', sort_list(1)
        stop
      endif
    endif

    do j=1,n
      i=sort_list(j)
      atom_id=booking_list_org(i)
      booking_list(j)=atom_id
    enddo
!
    if (debug_mode) then
      if ( booking_list(1) /= center_atom_id ) then
        write(*,*) 'ERROR(reorder_booked_atoms):booking_list(1) =', booking_list(1), center_atom_id
        stop
      endif
    endif
!
    if (debug_mode) then
      do j=1,n
        atom_id=booking_list(j)
        call get_vector_for_atom_pair(center_atom_id, atom_id, vec(1), vec(2), vec(3))
        loc_coord(1:3,j)=vec(1:3)
        distance(j)=dot_product(vec, vec)
        if (j /= 1) then
          diff=distance(j)-distance(j-1)
          if (diff < - eps) then
            write(*,*)'ERROR(reorder_booked_atoms):j,diff=',j, diff
            stop
          endif
        endif
      enddo
      do j=1,n
        write(*,*)'new order :', j, booking_list(j), distance(j)
      enddo
    endif
!
!   stop 'STOP MANUALLY'
!
  end subroutine reorder_booked_atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @ Sort routine (NOT efficient)
!
  subroutine sort_routine_test(value, sorted_order)
    implicit none
    real(8),                   intent(in)     :: value(:)
    integer,                   intent(out)    :: sorted_order(:)
    real(8), allocatable  :: value_wrk(:)
    integer :: n, ierr
    integer :: j,i
!
    n=size(value,1)
!
    allocate (value_wrk(n), stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ALloc. error.: value work'
      stop
    endif
    value_wrk(:)=value(:)
!
    do j=1,n
      i=minloc(value_wrk,1)
!     write(*,*)' new order:', j, i, value_wrk(i)
      sorted_order(j)=i
      value_wrk(i)=huge(0.0d0)
    enddo
!
    do j=1,n
      i=sorted_order(j)
      write(*,*)' sorted order:', j, i, value(i)
    enddo
!
!   stop
!
  end subroutine sort_routine_test
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_dstm_reordering


