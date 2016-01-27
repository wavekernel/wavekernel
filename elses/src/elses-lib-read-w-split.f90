!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_read_w_split
!
  implicit none
!
  private
  public :: read_w_split
  public :: test_read_w_split
  
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ Read a line with splitting
!
  subroutine read_w_split(chara_wrk, separator, number_of_items, max_word_length_set, chara_items)
    implicit none
    character(len=*),    intent(in) :: chara_wrk
    character(len=*),    intent(in) :: separator
    integer,          intent(inout) :: number_of_items
    integer,          intent(inout) :: max_word_length_set
    character(len=*),      optional :: chara_items(:)
    integer :: i_sta, i_end, i_end_word
    integer,             parameter  :: loop_max = 100
    integer :: i_try, j
    integer :: item_counter, max_word_length
!
    i_sta = 1
    i_end = len_trim(chara_wrk)
!
    item_counter=0
    max_word_length=0
    do i_try=1, loop_max
      if (i_try == loop_max) then
        write(*,*)'ERROR(read_w_split):i_try,loop_max=', i_try, loop_max
        stop
      endif
!     write(*,*)' i_sta, i_end=', i_sta, i_end
      if ( i_sta > i_end-1 ) exit
      j=index(chara_wrk(i_sta:i_end), separator)
!
      select case(j)
        case (0) 
          i_end_word = i_end
        case (1) 
          i_sta = i_sta+1
          cycle
        case default
          i_end_word = i_sta+j-2
      end select
!
!     write(*,*)'i_sta, i_end_word=', i_sta, i_end_word
      item_counter = item_counter +1
      max_word_length = max(max_word_length, i_end_word-i_sta+1)
!     write(*,'(a,a,a)') 'item = "', trim(chara_word), '"'
!
      if (present(chara_items)) then
        if (item_counter > size(chara_items,1)) then
          write(*,*)'ERROR(read_w_split):item_counter=',item_counter
          stop
        endif
        if (max_word_length > len(chara_items)) then
          write(*,*)'ERROR:max_word_length, len(chara_items)=', max_word_length, len(chara_items)
          stop
        endif
        chara_items(item_counter)=trim(adjustl(chara_wrk(i_sta:i_end_word)))
      endif
!
      i_sta = i_end_word+2
!
    enddo
!
    if (.not. present(chara_items)) then
      number_of_items=item_counter
      max_word_length_set = max_word_length
    endif
!
  end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ Test routine
!
  subroutine test_read_w_split(log_unit)
    implicit none
    integer, intent(in) :: log_unit
!
    character(len=*), parameter :: separator=' '
    integer :: number_of_items, max_word_length
    character(len=1024), allocatable :: chara_items(:)
    integer :: lu, ierr, i
!
    character(len=1024) :: chara_test
    character(len=1024) :: chara_items_result(4)
!
    chara_test='  Dog Book  Apple    Supercalifragilisticexpialidocious '
    chara_items_result(1)='Dog'
    chara_items_result(2)='Book'
    chara_items_result(3)='Apple'
    chara_items_result(4)='Supercalifragilisticexpialidocious'
!
    lu = log_unit
!
    if (lu > 0) write(lu,'(a)') '@@ test_read_w_split'
    if (lu > 0) write(lu,'(a,a,a)') '  input text = "', trim(chara_test), '"'
!
    number_of_items=0
    call read_w_split(chara_test, separator, number_of_items, max_word_length)
    if (lu > 0) then
      write(lu,*) ' number_of_items = ', number_of_items
      write(lu,*) ' max_word_length = ', max_word_length
    endif
!
    if (max_word_length > len(chara_items)) then
      write(*,*)'ERROR:max_word_length, len(chara_items)=', max_word_length, len(chara_items)
      stop
    endif

    allocate(chara_items(number_of_items), stat=ierr)
    if (ierr /= 0) stop 'Alloc. Error: chara_items'
!
    call read_w_split(chara_test, separator, number_of_items, max_word_length,chara_items)
!
    if (number_of_items /= 4) then
      write(*,*)'ERROR(test_read_w_split):number_of_items=', number_of_items
      stop
    endif
!
    if (max_word_length /= 34) then
      write(*,*)'ERROR(test_read_w_split):max_word_length=', max_word_length
      stop
    endif
!
    do i=1,number_of_items
      if (lu > 0) write(lu,'(a,a,a)') '  item = "', trim(chara_items(i)), '"'
      if (trim(chara_items(i)) /= trim(chara_items_result(i))) then
        write(*,*)'ERROR(test_read_w_split)'
        stop
      endif
    enddo
!
    if (lu > 0) write(lu,'(a)') '....end without error: test_read_w_split'
!
!   stop 'Stop manually'
!
  end subroutine
!
end module M_lib_read_w_split
