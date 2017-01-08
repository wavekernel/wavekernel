!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_cohp_dstm
!
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_file_is_set, log_unit   !(unchanged)
  implicit none
!  
  private
!
! Public routines
  public calc_cohp_dstm
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the partial trace
!       such as local energy for given atom
!
  subroutine calc_cohp_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&             dm_loc, mat_wrk, result_value, jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset, id_of_my_omp_thread)
!    
    use M_qm_domain,      only : jsv4jsd, noav, njsd, nval, atm_element !(unchanged)
    use M_cohp_dstm_plot, only : calc_cohp_dstm_plot  !(routine)
!
    implicit none
    integer,                   intent(in)  :: atm_index, dst_atm_index, orb_index, m_int
    integer,                   intent(in)  :: jsv4jsk(:), jjkset(:)
    real(DOUBLE_PRECISION),    intent(out) :: result_value 
    real(DOUBLE_PRECISION),    intent(in)  :: dm_loc(:)
    real(DOUBLE_PRECISION),    intent(in)  :: mat_wrk(:,:,:,:)         !! hamiltonian
    integer,                   intent(in)  :: booking_list_dstm(:,:)
    integer,                   intent(in)  :: booking_list_dstm_len(:)
    integer,                   intent(in)  :: id_of_my_omp_thread
!
    integer :: jsv2, ja2, ict4h, jj, jsd1, jsv1, jjkset1, ja1, jjk1
    integer :: alloc_size_dst, nval_max
    real(DOUBLE_PRECISION) :: ddsum
    integer :: jsk1, jsk2
!
    logical :: debug_mode
!
    integer :: ierr
    integer :: size1, size2, size3
    real(DOUBLE_PRECISION), allocatable :: dm_loc_mat(:,:,:)
    real(DOUBLE_PRECISION), allocatable :: cohp_loc(:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocate the arrays
!
    size1=size(mat_wrk,1)
    size2=size(mat_wrk,2)
    size3=size(mat_wrk,3)
!
    allocate(dm_loc_mat(size1,size2,size3), stat=ierr)
    if (ierr /= 0) then
      stop 'Alloc. Error (calc_partial_cohp_dstm):dm_loc_mat'
    endif
    dm_loc_mat(:,:,:)=0.0d0
!   
    allocate(cohp_loc(size3), stat=ierr)
    if (ierr /= 0) then
      stop 'Alloc. Error (calc_partial_cohp_dstm):cohp_loc'
    endif
    cohp_loc(:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set the local variables
!
    if (i_verbose >= 2) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif   
!
    jsv2=atm_index
    ja2=orb_index
    ict4h=1
!
    jsk2=1
    if (atm_index /= jsv4jsk(jsk2)) then
      write(*,*)'ERROR(calc_partial_trace_dstm)'
      write(*,*)'atm_index,jsv4jsk(jsk2) =',atm_index, jsv4jsk(jsk2)
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Trivial checking for the input quantities (non-essential)
!
    if (size(dm_loc,1) /= m_int) then
       write(*,*)'ERROR:convert_dm_and_edm:size(dm_loc,1)= ',size(dm_loc,1)
    endif
!
    if ((jsv2 <=0 ) .or. (jsv2 > noav)) then
       write(*,*)'ERROR:convert_dm_and_edm'
       write(*,*)' atm_index, orb_index=',jsv2, ja2
    endif   
!
    if ((ja2 <=0 ) .or. (ja2 > nval(atm_element(jsv2)))) then
       write(*,*)'ERROR:convert_dm_and_edm'
       write(*,*)' atm_index, orb_index=',jsv2, ja2
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Genarate DM as the matrix form of 'dm_loc_mat'
!
    jj=0
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      jsv1=jsv4jsk(jsk1)
      jjkset1=jjkset(jsk1)
      do ja1=1,nval(atm_element(jsv1))
        jjk1=jjkset1+ja1
        jj=jj+1
        if (jj > m_int) then
          write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
          write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
          stop
        endif
        dm_loc_mat(ja1,ja2,jsd1) =  dm_loc(jj)
      enddo
    enddo   
!
    if (m_int /= jj) then
       write(*,*)'ERROR(calc_partial_cohp_dstm):jj,m_int=',jj,m_int
       stop
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculate the partial COHP
!
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      jsv1=jsv4jsk(jsk1)
      jjkset1=jjkset(jsk1)
      ddsum=0.0d0
      do ja1=1,nval(atm_element(jsv1))
        ddsum=ddsum+dm_loc_mat(ja1,ja2,jsd1)*mat_wrk(ja1,ja2,jsd1,jsk2)
!       write(*,'(a,2i10,2f20.10)')'ja1,jsd1, dm, H=', ja1, jsd1, & 
!&                                   dm_loc_mat(ja1,ja2,jsd1), mat_wrk(ja1,ja2,jsd1,jsk2)
      enddo
      cohp_loc(jsd1)=2.0d0*ddsum
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Calculate the local energy
!
    ddsum=0.0d0
    do jsd1=1,booking_list_dstm_len(jsk2)
      ddsum=ddsum+cohp_loc(jsd1)
!     write(*,'(a,i10,f20.10)')' jsd1, cohp_loc=',jsd1, cohp_loc(jsd1)
    enddo   
    result_value=ddsum
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ plot the COHP into the file
!
    call calc_cohp_dstm_plot(atm_index,orb_index, cohp_loc,jsv4jsk,booking_list_dstm, & 
&                            booking_list_dstm_len,id_of_my_omp_thread)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    deallocate(dm_loc_mat, stat=ierr)
    if (ierr /= 0) then
      stop 'Dealloc. Error (calc_partial_cohp_dstm):dm_loc_mat'
    endif
!
    deallocate(cohp_loc, stat=ierr)
    if (ierr /= 0) then
      stop 'Dealloc. Error (calc_partial_cohp_dstm):cohp_loc'
    endif
!
  end subroutine calc_cohp_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_cohp_dstm

