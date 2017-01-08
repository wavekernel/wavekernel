!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_partial_trace_dstm
!
!
  use M_qm_domain,        only : i_verbose, DOUBLE_PRECISION !(unchanged)
  use M_io_dst_write_log, only : log_file_is_set, log_unit   !(unchanged)
  implicit none
!  
  private
!
! Public routines
  public set_interac_list_dstm
  public calc_partial_trace_dstm
  public set_booking_list_rev1_dstm
  public calc_partial_tb0_force_dstm
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the interaction basis list
!
  subroutine set_interac_list_dstm(jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset, m_int, interaction_list)
!    
    use M_qm_domain,   only :  noav, nval, atm_element !(unchanged)
!
    implicit none
    integer,                intent(inout)  :: m_int
    integer,                   intent(in)  :: jsv4jsk(:), jjkset(:)
    integer,                   intent(in)  :: booking_list_dstm(:,:)
    integer,                   intent(in)  :: booking_list_dstm_len(:)
    integer,                     optional  :: interaction_list(:)
!
    integer :: jj, jsd1, jsv1, jjkset1, ja1, jjk1
    integer :: jsk1, jsk2
!
!   logical, parameter :: debug_mode = .true.
!
    jsk2=1

!
    if ( present(interaction_list) ) then
      if (m_int /= size(interaction_list,1) ) then
        write(*,*) 'Stop:m_int, size(interaction_list,1) =', m_int, size(interaction_list,1)
        stop
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    jj=0
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      jsv1=jsv4jsk(jsk1)
      jjkset1=jjkset(jsk1)
      do ja1=1,nval(atm_element(jsv1))
        jjk1=jjkset1+ja1
        jj=jj+1
        if ( present(interaction_list) ) then
          if (jj > m_int) then
            write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
            write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
            stop
          endif
          interaction_list(jj)=jjk1
        endif
      enddo
    enddo   
!
    if ( present(interaction_list) ) then
      if (m_int /= jj) then
        write(*,*)'ERROR(convert_dm_and_edm):jj,m_int=',jj,m_int
        stop
     endif   
    else
     m_int=jj
    endif
!
  end subroutine set_interac_list_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the interaction basis list
!
  subroutine set_interac_list_dstm_old(dst_atm_index, atm_index, orb_index, m_int,  &
&                  jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset, interaction_list)
!    
    use M_qm_domain,   only :  noav, nval, atm_element !(unchanged)
!
    implicit none
    integer,                   intent(in)  :: atm_index, dst_atm_index, orb_index
    integer,                intent(inout)  :: m_int
    integer,                   intent(in)  :: jsv4jsk(:), jjkset(:)
    integer,                   intent(in)  :: booking_list_dstm(:,:)
    integer,                   intent(in)  :: booking_list_dstm_len(:)
    integer,                     optional  :: interaction_list(:)
!
    integer :: jsv2, ja2, ict4h, jj, jsd1, jsv1, jjkset1, ja1, jjk1
    integer :: alloc_size_dst, nval_max
    real(DOUBLE_PRECISION) :: ddd1, ddd2, ddsum
    integer :: jsk1, jsk2
!
    logical :: debug_mode
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
    if ( present(interaction_list) ) then
      if (m_int /= size(interaction_list,1) ) then
        write(*,*) 'Stop:m_int, size(interaction_list,1) =', m_int, size(interaction_list,1)
        stop
      endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    jj=0
    ddsum=0.0d0
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      jsv1=jsv4jsk(jsk1)
      jjkset1=jjkset(jsk1)
      do ja1=1,nval(atm_element(jsv1))
        jjk1=jjkset1+ja1
        jj=jj+1
        if ( present(interaction_list) ) then
          if (jj > m_int) then
            write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
            write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
            stop
          endif
          interaction_list(jj)=jjk1
        endif
      enddo
    enddo   
!
!   write(*,*)'partial_trace:r=',atm_index, dst_atm_index, orb_index, m_int, result_value
!
    if ( present(interaction_list) ) then
      if (m_int /= jj) then
        write(*,*)'ERROR(convert_dm_and_edm):jj,m_int=',jj,m_int
        stop
     endif   
    endif
!
  end subroutine set_interac_list_dstm_old
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the partial trace
!       such as local energy for given atom
!
  subroutine calc_partial_trace_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&                      dm_loc, mat_wrk, result_value, jsv4jsk, booking_list_dstm, booking_list_dstm_len, jjkset)
!    
    use M_qm_domain,   only : jsv4jsd, noav, njsd, nval, atm_element !(unchanged)
!
    implicit none
    integer,                   intent(in)  :: atm_index, dst_atm_index, orb_index, m_int
    integer,                   intent(in)  :: jsv4jsk(:), jjkset(:)
    real(DOUBLE_PRECISION),    intent(out) :: result_value 
    real(DOUBLE_PRECISION),    intent(in)  :: dm_loc(:)
    real(DOUBLE_PRECISION),    intent(in)  :: mat_wrk(:,:,:,:)         !! ham_tb0_dstm or overlap_dstm
    integer,                   intent(in)  :: booking_list_dstm(:,:)
    integer,                   intent(in)  :: booking_list_dstm_len(:)
!
    integer :: jsv2, ja2, ict4h, jj, jsd1, jsv1, jjkset1, ja1, jjk1
    integer :: alloc_size_dst, nval_max
    real(DOUBLE_PRECISION) :: ddd1, ddd2, ddsum
    integer :: jsk1, jsk2
!
    logical :: debug_mode
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
!
!   write(*,*)'partial_trace=',atm_index, dst_atm_index, orb_index, m_int
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
    jj=0
    ddsum=0.0d0
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      jsv1=jsv4jsk(jsk1)
!     jsk1=jsk4jsv(jsv1)
!     if (jsk1 == 0) then
!       write(*,*)'ERROR(convert_dm_and_edm):jsk1=',jsk1
!       stop
!     endif
      jjkset1=jjkset(jsk1)
      do ja1=1,nval(atm_element(jsv1))
        jjk1=jjkset1+ja1
        jj=jj+1
        if (jj > m_int) then
          write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
          write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
          stop
        endif
        if (debug_mode) then
          if ( (jj < lbound(dm_loc,1)) .or. (jj > ubound(dm_loc,1)) ) then
            write(*,*)'ERROR(dm_loc):jj, jsv2, ja2, jsv1, ja1=',jj, jsv2, ja2, jsv1, ja1
            stop
          endif
          if ( (ja1 < lbound(mat_wrk,1)) .or. (ja1 > ubound(mat_wrk,1)) ) then
            write(*,*)'ERROR(mat_loc,1): jsv2, ja2, jsv1, ja1=',jsv2, ja2, jsv1, ja1
            stop
          endif
          if ( (ja2 < lbound(mat_wrk,2)) .or. (ja2 > ubound(mat_wrk,2)) ) then
            write(*,*)'ERROR(mat_loc,2): jsv2, ja2, jsv1, ja1=',jsv2, ja2, jsv1, ja1
            stop
          endif
          if ( (jsd1 < lbound(mat_wrk,3)) .or. (jsd1 > ubound(mat_wrk,3)) ) then
            write(*,*)'ERROR(mat_loc,3): jsv2, ja2, jsv1, ja1=',jsv2, ja2, jsv1, ja1
            stop
          endif
        endif
        ddd1 =  dm_loc(jj)
        ddd2 = mat_wrk(ja1,ja2,jsd1,jsk2)
        ddsum = ddsum + ddd1*ddd2 
      enddo
    enddo   
    result_value=2.0d0*ddsum   ! Two : para-spin factor
!
!   write(*,*)'partial_trace:r=',atm_index, dst_atm_index, orb_index, m_int, result_value
!
    if (m_int /= jj) then
       write(*,*)'ERROR(convert_dm_and_edm):jj,m_int=',jj,m_int
       stop
    endif   
!
  end subroutine calc_partial_trace_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the reversed booking list
!
  subroutine set_booking_list_rev1_dstm(booking_list_rev1_dstm, booking_list_dstm, booking_list_dstm_len, num_atom_proj)
    implicit none
    integer, intent(out):: booking_list_rev1_dstm(:)
    integer, intent(in) :: booking_list_dstm(:,:)
    integer, intent(in) :: booking_list_dstm_len(:)
    integer, intent(in) :: num_atom_proj
    integer :: jsk1, jsk2, jsd2 , jsk_seed
!
    jsk_seed=1
!
    booking_list_rev1_dstm(:)=0
    do jsk1=1, num_atom_proj
      do jsd2=1,booking_list_dstm_len(jsk1)
        jsk2=booking_list_dstm(jsd2,jsk1)
        if ( jsk2 == jsk_seed ) then
          booking_list_rev1_dstm(jsk1)=jsd2
          exit
        endif
      enddo
    enddo

  end subroutine set_booking_list_rev1_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the TB force 
!
  subroutine calc_partial_tb0_force_dstm(atm_index, dst_atm_index, orb_index, m_int, & 
&                  dm_loc, edm_loc, d_ham_tb0_dstm, d_overlap_dstm,                 &
&                  atm_tb0_force_dstm, jsv4jsk, booking_list_dstm, booking_list_rev1_dstm, booking_list_dstm_len, jjkset)
!
!   use elses_mod_md_dat,   only : itemd
!   use M_md_dst,      only : myrank
!   use M_qm_dst_proj_cell,    only : dst_atm_list, dst_atm_list_rev, len_dst_atm_list !(unchnaged)
!   use M_qm_domain, only : noav, atm_element, nval, jsv4jsd, njsd !(unchanged)
!   use M_qm_geno_dst, only : d_ham_tb0_dst, d_overlap_dst
!   use M_qm_domain, only : ddhij, ddsij              !(unchanged)
!   use elses_mod_phys_const, only : para_spin_factor !(parameter)
    use M_qm_domain, only : atm_element, nval !(unchanged)
    implicit none
    real(DOUBLE_PRECISION),    intent(out) :: atm_tb0_force_dstm(:,:)
    integer,                   intent(in)  :: atm_index, dst_atm_index, orb_index, m_int
    integer,                   intent(in)  :: jsv4jsk(:), jjkset(:)
    real(DOUBLE_PRECISION),    intent(in)  ::  dm_loc(:)
    real(DOUBLE_PRECISION),    intent(in)  :: edm_loc(:)
    real(DOUBLE_PRECISION),    intent(in)  :: d_ham_tb0_dstm(:,:,:,:,:)
    real(DOUBLE_PRECISION),    intent(in)  :: d_overlap_dstm(:,:,:,:,:)
    integer,                   intent(in)  :: booking_list_dstm(:,:)
    integer,                   intent(in)  :: booking_list_rev1_dstm(:)
    integer,                   intent(in)  :: booking_list_dstm_len(:)
    logical, parameter :: debug_mode = .true.
!
    integer :: jsv2, ict4h, jj, jsv1, jjkset1, jjk1
    integer :: alloc_size_dst, nval_max
    real(DOUBLE_PRECISION) :: ddd1, ddd2, ddsum
    integer :: jsk1, jsk2

    integer :: jsd1, jsd2, nss1, nss2, ja1, ja2
!   integer :: dst_atm_index, dst_atm_index1
!   integer :: size1, size2
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set the local variables
!
    jsv2=atm_index
    ja2=orb_index
    jsk2=1
!
    if (atm_index /= jsv4jsk(jsk2)) then
      write(*,*)'ERROR(calc_partial_trace_dstm)'
      write(*,*)'ERROR(calc_partial_trace_dstm)'
      stop
    endif
!
    atm_tb0_force_dstm(:,:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Add the 'first-derivative' term
!
    jj=0
    ddsum=0.0d0
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      if (jsk1 == jsk2) then
        jj=jj+nval(atm_element(jsv2))
        cycle
      endif
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
        ddd1 =   dm_loc(jj)    !  dm(ja1, ja2, jsd1, jsk2)
        ddd2 =  edm_loc(jj)    ! edm(ja1, ja2, jsd1, jsk2)
        atm_tb0_force_dstm(:,jsk2)=atm_tb0_force_dstm(:,jsk2) &
&            - 2.0d0*ddd1*d_ham_tb0_dstm(:,ja1,ja2,jsd1,jsk2) &
&            + 2.0d0*ddd2*d_overlap_dstm(:,ja1,ja2,jsd1,jsk2) ! Two : para-spin factor
        if (i_verbose >= 1) then
          if (jsv2 < 3) then
            write(*,'(a,4i7,8f10.5)')'d_ham_tb=',ja1,ja2, jsv1, jsv2, &
&             d_ham_tb0_dstm(:,ja1,ja2,jsd1,jsk2),d_overlap_dstm(:,ja1,ja2,jsd1,jsk2),ddd1,ddd2
          endif
        endif  
      enddo
    enddo
    if (jj /= m_int) then
      write(*,*)'ERROR(alc_partial_tb0_force_dstm):jj,m_int=',jj,m_int
      stop
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Add the 'second-derivative' term
!
    jj=0
    ddsum=0.0d0
    do jsd1=1,booking_list_dstm_len(jsk2)
      jsk1=booking_list_dstm(jsd1, jsk2)
      if (jsk1 == jsk2) then
        jj=jj+nval(atm_element(jsv2))
        cycle
      endif
      jsd2=booking_list_rev1_dstm(jsk1)
      if (debug_mode) then 
        if (jsd2 < 1) then
          write(*,*)'ERROR(calc_partial_tb0_force_dstm:a):jsd2,jsk1,jsk2=',jsd2,jsk1,jsk2
          stop
        endif
        if (jsd2 > booking_list_dstm_len(jsk1)) then
          write(*,*)'ERROR(calc_partial_tb0_force_dstm:b):jsd2,jsk1,jsk2=',jsd2,jsk1,jsk2
          stop
        endif
        if (booking_list_dstm(jsd2, jsk1) /= jsk2) then
          write(*,*)'ERROR(calc_partial_tb0_force_dstm:c):jsd2,jsk1,jsk2=',jsd2,jsk1,jsk2
          stop
        endif
      endif
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
        ddd1 =   dm_loc(jj)    !  dm(ja1, ja2, jsd1, jsk2)
        ddd2 =  edm_loc(jj)    ! edm(ja1, ja2, jsd1, jsk2)
        atm_tb0_force_dstm(:,jsk1)=atm_tb0_force_dstm(:,jsk1) &
&            - 2.0d0*ddd1*d_ham_tb0_dstm(:,ja2,ja1,jsd2,jsk1) &
&            + 2.0d0*ddd2*d_overlap_dstm(:,ja2,ja1,jsd2,jsk1) ! Two : para-spin factor
      enddo
    enddo
    if (jj /= m_int) then
      write(*,*)'ERROR(alc_partial_tb0_force_dstm):jj,m_int=',jj,m_int
      stop
    endif
!
  end subroutine calc_partial_tb0_force_dstm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_partial_trace_dstm
