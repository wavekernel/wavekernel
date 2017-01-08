!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_set_s_as_i
!
  private
  public set_s_mat_as_i
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine set_s_mat_as_i(result_flag)
    use M_config,    only : config !(unchanged) 
    use M_qm_domain, only : dsij !(CHANGED!!!)
    use M_qm_domain, only : noav, njsd, jsv4jsd, nval, atm_element !(unchanged)
    implicit none
    logical, intent(out) :: result_flag
    integer :: ict4h
    integer :: jsv1, jsv2, jsd, ja1, ja2
!
    result_flag = .false.
    if (.not. config%calc%genoOption%set_S_as_I) return
!
    if (.not. allocated(dsij)) then 
      write(*,*)'INFO-XML: set_S_as_I tag is ignored, because DSIJ is not allocated'
      return 
    endif   
!
    result_flag = .true.
    write(*,*)'@@ Set S as I:INFO: S is set to be I'
!
    dsij(:,:,:,:)=0.0d0
    ict4h=1
    do jsv2=1,noav
      do jsd=1,njsd(jsv2,ict4h)
        jsv1=jsv4jsd(jsd,jsv2)
        do ja2=1,nval(atm_element(jsv2))
          do ja1=1,nval(atm_element(jsv1))
            if ((jsv2 == jsv1) .and. (ja1 == ja2)) then
              dsij(ja1,ja2,jsd,jsv2)=1.0d0
            endif   
          enddo
        enddo   
      enddo   
    enddo   
!
  end subroutine set_s_mat_as_i
!
end module M_qm_set_s_as_i


