!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_ini_dst
!
  use M_config,       only : config ! CHANGED ONLY FOR 
!
  private
  public ini_dst_set
!
contains
!
  subroutine ini_dst_set
    use M_qm_domain,        only : nval
    implicit none
    integer :: lu, i_v 
    integer :: j
    logical :: const_num_orbital
    integer :: nval_max
    character(len=20) :: mat_vec_type
!
    lu = config%calc%distributed%log_unit
    i_v= config%option%verbose
!
    if (i_v >= 1) then
      if (lu > 0) then
        write(lu,*)'@@ INI_SET_DST'
        write(lu,*)'size(nval,1) = ', size(nval,1)
      endif
    endif
!
    nval_max = maxval(nval)
!
    const_num_orbital=.true.
    do j=1, size(nval,1)
      if (nval(j) /= nval_max) const_num_orbital=.false.
    enddo
!
    config%calc%distributed%mat_vec_const_num_orbital = const_num_orbital
    config%calc%distributed%mat_vec_max_num_orbital   = maxval(nval)
!
    config%calc%distributed%mat_vec_switch_bcrs = .false. ! default setting
!
    mat_vec_type=trim(adjustl(config%calc%distributed%mat_vec_mode))
    if (mat_vec_type == 'default') mat_vec_type = 'dstm'
    if (mat_vec_type == 'CRS')     mat_vec_type = 'crs'
    if (mat_vec_type == 'DENS')    mat_vec_type = 'dens'
    if (mat_vec_type == 'DSTM')    mat_vec_type = 'dstm'
    if (mat_vec_type == 'BCRS')    mat_vec_type = 'bcrs'
!
    if (mat_vec_type == 'bcrs') then
      if (config%calc%distributed%mat_vec_const_num_orbital) then
        mat_vec_type = 'dstm'
        config%calc%distributed%mat_vec_switch_bcrs = .true.
      else
        mat_vec_type = 'crs'
        if (i_v >= 1) then
          if (lu > 0) then
            write(lu,*)'WARNING:dst mat-vec mode is changed from BCRS to CRS '
          endif
        endif

      endif
    endif
!
    config%calc%distributed%mat_vec_mode=trim(mat_vec_type)
!
    if (i_v >= 1) then
      if (lu > 0) then
        write(lu,*)'mat_vec_const_num_orbital =', const_num_orbital
        write(lu,*)'mat_vec_max_num_orbital   =', maxval(nval)
      endif
    endif
!
!
!   stop 'Stop manually'
!
  end subroutine ini_dst_set
!
!
end module M_ini_dst
