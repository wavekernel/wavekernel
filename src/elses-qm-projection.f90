!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_projection
!
  use elses_mod_ctrl,       only : i_verbose
!
  implicit none
!  
  private
  public proj_init
  public proj_init_end
  public proj_get_mat_size
  public proj_get_list
  public get_interac_list_num_proj_atom
  public get_interac_list_num_proj
  public set_interac_list_proj
  public matvec_mul_proj
  public inner_product_with_s_proj
  public convert_dm_and_edm
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_init_end(imode)
    use M_wall_clock_time, only : get_system_clock_time
    implicit none
    integer :: imode
    real(8) :: time_wrk, time_wrk_previous ! work variable for measuring the time
!
    write(*,*)'@@ qm_solver_projection'
!
    call get_system_clock_time(time_wrk)
    time_wrk_previous=time_wrk
!
    call elses_set_param_ctl_kr
!
    call get_system_clock_time(time_wrk)
    write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(1)  =',time_wrk-time_wrk_previous
    time_wrk_previous=time_wrk
!
    if (imode == 1) then
      call elses_alloc_kry_glo(imode)
!
      call get_system_clock_time(time_wrk)
      write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(2)  =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      call elses_set_rcut_kry(imode)
!      ----> Set rcut_kry(1:noav) : projection radius
!              noak_str(1:noav) : atoms in the projection radius
!
      call get_system_clock_time(time_wrk)
      write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(3)  =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      call elses_alloc_jsv4jsk_str(imode)
!
      call get_system_clock_time(time_wrk)
      write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(4)  =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
      call elses_set_jsv4jsk_str(imode)
!       ----> Set jsv4jsk_str(1:noak_max,noav) : booking list for projection
!           NOTE: The 'self' atom is set to be the first element
!                    ( jsv4jsk_str(1,atm_index) = atm_index )
!
      call get_system_clock_time(time_wrk)
      write(*,'(a,f10.5)')'TIME:qm_solver_gkrylov_dst:init_end:(5)  =',time_wrk-time_wrk_previous
      time_wrk_previous=time_wrk
!
    endif  
!
    if (imode == 2) then
      call elses_alloc_kry_glo(imode)
      call elses_alloc_jsv4jsk_str(imode)
    endif
!  
  end subroutine proj_init_end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_init
    implicit none
    integer :: imode
!
    write(*,*)'@@ qm_solver_projection'
!
    call elses_set_param_ctl_kr
!
    imode=1
    call elses_alloc_kry_glo(imode)
!
    imode=1
    call elses_set_rcut_kry(imode)
!    ----> Set rcut_kry(1:noav) : projection radius
!              noak_str(1:noav) : atoms in the projection radius
!
    imode=1
    call elses_alloc_jsv4jsk_str(imode)
    call elses_set_jsv4jsk_str(imode)
!    ----> Set jsv4jsk_str(1:noak_max,noav) : booking list for projection
!           NOTE: The 'self' atom is set to be the first element
!                    ( jsv4jsk_str(1,atm_index) = atm_index )
!
  end subroutine proj_init
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_get_mat_size(atm_index,mat_size,num_atom_proj)
!
    use elses_arr_kry_glo,    only : jsv4jsk_str, noak_str ! (unchanged)
    use elses_mod_js4jsv,     only : js4jsv                ! (unchanged)
    use elses_mod_tx,         only : jsei                  ! (unchanged)
    use M_qm_domain,          only : noav, nval            ! (unchanged)
!
    implicit none
    integer, intent(in)  :: atm_index
    integer, intent(out) :: mat_size
    integer, intent(out) :: num_atom_proj
!
    integer :: jj, jsk, jsv, js, nss, nval2, jsk_self
!
!   write(*,*)'@@ qm_proj_get_mat_size:atm_index=',atm_index
!
    num_atom_proj=noak_str(atm_index)
!       ---> number of atom for projection
!
    jj=0
    jsk_self=0
    do jsk=1,num_atom_proj
      jsv=jsv4jsk_str(jsk,atm_index)
      if ((jsv .le. 0) .or. (jsv .gt. noav)) then
        write(6,*)'ERROR!(qm_proj_get_mat_size):jsv,noav=',jsv,noav
        stop
      endif   
      if (jsv == atm_index) then
        if (jsk_self == 0) then
          jsk_self=jsk
        else
          write(6,*)'ERROR!(qm_proj_get_mat_size):jsv,jsk_self=',jsv,jsk_self
          stop
        endif   
      endif   
      js=js4jsv(jsv)
      nss=jsei(js)
      nval2=nval(nss)
      jj=jj+nval2
    enddo  
    mat_size=jj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Check that The 'self' atom is set to be the first element
!              ( jsv4jsk_str(1,atm_index) = atm_index )
!
    if (jsk_self /= 1) then
      write(*,*)'ERROR!(qm_proj_get_mat_size):jsk_self=',jsk_self
      stop
    endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Trivial checking
!
    if ((mat_size <= 0) .or. (mat_size > num_atom_proj*100)) then
      write(*,*)'ERROR!(qm_proj_get_mat_size)'
      write(*,*)'atm_index, mat_size=',atm_index,mat_size
      stop
    endif   
!
  end subroutine proj_get_mat_size
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine proj_get_list(atm_index,jsv4jsk,jsk4jsv,jjkset)
!
    use elses_arr_kry_glo,    only : jsv4jsk_str, noak_str ! (unchanged)
    use elses_mod_js4jsv,     only : js4jsv                ! (unchanged)
    use elses_mod_tx,         only : jsei                  ! (unchanged)
    use M_qm_domain,          only : noav, nval            ! (unchanged)
!
    implicit none
    integer, intent(in)      :: atm_index
    integer, intent(out)     :: jsv4jsk(:)
    integer, intent(out)     :: jsk4jsv(:)
    integer, intent(out)     :: jjkset(:)
!
    integer :: num_atom_proj
    integer :: jsk, jsv, jj, js, nss, nval2
    integer :: nval_upper_limit
!
    nval_upper_limit=100
!
    num_atom_proj=size(jsv4jsk,1)
    if (num_atom_proj /= noak_str(atm_index)) then
      write(*,*)'ERROR(proj_get_list):num_atom_proj=',num_atom_proj
      stop
    endif
!
    jsv4jsk(1:num_atom_proj)=jsv4jsk_str(1:num_atom_proj,atm_index)
    jsk4jsv(:)=0
!
    do jsk=1, num_atom_proj
      jsv=jsv4jsk(jsk)
      if ((jsv .le. 0) .or. (jsv .gt. noav)) then
         write(6,*)'ERROR!(SETKRY):jsk,jsv=',jsk,jsv
         stop
      endif   
      jsk4jsv(jsv)=jsk
    enddo
!
    jj=0
    do jsk=1, num_atom_proj
      jjkset(jsk)=jj
      jsv=jsv4jsk(jsk)
      if ((jsv .le. 0) .or. (jsv .gt. noav)) then
        write(6,*)'ERROR!(SETKRY:1100):JSK,JSV=',jsk,jsv
        stop
      endif   
      js=js4jsv(jsv)
      nss=jsei(js)
      nval2=nval(nss)
      if ((nval2 .le. 0) .or. (nval2 .gt. nval_upper_limit)) then
        write(6,*)'ERROR!(SETKRY:1100):JS,NVAL2=',js,nval2
        stop
      endif   
      jj=jj+nval2
     enddo   
!
  end subroutine proj_get_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_interac_list_num_proj_atom(jsv2,m_int)
!
    use M_qm_domain,          only : jsv4jsd, njsd, atm_element, nval !(unchanged)
    implicit none
    integer,          intent(in)   :: jsv2
    integer,          intent(out)  :: m_int
!
    integer :: jsv1, ict4h
    integer :: jsd, num_of_bases
!
    ict4h=1
!
    num_of_bases=0
    do jsd=1,njsd(jsv2,ict4h)
      jsv1=jsv4jsd(jsd,jsv2)
      num_of_bases=num_of_bases+nval(atm_element(jsv1))
    enddo   
    m_int=num_of_bases
!
    if (m_int <= 0) then
      write(*,*)'ERROR(get_interac_list_num_proj):jsv2, m_int=',jsv2, m_int
    endif
!
  end subroutine get_interac_list_num_proj_atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_interac_list_num_proj(j_src,m_int)
!
    use elses_mod_orb2,       only : j2js,j2ja                        !(unchanged)
    use elses_mod_js4jsv,     only : jsv4js                           !(unchanged)
    use M_qm_domain,          only : jsv4jsd, njsd, atm_element, nval !(unchanged)
    implicit none
    integer,          intent(in)   :: j_src
    integer,          intent(out)  :: m_int
!
    integer :: js, jsv2, jsv1, ja, ja1, ict4h
    integer :: jsd, num_of_bases
!
    js=j2js(j_src)
    ja=j2ja(j_src)
    jsv2=jsv4js(js)
    ict4h=1
!
    call get_interac_list_num_proj_atom(jsv2,m_int)
!
  end subroutine get_interac_list_num_proj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_interac_list_proj(j_src,list_member,jsk4jsv,jjkset)
!
    use elses_mod_orb2,       only : j2js,j2ja                        !(unchanged)
    use elses_mod_js4jsv,     only : jsv4js                           !(unchanged)
    use M_qm_domain,          only : jsv4jsd, njsd, atm_element, nval !(unchanged)
    implicit none
    integer,          intent(in)   :: j_src
    integer,          intent(in)   :: jsk4jsv(:)
    integer,          intent(in)   :: jjkset(:)
    integer,          intent(out)  :: list_member(:)
!
    integer :: js, ja, jsv2, ict4h, m_int
    integer :: jsd1, jsv1, jsk1, jjkset1, ja2, jjk2, ja1, jjk1, jj
!
    js=j2js(j_src)
    ja=j2ja(j_src)
    jsv2=jsv4js(js)
    ict4h=1
!
    m_int=size(list_member,1)
    list_member(:)=0
!
    jj=0
    do jsd1=1,njsd(jsv2,ict4h)
      jsv1=jsv4jsd(jsd1,jsv2)
      jsk1=jsk4jsv(jsv1)
      if (jsk1 .ne. 0) then
        jjkset1=jjkset(jsk1)
        do ja1=1,nval(atm_element(jsv1))
          jjk1=jjkset1+ja1
          jj=jj+1
          if (jj > m_int) then
            write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
            write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
            stop
          endif
          list_member(jj)=jjk1
        enddo
       endif  
    enddo   
    if (m_int /= jj) then
       write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
       stop
    endif   
!    
  end subroutine set_interac_list_proj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with H
!   mat_kind='H' : (vect_out) = H (vect_in)
!   mat_kind='S' : (vect_out) = S (vect_in)
!
  subroutine matvec_mul_proj(vect_in,vect_out,mat_kind,jsv4jsk,jsk4jsv,jjkset)
!
    use M_qm_domain, only : noav, nval, dhij, dsij, njsd, jsv4jsd, atm_element  ! (unchanged)
    implicit none
    complex(8),       intent(in)  :: vect_in(:)
    complex(8),       intent(out) :: vect_out(:)
    integer,          intent(in)  :: jsv4jsk(:)
    integer,          intent(in)  :: jsk4jsv(:)
    integer,          intent(in)  :: jjkset(:)
    character(len=*),  intent(in)  :: mat_kind
    integer :: m, noak, noav2
    integer :: ierr, ict4h
    integer :: jsk2, jsv2, jjkset2, jsd1, jsv1, jsk1 
    integer :: jjkset1, ja2, jjk2, ja1, jjk1
    real(8) :: dbigd
!
    ict4h=1
    m=size(vect_in,1)
    noav2=size(jsk4jsv,1)
    noak=size(jsv4jsk,1)
!
    if (m /= size(vect_out,1)) stop 'Abort:Size mismatch of m'
    if (noav2 /= noav) stop 'Abort:Size mismatch of noav'
    if (noak /= size(jjkset,1)) stop 'Abort:Size mismatch of noak'
    if ((mat_kind /='H') .and. (mat_kind /='S')) stop 'Abort:Wrong mode'
!
    vect_out(:)=0.0d0
!
    if (mat_kind == 'H') then
      do jsk2=1,noak
        jsv2=jsv4jsk(jsk2)
        jjkset2=jjkset(jsk2)
        do jsd1=1,njsd(jsv2,ict4h)
          jsv1=jsv4jsd(jsd1,jsv2)
          jsk1=jsk4jsv(jsv1)
          if ((jsk1 < 0) .or. (jsk1 > noak)) then
            write(*,*)'ERROR(matvec_mul_proj)'
            write(*,*)'jsv1, jsk1=',jsv1,jsk1
            stop
          endif   
          if (jsk1 .ne. 0) then
            jjkset1=jjkset(jsk1)
            do ja2=1,nval(atm_element(jsv2))
              jjk2=jjkset2+ja2
              do ja1=1,nval(atm_element(jsv1))
                jjk1=jjkset1+ja1
                dbigd=dhij(ja1,ja2,jsd1,jsv2)
                vect_out(jjk2)=vect_out(jjk2)+dbigd*vect_in(jjk1)
              enddo
            enddo
          endif  
        enddo   
      enddo   
    endif
!
    if (mat_kind == 'S') then
      do jsk2=1,noak
        jsv2=jsv4jsk(jsk2)
        jjkset2=jjkset(jsk2)
        do jsd1=1,njsd(jsv2,ict4h)
          jsv1=jsv4jsd(jsd1,jsv2)
          jsk1=jsk4jsv(jsv1)
          if ((jsk1 < 0) .or. (jsk1 > noak)) then
            write(*,*)'ERROR(matvec_mul_proj)'
            write(*,*)'jsv1, jsk1=',jsv1,jsk1
            stop
          endif   
          if (jsk1 .ne. 0) then
            jjkset1=jjkset(jsk1)
            do ja2=1,nval(atm_element(jsv2))
              jjk2=jjkset2+ja2
              do ja1=1,nval(atm_element(jsv1))
                jjk1=jjkset1+ja1
                dbigd=dsij(ja1,ja2,jsd1,jsv2)
                vect_out(jjk2)=vect_out(jjk2)+dbigd*vect_in(jjk1)
              enddo
            enddo
          endif  
        enddo   
      enddo   
    endif
!
  end subroutine matvec_mul_proj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inner product with the j-th column of the S matrix
!     within the 'compressed' format
!
  subroutine inner_product_with_s_proj(vect_in,result_value,j_src)
!
    use elses_mod_orb2,   only : j2js,j2ja                        !(unchanged)
    use elses_mod_js4jsv, only : jsv4js                           !(unchanged)
    use M_qm_domain,      only : jsv4jsd, njsd, atm_element, nval, &
&                                 dsij  !(unchanged)
    implicit none
    integer,          intent(in)   :: j_src
    real(8),          intent(in)   :: vect_in(:)
    real(8),          intent(out)  :: result_value
!
    integer :: m_in, js, jsv2, jsv1, ja2, ja1, ict4h, jsd, jj
    real(8) :: tmp
!
    m_in=size(vect_in,1)
!
    js=j2js(j_src)
    ja2=j2ja(j_src)
    jsv2=jsv4js(js)
    ict4h=1
!
    jj=0
    tmp=0.0d0
    do jsd=1,njsd(jsv2,ict4h)
      jsv1=jsv4jsd(jsd,jsv2)
      do ja1=1,nval(atm_element(jsv1))
        jj=jj+1
        tmp=tmp+dsij(ja1,ja2,jsd,jsv2)*vect_in(jj)
      enddo   
    enddo  
    result_value=tmp
!
    if (m_in /= jj) then
      write(*,*)'ERROR(inner_product_with_s_proj)'
      write(*,*)'m_in, jj=',m_in,jj
      stop
    endif
!
  end subroutine inner_product_with_s_proj
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert the density matrix and energy density matrix
!    into the arrays in the ELSES code 
!
  subroutine convert_dm_and_edm(atm_index,orb_index,dm_loc,edm_loc,jsv4jsk,jsk4jsv,jjkset)
!    
    use M_qm_domain,   only : dbij, dpij, jsv4jsd, noav, njsd, nval, atm_element
!           CHANGED   : dbij(:,:,orb_index,atom_index), dpij(:,:,orb_index,atom_index) 
!           Unchanged : jsv4jsd, njsd       
!
    implicit none
    integer, intent(in)      :: atm_index, orb_index
    real(8), intent(in)      ::  dm_loc(:) ! density matrix for a given basis
    real(8), intent(in)      :: edm_loc(:) ! energy density matrix for a given basis
    integer, intent(in)      :: jsv4jsk(:), jsk4jsv(:), jjkset(:) ! info. of the projection 
!
    integer :: m_int
    integer :: jsv2, ja2, ict4h, jj, jsd1, jsv1, jsk1, jjkset1, ja1, jjk1
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set the local variables
!
    jsv2=atm_index
    ja2=orb_index
    ict4h=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Get the number of basis in the interaction list
!     ----> m_int
!
    call get_interac_list_num_proj_atom(atm_index,m_int)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Trivial checking for the input quantities (non-essential)
!
    if (size(dm_loc,1) /= m_int) then
       write(*,*)'ERROR:convert_dm_and_edm:size(dm_loc,1)= ',size(dm_loc,1)
    endif
!
    if (size(edm_loc,1) /= m_int) then
       write(*,*)'ERROR:convert_dm_and_edm:size(edm_loc,1)= ',size(edm_loc,1)
    endif
!
    if ((jsv2 <=0 ).or. (jsv2 > noav)) then
       write(*,*)'ERROR:convert_dm_and_edm'
       write(*,*)' atm_index, orb_index=',jsv2, ja2
    endif   
!
    if ((ja2 <=0 ).or. (ja2 > nval(atm_element(jsv2)))) then
       write(*,*)'ERROR:convert_dm_and_edm'
       write(*,*)' atm_index, orb_index=',jsv2, ja2
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Convert the DM and EDM into dbij and dpij
!
!
!
    jj=0
    do jsd1=1,njsd(jsv2,ict4h)
      jsv1=jsv4jsd(jsd1,jsv2)
      jsk1=jsk4jsv(jsv1)
      if (jsk1 == 0) then
        write(*,*)'ERROR(convert_dm_and_edm):jsk1=',jsk1
        stop
      else
        jjkset1=jjkset(jsk1)
        do ja1=1,nval(atm_element(jsv1))
          jjk1=jjkset1+ja1
          jj=jj+1
          if (jj > m_int) then
            write(*,*)'ERROR(set_interac_list_proj):jj,m_int=',jj,m_int
            write(*,*)' jsd1, jsv1, jsk1, jjk1=',jsd1, jsv1, jsk1, jjk1
            stop
          endif
          dbij(ja1,ja2,jsd1,jsv2)=dm_loc(jj)
          dpij(ja1,ja2,jsd1,jsv2)=edm_loc(jj)
        enddo
       endif  
    enddo   
!
    if (m_int /= jj) then
       write(*,*)'ERROR(convert_dm_and_edm):jj,m_int=',jj,m_int
       stop
    endif   
!
  end subroutine convert_dm_and_edm
!
end module M_qm_projection
