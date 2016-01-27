!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_gkrylov_hist
!
  implicit none
!
  logical, allocatable :: proj_list_unchanged(:)
!     ----> proj_list_check(1:noav)
!            proj_list_check(jsv) = 1, 
!               if the proj list for the JSV-th atom is the same as in the previous MD step
!
  integer, allocatable :: proj_list_length_prv(:)
!     ----> proj_list_length_prv(1:noav)
!              The number of booked members for each atom.
!
  integer, allocatable :: proj_list_prv(:,:)
!     ----> booking list in the previous MD step, used as jsv1 = proj_list_prv(jsk,jsv2)
!
  integer, allocatable :: mat_size_str(:)
!     ---->  mat_size_str(1:noav) : matrix size for each projected subspace
!
  real(8), allocatable :: s_inv_e_j(:,:)
!     ---->   s_inv_e_j(matrix_size , n_tot_base) : S^{-1}_e_j for the j-th basis
!
  real(8), allocatable :: rho_e_j(:,:)
!     ---->   rho_e_j(matrix_size , n_tot_base) : rho e_j for the j-th basis
!
  real(8), allocatable ::  u_hst_str(:,:,:,:)
! real(8), allocatable :: su_hst_str(:,:,:,:)
! real(8), allocatable :: hu_hst_str(:,:,:,:)
  real(8), allocatable :: v_mat_kr_str(:,:,:,:)
  real(8), allocatable ::  u_b_hst_str(:,:,:)
!
  private
  public :: proj_list_unchanged
! public :: proj_list_length_prv
! public :: proj_list_prv
! public :: mat_size_str
  public :: s_inv_e_j
  public :: rho_e_j
  public :: u_hst_str, v_mat_kr_str
  public :: u_b_hst_str
!
  public :: check_proj_list
  public :: alloc_for_hist
  public :: set_s_inv_e_j_for_hist
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_proj_list
!
    use M_qm_domain,        only : i_verbose, noav                            !(unchanged)
    use elses_arr_kry_glo,  only : proj_list_length => noak_str, proj_list => jsv4jsk_str  !(unchanged)
!   use M_config,           only : config !(unchanged)
!
    implicit none
!
    logical :: flag_for_init
    integer :: ierr
    integer :: size1, size2
    integer :: size1b, size2b
    integer :: jsv, jsk
    integer :: n_count
!
    if (i_verbose >= 1) then
      write(*,*)'@@ check_proj_list'
    endif  
!
    if ( .not. allocated(proj_list_length_prv) ) then
      flag_for_init = .true. 
    else
      flag_for_init = .false.
    endif   
!
    if (i_verbose >= 1) then
      write(*,*)' flag_for_init=',flag_for_init
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial procedures 
!
    if (flag_for_init) then
!
      allocate (proj_list_length_prv(noav),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
!
      allocate (proj_list_unchanged(noav),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
      proj_list_unchanged(:)=.false.
!
      size1 = size(proj_list,1)
      size2 = size(proj_list,2)
      size1 = min(nint(size1*1.5),noav)
!        ---> size1 is set to be enough large
      allocate (proj_list_prv(size1,size2),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
!
      proj_list_length_prv(:)=proj_list_length(:)
      do jsv=1,noav
        do jsk=1,proj_list_length(jsv)
          proj_list_prv(jsk,jsv)=proj_list(jsk,jsv)
        enddo   
      enddo   
!
      return
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check if the projection list is the same as in the prev. step
!       ----> output : proj_list_unchanged
!
    do jsv=1,noav
!
      proj_list_unchanged(jsv) =.false.
      if ( proj_list_length_prv(jsv) == proj_list_length(jsv) ) then
        proj_list_unchanged(jsv) =.true.
        do jsk=1,proj_list_length(jsv)
          if ( proj_list(jsk,jsv) /= proj_list_prv(jsk,jsv)) then
            proj_list_unchanged(jsv) =.false.
          endif
        enddo   
      endif   
!
!     write(*,*)'check proj list=',jsv, proj_list_unchanged(jsv)
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Count up the number of atom for the same booking list (non-essential)
!       ----> n_count
!
    n_count=0
    do jsv=1,noav
      if ( proj_list_unchanged(jsv) ) n_count=n_count+1
    enddo   
    if (i_verbose >= 1) then
      write(*,*)'  number of atoms with unchnaged proj lsit =',n_count
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save the proj list for the next MD step
!       ----> output : proj_list_prv, proj_list_length_prv
!
    size1  = size(proj_list,1)
    size2  = size(proj_list,2)
    size1b = size(proj_list_prv,1)
    size2b = size(proj_list_prv,2)
    if (i_verbose >= 1) then
      write(*,*)'size1 , size2  =',size1, size2
      write(*,*)'size1b, size2b =',size1b,size2b
    endif  
!
    if (size1 > size1b) then
      if (i_verbose >= 1) then
        write(*,*)'  Reallocation of proj_list_prv'
      endif  
      deallocate (proj_list_prv,stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Dealloc. error (check_proj_list)'
        stop
      endif   
      allocate (proj_list_prv(size1,size2),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
    endif   
!
    proj_list_length_prv(:)=proj_list_length(:)
    do jsv=1,noav
      do jsk=1,proj_list_length(jsv)
        proj_list_prv(jsk,jsv)=proj_list(jsk,jsv)
      enddo   
    enddo   
!
  end subroutine check_proj_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine alloc_for_hist
!
    use M_config,           only : config !(unchanged)
    use M_qm_domain,        only : i_verbose, noav   !(unchanged)
    use M_qm_projection,    only : proj_get_mat_size !(routines)
    use elses_mod_orb2,     only : n_tot_base        !(unchanged)
    implicit none
    integer :: ierr
    integer :: atm_index, mat_size, num_atom_proj
    integer :: mat_size_max, num_atom_proj_max
    integer :: size_for_alloc
    logical :: flag_for_init
    logical :: flag_for_s_inv_e_j
    logical :: flag_for_rho_e_j
!
    if (i_verbose >= 1) then
      write(*,*)'@@ alloc_for_hist'
    endif  
!
    if ( .not. allocated(mat_size_str) ) then
      flag_for_init = .true. 
    else
      flag_for_init = .false.
    endif   
!
    if (flag_for_init) then
      allocate (mat_size_str(noav),stat=ierr)
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
    endif  
!
    do atm_index=1,noav
      call proj_get_mat_size(atm_index, mat_size, num_atom_proj) 
      mat_size_str(atm_index)=mat_size
    enddo
!   
    mat_size_max = maxval(mat_size_str)
    size_for_alloc=nint(1.5*mat_size_max)
!
    if (i_verbose >= 1) then
      write(*,*)' mat_size_max  = ',mat_size_max
      write(*,*)' size_for_alloc= ',size_for_alloc
    endif  
!
    if ( .not. allocated(s_inv_e_j) ) then
      flag_for_s_inv_e_j = .true.
    else
      flag_for_s_inv_e_j = .false.
    endif   
!
    if ( flag_for_s_inv_e_j ) then
      if (i_verbose >= 1) then
        write(*,'(a,f10.5)')'  alloc. of s_inv_e_j : memory size [GB] =', &
&            8.0d0*dble(size_for_alloc)*dble(n_tot_base)/1.0d9
      endif 
      allocate ( s_inv_e_j(size_for_alloc, n_tot_base), stat=ierr )
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
      s_inv_e_j(:,:)=0.0d0
    endif  
!
    if ( .not. allocated(rho_e_j) ) then
      flag_for_rho_e_j = .true.
    else
      flag_for_rho_e_j = .false.
    endif   
!
    if ( flag_for_rho_e_j ) then
      if (i_verbose >= 1) then
        write(*,'(a,f10.5)')'  alloc. of rho_e_j   : memory size [GB] =', &
&            8.0d0*dble(size_for_alloc)*dble(n_tot_base)/1.0d9
      endif 
      allocate ( rho_e_j(size_for_alloc, n_tot_base), stat=ierr )
      if (ierr /= 0) then
        write(*,*)'Alloc. error (check_proj_list)'
        stop
      endif   
      rho_e_j(:,:)=0.0d0
    endif  
!
    if (config%calc%solver%mode_for_large_memory >= 1) then
      call alloc_for_hist_str
    endif
!
  end subroutine alloc_for_hist
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine alloc_for_hist_str
!
    use M_config,           only : config !(unchanged)
    use M_qm_domain,        only : i_verbose, nval, noav,atm_element   !(unchanged)
    use M_qm_projection,    only : get_interac_list_num_proj !(routine)
    use elses_mod_orb2,     only : js2j !(unchanged)
!
    implicit none
    integer :: ierr
    integer :: nval_max
    integer :: kr_dim_max_input
!
    integer :: atm_index, orb_index,j
    integer :: mat_size_max, m_int
    integer :: size_for_alloc
!
    nval_max=maxval(nval)
    kr_dim_max_input=config%calc%solver%dimension
!
    if (allocated(u_hst_str)) then
      if (i_verbose >= 1) write(*,*)'@@ alloc_for_hist_str .. is skipped'
      return
    endif
!
    if (i_verbose >= 1) then
       write(*,*)'@@ alloc_for_hist_str'
       write(*,*)'  kr_dim_max_input =',kr_dim_max_input
       write(*,*)'  nval_max         =',nval_max
    endif
!   
    mat_size_max=0
    do atm_index=1,noav
      do orb_index=1, nval(atm_element(atm_index))
        j=js2j(orb_index,atm_index)
        call get_interac_list_num_proj(j,m_int)
        if (m_int > mat_size_max) mat_size_max=m_int
      enddo
    enddo
    size_for_alloc=nint(1.5*mat_size_max)
!
    if (i_verbose >= 1) write(*,*)'size_for_alloc=',size_for_alloc
!
    if (.not. allocated(u_hst_str)) then
      if  (i_verbose >= 1) then 
        write(*,'(a,f15.10)')' Hst_str size [GB] =',8.0d0*dble(size_for_alloc)   & 
&                                 *dble(kr_dim_max_input)*dble(nval_max)*dble(noav)/1.0d9
      endif  
      allocate ( u_hst_str(size_for_alloc, kr_dim_max_input, nval_max, noav), stat=ierr )
      if (ierr /= 0) then
        write(*,*)'Alloc. error (alloc_for_hist_str)'
        stop
      endif   
    endif  
!  
    if (.not. allocated(u_b_hst_str)) then
      if  (i_verbose >= 1) then 
        write(*,'(a,f15.10)')' Hst_str size [GB] =',8.0d0 &
&                                 *dble(kr_dim_max_input)*dble(nval_max)*dble(noav)/1.0d9
      endif  
      allocate ( u_b_hst_str(kr_dim_max_input, nval_max, noav), stat=ierr )
      if (ierr /= 0) then
        write(*,*)'Alloc. error (alloc_for_hist_str)'
        stop
      endif   
    endif  

!   if (.not. allocated(su_hst_str)) then
!     allocate ( su_hst_str(size_for_alloc, kr_dim_max_input, nval_max, noav), stat=ierr )
!     if (ierr /= 0) stop 'Alloc. error (alloc_for_hist_str)'
!   endif  
!
!   if (.not. allocated(hu_hst_str)) then
!     allocate ( hu_hst_str(size_for_alloc, kr_dim_max_input, nval_max, noav), stat=ierr )
!     if (ierr /= 0) stop 'Alloc. error (alloc_for_hist_str)'
!   endif  
!
    if (.not. allocated(v_mat_kr_str)) then
      allocate (v_mat_kr_str(kr_dim_max_input,kr_dim_max_input,nval_max,noav),stat=ierr)
      if (ierr /= 0) stop 'Alloc. error (alloc_for_hist_str)'
    endif  
!
  end subroutine alloc_for_hist_str
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_s_inv_e_j_for_hist(imode)
!
    use elses_mod_orb2,     only : n_tot_base        !(unchanged)
    use M_qm_domain,        only : i_verbose, nval, noav,atm_element   !(unchanged)
    use elses_mod_orb2,     only : js2j !(unchanged)
    implicit none
    integer, intent(in)  :: imode
    integer :: atm_index, orb_index
    integer :: i_show, j
    integer :: n_count
!
    i_show=4
    if (i_verbose >= 1) then
      write(*,*)'@@ set_s_inv_e_j_for_hist:imode=',imode
    endif  
!
    if (imode == 0) then
      write(*,*)'  imode=0:Reset s_inv_e_j'
      s_inv_e_j(:,:)=0.0d0
      return
    endif   
!
    write(*,*)'  imode=1:Set s_inv_e_j as historical one'
    n_count=0
    do atm_index=1,noav
      if ( .not. proj_list_unchanged(atm_index) ) then
        do orb_index=1,nval(atm_element(atm_index))
          j=js2j(orb_index,atm_index)
          s_inv_e_j(:,j)=0.0d0
          rho_e_j(:,j)=0.0d0
        enddo  
      else
        n_count=n_count+1
        if (atm_index <= i_show) then
          write(*,*)'  unchanged proj. list, atm_index=',atm_index
        endif  
      endif   
    enddo
!
    if (i_verbose >= 1) then
      write(*,*)' n_count=',n_count
    endif  
!
  end subroutine set_s_inv_e_j_for_hist
  
!
!
end module M_la_gkrylov_hist

