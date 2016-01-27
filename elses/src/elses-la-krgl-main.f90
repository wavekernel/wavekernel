!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_krgl_main
!
  implicit none
!
  real(8), allocatable :: wt_kr(:,:,:)
  real(8), allocatable :: eig_kr(:,:,:)
  integer, allocatable :: kr_dim_str(:,:)
!
  private
  public :: krgl_main
  public :: krgl_alloc
  public :: set_chemical_potential
  public :: plot_ldos_in_krgl
  public :: reset_dsij
  public :: sym_dbij_and_dpij
  public :: wt_kr, eig_kr, kr_dim_str
!
! The below routines are public only for elses-la-gkrylov-main.f90
  public :: calc_u_su_hu
  public :: get_num_of_nonzeros
  public :: calc_eig
  public :: check_matrix
  public :: schmidt_orth_w_s_mat
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine krgl_alloc(imode)
!
    use M_qm_domain, only : i_verbose, noav, nval !(unchanged)
    use M_config,  only : config !(unchaged)
    implicit none
    integer,          intent(in)  :: imode
    integer                       :: nval_max, ierr
    integer                       :: kr_dim_max_input
!
    kr_dim_max_input  = config%calc%solver%dimension
!
    if (i_verbose >= 1) then
      write(*,*)'@@ krgl_alloc:imode=',imode
      write(*,*)'   kr_dim_max_input=',kr_dim_max_input
    endif
!   
    nval_max=maxval(nval)
    if (i_verbose >= 1) then
      write(*,*)'   noav    =',noav
      write(*,*)'   nval_max=',nval_max
    endif
!
    if (imode == 1) then
!
      if (.not. allocated(kr_dim_str)) then
        allocate (kr_dim_str(nval_max,noav),stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in alloc (kr_dim_str)'
        kr_dim_str(:,:)=kr_dim_max_input
      endif   
!
      allocate (wt_kr(kr_dim_max_input,nval_max,noav),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (wt_kr)'
      wt_kr(:,:,:)=0.0d0
!
      allocate (eig_kr(kr_dim_max_input,nval_max,noav),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (eig_kr)'
      eig_kr(:,:,:)=0.0d0
!
    endif  
!
    if (imode == 2) then
      deallocate(wt_kr,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
      deallocate (eig_kr,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
    endif  
!
!
  end subroutine krgl_alloc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine krgl_main(j_src,b,jsv4jsk,jsk4jsv,jjkset,prc_index,scheme_mode)
!
!
!    INPUT 
!
    use M_qm_domain, only : i_verbose, chemical_potential,temp_for_electron, njsd !(unchanged)
    use M_config,  only : config !(unchaged)
    use M_lib_phys_const, only : ev4au !(unchaged)
    use elses_mod_orb2,   only : j2js,j2ja    ! (unchaged)
!
    use M_la_matvec_io,  only : get_interaction_list_number, & 
&       set_interaction_list, matvec_mul, inner_product_with_s          ! (routine)
    use M_la_matvec_routines, only : cg_s_mat_proj_r2, matvec_mul_proj_r ! (routine)
    use M_wall_clock_time,    only : get_system_clock_time     ! (routine)
    use M_qm_projection,      only : get_interac_list_num_proj, &
&                                    set_interac_list_proj, &
&                                    inner_product_with_s_proj, convert_dm_and_edm   ! (routine)
!
    use M_lib_math_func,  only : Fermi_Dirac_Func
!
    use M_la_gkrylov_hist, only : proj_list_unchanged, s_inv_e_j, rho_e_j !(CHANGED)
    use M_la_gkrylov_hist, only : u_hst_str, v_mat_kr_str !(CHANGED)
    use M_la_gkrylov_hist, only : u_b_hst_str             !(CHANGED)
!
    implicit none
    integer,          intent(in)  :: j_src
    integer,          intent(in)  :: jsv4jsk(:)
    integer,          intent(in)  :: jsk4jsv(:)
    integer,          intent(in)  :: jjkset(:)
    real(8),          intent(in)  :: b(:)
    integer,          intent(in)  :: prc_index
    character(len=32), intent(in) :: scheme_mode
!
    real(8), allocatable :: u(:)
    real(8), allocatable :: su(:)
    real(8), allocatable :: hu(:)
!
    real(8), allocatable ::  u_hst(:,:)
    real(8), allocatable :: su_hst(:,:)
    real(8), allocatable :: hu_hst(:,:)
!
    real(8), allocatable :: h_mat_kr(:,:)
    real(8), allocatable :: s_mat_kr(:,:)
    real(8), allocatable :: v_mat_kr(:,:)
    real(8), allocatable :: eig_wrk(:)
!
    real(8), allocatable ::  den_mat_kr(:,:)
    real(8), allocatable :: eden_mat_kr(:,:)
!
    real(8), allocatable ::  den_mat_kr_b(:)
    real(8), allocatable :: eden_mat_kr_b(:)
!
    real(8), allocatable ::   dm_loc(:)
    real(8), allocatable ::  edm_loc(:)
!
    integer :: ierr
!
    integer :: mat_dim
    integer :: m_int
    integer, allocatable     :: interaction_list(:)
!
    real(8) :: eps_c
    integer :: max_ite_for_cg_loop
!
    integer :: atm_index, orb_index, ini_elem
    integer :: i_show
!
    real(8) :: eps
    integer :: ite_cg
    integer :: i_check_mode
    integer :: imode_check_mat
!
    integer :: kr_dim, kr_dim_max, kr_dim_max_input, kr_dim_max_loop
    integer :: j, al
    integer :: m,n,m_pick
    real(8) :: ddd, ddd1, ddd2
    real(8) :: norm_factor, norm_factor1
!
    real(8) :: xmu, xbeta, f_occ
    integer :: jjk
!
    logical :: flag_for_s_inv
!
    integer :: i_kr_hst_str
    integer :: n_count, n_count2
    real(8) :: ddmax
!
    integer :: mArnoldi_q_wrk  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Initial procedure
!
    if (i_verbose >=1) then 
      i_check_mode=1
    else
      i_check_mode=0
    endif   
!
    kr_dim_max_input=config%calc%solver%dimension
    i_kr_hst_str=config%calc%solver%mode_for_large_memory
!
    mArnoldi_q_wrk      = config%calc%solver%mArnoldi_q
    eps_c               = config%calc%solver%inner_cg_loop%convergence_eps
    max_ite_for_cg_loop = config%calc%solver%inner_cg_loop%max_iteration
!
    if (mArnoldi_q_wrk > kr_dim_max_input-1) then
      write(*,*)'ERROR: mArnoldi_q =',mArnoldi_q_wrk
      stop
    endif
!
    atm_index=j2js(j_src)
    orb_index=j2ja(j_src)
    ini_elem=orb_index
!
    mat_dim=size(b,1)
!     ----> Matrix dimension
!
    i_show=0
    if (i_check_mode == 1) then
      if (atm_index <= 2) then
        i_show=1
      endif
    endif  
!
    if (i_show >= 1) then
      write(*,*)'atm_index=',atm_index
      write(*,*)'orb_index=',orb_index
      write(*,*)'the position of non-zero element in the initial vector=',ini_elem
      write(*,*)'mode for large memory=',i_kr_hst_str
      write(*,*)'mArnoldi_q           =',mArnoldi_q_wrk
      write(*,*)'inner_cg_loop%convergence_eps =', eps_c
      write(*,*)'inner_cg_loop%max_iteration   =', max_ite_for_cg_loop
    endif  
!
    flag_for_s_inv=.false.
    if ( allocated(s_inv_e_j) ) then
       flag_for_s_inv=.true.
    endif   
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Tribial checking (check mode only)
!
    if (i_check_mode == 1) then
      ddd=0.0d0
      do j=1,mat_dim
        if (j == ini_elem) then 
           ddd=ddd+dabs(b(j)-1.0d0)
        else
           ddd=ddd+dabs(b(j))
        endif   
      enddo
!     write(*,*)'ddd=',ddd
      if (dabs(ddd) > 1.0d-10) then
        write(*,*)'ERROR:ini_elem=',ini_elem
        stop
      endif   
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Construct the interaction list
!
    call get_interac_list_num_proj(j_src,m_int)
    if (i_show >= 1) then
      write(*,*)'m_int=',m_int
    endif  
!
    allocate (interaction_list(m_int),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc1'
    call set_interac_list_proj(j_src,interaction_list,jsk4jsv,jjkset)
!    
    if (i_show >= 1) then
      write(*,*)' ...set_interac_list_proj is ended'
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Allocation of working vectors
!
!
    allocate (u(mat_dim),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (u)'
    u(:)=0.0d0
!
    allocate (su(mat_dim),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (su)'
    su(:)=0.0d0
!
    allocate (hu(mat_dim),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (hu)'
    hu(:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Allocation of arraies with kr_dim_max_input
!
    allocate (u_hst(mat_dim,kr_dim_max_input),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (u_hst)'
    u_hst(:,:)=0.0d0
!
    allocate (su_hst(mat_dim,kr_dim_max_input),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (su_hst)'
    su_hst(:,:)=0.0d0
!
    allocate (hu_hst(mat_dim,kr_dim_max_input),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (hu_hst)'
    hu_hst(:,:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Flags for several techniques
!
    if (allocated(s_inv_e_j)) then
      flag_for_s_inv = .true.
    else
      flag_for_s_inv = .false.
    endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the initial basis : |l_1> = S^{-1} |b> ( non-normalized vector )
!
    ite_cg=max_ite_for_cg_loop
    eps=eps_c
    u(:)=b(:)
!
    if (i_show >= 1) then
      write(*,'(a,i10,F10.5,2I10)') 'ite_cg, eps, atm_index, orb_index=', & 
&                                   ite_cg, eps, atm_index, orb_index
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Krylov loop for  u_m, Su_m, Hu_m, ( m = 1,...., kr_dim_max )
!
    kr_dim_max_loop=kr_dim_max_input
!
    if ( (i_kr_hst_str == 1)  .and. (prc_index == 2) ) then
      kr_dim_max_loop=0
      kr_dim_max=kr_dim_str(orb_index,atm_index)
      if ( (kr_dim_max <= 0) .or. (kr_dim_max > kr_dim_max_input)) then
        write(*,*) 'ERROR(KR loop):kr_dim_max=',kr_dim_max
        stop
      endif   
       u_hst(:,:) =  0.0d0
      su_hst(:,:) =  0.0d0
      hu_hst(:,:) =  0.0d0
      do kr_dim=1,kr_dim_max
        do j=1,m_int
          jjk=interaction_list(j)
           u_hst(jjk,kr_dim) =  u_hst_str(j, kr_dim, orb_index, atm_index)
        enddo  ! NOTE: NOT ALL the elements are generated.
      enddo
    endif   
!
    do kr_dim=1,kr_dim_max_loop
!       u(:) : input |l_n> ( non-normalized vector) 
!
      kr_dim_max=kr_dim
      su(:)=0.0d0
      hu(:)=0.0d0
      norm_factor=0.0d0
      call calc_u_su_hu(u,su,hu,norm_factor,jsv4jsk,jsk4jsv,jjkset,ierr)
!        --->  u(:) :   |u_n> = |l_n> / |<l_n|S|l_n>|^{1/2}
!                     (S-normalized vector) 
!             su(:) :  S|u_n>
!             hu(:) :  H|u_n>
!            norm_factor : C = |<l_n|S|l_n>|^{1/2}
!
      if (ierr /= 0) then
        if (i_verbose >= 1) then 
          write(*,'(a,3i20)')'Info:terminated:',j_src,atm_index,kr_dim
          write(*,*)'  j_src    = ', j_src
          write(*,*)' mat_dim   = ', mat_dim
          write(*,*)' kr_dim    = ', kr_dim
          write(*,*)' prc_index = ', prc_index
          write(*,*)' atm_index = ', atm_index
          write(*,*)' orb_index = ', orb_index
          write(*,*)'  njsd     = ', njsd(atm_index,1)
        endif
        kr_dim_max=kr_dim-1
        exit
      endif   
!
!
      if (kr_dim == 1) norm_factor1=norm_factor
!
       u_hst(:,kr_dim) =  u(:)
      su_hst(:,kr_dim) = su(:)
      hu_hst(:,kr_dim) = hu(:)
!
      if ( (prc_index == 1) .and. (i_kr_hst_str == 1) ) then
        u_b_hst_str(kr_dim, orb_index, atm_index) =  u_hst(ini_elem,kr_dim)
        do j=1,m_int
          jjk=interaction_list(j)
           u_hst_str(j, kr_dim, orb_index, atm_index) =  u_hst(jjk,kr_dim)
        enddo  
      endif   
!
      if (kr_dim == kr_dim_max_input) cycle
!
      u(:)=hu(:)
      if (scheme_mode == 'ekrgl') then
        if (kr_dim == kr_dim_max_input - mArnoldi_q_wrk) then
          if (i_show >= 1) then
            write(*,*) ' Multiple KR space : kr_dim=',kr_dim
          endif  
          ite_cg=max_ite_for_cg_loop
          eps=eps_c
!
          if (flag_for_s_inv) then
            if (i_show >= 1) then
              write(*,*) ' S_inv_e_j is allocated'
              write(*,*) ' proj_list_unchanged=',atm_index, proj_list_unchanged(atm_index)
            endif   
            if (prc_index == 1) then
              u(1:mat_dim)=s_inv_e_j(1:mat_dim,j_src)
              call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
              s_inv_e_j(1:mat_dim,j_src)=u(1:mat_dim)
              if (i_show >= 1) then
                write(*,*)'Sinv(1): j_src, eps, ite_cg=',j_src,eps,ite_cg
              endif   
            else
              u(1:mat_dim)=s_inv_e_j(1:mat_dim,j_src)
              call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
              s_inv_e_j(1:mat_dim,j_src)=u(1:mat_dim)
              if (i_show >= 1) then
                write(*,*)'Sinv(2): j_src, eps, ite_cg=',j_src,eps,ite_cg
              endif  
            endif 
          else
            u(:)=0.0d0
            call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
          endif   
!
        endif
      endif  
!
      if (kr_dim /= kr_dim_max_input) then
         call schmidt_orth_w_s_mat(u,u_hst,su_hst,kr_dim)
!        --->  u(:) :   |l_{n+1}> 
      endif   
!
    enddo
    kr_dim_str(orb_index,atm_index)=kr_dim_max
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Shrink the reduced subspace, if S is deviated from I, because of insufficient bases
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
!
      call get_num_of_nonzeros(u_hst(:,kr_dim_max),n_count,n_count2)
      if (n_count < kr_dim_max*2) then
        do n=2,kr_dim_max
          ddmax=0.0d0
          do m=1,n
            ddd = dot_product(u_hst(:,m),su_hst(:,n)) ! ( u_m S u_n )
            if (m == n ) ddd = ddd-1.0d0
            ddd = dabs(ddd)
            if (ddmax < ddd) then 
              m_pick=m
              ddmax=ddd
            endif  
          enddo
          if (ddmax > 1.0d-10) then
            if (i_verbose >= 1) then
              call get_num_of_nonzeros(u_hst(:,n),n_count,n_count2)
              write(*,*)'INFO:Shrink kr_dim_max       = ', n-1, kr_dim_max
              write(*,*)'INFO:   atm_index,orb_index  = ', atm_index, orb_index
              write(*,*)'INFO:   n, m_pick            = ', n, m_pick
              write(*,*)'INFO:   ddmax                = ', ddmax
              write(*,*)'INFO:   n_count, n_count2    = ', n_count, n_count2
              write(*,*)'INFO:   mat_dim              = ', mat_dim
              write(*,*)'INFO:   njsd                 = ', njsd(atm_index,1)
            endif
            kr_dim_max=n-1
            kr_dim_str(orb_index,atm_index)=kr_dim_max
            exit   
          endif
        enddo
      endif  
!
    endif
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Construct the reduced matrices with kr_dim_max
!
    allocate (v_mat_kr(kr_dim_max,kr_dim_max),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (v_mat_kr)'
    v_mat_kr(:,:)=0.0d0
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
!
      allocate (h_mat_kr(kr_dim_max,kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (h_mat_kr)'
      h_mat_kr(:,:)=0.0d0
!
      allocate (s_mat_kr(kr_dim_max,kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (s_mat_kr)'
      s_mat_kr(:,:)=0.0d0
!
!
      allocate (eig_wrk(kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (eig_wrk)'
      eig_wrk(:)=0.0d0
!
      do n=1,kr_dim_max
        do m=1,kr_dim_max
          h_mat_kr(m,n)=dot_product(u_hst(:,m),hu_hst(:,n)) ! ( u_m H u_n )
          s_mat_kr(m,n)=dot_product(u_hst(:,m),su_hst(:,n)) ! ( u_m S u_n )
        enddo
      enddo
!
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Trivial check the matrices (check mode only)
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
!
    if ( i_check_mode ==1 ) then
      if (i_show >= 1) then
        do n=1,kr_dim_max
          call get_num_of_nonzeros(u_hst(:,n),n_count,n_count2)
          write(*,'(a,3i10)')'n, num of non-zero =',n, n_count, n_count2
        enddo  
      endif  
!
      imode_check_mat=1
      call check_matrix(s_mat_kr,ddd,imode_check_mat)
      if (i_show >= 1) then
         write(*,*)' Check (reduced S) = I : deviation=',ddd
      endif  
      if (ddd > 1.0d-9) then
        write(*,*)' Abort:Deviation  for (reduced S)=I:deviation=',ddd
        write(*,*)' j_src=',j_src
        write(*,*)' mat_dim   = ', mat_dim
        write(*,*)' kr_dim_max= ', kr_dim_max
        write(*,*)' prc_index = ', prc_index
        write(*,*)' atm_index = ', atm_index
        write(*,*)' orb_index = ', orb_index
        write(*,*)'  njsd     = ', njsd(atm_index,1)
        do n=1,kr_dim_max
          call get_num_of_nonzeros(u_hst(:,n),n_count,n_count2)
          write(*,'(a,3i10)')'DUMP:n, num of non-zero =',n, n_count, n_count2
          do m=1,kr_dim_max
            write(*,'(a,2i10,2f20.10)')'DUMP:m, n, H, S =',m, n, h_mat_kr(m,n), s_mat_kr(m,n)
          enddo
        enddo
        write(*,*)' Abort:Deviation  for (reduced S)=I:deviation=',ddd
        stop
      endif   
!
      imode_check_mat=2
      call check_matrix(h_mat_kr,ddd,imode_check_mat)
      if (i_show >= 1) then
        write(*,*)' Check symmetric matrix for (reduced H) : result=',ddd
      endif  
      if (ddd > 1.0d-10) then
        write(*,*)' Abort:Assymmetric for reduced S'
        write(*,*)' j_src=',j_src
        write(*,*)' mat_dim   = ', mat_dim
        write(*,*)' kr_dim_max= ', kr_dim_max
        write(*,*)' prc_index = ', prc_index
        write(*,*)' atm_index = ', atm_index
        write(*,*)' orb_index = ', orb_index
        write(*,*)'  njsd     = ', njsd(atm_index,1)
        stop
      endif   
!
      call check_matrix(s_mat_kr,ddd,imode_check_mat)
      if (i_show >= 1) then
        write(*,*)' Check symmetric matrix for (reduced S) : result=',ddd
      endif  
      if (ddd > 1.0d-10) then
        write(*,*)' Abort:Assymmetric for reduced S'
        write(*,*)' j_src=',j_src
        write(*,*)' mat_dim   = ', mat_dim
        write(*,*)' kr_dim_max= ', kr_dim_max
        write(*,*)' prc_index = ', prc_index
        write(*,*)' atm_index = ', atm_index
        write(*,*)' orb_index = ', orb_index
        write(*,*)'  njsd     = ', njsd(atm_index,1)
        stop
      endif   

!
    endif  
    endif  
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Monitor the matrix (check mode only)
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
    if ( i_check_mode == 1 ) then
!
      ddd=0.0d0
      do n=1,kr_dim_max
        ddd=ddd+h_mat_kr(n,n)
        if (i_show >= 1) then
          write(*,'(a,i10,f30.20)') 'H(n,n)=', n, h_mat_kr(n,n)*ev4au
        endif  
        if (dabs(h_mat_kr(n,n))*ev4au > 1000.0d0) then
          write(*,'(a,3i10,f30.20)') 'Warning:atm,orb,n,H(n,n)=', atm_index, orb_index, n, h_mat_kr(n,n)*ev4au
        endif
      enddo   
!
      if (i_show >= 1) then
        write(*,'(a,2i10,f30.20)')'atm, orb, Tr[H_kr] (eV) =',atm_index, orb_index, ddd*ev4au
      endif  
!
    endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Eigen states for the reduced matrix
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
      call calc_eig(v_mat_kr,eig_wrk,h_mat_kr,s_mat_kr)
!
      eig_kr(1:kr_dim_max,orb_index,atm_index)=eig_wrk(1:kr_dim_max)
!
      if ( i_kr_hst_str == 1 ) then
        if (.not. allocated(v_mat_kr_str)) stop 'Alloc. error (v_mat_kr_str)'
        do kr_dim=1,kr_dim_max
          v_mat_kr_str(1:kr_dim_max,kr_dim,orb_index,atm_index)=v_mat_kr(1:kr_dim_max,kr_dim)
        enddo  
      endif   
!
    else
!
      if (.not. allocated(v_mat_kr_str)) stop 'Alloc. error (v_mat_kr_str)'
      do kr_dim=1,kr_dim_max
        v_mat_kr(1:kr_dim_max,kr_dim)=v_mat_kr_str(1:kr_dim_max,kr_dim,orb_index,atm_index)
      enddo  
!
    endif
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot eigen levels (check mode only)
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
    if ( i_check_mode == 1 ) then
!
      ddd=0.0d0
      do n=1,kr_dim_max
        ddd=ddd+eig_wrk(n)*ev4au
!       if (i_show >= 1) then
!         write(*,'(a,i10,f30.20)') 'Eig [eV] =', n, eig_wrk(n)*ev4au
!       endif  
      enddo   
!
      if (i_show >= 1) then
        write(*,*)' Sum of eigen values (eV) =',ddd
      endif  
!
    endif  
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the weight with assumption of |j> = (norm_factor1)^{-1} S |u_1>
!              (check mode only)
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
    if ( i_check_mode == 1 ) then
      do al=1,kr_dim_max
        ddd=dot_product(v_mat_kr(1:kr_dim_max,al),su_hst(ini_elem,1:kr_dim_max))
        wt_kr(al,orb_index,atm_index)=ddd*v_mat_kr(1,al)/norm_factor1
      enddo
    endif  
    endif  
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot the weight (check mode only)
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
    if ( i_check_mode == 1 ) then
      if (i_show >= 1) then
        ddd=0.0d0
        do al=1,kr_dim_max
          ddd=ddd+wt_kr(al,orb_index,atm_index)
!         write(*,*)'wt1=',al,wt_kr(al,orb_index,atm_index)
        enddo   
        write(*,'(a,2i10,f30.20)')' sum of wt1 =',orb_index,atm_index,ddd
      endif
    endif  
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the weight 
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
    do al=1,kr_dim_max
!     ddd=0.0d0
!     do n=1,kr_dim_max
!       ddd=ddd+v_mat_kr(n,al)*su_hst(ini_elem,n)
!     enddo   
      ddd=dot_product(v_mat_kr(1:kr_dim_max,al),su_hst(ini_elem,1:kr_dim_max))
!     ddd1=0.0d0
!     do n=1,kr_dim_max
!       ddd1=ddd1+v_mat_kr(n,al)*u_hst(ini_elem,n)
!     enddo   
      ddd1=dot_product(v_mat_kr(1:kr_dim_max,al),u_hst(ini_elem,1:kr_dim_max))
      wt_kr(al,orb_index,atm_index)=ddd*ddd1
    enddo
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot the weight (check mode only)
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
    if (i_show >= 1) then
      ddd=0.0d0
      do al=1,kr_dim_max
        ddd=ddd+wt_kr(al,orb_index,atm_index)
        write(*,'(a,i10,2f30.20)') 'Eig [eV] =', al, eig_wrk(al)*ev4au,wt_kr(al,orb_index,atm_index)
      enddo   
      write(*,'(a,2i10,f30.20)')' sum of wt2 =',orb_index,atm_index,ddd
    endif
    endif  
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the density matrix when prc_index=2
!
    if (prc_index == 2) then
!
      allocate ( den_mat_kr(kr_dim_max,kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (den_mat_kr)'
      den_mat_kr(:,:)=0.0d0
!
      allocate (eden_mat_kr(kr_dim_max,kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (eden_mat_kr)'
      eden_mat_kr(:,:)=0.0d0
!
      allocate ( den_mat_kr_b(kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (den_mat_kr_b)'
      den_mat_kr_b(:)=0.0d0
!
      allocate ( eden_mat_kr_b(kr_dim_max),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (eden_mat_kr_b)'
      eden_mat_kr_b(:)=0.0d0
!
      allocate ( dm_loc(m_int),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (dm_loc)'
      dm_loc(:)=0.0d0
!
      allocate ( edm_loc(m_int),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc (edm_loc)'
      edm_loc(:)=0.0d0
!
      xmu=chemical_potential
      xbeta=1.0d0/temp_for_electron
!
      if (i_show >= 1) then
        write(*,*)'Generate DM, prc_index=',prc_index
        write(*,*)'chem_pot         =',chemical_potential
        write(*,*)'temp_for_electron=',temp_for_electron
      endif  
!
      do al=1,kr_dim_max
        ddd=(eig_kr(al,orb_index,atm_index)-xmu)*xbeta
        f_occ=Fermi_Dirac_Func(ddd)
!       if (i_show >= 1) then
!         write(*,*)'f_occ=',al,f_occ
!       endif  
        do m=1,kr_dim_max
          ddd=eig_kr(al,orb_index,atm_index)
           den_mat_kr(:,m) =  den_mat_kr(:,m) + v_mat_kr(:,al) * v_mat_kr(m,al) * f_occ
          eden_mat_kr(:,m) = eden_mat_kr(:,m) + v_mat_kr(:,al) * v_mat_kr(m,al) * f_occ * ddd
        enddo   
      enddo   
!
      do n=1,kr_dim_max
       den_mat_kr_b(n) =dot_product( den_mat_kr(n, 1:kr_dim_max), u_hst(ini_elem, 1:kr_dim_max))
       eden_mat_kr_b(n)=dot_product(eden_mat_kr(n, 1:kr_dim_max), u_hst(ini_elem, 1:kr_dim_max))
      enddo
!   
      do j=1,m_int
        jjk=interaction_list(j)
         dm_loc(j) = dot_product( u_hst(jjk,:),  den_mat_kr_b(:) )
        edm_loc(j) = dot_product( u_hst(jjk,:), eden_mat_kr_b(:) )
      enddo
!   
      call convert_dm_and_edm(atm_index,orb_index,dm_loc,edm_loc,jsv4jsk,jsk4jsv,jjkset)
!
      if ( allocated(rho_e_j) ) then
        if (i_show == 1) then
          write(*,*)'  rho_e_j is stored:j_src=', j_src, atm_index, orb_index
        endif   
        if (i_check_mode == 1) then
          if (size(rho_e_j,1) < mat_dim) then
            write(*,*)'ERROR: krgl_main (rho_e_j):size1=',size(rho_e_j,1)
            stop
          endif   
          if (size(rho_e_j,2) < j_src) then
            write(*,*)'ERROR: krgl_main (rho_e_j):size2=',size(rho_e_j,2)
            stop
          endif   
        endif  
        do al=1,kr_dim_max
          rho_e_j(1:mat_dim,j_src)= rho_e_j(1:mat_dim,j_src) + u_hst(1:mat_dim, al)*den_mat_kr_b(al)
        enddo
      endif
!
      deallocate ( den_mat_kr,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      deallocate (eden_mat_kr,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      deallocate ( den_mat_kr_b,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      deallocate ( eden_mat_kr_b,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      deallocate ( dm_loc,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
      deallocate ( edm_loc,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Deallocation
!
    deallocate (interaction_list,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (interaction_list)'
!
    deallocate (u, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (u)'
!
    deallocate (su, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (su)'
!
    deallocate (hu, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (hu)'
!
    if (allocated(h_mat_kr)) then
      deallocate (h_mat_kr, stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc (h_mat_kr)'
    endif  
!
    if (allocated(s_mat_kr)) then
      deallocate (s_mat_kr, stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc (s_mat_kr)'
    endif  
!
    deallocate (v_mat_kr ,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (v_mat_kr)'
!
    if (allocated(eig_wrk)) then
      deallocate (eig_wrk ,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc (eig_wrk)'
    endif  
!
!
  end subroutine krgl_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_u_su_hu(u,su,hu,norm_factor,jsv4jsk,jsk4jsv,jjkset,ierr)
!
    use M_la_matvec_routines, only : matvec_mul_proj_r ! (routine)
!
    implicit none
    integer,          intent(in)  :: jsv4jsk(:)
    integer,          intent(in)  :: jsk4jsv(:)
    integer,          intent(in)  :: jjkset(:)
    real(8),       intent(inout)  :: u(:)
    real(8),         intent(out)  :: su(:)
    real(8),         intent(out)  :: hu(:)
    real(8),         intent(out)  :: norm_factor
    integer,         intent(out)  :: ierr
!
    real(8) :: ddd
!
      call matvec_mul_proj_r(u,su,'S',jsv4jsk,jsk4jsv,jjkset)
!               ---> su : S |m_n> (non-normalized vector)
!
      call matvec_mul_proj_r(u,hu,'H',jsv4jsk,jsk4jsv,jjkset)
!               ---> hu : H |m_n> (non-normalized vector)
!
      ddd=dot_product(u(:),su(:))
!      
      if (ddd < 0.0d0) then
        write(*,*)'ERROR:(calc_u_su_hu): <u|S|u>=',ddd
        stop
      endif
!
      ddd=dsqrt(ddd)
!
      ierr=0
      if (dabs(ddd) < 1.0d-10) then
!        write(*,*)'ERROR:calc_u_su_hu'
         write(*,*)'Info(calc_u_su_hu):Too small normalization factor:=',ddd
         ierr=1
         return
      endif   
!
      if (dabs(ddd) > 1.0d4) then
!        write(*,*)'ERROR:calc_u_su_hu'
         write(*,*)'Warning(calc_u_su_hu):Too large norm ? : N = <u|S|u> =', ddd
      endif   
!
      norm_factor=1.0d0/ddd
!        ---> normalization factor
!
      u(:)   = u(:)/ddd
!        ---> u :   |u_n> (S-normalized vector)
!
      su(:) = su(:)/ddd
!        ---> su : S|u_n> 
      hu(:) = hu(:)/ddd
!        ---> hu : H|u_n> 

  end subroutine calc_u_su_hu
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!@@ Gram-Schmidt orthogonalization with S matrix
!   input v(:)
!   output : |v_new> = |v> - \sum_m^{n} q_{m} |u_m>
!                where  q_{m} = <u_m|S|v>
!
  subroutine schmidt_orth_w_s_mat(v,u_hst,su_hst,kr_dim)
!
!
    implicit none
    real(8),          intent(inout)  :: v(:)
    real(8),          intent(in)  :: u_hst(:,:)
    real(8),          intent(in)  :: su_hst(:,:)
    integer,          intent(in)  :: kr_dim
!
    integer :: m
    real(8) :: ddd
!
    do m=kr_dim,1,-1
      ddd=dot_product(su_hst(:,m),v(:))
!        -----> q_{m} = <u_m|S|v>
      v(:)=v(:)-ddd*u_hst(:,m)
!        -----> |v_new> = |v> - \sum_m^{n} q_{m} |u_m>
    enddo
!   
  end subroutine schmidt_orth_w_s_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_eig(v_mat_kr,eig_kr,h_mat_kr,s_mat_kr)
!
    implicit none
    real(8),         intent(out)  :: v_mat_kr(:,:)
    real(8),         intent(out)  :: eig_kr(:)
    real(8),         intent(in)   :: h_mat_kr(:,:)
    real(8),         intent(in)   :: s_mat_kr(:,:)
!
    integer :: n,ierr
    real(8), allocatable ::  a(:,:)
    real(8), allocatable ::  b(:,:)
    real(8), allocatable ::  eig(:)
    real(8), allocatable ::  wrk_vec(:)
!
    integer :: info, lda, lwork
    real(8), allocatable ::  work(:)
    integer :: i_check_eigen, i,j
    real(8) ::ddd
!
!   i_check_eigen=1
    i_check_eigen=0
!
    n=size(h_mat_kr,1)
!
    ierr=0
    if (size(v_mat_kr,1) /= n) ierr=1 
    if (size(v_mat_kr,2) /= n) ierr=1 
    if (size(eig_kr,  1) /= n) ierr=1 
    if (size(h_mat_kr,1) /= n) ierr=1 
    if (size(h_mat_kr,2) /= n) ierr=1 
    if (size(s_mat_kr,1) /= n) ierr=1
    if (size(s_mat_kr,2) /= n) ierr=1
!
    if (ierr /=0) then
      write(*,*) 'ERROR(calc_eig):n=', n
      write(*,*) 'ERROR(calc_eig) size(v_mat_kr)=', size(v_mat_kr,1), size(v_mat_kr,2)
      write(*,*) 'ERROR(calc_eig) size(eig_kr)  =', size(eig_kr,1)
      write(*,*) 'ERROR(calc_eig) size(h_mat_kr)=', size(h_mat_kr,1), size(h_mat_kr,2)
      write(*,*) 'ERROR(calc_eig) size(s_mat_kr)=', size(s_mat_kr,1), size(s_mat_kr,2)
      stop
    endif   
!
    lda=n
    lwork=3*n-1
    info=0
!
    allocate (a(n,n),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc1'
    allocate (eig(n),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc1'
    allocate (work(3*n-1),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc1'
!
    a(:,:)=h_mat_kr
!
    call dsyev("V","U",n,a,lda,eig,work,lwork,info)
    if (info /= 0) then
      write(*,*)'LAPACK ERROR:info=',info
      stop
    endif   
!
    v_mat_kr(:,:)=a(:,:)
    eig_kr(:)=eig(:)
!
    if (i_check_eigen == 1) then
      allocate (wrk_vec(n),stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in alloc'
      do j=1,n
        wrk_vec(:)=0.0d0
        wrk_vec(:)=matmul(h_mat_kr(:,:), v_mat_kr(:,j))
        wrk_vec(:)=wrk_vec(:)-eig_kr(j)*v_mat_kr(:,j)
        ddd=dot_product(wrk_vec(:),wrk_vec(:))
        if (dabs(ddd) >= 1.0d-10) then
          write(*,*)'ERROR(calc_eig): j, residual =',j, ddd
          stop
        endif
!       write(*,*)' residual =',j, ddd
      enddo   
      deallocate (wrk_vec,stat=ierr)
      if (ierr /= 0) stop 'Abort:ERROR in dealloc'
    endif
!
    if (i_check_eigen == 1) then
      do j=1,n
       do i=1,n
         ddd=dot_product(v_mat_kr(i,:),v_mat_kr(j,:))
         if (i == j) ddd=ddd-1.0d0
         if (dabs(ddd) >= 1.0d-10) then
           write(*,*)'ERROR(calc_eig): j, residual in completeness =',i,j, ddd
           stop
         endif
       enddo
      enddo   
    endif

    deallocate (a,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc(a)'
!
    deallocate (eig,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc(eig)'
!
    deallocate (work,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc(work)'

  end subroutine calc_eig
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_chemical_potential(rNelec)
!        -----> chemical_potential in (M_qm_domain)
!
    use M_config,  only : config !(unchaged)
    use M_lib_phys_const, only : ev4au !(unchanged)
    use M_qm_domain,      only : i_verbose, noav, nval, temp_for_electron, &
&                                atm_element, chemical_potential 
                     ! CHANGED : chemical potential, others: unchanged
    implicit none
    real(8),          intent(in)  :: rNelec
!
    integer :: kr_dim_max_input
    integer :: ipe, npe
!
    integer :: iloop, nloopmax, irec
    integer :: jsv, js, nss, nval2, ja, nrecl
    real(8)  xtemp, xbeta, ddd1, ddd2, ddemin, ddemax
    real(8)  xmu0, xmumn, xmumx, xidos0, xidosE, tdossum, tensum
    real(8)  rRd, rEd, yyy, xexp, xqinit, err
    integer :: n_count, n_count1
!
!
    write(*,*)'@@ set_chemical_potential'
    nloopmax=100
    xtemp=temp_for_electron
!        ----> Temperature in Fermi Distribution
    xbeta=1.0d0/xtemp
!
    kr_dim_max_input=config%calc%solver%dimension
!
    write(*,*)'xtemp [au,eV]=',xtemp,xtemp*ev4au
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Trivial checking
!
    if (.not. allocated(kr_dim_str)) then
      write(*,*)'ERROR(set_chemical_potential):not allocated'
      stop
    endif   
!
    if ( maxval(kr_dim_str) <=0 ) then
      write(*,*)'ERROR(set_chemical_potential)'
      write(*,*)' maxval(kr_dim_str)=',maxval(kr_dim_str)
      stop
    endif   
!
    if ( maxval(kr_dim_str) > config%calc%solver%dimension  ) then
      write(*,*)'ERROR(set_chemical_potential)'
      write(*,*)' maxval(kr_dim_str)=',maxval(kr_dim_str)
      stop
    endif   
!
    if (rNelec .lt. 0.99d0) then
      write(*,*)'@@ ERROR!:SET_KRY_CHEM_POT:Nelec=',rNelec
      stop
    endif   
!
    if (i_verbose >= 1) then
      n_count=0
      n_count1=0
      do jsv=1,noav
        do ja=1,nval(atm_element(jsv))
          n_count1=n_count1+1
          if (kr_dim_str(ja,jsv) == kr_dim_max_input) n_count=n_count+1
        enddo
      enddo   
      write(*,*)'INFO:# bases with terminated KR subspace=',n_count1-n_count
    endif  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Find the min and max of energy levels
!
      ddemin=minval(eig_kr)
      ddemax=maxval(eig_kr)
!
      write(6,*)' ddemin [au,eV]=',ddemin,ddemin*ev4au
      write(6,*)' ddemax [au,eV]=',ddemax,ddemax*ev4au
      if (dabs(ddemin-ddemax) .lt. 1.0d-10) then
        write(6,*)'ERROR!:SET_KRY_CHEM_POT'
        stop
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ The bisection method
!
      xmu0=ddemin
      xmumn=ddemin
      xmumx=ddemax
!
      do iloop=1,nloopmax
!
       if (iloop == nloopmax) then
         write(*,*)'Bisection was not converged ! '
         stop
       endif
!   
       xidos0= 0.0d0
       xidosE= 0.0d0
       tdossum=0.0d0
       tensum =0.0d0
!
!      OMP-NOTE : Shared scalors : xmu0,xbeta
!$omp  parallel &
!$omp& default(shared) & 
!$omp& private(ipe,npe) &
!$omp& private(jsv,js,nss,nval2,ja,irec,nrecl) &
!$omp& private(rRd,rEd,yyy,xexp) &
!$omp& reduction(+ : xidos0) &
!$omp& reduction(+ : xidosE) &
!$omp& reduction(+ : tdossum) &
!$omp& reduction(+ : tensum)
       ipe=0
       npe=0
!      ipe=omp_get_thread_num()+1
!      npe=omp_get_num_threads()
!      write(6,*)'ipe,npe=',ipe,npe
!$omp do schedule(static)
        do jsv=1,noav
          nss=atm_element(jsv)
          nval2=nval(nss)
          do ja=1,nval2
            nrecl=kr_dim_str(ja,jsv)
            do irec=nrecl,1,-1
              rRd=wt_kr(irec,ja,jsv)
              rEd=eig_kr(irec,ja,jsv)
              yyy=xbeta*(rEd-xmu0)
              if(yyy .gt. 100.d0) then
                xexp=0.0d0
              elseif(yyy.lt.-100.d0) then
                xexp=1.0d0
              else
                xexp=1.0d0/(dexp(yyy)+1.0d0)
              endif
              tdossum=tdossum +2.0d0*rRd
              tensum =tensum  +2.0d0*rRd*rEd
              xidos0=xidos0   +2.0d0*xexp*rRd
              xidosE=xidosE   +2.0d0*xexp*rRd*rEd
            enddo
          enddo
        enddo
!$omp end do
!$omp end parallel 
!
!
      xqinit=rNelec
      if (dabs(xqinit) .le. 1.0d-10) then
        write(6,*)'ERROR!(CHEM_POT):xqinit=',xqinit
        stop
      endif   
!
      err=(xqinit-xidos0)/xqinit
!     if (iloop .le. 3) then 
      if (i_verbose >= 1) write(*,'(a,i10,3f20.10)')'Bisec.',iloop,xmu0,xidos0, err
!     endif   
      if(abs(err) .le. 1.0d-12) exit
!
      if(err .gt.  0.0d0) then
        xmumn=xmu0 
        xmumx=xmumx
        xmu0=(xmumn+xmumx)*0.5d0
      else
        xmumn=xmumn
        xmumx=xmu0 
        xmu0=(xmumn+xmumx)*0.5d0
      endif
!       
      enddo
!
      write(6,'(a,i20,2f30.20)')'chem_pot [au,eV]=',iloop,xmu0,xmu0*ev4au
      chemical_potential=xmu0
!
!     write(*,*)'Set the chem pot as a quite large value'
!     xmu0=1.0d20
!     write(*,*)'    chem_pot=',xmu0
!     chemical_potential=xmu0
!
  end subroutine set_chemical_potential
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Plot LDOS in KRGL method (optional)
!
!
  subroutine plot_ldos_in_krgl
    use M_lib_phys_const,     only : pi !(constant)
    use M_config,     only : config !(unchanged )
    use elses_mod_phys_const, only : ev4au  !(unchanged )
!
!   use elses_arr_eig_leg, only : atmp, atmp2, eig2, f_occ  !(unchanged )
!   use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged )
!
    use elses_mod_file_io, only : vacant_unit
    use M_qm_domain,      only : i_verbose, noav, nval, temp_for_electron, &
&                                atm_element, chemical_potential 
!
    use M_lib_math_func,      only : Fermi_Dirac_Func !(function)
    use M_io_ctrl_dos_plot,   only : initial_preparation_for_dos

    implicit none
    integer :: i,j,k, step_count, ierr, jj
    integer :: number_energy_mesh
    integer :: atm_index, orb_index
    integer :: iunit1, iunit2
    real(8) :: sum_energy, sum_occupation
    real(8) :: width_of_peak
    real(8) :: x
    real(8) :: step_fn, weight
    real(8) :: energy_btm, energy_top, de
    real(8), allocatable  :: energy_mesh(:)
    real(8), allocatable  :: nos(:)
    real(8), allocatable  :: enos(:)
    real(8), allocatable  :: pnos(:,:)
    real(8), allocatable  :: lnos(:,:)
!
    integer :: imode
    integer :: kr_dim_max
    integer :: jsv, nss, nval2, ja, nrecl, irec
    real(8) :: rRd, rEd
    real(8) :: nos_prv, nos_tmp, dos_tmp
    real(8) :: enos_prv, enos_tmp, edos_tmp
!
    real(8) :: energy_inf, nos_inf, enos_inf
!
    kr_dim_max=config%calc%solver%dimension
!   width_of_peak=temp_for_electron
!
    if (i_verbose >= 1) then
      write(*,*)'@@ plot_ldos_in_krgl'
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call initial_preparation_for_dos(imode, number_energy_mesh, energy_btm, energy_top, width_of_peak)
    if (imode == 0) then
      if (i_verbose >= 1) write(*,*)'.... is skipped'
      return
    endif   
!
    if (i_verbose >= 1) then
       write(*,*)' number of energy mesh  =',number_energy_mesh
       write(*,*)' energy bottom   [eV]   =',energy_btm*ev4au
       write(*,*)' energy top      [eV]   =',energy_top*ev4au
       write(*,*)' width of peak      [eV] =',width_of_peak*ev4au
       write(*,*)' chemical potential [au] =',chemical_potential
       write(*,*)' chemical potential [eV] =',chemical_potential*ev4au
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate (energy_mesh(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!    
    allocate (nos(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    nos(:)=0.0d0
!
    allocate (enos(0:number_energy_mesh),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    enos(:)=0.0d0
!
!   allocate (pnos(0:number_energy_mesh,n_tot_base),stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in alloc:pnos'
!   pnos(:,:)=0.0d0
!
!   allocate (lnos(0:number_energy_mesh,noav),stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in alloc:lnos'
!   lnos(:,:)=0.0d0
!
    de=(energy_top-energy_btm)/dble(number_energy_mesh)
!      -----> de : energy interval [au]
!
    do jj=0,number_energy_mesh
      energy_mesh(jj)=energy_btm+de*dble(jj)
!      -----> energy mesh in au
    enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Total DOS
!
    iunit1=vacant_unit()
    open(iunit1, file='output_tdos.txt', status='unknown')
    iunit2=vacant_unit()
!   open(iunit2, file='output_edosinfo-edos-krgl.txt', status='unknown')
!
    energy_inf=10.0d60
!
    nos_inf=0.0d0
    enos_inf=0.0d0
    do jsv=1,noav
      nss=atm_element(jsv)
      nval2=nval(nss)
      do ja=1,nval2
        nrecl=kr_dim_max
        do irec=nrecl,1,-1
          rRd=wt_kr(irec,ja,jsv)*2.0d0
          rEd=eig_kr(irec,ja,jsv)
          do jj=0,number_energy_mesh
            x=(rEd-energy_mesh(jj))/width_of_peak
            step_fn=Fermi_Dirac_Func(x)
            nos(jj)  =  nos(jj)+step_fn*rRd
            enos(jj) = enos(jj)+step_fn*rRd*rEd
          enddo
          x=(rEd-energy_inf)/width_of_peak
          step_fn=Fermi_Dirac_Func(x)
          nos_inf  = nos_inf  + step_fn*rRd
          enos_inf = enos_inf + step_fn*rRd*rEd
        enddo
      enddo
    enddo
!
    write(*,*) 'nos_inf =',nos_inf, energy_inf
    write(*,*) 'enos_inf=',enos_inf, energy_inf
!
    nos_prv=nos(0)
    enos_prv=enos(0)
    do jj=0,number_energy_mesh
      nos_tmp =nos(jj)
      enos_tmp=enos(jj)*ev4au
      dos_tmp=(nos_tmp-nos_prv)/(de*ev4au)
!     edos_tmp=(enos_tmp-enos_prv)/de
      write(iunit1,'(a,i7,4f20.10)')'tdos=',jj,energy_mesh(jj)*ev4au,dos_tmp,nos_tmp,enos_tmp
!     write(iunit2,'(a,i7,3f20.10)')'edos=',jj,energy_mesh(jj),edos_tmp,enos_tmp
      nos_prv =nos_tmp
      enos_prv=enos_tmp
    enddo   
!    
    close(iunit1)
!   close(iunit2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Deallocation
!
    deallocate (nos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
    deallocate (enos,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
!   deallocate (pnos,stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
!   deallocate (lnos,stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
  end subroutine plot_ldos_in_krgl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine reset_dsij
     use M_qm_domain
     implicit none
     integer :: ict4h
     integer :: jsv1, jsv2, jsd, ja1, ja2
!
     write(*,*)'@@ RESET_DSIJ (for developper)'
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
  end subroutine reset_dsij
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_matrix(mat,result,imode_check_mat)
!
!    Check the input matrix, A, for 
!       (1) the deviation from the unit matrix ( A = I ) (imode=1)
!       (2) the deviation from the symmetrix matrix ( A^t = A ) (imode=2)
!
    implicit none
    real(8),         intent(in)  :: mat(:,:)
    real(8),         intent(out) :: result
    integer,         intent(in)  :: imode_check_mat
    real(8), allocatable         :: wrk_mat(:,:)
    integer :: m, ierr
    integer :: i,j
!    
    m=size(mat,1)
!
    allocate (wrk_mat(m,m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
    wrk_mat(:,:)=0.0d0
!
    if (imode_check_mat == 1) then
      wrk_mat(:,:)=dabs(mat(:,:))
      do i=1,m
        wrk_mat(i,i)=mat(i,i)-1.0d0
      enddo
      result=maxval(wrk_mat)
    endif  
!
    if (imode_check_mat == 2) then
      do j=1,m
        do i=1,m
          wrk_mat(i,j)=dabs((mat(i,j)-mat(j,i)))
        enddo
      enddo
      result=maxval(wrk_mat)
    endif  
!
    deallocate (wrk_mat)
!
  end subroutine check_matrix
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Symmetrize the DM and EDM matrices
!       ----> Note: An old routine 'elses_sym_dbij' is used.
!
  subroutine sym_dbij_and_dpij
    use M_qm_domain,        only : dbij, dpij
    use elses_arr_dbij2,    only : dbij2
    implicit none
    real(8), allocatable :: wrk_mat(:,:,:,:)
    integer :: mat_dim(4)
    integer :: j, ierr

!
    write(*,*)'@@ sym_dbij_and_dpij'
!
    if (allocated(dbij)  .eqv. .false.) then
      write(*,*)'ERROR:DBIJ is not allocated'
      stop
    endif   
!
    if (allocated(dpij)  .eqv. .false.) then
      write(*,*)'ERROR:DPIJ is not allocated'
      stop
    endif   
!
    if (allocated(dbij2)  .eqv. .false.) then
      write(*,*)'ERROR:DBIJ2 is not allocated'
      stop
    endif   
!
    do j=1,4
      mat_dim(j)=size(dbij,j)
      write(*,*)'matrix dimension =',j,mat_dim(j)
    enddo  
!
    allocate (wrk_mat(mat_dim(1),mat_dim(2),mat_dim(3),mat_dim(4)),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    wrk_mat(:,:,:,:)=dbij(:,:,:,:)
!      ----> dbij is copied into wrk_mat
!
    dbij(:,:,:,:)=dpij(:,:,:,:)
    call elses_sym_dbij
    dpij(:,:,:,:)=dbij(:,:,:,:)
!       -----> dpij is symmetrized
!
    dbij(:,:,:,:)=wrk_mat(:,:,:,:)
    call elses_sym_dbij
!
    deallocate (wrk_mat,stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
!   write(*,*)'Stop manually'
!   stop
!
  end subroutine sym_dbij_and_dpij
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Get the number of non-zero elements for a vector
!
  subroutine get_num_of_nonzeros(vect,num,num2)
    implicit none
    real(8), intent(in)  :: vect(:)
    integer, intent(out) :: num, num2
    real(8),  parameter  :: epsilon=1.0d-10 
    integer              :: vect_size, n_count, j
!
    vect_size=size(vect,1)
    n_count=0
!
    do j=1,vect_size
      if ( dabs(vect(j)) > epsilon ) then
        n_count=n_count+1
      endif
    enddo
!
    num  = n_count
    num2 = vect_size
!
  end subroutine get_num_of_nonzeros

!
end module M_la_krgl_main

