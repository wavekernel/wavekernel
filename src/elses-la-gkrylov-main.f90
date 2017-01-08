!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_gkrylov_main
!
  implicit none
!
  private
  public :: gkrylov_main  ! OUTDATED ROUTINE
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gkrylov_main(j_src,b,jsv4jsk,jsk4jsv,jjkset,prc_index,scheme_mode, & 
&         s_inv_e_j_wrk,dm_wrk,u_hst_wrk, v_mat_kr, eig_wrk, u_b_hst_wrk, wt_kr_wrk, kr_dim_max)
!
!
!
    use  M_la_krgl_main, only : calc_u_su_hu, get_num_of_nonzeros,  &
&                             check_matrix, schmidt_orth_w_s_mat !(routine)
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
    use M_lib_math_func,  only : Fermi_Dirac_Func  !(function)
!
!   use M_la_gkrylov_hist, only : proj_list_unchanged, s_inv_e_j, rho_e_j !(CHANGED)
!   use M_la_gkrylov_hist, only : u_hst_str, v_mat_kr_str !(CHANGED)
!   use M_la_gkrylov_hist, only : u_b_hst_str             !(CHANGED)
!
    use M_eig_standard,          only : eig_standard !(routine)
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
    real(8),        intent(inout) :: s_inv_e_j_wrk(:)
    real(8),        intent(inout) :: dm_wrk(:,:)
    real(8),        intent(inout) :: u_hst_wrk(:,:) 
    real(8),        intent(inout) :: v_mat_kr(:,:)
    real(8),        intent(inout) :: eig_wrk(:)
    real(8),        intent(inout) :: u_b_hst_wrk(:)
    real(8),        intent(inout) :: wt_kr_wrk(:)
    integer,        intent(inout) :: kr_dim_max
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
!
    real(8), allocatable ::  den_mat_kr(:,:)
    real(8), allocatable :: eden_mat_kr(:,:)
!
    real(8), allocatable ::  den_mat_kr_b(:)
    real(8), allocatable :: eden_mat_kr_b(:)
!
!   real(8), allocatable ::   dm_loc(:)
!   real(8), allocatable ::  edm_loc(:)
!

    integer :: ierr
!
    integer :: mat_dim
    integer :: m_int, m_int2
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
    integer :: kr_dim, kr_dim_max_input, kr_dim_max_loop
!   integer :: kr_dim, kr_dim_max, kr_dim_max_input, kr_dim_max_loop
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
    character(len=72) :: msg
!
    write(*,*)'ERROR:Outdated routine is called :gkrylov_main'
    stop
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
!     if ((atm_index >= 59) .and. (atm_index <= 60)) then
!       i_show=1
!     endif
    endif  
!
    if (i_show >= 1) then
      write(*,*)'atm_index=',atm_index
      write(*,*)'orb_index=',orb_index
      write(*,*)'the position of non-zero element in the initial vector=',ini_elem
      write(*,*)'mode for large memory=',i_kr_hst_str
    endif  
!
!
!   flag_for_s_inv=.false.
!   if ( allocated(s_inv_e_j) ) then
!      flag_for_s_inv=.true.
!   endif   
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
    m_int=size(dm_wrk,1)
    if (i_show >= 1) then
      write(*,*)'m_int=',m_int
    endif  
!
    if (i_check_mode == 1) then
      call get_interac_list_num_proj(j_src,m_int2)
      if (m_int /= m_int2) then
        write(*,*)' ERROR:m_int,m_int2=',m_int,m_int2
        stop
      endif   
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
!   if (allocated(s_inv_e_j)) then
!     flag_for_s_inv = .true.
!   else
!     flag_for_s_inv = .false.
!   endif
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Set the initial basis : |l_1> = S^{-1} |b> ( non-normalized vector )
!
    eps_c=-8.0d0
    max_ite_for_cg_loop=100
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
!     kr_dim_max=kr_dim_str(orb_index,atm_index)
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
          u_hst(jjk,kr_dim) =  u_hst_wrk(j, kr_dim)
!         u_hst(jjk,kr_dim) =  u_hst_str(j, kr_dim, orb_index, atm_index)
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
        u_b_hst_wrk(kr_dim) =  u_hst(ini_elem,kr_dim)
!       u_b_hst_str(kr_dim, orb_index, atm_index) =  u_hst(ini_elem,kr_dim)
        do j=1,m_int
          jjk=interaction_list(j)
           u_hst_wrk(j, kr_dim) =  u_hst(jjk,kr_dim)
!          u_hst_str(j, kr_dim, orb_index, atm_index) =  u_hst(jjk,kr_dim)
        enddo  
      endif   
!
      u(:)=hu(:)
      if ((scheme_mode == 'ekrgl') .and. (kr_dim == kr_dim_max_input/2)) then
        if (i_show >= 1) then
          write(*,*) ' Multiple KR space : kr_dim=',kr_dim
        endif  
        ite_cg=max_ite_for_cg_loop
        eps=eps_c
        u(1:mat_dim)=s_inv_e_j_wrk(1:mat_dim)
        call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
        s_inv_e_j_wrk(1:mat_dim)=u(1:mat_dim)
        if (i_show >= 1) then
          write(*,'(a,2i10,f30.20,i10)')'Sinv: prc_index, j_src, eps, ite_cg=', prc_index, j_src,eps,ite_cg
        endif   
      endif  
!
!
!         if (flag_for_s_inv) then
!           if (i_show >= 1) then
!             write(*,*) ' S_inv_e_j is allocated'
!             write(*,*) ' proj_list_unchanged=',atm_index, proj_list_unchanged(atm_index)
!           endif   
!           if (prc_index == 1) then
!             u(1:mat_dim)=s_inv_e_j(1:mat_dim,j_src)
!             call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
!             s_inv_e_j(1:mat_dim,j_src)=u(1:mat_dim)
!             if (i_show >= 1) then
!               write(*,*)'Sinv(1): j_src, eps, ite_cg=',j_src,eps,ite_cg
!             endif   
!           else
!             u(1:mat_dim)=s_inv_e_j(1:mat_dim,j_src)
!             call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
!             s_inv_e_j(1:mat_dim,j_src)=u(1:mat_dim)
!             if (i_show >= 1) then
!               write(*,*)'Sinv(2): j_src, eps, ite_cg=',j_src,eps,ite_cg
!             endif  
!           endif 
!
!         else
!           u(:)=0.0d0
!           call cg_s_mat_proj_r2(b, u, eps, ite_cg, jsv4jsk, jsk4jsv, jjkset)
!         endif   
!
!       endif
!     endif  
!
      if (kr_dim /= kr_dim_max_input) then
         call schmidt_orth_w_s_mat(u,u_hst,su_hst,kr_dim)
!        --->  u(:) :   |l_{n+1}> 
      endif   
!
    enddo
!   kr_dim_str(orb_index,atm_index)=kr_dim_max
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
!           kr_dim_str(orb_index,atm_index)=kr_dim_max
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
!   allocate (v_mat_kr(kr_dim_max,kr_dim_max),stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in alloc (v_mat_kr)'
!   v_mat_kr(:,:)=0.0d0
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
!     allocate (eig_wrk(kr_dim_max),stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in alloc (eig_wrk)'
!     eig_wrk(:)=0.0d0
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
      enddo   
!
      if (i_show >= 1) then
        write(*,*)' Tr[H_kr] (eV) =',ddd*ev4au
      endif  
!
    endif
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Eigen states for the reduced matrix
!
    if ( (prc_index == 1) .or. (i_kr_hst_str == 0) ) then
      msg='' 
      call eig_standard(v_mat_kr, eig_wrk, h_mat_kr, 1, size(h_mat_kr,1), msg)
      if (trim(msg) /= '') then
        write(*,*)' ERROR(eig_standard):',trim(msg)
        write(*,*)'  atm_index, orb_index =',atm_index, orb_index
        stop
      endif
!     call calc_eig(v_mat_kr,eig_wrk,h_mat_kr,s_mat_kr)
!
!     eig_kr(1:kr_dim_max,orb_index,atm_index)=eig_wrk(1:kr_dim_max)
!
!     if ( i_kr_hst_str == 1 ) then
!       if (.not. allocated(v_mat_kr_str)) stop 'Alloc. error (v_mat_kr_str)'
!       do kr_dim=1,kr_dim_max
!         v_mat_kr_str(1:kr_dim_max,kr_dim,orb_index,atm_index)=v_mat_kr(1:kr_dim_max,kr_dim)
!       enddo  
!     endif   
!
!   else
!
!     if (.not. allocated(v_mat_kr_str)) stop 'Alloc. error (v_mat_kr_str)'
!     do kr_dim=1,kr_dim_max
!       v_mat_kr(1:kr_dim_max,kr_dim)=v_mat_kr_str(1:kr_dim_max,kr_dim,orb_index,atm_index)
!     enddo  
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
        wt_kr_wrk(al)=ddd*v_mat_kr(1,al)/norm_factor1
!       wt_kr(al,orb_index,atm_index)=ddd*v_mat_kr(1,al)/norm_factor1
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
          ddd=ddd+wt_kr_wrk(al)
!         ddd=ddd+wt_kr(al,orb_index,atm_index)
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
      wt_kr_wrk(al)=ddd*ddd1
!     wt_kr(al,orb_index,atm_index)=ddd*ddd1
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
        ddd=ddd+wt_kr_wrk(al)
!       ddd=ddd+wt_kr(al,orb_index,atm_index)
        write(*,'(a,i10,2f30.20)') 'Eig [eV] =', al, eig_wrk(al)*ev4au,wt_kr_wrk(al)
!       write(*,'(a,i10,2f30.20)') 'Eig [eV] =', al, eig_wrk(al)*ev4au,wt_kr(al,orb_index,atm_index)
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
!     allocate ( dm_loc(m_int),stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in alloc (dm_loc)'
!     dm_loc(:)=0.0d0
!
!     allocate ( edm_loc(m_int),stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in alloc (edm_loc)'
!     edm_loc(:)=0.0d0
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
        ddd=(eig_wrk(al)-xmu)*xbeta
!       ddd=(eig_kr(al,orb_index,atm_index)-xmu)*xbeta
        f_occ=Fermi_Dirac_Func(ddd)
!       if (i_show >= 1) then
!         write(*,*)'f_occ=',al,f_occ
!       endif  
        do m=1,kr_dim_max
          ddd=eig_wrk(al)
!         ddd=eig_kr(al,orb_index,atm_index)
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
!     allocate ( dm_wrk(m_int,2),stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in alloc (dm_loc)'
!     dm_wrk(:,:)=0.0d0
!
      do j=1,m_int
        jjk=interaction_list(j)
         dm_wrk(j,1) = dot_product( u_hst(jjk,:),  den_mat_kr_b(:) )
         dm_wrk(j,2) = dot_product( u_hst(jjk,:), eden_mat_kr_b(:) )
      enddo
!   
!     call convert_dm_and_edm(atm_index,orb_index,dm_loc,edm_loc,jsv4jsk,jsk4jsv,jjkset)
!
!     if ( allocated(rho_e_j) ) then
!       if (i_show == 1) then
!         write(*,*)'  rho_e_j is stored:j_src=', j_src, atm_index, orb_index
!       endif   
!       if (i_check_mode == 1) then
!         if (size(rho_e_j,1) < mat_dim) then
!           write(*,*)'ERROR: krgl_main (rho_e_j):size1=',size(rho_e_j,1)
!           stop
!         endif   
!         if (size(rho_e_j,2) < j_src) then
!           write(*,*)'ERROR: krgl_main (rho_e_j):size2=',size(rho_e_j,2)
!           stop
!         endif   
!       endif  
!       do al=1,kr_dim_max
!         rho_e_j(1:mat_dim,j_src)= rho_e_j(1:mat_dim,j_src) + u_hst(1:mat_dim, al)*den_mat_kr_b(al)
!       enddo
!     endif
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
!     deallocate ( dm_loc,stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in dealloc'
!
!     deallocate ( edm_loc,stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in dealloc'
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
!   deallocate (v_mat_kr ,stat=ierr)
!   if (ierr /= 0) stop 'Abort:ERROR in dealloc (v_mat_kr)'
!
!   if (allocated(eig_wrk)) then
!     deallocate (eig_wrk ,stat=ierr)
!     if (ierr /= 0) stop 'Abort:ERROR in dealloc (eig_wrk)'
!   endif  
!
!
  end subroutine gkrylov_main
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_la_gkrylov_main


