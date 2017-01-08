!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_matvec_dens
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
  public :: make_mat_dens
  public :: calc_u_su_hu_dens
! public :: calc_u_su_hu_dens_dum
  public :: cg_s_mat_dens
!
! public :: matvec_mul_dens
! public :: matvec_mul_dens_dum
!
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_mat_dens(jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, mat_h, mat_s, &
&                            mat_dens_h_s )
!
    use M_qm_domain,     only : nval, atm_element  ! (unchanged)
!
    implicit none
!   real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
!   real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: jsv4jsk(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_h(:,:,:,:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_s(:,:,:,:)
    real(DOUBLE_PRECISION), intent(out) :: mat_dens_h_s(:,:,:) ! matrices H and S in the dense-matrix form
!
!   logical, parameter :: debug_mode = .false.
!   logical, parameter :: data_damp = .true.
!   real(DOUBLE_PRECISION), parameter  :: eps_value = 1.0d-10
!
    integer :: jsk2, jsk1, jsd, jsv2, jsv1, ja2, ja1
    integer :: num_atom_proj
    integer :: jjkset1, jjkset2, jjk1, jjk2
    integer :: nval1, nval2
    integer :: iloop, n_count, ierr, mat_dim
    integer :: i, j, n
    real(DOUBLE_PRECISION)  :: mat_value
!
!   write(*,*)'matvec_mul_dstm'
!
    num_atom_proj=size(mat_h, 4)
!
    mat_dens_h_s(:,:,:)=0.0d0
!
    do jsk2=1, num_atom_proj
      jjkset2=jjkset(jsk2)
      jsv2=jsv4jsk(jsk2)
      do ja2=1,nval(atm_element(jsv2))
        jjk2=jjkset2+ja2
        do jsd=1, booking_list_dstm_len(jsk2)
          jsk1=booking_list_dstm(jsd, jsk2)
!         if (debug_mode) then
!            if ( ( jsk1 <= 0 ) .or. (jsk1 > num_atom_proj)) then
!             stop 'ERROR(matvec_mul_dstm)'
!            endif
!         endif
          jsv1=jsv4jsk(jsk1)
          jjkset1=jjkset(jsk1)
          do ja1=1,nval(atm_element(jsv1))
            jjk1=jjkset1+ja1
            mat_dens_h_s(jjk1, jjk2,1)=mat_h(ja1,ja2,jsd,jsk2)
            mat_dens_h_s(jjk1, jjk2,2)=mat_s(ja1,ja2,jsd,jsk2)
          enddo
        enddo
      enddo
    enddo
!
!   if (data_damp) then
!     n=size(mat_dens_h_s, 1) 
!     do i=1, n
!       do j=1,n
!         if ( dabs(mat_dens_h_s(i, j, 1)) + dabs(mat_dens_h_s(i, j, 2)) > eps_value  ) then
!           write(*,*)'i,j, H(i,j), S(i,j)=', i,j, mat_dens_h_s(i, j, 1:2)
!         endif
!       enddo
!     enddo
!   endif
!
!   stop 'Stop manually'
!
  end subroutine make_mat_dens
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with DSTM matrix : A (= H or S)
!      (vect_out): = A (vect_in)
!
  subroutine matvec_mul_dens(vect_in, vect_out, mat_data, id_of_my_omp_thread)
!   use M_qm_domain,     only : nval, atm_element  ! (unchanged)
!
    use M_lib_mpi_wrapper,     only : mpi_wrapper_wtime        !(routine)
    use M_lib_timer_in_thread, only : matvec_timer_in_thread   !(CHANGED)
    implicit none
    logical, parameter                  :: use_blas = .true.
!   logical, parameter                  :: use_blas = .false.
    real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_data(:,:)
    real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: id_of_my_omp_thread
    integer n
    real(DOUBLE_PRECISION), parameter :: alpha=1.0d0
    real(DOUBLE_PRECISION), parameter :: beta=0.0d0
    real(DOUBLE_PRECISION)              :: time_data1, time_data2
!
    if (allocated(matvec_timer_in_thread)) then
      call mpi_wrapper_wtime(time_data1)
    endif
!
    vect_out(:)=0.0d0
    if (use_blas) then 
      n=size(vect_in,1)
      call dsymv('U', n, alpha, mat_data, n, vect_in, 1, beta, vect_out, 1)
    else
      vect_out=matmul(mat_data, vect_in)
    endif
!
    if (allocated(matvec_timer_in_thread)) then
      call mpi_wrapper_wtime(time_data2)
      matvec_timer_in_thread(id_of_my_omp_thread+1,1)=matvec_timer_in_thread(id_of_my_omp_thread+1,1)+(time_data2-time_data1)
      matvec_timer_in_thread(id_of_my_omp_thread+1,2)=matvec_timer_in_thread(id_of_my_omp_thread+1,2)+1.0d0
    endif
!
  end subroutine matvec_mul_dens
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_u_su_hu_dens(u,su,hu,norm_factor, mat_dens_h_s, ierr, id_of_my_omp_thread)
!
    implicit none
    real(8),       intent(inout)  :: u(:)
    real(8),         intent(out)  :: su(:)
    real(8),         intent(out)  :: hu(:)
    real(8),         intent(out)  :: norm_factor
    integer,         intent(out)  :: ierr
    real(DOUBLE_PRECISION), intent(in) :: mat_dens_h_s(:,:,:)
    integer,         intent(in)   :: id_of_my_omp_thread
!
!   integer,                intent(in) :: irp(:)
!   integer,                intent(in) :: icol(:)
!   real(DOUBLE_PRECISION), intent(in) :: val(:,:)
!
    real(8) :: ddd
!
!   stop 'Stop manually'
!
!     call matvec_mul_crs(u, su, irp, icol, val(:,2))
      call matvec_mul_dens(u, su, mat_dens_h_s(:,:,2), id_of_my_omp_thread)
!               ---> su : S |m_n> (non-normalized vector)
!
!
!     call matvec_mul_crs(u, hu, irp, icol, val(:,1))
      call matvec_mul_dens(u, hu, mat_dens_h_s(:,:,1), id_of_my_omp_thread)
!               ---> hu : H |m_n> (non-normalized vector)
!
      ddd=dot_product(u(:),su(:))
!
      if (ddd < 0.0d0) then
        write(*,*)'ERROR:(calc_u_su_hu_dstm): <u|S|u>=',ddd
        norm_factor=ddd  ! norm_factor is used as uSu in this case !!!
        ierr=2
        return
        stop
      endif
!
      ddd=dsqrt(ddd)
!      
      ierr=0
      if (dabs(ddd) < 1.0d-10) then
!        write(*,*)'ERROR:calc_u_su_hu'
         write(*,*)'Info(calc_u_su_hu_dstm):Too small normalization factor:=',ddd
         ierr=1
         return
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

  end subroutine calc_u_su_hu_dens
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cg_s_mat_dens(b,x,eps,kend, mat_dens, id_of_my_omp_thread)
!
!    ---> Solve linear eq.: S x = b (REAL VARIALE)
!          convergence criteria : 
!             log_10 | r / b | < EPS
!      kend : (in input ) the maximum iteration number
!           : (in output) the last iteration number that is executed
!
!   use M_la_matvec_io, only : matvec_mul
!   use M_qm_projection, only : matvec_mul_proj ! (routine)
!
    implicit none
    real(8),       intent(in)     :: b(:)
    real(8),       intent(inout)  :: x(:)
    real(8),       intent(inout)  :: eps
    integer,        intent(inout) :: kend
    real(DOUBLE_PRECISION), intent(in) :: mat_dens(:,:)
    integer,        intent(in)    :: id_of_my_omp_thread

!   integer,           intent(in) :: jsv4jsk(:)
!   integer,           intent(in) :: jsk4jsv(:)
!   integer,           intent(in) :: jjkset(:)
!
!   integer,                intent(in)  :: booking_list_dstm(:,:)
!   integer,                intent(in)  :: booking_list_dstm_len(:)
!   real(DOUBLE_PRECISION), intent(in)  :: overlap_dstm(:,:,:,:)
!
!   integer,                intent(in) :: irp(:)
!   integer,                intent(in) :: icol(:)
!   real(DOUBLE_PRECISION), intent(in) :: val(:,:)
!
    real(8), allocatable          :: r(:), p(:), ap(:)
    real(8)                       :: alpha, beta, rho0, rho1, tbs
    real(8)                       :: hg, hal, pap,hnor
    integer                       :: m, m2, i, ierr
    integer                       :: kk, kend_in
!
!
    kend_in=kend
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get matrix size : m
!
    m=size(b,1)
    m2=size(x,1)
    if (m /= m2) stop 'Parameter Mismatch:m,m2'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix allocation
!
    allocate (r(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (p(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (ap(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iniital value
!
    beta = 0.0d0
!
    call matvec_mul_dens(x, ap, mat_dens, id_of_my_omp_thread)
!        ---> Mat-vec multiplication : ap = S x
!
    r(:) = b(:)-ap(:)  ! r_0 = b - A x_0
    p(:) = r(:)        ! p_0 = r_0
!
    rho0 = dot_product(r(:),r(:))
    hnor = dot_product(b(:),b(:))
    hnor=dlog10(hnor)/2.0d0
!
    kk=-1
    hal = rho0
    hg = dlog10(hal)/2.0d0 - hnor
    if(hg .le. eps) then
      eps= dlog10(hal)/2.0d0 - hnor
      kend=kk
      return
    endif    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main CG loop
!
    do kk = 0, kend_in
!
      hal = rho0
      hg = dlog10(hal)/2.0d0 - hnor
!
      if(hg .le. eps) then
!        --  If converged, the true residual is calculated and exit---
!
         call matvec_mul_dens   (x, ap, mat_dens, id_of_my_omp_thread)
!          ---> Mat-vec multiplication : ap = S x
!
        r(:) = b(:)-ap(:)
        rho0 = dot_product(r(:),r(:))
        hal = rho0 + 1.0d-100
        eps= dlog10(hal)/2.0d0 - hnor
        kend=kk
        return
      end if
!
      p(:) = r(:) + beta*p(:)
!
      call matvec_mul_dens  (p, ap, mat_dens, id_of_my_omp_thread)
!        ---> Mat-vec multiplication : ap = S p
      pap = dot_product(p(:),ap(:))
!
      alpha = rho0/pap
      x(:) = x(:) +alpha*p(:)
      r(:) = r(:) -alpha*ap(:)
      rho1=dot_product(r(:),r(:))
      beta = rho1/rho0
      rho0 = rho1
!
    enddo
!
    return
!
!
  end subroutine cg_s_mat_dens
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module M_la_matvec_dens



