!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_matvec_crs
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
  public :: make_mat_crs
  public :: get_num_nonzero_elem
  public :: calc_u_su_hu_crs
  public :: cg_s_mat_crs
! public :: cg_s_mat_crs_dum
!
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make matrix in CSR format
!
!
  subroutine get_num_nonzero_elem(jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, mat_dstm, nnz)
!
    use M_qm_domain,     only : nval, atm_element  ! (unchanged)
!
    implicit none
    integer,                intent(in)  :: jjkset(:)
    integer,                intent(in)  :: jsv4jsk(:)
    integer,                intent(in)  :: booking_list_dstm(:,:)
    integer,                intent(in)  :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(in)  :: mat_dstm(:,:,:,:)
    integer,                intent(out) :: nnz
!
    logical, parameter :: debug_mode = .false.
!   logical, parameter :: debug_mode = .true.
!
    integer :: jsk2, jsk1, jsd, jsv2, jsv1, ja2, ja1
    integer :: num_atom_proj
    integer :: jjkset1, jjkset2, jjk1, jjk2
    integer :: nval1, nval2
    integer :: n_count, ierr, n_mat_size
!
!   write(*,*)'matvec_mul_dstm'
!
    num_atom_proj=size(mat_dstm, 4)
!
    n_count=0
    do jsk2=1, num_atom_proj
      jjkset2=jjkset(jsk2)
      jsv2=jsv4jsk(jsk2)
      nval2=nval(atm_element(jsv2))
      do jsd=1, booking_list_dstm_len(jsk2)
        jsk1=booking_list_dstm(jsd, jsk2)
        if (debug_mode) then
          if ( ( jsk1 <= 0 ) .or. (jsk1 > num_atom_proj)) then
            stop 'ERROR(matvec_mul_dstm)'
          endif
        endif
        jsv1=jsv4jsk(jsk1)
        nval1=nval(atm_element(jsv1))
        n_count=n_count+nval1*nval2
      enddo
    enddo
!   write(*,*)'n_count=', n_count
    nnz=n_count
!
  end subroutine get_num_nonzero_elem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_mat_crs(jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, mat_h, mat_s, &
&                            irp, icol, val )
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
    integer,                intent(out) :: irp(:)
    integer,                intent(out) :: icol(:)
    real(DOUBLE_PRECISION), intent(out) :: val(:,:)
!
    logical, parameter :: debug_mode = .false.
!
    integer :: jsk2, jsk1, jsd, jsv2, jsv1, ja2, ja1
    integer :: num_atom_proj
    integer :: jjkset1, jjkset2, jjk1, jjk2
    integer :: nval1, nval2
    integer :: iloop, n_count, ierr, mat_dim
    real(DOUBLE_PRECISION)  :: mat_value
!
!   write(*,*)'matvec_mul_dstm'
!
    num_atom_proj=size(mat_h, 4)
    mat_dim      =size(irp,1)-1
!
    n_count=0
    do jsk2=1, num_atom_proj
      jjkset2=jjkset(jsk2)
      jsv2=jsv4jsk(jsk2)
      do ja2=1,nval(atm_element(jsv2))
        jjk2=jjkset2+ja2
        irp(jjk2)=n_count+1
        do jsd=1, booking_list_dstm_len(jsk2)
          jsk1=booking_list_dstm(jsd, jsk2)
          if (debug_mode) then
            if ( ( jsk1 <= 0 ) .or. (jsk1 > num_atom_proj)) then
              stop 'ERROR(matvec_mul_dstm)'
            endif
          endif
          jsv1=jsv4jsk(jsk1)
          jjkset1=jjkset(jsk1)
          do ja1=1,nval(atm_element(jsv1))
            jjk1=jjkset1+ja1
            n_count=n_count+1
            val(n_count,1) =mat_h(ja1,ja2,jsd,jsk2)
            val(n_count,2) =mat_s(ja1,ja2,jsd,jsk2)
            icol(n_count)=jjk1
          enddo
        enddo
      enddo
    enddo
!   write(*,*)'jjk2, mat_dim, n_count=', jjk2, mat_dim, n_count
    irp(mat_dim+1)=n_count+1
!
!
!   stop 'Stop manually'
!
  end subroutine make_mat_crs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix-vector multiplication with CRS format
!      (vect_out): = A (vect_in)
!
  subroutine matvec_mul_crs(vect_in, vect_out, irp, icol, val, id_of_my_omp_thread)
    use M_lib_mpi_wrapper,     only : mpi_wrapper_wtime        !(routine)
    use M_lib_timer_in_thread, only : matvec_timer_in_thread   !(CHANGED)
    implicit none
    real(DOUBLE_PRECISION), intent(in)  :: vect_in(:)
    real(DOUBLE_PRECISION), intent(out) :: vect_out(:)
    integer,                intent(in)  :: irp(:)
    integer,                intent(in)  :: icol(:)
    real(DOUBLE_PRECISION), intent(in)  :: val(:)
    integer,                intent(in)  :: id_of_my_omp_thread
    integer                             :: n
    integer                             :: i, j_ptr
    real(DOUBLE_PRECISION)              :: s
!
    real(DOUBLE_PRECISION)              :: time_data1, time_data2
!
!   write(*,*)'INFO:MATVEC-CRS:id_of_my_omp_thread =', id_of_my_omp_thread
!
    if (allocated(matvec_timer_in_thread)) then
      call mpi_wrapper_wtime(time_data1)
    endif
!
    n=size(vect_in,1)
!
    do i=1,n
      s=0.0d0
      do j_ptr=irp(i), irp(i+1)-1
        s=s+val(j_ptr) * vect_in(icol(j_ptr))
      enddo
      vect_out(i)=s
    enddo
!
    if (allocated(matvec_timer_in_thread)) then
      call mpi_wrapper_wtime(time_data2)
      matvec_timer_in_thread(id_of_my_omp_thread+1,1)=matvec_timer_in_thread(id_of_my_omp_thread+1,1)+(time_data2-time_data1)
      matvec_timer_in_thread(id_of_my_omp_thread+1,2)=matvec_timer_in_thread(id_of_my_omp_thread+1,2)+1.0d0
    endif
!
  end subroutine matvec_mul_crs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_u_su_hu_crs(u,su,hu,norm_factor, irp, icol, val, ierr, id_of_my_omp_thread )
!
    implicit none
    real(8),       intent(inout)  :: u(:)
    real(8),         intent(out)  :: su(:)
    real(8),         intent(out)  :: hu(:)
    real(8),         intent(out)  :: norm_factor
    integer,         intent(out)  :: ierr
!
    integer,                intent(in) :: irp(:)
    integer,                intent(in) :: icol(:)
    real(DOUBLE_PRECISION), intent(in) :: val(:,:)
    integer,                intent(in) :: id_of_my_omp_thread
!
    real(8) :: ddd
!
!     call matvec_mul_crs_dum(u, su, jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, overlap_dstm)
      call matvec_mul_crs(u, su, irp, icol, val(:,2), id_of_my_omp_thread)
!               ---> su : S |m_n> (non-normalized vector)
!
!
!     call matvec_mul_crs_dum(u, hu, jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, ham_tot_dstm)
      call matvec_mul_crs(u, hu, irp, icol, val(:,1), id_of_my_omp_thread)
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

  end subroutine calc_u_su_hu_crs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cg_s_mat_crs(b,x,eps,kend, irp, icol, val, id_of_my_omp_thread)
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
!   integer,           intent(in) :: jsv4jsk(:)
!   integer,           intent(in) :: jsk4jsv(:)
!   integer,           intent(in) :: jjkset(:)
!
!   integer,                intent(in)  :: booking_list_dstm(:,:)
!   integer,                intent(in)  :: booking_list_dstm_len(:)
!   real(DOUBLE_PRECISION), intent(in)  :: overlap_dstm(:,:,:,:)
!
    integer,                intent(in) :: irp(:)
    integer,                intent(in) :: icol(:)
    real(DOUBLE_PRECISION), intent(in) :: val(:,:)
    integer,                intent(in) :: id_of_my_omp_thread
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
!   call matvec_mul_crs_dum(x, ap, jjkset, jsv4jsk, booking_list_dstm, booking_list_dstm_len, overlap_dstm)
    call matvec_mul_crs(x, ap, irp, icol, val(:,2), id_of_my_omp_thread)
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
         call matvec_mul_crs    (x, ap, irp, icol, val(:,2), id_of_my_omp_thread)
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
      call matvec_mul_crs    (p, ap, irp, icol, val(:,2), id_of_my_omp_thread)
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
  end subroutine cg_s_mat_crs
!
end module M_la_matvec_crs


