!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_la_gscocg_main
!
  implicit none
!
  private
  public :: gscocg_main
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gscocg_main(j_src,jsv4jsk,jsk4jsv,jjkset,b,z_seed,z_shift,nj,dos)
!
!           ( z_k S - H ) x(z_k) = b
!
!    INPUT 
!      j_src      : index 'j' in the outer loop
!      jsv4jsk(:) : projection list used in mat-vec multiplication
!      jsk4jsv(:) : projection list used in mat-vec multiplication
!      b(m)       : vector 'b' 
!      z_seed     : seed energy points
!      z_shift(n) : shift energy points
!
!
    use M_lib_phys_const, only : pi_math => pi ! (constant)
    use M_la_matvec_io,  only : get_interaction_list_number, & 
&       set_interaction_list, matvec_mul, inner_product_with_s,calc_nj ! (routine)
    use M_la_gscocg_cg_s_mat, only : cg_s_mat_proj             ! (routine)
    use M_wall_clock_time,    only : get_system_clock_time     ! (routine)
    use M_qm_projection,      only : get_interac_list_num_proj, &
&                                    set_interac_list_proj, matvec_mul_proj, &
&                                    inner_product_with_s_proj   ! (routine)
!
    implicit none
    integer,          intent(in)  :: j_src
    integer,          intent(in)  :: jsv4jsk(:)
    integer,          intent(in)  :: jsk4jsv(:)
    integer,          intent(in)  :: jjkset(:)
    complex(8),       intent(in)  :: b(:)
    complex(8),       intent(in)  :: z_seed
    complex(8),       intent(in)  :: z_shift(:)
!
    real(8) :: elapse_time1, elapse_time2, elapse_time3, elapse_time4,time1,time2
!
!  Arrays for 'seed' energy point (imcomplete)
!
    complex(8), allocatable       :: x(:)
!
!  Arrays for 'shift' energy point (imcomplete)
!
!    complex(8), allocatable       :: pi(:)

!MT
	complex(8), allocatable       :: r(:), p(:), q(:), pi_v(:), qs(:), rk(:)
	complex(8), allocatable	      :: sx(:,:),sp(:,:),sigma(:),salpha(:),sbeta(:)
	complex(8), allocatable	      :: pi_n(:),pi(:),pi_o(:)
	complex(8)                    :: alpha, alpha_old, beta, rho0, rho1, pap,salphak,sbetak
	real(8), allocatable	      :: hg2(:)
	real(8)                       :: hg, hal, hnor,eps,eps2,ccc
	integer, allocatable	      :: flg(:),iter(:)
	integer                       :: m, m2, n, i, j, ierr, ind
    integer                       :: kk, kend_in,kend

    complex(8), allocatable       ::sax(:),sbx(:)
    real(8)                       ::hes,tbs,tbs2,Rave1,Rave2,Rave3,Rave4,temp  
    real(8)                       :: temp2, dos(:,:), nj(:)
    real(8), allocatable          :: vect_in(:)
!MT
!
!  Arrays for the interaction list and related quantities
!
    integer :: m_int
    integer, allocatable     :: interaction_list(:)
    complex(8), allocatable  :: sx_int(:,:)
!    ---> A PARTIAL solution vector 'x' 
!                for the shift eqns :sx_int(1:m_int, 1:n)
!         The correspondence between the full solution vector sx(1:m,1:n)
!                 sx_int(ind, k) = sx(i,k)
!                      where i = interaction_list(ind)
!
!  Other variables
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get matrix size
!
    m=size(b,1)
!     ----> Matrix size
    n=size(z_shift,1)
!     ----> Number of energy points
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix allocation
!
    allocate (x(m),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'
!
    allocate (pi(n),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc'

!MT
	allocate (r(m), p(m), q(m), qs(m), rk(m))
	allocate (sx(m,n),sp(m,n),sigma(n),salpha(n),sbeta(n))
	allocate (pi_n(n),pi_o(n), pi_v(n))
	allocate (hg2(n))
	allocate (flg(n),iter(n))
        
        allocate (sax(m),sbx(m))
	

!   allocate every array defined before


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct the interaction list
!
    call get_interac_list_num_proj(j_src,m_int)
    write(*,*)'m_int=',m_int
!
    allocate (interaction_list(m_int),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc1'
    allocate (vect_in(m_int),stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc2'
!   
    call set_interac_list_proj(j_src,interaction_list,jsk4jsv,jjkset)
!    
    write(*,*)' ...set_interac_list_proj is ended'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The main loop of GSCOCG
!        write(*,*)'allocate finished, GSCOCG iteration start'

        kend=2000;eps=-8.0D0;eps2=-8.0D0
        hnor=0.0d0 
	do i=1,n
		pi(i)	=DCMPLX(1.0D0,0.0D0)
		pi_o(i)	=DCMPLX(1.0D0,0.0D0)
		pi_v(i)	=DCMPLX(1.0D0,0.0D0)
	enddo
	alpha_old=DCMPLX(1.0D0,0.0D0)
	do i=1,n
		do j=1,m
			sp(j,i)=DCMPLX(0.0D0,0.0D0)
			sx(j,i)=DCMPLX(0.0D0,0.0D0)
		enddo
	enddo
	do i=1,n
		sigma(i)=z_shift(i)-z_seed !(-0.50D0-DREAL(z_seed))+DBLE(I-1)/DCMPLX(2500.0D0,0.0D0)
		flg(i)=0
	enddo
!!!!!!!!!!!!!!!!

!	WRITE(*,*)'kend 1',kend

	beta=0.0D0;rho0=0.0D0
	do j=1,m
		r(j)=b(j)
		p(j)=DCMPLX(0.0D0,0.0D0)
	enddo

        kend_in=kend
	call cg_s_mat_proj(r,rk,eps,kend_in,jsv4jsk,jsk4jsv,jjkset)
!       stop 'inner CG loop ended'
!
	do j=1,m
		rho0=rho0+r(j)*rk(j)
	enddo

        time1=0.0d0;time2=0.0d0
loop1:	do kk=0,kend
!        exit loop1
        call get_system_clock_time(elapse_time1)
       

		pap=0.0D0;hal=0.0D0
		do j=1,m
			p(j)=rk(j)+beta*p(j)
			hal=hal+conjg(rk(j))*rk(j)
		enddo

		call matvec_mul_proj(p,q,'H',jsv4jsk,jsk4jsv,jjkset)
		call matvec_mul_proj(p,qs,'S',jsv4jsk,jsk4jsv,jjkset)
		do j=1,m
                   q(j)=z_seed*qs(j)-q(j)
		enddo
		do j=1,m
                   pap=pap+p(j)*q(j)
		enddo

                Rave1=0.0d0
                do ind=1,m_int
                   Rave1=Rave1+conjg(rk(interaction_list(ind)))*rk(interaction_list(ind))         
                enddo
                Rave2=0.0d0
                do ind=1,m
                   Rave2=Rave2+conjg(rk(ind))*rk(ind)             
                enddo
             
!Temp            
                temp=0.0d0
                do i=1,n
                   temp2=conjg(pi(i))*pi(i)
                   if ((dabs(temp2) < 1.0d-100) .and. (dabs(temp2) > 1.0d-100)) then
                     write(*,*)' temp2=',temp2
                   endif   
                   temp=temp+1.0d0/temp2
                enddo
                Rave1=Rave1*temp/n
                Rave2=Rave2*temp/n
                if(kk.eq.1)Rave3=Rave1
                if(kk.ge.1)write(48,*)kk,dlog10(Rave1/Rave3)
                if(kk.eq.1)Rave4=Rave2
                if(kk.ge.1)write(49,*)kk,dlog10(Rave2/Rave4)
	
!TEMP               
!                if(kk.ge.1)then
!                   temp=conjg(pi(501))*(pi(501))
!                   write(41,*)kk,dlog10(Rave1/(temp))
!                   temp=conjg(pi(1001))*(pi(1001))
!                   write(42,*)kk,dlog10(Rave1/(temp))
!                   temp=conjg(pi(1501))*(pi(1501))
!                   write(43,*)kk,dlog10(Rave1/(temp)) 
!                   temp=conjg(pi(2001))*(pi(2001))
!                   write(44,*)kk,dlog10(Rave1/(temp))
!                   temp=conjg(pi(501))*(pi(501))
!                   write(45,*)kk,dlog10(Rave2/(temp))
!                   temp=conjg(pi(1001))*(pi(1001))
!                   write(46,*)kk,dlog10(Rave2/(temp))
!                   temp=conjg(pi(1501))*(pi(1501))
!                   write(47,*)kk,dlog10(Rave2/(temp)) 
!                   temp=conjg(pi(2001))*(pi(2001))
!                   write(48,*)kk,dlog10(Rave2/(temp))
!                 endif
!TEMP

                
                hg=DLOG10(hal)/2.0D0  -hnor
                write(*,*)kk,hg    !,Rave1,Rave2
!                stop

!return the whole subroutine if possible
                
		if(hg.lt.eps2)then
!                   write(*,*)'calculate the true residual'
                   
                   do i=1,n
                      write(70+j_src,*)dreal(sigma(i)+dreal(z_seed)),iter(i)
                   enddo   
                   exit loop1
                   do i=1,n
                      
                      call matvec_mul_proj(sx(:,i),sax,'H',jsv4jsk,jsk4jsv,jjkset)
                      call matvec_mul_proj(sx(:,i),sbx,'S',jsv4jsk,jsk4jsv,jjkset)
                      
                      do j=1,m
                      sax(j)=b(j)-(z_shift(i)*sbx(j)-sax(j))
                      enddo

                      hes=0.0d0;hnor=0.0d0
                      do j=1,m
                         
                         hes=hes+conjg(sax(j))*sax(j)
                         
                         hnor=hnor+conjg(b(j))*b(j)
                      enddo 

                      hes=dlog10(hes)/2.0d0-dlog10(hnor)/2.0d0                     

                      write(71,'(2F12.5,I4)')dreal(sigma(i)+dreal(z_seed)),hes,iter(i)
                   enddo
!                   stop
                   exit loop1
                endif   
                   

		alpha=rho0/pap
!                write(*,*)'alpha=',alpha
		do j=1,m
			x(j)=x(j)+alpha*p(j)
		enddo

!!!!!!!!!for shifted systems	
                call get_system_clock_time(elapse_time3)
	
                loop2:	do i=1,n
			pi_n(i)=(DCMPLX(1.0D0,0.0D0)+alpha*sigma(i))*pi(i)+alpha*beta/alpha_old*(pi(i)-pi_o(i))
			sbeta(i)=(pi_o(i)/pi(i))*(pi_o(i)/pi(i))*beta
			salpha(i)=(pi(i)/pi_n(i))*alpha
			ccc=conjg(pi(i))*pi(i)

                        hg2(i)= DLOG10(hal)/2.0 -DLOG10(ccc)/2.0
			if(hg2(i).lt.eps2) flg(i)=1
                        if(flg(i).eq.1) cycle loop2
			iter(i)=kk  

			sbetak=sbeta(i)*pi(i)/pi_v(i)
			do j=1,m
				sp(j,i)=rk(j)+sbetak*sp(j,i)
                        enddo
                        salphak=salpha(i)/pi(i)
                        do j=1,m
				sx(j,i)=sx(j,i)+salphak*sp(j,i)
			enddo

			pi_v(i)=pi(i)
		enddo loop2
		
		do i=1,n
			pi_o(i)=pi(i)
			pi(i)=pi_n(i)
		enddo

                call get_system_clock_time(elapse_time4)

		alpha_old=alpha

		do j=1,m
			r(j)=r(j)-alpha*q(j)
		enddo	        
    
		rho1=0.0D0;kend_in=kend
		call  cg_s_mat_proj(r,rk,eps,kend_in,jsv4jsk,jsk4jsv,jjkset)	
 		do j=1,m
			rho1=rho1+rk(j)*r(j)
		enddo
		beta=rho1/rho0
		rho0=rho1
!                write(*,*)'hnor= ',hnor
!                stop

                call get_system_clock_time(elapse_time2)
!                write(*,*)'total_time=',elapse_time2-elapse_time1
!                write(*,*)'inner_time=',elapse_time4-elapse_time3
!                write(60,*)kk,(elapse_time4-elapse_time3),(elapse_time2-elapse_time1)
                time1=time1+(elapse_time2-elapse_time1)
                time2=time2+(elapse_time4-elapse_time3)
	enddo loop1
        write(60,*)time1,time1-time2,time2
!MT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!

    

!MT
    dos(:,:) = 0.0d0
    nj (:)   = 0.0d0
    do i=1,n
       do ind=1,m_int
          dos(ind,i)=-dimag(sx(interaction_list(ind),i))/pi_math
       enddo      
    enddo
!    call calc_nj(j_src,dos,nj)
    do i=1,n
       vect_in(1:m_int)=dos(1:m_int,i)
       call inner_product_with_s_proj(vect_in,nj(i),j_src)
       write(80+j_src,*)dreal(sigma(i)+dreal(z_seed)),nj(i)
    enddo
   
    write(80,*)j_src,sum(nj(1:n))*0.001d0 

!MT   
!
!    do ind=1,m_int
!      i=interaction_list(ind)
!      write(*,*)'ind, i=',ind,i
!      if ((i < 1) .or. (i > m)) then
!        write(*,*)'Abort(gscocg_main)'
!        stop
!      endif   
!    enddo
!
!    stop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate (interaction_list)
    deallocate (x,pi)
    deallocate (r, p, q, pi_v, qs, rk)
    deallocate (sx,sp,sigma,salpha,sbeta)
    deallocate (pi_n,pi_o)
    deallocate (hg2)
    deallocate (flg,iter)
    deallocate (sax,sbx)
! Final deallocation of matrices
!
  end subroutine gscocg_main
!
end module M_la_gscocg_main





