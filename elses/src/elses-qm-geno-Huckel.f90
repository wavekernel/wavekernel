!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_Huckel
  ! generalized non orthogoanal calculation, generalized Huckel
  use M_qm_domain, only: i_verbose, DOUBLE_PRECISION,&
       atm_element, e_num_on_atom
  use M_qm_geno_Huckel_atom_params, &
       only: AtomParameter, NUMBER_OF_ORBITALS, NUMBER_OF_CONF_VOIP, &
       UseVOIP, vanDerWaals, genoMethod, genoSTD, genoNOLIN, genoICON, &
       normDoubleZeta, sdelta=>small_delta, kappa, &
       constantK, Kval
  use elses_mod_phys_const, only : eV4au
  implicit none

  private
  public :: SetOverlap, SetDiagonalElements, SetNondiagonalElements
  public :: SetDiagonalShifts
  public :: RepulsiveEnergy, calc_SS_dSS

  integer, parameter :: MAXN=7, MAXL=2
  integer, allocatable :: cmb(:,:)
  real(DOUBLE_PRECISION), allocatable :: fct(:), fct2(:,:)
  !real(DOUBLE_PRECISION),allocatable :: cff(:,:,:,:,:,:)
  logical :: initFinished=.false.

contains

 subroutine initOVLCalc()
   implicit none
   integer :: j,k,p,q,n, ll, mm
!$OMP CRITICAL
   if (.not. initFinished) then ! this IF statement is necessary not to doubly allocate arrays
      allocate(cmb(0:MAXN,0:MAXN)) ! 2*MAXL < MAXN (MAXL=2, MAXN=7)
      allocate(fct(0:2*MAXN),fct2(0:2*MAXN,0:2*MAXN)) ! fct2(j,k)=fct(j)/fct(k) (j>=k)
!!$      allocate(cff(MAXN,0:MAXL,0:MAXL,0:MAXL,0:MAXN, 0:MAXL))
      !  (principle q.n,  L,     M,      p,     q,      k)
      cmb=0
      !cff=0d0
      do p=0, MAXN
         cmb(p,0)=1
         cmb(p,p)=1
         do q=1, p-1
            cmb(p,q)=cmb(p-1,q-1)+cmb(p-1,q)
         end do
      end do
      fct2=0d0
      fct2(0,0)=1d0
      do p=1, 2*MAXN
         fct2(p,p)=1d0
         do q=p-1, 0, -1
            fct2(p,q)=fct2(p,q+1)*(q+1)
         end do
      end do
      fct(:)=fct2(:,0)
!!$      do n=1, MAXN
!!$         do ll=0, min(n-1,2)
!!$            do mm=0, ll
!!$               do j=0, (ll-mm)/2
!!$                  do p=0, ll-mm-2*j
!!$                     do q=0, n-ll+2*j
!!$                        do k=0, mm
!!$                           cff(n,ll,mm,p,q,k)=&
!!$                                fct(ll+mm)*cmb(ll,j)*cmb(2*(ll-j),ll+mm)*&
!!$                                cmb(ll-mm-2*j,p)*cmb(n-ll+2*j,q)*cmb(mm,k)
!!$                        end do
!!$                     end do
!!$                  end do
!!$               end do
!!$            end do
!!$         end do
!!$      end do
      initFinished=.true.
   end if
!$OMP END CRITICAL
 end subroutine initOVLCalc

 subroutine SetOverlap(left_atom_spec,right_atom_spec,x,y,z,S,dS,bondMask)
  ! only for the non_diagonal part
  ! coded by m.ikeda, fla, Japan, 2008.09.05
  ! Ver-2008-09-03-ver1
  ! Ver-2008-09-11-ver2
  !   do loop maximum value from nval===> num_val_orb.
  !   maximum of num_val_orb nval_a ==> nval_orb_a
  !   maximum of num_val_orb nval_b ==> nval_orb_b
  !   added -verbose option
  !   overwrite B-function when lamda=0.
  ! Ver-2008-09-13-ver3
  !   correct s,p-d matrix element
  !   move inside endif
  !   correct Eq.(1.26) (2*l_a+1)!==>(2*l_a+1) in factor3 expression.
  ! ver-2008-09-18-ver4
  !   correct Eq.(1.56) upper half element {s,p} - {d}/* => +
  ! Ver-2008-09-20-ver5
  !   correct B(n,lamda) ==> B(n,-lamda)
  !           correct expressions according to the correct Eq.(1.26).
  !   correct TSKbond vector expression using Slater-Koster table.
  !           S(5,2), S(5,3) and S(5,4).
  ! Ver-2008-09-23-ver6
  !   move summation of two derivatives outside select statments.
  ! Ver-2008-09-27-ver7
  !   store TSK bond vector into Bond_matrix form.
  !   both for overlap and differential of overlap matrix elements.
  ! Ver-2008-09-27-ver8
  !   add TSKBond Bond_matrix and dBond_matrix
  ! Ver-2008-09-28-ver9
  !   move select statements into subroutine Rotation_Slater_Koster
  ! Ver-2008-10-01
  !   add print out S and dS matrix in formatted form.
  !       orbita_name(9)
  ! Ver-2008-10-02
  !   correct nval_orb_b ==>nval_orb_a in print out mode.
  !   improve the print format especially for the case which includes d-orbitals.
  ! Ver-2008-10-03
  !   correct s,p-d cross term: (npq_a+npq_b)==>(npq_a+npq_b+1)
  !   added dSS11=0.d0 satatement in s,p-d cross term in lower half elements.
  ! Ver-2008-10-16
  !   correct sign_factor: unnecessary factor
  !   remove sign_factor from Rotate_Slater_Koster routine
  !   sign_factor=1.do for all cases.
  ! Ver-2008-11-05
  !   sign_factor: necessary for calling skderive routine.
  !   correct : Rotation_Slater_Koster routine
  !   correct : summation of derivatives
  ! Ver-2008-11-25
  !   correct : Rotation_Slater_Koster routine
  !             In Vyy, Dyy, Vyy, sign_factor=-1 is irrelevant.
  ! Ver-2008-12-01
  !   correct : Vsd_sigma, dVsd_sigma, Vds_sigma, dVds_sigma
  !             were not defined for nval_orb=1 and  nvasl_orb=9.
  !             The correct treatment is included.
  ! Jul. 6, 2010 : Drastically rewritten by S.Y
    use MSlaterKosterTable
    use MSlaterKosterDerivative
    implicit none
    type(TSKBond)            ::  bond
    type(TSKBond)            :: dBond
    integer, intent(in)      :: left_atom_spec, right_atom_spec
    real(DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(DOUBLE_PRECISION), intent(out):: S(:,:), dS(:,:,:)
    real(DOUBLE_PRECISION)   :: R
    real(DOUBLE_PRECISION)   :: rl, rm, rn
    real(DOUBLE_PRECISION), dimension(0:2)  :: v_SS, v_dSS
    real(DOUBLE_PRECISION), dimension(0:2,2):: v2tmp, v2dtmp
    real(DOUBLE_PRECISION)   :: zeta_a(2), zeta_b(2), c_a(2), c_b(2), normalization
    integer                  :: nval_orb_a, nval_orb_b, max_nval, max_l_a, max_l_b, nzeta_a, nzeta_b
    integer                  :: n, II, JJ
    integer                  :: l_a, l_b, m_a, m_b
    integer                  :: npq_a, npq_b
    real(DOUBLE_PRECISION),optional         :: bondMask(:)

    if(.not. initFinished) call initOVLCalc() ! init combinatorial, factorial, etc

    nval_orb_b=AtomParameter( left_atom_spec)%num_val_orb
    nval_orb_a=AtomParameter(right_atom_spec)%num_val_orb
    max_nval = max(nval_orb_a, nval_orb_b)
    max_l_b  = nint(sqrt(real(nval_orb_b,DOUBLE_PRECISION)))-1
    max_l_a  = nint(sqrt(real(nval_orb_a,DOUBLE_PRECISION)))-1

    !Eq.(1.26)

    R = sqrt(x*x + y*y + z*z)
    rl= x/R
    rm= y/R
    rn= z/R

     bond%v(4)=0d0
    dBond%v(4)=0d0
    do l_b=0, max_l_b
       npq_b        = AtomParameter( left_atom_spec)%principal(l_b*l_b+1)
       c_b(2)       = AtomParameter( left_atom_spec)%c2(l_b*l_b+1)
       if( c_b(2) == 0d0 ) then
          nzeta_b   = 1
          c_b(1)    = 1d0
          zeta_b(1) = AtomParameter( left_atom_spec)%zeta(l_b*l_b+1)
       else
          nzeta_b   = 2
          c_b(1)    = AtomParameter( left_atom_spec)%c1(l_b*l_b+1)
          zeta_b(1)= AtomParameter( left_atom_spec)%zeta(l_b*l_b+1)
          zeta_b(2)= AtomParameter( left_atom_spec)%zeta2(l_b*l_b+1)
          normalization=normDoubleZeta(npq_b,zeta_b(1),zeta_b(2),c_b(1),c_b(2))
          if(abs(1-normalization) > 1d-6) c_b=c_b/normalization
       end if
       
       do l_a=0, max_l_a       ! The next line is an optimization for homogeneous case
          if(left_atom_spec==right_atom_spec .and. l_b > l_a ) cycle 
          
          npq_a         = AtomParameter(right_atom_spec)%principal(l_a*l_a+1)
          c_a(2)        = AtomParameter(right_atom_spec)%c2(l_a*l_a+1)
          if( c_a(2) == 0d0 ) then
             nzeta_a    = 1
             c_a(1)     = 1d0
             zeta_a(1)  = AtomParameter(right_atom_spec)%zeta(l_a*l_a+1)
          else
             nzeta_a    = 2
             c_a(1)     = AtomParameter(right_atom_spec)%c1(l_a*l_a+1)
             zeta_a(1)  = AtomParameter(right_atom_spec)%zeta(l_a*l_a+1)
             zeta_a(2)  = AtomParameter(right_atom_spec)%zeta2(l_a*l_a+1)
             normalization=normDoubleZeta(npq_a,zeta_a(1),zeta_a(2),c_a(1),c_a(2))
             if(abs(1-normalization) > 1d-6) c_a=c_a/normalization
          end if
          if(nzeta_a==1 .and. nzeta_b==1) then ! No double-zeta
             if(l_b<=l_a)then
                call calc_SS_dSS(npq_b,npq_a,l_b,l_a,zeta_b(1),zeta_a(1),R,v_SS,v_dSS)
             else
                call calc_SS_dSS(npq_a,npq_b,l_a,l_b,zeta_a(1),zeta_b(1),R,v_SS,v_dSS)
             end if
          else
              v2tmp=0d0
             v2dtmp=0d0
             do ii=1, nzeta_b
                do jj=1, nzeta_a
                   if(l_b<=l_a)then ! In l_b > l_a case,
                      ! calc_SS_dSS does not satisfy the requirement of
                      ! SlaterKosterCoefficient & SlaterKosterDerivative
                      call calc_SS_dSS(npq_b,npq_a,l_b,l_a,zeta_b(ii),zeta_a(jj),R,v_SS,v_dSS)
                   else
                      call calc_SS_dSS(npq_a,npq_b,l_a,l_b,zeta_a(jj),zeta_b(ii),R,v_SS,v_dSS)
                   end if
                    v2tmp(:,ii)= v2tmp(:,ii)+  v_SS(:)*c_a(jj)
                   v2dtmp(:,ii)=v2dtmp(:,ii)+ v_dSS(:)*c_a(jj)
                end do
             end do
              v_SS=0d0
             v_dSS=0d0
             do ii=1, nzeta_b
                 v_SS(:)= v_SS(:)+ v2tmp(:,ii)*c_b(ii)
                v_dSS(:)=v_dSS(:)+v2dtmp(:,ii)*c_b(ii)
             end do
          end if
           bond%v(1:3)= v_SS(0:2) ! v_SS(0:2) sigma, pi, delta ...
          dBond%v(1:3)=v_dSS(0:2) ! derivative (w.r.t r) of v_SS
          ! If l_b = 1, l_a = 2, then v_SS(0:1) = (/ pd\sigma, pd\pi /)
          if(l_a == l_b .and. present(bondMask))then
             ! For COHP analysis, p-p\pi component is set to be 0.
             ! We assume that each element of bondMask(:) is corresponding to
             ! ss-\sigma,n.a,pp-\sigma,pp-\pi,n.a,n.a,dd-\sigma,dd-\pi,dd-\delta
             ! n.a : not available
             do ii=0, l_a
                 bond%v(ii+1)= bond%v(ii+1)*bondMask(l_a*(l_a+1)+ii+1)
                dBond%v(ii+1)=dBond%v(ii+1)*bondMask(l_a*(l_a+1)+ii+1)
             end do
          end if
          
          do ii=l_b*l_b, l_b*(l_b+2)   ! s px py pz d_xy d_yz d_zx d_xx d_zz
             do jj=l_a*l_a, l_a*(l_a+2)! s px py pz d_xy d_yz d_zx d_xx d_zz
                if(l_b <= l_a)then ! In
                   ! calc_SS_dSS does not satisfy the requirement of
                   ! SlaterKosterCoefficient & SlaterKosterDerivative

                   S(ii+1,jj+1)=&
                        SlaterKosterCoefficient( bond,ii,jj,rl,rm,rn)
                   call SlaterKosterDerivative&
                        ( bond,ii,jj,rl,rm,rn,dS(:,ii+1,jj+1))
                   dS(:,ii+1,jj+1)= dS(:,ii+1,jj+1)/R+&
                        SlaterKosterCoefficient(dBond,ii,jj,rl,rm,rn)&
                        *(/rl,rm,rn/)
                else
                   ! see the correspondence between S.Y & M.I on Oct. 28, 2008.
                   S(ii+1,jj+1)=&
                        SlaterKosterCoefficient( bond,jj,ii,-rl,-rm,-rn)
                   call SlaterKosterDerivative&
                        ( bond,jj,ii,-rl,-rm,-rn,dS(:,ii+1,jj+1))
                   dS(:,ii+1,jj+1)=-dS(:,ii+1,jj+1)/R+&
                        SlaterKosterCoefficient(dBond,jj,ii,-rl,-rm,-rn)&
                        *(/rl,rm,rn/)
                end if
             end do
          end do
       end do
       if(left_atom_spec==right_atom_spec) then ! Optimization for homogeneous case
         do l_a=0, l_b-1
             do ii=l_b*l_b,l_b*(l_b+2)  !The order of orbitals follows Iguchi's notation
                do jj=l_a*l_a, l_a*(l_a+2) ! The situation is same as above line.
                    S(  ii+1,jj+1)=  S(  jj+1,ii+1)*(1-2*mod(l_a+l_b,2))
                   dS(:,ii+1,jj+1)= dS(:,jj+1,ii+1)*(1-2*mod(l_a+l_b,2))
                end do
             end do
          end do
       end if
    end do
  end subroutine SetOverlap

  subroutine calc_SS_dSS(npq_b,npq_a,l_b,l_a,zeta_b,zeta_a,R,SS,dSS)
    ! This routine is originally coded by M.Ikeda at FLA
    ! Then, SY modified.
    !----------------------------------------------------------
    ! Coding according to the document provided by Dr.Yamamoto.
    ! Eq.(1.26)
    ! added by m.ikeda, 2008/08/04, FLA, Japan.
    ! Store angular momentum 
    ! Orbital
    ! s       L=0 ===> num_val_orb(=nval)=1
    ! s,p     L=1 ===> num_val_orb(=nval)=4
    ! s,p,d   L=2 ===> num_val_orb(=nval)=9
    ! Jul. 9, 2010 : Drastically rewritten by S.Y
    implicit none
    integer,                intent(in) :: npq_a, npq_b, l_a, l_b
    real(DOUBLE_PRECISION), intent(in) :: zeta_a, zeta_b, R
    real(DOUBLE_PRECISION) :: SS(0:), dSS(0:)
    real(DOUBLE_PRECISION) :: R_h
    real(DOUBLE_PRECISION) :: factor, dfactor, facta, factb, tmp, suma, sumb, sumda, sumdb
    real(DOUBLE_PRECISION) :: AA(0:2*MAXN), BB(0:2*MAXN), dAA(0:2*MAXN), dBB(0:2*MAXN)
    integer :: j_a, ip_a, iq_a, iaa, ibb
    integer :: j_b, ip_b, iq_b
    integer :: k, k_prime, ii, mm

    R_h= R/2.d0
    
    do ii=0, npq_a+npq_b
       AA(ii)=A(ii, R_h*(zeta_a+zeta_b), dAA(ii))
       BB(ii)=B(ii, R_h*(zeta_b-zeta_a), dBB(ii))
    end do
    
     SS=0d0
    dSS=0d0
    do mm=0, min(l_a,l_b)
       ! Eq.(1.26)
       ! SS (0:2) : Overlap matrix representation of Eq.(1.26)
       ! dSS(3,0:2) : derivative of SS w.r.t position of right_atom.
       ! npq : Principal quantum numbers
       ! zeta: Screening parameters
       ! fct(n) = k! , fct2(n,m) = n! / m!, cmb(n,m)= n!/(m! * (n-m)!)
       factor&
            =(2.d0*zeta_a)**(2*npq_a+1)/fct(2*npq_a)&
            *(2.d0*zeta_b)**(2*npq_b+1)/fct(2*npq_b)
       factor=factor*real(2*l_a+1,DOUBLE_PRECISION)/2d0/fct2(l_a+mm,l_a-mm)&
                    *real(2*l_b+1,DOUBLE_PRECISION)/2d0/fct2(l_b+mm,l_b-mm)
       factor=sqrt(factor)*fct2(l_a+mm,l_a)*fct2(l_b+mm,l_b)/2**(l_a+l_b)
       dfactor= 0.5d0*(npq_a+npq_b+1)*(R_h)**(npq_a+npq_b)  *factor
        factor=                       (R_h)**(npq_a+npq_b+1)*factor
       do j_a=0, (l_a-mm)/2
          do ip_a=0, l_a-mm-2*j_a
             do iq_a=0, npq_a-l_a+2*j_a
                facta=cmb(l_a,j_a)*cmb(2*(l_a-j_a),l_a+mm)&
                     *cmb(l_a-mm-2*j_a,ip_a)*cmb(npq_a-l_a+2*j_a,iq_a)
                do j_b=0, (l_b-mm)/2
                   do ip_b=0, l_b-mm-2*j_b
                      do iq_b=0, npq_b-l_b+2*j_b
                         factb=cmb(l_b,j_b)*cmb(2*(l_b-j_b),l_b+mm)&
                              *cmb(l_b-mm-2*j_b,ip_b)*cmb(npq_b-l_b+2*j_b,iq_b)
                         iaa=npq_a+npq_b-ip_a-ip_b-iq_a-iq_b
                         ibb=l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b)
                         suma =0d0
                         sumda=0d0
                         do k=0, mm
                            suma = suma+cmb(mm,k)*(1-2*mod(k,2))* AA(iaa-2*k)
                            sumda=sumda+cmb(mm,k)*(1-2*mod(k,2))*dAA(iaa-2*k)
                         end do
                         sumb =0d0
                         sumdb=0d0
                         do k_prime=0, mm
                            sumb = sumb+cmb(mm,k_prime)*(1-2*mod(k_prime,2))* BB(ibb-2*k_prime)
                            sumdb=sumdb+cmb(mm,k_prime)*(1-2*mod(k_prime,2))*dBB(ibb-2*k_prime)
                         end do
                         tmp   = facta*factb*(1-2*mod(ip_a+iq_a+j_a+j_b+mm,2))
                          SS(mm)= SS(mm) + tmp *  suma * sumb
                         dSS(mm)=dSS(mm) + tmp *( sumda* sumb * 0.5d0 * (zeta_a+zeta_b)&
                              +                   suma * sumdb* 0.5d0 * (zeta_b-zeta_a))
                      enddo ! iq_b
                   enddo ! ip_b
                enddo ! j_b
             enddo ! iq_q
          enddo ! ip_a
       enddo ! j_a
       dSS(mm) = dfactor * SS(mm) + factor * dSS(mm)
        SS(mm) =  factor * SS(mm)
    end do
  contains
    function A(n,lambda,dA) !(also returns dA) for positive lambda
      real(DOUBLE_PRECISION)              :: A
      integer, intent(in)       :: n
      real(DOUBLE_PRECISION), intent(in)  :: lambda
      real(DOUBLE_PRECISION), intent(out) :: dA
      integer                   :: k
      real(DOUBLE_PRECISION)              :: e
!!$      if(lambda<0) stop ' lambda < 0'
      !Eq.(1.37), Eq(1.39)
      e=exp(-lambda)/lambda
      A=e
      dA=-e*(1d0+1d0/lambda)
      do k=1,n
         e=e*real(n-k+1,kind=DOUBLE_PRECISION)/lambda
         A=A+e
         dA=dA-e*(1d0+real(k+1,kind=DOUBLE_PRECISION)/lambda)
      enddo
    end function A
!!$    function A_N(n,lambda,dA) !(also returns dA) for negative lambda
!!$      real(DOUBLE_PRECISION)              :: A_N
!!$      integer, intent(in)       :: n
!!$      real(DOUBLE_PRECISION), intent(in)  :: lambda
!!$      real(DOUBLE_PRECISION), intent(out) :: dA
!!$      integer                   :: k
!!$      real(DOUBLE_PRECISION)              :: e
!!$      !if(lambda>0) stop ' lambda < 0'
!!$      A_N=1d0/lambda
!!$      do k=1, n
!!$         A_N=A_N*real(k,  kind=DOUBLE_PRECISION)/lambda
!!$      end do
!!$      dA   =-A_N*real(n+1,kind=DOUBLE_PRECISION)/lambda
!!$      e=1d0
!!$      k=0
!!$      do while ( abs(e) > 1d-16 )
!!$         A_N= A_N-e/real(k+n+1,kind=DOUBLE_PRECISION)
!!$         dA=dA+e/real(k+n+2,kind=DOUBLE_PRECISION)
!!$         k=k+1
!!$         e=e*(-lambda)/real(k,kind=DOUBLE_PRECISION)
!!$      end do
!!$    end function A_N
    function B(n,lambda,dB) ! (also returns dB)
      real(kind=DOUBLE_PRECISION):: B,dB,lambda,f,g
      integer     :: n,m,k
      intent(in)  :: n,lambda
      intent(out) :: dB
      real(DOUBLE_PRECISION), parameter :: CRITERION=5d0
      ! The more  CRITERION increase and the more precision increase.
      ! When CRITERION=0.5d0, it is comparable to single precision.

      ! Eq.(1.38), Eq(1.40)
      k=iand(n,1)
     if(abs(lambda) < CRITERION )then
         if(k == 0)then
            m = 0
            f = 2d0  ! 2*lamda^(2*m)/(2*m)!
            B = f/real(n+1,kind=DOUBLE_PRECISION)
            f = f*lambda/real(m+1,kind=DOUBLE_PRECISION)
            ! 2*lamda^(2*m+1)/(2*m+1)!
            dB= f/real(n+3,kind=DOUBLE_PRECISION)
            m = m+2
            do while( abs(f)>1d-16 )
               f = f*lambda/real(m,kind=DOUBLE_PRECISION)
               B = B+f/real(n+m+1,kind=DOUBLE_PRECISION)
               f = f*lambda/real(m+1,kind=DOUBLE_PRECISION)
               dB=dB+f/real(n+m+3,kind=DOUBLE_PRECISION)
               m = m+2
            end do
         else
            m = 0
            f =-2d0  !-2*lamda^(2*m)/(2*m)!
            dB= f/real(n+2,kind=DOUBLE_PRECISION)
            f = f*lambda/real(m+1,kind=DOUBLE_PRECISION)
            ! 2*lamda^(2*m+1)/(2*m+1)!
            B = f/real(n+2,kind=DOUBLE_PRECISION)
            m = m+2
            do while( abs(f)>1d-16 )
               f = f*lambda/real(m,kind=DOUBLE_PRECISION)
               dB=dB+f/real(n+m+2,kind=DOUBLE_PRECISION)
               f = f*lambda/real(m+1,kind=DOUBLE_PRECISION)
               B =B +f/real(n+m+2,kind=DOUBLE_PRECISION)
               m = m+2
            end do
         end if
      else

            B=-(1-2*k)*A(n,-lambda,dB)
            dB=(1-2*k)*dB
            f=A(n,lambda,g)
            B=B-f
            dB=dB-g
         
!!$         if (lambda > 0) then
!!$            B=-(1-2*k)*A_N(n,-lambda,dB)
!!$            dB=(1-2*k)*dB
!!$            f=A(n,lambda,g)
!!$            B=B-f
!!$            dB=dB-g
!!$         else
!!$            B=-(1-2*k)*A(n,-lambda,dB)
!!$            dB=(1-2*k)*dB
!!$            f=A_N(n,lambda,g)
!!$            B=B-f
!!$            dB=dB-g
!!$         end if
      end if
    end function B
  end subroutine calc_SS_dSS

  subroutine SetDiagonalElements(atom_spec,Occ,Ham)
    use elses_mod_phys_const
    integer, intent(in)      :: atom_spec
    real(DOUBLE_PRECISION), intent(in) :: Occ(:)
    real(DOUBLE_PRECISION), intent(out):: Ham(:,:)
    ! VOIP(inital_charge)
    ! If d-orbital exists, use VOIP_D & initial_occupation instead.
    integer :: i, j, num_val_orb
    real(DOUBLE_PRECISION) :: q, VOIP_spd(3,9), D_sigma, D_pi, D_delta
    q = AtomParameter(atom_spec)%num_val_elec-sum(Occ)+AtomParameter(atom_spec)%initial_charge
    num_val_orb = AtomParameter(atom_spec)%num_val_orb

         D_sigma = Occ(1)
         D_pi    = 0d0
         D_delta = 0d0
    if(num_val_orb .gt.1)&
         D_pi    = sum(Occ(2:4))
    if(num_val_orb .gt.4)&
         D_delta = sum(Occ(5:9))

    select case (num_val_orb)
    case(1,4)
       ! s-orbital or s- and p-orbital
       do i = 1, num_val_orb
          if (.not. UseVOIP) then
             Ham(i,i) = AtomParameter(atom_spec)%initial_diagonal_elements(i)
          else
             Ham(i,i) = -1.0d0 * ( AtomParameter(atom_spec)%VOIP(3,i) &
                                 + AtomParameter(atom_spec)%VOIP(2,i) * q &
                                 + AtomParameter(atom_spec)%VOIP(1,i) * q * q )
          endif
       end do

    case(9)
       ! s-, p- and d-orbital
       if(.not. UseVOIP) then
          do i = 1, num_val_orb
             Ham(i,i) = AtomParameter(atom_spec)%initial_diagonal_elements(i)
          enddo

       else
          do j = 1, num_val_orb
             do i = 1, NUMBER_OF_CONF_VOIP
                VOIP_spd(i, j) = AtomParameter(atom_spec)%VOIP_D(3,i,j) &
                               + AtomParameter(atom_spec)%VOIP_D(2,i,j) * q &
                               + AtomParameter(atom_spec)%VOIP_D(1,i,j) * q * q
             end do
          end do

          do i = 1, num_val_orb
             select case (i)
             case (1)
                Ham(i,i) = -1.0d0 * ( (2.0d0 - D_sigma - D_pi) * VOIP_spd(1, i) &
                         + (D_sigma - 1.0d0) * VOIP_spd(2, i) &
                         + D_pi * VOIP_spd(3, i) )
             case(2:4)
                Ham(i,i) = -1.0d0 * ( (2.0d0 - D_sigma - D_pi) * VOIP_spd(1, i) &
                         + (D_pi - 1.0d0) * VOIP_spd(2, i) &
                         + D_sigma * VOIP_spd(3, i) )
             case(5:9)
                Ham(i,i) = -1.0d0 * ( (1.0d0 - D_sigma - D_pi) * VOIP_spd(1, i) &
                         + D_sigma * VOIP_spd(2, i) &
                         + D_pi * VOIP_spd(3, i) )
             case default
                stop 'SetDiagonalElements: !!!Error!!! invalid orbital in VOIP'
             end select
          end do
       endif

    case default
        write(*,'("SetDiagonalElements: !!!ERROR!!! num_val_orb = ",I8)') num_val_orb
        stop 'Number of valence electron orbitals is not 1, 4 or 9!!'

    end select
!
    if(i_verbose >= 200) then
       write(*,'("SetDiagonalElements:")') 
       write(*,'("  Initial Num. Val. Elec. = ",ES21.14,", Num. Val. Elec. = ",ES21.14,&
            &"Initial Q= ",ES21.14,", Q = ",ES21.14)') &
            AtomParameter(atom_spec)%num_val_elec, sum(Occ), AtomParameter(atom_spec)%initial_charge, q
       if(UseVOIP) write(*,'("  D_sigma = ",ES21.14,", D_pi = ",ES21.14,", D_delta = ",ES21.14)') &
            D_sigma, D_pi, D_delta
       do i= 1, num_val_orb
          write(*,'(A,4(I1,A),ES13.6,A)') &
               '  H(', i, ',', i, ',', atom_spec,',', atom_spec, ') = ', Ham(i,i)*eV4au, ' eV'
       enddo
    endif
!
  end subroutine SetDiagonalElements

  subroutine SetDiagonalShifts(atom_spec,Ham)
    use elses_mod_phys_const
    integer, intent(in)      :: atom_spec
    real(DOUBLE_PRECISION), intent(inout):: Ham(:,:)
    integer :: i, num_val_orb
    num_val_orb = AtomParameter(atom_spec)%num_val_orb
    do i=1, num_val_orb
       Ham(i,i)=Ham(i,i)+AtomParameter(atom_spec)%extra_shifts(i)
    end do
    if(i_verbose >= 200) then
       write(*,'("SetDiagonalShifts:")') 
       do i= 1, num_val_orb
          write(*,'(A,4(I1,A),ES13.6,A)') &
               '  H(', i, ',', i, ',', atom_spec,',', atom_spec, ') = ', Ham(i,i)*eV4au, ' eV'
       enddo
    endif
!
  end subroutine SetDiagonalShifts

  subroutine SetNondiagonalElements(left_atom_spec, right_atom_spec,&
       x,y,z, S, left_H_diagonal, right_H_diagonal, Ham, dS, dHam)
    use M_qm_domain
    integer, intent(in)      :: left_atom_spec,right_atom_spec
    real(DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(DOUBLE_PRECISION), intent(in) :: S(:,:), dS(:,:,:)
    real(DOUBLE_PRECISION), intent(in) :: left_H_diagonal(:,:), right_H_diagonal(:,:)
    real(DOUBLE_PRECISION), intent(out):: Ham(:,:), dHam(:,:,:)
    integer :: nval1,nval2,ja1,ja2,w
    real(DOUBLE_PRECISION) :: r(3),lhd,rhd
    integer :: ja(2),atom_s(2)
    real(DOUBLE_PRECISION) :: rn(2),r0,delta
    real(DOUBLE_PRECISION) :: ff,ff2,cl1,cl2,zl1,zl2,hi1,hi2
    integer :: nl,i
    r(1:3)=(/ x,y,z /)
    atom_s(1:2)=(/ left_atom_spec,right_atom_spec /)
    nval1=AtomParameter(atom_s(1))%num_val_orb
    nval2=AtomParameter(atom_s(2))%num_val_orb

    do ja2=1,nval2
      do ja1=1,nval1
        ja(1:2)=(/ ja1,ja2 /)

        ! Eq.(2.2)
        hi1= left_H_diagonal(ja(1),ja(1))
        hi2=right_H_diagonal(ja(2),ja(2))

        if(abs(hi1+hi2) < 1e-5) then
           write(*,'("!!!!! ERROR !!!!!")')
           write(*,'("This is the case of Hii=-Hjj.")')
           write(*,'("The quantity Delta can not be evaluated.")')
           stop
        end if

        delta=(hi1-hi2)/(hi1+hi2)

        ! Correction of delta
        if(abs(delta)>1.0d0) then
           if( i_verbose >= 100 )then
              write(*,'("!!!!! WARNING !!!!!")')
              write(*,'("The value of Delta is corrected.")')
              write(*,'(ES14.5,"->",ES14.5)') delta, 1d0
              write(*,'("left atom=",I8,", right atom=",I8)') ja(:)
           end if
           delta=1.0d0
        end if

        ! Eq.(2.3)
        do i=1,2
          nl=AtomParameter(atom_s(i))%principal(ja(i))
          zl1=AtomParameter(atom_s(i))%zeta(ja(i))
          if(ja(i)<=4)then
            rn(i)=nl / zl1
          else
            cl1=AtomParameter(atom_s(i))%c1(ja(i))
            cl2=AtomParameter(atom_s(i))%c2(ja(i))
            zl2=AtomParameter(atom_s(i))%zeta2(ja(i))
!            ff=1.0D0
            ff=cl1*cl1 + cl2*cl2 + 2.0D0*cl1*cl2*&
              ((4.0D0*zl1*zl2/(zl1+zl2)/(zl1+zl2))**(nl+0.5D0))
            select case(genoMethod)
            case(genoICON)
               !JPC 1989 (93) 5366-5371 : but normalization is wrong ...
               ff2=cl1*cl1*zl1 + cl2*cl2*zl2 + 0.5D0*&
                    ((4.0D0*zl1*zl2)**(nl+0.5D0))/((zl1+zl2)**(2.0D0*nl))
            case default
               ! This is a standard form
               ff2=cl1*cl1*zl1 + cl2*cl2*zl2 + cl1*cl2*&
                    ((4.0D0*zl1*zl2)**(nl+0.5D0))/((zl1+zl2)**(2.0D0*nl))
            end select
            rn(i)=nl*ff/ff2
          endif
        enddo
        r0=rn(1)+rn(2)

        ! Eq.(2.7), dH_ij
        lhd= left_H_diagonal(ja1,ja1)
        rhd=right_H_diagonal(ja2,ja2)
        if(constantK) then
           Ham(ja1,ja2)=(lhd+rhd)/2.0D0 * Kval * S(ja1,ja2)
           do w=1,3
              dHam(w,ja1,ja2)=(lhd+rhd)/2.0D0 * Kval * dS(w,ja1,ja2)
           enddo
        else
           Ham(ja1,ja2)=K(r,r0,delta)*(lhd+rhd)/2.0D0*S(ja1,ja2)
           do w=1,3
              dHam(w,ja1,ja2)=(lhd+rhd)/2.0D0 * &
                   (dK(w,r,r0,delta)*S(ja1,ja2)+K(r,r0,delta)*dS(w,ja1,ja2))
           enddo
        end if

        if (i_verbose >= 50) then
           write(*,*)'  ja1,ja2     =',ja1,ja2
           write(*,*)'        K     =',K(r,r0,delta)
           write(*,*)'      H_ii[eV]=',lhd*ev4au
           write(*,*)'      H_jj[eV]=',rhd*ev4au
           write(*,*)'(H_ii+H_jj)/2 =',(lhd+rhd)/2.0D0*ev4au
           write(*,*)'      H_ij[eV]=',Ham(ja1,ja2)*ev4au
        endif        
        if(i_verbose >= 200) then
          write(*,'(a,4(i1,a),e13.6)') &
            'H(',ja1,',',ja2,',',left_atom_spec,',',right_atom_spec,&
            ')=',Ham(ja1,ja2)
        endif
      enddo
    enddo
  contains

    function K(r,r0,delta) result(kk)
      real(DOUBLE_PRECISION),intent(in) :: r(3),r0,delta
      real(DOUBLE_PRECISION) :: rr,q,Kp,kk
      rr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) !a.u.

      ! Eq.(2.4) - Eq.(2.6)
      q=((rr-r0)-dabs(rr-r0))*sdelta
      q=1.0D0+q*q
      Kp=kappa + delta*delta - delta*delta*delta*delta*kappa
      kk=1.0D0 + Kp*exp(-1.0D0*sdelta*(rr-r0))/q

    end function K

    function dK(j,r,r0,delta) result(dkk)
      integer,intent(in) :: j
      real(DOUBLE_PRECISION),intent(in) :: r(3),r0,delta
      real(DOUBLE_PRECISION) :: rr,q,Kp,dkk
      rr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) !a.u.

      ! Eq.(2.5) 
      Kp=kappa + delta*delta - delta*delta*delta*delta*kappa

      if (rr.ge.r0) then
        ! Eq.(2.4), Eq.(2.6)
        q=1.0D0
        dkk= -1.0D0*sdelta*Kp*exp(-1.0D0*sdelta*(rr-r0))*r(j)/rr
      else
        ! Eq.(2.4), Eq.(2.6)
        q=2.0D0*(rr-r0)*sdelta
        q=1.0D0+q*q
        dkk= -1.0D0*sdelta*Kp*exp(-1.0D0*sdelta*(rr-r0))*r(j)/rr*&
          (q+8.0D0*sdelta*(rr-r0))/q/q
      endif
    end function dK
  end subroutine SetNondiagonalElements

  function RepulsiveEnergy(left_atom, right_atom,&
       x,y,z, dRepulsiveEnergy, E_atompair_vdW_only)

    use M_config             ! unchanged
    implicit none
    real(kind=DOUBLE_PRECISION) :: RepulsiveEnergy
    integer, intent(in)         :: left_atom,right_atom
    integer                     :: left_atom_spec,right_atom_spec
    real(kind=DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(kind=DOUBLE_PRECISION), intent(out):: dRepulsiveEnergy(3)
    real(DOUBLE_PRECISION),  optional :: E_atompair_vdW_only ! Atom-pair Energy only for vdW interaction
!
    real(kind=DOUBLE_PRECISION) :: Za, Zb, R, zeta, s, t, occ, d, dt, ds
    real(kind=DOUBLE_PRECISION) :: z1, z2, w1, w2, d1, d2, zs
    real(kind=DOUBLE_PRECISION) :: EvanDerWaals, dEvanDerWaals, Ia, Ib, aa, ab
    ! EvanDerWaals is calculated by London dispersion formula -(3/2) ((Ia Ib)/(Ia+Ib)) aa ab /R^6
    real(kind=DOUBLE_PRECISION) :: fd, dfd, vdW_R_a, vdW_R_b, xab
    real(kind=DOUBLE_PRECISION) :: vdW_lambda
    real(kind=DOUBLE_PRECISION) :: vdW_lambda_default = 7.5d-4
    integer :: vdW_n = 8
    ! S.Y Sep 16, 2013
    ! The damping factor is introduced to EvanDerWaals. See the reference Phys. Rev. B 73, 205101 (2006).
    integer :: n, atomic_orbital
    !RepulsiveEnergy   Eq.(1.46)+(1.50)
    !dRepulsiveEnergy  Eq.(1.52)
!
    real(kind=DOUBLE_PRECISION) :: E_rest_non_vdW, dE_rest_non_vdW
!
    if (config%calc%genoOption%vanDerWaals_lambda > 0.0d0) then 
      vdW_lambda = config%calc%genoOption%vanDerWaals_lambda  ! value is set from file, if the dummy value
    else
      vdW_lambda = vdW_lambda_default                         ! default value
      config%calc%genoOption%vanDerWaals_lambda = vdW_lambda_default 
    endif   
!
     left_atom_spec=atm_element( left_atom)
    right_atom_spec=atm_element(right_atom)
    Za=AtomParameter( left_atom_spec)%num_val_elec+AtomParameter( left_atom_spec)%initial_charge
    Zb=AtomParameter(right_atom_spec)%num_val_elec+AtomParameter(right_atom_spec)%initial_charge
    ! Q = initial_charge (Not Mulliken charge)

    R=sqrt(x*x+y*y+z*z)

    if (R < 1.0d-10) then
      write(*,'(a,2i10,f30.20)')'ERROR(RepulsiveEnergy):i,j,R(i,j)=',left_atom, right_atom, R
      write(*,'(a)')'ERROR(RepulsiveEnergy):Two atoms are placed at the same position'
      stop
    endif

!
     E_rest_non_vdW=0.0d0
    dE_rest_non_vdW=0.0d0
!
!     
!   config%calc%interaction_range%cutoff_rest_non_vdW=huge(1.0d0)
!
    if (R < config%calc%interaction_range%cutoff_rest_non_vdW ) then
!
      t=0d0
     dt=0d0

     ! for A atom
      s=0d0
     ds=0d0
     do atomic_orbital=1, AtomParameter( left_atom_spec)%num_val_orb
       n   = AtomParameter( left_atom_spec)%principal(atomic_orbital)
       occ = AtomParameter( left_atom_spec)%initial_occupation(atomic_orbital)
       if(AtomParameter( left_atom_spec)%c2(atomic_orbital) == 0d0) then ! c2=0 => single zeta
          zeta= AtomParameter( left_atom_spec)%rescaled_zeta(atomic_orbital)
          s   =  s + occ*twobody2(n,zeta,R,d)
          ds  = ds + occ*d
       else
          z1 = AtomParameter( left_atom_spec)%rescaled_zeta(atomic_orbital)
          z2 = AtomParameter( left_atom_spec)%rescaled_zeta2(atomic_orbital)
          w1 = AtomParameter( left_atom_spec)%rescaled_c1(atomic_orbital)
          w2 = AtomParameter( left_atom_spec)%rescaled_c2(atomic_orbital)
          zeta = (z1+z2)/2
          zs = (z1*z2/(zeta*zeta))**(n+0.5d0)
          select case(genoMethod)
          case(genoICON)
             ! emulate erroneous ICON calculation
             s  =  s + occ*((twobody2(n,z1,R,d1)-1/R)*w1*w1 & 
                     + 2*(twobody2(n,zeta,R,d)-1/R)*w1*w2 + (twobody2(n,z2,R,d2)-1/R)*w2*w2+1/R)
             ds = ds + occ*((d1+1/R/R)*w1*w1 + 2*(d+1/R/R)*w1*w2*zs + (d2+1/R/R)*w2*w2-1/R/R)
          case default !! genoSTD
             s  =  s + occ*(twobody2(n,z1,R,d1)*w1*w1 + 2*twobody2(n,zeta,R,d)*w1*w2*zs + twobody2(n,z2,R,d2)*w2*w2)
             ds = ds + occ*(d1*w1*w1 + 2*d*w1*w2*zs + d2*w2*w2)
          end select
       end if
     end do
      t=Zb* s
     dt=Zb*ds
     ! for B atom
      s=0d0
     ds=0d0
     do atomic_orbital=1, AtomParameter(right_atom_spec)%num_val_orb
       n   = AtomParameter(right_atom_spec)%principal(atomic_orbital)
       occ = AtomParameter(right_atom_spec)%initial_occupation(atomic_orbital)
       if(AtomParameter(right_atom_spec)%c2(atomic_orbital) == 0d0) then ! c2=0 => single zeta
          zeta= AtomParameter(right_atom_spec)%rescaled_zeta(atomic_orbital)
          s   =  s + occ*twobody2(n,zeta,R,d)
          ds  = ds + occ*d
       else
          z1 = AtomParameter(right_atom_spec)%rescaled_zeta(atomic_orbital)
          z2 = AtomParameter(right_atom_spec)%rescaled_zeta2(atomic_orbital)
          w1 = AtomParameter(right_atom_spec)%rescaled_c1(atomic_orbital)
          w2 = AtomParameter(right_atom_spec)%rescaled_c2(atomic_orbital)
          zeta = (z1+z2)/2
          zs = (z1*z2/(zeta*zeta))**(n+0.5d0)
          select case(genoMethod)
          case(genoICON)
             ! emulate erroneous ICON calculation
             s  =  s + occ*((twobody2(n,z1,R,d1)-1/R)*w1*w1 & 
                     + 2*(twobody2(n,zeta,R,d)-1/R)*w1*w2 + (twobody2(n,z2,R,d2)-1/R)*w2*w2+1/R)
             ds = ds + occ*((d1+1/R/R)*w1*w1 + 2*(d+1/R/R)*w1*w2*zs + (d2+1/R/R)*w2*w2-1/R/R)
          case default
             s  =  s + occ*(twobody2(n,z1,R,d1)*w1*w1 + 2*twobody2(n,zeta,R,d)*w1*w2*zs + twobody2(n,z2,R,d2)*w2*w2)
             ds = ds + occ*(d1*w1*w1 + 2*d*w1*w2*zs + d2*w2*w2)
          end select
       end if
     end do
      t= t+Za* s
     dt=dt+Za*ds
!
     E_rest_non_vdW =   Za*Zb/R - 0.5d0 * t
    dE_rest_non_vdW = - Za*Zb/(R*R) - 0.5d0 * dt 

    endif  
!
!    write(*,'(a,2f20.10)')'R, t, cutoff (res) =', R, t
!

! Ref: Semiempirical van der Waals correction 
! to the density functional description of solids and molecular structures
! F. Ortmann, F. Bechstedt, W. G. Schmidt
! Phys. Rev. B 73, 205101 (2006)
! r_{cov} <-> van der Waals radius
!
!  Note : E_atompair_vdW_only = (1/2) * EvanDerWaals  
!            where the factor of (1/2) is one for the double counting
!
    EvanDerWaals  = 0.d0
    dEvanDerWaals = 0.d0
    if (config%calc%genoOption%vanDerWaals) then 
     if ( R < config%calc%interaction_range%cutoff_rest ) then
       Ia = AtomParameter(right_atom_spec)%ionization_energy
       Ib = AtomParameter( left_atom_spec)%ionization_energy
       aa = AtomParameter(right_atom_spec)%dipole_polarizability
       ab = AtomParameter( left_atom_spec)%dipole_polarizability
       vdW_R_a = AtomParameter(right_atom_spec)%vdW_radius
       vdW_R_b = AtomParameter( left_atom_spec)%vdW_radius
       !! London dispersion force
       EvanDerWaals  = - 1.5d0 * aa * ab * Ia * Ib / (Ia + Ib) / R**6
       dEvanDerWaals = - 6d0 * EvanDerWaals / R
       xab = R/(vdW_R_a + vdW_R_b)
       fd = fdamp(xab,dfd)
       dfd = dfd/(vdW_R_a + vdW_R_b)
       dEvanDerWaals = EvanDerWaals * dfd + dEvanderWaals * fd
       EvanDerWaals = EvanDerWaals * fd
       if(i_verbose >= 200)then
          write(*,*) 'vanDerWaals=', vanDerWaals
          write(*,*) "Ia,Ib,aa,ab=", Ia, Ib, aa, ab
          write(*,*) "EvanDerWaals =", EvanDerWaals
          write(*,*) "dEvanDerWaals =", dEvanDerWaals
       end if
     endif  
    end if

    RepulsiveEnergy     = 0.5d0 * (E_rest_non_vdW + EvanDerWaals) ! 0.5 for double counting
!
!   write(*,'(a,2f20.10)')'R, t, cutoff (res), E_Rep =', R, t
!
    if (present(E_atompair_vdW_only)) then
      E_atompair_vdW_only = 0.5d0 * EvanDerWaals  ! van der Waals only ( 0.5 for double counting )
    endif  
!
    dRepulsiveEnergy(1)=x
    dRepulsiveEnergy(2)=y
    dRepulsiveEnergy(3)=z
    dRepulsiveEnergy = 0.5d0 * (dE_rest_non_vdW +  dEvanDerWaals) * (dRepulsiveEnergy/R)
!
    if(i_verbose >= 200)then
       write(*,'("RepulsiveEnergy: genoMethod=",I3," (0:STD,100:NOLIN,200:ICON)")') genoMethod
       write(*,'("RepulsiveEnergy: R, Za, Zb = ",3ES22.14)') R, Za, Zb
       write(*,'("RepulsiveEnergy: val, dval = ",2ES22.14)') &
            RepulsiveEnergy, - 0.5d0 * (Za*Zb/(R*R) + 0.5d0 * dt - dEvanDerWaals)
       write(*,'("RepulsiveEnergy: dir cos   = ",3ES22.14)') x/R, y/R, z/R
    end if

  contains

    function fdamp(x,d)
      real(kind=DOUBLE_PRECISION)::fdamp,x,d
      intent(in) ::x
      intent(out)::d
      real(kind=DOUBLE_PRECISION) :: e, t, t1
      integer :: k
      t1 = x**(vdW_n-1)
      t  = x**(vdW_n)
      fdamp = 1d0 - dexp(-vdW_lambda * t)
      d = vdW_lambda * (vdW_n * t1) * dexp(-vdW_lambda * t)
      if(i_verbose >=200) then
         write(*,'("n=",I2," lambda=",ES21.14)') vdW_n, vdW_lambda
      end if
!       if(x>0.1)then
!          fdamp = 1d0 - dexp(-vdW_lambda * t)
!          d = vdW_lambda * (vdW_n * t1) * dexp(-vdW_lambda * t)
!       else
!          k = 1
!          fdamp = t
!          e = t
!          do while (abs(e) > 1d-17*abs(fdamp) )
!             fdamp = fdamp + e
!             k = k+1
!             e = e*t/k
!          end do
!          d = vdW_lambda * (vdW_n * t1) * (1d0+fdamp)
!      end if
         
    end function fdamp


    function twobody1(n,zeta,R,d)
      ! using infinite series expansion of incomplete Gamma Function
      real(kind=DOUBLE_PRECISION)::twobody1,zeta,R,e,s,zr2,g,t
      integer::n, k
      intent(in)::n,zeta,R
      real(kind=DOUBLE_PRECISION),intent(out) :: d ! derivative
      if(n<1) stop 'twobody1: out of range (n)'
      zr2=-2d0*zeta*R
      e=(zr2*zr2)/real((2*n)*(2*n+1),kind=DOUBLE_PRECISION)
      do k=2, n
         e=e*(zr2*zr2)/real((2*k-1)*(2*k),kind=DOUBLE_PRECISION)
      end do
      s=e
      g=-e*(2*n)*zeta/R   !  e* 2*n/zr2*(-2*zeta)
      d=g
      k=1
      do while (abs(e) > 1d-17*abs(s) )
         t=zr2/real(k,kind=DOUBLE_PRECISION)
         e=e*t*(k+2*n-1)/real(k+2*n+1,kind=DOUBLE_PRECISION)
         s=s+e
         g=g*t*(k+2*n)/real(k+2*n+1,kind=DOUBLE_PRECISION)
         d=d+g
         k=k+1
      end do
!      if(abs(zr2)<=1d0)     twobody=(1/real(n,kind=DOUBLE_PRECISION)-s)*zeta
      twobody1=(1/real(n,kind=DOUBLE_PRECISION)-s)*zeta
    end function twobody1

    function twobody2(n,zeta,R,d)
      ! using finite term expression of incomplete Gamma Function
      real(kind=DOUBLE_PRECISION)::twobody2,zeta,R,e,s,zr2,f,t
      integer::n, k
      intent(in)::n,zeta,R
      real(kind=DOUBLE_PRECISION),intent(out) :: d ! derivative
      if(n<1) stop 'twobody2: out of range (n)'
      zr2=2d0*zeta*R
      e=exp(-zr2)
      s=e
      f=e
      d=f
      do k=1, 2*n-1
         t=zr2/real(k,kind=DOUBLE_PRECISION)
         e=e*t*(2*n-k)/real(2*n-k+1,kind=DOUBLE_PRECISION)
         s=s+e
         f=f*t
         d=d+f
      end do
      d=d+f*zr2/real(k,kind=DOUBLE_PRECISION) ! k=2*n
      d=(d-1d0)/(R*R)
!     if(abs(zr2)>=1d0) twobody=(1d0-s)/R
      twobody2=(1d0-s)/R
    end function twobody2

  end function RepulsiveEnergy

end module M_qm_geno_Huckel
