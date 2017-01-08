module M_qm_geno_Huckel
  ! generalized non orthogoanal calculation, generalized Huckel
  use M_qm_domain, only: i_verbose, DOUBLE_PRECISION
  use M_qm_geno_Huckel_atom_params, only: AtomParameter, &
  NUMBER_OF_ORBITALS, NUMBER_OF_CONF_VOIP, UseVOIP, EmulateICON
  implicit none

  private
  public :: SetOverlap, SetDiagonalElements, SetNondiagonalElements
  public :: RepulsiveEnergy

contains

 subroutine SetOverlap(left_atom_spec,right_atom_spec,x,y,z,S,dS)
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
    use MSlaterKosterTable
    use MSlaterKosterDerivative
    implicit none
!   type(TSKBond)            :: bond ! unused
    type(TSKBond)            :: Bond_matrix(NUMBER_OF_ORBITALS,NUMBER_OF_ORBITALS)
    type(TSKBond)            :: dBond_matrix(NUMBER_OF_ORBITALS,NUMBER_OF_ORBITALS)
    type(TVector3d)          :: oresult
    integer, intent(in)      :: left_atom_spec, right_atom_spec
    real(kind=8), intent(in) :: x, y, z ! r*(direction cosines)
    real(kind=8), intent(out):: S(:,:), dS(:,:,:)
!   real(kind=8),            :: fact(NUMBER_OF_FACTORIALS)
    real(kind=8)             :: fact(0:50)
    real(kind=8)             :: R, R_h, factor, factor1, factor2, factor3, factor4     
    real(kind=8)             :: factjla, factjlpa, factjlqa                    
    real(kind=8)             :: factjlb, factjlpb, factjlqb                    
    real(kind=8)             :: factmk, factmk_p                 
    real(kind=8)             :: AA, BB, dAA, dBB
    real(kind=8)             :: SS(9,9), SS11, SS12, SS21, SS22 !, dummy(3) ! unused
    real(kind=8)             :: dSS(9,9), dSS11, dSS12, dSS21, dSS22            
    real(kind=8)             :: dSSV(3,9,9), sresult, deriv_mat, sign_factor
    real(kind=8)             :: AA12, AA21, AA22, dAA12, dAA21, dAA22
    real(kind=8)             :: BB12, BB21, BB22, dBB12, dBB21, dBB22   
    real(kind=8)             :: dfactor1, dfactor4, factor5                                   
    real(kind=8)             :: factor12, factor21, factor22   
    real(kind=8)             :: factor4_12, factor4_21, factor4_22   
    real(kind=8)             :: dfactor4_12, dfactor4_21, dfactor4_22   
    real(kind=8)             :: c1_a, c2_a, f_a        
    real(kind=8)             :: c1_b, c2_b, f_b        
    real(kind=8)             :: zeta_a, zeta_b , zeta2_a, zeta2_b
    real(kind=8)             :: rl, rm, rn
    real(kind=8)             :: Vss_sigma, Vsp_sigma, Vsd_sigma, Vps_sigma, Vds_sigma
    real(kind=8)             :: Vpp_sigma, Vpp_pi, Vpd_sigma, Vpd_pi, Vdp_sigma, Vdp_pi
    real(kind=8)             :: Vdd_sigma, Vdd_pi, Vdd_delta
    real(kind=8)             :: dVss_sigma, dVsp_sigma, dVsd_sigma, dVps_sigma, dVds_sigma
    real(kind=8)             :: dVpp_sigma, dVpp_pi, dVpd_sigma, dVpd_pi, dVdp_sigma, dVdp_pi
    real(kind=8)             :: dVdd_sigma, dVdd_pi, dVdd_delta
!   real(kind=8)             :: BB_exact, BB12_exact, BB21_exact, BB22_exact ! unused
!   integer                  :: nBB, nBB12, nBB21, nBB22 ! unused
    integer                  :: nval_orb_a, nval_orb_b
    integer                  :: n, II, JJ, lvec(9), mvec(9), max_fact, IIMAX        
    integer                  :: l_a, m_a 
    integer                  :: l_b !, m_b ! unused
    integer                  :: npq_a, npq_b
    integer                  :: j_a, ip_a, iq_a, k, mm
    integer                  :: j_a_max, ip_a_max, iq_a_max
    integer                  :: j_b, ip_b, iq_b, k_prime
    integer                  :: j_b_max, ip_b_max, iq_b_max
    character(len=8)        :: orbital_name(9)

    ! {s,p} <-> {s,p}
    !Eq.(1.26)

    ! Factorials are prepared and stored.
    ! fact(n)=n*(n-1)*(n-2)....2*1
    ! Definition of factorial is different from ICON-EDIT.
    ! ICON-EDIT: maximum of fact(n) is fact(25).
    ! added by m.ikeda, 2008/08/04, FLA, Japan.

    ! Maximum of factorials.
      max_fact=25
      fact(0)=1.d0
      do n=1, max_fact
      fact(n)=fact(n-1)*dfloat(n)
      enddo 
    !
      if(i_verbose.eq.100) then
        write(*,*) ' max_fact=', max_fact
        do n=0, max_fact
        write(*,*) ' n=', n, ' fact(n)=', fact(n)
      enddo
      end if
    ! Prepare and store lvec and mvec.
      II=0
      do l_a=0, 2
      do m_a= -l_a, l_a 
      II=II+1
      lvec(II)= l_a
      mvec(II)= m_a
      enddo
      enddo
    !
        if(i_verbose.eq.100) then
        write(*,*) ' print out list vector lvec(II) mvec(II)'
        do II=1, 9
        write(*,*) ' II=', II,' lvec(II)=', lvec(II),' mvec(II)=', mvec(II)
        enddo
        end if
    ! Coding according to the document provided by Dr.Yamamoto.
    ! Eq.(1.26)
    ! added by m.ikeda, 2008/08/04, FLA, Japan.

    ! Store angular momentum 
    ! Orbital
    ! s       L=0 ===> num_val_orb(=nval)=1
    ! s,p     L=1 ===> num_val_orb(=nval)=4
    ! s,p,d   L=2 ===> num_val_orb(=nval)=9
    !
    ! Matrix:SS(II,JJ) Eq.(1.26)
    ! II\JJ         1    2    3    4    5    6    7    8    9
    !              L=0  L=1  L=1  L=1  L=2  L=2  L=2  L=2  L=2
    !              m=0  m=-1 m=0  m=1  m=-2 m=-1 m=0  m=1  m=2
    ! 1  L=0/m= 0  S11   0   S13   0    0    0   S17   0    0
    ! 2  L=1/m=-1   0  (S22)  0    0    0  (S26)  0    0    0
    ! 3  L=1/m= 0  S31*  0   S33   0    0    0   S37   0    0
    ! 4  L=1/m= 1   0    0    0   S44   0    0    0   S48   0
    ! 5  L=2/m=-2   0    0    0    0  (S55)  0    0    0    0
    ! 6  L=2/m=-1   0  (S62*) 0    0    0  (S66)  0    0    0
    ! 7  L=2/m= 0  S71*  0   S73*  0    0    0   S77   0    0
    ! 8  L=2/m= 1   0    0    0   S84*  0    0    0   S88   0
    ! 9  L=2/m= 2   0    0    0    0    0    0    0    0   S99
    !
    ! Non-vanishing elements of matrix SS(II,JJ)
    ! Upper half elements are calculated as in 
    ! Vss_sigma=S11
    ! Vsp_sigma=S13
    ! Vsd_sigma=S17
    ! Vpp_pi   =S22=S44
    ! Vpd_pi   =S26=S48
    ! Vpp_sigma=S33
    ! Vpd_sigma=S37
    ! Vdd_delta=S55=S99
    ! Vdd_pi   =S66=S88
    ! Vdd_sigma=S77
    ! Lower half elements are calculated using some techniques.
    ! By changing left_atom and right_atom, lower elements are obtained.
    ! When one calls Sktable routine, the (l, m, n) direction cosines are 
    ! changed to (-l, -m, -n).
    ! Vps_sigma=S31*=S(1,3) :left_atom and right_atom are exchanged.
    ! Vdp_pi   =S62*=S(2,6) :left_atom and right_atom are exchanged.
    !               =S(4,8) :left_atom and right atom are exchanged.
    ! Vds_sigma=S71*=S(1,7) :left_atom and right_atom are exchanged.
    ! Vdp_sigma=S73*=S(3,7) :left_atom and right_atom are exchanged.
    ! Vdp_pi   =S84*=S(4,8) :left_atom and right_atom are exchanged.
    ! Using sktable routine, (-l, -m, -n) dierction consines are used.

    ! store parameters into simple arguments.
    ! l+m numbers

    ! R = dsqrt(x*2+y*2+z*2)
      R=dsqrt(x*x+y*y+z*z)
      R_h= R/2.d0
    ! Eq.(1.26)
    ! SS(II, JJ): Overlap matrix representation of Eq.(1.26)
    ! dSS(II,JJ): derivative of SS(II,JJ) with r of right_atom.
    ! Without the direction vector (l,m,n).
    ! Upper half elements of SS(II, JJ) are calculated first.
      nval_orb_b=AtomParameter(left_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(right_atom_spec)%num_val_orb
      do II=1, nval_orb_b
      do JJ=II, nval_orb_a
      SS(II,JJ)= 0.d0
      dSS(II,JJ)=0.d0
      if((mvec(II).eq.mvec(JJ)).and.mvec(II).ge.0) THEN
    ! write(*,*) ' upper half elements'
    ! write(*,*) ' II=',II,' JJ=',JJ
    ! Angular momentum
      l_b=AtomParameter(left_atom_spec)%angular(II)
      l_a=AtomParameter(right_atom_spec)%angular(JJ)
    ! Pricipal quantum numbers
      npq_b=AtomParameter(left_atom_spec)%principal(II)
      npq_a=AtomParameter(right_atom_spec)%principal(JJ)
    ! Screening parameters
      zeta_b=AtomParameter(left_atom_spec)%zeta(II)
      zeta_a=AtomParameter(right_atom_spec)%zeta(JJ)
    ! mm=iabs(mvec(II))
      mm=mvec(II)
      factor1=(R_h)**(npq_a+npq_b+1)
      dfactor1=0.5d0*(npq_a+npq_b+1)*(R_h)**(npq_a+npq_b)
      factor2=(2.d0*zeta_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta_b)**(2*npq_b+1)/fact(2*npq_b)
      factor3=&
              +fact(l_a-mm)*dfloat(2*l_a+1)/2.d0/fact(l_a+mm) &
              *fact(l_b-mm)*dfloat(2*l_b+1)/2.d0/fact(l_b+mm)
   !  factor3=&
   !          +fact(l_a-mm)*fact(2*l_a+1)/2.d0/fact(l_a+mm) &
   !          *fact(l_b-mm)*fact(2*l_b+1)/2.d0/fact(l_b+mm)
      factor4=factor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      dfactor4=dfactor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
   !  j_a_max=iabs((l_a-mm)/2)
   !  write(*,*) ' l_a=', l_a, ' mm=', mm
   !  write(*,*) ' j_a_max=', j_a_max
      j_a_max=(l_a-mm)/2
   !  write(*,*) ' j_a_max=', j_a_max
      do j_a=0, j_a_max 
      factjla=&
              +fact(l_a)/fact(j_a)/fact(l_a-j_a) &
              *fact(2*(l_a-j_a))/fact(l_a-mm-2*j_a)
      ip_a_max=l_a-mm-2*j_a
      do ip_a=0, ip_a_max
      factjlpa=&
               +fact(l_a-mm-2*j_a)/fact(ip_a)/fact(l_a-mm-2*j_a-ip_a)
      iq_a_max=npq_a-l_a+2*j_a
      do iq_a=0, iq_a_max
      factjlqa=&
               +fact(npq_a-l_a+2*j_a)/fact(iq_a)/fact(npq_a-l_a+2*j_a-iq_a)
   !  j_b_max=iabs((l_b-mm)/2)
      j_b_max=(l_b-mm)/2
      do j_b=0, j_b_max
      factjlb=&
              +fact(l_b)/fact(j_b)/fact(l_b-j_b) &
              *fact(2*(l_b-j_b))/fact(l_b-mm-2*j_b)
      ip_b_max=l_b-mm-2*j_b
      do ip_b=0, ip_b_max
      factjlpb=&
               +fact(l_b-mm-2*j_b)/fact(ip_b)/fact(l_b-mm-2*j_b-ip_b)
      iq_b_max=npq_b-l_b+2*j_b
      do iq_b=0, iq_b_max
      factjlqb=&
               +fact(npq_b-l_b+2*j_b)/fact(iq_b)/fact(npq_b-l_b+2*j_b-iq_b)
      do k=0, mm 
      factmk=fact(mm)/fact(k)/fact(mm-k)
      do k_prime=0, mm
      factmk_p=fact(mm)/fact(k_prime)/fact(mm-k_prime)
      factor5=&
              +(-1.d0)**(ip_a+iq_a+j_a+j_b+k+k_prime+mm)*factmk*factmk_p &
              *factjla*factjlpa*factjlqa &
              *factjlb*factjlpb*factjlqb  
      AA=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta_a+zeta_b), dAA)
      BB=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta_b-zeta_a), dBB)
   !  if(i_verbose.eq.1) then
   !  write(*,*) ' A(2, 2.45664494110666, dAA)=', A(2, 2.45664494110666d0, dAA)
   !  write(*,*) ' A=', A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta_a+zeta_b), dAA)
   !  write(*,*) ' npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k=',npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k
   !  write(*,*) ' R_h*(zeta_a+zeta_b)=', R_h*(zeta_a+zeta_b)
   !  write(*,*) ' A(2, 2.45664494110666, dAA)=', A(2, R_h*(zeta_a+zeta_b), dAA)
   !  write(*,*) ' AA=', AA, ' dAA=',dAA
   !  write(*,*) ' l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)=', &
   !               l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
   !  write(*,*) ' R_h*(zeta_a-zeta_b)=',R_h*(zeta_a-zeta_b)
   !  write(*,*) ' BB=', BB, ' dBB=',dBB
   !  write(*,*) ' factor4*factor5*AA*BB=', factor4*factor5*AA*BB
   !  end if
   !
   !         nBB=  l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
   !         BB_exact=2.d0/dfloat(nBB+1)
   !  if(((zeta_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.0) then
   !       BB=BB_exact
   !       write(*,*) ' nBB=', nBB
   !  end if
   !  if(((zeta_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.1) then
   !       BB=0.d0     
   !       write(*,*) ' nBB=', nBB
   !  end if
      SS(II,JJ)=SS(II,JJ)&
                +factor4*factor5*AA*BB
      dSS(II,JJ)=dSS(II,JJ)&
                 +dfactor4*factor5*AA*BB &
                 +factor4*factor5*0.5d0*(zeta_a+zeta_b)*dAA*BB &
                 +factor4*factor5*0.5d0*(zeta_b-zeta_a)*AA*dBB
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
      enddo
      enddo
    !
    ! Eq.(1.26)
    ! SS(II, JJ): Overlap matrix representation of Eq.(1.26)
    ! Lower half elements of SS(II, JJ) are calculated secondly.
    ! There are some tricks to calculate lower elements.
    ! These techniques are required in order to use SlaterKoster routine.
    ! left_atom and right=atom are exchanged.
    !
      nval_orb_b=AtomParameter(right_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(left_atom_spec)%num_val_orb
      do II=1, nval_orb_b
      do JJ=II+1, nval_orb_a
      SS(JJ,II)= 0.d0
      dSS(JJ,II)= 0.d0
      if((mvec(II).eq.mvec(JJ)).and.mvec(II).ge.0) THEN
    ! write(*,*) ' Lower half element'
    ! write(*,*) ' II=',JJ,' JJ=',II
    ! Angular momentum
      l_b=AtomParameter(right_atom_spec)%angular(II)
      l_a=AtomParameter(left_atom_spec)%angular(JJ)
    ! Pricipal quantum numbers
      npq_b=AtomParameter(right_atom_spec)%principal(II)
      npq_a=AtomParameter(left_atom_spec)%principal(JJ)
    ! Screening parameters
      zeta_b=AtomParameter(right_atom_spec)%zeta(II)
      zeta_a=AtomParameter(left_atom_spec)%zeta(JJ)
    ! mm=iabs(mvec(II))
      mm=mvec(II)
      factor1=(R_h)**(npq_a+npq_b+1)
      dfactor1=0.5d0*(npq_a+npq_b+1)*(R_h)**(npq_a+npq_b)
      factor2=(2.d0*zeta_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta_b)**(2*npq_b+1)/fact(2*npq_b)
      factor3=&
              +fact(l_a-mm)*dfloat(2*l_a+1)/2.d0/fact(l_a+mm) &
              *fact(l_b-mm)*dfloat(2*l_b+1)/2.d0/fact(l_b+mm)
    ! factor3=&
    !         +fact(l_a-mm)*fact(2*l_a+1)/2.d0/fact(l_a+mm) &
    !         *fact(l_b-mm)*fact(2*l_b+1)/2.d0/fact(l_b+mm)
      factor4=factor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      dfactor4=dfactor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
    ! j_a_max=iabs((l_a-mm)/2)
      j_a_max=(l_a-mm)/2
      do j_a=0, j_a_max 
      factjla=&
              +fact(l_a)/fact(j_a)/fact(l_a-j_a) &
              *fact(2*(l_a-j_a))/fact(l_a-mm-2*j_a)
      ip_a_max=l_a-mm-2*j_a
      do ip_a=0, ip_a_max
      factjlpa=&
               +fact(l_a-mm-2*j_a)/fact(ip_a)/fact(l_a-mm-2*j_a-ip_a)
      iq_a_max=npq_a-l_a+2*j_a
      do iq_a=0, iq_a_max
      factjlqa=&
               +fact(npq_a-l_a+2*j_a)/fact(iq_a)/fact(npq_a-l_a+2*j_a-iq_a)
    ! j_b_max=iabs((l_b-mm)/2)
      j_b_max=(l_b-mm)/2
      do j_b=0, j_b_max
      factjlb=&
              +fact(l_b)/fact(j_b)/fact(l_b-j_b) &
              *fact(2*(l_b-j_b))/fact(l_b-mm-2*j_b)
      ip_b_max=l_b-mm-2*j_b
      do ip_b=0, ip_b_max
      factjlpb=&
               +fact(l_b-mm-2*j_b)/fact(ip_b)/fact(l_b-mm-2*j_b-ip_b)
      iq_b_max=npq_b-l_b+2*j_b
      do iq_b=0, iq_b_max
      factjlqb=&
               +fact(npq_b-l_b+2*j_b)/fact(iq_b)/fact(npq_b-l_b+2*j_b-iq_b)
      do k=0, mm 
      factmk=fact(mm)/fact(k)/fact(mm-k)
      do k_prime=0, mm
      factmk_p=fact(mm)/fact(k_prime)/fact(mm-k_prime)
      factor5=&
              +(-1.d0)**(ip_a+iq_a+j_a+j_b+k+k_prime+mm)*factmk*factmk_p &
              *factjla*factjlpa*factjlqa &
              *factjlb*factjlpb*factjlqb  
      AA=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta_a+zeta_b), dAA)
      BB=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta_b-zeta_a), dBB)
   !
   !         nBB=  l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
   !         BB_exact=2.d0/dfloat(nBB+1)
   !  if(((zeta_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.0) then
   !       BB=BB_exact
   !       write(*,*) ' nBB=', nBB
   !  end if
   !  if(((zeta_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.1) then
   !       BB=0.d0     
   !       write(*,*) ' nBB=', nBB
   !  end if
   !
      SS(JJ,II)=SS(JJ,II)&
                +factor4*factor5*AA*BB
      dSS(JJ,II)=dSS(JJ,II)&
                 +dfactor4*factor5*AA*BB &
                 +factor4*factor5*0.5d0*(zeta_a+zeta_b)*dAA*BB &
                 +factor4*factor5*0.5d0*(zeta_b-zeta_a)*AA*dBB
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
      enddo
      enddo
    !
    ! {s,p} <-> {d}(double-zeta)
    !Eq.(1.56)
    !
    !
    ! do II=1,9
    ! write(*,*) ' c1=', AtomParameter(left_atom_spec)%c1(II)
    ! end do
    !
      nval_orb_b=AtomParameter(left_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(right_atom_spec)%num_val_orb
      if(nval_orb_a.eq.9) then
        if(nval_orb_b.le.4) then
           IIMAX=nval_orb_b
         else if(nval_orb_b.eq.9) then
           IIMAX=4
        end if
      do II=1, IIMAX
      do JJ=5, nval_orb_a
      SS11= 0.d0
      dSS11= 0.d0
    ! write(*,*) ' II=', II, ' JJ=', JJ
      if((mvec(II).eq.mvec(JJ)).and.mvec(II).ge.0) THEN
    ! Angular momentum
      l_b=AtomParameter(left_atom_spec)%angular(II)
      l_a=AtomParameter(right_atom_spec)%angular(JJ)
    ! Pricipal quantum numbers
      npq_b=AtomParameter(left_atom_spec)%principal(II)
      npq_a=AtomParameter(right_atom_spec)%principal(JJ)
    ! Screening parameters
      zeta_b=AtomParameter(left_atom_spec)%zeta(II)
      zeta_a=AtomParameter(right_atom_spec)%zeta(JJ)
      zeta2_a=AtomParameter(right_atom_spec)%zeta2(JJ)
    ! double zeta coefficients
    ! c1_a=AtomParameter(left_atom_spec)%c1(II)
    ! c2_a=AtomParameter(left_atom_spec)%c2(II)
      c1_a=AtomParameter(right_atom_spec)%c1(JJ)
      c2_a=AtomParameter(right_atom_spec)%c2(JJ)
    ! write(*,*) ' c1_a=', c1_a, ' c2_a=', c2_a
    ! normalization constant
      factor= (4.d0*zeta_a*zeta2_a/(zeta_a+zeta2_a)**2)**(dfloat(npq_a)+1.d0/2.d0)
      f_a=dsqrt(c1_a**2+c2_a**2+2.d0*c1_a*c2_a*factor)
    ! mm=iabs(mvec(II))
      mm=mvec(II)
      factor1=(R_h)**(npq_a+npq_b+1)
      dfactor1=0.5d0*(npq_a+npq_b+1)*(R_h)**(npq_a+npq_b)
    ! correct m.ikeda, 10.03.2008
    ! dfactor1=0.5d0*(npq_a+npq_b)*(R_h)**(npq_a+npq_b)
      factor2=(2.d0*zeta2_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta_b)**(2*npq_b+1)/fact(2*npq_b)
      factor3=&
              +fact(l_a-mm)*dfloat(2*l_a+1)/2.d0/fact(l_a+mm) &
              *fact(l_b-mm)*dfloat(2*l_b+1)/2.d0/fact(l_b+mm)
    ! factor3=&
    !         +fact(l_a-mm)*fact(2*l_a+1)/2.d0/fact(l_a+mm) &
    !         *fact(l_b-mm)*fact(2*l_b+1)/2.d0/fact(l_b+mm)
      factor4=factor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      dfactor4=dfactor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
    ! j_a_max=iabs((l_a-mm)/2)
      j_a_max=(l_a-mm)/2
      do j_a=0, j_a_max 
      factjla=&
              +fact(l_a)/fact(j_a)/fact(l_a-j_a) &
              *fact(2*(l_a-j_a))/fact(l_a-mm-2*j_a)
      ip_a_max=l_a-mm-2*j_a
      do ip_a=0, ip_a_max
      factjlpa=&
               +fact(l_a-mm-2*j_a)/fact(ip_a)/fact(l_a-mm-2*j_a-ip_a)
      iq_a_max=npq_a-l_a+2*j_a
      do iq_a=0, iq_a_max
      factjlqa=&
               +fact(npq_a-l_a+2*j_a)/fact(iq_a)/fact(npq_a-l_a+2*j_a-iq_a)
    ! j_b_max=iabs((l_b-mm)/2)
      j_b_max=(l_b-mm)/2
      do j_b=0, j_b_max
      factjlb=&
              +fact(l_b)/fact(j_b)/fact(l_b-j_b) &
              *fact(2*(l_b-j_b))/fact(l_b-mm-2*j_b)
      ip_b_max=l_b-mm-2*j_b
      do ip_b=0, ip_b_max
      factjlpb=&
               +fact(l_b-mm-2*j_b)/fact(ip_b)/fact(l_b-mm-2*j_b-ip_b)
      iq_b_max=npq_b-l_b+2*j_b
      do iq_b=0, iq_b_max
      factjlqb=&
               +fact(npq_b-l_b+2*j_b)/fact(iq_b)/fact(npq_b-l_b+2*j_b-iq_b)
      do k=0, mm 
      factmk=fact(mm)/fact(k)/fact(mm-k)
      do k_prime=0, mm
      factmk_p=fact(mm)/fact(k_prime)/fact(mm-k_prime)
      factor5=&
              +(-1.d0)**(ip_a+iq_a+j_a+j_b+k+k_prime+mm)*factmk*factmk_p &
              *factjla*factjlpa*factjlqa &
              *factjlb*factjlpb*factjlqb 
       AA=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta2_a+zeta_b), dAA)
       BB=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta_b-zeta2_a), dBB)
   !
   !         nBB=  l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
   !         BB_exact=2.d0/dfloat(nBB+1)
   !  if(((zeta2_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.0) then
   !       BB=BB_exact
   !       write(*,*) ' nBB=', nBB
   !  end if
   !  if(((zeta2_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.1) then
   !       BB=0.d0     
   !       write(*,*) ' nBB=', nBB
   !  end if
   !
         SS11=SS11&
              +factor4*factor5*AA*BB
         dSS11=dSS11&
               +dfactor4*factor5*AA*BB &
               +factor4*factor5*0.5d0*(zeta2_a+zeta_b)*dAA*BB &
               +factor4*factor5*0.5d0*(zeta_b-zeta2_a)*AA*dBB
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    !     write(*,*) ' s, p (left) X d (right)'
    !     write(*,*) ' c1_a=', c1_a, ' c2_a=' , c2_a, ' f_a=', f_a
          SS(II,JJ)=c1_a/f_a*SS(II,JJ)+c2_a/f_a*SS11
          dSS(II,JJ)=c1_a/f_a*dSS(II,JJ)+c2_a/f_a*dSS11
      endif
      enddo
      enddo
      end if

    ! {s,p} <-> {d}(double-zeta)
    !Eq.(1.56)
    ! Lower haf elements are calculated by changing
    ! left_atom amd right_atom.
      nval_orb_b=AtomParameter(right_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(left_atom_spec)%num_val_orb
      if(nval_orb_a.eq.9) then
         if(nval_orb_b.le.4) then
         IIMAX=nval_orb_b
         else if(nval_orb_b.eq.9) then
         IIMAX=4
         end if
      do II=1, IIMAX
      do JJ=5, nval_orb_a
      SS11= 0.d0
    ! added m.ikeda, 10.03.2008
      dSS11=0.d0
      if((mvec(II).eq.mvec(JJ)).and.mvec(II).ge.0) THEN
    ! Angular momentum
      l_b=AtomParameter(right_atom_spec)%angular(II)
      l_a=AtomParameter(left_atom_spec)%angular(JJ)
    ! Pricipal quantum numbers
      npq_b=AtomParameter(right_atom_spec)%principal(II)
      npq_a=AtomParameter(left_atom_spec)%principal(JJ)
    ! Screening parameters
      zeta_b=AtomParameter(right_atom_spec)%zeta(II)
      zeta_a=AtomParameter(left_atom_spec)%zeta(JJ)
      zeta2_a=AtomParameter(left_atom_spec)%zeta2(JJ)
    ! double zeta coefficients
      c1_a=AtomParameter(left_atom_spec)%c1(JJ)
      c2_a=AtomParameter(left_atom_spec)%c2(JJ)
    ! normalization constant
      factor= (4.d0*zeta_a*zeta2_a/(zeta_a+zeta2_a)**2)**(dfloat(npq_a)+1.d0/2.d0)
      f_a=dsqrt(c1_a**2+c2_a**2+2.d0*c1_a*c2_a*factor)
    ! mm=iabs(mvec(II))
      mm=mvec(II)
      factor1=(R_h)**(npq_a+npq_b+1)
      dfactor1=0.5d0*(npq_a+npq_b+1)*(R_h)**(npq_a+npq_b)
      factor2=(2.d0*zeta2_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta_b)**(2*npq_b+1)/fact(2*npq_b)
      factor3=&
              +fact(l_a-mm)*dfloat(2*l_a+1)/2.d0/fact(l_a+mm) &
              *fact(l_b-mm)*dfloat(2*l_b+1)/2.d0/fact(l_b+mm)
    ! factor3=&
    !         +fact(l_a-mm)*fact(2*l_a+1)/2.d0/fact(l_a+mm) &
    !         *fact(l_b-mm)*fact(2*l_b+1)/2.d0/fact(l_b+mm)
      factor4=factor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      dfactor4=dfactor1*dsqrt(factor2*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
    ! j_a_max=iabs((l_a-mm)/2)
      j_a_max=(l_a-mm)/2
      do j_a=0, j_a_max 
      factjla=&
              +fact(l_a)/fact(j_a)/fact(l_a-j_a) &
              *fact(2*(l_a-j_a))/fact(l_a-mm-2*j_a)
      ip_a_max=l_a-mm-2*j_a
      do ip_a=0, ip_a_max
      factjlpa=&
               +fact(l_a-mm-2*j_a)/fact(ip_a)/fact(l_a-mm-2*j_a-ip_a)
      iq_a_max=npq_a-l_a+2*j_a
      do iq_a=0, iq_a_max
      factjlqa=&
               +fact(npq_a-l_a+2*j_a)/fact(iq_a)/fact(npq_a-l_a+2*j_a-iq_a)
    ! j_b_max=iabs((l_b-mm)/2)
      j_b_max=(l_b-mm)/2
      do j_b=0, j_b_max
      factjlb=&
              +fact(l_b)/fact(j_b)/fact(l_b-j_b) &
              *fact(2*(l_b-j_b))/fact(l_b-mm-2*j_b)
      ip_b_max=l_b-mm-2*j_b
      do ip_b=0, ip_b_max
      factjlpb=&
               +fact(l_b-mm-2*j_b)/fact(ip_b)/fact(l_b-mm-2*j_b-ip_b)
      iq_b_max=npq_b-l_b+2*j_b
      do iq_b=0, iq_b_max
      factjlqb=&
               +fact(npq_b-l_b+2*j_b)/fact(iq_b)/fact(npq_b-l_b+2*j_b-iq_b)
      do k=0, mm 
      factmk=fact(mm)/fact(k)/fact(mm-k)
      do k_prime=0, mm
      factmk_p=fact(mm)/fact(k_prime)/fact(mm-k_prime)
       factor5=&
              +(-1.d0)**(ip_a+iq_a+j_a+j_b+k+k_prime+mm)*factmk*factmk_p &
              *factjla*factjlpa*factjlqa &
              *factjlb*factjlpb*factjlqb  
       AA=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta2_a+zeta_b), dAA)
       BB=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta_b-zeta2_a), dBB)
   !
   !         nBB=  l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
   !         BB_exact=2.d0/dfloat(nBB+1)
   !  if(((zeta2_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.0) then
   !       BB=BB_exact
   !       write(*,*) ' nBB=', nBB
   !  end if
   !  if(((zeta2_a-zeta_b).eq.0.d0).and.mod(nBB,2).eq.1) then
   !       BB=0.d0     
   !       write(*,*) ' nBB=', nBB
   !  end if
   !
          SS11=SS11&
               +factor4*factor5*AA*BB
          dSS11=dSS11&
                +dfactor4*factor5*AA*BB &
                +factor4*factor5*0.5d0*(zeta2_a+zeta_b)*dAA*BB &
                +factor4*factor5*0.5d0*(zeta_b-zeta2_a)*AA*dBB
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    !     write(*,*) ' c1_a=', c1_a, ' f_a=', f_a
          SS(JJ,II)=c1_a/f_a*SS(JJ,II)+c2_a/f_a*SS11
          dSS(JJ,II)=c1_a/f_a*dSS(JJ,II)+c2_a/f_a*dSS11
      endif
      enddo
      enddo
      end if

    ! {d}   <-> {d}(both sides are double-zeta)
    !Eq.(1.58)
      nval_orb_b=AtomParameter(left_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(right_atom_spec)%num_val_orb
      if(nval_orb_a.eq.9.and.nval_orb_b.eq.9) then
      do II=5, nval_orb_b
    ! mod by m.ikeda, fla, 12.12.2008
    ! do JJ=5, nval_orb_a
    ! mod end m.ikeda, fla, 12.12.2008
      do JJ=II, nval_orb_a
      SS12= 0.d0
      SS21= 0.d0
      SS22= 0.d0
      dSS12= 0.d0
      dSS21= 0.d0
      dSS22= 0.d0
      if((mvec(II).eq.mvec(JJ)).and.mvec(II).ge.0) THEN
    ! Angular momentum
      l_b=AtomParameter(left_atom_spec)%angular(II)
      l_a=AtomParameter(right_atom_spec)%angular(JJ)
    ! Pricipal quantum numbers
      npq_b=AtomParameter(left_atom_spec)%principal(II)
      npq_a=AtomParameter(right_atom_spec)%principal(JJ)
    ! Screening parameters
      zeta_b=AtomParameter(left_atom_spec)%zeta(II)
      zeta2_b=AtomParameter(left_atom_spec)%zeta2(II)
      zeta_a=AtomParameter(right_atom_spec)%zeta(JJ)
      zeta2_a=AtomParameter(right_atom_spec)%zeta2(JJ)
    ! double zeta coefficients
      c1_b=AtomParameter(left_atom_spec)%c1(II)
      c2_b=AtomParameter(left_atom_spec)%c2(II)
      c1_a=AtomParameter(right_atom_spec)%c1(II)
      c2_a=AtomParameter(right_atom_spec)%c2(II)
    ! normalization constant
      factor= (4.d0*zeta_a*zeta2_a/(zeta_a+zeta2_a)**2)**(dfloat(npq_a)+1.d0/2.d0)
      f_a=dsqrt(c1_a**2+c2_a**2+2.d0*c1_a*c2_a*factor)
      factor= (4.d0*zeta_b*zeta2_b/(zeta_b+zeta2_b)**2)**(dfloat(npq_b)+1.d0/2.d0)
      f_b=dsqrt(c1_b**2+c2_b**2+2.d0*c1_b*c2_b*factor)
    ! magnetic quantum number
    ! mm=iabs(mvec(II))
      mm=mvec(II)
      factor1=(R_h)**(npq_a+npq_b+1)
      dfactor1=0.5d0*(npq_a+npq_b+1)*(R_h)**(npq_a+npq_b)

      factor12=(2.d0*zeta2_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta_b)**(2*npq_b+1)/fact(2*npq_b)
      factor21=(2.d0*zeta_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta2_b)**(2*npq_b+1)/fact(2*npq_b)
      factor22=(2.d0*zeta2_a)**(2*npq_a+1)/fact(2*npq_a)*(2.d0*zeta2_b)**(2*npq_b+1)/fact(2*npq_b)

      factor3=&
              +fact(l_a-mm)*dfloat(2*l_a+1)/2.d0/fact(l_a+mm) &
              *fact(l_b-mm)*dfloat(2*l_b+1)/2.d0/fact(l_b+mm)
    ! factor3=&
    !         +fact(l_a-mm)*fact(2*l_a+1)/2.d0/fact(l_a+mm) &
    !         *fact(l_b-mm)*fact(2*l_b+1)/2.d0/fact(l_b+mm)
      factor4_12=factor1*dsqrt(factor12*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      factor4_21=factor1*dsqrt(factor21*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      factor4_22=factor1*dsqrt(factor22*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)

      dfactor4_12=dfactor1*dsqrt(factor12*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      dfactor4_21=dfactor1*dsqrt(factor21*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)
      dfactor4_22=dfactor1*dsqrt(factor22*factor3)/2**(l_a+l_b)/fact(l_a)/fact(l_b)

    ! j_a_max=iabs((l_a-mm)/2)
      j_a_max=(l_a-mm)/2
      do j_a=0, j_a_max 
      factjla=&
              +fact(l_a)/fact(j_a)/fact(l_a-j_a) &
              *fact(2*(l_a-j_a))/fact(l_a-mm-2*j_a)
      ip_a_max=l_a-mm-2*j_a
      do ip_a=0, ip_a_max
      factjlpa=&
               +fact(l_a-mm-2*j_a)/fact(ip_a)/fact(l_a-mm-2*j_a-ip_a)
      iq_a_max=npq_a-l_a+2*j_a
      do iq_a=0, iq_a_max
      factjlqa=&
               +fact(npq_a-l_a+2*j_a)/fact(iq_a)/fact(npq_a-l_a+2*j_a-iq_a)
    ! j_b_max=iabs((l_b-mm)/2)
      j_b_max=(l_b-mm)/2
      do j_b=0, j_b_max
      factjlb=&
              +fact(l_b)/fact(j_b)/fact(l_b-j_b) &
              *fact(2*(l_b-j_b))/fact(l_b-mm-2*j_b)
      ip_b_max=l_b-mm-2*j_b
      do ip_b=0, ip_b_max
      factjlpb=&
               +fact(l_b-mm-2*j_b)/fact(ip_b)/fact(l_b-mm-2*j_b-ip_b)
      iq_b_max=npq_b-l_b+2*j_b
      do iq_b=0, iq_b_max
      factjlqb=&
               +fact(npq_b-l_b+2*j_b)/fact(iq_b)/fact(npq_b-l_b+2*j_b-iq_b)
      do k=0, mm 
      factmk=fact(mm)/fact(k)/fact(mm-k)
      do k_prime=0, mm
      factmk_p=fact(mm)/fact(k_prime)/fact(mm-k_prime)
       factor5=&
              +(-1.d0)**(ip_a+iq_a+j_a+j_b+k+k_prime+mm)*factmk*factmk_p &
              *factjla*factjlpa*factjlqa &
              *factjlb*factjlpb*factjlqb  

    ! zeta_b zeta2_a term
       AA12=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta2_a+zeta_b), dAA12)
       BB12=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta_b-zeta2_a), dBB12)
    !
    !   nBB12=l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
    !   write(*,*) ' nBB12=', nBB12, ' BB12=', BB12
    !   BB12_exact=2.d0/dfloat(nBB12+1)
    !   if((zeta2_a-zeta_b).eq.0.0.and.mod(nBB12,2).eq.0) then
    !   BB12=BB12_exact
    !   end if
    !   if((zeta2_a-zeta_b).eq.0.0.and.mod(nBB12,2).eq.1) then
    !   BB12=0.d0
    !   end if
    !
        SS12=SS12&
             +factor4_12*factor5*AA12*BB12
        dSS12=dSS12&
              +dfactor4_12*factor5*AA12*BB12 &
              +factor4_12*factor5*0.5d0*(zeta2_a+zeta_b)*dAA12*BB12 &
              +factor4_12*factor5*0.5d0*(zeta_b-zeta2_a)*AA12*dBB12
    ! zeta2_b zeta_a term
       AA21=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta_a+zeta2_b), dAA21)
       BB21=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta2_b-zeta_a), dBB21)
    !
    !   nBB21=l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
    !   BB21_exact=2.d0/dfloat(nBB21+1)
    !   if((zeta_a-zeta2_b).eq.0.0.and.mod(nBB21,2).eq.0) then
    !   BB21=BB21_exact
    !   end if
    !   if((zeta_a-zeta2_b).eq.0.0.and.mod(nBB21,2).eq.1) then
    !   BB21=0.d0
    !   end if
    !
        SS21=SS21&
              +factor4_21*factor5*AA21*BB21
       dSS21=dSS21&
              +dfactor4_21*factor5*AA21*BB21 &
              +factor4_21*factor5*0.5d0*(zeta_a+zeta2_b)*dAA21*BB21 &
              +factor4_21*factor5*0.5d0*(zeta2_b-zeta_a)*AA21*dBB21
    ! zeta2_b zeta2_a term
       AA22=A(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k, R_h*(zeta2_a+zeta2_b), dAA22)
       BB22=B(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime), R_h*(zeta2_b-zeta2_a), dBB22)
    !
    !   nBB22=l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)
    !   BB22_exact=2.d0/dfloat(nBB22+1)
    !   if((zeta2_a-zeta2_b).eq.0.0.and.mod(nBB22,2).eq.0) then
    !   BB22=BB22_exact
    !   end if
    !   if((zeta2_a-zeta2_b).eq.0.0.and.mod(nBB22,2).eq.1) then
    !   BB22=0.d0
    !   end if
    !
        SS22=SS22&
              +factor4_22*factor5*AA22*BB22
        dSS22=dSS22&
              +dfactor4_22*factor5*AA22*BB22 &
              +factor4_22*factor5*0.5d0*(zeta2_a+zeta2_b)*dAA22*BB22 &
              +factor4_22*factor5*0.5d0*(zeta2_b-zeta2_a)*AA22*dBB22
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    ! Eq.(1.55)
    !     write(*,*) ' c1_a=', c1_a, ' c2_a=', c2_a
    !     write(*,*) ' c1_b=', c1_b, ' c2_b=', c2_b
    !     write(*,*) ' f_a=', f_a, ' f_b=', f_b
        SS(II,JJ)=&
                   +c1_b*c1_a/f_b/f_a*SS(II,JJ) &
                   +c1_b*c2_a/f_b/f_a*SS12 &
                   +c2_b*c1_a/f_b/f_a*SS21 &
                   +c2_b*c2_a/f_b/f_a*SS22
        dSS(II,JJ)=&
                   +c1_b*c1_a/f_b/f_a*dSS(II,JJ) &
                   +c1_b*c2_a/f_b/f_a*dSS12 &
                   +c2_b*c1_a/f_b/f_a*dSS21 &
                   +c2_b*c2_a/f_b/f_a*dSS22
      endif
      enddo
      enddo
      end if
    ! Store SS into Overlap matrix S 
    ! From Lm representation to Sss_sigma, Ssp_sigma, Spp_sigma
    ! Spp_/pi, etc.
    ! matrix elements
    ! Vss_sigma, Vsp_sigma, Vsd_sigma
    !            Vps_sigma, Vds_sigma
    ! Vpp_sigma, Vpd_sigma
    ! Vpp_sigma, Vpd_sigma
    ! Vdd_sigma, Vdd_pi,    Vdd_delta
    ! derivative of matrix elements with respect to r.
    ! dVss_sigma, dVsp_sigma, dVsd_sigma
    !             dVps_sigma, dVds_sigma
    ! dVpp_sigma, dVpd_sigma
    ! dVpp_sigma, dVpd_sigma
    ! dVdd_sigma, dVdd_pi,    dVdd_delta
    ! Euler angle rotation
    ! coded by m.ikeda, FLA, Japan. 09.01.2008
      rl= x/R
      rm= y/R
      rn= z/R
      nval_orb_b=AtomParameter(left_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(right_atom_spec)%num_val_orb
    ! case:nval_orb_b=1
    ! Vss_sigma/dVss_sigma
    if(nval_orb_b.ge.1) then
      if(nval_orb_a.ge.1) then
        Vss_sigma=SS(1,1)
        dVss_sigma=dSS(1,1)
         if(i_verbose >= 50) then
            write(*,*) ' Vss_sigma=', Vss_sigma
            write(*,*) ' dVss_sigma=', dVss_sigma
         end if
      end if
      if(nval_orb_a.ge.4) then
    ! Vsp_sigma/dVsp_sigma
        Vsp_sigma=SS(1,3)
        dVsp_sigma=dSS(1,3)
         if(i_verbose >= 50) then
            write(*,*) ' Vsp_sigma=', Vsp_sigma
            write(*,*) ' dVsp_sigma=', dVsp_sigma
         end if
      end if
    ! added by m.ikeda,fla, 12.01.2008
      if(nval_orb_a.eq.9) then
    ! Vsd_sigma/dVsd_sigma
        Vsd_sigma=SS(1,7)
    ! correct by m.ikeda, fla, 12.11.2008
    ! dVsd_sigma=SS(1,7)
        dVsd_sigma=dSS(1,7)
    ! mod end, m.ikeda, fla, 12.11.2008
         if(i_verbose >= 50) then
            write(*,*) ' Vsd_sigma=', Vsd_sigma
            write(*,*) ' dVsd_sigma=', dVsd_sigma
         end if
       end if
    end if
    ! mod end m.ikeda, fla, 12.01.2008

    if(nval_orb_b.ge.4) then
      if(nval_orb_a.ge.1) then
    ! Vps_sigma/dVps_sigma
        Vps_sigma=SS(3,1)
        dVps_sigma=dSS(3,1)
         if(i_verbose >= 50) then
            write(*,*) ' Vps_sigma=', Vps_sigma
            write(*,*) ' dVps_sigma=', dVps_sigma
         end if
      end if
      if(nval_orb_a.ge.4) then
    ! Vpp_sigma/dVpp_sigma
        Vpp_sigma=SS(3,3)
        dVpp_sigma=dSS(3,3)
    ! Vpp_pi/dVpp_pi    
        Vpp_pi=SS(4,4)
        dVpp_pi=dSS(4,4)
         if(i_verbose >= 50) then
            write(*,*) ' Vpp_sigma=', Vpp_sigma
            write(*,*) ' Vpp_pi   =', Vpp_pi
            write(*,*) ' dVpp_sigma=', dVpp_sigma
            write(*,*) ' dVpp_pi   =', dVpp_pi
         end if
      end if
      if(nval_orb_a.eq.9) then
    ! Vsd_sigma/dVsd_sigma
         Vsd_sigma=SS(1,7)
         dVsd_sigma=dSS(1,7)
    ! Vpd_sigma/dVpd_sigma
         Vpd_sigma=SS(3,7)
         dVpd_sigma=dSS(3,7)
    ! Vpd_pi/dVpd_pi    
         Vpd_pi=SS(4,8) 
         dVpd_pi=dSS(4,8) 
          if(i_verbose >= 50) then
             write(*,*) ' Vsd_sigma=', Vsd_sigma
             write(*,*) ' Vpd_sigma=', Vpd_sigma
             write(*,*) ' Vpd_pi   =', Vpd_pi
             write(*,*) ' dVsd_sigma=', dVsd_sigma
             write(*,*) ' dVpd_sigma=', dVpd_sigma
             write(*,*) ' dVpd_pi   =', dVpd_pi
          end if
      end if
    end if
    ! added by m.ikeda,fla, 12.01.2008
    if(nval_orb_b.eq.9) then 
      if(nval_orb_a.ge.1) then
    ! Vsd_sigma/dVsd_sigma
         Vds_sigma=SS(7,1)
    ! correct by m.ikeda, 12.11.2008
    ! dVds_sigma=SS(7,1)
         dVds_sigma=dSS(7,1)
    ! corrcetion end , m.ikeda, fla, 12.11.2008
           if(i_verbose >= 50) then
              write(*,*) ' Vds_sigma=', Vds_sigma
              write(*,*) ' dVds_sigma=', dVds_sigma
           end if
      end if
    ! mod end m.ikeda, fla, 12.01.2008
      if(nval_orb_a.ge.4) then
    ! Vds_sigma/dVds_sigma
         Vds_sigma=SS(7,1)
         dVds_sigma=dSS(7,1)
    ! Vdp_sigma/dVdp_sigma
         Vdp_sigma=SS(7,3)
         dVdp_sigma=dSS(7,3)
    ! Vdp_pi/dVdp_pi
         Vdp_pi=SS(8,4)
         dVdp_pi=dSS(8,4)
           if(i_verbose >= 50) then
              write(*,*) ' Vds_sigma=', Vds_sigma
              write(*,*) ' Vdp_sigma=', Vdp_sigma
              write(*,*) ' Vdp_pi   =', Vdp_pi
              write(*,*) ' dVds_sigma=', dVds_sigma
              write(*,*) ' dVdp_sigma=', dVdp_sigma
              write(*,*) ' dVdp_pi   =', dVdp_pi
           end if
      end if
      if(nval_orb_a.eq.9) then
    ! Sdd_sigma/dSdd_sigma
         Vdd_sigma=SS(7,7)
         dVdd_sigma=dSS(7,7)
    ! Sdd_pi/dSdd_pi
         Vdd_pi=SS(8,8)
         dVdd_pi=dSS(8,8)
    ! Sdd_delta/dSdd_delta
         Vdd_delta=SS(9,9)
         dVdd_delta=dSS(9,9)
           if(i_verbose >= 50) then
              write(*,*) ' Vdd_sigma=', Vdd_sigma
              write(*,*) ' Vdd_pi   =', Vdd_pi   
              write(*,*) ' Vdd_delta=', Vdd_delta
              write(*,*) ' dVdd_sigma=', dVdd_sigma
              write(*,*) ' dVdd_pi   =', dVdd_pi   
              write(*,*) ' dVdd_delta=', dVdd_delta
           end if
      end if
    end if
    ! test for TSKbond_vector
    !         dVdd_sigma=0.d0
    !         dVdd_pi=0.d0
    !         dVdd_delta=0.d0
    ! rotate by using Slater-Koster table
    !    input :
    !        (x,y,z)/sqrt(x*x+y*y+z*z) (direction cosines)
    !        AtomParameter(left_atom_spec)%angular(:)
    !        AtomParameter(right_atom_spec)%angular(:)
    !        NB. ( (1+3+5) x (1+3+5) block )
    ! dS :derivative of overlap matrix
    !
    ! Overlap matrix elements
    ! TSKbond vector Table
    !   ss:bond%v(1)=Vss_sigma, other components=0.
    !   sp:bond%v(1)=Vsp_sigma
    !   sd:bond%v(1)=Vsd_sigma
    !   ps:bond%v(1)=Vps_sigma
    !   pp:bond%v(1)=Vpp_sigma, bond%v(2)=Vpp_pi
    !   pd:bond%v(1)=Vpd_sigma, bond%v(2)=Vpd_pi
    !   ds:bond%v(1)=Vds_sigma
    !   dp:bond%v(1)=Vdp_sigma, bond%v(2)=Vdp_pi
    !   dd:bond%v(1)=Vdd_sigma, bond%v(2)=Vdd_pi, bond%v(3)=Vdd_delta
    !   Other components of bond vector are not used in Slater-Koster routines.
    ! II/JJ  1    2    3    4    5    6    7    8    9
    !        s    x    y    z    xy   yz   zx   xx   zz
    ! 1  s   ss | sp   sp   sp | sd   sd   sd   sd   sd
    !        ---|--------------|-----------------------
    ! 2  x   ps | pp   pp   pp | pd   pd   pd   pd   pd
    ! 3  y   ps | pp   pp   pp | pd   pd   pd   pd   pd
    ! 4  z   ps | pp   pp   pp | pd   pd   pd   pd   pd
    !        ---|--------------|-----------------------
    ! 5  xy  ds | dp   dp   dp | dd   dd   dd   dd   dd
    ! 6  yz  ds | dp   dp   dp | dd   dd   dd   dd   dd
    ! 7  zx  ds | dp   dp   dp | dd   dd   dd   dd   dd
    ! 8  xx  ds | dp   dp   dp | dd   dd   dd   dd   dd
    ! 9  zz  ds | dp   dp   dp | dd   dd   dd   dd   dd
    ! 
    !  xx=x^2-y^2
    !  zz=3*z^2-r^2
    ! Store into Overlap matrix elements.
    ! added by m.ikeda, FLA, 09.27.2008.  
      nval_orb_b=AtomParameter(left_atom_spec)%num_val_orb
      nval_orb_a=AtomParameter(right_atom_spec)%num_val_orb
    ! Store bond vectors into matrix forms.
    ! Bond_matrix(9, 9)%v(i), dBond_matrix(9, 9)%v(i)
    ! Type(TSKbond)  :: Bond_matrix(9,9), dBond_matrix(9,9)
      do II=1, nval_orb_b
      do JJ=1, nval_orb_a
      Bond_matrix(II,JJ)%v(1)=0.d0
      Bond_matrix(II,JJ)%v(2)=0.d0
      Bond_matrix(II,JJ)%v(3)=0.d0
      Bond_matrix(II,JJ)%v(4)=0.d0
      dBond_matrix(II,JJ)%v(1)=0.d0
      dBond_matrix(II,JJ)%v(2)=0.d0
      dBond_matrix(II,JJ)%v(3)=0.d0
      dBond_matrix(II,JJ)%v(4)=0.d0
      if(II.eq.1) then
        if(JJ.eq.1) then
         Bond_matrix(II, JJ)%v(1)=Vss_sigma
         dBond_matrix(II, JJ)%v(1)=dVss_sigma
    !   write(*,*)' Bond_matrix(1,1)%v(1)=', Bond_matrix(II,JJ)%v(1) 
        else if(JJ.ge.2.and.JJ.le.4) then
         Bond_matrix(II, JJ)%v(1)=Vsp_sigma
         dBond_matrix(II, JJ)%v(1)=dVsp_sigma
        else if(JJ.ge.5.and.JJ.le.9) then
         Bond_matrix(II, JJ)%v(1)=Vsd_sigma
         dBond_matrix(II, JJ)%v(1)=dVsd_sigma
        end if
      endif
    !  
      if(II.ge.2.and.II.le.4) then
        if(JJ.eq.1) then
         Bond_matrix(II, JJ)%v(1)=Vps_sigma
         dBond_matrix(II, JJ)%v(1)=dVps_sigma
        else if(JJ.ge.2.and.JJ.le.4) then
         Bond_matrix(II, JJ)%v(1)=Vpp_sigma
         Bond_matrix(II, JJ)%v(2)=Vpp_pi
         dBond_matrix(II, JJ)%v(1)=dVpp_sigma
         dBond_matrix(II, JJ)%v(2)=dVpp_pi
        else if(JJ.ge.5.and.JJ.le.9) then
         Bond_matrix(II, JJ)%v(1)=Vpd_sigma
         Bond_matrix(II, JJ)%v(2)=Vpd_pi
         dBond_matrix(II, JJ)%v(1)=dVpd_sigma
         dBond_matrix(II, JJ)%v(2)=dVpd_pi
        end if
      end if
    !  
      if(II.ge.5.and.II.le.9) then
        if(JJ.eq.1) then
         Bond_matrix(II, JJ)%v(1)=Vds_sigma
         dBond_matrix(II, JJ)%v(1)=dVds_sigma
        else if(JJ.ge.2.and.JJ.le.4) then
         Bond_matrix(II, JJ)%v(1)=Vdp_sigma
         Bond_matrix(II, JJ)%v(2)=Vdp_pi
         dBond_matrix(II, JJ)%v(1)=dVdp_sigma
         dBond_matrix(II, JJ)%v(2)=dVdp_pi
        else if(JJ.ge.5.and.JJ.le.9) then
         Bond_matrix(II, JJ)%v(1)=Vdd_sigma
         Bond_matrix(II, JJ)%v(2)=Vdd_pi
         Bond_matrix(II, JJ)%v(3)=Vdd_delta
         dBond_matrix(II, JJ)%v(1)=dVdd_sigma
         dBond_matrix(II, JJ)%v(2)=dVdd_pi
         dBond_matrix(II, JJ)%v(3)=dVdd_delta
        end if
      end if
    !    write(*,*) ' II=',II, ' JJ=', JJ
    !    write(*,*) ' Bond_matrix(II,JJ)%v(1)=', Bond_matrix(II,JJ)%v(1)
    !    write(*,*) ' Bond_matrix(II,JJ)%v(2)=', Bond_matrix(II,JJ)%v(2)
    !    write(*,*) ' Bond_matrix(II,JJ)%v(3)=', Bond_matrix(II,JJ)%v(3)
    !    write(*,*) ' dBond_matrix(II,JJ)%v(1)=', dBond_matrix(II,JJ)%v(1)
    !    write(*,*) ' dBond_matrix(II,JJ)%v(2)=', dBond_matrix(II,JJ)%v(2)
    !    write(*,*) ' dBond_matrix(II,JJ)%v(3)=', dBond_matrix(II,JJ)%v(3)
      end do
      end do
    !  
    ! rotation using Slater-Koster coefficients.
    ! subroutine Rotation_Slater_Koster
    ! rotation of Eulerian angles.
    !  
    ! call Rotation_Slater_Koster(Bond_matrix, dBond_matrix, rl, rm, rn, result, oresult, sign_factor)
      do II=1, nval_orb_b
      do JJ=1, nval_orb_a
      call Rotation_Slater_Koster(Bond_matrix(II,JJ), dBond_matrix(II,JJ), II, JJ, & 
                                  rl, rm, rn, sresult, oresult, deriv_mat, sign_factor)
      S(II,JJ)=sresult
    ! derivative
    !   sign_factor   = 1.d0
        dS(1,II,JJ)   = sign_factor*oresult%x/R
        dS(2,II,JJ)   = sign_factor*oresult%y/R
        dS(3,II,JJ)   = sign_factor*oresult%z/R
    ! derivative
    !   sign_factor=1.d0
        dSSV(1,II,JJ) = deriv_mat*rl
        dSSV(2,II,JJ) = deriv_mat*rm
        dSSV(3,II,JJ) = deriv_mat*rn
    !
    !   dSSV(1,II,JJ) = deriv_mat*sign_factor*rl
    !   dSSV(2,II,JJ) = deriv_mat*sign_factor*rm
    !   dSSV(3,II,JJ) = deriv_mat*sign_factor*rn
    ! total differential of overlap matrix
        dS(1,II,JJ)   = dS(1,II,JJ)+dSSV(1,II,JJ)
        dS(2,II,JJ)   = dS(2,II,JJ)+dSSV(2,II,JJ)
        dS(3,II,JJ)   = dS(3,II,JJ)+dSSV(3,II,JJ)
    !   dS(1,II,JJ)   = dSSV(1,II,JJ)
    !   dS(2,II,JJ)   = dSSV(2,II,JJ)
    !   dS(3,II,JJ)   = dSSV(3,II,JJ)
    !
    !   if(i_verbose.eq.1) then
    !   write(*,*) ' S(', II, ',', JJ,')=', S(II,JJ)
    !   write(*,*) '    dS(   1 ,', II, ',', JJ,')=', dS(1,II,JJ)
    !   write(*,*) '    dS(   2 ,', II, ',', JJ,')=', dS(2,II,JJ)
    !   write(*,*) '    dS(   3 ,', II, ',', JJ,')=', dS(3,II,JJ)
    !   end if
    !
      end do
      end do
         if (i_verbose >= 50) then
    ! prepare orbital name
    !  orbital_name(i):     s    x    y    z    xy   yz   zx   xx   zz
    !   xx=x^2-y^2
    !   zz=3*z^2-r^2
         orbital_name(1)='   s    '
         orbital_name(2)='   x    '
         orbital_name(3)='   y    '
         orbital_name(4)='   z    '
         orbital_name(5)='   xy   '
         orbital_name(6)='   yz   '
         orbital_name(7)='   zx   '
         orbital_name(8)='x^2-y^2 '
         orbital_name(9)='3z^2-r^2'
    !
            write(*,'(A)') '<<<set_overlap_geno off-diagonal verbose mode>>>'
            write(*,'("left_atom_spec=",I6,", right_atom_spec=",I6)') left_atom_spec, right_atom_spec
            write(*,'(A)') 'Overlap matrix S(II,JJ)'
         if(nval_orb_a.le.4) then
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,nval_orb_a)
            do II=1,nval_orb_b 
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(S(II,JJ),JJ=1,nval_orb_a)
            end do
            write(*,'(/,A)') 'differential of Overlap matrix'
            write(*,'(A)') '   x derivative dS(1,II,JJ)'
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,nval_orb_a)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(1,II,JJ),JJ=1,nval_orb_a)
            end do
            write(*,'(/,A)') '   y derivative dS(2,II,JJ)'
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,nval_orb_a)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(2,II,JJ),JJ=1,nval_orb_a)
            end do
            write(*,'(/,A)') '   z derivative dS(3,II,JJ)'
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,nval_orb_a)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(3,II,JJ),JJ=1,nval_orb_a)
            end do
         else if(nval_orb_a.ge.5) then
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,4)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,4)
            do II=1,nval_orb_b 
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(S(II,JJ),JJ=1,4)
            end do
            write(*,'(14x,9(I11,11x))') (JJ,JJ=5,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=5,nval_orb_a)
            do II=1,nval_orb_b 
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(S(II,JJ),JJ=5,nval_orb_a)
            end do
     !
            write(*,'(/,A)') 'differential of Overlap matrix'
            write(*,'(A)') '   x derivative dS(1,II,JJ)'
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,4)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,4)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(1,II,JJ),JJ=1,4)
            end do
            write(*,'(14x,9(I11,11x))') (JJ,JJ=5,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=5,nval_orb_a)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(1,II,JJ),JJ=5,nval_orb_a)
            end do
     !
            write(*,'(/,A)') '   y derivative dS(2,II,JJ)'
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,4)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,4)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(2,II,JJ),JJ=1,4)
            end do
            write(*,'(14x,9(I11,11x))') (JJ,JJ=5,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=5,nval_orb_a)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(2,II,JJ),JJ=5,nval_orb_a)
            end do
     !
            write(*,'(/,A)') '   z derivative dS(3,II,JJ)'
            write(*,'(14x,9(I11,11x))') (JJ,JJ=1,4)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=1,4)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(3,II,JJ),JJ=1,4)
            end do
            write(*,'(14x,9(I11,11x))') (JJ,JJ=5,nval_orb_a)
            write(*,'(14x,9(7x,A,7x))') (orbital_name(JJ),JJ=5,nval_orb_a)
            do II=1,nval_orb_b
               write(*,'(I4,2x,A,9ES22.14)') II,orbital_name(II),(dS(3,II,JJ),JJ=5,nval_orb_a)
            end do
         end if
     !
            write(*,'(A)') '<<<end set_overlap_geno off-diagonal verbose mode>>>'
         end if
    !  
    ! rotate by using Slater-Koster table
    !    input :
    !        (x,y,z)/sqrt(x*x+y*y+z*z) (direction cosines)
    !        AtomParameter(left_atom_spec)%angular(:)
    !        AtomParameter(right_atom_spec)%angular(:)
    !        NB. ( (1+3+5) x (1+3+5) block )
    ! dS :derivative of overlap matrix


  contains
    subroutine Rotation_Slater_Koster(bond, dbond, II, JJ, rl, rm, rn, & 
                                      sresult, oresult, deriv_mat, sign_factor)
  ! added by m.ikeda, fla, Japan.          
  ! Euler angle rotations using Slater-Koster coefficients.
    implicit none
    type(TSKBond), intent(in)    :: bond, dbond 
    type(TVector3d), intent(out) :: oresult
    real(kind=8), intent(in)     :: rl, rm, rn
    real(kind=8), intent(out)    :: sresult, deriv_mat, sign_factor
    integer,      intent(in)     :: II, JJ
  ! corrceted by m.ikeda,fla,2008.10.27
  ! sign for lower half elements of derivative
  !   p-p, s-d, d-d no change.
  !   s-p, p-d  need sign change
      oresult = TVector3d(0d0, 0d0, 0d0)
      sign_factor=1.d0
      select case(II)
      case(1)
      select case(JJ)
      case(1)
    ! matrix element
        sresult  = Vss(bond, rl, rm, rn)
        oresult  = TVector3d(0d0, 0d0, 0d0)
    ! derivative of matrix element
        deriv_mat=Vss(dbond, rl, rm, rn)
      case(2)
    ! matrix element
        sresult  = Vsx(bond, rl, rm, rn)
        oresult  = Dsx(bond, rl, rm, rn)
        deriv_mat= Vsx(dbond , rl, rm, rn)
      case(3)
        sresult  = Vsy(bond, rl, rm, rn)
        oresult  = Dsy(bond, rl, rm, rn)
        deriv_mat= Vsy(dbond, rl, rm, rn)
      case(4)
        sresult  = Vsz(bond, rl, rm, rn)
        oresult  = Dsz(bond, rl, rm, rn)
        deriv_mat= Vsz(dbond, rl, rm, rn)
      case(5)
        sresult  = Vs_xy(bond, rl, rm, rn)
        oresult  = Ds_xy(bond, rl, rm, rn)
        deriv_mat= Vs_xy(dbond, rl, rm, rn)
      case(6)
        sresult  = Vs_yz(bond, rl, rm, rn)
        oresult  = Ds_yz(bond, rl, rm, rn)
        deriv_mat= Vs_yz(dbond, rl, rm, rn)
      case(7)
        sresult  = Vs_zx(bond, rl, rm, rn)
        oresult  = Ds_zx(bond, rl, rm, rn)
        deriv_mat= Vs_zx(dbond, rl, rm, rn)
      case(8)
        sresult  = Vs_xx(bond, rl, rm, rn)
        oresult  = Ds_xx(bond, rl, rm, rn)
        deriv_mat= Vs_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vs_zz(bond, rl, rm, rn)
        oresult  = Ds_zz(bond, rl, rm, rn)
        deriv_mat= Vs_zz(dbond, rl, rm, rn)
      end select
    case(2)
      select case(JJ)
      case(1)
        sresult  = Vsx(bond, -rl, -rm, -rn)
        oresult  = Dsx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vsx(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vxx(bond, rl, rm, rn)
        oresult  = Dxx(bond, rl, rm, rn)
        deriv_mat= Vxx(dbond, rl, rm, rn)
      case(3)
        sresult  = Vxy(bond, rl, rm, rn)
        oresult  = Dxy(bond, rl, rm, rn)
        deriv_mat= Vxy(dbond, rl, rm, rn)
      case(4)
        sresult  = Vxz(bond, rl, rm, rn)
        oresult  = Dxz(bond, rl, rm, rn)
        deriv_mat= Vxz(dbond, rl, rm, rn)
      case(5)
        sresult  = Vx_xy(bond, rl, rm, rn)
        oresult  = Dx_xy(bond, rl, rm, rn)
        deriv_mat= Vx_xy(dbond, rl, rm, rn)
      case(6)
        sresult  = Vx_yz(bond, rl, rm, rn)
        oresult  = Dx_yz(bond, rl, rm, rn)
        deriv_mat= Vx_yz(dbond, rl, rm, rn)
      case(7)
        sresult  = Vx_zx(bond, rl, rm, rn)
        oresult  = Dx_zx(bond, rl, rm, rn)
        deriv_mat= Vx_zx(dbond, rl, rm, rn)
      case(8)
        sresult  = Vx_xx(bond, rl, rm, rn)
        oresult  = Dx_xx(bond, rl, rm, rn)
        deriv_mat= Vx_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vx_zz(bond, rl, rm, rn)
        oresult  = Dx_zz(bond, rl, rm, rn)
        deriv_mat= Vx_zz(dbond, rl, rm, rn)
      end select
    case(3)
      select case(JJ)
      case(1)
        sresult  = Vsy(bond, -rl, -rm, -rn)
        oresult  = Dsy(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vsy(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vxy(bond, -rl, -rm, -rn)
        oresult  = Dxy(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxy(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vyy(bond, rl, rm, rn)
        oresult  = Dyy(bond, rl, rm, rn)
! corrected by m.ikeda, fla, 2008.11.25.
!       sign_factor=-1.d0
! irrelevant factor
        deriv_mat= Vyy(dbond, rl, rm, rn)
      case(4)
        sresult  = Vyz(bond, rl, rm, rn)
        oresult  = Dyz(bond, rl, rm, rn)
        deriv_mat= Vyz(dbond, rl, rm, rn)
      case(5)
        sresult  = Vy_xy(bond, rl, rm, rn)
        oresult  = Dy_xy(bond, rl, rm, rn)
        deriv_mat= Vy_xy(dbond, rl, rm, rn)
      case(6)
        sresult  = Vy_yz(bond, rl, rm, rn)
        oresult  = Dy_yz(bond, rl, rm, rn)
        deriv_mat= Vy_yz(dbond, rl, rm, rn)
      case(7)
        sresult  = Vy_zx(bond, rl, rm, rn)
        oresult  = Dy_zx(bond, rl, rm, rn)
        deriv_mat= Vy_zx(dbond, rl, rm, rn)
      case(8)
        sresult  = Vy_xx(bond, rl, rm, rn)
        oresult  = Dy_xx(bond, rl, rm, rn)
        deriv_mat= Vy_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vy_zz(bond, rl, rm, rn)
        oresult  = Dy_zz(bond, rl, rm, rn)
        deriv_mat= Vy_zz(dbond, rl, rm, rn)
      end select
    case(4)
      select case(JJ)
      case(1)
        sresult  = Vsz(bond, -rl, -rm, -rn)
        oresult  = Dsz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vsz(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vxz(bond, -rl, -rm, -rn)
        oresult =  Dxz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxz(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vyz(bond, -rl, -rm, -rn)
        oresult  = Dyz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vyz(dbond, -rl, -rm, -rn)
      case(4)
        sresult  = Vzz(bond, rl, rm, rn)
        oresult  = Dzz(bond, rl, rm, rn)
        deriv_mat= Vzz(dbond, rl, rm, rn)
      case(5)
        sresult  = Vz_xy(bond, rl, rm, rn)
        oresult  = Dz_xy(bond, rl, rm, rn)
        deriv_mat= Vz_xy(dbond, rl, rm, rn)
      case(6)
        sresult  = Vz_yz(bond, rl, rm, rn)
        oresult  = Dz_yz(bond, rl, rm, rn)
        deriv_mat= Vz_yz(dbond, rl, rm, rn)
      case(7)
        sresult  = Vz_zx(bond, rl, rm, rn)
        oresult  = Dz_zx(bond, rl, rm, rn)
        deriv_mat= Vz_zx(dbond, rl, rm, rn)
      case(8)
        sresult  = Vz_xx(bond, rl, rm, rn)
        oresult  = Dz_xx(bond, rl, rm, rn)
        deriv_mat= Vz_xx(dbond , rl, rm, rn)
      case(9)
        sresult  = Vz_zz(bond, rl, rm, rn)
        oresult  = Dz_zz(bond, rl, rm, rn)
        deriv_mat= Vz_zz(dbond, rl, rm, rn)
      end select
    case(5)
      select case(JJ)
      case(1)
        sresult  = Vs_xy(bond, -rl, -rm, -rn)
        oresult  = Ds_xy(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vs_xy(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vx_xy(bond, -rl, -rm, -rn)
        oresult  = Dx_xy(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vx_xy(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vy_xy(bond, -rl, -rm, -rn)
        oresult  = Dy_xy(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vy_xy(dbond, -rl, -rm, -rn)
      case(4)
        sresult  = Vz_xy(bond, -rl, -rm, -rn)
        oresult  = Dz_xy(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vz_xy(dbond, -rl, -rm, -rn)
      case(5)
        sresult  = Vxy_xy(bond, rl, rm, rn)
        oresult  = Dxy_xy(bond, rl, rm, rn)
        deriv_mat= Vxy_xy(dbond, rl, rm, rn)
      case(6)
        sresult  = Vxy_yz(bond, rl, rm, rn)
        oresult  = Dxy_yz(bond, rl, rm, rn)
        deriv_mat= Vxy_yz(dbond, rl, rm, rn)
      case(7)
        sresult  = Vxy_zx(bond, rl, rm, rn)
        oresult  = Dxy_zx(bond, rl, rm, rn)
        deriv_mat= Vxy_zx(dbond, rl, rm, rn)
      case(8)
        sresult  = Vxy_xx(bond, rl, rm, rn)
        oresult  = Dxy_xx(bond, rl, rm, rn)
        deriv_mat= Vxy_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vxy_zz(bond, rl, rm, rn)
        oresult  = Dxy_zz(bond, rl, rm, rn)
        deriv_mat= Vxy_zz(dbond, rl, rm, rn)
      end select
    case(6)
      select case(JJ)
      case(1)
        sresult  = Vs_yz(bond, -rl, -rm, -rn)
        oresult  = Ds_yz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vs_yz(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vx_yz(bond, -rl, -rm, -rn)
        oresult  = Dx_yz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vx_yz(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vy_yz(bond, -rl, -rm, -rn)
        oresult  = Dy_yz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vy_yz(dbond, -rl, -rm, -rn)
      case(4)
        sresult  = Vz_yz(bond, -rl, -rm, -rn)
        oresult  = Dz_yz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vz_yz(dbond, -rl, -rm, -rn)
      case(5)
        sresult  = Vxy_yz(bond, -rl, -rm, -rn)
        oresult  = Dxy_yz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxy_yz(dbond, -rl, -rm, -rn)
      case(6)
        sresult  = Vyz_yz(bond, rl, rm, rn)
        oresult  = Dyz_yz(bond, rl, rm, rn)
        deriv_mat= Vyz_yz(dbond, rl, rm, rn)
      case(7)
        sresult  = Vyz_zx(bond, rl, rm, rn)
        oresult  = Dyz_zx(bond, rl, rm, rn)
   ! mod by m.ikeda, 12.12.2008
   !    deriv_mat= Vyz_yz(dbond, rl, rm, rn)
        deriv_mat= Vyz_zx(dbond, rl, rm, rn)
   ! mod end m.ikeda, 12.12.2008
      case(8)
        sresult  = Vyz_xx(bond, rl, rm, rn)
        oresult  = Dyz_xx(bond, rl, rm, rn)
        deriv_mat= Vyz_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vyz_zz(bond, rl, rm, rn)
        oresult  = Dyz_zz(bond, rl, rm, rn)
        deriv_mat= Vyz_zz(dbond, rl, rm, rn)
      end select
    case(7)
      select case(JJ)
      case(1)
        sresult  = Vs_zx(bond, -rl, -rm, -rn)
        oresult  = Ds_zx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vs_zx(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vx_zx(bond, -rl, -rm, -rn)
        oresult  = Dx_zx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vx_zx(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vy_zx(bond, -rl, -rm, -rn)
        oresult  = Dy_zx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vy_zx(dbond, -rl, -rm, -rn)
      case(4)
        sresult  = Vz_zx(bond, -rl, -rm, -rn)
        oresult  = Dz_zx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vz_zx(dbond, -rl, -rm, -rn)
      case(5)
        sresult  = Vxy_zx(bond, -rl, -rm, -rn)
        oresult  = Dxy_zx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxy_zx(dbond, -rl, -rm, -rn)
      case(6)
        sresult  = Vyz_zx(bond, -rl, -rm, -rn)
        oresult  = Dyz_zx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vyz_zx(dbond, -rl, -rm, -rn)
      case(7)
        sresult  = Vzx_zx(bond, rl, rm, rn)
        oresult  = Dzx_zx(bond, rl, rm, rn)
        deriv_mat= Vzx_zx(dbond, rl, rm, rn)
      case(8)
        sresult  = Vzx_xx(bond, rl, rm, rn)
        oresult  = Dzx_xx(bond, rl, rm, rn)
        deriv_mat= Vzx_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vzx_zz(bond, rl, rm, rn)
        oresult  = Dzx_zz(bond, rl, rm, rn)
        deriv_mat= Vzx_zz(dbond, rl, rm, rn)
      end select
    case(8)
      select case(JJ)
      case(1)
        sresult  = Vs_xx(bond, -rl, -rm, -rn)
        oresult  = Ds_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vs_xx(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vx_xx(bond, -rl, -rm, -rn)
        oresult  = Dx_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vx_xx(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vy_xx(bond, -rl, -rm, -rn)
        oresult  = Dy_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vy_xx(dbond, -rl, -rm, -rn)
      case(4)
        sresult  = Vz_xx(bond, -rl, -rm, -rn)
        oresult  = Dz_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vz_xx(dbond, -rl, -rm, -rn)
      case(5)
        sresult  = Vxy_xx(bond, -rl, -rm, -rn)
        oresult  = Dxy_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxy_xx(dbond, -rl, -rm, -rn)
      case(6)
        sresult  = Vyz_xx(bond, -rl, -rm, -rn)
        oresult  = Dyz_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vyz_xx(dbond, -rl, -rm, -rn)
      case(7)
        sresult  = Vzx_xx(bond, -rl, -rm, -rn)
        oresult  = Dzx_xx(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vzx_xx(dbond, -rl, -rm, -rn)
      case(8)
        sresult  = Vxx_xx(bond, rl, rm, rn)
        oresult  = Dxx_xx(bond, rl, rm, rn)
        deriv_mat= Vxx_xx(dbond, rl, rm, rn)
      case(9)
        sresult  = Vxx_zz(bond, rl, rm, rn)
        oresult  = Dxx_zz(bond, rl, rm, rn)
        deriv_mat= Vxx_zz(dbond, rl, rm, rn)
      end select
    case(9)
      select case(JJ)
      case(1)
        sresult  = Vs_zz(bond, -rl, -rm, -rn)
        oresult  = Ds_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vs_zz(dbond, -rl, -rm, -rn)
      case(2)
        sresult  = Vx_zz(bond, -rl, -rm, -rn)
        oresult  = Dx_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vx_zz(dbond, -rl, -rm, -rn)
      case(3)
        sresult  = Vy_zz(bond, -rl, -rm, -rn)
        oresult  = Dy_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vy_zz(dbond, -rl, -rm, -rn)
      case(4)
        sresult  = Vz_zz(bond, -rl, -rm, -rn)
        oresult  = Dz_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vz_zz(dbond, -rl, -rm, -rn)
      case(5)
        sresult  = Vxy_zz(bond, -rl, -rm, -rn)
        oresult  = Dxy_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxy_zz(dbond, -rl, -rm, -rn)
      case(6)
        sresult  = Vyz_zz(bond, -rl, -rm, -rn)
        oresult  = Dyz_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vyz_zz(dbond, -rl, -rm, -rn)
      case(7)
        sresult  = Vzx_zz(bond, -rl, -rm, -rn)
        oresult  = Dzx_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vzx_zz(dbond, -rl, -rm, -rn)
      case(8)
        sresult  = Vxx_zz(bond, -rl, -rm, -rn)
        oresult  = Dxx_zz(bond, -rl, -rm, -rn)
        sign_factor=-1.d0
        deriv_mat= Vxx_zz(dbond, -rl, -rm, -rn)
      case(9)
        sresult  = Vzz_zz(bond, rl, rm, rn)
        oresult  = Dzz_zz(bond, rl, rm, rn)
        deriv_mat= Vzz_zz(dbond, rl, rm, rn)
      end select
    end select 
    ! rotate by using Slater-Koster table
    !    input :
    !        (x,y,z)/sqrt(x*x+y*y+z*z) (direction cosines)
    !        AtomParameter(left_atom_spec)%angular(:)
    !        AtomParameter(right_atom_spec)%angular(:)
    !        NB. ( (1+3+5) x (1+3+5) block )
    end subroutine Rotation_Slater_Koster

    function A (n,lambda,dA) !(also returns dA)
      real(kind=8)              :: A
      integer, intent(in)       :: n
      real(kind=8), intent(in)  :: lambda
      real(kind=8), intent(out) :: dA
      integer                   :: k
      real(kind=8)              :: e

      !Eq.(1.37), Eq(1.39)
      e=exp(-lambda)/lambda
      A=e
      dA=-e*(lambda+1.d0)/lambda
      do k=1,n
         e=e*(n-k+1)/lambda
         A=A+e
         dA=dA-e*(lambda+dble(k)+1.d0)/lambda
      enddo
    end function A
    function B (n,lambda,dB) !(also returns dB)
      real(kind=8)              :: B
      integer, intent(in)       :: n
      real(kind=8), intent(in)  :: lambda
      real(kind=8), intent(out) :: dB
      logical                   :: even
      integer                   :: m
      real(kind=8)              :: e,s1,s2,eps,c1

      eps=1.d-16

      if(mod(n,2)==0)then
         even=.true.
         c1=2.d0
      else
         even=.false.
         c1=-2.d0
      end if

      !Eq.(1.38), Eq(1.40)
      if (even) then
         e=lambda/dble(n+3)
         s1=1.d0/dble(n+1)+e*lambda/dble(2)
         s2=e
      else
         e=1.d0/dble(n+2)
         s1=e*lambda
         s2=e
      endif
      m=1
      do while (abs(e)>eps*abs(s2))
         if (even) then
            e=e*dble(n+2*m+1)/dble(n+2*m+3)/dble(2*m)/dble(2*m+1)*lambda**2
            s1=s1+e*lambda/dble(2*m+2)
            s2=s2+e
         else
            e=e*dble(n+2*m)/dble(n+2*m+2)/dble(2*m-1)/dble(2*m)*lambda**2
            s1=s1+e*lambda/dble(2*m+1)
            s2=s2+e
         endif
         m=m+1
      enddo
      B=c1*s1
      dB=c1*s2 
    end function B
  end subroutine SetOverlap

  subroutine SetDiagonalElements(atom_spec,Occ,Ham)
    use elses_mod_phys_const
    integer, intent(in)      :: atom_spec
    real(kind=8), intent(in) :: Occ(:)
    real(kind=8), intent(out):: Ham(:,:)
    ! VOIP(inital_charge)
    ! If d-orbital exists, use VOIP_D & initial_occupation instead.
    integer :: i, j, num_val_orb
    real(8) :: q, VOIP_spd(3,9), D_sigma, D_pi, D_delta
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

  subroutine SetNondiagonalElements(left_atom_spec, right_atom_spec,&
       x,y,z, S, left_H_diagonal, right_H_diagonal, Ham, dS, dHam)
    use M_qm_domain
    use elses_mod_phys_const, only : angst,ev4au
    integer, intent(in)      :: left_atom_spec,right_atom_spec
    real(kind=8), intent(in) :: x, y, z ! r*(direction cosines)
    real(kind=8), intent(in) :: S(:,:), dS(:,:,:)
    real(kind=8), intent(in) :: left_H_diagonal(:,:), right_H_diagonal(:,:)
    real(kind=8), intent(out):: Ham(:,:), dHam(:,:,:)
    integer :: nval1,nval2,ja1,ja2,w
    real(8) :: r(3),lhd,rhd
    integer :: ja(2),atom_s(2)
    real(8) :: rn(2),r0,delta
    real(8) :: ff,ff2,cl1,cl2,zl1,zl2,hi1,hi2
    integer :: nl,i
    r(1:3)=(/ x,y,z /)
    atom_s(1:2)=(/ left_atom_spec,right_atom_spec /)
    nval1=AtomParameter(atom_s(1))%num_val_orb
    nval2=AtomParameter(atom_s(2))%num_val_orb

    do ja2=1,nval2
      do ja1=1,nval1
        ja(1:2)=(/ ja1,ja2 /)

        ! Eq.(2.2)
        hi1=left_H_diagonal(ja(1),ja(1))
        hi2=right_H_diagonal(ja(2),ja(2))
        delta=(hi1-hi2)/(hi1+hi2)

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
            if ( EmulateICON ) then
               !JPC 1989 (93) 5366-5371
               ff2=cl1*cl1*zl1 + cl2*cl2*zl2 + 0.5D0*&
                    ((4.0D0*zl1*zl2)**(nl+0.5D0))/((zl1+zl2)**(2.0D0*nl))
            else
               ff2=cl1*cl1*zl1 + cl2*cl2*zl2 + cl1*cl2*&
                    ((4.0D0*zl1*zl2)**(nl+0.5D0))/((zl1+zl2)**(2.0D0*nl))
            end if
            rn(i)=nl*ff/ff2
          endif
        enddo
        r0=rn(1)+rn(2)

        ! Eq.(2.7), dH_ij
        lhd=left_H_diagonal(ja1,ja1)
        rhd=right_H_diagonal(ja2,ja2)
        Ham(ja1,ja2)=K(r,r0,delta)*(lhd+rhd)/2.0D0*S(ja1,ja2)
!
        if (i_verbose >= 50) then
           write(*,*)'  ja1,ja2     =',ja1,ja2
           write(*,*)'        K     =',K(r,r0,delta)
           write(*,*)'      H_ii[eV]=',lhd*ev4au
           write(*,*)'      H_jj[eV]=',rhd*ev4au
           write(*,*)'(H_ii+H_jj)/2 =',(lhd+rhd)/2.0D0*ev4au
           write(*,*)'      H_ij[eV]=',Ham(ja1,ja2)*ev4au
        endif        
!
        do w=1,3
          dHam(w,ja1,ja2)=(lhd+rhd)/2.0D0 *&
           (dK(w,r,r0,delta)*S(ja1,ja2)+K(r,r0,delta)*dS(w,ja1,ja2))
        enddo

        if(i_verbose >= 200) then
          write(*,'(a,4(i1,a),e13.6)') &
            'H(',ja1,',',ja2,',',left_atom_spec,',',right_atom_spec,&
            ')=',Ham(ja1,ja2)
        endif
      enddo
    enddo
  contains

    function K(r,r0,delta) result(kk)
      real(8),intent(in) :: r(3),r0,delta
      real(8) :: rr,q,Kp,kk
      real(8),parameter :: delt=0.35D0*angst, kappa=1.0D0
      rr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) !a.u.

      ! Eq.(2.4) - Eq.(2.6)
      q=((rr-r0)-dabs(rr-r0))*delt
      q=1.0D0+q*q
      Kp=kappa + delta*delta - delta*delta*delta*delta*kappa
      kk=1.0D0 + Kp*exp(-1.0D0*delt*(rr-r0))/q

    end function K

    function dK(j,r,r0,delta) result(dkk)
      integer,intent(in) :: j
      real(8),intent(in) :: r(3),r0,delta
      real(8) :: rr,q,Kp,dkk
      real(8),parameter :: delt=0.35D0*angst, kappa=1.0D0
      rr=dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) !a.u.

      ! Eq.(2.5) 
      Kp=kappa + delta*delta - delta*delta*delta*delta*kappa

      if (rr.ge.r0) then
        ! Eq.(2.4), Eq.(2.6)
        q=1.0D0
        dkk= -1.0D0*delt*Kp*exp(-1.0D0*delt*(rr-r0))*r(j)/rr
      else
        ! Eq.(2.4), Eq.(2.6)
        q=2.0D0*(rr-r0)*delt
        q=1.0D0+q*q
        dkk= -1.0D0*delt*Kp*exp(-1.0D0*delt*(rr-r0))*r(j)/rr*&
          (q+8.0D0*delt*(rr-r0))/q/q
      endif
    end function dK
  end subroutine SetNondiagonalElements

  function RepulsiveEnergy(left_atom_spec, right_atom_spec,&
       x,y,z, dRepulsiveEnergy)
    implicit none
    integer, parameter :: DOUBLE_PRECISION=kind(1d0)
    real(kind=DOUBLE_PRECISION) :: RepulsiveEnergy
    integer, intent(in)      :: left_atom_spec,right_atom_spec
    real(kind=DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(kind=DOUBLE_PRECISION), intent(out):: dRepulsiveEnergy(3)
    real(kind=DOUBLE_PRECISION) :: Za, Zb, R, zeta, s, t, occ, d, dt, ds, z1, z2, w1, w2, d1, d2, zs
    integer :: n, atomic_orbital
    !RepulsiveEnergy   Eq.(1.46)+(1.50)
    !dRepulsiveEnergy  Eq.(1.52)

    Za=AtomParameter( left_atom_spec)%num_val_elec+AtomParameter( left_atom_spec)%initial_charge
    Zb=AtomParameter(right_atom_spec)%num_val_elec+AtomParameter(right_atom_spec)%initial_charge
    ! Q = initial_charge (Not Mulliken charge)

    R=sqrt(x*x+y*y+z*z)

    ! for A atom
     s=0d0
    ds=0d0
    do atomic_orbital=1, AtomParameter( left_atom_spec)%num_val_orb
       n   = AtomParameter( left_atom_spec)%principal(atomic_orbital)
       occ = AtomParameter( left_atom_spec)%initial_occupation(atomic_orbital)
       if(AtomParameter( left_atom_spec)%orbital(atomic_orbital) /= 'd') then
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
          if( EmulateICON ) then
             ! emulate erroneous ICON calculation
             s  =  s + occ*((twobody2(n,z1,R,d1)-1/R)*w1*w1 & 
&                    + 2*(twobody2(n,zeta,R,d)-1/R)*w1*w2 + (twobody2(n,z2,R,d2)-1/R)*w2*w2+1/R)
             ds = ds + occ*((d1+1/R/R)*w1*w1 + 2*(d+1/R/R)*w1*w2*zs + (d2+1/R/R)*w2*w2-1/R/R)
          else
             s  =  s + occ*(twobody2(n,z1,R,d1)*w1*w1 + 2*twobody2(n,zeta,R,d)*w1*w2*zs + twobody2(n,z2,R,d2)*w2*w2)
             ds = ds + occ*(d1*w1*w1 + 2*d*w1*w2*zs + d2*w2*w2)
          end if
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
       if(AtomParameter( left_atom_spec)%orbital(atomic_orbital) /= 'd') then
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
          if( EmulateICON ) then
             ! emulate erroneous ICON calculation
             s  =  s + occ*((twobody2(n,z1,R,d1)-1/R)*w1*w1 & 
&                    + 2*(twobody2(n,zeta,R,d)-1/R)*w1*w2 + (twobody2(n,z2,R,d2)-1/R)*w2*w2+1/R)
             ds = ds + occ*((d1+1/R/R)*w1*w1 + 2*(d+1/R/R)*w1*w2*zs + (d2+1/R/R)*w2*w2-1/R/R)
          else
             s  =  s + occ*(twobody2(n,z1,R,d1)*w1*w1 + 2*twobody2(n,zeta,R,d)*w1*w2*zs + twobody2(n,z2,R,d2)*w2*w2)
             ds = ds + occ*(d1*w1*w1 + 2*d*w1*w2*zs + d2*w2*w2)
          end if
       end if
    end do
     t= t+Za* s
    dt=dt+Za*ds

    RepulsiveEnergy  = 0.5d0 * (Za*Zb/R - 0.5d0 * t)

    dRepulsiveEnergy(1)=x
    dRepulsiveEnergy(2)=y
    dRepulsiveEnergy(3)=z
    dRepulsiveEnergy = - 0.5d0 * (Za*Zb/(R*R) + 0.5d0 * dt) * (dRepulsiveEnergy/R)

    if(i_verbose >= 200)then
       write(*,'("RepulsiveEnergy: EmulateICON:",L)') EmulateICON
       write(*,'("RepulsiveEnergy: R, Za, Zb = ",3ES22.14)') R, Za, Zb
       write(*,'("RepulsiveEnergy: val, dval = ",2ES22.14)') &
            RepulsiveEnergy, - 0.5d0 * (Za*Zb/(R*R) + 0.5d0 * dt)
       write(*,'("RepulsiveEnergy: dir cos   = ",3ES22.14)') x/R, y/R, z/R
    end if

  contains

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
