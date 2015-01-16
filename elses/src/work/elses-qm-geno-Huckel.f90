!================================================================
! ELSES version 0.02
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module M_qm_geno_Huckel
  ! generalized non orthogoanal calculation, generalized Huckel
  use M_qm_domain, only: i_verbose, DOUBLE_PRECISION
  use M_qm_geno_Huckel_atom_params, &
       only: AtomParameter, NUMBER_OF_ORBITALS, NUMBER_OF_CONF_VOIP, &
       UseVOIP, EmulateICON, normDoubleZeta, sdelta=>small_delta, kappa
  use elses_mod_phys_const, only : eV4au
  implicit none

  private
  public :: SetOverlap, SetDiagonalElements, SetNondiagonalElements
  public :: RepulsiveEnergy
  public :: calc_SS_dSS

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
    type(TSKBond)            ::  Bond_matrix(NUMBER_OF_ORBITALS,NUMBER_OF_ORBITALS)
    type(TSKBond)            :: dBond_matrix(NUMBER_OF_ORBITALS,NUMBER_OF_ORBITALS)
    type(TVector3d)          :: oresult
    integer, intent(in)      :: left_atom_spec, right_atom_spec
    real(DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(DOUBLE_PRECISION), intent(out):: S(:,:), dS(:,:,:)
    real(DOUBLE_PRECISION)   :: R 
    real(DOUBLE_PRECISION)   :: SS(9,9), dSS(9,9)
    real(DOUBLE_PRECISION)   :: dSSV(3,9,9), sresult, deriv_mat, sign_factor
    real(DOUBLE_PRECISION)   :: rl, rm, rn
    real(DOUBLE_PRECISION)   :: Vss_sigma, Vsp_sigma, Vsd_sigma, Vps_sigma, Vds_sigma
    real(DOUBLE_PRECISION)   :: Vpp_sigma, Vpp_pi, Vpd_sigma, Vpd_pi, Vdp_sigma, Vdp_pi
    real(DOUBLE_PRECISION)   :: Vdd_sigma, Vdd_pi, Vdd_delta
    real(DOUBLE_PRECISION)   :: dVss_sigma, dVsp_sigma, dVsd_sigma, dVps_sigma, dVds_sigma
    real(DOUBLE_PRECISION)   :: dVpp_sigma, dVpp_pi, dVpd_sigma, dVpd_pi, dVdp_sigma, dVdp_pi
    real(DOUBLE_PRECISION)   :: dVdd_sigma, dVdd_pi, dVdd_delta
    real(DOUBLE_PRECISION), dimension(0:2)  :: v_SS, v_dSS
    real(DOUBLE_PRECISION), dimension(0:2,2):: v2tmp, v2dtmp
    real(DOUBLE_PRECISION)   :: zeta_a, zeta_b, zeta2_a(2), zeta2_b(2), c_a(2), c_b(2), normalization
    integer                  :: nval_orb_a, nval_orb_b, max_nval, max_l_a, max_l_b, nzeta_a, nzeta_b
    integer                  :: n, II, JJ
    integer                  :: l_a, l_b, mm
    integer                  :: npq_a, npq_b
    character(len=8)         :: orbital_name(9)

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

     SS=0d0
    dSS=0d0
    do l_b=0, max_l_b
       npq_b    =AtomParameter( left_atom_spec)%principal(l_b*l_b+1)
       zeta_b   =AtomParameter( left_atom_spec)%zeta(l_b*l_b+1)
       do l_a=0, max_l_a       ! The next line is an optimization for homogeneous case
          if(left_atom_spec==right_atom_spec .and. l_a < l_b ) cycle 
          npq_a =AtomParameter(right_atom_spec)%principal(l_a*l_a+1)
          zeta_a=AtomParameter(right_atom_spec)%zeta(l_a*l_a+1)
          
          if(l_a /= 2 .and. l_b /= 2) then ! No double-zeta
             !nik!
             if(l_b <= l_a)then
                call calc_SS_dSS(npq_b,npq_a,l_b,l_a,zeta_b,zeta_a,R,v_SS,v_dSS)
             else
                call calc_SS_dSS(npq_a,npq_b,l_a,l_b,zeta_a,zeta_b,R,v_SS,v_dSS)
             end if
          else
             zeta2_a(1)=zeta_a             ! Double-zeta
             zeta2_b(1)=zeta_b
             if(l_a == 2) then
                zeta2_a(2)=AtomParameter(right_atom_spec)%zeta2(l_a*l_a+1)
                nzeta_a=2
                c_a(1)    =AtomParameter(right_atom_spec)%c1(l_a*l_a+1)
                c_a(2)    =AtomParameter(right_atom_spec)%c2(l_a*l_a+1)
                normalization=normDoubleZeta(npq_a,zeta2_a(1),zeta2_a(2),c_a(1),c_a(2))
                if(abs(1-normalization) > 1d-6) c_a=c_a/normalization
             else
                nzeta_a=1
                c_a(1)    =1d0
             end if
             if(l_b == 2) then
                zeta2_b(2)=AtomParameter( left_atom_spec)%zeta2(l_b*l_b+1)
                nzeta_b=2
                c_b(1)    =AtomParameter( left_atom_spec)%c1(l_b*l_b+1)
                c_b(2)    =AtomParameter( left_atom_spec)%c2(l_b*l_b+1)
                normalization=normDoubleZeta(npq_b,zeta2_b(1),zeta2_b(2),c_b(1),c_b(2))
                if(abs(1-normalization) > 1d-6) c_b=c_b/normalization
             else
                nzeta_b=1
                c_b(1)    =1d0
             end if
              v2tmp=0d0
             v2dtmp=0d0
             do ii=1, nzeta_b
                do jj=1, nzeta_a
                   !nik!
                   if(l_b <= l_a)then
                      call calc_SS_dSS(npq_b,npq_a,l_b,l_a,zeta2_b(ii),zeta2_a(jj),R,v_SS,v_dSS)
                   else
                      call calc_SS_dSS(npq_a,npq_b,l_a,l_b,zeta2_a(jj),zeta2_b(ii),R,v_SS,v_dSS)
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
          ii=l_b*(l_b+1)
          jj=l_a*(l_a+1)
          do mm=0, min(l_a,l_b)
              SS(ii+mm+1,jj+mm+1)= v_SS(mm)  ! In mm=0 case, memory copy occurs twice.
              SS(ii-mm+1,jj-mm+1)= v_SS(mm)  ! So, this is not a smart code.
             dSS(ii+mm+1,jj+mm+1)=v_dSS(mm)  ! But we prefer this code, for simplicity.
             dSS(ii-mm+1,jj-mm+1)=v_dSS(mm)
          end do
       end do
       
       if(left_atom_spec==right_atom_spec) then ! Optimization for homogeneous case
          do l_a=0, l_b-1
             do mm=-l_a, l_a
                ii=l_b*(l_b+1)+mm+1
                jj=l_a*(l_a+1)+mm+1
                 SS(ii,jj)= SS(jj,ii)!nik!*(1-2*mod(l_a+l_b,2))
                dSS(ii,jj)=dSS(jj,ii)!nik!*(1-2*mod(l_a+l_b,2))
             end do
          end do
       end if
    end do

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
    ! case:nval_orb_b=1
    ! Vss_sigma/dVss_sigma
    if(nval_orb_b.ge.1) then
      if(nval_orb_a.ge.1) then
        Vss_sigma=SS(1,1)
        dVss_sigma=dSS(1,1)
         if(i_verbose.eq.50) then
            write(*,*) ' Vss_sigma=', Vss_sigma
            write(*,*) ' dVss_sigma=', dVss_sigma
         end if
      end if
      if(nval_orb_a.ge.4) then
    ! Vsp_sigma/dVsp_sigma
        Vsp_sigma=SS(1,3)
        dVsp_sigma=dSS(1,3)
         if(i_verbose.eq.50) then
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
         if(i_verbose.eq.50) then
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
         if(i_verbose.eq.50) then
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
         if(i_verbose.eq.50) then
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
          if(i_verbose.eq.50) then
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
           if(i_verbose.eq.50) then
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
           if(i_verbose.eq.50) then
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
           if(i_verbose.eq.50) then
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
    real(DOUBLE_PRECISION), intent(in)     :: rl, rm, rn
    real(DOUBLE_PRECISION), intent(out)    :: sresult, deriv_mat, sign_factor
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
    implicit none
    integer, parameter       :: MAX_FACT=25
    integer,                intent(in) :: npq_a, npq_b, l_a, l_b
    real(DOUBLE_PRECISION), intent(in) :: zeta_a, zeta_b, R
    real(DOUBLE_PRECISION) :: SS(0:), dSS(0:)
    real(DOUBLE_PRECISION) :: R_h
    real(DOUBLE_PRECISION) :: factor1, factor2, factor3, factor4, factor5, dfactor1, dfactor4
    real(DOUBLE_PRECISION) :: factjla, factjlpa, factjlqa
    real(DOUBLE_PRECISION) :: factjlb, factjlpb, factjlqb
    real(DOUBLE_PRECISION) :: factmk, factmk_p
    real(DOUBLE_PRECISION) :: AA(0:14), BB(0:14), dAA(0:14), dBB(0:14)
    real(DOUBLE_PRECISION) :: fact(0:MAX_FACT)
    integer :: j_a, j_a_max, ip_a, ip_a_max, iq_a, iq_a_max
    integer :: j_b, j_b_max, ip_b, ip_b_max, iq_b, iq_b_max
    integer :: k, k_prime, ii, mm
    ! Factorials are prepared and stored.
    ! fact(n)=n*(n-1)*(n-2)....2*1
    ! Definition of factorial is different from ICON-EDIT.
    ! ICON-EDIT: maximum of fact(n) is fact(25).
    ! added by m.ikeda, 2008/08/04, FLA, Japan.
    fact(0)=1.d0
    do ii=1, MAX_FACT
       fact(ii)=fact(ii-1)*real(ii,DOUBLE_PRECISION)
    enddo

    R_h= R/2.d0
    
    do ii=0, npq_a+npq_b
       AA(ii)=A(ii, R_h*(zeta_a+zeta_b), dAA(ii))
       BB(ii)=B(ii, R_h*(zeta_b-zeta_a), dBB(ii))
    end do
    
     SS=0d0
    dSS=0d0
    do mm=0, min(l_a,l_b)
       ! Eq.(1.26)
       ! SS( -2:2): Overlap matrix representation of Eq.(1.26)
       ! dSS(-2:2): derivative of SS(II,JJ) with r of right_atom.
       ! Without the direction vector (l,m,n).
       ! Upper half elements of SS(II, JJ) are calculated first.
       ! Principal quantum numbers
       ! Screening parameters
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
                               !if(-ip_a-ip_b-iq_a-iq_b-2*k > 0) stop
                               !if(-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime) > npq_a+npq_b) stop
                               !if(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime) > npq_a+npq_b)stop
                               SS(mm)=SS(mm)&
                                    +factor4*factor5* &
                                    AA(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k)* &
                                    BB(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime))
                               dSS(mm)=dSS(mm)&
                                    +dfactor4*factor5*&
                                    AA(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k)*  &
                                    BB(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)) &
                                    +factor4*factor5*0.5d0*(zeta_a+zeta_b)*   &
                                    dAA(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k)* &
                                    BB(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime)) &
                                    +factor4*factor5*0.5d0*(zeta_b-zeta_a)*   &
                                    AA(npq_a+npq_b-ip_a-ip_b-iq_a-iq_b-2*k)*  &
                                    dBB(l_a+l_b-ip_a-ip_b+iq_a+iq_b-2*(j_a+j_b+k_prime))
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    end do
  contains
    function A (n,lambda,dA) !(also returns dA)
      real(DOUBLE_PRECISION)              :: A
      integer, intent(in)       :: n
      real(DOUBLE_PRECISION), intent(in)  :: lambda
      real(DOUBLE_PRECISION), intent(out) :: dA
      integer                   :: k
      real(DOUBLE_PRECISION)              :: e

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
    function B (n,lambda,dB) ! (also returns dB)
      real(kind=DOUBLE_PRECISION):: B,dB,lambda,f,g
      integer     :: n,m,k
      intent(in)  :: n,lambda
      intent(out) :: dB
      k=iand(n,1)
      if(abs(lambda) < .5d0)then
         ! The more cut off radius increase and the more precision increase.
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
      end if
    end function B
!     function B (n,lambda,dB) !(also returns dB)
!       real(DOUBLE_PRECISION)              :: B
!       integer, intent(in)                 :: n
!       real(DOUBLE_PRECISION), intent(in)  :: lambda
!       real(DOUBLE_PRECISION), intent(out) :: dB
!       logical                             :: even
!       integer                             :: m
!       real(DOUBLE_PRECISION)              :: e,s1,s2,eps,c1

!       eps=1.d-16

!       if(mod(n,2)==0)then
!          even=.true.
!          c1=2.d0
!       else
!          even=.false.
!          c1=-2.d0
!       end if

!       !Eq.(1.38), Eq(1.40)
!       if (even) then
!          e=lambda/dble(n+3)
!          s1=1.d0/dble(n+1)+e*lambda/dble(2)
!          s2=e
!       else
!          e=1.d0/dble(n+2)
!          s1=e*lambda
!          s2=e
!       endif
!       m=1
!       do while (abs(e)>eps*abs(s2))
!          if (even) then
!             e=e*dble(n+2*m+1)/dble(n+2*m+3)/dble(2*m)/dble(2*m+1)*lambda**2
!             s1=s1+e*lambda/dble(2*m+2)
!             s2=s2+e
!          else
!             e=e*dble(n+2*m)/dble(n+2*m+2)/dble(2*m-1)/dble(2*m)*lambda**2
!             s1=s1+e*lambda/dble(2*m+1)
!             s2=s2+e
!          endif
!          m=m+1
!       enddo
!       B=c1*s1
!       dB=c1*s2 
!     end function B
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

  function RepulsiveEnergy(left_atom_spec, right_atom_spec,&
       x,y,z, dRepulsiveEnergy)
    implicit none
    real(kind=DOUBLE_PRECISION) :: RepulsiveEnergy
    integer, intent(in)         :: left_atom_spec,right_atom_spec
    real(kind=DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(kind=DOUBLE_PRECISION), intent(out):: dRepulsiveEnergy(3)
    real(kind=DOUBLE_PRECISION) :: Za, Zb, R, zeta, s, t, occ, d, dt, ds
    real(kind=DOUBLE_PRECISION) :: z1, z2, w1, w2, d1, d2, zs
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
