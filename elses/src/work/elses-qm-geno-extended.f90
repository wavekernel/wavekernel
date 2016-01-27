!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2011 all rights reserved
!================================================================
module M_qm_geno_extended
!
! Module variables : changed 
!   --> dhij, dsij, atm_force
! Module variables : unchanged
!   --> all the other variables listed below
!
  use M_qm_domain
  use M_qm_geno_Huckel_atom_params, &
       only: AtomParameter, NUMBER_OF_ORBITALS, NUMBER_OF_CONF_VOIP, &
       UseVOIP, EmulateICON, normDoubleZeta, sdelta=>small_delta, kappa
  use elses_mod_phys_const, only : eV4au
  implicit none
!
   private
!
! Public routines
   public :: set_hamil_and_overlap_geno_ext
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_hamil_and_overlap_geno_ext
!    use M_qm_geno_Huckel, only : SetOverlap, SetDiagonalElements, SetNondiagonalElements
     use M_qm_geno_Huckel, only : SetDiagonalElements, SetNondiagonalElements
     use M_wall_clock_time,only : get_elapse_wall_clock_time
!
     implicit none
     integer, parameter :: ict=1
     integer jsv1, jsv2, nss1, nss2, nval1
     integer jsd1, ja1, ja2, jsd3, jsd4, j1, j2
     real(8) dvecx, dvecy, dvecz, w
     real(8) wclock_time1, wclock_time2
     real(8), allocatable :: S_wrk(:,:), dS_wrk(:,:,:)
     integer ierr
!
     if (i_verbose >= 1) then
       write(*,*)'@@ set_hamiltonian_geno_ext'
       call get_elapse_wall_clock_time(wclock_time1)
       wclock_time2=wclock_time1
     endif
!  
!------------------diagonal part----------
!$omp parallel do default(shared) private(nss1,nss2,jsv1,nval1)
     do jsv2=1,noav
       nss2=atm_element(jsv2)
       do jsd1=1,njsd(jsv2,ict)
         ham_tb0(:,:,jsd1,jsv2)=0.0d0
         dsij   (:,:,jsd1,jsv2)=0.0d0
         jsv1=jsv4jsd(jsd1,jsv2)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
!
         if (jsv1 .eq. jsv2) then
!            ---> on-site terms
            call SetDiagonalElements(nss1,e_num_on_basis(:,jsv2),ham_tb0(:,:,jsd1,jsv2))
            do ja1=1, nval1
               dsij(ja1,ja1,jsd1,jsv2)=1d0
            end do
         endif
!   
       enddo
     enddo
!$omp end parallel do
!------------------end diagonal part----------

     if (i_verbose >= 1) then
       call get_elapse_wall_clock_time(wclock_time1)
       write(*,*)'  time(set_hamiltonian_geno:diagonal) =', wclock_time1-wclock_time2
       wclock_time2=wclock_time1
     endif
!
!$omp parallel do default(shared) &
!$omp&private(nss1, nss2, jsv1, nval1, &
!$omp& dvecx, dvecy, dvecz, w, jsd3, jsd4, S_wrk, dS_wrk)
     do jsv2=1,noav
       nss2=atm_element(jsv2)
       do jsd1=1,njsd(jsv2,ict)
         jsv1=jsv4jsd(jsd1,jsv2)
         nss1=atm_element(jsv1)
!
         if (jsv1 .ne. jsv2) then
!            ---> off-site terms
!
           dvecx=atm_position(1,jsv2)-atm_position(1,jsv1)
           dvecy=atm_position(2,jsv2)-atm_position(2,jsv1)
           dvecz=atm_position(3,jsv2)-atm_position(3,jsv1)
!
           if (i_pbc_x == 1) dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
           if (i_pbc_y == 1) dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
           if (i_pbc_z == 1) dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0
!
           dvecx=dvecx*ax
           dvecy=dvecy*ay
           dvecz=dvecz*az
!             Vector R2 - R1 : (dvecx, dvecy, dvecz)
!
           w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!             Distance | R2 - R1 |
!
           if (allocated(S_wrk)) deallocate(S_wrk)
           allocate(S_wrk(size(dsij,1),size(dsij,2)),stat=ierr)
           if (ierr /= 0) then
             write(*,*)'Alloc. error:set_hamil_and_overlap_geno_ext'
             stop
           endif   
!
           if (allocated(dS_wrk)) deallocate(dS_wrk)
           allocate(dS_wrk( size(ddsij,1),size(ddsij,2),size(ddsij,3) ),stat=ierr)
           if (ierr /= 0) then
             write(*,*)'Alloc. error:set_hamil_and_overlap_geno_ext'
             stop
           endif   
!
           call SetOverlap_wrk(nss1,nss2,dvecx,dvecy,dvecz,S_wrk(:,:),dS_wrk(:,:,:))
           jsd3=jsd_diagonal(jsv1)
           jsd4=jsd_diagonal(jsv2)

           call SetNondiagonalElements(nss1,nss2,dvecx,dvecy,dvecz,&
                S_wrk(:,:), &
                ham_tb0(:,:,jsd3,jsv1),ham_tb0(:,:,jsd4,jsv2), &
                ham_tb0(:,:,jsd1,jsv2), dS_wrk(:,:,:),dham_tb0(:,:,:,jsd1,jsv2))

         endif
!
         if (i_verbose >= 200) then
            write(*,'(A)') '<<<set_hamiltonian_and_overlap_geno verbose mode>>>'
            write(*,'("jsd1=",I6,", jsv2=",I6)') jsd1, jsv2
            write(*,'(A)') 'Overlap matrix'
            do j1=1,size(dsij,1)
               write(*,'(9ES22.14)')(dsij(j1,j2,jsd1,jsv2),j2=1,size(dsij,2))
            end do
         end if
         if (i_verbose >= 100) then
            write(*,'(A)') 'Hamiltonian'
            do j1=1,size(dsij,1)
               write(*,'(9ES22.14)')(ham_tb0(j1,j2,jsd1,jsv2),j2=1,size(dsij,2))
            end do
            write(*,'(A)') '<<<end set_hamiltonian_and_overlap_geno verbose mode>>>'
         end if
!   
       enddo
     enddo
!$omp end parallel do
!
     if (i_verbose >= 15) call browsGreaterElementsOfOverlapMatrix()
     if (i_verbose >= 1) then
       call get_elapse_wall_clock_time(wclock_time1)
       write(*,*)'  time(set_hamiltonian_geno:offdiago) =', wclock_time1-wclock_time2
       wclock_time2=wclock_time1
     endif
!
     contains
       subroutine browsGreaterElementsOfOverlapMatrix()
         use M_lib_sort
         implicit none
         integer, allocatable :: ind(:)
         integer :: N, j1, j2, jsd1, jsv2, k, m, q, mOrb, njsd, noav
         if(size(dsij,1) /= size(dsij,2)) &
              stop 'browsGreaterElementsOfOverlapMatrix: illeagal size of S'
         mOrb=size(dsij,1)
         njsd=size(dsij,3)
         noav=size(dsij,4)
         N=mOrb*mOrb*njsd*noav
         allocate(ind(N))
         call dhsort(N, reshape(dsij,(/ N /) ),ind)
         write(*,'("browsGreaterElementsOfOverlapMatrix:")')
         write(*,'("N=",I9)') N
         write(*,'("S(j1,j2,jsd1,jsv2)")')
         write(*,'("Largest 100 overlaps")')
         
         ! largest element should be diagonal one ( & = 1.d0 )
         do k=N, 1, -1
            m=ind(k)-1
            q=m/mOrb
            j1=m-q*mOrb+1
            m=q
            q=m/mOrb
            j2=m-q*mOrb+1
            m=q
            q=m/njsd
            jsd1=m-q*njsd+1
            m=q
            q=m/noav
            jsv2=m-q*noav+1
            if(dsij(j1,j2,jsd1,jsv2) > 1d0) stop 'anomalous overlap > 1d0'
            if(dsij(j1,j2,jsd1,jsv2) < 1d0) exit ! loop termination
            if( j1 /= j2 )                  stop 'j1 /= j2'
            if( jsd1 /= 1 )                 stop 'not diagonal (jsd1 /= 1)'
         end do
         do k=k, max(k-100,1), -1
            m=ind(k)-1
            q=m/mOrb
            j1=m-q*mOrb+1
            m=q
            q=m/mOrb
            j2=m-q*mOrb+1
            m=q
            q=m/njsd
            jsd1=m-q*njsd+1
            m=q
            q=m/noav
            jsv2=m-q*noav+1
            write(*,'("S(",I1,",",I1,",",I3,",",I8,")=",ES14.7," :jsd1->",I8)')&
                           j1,    j2,   jsd1,  jsv2,  dsij(j1,j2,jsd1,jsv2),&
                           jsv4jsd(jsd1,jsv2)
         end do
         write(*,'("Smallest 100 overlaps")')
         do k=1, min(100,N)
            m=ind(k)-1
            q=m/mOrb
            j1=m-q*mOrb+1
            m=q
            q=m/mOrb
            j2=m-q*mOrb+1
            m=q
            q=m/njsd
            jsd1=m-q*njsd+1
            m=q
            q=m/noav
            jsv2=m-q*noav+1
            write(*,'("S(",I1,",",I1,",",I3,",",I8,")=",ES14.7," :jsd1->",I8)')&
                           j1,    j2,   jsd1,  jsv2,  dsij(j1,j2,jsd1,jsv2),&
                           jsv4jsd(jsd1,jsv2)
         end do
       end subroutine browsGreaterElementsOfOverlapMatrix

       integer function jsd_diagonal(jsv)
         integer,intent(in) :: jsv
         integer :: jsd
!!$         do jsd=1, njsd(jsv,ict)
!!$            if(jsv4jsd(jsd,jsv)==jsv)exit
!!$         end do
!!$         if(jsd>njsd(jsv,ict)) stop 'No diagonal element'
         jsd=1
         if(jsv4jsd(jsd,jsv) /= jsv) stop &
              'Unexpected order of the Hamiltonian data'
         jsd_diagonal=jsd
       end function jsd_diagonal

   end subroutine set_hamil_and_overlap_geno_ext

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine SetOverlap_wrk(left_atom_spec,right_atom_spec,x,y,z,S,dS)
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
    use M_qm_geno_Huckel, only : calc_SS_dSS
    implicit none
    type(TSKBond)            ::  Bond_matrix(NUMBER_OF_ORBITALS,NUMBER_OF_ORBITALS)
    type(TSKBond)            :: dBond_matrix(NUMBER_OF_ORBITALS,NUMBER_OF_ORBITALS)
    type(TVector3d)          :: oresult
    integer, intent(in)      :: left_atom_spec, right_atom_spec
    real(DOUBLE_PRECISION), intent(in) :: x, y, z ! r*(direction cosines)
    real(DOUBLE_PRECISION), intent(out):: S(:,:), dS(:,:,:)
    real(DOUBLE_PRECISION)   :: R, factor, factor1, factor2, factor3, factor4
    real(DOUBLE_PRECISION)   :: factjla, factjlpa, factjlqa
    real(DOUBLE_PRECISION)   :: factjlb, factjlpb, factjlqb
    real(DOUBLE_PRECISION)   :: factmk, factmk_p
    real(DOUBLE_PRECISION)   :: AA(0:14), BB(0:14), dAA(0:14), dBB(0:14)
    real(DOUBLE_PRECISION)   :: SS(9,9), dSS(9,9)
    real(DOUBLE_PRECISION)   :: dSSV(3,9,9), sresult, deriv_mat, sign_factor
    real(DOUBLE_PRECISION)   :: dfactor1, dfactor4, factor5
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
    integer                  :: l_a, l_b
    integer                  :: npq_a, npq_b
    integer                  :: j_a, ip_a, iq_a, k, mm
    integer                  :: j_a_max, ip_a_max, iq_a_max
    integer                  :: j_b, ip_b, iq_b, k_prime
    integer                  :: j_b_max, ip_b_max, iq_b_max
    character(len=8)         :: orbital_name(9)
!
    real(8) weight_factor
!
    weight_factor=0.0d0
!   weight_factor=1.0d0
!
    ! Factorials are prepared and stored.
    ! fact(n)=n*(n-1)*(n-2)....2*1
    ! Definition of factorial is different from ICON-EDIT.
    ! ICON-EDIT: maximum of fact(n) is fact(25).
    ! added by m.ikeda, 2008/08/04, FLA, Japan.

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
!       Vpp_pi=SS(4,4)
!       dVpp_pi=dSS(4,4)
        Vpp_pi=SS(4,4)*weight_factor
        dVpp_pi=dSS(4,4)*weight_factor
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
  end subroutine SetOverlap_wrk

end module M_qm_geno_extended
