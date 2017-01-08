!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno
!
! Module variables : changed 
!   --> dhij, dsij, atm_force
! Module variables : unchanged
!   --> all the other variables listed below
!
  use M_qm_domain
  implicit none
!
   private
!
! Public routines
   public :: set_hamiltonian_and_overlap_geno
   public :: set_rest_geno, set_qm_force_geno
   public :: initialize_charge,  renew_charge, rms_delta_q
   public :: write_hamiltonian_and_overlap_geno, write_Ecc_perAtom_geno
   public :: calc_e_rep_atom
!
   contains

     subroutine initialize_charge
       use M_qm_geno_Huckel_atom_params, only: GetInitialOcc, GetInitialENum
       call GetInitialOcc (e_num_on_basis)
       call GetInitialENum(e_num_on_atom)
       previous_e_num_on_basis=e_num_on_basis
       previous_e_num_on_atom =e_num_on_atom
     end subroutine initialize_charge

     subroutine renew_charge(mixing_ratio)
       use M_qm_geno_Huckel_atom_params, only: GetInitialOcc, GetInitialENum
       real(DOUBLE_PRECISION), intent(in)  :: mixing_ratio
       real(DOUBLE_PRECISION), allocatable :: delta(:,:)
       integer :: jsv2, nval2
       allocate(delta(size(e_num_on_basis,1),size(e_num_on_basis,2)))
       delta=e_num_on_basis-previous_e_num_on_basis
       e_num_on_basis = previous_e_num_on_basis + delta * mixing_ratio
       previous_e_num_on_basis = e_num_on_basis
       previous_e_num_on_atom  = e_num_on_atom
       ! renew e_num_on_atom
       do jsv2=1,noav
          nval2=nval(atm_element(jsv2))
          e_num_on_atom(jsv2)=sum(e_num_on_basis(1:nval2,jsv2))
       end do
       if( i_verbose >= 50 ) then
          write(*,'("renew_charge:")')
          write(*,'(" delta = ", ES21.14)') delta
          write(*,'(" new e_num_on_atom = ", ES21.14)') e_num_on_atom
       end if
       deallocate(delta)
     end subroutine renew_charge

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_hamiltonian_and_overlap_geno(bondMask)
     use M_config               ! (unchanged)
     use M_qm_geno_Huckel, only : SetOverlap, SetDiagonalElements,&
     SetNondiagonalElements, SetDiagonalShifts
     use M_wall_clock_time,only : get_elapse_wall_clock_time
!
     implicit none
     integer, parameter :: ict=1
     integer            :: jsv1, jsv2, nss1, nss2, nval1
     integer            :: jsd1, ja1, ja2, jsd3, jsd4, j1, j2
     real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz, w
     real(DOUBLE_PRECISION) :: wclock_time1, wclock_time2
     real(DOUBLE_PRECISION),optional   :: bondMask(:)
!
     if (i_verbose >= 1) then
       write(*,*)'@@ set_hamiltonian_geno'
       call get_elapse_wall_clock_time(wclock_time1)
       wclock_time2=wclock_time1
     endif
!  
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
!$omp& dvecx, dvecy, dvecz, w, jsd3, jsd4)
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
           if(present(bondMask))then
              call SetOverlap(nss1,nss2,dvecx,dvecy,dvecz,dsij(:,:,jsd1,jsv2),ddsij(:,:,:,jsd1,jsv2),bondMask)
           else
              call SetOverlap(nss1,nss2,dvecx,dvecy,dvecz,dsij(:,:,jsd1,jsv2),ddsij(:,:,:,jsd1,jsv2))
           end if
           jsd3=jsd_diagonal(jsv1)
           jsd4=jsd_diagonal(jsv2)

           call SetNondiagonalElements(nss1,nss2,dvecx,dvecy,dvecz,&
                dsij(:,:,jsd1,jsv2), &
                ham_tb0(:,:,jsd3,jsv1),ham_tb0(:,:,jsd4,jsv2), &
                ham_tb0(:,:,jsd1,jsv2), ddsij(:,:,:,jsd1,jsv2),dham_tb0(:,:,:,jsd1,jsv2))

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
!
!------------------extra diagonal shifts----------
!$omp parallel do default(shared) private(nss1,nss2,jsv1,nval1)
     do jsv2=1,noav
       nss2=atm_element(jsv2)
       do jsd1=1,njsd(jsv2,ict)
         jsv1=jsv4jsd(jsd1,jsv2)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
!
         if (jsv1 .eq. jsv2) then
!            ---> on-site terms
            call SetDiagonalShifts(nss1,ham_tb0(:,:,jsd1,jsv2))
         endif
!   
       enddo
     enddo
!$omp end parallel do
!------------------end extra diagonal shifts----------

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

   end subroutine set_hamiltonian_and_overlap_geno

   subroutine write_hamiltonian_and_overlap_geno(fileName,BinaryIO)
     use elses_mod_file_io
     integer, parameter :: ict=1
     integer :: jsv1, jsv2, njsd1, njsd2, nss1, nss2, nval1, nval2
     integer :: jsd1, ja1, ja2, jsd3, jsd4, j1, j2
     integer :: iUnit
     character(len=*) :: fileName
     logical :: BinaryIO

     iUNIT=vacant_unit()

     if(.not. BinaryIO)then
     open(iUnit, file=fileName, status='replace')
     write(iUnit,'(A)') &
          "# size of Overlap & Hamiltonian storage&
          & S(int1,int2,int3,int4), H(int1,int2,int3,int4)"
     write(iUnit,'(4(I8," "))')size(dsij,1),size(dsij,2),size(dsij,3),size(dsij,4)
     write(iUnit,'(A)') &
          "# List of neighbor atoms (int3-th neighbor of int4)"
     do jsv2=1, noav
        njsd2=njsd(jsv2,ict)
        write(iUnit,'(I8,"(=int4)-th atom : ",I8," neighbors")') jsv2, njsd2
        !read(iUnit,'(I8,A18,I8)') jsv2, string, njsd2
        write(iUnit,'(8(I8," "))')  (jsv4jsd(jsd1,jsv2) ,jsd1=1, njsd2)
     end do
     write(iUnit,'(A)') &
          "# Number of valence orbitals of int4-th atom"
     write(iUnit,'(8(I8," "))')  (nval(atm_element(jsv2)),jsv2=1,noav)
     write(iUnit,'(A)') &
          "# S(int1,int2,int3,int4) = double1, H(int1,int2,int3,int4) = double2"
     do jsv2=1,noav
        njsd2=njsd(jsv2,ict)
        nss2=atm_element(jsv2)
        nval2=nval(nss2)
        do jsd1=1,njsd2
           jsv1=jsv4jsd(jsd1,jsv2)
           nss1=atm_element(jsv1)
           nval1=nval(nss1)
           do j2=1, nval2
              do j1=1, nval1
                 write(iUnit,'(4(I8," "),2(ES23.16," "))')&
                      j1,j2,jsd1,jsv2,&
                      dsij(j1,j2,jsd1,jsv2), ham_tb0(j1,j2,jsd1,jsv2)
              end do
           end do
           select case (i_verbose)
           case (1:199)
              write(*,'("write_hamiltonian_and_overlap_geno:")')
              write(*,'(" verbose =",I4)')i_verbose
           case (200:)
              write(*,'(A)') '<<<set_hamiltonian_and_overlap_geno verbose mode>>>'
              write(*,'("jsd1=",I6,", jsv2=",I6)') jsd1, jsv2
              write(*,'(A)') 'Overlap matrix'
              do j1=1,size(dsij,1)
                 write(*,'(9ES22.14)')(dsij(j1,j2,jsd1,jsv2),j2=1,size(dsij,2))
              end do
              write(*,'(A)') 'Hamiltonian'
              do j1=1,size(dsij,1)
                 write(*,'(9ES22.14)')(ham_tb0(j1,j2,jsd1,jsv2),j2=1,size(dsij,2))
              end do
              write(*,'(A)') '<<<end set_hamiltonian_and_overlap_geno verbose mode>>>'
           case default
              ! do nothing
           end select
!   
        end do
     end do
     else
        open(iUnit, file=fileName, status='replace',FORM='UNFORMATTED')
        ! short cut for large (noav>=100) system        
        write(iUnit) size(dsij,1),size(dsij,2),size(dsij,3),size(dsij,4)
        write(iUnit) njsd(:,ict) !!! The range of 2nd argument of njsd : 0:nncut !!!
        write(iUnit) jsv4jsd
        write(iUnit) (nval(atm_element(jsv2)),jsv2=1,noav)
        write(iUnit) dsij
        write(iUnit) ham_tb0
     end if
     close(iUnit)
   end subroutine write_hamiltonian_and_overlap_geno

   subroutine write_Ecc_perAtom_geno(EccPerAtom,fileName,BinaryIO)
     use elses_mod_file_io
     integer, parameter :: ict=1
     integer :: jsv1, jsv2, njsd1, njsd2, nss1, nss2, nval1, nval2
     integer :: jsd1, ja1, ja2, jsd3, jsd4, j1, j2
     integer :: iUnit
     character(len=*) :: fileName
     logical :: BinaryIO
     real(kind=DOUBLE_PRECISION):: EccPerAtom(:)

     if(size(EccPerAtom,1) /= noav) then
        if(i_verbose >= 1) then
           write(*,*) '@@write_Ecc_perAtom_geno: size mismatch. size(EccPerAtom,1) /= noav'
           stop
        end if
     end if
     iUNIT=vacant_unit()
     write(*,*) fileName

     if(.not. BinaryIO)then
     open(iUnit, file=fileName, status='replace')
     write(iUnit,'(A)') &
          "# of atoms"
     write(iUnit,'(I8)') noav
     write(iUnit,'(A)') &
          "# Ecc/Atom"
     do jsv2=1,noav
        write(iUnit,'(ES22.14)') EccPerAtom(jsv2)
     end do
     else
        open(iUnit, file=fileName, status='replace',FORM='UNFORMATTED')
        ! short cut for large (noav>=100) system        
        write(iUnit) noav
        write(iUnit) EccPerAtom
     end if
     close(iUnit)
   end subroutine write_Ecc_perAtom_geno


!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_rest_geno(EccPerAtom)
!
     use M_config ! unchanged            
     use M_qm_geno_Huckel, only : RepulsiveEnergy ! function
!
     use elses_mod_ene, only : ecc           ! CHANGED
     use M_qm_domain,   only : E_vdW_only    ! CHANGED
!
     use elses_mod_phys_const, only : ev4au  ! parameter
     use M_qm_geno_rest_kernel, only : calc_geno_rest_kernel !routine
!
     implicit none
     integer, parameter :: ict=1
     integer            :: jsv1, jsv2, njsd1, njsd2, nss1, nss2, jsd1
     real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz, w, w_vdW
     real(DOUBLE_PRECISION) :: repulsive_force(3)
     real(DOUBLE_PRECISION),optional :: EccPerAtom(:)
     real(DOUBLE_PRECISION) :: E_atompair_vdW_only
     integer                :: lu, ierr
     integer, allocatable   :: local_booking_list(:)
     integer                :: len_local_booking_list
     logical                :: cutoff_rest_cell_max
!
     cutoff_rest_cell_max = config%calc%interaction_range%cutoff_rest_cellmax
     lu=config%calc%distributed%log_unit
!
     allocate (local_booking_list(noav), stat=ierr)
     if ( ierr /= 0 ) then 
       stop 'Alloc. Error (set_hamiltonian_and_overlap_geno:local_booking_list)'
     endif  
!
     if (i_verbose >= 1) then
       if (lu > 0) write(lu,*)'@@ set_rest_geno'
     endif
     if(present(EccPerAtom))then
        if (i_verbose >= 1) then
          if (lu > 0) write(lu,*)'   set_rest_geno: optionally calculate Ecc/Atom.'
        end if
        if(size(EccPerAtom,1) /= noav)then
           if (i_verbose >= 1) then
             if (lu > 0) write(lu,*)'   set_rest_geno: size mismatch.# of atom /= size(EccPerAtom,1)'
           end if
        end if
        EccPerAtom=0d0
     end if
!  
     ecc=0.0d0
     E_vdW_only = 0.0d0
     atm_force=0.0d0
!
!
     do jsv2=1,noav
       nss2=atm_element(jsv2)
!
       if (cutoff_rest_cell_max) then
         len_local_booking_list      = noav
         do jsv1=1,noav
           local_booking_list(jsv1) = jsv1
         enddo   
       else
         njsd2=njsd(jsv2,ict)
         len_local_booking_list=njsd2
         do jsd1=1,njsd2
           jsv1=jsv4jsd(jsd1,jsv2)
           local_booking_list(jsd1)=jsv1
         enddo
       endif   
!
       call calc_geno_rest_kernel(jsv2, local_booking_list, len_local_booking_list, w, w_vdw, atm_force)
!
       ecc = ecc + w
       E_vdW_only = E_vdW_only + w_vdW
       if(present(EccPerAtom)) EccPerAtom(jsv2)=w
!
     enddo
!
     deallocate (local_booking_list, stat=ierr)
     if ( ierr /= 0 ) then 
       stop 'Dealloc. Error (set_hamiltonian_and_overlap_geno:local_booking_list)'
     endif  
!
     if (i_verbose >= 1) then
       if (lu > 0) write(lu,*)'INFO:E_vdW_only [au, eV] =', E_vdW_only, E_vdW_only*ev4au
       if (config%calc%genoOption%vanDerWaals) then
         if (config%calc%genoOption%vanDerWaals_lambda > 0.0d0) then
           if (lu >0) write(lu,*)'INFO:vdW_lambda      =', config%calc%genoOption%vanDerWaals_lambda
         endif   
       endif   
     endif
!
     if (i_verbose >= 1) then
       if (lu >0) write(*,'("ecc [au], ecc[eV]=",ES22.14,",",ES22.14)') ecc, ecc*ev4au
     endif
!
     if (i_verbose >= 50) then
       do jsv1=1,noav
         write(*,'("Atom=",I6," :Force from repulsive part =",3ES23.14)') jsv1, atm_force(:,jsv1)
       end do
     endif
!
   end subroutine set_rest_geno


   subroutine set_qm_force_geno
     use elses_mod_phys_const, only : para_spin_factor
     implicit none
     integer, parameter :: ict=1
     integer :: jsv1, jsv2, njsd1, njsd2, nval1, nval2
     integer :: jsd1, jsd2, nss1, nss2, ja1, ja2
     real(DOUBLE_PRECISION),allocatable :: atm_force_tmp(:,:)
     allocate(atm_force_tmp(3,noav))
     atm_force_tmp=0.0d0
     do jsv2=1,noav
       njsd2=njsd(jsv2,ict)
       nss2=atm_element(jsv2)
       nval2=nval(nss2)
       do jsd1=1,njsd2
         jsv1=jsv4jsd(jsd1,jsv2)
!
         njsd1=njsd(jsv1,ict)
         do jsd2=1, njsd1
            if(jsv2 == jsv4jsd(jsd2,jsv1)) exit
         end do
         if( jsd2 == njsd1+1 ) stop 'set_qm_force_geno: no pair'
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         if (jsv1 /=  jsv2) then
            do ja2=1, nval2
               do ja1=1, nval1
                  atm_force_tmp(:,jsv2)=atm_force_tmp(:,jsv2)+&
                       ddhij(:,ja1,ja2,jsd1,jsv2)*dbij(ja1,ja2,jsd1,jsv2)-&
                       ddsij(:,ja1,ja2,jsd1,jsv2)*dpij(ja1,ja2,jsd1,jsv2)
                  atm_force_tmp(:,jsv1)=atm_force_tmp(:,jsv1)+&
                       ddhij(:,ja2,ja1,jsd2,jsv1)*dbij(ja1,ja2,jsd1,jsv2)-&
                       ddsij(:,ja2,ja1,jsd2,jsv1)*dpij(ja1,ja2,jsd1,jsv2)
               end do
            end do
         endif
       enddo
     enddo
     atm_force_tmp  = para_spin_factor * atm_force_tmp
     atm_force(:,:) = atm_force(:,:) - atm_force_tmp(:,:) ! F_{REP} + F_{TB}
!
     if (i_verbose >= 50) then
        write(*,'("set_qm_force_geno:")')
        do ja1=1, noav
           write(*,'("atom=",I8)')ja1
           write(*,'("    TB force (:,",I8,")",3ES22.14)')  ja1, (atm_force_tmp(ja2,ja1),ja2=1, 3)
        end do
     end if
     if (i_verbose >= 1) then
        do ja1=2, noav
           atm_force_tmp(:,1)=atm_force_tmp(:,1)+atm_force_tmp(:,ja1)
        end do
        write(*,'(" Sum of TB force =",3ES22.14)') -atm_force_tmp(:,1)
        atm_force_tmp=0d0
        do ja1=1, noav
           atm_force_tmp(:,1)=atm_force_tmp(:,1)+atm_force(:,ja1)
        end do
        write(*,'(" Sum of    force =",3ES22.14)')  atm_force_tmp(:,1)
     endif
!
     deallocate(atm_force_tmp)
   end subroutine set_qm_force_geno
!

   function rms_delta_q() result(dq)

     implicit none
     real(8) :: dq, nf
     integer :: i

     dq=0.0d0
     nf=0.0d0

     do i=1,noav
        dq=dq+(e_num_on_atom(i)-previous_e_num_on_atom(i))**2
        nf=nf+(previous_e_num_on_atom(i))**2
     end do

     dq=sqrt(dq/nf)
     
     if(i_verbose >= 200) then
        write(*,'("normalized dq=",F16.10)') dq
     end if
     
   end function rms_delta_q
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the repulsive energy for a given atom
   subroutine calc_e_rep_atom(atm_index, value_of_energy)
     use M_qm_geno_Huckel, only : RepulsiveEnergy
!
!
     implicit none
     integer,                intent(in)  :: atm_index
     real(DOUBLE_PRECISION), intent(out) :: value_of_energy
     integer, parameter :: ict=1
     integer            :: jsv1, jsv2, njsd1, njsd2, nss1, nss2, jsd1
     real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz
     real(DOUBLE_PRECISION) :: repulsive_force(3)
!
     value_of_energy=0.0d0
!
     jsv2=atm_index
     njsd2=njsd(jsv2,ict)
     nss2=atm_element(jsv2)
     do jsd1=1,njsd2
       jsv1=jsv4jsd(jsd1,jsv2)
       nss1=atm_element(jsv1)
!
       if (jsv1 /=  jsv2) then
         dvecx=atm_position(1,jsv2)-atm_position(1,jsv1)
         dvecy=atm_position(2,jsv2)-atm_position(2,jsv1)
         dvecz=atm_position(3,jsv2)-atm_position(3,jsv1)
         if (i_pbc_x == 1) dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
         if (i_pbc_y == 1) dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
         if (i_pbc_z == 1) dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0
         dvecx=dvecx*ax
         dvecy=dvecy*ay
         dvecz=dvecz*az
!             Vector R2 - R1 : (dvecx, dvecy, dvecz)
         value_of_energy=value_of_energy+RepulsiveEnergy(nss1, nss2,dvecx,dvecy,dvecz, repulsive_force)
       endif
     enddo
!
   end subroutine calc_e_rep_atom
!
end module M_qm_geno
