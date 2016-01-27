!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_geno_dst
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
  real(DOUBLE_PRECISION), allocatable :: mat_tb0_diag(:,:) !(calculated in set_mat_tb0_diag)
!     -----> Diagonal elements of the TB0 matrix
!
  real(DOUBLE_PRECISION), allocatable :: ham_tot_dst(:,:,:,:)     !(calculated in set_ham_tb0_and_overlap_dst0)
  real(DOUBLE_PRECISION), allocatable :: ham_tb0_dst(:,:,:,:)     !(calculated in set_ham_tb0_and_overlap_dst0)
  real(DOUBLE_PRECISION), allocatable :: d_ham_tb0_dst(:,:,:,:,:) !(calculated in set_ham_tb0_and_overlap_dst0)
  real(DOUBLE_PRECISION), allocatable :: overlap_dst(:,:,:,:)     !(calculated in set_ham_tb0_and_overlap_dst0)
  real(DOUBLE_PRECISION), allocatable :: d_overlap_dst(:,:,:,:,:) !(calculated in set_ham_tb0_and_overlap_dst0)
!
   private
!
! Public variables
   public :: ham_tot_dst, ham_tb0_dst, d_ham_tb0_dst, overlap_dst, d_overlap_dst
!
! Public routines
   public :: set_ham_tb0_and_overlap_dst
!

!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_ham_tb0_and_overlap_dst(bondMask)
!
     use M_qm_domain, only : atm_element, noav, nval, jsv4jsd, njsd, e_num_on_basis !(unchanged)
     use M_qm_domain, only : atm_position, i_pbc_x, i_pbc_y, i_pbc_z, ax, ay, az !(unchanged)
     use M_qm_geno_Huckel, only : SetOverlap, SetDiagonalElements, SetNonDiagonalElements, SetDiagonalShifts !(routines)
     use M_qm_dst_proj_cell, only : dst_atm_list, len_dst_atm_list !(unchanged)
     use elses_mod_noav,   only : noas !(unchanged)
     implicit none
     real(DOUBLE_PRECISION),optional   :: bondMask(:)
!
     real(DOUBLE_PRECISION), allocatable :: ham_tb0_onsite(:,:,:) 
     real(DOUBLE_PRECISION), allocatable :: ham_tb0_offsite(:,:) 
     real(DOUBLE_PRECISION), allocatable :: d_ham_tb0_offsite(:,:,:) 
     real(DOUBLE_PRECISION), allocatable :: overlap_offsite(:,:) 
     real(DOUBLE_PRECISION), allocatable :: d_overlap_offsite(:,:,:) 
     integer, parameter :: ict=1
     integer            :: ierr, nval_max
     integer            :: atm_index, dst_atm_index, orb_index
     integer            :: jsv2, nss2, jsd1, jsv1, nss1, jsd3, jsd4
     real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz, w
     logical :: set_ham_tb0_dst
     logical :: set_overlap_dst
     integer            :: len_dst_atm_list_alloc

!
     set_ham_tb0_dst=.true.
     set_overlap_dst=.true.
!
     len_dst_atm_list_alloc=min(noav,nint(1.5d0*len_dst_atm_list(2)))
!        ---> The factor of 1.5 is a tolerance factor.
!
     if (i_verbose >= 1) then
       write(*,*)'@@ set_ham_tb0_and_overlap_dst'
       write(*,*)'  noas =',noas
       write(*,*)'  len_dst_atm_list(1)(2) =',len_dst_atm_list(1), len_dst_atm_list(2)
       write(*,*)'  len_dst_atm_list_alloc =',len_dst_atm_list_alloc
       write(*,*)'  size(dst_atm_list,1)   =',size(dst_atm_list,1)
     endif
!  
     if (size(dst_atm_list,1) < len_dst_atm_list(2)) then
       write(*,*)'ERROR(set_ham_tb0_and_overlap_dst)'
       write(*,*)'  len_dst_atm_list(1) =',len_dst_atm_list(1)
       write(*,*)'  len_dst_atm_list(2) =',len_dst_atm_list(2)
       write(*,*)'  size(dst_atm_list,1)   =',size(dst_atm_list,1)
       stop
     endif
!   
     call set_mat_tb0_diag
!
     nval_max=maxval(nval)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of ham_tot_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     if ( allocated(ham_tot_dst)) then
       if (size(ham_tot_dst,4) < len_dst_atm_list(2)) then
         deallocate ( ham_tot_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (ham_tot_dst)'         
       endif
     endif
!
     if ( .not. allocated(ham_tot_dst)) then
       allocate ( ham_tot_dst( nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (ham_tot_dst)'
     endif  
     ham_tot_dst(:,:,:,:)=0.0d0
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of ham_tb0_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     len_dst_atm_list_alloc=size(ham_tot_dst,4)
     if ( allocated(ham_tb0_dst)) then
       if (size(ham_tb0_dst,4) /= len_dst_atm_list_alloc) then
         deallocate ( ham_tb0_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (ham_tb0_dst)'
       endif
     endif
!
     if ( .not. allocated(ham_tb0_dst)) then
       allocate ( ham_tb0_dst( nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (ham_tb0_dst)'
     endif  
     ham_tb0_dst(:,:,:,:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of overlap_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     len_dst_atm_list_alloc=size(ham_tot_dst,4)
     if ( allocated(overlap_dst)) then
       if (size(overlap_dst,4) /= len_dst_atm_list_alloc) then
         deallocate ( overlap_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (overlap_dst)'
       endif
     endif
!
     if ( .not. allocated(overlap_dst)) then
       allocate ( overlap_dst( nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (ham_tb0_dst)'
     endif  
     overlap_dst(:,:,:,:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of d_ham_tb0_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     len_dst_atm_list_alloc=size(ham_tot_dst,4)
     if ( allocated(d_ham_tb0_dst)) then
       if (size(d_ham_tb0_dst,5) /= len_dst_atm_list_alloc) then
         deallocate ( d_ham_tb0_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (d_ham_tb0_dst)'
       endif
     endif
!     
     if (.not. allocated(d_ham_tb0_dst)) then
       allocate ( d_ham_tb0_dst(3, nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (d_ham_tb0_dst)'
     endif
     d_ham_tb0_dst(:,:,:,:,:)=0.0d0
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Allocation of d_overlap_dst, if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     len_dst_atm_list_alloc=size(ham_tot_dst,4)
     if ( allocated(d_overlap_dst)) then
       if (size(d_overlap_dst,5) /= len_dst_atm_list_alloc) then
         deallocate ( d_overlap_dst,stat=ierr)
         if( ierr /= 0) stop 'ERROR in dealloc (d_ham_tb0_dst)'
       endif
     endif
!
     if (.not. allocated(d_overlap_dst)) then
       allocate ( d_overlap_dst( 3, nval_max, nval_max, noas, len_dst_atm_list_alloc ),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (d_overlap_dst)'
     endif  
     d_overlap_dst(:,:,:,:,:)=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!$omp  parallel default(shared) &
!$omp& private(ham_tb0_onsite, ham_tb0_offsite, d_ham_tb0_offsite, overlap_offsite, d_overlap_offsite) &
!$omp& private(dst_atm_index, atm_index, jsv2, nss2, jsd1, jsv1, nss1, orb_index) &
!$omp& private(dvecx, dvecy, dvecz, w)
!
     allocate (ham_tb0_onsite(nval_max, nval_max,2),stat=ierr)
     if( ierr /= 0) stop 'ERROR in alloc (ham_tb_onsite)'
!
     allocate (ham_tb0_offsite(nval_max, nval_max),stat=ierr)
     if( ierr /= 0) stop 'ERROR in alloc (ham_tb_offsite)'
!
     allocate (d_ham_tb0_offsite(3,nval_max, nval_max),stat=ierr)
     if( ierr /= 0) stop 'ERROR in alloc (d_ham_tb_offsite)'
!
     allocate (overlap_offsite(nval_max, nval_max),stat=ierr)
     if( ierr /= 0) stop 'ERROR in alloc (overlap_offsite)'
!
     allocate (d_overlap_offsite(3, nval_max, nval_max),stat=ierr)
     if( ierr /= 0) stop 'ERROR in alloc (d_overlap_offsite)'
!
!$omp  do schedule(static)
     do dst_atm_index=1,len_dst_atm_list(2)
      atm_index=dst_atm_list(dst_atm_index)
      jsv2=atm_index
      nss2=atm_element(jsv2)
      do jsd1=1,njsd(jsv2,ict)
        jsv1=jsv4jsd(jsd1,jsv2)
        nss1=atm_element(jsv1)
!
        if (jsv1 == jsv2) then
          if (set_ham_tb0_dst) then
            do orb_index=1,nval_max
              ham_tb0_dst(orb_index,orb_index,jsd1,dst_atm_index)=mat_tb0_diag(orb_index,jsv2)
            enddo
          endif   
          if (set_overlap_dst) then
            do orb_index=1,nval_max
              overlap_dst(orb_index,orb_index,jsd1,dst_atm_index)=1.0d0
            enddo
          endif   
          cycle
        endif  
!
        ham_tb0_onsite(:,:,:)=0.0d0
        ham_tb0_offsite(:,:)=0.0d0
        d_ham_tb0_offsite(:,:,:)=0.0d0
        overlap_offsite(:,:)=0.0d0
        d_overlap_offsite(:,:,:)=0.0d0
!
        do orb_index=1,maxval(nval)
          ham_tb0_onsite(orb_index,orb_index,1)=mat_tb0_diag(orb_index,jsv1)
          ham_tb0_onsite(orb_index,orb_index,2)=mat_tb0_diag(orb_index,jsv2)
        enddo   
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
           call SetOverlap(nss1,nss2,dvecx,dvecy,dvecz,overlap_offsite(:,:),d_overlap_offsite(:,:,:),bondMask)
        else
           call SetOverlap(nss1,nss2,dvecx,dvecy,dvecz,overlap_offsite(:,:),d_overlap_offsite(:,:,:))
        end if
!
        call SetNondiagonalElements(nss1,nss2,dvecx,dvecy,dvecz,&
&             overlap_offsite(:,:), &
&             ham_tb0_onsite(:,:,1), ham_tb0_onsite(:,:,2), &
&             ham_tb0_offsite(:,:), d_overlap_offsite(:,:,:), d_ham_tb0_offsite(:,:,:))
!
        if (set_ham_tb0_dst) then
            ham_tb0_dst(:,:,  jsd1,dst_atm_index) =   ham_tb0_offsite(:,:)
          d_ham_tb0_dst(:,:,:,jsd1,dst_atm_index) = d_ham_tb0_offsite(:,:,:)
        endif   
!
        if (set_overlap_dst) then
            overlap_dst(:,:,  jsd1,dst_atm_index) =   overlap_offsite(:,:)
          d_overlap_dst(:,:,:,jsd1,dst_atm_index) = d_overlap_offsite(:,:,:)
        endif   
!
      enddo
     enddo
!$omp end do
!
     ham_tot_dst(:,:,:,:)=ham_tb0_dst(:,:,:,:)
!
     deallocate (ham_tb0_onsite, stat=ierr)
     if( ierr /= 0) stop 'ERROR in dealloc (ham_tb_onsite)'
!
     deallocate (ham_tb0_offsite, stat=ierr)
     if( ierr /= 0) stop 'ERROR in dealloc (ham_tb_offsite)'
!
     deallocate (d_ham_tb0_offsite, stat=ierr)
     if( ierr /= 0) stop 'ERROR in dealloc (d_ham_tb_offsite)'
!
     deallocate (overlap_offsite, stat=ierr)
     if( ierr /= 0) stop 'ERROR in dealloc (overlap_offsite)'
!
     deallocate (d_overlap_offsite, stat=ierr)
     if( ierr /= 0) stop 'ERROR in dealloc (d_overlap_offsite)'
!
!$omp end parallel
!
   end subroutine set_ham_tb0_and_overlap_dst
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_mat_tb0_diag
!      OUTPUT : mat_tb0_diag(:,:)
!
     use M_qm_domain, only : atm_element, noav, nval, jsv4jsd, e_num_on_basis !(unchanged)
     use M_qm_geno_Huckel, only : SetDiagonalElements,SetDiagonalShifts !(routines)
     implicit none
     real(DOUBLE_PRECISION), allocatable :: wrk_mat(:,:) ! work array
     integer            :: ierr
     integer            :: jsv, nss1, jsd, orb_index
!
     if (i_verbose >= 1) then
       write(*,*)'@@ set_hamil_mat_diag'
     endif
!  
     if (.not. allocated(mat_tb0_diag)) then
       allocate (mat_tb0_diag(maxval(nval), noav),stat=ierr)
       if( ierr /= 0) stop 'ERROR in alloc (mat_tb0_diag)'
     endif
!   
!
!$omp  parallel default(shared) &
!$omp& private(jsv, nss1, jsd, wrk_mat,orb_index) 
!
     allocate (wrk_mat(maxval(nval), maxval(nval)),stat=ierr)
     if( ierr /= 0) stop 'ERROR in alloc (wrk_mat)'
!
!$omp  do schedule(static)
     do jsv=1,noav
!      write(*,*)'jsv(set_mat_tb0_diag)=',jsv
       nss1=atm_element(jsv)
       jsd=1
       if(jsv4jsd(jsd,jsv) /= jsv) stop 'Unexpected order of the Hamiltonian data'
       wrk_mat(:,:)=0.0d0
!      write(*,*)'size(wrk_mat,1)=',size(wrk_mat,1)
!      write(*,*)'size(wrk_mat,2)=',size(wrk_mat,2)
!
!      write(*,*)'Go Diag'
       call SetDiagonalElements(nss1,e_num_on_basis(:,jsv),wrk_mat(:,:))
!      write(*,*)'Go DiagShift'
       call SetDiagonalShifts(nss1,wrk_mat(:,:))
!
       mat_tb0_diag(:,jsv)=0.0d0
!      write(*,*)'Go orb_index_loop'
       do orb_index=1,nval(nss1)
!        write(*,*)'orb_index(set_mat_tb0_diag)=',orb_index
         mat_tb0_diag(orb_index,jsv)=wrk_mat(orb_index,orb_index)
       enddo
!
     enddo
!$omp end do
!
     deallocate (wrk_mat,stat=ierr)
     if( ierr /= 0) stop 'ERROR in dealloc (wrk_mat)'
!
!$omp end parallel
!
!
   end subroutine set_mat_tb0_diag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_geno_dst
