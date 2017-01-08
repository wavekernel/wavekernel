!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_dens_mat
!
!
  use M_qm_domain,     only : i_verbose !(unchanged)
!
  implicit none
  real(8), allocatable ::   dens_mat_dst(:,:,:,:) ! DIFFERENT VALUES AMONG NODES
  real(8), allocatable :: e_dens_mat_dst(:,:,:,:) ! DIFFERENT VALUES AMONG NODES
  real(8), allocatable :: e_num_on_atom_dst(:)    ! DIFFERENT VALUES AMONG NODES !(calculated in mulliken_dst_test)
  real(8), allocatable :: e_num_on_basis_dst(:,:) ! DIFFERENT VALUES AMONG NODES !(calculated in mulliken_dst_test)
!  
  private
!
! Public variables
  public dens_mat_dst, e_dens_mat_dst, e_num_on_atom_dst, e_num_on_basis_dst
!
! Public routines
  public alloc_dst_dens_mat
  public qm_calc_etb_dst
! public qm_calc_etb_dst_test
  public set_tb_force_dst
  public calc_mulliken_dst
  public calc_mulliken_dst_test

!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Allocate the (energy) density matrix for the distributed atoms
!       ----> dens_mat_dst, e_dens_mat_dst
!
  subroutine alloc_dst_dens_mat
!
    use M_qm_domain,        only : noav, nval       !(unchanged)
    use M_qm_dst_proj_cell, only : len_dst_atm_list !(routine)
    use M_qm_geno_dst,      only : ham_tot_dst      !(unchanged)  
!
    implicit none
    integer :: alloc_size_dst
!
    integer :: nval_max, matsize, ierr 
!
    if (i_verbose >= 1) then
      write(*,*)'@@ alloc_dst_dens_mat'
    endif
!  
    alloc_size_dst=size(ham_tot_dst,3)
    nval_max=maxval(nval)
!
    if (allocated(dens_mat_dst)) then
      matsize=size(dens_mat_dst,3)
      if (matsize >= alloc_size_dst) then
        write(*,*)' (e_)dens_mat_dst : already allocated and OK'
        return
      else
        deallocate (dens_mat_dst, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (dens_mat_dst)'
        deallocate (e_dens_mat_dst, stat=ierr)
        if (ierr /= 0) stop 'Abort:ERROR in dealloc (e_dens_mat_dst)'
      endif   
    endif   
!
    allocate (dens_mat_dst(nval_max, nval_max, alloc_size_dst, len_dst_atm_list(1)), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (dens_mat_dst)'
    dens_mat_dst(:,:,:,:)=0.0d0
!                                                                                                              
    allocate (e_dens_mat_dst(nval_max, nval_max, alloc_size_dst, len_dst_atm_list(1)), stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in alloc (e_dens_mat_dst)'
    e_dens_mat_dst(:,:,:,:)=0.0d0
!
  end subroutine alloc_dst_dens_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Deallocate the (energy) density matrix for the distributed atoms
!       ----> dens_mat_dst, e_dens_mat_dst
!
  subroutine dealloc_dst_dens_mat
!
    implicit none
    integer :: ierr
!
    if (i_verbose >= 1) then
      write(*,*)'@@ dealloc_dst_dens_mat'
    endif
!  
    deallocate (dens_mat_dst, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (dens_mat_dst)'
!                                                                                                              
    deallocate (e_dens_mat_dst, stat=ierr)
    if (ierr /= 0) stop 'Abort:ERROR in dealloc (e_dens_mat_dst)'
!
  end subroutine dealloc_dst_dens_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the TB energy by dens_mat_dst, e_dens_mat_dst
!
!           ETB = 2 Tr [rho H]        ( in imode = 1 )
!           ETB = 2 Tr [pi S]         ( in imode = 2 )
!           ETB = 2 Tr [rho H_{TB0}]  ( in imode = 3 )
!
  subroutine qm_calc_etb_dst(value_of_etb,imode)
!
!
    use M_qm_dst_proj_cell,    only : dst_atm_list, len_dst_atm_list          !(unchnaged)
    use M_qm_geno_dst, only : ham_tot_dst, ham_tb0_dst, overlap_dst !(unchanged)
    use M_qm_domain, only : noav, atm_element, nval, jsv4jsd, njsd  !(unchanged)
    implicit none
    integer      :: imode
    integer, parameter  :: ict4h=1
    real(kind=8) :: value_of_etb
    integer      :: jsv2, nss2, nval2, ja2, jsd1
    integer      :: jsv1, nss1, nval1, ja1
    real(kind=8) :: ddsum, ddd1, ddd2
    integer      :: dst_atm_index

!
!   if (i_verbose >= 1) then
!     write(*,*)'@@ qm_calc_etb:imode=',imode
!   endif
!  
    ddsum=0.0d0  
    do dst_atm_index=1,len_dst_atm_list(1)
      jsv2=dst_atm_list(dst_atm_index)
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         do ja2=1,nval2
           do ja1=1,nval1
             if (imode == 1) then
               ddd1=dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
               ddd2=ham_tot_dst(ja1,ja2,jsd1,dst_atm_index)
!              write(*,*)ja1,ja2,jsd1,jsv2,ddd1,ddd2
             endif  
             if (imode == 2) then
               ddd1=e_dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
               ddd2=overlap_dst(ja1,ja2,jsd1,dst_atm_index)
             endif
             if (imode == 3) then
               ddd1=dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
               ddd2=ham_tb0_dst(ja1,ja2,jsd1,dst_atm_index)
             endif  
             ddsum=ddsum+ddd1*ddd2
           enddo  
         enddo  
      enddo  
    enddo  
    ddsum=ddsum*2.0d0
    value_of_etb=ddsum
!
!     if (i_verbose >= 1) then
!       if (imode == 1) then
!         write(6,*)'ETB as 2 Tr[rho H] =',ddsum
!       endif   
!       if (imode == 2) then
!         write(6,*)'ETB as 2 Tr[pi S]  =',ddsum
!       endif   
!       if (imode == 3) then
!         write(6,*)'ETB as 2 Tr[rho H_{TB0}] =',ddsum
!       endif   
!     endif  
!
  end subroutine qm_calc_etb_dst
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Calculate the TB force by dens_mat_dst, e_dens_mat_dst
!
  subroutine set_tb_force_dst(atm_force_tmp)
!
    use elses_mod_md_dat,   only : itemd
    use M_md_dst,      only : myrank
    use M_qm_dst_proj_cell,    only : dst_atm_list, dst_atm_list_rev, len_dst_atm_list !(unchnaged)
    use M_qm_domain, only : noav, atm_element, nval, jsv4jsd, njsd !(unchanged)
    use M_qm_geno_dst, only : d_ham_tb0_dst, d_overlap_dst
!   use M_qm_domain, only : ddhij, ddsij              !(unchanged)
    use elses_mod_phys_const, only : para_spin_factor !(parameter)
    implicit none
    real(kind=8), intent(inout) :: atm_force_tmp(:,:)
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, njsd1, njsd2, nval1, nval2
    integer :: jsd1, jsd2, nss1, nss2, ja1, ja2
    integer :: dst_atm_index, dst_atm_index1
    integer :: size1, size2
!
    if (.not. allocated(dst_atm_list_rev)) then
      stop 'ERROR(set_tb_force_dst):dst_atm_list_rev is not allocated'
    endif
!
    size1=size(atm_force_tmp,1)
    size2=size(atm_force_tmp,2)
!    
    if (size1 /= 3) then
      write(*,*)'ERROR:atm_force_tmp: size1=',size1
      stop
    endif   
!
    if (size2 /= noav) then
      write(*,*)'ERROR:atm_force_tmp: size2=',size2
      stop
    endif   
!
!   write(*,'(a,2i10,f30.20)')'CHECK-IN:itemd, myrank, sum(ddhij) =',itemd, myrank, sum(ddhij)
!   write(*,'(a,2i10,f30.20)')'CHECK-IN:itemd, myrank, sum(ddsij) =',itemd, myrank, sum(ddsij)
!
    atm_force_tmp(:,:)=0.0d0
!
    do dst_atm_index=1,len_dst_atm_list(1)
      jsv2=dst_atm_list(dst_atm_index)
      njsd2=njsd(jsv2,ict)
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd2
        jsv1=jsv4jsd(jsd1,jsv2)
        dst_atm_index1=dst_atm_list_rev(jsv1)
        if ( (dst_atm_index1 < 1) .or. (dst_atm_index1 > len_dst_atm_list(2)) ) then
          write(*,*)'ERROR(set_tb_force_dst)'
          stop
        endif
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
                 atm_force_tmp(:,jsv2)=atm_force_tmp(:,jsv2) &
&                    +d_ham_tb0_dst(:,ja1,ja2,jsd1,dst_atm_index)*dens_mat_dst(ja1,ja2,jsd1,dst_atm_index) &
&                    -d_overlap_dst(:,ja1,ja2,jsd1,dst_atm_index)*e_dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
                 atm_force_tmp(:,jsv1)=atm_force_tmp(:,jsv1) &
&                    +d_ham_tb0_dst(:,ja2,ja1,jsd2,dst_atm_index1)*dens_mat_dst(ja1,ja2,jsd1,dst_atm_index) &
&                    -d_overlap_dst(:,ja2,ja1,jsd2,dst_atm_index1)*e_dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
              end do
           end do
        endif
      enddo
    enddo
    atm_force_tmp(:,:)  = - para_spin_factor * atm_force_tmp(:,:)
!
  end subroutine set_tb_force_dst
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calculate Mulliken charge :  q_i = sum_j S_ij Rho_ji   
!         OUTPUT : e_num_on_atom_dst, e_num_on_basis_dst
!
  subroutine calc_mulliken_dst
!
    use elses_mod_phys_const, only : para_spin_factor              !(parameter)
    use M_qm_domain, only : noav, atm_element, nval, jsv4jsd, njsd !(unchanged)
!   use M_qm_domain, only : dsij !(unchanged)
    use M_qm_domain, only : DOUBLE_PRECISION          !(parameter)
    use M_qm_dst_proj_cell,    only : dst_atm_list, len_dst_atm_list        !(unchnaged)
    use elses_mod_md_dat,     only : itemd !(unchanged)
    use M_qm_geno_dst, only : overlap_dst  !(unchanged)
    implicit none
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, njsd2, nss1, nss2, nval1, nval2
    integer :: jsd1, ja1, ja2
    real(DOUBLE_PRECISION) :: d_mulliken, dbijd, dsijd, e_num_temp
    integer :: ierr
    integer :: dst_atm_index
!
    ierr=0
    if (.not. allocated(overlap_dst)) ierr=1
    if (.not. allocated(dst_atm_list)) ierr=1
    if (.not. allocated(len_dst_atm_list)) ierr=1
    if (ierr /= 0) then
      write(*,*)'ERROR in allocateion :mulliken_dst_test'
      stop
    endif
!
    if (size(e_num_on_basis_dst,1) /= maxval(nval)) then
      write(*,*)'ERROR (mulliken_dst_test):size(e_num_on_basis_wrk,1)=',size(e_num_on_basis_dst,1)
      stop
    endif
!
    if ( .not. allocated(dens_mat_dst) ) then
      stop 'ERROR in alloc. of dens_mat_dst (mulliken_dst_test)'
    endif
!   
    if ( ( len_dst_atm_list(1) < 1 ) .or. ( len_dst_atm_list(1) > noav ) ) then
      write(*,*)'ERROR  :mulliken_dst_test: len_dst_atm_list(1)=',len_dst_atm_list(1)
      stop
    endif
!
    e_num_on_basis_dst(:,:) = 0.0d0
    e_num_on_atom_dst(:)    = 0.0d0
!
    do dst_atm_index=1,len_dst_atm_list(1)
      jsv2=dst_atm_list(dst_atm_index)
      njsd2=njsd(jsv2,ict)
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      e_num_temp=0d0
      do ja2=1,nval2
        d_mulliken=0.0d0
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          nss1=atm_element(jsv1)
          nval1=nval(nss1)
          do ja1=1,nval1
            dbijd=dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
            dsijd=overlap_dst(ja1,ja2,jsd1,dst_atm_index)
!           if (dabs(dsijd-dsij(ja1,ja2,jsd1,jsv2)) .gt. 1.0d-10) then
!             write(*,'(a,5i10,2f20.10)')'itemd,ja1,ja2,jav1,jsv2,dsij,dsij_dst=', &
!&                                        itemd,ja1,ja2,jsv1,jsv2,dsijd,dsij(ja1,ja2,jsd1,jsv2)
!           endif   
!           dsijd=dsij(ja1,ja2,jsd1,jsv2)
!           if (jsv2 < 3) then
!             if (i_verbose >=0) then
!                write(*,'(a,5i10,2f20.10)')'itemd,ja1,ja2,jsv1,jsv2,dbij,dsij=',itemd,ja1,ja2,jsv1,jsv2,dbijd,dsijd
!             endif  
!           endif  
!           if ((jsv2 < 3) .or. (jsv2 > 58)) then
!             write(*,'(a,5i10,2f20.10)')'itemd,ja1,ja2,jav1,jsv2,dsij,dsij_dst=', &
!&                                        itemd,ja1,ja2,jsv1,jsv2,dsijd,dsij(ja1,ja2,jsd1,jsv2)
!           endif
            d_mulliken=d_mulliken+para_spin_factor*dbijd*dsijd
!              : para_spin_factor(=2.0d0) is multiplied
          enddo
        enddo
        e_num_on_basis_dst(ja2,dst_atm_index) = d_mulliken 
        e_num_temp=e_num_temp+d_mulliken
      enddo
      e_num_on_atom_dst(dst_atm_index)=e_num_temp
    enddo
!
  end subroutine calc_mulliken_dst
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Calculate Mulliken charge :  q_i = sum_j S_ij Rho_ji   
!         OUTPUT : e_num_on_atom_dst, e_num_on_basis_dst
!
  subroutine calc_mulliken_dst_test
!
    use elses_mod_phys_const, only : para_spin_factor             !(parameter)
    use M_qm_domain, only : noav, atm_element, nval, jsv4jsd, njsd, dsij !(unchanged)
    use M_qm_domain, only : DOUBLE_PRECISION          !(parameter)
    use M_qm_dst_proj_cell,    only : dst_atm_list, len_dst_atm_list        !(unchnaged)
    use elses_mod_md_dat,     only : itemd !(unchanged)
    implicit none
    integer, parameter :: ict=1
    integer :: jsv1, jsv2, njsd2, nss1, nss2, nval1, nval2
    integer :: jsd1, ja1, ja2
    real(DOUBLE_PRECISION) :: d_mulliken, dbijd, dsijd, e_num_temp
    integer :: ierr
    integer :: dst_atm_index
!
    ierr=0
    if (.not. allocated(dsij)) ierr=1
    if (.not. allocated(dst_atm_list)) ierr=1
    if (.not. allocated(len_dst_atm_list)) ierr=1
    if (ierr /= 0) then
      write(*,*)'ERROR in allocateion :mulliken_dst_test'
      stop
    endif
!
    if (size(e_num_on_basis_dst,1) /= maxval(nval)) then
      write(*,*)'ERROR (mulliken_dst_test):size(e_num_on_basis_wrk,1)=',size(e_num_on_basis_dst,1)
      stop
    endif
!
    if ( .not. allocated(dens_mat_dst) ) then
      stop 'ERROR in alloc. of dens_mat_dst (mulliken_dst_test)'
    endif
!   
    if ( ( len_dst_atm_list(1) < 1 ) .or. ( len_dst_atm_list(1) > noav ) ) then
      write(*,*)'ERROR  :mulliken_dst_test: len_dst_atm_list(1)=',len_dst_atm_list(1)
      stop
    endif
!
    e_num_on_basis_dst(:,:) = 0.0d0
    e_num_on_atom_dst(:)    = 0.0d0
!
    do dst_atm_index=1,len_dst_atm_list(1)
      jsv2=dst_atm_list(dst_atm_index)
      njsd2=njsd(jsv2,ict)
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      e_num_temp=0d0
      do ja2=1,nval2
        d_mulliken=0.0d0
        do jsd1=1,njsd2
          jsv1=jsv4jsd(jsd1,jsv2)
          nss1=atm_element(jsv1)
          nval1=nval(nss1)
          do ja1=1,nval1
            dbijd=dens_mat_dst(ja1,ja2,jsd1,dst_atm_index)
            dsijd=dsij(ja1,ja2,jsd1,jsv2)
            if (jsv2 < 3) then
              if (i_verbose >=0) then
                 write(*,'(a,5i10,2f20.10)')'itemd,ja1,ja2,jsv1,jsv2,dbij,dsij=',itemd,ja1,ja2,jsv1,jsv2,dbijd,dsijd
              endif  
            endif  
            d_mulliken=d_mulliken+para_spin_factor*dbijd*dsijd
!              : para_spin_factor(=2.0d0) is multiplied
          enddo
        enddo
        e_num_on_basis_dst(ja2,dst_atm_index) = d_mulliken 
        e_num_temp=e_num_temp+d_mulliken
      enddo
      e_num_on_atom_dst(dst_atm_index)=e_num_temp
    enddo
!
  end subroutine calc_mulliken_dst_test
!
end module M_qm_dst_dens_mat
