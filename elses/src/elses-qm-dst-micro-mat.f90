!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_micro_mat
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
  implicit none
!
  private
  public set_ham_tb0_and_overlap_dstm
  public set_size_for_dstm_micro_mat
! public set_size_for_dst_micro_mat
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_ham_tb0_and_overlap_dstm(atm_index_seed, jsv4jsk, num_atom_proj, &
&                 booking_list_dstm, booking_list_dstm_len, overlap_dstm, ham_tb0_dstm, d_overlap_dstm, d_ham_tb0_dstm)
!               !THIS ROUTINE SHOULD BE CALLED WITHIN OMP LOOP
!
    use M_qm_domain,       only : nval, atm_element !(unchanged)
    use elses_mod_sel_sys, only : r_cut_book        !(unchanged)
    use M_md_get_distance, only : get_distance             !(routine)
    use M_md_get_distance, only : get_vector_for_atom_pair !(routine)
    use M_qm_geno_Huckel,  only : SetOverlap, SetNonDiagonalElements !(routines)
    use M_qm_domain,       only : c_system, S_is_I !(unchanged)
!
    implicit none
    integer, intent(in)     :: jsv4jsk(:)
    integer, intent(in)     :: num_atom_proj
    integer, intent(in)     :: atm_index_seed
    integer, intent(out)    :: booking_list_dstm(:,:)
    integer, intent(out)    :: booking_list_dstm_len(:)
    real(DOUBLE_PRECISION), intent(out)  ::   overlap_dstm(:,:,:,:)
    real(DOUBLE_PRECISION), intent(out)  ::   ham_tb0_dstm(:,:,:,:)
    real(DOUBLE_PRECISION), intent(out)  :: d_overlap_dstm(:,:,:,:,:)
    real(DOUBLE_PRECISION), intent(out)  :: d_ham_tb0_dstm(:,:,:,:,:)
!
    integer                 :: atm_index1, atm_index2
    integer                 :: ict4h, ishow, ierr
    integer                 :: size1, size2, size3, size4
    integer                 :: orb_index, orb_size
    integer                 :: jsk1, jsk2, jsd
    integer                 :: nval_max, nss1, nss2
    real(DOUBLE_PRECISION)  :: cutoff_radius, rcut1, dist
    real(DOUBLE_PRECISION)  :: dvecx, dvecy, dvecz
!
    real(DOUBLE_PRECISION), allocatable :: ham_tb0_dstm_diag(:,:)
    real(DOUBLE_PRECISION), allocatable :: ham_tb0_onsite(:,:,:) 
    real(DOUBLE_PRECISION), allocatable :: ham_tb0_offsite(:,:) 
    real(DOUBLE_PRECISION), allocatable :: d_ham_tb0_offsite(:,:,:) 
    real(DOUBLE_PRECISION), allocatable :: overlap_offsite(:,:) 
    real(DOUBLE_PRECISION), allocatable :: d_overlap_offsite(:,:,:) 
!
    if (c_system /= 'geno') then
      write(*,*)'INFO(set_ham_tb0_and_overlap_dstm):c_system=', c_system
      write(*,*)'  S_is_I=',S_is_I
    endif   
!
    size1=size(ham_tb0_dstm,1)
    size2=size(ham_tb0_dstm,2)
    size3=size(ham_tb0_dstm,3)
    size4=size(ham_tb0_dstm,4)
    nval_max=maxval(nval)
!
    if (num_atom_proj /= size4) then
      write(*,*)'ERROR(set_ham_and_overlap_dst_micro):size4=',size4
    endif
!
    if (nval_max /= size1) then
      write(*,*)'ERROR(set_ham_and_overlap_dst_micro):size1=',size1
    endif
!
    if (nval_max /= size2) then
      write(*,*)'ERROR(set_ham_and_overlap_dst_micro):size2=',size2
    endif
!
    allocate (ham_tb0_dstm_diag(nval_max,num_atom_proj),stat=ierr)
    if( ierr /= 0) stop 'ERROR in alloc (ham_tb0_dstm_diag)'
!
    booking_list_dstm(:,:)=0
    ham_tb0_dstm(:,:,:,:)=0.0d0
!
    if (.not. S_is_I) then
      overlap_dstm(:,:,:,:)=0.0d0
    endif  
!
    cutoff_radius=r_cut_book
    if( cutoff_radius .lt. huge(1d0) ) then
       rcut1=cutoff_radius*1.001d0
    else
       rcut1=cutoff_radius
    endif
!    
!   stop 'Stop manually:T. Hoshi(1)'
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
!   stop 'Stop manually:T. Hoshi(2)'
!
    if (.not. S_is_I) then
      allocate (overlap_offsite(nval_max, nval_max),stat=ierr)
      if( ierr /= 0) stop 'ERROR in alloc (overlap_offsite)'
!
      allocate (d_overlap_offsite(3, nval_max, nval_max),stat=ierr)
      if( ierr /= 0) stop 'ERROR in alloc (d_overlap_offsite)'
    endif  
!
!   stop 'Stop manually:T. Hoshi(3)'
!
    jsk1=1
    if (atm_index_seed /= jsv4jsk(jsk1)) then
      stop 'ERROR:set_ham_and_overlap_dst_micro'
    endif
!
!   stop 'Stop manually:T. Hoshi(4)'
!
    if (S_is_I) then
      write(*,*)'STOP manually:  S_is_I=',S_is_I
      stop
    endif   
!       
!   stop 'Stop manually:T. Hoshi(5)'
!
    do jsk1=1,num_atom_proj
      atm_index1=jsv4jsk(jsk1)
      orb_size=nval(atm_element(atm_index1))
      call set_mat_tb0_diag_loc(atm_index1, ham_tb0_dstm_diag(1:orb_size,jsk1))
      jsd=1
      booking_list_dstm(jsd,jsk1)=jsk1 ! booking the 'self' atom as jsd=1
      ham_tb0_dstm(:,:,jsd,jsk1)=0.0d0
      overlap_dstm(:,:,jsd,jsk1)=0.0d0
      do orb_index=1,nval_max
        ham_tb0_dstm(orb_index,orb_index,jsd,jsk1)=ham_tb0_dstm_diag(orb_index,jsk1)
        overlap_dstm(orb_index,orb_index,jsd,jsk1)=1.0d0
      enddo
    enddo
!
    jsd=1
    do jsk2=1,num_atom_proj
      atm_index2=jsv4jsk(jsk2)
      nss2=atm_element(atm_index2)
      do orb_index=1, nval(nss2)
        ham_tb0_onsite(orb_index,orb_index,2)= ham_tb0_dstm_diag(orb_index,jsk2)
      enddo
      jsd=1
      do jsk1=1,num_atom_proj
        atm_index1=jsv4jsk(jsk1)
        if ( atm_index1 == atm_index2 ) cycle
        call get_vector_for_atom_pair(atm_index1, atm_index2, dvecx, dvecy, dvecz)
        dist=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
        if (dist > rcut1) cycle
        jsd=jsd+1
        booking_list_dstm(jsd,jsk2)=jsk1 ! booking the 'self' atom as jsd=1
        nss1=atm_element(atm_index1)
        do orb_index=1, nval(nss1)
          ham_tb0_onsite(orb_index,orb_index,1)= ham_tb0_dstm_diag(orb_index,jsk1)
        enddo
        call SetOverlap(nss1,nss2,dvecx,dvecy,dvecz,overlap_offsite(:,:),d_overlap_offsite(:,:,:))
        call SetNondiagonalElements(nss1,nss2,dvecx,dvecy,dvecz,&
&             overlap_offsite(:,:), &
&             ham_tb0_onsite(:,:,1), ham_tb0_onsite(:,:,2), &
&             ham_tb0_offsite(:,:), d_overlap_offsite(:,:,:), d_ham_tb0_offsite(:,:,:))
          ham_tb0_dstm(:, :,    jsd, jsk2) =   ham_tb0_offsite(:,:)
        d_ham_tb0_dstm(:, :, :, jsd, jsk2) = d_ham_tb0_offsite(:,:,:)
          overlap_dstm(:, :,    jsd, jsk2) =   overlap_offsite(:,:)
        d_overlap_dstm(:, :, :, jsd, jsk2) = d_overlap_offsite(:,:,:)
      enddo
      booking_list_dstm_len(jsk2)=jsd  ! # of booked atoms
    enddo
!
!
  end subroutine set_ham_tb0_and_overlap_dstm
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_mat_tb0_diag_loc(atm_index, mat_tb0_diag_loc)
!      ----> Set diagonal elements for the given atom index
!
    use M_qm_domain,      only : atm_element, e_num_on_basis !(unchanged)
    use M_qm_geno_Huckel, only : SetDiagonalElements,SetDiagonalShifts !(routines)
    implicit none
    integer, intent(in)     :: atm_index
    real(DOUBLE_PRECISION), intent(out)  :: mat_tb0_diag_loc(:)
    real(DOUBLE_PRECISION), allocatable  :: wrk_mat(:,:)
    real(DOUBLE_PRECISION), allocatable  :: e_num_on_basis_wrk(:)
    integer            :: ierr
    integer            :: nss1, orb_index, orb_size
!
    nss1=atm_element(atm_index)
    orb_size=size(mat_tb0_diag_loc,1)
!
    allocate (wrk_mat(orb_size, orb_size),stat=ierr)
    if( ierr /= 0) stop 'ERROR in alloc (wrk_mat)'
!
    allocate (e_num_on_basis_wrk(orb_size),stat=ierr)
    if( ierr /= 0) stop 'ERROR in alloc (e_num_on_basis_wrk)'
!
    e_num_on_basis_wrk(:)=0.0d0
    if (allocated(e_num_on_basis)) then
      if (size(e_num_on_basis,1) < size(e_num_on_basis_wrk,1)) then
        write(*,*)'ERROR!:size of e_num_on_basis_wrk =',size(e_num_on_basis_wrk,1)
        write(*,*)'   size of e_num_on_basis         =',size(e_num_on_basis,1)
        write(*,*)'   atm_index =', atm_index
        write(*,*)'   nss1      =', nss1
        write(*,*)'  orb_size   =', orb_size
        stop
      else
        e_num_on_basis_wrk(1:orb_size)=e_num_on_basis(1:orb_size,atm_index) 
      endif   
    endif   
!
    call SetDiagonalElements(nss1,e_num_on_basis_wrk,wrk_mat(:,:))
!     Note : e_num_on_basis is not needed (for non VOIP method)
!
    call SetDiagonalShifts(nss1,wrk_mat(:,:))
!
    do orb_index=1,orb_size
      mat_tb0_diag_loc(orb_index)=wrk_mat(orb_index,orb_index)
    enddo
!
    deallocate (wrk_mat,stat=ierr)
    if( ierr /= 0) stop 'ERROR in dealloc (wrk_mat)'
!
   end subroutine set_mat_tb0_diag_loc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_size_for_dstm_micro_mat(atm_index, jsv4jsk, num_atom_proj, size_for_interacted_atoms)  
!               !THIS ROUTINE SHOULD BE CALLED WITHIN OMP LOOP
!
!   use M_qm_domain,   only : njsd,noav !(unchanged)
    use elses_mod_sel_sys, only : r_cut_book        !(unchanged)
    use M_md_get_distance, only : get_vector_for_atom_pair !(routine)
!
    implicit none
    integer, intent(in)     :: jsv4jsk(:)
    integer, intent(in)     :: num_atom_proj
    integer, intent(in)     :: atm_index
    integer, intent(out)    :: size_for_interacted_atoms
    integer                 :: atm_index1, atm_index2
    integer                 :: jsk1, jsk2, jsd
!   integer                 :: ict4h, ishow
!   logical, parameter      :: check_mode=.true.
    real(DOUBLE_PRECISION)  :: cutoff_radius, rcut1
    real(DOUBLE_PRECISION)  :: dist, dvecx, dvecy, dvecz
    logical                 :: central_atom_is_included
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the interaction radius
!
    cutoff_radius=r_cut_book
    if( cutoff_radius .lt. huge(1d0) ) then
       rcut1=cutoff_radius*1.001d0
    else
       rcut1=cutoff_radius
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    central_atom_is_included = .false.
    size_for_interacted_atoms=0
!
    do jsk2=1,num_atom_proj
      atm_index2=jsv4jsk(jsk2)
      if (atm_index2 == atm_index) central_atom_is_included = .true.
      jsd=1
      do jsk1=1,num_atom_proj
        atm_index1=jsv4jsk(jsk1)
        if ( atm_index1 == atm_index2 ) cycle
        call get_vector_for_atom_pair(atm_index1, atm_index2, dvecx, dvecy, dvecz)
        dist=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
        if (dist > rcut1) cycle
        jsd=jsd+1
      enddo
      size_for_interacted_atoms=max(size_for_interacted_atoms,jsd)
    enddo
!
    if (.not. central_atom_is_included) then
      write(*,*)'ERROR(set_size_for_dstm_micro_mat)'
      write(*,*)' central_atom_is_included =',central_atom_is_included
      stop
    endif
!
  end subroutine set_size_for_dstm_micro_mat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_size_for_dst_micro_mat(atm_index2, jsv4jsk, num_atom_proj, size_for_interacted_atoms)  
!               !THIS ROUTINE SHOULD BE CALLED WITHIN OMP LOOP
!
    use M_qm_domain,   only : njsd,noav !(unchanged)
    implicit none
    integer, intent(in)     :: jsv4jsk(:)
    integer, intent(in)     :: num_atom_proj
    integer, intent(in)     :: atm_index2
    integer, intent(out)    :: size_for_interacted_atoms
    integer                 :: atm_index1
    integer                 :: jsk1, jsk2
    integer                 :: ict4h, ishow
    logical, parameter      :: check_mode=.true.
    real(DOUBLE_PRECISION)  :: dist
!
    ict4h=1
    ishow=noav
!
    size_for_interacted_atoms=0
    do jsk1=1,num_atom_proj
      atm_index1=jsv4jsk(jsk1)
      if (check_mode) then 
        if ( (njsd(atm_index1, ict4h) < 1) .or. ( njsd(atm_index1, ict4h) > noav )) then
          write(*,*)'ERROR(set_size_for_dst_micro_mat):njsd=',njsd(atm_index1, ict4h)
        endif
      endif
      size_for_interacted_atoms=max(njsd(atm_index1, ict4h), size_for_interacted_atoms)
    enddo
!
!   if (i_verbose >= 1) then
!     if (atm_index2 <= ishow) then
!       write(*,*)'atm_index, mat size for int atoms=',atm_index2, size_for_interacted_atoms
!     endif   
!   endif
!
!   stop 'Stop manually'
!
  end subroutine set_size_for_dst_micro_mat
!
end module M_qm_dst_micro_mat


