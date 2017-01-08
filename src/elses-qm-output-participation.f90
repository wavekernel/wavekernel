!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_output_participation

  use M_qm_domain,        only : i_verbose !(unchanged)
  implicit none
!
  private
  public :: calc_participation
  public :: gene_participation
!
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine gene_participation(part_ratio,log_unit)
!
    use elses_mod_phys_const, only : ev4au  !(unchanged)
    use elses_arr_eig_leg, only : atmp, eig2, f_occ  !(unchanged)
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged)
    use  M_qm_domain,   only : noav !(unchanged)
    use elses_mod_file_io,      only : vacant_unit  !(function)
!
    implicit none
    real(kind=8), intent(inout) :: part_ratio(:)
    integer,      intent(in)    :: log_unit
    integer :: k, j
    integer :: atm_index, orb_index 
    integer :: ierr
    real(kind=8), allocatable ::  weight(:,:)
    real(kind=8) :: wfn, p_ratio, wt_sum1, wt_sum2
    integer :: iunit
!
    if (log_unit >0) then 
      if (i_verbose >= 1) write(log_unit,'(a)')'@@ gene_participation'
    endif
!
    part_ratio(:)=0.0d0
!
    allocate (weight(noav,2), stat=ierr)
    weight(:,:)=0.0d0
!
    iunit=vacant_unit()
!   open(iunit, file=filename)
!
    do k=1, n_tot_base
      weight(:,:)=0.0d0
      do j=1, n_tot_base
        atm_index=j2js(j) 
        orb_index=j2ja(j) 
        wfn=atmp(j,k)
!       write(99, '(4i10,f20.10)') k, j, atm_index, orb_index, wfn
        weight(atm_index,1)=weight(atm_index,1)+wfn*wfn
        weight(atm_index,2)=weight(atm_index,2)+wfn*wfn*wfn*wfn
      enddo   
!     do atm_index=1, noav
!       write(98, '(2i10,f20.10)') k, atm_index, weight(atm_index,1)
!     enddo   
      wt_sum1=sum(weight(:,1))
      wt_sum2=sum(weight(:,2))
      if (wt_sum2 < -1.0d-10) then
        write(*,*)'Error(calc_participation):k, weight ', k, wt_sum2
        stop
      endif   
      p_ratio=wt_sum1*wt_sum1/wt_sum2
      part_ratio(k)=p_ratio
!     write(iunit, '(i10,3f20.10)') k, eig2(k)*ev4au, p_ratio, f_occ(k)
    enddo   
!
    deallocate(weight, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR(calc_participation)' 
      stop
    endif   
!
!   close(iunit)
!
!   stop 'Stop manually in calc_participation'
!
  end subroutine gene_participation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine calc_participation
!
    use M_config                            !(unchanged)
    use elses_mod_phys_const, only : ev4au  !(unchanged)
    use elses_arr_eig_leg, only : atmp, eig2, f_occ, eigenvectors, desc_eigenvectors  !(unchanged)
    use elses_mod_orb2,  only : j2js,j2ja,js2j,n_tot_base  !(unchanged)
    use  M_qm_domain,   only : noav !(unchanged)
    use elses_mod_file_io,      only : vacant_unit  !(function)
    use M_lib_dst_info, only : root_node  !(unchanged)
    use M_eig_solver_center, only : gather_vector_in_matrix  !(function)
!
    implicit none
    character(len=*), parameter :: filename_header="output_participation.txt"
    character(len=128)          :: filename
    integer :: k, j
    integer :: atm_index, orb_index 
    integer :: ierr
    real(kind=8), allocatable ::  weight(:,:)
    real(kind=8) :: wfn, p_ratio, wt_sum1, wt_sum2
    integer :: iunit
    integer :: lu, i_v
    real(8), allocatable :: eigenvector(:)
!
    lu  = config%calc%distributed%log_unit
    i_v = config%option%verbose
!
    if (i_v >= 1) then
      if (lu > 0) write(lu,*)'@@ calc_participation'
    endif
!
!
    allocate (weight(noav,2), eigenvector(n_tot_base), stat=ierr)
    if (root_node) then
      weight(:,:)=0.0d0
      iunit=vacant_unit()
      filename=trim(config%option%output_dir)//filename_header
      open(iunit, file=filename)
    end if
!
    do k=1, n_tot_base
      if (trim(config%calc%solver%scheme) == 'eigen_mpi') then
        call gather_vector_in_matrix(eigenvectors, desc_eigenvectors, k, eigenvector)
      else
        eigenvector(1:n_tot_base) = atmp(1:n_tot_base,k)
      end if
      if (root_node) then
        weight(:,:)=0.0d0
        do j=1, n_tot_base
          atm_index=j2js(j)
          orb_index=j2ja(j)
          wfn=eigenvector(j)
!         write(99, '(4i10,f20.10)') k, j, atm_index, orb_index, wfn
          weight(atm_index,1)=weight(atm_index,1)+wfn*wfn
          weight(atm_index,2)=weight(atm_index,2)+wfn*wfn*wfn*wfn
        enddo
!       do atm_index=1, noav
!         write(98, '(2i10,f20.10)') k, atm_index, weight(atm_index,1)
!       enddo
        wt_sum1=sum(weight(:,1))
        wt_sum2=sum(weight(:,2))
        if (wt_sum2 < -1.0d-10) then
          write(*,*)'Error(calc_participation):k, weight ', k, wt_sum2
          stop
        endif
        p_ratio=wt_sum1*wt_sum1/wt_sum2
        write(iunit, '(i10,3f20.10)') k, eig2(k)*ev4au, p_ratio, f_occ(k)
      end if
    enddo   
!
    deallocate(weight, eigenvector, stat=ierr)
    if (ierr /= 0) then
      write(*,*)'ERROR(calc_participation)' 
      stop
    endif   
!
    if (root_node) then
      close(iunit)
    end if
!
!   stop 'Stop manually in calc_participation'
!
  end subroutine calc_participation
!
end module M_output_participation

