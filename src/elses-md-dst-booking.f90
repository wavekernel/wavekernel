!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_dst_booking
!
  implicit none
!
  integer, allocatable :: len_bk_list_dst(:,:) !  (DIFFERENT VALUE AMONG THE NODES)
!        ! len_bk_list_dst(:,1) : length of the booking list for the distributed atoms
!
  integer, allocatable :: bk_list_dst(:,:,:) !  (DIFFERENT VALUE AMONG THE NODES)
!        ! bk_list_dst(:,:,1) : lthe booking list for the distributed atoms
!
  private
!
  public len_bk_list_dst
!
  public set_len_bk_list_dst
  public dealloc_bk_list_dst
!
  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Get the length of the booking list for a given radius
!       OUTPUT : len_bk_list_dst(:,1) in the module variable
!       
!         non-order-N calculation (temporaly) !!!
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
   subroutine set_len_bk_list_dst(cutoff_radius)
!
     use M_qm_dst_proj_cell,   only  : len_dst_atm_list, dst_atm_list !(unchanged)
!
     use M_qm_domain, only : i_verbose, noav, ax, ay, az, & 
&                  i_pbc_x, i_pbc_y, i_pbc_z, atm_position   !(unchanged)
!
     implicit none
     real(8), intent(in)  :: cutoff_radius
     integer i_show, ierr
     integer dst_atm_index, jsv1, jsv2, jsd
     real(8) ::  rcut1
     real(8) ::  dvecx, dvecy, dvecz, w
!
     if (i_verbose >= 1) then
       write(6,*)'@@@ set_len_bk_list_dst'
       write(6,*)' cutoff radius =',cutoff_radius
     endif  
!
     i_show=5
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Check the parapmeters
!
     if (cutoff_radius < 0.0d0) then 
       write(6,*)'ERROR:set_len_bk_list_dst'
       stop
     endif
!
     if ( ( len_dst_atm_list(1) < 1 ) .or. (len_dst_atm_list(1) > noav ) ) then 
       write(6,*)'ERROR:set_len_bk_list_dst'
       stop
     endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     if ( allocated(len_bk_list_dst) ) then
       deallocate(len_bk_list_dst, stat=ierr)
       if (ierr /= 0) stop 'ERROR in dealloc. len_bk_list_dst'
     endif   
!
     allocate (len_bk_list_dst(len_dst_atm_list(1),2), stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. len_bk_list_dst'
     len_bk_list_dst(:,:)=noav
!
     if( cutoff_radius .lt. huge(1d0) ) then
        rcut1=cutoff_radius*1.001d0
     else
        return
     end if
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!$omp  parallel & 
!$omp& default(shared) &
!$omp& private(dvecx,dvecy,dvecz,w) & 
!$omp& private(jsv1,jsv2,jsd, dst_atm_index)
!$omp  do schedule(static)
     do dst_atm_index=1,len_dst_atm_list(1)
       jsv2=dst_atm_list(dst_atm_index)
!
       jsd=1
!   
       do jsv1=1,noav
         if (jsv1 == jsv2) cycle
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
         w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
         if (w .lt. rcut1) then
           jsd=jsd+1
         endif
       enddo  
       len_bk_list_dst(dst_atm_index,1)=jsd
!
     enddo
!$omp end do
!$omp end parallel 
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   end subroutine set_len_bk_list_dst
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Set the booking list for the distributed atoms
!       OUTPUT : bk_list_dst(:,:,1) in the module variable
!       
!         non-order-N calculation (temporaly) !!!
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
   subroutine set_bk_list_dst(cutoff_radius,alloc_size)
!
     use M_qm_dst_proj_cell,   only  : len_dst_atm_list, dst_atm_list !(unchanged)
!
     use M_qm_domain, only : i_verbose, noav, ax, ay, az, & 
&                  i_pbc_x, i_pbc_y, i_pbc_z, atm_position   !(unchanged)
!
     implicit none
     real(8), intent(in)  :: cutoff_radius
     integer, intent(in)  :: alloc_size
     integer i_show, ierr
     integer dst_atm_index, jsv1, jsv2, jsd
     real(8) ::  rcut1
     real(8) ::  dvecx, dvecy, dvecz, w
!
     if (i_verbose >= 1) then
       write(6,*)'@@@ get_length_of_list'
       write(6,*)' cutoff radius =',cutoff_radius
     endif  
!
     i_show=5
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Check the parapmeters
!
     if (cutoff_radius < 0.0d0) then 
       write(6,*)'ERROR:set_len_bk_list_dst'
       stop
     endif
!
     if ( ( len_dst_atm_list(1) < 1 ) .or. (len_dst_atm_list(1) > noav ) ) then 
       write(6,*)'ERROR:set_len_bk_list_dst'
       stop
     endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
     allocate (bk_list_dst(alloc_size, len_dst_atm_list(1), 2), stat=ierr)
     if (ierr /= 0) stop 'ERROR in alloc. len_bk_list_dst'
     bk_list_dst(:,:,:)=0
!
     if( cutoff_radius .lt. huge(1d0) ) then
        rcut1=cutoff_radius*1.001d0
     else
        rcut1=cutoff_radius
     end if
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!$omp  parallel & 
!$omp& default(shared) &
!$omp& private(dvecx,dvecy,dvecz,w) & 
!$omp& private(jsv1,jsv2,jsd, dst_atm_index)
!$omp  do schedule(static)
     do dst_atm_index=1,len_dst_atm_list(1)
       jsv2=dst_atm_list(dst_atm_index)
!
       jsd=1
!   
       do jsv1=1,noav
         if (jsv1 == jsv2) cycle
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
         w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
         if (w .lt. rcut1) then
           jsd=jsd+1
         endif
       enddo  
       len_bk_list_dst(dst_atm_index,1)=jsd
!
     enddo
!$omp end do
!$omp end parallel 
!     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   end subroutine set_bk_list_dst
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  @@ Deallocate the array in the module variable
!
      !! Copyright (C) ELSES. 2007-2016 all rights reserved
   subroutine dealloc_bk_list_dst
!
     implicit none
     integer :: ierr
!
     deallocate (len_bk_list_dst, stat=ierr)
     if (ierr /= 0) stop 'ERROR in dealloc. len_bk_list_dst'
!
     deallocate (bk_list_dst, stat=ierr)
     if (ierr /= 0) stop 'ERROR in dealloc. bk_list_dst'
!
   end subroutine dealloc_bk_list_dst
!
end module M_md_dst_booking
