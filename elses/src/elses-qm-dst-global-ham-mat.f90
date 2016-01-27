!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_dst_global_ham_mat
!
  use M_qm_domain,   only : i_verbose, DOUBLE_PRECISION !(unchanged)
!
  implicit none
!
   private
!
! Public routines
   public :: set_ham_tb0_and_overlap_global
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine set_ham_tb0_and_overlap_global
!
     use M_qm_domain,   only : atm_element, noav, nval, jsv4jsd, njsd !(unchanged)
     use M_qm_domain,   only : ham_tb0, dhij, dsij, ddhij, ddsij, dham_tb0  !(CHANGED)
     use M_qm_geno_dst, only : ham_tb0_dst, d_ham_tb0_dst, overlap_dst, d_overlap_dst !(unchanged)
     use M_qm_dst_proj_cell,  only : dst_atm_list, len_dst_atm_list  !(unchanged)
     use M_md_dst,      only : set_dst_final                   !(unchanged)
     use M_lib_mpi_wrapper, only : mpi_wrapper_allreduce_r4, mpi_wrapper_allreduce_r5 !(unchanged)
     use M_md_dst,      only : myrank
     use elses_mod_md_dat,     only : itemd
!
     implicit none
     integer, parameter :: ict=1
     integer            :: ierr
     integer            :: dst_atm_index, atm_index, jsv2, jsd1, jsv1, ja1, ja2
!
     if (i_verbose >= 1) then
       write(*,*)'@@ set_ham_tb0_and_overlap_global'
     endif
!  
     if (.not. allocated(ham_tb0)) then
       stop 'ERROR:Not allocated : ham_tb0'
     endif
!   
     if (i_verbose >= 0) then
       do dst_atm_index=1,len_dst_atm_list(1)
         atm_index=dst_atm_list(dst_atm_index)
         jsv2=atm_index
         do jsd1=1,njsd(jsv2,ict)
           jsv1=jsv4jsd(jsd1,jsv2)
           do ja2=1,nval(atm_element(jsv2))
             do ja1=1,nval(atm_element(jsv1))
               if ((atm_index < 5) .and. (atm_index > 59)) then 
                 write(*,'(a,5i10,2f20.10)')'s_dst_data=',myrank, &
&                 ja1,ja2,jsv1,jsv2, overlap_dst(ja1,ja2,jsd1,dst_atm_index)
               endif  
             enddo
           enddo  
         enddo
       enddo
     endif
!
     ham_tb0(:,:,:,:)=0.0d0
     dsij(:,:,:,:)=0.0d0
     dham_tb0(:,:,:,:,:)=0.0d0
     ddsij(:,:,:,:,:)=0.0d0
!
!    call plot_ham_tb0_and_overlap_global
!    stop
!
!$omp  parallel default(shared) &
!$omp& private(dst_atm_index, atm_index, jsv2, jsd1)
!$omp  do schedule(static)
     do dst_atm_index=1,len_dst_atm_list(1)
       atm_index=dst_atm_list(dst_atm_index)
       jsv2=atm_index
       do jsd1=1,njsd(jsv2,ict)
         ham_tb0(:,:,jsd1,jsv2)    =   ham_tb0_dst(:,:,  jsd1,dst_atm_index)
         dsij(:,:,jsd1,jsv2)       =   overlap_dst(:,:,  jsd1,dst_atm_index)
         dham_tb0(:,:,:,jsd1,jsv2) = d_ham_tb0_dst(:,:,:,jsd1,dst_atm_index) 
         ddsij(:,:,:,jsd1,jsv2)    = d_overlap_dst(:,:,:,jsd1,dst_atm_index) 
       enddo
     enddo
!$omp end do
!$omp end parallel
!
     call mpi_wrapper_allreduce_r4(ham_tb0)
     call mpi_wrapper_allreduce_r4(dsij)
     call mpi_wrapper_allreduce_r5(dham_tb0)
     call mpi_wrapper_allreduce_r5(ddsij)
!
     dhij(:,:,:,:)  = ham_tb0(:,:,:,:)
     ddhij(:,:,:,:,:)=dham_tb0(:,:,:,:,:)
!
     write(*,*)'myrank,sum(dhij)  =',myrank, itemd, sum(dhij)
     write(*,*)'myrank,sum(dsij)  =',myrank, itemd, sum(dsij)
     write(*,*)'myrank,sum(ddhij) =',myrank, itemd, sum(ddhij)
     write(*,*)'myrank,sum(ddsij) =',myrank, itemd, sum(ddsij)
!
     if (myrank == 0) then
       call plot_ham_tb0_and_overlap_global
     endif  
!    call set_dst_final
!    stop
!
   end subroutine set_ham_tb0_and_overlap_global
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine plot_ham_tb0_and_overlap_global
!
     use M_qm_domain, only : atm_element, noav, nval, jsv4jsd, njsd !(unchanged)
     use M_qm_domain, only : ham_tb0, dsij !(unchanged)
     implicit none
     integer :: atm_index, jsv2, jsd1, jsv1, ja2, ja1
     integer, parameter :: ict=1
!
     if (i_verbose >= 1) return
!
     do atm_index=1,noav
       jsv2=atm_index
       do jsd1=1,njsd(jsv2,ict)
         jsv1=jsv4jsd(jsd1,jsv2)
         do ja2=1,nval(atm_element(jsv2))
           do ja1=1,nval(atm_element(jsv1))
             if (jsv2 < 5) then
               write(*,'(a,4i10,2f20.10)')'ham_data=',ja1,ja2,jsv1,jsv2, &
&                                        ham_tb0(ja1,ja2,jsd1,jsv2),dsij(ja1,ja2,jsd1,jsv2)
             endif
           enddo
         enddo  
       enddo
     enddo
!
   end subroutine plot_ham_tb0_and_overlap_global
!
end module M_qm_dst_global_ham_mat
