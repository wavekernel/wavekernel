!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_get_distance
!
   use  M_qm_domain, only : i_verbose, DOUBLE_PRECISION
   implicit none
!
   private
!
! Public routine
   public  :: get_distance
   public  :: get_vector_for_atom_pair
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Get the distance for an atom pair
!
   subroutine get_vector_for_atom_pair(jsv1, jsv2, dvecx, dvecy, dvecz)
!
     use M_qm_domain, only : noav, ax, ay, az             !(unchanged)
     use M_qm_domain, only : i_pbc_x, i_pbc_y, i_pbc_z    !(unchanged)
     use M_qm_domain, only : atm_position                 !(unchanged)
!
     implicit none
     integer, intent(in)                  :: jsv1, jsv2
     real(DOUBLE_PRECISION), intent(out)  :: dvecx, dvecy, dvecz
!
     if (jsv1 == jsv2) then
      dvecx=0.0d0
      dvecy=0.0d0
      dvecz=0.0d0
      return
     endif
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
!
   end subroutine get_vector_for_atom_pair
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  @@ Get the distance for an atom pair
!
   subroutine get_distance(jsv1, jsv2, dist)
!
     use M_qm_domain, only : noav, ax, ay, az             !(unchanged)
     use M_qm_domain, only : i_pbc_x, i_pbc_y, i_pbc_z    !(unchanged)
     use M_qm_domain, only : atm_position                 !(unchanged)
!
     implicit none
     integer, intent(in)                  :: jsv1, jsv2
     real(DOUBLE_PRECISION), intent(out)  :: dist
!
     real(DOUBLE_PRECISION) :: dvecx, dvecy, dvecz
!
     if (jsv1 == jsv2) then
      dist=0.0d0
      return
     endif
!   
     call get_vector_for_atom_pair(jsv1, jsv2, dvecx, dvecy, dvecz)
!
     dist=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
   end subroutine get_distance
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_md_get_distance


