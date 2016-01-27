!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_io_ctrl_dos_plot
!
! use M_qm_domain,          only : i_verbose !(unchanged)
  use elses_mod_phys_const, only : ev4au     !(parameter)
!
   private
!
! Public routines
   public initial_preparation_for_dos
!
   contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initial_preparation_for_dos(imode, number_energy_mesh, energy_btm, energy_top, width_of_peak, i_verbose_in)
!
    use elses_mod_file_io, only : vacant_unit
!
    implicit none
    integer,          intent(out)  :: imode, number_energy_mesh
    real(8),          intent(out)  :: energy_btm, energy_top, width_of_peak
    integer,          optional     :: i_verbose_in
    integer :: iunit
    integer :: i_verbose
!
!
    character(len=100) :: file_name
    logical :: file_exist
!
    file_name='input_energy_mesh.txt'
!
    if (present(i_verbose_in)) then 
      i_verbose = i_verbose_in
    else
      i_verbose = 1
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Dummy values
!
    imode=0
    number_energy_mesh=0
    energy_btm=0.0d0
    energy_top=0.0d0
    width_of_peak=0.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    inquire (file=trim(file_name), exist=file_exist)
!
    if (file_exist .eqv. .false.) then
      if (i_verbose >= 1) write(*,*)'The file does not exist:',file_name
      return
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    iunit=vacant_unit()
!
    open(iunit, file=file_name, status='unknown')
!
    if (i_verbose >= 1) then 
      write(*,*)' open file for energy mesh:',file_name
    endif  
!
    read(iunit,*) imode
    read(iunit,*) number_energy_mesh
    read(iunit,*) energy_btm
    read(iunit,*) energy_top
!      ----> enery botom and top for DOS plotting [eV]
    read(iunit,*) width_of_peak
!
    energy_btm=energy_btm/ev4au
    energy_top=energy_top/ev4au
    width_of_peak=width_of_peak/ev4au
!
    close(iunit)
!
   end subroutine initial_preparation_for_dos
!
end module M_io_ctrl_dos_plot



