!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_output_rest_pair
!
  use M_qm_domain, only : i_verbose, DOUBLE_PRECISION !(unchanged)
!
  implicit none
  integer, allocatable :: called_first(:)
!
  private
!
! Public routines
  public plot_rest_energy_pair
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! @@ Plot the rest-part energy for each atom pair
!
  subroutine plot_rest_energy_pair
!
    use elses_mod_file_io, only : vacant_unit !(function)
    use M_lib_phys_const, only : ev4au               !(unchanged) 
    use M_qm_domain ,  only : i_verbose, atm_element, nval, &
&                             jsv4jsd, njsd, noav, atm_element !(unchanged) 
    use elses_mod_js4jsv,    only : js4jsv !(unchanged)
    use M_config,            only : config !(unchanged)
    use elses_mod_md_dat,    only : itemd  !(unchanged)
    use M_qm_geno_Huckel,    only : RepulsiveEnergy !(function)
    use M_md_get_distance,   only : get_vector_for_atom_pair !(routine)
!
!
    implicit none
!
!   integer      :: imode
!   integer, parameter  :: ict4h=1
!   real(kind=8) :: value_of_etb
    integer      :: jsv2, nss2, nval2, ja2, jsd1
    integer      :: jsv1, nss1, nval1, ja1
    real(kind=8) :: ddd, ddsum
    real(kind=8) :: dvecx, dvecy, dvecz
    integer      :: n_count
    integer      :: js1, js2
    integer      :: iunit, ierr
    logical      :: flag_for_first_write
    integer, parameter  :: ict4h = 1
    character(len=32) :: filename_wrk
    real(DOUBLE_PRECISION) :: repulsive_force(3)
!
    if ( .not. config%output%bond_list%set ) return
    if (mod(itemd-1,config%output%bond_list%interval) /= 0 ) return
    filename_wrk='output_rest_energy_pair.txt'
!
    if (.not. allocated(called_first)) then
      flag_for_first_write = .true.
      allocate(called_first(1), stat=ierr)
      if (ierr /=0 ) stop 'Alloc. error (plot_rest_energy_pair):called_first'
    else
      flag_for_first_write = .false.
    endif   
!
    if (i_verbose >= 1) then
      write(*,*)'@@ plot_rest_energy_pair'
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iunit=vacant_unit()
!
    if ( flag_for_first_write ) then
      if (i_verbose >= 1) then
        write(*,*)' The output rest-energy-pair file is initialized : ',trim(filename_wrk)
      endif
      open(iunit, file=filename_wrk, status='unknown')
    else
      open(iunit, file=filename_wrk, status='unknown', position='append')
    endif   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Count the number of the atom pairs
!  
    n_count=0
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      do jsd1=1,njsd(jsv2,ict4h)
        jsv1=jsv4jsd(jsd1,jsv2)
        n_count=n_count+1
!       if (jsv1 /= jsv2) n_count=n_count+1
      enddo
    enddo
    if (i_verbose >= 1) then
      write(*,*)'   n_count=',n_count
    endif
!
    write(iunit,'(a)')'%%MatrixMarket matrix coordinate real general'
    write(iunit,'(a,i10)')'% step_count=', config%system%structure%mdstep
    write(iunit,*)noav, noav, n_count
!
    ddsum=0.0d0
    do jsv2=1,noav
      nss2=atm_element(jsv2)
      nval2=nval(nss2)
      js2=js4jsv(jsv2)
      do jsd1=1,njsd(jsv2,ict4h)
         jsv1=jsv4jsd(jsd1,jsv2)
         js1=js4jsv(jsv1)
         nss1=atm_element(jsv1)
         nval1=nval(nss1)
         call get_vector_for_atom_pair(jsv1, jsv2, dvecx, dvecy, dvecz)
!          ----> (dvecx, dvecy, dvecz) : Vector r = R_2 - R_1
!
         if (jsv1 == jsv2) then
           ddd=0.0d0
         else
           ddd=RepulsiveEnergy(jsv1, jsv2, dvecx, dvecy, dvecz, repulsive_force)
         endif   
         ddsum=ddsum+ddd
!
         write(iunit,'(2i10,f20.10)')jsv1, jsv2, ddd*ev4au
      enddo  
    enddo  
!
    close(iunit)
!
    if (i_verbose >= 1) then
      write(*,'(a,2f20.10)')'INFO:Rest energy [au,eV] =', ddsum, ddsum*ev4au
    endif
!
  end subroutine plot_rest_energy_pair
!
end module M_qm_output_rest_pair


