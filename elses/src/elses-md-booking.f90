!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_md_booking
!
  public :: get_max_size_of_booking_list
  public :: plot_booking_list 
!
  contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Get the maximum size of booking list for a given radius
!       (Non order-N routine)
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine get_max_size_of_booking_list(rcut,max_size)
    use elses_mod_ctrl,       only : i_verbose
    use elses_mod_noav,       only : noav
    use elses_mod_sim_cell,   only : noa, ax, ay, az, &
&                                    i_pbc_x, i_pbc_y,i_pbc_z 
    use elses_mod_tx,         only : tx, ty, tz
    use elses_mod_js4jsv,     only : js4jsv
!
    implicit none
    real(kind=8), intent(in)  ::  rcut
    integer,      intent(out) ::  max_size

!
    integer ishow,ict
    real(kind=8) :: rcut1
    real(kind=8) ::  ddd,dvecx,dvecy,dvecz,dxc,dyc,dzc,w
!
    integer ipe,npe
    integer js1,js2,jsv1,jsv2,js1d,js1b,jsd
    integer count, max_count
!
!
    if (i_verbose >= 1) then
      write(*,*)'@@ get_max_size_of_booking_list:rcut=',rcut
    endif  
!
    rcut1=rcut*1.001d0
!      radius  for booking
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    max_count=0
!$omp  parallel &
!$omp& default(shared) &
!$omp& private(ipe,npe) & 
!$omp& private(dvecx,dvecy,dvecz,w) &
!$omp& private(jsv1,jsv2,js1,js2,js1d,js1b,count) &
!$omp& reduction(max : max_count)
!     ipe=omp_get_thread_num()+1
!     npe=omp_get_num_threads()
!     write(6,*)'ipe,npe=',ipe,npe
!$omp  do schedule(static)
      do jsv2=1,noav
        js2=js4jsv(jsv2)
        count=0
        do jsv1=1,noav
          js1=js4jsv(jsv1)
          dvecx=tx(js2)-tx(js1)
          dvecy=ty(js2)-ty(js1)
          dvecz=tz(js2)-tz(js1)
!
          if (i_pbc_x == 1)dvecx = modulo(dvecx+0.5d0,1.0d0) - 0.5d0
          if (i_pbc_y == 1)dvecy = modulo(dvecy+0.5d0,1.0d0) - 0.5d0
          if (i_pbc_z == 1)dvecz = modulo(dvecz+0.5d0,1.0d0) - 0.5d0
!
          dvecx=dvecx*ax
          dvecy=dvecy*ay
          dvecz=dvecz*az
          w=dsqrt(dvecx*dvecx+dvecy*dvecy+dvecz*dvecz)
!
          if (w .le. rcut1) then
            count=count+1
          endif  
!
        enddo
!       write(*,*)'js,count=',jsv2,count
        max_count=max(max_count,count)
      enddo
!$omp end do
!$omp end parallel 
!     
      max_size=max_count

  end subroutine get_max_size_of_booking_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot the booking list into the file 'Booking_info.txt'
!
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine plot_booking_list 
    use elses_mod_ctrl,       only : i_verbose
    use elses_mod_phys_const, only : angst
    use elses_mod_noav,       only : noav
    use elses_mod_sim_cell,   only : noa, ax, ay, az, &
&                                    i_pbc_x, i_pbc_y,i_pbc_z 
    use elses_mod_tx,         only : tx, ty, tz
    use elses_mod_js4jsv,     only : js4jsv
    use elses_mod_jsv4jsd,    only : jsv4jsd, njsd
    use elses_mod_file_io,    only : vacant_unit
    use elses_mod_multi,      only : ict4h
!
    implicit none
    character(len=*), parameter :: FNAME="info-booking-list.txt"
    integer ishow,iunit, njsd2
    real(kind=8) :: dxc,dyc,dzc,drr
!
    integer js1,js2,jsv1,jsv2,js1d,js1b,jsd1
!
    iunit=vacant_unit()
    if (i_verbose >= 1) then 
      write(*,*)' @@ plot_booking_list'
      write(*,*)'   ---> file: info-booking-list.txt'
    endif   
    open(iunit,file=FNAME,status='replace')
!
    do jsv2=1,noav
      js2=js4jsv(jsv2)
      njsd2=njsd(jsv2,ict4h)
      write(iunit,'(A4,I20,I20)') 'atom', jsv2, njsd2
      do jsd1=1,njsd2
        jsv1=jsv4jsd(jsd1,jsv2)
        js1=js4jsv(jsv1)
        dxc=tx(js2)-tx(js1)
        dyc=ty(js2)-ty(js1)
        dzc=tz(js2)-tz(js1)
        if (i_pbc_x == 1) dxc = modulo(dxc + 0.5d0, 1.0d0) - 0.5d0
        if (i_pbc_y == 1) dyc = modulo(dyc + 0.5d0, 1.0d0) - 0.5d0
        if (i_pbc_z == 1) dzc = modulo(dzc + 0.5d0, 1.0d0) - 0.5d0
        dxc=dxc*ax
        dyc=dyc*ay
        dzc=dzc*az
        drr=dsqrt(dxc*dxc+dyc*dyc+dzc*dzc)*angst
!           ---> Distance | R_1 - R_2 | in Angstrom
       write(iunit,'(3I10,F15.8)') jsv2,jsv1,jsd1,drr
      enddo
    enddo   
!
  end subroutine plot_booking_list
!
end module M_md_booking

