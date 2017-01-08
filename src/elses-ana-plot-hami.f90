!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!zzz  @@@@ elses-ana-plot-hami.f @@@@@
!zzz  @@@@@ 2008/03/31 @@@@@
!2008!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!0101T.Hoshi; A modi. (v.0.0.0a_2008_01_01_hoshi)
!0331T.Hoshi; Plot the mapping 
!               ( basis number ) <-----> ( atom, orbital )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @@ Plot the hamiltonian into file
!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_plot_hami
   use elses_mod_tx,         only : jsei
   use elses_mod_js4jsv,     only : js4jsv 
   use elses_mod_jsv4jsd,    only : jsv4jsd, njsd
   use elses_mod_orb1,       only : nval
   use elses_mod_orb2,       only : j2js,j2ja,js2j,n_tot_base
!  use elses_arr_dbij,       only : dbij
   use elses_arr_dhij,       only : dhij
   use elses_mod_multi,      only : ict4h
   use elses_mod_orb2,       only : js2j
   use elses_mod_noav,       only : noav
   use elses_mod_ctrl,       only : i_verbose
   use elses_mod_phys_const, only : ev4au
   use elses_mod_file_io
!
!  use elses_param_ctl_kr, only : noav_kr_def
!
   implicit none
!
   integer :: ii, iisum, iisum2
   integer :: jsv2,      js2, ja2, njsd2, nss, nval2, jg
   integer :: jsv1,jsd1, js1, ja1,             nval1, ig
   real(8) :: dlijd, ahijd, dtmp
   integer :: iunit
   integer :: jj, js, ja
   character(len=*), parameter :: FNAME="Hamiltonian.txt"
!
   if (i_verbose >= 1) then 
     write(*,*)'@@ Plot the Hamiltonian into file'
   endif   
!
   iunit=vacant_unit()
   if (i_verbose >= 1) then 
     write(*,*)' get vacant unit:iunit=',iunit
   endif   
   open(iunit,file=FNAME,status='replace')

!!!!!!!!!!!!!!!!!
! Count up the number of non-zero elements
!
   iisum=0
   iisum2=0
   do jsv2=1,noav
      js2=js4jsv(jsv2)
      njsd2=njsd(jsv2,ict4h)
      nss=jsei(js2)
      nval2=nval(nss)
      iisum2=iisum2+nval2
      do jsd1=1,njsd2
        jsv1=jsv4jsd(jsd1,jsv2)
        if ((jsv1 .le. 0) .or. (jsv1 .gt. noav)) then
          write(6,*)'error:(plot_hami):jsv1,noav=',jsv1,noav
          stop
        endif
        js1=js4jsv(jsv1)
        if ((js1 .le. 0) .or. (js1 .gt. noav)) then
          write(6,*)'error!(enetb):jsv1,js1=',jsv1,js1
          stop
        endif
        nss=jsei(js1)
        nval1=nval(nss)
        do ja2=1,nval2
          do ja1=1,nval1
            iisum=iisum+1
          enddo
        enddo
     enddo
   enddo
!
   if (i_verbose >= 1) then 
     write(*,*)' noav,iisum,iisum2=',noav,iisum,iisum2
   endif   
   write(iunit,"(i15)") iisum2
!      ---> number of matrix dimension
   write(iunit,"(i15)") iisum
!      ---> number of non-zero elements
!
!
   ii=0
   do jsv2=1,noav
     js2=js4jsv(jsv2)
     njsd2=njsd(jsv2,ict4h)
     nss=jsei(js2)
     nval2=nval(nss)
     do ja2=1,nval2
         jg=js2j(ja2,js2)
      do jsd1=1,njsd2
        jsv1=jsv4jsd(jsd1,jsv2)
        if ((jsv1 .le. 0) .or. (jsv1 .gt. noav)) then
          write(6,*)'error!(enetb):jsv1,noav=',jsv1,noav
          stop
        endif
        js1=js4jsv(jsv1)
        if ((js1 .le. 0) .or. (js1 .gt. noav)) then
          write(6,*)'error!(enetb):jsv1,js1=',jsv1,js1
          stop
        endif
        nss=jsei(js1)
        nval1=nval(nss)
        do ja1=1,nval1
          ig=js2j(ja1,js1)
!         dlijd=dbij(ja1,ja2,jsd1,jsv2)
          ahijd=dhij(ja1,ja2,jsd1,jsv2)*ev4au
          dtmp=0.0d0
          ii=ii+1
          write(iunit,"(3i15,2f24.18)") ii,jg,ig,ahijd,dtmp
        enddo
      enddo   
     enddo   
   enddo   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plot the mapping:  ( basis number ) <-----> ( atom, orbital )
!
   if (i_verbose >= 1) then 
     write(*,*)' n_tot_base=',n_tot_base
   endif   
!
   if ((n_tot_base <= 1) .or. (n_tot_base > noav*1000)) then 
     write(*,*)'ERROR:n_tot_base=',n_tot_base
     stop
   endif   
!
   do jj=1,n_tot_base
      js=j2js(jj)
      ja=j2ja(jj)
      write(iunit,"(3i15)") jj,js,ja
   enddo
!
   close(iunit)
!
end subroutine elses_plot_hami
