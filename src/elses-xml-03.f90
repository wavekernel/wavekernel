!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!zzz  @@@@ elses-xml-03.f @@@@@
!zzz  @@@@@ 2008/01/07 @@@@@
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!@@ Read the temperature (level-broadening) 
!                parameter for elctron system
!
!! Copyright (C) ELSES. 2007-2016 all rights reserved
subroutine elses_xml_load_elec_temperature
   use elses_xml_config
   use elses_mod_ctrl,       only : i_verbose
!  use elses_mod_phys_const, only : ev2kel, ev4au
   use elses_mod_phys_const, only : ev4au
   use elses_mod_elec_cond, only : temp_for_electron
!
   implicit none 
   real(8) :: tempd
!
   if (i_verbose >= 1) then
     write(*,*)'@@ LSES_LOAD_ELEC_TEMPERATURE'
   endif  
!
   tempd=1.0d0/200d0
!        ---> temperature in au
!
   temp_for_electron=tempd
!
   if (i_verbose >= 1) then
     write(*,*)' temperature for elec (au,eV)=',        &
&          temp_for_electron,temp_for_electron*ev4au
   endif  

!
end subroutine elses_xml_load_elec_temperature
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
