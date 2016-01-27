!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_db_polarization
!
! Database of dipole polarization 
!   Ref. CRC Handbook of Chemistry and Physics, 94th Edition
!        ed. William M. Haynes, CRC Press, Boca Raton (2013), p. 10-188-189
!          ISBN:978-1466571143
!
!
   implicit none
!
   private
   public :: get_polarization_db_unit
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the dipole polarization in the original unit
!
  subroutine get_polarization_db_unit(atomic_num, result_value, unit_name)
   implicit none
   integer,           intent(in)    :: atomic_num
   real(8),           intent(out)   :: result_value
   character(len=*),  intent(inout) :: unit_name
!
   unit_name='cm3'
   result_value = -1.0d0 ! dummy value
!
   select case(atomic_num)
     case (1) ! H
       result_value = 0.666793d-24
     case (2) ! He
       result_value = 0.2050522d-24
     case (3) ! Li
       result_value = 24.33d-24
     case (4) ! Be
       result_value = 5.60d-24
     case (5) ! B
       result_value = 3.03d-24
     case (6) ! C
       result_value = 1.67d-24
     case (7) ! N
       result_value = 1.10d-24
     case (8) ! O
       result_value = 0.802d-24
     case (9) ! F
       result_value = 0.557d-24
     case (10) ! Ne
       result_value = 0.39432d-24
     case (11) ! Na
       result_value = 24.11d-24
     case (12) ! Mg
       result_value = 10.6d-24
     case (13) ! Al
       result_value = 6.8d-24
     case (14) ! Si
       result_value = 5.53d-24
     case (15) ! P
       result_value = 3.63d-24
     case (16) ! S
       result_value = 2.90d-24
     case (17) ! Cl
       result_value = 2.18d-24
     case (18) ! Ar
       result_value = 1.6411d-24
     case (19) ! K
       result_value = 43.06d-24
     case (20) ! Ca
       result_value = 22.8d-24
     case (21) ! Sc
       result_value = 17.8d-24
     case (22) ! Ti
       result_value = 14.6d-24
     case (23) ! V
       result_value = 12.4d-24
     case (24) ! Cr
       result_value = 11.6d-24
     case (25) ! Mn
       result_value = 9.4d-24
     case (26) ! Fe
       result_value = 8.4d-24
     case (27) ! Co
       result_value = 7.5d-24
     case (28) ! Ni
       result_value = 6.8d-24
     case (29) ! Cu
       result_value = 6.2d-24
     case (30) ! Zn
       result_value = 5.75d-24
     case (31) ! Ga
       result_value = 8.12d-24
     case (32) ! Ge
       result_value = 5.84d-24
     case (33) ! As
       result_value = 4.31d-24
     case (34) ! Se
       result_value = 3.77d-24
     case (35) ! Br
       result_value = 3.05d-24
     case (36) ! Kr
       result_value = 2.4844d-24
     case (37) ! Rb
       result_value = 47.24d-24
     case (38) ! Sr
       result_value = 27.6d-24
     case (39) ! Y
       result_value = 22.7d-24
     case (40) ! Zr
       result_value = 17.9d-24
     case (41) ! Nb
       result_value = 15.7d-24
     case (42) ! Mo
       result_value = 12.8d-24
     case (43) ! Tc
       result_value = 11.4d-24
     case (44) ! Ru
       result_value = 9.6d-24
     case (45) ! Rh
       result_value = 8.6d-24
     case (46) ! Pd
       result_value = 4.8d-24
     case (47) ! Ag
       result_value = 6.78d-24
     case (48) ! Cd
       result_value = 7.36d-24
     case (49) ! In
       result_value = 10.2d-24
     case (50) ! Sn
       result_value = 6.28d-24
     case (51) ! Sb
       result_value = 6.6d-24
     case (52) ! Te
       result_value = 5.5d-24
     case (53) ! I
       result_value = 5.35d-24
     case (54) ! Xe
       result_value = 4.044d-24
     case (55) ! Cs
       result_value = 59.42d-24
     case (56) ! Ba
       result_value = 39.7d-24
     case (57) ! La
       result_value = 31.1d-24
     case (58) ! Ce
       result_value = 29.6d-24
     case (59) ! Pr
       result_value = 28.2d-24
     case (60) ! Nd
       result_value = 31.4d-24
     case (61) ! Pm
       result_value = 30.1d-24
     case (62) ! Sm
       result_value = 28.8d-24
     case (63) ! Eu
       result_value = 27.7d-24
     case (64) ! Gd
       result_value = 23.5d-24
     case (65) ! Tb
       result_value = 25.5d-24
     case (66) ! Dy
       result_value = 24.5d-24
     case (67) ! Ho
       result_value = 23.6d-24
     case (68) ! Er
       result_value = 22.7d-24
     case (69) ! Tm
       result_value = 21.8d-24
     case (70) ! Yb
       result_value = 20.9d-24
     case (71) ! Lu
       result_value = 21.9d-24
     case (72) ! Hf
       result_value = 16.2d-24
     case (73) ! Ta
       result_value = 13.1d-24
     case (74) ! W
       result_value = 11.1d-24
     case (75) ! Re
       result_value = 9.7d-24
     case (76) ! Os
       result_value = 8.5d-24
     case (77) ! Ir
       result_value = 7.6d-24
     case (78) ! Pt
       result_value = 6.5d-24
     case (79) ! Au
       result_value = 5.8d-24
     case (80) ! Hg
       result_value = 5.02d-24
     case (81) ! Tl
       result_value = 7.6d-24
     case (82) ! Pb
       result_value = 6.98d-24
     case (83) ! Bi
       result_value = 7.4d-24
     case (84) ! Po
       result_value = 6.8d-24
     case (85) ! At
       result_value = 6.0d-24
     case (86) ! Rn
       result_value = 5.3d-24
   end select   
!
  end subroutine get_polarization_db_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_db_polarization
