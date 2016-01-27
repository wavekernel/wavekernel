!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_db_covalent_radius
!
! Database of covalent radius
!   Ref. CRC Handbook of Chemistry and Physics, 94th Edition
!        ed. William M. Haynes, CRC Press, Boca Raton (2013), p. 9-49-50
!          ISBN:978-1466571143
!
   implicit none
!
   private
   public :: get_covalent_radius_db_unit
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the covalent radius in the original unit
!
  subroutine get_covalent_radius_db_unit(atomic_num, result_value, unit_name)
   implicit none
   integer,           intent(in)    :: atomic_num
   real(8),           intent(out)   :: result_value
   character(len=*),  intent(inout) :: unit_name
!
   unit_name='Angstrom'
   result_value = -1.0d0 ! dummy value
!
   select case(atomic_num)
     case (1) ! H
       result_value = 0.32d0
     case (2) ! He
       result_value = 0.37d0
     case (3) ! Li
       result_value = 1.30d0
     case (4) ! Be
       result_value = 0.99d0
     case (5) ! B
       result_value = 0.84d0
     case (6) ! C
       result_value = 0.75d0
     case (7) ! N
       result_value = 0.71d0
     case (8) ! O
       result_value = 0.64d0
     case (9) ! F
       result_value = 0.60d0
     case (10) ! Ne
       result_value = 0.62d0
     case (11) ! Na
       result_value = 1.60d0
     case (12) ! Mg
       result_value = 1.40d0
     case (13) ! Al
       result_value = 1.24d0
     case (14) ! Si
       result_value = 1.14d0
     case (15) ! P
       result_value = 1.09d0
     case (16) ! S
       result_value = 1.04d0
     case (17) ! Cl
       result_value = 1.00d0
     case (18) ! Ar
       result_value = 1.01d0
     case (19) ! K
       result_value = 2.00d0
     case (20) ! Ca
       result_value = 1.74d0
     case (21) ! Sc
       result_value = 1.59d0
     case (22) ! Ti
       result_value = 1.48d0
     case (23) ! V
       result_value = 1.44d0
     case (24) ! Cr
       result_value = 1.30d0
     case (25) ! Mn
       result_value = 1.29d0
     case (26) ! Fe
       result_value = 1.24d0
     case (27) ! Co
       result_value = 1.18d0
     case (28) ! Ni
       result_value = 1.17d0
     case (29) ! Cu
       result_value = 1.22d0
     case (30) ! Zn
       result_value = 1.20d0
     case (31) ! Ga
       result_value = 1.23d0
     case (32) ! Ge
       result_value = 1.20d0
     case (33) ! As
       result_value = 1.20d0
     case (34) ! Se
       result_value = 1.18d0
     case (35) ! Br
       result_value = 1.17d0
     case (36) ! Kr
       result_value = 1.16d0
     case (37) ! Rb
       result_value = 2.15d0
     case (38) ! Sr
       result_value = 1.90d0
     case (39) ! Y
       result_value = 1.76d0
     case (40) ! Zr
       result_value = 1.64d0
     case (41) ! Nb
       result_value = 1.56d0
     case (42) ! Mo
       result_value = 1.46d0
     case (43) ! Tc
       result_value = 1.38d0
     case (44) ! Ru
       result_value = 1.36d0
     case (45) ! Rh
       result_value = 1.34d0
     case (46) ! Pd
       result_value = 1.30d0
     case (47) ! Ag
       result_value = 1.36d0
     case (48) ! Cd
       result_value = 1.40d0
     case (49) ! In
       result_value = 1.42d0
     case (50) ! Sn
       result_value = 1.40d0
     case (51) ! Sb
       result_value = 1.40d0
     case (52) ! Te
       result_value = 1.37d0
     case (53) ! I
       result_value = 1.36d0
     case (54) ! Xe
       result_value = 1.36d0
     case (55) ! Cs
       result_value = 2.38d0
     case (56) ! Ba
       result_value = 2.06d0
     case (57) ! La
       result_value = 1.94d0
     case (58) ! Ce
       result_value = 1.84d0
     case (59) ! Pr
       result_value = 1.90d0
     case (60) ! Nd
       result_value = 1.88d0
     case (61) ! Pm
       result_value = 1.86d0
     case (62) ! Sm
       result_value = 1.85d0
     case (63) ! Eu
       result_value = 1.83d0
     case (64) ! Gd
       result_value = 1.82d0
     case (65) ! Tb
       result_value = 1.81d0
     case (66) ! Dy
       result_value = 1.80d0
     case (67) ! Ho
       result_value = 1.79d0
     case (68) ! Er
       result_value = 1.77d0
     case (69) ! Tm
       result_value = 1.77d0
     case (70) ! Yb
       result_value = 1.78d0
     case (71) ! Lu
       result_value = 1.74d0
     case (72) ! Hf
       result_value = 1.64d0
     case (73) ! Ta
       result_value = 1.58d0
     case (74) ! W
       result_value = 1.50d0
     case (75) ! Re
       result_value = 1.41d0
     case (76) ! Os
       result_value = 1.36d0
     case (77) ! Ir
       result_value = 1.32d0
     case (78) ! Pt
       result_value = 1.30d0
     case (79) ! Au
       result_value = 1.30d0
     case (80) ! Hg
       result_value = 1.32d0
     case (81) ! Tl
       result_value = 1.44d0
     case (82) ! Pb
       result_value = 1.45d0
     case (83) ! Bi
       result_value = 1.50d0
     case (84) ! Po
       result_value = 1.42d0
     case (85) ! At
       result_value = 1.48d0
     case (86) ! Rn
       result_value = 1.46d0
   end select   
!
  end subroutine get_covalent_radius_db_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_db_covalent_radius

