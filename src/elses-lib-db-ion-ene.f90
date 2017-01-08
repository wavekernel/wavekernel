!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_db_ion_ene
!
! Database of ionization energy
!   Ref. CRC Handbook of Chemistry and Physics, 94th Edition, p. 10-197-198
!        ed. William M. Haynes, CRC Press, Boca Raton (2013)
!          ISBN:978-1466571143
!
   implicit none
!
   private
   public :: get_ion_ene_db_unit
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the ionization energy in the original unit
!
  subroutine get_ion_ene_db_unit(atomic_num, result_value, unit_name)
   implicit none
   integer,           intent(in)    :: atomic_num
   real(8),           intent(out)   :: result_value
   character(len=*), intent(inout)  :: unit_name
!
   unit_name='eV'
   result_value = -1.0d0 ! dummy value
!
   select case(atomic_num)
     case (1) ! H
       result_value = 13.598433d0
     case (2) ! He
       result_value = 24.587387d0
     case (3) ! Li
       result_value = 5.391719d0
     case (4) ! Be
       result_value = 9.32270d0
     case (5) ! B
       result_value = 8.29802d0
     case (6) ! C
       result_value = 11.26030d0
     case (7) ! N
       result_value = 14.5341d0
     case (8) ! O
       result_value = 13.61805d0
     case (9) ! F
       result_value = 17.4228d0
     case (10) ! Ne
       result_value = 21.56454d0
     case (11) ! Na
       result_value = 5.139076d0
     case (12) ! Mg
       result_value = 7.646235d0
     case (13) ! Al
       result_value = 5.985768d0
     case (14) ! Si
       result_value = 8.15168d0
     case (15) ! P
       result_value = 10.48669d0
     case (16) ! S
       result_value = 10.36001d0
     case (17) ! Cl
       result_value = 12.96763d0
     case (18) ! Ar
       result_value = 15.759610d0
     case (19) ! K
       result_value = 4.3406633d0
     case (20) ! Ca
       result_value = 6.11316d0
     case (21) ! Sc
       result_value = 6.56149d0
     case (22) ! Ti
       result_value = 6.82812d0
     case (23) ! V
       result_value = 6.74619d0
     case (24) ! Cr
       result_value = 6.76651d0
     case (25) ! Mn
       result_value = 7.43402d0
     case (26) ! Fe
       result_value = 7.9024d0
     case (27) ! Co
       result_value = 7.88101d0
     case (28) ! Ni
       result_value = 7.6398d0
     case (29) ! Cu
       result_value = 7.72638d0
     case (30) ! Zn
       result_value = 9.394199d0
     case (31) ! Ga
       result_value = 5.999301d0
     case (32) ! Ge
       result_value = 7.89943d0
     case (33) ! As
       result_value = 9.7886d0
     case (34) ! Se
       result_value = 9.75239d0
     case (35) ! Br
       result_value = 11.8138d0
     case (36) ! Kr
       result_value = 13.99961d0
     case (37) ! Rb
       result_value = 4.177128d0
     case (38) ! Sr
       result_value = 5.69485d0
     case (39) ! Y
       result_value = 6.2173d0
     case (40) ! Zr
       result_value = 6.63390d0
     case (41) ! Nb
       result_value = 6.75885d0
     case (42) ! Mo
       result_value = 7.09243d0
     case (43) ! Tc
       result_value = 7.28d0
     case (44) ! Ru
       result_value = 7.36050d0
     case (45) ! Rh
       result_value = 7.45890d0
     case (46) ! Pd
       result_value = 8.3369d0
     case (47) ! Ag
       result_value = 7.57623d0
     case (48) ! Cd
       result_value = 8.99382d0
     case (49) ! In
       result_value = 5.78636d0
     case (50) ! Sn
       result_value = 7.34392d0
     case (51) ! Sb
       result_value = 8.60839d0
     case (52) ! Te
       result_value = 9.0096d0
     case (53) ! I
       result_value = 10.45126d0
     case (54) ! Xe
       result_value = 12.12984d0
     case (55) ! Cs
       result_value = 3.893905d0
     case (56) ! Ba
       result_value = 5.211664d0
     case (57) ! La
       result_value = 5.5769d0
     case (58) ! Ce
       result_value = 5.5387d0
     case (59) ! Pr
       result_value = 5.473d0
     case (60) ! Nd
       result_value = 5.5250d0
     case (61) ! Pm
       result_value = 5.582d0
     case (62) ! Sm
       result_value = 5.6437d0 
     case (63) ! Eu
       result_value = 5.67038d0
     case (64) ! Gd
       result_value = 6.14980d0
     case (65) ! Tb
       result_value = 5.8638d0
     case (66) ! Dy
       result_value = 5.9389d0
     case (67) ! Ho
       result_value = 6.0215d0
     case (68) ! Er
       result_value = 6.1077d0
     case (69) ! Tm
       result_value = 6.18431d0
     case (70) ! Yb
       result_value = 6.25416d0
     case (71) ! Lu
       result_value = 5.42586d0
     case (72) ! Hf
       result_value = 6.82507d0
     case (73) ! Ta
       result_value = 7.54957d0
     case (74) ! W
       result_value = 7.86403d0
     case (75) ! Re
       result_value = 7.83352d0
     case (76) ! Os
       result_value = 8.43823d0
     case (77) ! Ir
       result_value = 8.96702d0
     case (78) ! Pt
       result_value = 8.9588d0
     case (79) ! Au
       result_value = 9.22553d0
     case (80) ! Hg
       result_value = 10.4375d0
     case (81) ! Tl
       result_value = 6.108194d0
     case (82) ! Pb
       result_value = 7.41663d0
     case (83) ! Bi
       result_value = 7.2855d0
     case (84) ! Po
       result_value = 8.414d0
     case (85) ! At
       result_value = -1.0d0
     case (86) ! Rn
       result_value = 10.7485d0
   end select   
!
  end subroutine get_ion_ene_db_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_db_ion_ene

