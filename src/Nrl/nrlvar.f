!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
!zzz  @@@@ nrlvar.f @@@@@
!zzz  @@@@@ 2007/04/01 @@@@@
!ccc2007cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0316: Modified for the compatibility of LSES code
!ccc2008cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!0401: Prepared (NT08A-123p23); i_pbc_x, i_pbc_y, i_pbc_z
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module MNrlVariables
        use elses_mod_phys_const, only : aumass => au_mass
        use elses_mod_mass,     only : awt, amm
        use elses_mod_js4jsv,   only : js4jsv 
        use elses_mod_jsv4jsd,  only : jsv4jsd, njsd
        use elses_mod_tx,       only : tx, ty, tz, jsei
        use elses_mod_orb1,     only : nval
        use elses_mod_noav,     only : noav
!       use elses_mod_sim_cell, only : iperiodic, ax, ay, az
        use elses_mod_sim_cell, only : ax, ay, az,
     +                                 i_pbc_x, i_pbc_y, i_pbc_z
!       target :: js4jsv, jsv4jsd, jsei, nval, njsd
      end module
