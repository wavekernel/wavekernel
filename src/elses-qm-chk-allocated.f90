!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_qm_chk_allocated
!
!
   use M_config,           only : config !(unchanged)
!
   implicit none
!
   private
   public chk_allocated_dst
!
   contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine chk_allocated_dst
!
     use M_qm_domain
     use elses_mod_tx,         only : tx, ty, tz, jsei
     use elses_mod_txp,        only : txp, typ, tzp
     use elses_mod_foi,        only : foi
     use elses_mod_foiold,     only : foiold
     use elses_mod_iflag,      only : iflag
     use elses_mod_mass,       only : amm
     use elses_mod_elem_name,  only : elem_name
     use elses_mod_vel,        only : velx, vely, velz
     use elses_mod_js4jsv,     only : js4jsv , jsv4js
     use elses_mod_jsv4jsd,    only : jsv4jsd, njsd
     use elses_mod_ene,        only : enea
     use elses_mod_multi,      only : intf
!
     use elses_mod_orb2,       only : j2js,j2ja,js2j, dbx, dby, dbz
     use elses_mod_r_base,     only : r_base
     use elses_arr_dbij,       only : dbij 
     use elses_arr_dhij,       only : dhij 
     use elses_arr_dfij,       only : dfij 
     use elses_arr_dbij2,      only : dbij2
     use arr_kry_ldos,         only : dldos
     use elses_arr_dhij_cohp,  only : dhij_cohp
     use elses_arr_kry_glo,    only : rA, rB, rE, rR, rW, rxmu, rRho, rFocc, rRNave, rRNdia, nrecg 
     use elses_arr_kry_glo,    only : rcut_kry, noak_str, noak_min, rKnstr, jsv4jsk_str
!
     implicit none
     integer :: iv, lu
     character(len=5) chara
!
     iv = config%option%verbose
     lu = config%calc%distributed%log_unit
! 
     if (lu <= 0) return
     if (iv <= 0) return
!
     write(lu,*)'@@ chk_allocated_dst (check allocate status)'
     write(lu,*)'INFO:use_matom                           =', config%system%structure%use_matom
     write(lu,*)'INFO:calc_force_mode                     =', trim(config%calc%calc_force_mode)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables in elses-mod-01.f
!
     write(lu,*)'INFO:allocated(tx)                       =', allocated(tx)
     write(lu,*)'INFO:allocated(ty)                       =', allocated(ty)
     write(lu,*)'INFO:allocated(tz)                       =', allocated(tz)
     write(lu,*)'INFO:allocated(jsei)                     =', allocated(jsei)
     write(lu,*)'INFO:allocated(txp)                      =', allocated(txp)
     write(lu,*)'INFO:allocated(typ)                      =', allocated(typ)
     write(lu,*)'INFO:allocated(tzp)                      =', allocated(tzp)
     write(lu,*)'INFO:allocated(foi)                      =', allocated(foi)
     write(lu,*)'INFO:allocated(foiold)                   =', allocated(foiold)
     write(lu,*)'INFO:allocated(iflag)                    =', allocated(iflag)
     write(lu,*)'INFO:allocated(elem_name)                =', allocated(elem_name)
     write(lu,*)'INFO:allocated(amm)                      =', allocated(amm)
     write(lu,*)'INFO:allocated(velx)                     =', allocated(velx)
     write(lu,*)'INFO:allocated(vely)                     =', allocated(vely)
     write(lu,*)'INFO:allocated(velz)                     =', allocated(velz)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables in elses-mod-02.f
!
     write(lu,*)'INFO:allocated(js4jsv)                   =', allocated(js4jsv)
     write(lu,*)'INFO:allocated(jsv4js)                   =', allocated(jsv4js)
     write(lu,*)'INFO:allocated(jsv4jsd)                  =', allocated(jsv4jsd)
     write(lu,*)'INFO:allocated(njsd)                     =', allocated(njsd)
     write(lu,*)'INFO:allocated(enea)                     =', allocated(enea)
     write(lu,*)'INFO:allocated(intf)                     =', allocated(intf)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables in elses-mod-03.f
!
     write(lu,*)'INFO:allocated(nval)                     =', allocated(nval)
     write(lu,*)'INFO:allocated(j2js)                     =', allocated(j2js)
     write(lu,*)'INFO:allocated(j2ja)                     =', allocated(j2ja)
     write(lu,*)'INFO:allocated(js2j)                     =', allocated(js2j)
     write(lu,*)'INFO:allocated(dbx)                      =', allocated(dbx)
     write(lu,*)'INFO:allocated(dby)                      =', allocated(dby)
     write(lu,*)'INFO:allocated(dbz)                      =', allocated(dbz)
     write(lu,*)'INFO:allocated(r_base)                   =', allocated(r_base)
     write(lu,*)'INFO:allocated(dbij)                     =', allocated(dbij)
     write(lu,*)'INFO:allocated(dhij)                     =', allocated(dhij)
     write(lu,*)'INFO:allocated(dfij)                     =', allocated(dfij)
     write(lu,*)'INFO:allocated(dbij2)                    =', allocated(dbij2)
     write(lu,*)'INFO:allocated(dldos)                    =', allocated(dldos)
     write(lu,*)'INFO:allocated(dhij_cohp)                =', allocated(dhij_cohp)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables in elses-mod-04.f
!
     write(lu,*)'INFO:allocated(rA)                       =', allocated(rA)
     write(lu,*)'INFO:allocated(rB)                       =', allocated(rB)
     write(lu,*)'INFO:allocated(rE)                       =', allocated(rE)
     write(lu,*)'INFO:allocated(rR)                       =', allocated(rR)
     write(lu,*)'INFO:allocated(rW)                       =', allocated(rW)
     write(lu,*)'INFO:allocated(rxmu)                     =', allocated(rxmu)
     write(lu,*)'INFO:allocated(rRho)                     =', allocated(rRho)
     write(lu,*)'INFO:allocated(rFocc)                    =', allocated(rFocc)
     write(lu,*)'INFO:allocated(rRNave)                   =', allocated(rRNave)
     write(lu,*)'INFO:allocated(rRNdia)                   =', allocated(rRNdia)
     write(lu,*)'INFO:allocated(nrecg)                    =', allocated(nrecg)
     write(lu,*)'INFO:allocated(rcut_kry)                 =', allocated(rcut_kry)
     write(lu,*)'INFO:allocated(noak_str)                 =', allocated(noak_str)
     write(lu,*)'INFO:allocated(noak_min)                 =', allocated(noak_min)
     write(lu,*)'INFO:allocated(jsv4jsk_str)              =', allocated(jsv4jsk_str)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables in elses-qm-domain.f90
!
     write(lu,*)'INFO:alloced(atm_position)               =', allocated(atm_position)
     write(lu,*)'INFO:alloced(atm_force)                  =', allocated(atm_force)
     write(lu,*)'INFO:alloced(atm_force_csc)              =', allocated(atm_force_csc)
     write(lu,*)'INFO:alloced(atm_force_tb0)              =', allocated(atm_force_tb0)
     write(lu,*)'INFO:alloced(atm_element)                =', allocated(atm_element)
     write(lu,*)'INFO:alloced(ham_tb0)                    =', allocated(ham_tb0)
     write(lu,*)'INFO:alloced(ham_csc)                    =', allocated(ham_csc)
     write(lu,*)'INFO:alloced(gamma_csc)                  =', allocated(gamma_csc)
     write(lu,*)'INFO:alloced(dgamma_csc)                 =', allocated(dgamma_csc)
     write(lu,*)'INFO:alloced(fgamma_csc)                 =', allocated(fgamma_csc)
!
     write(lu,*)'INFO:alloced(e_num_on_basis)             =', allocated(e_num_on_basis)
     write(lu,*)'INFO:alloced(previous_e_num_on_basis)    =', allocated(previous_e_num_on_basis)
     write(lu,*)'INFO:alloced(delta_e_numm)               =', allocated(delta_e_num)
     write(lu,*)'INFO:alloced(e_num_on_atom)              =', allocated(e_num_on_atom)
     write(lu,*)'INFO:alloced(tau_csc)                    =', allocated(tau_csc)
     write(lu,*)'INFO:alloced(previous_e_num_on_atom)     =', allocated(previous_e_num_on_atom)
     write(lu,*)'INFO:alloced(ddhij)                      =', allocated(ddhij)
     write(lu,*)'INFO:alloced(ddsij)                      =', allocated(ddsij)
     write(lu,*)'INFO:alloced(dham_tb0)                   =', allocated(dham_tb0)
!
!    write(*,*)'Press any key'
!    read(5,*) chara
!
   end subroutine chk_allocated_dst
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_qm_chk_allocated

