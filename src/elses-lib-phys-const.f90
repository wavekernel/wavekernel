!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module M_lib_phys_const
!
   real(8), parameter :: pi=.314159265358979323D1
!
   real(8), parameter :: ev4au=2.0d0*13.6058d0
!                        --->   1 au = (ev4au) eV
!
   real(8), parameter :: angst=0.529177d0
!                        --->   1 au = 0.529177 A
!
   real(8), parameter :: au_fsec=124.06948d0/3.0d0
!                        --->   1 fsec = (AUFSEC) au
!
   real(8), parameter :: au_mass=1.6605655d-24/9.109534d-28
!
   real(8), parameter :: ev2kel=1.60217733d0/1.380658d0*10000.0d0
!                        --->   1 eV  = (EV2KEL) kelvin
!
   real(8), parameter :: ene_j4ev=1.60217733d-19
!                        ---> 1 [eV] = 1.60217733 x 10^{-19} [J]
!
!  real(8), parameter :: ene_j4ev=1.60219d-19 (OLD)
!                        ---> 1 [eV] = 1.60219 x 10^{-19} [J] (OLD)
!
   real(8), parameter :: kcal_per_mol = 627.50955d0
!                        --->  1 [au] = 627.50955 [kcal/mol]
!
   real(8), parameter :: para_spin_factor=2.0d0
!
   public convert_unit      ! (function)
   public test_convert_unit ! (routine)
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function convert_unit( mode_name, unit_name) result(unit_value)
    implicit none
    character(len=*), intent(in) :: unit_name
    character(len=*), intent(in) :: mode_name
    real(8) :: c, unit_value
!
    c=1.0d0
    select case(unit_name)
    case("au", "AU", "a.u.", "A.U.")
       c = 1.d0
!
    case("hartree", "Hartree", "HARTREE")
       c = 1.d0
    case("eV", "ev", "EV")
       c = ev4au
    case("J", "Joule")
       c = ev4au*ene_j4ev 
    case("ry", "Ry", "RY", "rydberg", "Rydberg")
       c = 2.d0
    case("K","kelvin","Kelvin")
       c = ev4au*ev2kel
!
    case("ang","Ang","angstrom","Angstrom")
       c = angst
    case("m","meter")
       c = angst * 1.0d-10
    case("cm","centimeter")
       c = angst * 1.0d-8
!
    case("cm3","centimeter3")
       c = ( angst * 1.0d-8 )**3
!
    case("Pa")
       c = (ev4au * ene_j4ev) / (angst * 1.0d-10)**3 ! Pa = J / m^3
!
    case("amu")
       c = 1.0d0/au_mass
    case("fsec")
       c = 1.0d0/au_fsec
    case default
       write(6,*)'# Error (Convert_unit): unit= ', trim(unit_name)
       stop
    end select
!
    select case(mode_name)
    case("into")
       unit_value=c
    case("from")
       unit_value=1.0d0/c
    case default
       write(6,*)'# Error (Convert_unit): mode= ', trim(mode_name)
       stop
    end select   

  end function convert_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine test_convert_unit(log_unit)
    implicit none
    real(8) :: value
    integer, intent(in) :: log_unit
!
    if (log_unit <= 0) return
!
    write(log_unit,'(a)')'@@@ test_convert_unit' 
!
    value=1.0d0
    write(log_unit,*)'energy: 1 au = X eV       : X=', value*convert_unit('into','eV')
    write(log_unit,*)'energy: 1 au = X J        : X=', value*convert_unit('into','J')
!
    value=1.0d0
    write(log_unit,*)'length: X au = 1 Angstrom : X=', value*convert_unit('from','Angstrom')
!
    value=1.0d0
    write(log_unit,*)'length: 1 au = X m        : X=', value*convert_unit('into','m')
!
    value=1.0d0
    write(log_unit,*)'pressure: 1 au = X GPa    : X=', value*convert_unit('into','Pa')*1.0d-9
!
!
  end subroutine test_convert_unit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
end module M_lib_phys_const

