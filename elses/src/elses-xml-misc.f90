!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module elses_xml_misc
  implicit none

contains

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine XML_error(str)
    character(len=*), intent(in), optional   :: str
    write(unit=0,fmt="(a)") "XML error"
    if( present(str) ) then
       write(unit=0,fmt="(a)") trim(str)
    endif
    stop
  end subroutine XML_error


  function XML_TO_AU( unit ) result(c)
    character(len=*), intent(in) :: unit
    real*8 :: c

    c = 1.d0
    select case(unit)
    case("au", "AU", "a.u.", "A.U.")
       c = 1.d0
    case("hartree", "Hartree", "HARTREE")
       c = 1.d0
    case("eV", "ev", "EV")
       c = 1.d0/(2.0d0*13.6058d0)
    case("ry", "Ry", "RY", "rydberg", "Rydberg")
       c = 1.d0/2.d0
    case("K","kelvin","Kelvin")
       c = 1.d0/((2.0d0*13.6058d0)*(1.60217733d0/1.380658d0*10000.0d0))
    case("ang","Ang","angstrom","Angstrom")
       c = 1.d0/0.529177d0
    case("amu")
       c = 1.6605655d-24/9.109534d-28
    case("fsec")
       c = 124.06948d0/3.0d0
    case default
       write(*,'(a,a)') '# Error!(XML_TO_AU): unknown unit ', trim(unit)
       stop
    end select
  end function XML_TO_AU

  function XML_AU_TO( unit ) result(c)
    use M_lib_phys_const
    character(len=*), intent(in) :: unit
    real(kind(1d0)) :: c

    c = 1.d0
    select case(unit)
    case("au", "AU", "a.u.", "A.U.")
       c = 1.d0
    case("hartree", "Hartree", "HARTREE")
       c = 1.d0
    case("eV", "ev", "EV")
       !c = 2.0d0*13.6058d0
       c = ev4au
    case("ry", "Ry", "RY", "rydberg", "Rydberg")
       c = 2.d0
    case("K","kelvin","Kelvin")
       !c = (2.0d0*13.6058d0)*(1.60217733d0/1.380658d0*10000.0d0)
       c = ev4au*ev2kel
    case("ang","Ang","angstrom","Angstrom")
       !c = 0.529177d0
       c = angst
    case("amu")
       !c = 9.109534d-28/1.6605655d-24
       c = 1d0/au_mass
    case("fsec")
       !c = 3.0d0/124.06948d0
       c = 1d0/au_fsec
    case default
       write(*,'(a,a)') '# Error!: unknown unit ', trim(unit)
       stop
    end select
  end function XML_AU_TO

  
end module elses_xml_misc
