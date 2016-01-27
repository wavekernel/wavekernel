!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
program elses_xml_selftest
  use elses_xml_config
  implicit none

  call main
  stop

contains

  ! the main function of this program
  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine main
    use m_dom_debug, only : dom_debug
    use M_options
    implicit none

    call option_default( config%option )
    call elses_process_options

    if( config%option%verbose>0 ) then
       dom_debug = .true.
    else
       dom_debug = .false.
    end if

    if( config%option%filename == "" ) then
       write(*,*) "==================================================="
       write(*,*) "Error! : command_argument function is not working. "
       write(*,*) "==================================================="
       return
    end if

    call config_load( config, config%option%filename )
!    call config_print( config )

    if( &
         config%name /= "selftest" .or. &
         config%system%temperature /= 1.d0 .or. &
         config%system%structure%name /= "selftest" .or. &
         config%system%structure%mdstep /= 0 .or. &
         config%system%structure%unitcell%vectorA(1) /= 10.0 .or. &
         config%system%structure%unitcell%vectorB(2) /= 10.0 .or. &
         config%system%structure%unitcell%vectorC(3) /= 10.0 .or. &
         config%system%structure%heatbath%massperatom /= 25.0 .or. &
         config%system%structure%nelement /= 1 .or. &
         config%system%structure%velement(1)%name /= "C" .or. &
         config%system%structure%velement(1)%classic%charge /= 4.00 .or. &
         config%system%structure%velement(1)%quantum%orbital /= 4 .or. &
         config%system%structure%natom /= 4 .or. &
         config%system%structure%vatom(1)%position(1) /= 0.0 .or. &
         config%system%structure%vatom(1)%velocity(1) /= 0.0 .or. &
         config%system%structure%vatom(2)%position(1) /= 10.0 .or. &
         config%system%structure%vatom(2)%velocity(1) /= 1.0 ) then
       goto 1000 ! error handling block
    end if

    write(*,*) "==================================================="
    write(*,*) "Success! : the XML library is working normally.    "
    write(*,*) "==================================================="
    return

1000 continue ! error handling block

    write(*,*) "==================================================="
    write(*,*) "Error! : the XML libirary is not working correctly."
    write(*,*) "==================================================="
    return
  end subroutine main

end program elses_xml_selftest
