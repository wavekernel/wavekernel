module dos_module
  real(8) :: Edos_min, Edos_max
  integer :: Ndos
  real(8), allocatable :: vTDOS(:) ! total dos, index:idos
  real(8), allocatable :: vPDOS(:,:) ! partial dos, index:idos,iao

  real(8), allocatable :: vITDOS(:) ! integrated TDOS, index:idos
  real(8), allocatable :: vIPDOS(:,:) ! integrated PDOS, index:idos,iao
  real(8), allocatable :: vIETDOS(:) ! integrated E*TDOS, index:idos
  real(8), allocatable :: vIEPDOS(:,:) ! integrated E*PDOS, index:idos,iao
  real(8) :: dEdos

  integer :: Nband, Nao, Nka, Nkb, Nkc
  real(8), allocatable :: vEband(:,:,:,:) ! iband,ika,ikb,ikc
  real(8), allocatable :: vPband(:,:,:,:,:) ! iao,iband,ika,ikb,ikc

  integer, parameter :: ID_METHOD_HISTOGRAM = 1
  integer, parameter :: ID_METHOD_GAUSSIAN  = 2
  integer, parameter :: ID_METHOD_TETRAHEDRON = 3

contains
  integer function dosIndex( E )
    implicit none
    real(8), intent(in) :: E

    dosIndex = int((E-(Edos_min-0.5d0*dEdos))/dEdos) + 1

    return
  end function dosIndex

  real(8) function dosEnergy( idos )
    implicit none
    integer, intent(in) :: idos

    dosEnergy = Edos_min + dEdos *(idos-1)

    return
  end function dosEnergy

  ! return [1:m]
  function modp( n, m ) result(r)
    implicit none
    integer, intent(in) :: n, m
    integer :: r

    r=mod(n-1,m)+1

    return
  end function modp

  
  function histogram( Eband, Edos ) result(dos)
    implicit none
    real(8), intent(in) :: Eband, Edos
    real(8) :: dos
    real(8) :: El, Eu

    El = Eband-0.5d0*dEdos
    Eu = Eband+0.5d0*dEdos

    if( Edos<El ) then
       dos = 0.0d0
    else if( El<=Edos .and. Edos<Eu ) then
       dos = 1.0d0/dEdos
    else
       dos = 0.0d0
    end if

    return
  end function histogram

  function integ_histogram( Eband, Edos ) result(integ_dos)
    implicit none
    real(8), intent(in) :: Eband, Edos
    real(8) :: integ_dos
    real(8) :: El, Eu

    El = Eband-0.5d0*dEdos
    Eu = Eband+0.5d0*dEdos

    if( Edos<El ) then
       integ_dos = 0.0d0
    else if( El<=Edos .and. Edos<Eu ) then
       integ_dos = 1.0d0/dEdos*(Edos-El)
    else
       integ_dos = 1.0d0
    end if

    return
  end function integ_histogram

  function gaussian( Eband, Edos ) result(dos)
    implicit none
    real(8), intent(in) :: Eband, Edos
    real(8) :: dos
    real(8), parameter :: M_SQRT1_2  = 0.70710678118654752440d0 !< sqrt(1/2)
    real(8), parameter :: M_2_SQRTPI = 1.12837916709551257390d0 !< 2/sqrt(Pi)

    dos = 0.5d0*M_SQRT1_2*M_2_SQRTPI/dEdos*exp(-0.5d0*(Edos-Eband)**2/dEdos**2)

    return
  end function gaussian


  function integ_gaussian( Eband, Edos ) result(integ_dos)
    implicit none
    real(8), intent(in) :: Eband, Edos
    real(8) :: integ_dos
    real(8), parameter :: M_SQRT1_2  = 0.70710678118654752440d0 !< sqrt(1/2)
    real(8), parameter :: M_2_SQRTPI = 1.12837916709551257390d0 !< 2/sqrt(Pi)
    real(8) :: derf

    integ_dos = 0.5d0*derf( M_SQRT1_2*(Edos-Eband)/dEdos ) + 0.5d0

    return
  end function integ_gaussian


  function tetrahedron( Eband, Edos ) result(dos)
    implicit none
    real(8), intent(in) :: Eband(4), Edos
    real(8) :: dos

    if( .false. ) then
    else if( Edos < Eband(1) ) then
       dos = 0.0d0
    else if( Edos < Eband(2) ) then
       dos = &
            + (1.0d0/2.0d0)*(Edos-Eband(1))**2 &
            /((Eband(2)-Eband(1))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1)))
    else if( Edos < Eband(3) ) then
       dos = &
            + (1.0d0/2.0d0)*(Eband(4)-Edos)*(Edos-Eband(2)) &
            /((Eband(3)-Eband(2))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1))) &
            + (1.0d0/2.0d0)*(Eband(3)-Edos)*(Edos-Eband(1)) &
            /((Eband(3)-Eband(2))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1)))
    else if( Edos < Eband(4) ) then
       dos = &
            + (1.0d0/2.0d0)*(Eband(4)-Edos)**2 &
            /((Eband(4)-Eband(3))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1)))
    else
       dos = 0.0d0
    end if

    return
  end function tetrahedron

  function integ_tetrahedron( Eband, Edos ) result(integ_dos)
    implicit none
    real(8), intent(in) :: Eband(4), Edos
    real(8) :: integ_dos

    if( .false. ) then
    else if( Edos < Eband(1) ) then
       integ_dos = 0.0d0
    else if( Edos < Eband(2) ) then
       integ_dos = &
            + (1.0d0/6.0d0)*(Edos-Eband(1))**3 &
            /((Eband(2)-Eband(1))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1)))
    else if( Edos < Eband(3) ) then
       integ_dos = &
            + (1.0d0/6.0d0)*(Eband(4)-Edos)*(Edos-Eband(2))**2 &
            /((Eband(3)-Eband(2))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1))) &
            + (1.0d0/6.0d0)*(Eband(3)-Edos)*(Edos-Eband(1))**2 &
            /((Eband(3)-Eband(2))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1))) &
            + (1.0d0/6.0d0)*(Edos-Eband(1))*(Edos-Eband(2)) &
            /((Eband(3)-Eband(2))*(Eband(4)-Eband(1)))
    else if( Edos < Eband(4) ) then
       integ_dos = &
            + (1.0d0/6.0d0)*( 1.0d0 - (Eband(4)-Edos)**3 &
            /((Eband(4)-Eband(3))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1))) )
    else
       integ_dos = 1.0d0/6.0d0
    end if

    return
  end function integ_tetrahedron


  function integ_E3( Ea, Eb, Ec, Ed ) result(integ)
    implicit none
    real(8), intent(in) :: Ea, Eb, Ec, Ed
    real(8) :: integ

    integ = 1.0d0/12.0d0 * &
         (Ed-Ec)*( 3.0d0*(Ec+Ed)**3 - 4.0d0*(Ec*Ec+Ed*Ed+Ec*Ed)*(Ea+Eb) &
         + 6.0d0*(Ea*Eb-Ec*Ed)*(Ec+Ed) )

    return
  end function integ_E3

  function integ_energy_histogram( Eband, Ef ) result(IETDOS)
    implicit none
    real(8), intent(in) :: Eband, Ef
    real(8) :: IETDOS

    real(8) :: El, Eu

    El = Eband-0.5d0*dEdos
    Eu = Eband+0.5d0*dEdos

    if( Ef<El ) then
       IETDOS = 0.0d0
    else if( El<=Ef .and. Ef<Eu ) then
       IETDOS = 1.0d0/dEdos*( Ef**2 - El**2 )
    else
       IETDOS = Eband
    end if

    return
  end function integ_energy_histogram

  function integ_energy_gaussian( Eband, Ef ) result(IETDOS)
    implicit none
    real(8), intent(in) :: Eband, Ef
    real(8) :: IETDOS
    real(8), parameter :: M_SQRT1_2  = 0.70710678118654752440d0 !< sqrt(1/2)
    real(8), parameter :: M_2_SQRTPI = 1.12837916709551257390d0 !< 2/sqrt(Pi)
    real(8) :: derf

    IETDOS = -dEdos*0.5d0*M_SQRT1_2*M_2_SQRTPI*exp(-0.5d0*(Ef-Eband)**2/dEdos**2) &
         + Eband*(0.5d0*derf( M_SQRT1_2*(Ef-Eband)/dEdos ) + 0.5d0)

    return
  end function integ_energy_gaussian

  function integ_energy_tetrahedron( Eband, Ef ) result(IETDOS)
    implicit none
    real(8), intent(in) :: Eband(4), Ef
    real(8) :: IETDOS

    if( .false. ) then
    else if( Ef < Eband(1) ) then
       IETDOS = 0.0d0
    else if( Ef < Eband(2) ) then
       IETDOS = &
            + 0.5d0* ( integ_E3(Eband(1),Eband(1),Eband(1),Ef) ) &
            /((Eband(2)-Eband(1))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1)))

    else if( Ef < Eband(3) ) then
       IETDOS = &
            + 0.5d0* ( integ_E3(Eband(1),Eband(1),Eband(1),Eband(2) ) ) &
            /((Eband(2)-Eband(1))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1))) &
            - 0.5d0* ( integ_E3(Eband(2),Eband(4),Eband(2),Ef) ) &
            /((Eband(3)-Eband(2))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1))) &
            - 0.5d0* ( integ_E3(Eband(1),Eband(3),Eband(2),Ef) ) &
            /((Eband(3)-Eband(2))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1)))

    else if( Ef < Eband(4) ) then
       IETDOS = &
            + 0.5d0* ( integ_E3(Eband(1),Eband(1),Eband(1),Eband(2)) ) &
            /((Eband(2)-Eband(1))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1))) &
            - 0.5d0* ( integ_E3(Eband(2),Eband(4),Eband(2),Eband(3)) ) &
            /((Eband(3)-Eband(2))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1))) &
            - 0.5d0* ( integ_E3(Eband(1),Eband(3),Eband(2),Eband(3)) ) &
            /((Eband(3)-Eband(2))*(Eband(3)-Eband(1))*(Eband(4)-Eband(1))) &
            + 0.5d0* ( integ_E3(Eband(4),Eband(4),Eband(3),Ef) ) &
            /((Eband(4)-Eband(3))*(Eband(4)-Eband(2))*(Eband(4)-Eband(1)))
    else
       IETDOS = &
            + (1.0d0/24.0d0)* ( Eband(1) + Eband(2) + Eband(3) + Eband(4) )
    end if

    return
  end function integ_energy_tetrahedron




  logical function match(fa,fb)
    implicit none
    real(8) :: fa, fb

    if( abs(fa-fb)<=1.0d-13 ) then
       match = .true.
    else
       match = .false.
    end if
  end function match

  logical function unmatch(fa,fb)
    implicit none
    real(8) ::  fa, fb

    unmatch = .not. match(fa,fb)

  end function unmatch


  real(8) function temp_log( a )
    implicit none
    real(8) :: a
    
    if( a<0.0d0 ) then
       temp_log = 1.0d0
    else
       temp_log = 0.0d0
    end if
    
    return
  end function temp_log

  subroutine NOCULC(dos_tetra,Eband,Edos)
    implicit none
    real(8), intent(out) :: dos_tetra(4)
    real(8), intent(in) :: Eband(4), Edos
    real(8) :: dos

    real(8) :: CR(4), f1,f2,f3,f4

    f1 = Edos-Eband(1)
    f2 = Edos-Eband(2)
    f3 = Edos-Eband(3)
    f4 = Edos-Eband(4)

    if( .false. ) then
    else if( unmatch(f1,f2) .and. unmatch(f1,f3) .and. unmatch(f1,f4) &
         .and. unmatch(f2,f3) .and. unmatch(f2,f4) .and. unmatch(f3,f4) ) then
       call NOCULC1(CR,f1,f2,f3,f4)
       ! check
       dos_tetra(1) = -CR(1)
       dos_tetra(2) = -CR(2)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(4)
    else if( match(f1,f2) .and. unmatch(f1,f3) .and. unmatch(f1,f4) .and. unmatch(f3,f4) ) then
       call NOCULC2(CR,f4,f3,f2,f1)
       ! check
       dos_tetra(1) = -CR(4)
       dos_tetra(2) = -CR(3)
       dos_tetra(3) = -CR(2)
       dos_tetra(4) = -CR(1)
    else if( match(f1,f3) .and. unmatch(f1,f2) .and. unmatch(f1,f4) .and. unmatch(f2,f4) ) then
       call NOCULC2(CR,f4,f2,f3,f1)
       ! check
       dos_tetra(1) = -CR(4)
       dos_tetra(2) = -CR(2)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(1)
    else if( match(f1,f4) .and. unmatch(f1,f2) .and. unmatch(f1,f3) .and. unmatch(f2,f3) ) then
       call NOCULC2(CR,f3,f2,f4,f1)
       ! check
       dos_tetra(1) = -CR(3)
       dos_tetra(2) = -CR(2)
       dos_tetra(3) = -CR(1)
       dos_tetra(4) = -CR(4)
    else if( match(f2,f3) .and. unmatch(f1,f2) .and. unmatch(f2,f4) .and. unmatch(f1,f4) ) then
       call NOCULC2(CR,f4,f1,f3,f2)
       ! check
       dos_tetra(1) = -CR(2)
       dos_tetra(2) = -CR(4)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(1)
    else if( match(f2,f4) .and. unmatch(f1,f2) .and. unmatch(f2,f3) .and. unmatch(f1,f3) ) then
       call NOCULC2(CR,f3,f1,f4,f2)
       ! check
       dos_tetra(1) = -CR(2)
       dos_tetra(2) = -CR(4)
       dos_tetra(3) = -CR(1)
       dos_tetra(4) = -CR(3)
    else if( match(f3,f4) .and. unmatch(f1,f3) .and. unmatch(f2,f3) .and. unmatch(f1,f2) ) then
       call NOCULC2(CR,f2,f1,f4,f3)
       ! check
       dos_tetra(1) = -CR(2)
       dos_tetra(2) = -CR(1)
       dos_tetra(3) = -CR(4)
       dos_tetra(4) = -CR(3)
    else if( match(f1,f2) .and. match(f3,f4) .and. unmatch(f1,f3) ) then
       call NOCULC3(CR,f4,f3,f2,f1)
       ! check
       dos_tetra(1) = -CR(4)
       dos_tetra(2) = -CR(3)
       dos_tetra(3) = -CR(2)
       dos_tetra(4) = -CR(1)
    else if( match(f1,f3) .and. match(f2,f4) .and. unmatch(f1,f2) ) then
       call NOCULC3(CR,f4,f2,f3,f1)
       ! check
       dos_tetra(1) = -CR(4)
       dos_tetra(2) = -CR(2)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(1)
    else if( match(f1,f4) .and. match(f2,f3) .and. unmatch(f1,f2) ) then
       call NOCULC3(CR,f3,f2,f4,f1)
       ! check
       dos_tetra(1) = -CR(3)
       dos_tetra(2) = -CR(1)
       dos_tetra(3) = -CR(2)
       dos_tetra(4) = -CR(4)
    else if( match(f1,f2) .and. match(f1,f3) .and. unmatch(f1,f4) ) then
       call NOCULC4(CR,f4,f3,f2,f1)
       ! check
       dos_tetra(1) = -CR(4)
       dos_tetra(2) = -CR(3)
       dos_tetra(3) = -CR(2)
       dos_tetra(4) = -CR(1)
    else if( match(f1,f2) .and. match(f1,f4) .and. unmatch(f1,f3) ) then
       call NOCULC4(CR,f3,f4,f2,f1)
       ! check
       dos_tetra(1) = -CR(3)
       dos_tetra(2) = -CR(4)
       dos_tetra(3) = -CR(1)
       dos_tetra(4) = -CR(2)
    else if( match(f1,f3) .and. match(f1,f4) .and. unmatch(f1,f2) ) then
       call NOCULC4(CR,f2,f4,f3,f1)
       ! check
       dos_tetra(1) = -CR(2)
       dos_tetra(2) = -CR(1)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(4)
    else if( match(f2,f3) .and. match(f2,f4) .and. unmatch(f1,f2) ) then
       call NOCULC4(CR,f1,f4,f3,f2)
       ! check
       dos_tetra(1) = -CR(1)
       dos_tetra(2) = -CR(2)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(4)
    else if( match(f1,f2) .and. match(f1,f3) .and. match(f1,f4) ) then
       call NOCULC5(CR,f1,f2,f3,f4)
       ! check
       dos_tetra(1) = -CR(1)
       dos_tetra(2) = -CR(2)
       dos_tetra(3) = -CR(3)
       dos_tetra(4) = -CR(4)
    end if

  end subroutine NOCULC

  subroutine NOCULC1(CR,f1,f2,f3,f4)
    implicit none
    real(8) :: CR(4), f1,f2,f3,f4
    real(8) :: Ci1,Ci2,Ci3,Ci4
    real(8) :: AL1,AL2,AL3,AL4
    real(8) :: cc1,cc2,cc3,cc4
    real(8) :: cd1,cd2,cd3,cd4,ce1,ce2,ce3,ce4

    !===== for each band (kb)
    if(f1*f2*f3*f4.eq.0.0d0) return
    AL1=temp_log(f1)
    AL2=temp_log(f2)
    AL3=temp_log(f3)
    AL4=temp_log(f4)

    cc1=(f2-f1)*(f3-f1)*(f4-f1)
    cc2=(f1-f2)*(f3-f2)*(f4-f2)
    cc3=(f1-f3)*(f2-f3)*(f4-f3)
    cc4=(f1-f4)*(f2-f4)*(f3-f4)

    cd1=AL1*f1**3
    cd2=AL2*f2**3
    cd3=AL3*f3**3
    cd4=AL4*f4**3

    ce1 = AL1*(f2/(f2-f1)+f3/(f3-f1)+f4/(f4-f1))
    ce2 = AL2*(f3/(f3-f2)+f4/(f4-f2)+f1/(f1-f2))
    ce3 = AL3*(f4/(f4-f3)+f1/(f1-f3)+f2/(f2-f3))
    ce4 = AL4*(f1/(f1-f4)+f2/(f2-f4)+f3/(f3-f4))

    !----- for each vertex (iv)
    !      for vertex 1
    Ci1 = - f1**2*ce1/cc1 &
         - cd2/(cc2*(f2-f1)) - cd3/(cc3*(f3-f1)) - cd4/(cc4*(f4-f1))
    !      for vertex 2
    Ci2 = - f2**2*ce2/cc2 &
         - cd3/(cc3*(f3-f2)) - cd4/(cc4*(f4-f2)) - cd1/(cc1*(f1-f2))
    !      for vertex 3
    Ci3 = - f3**2*ce3/cc3 &
         - cd4/(cc4*(f4-f3)) - cd1/(cc1*(f1-f3)) - cd2/(cc2*(f2-f3))
    !      for vertex 4
    Ci4 = - f4**2*ce4/cc4 &
         - cd1/(cc1*(f1-f4)) - cd2/(cc2*(f2-f4)) - cd3/(cc3*(f3-f4))
    !====   
    CR(1) =  Ci1/6.d0
    CR(2) =  Ci2/6.d0
    CR(3) =  Ci3/6.d0
    CR(4) =  Ci4/6.d0

    return
  end subroutine NOCULC1

  ! f4=f3, ne. f1,f2
  subroutine NOCULC2(CR,f1,f2,f3,f4)
    implicit none
    real(8) :: CR(4), f1,f2,f3,f4

    real(8) :: Ci1,Ci2,Ci3,Ci4
    real(8) :: cf1,cf2,cf3,cf4,cf

    !----- for each vertex (iv)
    ! ---- Factorization of terms temp_log(f4) is essential.
    if(f1*f2*f3*f4.eq.0.0d0) return

    !      for vertex 4 and vertex 3 
    cf1 = 0.0d0
    cf2 = +temp_log(f1)*f1**3/(f1-f2)/((f1-f4)**3)
    cf3 = -temp_log(f2)*f2**3/(f1-f2)/((f2-f4)**3)
    cf4 = +temp_log(f4)*f4*( f1**2*(f2-f4)**2+f2**2*(f1-f4)**2 &
         +f1*f2*(f1-f4)*(f2-f4) ) &
         /((f1-f4)**3*(f2-f4)**3)       
    cf=cf1+cf2+cf3+cf4
    Ci4 = cf
    Ci3 = Ci4

    !      for vertex 2
    Ci2  = temp_log(f1)*f1**3/((f1-f2)**2*(f1-f4)**2)
    cf   = f1*(f2-f4) -2.d0*f4*(f1-f2)
    Ci2  = Ci2 - cf*f2**2*temp_log(f2)/((f2-f4)**3*(f1-f2)**2)
    cf   = 2.d0*f2*(f4-f1)+f1*(f4-f2)
    Ci2  = Ci2 + cf*f4**2*temp_log(f4)/((f2-f4)**3*(f1-f4)**2)

    !      for vertex 1
    Ci1  = temp_log(f2)*f2**3/((f2-f1)**2*(f2-f4)**2)
    cf   = f2*(f1-f4) -2.d0*f4*(f2-f1)
    Ci1  = Ci1 - cf*f1**2*temp_log(f1)/((f1-f4)**3*(f2-f1)**2)
    cf   = 2.d0*f1*(f4-f2)+f2*(f4-f1)
    Ci1  = Ci1 + cf*f4**2*temp_log(f4)/((f1-f4)**3*(f2-f4)**2)
    !====   
    CR(1) =  Ci1/6.d0
    CR(2) =  Ci2/6.d0
    CR(3) =  Ci3/6.d0
    CR(4) =  Ci4/6.d0

    return
  end subroutine NOCULC2

  ! f1=f2 .ne. f3=f4
  subroutine NOCULC3(CR,f1,f2,f3,f4)
    implicit none
    real(8) :: CR(4), f1,f2,f3,f4

    real(8) :: Ci1,Ci2,Ci3,Ci4
    real(8) :: cf

    !----- for each vertex (iv)
    if(f1*f2*f3*f4.eq.0.0d0) return

    !      for vertex 3, vertex4
    Ci3 = - f1**2*f4*temp_log(f1)/(f1-f4)+f1**2*f4*temp_log(f4)/(f1-f4)
    Ci4 = Ci3
    !      for vertex 2, vertex2
    Ci1 = + f1*f4**2*temp_log(f1)/(f1-f4)-f1*f4**2*temp_log(f4)/(f1-f4)
    Ci2=Ci1
    !====   
    CR(1) =  Ci1/(2.d0*(f1-f4)**3)
    CR(2) =  Ci2/(2.d0*(f1-f4)**3)
    CR(3) =  Ci3/(2.d0*(f1-f4)**3)
    CR(4) =  Ci4/(2.d0*(f1-f4)**3)

    return
  end subroutine NOCULC3

  ! f4=f3=f2 .ne. f1
  subroutine NOCULC4(CR,f1,f2,f3,f4)
    implicit none
    real(8) :: CR(4), f1,f2,f3,f4

    real(8) :: Ci1,Ci2,Ci3,Ci4
    real(8) :: cf

    !----- for each vertex (iv)
    if(f1*f2*f3*f4.eq.0.0d0) return
    !      for vertex4, vertex3, and vertex f2
    Ci4 = + f1**3*temp_log(f1)/(f1-f4)-f1**3*temp_log(f4)/(f1-f4)
    Ci3 = Ci4
    Ci2 = Ci4
    !      for vertex 1
    Ci1 = - 3.d0*f1**2*f4*temp_log(f1)/(f1-f4)+3.d0*f1**2*f4*temp_log(f4)/(f1-f4)
    !====   
    CR(1) =  Ci1/(6.d0*(f1-f4)**3)
    CR(2) =  Ci2/(6.d0*(f1-f4)**3)
    CR(3) =  Ci3/(6.d0*(f1-f4)**3)
    CR(4) =  Ci4/(6.d0*(f1-f4)**3)

    return
  end subroutine NOCULC4

  ! f1=f2=f3=f4
  subroutine NOCULC5(CR,f1,f2,f3,f4)
    implicit none
    real(8) :: CR(4), f1,f2,f3,f4

    real(8) :: Ci1,Ci2,Ci3,Ci4

    !----- for each vertex (iv)
    if(f1*f2*f3*f4.eq.0.0d0) return

    !      for vertex 1, 2, 3, 4
    Ci1 = 0.0d0
    Ci2 = 0.0d0
    Ci3 = 0.0d0
    Ci4 = 0.0d0
    !====   
    CR(1) =  Ci1
    CR(2) =  Ci2
    CR(3) =  Ci3
    CR(4) =  Ci4

    return
  end subroutine NOCULC5

end module dos_module

program dos
  use dos_module
  implicit none


  ! ===== Variables for command line options
  character(len=32) :: argv
  integer :: id_method
  real(8) :: Nval, Efer, Eorb

  ! ===== Variables read from "ForDosCalc.txt"
  character(len=*),parameter :: fileNameParms  = "ForDosCalc.txt"

  ! ===== Variables read from "EigenEnergy.txt"
  character(len=*),parameter :: fileNameEnergy = "EigenEnergy.txt"

  ! ===== Variables read from "EigenWave.txt"
  character(len=*),parameter :: fileNameWave   = "EigenWave.txt"


  ! ===== output file 
  character(len=*),parameter :: fileNameOutput= "DOS.txt"

  Nval = 0.0d0
  if( command_argument_count() >= 1 ) then
     call get_command_argument(1,argv)
     select case(argv)
     case("-histogram","-h")
        id_method = ID_METHOD_HISTOGRAM
     case("-gaussian","-g")
        id_method = ID_METHOD_GAUSSIAN
     case("-tetrahedron","-t")
        id_method = ID_METHOD_TETRAHEDRON
     case default
        write(*,*) "unknown argument:", trim(argv)
        write(*,*) "use argument -histogram or -gaussian or -tetrahedron"
        stop
     end select

     if( command_argument_count() == 2 ) then
        call get_command_argument(2,argv)
        read(argv,*) Nval
     end if
  else
     id_method = ID_METHOD_HISTOGRAM
  end if

  call readParam(fileNameParms)

  call readEnergy(fileNameEnergy)
  call readWave(fileNameWave)

  call calcDOS(id_method)
  call calcIDOS

  if( Nval /= 0.0d0 ) then
     call calcFermi(Efer,Nval,id_method)
     call calcEorbital(Eorb,Efer,id_method)
  end if

  write(*,'(a,f12.6)') "#  Nelec    = ", Nval
  write(*,'(a,f12.6)') "#  Efermi   = ", Efer
  write(*,'(a,f12.6)') "#  Eorbital = ", Eorb


  call outputDOS(fileNameOutput)

  call clearDOS

end program dos


subroutine readParam(fileName)
  use dos_module
  implicit none
  character(len=*), intent(in) :: fileName
  integer, parameter    :: unit=10

  !----- read some parameters from "ForBandCalc.txt" ----------------------
  open(unit,file=fileName)
  read(unit,*)                         ! skip comment
  read(unit,*) Edos_min, Edos_max
  read(unit,*)                         ! skip comment
  read(unit,*) Ndos
  close(unit)

  dEdos = (Edos_max-Edos_min)/(Ndos-1)
  allocate( vTDOS(Ndos) )
  vTDOS(:) = 0.0d0

  return
end subroutine readParam

subroutine readEnergy(fileName)
  use dos_module
  implicit none
  character(len=*), intent(in) :: fileName
  integer, parameter    :: unit=10

  integer :: ika, ikb, ikc
  integer :: dumi1,dumi2,dumi3
  real(8) :: dumr1,dumr2,dumr3

  open(unit, file=fileName, status='old')
  read(unit,*) Nka, Nkb, Nkc, Nband
  allocate( vEband(Nband,Nka,Nkb,Nkc) )
  do ika=1, Nka
     do ikb=1, Nkb
        do ikc=1, Nkc
           read(unit,*) dumi1,dumi2,dumi3,dumr1,dumr2,dumr3,vEband(:,ika,ikb,ikc)
        end do
     end do
  end do
  close(unit)

  return
end subroutine readEnergy

subroutine readWave(fileName)
  use dos_module
  implicit none
  character(len=*), intent(in) :: fileName
  integer, parameter    :: unit=10

  integer :: ika, ikb, ikc, m
  integer :: dumi1,dumi2,dumi3
  real(8) :: dumr1

  open(unit, file=fileName, status='old')
  read(unit,*) Nka, Nkb, Nkc, Nband, Nao
  allocate( vPband(Nao,Nband,Nka,Nkb,Nkc) )

  allocate( vPDOS(Ndos,Nao) )
  vPDOS(:,:) = 0.0d0

  do ika=1, Nka
     do ikb=1, Nkb
        do ikc=1, Nkc
           do m=1, Nband
              read(unit,*) dumi1,dumi2,dumi3,dumr1,vPband(:,m,ika,ikb,ikc)
           end do
        end do
     end do
  end do
  close(unit)

  return
end subroutine readWave


subroutine calcDOS( id_method )
  use dos_module
  implicit none
  integer, intent(in) :: id_method

  select case(id_method)
  case(ID_METHOD_HISTOGRAM)
     call calcDOS_Histogram
  case(ID_METHOD_GAUSSIAN)
     call calcDOS_Gaussian
  case(ID_METHOD_TETRAHEDRON)
     call calcDOS_Tetrahedron
  end select

  return
end subroutine calcDOS

subroutine calcFermi( Efer, Nval, id_method )
  use dos_module
  implicit none
  real(8), intent(out) :: Efer
  real(8), intent(in) :: Nval
  integer, intent(in) :: id_method

  select case(id_method)
  case(ID_METHOD_HISTOGRAM)
     call calcFermi_Histogram(Efer,Nval)
  case(ID_METHOD_GAUSSIAN)
     call calcFermi_Gaussian(Efer,Nval)
  case(ID_METHOD_TETRAHEDRON)
     call calcFermi_Tetrahedron(Efer,Nval)
  end select

  return
end subroutine calcFermi

subroutine calcEorbital( Eorb, Efer, id_method )
  use dos_module
  implicit none
  real(8), intent(out) :: Eorb
  real(8), intent(in) :: Efer
  integer, intent(in) :: id_method

  select case(id_method)
  case(ID_METHOD_HISTOGRAM)
     call calcEorbital_Histogram(Eorb,Efer)
  case(ID_METHOD_GAUSSIAN)
     call calcEorbital_Gaussian(Eorb,Efer)
  case(ID_METHOD_TETRAHEDRON)
     call calcEorbital_Tetrahedron(Eorb,Efer)
  end select

  return
end subroutine calcEorbital

subroutine calcDOS_Histogram
  use dos_module
  implicit none

  integer :: ika, ikb, ikc, iband
  integer :: idos
  real(8) :: Eband, dos
  real(8), allocatable :: Pwave(:)
  allocate( Pwave(Nao) )

  do iband=1, Nband
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              Eband = vEband(iband,ika,ikb,ikc)
              Pwave = vPband(:,iband,ika,ikb,ikc)
              idos  = dosIndex(Eband)

              if( 0<idos .and. idos<=Ndos ) then
                 dos = histogram(Eband,dosEnergy(idos))
                 vTDOS(idos) = vTDOS(idos) + dos
                 vPDOS(idos,:) = vPDOS(idos,:) + dos*Pwave(:)
              end if
           end do ! ika
        end do ! ikb
     end do ! ikc
  end do ! iband

  vTDOS(:) = vTDOS(:)*(1.0d0/(Nka*Nkb*Nkc))
  vPDOS(:,:) = vPDOS(:,:)*(1.0d0/(Nka*Nkb*Nkc))

  deallocate(Pwave)

  return
end subroutine calcDOS_Histogram


subroutine calcFermi_Histogram( Ef, Nele )
  use dos_module
  implicit none

  real(8), intent(out) :: Ef
  real(8), intent(in) :: Nele

  integer :: ika, ikb, ikc, iband
  real(8) :: Eband, integ_tdos
  real(8) :: Na, Nb
  real(8) :: Ea, Eb

  ! set the first trial range of energy
  Ea = minval(vEband)
  Eb = maxval(vEband)
  Ef = (Ea+Eb)*0.5d0

  do
     integ_tdos = 0.0d0
     do iband=1, Nband
        do ika=1, Nka
           do ikb=1, Nkb
              do ikc=1, Nkc
                 Eband = vEband(iband,ika,ikb,ikc)

                 integ_tdos = integ_tdos + integ_histogram(Eband,Ef)
              end do ! ika
           end do ! ikb
        end do ! ikc
     end do ! iband
     integ_tdos = integ_tdos*(1.0d0/(Nka*Nkb*Nkc))

     if( integ_tdos > Nele ) then
        Eb = Ef
        Ef = (Ea+Eb)*0.5d0
     else
        Ea = Ef
        Ef = (Ea+Eb)*0.5d0
     end if
     if( Eb-Ea<1.0d-8 ) exit
  end do

  return
end subroutine calcFermi_Histogram

subroutine calcEorbital_Histogram( Eorb, Ef )
  use dos_module
  implicit none

  real(8), intent(out) :: Eorb
  real(8), intent(in) :: Ef

  integer :: ika, ikb, ikc, iband
  real(8) :: Eband, integ_etdos

  integ_etdos = 0.0d0
  do iband=1, Nband
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              Eband = vEband(iband,ika,ikb,ikc)

              integ_etdos = integ_etdos + integ_energy_histogram(Eband,Ef)
           end do ! ika
        end do ! ikb
     end do ! ikc
  end do ! iband
  integ_etdos = integ_etdos*(1.0d0/(Nka*Nkb*Nkc))

  Eorb = integ_etdos

  return
end subroutine calcEorbital_Histogram

subroutine calcDOS_Gaussian
  use dos_module
  implicit none

  integer :: ika, ikb, ikc, iband
  integer :: idos, idos_s, idos_e
  real(8) :: Eband, dos
  real(8), allocatable :: Pwave(:)
  allocate( Pwave(Nao) )

  do iband=1, Nband
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              Eband = vEband(iband,ika,ikb,ikc)
              Pwave = vPband(:,iband,ika,ikb,ikc)

              idos_s = dosIndex(Eband-dEdos)
              idos_e = dosIndex(Eband+dEdos)

              do idos=idos_s, idos_e
                 if( 0<idos .and. idos<=Ndos ) then
                    dos = gaussian(Eband,dosEnergy(idos))
                    vTDOS(idos) = vTDOS(idos) + dos
                    vPDOS(idos,:) = vPDOS(idos,:) + dos*Pwave(:)
                 end if
              end do ! idos
           end do ! ika
        end do ! ikb
     end do ! ikc
  end do ! iband

  vTDOS(:) = vTDOS(:)*(1.0d0/(Nka*Nkb*Nkc))
  vPDOS(:,:) = vPDOS(:,:)*(1.0d0/(Nka*Nkb*Nkc))

  deallocate(Pwave)

  return
end subroutine calcDOS_Gaussian



subroutine calcFermi_Gaussian( Ef, Nele )
  use dos_module
  implicit none

  real(8), intent(out) :: Ef
  real(8), intent(in) :: Nele

  integer :: ika, ikb, ikc, iband
  real(8) :: Eband, integ_tdos
  real(8) :: Na, Nb
  real(8) :: Ea, Eb

  ! set the first trial range of energy
  Ea = minval(vEband)
  Eb = maxval(vEband)
  Ef = (Ea+Eb)*0.5d0

  do
     integ_tdos = 0.0d0
     do iband=1, Nband
        do ika=1, Nka
           do ikb=1, Nkb
              do ikc=1, Nkc
                 Eband = vEband(iband,ika,ikb,ikc)

                 integ_tdos = integ_tdos + integ_gaussian(Eband,Ef)
              end do ! ika
           end do ! ikb
        end do ! ikc
     end do ! iband
     integ_tdos = integ_tdos*(1.0d0/(Nka*Nkb*Nkc))

     if( integ_tdos > Nele ) then
        Eb = Ef
        Ef = (Ea+Eb)*0.5d0
     else
        Ea = Ef
        Ef = (Ea+Eb)*0.5d0
     end if
     if( Eb-Ea<1.0d-8 ) exit
  end do

  return
end subroutine calcFermi_Gaussian


subroutine calcEorbital_Gaussian( Eorb, Ef )
  use dos_module
  implicit none

  real(8), intent(out) :: Eorb
  real(8), intent(in) :: Ef

  integer :: ika, ikb, ikc, iband
  real(8) :: Eband, integ_etdos

  integ_etdos = 0.0d0
  do iband=1, Nband
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              Eband = vEband(iband,ika,ikb,ikc)

              integ_etdos = integ_etdos + integ_energy_gaussian(Eband,Ef)
           end do ! ika
        end do ! ikb
     end do ! ikc
  end do ! iband
  integ_etdos = integ_etdos*(1.0d0/(Nka*Nkb*Nkc))

  Eorb = integ_etdos

  return
end subroutine calcEorbital_Gaussian





subroutine calcDOS_Tetrahedron
  use dos_module
  implicit none

  integer :: ika, ikb, ikc, iband
  integer :: idos, idos_s, idos_e
  integer :: t, n, m, j
  integer :: jka, jkb, jkc
  real(8) :: Eband(4), tdos, integ_tdos, work, dos_tetra(4)
  real(8), allocatable :: Pwave(:,:)
  integer :: index(4,6)

  index(:,1) = (/ 0,1,3,7 /)
  index(:,2) = (/ 0,1,5,7 /)
  index(:,3) = (/ 0,4,5,7 /)
  index(:,4) = (/ 0,4,6,7 /)
  index(:,5) = (/ 0,2,6,7 /)
  index(:,6) = (/ 0,2,3,7 /)

  allocate( Pwave(Nao,4) )


  do iband=1, Nband
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              do t=1, 6
                 do n=1, 4
                    j = index(n,t)
                    jka = modp(ika+mod(j/1,2),Nka)
                    jkb = modp(ikb+mod(j/2,2),Nkb)
                    jkc = modp(ikc+mod(j/4,2),Nkc)

                    Eband(n)  = vEband(iband,jka,jkb,jkc)
                    Pwave(:,n)  = vPband(:,iband,jka,jkb,jkc)
                 end do

                 do n=1, 4
                    do m=n+1, 4
                       if( Eband(n) > Eband(m) ) then
                          work = Eband(n)
                          Eband(n) = Eband(m)
                          Eband(m) = work
                       end if
                    end do
                 end do

                 idos_s = dosIndex(Eband(1))
                 idos_e = dosIndex(Eband(4))

                 do idos=idos_s, idos_e
                    if( 0<idos .and. idos<=Ndos ) then
                       tdos = tetrahedron(Eband,dosEnergy(idos))
                       vTDOS(idos) = vTDOS(idos) + tdos

                       call NOCULC(dos_tetra,Eband,dosEnergy(idos))

                       do n=1, 4
                          vPDOS(idos,:) = vPDOS(idos,:) + dos_tetra(n)*Pwave(:,n)
                       end do
                    end if
                 end do ! idos

              end do ! t
           end do ! ika
        end do ! ikb
     end do ! ikc
  end do ! iband

!!$  do idos=1, Ndos
!!$     vTDOS(idos) = sum(vPDOS(idos,:))
!!$  end do

  vTDOS(:) = vTDOS(:)*(1.0d0/(Nka*Nkb*Nkc))
  vPDOS(:,:) = vPDOS(:,:)*(1.0d0/(Nka*Nkb*Nkc))

  deallocate(Pwave)

  return
end subroutine calcDOS_Tetrahedron


subroutine calcFermi_Tetrahedron( Ef, Nele )
  use dos_module
  implicit none

  real(8), intent(out) :: Ef
  real(8), intent(in) :: Nele

  integer :: ika, ikb, ikc, iband
  integer :: t, n, m, j
  integer :: jka, jkb, jkc
  real(8) :: Eband(4), work, integ_tdos
  integer :: index(4,6)
  real(8) :: Na, Nb
  real(8) :: Ea, Eb

  index(:,1) = (/ 0,1,3,7 /)
  index(:,2) = (/ 0,1,5,7 /)
  index(:,3) = (/ 0,4,5,7 /)
  index(:,4) = (/ 0,4,6,7 /)
  index(:,5) = (/ 0,2,6,7 /)
  index(:,6) = (/ 0,2,3,7 /)

  ! set the first trial range of energy
  Ea = minval(vEband)
  Eb = maxval(vEband)
  Ef = (Ea+Eb)*0.5d0

!!$  write(*,'(a)') "# Searching fermi energy"
!!$  write(*,'(a)') "#   Ea,         Ef,         Eb,         ITDOS(Ef),  Nval"
  do
     integ_tdos = 0.0d0
     do iband=1, Nband
        do ika=1, Nka
           do ikb=1, Nkb
              do ikc=1, Nkc
                 do t=1, 6
                    do n=1, 4
                       j = index(n,t)
                       jka = modp(ika+mod(j/1,2),Nka)
                       jkb = modp(ikb+mod(j/2,2),Nkb)
                       jkc = modp(ikc+mod(j/4,2),Nkc)

                       Eband(n)  = vEband(iband,jka,jkb,jkc)
                    end do

                    do n=1, 4
                       do m=n+1, 4
                          if( Eband(n) > Eband(m) ) then
                             work = Eband(n)
                             Eband(n) = Eband(m)
                             Eband(m) = work
                          end if
                       end do
                    end do

                    integ_tdos = integ_tdos + integ_tetrahedron(Eband,Ef)
                 end do ! t
              end do ! ika
           end do ! ikb
        end do ! ikc
     end do ! iband
     integ_tdos = integ_tdos*(1.0d0/(Nka*Nkb*Nkc))

!!$     write(*,"(5f12.6)") Ea, Ef, Eb, integ_tdos, Nele
     if( integ_tdos > Nele ) then
        Eb = Ef
        Ef = (Ea+Eb)*0.5d0
     else
        Ea = Ef
        Ef = (Ea+Eb)*0.5d0
     end if
     if( Eb-Ea<1.0d-6 ) exit
  end do


  return
end subroutine calcFermi_Tetrahedron




subroutine calcEorbital_Tetrahedron( Eorb, Ef )
  use dos_module
  implicit none

  real(8), intent(out) :: Eorb
  real(8), intent(in) :: Ef

  integer :: ika, ikb, ikc, iband
  integer :: t, n, m, j
  integer :: jka, jkb, jkc
  real(8) :: Eband(4), work, integ_etdos
  integer :: index(4,6)


  index(:,1) = (/ 0,1,3,7 /)
  index(:,2) = (/ 0,1,5,7 /)
  index(:,3) = (/ 0,4,5,7 /)
  index(:,4) = (/ 0,4,6,7 /)
  index(:,5) = (/ 0,2,6,7 /)
  index(:,6) = (/ 0,2,3,7 /)

  integ_etdos = 0.0d0
  do iband=1, Nband
     do ika=1, Nka
        do ikb=1, Nkb
           do ikc=1, Nkc
              do t=1, 6
                 do n=1, 4
                    j = index(n,t)
                    jka = modp(ika+mod(j/1,2),Nka)
                    jkb = modp(ikb+mod(j/2,2),Nkb)
                    jkc = modp(ikc+mod(j/4,2),Nkc)

                    Eband(n)  = vEband(iband,jka,jkb,jkc)
                 end do

                 do n=1, 4
                    do m=n+1, 4
                       if( Eband(n) > Eband(m) ) then
                          work = Eband(n)
                          Eband(n) = Eband(m)
                          Eband(m) = work
                       end if
                    end do
                 end do

                 integ_etdos = integ_etdos + integ_energy_tetrahedron(Eband,Ef)
              end do ! t
           end do ! ika
        end do ! ikb
     end do ! ikc
  end do ! iband
  integ_etdos = integ_etdos*(1.0d0/(Nka*Nkb*Nkc))

  Eorb = integ_etdos

  return
end subroutine calcEorbital_Tetrahedron



subroutine calcIDOS
  use dos_module
  implicit none

  real(8) :: sum1, sumE, E
  integer :: idos, iao

  allocate( vITDOS(Ndos) )
  allocate( vIPDOS(Ndos,Nao) )
  allocate( vIETDOS(Ndos) )
  allocate( vIEPDOS(Ndos,Nao) )

  sum1 = 0.0d0
  sumE = 0.0d0
  do idos=1, Ndos
     E = dosEnergy(idos)
     sum1 = sum1 +   vTDOS(idos)*dEdos
     sumE = sumE + E*vTDOS(idos)*dEdos
     vITDOS(idos)  = sum1
     vIETDOS(idos) = sumE
  end do

  do iao=1, Nao
     sum1 = 0.0d0
     sumE = 0.0d0
     do idos=1, Ndos
        E = dosEnergy(idos)
        sum1 = sum1 +   vPDOS(idos,iao)*dEdos
        sumE = sumE + E*vPDOS(idos,iao)*dEdos
        vIPDOS(idos,iao)  = sum1
        vIEPDOS(idos,iao) = sumE
     end do
  end do

  return
end subroutine calcIDOS


subroutine outputDOS(fileName)
  use dos_module
  implicit none
  character(len=*), intent(in) :: fileName
  integer, parameter    :: unit=10
  integer :: idos, iao

  open(unit,file=fileName)

  write(unit,*) "index, energy, ", &
       "TDOS, PDOS..., ", &
       "integrated TDOS,integrated PDOS...,", &
       "integrated E*DOS,integrated E*PDOS...,"
       
  do idos=1, Ndos
     write(unit,'(1I5,1E12.3$)') idos, dosEnergy(idos)

     write(unit,'(E15.5$)') vTDOS(idos)
     do iao=1, Nao
        write(unit,'(E15.5$)') vPDOS(idos,iao)
     end do

     write(unit,'(E15.5$)') vITDOS(idos)
     do iao=1, Nao
        write(unit,'(E15.5$)') vIPDOS(idos,iao)
     end do

     write(unit,'(E15.5$)') vIETDOS(idos)
     do iao=1, Nao
        write(unit,'(E15.5$)') vIEPDOS(idos,iao)
     end do

     write(unit,*)
  end do

  close(unit)
  return
end subroutine outputDOS



subroutine clearDOS
  use dos_module
  implicit none

  deallocate( vTDOS )
  deallocate( vPDOS )
  deallocate( vEband )
  deallocate( vPband )

  deallocate( vITDOS )
  deallocate( vIPDOS )
  deallocate( vIETDOS )
  deallocate( vIEPDOS )

  return
end subroutine clearDOS
