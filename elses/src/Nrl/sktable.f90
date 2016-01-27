!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module MSlaterKosterTable
  implicit none

!  private
  public TSKBond
  public SlaterKosterCoefficient

  !!
  ! v(1): sigma bonding
  ! v(2): pi bonding
  ! v(3): delta bonding
  ! v(4): phi bonding
  !
  type TSKBond
    real(8) :: v(4)
  end type

contains

  double precision function Vss(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vss = bond%v(1)
  end function

  double precision function Vsx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vsx = l*bond%v(1)
  end function

  double precision function Vsy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision l, m, n

    Vsy = Vsx(bond, m, n, l)
  end function

  double precision function Vsz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision l, m, n

    Vsz = Vsx(bond, n, l, m)
  end function

  double precision function Vxx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vxx = l*l*bond%v(1) + (1 - l*l)*bond%v(2)
  end function

  double precision function Vxy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vxy = l*m*(bond%v(1) - bond%v(2))
  end function

  double precision function Vxz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vxz = l*n*(bond%v(1) - bond%v(2))
  end function

  double precision function Vyy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vyy = Vxx(bond, m, n, l)
  end function

  double precision function Vyz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vyz = Vxy(bond, m, n, l)
  end function

  double precision function Vzz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vzz = Vxx(bond, n, l, m)
  end function

  double precision function Vs_xy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n
    
    Vs_xy = dsqrt(3.0d0)*l*m*bond%v(1)
  end function

  double precision function Vs_yz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vs_yz = Vs_xy(bond, m, n, l)
  end function

  double precision function Vs_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vs_zx = Vs_xy(bond, n, l, m)
  end function

  double precision function Vs_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vs_xx = dsqrt(3d0)*0.5d0*(l**2-m**2)*bond%v(1)
  end function

  double precision function Vs_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vs_zz = (n**2 - 0.5d0*(l**2+m**2))*bond%v(1)
  end function

  double precision function Vx_xy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vx_xy = dsqrt(3d0)*l*l*m*bond%v(1) + &
              m*(1-2*l**2)*bond%v(2)
  end function

  double precision function Vx_yz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vx_yz = dsqrt(3d0)*l*m*n*bond%v(1) - &
                   2*l*m*n*bond%v(2)
  end function

  double precision function Vx_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vx_zx = dsqrt(3d0)*l**2*n*bond%v(1) + &
               n*(1-2*l**2)*bond%v(2)
  end function

  double precision function Vy_xy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vy_xy = Vx_zx(bond, m, n, l)
  end function

  double precision function Vy_yz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vy_yz = Vx_xy(bond, m, n, l)
  end function

  double precision function Vy_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vy_zx = Vx_yz(bond, m, n, l)
  end function

  double precision function Vz_xy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vz_xy = Vx_yz(bond, n, l, m)
  end function

  double precision function Vz_yz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vz_yz = Vx_zx(bond, n, l, m)
  end function

  double precision function Vz_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l, m, n

    Vz_zx = Vx_xy(bond, n, l, m)
  end function

  double precision function Vx_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vx_xx = dsqrt(3d0)*0.5d0*l*(l**2-m**2)*bond%v(1) + &
                          l*(1-l**2+m**2)*bond%v(2)
  end function

  double precision function Vy_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vy_xx = dsqrt(3d0)*0.5d0*m*(l**2-m**2)*bond%v(1) &
                    -m*(1 + l**2 - m**2)*bond%v(2)
  end function

  double precision function Vz_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vz_xx = dsqrt(3d0)*0.5d0*n*(l**2-m**2)*bond%v(1) &
                          -n*(l**2-m**2)*bond%v(2)
  end function

  double precision function Vx_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vx_zz = l*(n**2 - 0.5d0*(l**2+m**2))*bond%v(1) &
                        - dsqrt(3d0)*l*n**2*bond%v(2)
  end function 

  double precision function Vy_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vy_zz = m*(n**2 - 0.5d0*(l**2+m**2))*bond%v(1) &
                        - dsqrt(3d0)*m*n**2*bond%v(2)
  end function 

  double precision function Vz_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision :: l,m,n

    Vz_zz = n*(n**2 - 0.5d0*(l**2+m**2))*bond%v(1) &
                + dsqrt(3d0)*n*(l**2+m**2)*bond%v(2)
  end function

  double precision function Vxy_xy(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxy_xy = 3*l**2*m**2*bond%v(1) &
            +(l**2+m**2-4*l**2*m**2)*bond%v(2) &
            +(n**2+l**2*m**2)*bond%v(3)
  end function

  double precision function Vxy_yz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxy_yz = 3*l*m**2*n*bond%v(1) &
            +l*n*(1-4*m**2)*bond%v(2) &
            +l*n*(m**2-1)*bond%v(3)
  end function

  double precision function Vxy_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxy_zx = 3*l**2*m*n*bond%v(1) &
            +m*n*(1-4*l**2)*bond%v(2) &
            +m*n*(l**2-1)*bond%v(3)
  end function

  double precision function Vyz_yz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vyz_yz = Vxy_xy(bond, m, n, l)
  end function

  double precision function Vyz_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vyz_zx = Vxy_yz(bond, m, n, l)
  end function

  double precision function Vzx_zx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vzx_zx = Vxy_xy(bond, n, l, m)
  end function

  double precision function Vxy_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxy_xx = 3.0d0/2.0d0*l*m*(l**2-m**2)*bond%v(1) &
            +2*l*m*(m**2-l**2)*bond%v(2) &
            +1/2.0d0*l*m*(l**2-m**2)*bond%v(3)
  end function

  double precision function Vyz_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vyz_xx = 3.0d0/2.0d0*m*n*(l**2-m**2)*bond%v(1) &
            -m*n*(1+2*(l**2-m**2))*bond%v(2) &
            +m*n*(1+1/2.0d0*(l**2-m**2))*bond%v(3)
  end function

  double precision function Vzx_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vzx_xx = 3.0d0/2.0d0*n*l*(l**2-m**2)*bond%v(1) &
            +n*l*(1-2*(l**2-m**2))*bond%v(2) &
            -n*l*(1-1/2.0d0*(l**2-m**2))*bond%v(3)
  end function

  double precision function Vxy_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxy_zz = dsqrt(3.0d0)*l*m*(n**2-1.0d0/2*(l**2+m**2))*bond%v(1) &
            -2*dsqrt(3d0)*l*m*n**2*bond%v(2) &
            +dsqrt(3.0d0)/2*l*m*(1+n**2)*bond%v(3)
  end function

  double precision function Vyz_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vyz_zz = dsqrt(3.0d0)*m*n*(n**2-1.0d0/2*(l**2+m**2))*bond%v(1) &
            +dsqrt(3.0d0)*m*n*(l**2+m**2-n**2)*bond%v(2) &
            -dsqrt(3.0d0)/2*m*n*(l**2+m**2)*bond%v(3)
  end function

  double precision function Vzx_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vzx_zz = dsqrt(3.0d0)*n*l*(n**2-1.0d0/2*(l**2+m**2))*bond%v(1) &
            +dsqrt(3.0d0)*n*l*(l**2+m**2-n**2)*bond%v(2) &
            -dsqrt(3.0d0)/2*n*l*(l**2+m**2)*bond%v(3)
  end function

  double precision function Vxx_xx(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxx_xx = 3.0d0/4*(l**2-m**2)**2*bond%v(1) &
            +(l**2+m**2-(l**2-m**2)**2)*bond%v(2) &
            +(n**2+1.0d0/4*(l**2-m**2)**2)*bond%v(3)
  end function

  double precision function Vxx_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vxx_zz = dsqrt(3.0d0)/2*(l**2-m**2)*(n**2-1.0d0/2*(l**2+m**2)) &
                                       *bond%v(1) &
            +dsqrt(3.0d0)*n**2*(m**2-l**2)*bond%v(2) &
            +dsqrt(3.0d0)/4*(1+n**2)*(l**2-m**2)*bond%v(3)
  end function

  double precision function Vzz_zz(bond, l, m, n)
    type(TSKBond) :: bond
    double precision     :: l, m, n

    Vzz_zz = (n**2-1.0d0/2*(l**2+m**2))**2*bond%v(1) &
            +3*n**2*(l**2+m**2)*bond%v(2) &
            +3.0d0/4*(l**2+m**2)**2*bond%v(3)
  end function

  !!
  ! Return Slater Koster coefficient.
  ! Note that Vsx = Vxs and so on.
  !
  real(8) function SlaterKosterCoefficient(ibond, lm1, lm2, l, m, n) &
       result (oresult)
    type(TSKBond), intent(in) :: ibond
    integer, intent(in) :: lm1, lm2
    real(8), intent(in) :: l, m, n

    select case(lm1)
    case(0)
      select case(lm2)
      case(0)
        oresult = Vss(ibond, l, m, n)
      case(1)
        oresult = Vsx(ibond, l, m, n)
      case(2)
        oresult = Vsy(ibond, l, m, n)
      case(3)
        oresult = Vsz(ibond, l, m, n)
      case(4)
        oresult = Vs_xy(ibond, l, m, n)
      case(5)
        oresult = Vs_yz(ibond, l, m, n)
      case(6)
        oresult = Vs_zx(ibond, l, m, n)
      case(7)
        oresult = Vs_xx(ibond, l, m, n)
      case(8)
        oresult = Vs_zz(ibond, l, m, n)
      end select
    case(1)
      select case(lm2)
      case(0)
        oresult = Vsx(ibond, l, m, n)
      case(1)
        oresult = Vxx(ibond, l, m, n)
      case(2)
        oresult = Vxy(ibond, l, m, n)
      case(3)
        oresult = Vxz(ibond, l, m, n)
      case(4)
        oresult = Vx_xy(ibond, l, m, n)
      case(5)
        oresult = Vx_yz(ibond, l, m, n)
      case(6)
        oresult = Vx_zx(ibond, l, m, n)
      case(7)
        oresult = Vx_xx(ibond, l, m, n)
      case(8)
        oresult = Vx_zz(ibond, l, m, n)
      end select
    case(2)
      select case(lm2)
      case(0)
        oresult = Vsy(ibond, l, m, n)
      case(1)
        oresult = Vxy(ibond, l, m, n)
      case(2)
        oresult = Vyy(ibond, l, m, n)
      case(3)
        oresult = Vyz(ibond, l, m, n)
      case(4)
        oresult = Vy_xy(ibond, l, m, n)
      case(5)
        oresult = Vy_yz(ibond, l, m, n)
      case(6)
        oresult = Vy_zx(ibond, l, m, n)
      case(7)
        oresult = Vy_xx(ibond, l, m, n)
      case(8)
        oresult = Vy_zz(ibond, l, m, n)
      end select
    case(3)
      select case(lm2)
      case(0)
        oresult = Vsz(ibond, l, m, n)
      case(1)
        oresult = Vxz(ibond, l, m, n)
      case(2)
        oresult = Vyz(ibond, l, m, n)
      case(3)
        oresult = Vzz(ibond, l, m, n)
      case(4)
        oresult = Vz_xy(ibond, l, m, n)
      case(5)
        oresult = Vz_yz(ibond, l, m, n)
      case(6)
        oresult = Vz_zx(ibond, l, m, n)
      case(7)
        oresult = Vz_xx(ibond, l, m, n)
      case(8)
        oresult = Vz_zz(ibond, l, m, n)
      end select
    case(4)
      select case(lm2)
      case(0)
        oresult = Vs_xy(ibond, l, m, n)
      case(1)
        oresult = Vx_xy(ibond, l, m, n)
      case(2)
        oresult = Vy_xy(ibond, l, m, n)
      case(3)
        oresult = Vz_xy(ibond, l, m, n)
      case(4)
        oresult = Vxy_xy(ibond, l, m, n)
      case(5)
        oresult = Vxy_yz(ibond, l, m, n)
      case(6)
        oresult = Vxy_zx(ibond, l, m, n)
      case(7)
        oresult = Vxy_xx(ibond, l, m, n)
      case(8)
        oresult = Vxy_zz(ibond, l, m, n)
      end select
    case(5)
      select case(lm2)
      case(0)
        oresult = Vs_yz(ibond, l, m, n)
      case(1)
        oresult = Vx_yz(ibond, l, m, n)
      case(2)
        oresult = Vy_yz(ibond, l, m, n)
      case(3)
        oresult = Vz_yz(ibond, l, m, n)
      case(4)
        oresult = Vxy_yz(ibond, l, m, n)
      case(5)
        oresult = Vyz_yz(ibond, l, m, n)
      case(6)
        oresult = Vyz_zx(ibond, l, m, n)
      case(7)
        oresult = Vyz_xx(ibond, l, m, n)
      case(8)
        oresult = Vyz_zz(ibond, l, m, n)
      end select
    case(6)
      select case(lm2)
      case(0)
        oresult = Vs_zx(ibond, l, m, n)
      case(1)
        oresult = Vx_zx(ibond, l, m, n)
      case(2)
        oresult = Vy_zx(ibond, l, m, n)
      case(3)
        oresult = Vz_zx(ibond, l, m, n)
      case(4)
        oresult = Vxy_zx(ibond, l, m, n)
      case(5)
        oresult = Vyz_zx(ibond, l, m, n)
      case(6)
        oresult = Vzx_zx(ibond, l, m, n)
      case(7)
        oresult = Vzx_xx(ibond, l, m, n)
      case(8)
        oresult = Vzx_zz(ibond, l, m, n)
      end select
    case(7)
      select case(lm2)
      case(0)
        oresult = Vs_xx(ibond, l, m, n)
      case(1)
        oresult = Vx_xx(ibond, l, m, n)
      case(2)
        oresult = Vy_xx(ibond, l, m, n)
      case(3)
        oresult = Vz_xx(ibond, l, m, n)
      case(4)
        oresult = Vxy_xx(ibond, l, m, n)
      case(5)
        oresult = Vyz_xx(ibond, l, m, n)
      case(6)
        oresult = Vzx_xx(ibond, l, m, n)
      case(7)
        oresult = Vxx_xx(ibond, l, m, n)
      case(8)
        oresult = Vxx_zz(ibond, l, m, n)
      end select
    case(8)
      select case(lm2)
      case(0)
        oresult = Vs_zz(ibond, l, m, n)
      case(1)
        oresult = Vx_zz(ibond, l, m, n)
      case(2)
        oresult = Vy_zz(ibond, l, m, n)
      case(3)
        oresult = Vz_zz(ibond, l, m, n)
      case(4)
        oresult = Vxy_zz(ibond, l, m, n)
      case(5)
        oresult = Vyz_zz(ibond, l, m, n)
      case(6)
        oresult = Vzx_zz(ibond, l, m, n)
      case(7)
        oresult = Vxx_zz(ibond, l, m, n)
      case(8)
        oresult = Vzz_zz(ibond, l, m, n)
      end select
    end select 
  end function
end module

