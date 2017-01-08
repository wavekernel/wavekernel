!================================================================
! ELSES version 0.06
! Copyright (C) ELSES. 2007-2016 all rights reserved
!================================================================
module MSlaterKosterDerivative
  use MSlaterKosterTable

  implicit none

!  private
  public SlaterKosterDerivative

  type TVector3d
    real(8) :: x
    real(8) :: y
    real(8) :: z
  end type

  integer, parameter :: TB_S    = 0, &
                        TB_PX   = 1, &
                        TB_PY   = 2, &
                        TB_PZ   = 3, &
                        TB_DXY  = 4, &
                        TB_DYZ  = 5, &
                        TB_DZX  = 6, &
                        TB_DXX  = 7, &
                        TB_DZZ  = 8

contains

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine Swap(vector, x, y, z)
    type(TVector3d), intent(inout) :: vector
    integer, intent(in)    :: x, y, z

    type(TVector3d) :: the_buffer

    the_buffer = vector

    call SwapElement(the_buffer%x, x)
    call SwapElement(the_buffer%y, y)
    call SwapElement(the_buffer%z, z)
  contains
    subroutine SwapElement(element, axis)
      real(8),  intent(inout) :: element
      integer, intent(in)    :: axis 
 
      select case(axis)
      case(1)
        vector%x = element
      case(2)
        vector%y = element
      case(3)
        vector%z = element
      end select
    end subroutine 
  end subroutine

  type(TVector3d) function Dsx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult%x = (- l * l + 1) * bond%v(1)
    oresult%y = (- m * l)     * bond%v(1)
    oresult%z = (- n * l)     * bond%v(1)
  end function

  type(TVector3d) function Dsy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dsx(bond, m, l, n)
    call Swap(oresult, 2, 1, 3)
  end function

  type(TVector3d) function Dsz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dsx(bond, n, l, m)
    call Swap(oresult, 3, 1, 2)
  end function

  type(TVector3d) function Dxx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    real(8) :: the_Vpp

    the_Vpp = bond%v(1) - bond%v(2)

    oresult%x = (- 2 * l * l * l + 2 * l) * the_Vpp
    oresult%y = (- 2 * m * l * l)         * the_Vpp
    oresult%z = (- 2 * n * l * l)         * the_Vpp
  end function

  type(TVector3d) function Dxy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    real(8) :: the_Vpp

    the_Vpp = bond%v(1) - bond%v(2)

    oresult%x = (- 2 * l * l * m + m) * the_Vpp
    oresult%y = (- 2 * m * l * m + l) * the_Vpp
    oresult%z = (- 2 * n * l * m)     * the_Vpp
  end function

  type(TVector3d) function Dxz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxy(bond, l, n, m)
    call Swap(oresult, 1, 3, 2)
  end function

  type(TVector3d) function Dyy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxx(bond, m, l, n)
    call Swap(oresult, 2, 1, 3) 
  end function   

  type(TVector3d) function Dyz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxy(bond, m, n, l)
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Dzz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxx(bond, n, m, l)
    call Swap(oresult, 3, 2, 1)
  end function 

  type(TVector3d) function Ds_xy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult%x = sqrt(3d0) * m * (-2 * l**2 + 1) * bond%v(1)
    oresult%y = sqrt(3d0) * l * (-2 * m**2 + 1) * bond%v(1)
    oresult%z = sqrt(3d0) * (-2 * l * m * n)    * bond%v(1)
  end function

  type(TVector3d) function Ds_yz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Ds_xy(bond, m, n, l)
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Ds_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Ds_xy(bond, n, l, m)
    call Swap(oresult, 3, 1, 2)
  end function

  type(TVector3d) function Ds_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult%x = sqrt(3d0) * l * (-l**2 + m**2 + 1) * bond%v(1)
    oresult%y = sqrt(3d0) * m * (-l**2 + m**2 - 1) * bond%v(1)
    oresult%z = sqrt(3d0) * n * (-l**2 + m**2)     * bond%v(1)
  end function

  type(TVector3d) function Ds_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult%x = -3 * l * n**2        * bond%v(1)
    oresult%y = -3 * m * n**2        * bond%v(1)
    oresult%z =  3 * n * (-n**2 + 1) * bond%v(1)
  end function

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ConvertPDParameter(bond, v1, v2)
    type(TSKBond),    intent(in)  :: bond
    real(8), intent(out) :: v1, v2

    v1 = sqrt(3d0) * bond%v(1) - 2 * bond%v(2)
    v2 = bond%v(2)
  end subroutine

  type(TVector3d) function Dx_xy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertPDParameter(bond, v1, v2)

    oresult%x = l * m * ((-3 * l**2 + 2) * v1 -              v2)
    oresult%y = l**2 *   (-3 * m**2 + 1) * v1 + (1 - m**2) * v2
    oresult%z = m * n * ( -3 * l**2      * v1 -              v2)
  end function

  type(TVector3d) function Dx_yz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult =  Dz_xy(bond, m, n, l)
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Dx_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dy_xy(bond, n, l, m)
    call Swap(oresult, 3, 1, 2)
  end function

  type(TVector3d) function Dy_xy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertPDParameter(bond, v1, v2)

    oresult%x = m**2  *  (-3 * l**2 + 1) * v1 + (-l**2 + 1) * v2
    oresult%y = l * m * ((-3 * m**2 + 2) * v1 -               v2)
    oresult%z = l * n *  (-3 * m**2      * v1 -               v2)
  end function

  type(TVector3d) function Dy_yz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dx_xy(bond, m, n, l) 
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Dy_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dz_xy(bond, n, l, m)
    call Swap(oresult, 3, 1, 2)
  end function

  type(TVector3d) function Dz_xy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertPDParameter(bond, v1, v2)

    oresult%x = m * n * (-3 * l**2 + 1) * v1
    oresult%y = l * n * (-3 * m**2 + 1) * v1
    oresult%z = l * m * (-3 * n**2 + 1) * v1 
  end function

  type(TVector3d) function Dz_yz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dy_xy(bond, m, n, l)
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Dz_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dx_xy(bond, n, l, m)
    call Swap(oresult, 3, 1, 2)
  end function  
  
  type(TVector3d) function Dx_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertPDParameter(bond, v1, v2)

    oresult%x = 0.5d0 * ( 3 * l**2 * (-l**2 + m**2 + 1) - m**2) * v1 &
                          + (-l**2 + 1) * v2
    oresult%y = l * m * (0.5d0 *     (-3 * l**2 + 3 * m**2 - 2) * v1 &
                          -               v2)
    oresult%z = l * n * (1.5d0 *     (-l**2 + m**2)             * v1 &
                          -               v2)
  end function

  type(TVector3d) function Dy_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m ,n
    real(8)             :: v1, v2

    call ConvertPDParameter(bond, v1, v2)

    oresult%x = l * m * (0.5d0 *    (-3 * l**2 + 3 * m**2 + 2) * v1 &
                          +              v2)
    oresult%y = 0.5d0 * (3 * m**2 * (-l**2 + m**2 - 1) + l**2) * v1 &
                          + (m**2 - 1) * v2
    oresult%z = m * n * (1.5d0 *    (-l**2 + m**2)             * v1 &
                          +              v2)
  end function

  type(TVector3d) function Dz_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m ,n
    real(8)             :: v1, v2

    call ConvertPDParameter(bond, v1, v2)

    oresult%x = 0.5d0 * l * n * (-3 * l**2 + 3 * m**2 + 2) * v1
    oresult%y = 0.5d0 * m * n * (-3 * l**2 + 3 * m**2 - 2) * v1
    oresult%z = 0.5d0 * (l**2 - m**2) * (-3 * n**2 + 1)    * v1
  end function

  type(TVector3d) function Dx_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c

    call ConvertPDParameter(bond, v1, v2)

    c = 1d0 / (2d0 * sqrt(3d0))

    oresult%x = c * (((-3 * l**2 + 1) * (3 * n**2 -1) - 2 * l**2) * v1 &
                     + 2 * (l**2 - 1) * v2)
    oresult%y = c * l * m * ((-9 * n**2 + 1) * v1 &
                             + 2 * v2)
    oresult%z = c * l * n * ((-9 * n**2 + 7) * v1 &
                             + 2 * v2)
  end function

  type(TVector3d) function Dy_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c

    call ConvertPDParameter(bond, v1, v2)

    c = 1d0 / (2d0 * sqrt(3d0))

    oresult%x = c * l * m * ((-9 * n**2 + 1) * v1 &
                             + 2 * v2)
    oresult%y = c * (((-3 * m**2 + 1) * (3 * n**2 -1) - 2 * m**2) * v1 &
                     + 2 * (m**2 - 1) * v2)
    oresult%z = c * m * n * ((-9 * n**2 + 7) * v1 &
                             + 2 * v2)
  end function
   
  type(TVector3d) function Dz_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c

    call ConvertPDParameter(bond, v1, v2)

    c = 1d0 / (2d0 * sqrt(3d0))

    oresult%x = c * l * n * ((-9 * n**2 + 1) * v1 &
                             - 4 * v2)
    oresult%y = c * m * n * ((-9 * n**2 + 1) * v1 &
                             - 4 * v2)
    oresult%z = c * ((-(-3 * n**2 + 1)**2 + 4 * n**2) * v1 &
                     + 4 * (-n**2 + 1) * v2)
  end function

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine ConvertDDParameter(bond, v1, v2)
    type(TSKBond),     intent(in) :: bond
    real(8) , intent(out) :: v1, v2

    v1 = 3 * bond%v(1) - 4 * bond%v(2) + bond%v(3)
    v2 = bond%v(2) - bond%v(3)
  end subroutine

  type(TVector3d) function Dxy_xy(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertDDParameter(bond, v1, v2)

    oresult%x = 2 * l * (m**2 * (-2 * l**2 + 1) * v1 +  n**2      * v2)
    oresult%y = 2 * m * (l**2 * (-2 * m**2 + 1) * v1 +  n**2      * v2)
    oresult%z = 2 * n * (- 2 * l**2 * m**2      * v1 + (n**2 - 1) * v2)
  end function

  type(TVector3d) function Dxy_yz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertDDParameter(bond, v1, v2)

    oresult%x = n * (m**2 * (-4 * l**2 + 1)      * v1 + (-2 * l**2 + 1) * v2)
    oresult%y = 2 * l * m * n * ((-2 * m**2 + 1) * v1 - v2)
    oresult%z = l * (m**2 * (-4 * n**2 + 1)      * v1 + (-2 * n**2 + 1) * v2)
  end function

  type(TVector3d) function Dxy_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertDDParameter(bond, v1, v2)

    oresult%x = 2 * l * m * n * ((-2 * l**2 + 1) * v1 - v2)
    oresult%y = n * (l**2 * (-4 * m**2 + 1)      * v1 + (-2 * m**2 + 1) * v2)
    oresult%z = m * (l**2 * (-4 * n**2 + 1)      * v1 + (-2 * n**2 + 1) * v2)
  end function

  type(TVector3d) function Dyz_yz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxy_xy(bond, m, n, l)
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Dyz_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxy_yz(bond, m, n, l)
    call Swap(oresult, 2, 3, 1)
  end function

  type(TVector3d) function Dzx_zx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n

    oresult = Dxy_xy(bond, n, l, m)
    call Swap(oresult, 3, 1, 2)
  end function

  type(TVector3d) function Dxy_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertDDParameter(bond, v1, v2)

    oresult%x = 0.5d0 * m &
            * ((l**2 - m**2) * (-4 * l**2 + 1) + 2 * l**2) * v1
    oresult%y = 0.5d0 * l &
            * ((l**2 - m**2) * (-4 * m**2 + 1) - 2 * m**2) * v1
    oresult%z = 2 * l * m * n * (-l**2 + m**2) * v1
  end function

  type(TVector3d) function Dyz_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertDDParameter(bond, v1, v2)

    oresult%x = l * m * n * ((2 * (-l**2 + m**2) + 1) * v1 + 2 * v2)
    oresult%y = 0.5d0 * n  * ( &
                  ((l**2 - m**2) * (1 - 4 * m**2) - 2 * m**2) * v1   &
                + (4 * m**2 - 2)                              * v2)
    oresult%z = 0.5d0 * m &
            * ((l**2 - m**2) * (-4 * n**2 + 1)   * v1 + (4 * n**2 - 2) * v2)
  end function

  type(TVector3d) function Dzx_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2

    call ConvertDDParameter(bond, v1, v2)

    oresult%x = 0.5d0 * n                                      &
         * (((- 4 * l**2 + 1) * (l**2 - m**2) + 2 * l**2) * v1 &
            + (-4 * l**2 + 2) * v2)
    oresult%y = l * m * n * ((2 * (-l**2 + m**2) - 1) * v1 - 2 * v2)
    oresult%z = 0.5d0 * l &
            * ((l**2 - m**2) * (-4 * n**2 + 1)   * v1 + (-4 * n**2 + 2) * v2)
  end function

  type(TVector3d) function Dxy_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c 

    call ConvertDDParameter(bond, v1, v2)
    c = 1d0 / (2d0 * sqrt(3d0))

    oresult%x = c * m * (((3 * n**2 - 1) * (-4 * l**2 + 1) - 2 * l**2) * v1 &
                         + 4 * (2 * l**2 - 1) * v2)
    oresult%y = c * l * (((3 * n**2 - 1) * (-4 * m**2 + 1) - 2 * m**2) * v1 &
                         + 4 * (2 * m**2 - 1) * v2)
    oresult%z = c * 4 * l * m * n * ((-3 * n**2 + 2) * v1 + 2 * v2)
  end function

  type(TVector3d) function Dyz_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c 

    call ConvertDDParameter(bond, v1, v2)
    c = 1d0 / (2d0 * sqrt(3d0))

    oresult%x = c * 2 * l * m * n * ((-6 * n**2 + 1) * v1 - 2 * v2)
    oresult%y = c * n * (((-4 * m**2 + 1) * (3 * n**2 - 1) - 2 * m**2) * v1 & 
                         + 2 * (-2 * m**2 + 1) * v2)
    oresult%z = c * m * (((-4 * n**2 + 1) * (3 * n**2 - 1) + 4 * n**2) * v1 &
                         + 2 * (-2 * n**2 + 1) * v2)
  end function

  type(TVector3d) function Dzx_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c 

    call ConvertDDParameter(bond, v1, v2)
    c = 1d0 / (2d0 * sqrt(3d0))

    oresult%x = c * n * (((-4 * l**2 + 1) * (3 * n**2 - 1) - 2 * l**2) * v1 &
                         + 2 * (-2 * l**2 + 1) * v2)
    oresult%y = c * 2 * l * m * n * ((-6 * n**2 + 1) * v1 - 2 * v2)
    oresult%z = c * l * (((-4 * n**2 + 1) * (3 * n**2 - 1) + 4 * n**2) * v1 &
                         + 2 * (-2 * n**2 + 1) * v2)
  end function

  type(TVector3d) function Dxx_xx(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, t

    call ConvertDDParameter(bond, v1, v2)

    t = l**2 - m**2

    oresult%x = l * (t * (-t + 1) * v1 + 2d0 *  n**2      * v2)
    oresult%y = m * (t * (-t - 1) * v1 + 2d0 *  n**2      * v2)
    oresult%z = n * (-t**2        * v1 + 2d0 * (n**2 - 1) * v2)
  end function

  type(TVector3d) function Dxx_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c, t

    call ConvertDDParameter(bond, v1, v2)

    c = 1d0/sqrt(3d0)
    t = l**2 - m**2

    oresult%x = c * l * ((-t * (3 * n**2 - 1) + n**2 - l**2) * v1 &
                         + 2 * (t - 1) * v2)
    oresult%y = c * m * ((-t * (3 * n**2 - 1) - n**2 + m**2) * v1 &
                         + 2 * (t + 1) * v2)
    oresult%z = c * n * t * ((-3 * n**2 + 2) * v1 + 2 * v2)
  end function

  type(TVector3d) function Dzz_zz(bond, l, m, n) result (oresult)
    type(TSKBond),    intent(in) :: bond
    real(8), intent(in) :: l, m, n
    real(8)             :: v1, v2, c, t

    call ConvertDDParameter(bond, v1, v2)

    c = 1d0/sqrt(3d0)
    t = 3 * n**2 - 1

    oresult%x = l * n**2 * (-t * v1 - 2 * v2)
    oresult%y = m * n**2 * (-t * v1 - 2 * v2)
    oresult%z = n * (-n**2 + 1) * (t * v1 + 2 * v2)
  end function

  !! Copyright (C) ELSES. 2007-2016 all rights reserved
  subroutine SlaterKosterDerivative(ibond, lm1, lm2, l, m, n, deriv)
    type(TSKBond),     intent(in) :: ibond
    integer, intent(in) :: lm1, lm2
    real(8),  intent(in) :: l, m, n
    real(8), intent(out) :: deriv(:)
    type(TVector3d) :: oresult

    oresult = TVector3d(0d0, 0d0, 0d0)

    select case(lm1)
    case(TB_S)
      select case(lm2)
      case(TB_S)
        oresult = TVector3d(0d0, 0d0, 0d0)
      case(TB_PX)
        oresult = Dsx(ibond, l, m, n)
      case(TB_PY)
        oresult = Dsy(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dsz(ibond, l, m, n)
      case(TB_DXY)
        oresult = Ds_xy(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Ds_yz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Ds_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Ds_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Ds_zz(ibond, l, m, n)
      end select
    case(TB_PX)
      select case(lm2)
      case(TB_S)
        oresult = Dsx(ibond, l, m, n)
      case(TB_PX)
        oresult = Dxx(ibond, l, m, n)
      case(TB_PY)
        oresult = Dxy(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dxz(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dx_xy(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dx_yz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dx_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dx_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dx_zz(ibond, l, m, n)
      end select
    case(TB_PY)
      select case(lm2)
      case(TB_S)
        oresult = Dsy(ibond, l, m, n)
      case(TB_PX)
        oresult = Dxy(ibond, l, m, n)
      case(TB_PY)
        oresult = Dyy(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dyz(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dy_xy(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dy_yz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dy_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dy_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dy_zz(ibond, l, m, n)
      end select
    case(TB_PZ)
      select case(lm2)
      case(TB_S)
        oresult = Dsz(ibond, l, m, n)
      case(TB_PX)
        oresult = Dxz(ibond, l, m, n)
      case(TB_PY)
        oresult = Dyz(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dzz(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dz_xy(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dz_yz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dz_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dz_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dz_zz(ibond, l, m, n)
      end select
    case(TB_DXY)
      select case(lm2)
      case(TB_S)
        oresult = Ds_xy(ibond, l, m, n)
      case(TB_PX)
        oresult = Dx_xy(ibond, l, m, n)
      case(TB_PY)
        oresult = Dy_xy(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dz_xy(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dxy_xy(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dxy_yz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dxy_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dxy_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dxy_zz(ibond, l, m, n)
      end select
    case(TB_DYZ)
      select case(lm2)
      case(TB_S)
        oresult = Ds_yz(ibond, l, m, n)
      case(TB_PX)
        oresult = Dx_yz(ibond, l, m, n)
      case(TB_PY)
        oresult = Dy_yz(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dz_yz(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dxy_yz(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dyz_yz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dyz_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dyz_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dyz_zz(ibond, l, m, n)
      end select
    case(TB_DZX)
      select case(lm2)
      case(TB_S)
        oresult = Ds_zx(ibond, l, m, n)
      case(TB_PX)
        oresult = Dx_zx(ibond, l, m, n)
      case(TB_PY)
        oresult = Dy_zx(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dz_zx(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dxy_zx(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dyz_zx(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dzx_zx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dzx_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dzx_zz(ibond, l, m, n)
      end select
    case(TB_DXX)
      select case(lm2)
      case(TB_S)
        oresult = Ds_xx(ibond, l, m, n)
      case(TB_PX)
        oresult = Dx_xx(ibond, l, m, n)
      case(TB_PY)
        oresult = Dy_xx(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dz_xx(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dxy_xx(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dyz_xx(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dzx_xx(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dxx_xx(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dxx_zz(ibond, l, m, n)
      end select
    case(TB_DZZ)
      select case(lm2)
      case(TB_S)
        oresult = Ds_zz(ibond, l, m, n)
      case(TB_PX)
        oresult = Dx_zz(ibond, l, m, n)
      case(TB_PY)
        oresult = Dy_zz(ibond, l, m, n)
      case(TB_PZ)
        oresult = Dz_zz(ibond, l, m, n)
      case(TB_DXY)
        oresult = Dxy_zz(ibond, l, m, n)
      case(TB_DYZ)
        oresult = Dyz_zz(ibond, l, m, n)
      case(TB_DZX)
        oresult = Dzx_zz(ibond, l, m, n)
      case(TB_DXX)
        oresult = Dxx_zz(ibond, l, m, n)
      case(TB_DZZ)
        oresult = Dzz_zz(ibond, l, m, n)
      end select
    end select

    deriv(1) = oresult%x
    deriv(2) = oresult%y
    deriv(3) = oresult%z
  end subroutine
end module
