module bspline
  use utils, only : wp

  implicit none
  !---------------
  ! Variables
  !---------------
  !  integer :: i,j

  !---------------
  ! Logic
  !---------------
  contains
    function linear_spline(y,x)
      !---------------
      ! Variables
      !---------------
      ! y are the three points on top of which
      ! the linear spline sits. Everywhere else
      ! it is zero.
      ! x is used for computing the desired value
      ! this spline assumes in this given point.
      real(wp), intent(in) :: y(3),x
      real(wp) :: h
      real(wp) :: linear_spline
      h = y(2) - y(1)
      linear_spline = 0

      if ( x < y(1) ) then
        linear_spline = 0
      elseif ( x <= y(2) ) then
        linear_spline = (x - y(1))/h
      elseif ( x <= y(3) ) then
        linear_spline =  (y(3) - x)/h
      end if
    end function linear_spline

    function linear_spline_deriv(y,x)
      !---------------
      ! Variables
      !---------------
      ! y are the three points on top of which
      ! the linear spline sits. Everywhere else
      ! it is zero.
      ! x is used for computing the desired value
      ! this spline assumes in this given point.
      real(wp), intent(in) :: y(3),x
      real(wp) :: h
      real(wp) :: linear_spline_deriv
      h = y(2) - y(1)
      linear_spline_deriv = 0.0_wp

      if ( x < y(1) ) then
        linear_spline_deriv = 0
      elseif ( x <= y(2) ) then
        linear_spline_deriv = 1.0_wp/h
      elseif ( x <= y(3) ) then
        linear_spline_deriv =  -1.0_wp/h
      end if
    end function linear_spline_deriv

    function cubic_spline(y, x)
      !---------------
      ! Variables
      !---------------
      ! y are the five points on top of which
      ! the cubi spline sits. Everywhere else
      ! it is zero.
      ! x is used for computing the desired value
      ! this spline assumes in this given point.
      real(wp), intent(in) :: y(5), x
      real(wp) :: h
      integer :: i,j
      real(wp) :: cubic_spline

      h = y(2) - y(1)
      if ( x < y(1) .or. x > y(5)) then
        cubic_spline = 0

      elseif ( x >= y(1) .and. x < y(2) ) then
        cubic_spline = ((x-y(1))**3)/(6*h**3)

      elseif ( x >= y(2) .and. x < y(3) ) then
        cubic_spline = &
        ( ((x-y(1))**2)*(y(3)-x) &
            + (x-y(1))*(y(4)-x)*(x-y(2))   &
        + ( (y(5)-x)*(x-y(2))**2 ) )/(6*h**3)

      elseif ( x >= y(3) .and. x < y(4) ) then
        cubic_spline = &
          ( (y(5)-x)*(x-y(2))*(y(4)-x) &
            + ((y(5)-x)**2)*(x-y(3))  &
        + (x-y(1))*(y(4)-x)**2 )/(6*h**3)

      elseif ( x >= y(4) .and. x < y(5) ) then
        cubic_spline = ((y(5)-x)**3)/(6*h**3)

      end if

    end function cubic_spline

    function cubic_spline_deriv(y, x)
      !---------------
      ! Variables
      !---------------
      ! y are the five points on top of which
      ! the cubi spline sits. Everywhere else
      ! it is zero.
      ! x is used for computing the desired value
      ! this spline assumes in this given point.
      real(wp), intent(in) :: y(5), x
      real(wp) :: h
      integer :: i,j
      real(wp) :: cubic_spline_deriv

      h = y(2) - y(1)
      if ( x < y(1) .or. x > y(5)) then
        cubic_spline_deriv = 0.0_wp

      elseif ( x >= y(1) .and. x < y(2) ) then
        cubic_spline_deriv = ((x-y(1))**2)/(2*h**3)

      elseif ( x >= y(2) .and. x < y(3) ) then
        cubic_spline_deriv = &
          - (9*x**2 - ( 2*y(5) + 2*y(4) + 2*y(3) + 6*y(2) + 6*y(1))*x &
          + 2*y(2)*y(5) &
          +( y(2) + y(1) )*y(4) + 2*y(1)*y(3) &
          + y(2)**2 + y(1)*y(2) + y(1)**2) / (6*h**3)

      elseif ( x >= y(3) .and. x < y(4) ) then
        cubic_spline_deriv = &
          (9*x**2 - ( 6*y(5) + 6*y(4) + 2*y(3) + 2*y(2) + 2*y(1) )*x &
          + y(5)**2 &
          + ( y(4)+2*y(3)+y(2) )*y(5) + y(4)**2 &
          + ( y(2)+2*y(1) )*y(4)) / (6*h**3)

      elseif ( x >= y(4) .and. x < y(5) ) then
        cubic_spline_deriv = (-(y(5)-x)**2)/(2*h**3)

      end if
    end function cubic_spline_deriv


    !--------------------
    ! New implementation
    !--------------------
    function cubic_spline_basis(x)
      real(wp) :: cubic_spline_basis
      real(wp), intent(in) :: x

      if ( x <= -2.0 ) then
        cubic_spline_basis = 0
      elseif ( -2.0 < x .and. x <= -1.0 ) then
        cubic_spline_basis = (1.0/4.0)*(2.0+x)**3
      elseif ( -1.0 < x .and. x <= 0.0 ) then
        cubic_spline_basis = (1.0/4.0)*((2.0+x)**3 - 4*(1.0+x)**3)
      elseif ( 0.0 < x .and. x <= 1.0 ) then
        cubic_spline_basis = (1.0/4.0)*((2.0-x)**3 - 4.0*(1-x)**3)
      elseif ( 1.0 < x .and. x <= 2.0 ) then
        cubic_spline_basis = (1.0/4.0)*(2.0-x)**3
      elseif( 2.0 < x ) then
        cubic_spline_basis = 0
      end if

    end function cubic_spline_basis

    function cubic_spline_basis_deriv(x)
      real(wp) :: cubic_spline_basis_deriv
      real(wp), intent(in) :: x

      if ( x <= -2.0 ) then
        cubic_spline_basis_deriv = 0.0
      elseif ( -2.0 < x .and. x <= -1.0 ) then
        cubic_spline_basis_deriv = (1.0/4.0)*3.0*(2.0+x)**2
      elseif ( -1.0 < x .and. x <= 0.0 ) then
        cubic_spline_basis_deriv = -(1.0/4.0)*((9*x**2)+(12*x))
      elseif ( 0.0 < x .and. x <= 1.0 ) then
        cubic_spline_basis_deriv = (1.0/4.0)*(3*x*(3*x-4.0))
      elseif ( 1.0 < x .and. x <= 2.0 ) then
        cubic_spline_basis_deriv = -(1.0/4.0)*3.0*(2.0-x)**2
      elseif( 2.0 < x ) then
        cubic_spline_basis_deriv = 0.0
      end if

    end function cubic_spline_basis_deriv

    function linear_spline_basis(x)
      real(wp) :: linear_spline_basis
      real(wp), intent(in) :: x

      if ( x <= -1.0 .or. x >= 1.0 ) then
        linear_spline_basis = 0.0
      elseif ( x < 0.0 ) then
        linear_spline_basis = (x+1.0)
      elseif ( x < 1.0 ) then
        linear_spline_basis = (1.0-x)
      end if

    end function linear_spline_basis

    function linear_spline_basis_deriv(x)
      real(wp) :: linear_spline_basis_deriv
      real(wp), intent(in) :: x

      if ( x <= -1.0 .or. x >= 1.0 ) then
        linear_spline_basis_deriv = 0.0
      elseif ( x < 0.0 ) then
        linear_spline_basis_deriv = x
      elseif ( x < 1.0 ) then
        linear_spline_basis_deriv = -x
      end if

    end function linear_spline_basis_deriv

end module bspline
