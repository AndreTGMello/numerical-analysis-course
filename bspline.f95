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

end module bspline
