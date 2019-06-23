module user_fun
  use utils

  implicit none

contains
  function k(x)
    real(wp), intent(in) :: x
    real(wp) :: k
    k = 1
  end function k

  function q(x)
    real(wp), intent(in) :: x
    real(wp) :: q
    q = 0
  end function q

  function f(x)
    real(wp), intent(in) :: x
    real(wp) :: f
    f = (pi**2)*(sin(pi*x)-9*sin(3*pi*x))
    !f = -2
  end function f

  function u(x)
    real(wp), intent(in) :: x
    real(wp) :: u
    u = sin(pi*x) - sin(3*pi*x)
    !u = (x-1.0)*x
  end function u

end module user_fun
