module num_integration

  implicit none

  contains

    function trapezoid(f,a,b,h)
      real, intent(in) :: a,b,h
      real, external :: f
      integer :: i,j,steps
      real :: trapezoid

      steps = (b-a)/h
      trapezoid = 0
      do i = 1, steps
        trapezoid = trapezoid + (f(a+(i-1)*h) + f(a+i*h))*h/2
      end do
    end function trapezoid

    function romberg(f,a,b,k)
      real, intent(in) :: a,b
      real, allocatable :: r(:)
      real, external :: f
      real :: h
      integer, intent(in) :: k
      integer :: i
      real :: romberg

      h = b - a
      allocate(r(k))

      do i = 1,k
        r(i) = trapezoid(f, a, b, h/2**(i-1))
      end do
      romberg = r(size(r))
    end function romberg

end module num_integration
