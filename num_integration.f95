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
      integer :: i,j,m,np
      real :: romberg

      h = b - a
      allocate(r(k))

      r(1) = h*(f(a)+f(b))/2
      np = 1
      do i = 2,k
        h = h/2
        r(i) = r(i-1)/2
        do j = 1,np
          r(i) = r(i) + h*f(a+(2*j-1)*h)
        end do
        m = 1
        do j = 2,i
          m = 4*m
          r(i-j+1) = (m*r(i-j+2)-r(i-j+1))/(m-1)
        end do
        np = np*2
      end do
      romberg = r(1)
    end function romberg

end module num_integration
