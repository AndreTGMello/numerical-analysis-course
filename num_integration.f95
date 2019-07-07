module num_integration
  use utils

  implicit none

  contains

    function trapezoid(f,a,b,h)
      real(wp), intent(in) :: a,b,h
      integer :: i,j,steps
      real(wp) :: trapezoid
      class(fun), intent(in) :: f

      steps = (b-a)/h
      trapezoid = 0
      do i = 1, steps
        trapezoid = trapezoid + (f%eval(a+(i-1)*h) + f%eval(a+i*h))*h/2
      end do
    end function trapezoid

    function romberg(f,a,b,k)
      real(wp), intent(in) :: a,b
      real(wp), allocatable :: r(:)
      class(fun), intent(in) :: f
      real(wp) :: h
      integer, intent(in) :: k
      integer :: i,j,m,np
      real(wp) :: romberg

      h = b - a
      allocate(r(k))

      r(1) = h*(f%eval(a)+f%eval(b))/2
      write(*,*) "h ", h
      write(*,*) "a ", a
      write(*,*) "b ", b
      write(*,*) "f(a) ", f%eval(a)
      write(*,*) "f(b) ", f%eval(b)
      write(*,*) "r(1) ", r(1)

      np = 1
      do i = 2,k
        h = h/2
        r(i) = r(i-1)/2
        do j = 1,np
          r(i) = r(i) + h*f%eval(a+(2*j-1)*h)
        end do
        m = 1
        do j = 2,i
          m = 4*m
          r(i-j+1) = (m*r(i-j+2)-r(i-j+1))/(m-1)
        end do
        np = np*2
      end do
      romberg = r(1)
      write(*,*) "r(end) ", r(1)
      write(*,*)
    end function romberg

end module num_integration
