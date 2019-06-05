module bspline

  implicit none
  !---------------
  ! Variables
  !---------------
!  integer :: i,j

  !---------------
  ! Logic
  !---------------
  contains
    function cubic_spline(y, x)
      !---------------
      ! Variables
      !---------------
      ! y is the partition
      ! x are the function values at y
      ! b is the values for the spline
      ! y_x is the extended partition
      ! h is the step size, fixed in this case
      real, allocatable, intent(in) :: y(:), x(:)
      real, allocatable, intent(out) :: b(:)
      real, allocatable :: y_x(:)
      real :: h
      integer :: i,j

      h = y(2) - y(1)

      ! Build extended partition y_x (y_extended)
      allocate(y_x(size(y)+2*3))
      y_x(1) = y(1) - 3*h
      do i = 2, 2*3+size(y)
        y_x(i) = y_x(i) + y_x(1)*i*h
      end do

      allocate(b(size(y_x)))
      do j = 1, 4
        do i = 1, size(b)

        end do
      end do

    end function cubic_spline


end module bspline
