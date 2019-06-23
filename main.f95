program main
use utils
use g_elimination
use num_integration
use bspline
use finite_elements

implicit none

!---------------
! Variables
!---------------

integer :: i,j,romberg_iterations
integer, allocatable :: n(:)
real(wp) :: x,error
real(wp), allocatable :: error_array_linear(:),error_array_cubic(:)
type(gpf):: gp

!---------------
! Logic
!---------------

write(*,*) "Program started!"

! Array of number of internal nodes for the integration
allocate(n(7))

! Array for plotting errors
allocate(error_array_linear(7))
allocate(error_array_cubic(7))

! Number of internal nodes for the integration
n(1:) = [7, 15, 31, 63, 127, 255, 511]

romberg_iterations = 12
! ---------------------------
! Solving with linear splines
! ---------------------------
do i = 1, 7
  call solve_finite_elements(n(i),linear_spline_basis,&
    linear_spline_basis,.false.,romberg_iterations,error)
  error_array_linear(i) = error
end do

call gp%plot(linspace(1.0_wp,7.0_wp,7),error_array_linear(1:7))

! ---------------------------
! Solving with cubic splines
! ---------------------------
do i = 1, 7
!  call solve_finite_elements(n(i),cubic_spline_basis,&
!    cubic_spline_basis_deriv,.true.,romberg_iterations,error)
!  error_array_cubic(i) = error
end do

!call gp%plot(linspace(1.0_wp,7.0_wp,7), error_array_cubic(1:7))

contains

end program main
