module test_functions

implicit none

real, parameter :: pi = 4.D0*datan(1.D0)

contains

  function poly0(x)
    real, intent(in) :: x
    real :: poly0

    poly0 = 1
  end function poly0

  function poly1(x)
    real, intent(in) :: x
    real :: poly1

    poly1 = x + 1
  end function poly1

  function poly2(x)
    real, intent(in) :: x
    real :: poly2

    poly2 = 3*x**2 + 2
  end function poly2

  function poly10(x)
    real, intent(in) :: x
    real :: poly10

    poly10 = 3*x**2 + 2
  end function poly10

  function ftrig(x)
    real, intent(in) :: x
    real :: ftrig

    ftrig = 8*x**7/77 + (4*sin(2*x*pi) + 9*cos(11*x*pi))/5
  end function ftrig

end module test_functions


program romgerg_tests
use utils
use num_integration
use test_functions

implicit none

!---------------
! Variables
!---------------

integer :: i,j,k,l,m,n,total_tests,passed_tests
real :: apoly1_0_1, apoly2_0_1, apoly10_0_1, aftrig_0_1
real :: apoly1_8_20, apoly2_8_20, apoly10_8_20, aftrig_8_20

real :: a,b,true,result
real ( kind = 4 ) abserr
real ( kind = 4 ), parameter :: epsabs = 0.0E+00
real ( kind = 4 ), parameter :: epsrel = 0.001E+00
integer ( kind = 4 ) ier
integer ( kind = 4 ), parameter :: key = 6
integer ( kind = 4 ) neval

!---------------
! Logic
!---------------

total_tests = 1
passed_tests = 0

write(*,*)
write(*,*) "!--------------------------------------"
write(*,*) "! Testing integration by Trapezoid Rule"
write(*,*) "!--------------------------------------"
write(*,*)

write(*,*)
write(*,*) "Test number ", total_tests
a = -1.0
b = 1.0
result = trapezoid(poly0, a, b, (b-a)/1e5)
true = 2.0
write (*, '(a)') '  Integrand is f(x)=1'
write ( *, '(a,g14.6)' ) '  Integral left endpoint A     ', a
write ( *, '(a,g14.6)' ) '  Integral right endpoint B    ', b
write ( *, '(a,g14.6)' ) '  Integral is                  ', true
write ( *, '(a,g14.6)' ) '  Estimated integral is        ', result
if ( close(result,true,1e-1) ) then
  passed_tests = passed_tests + 1
  write(*,*) "Test passed"
end if
write(*,*)

write(*,*)
total_tests = total_tests + 1
write(*,*) "Test number ", total_tests
a = -1.0
b = 1.0
result = trapezoid(poly1, a, b, (b-a)/1e5)
call qag ( poly1, a, b, epsabs, epsrel, key, true, abserr, neval, ier )
write ( *, '(a)' ) '  Integrand is f(x)=x+1'
write ( *, '(a,g14.6)' ) '  Integral left endpoint A     ', a
write ( *, '(a,g14.6)' ) '  Integral right endpoint B    ', b
write ( *, '(a,g14.6)' ) '  Integral is                  ', true
write ( *, '(a,g14.6)' ) '  Estimated integral is        ', result
if ( close(result,true,1e-1) ) then
  passed_tests = passed_tests + 1
  write(*,*) "Test passed"
end if
write(*,*)


write(*,*) "!--------------------------"
write(*,'(a,g14.6)') "  Total tests", total_tests
write(*,'(a,g14.6)') "  Tests passed", passed_tests
write(*,*) "!--------------------------"
end program romgerg_tests
