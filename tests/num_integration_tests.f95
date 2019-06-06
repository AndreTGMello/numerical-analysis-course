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

  function poly5(x)
    real, intent(in) :: x
    real :: poly5

    poly5 = 33.33*x**5 + 8.8*x**4 + x**3 + 4.5*x*2 + 2.123
  end function poly5

  function ftrig(x)
    real, intent(in) :: x
    real :: ftrig

    ftrig = 8*x**7/77 + (4*sin(2*x*pi) + 9*cos(11*x*pi))/5
  end function ftrig

  function fe(x)
    real, intent(in) :: x
    real :: fe

    fe = exp(x)
  end function fe

  function fln(x)
    real, intent(in) :: x
    real :: fln

    fln = log(x)
  end function fln

  subroutine test_integration(f,a,b,k,expression,passed_tests,total_tests,romberg,close)
    real, intent(in) :: a,b
    real :: true,result
    real ::  abserr
    real, parameter :: epsabs = 0.0E+00
    real, parameter :: epsrel = 0.001E+00
    integer, parameter :: key = 6
    integer :: ier
    integer, intent(inout) :: passed_tests, total_tests
    integer :: neval
    integer, intent(in) :: k
    real, external :: f
    real, external :: romberg
    logical, external :: close
    character (len=:), allocatable, intent(in) :: expression

    write(*,*)
    write(*,*) "Test number ", total_tests
    total_tests = total_tests + 1
    result = romberg(f, a, b, k)
    call qag ( f, a, b, epsabs, epsrel, key, true, abserr, neval, ier )
    write (*, '(a,a)') '  Integrand is ', expression
    write ( *, '(a,g14.6)' ) '  Integral left endpoint A     ', a
    write ( *, '(a,g14.6)' ) '  Integral right endpoint B    ', b
    write (*, '(a,g14.6)') '  #Romberg iterations(k) is', k
    write ( *, '(a,g14.6)' ) '  Integral is                  ', true
    write ( *, '(a,g14.6)' ) '  Estimated integral is        ', result
    if ( close(result,true,1e-1) ) then
      passed_tests = passed_tests + 1
      write(*,*) "Test passed"
    else
      write(*,*) "!-----------!"
      write(*,*) "!Test failed!"
      write(*,*) "!-----------!"
    end if
    write(*,*)

  end subroutine test_integration

end module test_functions


program romgerg_tests
use utils
use num_integration
use test_functions

implicit none

!---------------
! Variables
!---------------

integer :: i,j,l,m,n,total_tests,passed_tests
character (len=:), allocatable :: expression
real :: a(3),b(3)
integer :: k(3)

!---------------
! Logic
!---------------

total_tests = 1
passed_tests = 0
k = [2, 5, 10]
a = [-2.0, 0.0, 1.0]
b = [1.0, 1.0, 5.0]


write(*,*)
write(*,*) "!-------------------------------"
write(*,*) "! Testing integration by Romberg"
write(*,*) "!-------------------------------"
write(*,*)

do i = 1, size(k)
  write(*,*)
  write(*,*) "!----------------------------------"
  write(*,*) "! With k = ",k(i)
  write(*,*) "!----------------------------------"
  do j = 1, size(a)
    expression = trim("f(x)=1")
    call test_integration(poly0,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=x+1")
    call test_integration(poly1,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=3*x**2 + 2")
    call test_integration(poly2,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=33.33*x**5 + 8.8*x**4 + x**3 + 4.5*x*2 + 2.123")
    call test_integration(poly5,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=8*x**7/77 + (4*sin(2*x*pi) + 9*cos(11*x*pi))/5")
    call test_integration(ftrig,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=e^x")
    call test_integration(fe,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    if ( a(j) > 0 ) then
      expression = trim("f(x)=ln(x)")
      call test_integration(fln,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)
    end if
  end do
  write(*,*)
end do

write(*,*) "!---------------------------"
write(*,'(a,g14.6)') " !Total tests ", total_tests-1
write(*,'(a,g14.6)') " !Tests passed", passed_tests
write(*,*) "!---------------------------"
end program romgerg_tests
