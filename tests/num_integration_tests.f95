module test_function_object
  ! Object oriented approach for function
  ! evaluation.

  use utils, only : wp
  use num_integration

  implicit none

  abstract interface
    function f_interface(x)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: f_interface
    end function f_interface
  end interface

  ! function_object is an object that can point
  ! to any function that follow the f_interface.
  ! This allows for the object to change its behavior
  ! in runtime.
  type, extends(fun) :: function_object
    procedure(f_interface), pointer, nopass :: f
  contains
    procedure :: eval => specialized_function
  end type function_object

  contains

  function specialized_function(self, x)
    class(function_object), intent(in) :: self
    real(wp), intent(in) :: x
    real(wp) :: specialized_function
    specialized_function = self%f(x)
  end function specialized_function

end module test_function_object

program romgerg_tests
use utils
use num_integration
use test_function_object

! Note that not every test is meant to pass.
! Estimations by Romberg's using a low number
! of iterations should not give a decent
! approximation for higher order polynomials
! or complex functions.

implicit none

!---------------
! Variables
!---------------

integer :: i,j,l,m,n,total_tests,passed_tests
character (len=:), allocatable :: expression
real(wp) :: a(3),b(3)
integer :: k(3)
type(function_object) :: fobj

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
    fobj%f => poly0
    call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=x+1")
    fobj%f => poly1
    call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=3*x**2 + 2")
    fobj%f => poly2
    call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=33.33*x**5 + 8.8*x**4 + x**3 + 4.5*x*2 + 2.123")
    fobj%f => poly5
    call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=8*x**7/77 + (4*sin(2*x*pi) + 9*cos(11*x*pi))/5")
    fobj%f => ftrig
    call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    expression = trim("f(x)=e^x")
    fobj%f => fe
    call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)

    if ( a(j) > 0 ) then
      expression = trim("f(x)=ln(x)")
      fobj%f => fln
      call test_integration(fobj,a(j),b(j),k(i),expression,passed_tests,total_tests,romberg,close)
    end if
  end do
  write(*,*)
end do

write(*,*) "!---------------------------"
write(*,'(a,g14.6)') " !Total tests ", total_tests-1
write(*,'(a,g14.6)') " !Tests passed", passed_tests
write(*,*) "!---------------------------"


contains

  function poly0(x)
    real(wp), intent(in) :: x
    real(wp) :: poly0

    poly0 = 1.0_wp
  end function poly0

  function poly1(x)
    real(wp), intent(in) :: x
    real(wp) :: poly1

    poly1 = x + 1
  end function poly1

  function poly2(x)
    real(wp), intent(in) :: x
    real(wp) :: poly2

    poly2 = 3*x**2 + 2
  end function poly2

  function poly5(x)
    real(wp), intent(in) :: x
    real(wp) :: poly5

    poly5 = 33.33*x**5 + 8.8*x**4 + x**3 + 4.5*x*2 + 2.123
  end function poly5

  function ftrig(x)
    real(wp), intent(in) :: x
    real(wp) :: ftrig

    ftrig = 8*x**7/77 + (4*sin(2*x*pi) + 9*cos(11*x*pi))/5
  end function ftrig

  function fe(x)
    real(wp), intent(in) :: x
    real(wp) :: fe

    fe = exp(x)
  end function fe

  function fln(x)
    real(wp), intent(in) :: x
    real(wp) :: fln

    fln = log(x)
  end function fln

  subroutine test_integration(fobj,a,b,k,expression,passed_tests,total_tests,romberg,close)
    real(wp), intent(in) :: a,b
    real(wp) :: true,result,abserr
    real(wp), parameter :: epsabs = 0.0D+00
    real(wp), parameter :: epsrel = 0.001D+00
    integer, intent(inout) :: passed_tests, total_tests
    integer, intent(in) :: k
    integer, parameter :: key = 6
    integer, parameter :: limit = 500
    integer, parameter :: lenw = 4 * limit
    integer :: ier
    integer :: iwork(limit)
    integer :: last
    integer :: neval
    real(wp) :: work(lenw)

!    real(wp), external :: f
!    real(wp), external :: fcomp
    class(function_object), intent(in) :: fobj
    real(wp), external :: romberg
    logical, external :: close

    character (len=:), allocatable, intent(in) :: expression

    write(*,*)
    write(*,*) "Test number ", total_tests
    total_tests = total_tests + 1
    result = romberg(fobj, a, b, k)
    call dqag ( fobj%f , a, b, epsabs, epsrel, key, true, abserr, neval, ier,&
    limit, lenw, last, iwork, work )

    write (*, '(A,A)') '  Integrand is ', expression
    write ( *, '(A,F15.5)' ) '  Integral left endpoint A     ', a
    write ( *, '(A,F15.5)' ) '  Integral right endpoint B    ', b
    write (*, '(A,I15)') '  #Romberg iterations(k) is', k
    write ( *, '(A,F15.5)' ) '  Integral is                  ', true
    write ( *, '(A,F15.5)' ) '  Estimated integral is        ', result
    if ( close(result,true,1.0_wp*1e-5) ) then
      passed_tests = passed_tests + 1
      write(*,*) "!Test passed!"
    else
      write(*,*) "!-----------!"
      write(*,*) "!Test failed!"
      write(*,*) "!-----------!"
    end if
    write(*,*)

  end subroutine test_integration
end program romgerg_tests
