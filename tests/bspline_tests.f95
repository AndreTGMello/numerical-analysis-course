program bspline_tests
use utils
use bspline
use ogpf

implicit none
!---------------
! Variables
!---------------
integer :: i,j,total_tests, passed_tests
real, allocatable :: y(:),x(:)
real :: result,true

!---------------
! Logic
!---------------
total_tests = 1
passed_tests = 0

write(*,*)
write(*,*) "!-----------------------"
write(*,*) "! Testing linear splines"
write(*,*) "!-----------------------"
write(*,*)

allocate(y(3))
allocate(x(4))

y = [1,2,3]
x = [0,1,2,3,4]
call test_linear_spline(y,x,total_tests,passed_tests)

y = [0.5,2.5,4.5]
x = [0.4,0.5,2.5,4.5,4.6]
call test_linear_spline(y,x,total_tests,passed_tests)

deallocate(y)
deallocate(x)

write(*,*)
write(*,*) "!----------------------"
write(*,*) "! Testing cubic splines"
write(*,*) "!----------------------"
write(*,*)

allocate(y(5))
allocate(x(7))

y = [1,2,3,4,5]
x = [0,1,2,3,4,5,6]
call test_cubic_spline(y,x,total_tests,passed_tests)


y = [10.55,40.55,70.55,100.55,130.55]
x = [10.00,10.55,40.55,70.55,100.55,130.55,131.00]
call test_cubic_spline(y,x,total_tests,passed_tests)

write(*,*) "!---------------------------"
write(*,'(a,g14.6)') " !Total tests ", total_tests-1
write(*,'(a,g14.6)') " !Tests passed", passed_tests
write(*,*) "!---------------------------"

contains
subroutine test_linear_spline(y,x,total_tests,passed_tests)
  real, allocatable, intent(in) :: y(:),x(:)
  integer, intent(inout) :: passed_tests,total_tests
  real :: true,result
  integer :: i

  write(*,*)
  write(*,*) "Test number ", total_tests
  write(*,*) "Grid is "
  write(*,*) y

  do i = 1, size(x)
    result = linear_spline(y,x(i))
    if ( x(i) <= y(1) .or. x(i) >= y(3) ) then
      true = 0.0
    else
      true = 1.0
    end if

    write(*,*) "Linear spline"
    write(*,*) "Evaluating at ", x(i)
    write(*,*) "Expected value is ", true
    write(*,*) "Calculated value is ", result
    if ( close(result,true,1e-5) ) then
      passed_tests = passed_tests + 1
      write(*,*) "Test passed"
    else
      write(*,*) "!-----------!"
      write(*,*) "!Test failed!"
      write(*,*) "!-----------!"
    end if
    total_tests = total_tests + 1
    write(*,*)
  end do

  call plot_spline(y, linear_spline)
end subroutine test_linear_spline

subroutine test_cubic_spline(y,x,total_tests,passed_tests)
  real, allocatable, intent(in) :: y(:),x(:)
  integer, intent(inout) :: passed_tests,total_tests
  real :: true,result
  integer :: i

  write(*,*)
  write(*,*) "Test number ", total_tests
  write(*,*) "Grid is "
  write(*,*) y

  do i = 1, size(x)
    true = 10.0
    result = cubic_spline(y,x(i))
    if ( x(i) <= y(1) .or. x(i) >= y(5) ) then
      true = 0.0
    elseif ( x(i) == y(2) .or. x(i) == y(4) ) then
      true = 1/6
    elseif ( x(i) == y(3) ) then
      true = 2/3
    end if

    write(*,*) "Cubic spline"
    write(*,*) "Evaluating at ", x(i)
    write(*,*) "Expected value is ", true
    write(*,*) "Calculated value is ", result
    if ( close(result,true,1e-5) ) then
      passed_tests = passed_tests + 1
      write(*,*) "Test passed"
    else
      write(*,*) "!-----------!"
      write(*,*) "!Test failed!"
      write(*,*) "!-----------!"
    end if
    total_tests = total_tests + 1
    write(*,*)
  end do
  call plot_spline(y, cubic_spline)
end subroutine test_cubic_spline

subroutine plot_spline(y,spline)
  real, external :: spline
  integer, parameter :: i = 101
  integer :: j
  real, intent(in) :: y(:)
  real :: x(i), fx(i)
  real(wp) :: a, b, dble_x(i), dble_fx(i)
  type(gpf):: gp
  a = dble(y(1))
  b = dble(y(size(y)))
  x = linspace(a,b,i)
  dble_x = dble(x)

  do j = 1,i
    fx(j) = spline(y,x(j))
  end do
  dble_fx = dble(fx)

  call gp%title('B-Spline')
  call gp%options('set key top right; set grid')
  call gp%plot(dble_x,dble_fx,'title "B(x)" with lines lt 1 lw 1')
end subroutine plot_spline

end program bspline_tests
