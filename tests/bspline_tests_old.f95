program bspline_tests
use utils
use bspline
use ogpf

implicit none
!---------------
! Variables
!---------------
integer :: i,j,total_tests, passed_tests
real(wp), allocatable :: y(:),x(:)
real(wp) :: result,true

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

y = linspace(1.0_wp,3.0_wp,3)
x = [y(1)-0.1,y(1),y(2),y(3),y(3)+0.1]
!call test_linear_spline(y,x,linear_spline,total_tests,passed_tests)

y = linspace(200.42_wp,799.99_wp,3)
x = [y(1)-0.1,y(1),y(2),y(3),y(3)+0.1]
!call test_linear_spline(y,x,linear_spline,total_tests,passed_tests)

write(*,*)
write(*,*) "!----------------------------------"
write(*,*) "! Testing linear spline derivatives"
write(*,*) "!----------------------------------"
write(*,*)

y = linspace(1.0_wp,3.0_wp,3)
x = [y(1)-0.1,y(1),y(2),y(3),y(3)+0.1]
!call test_linear_spline(y,x,linear_spline_deriv,total_tests,passed_tests)
!call plot_spline(y,linear_spline_deriv)

y = linspace(200.42_wp,799.99_wp,3)
x = [y(1)-0.1,y(1),y(2),y(3),y(3)+0.1]
!call test_linear_spline(y,x,linear_spline_deriv,total_tests,passed_tests)
!call plot_spline(y,linear_spline_deriv)

deallocate(y)
deallocate(x)

write(*,*)
write(*,*) "!----------------------"
write(*,*) "! Testing cubic splines"
write(*,*) "!----------------------"
write(*,*)

allocate(y(5))
allocate(x(7))

y = linspace(1.0_wp,5.0_wp,5)
x = [y(1)-0.1,y(1),y(2),y(3),y(4),y(5),y(5)+0.1]
!call test_cubic_spline(y,x,cubic_spline,total_tests,passed_tests)


y = linspace(234.432_wp,987.789_wp,5)
x = [y(1)-0.1,y(1),y(2),y(3),y(4),y(5),y(5)+0.1]
!call test_cubic_spline(y,x,cubic_spline,total_tests,passed_tests)

write(*,*)
write(*,*) "!--------------------------------"
write(*,*) "! Testing cubic spline derivative"
write(*,*) "!--------------------------------"
write(*,*)

y = linspace(1.0_wp,5.0_wp,5)
x = [y(1)-0.1,y(1),y(2),y(3),y(4),y(5),y(5)+0.1]
!call test_cubic_spline(y,x,cubic_spline_deriv,total_tests,passed_tests)
!call plot_spline(y,cubic_spline_deriv)

y = linspace(234.432_wp,987.789_wp,5)
x = [y(1)-0.1,y(1),y(2),y(3),y(4),y(5),y(5)+0.1]
!call test_cubic_spline(y,x,cubic_spline_deriv,total_tests,passed_tests)
!call plot_spline(y,cubic_spline_deriv)

y = linspace(0.0_wp,1.0_wp,5)
y(2) = 0.0_wp
y(3) = 0.0_wp
x = [y(1)-0.1,y(1),y(2),y(3),y(4),y(5),y(5)+0.1]
call test_cubic_spline(y,x,cubic_spline,total_tests,passed_tests)
!call plot_spline(y,cubic_spline_deriv)

!call plot_new_cubic_spline(cubic_spline_basis)

write(*,*) "!---------------------------"
write(*,'(a,g14.6)') " !Total tests ", total_tests-1
write(*,'(a,g14.6)') " !Tests passed", passed_tests
write(*,*) "!---------------------------"

contains
subroutine test_linear_spline(y,x,f,total_tests,passed_tests)
  real(wp), allocatable, intent(in) :: y(:),x(:)
  integer, intent(inout) :: passed_tests,total_tests
  real(wp) :: true,result
  real(wp), external :: f
  integer :: i

  write(*,*)
  write(*,*) "Test number ", total_tests
  write(*,*) "Grid is "
  write(*,*) y

  do i = 1, size(x)
    result = f(y,x(i))
    if ( x(i) <= y(1) .or. x(i) >= y(3) ) then
      true = 0.0
    else
      true = 1.0
    end if

    write(*,"(A,TR8,F15.5)") "Evaluating at", x(i)
    write(*,"(A,TR4,F15.5)") "Expected value is", true
    write(*,"(A,TR2,F15.5)") "Calculated value is", result
    if ( close(result,true,1.0_wp*1e-5) ) then
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

  call plot_spline(y, f)
end subroutine test_linear_spline

subroutine test_cubic_spline(y,x,f,total_tests,passed_tests)
  real(wp), intent(in) :: y(5),x(7)
  integer, intent(inout) :: passed_tests,total_tests
  real(wp) :: true,result
  real(wp), external :: f
  integer :: i

  write(*,*)
  write(*,"(A,I10)") "Test number ", total_tests
  write(*,*) "Grid is "
  write(*,*) y
  write(*,*)

  do i = 1, size(x)
    result = f(y,x(i))
    true = 0.0
    if ( close(x(i),y(2),1.0_wp*1e-1) .or. close(x(i),y(4),1.0_wp*1e-1) ) then
      true = 0.16666666_wp
    elseif ( close(x(i),y(3),1.0_wp*1e-1) ) then
      true = 0.66666666_wp
    end if

    write(*,"(A,TR8,F15.5)") "Evaluating at", x(i)
    write(*,"(A,TR4,F15.5)") "Expected value is", true
    write(*,"(A,TR2,F15.5)") "Calculated value is", result
    if ( close(result,true,1.0_wp*1e-5) ) then
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
  call plot_spline(y, f)
end subroutine test_cubic_spline

subroutine plot_new_cubic_spline(spline)
  real(wp), external :: spline
  integer, parameter :: i = 101
  integer, parameter :: n = 5
  integer, parameter :: m = n+2
  integer :: j
  real(wp) :: y(m)
  real(wp) :: x(i),fx(i),h
  type(gpf):: gp

  y = linspace(0.0_wp,1.0_wp,m)
  h = y(2) - y(1)

  x = linspace(y(1),y(3),i)
  do j = 1,i
    fx(j) = spline(x(j)/h) - 4*spline((x(j)+h)/h)
  end do
  call gp%title('B-Spline Phi_0')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  x = linspace(y(1),y(4),i)
  do j = 1,i
    fx(j) = spline((x(j)-y(2))/h) - spline((x(j)+h)/h)
  end do
  call gp%title('B-Spline Phi_1')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  x = linspace(y(1),y(5),i)
  do j = 1,i
    fx(j) = spline((x(j)-y(3))/h)
  end do
  call gp%title('B-Spline Phi_2')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  x = linspace(y(m-3),y(m),i)
  do j = 1,i
    fx(j) = spline((x(j)-y(n+1))/h) - spline((x(j)-(n+2)*h)/h)
  end do
  call gp%title('B-Spline Phi_6')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  x = linspace(y(m-2),y(m),i)
  do j = 1,i
    fx(j) = spline((x(j)-y(n+2))/h) - 4*spline((x(j)-(n+2)*h)/h)
  end do
  call gp%title('B-Spline Phi_7')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

end subroutine plot_new_cubic_spline

subroutine plot_spline(y,spline)
  real(wp), external :: spline
  integer, parameter :: i = 101
  integer :: j
  real(wp), intent(in) :: y(:)
  real(wp) :: x(i),fx(i)
  type(gpf):: gp
  x = linspace(y(1),y(size(y)),i)

  do j = 1,i
    fx(j) = spline(y,x(j))
  end do

  call gp%title('B-Spline')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')
end subroutine plot_spline

end program bspline_tests
