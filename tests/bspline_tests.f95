program bspline_tests
use utils
use bspline
use finite_elements
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

!call plot_linear_spline()

write(*,*)
write(*,*) "!----------------------------------"
write(*,*) "! Testing linear spline derivatives"
write(*,*) "!----------------------------------"
write(*,*)

!call plot_linear_spline_deriv()

write(*,*)
write(*,*) "!----------------------"
write(*,*) "! Testing cubic splines"
write(*,*) "!----------------------"
write(*,*)

!call plot_cubic_spline(cubic_spline_basis)

write(*,*)
write(*,*) "!---------------------------------"
write(*,*) "! Testing cubic spline derivatives"
write(*,*) "!---------------------------------"
write(*,*)

call plot_cubic_spline(cubic_spline_basis_deriv)

contains

subroutine plot_cubic_spline(s)
  real(wp), external :: s
  integer, parameter :: i = 101
  integer, parameter :: n = 7
  integer, parameter :: m = n+2
  integer :: j,k
  real(wp) :: y(m)
  real(wp) :: x(i),fx(i),h
  type(gpf):: gp
  type(inner_product_obj) :: ipo

  y = linspace(0.0_wp,1.0_wp,m)
  h = y(2) - y(1)
  write(*,*) "Node grid:"
  write(*,'(F10.5)', advance="no") y
  write(*,*)

  k = 0
  call set_inner_product(ipo,y,k,n,constant_fun,s,.true.)
  x = linspace(y(1),y(3),i)
  do j = 1,i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Cubic B-Spline Phi_0')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 1
  call set_inner_product(ipo,y,k,n,constant_fun,s,.true.)
  x = linspace(y(1),y(4),i)
  do j = 1,i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Cubic B-Spline Phi_1')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 2
  call set_inner_product(ipo,y,k,n,constant_fun,s,.true.)
  x = linspace(y(1),y(5),i)
  do j = 1,i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Cubic B-Spline Phi_2')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 7
  call set_inner_product(ipo,y,k,n,constant_fun,s,.true.)
  x = linspace(y(m-3),y(m),i)
  do j = 1,i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Cubic B-Spline Phi_7')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 8
  call set_inner_product(ipo,y,k,n,constant_fun,s,.true.)
  x = linspace(y(m-2),y(m),i)
  do j = 1,i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Cubic B-Spline Phi_8')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

end subroutine plot_cubic_spline

subroutine plot_linear_spline()
  integer, parameter :: i = 101
  integer, parameter :: n = 7
  integer, parameter :: m = n+2
  integer :: j,k
  real(wp) :: x(i),fx(i),h
  type(gpf):: gp
  type(inner_product_obj) :: ipo

  y = linspace(0.0_wp,1.0_wp,m)
  h = y(2) - y(1)

  ! --------------
  ! B-Spline Tests
  ! --------------
  ! ---------
  ! Test One
  ! ---------
  write(*,*) "Node grid:"
  write(*,'(F10.5)', advance="no") y
  write(*,*)

  k = 0
  x = linspace(y(1),y(3),i)
  call set_inner_product(ipo,y,k,n,&
                        constant_fun,linear_spline_basis,.false.)
  do j = 1, i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Linear B-Spline Phi_0')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 1
  x = linspace(y(2),y(4),i)
  call set_inner_product(ipo,y,k,n,&
                        constant_fun,linear_spline_basis,.false.)
  do j = 1, i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Linear B-Spline Phi_1')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 6
  x = linspace(y(m-2),y(m),i)
  call set_inner_product(ipo,y,k,n,&
                        constant_fun,linear_spline_basis,.false.)
  do j = 1, i
    fx(j) = ipo%eval(x(j))
  end do
  call gp%title('Linear B-Spline Phi_7')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

end subroutine plot_linear_spline


subroutine plot_linear_spline_deriv()
  integer, parameter :: i = 101
  integer, parameter :: n = 7
  integer, parameter :: m = n+2
  integer :: j,k
  real(wp) :: x(i),fx(i),h
  type(gpf):: gp

  y = linspace(0.0_wp,1.0_wp,m)
  h = y(2) - y(1)

  write(*,*) "Node grid:"
  write(*,'(F10.5)', advance="no") y
  write(*,*)

  k = 0
  x = linspace(y(1),y(3),i)
  do j = 1, i
    if ( x(j) == 0 ) then
      fx(j) = 1.0/h
    elseif ( x(j) == size(y) ) then
      fx(j) = -1.0/h
    elseif ( x(j) <= y( k+2 ) ) then
      fx(j) = 1.0/h
    elseif ( x(j) <= y( k+3 ) ) then
      fx(j) = -1.0/h
    else
      fx(j) = 0.0/h
      !write(*,*) "Error ", x
    end if
  end do
  call gp%title('Linear B-Spline Phi_0 Deriv')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 1
  x = linspace(y(2),y(4),i)
  do j = 1, i
    if ( x(j) == 0 ) then
      fx(j) = 1.0/h
    elseif ( x(j) == size(y) ) then
      fx(j) = -1.0/h
    elseif ( x(j) <= y( k+2 ) ) then
      fx(j) = 1.0/h
    elseif ( x(j) <= y( k+3 ) ) then
      fx(j) = -1.0/h
    else
      fx(j) = 0.0/h
      !write(*,*) "Error ", x
    end if
  end do
  call gp%title('Linear B-Spline Phi_1 Deriv')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

  k = 6
  x = linspace(y(m-2),y(m),i)
  do j = 1, i
    if ( x(j) == 0 ) then
      fx(j) = 1.0/h
    elseif ( x(j) == size(y) ) then
      fx(j) = -1.0/h
    elseif ( x(j) <= y( k+2 ) ) then
      fx(j) = 1.0/h
    elseif ( x(j) <= y( k+3 ) ) then
      fx(j) = -1.0/h
    else
      fx(j) = 0.0/h
      !write(*,*) "Error ", x
    end if
  end do
  call gp%title('Linear B-Spline Phi_7 Deriv')
  call gp%options('set key top right; set grid')
  call gp%plot(x,fx,'title "B(x)" with lines lt 1 lw 1')

end subroutine plot_linear_spline_deriv

function constant_fun(x)
  real(wp) :: constant_fun
  real(wp), intent(in) :: x
  constant_fun = 1.0_wp
end function constant_fun

end program bspline_tests
