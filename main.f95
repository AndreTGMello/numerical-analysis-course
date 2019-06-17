module user_defined_functions
  use num_integration
  use utils, only : wp
  implicit none

  abstract interface
    function f_interface(x)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: f_interface
    end function f_interface
    function linear_spline_interface(y,x)
      import :: wp
      real(wp), intent(in) :: y(3),x
      real(wp) :: linear_spline_interface
    end function linear_spline_interface
    function cubic_spline_interface(y,x)
      import :: wp
      real(wp), intent(in) :: y(5),x
      real(wp) :: cubic_spline_interface
    end function cubic_spline_interface
  end interface

  type, extends(fun) :: inner_product_l_obj
    ! Support nodes for base i and base j of the spline
    real(wp), allocatable :: yi(:), yj(:)

    ! Pointers to f, q and k function, as described by
    ! the formula.
    procedure(f_interface), pointer, nopass :: f,q,k

    ! Pointers to linear_spline functions
    procedure(linear_spline_interface), pointer, nopass :: l,dl

    ! Pointers to cubic_spline functions
    procedure(cubic_spline_interface), pointer, nopass :: c,dc

    ! Only cubic or linear splines pointers will be
    ! assigned, depending on the size of the y array.

  contains
    procedure :: eval => inner_product_l_eval
  end type inner_product_l_obj

  type, extends(fun) :: inner_product_obj
    ! Support nodes for base i and base j of the spline
    real(wp), allocatable :: y(:)

    ! Pointers to f, q and k function, as described by
    ! the formula.
    procedure(f_interface), pointer, nopass :: f

    ! Pointers to linear_spline functions
    procedure(linear_spline_interface), pointer, nopass :: l

    ! Pointers to cubic_spline functions
    procedure(cubic_spline_interface), pointer, nopass :: c

    ! Only cubic or linear splines pointers will be
    ! assigned, depending on the size of the y array.

  contains
    procedure :: eval => inner_product_eval
  end type inner_product_obj

contains
  function inner_product_l_eval(self, x)
    class(inner_product_l_obj), intent(in) :: self
    real(wp), intent(in) :: x
    real(wp) :: inner_product_l_eval

    if ( size(self%yi) == 3 ) then
      ! Inner product for linear_splines
      inner_product_l_eval = &
      self%k(x) * self%dl(self%yi,x) * self%dl(self%yj,x) &
      + self%q(x) * self%l(self%yi,x) * self%l(self%yj,x)

    elseif ( size(self%yi) == 5 ) then
      ! Inner product for cubic_splines
      inner_product_l_eval =  &
      self%k(x) * self%dc(self%yi,x) * self%dc(self%yj,x) &
      + self%q(x) * self%c(self%yi,x) * self%c(self%yj,x)

    else
      ! Something went wrong
      write(*,*) "Incorrect usage."
      inner_product_l_eval = 0.0_wp
    end if
  end function inner_product_l_eval

  function inner_product_eval(self, x)
    class(inner_product_obj), intent(in) :: self
    real(wp), intent(in) :: x
    real(wp) :: inner_product_eval

    if ( size(self%y) == 3 ) then
      ! Inner product for linear_splines
      inner_product_eval = &
      self%f(x) * self%l(self%y,x)

    elseif ( size(self%y) == 5 ) then
      ! Inner product for cubic_splines
      inner_product_eval =  &
      self%f(x) * self%c(self%y,x)
    else
      ! Something went wrong
      write(*,*) "Incorrect usage."
      inner_product_eval = 0.0_wp
    end if
  end function inner_product_eval

  subroutine set_inner_product_l(ipo,yi,yj,q,k,f,phi,dphi)
    class(inner_product_l_obj), intent(inout) :: ipo
    real(wp), external :: q,k,f,phi,dphi
    real(wp), intent(in) :: yi(:),yj(:)
!    allocate(ipo%yi(size(yi)))
!    allocate(ipo%yj(size(yj)))
    ipo%yi = yi
    ipo%yj = yj
    ipo%q => q
    ipo%k => k
    ipo%f => f
    if ( size(yi) == 3 ) then
      ! Pointers to linear_spline functions
      ipo%l => phi
      ipo%dl => dphi
    elseif ( size(yi) == 5 ) then
      ! Pointers to cubic_spline functions
      ipo%c => phi
      ipo%dc => dphi
    end if
  end subroutine set_inner_product_l

  subroutine set_inner_product(ipo,y,f,phi)
    class(inner_product_obj), intent(inout) :: ipo
    real(wp), external :: f,phi
    real(wp), intent(in) :: y(:)
    ipo%y = y
    ipo%f => f
    if ( size(y) == 3 ) then
      ! Pointers to linear_spline functions
      ipo%l => phi
    elseif ( size(y) == 5 ) then
      ! Pointers to cubic_spline functions
      ipo%c => phi
    end if
  end subroutine set_inner_product

end module user_defined_functions

program main
use utils
!use bm_ge
use g_elimination
use num_integration
use bspline
use user_defined_functions
use ogpf, only : linspace

implicit none


!---------------
! Variables
!---------------

integer :: i,j,l,o,m,ub,lb,ms
integer :: num_nodes
integer :: num_support_nodes
integer, allocatable :: n(:)
real(wp), parameter :: alpha = 0.0_wp
real(wp), parameter :: omega = 1.0_wp
real(wp) :: a,b,h,x_0,x_1,x,result
real(wp), allocatable :: integration_range(:),support_node_i(:),support_node_j(:), f_phi(:)
type(inner_product_l_obj) :: ipol
type(inner_product_obj) :: ipo
type(col), allocatable :: Row(:)

!---------------
! Logic
!---------------

write(*,*) "Main"

! ---------------------------
! Solving with linear splines
! ---------------------------

! Romberg iterations
l = 20

! Number of support nodes
! 3 = Linear Spline
! 5 = Cubic Spline
num_support_nodes = 3

! Array of number of internal nodes for the integration
allocate(n(2))

! Number of internal nodes for the integration
n(1) = 7

! Banded matrix structure
allocate(Row(n(1)))

! Inner product between f and the splines
allocate(f_phi(n(1)))

! Total number of integration nodes
! Total nodes = internal nodes + 2
num_nodes = n(1)+2

! Nodes in the integration range
! Range is given by alpha and omegra params
allocate(integration_range(num_nodes))
integration_range = linspace(alpha,omega,num_nodes)

! Stepsize
h = integration_range(2) - integration_range(1)

! Support nodes for the splines
allocate(support_node_i(num_support_nodes))
allocate(support_node_j(num_support_nodes))

! In the case of linear splines,
! first and last column will have one less
! element.
allocate(Row(1)%col(num_support_nodes-1))
allocate(Row(size(Row))%col(num_support_nodes-1))

do m = 3,num_nodes-2
  allocate(Row(m-1)%col(num_support_nodes))

  ! Construct support nodes for the spline i
  do i = 1,num_support_nodes
    support_node_i(i) = integration_range(m) + (i-1)*h
  end do

  call set_inner_product(ipo,support_node_i,&
                            f,&
                            linear_spline&
                            )

  a = support_node_i(1)
  b = support_node_i(size(support_node_i))
  result = romberg(ipo, a, b, l)
  f_phi(m-1) = result

  ! Variable o helps with indexing
  o = 1
  do j = m-1,m+1
    ! Construct support nodes for the spline j
    do i = 1,num_support_nodes
      support_node_j(i) = integration_range(j) + (i-1)*h
    end do

    write(*,*) support_node_i
    write(*,*) support_node_j
    write(*,*)

    ! Define integration limits
    a = support_node_i(1)
    b = support_node_j(size(support_node_j))

    ! Construct the inner product, according to the
    ! given formula
    call set_inner_product_l(ipol,support_node_i,support_node_j,&
                            q,k,f,&
                            linear_spline,linear_spline_deriv&
                            )

    ! Integrate inner product using Romberg's
    result = romberg(ipol, a, b, l)
    ! Save result in the banded matrix structure
    Row(m-1)%col(o) = result

    o = o + 1
  end do
end do

! Test
call pprint_band(Row,1,1)
write(*,*) f_phi

contains

! User defined functions for the problem
! to be solved go here

function q(x)
    real(wp), intent(in) :: x
    real(wp) :: q
    q = 0
  end function q

  function k(x)
    real(wp), intent(in) :: x
    real(wp) :: k
    k = 1
  end function k

  function f(x)
    real(wp), intent(in) :: x
    real(wp) :: f
    f = (pi**2)*(sin(pi*x)-9*sin(3*pi*x))
  end function f

end program main
