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

  function estimated_function(alpha,x,support_nodes,spline_base)
    real(wp), external :: spline_base
    real(wp), intent(in) :: alpha(:),support_nodes(:,:)
    real(wp), intent(in) :: x
    real(wp) :: estimated_function
    integer :: i
    estimated_function = 0
    do i = 1, size(alpha)
      estimated_function = estimated_function + alpha(i)*spline_base(support_nodes(i,1:),x)
    end do
  end function estimated_function

end module user_defined_functions

program main
use utils
!use bm_ge
use g_elimination
use num_integration
use bspline
use user_defined_functions
use ogpf

implicit none


!---------------
! Variables
!---------------

integer :: i,j,l,o,m,ub,lb,ms,lini,lend,c,romberg_iterations
integer :: num_nodes
integer :: num_support_nodes
integer, allocatable :: n(:)
real(wp), parameter :: alpha = 0.0_wp
real(wp), parameter :: omega = 1.0_wp
real(wp) :: a,b,h,x_0,x_1,x,result
real(wp), allocatable :: integration_range(:),support_node_i(:)&
                          ,support_node_j(:), f_phi(:), coef(:)&
                          ,support_nodes_array(:,:), plot_expected(:)&
                          ,plot_calculated(:)
type(inner_product_l_obj) :: ipol
type(inner_product_obj) :: ipo
type(col), allocatable :: Row(:)
type(gpf):: gp

!---------------
! Logic
!---------------

write(*,*) "Main"

! Array of number of internal nodes for the integration
allocate(n(4))

! Number of internal nodes for the integration
n(1:) = [7, 15, 31, 63]

romberg_iterations = 11
! ---------------------------
! Solving with linear splines
! ---------------------------
do c = 1, size(n)
  call solve_finite_elements(n(c),linear_spline,linear_spline_deriv,.false.,3)
end do

! ---------------------------
! Solving with cubic splines
! ---------------------------
do c = 1, size(n)
  call solve_finite_elements(n(c),cubic_spline,cubic_spline_deriv,.true.,romberg_iterations)
end do

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

  function u(x)
    real(wp), intent(in) :: x
    real(wp) :: u
    u = sin(pi*x) - sin(3*pi*x)
  end function u

  subroutine solve_finite_elements(n,base_spline,base_spline_deriv,type,romberg_iterations)
    integer, intent(in) :: n,romberg_iterations
    real(wp), external :: base_spline,base_spline_deriv
    ! Flag:
    ! If .true. cubic, if .false. linear
    logical, intent(in) :: type
    integer :: i,j,l,o,m,ub,lb,ms,lini,lend,c
    integer :: total_num_nodes,num_support_nodes,total_num_support_nodes
    real(wp), parameter :: alpha = 0.0_wp
    real(wp), parameter :: omega = 1.0_wp
    real(wp) :: a,b,h,x_0,x_1,x,result
    real(wp), allocatable :: integration_range(:),support_node_i(:)&
                              ,support_node_j(:), f_phi(:), coef(:)&
                              ,support_nodes_array(:,:), plot_expected(:)&
                              ,plot_calculated(:)
    type(inner_product_l_obj) :: ipol
    type(inner_product_obj) :: ipo
    type(col), allocatable :: Row(:)
    type(gpf):: gp

    ! Number of upper and lower bands
    ! ub = lb = 1 iff linear splines
    ! ub = lb = 3 iff cubic splines
    ! Number of support nodes
    ! 3 = Linear Spline
    ! 5 = Cubic Spline
    if ( type ) then
      write(*,*) "Solving with cubic splines"
      ub = 3
      lb = 3
      ! Number of internal nodes
      num_support_nodes = 3
      ! Total number of nodes
      total_num_support_nodes = num_support_nodes+2
    else
      write(*,*) "Solving with linear splines"
      ub = 1
      lb = 1
      ! Number of internal nodes
      num_support_nodes = 1
      ! Total number of nodes
      total_num_support_nodes = num_support_nodes+2
    end if

    ! Banded matrix structure
    allocate(Row(n))

    ! Array for the result of the
    ! inner product between f and the splines
    allocate(f_phi(n))

    ! Coeficients alpha
    allocate(coef(n))

    ! Total number of integration nodes
    ! total_num_nodes = internal nodes + 2
    total_num_nodes = n+2

    ! Saving results for plots
    allocate(plot_expected(total_num_nodes))
    allocate(plot_calculated(total_num_nodes))

    ! Nodes in the integration range
    ! Range is given by alpha and omegra params
    allocate(integration_range(total_num_nodes))
    integration_range = linspace(alpha,omega,total_num_nodes)

    ! Stepsize
    h = integration_range(2) - integration_range(1)

    ! Support nodes for the splines
    allocate(support_node_i(total_num_support_nodes))
    allocate(support_node_j(total_num_support_nodes))
    allocate(support_nodes_array(n,total_num_support_nodes))

    do m = 1,n
      ! Construct support nodes for the spline i
      do i = 1,total_num_support_nodes
        support_node_i(i) = integration_range(m) + (i-1)*h
      end do

      ! Save this support node array for later use
      ! at the `estimated_function` method
      support_nodes_array(m,1:) = support_node_i(1:)

      call set_inner_product(ipo,support_node_i,&
                            f,&
                            base_spline&
                            )

      a = support_node_i(1)
      b = support_node_i(size(support_node_i))
      result = romberg(ipo, a, b, romberg_iterations)
      f_phi(m) = result

      ! Calculatin how many support nodes
      ! phi_j shall be used in the L inner
      ! product agains phi_i
      lini = -num_support_nodes
      lend = num_support_nodes
      if ( m < lb+1 ) then
        lini = 0+(1-m)
        lend = num_support_nodes
      elseif ( m > n-lb ) then
        lini = -num_support_nodes
        lend = num_support_nodes+(lb+1-m)
      end if
      allocate(Row(m)%col(2*(num_support_nodes)+1-(abs(lini+lend))))

      ! Variable o helps with indexing
      o = 1
      do j = lini,lend
        ! Construct support nodes for the spline j
        do i = 1,total_num_support_nodes
          support_node_j(i) = integration_range(m+j) + (i-1)*h
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
                                base_spline,base_spline_deriv&
                                )

        ! Integrate inner product using Romberg's
        result = romberg(ipol, a, b, romberg_iterations)
        ! Save result in the banded matrix structure
        Row(m)%col(o) = result

        o = o + 1
      end do
    end do

    ! Test
    call pprint_band(Row,lb,ub)

    call gaussian_elimination_banded(Row,coef,f_phi,1,1)

    call pprint_band(Row,lb,ub)


    do i = 1, size(integration_range)
      plot_calculated(i) = estimated_function(f_phi,integration_range(i),support_nodes_array,base_spline)
      plot_expected(i) = u(integration_range(i))
    end do

    call gp%title('Test')
    call gp%options('set key top right; set grid')
    call gp%plot(x1=integration_range,y1=plot_expected,ls1='title "Expected" with lines lt 1 lw 1',&
    x2=integration_range,y2=plot_calculated,ls2='title "Calculated" with lines lt 2 lw 1')

  end subroutine solve_finite_elements
end program main
