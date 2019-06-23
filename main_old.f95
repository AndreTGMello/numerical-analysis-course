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
    procedure(f_interface), pointer, nopass :: k,q

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

  subroutine set_inner_product_l(ipo,yi,yj,k,q,phi,dphi)
    class(inner_product_l_obj), intent(inout) :: ipo
    real(wp), external :: q,k,f,phi,dphi
    real(wp), intent(in) :: yi(:),yj(:)
!    allocate(ipo%yi(size(yi)))
!    allocate(ipo%yj(size(yj)))
    ipo%yi = yi
    ipo%yj = yj
    ipo%k => k
    ipo%q => q
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

romberg_iterations = 7
! ---------------------------
! Solving with linear splines
! ---------------------------
do i = 1, 1
!  call solve_finite_elements(n(i),linear_spline,linear_spline_deriv,.false.,romberg_iterations,error)
!  error_array_linear(i) = error
end do

!call gp%plot(linspace(1.0_wp,4.0_wp,4),error_array_linear(1:4))

! ---------------------------
! Solving with cubic splines
! ---------------------------
do i = 1, 2
  call solve_finite_elements(n(i),cubic_spline,cubic_spline_deriv,.true.,romberg_iterations,error)
  error_array_cubic(i) = error
end do

!call gp%plot(linspace(1.0_wp,4.0_wp,4), error_array_cubic(1:4))

contains

! User defined functions for the problem
! to be solved go here

  function k(x)
    real(wp), intent(in) :: x
    real(wp) :: k
    k = 1
  end function k

  function q(x)
    real(wp), intent(in) :: x
    real(wp) :: q
    q = 0
  end function q

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

  subroutine solve_finite_elements(n,base_spline,base_spline_deriv,type,romberg_iterations,error)
    integer, intent(in) :: n,romberg_iterations
    real(wp), intent(inout) :: error
    real(wp), external :: base_spline,base_spline_deriv
    ! Flag:
    ! If .true. cubic, if .false. linear
    logical, intent(in) :: type
    integer :: i,j,l,o,m,ub,lb,ms,lini,lend,c
    integer :: num_basis,total_num_nodes,num_nodes_ep,num_support_nodes,total_num_support_nodes
    real(wp), parameter :: alpha = 0.0_wp
    real(wp), parameter :: omega = 1.0_wp
    real(wp) :: a,b,h,x,result
    real(wp), allocatable :: integration_range(:),support_node_i(:)&
                              ,support_node_j(:), f_phi(:), coef(:)&
                              ,support_nodes_array(:,:), plot_expected(:)&
                              ,plot_calculated(:)
    type(inner_product_l_obj) :: ipol
    type(inner_product_obj) :: ipo
    type(col), allocatable :: Row(:)
    type(gpf):: gp

    write(*,*) "-------------------------"
    write(*,*) "Number of internal nodes:"
    write(*,*) n
    write(*,*)
    write(*,*) "Solving with:"

    ! Number of upper and lower bands
    ! ub = lb = 1 iff linear splines
    ! ub = lb = 3 iff cubic splines
    ! Number of support nodes
    ! 3 = Linear Spline
    ! 5 = Cubic Spline
    if ( type ) then
      write(*,*) "Cubic Splines"
      write(*,*)
      ub = 3
      lb = 3
      ! Number of internal nodes
      num_support_nodes = 3
      ! Total number of nodes
      total_num_support_nodes = num_support_nodes+2
    else
      write(*,*) "Linear Splines"
      write(*,*)
      ub = 1
      lb = 1
      ! Number of internal nodes
      num_support_nodes = 1
      ! Total number of nodes
      total_num_support_nodes = num_support_nodes+2
    end if

    ! Total number of integration nodes
    ! total_num_nodes = internal nodes + 2
    total_num_nodes = n+2

    ! Saving results for plots
    allocate(plot_expected(total_num_nodes))
    allocate(plot_calculated(total_num_nodes))

    ! Number of b-spline basis
    num_basis = total_num_nodes+(num_support_nodes-1)

    ! Extended partition
    num_nodes_ep = total_num_nodes+2*num_support_nodes

    ! Banded matrix structure
    allocate(Row(num_basis))

    ! Array for the result of the
    ! inner product between f and the splines
    allocate(f_phi(num_basis))

    ! Coefficients alpha
    allocate(coef(num_basis))

    ! Stepsize
    h = 1.0_wp/(n+1.0_wp)
    write(*,*) "Step size"
    write(*,*) h
    write(*,*)

    ! Nodes in the integration range
    allocate(integration_range(num_nodes_ep))

    do i = 1, num_nodes_ep
      integration_range(i) = (i-total_num_support_nodes+1)*h
    end do

    ! Support nodes for the splines
    allocate(support_node_i(total_num_support_nodes))
    allocate(support_node_j(total_num_support_nodes))
    allocate(support_nodes_array(num_basis,total_num_support_nodes))

    write(*,*) "Creating banded matrix"
    write(*,*) ". . ."
    do m = 1,num_basis
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

      a = max(support_node_i(1),0.0_wp)
      b = min(support_node_i(size(support_node_i)),1.0_wp)
      f_phi(m) = romberg(ipo, a, b, romberg_iterations)

      ! Calculating how many support nodes
      ! phi_j shall be used in the L inner
      ! product agains phi_i
      lini = -num_support_nodes
      lend = num_support_nodes
      if ( m < lb+1 ) then
        lini = 0+(1-m)
        lend = num_support_nodes
      elseif ( m > num_basis-ub ) then
        lini = -num_support_nodes
        lend = num_support_nodes-(ub+m-num_basis)
      end if
      allocate(Row(m)%col(2*(num_support_nodes)+1-(abs(lini+lend))))

      ! Variable o helps with indexing
      o = 1
      do j = lini,lend
        ! Set support node values for the spline j
        do i = 1,total_num_support_nodes
          support_node_j(i) = integration_range(m+j) + (i-1)*h
        end do

        ! Test
        !write(*,*) support_node_i
        !write(*,*) support_node_j
        !write(*,*)

        ! Set values the inner product,
        ! according to the given formula
        call set_inner_product_l(ipol,support_node_i,support_node_j,&
                            k,q,&
                            base_spline,base_spline_deriv&
                            )

        ! Define integration limits
        a = min(support_node_i(1), support_node_j(1))
        b = max(support_node_i(total_num_support_nodes),&
                support_node_j(total_num_support_nodes))

        ! Integrate inner product using Romberg's
        ! Save result in the banded matrix structure
        Row(m)%col(o) = romberg(ipol, a, b, romberg_iterations)

        o = o + 1
      end do
    end do
    write(*,*) "Banded matrix creation complete"
    write(*,*)

    ! Test
    call pprint_band(Row,lb,ub)
    write(*,*)
    write(*,'(F10.5)') coef
    write(*,*)
    write(*,'(F10.5)') f_phi
    write(*,*)

    write(*,*) "Performing gaussian elimination"
    write(*,*) ". . ."
    call gaussian_elimination_banded(Row,coef,f_phi,lb,ub)
    write(*,*) "Gaussian elimination complete"
    write(*,*)

    ! Test
    call pprint_band(Row,lb,ub)
    write(*,*)
    write(*,'(F10.5)') coef
    write(*,*)
    write(*,'(F10.5)') f_phi
    write(*,*)

    write(*,*) "Calculating Error`s Infinity Norm"
    write(*,*) ". . ."

    error = 0.0_wp
    do i = 1, total_num_nodes

      x = integration_range(i+num_support_nodes)
      ! Test
      !write(*,*) x
      plot_expected(i) = u(x)

      plot_calculated(i) = estimated_function(coef,x,support_nodes_array,base_spline)

      if ( abs(plot_calculated(i)-plot_expected(i)) > error ) then
        error = abs(plot_calculated(i)-plot_expected(i))
      end if

    end do
    write(*,*) "Error: ",error
    write(*,*)

    write(*,*) "Now plotting"
    write(*,*) ". . ."
    if ( type ) then
      call gp%title('Cubic Spline')
    else
      call gp%title('Linear Spline')
    end if
    call gp%options('set key top right; set grid')
    call gp%plot(x1=integration_range(1+num_support_nodes:num_nodes_ep-num_support_nodes),&
            y1=plot_expected,&
            ls1='title "Expected" with lines lt 1 lw 1',&
            x2=integration_range(1+num_support_nodes:num_nodes_ep-num_support_nodes),&
            y2=plot_calculated,&
            ls2='title "Calculated" with lines lt 2 lw 1')

    write(*,*) "Plotting done"
    write(*,*)
    write(*,*) "Problem solved"
    write(*,*) "---------------"
    write(*,*)
  end subroutine solve_finite_elements
end program main
