module bspline
  use num_integration
  use utils, only : wp

  implicit none
  
  abstract interface
    function f_interface(x)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: f_interface
    end function f_interface
  end interface

  type, extends(fun) :: inner_product_l_obj
    integer :: i,j,n
    logical :: flag ! To decide weather linear or cubic spline
    real(wp), allocatable :: y(:) ! Nodes on the integration range

    ! Pointers to k and q function, as described by
    ! the formula.
    procedure(f_interface), pointer, nopass :: k,q

    ! Pointers to spline basis function (and derivative)
    procedure(f_interface), pointer, nopass :: s,ds

  contains
    procedure :: eval => inner_product_l_eval
  end type inner_product_l_obj

  type, extends(fun) :: inner_product_obj
    integer :: i,n
    logical :: flag ! To decide weather linear or cubic spline
    real(wp), allocatable :: y(:) ! Nodes on the integration range

    procedure(f_interface), pointer, nopass :: f

    ! Pointers to linear_spline functions
    procedure(f_interface), pointer, nopass :: s

  contains
    procedure :: eval => inner_product_eval
  end type inner_product_obj

contains
  function s_eval(x,y,i,s)
    real(wp) :: s_eval
    real(wp) :: h
    real(wp), intent(in) :: x,y(:)
    real(wp), external :: s
    integer, intent(in) :: i
    integer :: n

    n = size(y) - 2
    h = y(2) - y(1)

    if ( i == 0 ) then
      s_eval = s( x/h ) - 4*s( (x + h)/h )
    elseif ( i == 1 ) then
      s_eval = s( (x - y(2))/h ) - s( (x + h)/h )
    elseif ( i == n ) then
      s_eval = s( (x - y( n + 1 ))/h )&
              - s( (x - (n + 2)*h)/h )
    elseif ( i == n + 1 ) then
      s_eval = s( (x-y( n + 2 ))/h )&
              - 4*s( (x- (n + 2)*h)/h )
    else
      s_eval = s( (x - y(i + 1) )/h )
    end if

  end function s_eval

  function inner_product_l_eval(self, x)
    class(inner_product_l_obj), intent(in) :: self
    real(wp), intent(in) :: x
    real(wp) :: inner_product_l_eval
    real(wp) :: phi_i,phi_j,dphi_i,dphi_j,h
    integer :: k

    h = self%y(2) - self%y(1)
    if ( self%flag ) then
      ! Basis i
      phi_i = s_eval(x,self%y,self%i,self%s)
      ! Basis i
      phi_j = s_eval(x,self%y,self%j,self%s)
      ! Basis i derivative
      dphi_i = s_eval(x,self%y,self%i,self%ds)
      ! Basis j derivative
      dphi_j = s_eval(x,self%y,self%j,self%ds)
    else
      ! Basis i
      phi_i = self%s( (x-self%y( self%i+2 ) )/h )
      ! Basis j
      phi_j = self%s( (x-self%y( self%j+2 ) )/h )

      ! Basis i deriv
      if ( x == 0 ) then
        dphi_i = 1.0/h
      elseif ( x == size(self%y) ) then
        dphi_i = -1.0/h
      elseif ( x <= self%y( self%i+2 ) ) then
        dphi_i = 1.0/h
      elseif ( x <= self%y( self%i+3 ) ) then
        dphi_i = -1.0/h
      else
        dphi_i = 0.0/h
        !write(*,*) "Error ", x
      end if

      ! Basis j derivative
      if ( x == 0 ) then
        dphi_j = 1.0/h
      elseif ( x == size(self%y) ) then
        dphi_j = -1.0/h
      elseif ( x <= self%y( self%j+2 ) ) then
        dphi_j = 1.0/h
      elseif ( x <= self%y( self%j+3 ) ) then
        dphi_j = -1.0/h
      else
        dphi_j = 0.0/h
        !write(*,*) "Error ", x
      end if

    end if
    !write(*,*) "x ", x
    !write(*,*) "phi_i ", phi_i
    !write(*,*) "phi_j ", phi_j
    !write(*,*) "dphi_j ", dphi_i
    !write(*,*) "dphi_j ", dphi_j

    inner_product_l_eval = &
      self%k(x) * dphi_i * dphi_j &
      + self%q(x) * phi_i * phi_j

  end function inner_product_l_eval

  function inner_product_eval(self, x)
    class(inner_product_obj), intent(in) :: self
    real(wp), intent(in) :: x
    real(wp) :: inner_product_eval,h,phi_i
    integer :: k

    h = self%y(2) - self%y(1)
    if ( self%flag ) then
      phi_i = s_eval(x,self%y,self%i,self%s)

    else
      phi_i = self%s( ( x-self%y( self%i+2 ) )/h )

    end if

    inner_product_eval =  &
      self%f(x)*phi_i
  end function inner_product_eval

  subroutine set_inner_product_l(ipo,y,i,j,n,k,q,s,ds,flag)
    class(inner_product_l_obj), intent(inout) :: ipo
    real(wp), external :: q,k,s,ds
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: i,j,n
    logical, intent(in) :: flag

    ipo%y = y
    ipo%i = i
    ipo%j = j
    ipo%n = n
    ipo%k => k
    ipo%q => q
    ipo%s => s
    ipo%ds => ds
    ipo%flag = flag
  end subroutine set_inner_product_l

  subroutine set_inner_product(ipo,y,i,n,f,s,flag)
    class(inner_product_obj), intent(inout) :: ipo
    real(wp), external :: f,s
    real(wp), intent(in) :: y(:)
    integer, intent(in) :: i,n
    logical, intent(in) :: flag

    ipo%y = y
    ipo%i = i
    ipo%n = n
    ipo%f => f
    ipo%s => s
    ipo%flag = flag
  end subroutine set_inner_product

  function estimated_function(a,y,x,s,flag)
    real(wp), intent(in) :: a(:) ! Coefficient vector
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y(:) ! Integration nodes
    real(wp), external :: s ! Spline basis function
    logical, intent(in) :: flag
    real(wp) :: h
    real(wp) :: estimated_function,phi_i
    integer :: i,n

    h = y(2) - y(1)
    estimated_function = 0.0

    if ( flag ) then
      do i = 0, size(a)-1
        phi_i = s_eval(x,y,i,s)

        ! Test
        !write(*,*) i
        !write(*,*) "phi_i ", phi_i
        !write(*,*) "a_i ", a(i+1)
        estimated_function = estimated_function + a(i+1)*phi_i
      end do
      !write(*,*)
    else
      do i = 0, size(a)-1
        phi_i = s( ( x-y( i+2 ) )/h )
        !write(*,*) i
        !write(*,*) "phi_i ", phi_i
        !write(*,*) "a_i ", a(i+1)
        estimated_function = estimated_function + phi_i*a(i+1)
      end do
    end if
    !write(*,*)
  end function estimated_function

  function linear_spline(y,x)
    !---------------
    ! Variables
    !---------------
    ! y are the three points on top of which
    ! the linear spline sits. Everywhere else
    ! it is zero.
    ! x is used for computing the desired value
    ! this spline assumes in this given point.
    real(wp), intent(in) :: y(3),x
    real(wp) :: h
    real(wp) :: linear_spline
    h = y(2) - y(1)
    linear_spline = 0

    if ( x < y(1) ) then
      linear_spline = 0
    elseif ( x <= y(2) ) then
      linear_spline = (x - y(1))/h
    elseif ( x <= y(3) ) then
      linear_spline =  (y(3) - x)/h
    end if
  end function linear_spline

  function linear_spline_deriv(y,x)
    !---------------
    ! Variables
    !---------------
    ! y are the three points on top of which
    ! the linear spline sits. Everywhere else
    ! it is zero.
    ! x is used for computing the desired value
    ! this spline assumes in this given point.
    real(wp), intent(in) :: y(3),x
    real(wp) :: h
    real(wp) :: linear_spline_deriv
    h = y(2) - y(1)
    linear_spline_deriv = 0.0_wp

    if ( x < y(1) ) then
      linear_spline_deriv = 0
    elseif ( x <= y(2) ) then
      linear_spline_deriv = 1.0_wp/h
    elseif ( x <= y(3) ) then
      linear_spline_deriv =  -1.0_wp/h
    end if
  end function linear_spline_deriv

  function cubic_spline(y, x)
    !---------------
    ! Variables
    !---------------
    ! y are the five points on top of which
    ! the cubi spline sits. Everywhere else
    ! it is zero.
    ! x is used for computing the desired value
    ! this spline assumes in this given point.
    real(wp), intent(in) :: y(5), x
    real(wp) :: h
    integer :: i,j
    real(wp) :: cubic_spline

    h = y(2) - y(1)
    if ( x < y(1) .or. x > y(5)) then
      cubic_spline = 0

    elseif ( x >= y(1) .and. x < y(2) ) then
      cubic_spline = ((x-y(1))**3)/(6*h**3)

    elseif ( x >= y(2) .and. x < y(3) ) then
      cubic_spline = &
      ( ((x-y(1))**2)*(y(3)-x) &
          + (x-y(1))*(y(4)-x)*(x-y(2))   &
      + ( (y(5)-x)*(x-y(2))**2 ) )/(6*h**3)

    elseif ( x >= y(3) .and. x < y(4) ) then
      cubic_spline = &
        ( (y(5)-x)*(x-y(2))*(y(4)-x) &
          + ((y(5)-x)**2)*(x-y(3))  &
      + (x-y(1))*(y(4)-x)**2 )/(6*h**3)

    elseif ( x >= y(4) .and. x < y(5) ) then
      cubic_spline = ((y(5)-x)**3)/(6*h**3)

    end if

  end function cubic_spline

  function cubic_spline_deriv(y, x)
    !---------------
    ! Variables
    !---------------
    ! y are the five points on top of which
    ! the cubi spline sits. Everywhere else
    ! it is zero.
    ! x is used for computing the desired value
    ! this spline assumes in this given point.
    real(wp), intent(in) :: y(5), x
    real(wp) :: h
    integer :: i,j
    real(wp) :: cubic_spline_deriv

    h = y(2) - y(1)
    if ( x < y(1) .or. x > y(5)) then
      cubic_spline_deriv = 0.0_wp

    elseif ( x >= y(1) .and. x < y(2) ) then
      cubic_spline_deriv = ((x-y(1))**2)/(2*h**3)

    elseif ( x >= y(2) .and. x < y(3) ) then
      cubic_spline_deriv = &
        - (9*x**2 - ( 2*y(5) + 2*y(4) + 2*y(3) + 6*y(2) + 6*y(1))*x &
        + 2*y(2)*y(5) &
        +( y(2) + y(1) )*y(4) + 2*y(1)*y(3) &
        + y(2)**2 + y(1)*y(2) + y(1)**2) / (6*h**3)

    elseif ( x >= y(3) .and. x < y(4) ) then
      cubic_spline_deriv = &
        (9*x**2 - ( 6*y(5) + 6*y(4) + 2*y(3) + 2*y(2) + 2*y(1) )*x &
        + y(5)**2 &
        + ( y(4)+2*y(3)+y(2) )*y(5) + y(4)**2 &
        + ( y(2)+2*y(1) )*y(4)) / (6*h**3)

    elseif ( x >= y(4) .and. x < y(5) ) then
      cubic_spline_deriv = (-(y(5)-x)**2)/(2*h**3)

    end if
  end function cubic_spline_deriv


  !--------------------
  ! New implementation
  !--------------------
  function cubic_spline_basis(x)
    real(wp) :: cubic_spline_basis
    real(wp), intent(in) :: x

    if ( x <= -2.0 ) then
      cubic_spline_basis = 0
    elseif ( -2.0 < x .and. x <= -1.0 ) then
      cubic_spline_basis = (1.0/4.0)*(2.0+x)**3
    elseif ( -1.0 < x .and. x <= 0.0 ) then
      cubic_spline_basis = (1.0/4.0)*((2.0+x)**3 - 4*(1.0+x)**3)
    elseif ( 0.0 < x .and. x <= 1.0 ) then
      cubic_spline_basis = (1.0/4.0)*((2.0-x)**3 - 4.0*(1-x)**3)
    elseif ( 1.0 < x .and. x <= 2.0 ) then
      cubic_spline_basis = (1.0/4.0)*(2.0-x)**3
    elseif( 2.0 < x ) then
      cubic_spline_basis = 0
    end if

  end function cubic_spline_basis

  function cubic_spline_basis_deriv(x)
    real(wp) :: cubic_spline_basis_deriv
    real(wp), intent(in) :: x

    if ( x <= -2.0 ) then
      cubic_spline_basis_deriv = 0.0
    elseif ( -2.0 < x .and. x <= -1.0 ) then
      cubic_spline_basis_deriv = (1.0/4.0)*3.0*(2.0+x)**2
    elseif ( -1.0 < x .and. x <= 0.0 ) then
      cubic_spline_basis_deriv = -(1.0/4.0)*((9*x**2)+(12*x))
    elseif ( 0.0 < x .and. x <= 1.0 ) then
      cubic_spline_basis_deriv = (1.0/4.0)*(3*x*(3*x-4.0))
    elseif ( 1.0 < x .and. x <= 2.0 ) then
      cubic_spline_basis_deriv = -(1.0/4.0)*3.0*(2.0-x)**2
    elseif( 2.0 < x ) then
      cubic_spline_basis_deriv = 0.0
    end if

  end function cubic_spline_basis_deriv

  function linear_spline_basis(x)
    real(wp) :: linear_spline_basis
    real(wp), intent(in) :: x

    if ( x <= -1.0 .or. x >= 1.0 ) then
      linear_spline_basis = 0.0
    elseif ( x < 0.0 ) then
      linear_spline_basis = (x+1.0)
    elseif ( x < 1.0 ) then
      linear_spline_basis = (1.0-x)
    end if

  end function linear_spline_basis

  function linear_spline_basis_deriv(x)
    real(wp) :: linear_spline_basis_deriv
    real(wp), intent(in) :: x

    if ( x <= -1.0 .or. x >= 1.0 ) then
      linear_spline_basis_deriv = 0.0
    elseif ( x < 0.0 ) then
      linear_spline_basis_deriv = x
    elseif ( x < 1.0 ) then
      linear_spline_basis_deriv = -x
    end if

  end function linear_spline_basis_deriv

end module bspline
