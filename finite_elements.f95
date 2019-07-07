module finite_elements
  use utils
  use user_func
  use g_elimination
  use num_integration
  use ogpf

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

  ! Cubic spline calculation
  ! First two and last two basis
  ! need to be calculated in a
  ! different fashion from the
  ! regular ones
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
      dphi_i = dphi_i/h
      ! Basis j derivative
      dphi_j = s_eval(x,self%y,self%j,self%ds)
      dphi_j = dphi_j/h
    else
      ! Basis i
      phi_i = self%s( (x-self%y( self%i+2 ) )/h )
      ! Basis j
      phi_j = self%s( (x-self%y( self%j+2 ) )/h )

      ! Basis i deriv
      if ( close(x,self%y(1),1.0_wp*1e-3) ) then
        dphi_i = 1.0/h
      elseif ( close(x,self%y(size(self%y)),1.0_wp*1e-3) ) then
        dphi_i = -1.0/h
      elseif ( x <= self%y( self%i+2 ) ) then
        dphi_i = 1.0/h
      elseif ( x <= self%y( self%i+3 ) ) then
        dphi_i = -1.0/h
      else
        dphi_i = 0.0
        !write(*,*) "Error ", x
      end if

      ! Basis j derivative
      if ( close(x,self%y(1),1.0_wp*1e-3) ) then
        dphi_j = 1.0/h
      elseif ( close(x,self%y(size(self%y)),1.0_wp*1e-3) ) then
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
    logical, intent(in) :: flag ! True = Cubic; False = Linear
    real(wp) :: h ! Stepsize
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
        estimated_function = estimated_function + phi_i*a(i+1)
      end do
      !write(*,*)
    else
      do i = 0, size(a)-1
        phi_i = s( ( x-y( i+2 ) )/h )
        !write(*,*) "i ", i
        !write(*,*) "x ", y(i+2)
        !write(*,*) "phi_i ", phi_i
        !write(*,*) "a_i ", a(i+1)
        !write(*,*) "u_n ", estimated_function
        !write(*,*)
        estimated_function = estimated_function + phi_i*a(i+1)
      end do
      !write(*,*)
    end if
    !write(*,*)
  end function estimated_function

  subroutine solve_finite_elements(n,base_spline,base_spline_deriv,flag,romberg_iterations,error)
      integer, intent(in) :: n,romberg_iterations
      real(wp), intent(inout) :: error
      real(wp), external :: base_spline,base_spline_deriv
      ! Flag:
      ! If .true. cubic, if .false. linear
      logical, intent(in) :: flag
      integer :: i,j,l,o,m,ub,lb,ms,c
      integer :: sp_dim,total_num_nodes,num_support_nodes,total_num_support_nodes,&
      i_padding,j_padding,node_index,fix_range
      real(wp), parameter :: alpha = 0.0_wp
      real(wp), parameter :: omega = 1.0_wp
      real(wp) :: a,b,h,x,result,aux
      real(wp), allocatable :: integration_range(:),&
                                f_phi(:), coef(:),&
                                plot_expected(:),&
                                plot_calculated(:)
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
      ! Linear Spline = 3
      ! Cubic Spline = 5
      if ( flag ) then
        write(*,*) "Cubic Splines"
        write(*,*)
        ub = 3
        lb = 3
        ! Number of internal nodes
        num_support_nodes = 3
        ! Total number of nodes
        total_num_support_nodes = num_support_nodes+2
        ! Fix nodes range
        fix_range = 2
      else
        write(*,*) "Linear Splines"
        write(*,*)
        ub = 1
        lb = 1
        ! Number of internal nodes
        num_support_nodes = 1
        ! Total number of nodes
        total_num_support_nodes = num_support_nodes+2
        ! Fix nodes range
        fix_range = 1
      end if

      ! Total number of integration nodes
      ! total_num_nodes = internal nodes + 2
      total_num_nodes = n+2

      ! Dimension of the spline
      if ( flag ) then
        sp_dim = n+2
      else
        sp_dim = n
      end if

      ! Saving results for plots
      allocate(plot_expected(total_num_nodes))
      allocate(plot_calculated(total_num_nodes))

      ! Banded matrix structure
      allocate(Row(sp_dim))
      m = ub
      o = -1
      do i = 1, sp_dim
        if ( i <= lb+1 ) then
          o = o + 1
        elseif ( i > sp_dim-ub ) then
          m = m - 1
        end if
        allocate(Row(i)%col(m+o+1))
        Row(i)%col(m+o+1) = 0.0_wp
      end do

      ! Array for the result of the
      ! inner product between f and the splines
      allocate(f_phi(sp_dim))

      ! Coefficients alpha
      allocate(coef(sp_dim))
      coef = 0.0

      ! Stepsize
      h = 1.0_wp/(n+1.0_wp)
      write(*,*) "Step size"
      write(*,*) h
      write(*,*)

      ! Nodes in the integration range
      allocate(integration_range(total_num_nodes))

      do i = 0, total_num_nodes-1
        integration_range(i+1) = i*h
      end do
      write(*,*) "Integration range:"
      write(*,'(F10.5)') integration_range(1:3)
      write(*,*) "..."
      write(*,'(F10.5)') integration_range(total_num_nodes-2:total_num_nodes)

      write(*,*) "Creating banded matrix"
      write(*,*) ". . ."

      ! Index fixing
      i_padding = 0
      j_padding = 0

      do i = 0,sp_dim-1
        ! Define integration limits
        if ( flag ) then
          node_index = (i+1) - fix_range
          if ( node_index <= 0 ) then
            aux = alpha
          else
            aux = integration_range(node_index)
          end if
          a = max(aux, alpha)

          node_index = (i+1) + fix_range
          if ( node_index > total_num_nodes ) then
            aux = omega
          else
            aux = integration_range(node_index)
          end if
          b = min(aux, omega)
        else
          a = integration_range(i+1)
          b = integration_range(i+3)
        end if

        call set_inner_product(ipo,integration_range,i,n,&
                              f,base_spline,flag)
        ! Test
        !write(*,*) "i ", i
        !write(*,*) "node ", node_index
        !write(*,'(A,F10.5)') "a ", a
        !write(*,*) "node ", node_index
        !write(*,'(A,F10.5)') "b ", b
        !write(*,*)

        f_phi(i+1) = romberg(ipo, a, b, romberg_iterations)

        if ( i >= lb ) then
          i_padding = i-lb
        end if
        do j = i,min(i+num_support_nodes,sp_dim-1)

          ! Define integration limits
          if ( flag ) then
            node_index = (j+1) - fix_range
            if ( node_index <= 0 ) then
              aux = alpha
            else
              aux = integration_range(node_index)
            end if
            a = max(aux, alpha)

            node_index = (i+1) + fix_range
            if ( node_index > total_num_nodes ) then
              aux = omega
            else
              aux = integration_range(node_index)
            end if
            b = min(aux, omega)
          else
            a = integration_range(j+1)
            b = integration_range(i+3)
          end if

          ! Set values the inner product L,
          ! according to the given formula
          call set_inner_product_l(ipol,&
                              integration_range,i,j,n,&
                              k,q,&
                              base_spline,base_spline_deriv,&
                              flag)
          ! Test
          write(*,*) "j ", j
          write(*,*) "node ", node_index
          write(*,'(A,F10.5)') "a ", a
          write(*,*) "i ", i
          write(*,*) "node ", node_index
          write(*,'(A,F10.5)') "b ", b
          write(*,*)

          ! Integrate inner product using Romberg's
          ! Save result in the banded matrix structure
          aux = romberg(ipol, a, b, romberg_iterations)

          Row(i+1)%col(j+1-i_padding) = aux
          ! Because of symmetry
          if ( i /= j ) then
            if ( j >= lb ) then
              j_padding = j-lb
            end if
            Row(j+1)%col(i+1-j_padding) = aux
          end if
        end do
        ! Test
        !call pprint_band(Row,lb,ub)
        !call print_band(Row,lb,ub)
      end do
      write(*,*) "Banded matrix creation complete"
      write(*,*)

      ! Test
      !call pprint_band(Row,lb,ub)
      write(*,*) "A"
      call print_band(Row,lb,ub)
      write(*,*)
      !write(*,'(F10.5)') coef
      !write(*,*)
      write(*,*) "b"
      write(*,*)
      write(*,'(F10.5)') f_phi

      write(*,*) "Performing gaussian elimination"
      write(*,*) ". . ."
      call gaussian_elimination_banded(Row,coef,f_phi,lb,ub)
      write(*,*) "Gaussian elimination complete"
      write(*,*)

      ! Test
      !call pprint_band(Row,lb,ub)
      !call print_band(Row,lb,ub)
      !write(*,*)
      !write(*,'(F10.5)') coef
      !write(*,*)
      !write(*,'(F10.5)') f_phi
      !write(*,*)

      write(*,*)
      write(*,*) "Calculating Error`s Infinity Norm"
      write(*,*) ". . ."

      error = 0.0_wp
      do i = 1, total_num_nodes

        x = integration_range(i)
        ! Test
        !write(*,*) x
        plot_expected(i) = u(x)

        plot_calculated(i) = estimated_function(coef,integration_range,x,base_spline,flag)

        if ( abs(plot_calculated(i)-plot_expected(i)) > error ) then
          error = abs(plot_calculated(i)-plot_expected(i))
        end if

      end do
      write(*,*) "Error: ",error
      write(*,*)

      write(*,*) "Now plotting"
      write(*,*) ". . ."
      if ( flag ) then
        call gp%title('Cubic Spline')
      else
        call gp%title('Linear Spline')
      end if
      call gp%options('set key top right; set grid')
      call gp%plot(x1=integration_range,&
              y1=plot_expected,&
              ls1='title "Expected" with lines lt 1 lw 1',&
              x2=integration_range,&
              y2=plot_calculated,&
              ls2='title "Calculated" with lines lt 2 lw 1')

      write(*,*) "Plotting done"
      write(*,*)
      write(*,*) "Problem solved"
      write(*,*) "---------------"
      write(*,*)
    end subroutine solve_finite_elements
end module finite_elements
