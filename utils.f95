module utils
implicit none

! wp = working precision
! wp is set to double precision
integer, parameter :: wp = kind(1.0d0)

! Data structure for storing banded matrices
type col
  real(wp), allocatable :: Col(:)
endtype col

contains

	! Identity (==) operations won't work with floating point
	! arithmetic. In this case, the close function verifies
	! wether two real numbers are close enough (according to tol)
	function close(a,b,tol)
		real(wp), intent(in) :: a,b,tol
		logical :: close
		close = .false.
		if ( abs(a-b) < tol ) then
			close = .true.
		end if
	end function close

	! Subroutine for pretty printing matrices
	subroutine pprint_mat(matrix)
	  integer :: i,j,m,n
		real(wp), intent(in) :: matrix(:,:)
		m = size(matrix, dim=1)
		n = size(matrix, dim=2)

		! n = shape(M)
		write(*,*) "!-------------------"
		do i = 1,size(matrix,1)
			do j = 1,size(matrix,2)
				write(*,fmt="(F10.5,3X)", advance="no") matrix(i,j)
			end do
			write(*,*)
		end do
		write(*,*) "!-------------------"
	end subroutine

	! Subroutine for pretty printing banded matrices
	subroutine pprint_band(Row,lb,ub)
		integer :: i,j,ms,lb,ub
		type(col), allocatable, intent(in) :: Row(:)
		ms = size(Row)
		write(*,*)
		write(*,*) "!-------------------"
		write(*,*)
		do i = 1, ms
			if ( i < lb+1 ) then
				do j = 1, size(Row(i)%Col, dim=1)
					write(*, fmt="(F10.5, 3X)", advance="no") Row(i)%Col(j)
				end do
				do j = 1, ms - size(Row(i)%Col, dim=1)
					write(*, fmt="(F10.5, 3X)", advance="no") 0.0
				end do
			else
				do j = 1, ms - size(Row(i)%Col, dim=1)
					write(*, fmt="(F10.5, 3X)", advance="no") 0.0
				end do
				do j = 1, size(Row(i)%Col, dim=1)
					write(*, fmt="(F10.5, 3X)", advance="no") Row(i)%Col(j)
				end do
			end if
			write(*,*)
		end do
		write(*,*) "!-------------------"
		write(*,*)
	end subroutine

!TODO inner_product_wrapper

!  function inner_product_wrapper(x)
!    real, intent(out) :: inner_product_wrapper
!    inner_product_wrapper = inner_product()
!  end function inner_product_wrapper

  function inner_product(phi_i, phi_j, dphi_i, dphi_j, q, k, x)
		real(wp) :: inner_product
    real(wp), intent(in) :: x
    real(wp), external :: phi_i, phi_j, dphi_i, dphi_j, q, k

!		interface
!			function phi_i(x)
!				real, intent(in) :: x
!				real, intent(out) :: phi_i
!			end function
!		end interface
!		interface
!			function phi_j(x)
!				real, intent(in) :: x
!				real, intent(out) :: phi_j
!			end function
!		end interface


		inner_product = q(x) * dphi_i(x) * dphi_j(x) + k(x) + phi_i(x) + phi_j(x)
  end function inner_product
end module utils
