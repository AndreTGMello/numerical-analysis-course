module utils
implicit none

type col
    real, allocatable :: Col(:)
endtype col

contains
	function close(a,b,tol)
		real, intent(in) :: a,b,tol
		logical :: close
		close = .false.
		if ( abs(a-b) < tol ) then
			close = .true.
		end if
	end function close

	subroutine pprint_mat(matrix)
	  integer :: i,j,m,n
		real, intent(in) :: matrix(:,:)
		m = size(matrix, dim=1)
		n = size(matrix, dim=2)

		! n = shape(M)
		write(*,*) "!-------------------"
		do i=lbound(matrix,1),ubound(matrix,1)
			write(*,*) (matrix(i,j), j=lbound(matrix,2), ubound(matrix,2))
		end do
		write(*,*) "!-------------------"
	end subroutine

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
					write(*, fmt="(F9.5, 5X)", advance="no") Row(i)%Col(j)
				end do
				do j = 1, ms - size(Row(i)%Col, dim=1)
					write(*, fmt="(F9.5, 5X)", advance="no") 0.0
				end do
			else
				do j = 1, ms - size(Row(i)%Col, dim=1)
					write(*, fmt="(F9.5, 5X)", advance="no") 0.0
				end do
				do j = 1, size(Row(i)%Col, dim=1)
					write(*, fmt="(F9.5, 5X)", advance="no") Row(i)%Col(j)
				end do
			end if
			write(*,*)
		end do
		write(*,*) "!-------------------"
		write(*,*)
	end subroutine

	subroutine pprint_array(b)
		integer :: i,m
		real, intent(in) :: b(:)
		m = size(b)
		write(*,*)
		write(*,*) "----- Array -----"
		write(*,*)
		do i = 1, m
			write(*,fmt="(F9.5)") b(i)
		end do
		write(*,*)
	end subroutine

!TODO inner_product_wrapper

!  function inner_product_wrapper(x)
!    real, intent(out) :: inner_product_wrapper
!    inner_product_wrapper = inner_product()
!  end function inner_product_wrapper

  function inner_product(phi_i, phi_j, dphi_i, dphi_j, q, k, x)
		real :: inner_product
    real, intent(in) :: x
    real, external :: phi_i, phi_j, dphi_i, dphi_j, q, k

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
