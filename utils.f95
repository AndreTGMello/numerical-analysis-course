module utils
implicit none

! wp = working precision
! wp is set to double precision
integer, parameter :: wp = kind(1.0d0)
real(wp), parameter :: pi = 4.D0*datan(1.D0)

! Data structure for storing banded matrices
type col
  real(wp), allocatable :: Col(:)
endtype col

! Abstract `function` object
! This will be used to make the numerical
! Integration procedure more flexible
type, abstract :: fun
contains
  procedure(eval_iface), deferred :: eval
end type fun

interface
  function eval_iface(self, x)
    import :: fun
		import :: wp
    class(fun), intent(in) :: self
    real(wp), intent(in) :: x
    real(wp) :: eval_iface
  end function eval_iface
end interface

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
		integer :: i,j,pad,ms,lb,ub
		type(col), allocatable, intent(in) :: Row(:)
		ms = size(Row)
		pad = 0
		write(*,*)
		write(*,*) "!-------------------"
		write(*,*)
		do i = 1, ms
			if ( i <= lb+1 ) then
				do j = 1, size(Row(i)%col)
					write(*, fmt="(F10.5, 3X)", advance="no") Row(i)%Col(j)
				end do
				do j = 1, ms - size(Row(i)%col)
					write(*, fmt="(F10.5, 3X)", advance="no") 0.0
				end do
			else
				pad = pad + 1
				do j = 1, pad !ms - size(Row(i)%Col, dim=1)
					write(*, fmt="(F10.5, 3X)", advance="no") 0.0
				end do
				do j = 1, size(Row(i)%col)
					write(*, fmt="(F10.5, 3X)", advance="no") Row(i)%Col(j)
				end do
				do j = pad+size(Row(i)%col), ms-1
					write(*, fmt="(F10.5, 3X)", advance="no") 0.0
				end do
			end if
			write(*,*)
		end do
		write(*,*) "!-------------------"
		write(*,*)
	end subroutine

end module utils
