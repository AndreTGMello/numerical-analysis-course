module utils
implicit none

type col
    real, allocatable :: Col(:)
endtype col

contains
	subroutine pprint_mat(matrix,m,n)
	  integer :: i,j,m,n
		real, dimension(m,n) :: matrix

		! n = shape(M)
		print *, "!-------------------"
		do i=lbound(matrix,1),ubound(matrix,1)
			write(*,*) (matrix(i,j), j=lbound(matrix,2), ubound(matrix,2))
		end do
		print *, "!-------------------",char(10)
	end subroutine

	subroutine pprint_band(Row,ms,lb,ub)
		integer :: i,j,ms,lb,ub
		type(col) :: Row(ms)
		write(*,*)
		write(*,*) "----- Matrix -----"
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
		write(*,*)
	end subroutine

	subroutine pprint_array(b,m)
		integer :: i,m
		real :: b(m)

		write(*,*)
		write(*,*) "----- Array -----"
		write(*,*)
		do i = 1, m
			write(*,fmt="(F9.5)") b(i)
		end do
		write(*,*)
	end subroutine
end module utils
