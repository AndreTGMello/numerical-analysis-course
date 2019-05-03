module utils
implicit none

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

end module utils
