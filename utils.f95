module utils
implicit none

contains
	subroutine pprint_mat(matrix,m,n)
	  integer :: i,j,m,n
		real, dimension(m,n) :: matrix

		! n = shape(M)

		do i=lbound(matrix,1),ubound(matrix,1)
			write(*,*) (matrix(i,j), j=lbound(matrix,2), ubound(matrix,2))
		end do
	end subroutine

end module utils
