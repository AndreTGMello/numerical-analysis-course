program bm_ge
  ! Program for banded matrix gauss elimination
  use utils

  implicit none

  !---------------
  ! Variables
  !---------------
  real, dimension(:,:), allocatable :: matrix
  integer :: i,j,m,n,k

  !---------------
  ! Procedures
  !---------------
  print *, "Enter square matrix number of rows/columns:"
  read *, m
  print *, "Enter number of diagonals:"
  read *, n
  allocate(matrix(m,n))

  do i=1,m
    do j=1,n
      matrix(i,j)=i+j
    end do
  end do

  call pprint_mat(matrix,m,n)

  deallocate(matrix)
end program bm_ge
