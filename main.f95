program main
use bm_ge
use utils

implicit none

!---------------
! Variables
!---------------
integer :: i,j,k,l,m,n
real, dimension(5) :: b
real, dimension(:,:), allocatable :: A

!---------------
! Logic
!---------------
k = 5
m = 6
allocate(A(m,k))
A = reshape((/0,0,8,13,18,22,0,4,9,14,19,23,1,5,10,15,20,24,2,6,11,16,21,0,3,7,12,17,0,0/), shape(A))

call pprint_mat(A,m,k)
deallocate(A)
end program main
