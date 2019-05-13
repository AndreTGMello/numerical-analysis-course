program main
use bm_ge
use utils

implicit none

type diagonal
    real, allocatable :: d(:)
endtype diagonal

!---------------
! Variables
!---------------

integer :: i,j,k,l,m,n,ub,lb
real :: b(5)
type(diagonal), allocatable :: D(:)

!---------------
! Logic
!---------------

! ub -> upper band
! lb -> lower band
ub = 2
lb = 2
k = ub+lb+1
m = 6

allocate(D(k))

! Allocate upper band diagonals
do i = ub, 1, -1
    allocate(D(i)%d(m-i))
end do
! Allocate main diagonal (i = j)
allocate(D(ub+1)%d(m))
! Allocate lower band diagonals
do i = 1, lb
    allocate(D(ub+1+i)%d(m-i))
end do

do i = 1, size(D, dim=1)
	print *, "!-------------------",char(10)
    do j = 1, size(D(i)%d(:))
        print *, D(i)%d(j)
    end do
	print *, "!-------------------",char(10)
end do

!call pprint_mat(A,m,k)
deallocate(D)
end program main
