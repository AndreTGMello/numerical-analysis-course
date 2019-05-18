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
real :: b(5), multiplier
type(diagonal), allocatable :: D(:)

!---------------
! Logic
!---------------

! ub -> upper band
! lb -> lower band
ub = 2
lb = 2
! Square matrix size
m = 4

allocate(D(1+lb+ub))

! Allocate main diagonal (i = j)
allocate(D(1)%d(m))
! Allocate lower band diagonals
do i = 1, lb
    allocate(D(i+1)%d(m-i))
end do
! Allocate upper band diagonals
do i = 1, ub
    allocate(D(lb+1+i)%d(m-i))
end do

do i = 1, size(D, dim=1)
	print *, "!-------------------",char(10)
    do j = 1, size(D(i)%d(:))
        print *, D(i)%d(j)
    end do
	print *, "!-------------------",char(10)
end do

! Gauss Elimination
do i = 1, m-1
    ! For each lower diagonal, do:
    do j = 2, lb+1
        multiplier = D(j)%d(i)/D(1)%d(i)
        D(j)%d(i) = 0
        ! For each upper diagonal, subtract from the diagonal below:
        do k = 1+lb+1, 1+lb+ub
            D(mod(k-1,lb+1))%(i+1) = D(mod())%d(i+1) - multiplier*D(k)%d(i)
        end do
    end do
end do

!call pprint_mat(A,m,k)
do i=1, size(D)
    deallocate(D(i)%d)
end do
deallocate(D)

end program main
