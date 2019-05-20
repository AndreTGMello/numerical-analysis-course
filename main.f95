program main
use bm_ge
use utils

implicit none

type row
    real, allocatable :: r(:)
endtype row

!---------------
! Variables
!---------------

integer :: i,j,k,l,m,n,ub,lb
real :: b(5), multiplier
type(row), allocatable :: A(:)

!---------------
! Logic
!---------------

! ub -> upper band
! lb -> lower band
ub = 2
lb = 2
! Square matrix size
m = 4

allocate(A(m))

do i = 1, m
    if ( i == 1 .or. i == m ) then
        allocate(A(i)%r(ub+1))
    else
        allocate(A(i)%r(lb+ub+1))
    end if
end do

do i = 1, size(A, dim=1)
	print *, "!-------------------",char(10)
    do j = 1, size(A(i)%r(:))
        print *, A(i)%r(j)
    end do
	print *, "!-------------------",char(10)
end do

!call pprint_mat(A,m,k)
do i=1, size(A)
    deallocate(A(i)%r)
end do
deallocate(A)

end program main
