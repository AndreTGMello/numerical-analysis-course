program main
use bm_ge
use g_elimination
use utils

implicit none


!---------------
! Variables
!---------------

integer :: i,j,k,l,m,n,ub,lb,ms
real :: multiplier
real, allocatable :: b(:), x(:), y(:), A(:,:)
type(col), allocatable :: Row(:)

!---------------
! Logic
!---------------

! ub -> upper band
! lb -> lower band
ub = 2
lb = 4
! Square matrix size
ms = 6

allocate(b(ms))
allocate(x(ms))
allocate(Row(ms))

do i = 1, ms
    if ( i == 1 ) then
        allocate(Row(i)%Col(ub+1))
    elseif ( i == ms ) then
        allocate(Row(i)%Col(lb+1))
    else
        allocate(Row(i)%Col(lb+ub+1))
    end if
end do

!Row(1)%Col = [8, 1]
!Row(2)%Col = [2, 10, 4]
!Row(3)%Col = [5, 10, 1]
!Row(4)%Col = [2, 3]
!b = [1, 3, 6, 20]

Row(1)%Col = [15, 7, 7]
Row(2)%Col = [2, 7, 1, 2]
Row(3)%Col = [5, 1, 15, 2, 1]
Row(4)%Col = [3, 2, 3, 14, 1, 2]
Row(5)%Col = [1, 3, 2, 4, 12, 3]
Row(6)%Col = [2, 3, 1, 5, 13]
b = [2, 1, 2, 4, 2, 6]

call pprint_band(Row,ms,lb,ub)
call pprint_array(b,size(b))

! Gaussian Elimination
do i = 1, ms-1
    ! For each row:
    j = i + 1
    do while ( j <= i+lb .and. j <= ms)
        l = 0
        m = 0
        ! If row is below the number of lower + uper bands,
        ! adjust index
        if ( i > lb+1 ) then
            l = ms - size(Row(i)%Col)
            m = ms - size(Row(j)%Col)
        elseif ( j > lb+1 ) then
            m = ms - size(Row(j)%Col)
        end if
        multiplier = Row(j)%Col(i-m)/Row(i)%Col(i-l)
        Row(j)%Col(i-m) = 0
        b(j) = b(j) - multiplier*b(i)
        ! For each column:
        k = i + 1
        do while ( k <= i+ub .and. k <= ms )
            Row(j)%Col(k-m) = Row(j)%Col(k-m) - multiplier*Row(i)%Col(k-l)
            k = k + 1
        end do
        j = j + 1
    end do
end do

! Solving the system
! k and l are auxiliary variables for index fixing/shifting
! due to the way the band matrix is stored.
l = ub + 1
do i = ms, ms-ub, -1
    x(i) = b(i)
    k = ms
    do j = size(Row(i)%Col, dim=1), size(Row(i)%Col, dim=1) - ub + l, -1
        x(i) = x(i) - Row(i)%Col(j)*x(k)
        k = k-1
        write(*,*)
    end do
    x(i) = x(i)/Row(i)%Col(j)
    l = l - 1
end do

do i = ms - ub - 1, 1, -1
    x(i) = b(i)
    k =  size(Row(i)%Col)
    do j = size(Row(i)%Col), size(Row(i)%Col) - ub + 1, -1
        x(i) = x(i) - Row(i)%Col(j)*x(k)
        k = k-1
        write(*,*)
    end do
    x(i) = x(i)/Row(i)%Col(j)
end do

allocate(A(ms,ms))
allocate(y(ms))

A(1,1:) = [15, 7, 7, 0, 0, 0]
A(2,1:) = [2, 7, 1, 2, 0, 0]
A(3,1:) = [5, 1, 15, 2, 1, 0]
A(4,1:) = [3, 2, 3, 14, 1, 2]
A(5,1:) = [1, 3, 2, 4, 12, 3]
A(6,1:) = [0, 2, 3, 1, 5, 13]
b = [2, 1, 2, 4, 2, 6]

print *, "---- AFTER GE ----",char(10)

call pprint_band(Row,ms,lb,ub)
call pprint_array(b,size(b))
call pprint_array(x,size(x))
write(*,*) "---TEST---"
call gauss_elimination(A, ms, x, b)
call pprint_mat(A,ms,ms)

do i=1, size(Row)
    deallocate(Row(i)%Col)
end do
deallocate(Row)
deallocate(A)
deallocate(b)
deallocate(y)
deallocate(x)

end program main
