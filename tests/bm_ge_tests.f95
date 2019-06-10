program bm_ge_tests
use bm_ge
use g_elimination
use utils

implicit none

!---------------
! Variables
!---------------

integer :: i,j,k,l,m,n,ub,lb,ms
real(wp) :: multiplier
real(wp), allocatable :: b(:), x(:), y(:), A(:,:)
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

call pprint_band(Row,lb,ub)
call pprint_array(b)
call gaussian_elimination_banded(Row,x,b,lb,ub)

print *, "---- AFTER GE ----",char(10)

call pprint_band(Row,lb,ub)
call pprint_array(b)
call pprint_array(x)
write(*,*) "---TEST---"
allocate(A(ms,ms))
allocate(y(ms))
A(1,1:) = [15, 7, 7, 0, 0, 0]
A(2,1:) = [2, 7, 1, 2, 0, 0]
A(3,1:) = [5, 1, 15, 2, 1, 0]
A(4,1:) = [3, 2, 3, 14, 1, 2]
A(5,1:) = [1, 3, 2, 4, 12, 3]
A(6,1:) = [0, 2, 3, 1, 5, 13]
b = [2, 1, 2, 4, 2, 6]
call gaussian_elimination(A, y, b)
call pprint_array(b)
call pprint_array(x)
call pprint_mat(A)

do i=1, size(Row)
    deallocate(Row(i)%Col)
end do
deallocate(Row)
deallocate(A)
deallocate(b)
deallocate(y)
deallocate(x)

end program bm_ge_tests
