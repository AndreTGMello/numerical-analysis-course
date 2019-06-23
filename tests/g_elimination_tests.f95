program bm_ge_tests
use bm_ge
use g_elimination
use utils

implicit none

!---------------
! Variables
!---------------

integer :: i,j,k,l,m,n,o,ub,lb,ms,pad
real(wp) :: multiplier
real(wp), allocatable :: b(:), x(:), y(:), A(:,:)
type(col), allocatable :: Row(:)

!---------------
! Logic
!---------------

! -------
! Test 1
! -------
! ub -> upper band
! lb -> lower band
ub = 2
lb = 4
! Square matrix size
ms = 6

allocate(b(ms))
allocate(x(ms))
allocate(Row(ms))

m = ub
o = -1
do i = 1, ms
  if ( i <= lb+1 ) then
    o = o + 1
  elseif ( i > ms-ub ) then
    m = m - 1
  end if
  allocate(Row(i)%col(m+o+1))
end do

Row(1)%Col = [15, 7, 7]
Row(2)%Col = [2, 7, 1, 2]
Row(3)%Col = [5, 1, 15, 2, 1]
Row(4)%Col = [3, 2, 3, 14, 1, 2]
Row(5)%Col = [1, 3, 2, 4, 12, 3]
Row(6)%Col = [2, 3, 1, 5, 13]
b = [2, 1, 2, 4, 2, 6]

call test_elimination(Row, b, ms, ub, lb)

deallocate(b)
deallocate(x)
m = ub
o = -1
do i = 1, ms
  if ( i <= lb+1 ) then
    o = o + 1
  elseif ( i > ms-ub ) then
    m = m - 1
  end if
  deallocate(Row(i)%col)
end do
deallocate(Row)

! -------
! Test 2
! -------
! ub -> upper band
! lb -> lower band
ub = 4
lb = 1
! Square matrix size
ms = 8

allocate(b(ms))
allocate(x(ms))
allocate(Row(ms))

m = ub
o = -1
do i = 1, ms
  if ( i <= lb+1 ) then
    o = o + 1
  elseif ( i > ms-ub ) then
    m = m - 1
  end if
  allocate(Row(i)%col(m+o+1))
end do

Row(1)%Col = [25, 7, 7, 1, 9]
Row(2)%Col = [2, 20, 1, 2, 4, 1]
Row(3)%Col = [5, 15, 2, 1, 3, 1]
Row(4)%Col = [3, 14, 1, 2, 2, 4]
Row(5)%Col = [4, 22, 3, 5, 1]
Row(6)%Col = [2, 33, 1, 5]
Row(7)%Col = [2, 11, 1]
Row(8)%Col = [2, 10]
b = [2, 2, 2, 4, 2, 6, 7, 10]

call test_elimination(Row, b, ms, ub, lb)

deallocate(b)
deallocate(x)
m = ub
o = -1
do i = 1, ms
  if ( i <= lb+1 ) then
    o = o + 1
  elseif ( i > ms-ub ) then
    m = m - 1
  end if
  deallocate(Row(i)%col)
end do
deallocate(Row)

! -------
! Test 3
! -------
! ub -> upper band
! lb -> lower band
ub = 3
lb = 3
! Square matrix size
ms = 9

allocate(b(ms))
allocate(x(ms))
allocate(Row(ms))

m = ub
o = -1
do i = 1, ms
  if ( i <= lb+1 ) then
    o = o + 1
  elseif ( i > ms-ub ) then
    m = m - 1
  end if
  allocate(Row(i)%col(m+o+1))
end do

Row(1)%Col = [0.18750, 0.10078,-0.04688,-0.00234]
Row(2)%Col = [0.10078, 0.24380,-0.03287,-0.05625,-0.00234]
Row(3)%Col = [-0.04688,-0.03287, 0.18750,-0.03522,-0.05625,-0.00234]
Row(4)%Col = [-0.00234,-0.05625,-0.03522, 0.18750,-0.03522,-0.05625,-0.00234]
Row(5)%Col = [-0.00234,-0.05625,-0.03522, 0.18750,-0.03522,-0.05625,-0.00234]
Row(6)%Col = [-0.00234,-0.05625,-0.03522, 0.18750,-0.03522,-0.05625,-0.00234]
Row(7)%Col = [-0.00234,-0.05625,-0.03522, 0.18750,-0.03287,-0.04688]
Row(8)%Col = [-0.00234,-0.05625,-0.03287, 0.24380, 0.10078]
Row(9)%Col = [-0.00234,-0.04688, 0.10078, 0.18750]

b = [-3.18305,-11.48596,-8.04394,6.70978,14.98292,6.70978,-8.04394,-11.48596,-3.18305]

call test_elimination(Row, b, ms, ub, lb)

deallocate(b)
deallocate(x)
do i = 1, ms
  if ( i <= lb+1 ) then
    o = o + 1
  elseif ( i > ms-ub ) then
    m = m - 1
  end if
  deallocate(Row(i)%col)
end do
deallocate(Row)


contains

	subroutine test_elimination(Row,b,ms,ub,lb)
		integer, intent(in) :: ub,lb,ms
		type(col), allocatable, intent(inout) :: Row(:)
		real(wp), allocatable, intent(inout) :: b(:)
		integer :: i,j,k,l,m,n,pad
		real(wp), allocatable :: c(:), x(:), y(:), A(:,:)

		allocate(x(size(b)))
		allocate(c(size(b)))
		c(1:) = b(1:)
		allocate(A(ms,ms))
		allocate(y(ms))
		pad = 0
		k = 0
		do i = 1, ms
			k = 0
			if ( i <= lb+1 ) then
				do j = 1, size(Row(i)%col)
					k = k + 1
					A(i,k) = Row(i)%Col(j)
				end do
				do j = 1, ms - size(Row(i)%col)
					k = k + 1
					A(i,k) =  0.0
				end do
			else
				pad = pad + 1
				do j = 1, pad
					k = k + 1
					A(i,k) =  0.0
				end do
				do j = 1, size(Row(i)%col)
					k = k + 1
					A(i,k) =  Row(i)%Col(j)
				end do
				do j = pad+size(Row(i)%col), ms-1
					k = k + 1
					A(i,k) =  0.0
				end do
			end if
		end do

		write(*,*) "--- TEST BEGIN ---"
		write(*,*)
		write(*,*) "--- BANDED Matrix ---"
		write(*,*) "--- Before Gaussian ---"
		call pprint_band(Row,lb,ub)
		write(*,*) "--- b ---"
		write(*,fmt="(F10.5, 3X)") b
		call gaussian_elimination_banded(Row,x,b,lb,ub)
		write(*,*)
		write(*,*) "--- After Gaussian ---"
		call pprint_band(Row,lb,ub)
		write(*,*) "--- b ---"
		write(*,fmt="(F10.5, 3X)") b
		write(*,*) "--- x ---"
		write(*,fmt="(F10.5, 3X)") x
		write(*,*)
		write(*,*) "--- REGULAR Matrix ---"
		write(*,*) "--- Before Gaussian ---"

		call pprint_mat(A)
		write(*,*) "--- b ---"
		write(*,fmt="(F10.5, 3X)") c
		call gaussian_elimination(A, y, c)
		write(*,*) "--- After Gaussian ---"
		call pprint_mat(A)
		write(*,*) "--- b ---"
		write(*,fmt="(F10.5, 3X)") c
		write(*,*) "--- x ---"
		write(*,fmt="(F10.5, 3X)") y
		write(*,*)
		write(*,*)
		write(*,*)
	end subroutine test_elimination

end program bm_ge_tests
