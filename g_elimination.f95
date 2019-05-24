module g_elimination
use utils

implicit none
contains
  subroutine gauss_elimination(A, m, x, b)
    integer :: i, j, k
    integer, intent(in) :: m
    real, intent(inout) :: A(m,m), x(m), b(m)
    ! Gauss elimination (on regular matrix)
    do i = 1, m-1
        do j = i+1, m
            A(j,i) = A(j,i)/A(i,i)
            b(j) = b(j) - b(i)*A(j,i)
            do k = i+1, m
                A(j,k) = A(j,k) - A(i,k)*A(j,i)
            end do
            A(j,i) = 0
        end do
    end do
    ! Solving the system (on regular matrix)
    do i = m, 1, -1
        x(i) = b(i)
        do j = m, i+1, -1
            x(i) = x(i) - A(i,j)*x(j)
        end do
        x(i) = x(i)/A(i,i)
    end do
    call pprint_array(x, size(x))
  end subroutine

end module g_elimination
