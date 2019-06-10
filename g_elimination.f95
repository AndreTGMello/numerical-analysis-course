module g_elimination
use utils

implicit none
contains
  subroutine gaussian_elimination(A, x, b)
    integer :: i, j, k, m
    real(wp), allocatable, intent(inout) :: A(:,:), x(:), b(:)
    m = size(A, dim=1)
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
  end subroutine

  subroutine gaussian_elimination_banded(Row, x, b, lb, ub)
    integer :: i,j,k,l,m,n,ms
    integer, intent(in) :: lb, ub
    real(wp), allocatable, intent(inout) :: x(:), b(:)
    type(col), allocatable, intent(inout) :: Row(:)

    ! Matrix Size:
    ms = size(Row)

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
            Row(j)%Col(i-m) = Row(j)%Col(i-m)/Row(i)%Col(i-l)
            b(j) = b(j) - Row(j)%Col(i-m)*b(i)
            ! For each column:
            k = i + 1
            do while ( k <= i+ub .and. k <= ms )
                Row(j)%Col(k-m) = Row(j)%Col(k-m) - Row(j)%Col(i-m)*Row(i)%Col(k-l)
                k = k + 1
            end do
            Row(j)%Col(i-m) = 0
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
  end subroutine gaussian_elimination_banded

end module g_elimination
