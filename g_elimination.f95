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
    integer :: i,j,k,l,m,n,ms,fix_col_index,fix_diag_col_index,temp_fix_col_index,diagonal,padding
    integer, intent(in) :: lb, ub
    real(wp), allocatable, intent(inout) :: x(:), b(:)
    type(col), allocatable, intent(inout) :: Row(:)

    ! Matrix Size:
    ms = size(Row)

    fix_col_index = 0
    fix_diag_col_index = 0
    ! Gaussian Elimination
    do i = 1, ms-1
        ! For each row:
        j = i + 1

        ! Fix pivot index
        if ( i > (lb+1) ) then
          fix_diag_col_index = fix_diag_col_index + 1
        end if

        do while ( j <= i+lb .and. j <= ms)

          ! Fix column index
          fix_col_index = 0
          if ( j > (lb+1) ) then
            fix_col_index = j - (lb+1)
          end if

          ! Multiplier
          Row(j)%Col(i-fix_col_index) = Row(j)%Col(i-fix_col_index)&
                        /Row(i)%Col(i-fix_diag_col_index)

          ! Apply multiplier to the systems solution (b)
          b(j) = b(j) - Row(j)%Col(i-fix_col_index)*b(i)

          ! For each column:
          k = i + 1
          do while ( k <= i+ub .and. k <= ms )
              Row(j)%Col(k-fix_col_index) = Row(j)%Col(k-fix_col_index)&
              - Row(i)%Col(k-fix_diag_col_index)*Row(j)%Col(i-fix_col_index)
              k = k + 1
          end do

          ! Set multiplier to zero
          Row(j)%Col(i-fix_col_index) = 0

          j = j + 1
        end do
!        call pprint_band(Row,lb,ub)
    end do


    ! Solving the system
    ! k, diagonal and padding are auxiliary variables
    ! for index fixing/shifting
    ! due to the way the band matrix is stored.

    padding = 0
    diagonal = size(Row(ms)%Col)
    do i = ms, 1, -1
      ! Fix index
      if ( ms-i > ub ) then
        padding = padding + 1
      end if
      if ( i <= lb ) then
        diagonal = diagonal - 1
      end if
      k = ms - padding

      ! Regular algorithm
      x(i) = b(i)
      do j = size(Row(i)%Col), diagonal + 1, -1
        x(i) = x(i) - Row(i)%Col(j)*x(k)
        k = k -1
      end do
      x(i) = x(i)/Row(i)%Col(j)
    end do

  end subroutine gaussian_elimination_banded

end module g_elimination
