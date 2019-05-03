module bm_ge
  ! Program for banded matrix gauss elimination
  use utils

  implicit none

  !---------------
  ! Variables
  !---------------
!  integer :: i,j

  !---------------
  ! Logic
  !---------------
  contains
    function solve_banded_matrix(A,b,m,n) result(x)
      integer, intent(in) :: m,n
      real, dimension(m,m), intent(in) :: A
      real, dimension(m), intent(in) :: b
      real, dimension(m) :: x
      ! m = rows of A = square matrix row/col numbers
      ! n = cols of A = number of diagonals
      x = b
    end function

end module bm_ge
