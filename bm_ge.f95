program bm_ge
  ! Program for banded matrix gauss elimination
  use Utils

  implicit none

  real, dimension(3,3) :: matrix
  integer :: i,j

  do i=1,3
    do j=1,3
      matrix(i,j)=i+j
    end do
  end do

  pprint_mat(matrix)

end program bm_ge
