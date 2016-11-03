module kernels
  public :: simpleLaplaceIter,&
            error
  contains
  subroutine simpleLaplaceIter(nrows, ncols, ain, aout)
    implicit none
    integer :: nrows, ncols
    real(8), dimension(0:nrows+1,0:ncols+1),intent(in) :: ain
    real(8), dimension(0:nrows+1,0:ncols+1),intent(out) :: aout
    integer :: col,row

    do col=1,ncols
      do row=1,nrows
        aout(col,row) = 0.25 * (&
                       ain(col-1,row) + &
                       ain(col+1,row) + &
                       ain(col,row-1) + &
                       ain(col,row+1) )
      enddo
    enddo

  end subroutine simpleLaplaceIter
  function error(nrows, ncols, a, b)
    implicit none
    integer :: nrows, ncols
    real(8), dimension(0:nrows+1,0:ncols+1),intent(in)  :: a,b
    integer :: row, col
    real(8) :: error, diff, tot

    tot = 0
    do col=1,ncols
      do row=1,nrows
        diff = a(col,row) - b(col,row)
        tot = tot + diff*diff
      enddo
    enddo

    error = sqrt(tot)
  end function error
end module kernels
