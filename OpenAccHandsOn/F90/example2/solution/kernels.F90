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

    !$acc parallel loop present(ain,aout) gang worker
    do col=1,ncols
      !$acc loop vector
      do row=1,nrows
        aout(col,row) = 0.25 * (&
                       ain(col-1,row) + &
                       ain(col+1,row) + &
                       ain(col,row-1) + &
                       ain(col,row+1) )
      enddo
    enddo
    !$acc end parallel loop

  end subroutine simpleLaplaceIter
  function error(nrows, ncols, a, b)
    implicit none
    integer :: nrows, ncols
    real(8), dimension(0:nrows+1,0:ncols+1),intent(in)  :: a,b
    integer :: row, col
    real(8) :: error, diff, tot

    tot = 0
    !$acc parallel loop reduction(+:tot) present(a,b) gang worker
    do col=1,ncols
      !$acc loop vector
      do row=1,nrows
        diff = a(col,row) - b(col,row)
        tot = tot + diff*diff
      enddo
    enddo
    !$acc end parallel loop

    error = sqrt(tot)
  end function error
end module kernels
