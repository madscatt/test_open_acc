program example2
  use kernels
  integer, parameter              :: N=1000, ITERS=1000
  real(8), parameter              :: MAXTOL= 2e-2
  real(8), dimension(0:N+1,0:N+1) :: inarr, outarr
  integer                         :: i, iter=0
  real(8)                         :: err

  inarr  = 0.0
  outarr = 0.0
  inarr(N/2:N/2+1,N/2:N/2+1) = real(N)
  err = MAXTOL+1

  do while ((err.gt.Maxtol).and.(iter.lt.ITERS))
    call simpleLaplaceIter(N,N,inarr,outarr)
    call simpleLaplaceIter(N,N,outarr,inarr)
    err = error(N,N,inarr,outarr)
    iter = iter + 2
    if(mod(iter,100).eq.0) then
      print *, "Iters:",iter,"Error:",err
    endif
  enddo

  print *, "Final Iters:",iter,"Error:",err

end program example2
