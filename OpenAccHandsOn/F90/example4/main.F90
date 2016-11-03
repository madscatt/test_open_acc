program example4
  integer, parameter :: N=10000000, HN=10000
  integer            :: a(N), h(HN), i
  real               :: tmp(HN)

  do i = 1, N, HN
    call random_number(tmp)
    a(i:i+HN-1) = tmp * HN + 1
  end do

  h(:) = 0

  do i=1,N
    h(a(i)) = h(a(i)) + 1
  enddo

  print *, h(1:5),maxval(h),minval(h)

end program example4
