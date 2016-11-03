program example1
  integer, parameter :: N=1000000,&
                        ITERS=1000
  real, dimension(N) :: a, b
  integer            :: i, iter

  a = 1.0
  b = 0.5

  do iter=1,ITERS
    do i=1,N
      b(i) = b(i) * b(i) + 1.0
    enddo
    do i=1,N
      a(i) = b(i) + a(i)
    enddo
    do i=1,N
      b(i) = b(i) / a(i)
      if (a(i).gt.10.0) a(i) = a(i) - 10.0
    enddo
  enddo

  ! Print from arrays to avoid dead-code elimination
  print *, a(1),a(N),b(1),b(N)

end program example1
