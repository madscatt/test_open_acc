program example4
  use curand
  integer, parameter :: N=10000000, HN=10000
  integer            :: a(N), h(HN), i
  type(curandGenerator) :: g

  istat = curandCreateGenerator(g,CURAND_RNG_PSEUDO_XORWOW)

  !$acc data create(a) copyout(h)

  !$acc host_data use_device(a)
  istat = curandGenerate(g, a, N)
  !$acc end host_data
  
  !$acc kernels
  a = mod(abs(a),HN) + 1
  !$acc end kernels

  !$acc kernels
  h(:) = 0
  !$acc end kernels

  !$acc parallel loop
  do i=1,N
    !$acc atomic 
    h(a(i)) = h(a(i)) + 1
  enddo
  !$acc end parallel loop
  !$acc end data

  print *, h(1:5),maxval(h),minval(h)

end program example4
