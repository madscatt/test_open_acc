#include <cstdio>
#include <cstdlib>
#include "curand.h"

using namespace std;

int main() {
  int N=100000000;
  int HN=10000;
  int maxv, minv, istat;
  curandGenerator_t g;
  
  int *a, *h;
  a=(int*)malloc(N*sizeof(int));
  h=(int*)malloc(HN*sizeof(int));

  #pragma acc data create(a[0:N]) copyout(h[0:HN])
  {
    istat = curandCreateGenerator(&g, CURAND_RNG_PSEUDO_DEFAULT);
    if (istat != CURAND_STATUS_SUCCESS) printf("Error in curand\n");
    #pragma acc host_data use_device(a)
    {
    istat = curandGenerate(g, (unsigned int*)a, N);
    if (istat != CURAND_STATUS_SUCCESS) printf("Error in curand\n");
    }
    istat = curandDestroyGenerator(g);
    if (istat != CURAND_STATUS_SUCCESS) printf("Error in curand\n");

    #pragma acc parallel loop
    for(int i=0;i<N;i++)
      a[i]=(0x7fffffff & a[i])%HN;

    #pragma acc parallel loop
    for(int i=0;i<HN;i++) 
      h[i]=0;
    
    #pragma acc parallel loop
    for(int i=0;i<N;i++) {
      #pragma acc atomic
      h[a[i]]+=1;
    }
  }

  maxv = 0;
  minv = N;
  for(int i=0;i<HN;i++) {
    if (h[i]>maxv) maxv = h[i];
    if (h[i]<minv) minv = h[i];
  }
  printf("%i %i %i %i %i %i %i %i\n",h[0],h[1],h[2],h[3],h[4],maxv,minv,N/HN);
  return 0;
}
