#include <cstdio>
#include <cstdlib>
using namespace std;

int main() {
  int N=1000000;
  int ITERS=1000;

  float *a,*b;
  a=(float*)malloc(N*sizeof(float));
  b=(float*)malloc(N*sizeof(float));

  for(int i=0;i<N;i++)
    a[i]=1.0;
  for(int i=0;i<N;i++)
    b[i]=0.5;

  #pragma acc kernels
  for(int iter=0;iter<ITERS;iter++) {
    for(int i=0;i<N;i++)
      b[i]=b[i]*b[i]+1.0f;
    for(int i=0;i<N;i++)
      a[i]=b[i]+a[i];
    for(int i=0;i<N;i++) {
      b[i]=b[i]/a[i];
      if (a[i]>10.0f) a[i]=a[i]-10.0;
    }
  }
  printf("%f %f %f %f\n",a[0],a[N-1],b[0],b[N-1]);
  return 0;
}
