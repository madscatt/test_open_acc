#include <cstdio>
#include <cstdlib>
using namespace std;

int main() {
  int N=100000000;
  int HN=10000;
  int maxv, minv;

  unsigned int *a, *h;
  a=(unsigned int*)malloc(N*sizeof(unsigned int));
  h=(unsigned int*)malloc(HN*sizeof(unsigned int));


  for(int i=0;i<N;i++)
    a[i]=lrand48()%HN;;

  for(int i=0;i<HN;i++)
    h[i]=0;

  for(int i=0;i<N;i++) {
    h[a[i]]++;
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
