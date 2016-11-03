#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;

#include "kernels.h"

int N=1000;
int MAX_ITERS=1000;
double TOL=2e-2;


int main() {
  double *in,*out;
  double tol;
  size_t size=(N+2)*(N+2);
  int k,iter=0;
  clock_t begin = clock();  
  in=(double*)malloc(size*sizeof(double));
  out=(double*)malloc(size*sizeof(double));

  for(int i=0;i<size;i++) {
    in[i]=0.0;
    out[i]=0.0;
  }

  k = (N+3)*(N/2)+1;
  in[k]     = 1000.0;
  in[k+1]   = 1000.0;
  in[k+N+2] = 1000.0;
  in[k+N+3] = 1000.0;
  do {
    simpleLaplaceIter(N,N,in,out);
    simpleLaplaceIter(N,N,out,in);
    tol=error(N,N,in,out);
    iter+=2;
    if(iter%100==0) {
      printf("Iters: %d, Tol: %lg\n", iter,tol);
    }
  } while (iter<MAX_ITERS && tol>TOL);

  printf("Final Iters: %d, Tol: %lg\n", iter,tol);
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  printf("TIME = %f\n",elapsed_secs) ;
  return 0;
}
