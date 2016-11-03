#include <cstdio>
#include <cstdlib>
using namespace std;

#include "kernels.h"

int N=1000;
int MAX_ITERS=1000;
double TOL=1e-9;


int main() {
  double *in,*out;
  double tol;
  size_t size=(N+2)*(N+2);
  int iter=0;
  
  in=(double*)malloc(size*sizeof(double));
  out=(double*)malloc(size*sizeof(double));

  for(int i=0;i<N;i++)
    in[i]=drand48();
  for(int i=0;i<N;i++)
    out[i]=in[i];

  do {
    simpleLaplaceIter(N,N,in,out);
    simpleLaplaceIter(N,N,out,in);
    tol=error(N,N,in,out);
    iter+=2;
    if(iter%100==0) {
      printf("Iters: %d, Tol: %lg\n", iter,tol);
    }
  } while (iter<MAX_ITERS && tol>TOL);

  printf("Iters: %d, Tol: %lg\n", iter,tol);
  return 0;
}
