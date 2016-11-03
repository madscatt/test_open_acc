#pragma once
#include<cmath>
#include "matrix.h"

struct vector {
  unsigned int n;
  double *coefs;
};

void allocate_vector(vector &v, unsigned int n) {
  v.n=n;
  v.coefs=(double*)malloc(n*sizeof(double));
  #pragma acc enter data create(v.coefs[0:n]) 
}

void free_vector(vector &v) {
  double *vcoefs=v.coefs;
  #pragma acc exit data delete(vcoefs)
  free(v.coefs);

}

void initialize_vector(vector &v,double val) {

  for(int i=0;i<v.n;i++)
    v.coefs[i]=val;

  #pragma acc update device(v.coefs[0:v.n])
  
}

