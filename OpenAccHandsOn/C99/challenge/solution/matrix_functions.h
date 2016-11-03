#pragma once
#include "vector.h"
#include "matrix.h"

void matvec(const matrix& A, const vector& x, const vector &y) {

  unsigned int num_rows=A.num_rows;
  unsigned int *row_offsets=A.row_offsets;
  unsigned int *cols=A.cols;
  double *Acoefs=A.coefs;
  double *xcoefs=x.coefs;
  double *ycoefs=y.coefs;

  #pragma acc parallel loop gang worker num_workers(4) vector_length(32) present(row_offsets,cols,Acoefs,xcoefs,ycoefs) async(1)
  for(int i=0;i<num_rows;i++) {
    double sum=0;
    int row_start=row_offsets[i];
    int row_end=row_offsets[i+1];
    #pragma acc loop vector reduction(+:sum)
    for(int j=row_start;j<row_end;j++) {
      unsigned int Acol=cols[j];
      double Acoef=Acoefs[j];
      double xcoef=xcoefs[Acol];
      sum+=Acoef*xcoef;
    }
    ycoefs[i]=sum;
  }
}
