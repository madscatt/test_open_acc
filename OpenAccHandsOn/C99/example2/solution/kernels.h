#pragma once
#include <cmath>

#pragma acc routine seq
int IDX(int row, int col, int LDA) {
  return row*LDA+col;
}

void simpleLaplaceIter(int ROWS, int COLS, double * Ain, double * Aout)
{
  int lda=COLS+2;

  #pragma acc parallel loop present(Ain,Aout) gang worker
  for(int row=1;row<ROWS;row++) {
    #pragma acc loop vector
    for(int col=1;col<COLS;col++) {
      Aout[IDX(row,col,lda)]= .25 * 
        ( 
         + Ain[IDX(row  ,col-1,lda)] 
         + Ain[IDX(row  ,col+1,lda)] 
         + Ain[IDX(row-1,col  ,lda)] 
         + Ain[IDX(row+1,col  ,lda)] 
        );
    }
  }
}

double error(int ROWS, int COLS, double *A, double *B) {
  int lda=COLS+2;

  double sum=0.0;

  #pragma acc parallel loop reduction(+:sum) present(A,B) gang worker
  for(int row=1;row<ROWS;row++) {
    #pragma acc loop vector
    for(int col=1;col<COLS;col++) {
      double diff=A[IDX(row,col,lda)]-B[IDX(row,col,lda)];
      sum+=diff*diff;

    }
  }
  return sqrt(sum);
}
