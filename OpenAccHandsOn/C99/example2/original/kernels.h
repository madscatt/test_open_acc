#pragma once
#include <cmath>

int IDX(int row, int col, int LDA) {
  return row*LDA+col;
}

void simpleLaplaceIter(int ROWS, int COLS, double * Ain, double * Aout)
{
  int lda=COLS+2;

  for(int col=1;col<COLS;col++) {
    for(int row=1;row<ROWS;row++) {
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

  for(int col=1;col<COLS;col++) {
    for(int row=1;row<ROWS;row++) {
      double diff=A[IDX(row,col,lda)]-B[IDX(row,col,lda)];
      sum+=diff*diff;

    }
  }
  return sqrt(sum);
}
