#include <stdlib.h>
int random_crs_matprod__
(double *v,int *idx,int *nnz,int *ii,
 double *x,double *y,int *size)
{
  double sum;
  int m = *size,n = *size, shift = -1,i,j,jrow;

  x    = x + shift;    /* shift for Fortran start by 1 indexing */
  v    += shift; /* shift for Fortran start by 1 indexing */
  idx  += shift;
  for (i=0; i<m; i++) {
    jrow = ii[i];
    n    = ii[i+1] - jrow;
    sum  = 0.0;
    for (j=0; j<n; j++) {
      sum += v[jrow]*x[idx[jrow]]; jrow++;
     }
    y[i] = sum;
  }
  return 0;
}
int random_crs_matprod_t__
(double *v,int *idx,int *nnz,int *ii,
 double *x,double *y,int *size)
{
  return 0;
}
