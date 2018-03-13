#include <mkl.h>

#include <stdlib.h>
#include <string.h>

double starttimer_(void);
double stoptimer_(void);
void add_mult_time_(double*);

int random_crs_matprod_ (
    double *v,int *idx,int *nnz,int *ii,
    double *x,double *y,int *size)
{
  char *version = getenv("SPMV_VERSION");
  double start = starttimer_();

  if(version && strcmp(version, "MKL") == 0) {
    mkl_dcsrgemv("N", size, v, ii, idx, x, y);
  } else {
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
  }

  double stop = stoptimer_() - start;
  add_mult_time_(&stop);
  return 0;
}
int random_crs_matprod_t_
(double *v,int *idx,int *nnz,int *ii,
 double *x,double *y,int *size)
{
  return random_crs_matprod_(v, idx, nnz, ii, x, y, size);
}
