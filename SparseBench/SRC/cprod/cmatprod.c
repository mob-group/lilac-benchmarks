#include <mkl.h>

#include <stdlib.h>
#include <string.h>

double starttimer_(void);
double stoptimer_(void);
void add_mult_time_(double*);
void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows);

int random_crs_matprod_ (
    double *v,int *idx,int *nnz,int *ii,
    double *x,double *y,int *size)
{
  double start = starttimer_();

  spmv_harness_(y, v, x, ii, idx, size);

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
