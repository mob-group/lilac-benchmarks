#include <mkl.h>

#include <stdint.h>
#include <stdio.h>

struct {
  int32_t naa;
  int32_t nzz;
  int32_t firstrow;
  int32_t lastrow;
  int32_t firstcol;
  int32_t lastcol;
} common;

void c_spmv_(
    int *restrict rowstr,
    double *restrict a,
    double *restrict vec,
    int *restrict colidx,
    double *restrict outmat)
{
  int upper = common.lastrow - common.firstrow + 1;
  for(int j = 0; j < upper; ++j) {
    double sum = 0.0;
    for(int k = rowstr[j] - 1; k < rowstr[j+1] - 1; ++k) {
      sum += a[k] * vec[colidx[k] - 1];
    }
    outmat[j] = sum;
  }
}

static long long int *ia;
static long long int *ja;

void mkl_spmv_(
    int *restrict rowstr,     // ia
    double *restrict a,       // a
    double *restrict vec,     // x
    int *restrict colidx,     // ja
    double *restrict outmat)  // y
{
  MKL_INT n_rows = common.lastrow - common.firstrow + 1;

  mkl_dcsrgemv(
    "n", &n_rows,
    a, rowstr, colidx, 
    vec, outmat
  );
}
