#include <mkl.h>

#include <stdio.h>

void c_spmv_(
  double *out,
  double *vec,
  double *A, int *IA, int *JA,
  int *n_min, int *non_zero,
  int *n_col
)
{
  for(int J3 = 0; J3 < *n_min; ++J3) {
    if(n_col[J3] == 0) { continue; }
    double sum = 0.0;

    for(int J2 = IA[J3] - 1; J2 < IA[J3+1]-1; ++J2) {
      sum += vec[J2] * A[JA[J2] - 1];
    }

    out[J3] = sum;
  }
}

void mkl_spmv_(
  double *out,
  double *vec,
  double *A, int *IA, int *JA,
  int *n_min, int *non_zero,
  int *n_col
)
{
  mkl_dcsrgemv("n", n_min, A, IA, JA, vec, out);
}
