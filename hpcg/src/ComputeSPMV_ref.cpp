
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file ComputeSPMV_ref.cpp

 HPCG routine
 */

#include "ComputeSPMV_ref.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

#include <algorithm>
#include <cstring>
#include <iostream>

/*!
  Routine to compute matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference SPMV implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV
*/
int ComputeSPMV_ref( const SparseMatrix & A, Vector & x, Vector & y) {

  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
#endif
  const double * const xv = x.values;
  double * const yv = y.values;
  const local_int_t nrow = A.localNumberOfRows;
#ifndef HPCG_NO_OPENMP
  #pragma omp parallel for
#endif
  for (local_int_t i=0; i< nrow; i++)  {
    double sum = 0.0;
    const double * const cur_vals = A.matrixValues[i];
    const local_int_t * const cur_inds = A.mtxIndL[i];
    const int cur_nnz = A.nonzerosInRow[i];

    for (int j=0; j< cur_nnz; j++)
      sum += cur_vals[j]*xv[cur_inds[j]];
    yv[i] = sum;
  }
  return 0;
}

int ComputeSPMV_mkl( const SparseMatrix & A, Vector  & x, Vector & y) {
  int nnz = A.totalNumberOfNonzeros;
  MKL_INT n_rows = A.localNumberOfRows;
  size_t a_off = 0;
  size_t j_off = 0;
  double *a = (double *)malloc(sizeof(double) * nnz);
  int *ia = (int *)malloc(sizeof(int) * (n_rows + 1));
  int *ja = (int *)malloc(sizeof(int) * nnz);
  ia[0] = 0;

  for(auto i = 0; i < n_rows; ++i) {
    auto len = A.nonzerosInRow[i];
    ia[i+1] = ia[i] + len;
    std::copy(A.matrixValues[i], A.matrixValues[i] + len, a + a_off);
    a_off += len;

    for(auto j = 0; j < len; ++j) {
      ja[j_off++] = A.mtxIndL[i][j];
    }
  }

  mkl_cspblas_dcsrgemv("N", &n_rows, a, ia, ja, x.values, y.values);
  free(a); free(ia); free(ja);
  return 0;
}
