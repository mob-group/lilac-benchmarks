#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static cublasHandle_t cublasH = NULL;
static cusparseHandle_t cusparseH = NULL;
static cudaStream_t stream = NULL;
static cusparseMatDescr_t descrA = NULL;

static cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
static cusparseStatus_t cusparseStat = CUSPARSE_STATUS_SUCCESS;

static cudaError_t cudaStat1 = cudaSuccess;
static cudaError_t cudaStat2 = cudaSuccess;
static cudaError_t cudaStat3 = cudaSuccess;
static cudaError_t cudaStat4 = cudaSuccess;
static cudaError_t cudaStat5 = cudaSuccess;

static double *d_csrValA = NULL;
static double *d_x = NULL;
static double *d_y = NULL;

static int *d_index_ptr = NULL;
static int *h_index_ptr = NULL;
static double *d_data_ptr = NULL;
static double *h_data_ptr = NULL;

static double h_one = 1.0;
static double h_zero = 0.0;

static int last_n = -1;
static int last_nnzA = -1;

void setup(int n, int nnzA)
{
  static bool ready = false;
  if(!ready) {
    // Set up CUDA libraries
    cudaStat1 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    assert(cudaSuccess == cudaStat1);

    cublasStat = cublasCreate(&cublasH);
    assert(CUBLAS_STATUS_SUCCESS == cublasStat);

    cublasStat = cublasSetStream(cublasH, stream);
    assert(CUBLAS_STATUS_SUCCESS == cublasStat);

    cusparseStat = cusparseCreate(&cusparseH);
    assert(CUSPARSE_STATUS_SUCCESS == cusparseStat);

    cusparseStat = cusparseSetStream(cusparseH, stream);
    assert(CUSPARSE_STATUS_SUCCESS == cusparseStat);

    // Set up matrix
    cusparseStat = cusparseCreateMatDescr(&descrA);
    assert(CUSPARSE_STATUS_SUCCESS == cusparseStat);

    cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ONE);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL );

    ready = true;
  }

  if(n != last_n || nnzA != last_nnzA) {
    cudaStat1 = cudaMalloc ((void**)&d_index_ptr , sizeof(int) * (n + 1 + nnzA));
    h_index_ptr = malloc(sizeof(int) * (n + 1 + nnzA));

    cudaStat2 = cudaMalloc ((void**)&d_data_ptr, sizeof(double) * (nnzA + n));
    h_data_ptr = malloc(sizeof(double) * (nnzA + n));

    cudaStat3 = cudaMalloc ((void**)&d_y, sizeof(double) * n);

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(h_index_ptr && h_data_ptr);

    last_n = n;
    last_nnzA = nnzA;
  }
}

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  int n = *rows + 1;
  int nnzA = rowstr[n] - rowstr[0];

  setup(n, nnzA);

  // Do device copy
  memcpy(h_index_ptr, rowstr, sizeof(int) * (n+1));
  memcpy(h_index_ptr + n+1, colidx, sizeof(int) * nnzA);
  cudaStat1 = cudaMemcpy(d_index_ptr, h_index_ptr, sizeof(int) * (n + 1 + nnzA), cudaMemcpyHostToDevice);

  memcpy(h_data_ptr, a, sizeof(double) * (nnzA));
  memcpy(h_data_ptr + nnzA, iv, sizeof(double) * n);
  cudaStat2 = cudaMemcpy(d_data_ptr, h_data_ptr, sizeof(double) * (nnzA + n), cudaMemcpyHostToDevice);

  assert(cudaSuccess == cudaStat1);
  assert(cudaSuccess == cudaStat2);
  
  // Do the SPMV
  cusparseStat = cusparseDcsrmv_mp(cusparseH,
                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                   n,
                                   n,
                                   nnzA,
                                   &h_one,
                                   descrA,
                                   d_data_ptr,
                                   d_index_ptr,
                                   d_index_ptr + (n+1),
                                   d_data_ptr + (nnzA),
                                   &h_zero,
                                   d_y);
  assert(CUSPARSE_STATUS_SUCCESS == cusparseStat);

  cudaStat1 = cudaMemcpy(ov, d_y, sizeof(double) * n, cudaMemcpyDeviceToHost);
  assert(cudaSuccess == cudaStat1);

  return NULL;
}
