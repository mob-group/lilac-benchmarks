#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <unistd.h>
#include <signal.h>
#include <sys/mman.h>

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

static int *d_csrRowPtrA = NULL;
static int *d_csrColIndA = NULL;
static double *d_csrValA = NULL;
static double *d_x = NULL;
static double *d_y = NULL;

void setup(void)
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
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);

    ready = true;
  }
}

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  setup();

  static double *last_a = NULL;
  static int *last_rowstr = NULL;
  static int *last_colidx = NULL;
  static int last_rows = -1;
  static double *last_iv = NULL;

  static int cols = 0;

  int n = *rows;
  int nnzA = rowstr[n] - rowstr[0];
  if(cols == 0) {
    for(int i = rowstr[0]; i < rowstr[n]; ++i) {
      if(colidx[i] >= cols) {
        cols = colidx[i];
      }
    }
  }

  if( a != last_a ||
      rowstr != last_rowstr ||
      colidx != last_colidx ||
      *rows != last_rows ||
      iv != last_iv)
  {
    #define free_if_needed(mem) do { if(mem) { cudaFree(mem); } } while(0);
    free_if_needed(d_csrRowPtrA);
    free_if_needed(d_csrColIndA);
    free_if_needed(d_csrValA);
    free_if_needed(d_x);
    free_if_needed(d_y);
    #undef free_if_needed

    cudaStat1 = cudaMalloc ((void**)&d_csrRowPtrA, sizeof(int) * (n+1) );
    cudaStat2 = cudaMalloc ((void**)&d_csrColIndA, sizeof(int) * nnzA );
    cudaStat3 = cudaMalloc ((void**)&d_csrValA   , sizeof(double) * nnzA );
    cudaStat4 = cudaMalloc ((void**)&d_x         , sizeof(double) * cols );
    cudaStat5 = cudaMalloc ((void**)&d_y         , sizeof(double) * n );

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);
    assert(cudaSuccess == cudaStat5);

    cudaStat1 = cudaMemcpy(d_csrRowPtrA, rowstr, sizeof(int) * (n+1), cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(d_csrColIndA, colidx, sizeof(int) * nnzA, cudaMemcpyHostToDevice);
    cudaStat3 = cudaMemcpy(d_csrValA, a, sizeof(double) * nnzA, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);

    cudaStat1 = cudaMemcpy(d_x, iv, sizeof(double) * cols, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cudaStat1);

    last_a = a;
    last_rowstr = rowstr;
    last_colidx = colidx;
    last_rows = *rows;
    last_iv = iv;
  }

  double one = 1.0;
  double zero = 0.0;
  cusparseStat = cusparseDcsrmv(cusparseH,
                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                   n,
                                   cols,
                                   nnzA,
                                   &one,
                                   descrA,
                                   d_csrValA,
                                   d_csrRowPtrA,
                                   d_csrColIndA,
                                   d_x,
                                   &zero,
                                   d_y);
  assert(CUSPARSE_STATUS_SUCCESS == cusparseStat);

  cudaStat1 = cudaMemcpy(ov, d_y, sizeof(double) * n, cudaMemcpyDeviceToHost);
  assert(cudaSuccess == cudaStat1);

  return NULL;
}
