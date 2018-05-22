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

void setup(int rows, int cols, int nnzA)
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

  static int last_rows = -1;
  static int last_cols = -1;
  static int last_nnzA = -1;
  if(cols != last_cols || rows != last_rows || nnzA != last_nnzA) {
    cudaStat1 = cudaMalloc ((void**)&d_csrRowPtrA, sizeof(int) * (rows+1) );
    cudaStat2 = cudaMalloc ((void**)&d_csrColIndA, sizeof(int) * nnzA );
    cudaStat3 = cudaMalloc ((void**)&d_csrValA   , sizeof(double) * nnzA );
    cudaStat4 = cudaMalloc ((void**)&d_x         , sizeof(double) * cols );
    cudaStat5 = cudaMalloc ((void**)&d_y         , sizeof(double) * rows );

    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);
    assert(cudaSuccess == cudaStat5);

    last_rows = rows;
    last_cols = cols;
    last_nnzA = nnzA;
  }
}

struct sigaction old_sigaction;

static double *a_begin = NULL;
static double *a_end = NULL;
static void *aligned_a_begin = NULL;
static void *aligned_a_end = NULL;

static int *rowstr_begin = NULL;
static int *rowstr_end = NULL;
static void *aligned_rowstr_begin = NULL;
static void *aligned_rowstr_end = NULL;

static int *colidx_begin = NULL;
static int *colidx_end = NULL;
static void *aligned_colidx_begin = NULL;
static void *aligned_colidx_end = NULL;

static void handler(int sig, siginfo_t *si, void *unused)
{
  void *addr = si->si_addr;
  if((addr >= aligned_a_begin && addr < aligned_a_end) ||
     (addr >= aligned_rowstr_begin && addr < aligned_rowstr_end) ||
     (addr >= aligned_colidx_begin && addr < aligned_colidx_end)) 
  {
    mprotect(aligned_a_begin, aligned_a_end - aligned_a_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
    mprotect(aligned_rowstr_begin, aligned_rowstr_end - aligned_rowstr_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
    mprotect(aligned_colidx_begin, aligned_colidx_end - aligned_colidx_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
    
    a_begin = NULL;
    rowstr_begin = NULL;
    colidx_begin = NULL;
  } else if(old_sigaction.sa_sigaction) {
    old_sigaction.sa_sigaction(sig, si, unused);
  }
}

#define ALIGN(name) \
  aligned_##name##_begin = name##_begin; \
  aligned_##name##_begin -= (size_t)aligned_##name##_begin % sysconf(_SC_PAGE_SIZE); \
  aligned_##name##_end = name##_end + sysconf(_SC_PAGE_SIZE) - 1; \
  aligned_##name##_end -= (size_t)aligned_##name##_end % sysconf(_SC_PAGE_SIZE); \
  mprotect(aligned_##name##_begin, aligned_##name##_end - aligned_##name##_begin, PROT_READ | PROT_EXEC);

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  static int cols = 0;

  int n = *rows + 1;
  int nnzA = rowstr[n] - rowstr[0];
  if(cols == 0) {
    for(int i = rowstr[0]; i < rowstr[n]; ++i) {
      if(colidx[i] >= cols) {
        cols = colidx[i];
      }
    }
  }

  setup(n, cols, nnzA);

  if((a_begin != NULL && a_begin != a) ||
     (rowstr_begin != NULL && rowstr_begin != rowstr) ||
     (colidx_begin != NULL && colidx_begin != colidx))
  {
    a_begin = NULL;
    rowstr_begin = NULL;
    colidx_begin = NULL;
  }

  if(a_begin == NULL || rowstr_begin == NULL || colidx_begin == NULL) {
    // Do device copy
    cudaStat1 = cudaMemcpy(d_csrRowPtrA, rowstr, sizeof(int) * (n+1), cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(d_csrColIndA, colidx, sizeof(int) * nnzA, cudaMemcpyHostToDevice);
    cudaStat3 = cudaMemcpy(d_csrValA, a, sizeof(double) * nnzA, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);

    struct sigaction sa;
    sa.sa_flags = SA_SIGINFO;
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = handler;
    sigaction(SIGSEGV, &sa, &old_sigaction);

    a_begin = a;
    a_end = a_begin + nnzA;
    ALIGN(a);

    rowstr_begin = rowstr;
    rowstr_end = rowstr_begin + n + 1;
    ALIGN(rowstr);

    colidx_begin = colidx;
    colidx_end = colidx_begin + nnzA;
    ALIGN(colidx);
  }

  cudaStat1 = cudaMemcpy(d_x, iv, sizeof(double) * cols, cudaMemcpyHostToDevice);
  assert(cudaSuccess == cudaStat1);

  // Do the SPMV
  double one = 1.0;
  double zero = 0.0;
  cusparseStat = cusparseDcsrmv_mp(cusparseH,
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
