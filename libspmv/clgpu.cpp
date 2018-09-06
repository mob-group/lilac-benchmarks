#define CL_HPP_ENABLE_EXCEPTION
#define CL_HPP_MINIMUM_OPENCL_VERSION BUILD_CLVERSION
#define CL_HPP_TARGET_OPENCL_VERSION BUILD_CLVERSION

#include "clgpu-model.h"
#include "native-impl.h"

#include <CL/cl.hpp>
#include <clSPARSE.h>
#include <clSPARSE-error.h>

#include <algorithm>
#include <csignal>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <sys/mman.h>

namespace {
  cl_int cl_status;
  clsparseStatus status;

  std::vector<cl::Device> devices{};
  std::vector<cl::Platform> platforms{};

  cl::Device* device;
  cl::Platform* platform;

  cl::Context context;
  cl::CommandQueue queue;
  clsparseControl control;

  clsparseScalar alpha;
  clsparseScalar beta;

  cldenseVector x;
  cldenseVector y;
  clsparseCsrMatrix A;
}

void set_platform_device(int p, int d)
{
  cl_status = cl::Platform::get(&platforms);
  if(cl_status != CL_SUCCESS) {
    std::cout << "Problem getting OpenCL platforms"
              << " [" << cl_status << "]" << '\n';
    std::exit(2);
  }

  platform = &platforms[p];

  cl_status = platform->getDevices(CL_DEVICE_TYPE_GPU, &devices);
  if(cl_status != CL_SUCCESS) {
    std::cout << "Problem getting devices from platform"
              << platform->getInfo<CL_PLATFORM_NAME>()
              << " error: [" << cl_status << "]" << '\n';
  }

  device = &devices[d];
}

void init_alpha_beta()
{
  clsparseInitScalar(&alpha);
  alpha.value = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                               sizeof(double), nullptr, 
                               &cl_status);

  clsparseInitScalar(&beta);
  beta.value = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              sizeof(double), nullptr, 
                              &cl_status);

  auto halpha = static_cast<double *>(clEnqueueMapBuffer(
      queue(), alpha.value, CL_TRUE, CL_MAP_WRITE,
      0, sizeof(double), 0, nullptr, nullptr, &cl_status));
  *halpha = 1.0f;
  cl_status = clEnqueueUnmapMemObject(queue(), alpha.value, halpha,
                                      0, nullptr, nullptr);

  auto hbeta = static_cast<double *>(clEnqueueMapBuffer(
        queue(), beta.value, CL_TRUE, CL_MAP_WRITE,
        0, sizeof(double), 0, nullptr, nullptr, &cl_status));
  *hbeta = 0.0f;
  cl_status = clEnqueueUnmapMemObject(queue(), beta.value, hbeta,
                                      0, nullptr, nullptr);
}

void f_init_alpha_beta()
{
  clsparseInitScalar(&alpha);
  alpha.value = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                               sizeof(float), nullptr, 
                               &cl_status);

  clsparseInitScalar(&beta);
  beta.value = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              sizeof(float), nullptr, 
                              &cl_status);

  auto halpha = static_cast<float *>(clEnqueueMapBuffer(
      queue(), alpha.value, CL_TRUE, CL_MAP_WRITE,
      0, sizeof(float), 0, nullptr, nullptr, &cl_status));
  *halpha = 1.0f;
  cl_status = clEnqueueUnmapMemObject(queue(), alpha.value, halpha,
                                      0, nullptr, nullptr);

  auto hbeta = static_cast<float *>(clEnqueueMapBuffer(
        queue(), beta.value, CL_TRUE, CL_MAP_WRITE,
        0, sizeof(float), 0, nullptr, nullptr, &cl_status));
  *hbeta = 0.0f;
  cl_status = clEnqueueUnmapMemObject(queue(), beta.value, hbeta,
                                      0, nullptr, nullptr);
}

void setup(int rows, int cols, int nnz)
{
  static bool ready = false;
  if(!ready) {
    set_platform_device(0, 0);

    context = cl::Context(*device);
    queue = cl::CommandQueue(context, *device);
    
    // Setup clsparse
    status = clsparseSetup();
    if(status != clsparseSuccess) {
      std::cout << "Problem setting up clSPARSE\n";
      std::exit(3);
    }

    auto createResult = clsparseCreateControl(queue());
    CLSPARSE_V(createResult.status, "Failed to create status control");
    control = createResult.control;

    init_alpha_beta();

    clsparseInitVector(&x);
    clsparseInitVector(&y);
    clsparseInitCsrMatrix(&A);

    // Setup GPU buffers
    ready = true;
  }

  static int last_rows = -1;
  static int last_cols = -1;
  static int last_nnzA = -1;
  if(cols != last_cols || rows != last_rows || nnz != last_nnzA) {
    A.num_rows = rows;
    A.num_cols = cols;
    A.num_nonzeros = nnz;

    x.num_values = A.num_cols;
    y.num_values = A.num_rows;

    if(A.values) { clReleaseMemObject(A.values); }
    if(A.col_indices) { clReleaseMemObject(A.col_indices); }
    if(A.row_pointer) { clReleaseMemObject(A.row_pointer); }
    if(x.values) { clReleaseMemObject(x.values); }
    if(y.values) { clReleaseMemObject(y.values); }

    A.values = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              sizeof(double) * A.num_nonzeros, 
                              nullptr, &cl_status);

    A.col_indices = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                                   sizeof(int) * A.num_nonzeros, 
                                   nullptr, &cl_status);

    A.row_pointer = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                                   sizeof(int) * (A.num_rows + 1), 
                                   nullptr, &cl_status);

    x.values = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              x.num_values * sizeof(double),
                              nullptr, &cl_status);

    y.values = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              y.num_values * sizeof(double),
                              nullptr, &cl_status);

    last_cols = cols;
    last_rows = rows;
    last_nnzA = nnz;
  }
}

void f_setup(int rows, int cols, int nnz)
{
  static bool ready = false;
  if(!ready) {
    set_platform_device(0, 0);

    context = cl::Context(*device);
    queue = cl::CommandQueue(context, *device);
    
    // Setup clsparse
    status = clsparseSetup();
    if(status != clsparseSuccess) {
      std::cout << "Problem setting up clSPARSE\n";
      std::exit(3);
    }

    auto createResult = clsparseCreateControl(queue());
    CLSPARSE_V(createResult.status, "Failed to create status control");
    control = createResult.control;

    f_init_alpha_beta();

    clsparseInitVector(&x);
    clsparseInitVector(&y);
    clsparseInitCsrMatrix(&A);

    // Setup GPU buffers
    ready = true;
  }

  static int last_rows = -1;
  static int last_cols = -1;
  static int last_nnzA = -1;
  if(cols != last_cols || rows != last_rows || nnz != last_nnzA) {
    A.num_rows = rows;
    A.num_cols = cols;
    A.num_nonzeros = nnz;

    x.num_values = A.num_cols;
    y.num_values = A.num_rows;

    if(A.values) { clReleaseMemObject(A.values); }
    if(A.col_indices) { clReleaseMemObject(A.col_indices); }
    if(A.row_pointer) { clReleaseMemObject(A.row_pointer); }
    if(x.values) { clReleaseMemObject(x.values); }
    if(y.values) { clReleaseMemObject(y.values); }

    A.values = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              sizeof(float) * A.num_nonzeros, 
                              nullptr, &cl_status);

    A.col_indices = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                                   sizeof(int) * A.num_nonzeros, 
                                   nullptr, &cl_status);

    A.row_pointer = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                                   sizeof(int) * (A.num_rows + 1), 
                                   nullptr, &cl_status);

    x.values = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              x.num_values * sizeof(float),
                              nullptr, &cl_status);

    y.values = clCreateBuffer(context(), CL_MEM_READ_ONLY, 
                              y.num_values * sizeof(float),
                              nullptr, &cl_status);

    last_cols = cols;
    last_rows = rows;
    last_nnzA = nnz;
  }
}

struct sigaction old_sigaction;

static double *a_begin = NULL;
static double *a_end = NULL;
static double *aligned_a_begin = NULL;
static double *aligned_a_end = NULL;

static float *af_begin = NULL;
static float *af_end = NULL;
static float *aligned_af_begin = NULL;
static float *aligned_af_end = NULL;

static int *rowstr_begin = NULL;
static int *rowstr_end = NULL;
static int *aligned_rowstr_begin = NULL;
static int *aligned_rowstr_end = NULL;

static int *colidx_begin = NULL;
static int *colidx_end = NULL;
static int *aligned_colidx_begin = NULL;
static int *aligned_colidx_end = NULL;

static void handler(int sig, siginfo_t *si, void *unused)
{
  void *addr = si->si_addr;
  if((addr >= aligned_a_begin && addr < aligned_a_end) ||
     (addr >= aligned_rowstr_begin && addr < aligned_rowstr_end) ||
     (addr >= aligned_colidx_begin && addr < aligned_colidx_end)) 
  {
    mprotect(aligned_a_begin, aligned_a_end - aligned_a_begin, 
             PROT_READ | PROT_WRITE | PROT_EXEC);
    mprotect(aligned_rowstr_begin, aligned_rowstr_end - aligned_rowstr_begin, 
             PROT_READ | PROT_WRITE | PROT_EXEC);
    mprotect(aligned_colidx_begin, aligned_colidx_end - aligned_colidx_begin, 
             PROT_READ | PROT_WRITE | PROT_EXEC);
    
    a_begin = NULL;
    rowstr_begin = NULL;
    colidx_begin = NULL;
  } else if(old_sigaction.sa_sigaction) {
    old_sigaction.sa_sigaction(sig, si, unused);
  }
}

static void f_handler(int sig, siginfo_t *si, void *unused)
{
  void *addr = si->si_addr;
  if((addr >= aligned_af_begin && addr < aligned_af_end) ||
     (addr >= aligned_rowstr_begin && addr < aligned_rowstr_end) ||
     (addr >= aligned_colidx_begin && addr < aligned_colidx_end)) 
  {
    mprotect(aligned_af_begin, aligned_af_end - aligned_af_begin, 
             PROT_READ | PROT_WRITE | PROT_EXEC);
    mprotect(aligned_rowstr_begin, aligned_rowstr_end - aligned_rowstr_begin, 
             PROT_READ | PROT_WRITE | PROT_EXEC);
    mprotect(aligned_colidx_begin, aligned_colidx_end - aligned_colidx_begin, 
             PROT_READ | PROT_WRITE | PROT_EXEC);
    
    af_begin = NULL;
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

extern "C" {

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  int n_rows = *rows;
  int nnzA = rowstr[n_rows] - rowstr[0];

  auto impl = predict(n_rows, nnzA);
  if(impl == CL_NATIVE_IMPL) {
    native_spmv(ov, a, iv, rowstr, colidx, rows);
  } else if(impl == CL_GPU_IMPL) {
    int n_cols = 0;
    for(int i = rowstr[0] - 1; i < rowstr[n_rows] - 1; i++) {
        if(colidx[i] >= n_cols) {
          n_cols = colidx[i];
        }
    }

    setup(n_rows, n_cols, nnzA);

    if((a_begin != nullptr && a_begin != a) ||
       (rowstr_begin != nullptr && rowstr_begin != rowstr) ||
       (colidx_begin != nullptr && colidx_begin != colidx))
    {
      a_begin = nullptr;
      rowstr_begin = nullptr;
      colidx_begin = nullptr;
    }

    if(a_begin == nullptr || rowstr_begin == nullptr || colidx_begin == nullptr) {
      cl_status = clEnqueueWriteBuffer(queue(), A.values, true, 0, 
                                       sizeof(double) * A.num_nonzeros, a, 
                                       0, nullptr, nullptr);

      auto one_based = new int[A.num_nonzeros + A.num_rows + 1];
      const auto sub_one = [](int v) { return v - 1; };

      auto it = std::transform(colidx, colidx + A.num_nonzeros, one_based, sub_one);
      std::transform(rowstr, rowstr + A.num_rows + 1, it, sub_one);

      cl_status = clEnqueueWriteBuffer(queue(), A.col_indices, true, 0,
                                       sizeof(int) * A.num_nonzeros, one_based,
                                       0, nullptr, nullptr);

      cl_status = clEnqueueWriteBuffer(queue(), A.row_pointer, true, 0,
                                       sizeof(int) * (A.num_rows + 1), 
                                       one_based + A.num_nonzeros,
                                       0, nullptr, nullptr);

      clsparseCsrMetaCreate(&A, control);

      delete[] one_based;

      cl_status = clEnqueueWriteBuffer(queue(), y.values, true, 0,
                                       sizeof(double) * y.num_values, ov,
                                       0, nullptr, nullptr);

      struct sigaction sa;
      sa.sa_flags = SA_SIGINFO;
      sigemptyset(&sa.sa_mask);
      sa.sa_sigaction = handler;
      sigaction(SIGSEGV, &sa, &old_sigaction);

      a_begin = a;
      a_end = a_begin + nnzA;
      ALIGN(a);

      rowstr_begin = rowstr;
      rowstr_end = rowstr_begin + n_rows + 1;
      ALIGN(rowstr);

      colidx_begin = colidx;
      colidx_end = colidx_begin + nnzA;
      ALIGN(colidx);
    }

    cl_status = clEnqueueWriteBuffer(queue(), x.values, true, 0,
                                     sizeof(double) * x.num_values, iv,
                                     0, nullptr, nullptr);

    status = clsparseDcsrmv(&alpha, &A, &x, &beta, &y, control);
    if(status != clsparseSuccess) {
      std::cout << "Problem performing SPMV!\n";
    }

    cl_status = clEnqueueReadBuffer(queue(), y.values, true, 0, sizeof(double) * y.num_values, ov, 0, nullptr, nullptr);
  }
}

void* f_spmv_harness_(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
  int n_rows = *rows;
  int nnzA = rowstr[n_rows] - rowstr[0];

  auto impl = predict(n_rows, nnzA);
  if(impl == CL_NATIVE_IMPL) {
    f_native_spmv(ov, a, iv, rowstr, colidx, rows);
  } else if(impl == CL_GPU_IMPL) {
    int n_cols = 0;
    for(int i = rowstr[0] - 1; i < rowstr[n_rows] - 1; i++) {
        if(colidx[i] >= n_cols) {
          n_cols = colidx[i];
        }
    }

    f_setup(n_rows, n_cols, nnzA);

    if((af_begin != nullptr && af_begin != a) ||
       (rowstr_begin != nullptr && rowstr_begin != rowstr) ||
       (colidx_begin != nullptr && colidx_begin != colidx))
    {
      af_begin = nullptr;
      rowstr_begin = nullptr;
      colidx_begin = nullptr;
    }

    if(af_begin == nullptr || rowstr_begin == nullptr || colidx_begin == nullptr) {
      cl_status = clEnqueueWriteBuffer(queue(), A.values, true, 0, 
                                       sizeof(float) * A.num_nonzeros, a, 
                                       0, nullptr, nullptr);

      auto one_based = new int[A.num_nonzeros + A.num_rows + 1];
      const auto sub_one = [](int v) { return v - 1; };

      auto it = std::transform(colidx, colidx + A.num_nonzeros, one_based, sub_one);
      std::transform(rowstr, rowstr + A.num_rows + 1, it, sub_one);

      cl_status = clEnqueueWriteBuffer(queue(), A.col_indices, true, 0,
                                       sizeof(int) * A.num_nonzeros, one_based,
                                       0, nullptr, nullptr);

      cl_status = clEnqueueWriteBuffer(queue(), A.row_pointer, true, 0,
                                       sizeof(int) * (A.num_rows + 1), 
                                       one_based + A.num_nonzeros,
                                       0, nullptr, nullptr);

      clsparseCsrMetaCreate(&A, control);

      delete[] one_based;

      cl_status = clEnqueueWriteBuffer(queue(), y.values, true, 0,
                                       sizeof(float) * y.num_values, ov,
                                       0, nullptr, nullptr);

      struct sigaction sa;
      sa.sa_flags = SA_SIGINFO;
      sigemptyset(&sa.sa_mask);
      sa.sa_sigaction = handler;
      sigaction(SIGSEGV, &sa, &old_sigaction);

      af_begin = a;
      af_end = af_begin + nnzA;
      ALIGN(af);

      rowstr_begin = rowstr;
      rowstr_end = rowstr_begin + n_rows + 1;
      ALIGN(rowstr);

      colidx_begin = colidx;
      colidx_end = colidx_begin + nnzA;
      ALIGN(colidx);
    }

    cl_status = clEnqueueWriteBuffer(queue(), x.values, true, 0,
                                     sizeof(float) * x.num_values, iv,
                                     0, nullptr, nullptr);

    status = clsparseDcsrmv(&alpha, &A, &x, &beta, &y, control);
    if(status != clsparseSuccess) {
      std::cout << "Problem performing SPMV!\n";
    }

    cl_status = clEnqueueReadBuffer(queue(), y.values, true, 0, sizeof(float) * y.num_values, ov, 0, nullptr, nullptr);
  }
}

}
