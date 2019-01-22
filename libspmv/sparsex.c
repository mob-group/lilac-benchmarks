#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <sys/mman.h>
#include <unistd.h>

#include <sparsex/sparsex.h>

struct sigaction old_sigaction;
static void* data_begin    = 0;
static void* data_end      = 0;
static void* aligned_begin = 0;
static void* aligned_end   = 0;

static void handler(int sig, siginfo_t *si, void *unused)
{
    if(si->si_addr >= aligned_begin && si->si_addr < aligned_end)
    {
        mprotect(aligned_begin, aligned_end-aligned_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
        data_begin = 0;
    }
    else if(old_sigaction.sa_sigaction)
    {
        old_sigaction.sa_sigaction(sig, si, unused);
    }
}

void setup(void)
{
  static bool setup = false;
  if(!setup) {
    spx_init();
    spx_option_set("spx.rt.nr_threads", "24");
    setup = true;
  }
}

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  static int cols = 0;
  static spx_matrix_t *mat = NULL;
  static spx_partition_t *parts = NULL;

  setup();

  if(data_begin != 0 && data_begin != colidx)
  {
      mprotect(aligned_begin, aligned_end-aligned_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
      data_begin = 0;
  }
  if(data_begin == 0)
  {
      cols = 0;
      for(int i = rowstr[0] - 1; i < rowstr[*rows] - 1; i++) {
          if(colidx[i] >= cols) {
            cols = colidx[i];
          }
      }
      data_begin = colidx;
      data_end   = colidx + rowstr[*rows];

      spx_input_t *inp = spx_input_load_csr(rowstr, colidx, a, *rows, cols, SPX_INDEX_ONE_BASED);
      mat = spx_mat_tune(inp);
      parts = spx_mat_get_partition(mat);

      struct sigaction sa;
      sa.sa_flags = SA_SIGINFO;
      sigemptyset(&sa.sa_mask);
      sa.sa_sigaction = handler;

      sigaction(SIGSEGV, &sa, &old_sigaction);

      aligned_begin = data_begin;
      aligned_end   = data_end + sysconf(_SC_PAGE_SIZE) - 1;
      aligned_begin -= (size_t)aligned_begin % sysconf(_SC_PAGE_SIZE);
      aligned_end   -= (size_t)aligned_end   % sysconf(_SC_PAGE_SIZE);

      mprotect(aligned_begin, aligned_end-aligned_begin, PROT_READ | PROT_EXEC);
  }

  spx_vector_t *x = spx_vec_create_from_buff(iv, NULL, cols, parts, SPX_VEC_AS_IS);
  spx_vector_t *y = spx_vec_create_from_buff(ov, NULL, *rows, parts, SPX_VEC_AS_IS);

  spx_matvec_mult(1.0, mat, x, y);

  return NULL;
}

void* f_spmv_harness_(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
  return NULL;
}
