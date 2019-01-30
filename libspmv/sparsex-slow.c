#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/sysinfo.h>
#include <unistd.h>

#include <sparsex/sparsex.h>

void setup(void)
{
  static bool setup = false;
  if(!setup) {
    char *cores;
    asprintf(&cores, "%d", get_nprocs());

    spx_init();
    spx_option_set("spx.rt.nr_threads", cores);
    setup = true;

    free(cores);
  }
}

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  int cols = 0;
  spx_matrix_t *mat = NULL;
  spx_partition_t *parts = NULL;

  setup();

  for(int i = rowstr[0] - 1; i < rowstr[*rows] - 1; i++) {
      if(colidx[i] >= cols) {
        cols = colidx[i];
      }
  }

  spx_input_t *inp = spx_input_load_csr(rowstr, colidx, a, *rows, cols, SPX_INDEX_ONE_BASED);
  mat = spx_mat_tune(inp);
  parts = spx_mat_get_partition(mat);

  spx_vector_t *x = spx_vec_create_from_buff(iv, NULL, cols, parts, SPX_VEC_AS_IS);
  spx_vector_t *y = spx_vec_create_from_buff(ov, NULL, *rows, parts, SPX_VEC_AS_IS);

  spx_matvec_mult(1.0, mat, x, y);

  return NULL;
}

void* f_spmv_harness_(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
  return NULL;
}
