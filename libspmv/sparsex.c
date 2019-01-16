#include <sparsex/sparsex.h>

#include <stdio.h>

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  spx_init();

  int cols = 0;
  for(int i = rowstr[0] - 1; i < rowstr[*rows] - 1; i++) {
    if(colidx[i] >= cols) {
      cols = colidx[i];
    }
  }

  spx_input_t *inp = spx_input_load_csr(rowstr, colidx, a, *rows, cols, SPX_INDEX_ONE_BASED);
  spx_matrix_t *mat = spx_mat_tune(inp);
  spx_partition_t *parts = spx_mat_get_partition(mat);

  spx_vector_t *x = spx_vec_create_from_buff(iv, NULL, cols, parts, SPX_VEC_AS_IS);
  spx_vector_t *y = spx_vec_create_from_buff(ov, NULL, *rows, parts, SPX_VEC_AS_IS);

  spx_matvec_mult(1.0, mat, x, y);
  
  return NULL;
}

void* f_spmv_harness_(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
  return NULL;
}
