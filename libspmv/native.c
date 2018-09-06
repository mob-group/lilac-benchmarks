#include "native-impl.h"

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  native_spmv(ov, a, iv, rowstr, colidx, rows);
}

void* f_spmv_harness_(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
  f_native_spmv(ov, a, iv, rowstr, colidx, rows);
}
