#include "native-impl.h"

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  native_spmv(ov, a, iv, rowstr, colidx, rows);
}
