#ifndef NATIVE_IMPL_H
#define NATIVE_IMPL_H

void* native_spmv(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows);
void* f_native_spmv(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows);

#endif
