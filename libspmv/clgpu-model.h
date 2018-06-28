#ifndef CLGPU_MODEL_H
#define CLGPU_MODEL_H

#define CL_NATIVE_IMPL 0
#define CL_GPU_IMPL 1

int predict(int rows, int nnz);

#endif
