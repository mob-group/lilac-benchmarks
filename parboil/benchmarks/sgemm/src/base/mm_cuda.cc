#include "harness.hpp"
#include <cuda_runtime.h>
#include "cublas_v2.h"

namespace {
template<typename type_in, typename type_out>
void CudaRead_update(type_in* in, int size, type_out& out) {
    cudaMemcpy(out, in, sizeof(type_in)*size,
               cudaMemcpyHostToDevice);
}

template<typename type_in, typename type_out>
void CudaRead_construct(int size, type_out& out) {
    cudaMalloc(&out, sizeof(type_in)*size);
}

template<typename type_in, typename type_out>
void CudaRead_destruct(int size, type_out& out) {
    cudaFree(out);
}

template<typename type_in, typename type_out>
using CudaRead = ReadObject<type_in, type_out,
    CudaRead_update<type_in,type_out>,
    CudaRead_construct<type_in,type_out>,
    CudaRead_destruct<type_in,type_out>>;

template<typename type_in, typename type_out>
void CudaWrite_update(type_in* in, int size, type_out& out) {
    cudaMemcpy(in, out, sizeof(type_in)*size,
               cudaMemcpyDeviceToHost);
}

template<typename type_in, typename type_out>
void CudaWrite_construct(int size, type_out& out) {
    cudaMalloc(&out, sizeof(type_in)*size);
}

template<typename type_in, typename type_out>
void CudaWrite_destruct(int size, type_out& out) {
    cudaFree(out);
}

template<typename type_in, typename type_out>
using CudaWrite = WriteObject<type_in, type_out,
    CudaWrite_update<type_in,type_out>,
    CudaWrite_construct<type_in,type_out>,
    CudaWrite_destruct<type_in,type_out>>;

struct Functor
{
    Functor() {
        cublasCreate(&handle);
    }

    ~Functor() {
        cublasDestroy(handle);
    }

    cublasHandle_t handle;
    CudaRead<float,float*> shadow_d_A;
    CudaRead<float,float*> shadow_d_B;
    CudaWrite<float,float*> shadow_d_C;

    void operator()(float* output, float* left, float* right, int N, int M, int K) {
        auto d_A = shadow_d_A(left, M*K);
        auto d_B = shadow_d_B(right, K*N);
        auto d_C = shadow_d_C(output, M*N);
        float alpha = 1.0;
        float beta  = 0.0;
        cublasSgemm_v2(handle, CUBLAS_OP_N, CUBLAS_OP_T,
                       M, N, K, &alpha, d_A, M, d_B, N, &beta, d_C, M);
    }
};
}

extern "C"
void mm_harness(float* output, float* left, float* right, int N, int M, int K) {
    static Functor functor;
    functor(output, left, right, N, M, K);
}
