#include "harness.hpp"
#include <cblas.h>

#include <iostream>

namespace {
struct Functor
{
    void operator()(float* output, float* left, float* right, int N, int M, int K) {
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                     M, N, K, 1.0, left, M, right, N, 0.0, output, M);
    }
};
}

extern "C"
void mm_harness(float* output, float* left, float* right, int N, int M, int K) {
    static Functor functor;
    functor(output, left, right, N, M, K);
}
