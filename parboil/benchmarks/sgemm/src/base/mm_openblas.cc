#include "llvm/IDL/harness.hpp"
#include <cblas.h>

namespace {
struct Functor
{
    void operator()(double* output, double* left, double* right, int N, int M, int K) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                     N, M, K, 1.0, left, K, right, N, 0.0, output, N);
    }
};
}

extern "C"
void mm_harness(double* output, double* left, double* right, int N, int M, int K) {
    static Functor functor;
    functor(output, left, right, N, M, K);
}
