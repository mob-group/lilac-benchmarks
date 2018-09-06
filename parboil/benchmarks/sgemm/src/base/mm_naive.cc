#include "harness.hpp"

namespace {
struct Functor
{
    void operator()(double* output, double* left, double* right, int N, int M, int K) {
        int i, j, k;
        for(i = 0; i < N; i++) {
          for(j = 0; j < M; j++) {
            double value = 0.0;
            for(k = 0; k < K; k++)
              value += left[i*K+k] * right[k*M+j];
            output[i*M+j] = value;
          }
        }
    }
};
}

extern "C"
void mm_harness(double* output, double* left, double* right, int N, int M, int K) {
    static Functor functor;
    functor(output, left, right, N, M, K);
}
