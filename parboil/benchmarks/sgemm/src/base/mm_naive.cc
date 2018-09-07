#include "harness.hpp"

namespace {
struct Functor
{
    void operator()(float* output, const float* left, const float* right, int N, int M, int K) {
      int lda = M;
      int ldb = N;
      int ldc = M;

        for(int nn = 0; nn < N; nn++) {
          for(int mm = 0; mm < M; mm++) {
            float value = 0.0;
            for(int i = 0; i < K; i++)
              value += left[mm + i * lda] * right[nn + i * ldb];
            output[nn*ldc+mm] = value;
          }
        }
    }
};
}

extern "C"
void mm_harness(float* output, const float* left, const float* right, int N, int M, int K) {
    static Functor functor;
    functor(output, left, right, N, M, K);
}
