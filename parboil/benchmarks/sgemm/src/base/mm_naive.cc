#include "harness.hpp"

namespace {
struct Functor
{
    void operator()(float* output, const float* left, const float* right, int N, int M, int K) {
      int lda = M;
      int ldb = N;
      int ldc = M;

      for(int i = 0; i < K; i++) {
        for(int nn = 0; nn < N; nn++) {
          for(int mm = 0; mm < M; mm++) {
              output[nn*ldc + mm] += left[mm + i * lda] * right[nn + i * ldb];
            }
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
