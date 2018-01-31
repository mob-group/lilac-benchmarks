MODULE FFTW3

  USE PREC, ONLY: INT32

  integer ( kind = INT32 ), parameter :: fftw_r2hc = 0
  integer ( kind = INT32 ), parameter :: fftw_hc2r = 1
  integer ( kind = INT32 ), parameter :: fftw_dht = 2
  integer ( kind = INT32 ), parameter :: fftw_redft00 = 3
  integer ( kind = INT32 ), parameter :: fftw_redft01 = 4
  integer ( kind = INT32 ), parameter :: fftw_redft10 = 5
  integer ( kind = INT32 ), parameter :: fftw_redft11 = 6
  integer ( kind = INT32 ), parameter :: fftw_rodft00 = 7
  integer ( kind = INT32 ), parameter :: fftw_rodft01 = 8
  integer ( kind = INT32 ), parameter :: fftw_rodft10 = 9
  integer ( kind = INT32 ), parameter :: fftw_rodft11 = 10
  integer ( kind = INT32 ), parameter :: fftw_forward = -1
  integer ( kind = INT32 ), parameter :: fftw_backward = +1
  integer ( kind = INT32 ), parameter :: fftw_measure = 0
  integer ( kind = INT32 ), parameter :: fftw_destroy_input = 1
  integer ( kind = INT32 ), parameter :: fftw_unaligned = 2
  integer ( kind = INT32 ), parameter :: fftw_conserve_memory = 4
  integer ( kind = INT32 ), parameter :: fftw_exhaustive = 8
  integer ( kind = INT32 ), parameter :: fftw_preserve_input = 16
  integer ( kind = INT32 ), parameter :: fftw_patient = 32
  integer ( kind = INT32 ), parameter :: fftw_estimate = 64
  integer ( kind = INT32 ), parameter :: fftw_estimate_patient = 128
  integer ( kind = INT32 ), parameter :: fftw_believe_pcost = 256
  integer ( kind = INT32 ), parameter :: fftw_dft_r2hc_icky = 512
  integer ( kind = INT32 ), parameter :: fftw_nonthreaded_icky = 1024
  integer ( kind = INT32 ), parameter :: fftw_no_buffering = 2048
  integer ( kind = INT32 ), parameter :: fftw_no_indirect_op = 4096
  integer ( kind = INT32 ), parameter :: fftw_allow_large_generic = 8192
  integer ( kind = INT32 ), parameter :: fftw_no_rank_splits = 16384
  integer ( kind = INT32 ), parameter :: fftw_no_vrank_splits = 32768
  integer ( kind = INT32 ), parameter :: fftw_no_vrecurse = 65536
  integer ( kind = INT32 ), parameter :: fftw_no_simd = 131072

END MODULE FFTW3
