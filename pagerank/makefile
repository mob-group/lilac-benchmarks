ifeq ($(SPMV_VERSION),mkl)
  C_LIB  = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lmkl-spmv
else ifeq ($(SPMV_VERSION),gpu)
  C_LIB = -L${CUDA_ROOT}/lib64 -lcublas -lcusparse -lcudart -lgpu-spmv
else ifeq ($(SPMV_VERSION),opencl)
  C_LIB = -lopencl-spmv
else ifeq ($(SPMV_VERSION),clgpu)
  C_LIB = -lclgpu-spmv
else
  C_LIB = -lnative-spmv
endif

C_LIB := -L$(SPMV_ROOT)/lib $(C_LIB)

pagerank: pagerank.cc library.cc
	$(CXX) -o $@ $^ -O2 $(C_LIB) -lm

clean:
	rm -f pagerank
