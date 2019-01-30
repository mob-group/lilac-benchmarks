ifeq ($(SPMV_VERSION),mkl)
  C_LIB  = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lmkl-spmv
else ifeq ($(SPMV_VERSION),mkl-slow)
  C_LIB  = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lmkl-slow-spmv
else ifeq ($(SPMV_VERSION),gpu)
  C_LIB = -L${CUDA_ROOT}/lib64 -lcublas -lcusparse -lcudart -lgpu-spmv
else ifeq ($(SPMV_VERSION),opencl00)
  C_LIB = -lopencl00-spmv
else ifeq ($(SPMV_VERSION),opencl01)
  C_LIB = -lopencl01-spmv
else ifeq ($(SPMV_VERSION),opencl10)
  C_LIB = -lopencl10-spmv
else ifeq ($(SPMV_VERSION),opencl11)
  C_LIB = -lopencl11-spmv
else ifeq ($(SPMV_VERSION),sparsex)
  C_LIB = -lsparsex-spmv
else ifeq ($(SPMV_VERSION),sparsex-slow)
  C_LIB = -lsparsex-slow-spmv
else
  C_LIB = -lnative-spmv
endif

C_LIB := -L$(SPMV_ROOT)/lib $(C_LIB)
CXXFLAGS := -std=c++11

bfs: bfs.cc library.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ -O2 -lm $(C_LIB) 

clean:
	rm -f bfs pagerank
