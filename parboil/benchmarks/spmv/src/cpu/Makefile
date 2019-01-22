# (c) 2010 The Board of Trustees of the University of Illinois

ifeq ($(SPMV_VERSION),mkl)
  C_LIB  = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lmkl-spmv
else ifeq ($(SPMV_VERSION),gpu)
  C_LIB = -L${CUDA_ROOT}/lib64 -lcublas -lcusparse -lcudart -lgpu-spmv
else ifeq ($(SPMV_VERSION),opencl)
  C_LIB = -lopencl-spmv
else ifeq ($(SPMV_VERSION),clgpu)
  C_LIB = -lclgpu-spmv
else ifeq ($(SPMV_VERSION),sparsex)
	C_LIB = -lfsparsex-spmv
else
  C_LIB = -lnative-spmv
endif

C_LIB := -L$(SPMV_ROOT)/lib $(C_LIB)

LANGUAGE=c
TOOLS_SRC=common_src/convert-dataset
SRCDIR_OBJS=main.o file.o 
APP_CFLAGS=-std=c99 -I$(TOOLS_SRC) $(C_LIB) -D_POSIX_C_SOURCE=200809L

include $(TOOLS_SRC)/commontools.mk
