CFLAGS=-std=c99 -D_GNU_SOURCE -O3
CXXFLAGS=-I${CL_ROOT}/include -L${CL_ROOT}/lib -L${CL_ROOT}/lib64 -std=c++11 -O3
LDFLAGS=-ldl

SPARSEX_LLVM_LIBDIR=$(shell $(LLVM_CONFIG) --libdir)
SPARSEX_CFLAGS=$(shell $(SPARSEX_CONFIG) --cppflags)
SPARSEX_LDFLAGS=$(shell $(SPARSEX_CONFIG) --ldflags)

all: test

native.so: native.c native-impl.c
	${CC} ${CFLAGS} -fPIC -shared $^ -o $@

mkl.so:	mkl.c mkl-model.c native-impl.c
	${CC} ${CFLAGS} -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl -fPIC -shared $^ -o $@

gpu.so:	gpu.c
	${CC} ${CFLAGS} -L${CUDA_ROOT}/lib64 -I${CUDA_ROOT}/include -lcublas -lcusparse -lcudart -fPIC -shared $< -o $@

unsafe.so: gpu-unsafe.c
	${CC} ${CFLAGS} -L${CUDA_ROOT}/lib64 -I${CUDA_ROOT}/include -lcublas -lcusparse -lcudart -fPIC -shared $< -o $@

opencl.so: opencl.cpp
	${CXX} ${CXXFLAGS} -fPIC -shared $< -lOpenCL -lclSPARSE -lclsparseTimer -o $@

clgpu.so: clgpu.cpp clgpu-model.c native-impl.c
	${CXX} ${CXXFLAGS} -fPIC -shared $^ -lOpenCL -lclSPARSE -lclsparseTimer -o $@

sparsex.so: sparsex.c
	$(CC) $(CFLAGS) -fPIC -shared $(SPARSEX_CFLAGS) $^ -o $@ $(SPARSEX_LDFLAGS) -L$(SPARSEX_LLVM_LIBDIR)

test: test.cpp
	${CXX} ${CXXFLAGS} $< -o $@ ${LDFLAGS}

install: $(platform).so
	cp $< ${SPMV_ROOT}/lib/lib$(platform)-spmv.so

clean:
	rm -f *.so
	rm -f test
