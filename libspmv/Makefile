CFLAGS=-std=c99 -D_GNU_SOURCE -O3
CXXFLAGS=-I${CL_ROOT}/include -L${CL_ROOT}/lib -L${CL_ROOT}/lib64 -std=c++11 -O3
LDFLAGS=-ldl

all: test

native.so: native.c
	${CC} ${CFLAGS} -fPIC -shared $< -o $@

mkl.so:	mkl.c
	${CC} ${CFLAGS} -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl -fPIC -shared $< -o $@

gpu.so:	gpu.c
	${CC} ${CFLAGS} -I${CUDA_ROOT}/include -Wl,-rpath,${CUDA_ROOT}/lib64 -fPIC -shared $< -o $@

opencl.so: opencl.cpp
	${CXX} ${CXXFLAGS} -fPIC -shared $< -lOpenCL -lclSPARSE -lclsparseTimer -o $@

clgpu.so: clgpu.cpp
	${CXX} ${CXXFLAGS} -fPIC -shared $< -lOpenCL -lclSPARSE -lclsparseTimer -o $@

test: test.cpp
	${CXX} ${CXXFLAGS} $< -o $@ ${LDFLAGS}

install: $(platform).so
	cp $< ${SPMV_ROOT}/lib/lib$(platform)-spmv.so

clean:
	rm -f *.so
	rm -f test
