CFLAGS=-std=c99 -D_GNU_SOURCE -g
CXXFLAGS=
LDFLAGS=-ldl

all: \
	native.so \
	mkl.so \
	gpu.so \
	test

native.so:	native.c
	${CC} ${CFLAGS} -fPIC -shared $< -o $@

mkl.so:	mkl.c
	${CC} ${CFLAGS} -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl -fPIC -shared $< -o $@

gpu.so:	gpu.c
	${CC} ${CFLAGS} -I${CUDA_ROOT}/include -Wl,-rpath,${CUDA_ROOT}/lib64 -fPIC -shared $< -o $@

test: test.cpp
	${CXX} ${CXXFLAGS} $< -o $@ ${LDFLAGS}

install: all
	cp mkl.so ${SPMV_INSTALL_PREFIX}/lib/libmkl-spmv.so
	cp gpu.so ${SPMV_INSTALL_PREFIX}/lib/libgpu-spmv.so
	cp native.so ${SPMV_INSTALL_PREFIX}/lib/libnative-spmv.so

clean:
	rm -f *.so
	rm -f test
