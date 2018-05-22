CFLAGS=-std=c99 -D_GNU_SOURCE -g

all: \
	native.so \
	mkl.so \
	gpu.so

native.so:	native.c
	${CC} ${CFLAGS} -fPIC -shared $< -o $@

mkl.so:	mkl.c
	${CC} ${CFLAGS} -I${MKLROOT}/include -Wl,-rpath,${MKLROOT}/lib/intel64 -fPIC -shared $< -o $@

gpu.so:	gpu.c
	${CC} ${CFLAGS} -I${CUDA_ROOT}/include -Wl,-rpath,${CUDA_ROOT}/lib64 -fPIC -shared $< -o $@

clean:
	rm *.so
