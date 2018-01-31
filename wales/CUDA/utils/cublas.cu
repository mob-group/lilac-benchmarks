/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/) 
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File cublas.cu: Implementation of class Cublas.
 *
 **/

#include "cublas.h"

Cublas::Cublas()
{
	CublasSafeCall( cublasCreate(&m_cublasHandle) );
}

Cublas::~Cublas()
{
	CublasSafeCall( cublasDestroy(m_cublasHandle) );
}

// Vector operations

void Cublas::dispatchAxpy(const size_t n, double *d_dst, const double *d_y, const double *d_x, const double *a, bool isDevicePointer) const
{
	const cublasPointerMode_t mode = isDevicePointer ? CUBLAS_POINTER_MODE_DEVICE
		: CUBLAS_POINTER_MODE_HOST;

	CublasSafeCall( cublasSetPointerMode(m_cublasHandle, mode) );

	if (d_dst != d_y)
		CudaSafeCall( cudaMemcpy(d_dst, d_y, n * sizeof(double), cudaMemcpyDeviceToDevice) );

	CublasSafeCall( cublasDaxpy(m_cublasHandle, int(n), a, d_x, 1, d_dst, 1) );
}

void Cublas::dispatchScale(const size_t n, double *d_dst, const double *d_x, const double *a, bool isDevicePointer) const
{
	const cublasPointerMode_t mode = isDevicePointer ? CUBLAS_POINTER_MODE_DEVICE
		: CUBLAS_POINTER_MODE_HOST;

	CublasSafeCall( cublasSetPointerMode(m_cublasHandle, mode) );

	if (d_dst != d_x)
		CudaSafeCall( cudaMemcpy(d_dst, d_x, n * sizeof(double), cudaMemcpyDeviceToDevice) );

	CublasSafeCall( cublasDscal(m_cublasHandle, int(n), a, d_dst, 1) );
}

void Cublas::dispatchDot(const size_t n, double *dst, const double *d_x, const double *d_y, bool isDstDevicePointer) const
{
	const cublasPointerMode_t mode = isDstDevicePointer ? CUBLAS_POINTER_MODE_DEVICE
		: CUBLAS_POINTER_MODE_HOST;

	CublasSafeCall( cublasSetPointerMode(m_cublasHandle, mode) );

	CublasSafeCall( cublasDdot(m_cublasHandle, int(n), d_x, 1, d_y, 1, dst) );
}

void Cublas::dispatchNrm2(const size_t n, double *dst, const double *d_x, bool isDevicePointer) const
{
	const cublasPointerMode_t mode = isDevicePointer ? CUBLAS_POINTER_MODE_DEVICE
		: CUBLAS_POINTER_MODE_HOST;

	CublasSafeCall( cublasSetPointerMode(m_cublasHandle, mode) );

	CublasSafeCall( cublasDnrm2(m_cublasHandle, int(n), d_x, 1, dst) );
}

