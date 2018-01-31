/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/) 
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File cublas.h: cuBLAS functionality. 
 *
 **/

#ifndef CUBLAS_H
#define CUBLAS_H

#include "error_checking.h"

class Cublas 
{
	public:
		Cublas();
		~Cublas();

		// axpy  computes  dst = a * x + y
		// scale computes  dst = a * x
		// dot   computes  dst = x^T y
		// nrm2  computes  dst = sqrt(x^T y)
		//
		// x, y, dest (for axpy / scale) are n-vectors,
		// a,    dest (for dot)          are scalars.
		//
		// aDevicePointer / dstDevicePointer indicate whether
		// dst and a point to memory on the device or host.
		// All other pointers (marked with a d_) must point to device memory.

		void dispatchAxpy (const size_t n, double *d_dst, const double *d_y, const double *d_x, const double *a, 
				bool isDevicePointer    = true) const;
		void dispatchScale(const size_t n, double *d_dst, const double *d_x,                    const double *a, 
				bool isDevicePointer    = true) const;
		void dispatchDot  (const size_t n, double *dst,   const double *d_x, const double *d_y,                 
				bool isDstDevicePointer = true) const;
		void dispatchNrm2 (const size_t n, double *dst,   const double *d_x, 
				bool isDevicePointer   = true) const;

	private:
		mutable cublasHandle_t m_cublasHandle;

};

#endif /* end of include guard: CUBLAS_H */
