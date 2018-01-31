/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File linesearch.h: Linesearch for Lbfgs.
 * NOTE: Included from lbfgs.cu, not to be used on its own!
 **/

#ifndef LBFGS_LINESEARCH_H
#define LBFGS_LINESEARCH_H

#include "lbfgs.h"

namespace gpu_lbfgs
{
	// Variables on the GPU. 

	__device__ double fNew;
	__device__ double factor;
	__device__ double stepSize;
	__device__ double evPercent;

	__device__ bool status;

	__constant__ double maxStep;
	__constant__ double maxFkRise;

	//Kernels.

	__global__ void adjustStepSize();
	__global__ void checkMaxFkRiseAndEvPc(const double *d_fk, const bool isRayleighRitz);
	__global__ void reduceStepSize();
}



bool Lbfgs::linesearch(double *d_x, double *d_z, double *d_fk, double *d_gk, size_t it, Lbfgs::lStatus &lbfgsStatus, double *d_step,
		double *d_tmp, const double *m_d_zWork, double &evPercent, double *outRms, double *d_gkSave)
{
	using namespace gpu_lbfgs;

	bool printingOn = m_debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = m_debugPrinting.getFileHandle();

	bool isTimingLinesearch = m_timer_linesearch.getTimingOn();

	const size_t numDimensions = m_costFunction.getNumDimensions();

	const bool isRayleighRitz  = m_costFunction.getIsRayleighRitz();

	double *d_xNew, *d_gNew, *d_fNew, *d_factor, *d_stepSize;

	const double one      =  1.0;
	const double minusOne = -1.0;

	double gNorm, dotted;

	CudaSafeCall( cudaMalloc((void**)&d_xNew, numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc((void**)&d_gNew, numDimensions * sizeof(double)) );

	CudaSafeCall( cudaGetSymbolAddress((void**)&d_fNew,     gpu_lbfgs::fNew) );
	CudaSafeCall( cudaGetSymbolAddress((void**)&d_factor,   gpu_lbfgs::factor) );
	CudaSafeCall( cudaGetSymbolAddress((void**)&d_stepSize, gpu_lbfgs::stepSize) );

	bool isFirstStep = (m_cumulativeIter == 0);
	if (isFirstStep) {
		// Make the first guess for the step length cautious. 
		m_cublas.dispatchNrm2(numDimensions, &gNorm, d_gk, false); // gNorm = sqrt(gk Dot gk)	
		double firstFactor = std::min(1/gNorm, gNorm);
		CudaSafeCall( cudaMemcpy(d_factor, &firstFactor, sizeof(double), cudaMemcpyHostToDevice) );
	}
	else {
		CudaSafeCall( cudaMemcpy(d_factor, &one, sizeof(double), cudaMemcpyHostToDevice) );
	}

	// If the step is pointing uphill, invert it. 
	m_cublas.dispatchDot(numDimensions, &dotted, d_z, d_gk, false); // dotted = z Dot gk
	bool isStepUphill = (dotted > 0.0);
	if (isStepUphill) {
		if (printingOn) {
			fileHandle << " Warning: step direction was uphill (positive projection onto gradient) - reversing" 
				<< std::endl;
		}
		m_cublas.dispatchScale(numDimensions, d_z, d_z, &minusOne, false); // z = -z
	}

	m_cublas.dispatchNrm2(numDimensions, d_stepSize, d_z); // stepSize = sqrt(z Dot z)
	// Make sure the step is no larger than maxStep. 
	adjustStepSize<<<1,1>>>();		
	CudaCheckError();
	cudaDeviceSynchronize();

	int nRed = 0;
	for (nRed = 0; nRed < m_nRedMax; ++nRed) {
		m_cublas.dispatchAxpy(numDimensions, d_xNew, d_x, d_z, d_factor); // xNew = x + z*factor

		if (isTimingLinesearch) {
			m_timer_linesearch.stop();
		}

		// Calculate energy and gradient for proposed step. 
		m_costFunction.basePotential(d_xNew, d_fNew, d_gNew);

		bool isColdFusion = m_costFunction.getIsColdFusion();
		if (isColdFusion) {
			lbfgsStatus = LBFGS_COLD_FUSION_DIAGNOSED;
			return false;
		}

		if (isTimingLinesearch) {
			m_timer_linesearch.start();
		}

		if (m_projectGrad) {
			// Save the true gradient. 
			CudaSafeCall( cudaMemcpy(d_gkSave, d_gNew, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );

			m_cublas.dispatchDot(numDimensions, d_tmp, m_d_zWork, d_gNew); // tmp = zWork Dot gNew

			// tmp = -tmp
			updateNegative<<<1, 1>>>(d_tmp);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy(numDimensions, d_gNew, d_gNew, m_d_zWork, d_tmp); // gNew += zWork*tmp
		}
		m_costFunction.calculateRms(d_gNew, outRms);

		// Check whether energy has risen too much and % change in eigenvalue if R-R minimization. 
		checkMaxFkRiseAndEvPc<<<1,1>>>(d_fk, isRayleighRitz);
		CudaCheckError();
		cudaDeviceSynchronize();

		if (isRayleighRitz) {
			CudaSafeCall( cudaMemcpyFromSymbol(&evPercent, gpu_lbfgs::evPercent, sizeof(double)) );
		}

		int isBelowMaxFkRise=false;
		CudaSafeCall( cudaMemcpyFromSymbol(&isBelowMaxFkRise, gpu_lbfgs::status, sizeof(bool)) );

		double b, testFactor, testStepSize;
		if (isBelowMaxFkRise) {
			if (printingOn) {
				CudaSafeCall( cudaMemcpy(&testFactor,   d_factor,   sizeof(double), cudaMemcpyDeviceToHost) );
				CudaSafeCall( cudaMemcpy(&testStepSize, d_stepSize, sizeof(double), cudaMemcpyDeviceToHost) );
				b = testFactor * testStepSize;
				fileHandle << " LBFGS step = " << std::setw(13) << std::fixed << std::setprecision(9) << b 
					<< std::endl;
			}
			break;
		}
		else {
			// If energy rose too much, reduce stepSize and continue the loop. 
			reduceStepSize<<<1,1>>>();
			CudaCheckError();
			cudaDeviceSynchronize();
			if (printingOn) {
				double testfk, testfNew;
				CudaSafeCall( cudaMemcpy(&testFactor,   d_factor,   sizeof(double), cudaMemcpyDeviceToHost) );
				CudaSafeCall( cudaMemcpy(&testStepSize, d_stepSize, sizeof(double), cudaMemcpyDeviceToHost) );
				CudaSafeCall( cudaMemcpy(&testfk,       d_fk,       sizeof(double), cudaMemcpyDeviceToHost) );
				CudaSafeCall( cudaMemcpy(&testfNew,     d_fNew,     sizeof(double), cudaMemcpyDeviceToHost) );
				b = testFactor * testStepSize;

				fileHandle << " Function increased from " << testfk << " to " << testfNew 
					<< " - reducing step size to " << b << std::endl;
			}
		}
	}

	bool isAboveMaxIt = (nRed >= m_nRedMax);
	if (isAboveMaxIt) {
		if (printingOn) {
			fileHandle << "Warning: a step size cannot be found where the maximum allowed function increase is not exceeded" 
				<< std::endl;
		}
	}

	CudaSafeCall( cudaMemcpy(d_x,  d_xNew, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );
	CudaSafeCall( cudaMemcpy(d_gk, d_gNew, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );
	CudaSafeCall( cudaMemcpy(d_fk, d_fNew,                 sizeof(double), cudaMemcpyDeviceToDevice) );

	CudaSafeCall( cudaMemcpy(d_step, d_factor, sizeof(double), cudaMemcpyDeviceToDevice) );

	if (printingOn) {
		double printEne = 0.0;
		CudaSafeCall( cudaMemcpy(&printEne, d_fk, sizeof(double), cudaMemcpyDeviceToHost) );
		fileHandle << " Function value and RMS = " << std::setw(20) << std::fixed << std::setprecision(10) << printEne 
			<< std::setw(20) << std::fixed << std::setprecision(10) << *outRms << " after " 
			<< std::setw(7) << (it + 1) << " LBFGS steps" << std::endl;
		if (isRayleighRitz) {
			fileHandle << " Eigenvalue % change = " << std::setw(13) << std::fixed << std::setprecision(9) << evPercent 
				<< std::endl;
		}
	}

	CudaSafeCall( cudaFree(d_xNew) );
	CudaSafeCall( cudaFree(d_gNew) );
	return true;
}



namespace gpu_lbfgs
{
	__global__ void adjustStepSize()
	{
		double a = factor * stepSize;
		if (a > maxStep) {
			factor = maxStep / stepSize;
		}
	}

	__global__ void checkMaxFkRiseAndEvPc(const double *d_fk, const bool isRayleighRitz)
	{
		double df;
		df = fNew - *d_fk;
		if (isRayleighRitz) {
			df = df/fmax(abs(fNew), 1.0e-100);
			evPercent = 100.0*abs(df);
		}

		if (df < maxFkRise) {
			status = true;
		}
		else {
			status = false;
		}
	}

	__global__ void reduceStepSize()
	{
		factor /= 10.0;
	}
}

#endif // LBFGS_LINESEARCH_H

