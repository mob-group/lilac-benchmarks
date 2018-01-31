/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/) 
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File lbfgs.cu: Implementation of class Lbfgs. 
 *
 **/

#include <iomanip>

#include "lbfgs.h"

namespace gpu_lbfgs
{

	// Variables on the GPU. 

	__device__ double fkm1;
	__device__ double tmp;
	__device__ double step;
	__device__ double tmp2;

	// Kernels. 

	// Small helper kernels for scalar operations in device memory needed during updates. *** Use with a single thread only! ***
	__global__ void update1   (double *alphaOut, const double *sDotZ, const double *rho, double *minusAlphaOut); // First update loop. 
	__global__ void update2   (double *alphaMinusBetaOut, const double *rho, const double *yDotZ, const double *alpha); // Second update loop. 
	__global__ void update3   (double *rhoOut, double *H0Out, double *yDotS, double *yDotY); // After linesearch.

	__global__ void updateNegative(double *tmpOut);
}



// File linesearch.h is no real header, it contains
// part of the implementation and must be included
// after the variables above have been declared.
#include "linesearch.h" 

Lbfgs::Lbfgs(CostFunction &costFunc, Printing &debugPrinting, Timer &timer_total, Timer &timer_updates, 
		Timer &timer_linesearch, Cublas &cublas, size_t numDimensions, double H0, int updates, int maxIter, 
		double gradEps, double maxStep, double maxFkRise, int nRedMax)
	: m_costFunction(costFunc)
	, m_debugPrinting(debugPrinting)
	, m_timer_total(timer_total)
	, m_timer_updates(timer_updates)
	, m_timer_linesearch(timer_linesearch)
	, m_cublas(cublas)
	, m_maxIter(maxIter)
	, m_cumulativeIter(0)
	, m_maxStep(maxStep)
	, m_gradientEps(gradEps)
	, m_maxFkRise(maxFkRise)
	, m_updates(updates)
	, m_evalPercentEps(1.0)
	, m_nRedMax(nRedMax)
	, m_H0(H0)
	  , m_projectGrad(false)
{
	CudaSafeCall( cudaMalloc(&m_d_H0, sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_s, m_updates * numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_y, m_updates * numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_alpha, m_updates * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_rho, m_updates * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_zWork, numDimensions * sizeof(double)) );

	CudaSafeCall( cudaMemcpy(m_d_H0, &m_H0, sizeof(double), cudaMemcpyHostToDevice) );

	CudaSafeCall( cudaMemcpyToSymbol(gpu_lbfgs::maxStep, &m_maxStep,  sizeof(double)) );
	CudaSafeCall( cudaMemcpyToSymbol(gpu_lbfgs::maxFkRise, &m_maxFkRise,  sizeof(double)) );
}



Lbfgs::~Lbfgs()
{
	CudaSafeCall( cudaFree(m_d_H0) );
	CudaSafeCall( cudaFree(m_d_s) );
	CudaSafeCall( cudaFree(m_d_y) );
	CudaSafeCall( cudaFree(m_d_alpha));
	CudaSafeCall( cudaFree(m_d_rho) );
	CudaSafeCall( cudaFree(m_d_zWork) );
}



// Get/set functions. 

void   Lbfgs::setMaxIterations(size_t maxIter)
{
	m_maxIter = maxIter;
}

void   Lbfgs::setEvalPercentEps(double evalPercentEps)
{
	m_evalPercentEps = evalPercentEps;
}

void Lbfgs::setProjectGrad(bool projectGrad)
{
	m_projectGrad = projectGrad;
}

void Lbfgs::setDevzWork(double *d_zWork, const size_t arraySize)
{
	CudaSafeCall( cudaMemcpy(m_d_zWork, d_zWork, arraySize * sizeof(double), cudaMemcpyDeviceToDevice) );
}



std::string Lbfgs::statusToString(Lbfgs::lStatus lbfgsStatus)
{
	switch (lbfgsStatus) {
		case LBFGS_CONVERGED :
			return "Convergence achieved";
		case LBFGS_REACHED_MAX_ITER :
			return "Reached maximum number of iterations";
		case LBFGS_COLD_FUSION_DIAGNOSED :
			return "Cold fusion diagnosed";
		case LBFGS_RMS_IS_NAN :
			return "RMS force was NaN - now reset to 1.0e-100 - if using AMBER igb=1 this can be due to negative Born radii";
		default :
			return "Unknown status";
	}
}



Lbfgs::lStatus Lbfgs::minimize(double *d_x, double *d_fk, double *d_gk, double *outRms, int *itDone, bool resetHistory)
{
	bool printingOn = m_debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = m_debugPrinting.getFileHandle();

	if (resetHistory) {
		m_cumulativeIter = 0;
		CudaSafeCall( cudaMemcpy(m_d_H0, &m_H0, sizeof(double), cudaMemcpyHostToDevice) );
	}

	if (printingOn) {
		if (resetHistory) {
			fileHandle << " Resetting LBFGS minimizer" << std::endl;
		}
		else {
			fileHandle << " Not resetting LBFGS minimizer" << std::endl;
		}
	}

	bool isTimingTotal = m_timer_total.getTimingOn();
	bool isTimingUpdates = m_timer_updates.getTimingOn();
	bool isTimingLinesearch = m_timer_linesearch.getTimingOn();

	if (isTimingTotal) {
		m_timer_total.start();
	}

	using namespace gpu_lbfgs;
	const size_t numDimensions = m_costFunction.getNumDimensions();
	const bool isRayleighRitz = m_costFunction.getIsRayleighRitz();

	double *d_fkm1;   // f_(k-1), function value at x_(k-1). 
	double *d_gkm1;   // g_(k-1), gradient at x_(k-1). 
	double *d_gkSave; // True gradient with no projection. 
	double *d_z;      // Search direction. 

	double *d_step;         // Current step length. 
	double *d_tmp, *d_tmp2; // Temporary storage for intermediate results. 

	// Allocations. 

	CudaSafeCall( cudaMalloc(&d_gkm1,   numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_gkSave, numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_z,      numDimensions * sizeof(double)) );

	// Addresses of global symbols. 

	CudaSafeCall( cudaGetSymbolAddress((void**)&d_fkm1, gpu_lbfgs::fkm1) );
	CudaSafeCall( cudaGetSymbolAddress((void**)&d_tmp,  gpu_lbfgs::tmp ) );
	CudaSafeCall( cudaGetSymbolAddress((void**)&d_tmp2, gpu_lbfgs::tmp2) );
	CudaSafeCall( cudaGetSymbolAddress((void**)&d_step, gpu_lbfgs::step) );

	if (printingOn) {
		fileHandle << " CUDA LBFGS" << std::endl;
	}

	if (isRayleighRitz) {
		bool shouldNormalise = true;
		m_costFunction.orthogopt(d_x, shouldNormalise);
	}

	// Initialize. 

	m_costFunction.basePotential(d_x, d_fk, d_gk);

	if (m_projectGrad) {
		// Save the true gradient. 
		CudaSafeCall( cudaMemcpy(d_gkSave, d_gk, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );

		double thisDot;
		m_cublas.dispatchDot(numDimensions, &thisDot, m_d_zWork, m_d_zWork, false); // thisDot = zWork Dot zWork
		bool iszWorkTooLarge = (abs(thisDot - 1.0) > 1.0e-10); 
		if (iszWorkTooLarge) {
			thisDot = 1.0/sqrt(thisDot);
			if (printingOn) {
				fileHandle << " Renormalising uphill vector by factor of " << thisDot << std::endl;
			}
			m_cublas.dispatchScale(numDimensions, m_d_zWork, m_d_zWork, &thisDot, false); // zWork = zWork*thisDot
		}
		m_cublas.dispatchDot(numDimensions, d_tmp, m_d_zWork, d_gk); // tmp = zWork Dot gk

		// tmp = -tmp
		updateNegative<<<1, 1>>>(d_tmp);
		CudaCheckError();
		cudaDeviceSynchronize();

		m_cublas.dispatchAxpy(numDimensions, d_gk, d_gk, m_d_zWork, d_tmp); // gk += zWork*tmp
	}	

	m_costFunction.calculateRms(d_gk, outRms);

	double evPercent = 0.0;

	lStatus lbfgsStatus = LBFGS_REACHED_MAX_ITER;
	size_t it = 0;

	if (printingOn) {
		double printEne = 0.0;
		CudaSafeCall( cudaMemcpy(&printEne, d_fk, sizeof(double), cudaMemcpyDeviceToHost) );
		fileHandle << " Function value and RMS = " 
			<< std::setw(20) << std::fixed << std::setprecision(10) << printEne 
			<< std::setw(20) << std::fixed << std::setprecision(10) << *outRms << " after " 
			<< std::setw(7) << it << " LBFGS steps" << std::endl;
	}

	bool isColdFusion = m_costFunction.getIsColdFusion();
	if (isColdFusion) {
		lbfgsStatus = LBFGS_COLD_FUSION_DIAGNOSED;
		it = m_maxIter;
	}

	for (it; it < m_maxIter; ++it) {

		if  (isnan(*outRms)){
			lbfgsStatus = LBFGS_RMS_IS_NAN;
			break;
		}

		// Check for convergence. 
		bool isRmsBelowGradEps = (*outRms < m_gradientEps);
		bool isEvPercentageChangeBelowTol = (evPercent < m_evalPercentEps);
		if (isRmsBelowGradEps && isEvPercentageChangeBelowTol) {
			lbfgsStatus = LBFGS_CONVERGED;
			break;
		}

		// Find search direction. 

		if (isTimingUpdates) {
			m_timer_updates.start();
		}

		const double minusOne = -1.0;
		m_cublas.dispatchScale(numDimensions, d_z, d_gk, &minusOne, false); // z = -gk

		const size_t MAX_IDX = std::min<size_t>(m_cumulativeIter, m_updates);

		for (size_t i = 1; i <= MAX_IDX; ++i) {
			size_t idx = ((m_cumulativeIter - i) % m_updates);

			m_cublas.dispatchDot(numDimensions, d_tmp, m_d_s + idx * numDimensions, d_z); // tmp = s Dot z

			// alpha = tmp * rho
			// tmp = -alpha
			update1<<<1, 1>>>(m_d_alpha + idx, d_tmp, m_d_rho + idx, d_tmp);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy(numDimensions, d_z, d_z, m_d_y + idx * numDimensions, d_tmp); // z += tmp*y
		}

		m_cublas.dispatchScale(numDimensions, d_z, d_z, m_d_H0); // z = H0*z

		for (size_t i = MAX_IDX; i > 0; --i) {
			size_t idx = ((m_cumulativeIter - i) % m_updates);

			m_cublas.dispatchDot(numDimensions, d_tmp, m_d_y + idx * numDimensions, d_z); // tmp = y Dot z

			// beta = rho * tmp
			// tmp = alpha - beta
			update2<<<1, 1>>>(d_tmp, m_d_rho + idx, d_tmp, m_d_alpha + idx);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy(numDimensions, d_z, d_z, m_d_s + idx * numDimensions, d_tmp); // z += tmp*s
		}

		if (isTimingUpdates) {
			m_timer_updates.stop();
		}

		if (m_projectGrad) {
			m_cublas.dispatchDot(numDimensions, d_tmp, m_d_zWork, d_z); // tmp = zWork Dot z

			// tmp = -tmp
			updateNegative<<<1, 1>>>(d_tmp);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy(numDimensions, d_z, d_z, m_d_zWork, d_tmp); // z += zWork*tmp
		}

		CudaSafeCall( cudaMemcpy(d_fkm1, d_fk,                 sizeof(double), cudaMemcpyDeviceToDevice) ); // fkm1 = fk;
		CudaSafeCall( cudaMemcpy(d_gkm1, d_gk, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) ); // gkm1 = gk;

		if (isTimingLinesearch) {
			m_timer_linesearch.start();
		}

		// The linesearch is defined in linesearch.h. 
		if (!linesearch(d_x, d_z, d_fk, d_gk, it, lbfgsStatus, d_step,
					d_tmp, m_d_zWork, evPercent, outRms, d_gkSave)) {
			break;
		}

		if (isTimingLinesearch) {
			m_timer_linesearch.stop();
		}
		if (isTimingUpdates) {
			m_timer_updates.start();
		}

		// Update s, y, rho and H_0. 

		// s   = x_k - x_{k-1} = step * z
		// y   = g_k - g_{k-1}
		// rho = 1 / (y^T s)
		// H_0 = (y^T s) / (y^T y)

		double *d_curS = m_d_s + (m_cumulativeIter % m_updates) * numDimensions;
		double *d_curY = m_d_y + (m_cumulativeIter % m_updates) * numDimensions;

		m_cublas.dispatchScale(numDimensions, d_curS, d_z,  d_step); // s = step*z
		m_cublas.dispatchAxpy (numDimensions, d_curY, d_gk, d_gkm1, &minusOne, false); // y = gk - gkm1

		if (isRayleighRitz) {
			bool shouldNormalise = true;
			m_costFunction.orthogopt(d_x, shouldNormalise);
		}

		m_cublas.dispatchDot(numDimensions, d_tmp,  d_curY, d_curS); // tmp  = y Dot s
		m_cublas.dispatchDot(numDimensions, d_tmp2, d_curY, d_curY); // tmp2 = y Dot y

		// rho = 1 / tmp
		// H0 = tmp / tmp2
		update3<<<1, 1>>>(m_d_rho + (m_cumulativeIter % m_updates), m_d_H0, d_tmp, d_tmp2);
		CudaCheckError();
		cudaDeviceSynchronize();

		if (isTimingUpdates) {
			m_timer_updates.stop();
		}

		++m_cumulativeIter;
	}

	*itDone = it;

	if (isRayleighRitz) {
		// Normalise. 
		double thisFactor;
		m_cublas.dispatchNrm2(numDimensions, &thisFactor, d_x, false); // thisFactor = sqrt(x Dot x)
		thisFactor = 1.0/thisFactor;
		m_cublas.dispatchScale(numDimensions, d_x, d_x, &thisFactor, false); // x = thisfactor*x
	}

	if (printingOn) {
		double diagOut;
		CudaSafeCall( cudaMemcpy(&diagOut, m_d_H0, sizeof(double), cudaMemcpyDeviceToHost) );
		fileHandle << " Diagonal inverse Hessian elements = " << std::setw(20) << std::setprecision(10) << diagOut 
			<< std::endl;
		fileHandle << " Reason for termination: " << statusToString(lbfgsStatus) << std::endl;

		fileHandle << std::endl;
	}

	if (m_projectGrad) {
		CudaSafeCall( cudaMemcpy(d_gk, d_gkSave, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );
		m_costFunction.calculateRms(d_gkSave, outRms);
		m_costFunction.setDevCoords(d_x);
		m_costFunction.setDevEnergy(d_fk);
		m_costFunction.setDevGradient(d_gk);
	}

	// Deallocations. 
	CudaSafeCall( cudaFree(d_gkm1) );
	CudaSafeCall( cudaFree(d_gkSave) );
	CudaSafeCall( cudaFree(d_z) );

	if (isTimingTotal) {
		m_timer_total.stop();
	}

	return lbfgsStatus;
}



// Device / kernel functions
namespace gpu_lbfgs
{
	__global__ void update1(double *alphaOut, const double *sDotZ, const double *rho, double *minusAlphaOut)
	{
		*alphaOut      = *sDotZ * (*rho);
		*minusAlphaOut = -*alphaOut;
	}

	__global__ void update2(double *alphaMinusBetaOut, const double *rho, const double *yDotZ, const double *alpha)
	{
		const double beta = *rho * *yDotZ;
		*alphaMinusBetaOut = *alpha - beta;
	}

	__global__ void update3(double *rhoOut, double *H0Out, double *yDotS, double *yDotY)
	{
		if (abs(*yDotS) < 1.0e-20) {
			if (*yDotS >= 0.0) {
				*yDotS = 1.0e-20;
			}
			else {
				*yDotS = -1.0e-20;
			}
		}

		*rhoOut = 1.0 / *yDotS;

		if (abs(*yDotY) < 1.0e-20) {
			if (*yDotY >= 0.0) {
				*yDotY = 1.0e-20;
			}
			else {
				*yDotY = -1.0e-20;
			}
		}


		*H0Out = *yDotS / *yDotY;
	}

	__global__ void updateNegative(double *tmpOut)
	{
		*tmpOut = -*tmpOut;
	}

}
