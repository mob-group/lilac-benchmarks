/**
 *
 * File bfgsts.cu: Implementation of class Bfgsts. 
 *
 **/

#include <iomanip>

#include "bfgsts.h"

Bfgsts::Bfgsts(CostFunction &costFunc, Lbfgs &xLbfgs, Lbfgs &pLbfgs, Printing &debugPrinting, Timer &timer_total, 
		Timer &timer_eFolSteps, Cublas &cublas, size_t numDimensions, double pushOffCut, double pushOffMag, 
		double maxEvecStep, double maxMaxStep, double minMaxStep, double trustRadius, int pMaxIter1, int pMaxIter2, 
		int bItMax, double evecOverlapTol, double evStepMagTol)
	: m_costFunction(costFunc)
	, m_xLbfgs(xLbfgs)
	, m_pLbfgs(pLbfgs)
	, m_debugPrinting(debugPrinting)
	, m_cublas(cublas)
	, m_timer_total(timer_total)
	, m_timer_eFolSteps(timer_eFolSteps)
	, m_pushOffCut(pushOffCut)
	, m_pushOffMag(pushOffMag)
	, m_maxEvecStep(maxEvecStep)
	, m_maxMaxStep(maxMaxStep)
	, m_minMaxStep(minMaxStep)
	, m_trustRadius(trustRadius)
	, m_pMaxIter1(pMaxIter1)
	, m_pMaxIter2(pMaxIter2)
	, m_bItMax(bItMax)
	  , m_evecOverlapTol(evecOverlapTol)
	  , m_evStepMagTol(evStepMagTol) {}



std::string Bfgsts::statusToString(Bfgsts::bStatus bfgstsStatus)
{
	switch (bfgstsStatus) {
		case BFGSTS_CONVERGED_TO_TS :
			return "Converged to a transition state of Hessian index 1";
		case BFGSTS_REACHED_MAX_ITER :
			return "Reached maximum number of iterations";
		case BFGSTS_COLD_FUSION_DIAGNOSED :
			return "Cold fusion diagnosed";
		default :
			return "Unknown status";
	}
}



Bfgsts::bStatus Bfgsts::findTs(double *d_coords, double *d_evec, double *energy, double *evalMin, double *outRms, int *itDone)
{
	bool printingOn = m_debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = m_debugPrinting.getFileHandle();

	bool isTimingTotal = m_timer_total.getTimingOn();
	bool isTimingSteps = m_timer_eFolSteps.getTimingOn();

	if (isTimingTotal) {
		m_timer_total.start();
	}

	const size_t numDimensions = m_costFunction.getNumDimensions();

	// Status returned by Lbfgs which indicates whether the minimization was successful. 
	Lbfgs::lStatus xLbfgsStatus, pLbfgsStatus;

	double *d_evecOld; // Copy of d_evec for calculation of overlap.
	double *d_energy;
	double *d_evalMin;
	double *d_gradient;
	double *d_xGrad;

	CudaSafeCall( cudaMalloc(&d_evecOld, numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_energy, sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_evalMin, sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_gradient, numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_xGrad, numDimensions * sizeof(double)) );

	bStatus bfgstsStatus = BFGSTS_REACHED_MAX_ITER;

	double *xRms = new double;
	int *xItDone = new int;
	int *pItDone = new int;

	double efStepMag   = 0.0;
	double evecOverlap = 0.0;
	bool isNegativeEv  = false;

	m_costFunction.setIsRayleighRitz(false);
	m_costFunction.basePotential(d_coords, d_energy, d_gradient);

	m_costFunction.calculateRms(d_gradient, outRms);

	m_costFunction.setDevEnergy(d_energy);
	m_costFunction.setDevGradient(d_gradient);

	m_costFunction.freeze(d_evec);

	*xItDone = 0;
	*pItDone = 0;

	size_t it = 1;
	for (it = 1; it <= m_bItMax; ++it) {

		fileHandle << std::endl;
		fileHandle << " Beginning of optimization cycle " << it << "." << std::endl;
		fileHandle << " -------------------------------" << std::endl;

		bool isFirstIter = (it == 1);

		CudaSafeCall( cudaMemcpy(d_evecOld, d_evec, numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );

		// True if we are minimizing Rayleigh-Ritz ratio and want to calculate secdiag. 
		m_costFunction.setIsRayleighRitz(true);
		fileHandle << " Lowest eigenvector search (Rayleigh-Ritz):" << std::endl;
		// Perform the minimization. 
		xLbfgsStatus = m_xLbfgs.minimize(d_evec, d_evalMin, d_xGrad, xRms, xItDone);

		CudaSafeCall( cudaMemcpy(evalMin, d_evalMin, sizeof(double), cudaMemcpyDeviceToHost) );

		bool isEvTooSmall = (abs(*evalMin) < 1.0e-100);
		if (isEvTooSmall) {
			if (*evalMin >= 0) {
				*evalMin = 1.0e-20;
			}
			else {
				*evalMin = -1.0e-20;
			}
		}

		isNegativeEv = (*evalMin < 0.0);

		if (isFirstIter && !isNegativeEv) {
			fileHandle << " Warning: initial eigenvalue is positive - increase NEBK?" << std::endl;
		}

		bool isxLbfgsConverged = (xLbfgsStatus == Lbfgs::LBFGS_CONVERGED);

		m_cublas.dispatchDot(numDimensions, &evecOverlap, d_evec, d_evecOld, false); // evecOverlap = evec Dot evecOld
		if (isxLbfgsConverged) {
			fileHandle << " Smallest eigenvalue converged in " << std::setw(7) << *xItDone << " steps" 
				<< std::endl; 
		}
		else {
			fileHandle << " Warning: smallest eigenvalue did not converge" << std::endl;
		}

		m_costFunction.freeze(d_evec);

		if (isTimingSteps) {
			m_timer_eFolSteps.start();
		}

		m_pLbfgs.setDevzWork(d_evec, numDimensions);

		// Step uphill parallel to the lowest eigenvector. 
		stepUphill(d_coords, d_evec, evalMin, d_energy, d_gradient, isxLbfgsConverged, isNegativeEv, efStepMag, outRms, it);

		if (isTimingSteps) {
			m_timer_eFolSteps.stop();
		}

		bool isOverlap = ((1.0 - abs(evecOverlap)) <= m_evecOverlapTol);
		if (!isFirstIter && isOverlap && isNegativeEv) {
			m_pLbfgs.setMaxIterations(m_pMaxIter2);
		}
		else {
			m_pLbfgs.setMaxIterations(m_pMaxIter1);
		}

		m_costFunction.setIsRayleighRitz(false);
		fileHandle << " Minimization perpendicular to lowest eigenvector:" << std::endl;
		bool shouldReset = false;
		// Perform minimization, projecting out uphill components of gradient along eigenvector.
		pLbfgsStatus = m_pLbfgs.minimize(d_coords, d_energy, d_gradient, outRms, pItDone, shouldReset);

		fileHandle << " True RMS grad = " << std::setw(15) << std::setprecision(7) << *outRms << std::endl;

		bool ispLbfgsConverged  = (pLbfgsStatus == Lbfgs::LBFGS_CONVERGED);
		bool ispLbfgsColdFusion = (pLbfgsStatus == Lbfgs::LBFGS_COLD_FUSION_DIAGNOSED);
		bool isEfStepConverged  = (efStepMag <= m_evStepMagTol);
		if (ispLbfgsConverged && isNegativeEv && isEfStepConverged) {
			bfgstsStatus = BFGSTS_CONVERGED_TO_TS;
			break;
		}
		else if (ispLbfgsColdFusion) {
			bfgstsStatus = BFGSTS_COLD_FUSION_DIAGNOSED;
			break;
		}
	}

	*itDone = it;

	CudaSafeCall( cudaMemcpy(energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost) );

	fileHandle << std::endl;
	fileHandle << " Reason for termination: " << statusToString(bfgstsStatus) << std::endl;
	fileHandle << std::endl;

	delete xItDone;
	delete pItDone;
	delete xRms;

	CudaSafeCall( cudaFree(d_evecOld) );
	CudaSafeCall( cudaFree(d_energy) );
	CudaSafeCall( cudaFree(d_evalMin) );
	CudaSafeCall( cudaFree(d_gradient) );
	CudaSafeCall( cudaFree(d_xGrad) );

	if (isTimingTotal) {
		m_timer_total.stop();
	}

	return bfgstsStatus;
}



void Bfgsts::stepUphill(double *d_coords, const double *d_evec, const double *evalMin, double *d_energy, double *d_gradient, 
		const bool isxLbfgsConverged, const bool &isNegativeEv, double &efStepMag, double *outRms, size_t it)
{
	bool printingOn = m_debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = m_debugPrinting.getFileHandle();

	const size_t numDimensions = m_costFunction.getNumDimensions();

	double fob;
	double efStep;
	double trustRatio;

	m_cublas.dispatchDot(numDimensions, &fob, d_evec, d_gradient, false); // fob = evec Dot gradient

	bool shouldStepAway = false;
	bool couldPushoff = (*outRms < m_pushOffCut);
	if (couldPushoff && !isNegativeEv && isxLbfgsConverged) {
		bool isFirstIter = (it == 1);
		bool isIterMultipleOfFour = ((it - 1) % 4 == 0);
		if (isFirstIter) {
			fileHandle << " Stepping away from minimum along softest mode + direction" << std::endl;
			shouldStepAway = true;
		}
		else if (isIterMultipleOfFour) {
			fileHandle << " Stepping away from solution of wrong index" << std::endl;
			if (m_pushOffMag != 0.0) {
				fob = m_pushOffMag;
			}
			else {
				fob = m_maxEvecStep;
			}
		}
	}

	double xp1 = abs(*evalMin)/2.0;
	double xp2 = 1.0 + 4.0*(fob/(*evalMin))*(fob/(*evalMin));
	xp1 = -xp1*(1.0 + sqrt(xp2));
	efStep = -fob/xp1;

	if (shouldStepAway) {
		if (m_pushOffMag != 0.0) {
			efStep = m_pushOffMag;
		}
		else {
			efStep = m_maxEvecStep/10.0;
		}
	}

	efStepMag = abs(efStep);

	if (!shouldStepAway) {
		double scale = std::min(m_maxEvecStep/std::max(abs(efStep),1.0e-10),1.0);
		efStep = scale*efStep;
	}

	if (isxLbfgsConverged && !isNegativeEv) {
		efStep = fob*m_maxMaxStep/abs(fob);
	}

	fileHandle << " EF step uphill of magnitude " << efStepMag << " is taken" << std::endl;
	// Take step uphill along the lowest eigenvector. 
	m_cublas.dispatchAxpy(numDimensions, d_coords, d_coords, d_evec, &efStep, false); // coords += evec*efStep

	double eOld;
	CudaSafeCall( cudaMemcpy(&eOld, d_energy, sizeof(double), cudaMemcpyDeviceToHost) );

	m_costFunction.setIsRayleighRitz(false);
	m_costFunction.basePotential(d_coords, d_energy, d_gradient);

	m_costFunction.calculateRms(d_gradient, outRms);

	double fobNew;
	m_cublas.dispatchDot(numDimensions, &fobNew, d_evec, d_gradient, false); // fobNew = evec Dot gradient

	bool isZeroEnergy = (eOld == 0.0);
	if (!isZeroEnergy && isNegativeEv) {
		double predictedEv = fobNew - fob;
		double actualEv = (efStep*(*evalMin));
		trustRatio = abs(1.0 - predictedEv/actualEv) ;
		if (trustRatio > m_trustRadius) {
			m_maxEvecStep = std::max(m_maxEvecStep/1.1, m_minMaxStep);
		}
		else {
			m_maxEvecStep = std::min(m_maxEvecStep*1.1, m_maxMaxStep);
		}
	}

	fileHandle << std::endl;
	fileHandle << " ----------------------------------------------------------------------------------------" 
		<< std::endl;
	fileHandle << "  Vector        Gradient          Secder            Step        Max step     Trust ratio " 
		<< std::endl;
	fileHandle << " ----------------------------------------------------------------------------------------" 
		<< std::endl;
	fileHandle << "     1   " << std::setw(15) << std::setprecision(6) << fob << " " 
		<< std::setw(15) << std::setprecision(6) << *evalMin << " " 
		<< std::setw(15) << std::setprecision(6) << efStep << " " 
		<< std::setw(15) << std::setprecision(6) << m_maxEvecStep << " " 
		<< std::setw(15) << std::setprecision(6) << trustRatio << " " << std::endl;
	fileHandle << " ----------------------------------------------------------------------------------------" 
		<< std::endl;
	fileHandle << std::endl;
}
