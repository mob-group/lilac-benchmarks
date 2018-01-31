/**
 *
 * File setup_bfgsts.cu: Communicates with OPTIM and sets up a BFGSTS calculation. 
 *
 **/

#include <stdbool.h>
#include <iostream>

#include "bfgsts.h"
#include "potential.h"

// Can call this function with any potential class. 
	template <class P>
void setup(Printing &debugPrinting, Cublas &cublas, const size_t numDimensions, double *xGradEps, _Bool *isConverged, 
		double *energy, int *xMaxIter, int *bIter, double *xMaxStep, double *xMaxFkRise, double *outRms, P &potential, 
		_Bool *printDebug, _Bool *timeCuda, int *nPotCalls, _Bool *isColdFusion, double *coldFusionLim, double *xH0, 
		int *xUpdates, double *evec, double *evalPercentEps, double *evalMin, double *coords, double *pushOffCut, 
		double *pushOffMag, double *maxEvecStep, double *maxMaxStep, double *minMaxStep, double *trustRadius, int *pMaxIter1, 
		int *pMaxIter2, int *bItMax, double *evecOverlapTol, int *pUpdates, double *gradEps, double *pMaxStep, 
		double *pMaxFkRise, double *pH0, double *evStepMagTol)
{
	bool isBfgsts = true;
	potential.setIsBfgsts(isBfgsts);

	bool printingOn = debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = debugPrinting.getFileHandle();

	fileHandle << "CUDA BFGSTS" << std::endl;

	Timer timer_xLbfgsTotal     ("GPU_BFGSTS_Rayleigh_Ritz_total"     );
	// Timing always turned off for updates and linesearch - setTimingOn() to turn on. 
	Timer timer_xLbfgsUpdates   ("GPU_BFGSTS_Rayleigh_Ritz_updates"   );
	Timer timer_xLbfgsLinesearch("GPU_BFGSTS_Rayleigh_Ritz_linesearch");

	Timer timer_pLbfgsTotal     ("GPU_BFGSTS_subspace_min_total"      );
	// Timing always turned off for updates and linesearch - setTimingOn() to turn on. 
	Timer timer_pLbfgsUpdates   ("GPU_BFGSTS_subspace_min_updates"    );
	Timer timer_pLbfgsLinesearch("GPU_BFGSTS_subspace_min_linesearch" );

	Timer timer_total           ("GPU_BFGSTS_total"   );
	Timer timer_eFolSteps       ("GPU_BFGSTS_EF_steps");

	if (*timeCuda) {
		timer_xLbfgsTotal.setTimingOn();
		timer_pLbfgsTotal.setTimingOn();
		timer_total.setTimingOn();
		timer_eFolSteps.setTimingOn();
	}

	int xnRedMax = 3;
	// Lbfgs constructor. 
	Lbfgs xLbfgs(potential, debugPrinting, timer_xLbfgsTotal, timer_xLbfgsUpdates, timer_xLbfgsLinesearch, cublas, 
			numDimensions, *xH0, *xUpdates, *xMaxIter, *xGradEps, *xMaxStep, *xMaxFkRise, xnRedMax);
	xLbfgs.setEvalPercentEps(*evalPercentEps);

	int pnRedMax = 10;
	// Second Lbfgs constructor. 
	Lbfgs pLbfgs(potential, debugPrinting, timer_pLbfgsTotal, timer_pLbfgsUpdates, timer_pLbfgsLinesearch, cublas, 
			numDimensions, *pH0, *pUpdates, *pMaxIter1, *gradEps, *pMaxStep, *pMaxFkRise, pnRedMax);
	pLbfgs.setProjectGrad(true);

	// Bfgsts constructor. 
	Bfgsts bfgsts(potential, xLbfgs, pLbfgs, debugPrinting, timer_total, timer_eFolSteps, cublas, numDimensions, *pushOffCut, 
			*pushOffMag, *maxEvecStep, *maxMaxStep, *minMaxStep, *trustRadius, *pMaxIter1, *pMaxIter2, *bItMax, 
			*evecOverlapTol, *evStepMagTol);

	// Status returned which indicates whether the search was successful (see bfgsts.h). 
	Bfgsts::bStatus bfgstsStatus;

	double *d_evec; // Eigenvector on the GPU.
	double *d_coords; // Coordinates on the GPU.

	CudaSafeCall( cudaMalloc(&d_evec,   numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_coords, numDimensions * sizeof(double)) );

	CudaSafeCall( cudaMemcpy(d_evec,   evec,   numDimensions * sizeof(double), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(d_coords, coords, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );

	// Main BFGSTS calculation. 
	bfgstsStatus = bfgsts.findTs(d_coords, d_evec, energy, evalMin, outRms, bIter);

	if (bfgstsStatus == Bfgsts::BFGSTS_CONVERGED_TO_TS) {
		*isConverged = 1;
		*isColdFusion = 0;
	}
	else {
		*isConverged = 0;
		if (bfgstsStatus == Bfgsts::BFGSTS_COLD_FUSION_DIAGNOSED) {
			*isColdFusion = 1;
		}
	}

	CudaSafeCall( cudaMemcpy(evec,   d_evec,   numDimensions * sizeof(double), cudaMemcpyDeviceToHost) );
	CudaSafeCall( cudaMemcpy(coords, d_coords, numDimensions * sizeof(double), cudaMemcpyDeviceToHost) );

	CudaSafeCall( cudaFree(d_evec) );
	CudaSafeCall( cudaFree(d_coords) );

	*nPotCalls = potential.getnCalls();

	if (*timeCuda) {
		timer_xLbfgsTotal.saveMeasurement();
		timer_pLbfgsTotal.saveMeasurement();
		timer_total.saveMeasurement();
		timer_eFolSteps.saveMeasurement();
	}
}



extern "C" void setup_bfgsts(int *nAtoms, double *coords, double *xGradEps, _Bool *isConverged, double *energy, int *xMaxIter, 
		int *bItMax, int *bIter, double *xMaxStep, double *xMaxFkRise, double *outRms, char *cudaPotential, 
		_Bool *printDebug, _Bool *timeCuda, int *nPotCalls, _Bool *isColdFusion, double *coldFusionLim, double *xH0, 
		int *xUpdates, double *evalMin, double *evec, _Bool *isAtomisticNotRigid, int *nDegFreedomF, int *nRigidBodyF, 
		int *nRigidSitesPerBodyF, int *rigidGroupsF, int *rigidMaxSiteF, double *sitesRigidBodyF, int *rigidSinglesF, 
		double *rigidInverseF, int *nSecDiag, double *evalPercentEps, double *pushOffCut, double *pushOffMag, 
		double *maxEvecStep, double *maxMaxStep, double *minMaxStep, double *trustRadius, int *pMaxIter1, int *pMaxIter2, 
		double *evecOverlapTol, int *pUpdates, double *gradEps, double *pMaxStep, double *pMaxFkRise, double *pH0, 
		double *evStepMagTol, _Bool *shouldFreeze, _Bool *isAtomFrozen, int *nFreeze, double *potTimeElapsed, 
		bool *isAaConvergence, int *nAddTarget, double *ljAddRep, double *ljAddAtt)
{
	const size_t numDimensions = 3*(*nAtoms);
	const double aaConvThreshold = 5*(*gradEps);

	// Rigid body parameters may not be intialised in Fortran code if rigid body framework not being used. 
	int nDegFreedom  = 0;
	int nRigidBody   = 0;
	int rigidMaxSite = 0;

	int *nRigidSitesPerBody = NULL;
	int *rigidGroups        = NULL;
	double *sitesRigidBody  = NULL;
	int *rigidSingles       = NULL;
	double *rigidInverse    = NULL;

	if (!(*isAtomisticNotRigid)) {
		nDegFreedom  = *nDegFreedomF;
		nRigidBody   = *nRigidBodyF;
		rigidMaxSite = *rigidMaxSiteF;

		nRigidSitesPerBody = nRigidSitesPerBodyF;
		rigidGroups        = rigidGroupsF;
		sitesRigidBody     = sitesRigidBodyF;
		rigidSingles       = rigidSinglesF;
		rigidInverse       = rigidInverseF;
	}

	// This is first opened in OPTIM fetchz to ensure file contents overwritten only at beginning of new run. 
	Printing debugPrinting("GPU_debug_out");
	if (*printDebug) {
		debugPrinting.setPrintingOn();
	}

	Timer timer_potential     ("GPU_potential");
	timer_potential.setTimingOn();

	// Set up cuBLAS. 
	Cublas cublas;

	// L specifies the Lennard-Jones potential. 
	if (*cudaPotential == 'L') {
		// Create an instance of the appropriate class for the potential, LjPotential. 
		LjPotential potential(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, rigidMaxSite, 
				nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, coords, 
				*nSecDiag, *isAtomisticNotRigid, aaConvThreshold, *coldFusionLim, *shouldFreeze, isAtomFrozen, 
				*nFreeze, *isAaConvergence, *nAddTarget, ljAddRep, ljAddAtt);

		setup<LjPotential>(debugPrinting, cublas, numDimensions, xGradEps, isConverged, energy, xMaxIter, bIter, xMaxStep, xMaxFkRise, 
				outRms, potential, printDebug, timeCuda, nPotCalls, isColdFusion, coldFusionLim, xH0, xUpdates, 
				evec, evalPercentEps, evalMin, coords, pushOffCut, pushOffMag, maxEvecStep, maxMaxStep, 
				minMaxStep, trustRadius, pMaxIter1, pMaxIter2, bItMax, evecOverlapTol, pUpdates, gradEps, 
				pMaxStep, pMaxFkRise, pH0, evStepMagTol);
	}
	// A specifies the AMBER potential. 
	else if (*cudaPotential == 'A') {
		AmberPotential potential(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, rigidMaxSite, 
				nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, coords, 
				*nSecDiag, *isAtomisticNotRigid, aaConvThreshold, *coldFusionLim, *shouldFreeze, isAtomFrozen, 
				*nFreeze, *isAaConvergence, *nAddTarget, ljAddRep, ljAddAtt);

		setup<AmberPotential>(debugPrinting, cublas, numDimensions, xGradEps, isConverged, energy, xMaxIter, bIter, xMaxStep, xMaxFkRise, 
				outRms, potential, printDebug, timeCuda, nPotCalls, isColdFusion, coldFusionLim, xH0, 
				xUpdates, evec, evalPercentEps, evalMin, coords, pushOffCut, pushOffMag, maxEvecStep, 
				maxMaxStep, minMaxStep, trustRadius, pMaxIter1, pMaxIter2, bItMax, evecOverlapTol, pUpdates, 
				gradEps, pMaxStep, pMaxFkRise, pH0, evStepMagTol);
	}
	else {
		std::cerr << "The specified potential has not been recognised. " << std::endl;
		exit(EXIT_FAILURE);
	}

	if (*timeCuda) {
		timer_potential.saveMeasurement();
	}

	*potTimeElapsed = timer_potential.elapsed()/1000.0;

}

