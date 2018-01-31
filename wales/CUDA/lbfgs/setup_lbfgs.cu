/**
 *
 * File setup_lbfgs.cu: Communicates with GMIN/OPTIM and sets up and starts an LBFGS minmization. 
 *
 **/

#include <stdbool.h>
#include <iostream>

#include "lbfgs.h"
#include "potential.h"

// Can call this function with any potential class. 
	template <class P>
void setup(Printing &debugPrinting, Cublas &cublas, const size_t numDimensions, double *x, double *gradEps, _Bool *isConverged, 
		double *energy, int *maxIter, int *itDone, double *maxStep, double *maxFkRise, double *outRms, P &potential, 
		_Bool *printDebug, _Bool *timeCuda, int *nPotCalls, double *H0, int *updates, _Bool *isColdFusion, _Bool *projectGrad)
{
	Timer timer_total     ("GPU_LBFGS_total"     );
	Timer timer_updates   ("GPU_LBFGS_updates"   );
	Timer timer_linesearch("GPU_LBFGS_linesearch");

	if (*timeCuda) {
		timer_total.setTimingOn();
		timer_updates.setTimingOn();
		timer_linesearch.setTimingOn();
	}

	double *d_x;  // Coordinates on the GPU. 
	double *d_fk; // Energy on the GPU. 
	double *d_gk; // Gradient on the GPU. 

	// Status returned which indicates whether the minimization was successful (see lbfgs.h).
	Lbfgs::lStatus lbfgsStatus;

	// Lbfgs constructor. 
	Lbfgs lbfgs(potential, debugPrinting, timer_total, timer_updates, timer_linesearch, cublas, numDimensions, *H0, *updates, *maxIter, *gradEps, *maxStep, *maxFkRise, 10);
	lbfgs.setProjectGrad(*projectGrad);

	CudaSafeCall( cudaMalloc(&d_x, numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_fk, sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_gk, numDimensions * sizeof(double)) );

	CudaSafeCall( cudaMemcpy(d_x, x, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );

	// Perform the minimization. 
	lbfgsStatus = lbfgs.minimize(d_x, d_fk, d_gk, outRms, itDone);

	CudaSafeCall( cudaMemcpy(x, d_x, numDimensions * sizeof(double), cudaMemcpyDeviceToHost) );
	CudaSafeCall( cudaMemcpy(energy, d_fk, sizeof(double), cudaMemcpyDeviceToHost) );

	if (lbfgsStatus == Lbfgs::LBFGS_CONVERGED) {
		*isConverged = 1;
		*isColdFusion = 0;
	}
	else {
		*isConverged = 0;
		if (lbfgsStatus == Lbfgs::LBFGS_COLD_FUSION_DIAGNOSED) {
			*isColdFusion = 1;
		}
	}

	CudaSafeCall( cudaFree(d_x) );	
	CudaSafeCall( cudaFree(d_fk) );
	CudaSafeCall( cudaFree(d_gk) );

	*nPotCalls = potential.getnCalls();

	if (*timeCuda) {
		timer_total.saveMeasurement();
		timer_updates.saveMeasurement();
		timer_linesearch.saveMeasurement();
	}
}



extern "C" void setup_lbfgs(int *n, double *x, double *gradEps, _Bool *isConverged, double *energy, int *maxIter, int *itDone, 
		double *maxStep, double *maxFkRise, double *outRms, char *cudaPot, _Bool *printDebug, _Bool *timeCuda, 
		int *nPotCalls, _Bool *isColdFusion, double *coldFusionLim, double *H0, int *updates, 
		_Bool *isAtomisticNotRigid, int *nDegFreedomF, int *nRigidBodyF, int *nRigidSitesPerBodyF, 
		int *rigidGroupsF, int *rigidMaxSiteF, double *sitesRigidBodyF, int *rigidSinglesF, 
		double *aaConvThresholdF, double *rigidInverseF, _Bool *projectGrad, _Bool *shouldFreeze, 
		_Bool *isAtomFrozen, int *nFreeze, double *potTimeElapsed, bool *isAaConvergence, int *nAddTarget, double *ljAddRep, 
		double *ljAddAtt)
{
	const size_t numDimensions = *n;

	const int nSecDiag = 0;

	const double aaConvThreshold = 5*(*aaConvThresholdF);

	// Rigid body parameters may not be intialised if rigid body framework not being used.
	int nDegFreedom  = 0;
	int nRigidBody   = 0;
	int rigidMaxSite = 0;

	int *nRigidSitesPerBody = NULL; 
	int *rigidGroups        = NULL; 
	double *sitesRigidBody  = NULL; 
	int *rigidSingles       = NULL; 
	double *rigidInverse    = NULL;

	double *coords = NULL;

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

	// This is first opened in GMIN main to ensure file contents overwritten only at beginning of new run.
	Printing debugPrinting("GPU_debug_out");
	if (*printDebug) {
		debugPrinting.setPrintingOn();
	}

	Timer timer_potential     ("GPU_potential");
	timer_potential.setTimingOn();

	// Set up cuBLAS. 
	Cublas cublas;

	// L specifies the Lennard-Jones potential.
	if (*cudaPot == 'L') {
		// Create an instance of the appropriate class for the potential, LjPotential.
		LjPotential potential(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, 
				rigidMaxSite, nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, 
				coords, nSecDiag, *isAtomisticNotRigid, aaConvThreshold, *coldFusionLim, *shouldFreeze, 
				isAtomFrozen, *nFreeze, *isAaConvergence, *nAddTarget, ljAddRep, ljAddAtt);

		setup<LjPotential>(debugPrinting, cublas, numDimensions, x, gradEps, isConverged, energy, maxIter, itDone, 
				maxStep, maxFkRise, outRms, potential, printDebug, timeCuda, nPotCalls, H0, updates, isColdFusion, 
				projectGrad);
	}
	// A specifies the AMBER potential.
	else if (*cudaPot == 'A') {
		AmberPotential potential(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, 
				rigidMaxSite, nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, 
				coords, nSecDiag, *isAtomisticNotRigid, aaConvThreshold, *coldFusionLim, *shouldFreeze, 
				isAtomFrozen, *nFreeze, *isAaConvergence, *nAddTarget, ljAddRep, ljAddAtt);

		setup<AmberPotential>(debugPrinting, cublas, numDimensions, x, gradEps, isConverged, energy, maxIter, itDone, 
				maxStep, maxFkRise, outRms, potential, printDebug, timeCuda, nPotCalls, H0, updates, isColdFusion, 
				projectGrad); 
	}
	else {
		std::cerr << "The specified potential has not been recognised" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (*timeCuda) {
		timer_potential.saveMeasurement();
	}

	*potTimeElapsed = timer_potential.elapsed()/1000.0;
}
