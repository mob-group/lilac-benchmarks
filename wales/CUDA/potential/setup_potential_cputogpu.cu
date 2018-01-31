/**
 *
 * File setup_potential_cputogpu.cu: Communicates with GMIN/OPTIM and calls the potential. 
 *
 **/

#include <stdbool.h>
#include <iostream>

#include "cost_function.h"
#include "potential.h"

// Can call this function with any potential class. 
	template <class P>
void setup(const size_t numDimensions, double *x, double *energy, double *gradients, P &potential)
{
	double *d_x;  // Coordinates on the GPU.
	double *d_fk; // Energy on the GPU.
	double *d_gk; // Gradient on the GPU.

	CudaSafeCall( cudaMalloc(&d_x, numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_fk, sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_gk, numDimensions * sizeof(double)) );

	CudaSafeCall( cudaMemcpy(d_x, x, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );

	// Perform the potential calculation.
	potential.computeEnergyAndGradient(d_x, d_fk, d_gk);

	CudaSafeCall( cudaMemcpy(energy, d_fk, sizeof(double), cudaMemcpyDeviceToHost) );
	CudaSafeCall( cudaMemcpy(gradients, d_gk, numDimensions * sizeof(double), cudaMemcpyDeviceToHost) );

	CudaSafeCall( cudaFree(d_x) );
	CudaSafeCall( cudaFree(d_fk) );
	CudaSafeCall( cudaFree(d_gk) );
}



extern "C" void setup_potential_cputogpu(int *n, double *x, double *energy, double *gradients, int *nAddTarget, double *ljAddRep, 
		double *ljAddAtt, char *cudaPot, _Bool *timeCuda, double *potTimeElapsed)
{
	const size_t numDimensions = *n;

	const int nSecDiag = 0;

	int nFreeze = 0;

	const double aaConvThreshold = 0;

	double coldFusionLim = -1.0e6;

	bool shouldFreeze = false;
	bool isAaConvergence = false;

	bool *isAtomFrozen = NULL;

	// Rigid body framework not being used.
	bool isAtomisticNotRigid = true;

	int nDegFreedom  = 0;
	int nRigidBody   = 0;
	int rigidMaxSite = 0;

	int *nRigidSitesPerBody = NULL; 
	int *rigidGroups        = NULL; 
	double *sitesRigidBody  = NULL; 
	int *rigidSingles       = NULL; 
	double *rigidInverse    = NULL;

	double *coords = NULL;

	// This is first opened in GMIN main to ensure file contents overwritten only at beginning of new run.
	Printing debugPrinting("GPU_debug_out");
	bool printDebug = false;
	if (printDebug) {
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
				coords, nSecDiag, isAtomisticNotRigid, aaConvThreshold, coldFusionLim, shouldFreeze, 
				isAtomFrozen, nFreeze, isAaConvergence, *nAddTarget, ljAddRep, ljAddAtt);

		setup<LjPotential>(numDimensions, x, energy, gradients, potential);
	}
	// A specifies the AMBER potential.
	else if (*cudaPot == 'A') {
		AmberPotential potential(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, 
				rigidMaxSite, nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, 
				coords, nSecDiag, isAtomisticNotRigid, aaConvThreshold, coldFusionLim, shouldFreeze, 
				isAtomFrozen, nFreeze, isAaConvergence, *nAddTarget, ljAddRep, ljAddAtt);

		setup<AmberPotential>(numDimensions, x, energy, gradients, potential);
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
