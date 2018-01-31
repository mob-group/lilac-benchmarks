/**
 *
 * File amber_potential.cu: Implementation of class AmberPotential, inherited from class CostFunction. 
 *
 **/

#include "cost_function.h"
#include "potential.h"
#include "gpu.h"

AmberPotential::AmberPotential(Printing &debugPrinting, Timer &timer_potential, Cublas &cublas, 
		size_t numDimensions, int nDegFreedom, int nRigidBody, int rigidMaxSite, 
		int *nRigidSitesPerBody, int *rigidGroups, double *sitesRigidBody, int *rigidSingles, 
		double *rigidInverse, double *x, int nSecDiag, bool isAtomisticNotRigid, 
		double aaConvThreshold, double coldFusionLim, bool shouldFreeze, bool *isAtomFrozen, 
		int nFreeze, bool isAaConvergence, int nAddTarget, double *ljAddRep, double *ljAddAtt)
	: CostFunction(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, rigidMaxSite, 
			nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, x, nSecDiag, 
			isAtomisticNotRigid, aaConvThreshold, coldFusionLim, shouldFreeze, isAtomFrozen, nFreeze, 
			isAaConvergence, nAddTarget, ljAddRep, ljAddAtt) {}



void AmberPotential::computeEnergyAndGradient(const double *d_x, double *d_f, double *d_gradf)
{
	double testEnergy;

	// Copy coordinates to correct location in device memory. 
	gminoptim_copy_crd(d_x);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Energy calculation. 
	gminoptim_gpu_energy(d_f);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Copy forces back. 
	gminoptim_copy_frc(d_gradf);
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaMemcpy(&testEnergy, d_f, sizeof(double), cudaMemcpyDeviceToHost) );
	// Test for cold fusion. 
	if (testEnergy < m_coldFusionLim) {
		m_isColdFusion = true;
	}
}
