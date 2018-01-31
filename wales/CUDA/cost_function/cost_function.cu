/**
 *
 * File cost_function.cu: Implementation of class CostFunction, not including rigid bodies.
 *
 **/

#include<iomanip>

#include "cost_function.h"

CostFunction::CostFunction(Printing &debugPrinting, Timer &timer_potential, Cublas &cublas, size_t numDimensions, int nDegFreedom, 
		int nRigidBody, int rigidMaxSite, int *nRigidSitesPerBody, int *rigidGroups, double *sitesRigidBody, 
		int *rigidSingles, double *rigidInverse, double *coords, int nSecDiag, bool isAtomisticNotRigid, 
		double aaConvThreshold, double coldFusionLim, bool shouldFreeze, bool *isAtomFrozen, int nFreeze, 
		bool isAaConvergence, int nAddTarget, double *ljAddRep, double *ljAddAtt)
	: m_timer_potential(timer_potential)
	, m_numDimensions(numDimensions)
	, m_nDegFreedom(nDegFreedom)
	, m_nRigidBody(nRigidBody)
	, m_rigidMaxSite(rigidMaxSite)
	, m_nRigidSitesPerBody(nRigidSitesPerBody)
	, m_rigidGroups(rigidGroups)
	, m_sitesRigidBody(sitesRigidBody)
	, m_rigidSingles(rigidSingles)
	, m_rigidInverse(rigidInverse)
	, m_isColdFusion(false)
	, m_isRayleighRitz(false)
	, m_coords(coords)
	, m_d_nRigidSitesPerBody(0)
	, m_d_sitesRigidBody(0)
	, m_d_rigidGroups(0)
	, m_d_rigidSingles(0) 
	, m_d_rigidInverse(0)
	, m_d_xRigid(0)
	, m_d_nRigidBody(0)
	, m_d_rigidMaxSite(0)
	, m_d_coords(0)
	, m_nSecDiag(nSecDiag)
	, m_isAtomisticNotRigid(isAtomisticNotRigid)
	, m_aaConvThreshold(aaConvThreshold)
	, m_nCalls(0)
	, m_coldFusionLim(coldFusionLim)
	, m_shouldFreeze(shouldFreeze)
	, m_isAtomFrozen(isAtomFrozen)
	, m_nFreeze(nFreeze)
	, m_nLargeRigidBody(0)
	, m_debugPrinting(debugPrinting)
	, m_cublas(cublas)
	, m_isBfgsts(false)
	, m_isAaConvergence(isAaConvergence)
	, m_nAddTarget(nAddTarget)
	, m_ljAddRep(ljAddRep)
	  , m_ljAddAtt(ljAddAtt)
{
	for (int i = 0; i < m_nRigidBody; ++i) {
		if (m_nRigidSitesPerBody[i] > 32) {
			m_nLargeRigidBody += 1;
		}
	}

	m_largeRigidIndices = new int[m_nLargeRigidBody];

	double *zerosx;
	double *zerosy;
	double *zerosz;

	zerosx = new double[numDimensions];
	zerosy = new double[numDimensions];
	zerosz = new double[numDimensions];

	CudaSafeCall( cudaMalloc(&m_d_zerosx, m_numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_zerosy, m_numDimensions * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_zerosz, m_numDimensions * sizeof(double)) );

	if (m_rigidInverse != NULL) {
		CudaSafeCall( cudaMalloc(&m_d_gkAtoms, numDimensions * sizeof(double)) );
	}

	CudaSafeCall( cudaMalloc(&m_d_nRigidSitesPerBody, m_nRigidBody * sizeof(int)) );
	CudaSafeCall( cudaMalloc(&m_d_sitesRigidBody, m_nRigidBody * 3 * m_rigidMaxSite * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_rigidGroups, m_nRigidBody * m_rigidMaxSite * sizeof(int)) );
	CudaSafeCall( cudaMalloc(&m_d_rigidSingles, (m_nDegFreedom/3 - 2*m_nRigidBody) * sizeof(int)) );
	CudaSafeCall( cudaMalloc(&m_d_rigidInverse, m_nRigidBody * 3 * 3 * sizeof(double)) );

	CudaSafeCall( cudaMalloc(&m_d_nRigidBody,   sizeof(int)) );
	CudaSafeCall( cudaMalloc(&m_d_rigidMaxSite, sizeof(int)) );

	CudaSafeCall( cudaMalloc(&m_d_xRigid,  m_nDegFreedom * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_gkRigid, m_nDegFreedom * sizeof(double)) );

	if (m_coords != NULL) {
		CudaSafeCall( cudaMalloc(&m_d_coords,   numDimensions * sizeof(double)) );
		CudaSafeCall( cudaMalloc(&m_d_energy, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&m_d_gradient, numDimensions * sizeof(double)) );
	}

	if (m_isAtomFrozen != NULL) {
		CudaSafeCall( cudaMalloc(&m_d_isAtomFrozen, numDimensions/3 * sizeof(bool)) );
	}

	CudaSafeCall( cudaMalloc(&m_d_nLargeRigidBody, sizeof(int)) );
	CudaSafeCall( cudaMalloc(&m_d_largeRigidIndices, m_nLargeRigidBody * sizeof(int)) );

	CudaSafeCall( cudaMalloc(&m_d_ljAddRep, m_nAddTarget * m_nAddTarget * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&m_d_ljAddAtt, m_nAddTarget * m_nAddTarget * sizeof(double)) );

	int counter = 0;
	while (counter < m_nLargeRigidBody) {
		for (int i = 0; i < m_nRigidBody; ++i) {
			if (m_nRigidSitesPerBody[i] > 32) {
				m_largeRigidIndices[counter] = i;
				++counter;
			}
		}
	}

	for (size_t i = 1; i < (numDimensions + 1); ++i) {
		if ((i + 2)%3 == 0) {
			zerosx[i-1] = 1.0;
		}
		else {
			zerosx[i-1] = 0.0;
		}
	}

	for (size_t i = 1; i < (numDimensions + 1); ++i) {
		if ((i + 1)%3 == 0) {
			zerosy[i-1] = 1.0;
		}
		else {
			zerosy[i-1] = 0.0;
		}
	}

	for (size_t i = 1; i < (numDimensions + 1); ++i) {
		if (i%3 == 0) {
			zerosz[i-1] = 1.0;
		}
		else {
			zerosz[i-1] = 0.0;
		}
	}

	CudaSafeCall( cudaMemcpy(m_d_nRigidSitesPerBody, m_nRigidSitesPerBody, m_nRigidBody * sizeof(int), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_sitesRigidBody, m_sitesRigidBody, m_nRigidBody * 3 * m_rigidMaxSite * sizeof(double), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_rigidGroups, m_rigidGroups, m_nRigidBody * m_rigidMaxSite * sizeof(int), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_rigidSingles, m_rigidSingles, (m_nDegFreedom/3 - 2*m_nRigidBody) * sizeof(int), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_rigidInverse, m_rigidInverse, m_nRigidBody * 3 * 3 * sizeof(double), cudaMemcpyHostToDevice) );

	CudaSafeCall( cudaMemcpy(m_d_nRigidBody, &m_nRigidBody, sizeof(int), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_rigidMaxSite, &m_rigidMaxSite, sizeof(int), cudaMemcpyHostToDevice) );

	if (m_coords != NULL) {
		CudaSafeCall( cudaMemcpy(m_d_coords, m_coords, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );
	}

	if (m_isAtomFrozen != NULL) {
		CudaSafeCall( cudaMemcpy(m_d_isAtomFrozen, m_isAtomFrozen, numDimensions/3 * sizeof(bool), cudaMemcpyHostToDevice) );
	}

	CudaSafeCall( cudaMemcpy(m_d_zerosx, zerosx, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_zerosy, zerosy, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_zerosz, zerosz, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );

	CudaSafeCall( cudaMemcpy(m_d_nLargeRigidBody, &m_nLargeRigidBody, sizeof(int), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_largeRigidIndices, m_largeRigidIndices, m_nLargeRigidBody * sizeof(int), cudaMemcpyHostToDevice) );

	CudaSafeCall( cudaMemcpy(m_d_ljAddRep, m_ljAddRep, m_nAddTarget * m_nAddTarget * sizeof(double), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(m_d_ljAddAtt, m_ljAddAtt, m_nAddTarget * m_nAddTarget * sizeof(double), cudaMemcpyHostToDevice) );

	delete [] zerosx;
	delete [] zerosy;
	delete [] zerosz;
}



CostFunction::~CostFunction() 
{
	CudaSafeCall( cudaFree(m_d_zerosx) );
	CudaSafeCall( cudaFree(m_d_zerosy) );
	CudaSafeCall( cudaFree(m_d_zerosz) );

	if (m_rigidInverse != NULL) {
		CudaSafeCall( cudaFree(m_d_gkAtoms));
	}

	CudaSafeCall( cudaFree(m_d_nRigidSitesPerBody));
	CudaSafeCall( cudaFree(m_d_sitesRigidBody));
	CudaSafeCall( cudaFree(m_d_rigidGroups));
	CudaSafeCall( cudaFree(m_d_rigidSingles));
	CudaSafeCall( cudaFree(m_d_rigidInverse));

	CudaSafeCall( cudaFree(m_d_xRigid) );
	CudaSafeCall( cudaFree(m_d_gkRigid) );

	CudaSafeCall( cudaFree(m_d_nRigidBody) );
	CudaSafeCall( cudaFree(m_d_rigidMaxSite) );

	if (m_coords != NULL) {
		CudaSafeCall( cudaFree(m_d_coords) );
		CudaSafeCall( cudaFree(m_d_energy) );
		CudaSafeCall( cudaFree(m_d_gradient) );
	}

	if (m_isAtomFrozen != NULL) {
		CudaSafeCall( cudaFree(m_d_isAtomFrozen) );
	}

	CudaSafeCall( cudaFree(m_d_nLargeRigidBody) );
	CudaSafeCall( cudaFree(m_d_largeRigidIndices) );

	CudaSafeCall( cudaFree(m_d_ljAddRep) );
	CudaSafeCall( cudaFree(m_d_ljAddAtt) );

	delete [] m_largeRigidIndices;
}



// Get/set functions. 

size_t CostFunction::getNumDimensions() const
{
	return m_numDimensions;
}

int CostFunction::getnCalls() const
{
	return m_nCalls;
}

bool CostFunction::getIsColdFusion() const
{
	return m_isColdFusion;
}

bool CostFunction::getIsRayleighRitz() const
{
	return m_isRayleighRitz;
}

void CostFunction::setIsRayleighRitz(bool isRayleighRitz)
{
	m_isRayleighRitz = isRayleighRitz;
}

void CostFunction::setDevCoords(double *d_coords)
{
	CudaSafeCall( cudaMemcpy(m_d_coords, d_coords, m_numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );
}

void CostFunction::setDevEnergy(double *d_energy)
{
	CudaSafeCall( cudaMemcpy(m_d_energy, d_energy, sizeof(double), cudaMemcpyDeviceToDevice) );
}

void CostFunction::setDevGradient(double *d_gradient)
{
	CudaSafeCall( cudaMemcpy(m_d_gradient, d_gradient, m_numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );
}

void CostFunction::setIsBfgsts(bool isBfgsts)
{
	m_isBfgsts = isBfgsts;
}



namespace gpu_secdiag
{
	// Variables on the GPU. 
	__constant__ double zeta;
	__constant__ double minusZeta;

	// Kernels. 
	__global__ void updateDiag(const double *d_fPlus, const double *d_fMinus, const double *m_d_energy, double *d_f, 
			double *d_eigen2, double *d_eigen3, double *d_twoEigen, double *d_twoEigen2);
	__global__ void updateDiag4(double *d_eigen4, double *d_twoEigen4);
	__global__ void calculateDebug(const double *d_fPlus, const double *d_fMinus, double *d_temp);
}



namespace gpu_freeze
{
	// Variables on the GPU. 
	__constant__ int numDimensions;

	//Kernels. 
	__global__ void freezeGrad(double *d_gradf, const bool *m_d_isAtomFrozen);
}



namespace gpu_orthogopt
{
	// Variables on the GPU. 
	__constant__ int numDimensions;

	// Kernels. 
	__global__ void updateCm(double *d_cmx, double *d_cmy, double *d_cmz);
	__global__ void updatevDotTrans(double *d_tmp, double *d_vDot);
	__global__ void updateTmp(double *d_tmp);
	__global__ void updatevDotRot(const double *d_tmp, double *d_tmp2, double *d_vDot);

	__global__ void rotationalx1(const double *d_x, const double *m_d_coords, double *d_tempArray, const double *d_cmy, 
			const double *d_cmz);
	__global__ void rotationalx2(double *d_x, const double *m_d_coords, const double *d_cmy, const double *d_cmz, 
			const double *d_tmp2);
	__global__ void rotationaly1(const double *d_x, const double *m_d_coords, double *d_tempArray, const double *d_cmx, 
			const double *d_cmz);
	__global__ void rotationaly2(double *d_x, const double *m_d_coords, const double *d_cmx, const double *d_cmz, 
			const double *d_tmp2);
	__global__ void rotationalz1(const double *d_x, const double *m_d_coords, double *d_tempArray, const double *d_cmx, 
			const double *d_cmz);
	__global__ void rotationalz2(double *d_x, const double *m_d_coords, const double *d_cmx, const double *d_cmy, 
			const double *d_tmp2);
}



// d_x is the estimated eigenvector and direction along which to compute the curvature. 
// d_f is the curvature. 
// d_gradf is the gradient of the curvature. 
void CostFunction::basePotential(double *d_x, double *d_f, double *d_gradf)
{
	bool printingOn = m_debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = m_debugPrinting.getFileHandle();

	bool isTiming = m_timer_potential.getTimingOn();
	if (isTiming) {
		m_timer_potential.start();
	}

	// This computes the curvature of the potential energy surface along the direction d_x at point m_d_coords. 
	if (m_isRayleighRitz) {
		const double zeta        =  1.0e-3;
		const double minusZeta   = -zeta;
		const double minusOne    = -1.0;
		const double oneOverZeta =  1.0/zeta;
		const double twoOverZeta =  2.0/zeta;

		double proj = 0.0;

		double *d_xCopy;
		double *d_nDimTmp;
		double *d_fPlus;
		double *d_fMinus;
		double *d_grad1;
		double *d_grad2;
		double *d_zeta;
		double *d_minusZeta;
		double *d_eigen2;
		double *d_eigen3;
		double *d_eigen4;
		double *d_twoEigen;
		double *d_twoEigen2;
		double *d_twoEigen4;
		double *d_tmp;

		CudaSafeCall( cudaMalloc(&d_xCopy, m_numDimensions * sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_nDimTmp, m_numDimensions * sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_fPlus, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_fMinus, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_grad1, m_numDimensions * sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_grad2, m_numDimensions * sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_eigen2, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_eigen3, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_eigen4, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_twoEigen, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_twoEigen2, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_twoEigen4, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_tmp, sizeof(double)) );

		CudaSafeCall( cudaGetSymbolAddress((void**)&d_zeta, gpu_secdiag::zeta) );
		CudaSafeCall( cudaGetSymbolAddress((void**)&d_minusZeta, gpu_secdiag::minusZeta) );

		CudaSafeCall( cudaMemcpy(d_xCopy, d_x, m_numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );

		CudaSafeCall( cudaMemcpyToSymbol(gpu_secdiag::zeta, &zeta, sizeof(double)) );
		CudaSafeCall( cudaMemcpyToSymbol(gpu_secdiag::minusZeta, &minusZeta, sizeof(double)) );

		if ((m_nSecDiag > 3) || (m_nSecDiag < 1)) {
			std::cerr << "This method for calculating secdiag is not recognised. " << std::endl;
			exit(EXIT_FAILURE);
		}

		bool shouldNormalise = true;
		orthogopt(d_xCopy, shouldNormalise);

		freeze(d_xCopy);

		// Compute the energy and gradient at coords + zeta*xCopy. 
		m_cublas.dispatchAxpy(m_numDimensions, d_nDimTmp, m_d_coords, d_xCopy, d_zeta); // nDimTmp = coords + zeta*xCopy
		rigidTransforms(d_nDimTmp, d_fPlus, d_grad1);

		// Second order central differences expansion.
		if (m_nSecDiag <= 2) {
			// Compute the energy and gradient at coords - zeta*xCopy. 
			m_cublas.dispatchAxpy(m_numDimensions, d_nDimTmp, m_d_coords, d_xCopy, d_minusZeta); // nDimTmp = coords - zeta*xCopy
			rigidTransforms(d_nDimTmp, d_fMinus, d_grad2);

			m_cublas.dispatchAxpy(m_numDimensions, d_nDimTmp, d_grad1, d_grad2, &minusOne, false); // nDimTmp = grad1 - grad2
			m_cublas.dispatchDot(m_numDimensions, d_eigen2, d_nDimTmp, d_xCopy); // eigen2 = nDimTmp Dot xCopy 

			gpu_secdiag::updateDiag<<<1,1>>>(d_fPlus, d_fMinus, m_d_energy, d_f, d_eigen2, d_eigen3, d_twoEigen, 
					d_twoEigen2);
			CudaCheckError();
			cudaDeviceSynchronize();

			if (printingOn) {
				gpu_secdiag::calculateDebug<<<1,1>>>(d_fPlus, d_fMinus, d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();
			}

			m_cublas.dispatchScale(m_numDimensions, d_nDimTmp, d_nDimTmp, &oneOverZeta, false); // nDimTmp = (1/zeta)*nDimTmp
			if (m_nSecDiag == 2) {
				m_cublas.dispatchAxpy(m_numDimensions, d_gradf, d_nDimTmp, d_xCopy, d_twoEigen2); // gradf = nDimTmp + twoEigen2*xCopy
			}
			else {
				m_cublas.dispatchAxpy(m_numDimensions, d_gradf, d_nDimTmp, d_xCopy, d_twoEigen); // gradf = nDimTmp + twoEigen*xCopy
			}
		}
		// First order forward finite differences (m_nSecDiag == 3). 
		else {
			m_cublas.dispatchAxpy(m_numDimensions, d_nDimTmp, d_grad1, m_d_gradient, &minusOne, false); // nDimTmp = grad1 - gradient
			m_cublas.dispatchDot(m_numDimensions, d_eigen4, d_nDimTmp, d_xCopy); // eigen4 = nDimTmp Dot xCopy
			gpu_secdiag::updateDiag4<<<1,1>>>(d_eigen4, d_twoEigen4);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchScale(m_numDimensions, d_nDimTmp, d_nDimTmp, &twoOverZeta, false); // nDimTmp = (2/zeta)*nDimTmp
			m_cublas.dispatchAxpy(m_numDimensions, d_gradf, d_nDimTmp, d_xCopy, d_twoEigen4); // gradf = nDimTmp + twoEigen4*xCopy
		}

		shouldNormalise = false;
		orthogopt(d_gradf, shouldNormalise);

		m_cublas.dispatchDot(m_numDimensions, &proj, d_gradf, d_xCopy, false); // proj = gradf Dot xCopy
		proj = -proj;
		m_cublas.dispatchAxpy(m_numDimensions, d_gradf, d_gradf, d_xCopy, &proj, false); // gradf += proj*xCopy

		if (printingOn) {
			double printDiag, printDiag2, printDiag3, printEnePlus, printEneMinus, printEnergy;

			CudaSafeCall( cudaMemcpy(&printDiag,     d_f,        sizeof(double), cudaMemcpyDeviceToHost) );
			CudaSafeCall( cudaMemcpy(&printDiag2,    d_eigen2,   sizeof(double), cudaMemcpyDeviceToHost) );
			CudaSafeCall( cudaMemcpy(&printDiag3,    d_eigen3,   sizeof(double), cudaMemcpyDeviceToHost) );
			CudaSafeCall( cudaMemcpy(&printEnePlus,  d_fPlus,    sizeof(double), cudaMemcpyDeviceToHost) );
			CudaSafeCall( cudaMemcpy(&printEneMinus, d_fMinus,   sizeof(double), cudaMemcpyDeviceToHost) );
			CudaSafeCall( cudaMemcpy(&printEnergy,   m_d_energy, sizeof(double), cudaMemcpyDeviceToHost) );

			fileHandle << " D,D2,D3,E+,E-,E = " << std::setw(15) << std::setprecision(5) << printDiag 
				<< std::setw(15) << std::setprecision(5) << printDiag2 
				<< std::setw(15) << std::setprecision(5) << printDiag3 
				<< std::setw(20) << std::setprecision(12) << printEnePlus 
				<< std::setw(20) << std::setprecision(12) << printEneMinus 
				<< std::setw(20) << std::setprecision(12) << printEnergy << std::endl;

			if (m_nSecDiag <= 2) {
				double printTmp;
				CudaSafeCall( cudaMemcpy(&printTmp, d_tmp, sizeof(double), cudaMemcpyDeviceToHost) );
				fileHandle << " Predicted gradient component = " 
					<< std::setw(20) << std::setprecision(10) << printTmp << std::endl;
			}
		}

		if (m_nSecDiag == 2) {
			CudaSafeCall( cudaMemcpy(d_f, d_eigen2, sizeof(double), cudaMemcpyDeviceToDevice) );
		}
		else if (m_nSecDiag == 3) {
			CudaSafeCall( cudaMemcpy(d_f, d_eigen4, sizeof(double), cudaMemcpyDeviceToDevice) );
		}

		CudaSafeCall( cudaFree(d_xCopy) );
		CudaSafeCall( cudaFree(d_nDimTmp) );
		CudaSafeCall( cudaFree(d_fPlus) );
		CudaSafeCall( cudaFree(d_fMinus) );
		CudaSafeCall( cudaFree(d_grad1) );
		CudaSafeCall( cudaFree(d_grad2) );
		CudaSafeCall( cudaFree(d_eigen2) );
		CudaSafeCall( cudaFree(d_eigen3) );
		CudaSafeCall( cudaFree(d_eigen4) );
		CudaSafeCall( cudaFree(d_twoEigen) );
		CudaSafeCall( cudaFree(d_twoEigen2) );
		CudaSafeCall( cudaFree(d_twoEigen4) );
		CudaSafeCall( cudaFree(d_tmp) );

	}
	else {
		rigidTransforms(d_x, d_f, d_gradf);
	}

	if (isTiming) {
		m_timer_potential.stop();
	}

}



void CostFunction::rigidTransforms(double *d_x, double *d_f, double *d_gradf)
{
	if (!m_isAtomisticNotRigid) {
		// Rigid bodies coordinate transformation. 
		transformRigidToC(d_x);
	}

	computeEnergyAndGradient(d_x, d_f, d_gradf);

	if (!m_isAtomisticNotRigid) {
		CudaSafeCall( cudaMemcpy(m_d_gkAtoms, d_gradf, m_numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );
		// Rigid bodies gradient transformation. 
		transformGrad(d_gradf, d_x);
	}

	freeze(d_gradf);

	m_nCalls += 1;
}



void CostFunction::calculateRms(const double *d_gradf, double *outRms)
{
	double gkNormSquared;
	m_cublas.dispatchDot(m_numDimensions, &gkNormSquared, d_gradf, d_gradf, false); // gkNormSquared = gradf Dot gradf
	if (!m_isAtomisticNotRigid) {
		*outRms = fmax(sqrt(gkNormSquared/m_nDegFreedom), 1.0e-100);
		if (m_isAaConvergence && (*outRms < m_aaConvThreshold) && !m_isBfgsts) {
			aaConvergence(m_d_gkAtoms, outRms);
			*outRms = fmax(sqrt(*outRms/m_numDimensions), 1.0e-100);
		}
	}
	else {
		*outRms = fmax(sqrt(gkNormSquared/m_numDimensions), 1.0e-100);
	}
}



void CostFunction::freeze(double *d_gradf)
{
	if (m_shouldFreeze) {
		CudaSafeCall( cudaMemcpyToSymbol(gpu_freeze::numDimensions, &m_numDimensions,  sizeof(int)) );

		dim3 blockDim;
		dim3 gridDim;

		blockDim.x = 1024;
		gridDim.x = ((m_numDimensions/3) + blockDim.x - 1)/blockDim.x;

		gpu_freeze::freezeGrad<<<gridDim, blockDim>>>(d_gradf, m_d_isAtomFrozen);
		CudaCheckError();
		cudaDeviceSynchronize();

	}
}



void CostFunction::orthogopt(double *d_x, const bool shouldNormalise)
{
	bool printingOn = m_debugPrinting.getPrintingOn();
	std::ofstream &fileHandle = m_debugPrinting.getFileHandle();

	CudaSafeCall( cudaMemcpyToSymbol(gpu_orthogopt::numDimensions, &m_numDimensions,  sizeof(int)) );
	double *d_tmp;
	CudaSafeCall( cudaMalloc(&d_tmp, sizeof(double)) );
	if (m_nFreeze < 3) {
		dim3 blockDim;
		dim3 gridDim;

		double *d_tempArray;

		double *d_cmx;
		double *d_cmy;
		double *d_cmz;

		double *d_tmp2;
		double *d_vDot;

		CudaSafeCall( cudaMalloc(&d_cmx, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_cmy, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_cmz, sizeof(double)) );

		CudaSafeCall( cudaMalloc(&d_tmp2, sizeof(double)) );
		CudaSafeCall( cudaMalloc(&d_vDot, sizeof(double)) );


		CudaSafeCall( cudaMalloc(&d_tempArray, m_numDimensions * sizeof(double)) );

		CudaSafeCall( cudaMemcpy(d_tempArray, m_d_zerosx, m_numDimensions * sizeof(double), cudaMemcpyDeviceToDevice) );

		double vDot;

		m_cublas.dispatchDot(m_numDimensions, d_cmx, m_d_zerosx, m_d_coords); // cmx = zerosx Dot coords
		m_cublas.dispatchDot(m_numDimensions, d_cmy, m_d_zerosy, m_d_coords); // cmy = zerosy Dot coords
		m_cublas.dispatchDot(m_numDimensions, d_cmz, m_d_zerosz, m_d_coords); // cmz = zerosz Dot coords

		gpu_orthogopt::updateCm<<<1,1>>>(d_cmx, d_cmy, d_cmz);
		CudaCheckError();
		cudaDeviceSynchronize();

		size_t nCheck = 0;
		for (nCheck =0; nCheck < 100; ++nCheck) {
			vDot = 0.0;
			CudaSafeCall( cudaMemcpy(d_vDot, &vDot, sizeof(double), cudaMemcpyHostToDevice) );

			m_cublas.dispatchDot(m_numDimensions, d_tmp, m_d_zerosx, d_x); // tmp = zerosx Dot x

			gpu_orthogopt::updatevDotTrans<<<1,1>>>(d_tmp, d_vDot);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy (m_numDimensions, d_x, d_x, m_d_zerosx, d_tmp); // x += zerosx * tmp

			if (shouldNormalise) {
				m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

				gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();

				m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp);  // x = tmp*x
			}

			m_cublas.dispatchDot(m_numDimensions, d_tmp, m_d_zerosy, d_x); // tmp = zerosy Dot x

			gpu_orthogopt::updatevDotTrans<<<1,1>>>(d_tmp, d_vDot);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy (m_numDimensions, d_x, d_x, m_d_zerosy, d_tmp); // x += zerosy * tmp

			if (shouldNormalise) {
				m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

				gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();

				m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp); // x = tmp*x
			}

			m_cublas.dispatchDot(m_numDimensions, d_tmp, m_d_zerosz, d_x); // tmp = zerosz Dot x

			gpu_orthogopt::updatevDotTrans<<<1,1>>>(d_tmp, d_vDot);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchAxpy (m_numDimensions, d_x, d_x, m_d_zerosz, d_tmp); // x += zerosz * tmp

			if (shouldNormalise) {
				m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

				gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();

				m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp); // x = tmp*x
			}


			blockDim.x = 1024;
			gridDim.x = ((m_numDimensions) + blockDim.x - 1)/blockDim.x;

			gpu_orthogopt::rotationalx1<<<gridDim, blockDim>>>(d_x, m_d_coords, d_tempArray, d_cmy, d_cmz);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchDot(m_numDimensions, d_tmp, m_d_zerosy, d_tempArray); // tmp = zerosy Dot tempArray
			m_cublas.dispatchDot(m_numDimensions, d_tmp2, m_d_zerosz, d_tempArray); // tmp2 = zerosz Dot tempArray

			gpu_orthogopt::updatevDotRot<<<1,1>>>(d_tmp, d_tmp2, d_vDot);
			CudaCheckError();
			cudaDeviceSynchronize();

			gpu_orthogopt::rotationalx2<<<gridDim, blockDim>>>(d_x, m_d_coords, d_cmy, d_cmz, d_tmp2);
			CudaCheckError();
			cudaDeviceSynchronize();

			if (shouldNormalise) {
				m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

				gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();

				m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp); // x = tmp*x
			}

			gpu_orthogopt::rotationaly1<<<gridDim, blockDim>>>(d_x, m_d_coords, d_tempArray, d_cmx, d_cmz);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchDot(m_numDimensions, d_tmp, m_d_zerosx, d_tempArray); // tmp = zerosx Dot tempArray
			m_cublas.dispatchDot(m_numDimensions, d_tmp2, m_d_zerosz, d_tempArray); // tmp2 = zerosz Dot tempArray

			gpu_orthogopt::updatevDotRot<<<1,1>>>(d_tmp, d_tmp2, d_vDot);
			CudaCheckError();
			cudaDeviceSynchronize();

			gpu_orthogopt::rotationaly2<<<gridDim, blockDim>>>(d_x, m_d_coords, d_cmx, d_cmz, d_tmp2);
			CudaCheckError();
			cudaDeviceSynchronize();

			if (shouldNormalise) {
				m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

				gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();

				m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp);  // x = tmp*x
			}

			gpu_orthogopt::rotationalz1<<<gridDim, blockDim>>>(d_x, m_d_coords, d_tempArray, d_cmx, d_cmy);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchDot(m_numDimensions, d_tmp, m_d_zerosx, d_tempArray); // tmp = zerosx Dot tempArray
			m_cublas.dispatchDot(m_numDimensions, d_tmp2, m_d_zerosy, d_tempArray); // tmp2 = zerosz Dot tempArray

			gpu_orthogopt::updatevDotRot<<<1,1>>>(d_tmp, d_tmp2, d_vDot);
			CudaCheckError();
			cudaDeviceSynchronize();

			gpu_orthogopt::rotationalz2<<<gridDim, blockDim>>>(d_x, m_d_coords, d_cmx, d_cmy, d_tmp2);
			CudaCheckError();
			cudaDeviceSynchronize();

			if (shouldNormalise) {
				m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

				gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
				CudaCheckError();
				cudaDeviceSynchronize();

				m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp); // x = tmp*x
			}

			CudaSafeCall( cudaMemcpy(&vDot, d_vDot, sizeof(double), cudaMemcpyDeviceToHost) );
			if (vDot <= 1.0e-6) {
				break;
			}
		}

		if (nCheck==100) {
			fileHandle << " Warning: cannot orthogonalise to known eigenvectors using orthogopt" << std::endl;
		}

		CudaSafeCall( cudaFree(d_cmx) );
		CudaSafeCall( cudaFree(d_cmy) );
		CudaSafeCall( cudaFree(d_cmz) );

		CudaSafeCall( cudaFree(d_tmp2) );
		CudaSafeCall( cudaFree(d_vDot) );

		CudaSafeCall( cudaFree(d_tempArray) );
	}
	else {
		// Simply normalise. 
		if (shouldNormalise) {
			m_cublas.dispatchNrm2(m_numDimensions, d_tmp, d_x); // tmp = sqrt(x Dot x)

			gpu_orthogopt::updateTmp<<<1,1>>>(d_tmp);
			CudaCheckError();
			cudaDeviceSynchronize();

			m_cublas.dispatchScale(m_numDimensions, d_x, d_x, d_tmp); // x = tmp*x
		}
	}

	CudaSafeCall( cudaFree(d_tmp) );
}



namespace gpu_secdiag
{
	__global__ void updateDiag(const double *d_fPlus, const double *d_fMinus, const double *m_d_energy, double *d_f, 
			double *d_eigen2, double *d_eigen3, double *d_twoEigen, double *d_twoEigen2)
	{
		// d_f is a second order central differences expansion of the curvature using energies. 
		*d_f = (*d_fPlus + *d_fMinus - 2.0*(*m_d_energy))/(zeta*zeta);
		// d_eigen2 is a second order central differences expansion of the curvature using gradients. 
		*d_eigen2 = *d_eigen2 / (2.0*zeta);
		// d_eigen3 is a more accurate estimate of the diagonal second derivative, 
		// but it cannot be differentiated analytically.
		*d_eigen3 = 2.0*(*d_f - *d_eigen2/2.0);

		*d_twoEigen  = -2.0*(*d_f);
		*d_twoEigen2 = -2.0*(*d_eigen2);
	}

	__global__ void updateDiag4(double *d_eigen4, double *d_twoEigen4)
	{
		*d_eigen4 = *d_eigen4/zeta;
		*d_twoEigen4 = -2.0*(*d_eigen4);
	}

	__global__ void calculateDebug(const double *d_fPlus, const double *d_fMinus, double *d_tmp)
	{
		*d_tmp = (*d_fPlus - *d_fMinus)/(2*zeta);
	}
}



namespace gpu_freeze
{

	__global__ void freezeGrad(double *d_gradient, const bool *m_d_isAtomFrozen)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;
		while (tid < (numDimensions/3)) {
			if (m_d_isAtomFrozen[tid]) {
				d_gradient[3*tid+0] = 0.0;
				d_gradient[3*tid+1] = 0.0;
				d_gradient[3*tid+2] = 0.0;
			}

			tid += blockDim.x * gridDim.x;
		}
	}
}



namespace gpu_orthogopt
{

	__global__ void rotationalx1(const double *d_x, const double *m_d_coords, double *d_tempArray, const double *d_cmy, 
			const double *d_cmz)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)) {
			d_tempArray[3*tid+1] = d_x[3*tid+1]*(m_d_coords[3*tid+2] - *d_cmz) - 
				(d_x[3*tid+2]*(m_d_coords[3*tid+1] - *d_cmy));
			d_tempArray[3*tid+2] = (m_d_coords[3*tid+2] - *d_cmz)*(m_d_coords[3*tid+2] - *d_cmz) + 
				(m_d_coords[3*tid+1] - *d_cmy)*(m_d_coords[3*tid+1] - *d_cmy);

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void rotationalx2(double *d_x, const double *m_d_coords, const double *d_cmy, const double *d_cmz, 
			const double *d_tmp2)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)) {
			d_x[3*tid+1] -= *d_tmp2*(m_d_coords[3*tid+2] - *d_cmz);
			d_x[3*tid+2] += *d_tmp2*(m_d_coords[3*tid+1] - *d_cmy);

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void rotationaly1(const double *d_x, const double *m_d_coords, double *d_tempArray, const double *d_cmx, 
			const double *d_cmz)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)) {
			d_tempArray[3*tid] = -d_x[3*tid+0]*(m_d_coords[3*tid+2] - *d_cmz) + 
				(d_x[3*tid+2]*(m_d_coords[3*tid+0] - *d_cmx));
			d_tempArray[3*tid+2] = (m_d_coords[3*tid+2] - *d_cmz)*(m_d_coords[3*tid+2] - *d_cmz) + 
				(m_d_coords[3*tid+0] - *d_cmx)*(m_d_coords[3*tid+0] - *d_cmx);

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void rotationaly2(double *d_x, const double *m_d_coords, const double *d_cmx, const double *d_cmz, 
			const double *d_tmp2)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)) {
			d_x[3*tid+0] += *d_tmp2*(m_d_coords[3*tid+2] - *d_cmz);
			d_x[3*tid+2] -= *d_tmp2*(m_d_coords[3*tid+0] - *d_cmx);

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void rotationalz1(const double *d_x, const double *m_d_coords, double *d_tempArray, const double *d_cmx, 
			const double *d_cmy)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)) {
			d_tempArray[3*tid+0] = d_x[3*tid+0]*(m_d_coords[3*tid+1] - *d_cmy) - 
				(d_x[3*tid+1]*(m_d_coords[3*tid+0] - *d_cmx));
			d_tempArray[3*tid+1] = (m_d_coords[3*tid+1] - *d_cmy)*(m_d_coords[3*tid+1] - *d_cmy) + 
				(m_d_coords[3*tid+0] - *d_cmx)*(m_d_coords[3*tid+0] - *d_cmx);

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void rotationalz2(double *d_x, const double *m_d_coords, const double *d_cmx, const double *d_cmy, 
			const double *d_tmp2)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)) {
			d_x[3*tid+0] -= *d_tmp2*(m_d_coords[3*tid+1] - *d_cmy);
			d_x[3*tid+1] += *d_tmp2*(m_d_coords[3*tid+0] - *d_cmx);

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void updateCm(double *d_cmx, double *d_cmy, double *d_cmz)
	{
		*d_cmx = *d_cmx/(numDimensions/3.0);
		*d_cmy = *d_cmy/(numDimensions/3.0);
		*d_cmz = *d_cmz/(numDimensions/3.0);
	}

	__global__ void updatevDotTrans(double *d_tmp, double *d_vDot)
	{
		double rootN = sqrt(numDimensions/3.0);
		*d_tmp = *d_tmp/rootN;
		*d_vDot = fmax(*d_vDot, abs(*d_tmp));
		*d_tmp = -(*d_tmp)/rootN;
	}

	__global__ void updateTmp(double *d_tmp)
	{
		*d_tmp = 1.0/(*d_tmp);
	}

	__global__ void updatevDotRot(const double *d_tmp, double *d_tmp2, double *d_vDot)
	{
		*d_vDot = fmax(*d_vDot, abs(*d_tmp)/sqrt(*d_tmp2));
		*d_tmp2 = *d_tmp/(*d_tmp2);
	}

}
