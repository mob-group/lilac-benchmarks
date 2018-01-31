/**
 *
 * File rigid_bodies.cu: Implementation of rigid body transformation functions from class CostFunction.
 *
 **/

#include "cost_function.h"

namespace gpu_rigid_bodies
{
	// Variables on the GPU. 

	__device__ bool shouldFindDeriv;


	// Kernels. 

	__global__ void transform(const int *m_d_nRigidBody, const double *m_d_xRigid, double *d_grmi);

	__global__ void transform2(double *d_x, const double *m_d_xRigid, const int *m_d_nRigidBody, 
			const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const double *m_d_sitesRigidBody, const double *d_grmi, const int *m_d_rigidMaxSite);

	__global__ void singleAtoms(double *d_x, const double *m_d_xRigid, const int *m_d_rigidSingles, const int *m_d_nRigidBody, 
			const int singlesThreads);

	__global__ void gradTransform1a(const double *m_d_xRigid, double *d_grmi1, double *d_grmi2, double *d_grmi3, 
			const int *m_d_nRigidBody);

	__global__ void gradTransform1b(const int *m_d_nRigidBody, const double *m_d_xRigid, double *d_grmi0);

	__global__ void intermediate(const double *d_gk, const int *m_d_nRigidSitesPerBody, const double *d_grmi1, 
			const double *d_grmi2, const double *d_grmi3, const double *m_d_sitesRigidBody, 
			const int *m_d_rigidGroups, double *d_tempArray, const int *m_d_nRigidBody, 
			const int *m_d_rigidMaxSite);

	__global__ void gradSingleAtoms(const double *d_gk, double *m_d_gkRigid, const int *m_d_rigidSingles, 
			const int *m_d_nRigidBody, const int singlesThreads);

	__global__ void aaConvRmdrvt(double *d_rmi10, double *d_rmi20, double *d_rmi30);

	__global__ void intermediate2(const double *d_gk, const int *m_d_nRigidSitesPerBody, const double *d_grmi0, 
			const double *d_grmi10, const double *d_grmi20, const double *d_grmi30, 
			const double *m_d_sitesRigidBody, const int *m_d_rigidGroups, double *d_tempArray, 
			const int *m_d_nRigidBody, const int *m_d_rigidMaxSite);

	__global__ void aaConvTorque(const double *d_torques, const double *m_d_gkRigid, const double *d_grmi0, 
			const int *m_d_nRigidSitesPerBody, const double *m_d_rigidInverse, double *d_rmsArray, 
			const int *m_d_nRigidBody);

	__global__ void warpReduce1(const int *m_d_nRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const int *m_d_rigidMaxSite, const double *d_gk, double *m_d_gkRigid);

	__global__ void warpReduce2(const int *m_d_nRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidMaxSite, 
			double *m_d_gkRigid, const double *d_tempArray);

	__global__ void warpReduce3(const int *m_d_nRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidMaxSite, 
			double *d_torques, const double *d_tempArray);

	__global__ void fullReduce1a(const int *m_d_nLargeRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const int *m_d_rigidMaxSite, const double *d_gk, const int *m_d_largeRigidIndices, double *d_outArray, 
			const int roundedMaxSite);

	__global__ void fullReduce2a(const int *m_d_nLargeRigidBody, double *m_d_gkRigid, const int *m_d_largeRigidIndices, 
			double *d_outArray, const int roundedMaxSite, const int outSize, const int blocks);

	__global__ void fullReduce1b(const int *m_d_nLargeRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups,
			const int *m_d_rigidMaxSite, const int *m_d_largeRigidIndices, double *d_outArray,
			const double *d_tempArray, const int roundedMaxSite);

	__global__ void fullReduce2b(const int *m_d_nLargeRigidBody, double *m_d_gkRigid, const int *m_d_largeRigidIndices, 
			double *d_outArray, const int *m_d_nRigidBody, const int roundedMaxSite, const int outSize, 
			const int blocks);

	__global__ void fullReduce2c(const int *m_d_nLargeRigidBody, double *torques, const int *m_d_largeRigidIndices, 
			double *d_outArray, const int *m_d_nRigidBody, const int roundedMaxSite, const int outSize, const int blocks);

	__global__ void fullReduce(const double *in, double* out, const int N);

	__global__ void fullDotReduce(const double *in, double* out, const int N);

}



void CostFunction::transformRigidToC(double *d_x)
{
	using namespace gpu_rigid_bodies;

	// Copy atomistic coordinates to rigid coordinates. 
	CudaSafeCall( cudaMemcpy(m_d_xRigid, d_x, m_nDegFreedom * sizeof(double), cudaMemcpyDeviceToDevice) );

	// Temporary global memory storage for rotation matrices. 
	double *d_grmi;
	CudaSafeCall( cudaMalloc(&d_grmi, 9 * m_nRigidBody * sizeof(double)) );

	bool shouldFindDerivatives = false;
	CudaSafeCall( cudaMemcpyToSymbol(gpu_rigid_bodies::shouldFindDeriv, &shouldFindDerivatives,  sizeof(bool)) );

	dim3 blockDim;
	blockDim.x = 1024;
	dim3 gridDim;
	gridDim.x = (m_nRigidBody + blockDim.x - 1)/blockDim.x;

	// Calculation of rotation matrix for each rigid body. 
	gpu_rigid_bodies::transform<<<gridDim, blockDim>>>(m_d_nRigidBody, m_d_xRigid, d_grmi);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Launch no. of threads equal to max. no. of sites per rigid body * number of rigid bodies for easy indexing. 
	// Not all threads wil be used as number of sites per rigid body can be less than maximum. 
	int noThreads = m_rigidMaxSite*m_nRigidBody; 
	blockDim.x = 1024;
	gridDim.x = (noThreads + blockDim.x - 1)/blockDim.x;

	// Rigid body to atomistic mapping. 
	gpu_rigid_bodies::transform2<<<gridDim, blockDim>>>(d_x, m_d_xRigid, m_d_nRigidBody, m_d_nRigidSitesPerBody, 
			m_d_rigidGroups, m_d_sitesRigidBody, d_grmi, m_d_rigidMaxSite); 
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaFree(d_grmi) );

	if (m_nDegFreedom > 6*m_nRigidBody) {
		int singlesThreads = (m_nDegFreedom - 6*m_nRigidBody)/3;

		blockDim.x = 1024;
		gridDim.x = (singlesThreads + blockDim.x - 1)/blockDim.x;

		// Treatment of single atoms not in rigid bodies. 
		gpu_rigid_bodies::singleAtoms<<<gridDim, blockDim>>>(d_x, m_d_xRigid, m_d_rigidSingles, m_d_nRigidBody, singlesThreads);
		CudaCheckError();
		cudaDeviceSynchronize();
	}
}



void CostFunction::transformGrad(double *d_gk, double *d_x)
{
	using namespace gpu_rigid_bodies;

	const size_t numDimensions = m_numDimensions;

	double *zeros;
	zeros = new double[numDimensions];
	for (size_t i = 0; i < numDimensions; ++i) {
		zeros[i] = 0.0;
	}

	// Initialise rigid gradient with zeros. 
	CudaSafeCall( cudaMemcpy(m_d_gkRigid, zeros, m_nDegFreedom * sizeof(double), cudaMemcpyHostToDevice) );

	// Temporary global memory storage for rotation matrices. 
	double *d_grmi1;
	double *d_grmi2;
	double *d_grmi3;

	CudaSafeCall( cudaMalloc(&d_grmi1, 9 * m_nRigidBody * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_grmi2, 9 * m_nRigidBody * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_grmi3, 9 * m_nRigidBody * sizeof(double)) );

	bool shouldFindDerivatives = true;
	CudaSafeCall( cudaMemcpyToSymbol(gpu_rigid_bodies::shouldFindDeriv, &shouldFindDerivatives,  sizeof(bool)) );

	dim3 blockDim;
	blockDim.x = 256;
	dim3 gridDim;
	gridDim.x = (m_nRigidBody + blockDim.x - 1)/blockDim.x;

	// Calculation of derivates of rotation matrices - one thread per rigid body. 
	gpu_rigid_bodies::gradTransform1a<<<gridDim, blockDim>>>(m_d_xRigid, d_grmi1, d_grmi2, d_grmi3, m_d_nRigidBody);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Reduction - projection of atomistic forces onto translational degrees of freedom of each rigid body.
	blockDim.x = 1024;
	gridDim.x = (32*m_nRigidBody + blockDim.x - 1)/blockDim.x; // 32 is warpSize
	// Reduction only takes place for rigid bodies with 32 or fewer sites. 
	gpu_rigid_bodies::warpReduce1<<<gridDim, blockDim>>>(m_d_nRigidBody, m_d_nRigidSitesPerBody, m_d_rigidGroups, m_d_rigidMaxSite, 
			d_gk, m_d_gkRigid);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Reduction for rigid bodies with more than 32 sites. 
	if (m_nLargeRigidBody != 0){
		int blockSize = 1024;
		blockDim.x = blockSize;

		// Round rigidMaxSite to nearest block size so that per block reduction doesn't sum components from different bodies. 
		int roundedMaxSite = ((m_rigidMaxSite + blockSize - 1)/blockSize)*blockSize;
		int numThreads = roundedMaxSite*m_nLargeRigidBody;

		int blocks = (numThreads + blockDim.x - 1)/blockDim.x;
		gridDim.x = blocks;

		double *d_outArray;
		CudaSafeCall( cudaMalloc(&d_outArray, 3 * blocks * sizeof(double)) );

		gpu_rigid_bodies::fullReduce1a<<<gridDim, blockDim>>>(m_d_nLargeRigidBody, m_d_nRigidSitesPerBody, 
				m_d_rigidGroups, m_d_rigidMaxSite, d_gk, m_d_largeRigidIndices, d_outArray, roundedMaxSite);
		CudaCheckError();
		cudaDeviceSynchronize();

		while (blocks > m_nLargeRigidBody) {
			int outSize = blocks;

			roundedMaxSite = ((blocks/m_nLargeRigidBody + blockSize - 1)/blockSize)*blockSize;
			numThreads = roundedMaxSite*m_nLargeRigidBody;

			blocks = (numThreads + blockDim.x - 1)/blockDim.x;
			gridDim.x = blocks;

			gpu_rigid_bodies::fullReduce2a<<<gridDim, blockDim>>>(m_d_nLargeRigidBody, m_d_gkRigid, 
					m_d_largeRigidIndices, d_outArray, roundedMaxSite, outSize, blocks);
			CudaCheckError();
			cudaDeviceSynchronize();
		}

		CudaSafeCall( cudaFree(d_outArray) );
	}

	// Temporary global memory storage for atomistic components of rigid body forces. 
	double *d_tempArray;
	int tempArraySize = 3 * m_nRigidBody * m_rigidMaxSite;
	CudaSafeCall( cudaMalloc(&d_tempArray, tempArraySize * sizeof(double)) );

	// Launch no. of threads equal to max. no. of sites per rigid body * number of rigid bodies for easy indexing. 
	// Not all threads wil be used as number of sites per rigid body can be less than maximum. 
	int noThreads = m_rigidMaxSite * m_nRigidBody;
	blockDim.x = 1024;
	gridDim.x = (noThreads + blockDim.x - 1)/blockDim.x;

	// Projection of forces for each atom onto rotational degrees of freedom of rigid bodies (not summed for each body). 
	gpu_rigid_bodies::intermediate<<<gridDim, blockDim>>>(d_gk, m_d_nRigidSitesPerBody, d_grmi1, d_grmi2, d_grmi3, 
			m_d_sitesRigidBody, m_d_rigidGroups, d_tempArray, 
			m_d_nRigidBody, m_d_rigidMaxSite);
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaFree(d_grmi1) );
	CudaSafeCall( cudaFree(d_grmi2) );
	CudaSafeCall( cudaFree(d_grmi3) );

	// Reduction of atomistic components (d_tempArray) gives forces projected onto rotational d.o.f. of each rigid body.
	blockDim.x = 1024;
	gridDim.x = (32*m_nRigidBody + blockDim.x - 1)/blockDim.x; // 32 is warpSize
	// Reduction only takes place for rigid bodies with 32 or fewer sites.
	gpu_rigid_bodies::warpReduce2<<<gridDim, blockDim>>>(m_d_nRigidBody, m_d_nRigidSitesPerBody, m_d_rigidMaxSite, 
			m_d_gkRigid, d_tempArray);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Reduction for rigid bodies with more than 32 sites.
	if (m_nLargeRigidBody != 0){
		int blockSize = 1024;
		blockDim.x = blockSize;

		// Round rigidMaxSite to nearest block size so that per block reduction doesn't sum components from different bodies.
		int roundedMaxSite = ((m_rigidMaxSite + blockSize - 1)/blockSize)*blockSize;
		int numThreads = roundedMaxSite*m_nLargeRigidBody;

		int blocks = (numThreads + blockDim.x - 1)/blockDim.x;
		gridDim.x = blocks;

		double *d_outArray;
		CudaSafeCall( cudaMalloc(&d_outArray, 3 * blocks * sizeof(double)) );

		gpu_rigid_bodies::fullReduce1b<<<gridDim, blockDim>>>(m_d_nLargeRigidBody, m_d_nRigidSitesPerBody, m_d_rigidGroups,
				m_d_rigidMaxSite, m_d_largeRigidIndices, d_outArray, d_tempArray, roundedMaxSite);
		CudaCheckError();
		cudaDeviceSynchronize();

		while (blocks > m_nLargeRigidBody) {
			int outSize = blocks;

			roundedMaxSite = ((blocks/m_nLargeRigidBody + blockSize - 1)/blockSize)*blockSize;
			numThreads = roundedMaxSite*m_nLargeRigidBody;

			blocks = (numThreads + blockDim.x - 1)/blockDim.x;
			gridDim.x = blocks;

			gpu_rigid_bodies::fullReduce2b<<<gridDim, blockDim>>>(m_d_nLargeRigidBody, m_d_gkRigid, 
					m_d_largeRigidIndices, d_outArray, m_d_nRigidBody, roundedMaxSite, outSize, blocks);
			CudaCheckError();
			cudaDeviceSynchronize();
		}

		CudaSafeCall( cudaFree(d_outArray) );
	}

	CudaSafeCall( cudaFree(d_tempArray) );

	if (m_nDegFreedom > 6*m_nRigidBody) {
		int singlesThreads = (m_nDegFreedom - 6*m_nRigidBody)/3;

		blockDim.x = 1024;
		gridDim.x = (singlesThreads + blockDim.x - 1)/blockDim.x;

		// Treatment of single atoms not in rigid bodies. 
		gpu_rigid_bodies::gradSingleAtoms<<<gridDim, blockDim>>>(d_gk, m_d_gkRigid, m_d_rigidSingles, 
				m_d_nRigidBody, singlesThreads);
		CudaCheckError();
		cudaDeviceSynchronize();
	}

	// Copy rigid coordinates/gradient to atomistic coordinates/gradient and pad with zeros. 
	CudaSafeCall( cudaMemcpy(d_x,  zeros, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );
	CudaSafeCall( cudaMemcpy(d_gk, zeros, numDimensions * sizeof(double), cudaMemcpyHostToDevice) );

	CudaSafeCall( cudaMemcpy(d_x,  m_d_xRigid,  m_nDegFreedom * sizeof(double), cudaMemcpyDeviceToDevice) );
	CudaSafeCall( cudaMemcpy(d_gk, m_d_gkRigid, m_nDegFreedom * sizeof(double), cudaMemcpyDeviceToDevice) );

	delete [] zeros;
}



void CostFunction::aaConvergence(const double *d_gk, double *outRms)
{
	using namespace gpu_rigid_bodies;

	*outRms = 0.0;

	double *d_rmi10;
	double *d_rmi20;
	double *d_rmi30;

	CudaSafeCall( cudaMalloc(&d_rmi10, 9 * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_rmi20, 9 * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_rmi30, 9 * sizeof(double)) );

	bool shouldFindDerivatives = true;
	CudaSafeCall( cudaMemcpyToSymbol(gpu_rigid_bodies::shouldFindDeriv, &shouldFindDerivatives,  sizeof(bool)) );

	gpu_rigid_bodies::aaConvRmdrvt<<<1, 1>>>(d_rmi10, d_rmi20, d_rmi30);
	CudaCheckError();
	cudaDeviceSynchronize();

	double *d_grmi0;

	CudaSafeCall( cudaMalloc(&d_grmi0, 9 * m_nRigidBody * sizeof(double)) );

	shouldFindDerivatives = false;
	CudaSafeCall( cudaMemcpyToSymbol(gpu_rigid_bodies::shouldFindDeriv, &shouldFindDerivatives,  sizeof(bool)) );

	dim3 blockDim;
	blockDim.x = 1024;
	dim3 gridDim;
	gridDim.x = (m_nRigidBody + blockDim.x - 1)/blockDim.x;

	gpu_rigid_bodies::gradTransform1b<<<gridDim, blockDim>>>(m_d_nRigidBody, m_d_xRigid, d_grmi0);
	CudaCheckError();
	cudaDeviceSynchronize();

	double *d_tempArray;
	int tempArraySize = 3 * m_nRigidBody * m_rigidMaxSite;
	CudaSafeCall( cudaMalloc(&d_tempArray, tempArraySize * sizeof(double)) );

	int noThreads = m_rigidMaxSite * m_nRigidBody;
	blockDim.x = 1024;
	gridDim.x = (noThreads + blockDim.x - 1)/blockDim.x;

	gpu_rigid_bodies::intermediate2<<<gridDim, blockDim>>>(d_gk, m_d_nRigidSitesPerBody, d_grmi0, d_rmi10, d_rmi20, 
			d_rmi30, m_d_sitesRigidBody, m_d_rigidGroups, d_tempArray, 
			m_d_nRigidBody, m_d_rigidMaxSite);
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaFree(d_rmi10) );
	CudaSafeCall( cudaFree(d_rmi20) );
	CudaSafeCall( cudaFree(d_rmi30) );

	double *d_torques;
	CudaSafeCall( cudaMalloc(&d_torques, 3 * m_nRigidBody * sizeof(double)) );

	// Reduction. 
	blockDim.x = 1024;
	gridDim.x = (32*m_nRigidBody + blockDim.x - 1)/blockDim.x; // 32 is warpSize
	// Reduction only takes place for rigid bodies with 32 or fewer sites.
	gpu_rigid_bodies::warpReduce3<<<gridDim, blockDim>>>(m_d_nRigidBody, m_d_nRigidSitesPerBody, m_d_rigidMaxSite, d_torques, 
			d_tempArray);
	CudaCheckError();
	cudaDeviceSynchronize();

	// Reduction for rigid bodies with more than 32 sites.
	if (m_nLargeRigidBody != 0){
		int blockSize = 1024;
		blockDim.x = blockSize;

		// Round rigidMaxSite to nearest block size so that per block reduction doesn't sum components from different bodies.
		int roundedMaxSite = ((m_rigidMaxSite + blockSize - 1)/blockSize)*blockSize;
		int numThreads = roundedMaxSite*m_nLargeRigidBody;

		int blocks = (numThreads + blockDim.x - 1)/blockDim.x;
		gridDim.x = blocks;

		double *d_outArray;
		CudaSafeCall( cudaMalloc(&d_outArray, 3 * blocks * sizeof(double)) );

		gpu_rigid_bodies::fullReduce1b<<<gridDim, blockDim>>>(m_d_nLargeRigidBody, m_d_nRigidSitesPerBody, m_d_rigidGroups,
				m_d_rigidMaxSite, m_d_largeRigidIndices, d_outArray, d_tempArray, roundedMaxSite);
		CudaCheckError();
		cudaDeviceSynchronize();

		while (blocks > m_nLargeRigidBody) {
			int outSize = blocks;

			roundedMaxSite = ((blocks/m_nLargeRigidBody + blockSize - 1)/blockSize)*blockSize;
			numThreads = roundedMaxSite*m_nLargeRigidBody;

			blocks = (numThreads + blockDim.x - 1)/blockDim.x;
			gridDim.x = blocks;

			gpu_rigid_bodies::fullReduce2c<<<gridDim, blockDim>>>(m_d_nLargeRigidBody, d_torques, 
					m_d_largeRigidIndices, d_outArray, m_d_nRigidBody, roundedMaxSite, outSize, blocks);
			CudaCheckError();
			cudaDeviceSynchronize();
		}

		CudaSafeCall( cudaFree(d_outArray) );
	}

	CudaSafeCall( cudaFree(d_tempArray) );

	double *d_rmsArray;
	CudaSafeCall( cudaMalloc(&d_rmsArray, m_nRigidBody * sizeof(double)) );

	blockDim.x = 1024;
	gridDim.x = (m_nRigidBody + blockDim.x - 1)/blockDim.x;

	gpu_rigid_bodies::aaConvTorque<<<gridDim, blockDim>>>(d_torques, m_d_gkRigid, d_grmi0, m_d_nRigidSitesPerBody, 
			m_d_rigidInverse, d_rmsArray, m_d_nRigidBody);
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaFree(d_grmi0) );
	CudaSafeCall( cudaFree(d_torques) );

	double temp1 = 0.0;
	double temp2 = 0.0;

	// Reduction. 
	blockDim.x = 512;
	gridDim.x = min((m_nRigidBody + blockDim.x - 1) / blockDim.x, 1024);

	double *d_outArray;
	CudaSafeCall( cudaMalloc(&d_outArray, gridDim.x * sizeof(double)) );

	gpu_rigid_bodies::fullReduce<<<gridDim, blockDim>>>(d_rmsArray, d_outArray, m_nRigidBody);
	CudaCheckError();
	cudaDeviceSynchronize();

	gpu_rigid_bodies::fullReduce<<<1, 1024>>>(d_outArray, d_outArray, gridDim.x);
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaMemcpy(&temp1, d_outArray, sizeof(double), cudaMemcpyDeviceToHost) );

	CudaSafeCall( cudaFree(d_outArray) );
	CudaSafeCall( cudaFree(d_rmsArray) );

	if (m_nDegFreedom > 6*m_nRigidBody) {
		int arraysize = m_nDegFreedom - 6*m_nRigidBody;

		// Dot product.  
		blockDim.x = 512;
		gridDim.x = min((arraysize + blockDim.x - 1) / blockDim.x, 1024);

		CudaSafeCall( cudaMalloc(&d_outArray, gridDim.x * sizeof(double)) );

		gpu_rigid_bodies::fullDotReduce<<<gridDim, blockDim>>>((m_d_gkRigid + 6*m_nRigidBody), d_outArray, arraysize);
		CudaCheckError();
		cudaDeviceSynchronize();

		gpu_rigid_bodies::fullReduce<<<1, 1024>>>(d_outArray, d_outArray, gridDim.x);
		CudaCheckError();
		cudaDeviceSynchronize();

		CudaSafeCall( cudaMemcpy(&temp2, d_outArray, sizeof(double), cudaMemcpyDeviceToHost) );

		CudaSafeCall( cudaFree(d_outArray) );
	}

	*outRms = temp1 + temp2;
}



namespace gpu_rigid_bodies
{
	__device__ void rmdrvt(double3 p, double rmi[9], double drmi1[9], double drmi2[9], double drmi3[9], 
			bool shouldFindDeriv)
	{
		double theta2 = p.x*p.x + p.y*p.y + p.z*p.z;
		if (theta2 < 1.0e-12) {
			// Rotation matrix if magnitude of rotation is zero. 
			rmi[0] =  1.0; // RM(1,1)
			rmi[1] =  p.z; // RM(2,1)
			rmi[2] = -p.y; // RM(3,1)
			rmi[3] = -p.z; // RM(1,2)
			rmi[4] =  1.0; // RM(2,2)
			rmi[5] =  p.x; // RM(3,2)
			rmi[6] =  p.y; // RM(1,3)
			rmi[7] = -p.x; // RM(2,3)
			rmi[8] =  1.0; // RM(3,3)

			if (shouldFindDeriv) {
				drmi1[0] =  0.0;
				drmi1[1] =  0.5*p.y;
				drmi1[2] =  0.5*p.z;
				drmi1[3] =  0.5*p.y;
				drmi1[4] = -1.0*p.x;
				drmi1[5] =  1.0;
				drmi1[6] =  0.5*p.z;
				drmi1[7] = -1.0;
				drmi1[8] = -1.0*p.x;

				drmi2[0] = -1.0*p.y;
				drmi2[1] =  0.5*p.x;
				drmi2[2] = -1.0;
				drmi2[3] =  0.5*p.x;
				drmi2[4] =  0.0;
				drmi2[5] =  0.5*p.z;
				drmi2[6] =  1.0;
				drmi2[7] =  0.5*p.z;
				drmi2[8] = -1.0*p.y;

				drmi3[0] = -1.0*p.z;
				drmi3[1] =  1.0;
				drmi3[2] =  0.5*p.x;
				drmi3[3] = -1.0;
				drmi3[4] = -1.0*p.z;
				drmi3[5] =  0.5*p.y;
				drmi3[6] =  0.5*p.x;
				drmi3[7] =  0.5*p.y;
				drmi3[8] =  0.0;
			}
		}

		else {
			double theta = sqrt(theta2); // Magnitude of rotation vector. 
			double ct = cos(theta);
			double st = sin(theta);
			double theta3 = 1.0/(theta2*theta);

			theta = 1.0 / theta;

			double3 pn;
			pn.x = theta*p.x;
			pn.y = theta*p.y;
			pn.z = theta*p.z; 

			// Skew-symmetric matrix. 
			double e[9];
			e[0] =  0.0; // E(1,1)
			e[1] =  pn.z; // E(2,1)
			e[2] = -pn.y; // E(3,1)
			e[3] = -pn.z; // E(1,2)
			e[4] =  0.0; // E(2,2)
			e[5] =  pn.x; // E(3,2)
			e[6] =  pn.y; // E(1,3)
			e[7] = -pn.x; // E(2,3)
			e[8] =  0.0; // E(3,3)

			double esq[9];
			esq[0] = e[0]*e[0] + e[3]*e[1] + e[6]*e[2];
			esq[1] = e[1]*e[0] + e[4]*e[1] + e[7]*e[2];
			esq[2] = e[2]*e[0] + e[5]*e[1] + e[8]*e[2];
			esq[3] = e[0]*e[3] + e[3]*e[4] + e[6]*e[5];
			esq[4] = e[1]*e[3] + e[4]*e[4] + e[7]*e[5];
			esq[5] = e[2]*e[3] + e[5]*e[4] + e[8]*e[5];
			esq[6] = e[0]*e[6] + e[3]*e[7] + e[6]*e[8];
			esq[7] = e[1]*e[6] + e[4]*e[7] + e[7]*e[8];
			esq[8] = e[2]*e[6] + e[5]*e[7] + e[8]*e[8];

			// Rotation matrix. 
			rmi[0] = 1.0 + (1.0 - ct)*esq[0] + st*e[0];
			rmi[1] = 0.0 + (1.0 - ct)*esq[1] + st*e[1];
			rmi[2] = 0.0 + (1.0 - ct)*esq[2] + st*e[2];
			rmi[3] = 0.0 + (1.0 - ct)*esq[3] + st*e[3];
			rmi[4] = 1.0 + (1.0 - ct)*esq[4] + st*e[4];
			rmi[5] = 0.0 + (1.0 - ct)*esq[5] + st*e[5];
			rmi[6] = 0.0 + (1.0 - ct)*esq[6] + st*e[6];
			rmi[7] = 0.0 + (1.0 - ct)*esq[7] + st*e[7];
			rmi[8] = 1.0 + (1.0 - ct)*esq[8] + st*e[8];

			if (shouldFindDeriv) {
				double de1[9];
				de1[0] =  0.0;
				de1[1] = -p.z*p.x*theta3;
				de1[2] =  p.y*p.x*theta3;
				de1[3] =  p.z*p.x*theta3;
				de1[4] =  0.0;
				de1[5] =  (theta - p.x*p.x*theta3);
				de1[6] = -p.y*p.x*theta3;
				de1[7] = -(theta - p.x*p.x*theta3);
				de1[8] =  0.0;

				double de2[9];
				de2[0] =  0.0;
				de2[1] = -p.z*p.y*theta3;
				de2[2] = -(theta - p.y*p.y*theta3);
				de2[3] =  p.z*p.y*theta3;
				de2[4] =  0.0;
				de2[5] = -p.x*p.y*theta3;
				de2[6] =  (theta - p.y*p.y*theta3);
				de2[7] =  p.x*p.y*theta3;
				de2[8] =  0.0;

				double de3[9];
				de3[0] =  0.0;
				de3[1] =  (theta - p.z*p.z*theta3);
				de3[2] =  p.y*p.z*theta3;
				de3[3] = -(theta - p.z*p.z*theta3);
				de3[4] =  0.0;
				de3[5] = -p.x*p.z*theta3;
				de3[6] = -p.y*p.z*theta3;
				de3[7] =  p.x*p.z*theta3;
				de3[8] =  0.0;

				drmi1[0] = st*pn.x*esq[0] + (1.0 - ct) * 
					((de1[0]*e[0] + de1[3]*e[1] + de1[6]*e[2]) + (e[0]*de1[0] + e[3]*de1[1] + e[6]*de1[2])) + 
					ct*pn.x*e[0] + st*de1[0];
				drmi1[1] = st*pn.x*esq[1] + (1.0 - ct) * 
					((de1[1]*e[0] + de1[4]*e[1] + de1[7]*e[2]) + (e[1]*de1[0] + e[4]*de1[1] + e[7]*de1[2])) + 
					ct*pn.x*e[1] + st*de1[1];
				drmi1[2] = st*pn.x*esq[2] + (1.0 - ct) * 
					((de1[2]*e[0] + de1[5]*e[1] + de1[8]*e[2]) + (e[2]*de1[0] + e[5]*de1[1] + e[8]*de1[2])) + 
					ct*pn.x*e[2] + st*de1[2];
				drmi1[3] = st*pn.x*esq[3] + (1.0 - ct) * 
					((de1[0]*e[3] + de1[3]*e[4] + de1[6]*e[5]) + (e[0]*de1[3] + e[3]*de1[4] + e[6]*de1[5])) + 
					ct*pn.x*e[3] + st*de1[3];
				drmi1[4] = st*pn.x*esq[4] + (1.0 - ct) * 
					((de1[1]*e[3] + de1[4]*e[4] + de1[7]*e[5]) + (e[1]*de1[3] + e[4]*de1[4] + e[7]*de1[5])) + 
					ct*pn.x*e[4] + st*de1[4];
				drmi1[5] = st*pn.x*esq[5] + (1.0 - ct) * 
					((de1[2]*e[3] + de1[5]*e[4] + de1[8]*e[5]) + (e[2]*de1[3] + e[5]*de1[4] + e[8]*de1[5])) + 
					ct*pn.x*e[5] + st*de1[5];
				drmi1[6] = st*pn.x*esq[6] + (1.0 - ct) * 
					((de1[0]*e[6] + de1[3]*e[7] + de1[6]*e[8]) + (e[0]*de1[6] + e[3]*de1[7] + e[6]*de1[8])) + 
					ct*pn.x*e[6] + st*de1[6];
				drmi1[7] = st*pn.x*esq[7] + (1.0 - ct) * 
					((de1[1]*e[6] + de1[4]*e[7] + de1[7]*e[8]) + (e[1]*de1[6] + e[4]*de1[7] + e[7]*de1[8])) + 
					ct*pn.x*e[7] + st*de1[7];
				drmi1[8] = st*pn.x*esq[8] + (1.0 - ct) * 
					((de1[2]*e[6] + de1[5]*e[7] + de1[8]*e[8]) + (e[2]*de1[6] + e[5]*de1[7] + e[8]*de1[8])) + 
					ct*pn.x*e[8] + st*de1[8];

				drmi2[0] = st*pn.y*esq[0] + (1.0 - ct) * 
					((de2[0]*e[0] + de2[3]*e[1] + de2[6]*e[2]) + (e[0]*de2[0] + e[3]*de2[1] + e[6]*de2[2])) + 
					ct*pn.y*e[0] + st*de2[0];
				drmi2[1] = st*pn.y*esq[1] + (1.0 - ct) * 
					((de2[1]*e[0] + de2[4]*e[1] + de2[7]*e[2]) + (e[1]*de2[0] + e[4]*de2[1] + e[7]*de2[2])) + 
					ct*pn.y*e[1] + st*de2[1];
				drmi2[2] = st*pn.y*esq[2] + (1.0 - ct) * 
					((de2[2]*e[0] + de2[5]*e[1] + de2[8]*e[2]) + (e[2]*de2[0] + e[5]*de2[1] + e[8]*de2[2])) + 
					ct*pn.y*e[2] + st*de2[2];
				drmi2[3] = st*pn.y*esq[3] + (1.0 - ct) * 
					((de2[0]*e[3] + de2[3]*e[4] + de2[6]*e[5]) + (e[0]*de2[3] + e[3]*de2[4] + e[6]*de2[5])) + 
					ct*pn.y*e[3] + st*de2[3];
				drmi2[4] = st*pn.y*esq[4] + (1.0 - ct) * 
					((de2[1]*e[3] + de2[4]*e[4] + de2[7]*e[5]) + (e[1]*de2[3] + e[4]*de2[4] + e[7]*de2[5])) + 
					ct*pn.y*e[4] + st*de2[4];
				drmi2[5] = st*pn.y*esq[5] + (1.0 - ct) * 
					((de2[2]*e[3] + de2[5]*e[4] + de2[8]*e[5]) + (e[2]*de2[3] + e[5]*de2[4] + e[8]*de2[5])) + 
					ct*pn.y*e[5] + st*de2[5];
				drmi2[6] = st*pn.y*esq[6] + (1.0 - ct) * 
					((de2[0]*e[6] + de2[3]*e[7] + de2[6]*e[8]) + (e[0]*de2[6] + e[3]*de2[7] + e[6]*de2[8])) + 
					ct*pn.y*e[6] + st*de2[6];
				drmi2[7] = st*pn.y*esq[7] + (1.0 - ct) * 
					((de2[1]*e[6] + de2[4]*e[7] + de2[7]*e[8]) + (e[1]*de2[6] + e[4]*de2[7] + e[7]*de2[8])) + 
					ct*pn.y*e[7] + st*de2[7];
				drmi2[8] = st*pn.y*esq[8] + (1.0 - ct) * 
					((de2[2]*e[6] + de2[5]*e[7] + de2[8]*e[8]) + (e[2]*de2[6] + e[5]*de2[7] + e[8]*de2[8])) + 
					ct*pn.y*e[8] + st*de2[8];

				drmi3[0] = st*pn.z*esq[0] + (1.0 - ct) * 
					((de3[0]*e[0] + de3[3]*e[1] + de3[6]*e[2]) + (e[0]*de3[0] + e[3]*de3[1] + e[6]*de3[2])) + 
					ct*pn.z*e[0] + st*de3[0];
				drmi3[1] = st*pn.z*esq[1] + (1.0 - ct) * 
					((de3[1]*e[0] + de3[4]*e[1] + de3[7]*e[2]) + (e[1]*de3[0] + e[4]*de3[1] + e[7]*de3[2])) + 
					ct*pn.z*e[1] + st*de3[1];
				drmi3[2] = st*pn.z*esq[2] + (1.0 - ct) * 
					((de3[2]*e[0] + de3[5]*e[1] + de3[8]*e[2]) + (e[2]*de3[0] + e[5]*de3[1] + e[8]*de3[2])) + 
					ct*pn.z*e[2] + st*de3[2];
				drmi3[3] = st*pn.z*esq[3] + (1.0 - ct) * 
					((de3[0]*e[3] + de3[3]*e[4] + de3[6]*e[5]) + (e[0]*de3[3] + e[3]*de3[4] + e[6]*de3[5])) + 
					ct*pn.z*e[3] + st*de3[3];
				drmi3[4] = st*pn.z*esq[4] + (1.0 - ct) * 
					((de3[1]*e[3] + de3[4]*e[4] + de3[7]*e[5]) + (e[1]*de3[3] + e[4]*de3[4] + e[7]*de3[5])) + 
					ct*pn.z*e[4] + st*de3[4];
				drmi3[5] = st*pn.z*esq[5] + (1.0 - ct) * 
					((de3[2]*e[3] + de3[5]*e[4] + de3[8]*e[5]) + (e[2]*de3[3] + e[5]*de3[4] + e[8]*de3[5])) + 
					ct*pn.z*e[5] + st*de3[5];
				drmi3[6] = st*pn.z*esq[6] + (1.0 - ct) * 
					((de3[0]*e[6] + de3[3]*e[7] + de3[6]*e[8]) + (e[0]*de3[6] + e[3]*de3[7] + e[6]*de3[8])) + 
					ct*pn.z*e[6] + st*de3[6];
				drmi3[7] = st*pn.z*esq[7] + (1.0 - ct) * 
					((de3[1]*e[6] + de3[4]*e[7] + de3[7]*e[8]) + (e[1]*de3[6] + e[4]*de3[7] + e[7]*de3[8])) + 
					ct*pn.z*e[7] + st*de3[7];
				drmi3[8] = st*pn.z*esq[8] + (1.0 - ct) * 
					((de3[2]*e[6] + de3[5]*e[7] + de3[8]*e[8]) + (e[2]*de3[6] + e[5]*de3[7] + e[8]*de3[8])) + 
					ct*pn.z*e[8] + st*de3[8];

			}
		}
	}

	__global__ void transform(const int *m_d_nRigidBody, const double *m_d_xRigid, double *d_grmi)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (*m_d_nRigidBody)) {

			double3 p; // Rotation vector. 

			p.x = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+0];
			p.y = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+1];
			p.z = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+2];

			double rmi[9]; // Rotation matrix. 

			double drmi1[9];
			double drmi2[9];
			double drmi3[9]; 

			// Calculate rotation matrices. 
			gpu_rigid_bodies::rmdrvt(p, rmi, drmi1, drmi2, drmi3, shouldFindDeriv);


			// Store in temporary array in global memory. 
			d_grmi[0+9*tid]= rmi[0];
			d_grmi[1+9*tid]= rmi[1];
			d_grmi[2+9*tid]= rmi[2];

			d_grmi[3+9*tid]= rmi[3];
			d_grmi[4+9*tid]= rmi[4];
			d_grmi[5+9*tid]= rmi[5];

			d_grmi[6+9*tid]= rmi[6];
			d_grmi[7+9*tid]= rmi[7];
			d_grmi[8+9*tid]= rmi[8];

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void transform2(double *d_x, const double *m_d_xRigid, const int *m_d_nRigidBody, 
			const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const double *m_d_sitesRigidBody, const double *d_grmi, const int *m_d_rigidMaxSite)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		int thisRigidBody = tid / (*m_d_rigidMaxSite);

		while (tid < (*m_d_nRigidBody) * (*m_d_rigidMaxSite)) {
			int i = tid - (thisRigidBody * (*m_d_rigidMaxSite));

			// Only use threads <= to *actual* number of sites. 
			if (i < m_d_nRigidSitesPerBody[thisRigidBody]) {
				int myAtom = m_d_rigidGroups[tid];

				d_x[3*myAtom-3] = m_d_xRigid[3*thisRigidBody] +  
					d_grmi[0+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi[3+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi[6+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];

				d_x[3*myAtom-2] = m_d_xRigid[3*thisRigidBody+1] + 
					d_grmi[1+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi[4+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi[7+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];

				d_x[3*myAtom-1] = m_d_xRigid[3*thisRigidBody+2] + 
					d_grmi[2+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi[5+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi[8+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void singleAtoms(double *d_x, const double *m_d_xRigid, const int *m_d_rigidSingles, const int *m_d_nRigidBody, 
			const int singlesThreads)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < singlesThreads) {
			int thisAtom = m_d_rigidSingles[tid];

			d_x[3*thisAtom-3] = m_d_xRigid[6*(*m_d_nRigidBody)+3*tid+0];
			d_x[3*thisAtom-2] = m_d_xRigid[6*(*m_d_nRigidBody)+3*tid+1];
			d_x[3*thisAtom-1] = m_d_xRigid[6*(*m_d_nRigidBody)+3*tid+2];

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void gradTransform1a(const double *m_d_xRigid, double *d_grmi1, double *d_grmi2, double *d_grmi3, 
			const int *m_d_nRigidBody)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (*m_d_nRigidBody)) {

			double3 pi; // Rotation vector. 

			pi.x = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+0];
			pi.y = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+1];
			pi.z = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+2];

			double drmi0[9];
			// Rotation matrices. 
			double drmi1[9];
			double drmi2[9];
			double drmi3[9];

			// Calculate rotation matrices. 
			gpu_rigid_bodies::rmdrvt(pi, drmi0, drmi1, drmi2, drmi3, shouldFindDeriv);

			// Store in temporary arrays in global memory. 

			d_grmi1[0+9*tid]= drmi1[0];
			d_grmi1[1+9*tid]= drmi1[1];
			d_grmi1[2+9*tid]= drmi1[2];

			d_grmi1[3+9*tid]= drmi1[3];
			d_grmi1[4+9*tid]= drmi1[4];
			d_grmi1[5+9*tid]= drmi1[5];

			d_grmi1[6+9*tid]= drmi1[6];
			d_grmi1[7+9*tid]= drmi1[7];
			d_grmi1[8+9*tid]= drmi1[8];


			d_grmi2[0+9*tid]= drmi2[0];
			d_grmi2[1+9*tid]= drmi2[1];
			d_grmi2[2+9*tid]= drmi2[2];

			d_grmi2[3+9*tid]= drmi2[3];
			d_grmi2[4+9*tid]= drmi2[4];
			d_grmi2[5+9*tid]= drmi2[5];

			d_grmi2[6+9*tid]= drmi2[6];
			d_grmi2[7+9*tid]= drmi2[7];
			d_grmi2[8+9*tid]= drmi2[8];


			d_grmi3[0+9*tid]= drmi3[0];
			d_grmi3[1+9*tid]= drmi3[1];
			d_grmi3[2+9*tid]= drmi3[2];

			d_grmi3[3+9*tid]= drmi3[3];
			d_grmi3[4+9*tid]= drmi3[4];
			d_grmi3[5+9*tid]= drmi3[5];

			d_grmi3[6+9*tid]= drmi3[6];
			d_grmi3[7+9*tid]= drmi3[7];
			d_grmi3[8+9*tid]= drmi3[8];

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void gradTransform1b(const int *m_d_nRigidBody, const double *m_d_xRigid, double *d_grmi0)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (*m_d_nRigidBody)) {

			double3 pi;
			pi.x = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+0];
			pi.y = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+1];
			pi.z = m_d_xRigid[3*(*m_d_nRigidBody)+3*tid+2];

			double drmi0[9];
			double drmi1[9];
			double drmi2[9];
			double drmi3[9];

			gpu_rigid_bodies::rmdrvt(pi, drmi0, drmi1, drmi2, drmi3, shouldFindDeriv);

			d_grmi0[0+9*tid]= drmi0[0];
			d_grmi0[1+9*tid]= drmi0[1];
			d_grmi0[2+9*tid]= drmi0[2];

			d_grmi0[3+9*tid]= drmi0[3];
			d_grmi0[4+9*tid]= drmi0[4];
			d_grmi0[5+9*tid]= drmi0[5];

			d_grmi0[6+9*tid]= drmi0[6];
			d_grmi0[7+9*tid]= drmi0[7];
			d_grmi0[8+9*tid]= drmi0[8];

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void intermediate(const double *d_gk, const int *m_d_nRigidSitesPerBody, const double *d_grmi1, 
			const double *d_grmi2, const double *d_grmi3, const double *m_d_sitesRigidBody, 
			const int *m_d_rigidGroups, double *d_tempArray, const int *m_d_nRigidBody, 
			const int *m_d_rigidMaxSite)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		int thisRigidBody = tid / (*m_d_rigidMaxSite);

		while (tid < (*m_d_nRigidBody) * (*m_d_rigidMaxSite)) {

			// Initialise entire array to zeros. 
			d_tempArray[3*tid+2] = 0.0;
			d_tempArray[3*tid+1] = 0.0;
			d_tempArray[3*tid+0] = 0.0;

			int i = tid - (thisRigidBody * (*m_d_rigidMaxSite));

			// Only use threads <= to *actual* number of sites. 
			if (i < m_d_nRigidSitesPerBody[thisRigidBody]) {
				int myAtom = m_d_rigidGroups[tid];

				double3 dr1;
				double3 dr2;
				double3 dr3;

				// Matrix multiplications between rotation matrices and rigid body sites. 
				dr1.x = d_grmi1[0+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi1[3+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi1[6+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr1.y = d_grmi1[1+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi1[4+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi1[7+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr1.z = d_grmi1[2+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi1[5+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi1[8+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];

				dr2.x = d_grmi2[0+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi2[3+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi2[6+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr2.y = d_grmi2[1+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi2[4+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi2[7+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr2.z = d_grmi2[2+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi2[5+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi2[8+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];

				dr3.x = d_grmi3[0+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi3[3+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi3[6+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr3.y = d_grmi3[1+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi3[4+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi3[7+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr3.z = d_grmi3[2+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi3[5+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi3[8+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];

				// Store dot product with gradient in temporary array in global memory. 
				d_tempArray[3*i+(*m_d_rigidMaxSite)*3*thisRigidBody] = d_gk[3*myAtom-3]*dr1.x + 
					d_gk[3*myAtom-2]*dr1.y + 
					d_gk[3*myAtom-1]*dr1.z;
				d_tempArray[1+3*i+(*m_d_rigidMaxSite)*3*thisRigidBody] = d_gk[3*myAtom-3]*dr2.x + 
					d_gk[3*myAtom-2]*dr2.y + 
					d_gk[3*myAtom-1]*dr2.z;
				d_tempArray[2+3*i+(*m_d_rigidMaxSite)*3*thisRigidBody] = d_gk[3*myAtom-3]*dr3.x + 
					d_gk[3*myAtom-2]*dr3.y + 
					d_gk[3*myAtom-1]*dr3.z;

			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void gradSingleAtoms(const double *d_gk, double *m_d_gkRigid, const int *m_d_rigidSingles, 
			const int *m_d_nRigidBody, const int singlesThreads)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < singlesThreads) {
			int thisAtom = m_d_rigidSingles[tid];

			m_d_gkRigid[6*(*m_d_nRigidBody)+3*tid+0] = d_gk[3*thisAtom-3];
			m_d_gkRigid[6*(*m_d_nRigidBody)+3*tid+1] = d_gk[3*thisAtom-2];
			m_d_gkRigid[6*(*m_d_nRigidBody)+3*tid+2] = d_gk[3*thisAtom-1];

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void aaConvRmdrvt(double *d_rmi10, double *d_rmi20, double *d_rmi30)
	{
		double3 p; // Rotation vector. 

		p.x = 0.0;
		p.y = 0.0;
		p.z = 0.0;

		// Rotation matrices. 
		double rmi0[9];
		double drmi10[9];
		double drmi20[9];
		double drmi30[9];

		// Calculate rotation matrices. 
		gpu_rigid_bodies::rmdrvt(p, rmi0, drmi10, drmi20, drmi30, shouldFindDeriv);

		// Store in temporary arrays in global memory. 

		d_rmi10[0] = drmi10[0];
		d_rmi10[1] = drmi10[1];
		d_rmi10[2] = drmi10[2];
		d_rmi10[3] = drmi10[3];
		d_rmi10[4] = drmi10[4];
		d_rmi10[5] = drmi10[5];
		d_rmi10[6] = drmi10[6];
		d_rmi10[7] = drmi10[7];
		d_rmi10[8] = drmi10[8];

		d_rmi20[0] = drmi20[0];
		d_rmi20[1] = drmi20[1];
		d_rmi20[2] = drmi20[2];
		d_rmi20[3] = drmi20[3];
		d_rmi20[4] = drmi20[4];
		d_rmi20[5] = drmi20[5];
		d_rmi20[6] = drmi20[6];
		d_rmi20[7] = drmi20[7];
		d_rmi20[8] = drmi20[8];

		d_rmi30[0] = drmi30[0];
		d_rmi30[1] = drmi30[1];
		d_rmi30[2] = drmi30[2];
		d_rmi30[3] = drmi30[3];
		d_rmi30[4] = drmi30[4];
		d_rmi30[5] = drmi30[5];
		d_rmi30[6] = drmi30[6];
		d_rmi30[7] = drmi30[7];
		d_rmi30[8] = drmi30[8];

	}

	__global__ void intermediate2(const double *d_gk, const int *m_d_nRigidSitesPerBody, const double *d_grmi0, 
			const double *d_grmi10, const double *d_grmi20, const double *d_grmi30, 
			const double *m_d_sitesRigidBody, const int *m_d_rigidGroups, double *d_tempArray, 
			const int *m_d_nRigidBody, const int *m_d_rigidMaxSite)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		int thisRigidBody = tid / (*m_d_rigidMaxSite);

		while (tid < (*m_d_nRigidBody) * (*m_d_rigidMaxSite)) {

			d_tempArray[3*tid+2] = 0.0;
			d_tempArray[3*tid+1] = 0.0;
			d_tempArray[3*tid+0] = 0.0;

			int i = tid - (thisRigidBody * (*m_d_rigidMaxSite));

			if (i < m_d_nRigidSitesPerBody[thisRigidBody]) {
				int myAtom = m_d_rigidGroups[tid];

				double3 dr0;
				double3 dr1;
				double3 dr2;
				double3 dr3;

				dr0.x = d_grmi0[0+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi0[3+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi0[6+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr0.y = d_grmi0[1+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi0[4+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi0[7+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];
				dr0.z = d_grmi0[2+9*thisRigidBody]*
					m_d_sitesRigidBody[i+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi0[5+9*thisRigidBody]*
					m_d_sitesRigidBody[i+(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)] + 
					d_grmi0[8+9*thisRigidBody]*
					m_d_sitesRigidBody[i+2*(*m_d_rigidMaxSite)+thisRigidBody*3*(*m_d_rigidMaxSite)];

				dr1.x = d_grmi10[0]*dr0.x + d_grmi10[3]*dr0.y + d_grmi10[6]*dr0.z;
				dr1.y = d_grmi10[1]*dr0.x + d_grmi10[4]*dr0.y + d_grmi10[7]*dr0.z;
				dr1.z = d_grmi10[2]*dr0.x + d_grmi10[5]*dr0.y + d_grmi10[8]*dr0.z;

				dr2.x = d_grmi20[0]*dr0.x + d_grmi20[3]*dr0.y + d_grmi20[6]*dr0.z;
				dr2.y = d_grmi20[1]*dr0.x + d_grmi20[4]*dr0.y + d_grmi20[7]*dr0.z;
				dr2.z = d_grmi20[2]*dr0.x + d_grmi20[5]*dr0.y + d_grmi20[8]*dr0.z;

				dr3.x = d_grmi30[0]*dr0.x + d_grmi30[3]*dr0.y + d_grmi30[6]*dr0.z;
				dr3.y = d_grmi30[1]*dr0.x + d_grmi30[4]*dr0.y + d_grmi30[7]*dr0.z;
				dr3.z = d_grmi30[2]*dr0.x + d_grmi30[5]*dr0.y + d_grmi30[8]*dr0.z;


				d_tempArray[0+3*i+(*m_d_rigidMaxSite)*3*thisRigidBody] = d_gk[3*myAtom-3]*dr1.x + 
					d_gk[3*myAtom-2]*dr1.y + 
					d_gk[3*myAtom-1]*dr1.z;
				d_tempArray[1+3*i+(*m_d_rigidMaxSite)*3*thisRigidBody] = d_gk[3*myAtom-3]*dr2.x + 
					d_gk[3*myAtom-2]*dr2.y + 
					d_gk[3*myAtom-1]*dr2.z;
				d_tempArray[2+3*i+(*m_d_rigidMaxSite)*3*thisRigidBody] = d_gk[3*myAtom-3]*dr3.x + 
					d_gk[3*myAtom-2]*dr3.y + 
					d_gk[3*myAtom-1]*dr3.z;

			}

			tid += blockDim.x * gridDim.x;
		}
	}


	__global__ void aaConvTorque(const double *d_torques, const double *m_d_gkRigid, const double *d_grmi0, 
			const int *m_d_nRigidSitesPerBody, const double *m_d_rigidInverse, double *d_rmsArray, 
			const int *m_d_nRigidBody)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (*m_d_nRigidBody)) {

			double3 torque;

			torque.x = d_grmi0[0+9*tid]*d_torques[3*tid+0] + d_grmi0[1+9*tid]*d_torques[3*tid+1] + 
				d_grmi0[2+9*tid]*d_torques[3*tid+2];
			torque.y = d_grmi0[3+9*tid]*d_torques[3*tid+0] + d_grmi0[4+9*tid]*d_torques[3*tid+1] + 
				d_grmi0[5+9*tid]*d_torques[3*tid+2];
			torque.z = d_grmi0[6+9*tid]*d_torques[3*tid+0] + d_grmi0[7+9*tid]*d_torques[3*tid+1] + 
				d_grmi0[8+9*tid]*d_torques[3*tid+2];

			double3 matmul;

			matmul.x =  m_d_rigidInverse[tid+0*(*m_d_nRigidBody)]*torque.x + 
				m_d_rigidInverse[tid+1*(*m_d_nRigidBody)]*torque.y + 
				m_d_rigidInverse[tid+2*(*m_d_nRigidBody)]*torque.z;
			matmul.y =  m_d_rigidInverse[tid+3*(*m_d_nRigidBody)]*torque.x + 
				m_d_rigidInverse[tid+4*(*m_d_nRigidBody)]*torque.y + 
				m_d_rigidInverse[tid+5*(*m_d_nRigidBody)]*torque.z;
			matmul.z =  m_d_rigidInverse[tid+6*(*m_d_nRigidBody)]*torque.x + 
				m_d_rigidInverse[tid+7*(*m_d_nRigidBody)]*torque.y + 
				m_d_rigidInverse[tid+8*(*m_d_nRigidBody)]*torque.z;

			double rmscontrib1 = matmul.x*torque.x + matmul.y*torque.y + matmul.z*torque.z;

			double rmscontrib2 = (1.0/m_d_nRigidSitesPerBody[tid])*(m_d_gkRigid[3*tid+0]*m_d_gkRigid[3*tid+0] + 
					m_d_gkRigid[3*tid+1]*m_d_gkRigid[3*tid+1] + 
					m_d_gkRigid[3*tid+2]*m_d_gkRigid[3*tid+2]);


			d_rmsArray[tid] = rmscontrib1 + rmscontrib2;

			tid += blockDim.x * gridDim.x;
		}
	}

	// Undocumented CUDA 6.5 function copied and pasted here as not present in CUDA 6.0 and below. 
	static __device__ __inline__ 
		double __myshfl_down(double var, unsigned int delta, int width=warpSize) {
			float lo, hi;
			asm volatile("mov.b64 {%0,%1}, %2;" : "=f"(lo), "=f"(hi) : "d"(var));
			hi = __shfl_down(hi, delta, width);
			lo = __shfl_down(lo, delta, width);
			asm volatile("mov.b64 %0, {%1,%2};" : "=d"(var) : "f"(lo), "f"(hi));
			return var;
		}

	// See 'Faster Parallel Reductions on Kepler' from the NVIDIA dev blog. 
	__inline__ __device__
		double warpReduceSum(double val) {
			for (int offset = warpSize/2; offset > 0; offset /= 2) {
				val += __myshfl_down(val, offset);
			}
			return val;
		}

	__global__ void warpReduce1(const int *m_d_nRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const int *m_d_rigidMaxSite, const double *d_gk, double *m_d_gkRigid)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (warpSize*(*m_d_nRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + warpSize) % warpSize);
			int currentBody = tid/warpSize;

			if ((index < m_d_nRigidSitesPerBody[currentBody]) && (m_d_nRigidSitesPerBody[currentBody] <= warpSize)) {
				int myAtom = m_d_rigidGroups[index+(*m_d_rigidMaxSite)*currentBody];

				elements.x = d_gk[3*myAtom-3];
				elements.y = d_gk[3*myAtom-2];
				elements.z = d_gk[3*myAtom-1];
			}

			elements.x = warpReduceSum(elements.x);
			elements.y = warpReduceSum(elements.y);
			elements.z = warpReduceSum(elements.z);

			if (tid % warpSize == 0) {
				m_d_gkRigid[3*currentBody] = elements.x;
				m_d_gkRigid[3*currentBody+1] = elements.y;
				m_d_gkRigid[3*currentBody+2] = elements.z;
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void warpReduce2(const int *m_d_nRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidMaxSite, 
			double *m_d_gkRigid, const double *d_tempArray)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (warpSize*(*m_d_nRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + warpSize) % warpSize);
			int currentBody = tid/warpSize;

			if ((index < m_d_nRigidSitesPerBody[currentBody]) && (m_d_nRigidSitesPerBody[currentBody] <= warpSize)) {
				elements.x = d_tempArray[0+3*index+3*(*m_d_rigidMaxSite)*currentBody];
				elements.y = d_tempArray[1+3*index+3*(*m_d_rigidMaxSite)*currentBody];
				elements.z = d_tempArray[2+3*index+3*(*m_d_rigidMaxSite)*currentBody];
			}

			elements.x = warpReduceSum(elements.x);
			elements.y = warpReduceSum(elements.y);
			elements.z = warpReduceSum(elements.z);

			if (tid % warpSize == 0) {
				m_d_gkRigid[3*(*m_d_nRigidBody)+3*currentBody] = elements.x;
				m_d_gkRigid[3*(*m_d_nRigidBody)+3*currentBody+1] = elements.y;
				m_d_gkRigid[3*(*m_d_nRigidBody)+3*currentBody+2] = elements.z;
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void warpReduce3(const int *m_d_nRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidMaxSite, 
			double *d_torques, const double *d_tempArray)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (warpSize*(*m_d_nRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + warpSize) % warpSize);
			int currentBody = tid/warpSize;

			if ((index < m_d_nRigidSitesPerBody[currentBody]) && (m_d_nRigidSitesPerBody[currentBody] <= warpSize)) {
				elements.x = d_tempArray[0+3*index+3*(*m_d_rigidMaxSite)*currentBody];
				elements.y = d_tempArray[1+3*index+3*(*m_d_rigidMaxSite)*currentBody];
				elements.z = d_tempArray[2+3*index+3*(*m_d_rigidMaxSite)*currentBody];
			}

			elements.x = warpReduceSum(elements.x);
			elements.y = warpReduceSum(elements.y);
			elements.z = warpReduceSum(elements.z);

			if (tid % warpSize == 0) {
				d_torques[3*currentBody] = elements.x;
				d_torques[3*currentBody+1] = elements.y;
				d_torques[3*currentBody+2] = elements.z;
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__inline__ __device__
		double blockReduceSum(double val) {
			static __shared__ double shared[32]; // Shared mem for 32 partial sums
			int lane = threadIdx.x % warpSize;
			int wid = threadIdx.x / warpSize;

			val = warpReduceSum(val);     // Each warp performs partial reduction

			// Write reduced value to shared memory
			if (lane==0) {
				shared[wid]=val;
			} 

			__syncthreads();              // Wait for all partial reductions

			//read from shared memory only if that warp existed
			val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

			//Final reduce within first warp
			if (wid==0) {
				val = warpReduceSum(val);
			} 

			return val;
		}

	__global__ void fullReduce1a(const int *m_d_nLargeRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const int *m_d_rigidMaxSite, const double *d_gk, const int *m_d_largeRigidIndices, double *d_outArray, 
			const int roundedMaxSite) 
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (roundedMaxSite*(*m_d_nLargeRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + roundedMaxSite) % roundedMaxSite);
			int thisBody = tid/roundedMaxSite;
			int currentBody = m_d_largeRigidIndices[thisBody];

			if (index < m_d_nRigidSitesPerBody[currentBody]) {
				int myAtom = m_d_rigidGroups[index+(*m_d_rigidMaxSite)*currentBody];
				elements.x = d_gk[3*myAtom-3];
				elements.y = d_gk[3*myAtom-2];
				elements.z = d_gk[3*myAtom-1];
			}

			elements.x = blockReduceSum(elements.x);
			__syncthreads();
			elements.y = blockReduceSum(elements.y);
			__syncthreads();
			elements.z = blockReduceSum(elements.z);

			if (threadIdx.x==0) {
				d_outArray[3*blockIdx.x+0] = elements.x;
				d_outArray[3*blockIdx.x+1] = elements.y;
				d_outArray[3*blockIdx.x+2] = elements.z;
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void fullReduce1b(const int *m_d_nLargeRigidBody, const int *m_d_nRigidSitesPerBody, const int *m_d_rigidGroups, 
			const int *m_d_rigidMaxSite, const int *m_d_largeRigidIndices, double *d_outArray, 
			const double *d_tempArray, const int roundedMaxSite) 
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (roundedMaxSite*(*m_d_nLargeRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + roundedMaxSite) % roundedMaxSite);
			int thisBody = tid/roundedMaxSite;
			int currentBody = m_d_largeRigidIndices[thisBody];

			if (index < m_d_nRigidSitesPerBody[currentBody]) {
				elements.x = d_tempArray[0+3*index+3*(*m_d_rigidMaxSite)*currentBody];
				elements.y = d_tempArray[1+3*index+3*(*m_d_rigidMaxSite)*currentBody];
				elements.z = d_tempArray[2+3*index+3*(*m_d_rigidMaxSite)*currentBody];
			}

			elements.x = blockReduceSum(elements.x);
			__syncthreads();
			elements.y = blockReduceSum(elements.y);
			__syncthreads();
			elements.z = blockReduceSum(elements.z);

			if (threadIdx.x==0) {
				d_outArray[3*blockIdx.x+0] = elements.x;
				d_outArray[3*blockIdx.x+1] = elements.y;
				d_outArray[3*blockIdx.x+2] = elements.z;
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void fullReduce2a(const int *m_d_nLargeRigidBody, double *m_d_gkRigid, const int *m_d_largeRigidIndices, 
			double *d_outArray, const int roundedMaxSite, const int outSize, const int blocks) 
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (roundedMaxSite*(*m_d_nLargeRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + roundedMaxSite) % roundedMaxSite);
			int thisBody = tid/roundedMaxSite;
			int currentBody = m_d_largeRigidIndices[thisBody];

			int sectionSize = outSize/(*m_d_nLargeRigidBody);
			if (index < sectionSize) {
				int var = index+thisBody*sectionSize;
				elements.x = d_outArray[3*var+0];
				elements.y = d_outArray[3*var+1];
				elements.z = d_outArray[3*var+2];
			}

			elements.x = blockReduceSum(elements.x);
			__syncthreads();
			elements.y = blockReduceSum(elements.y);
			__syncthreads();
			elements.z = blockReduceSum(elements.z);

			if (threadIdx.x==0) {
				if (blocks == (*m_d_nLargeRigidBody)) {
					m_d_gkRigid[3*currentBody] = elements.x;
					m_d_gkRigid[3*currentBody+1] = elements.y;
					m_d_gkRigid[3*currentBody+2] = elements.z;
				}
				else {
					d_outArray[3*blockIdx.x+0] = elements.x;
					d_outArray[3*blockIdx.x+1] = elements.y;
					d_outArray[3*blockIdx.x+2] = elements.z;
				}
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void fullReduce2b(const int *m_d_nLargeRigidBody, double *m_d_gkRigid, const int *m_d_largeRigidIndices, 
			double *d_outArray, const int *m_d_nRigidBody, const int roundedMaxSite, const int outSize, const int blocks) 
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (roundedMaxSite*(*m_d_nLargeRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + roundedMaxSite) % roundedMaxSite);
			int thisBody = tid/roundedMaxSite;
			int currentBody = m_d_largeRigidIndices[thisBody];

			int sectionSize = outSize/(*m_d_nLargeRigidBody);
			if (index < sectionSize) {
				int var = index+thisBody*sectionSize;
				elements.x = d_outArray[3*var+0];
				elements.y = d_outArray[3*var+1];
				elements.z = d_outArray[3*var+2];
			}

			elements.x = blockReduceSum(elements.x);
			__syncthreads();
			elements.y = blockReduceSum(elements.y);
			__syncthreads();
			elements.z = blockReduceSum(elements.z);

			if (threadIdx.x==0) {
				if (blocks == (*m_d_nLargeRigidBody)) {
					m_d_gkRigid[3*(*m_d_nRigidBody)+3*currentBody] = elements.x;
					m_d_gkRigid[3*(*m_d_nRigidBody)+3*currentBody+1] = elements.y;
					m_d_gkRigid[3*(*m_d_nRigidBody)+3*currentBody+2] = elements.z;
				}
				else {
					d_outArray[3*blockIdx.x+0] = elements.x;
					d_outArray[3*blockIdx.x+1] = elements.y;
					d_outArray[3*blockIdx.x+2] = elements.z;
				}
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void fullReduce2c(const int *m_d_nLargeRigidBody, double *d_torques, const int *m_d_largeRigidIndices, 
			double *d_outArray, const int *m_d_nRigidBody, const int roundedMaxSite, const int outSize, 
			const int blocks) 
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (roundedMaxSite*(*m_d_nLargeRigidBody))) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + roundedMaxSite) % roundedMaxSite); // Technically unnecessary with one body at a time - equivalent to tid
			int thisBody = tid/roundedMaxSite; // Also unnecessary with one body - could just pass as argument
			int currentBody = m_d_largeRigidIndices[thisBody];

			int sectionSize = outSize/(*m_d_nLargeRigidBody);
			if (index < sectionSize) {
				int var = index+thisBody*sectionSize;
				elements.x = d_outArray[3*var+0];
				elements.y = d_outArray[3*var+1];
				elements.z = d_outArray[3*var+2];
			}

			elements.x = blockReduceSum(elements.x);
			__syncthreads();
			elements.y = blockReduceSum(elements.y);
			__syncthreads();
			elements.z = blockReduceSum(elements.z);

			if (threadIdx.x==0) {
				if (blocks == (*m_d_nLargeRigidBody)) {
					d_torques[3*currentBody] = elements.x;
					d_torques[3*currentBody+1] = elements.y;
					d_torques[3*currentBody+2] = elements.z;
				}
				else {
					d_outArray[3*blockIdx.x+0] = elements.x;
					d_outArray[3*blockIdx.x+1] = elements.y;
					d_outArray[3*blockIdx.x+2] = elements.z;
				}
			}

			tid += blockDim.x * gridDim.x;
		}

	}

	__global__ void fullReduce(const double *in, double* out, const int N) 
	{
		double sum = 0;
		//reduce multiple elements per thread
		for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x * gridDim.x) {
			sum += in[i];
		}
		sum = blockReduceSum(sum);
		if (threadIdx.x==0) {
			out[blockIdx.x]=sum;
		}
	}

	__global__ void fullDotReduce(const double *in, double* out, const int N)
	{
		double sum = 0;
		//reduce multiple elements per thread
		for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x * gridDim.x) {
			sum += in[i]*in[i];
		}
		sum = blockReduceSum(sum);
		if (threadIdx.x==0) {
			out[blockIdx.x]=sum;
		}
	}

}
