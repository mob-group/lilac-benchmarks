/**
 *
 * File lj_potential.cu: Implementation of class LjPotential, inherited from class CostFunction. 
 *
 **/

#include "cost_function.h"
#include "potential.h"

namespace gpu_lj
{
	// Variables on the GPU. 

	__constant__ size_t numDimensions;
	__constant__ int nAddTarget;


	// Kernels.

	__global__ void fullReduce1(const double *in, double *out, double *d_gradf, const int N, const int roundedMaxSite, 
			const int blocks);
	__global__ void kernelComputeEneGrad(const double *d_x, const double *m_d_ljAddRep, const double *m_d_ljAddAtt, 
			double *d_energyContrib, double *d_gradientContrib);
	__global__ void deviceReduceKernel(const double *in, double* out, const int N);
}



LjPotential::LjPotential(Printing &debugPrinting, Timer &timer_potential, Cublas &cublas, size_t numDimensions,
		int nDegFreedom, int nRigidBody, int rigidMaxSite, int *nRigidSitesPerBody,
		int *rigidGroups, double *sitesRigidBody, int *rigidSingles, double *rigidInverse,
		double *coords, int nSecDiag, bool isAtomisticNotRigid, double aaConvThreshold,
		double coldFusionLim, bool shouldFreeze, bool *isAtomFrozen, int nFreeze, bool isAaConvergence, int nAddTarget,
		double *ljAddRep, double *ljAddAtt)
: CostFunction(debugPrinting, timer_potential, cublas, numDimensions, nDegFreedom, nRigidBody, rigidMaxSite,
		nRigidSitesPerBody, rigidGroups, sitesRigidBody, rigidSingles, rigidInverse, coords,
		nSecDiag, isAtomisticNotRigid, aaConvThreshold, coldFusionLim, shouldFreeze, isAtomFrozen,
		nFreeze, isAaConvergence, nAddTarget, ljAddRep, ljAddAtt)
{
	CudaSafeCall( cudaMemcpyToSymbol(gpu_lj::numDimensions, &m_numDimensions, sizeof(size_t)) );
	CudaSafeCall( cudaMemcpyToSymbol(gpu_lj::nAddTarget, &m_nAddTarget, sizeof(int)) );
}



void LjPotential::computeEnergyAndGradient(const double *d_x, double *d_f, double *d_gradf)
{
	using namespace gpu_lj;

	dim3 blockDim;
	dim3 gridDim;

	// Set launch parameters for main potential kernel. 
	int numThreads = (m_numDimensions/3)*(m_numDimensions/3);
	blockDim.x = 512;
	int numBlocks = (numThreads + blockDim.x - 1)/blockDim.x;
	gridDim.x = numBlocks;

	double *d_energyContrib;
	double *d_gradientContrib;

	CudaSafeCall( cudaMalloc(&d_energyContrib, numThreads * sizeof(double)) );
	CudaSafeCall( cudaMalloc(&d_gradientContrib, 3 * numThreads * sizeof(double)) );

	gpu_lj::kernelComputeEneGrad<<<gridDim, blockDim>>>(d_x, m_d_ljAddRep, m_d_ljAddAtt, d_energyContrib, d_gradientContrib);
	CudaCheckError();
	cudaDeviceSynchronize();


	// Summation of energy contributions
	int threads = 512;
	int blocks = min((numThreads + threads - 1) / threads, 1024);

	double *d_out;

	CudaSafeCall( cudaMalloc(&d_out, blocks * sizeof(double)) );

	gpu_lj::deviceReduceKernel<<<blocks, threads>>>(d_energyContrib, d_out, numThreads);
	CudaCheckError();
	cudaDeviceSynchronize();

	gpu_lj::deviceReduceKernel<<<1, 1024>>>(d_out, d_out, blocks);
	CudaCheckError();
	cudaDeviceSynchronize();

	CudaSafeCall( cudaMemcpy(d_f, d_out, sizeof(double), cudaMemcpyDeviceToDevice) );


	// Summation of gradient contributions
	int blockSize = 1024;
	blockDim.x = blockSize;

	// Round m_numDimensions to nearest block size so that per block reduction doesn't sum components from different atoms. 
	int roundedMaxSite = ((m_numDimensions/3 + blockSize - 1)/blockSize)*blockSize;

	numThreads = roundedMaxSite*m_numDimensions/3;

	blocks = (numThreads + blockDim.x - 1)/blockDim.x;
	gridDim.x = blocks;

	double *d_outArray;

	CudaSafeCall( cudaMalloc(&d_outArray, 3 * blocks * sizeof(double)) );

	gpu_lj::fullReduce1<<<gridDim, blockDim>>>(d_gradientContrib, d_outArray, d_gradf, m_numDimensions/3, roundedMaxSite, 
			blocks);
	CudaCheckError();
	cudaDeviceSynchronize();

	while (blocks > m_numDimensions/3) {
		int oldBlocks = blocks;

		roundedMaxSite = ((blocks/(m_numDimensions/3) + blockSize - 1)/blockSize)*blockSize;
		numThreads = roundedMaxSite*m_numDimensions/3;

		blocks = (numThreads + blockDim.x - 1)/blockDim.x;
		gridDim.x = blocks;

		gpu_lj::fullReduce1<<<gridDim, blockDim>>>(d_outArray, d_outArray, d_gradf, oldBlocks/(m_numDimensions/3), 
				roundedMaxSite, blocks);
		CudaCheckError();
		cudaDeviceSynchronize();
	}

	CudaSafeCall( cudaFree(d_energyContrib) );
	CudaSafeCall( cudaFree(d_gradientContrib) );
	CudaSafeCall( cudaFree(d_out) );
	CudaSafeCall( cudaFree(d_outArray) );
}



namespace gpu_lj 
{
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

	__inline__ __device__
		double warpReduceSum(double val) {
			for (int offset = warpSize/2; offset > 0; offset /= 2) 
				val += __myshfl_down(val, offset);
			return val;
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

	__global__ void deviceReduceKernel(const double *in, double* out, const int N) {
		double sum = 0;
		//reduce multiple elements per thread.
		for (int i = blockIdx.x * blockDim.x + threadIdx.x;
				i < N;
				i += blockDim.x * gridDim.x) {
			sum += in[i];
		}
		sum = blockReduceSum(sum);
		if (threadIdx.x==0)
			out[blockIdx.x]=sum;
	}

	__global__ void fullReduce1(const double *in, double* out, double *d_gradf, const int N, const int roundedMaxSite, 
			const int blocks)
	{
		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (roundedMaxSite*numDimensions/3)) {

			double3 elements;
			elements.x = 0.0;
			elements.y = 0.0;
			elements.z = 0.0;

			int index = ((tid + roundedMaxSite) % roundedMaxSite);
			int thisAtom = tid/roundedMaxSite;

			if (index < N) {
				int var = index + thisAtom*N;
				elements.x = in[3*var+0];
				elements.y = in[3*var+1];
				elements.z = in[3*var+2];
			}

			elements.x = blockReduceSum(elements.x);
			__syncthreads();
			elements.y = blockReduceSum(elements.y);
			__syncthreads();
			elements.z = blockReduceSum(elements.z);

			if (threadIdx.x==0) {
				if (blocks == numDimensions/3) {
					d_gradf[3*thisAtom+0] = elements.x;
					d_gradf[3*thisAtom+1] = elements.y;
					d_gradf[3*thisAtom+2] = elements.z;
				}
				else {
					out[3*blockIdx.x+0] = elements.x;
					out[3*blockIdx.x+1] = elements.y;
					out[3*blockIdx.x+2] = elements.z;
				}
			}

			tid += blockDim.x * gridDim.x;
		}
	}

	__global__ void kernelComputeEneGrad(const double *d_x, const double *m_d_ljAddRep, const double *m_d_ljAddAtt, 
			double *d_energyContrib, double *d_gradientContrib)
	{
		// Each thread deals with one atom interacting with one other atom.

		int tid = blockIdx.x * blockDim.x + threadIdx.x;

		while (tid < (numDimensions/3)*(numDimensions/3)) {

			int refAtom = tid / (numDimensions/3); // Integer division rounds down.
			int myAtom = tid % (numDimensions/3);

			// Read the coordinates from global memory into memory local to each thread. 
			// Could reorder d_x to xxxyyyzzz pattern for coalesced access (in setup routine before L-BFGS)
			// Extra reads probably covered by L2 cache on Kepler though - would need to test
			double myPositionX = d_x[3*myAtom+0];
			double myPositionY = d_x[3*myAtom+1];
			double myPositionZ = d_x[3*myAtom+2];

			double refPositionX = d_x[3*refAtom+0];
			double refPositionY = d_x[3*refAtom+1];
			double refPositionZ = d_x[3*refAtom+2];

			int mj1 = refAtom % nAddTarget;
			int mj2 = myAtom % nAddTarget;

			// Carry out energy and gradient calculations. 

			// Calculate the distance vector. 
			double3 r;
			r.x = myPositionX - refPositionX;
			r.y = myPositionY - refPositionY;
			r.z = myPositionZ - refPositionZ;

			// Distance squared. 
			double rsq = r.x * r.x + r.y * r.y + r.z * r.z;

			double attractive = m_d_ljAddAtt[mj2+mj1*nAddTarget];
			double repulsive = m_d_ljAddRep[mj2+mj1*nAddTarget];

			// Ignore pair if distance is too small, this is usually the self-interaction of the atom. 
			double energyContrib = 0.0;
			double gradContrib = 0.0;

			if (rsq >= 1.0e-6) {
				// Calculate 1/r**2, 1/r**6, 1/r**12. 
				double ir2  = 1.0 / rsq;
				double ir6  = ir2*ir2*ir2;
				double ir12 = ir6*ir6;

				// Calculate the energy. 
				energyContrib = 2.0 * (ir12*repulsive-ir6*attractive);
				// Calculate the gradient. 
				gradContrib = 4.0 * (12.0 * ir12 * repulsive -
						6.0 * ir6 * attractive) * ir2;
			}

			d_energyContrib[tid] = energyContrib;
			d_gradientContrib[3*tid+0] = gradContrib * r.x;
			d_gradientContrib[3*tid+1] = gradContrib * r.y;
			d_gradientContrib[3*tid+2] = gradContrib * r.z;

			tid += blockDim.x * gridDim.x;
		}
	}

}
