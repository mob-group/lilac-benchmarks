/**
 *
 * File bfgsts.h: Definition of bfgsts. 
 *
 **/

#ifndef BFGSTS_H
#define BFGSTS_H

#include "lbfgs.h"

class Bfgsts 
{
	public:
		Bfgsts(CostFunction &costFunc, Lbfgs &xLbfgs, Lbfgs &pLbfgs, Printing &debugPrinting, Timer &timer_total, 
				Timer &timer_eFolSteps, Cublas &cublas, size_t numDimensions, double pushOffCut, double pushOffMag, 
				double maxEvecStep, double maxMaxStep, double minMaxStep, double trustRadius, int pMaxIter1, int pMaxIter2, 
				int bItMax, double evecOverlapTol, double evStepMagTol);

		enum bStatus {
			BFGSTS_CONVERGED_TO_TS,
			BFGSTS_REACHED_MAX_ITER,
			BFGSTS_COLD_FUSION_DIAGNOSED
		};

		// Returns a string describing the status indicated by the value of bfgsStatus.
		static std::string statusToString(bStatus bfgsStatus);

		// Main Bfgsts calculation. 
		bStatus findTs(double *d_coords, double *d_evec, double *energy, double *evalMin, double *outRms, int *itDone);

	private:
		// Eigenvector following step uphill.  
		void stepUphill(double *d_coords, const double *d_evec, const double *evalMin, double *d_energy, 
				double *d_gradient, const bool isxLbfgsConverged, const bool &isNegativeEv, double &stpMag, 
				double *outRms, size_t it);

		CostFunction &m_costFunction;

		Lbfgs &m_xLbfgs; // Rayleigh-Ritz minimizer. 
		Lbfgs &m_pLbfgs; // Tangent space minimizer. 

		Printing &m_debugPrinting;

		Timer &m_timer_total;
		Timer &m_timer_eFolSteps;

		Cublas &m_cublas;


		const double m_pushOffCut;     // RMS threshold under which a pushoff will be applied. 
		const double m_pushOffMag;     // Magnitude of step away from transition state. 
		const double m_maxMaxStep;     // Max. allowed uphill step size. 
		const double m_minMaxStep;     // Min. value that max. uphill step size is allowed to fall to. 
		const double m_trustRadius;    // Trust radius for adjusting max. step along eigenvector. 
		const double m_evecOverlapTol; // Overlap between current and previous eigenvector estimate. 
		const double m_evStepMagTol;   // Tolerance for eigenvector overlap (no. tangent space steps small->large). 

		double       m_maxEvecStep;    // Max. step along the eigenvector.

		const int    m_pMaxIter1;      //  No. of LBFGS steps allowed in restricted minimization (grad projected). 
		const int    m_pMaxIter2;      // No. of LBFGS steps allowed in tangent space after eigenvalue converged. 
		const int    m_bItMax;         // Max. no. of bfgsts iterations allowed. 

};
#endif /* end of include guard: BFGSTS_H */
