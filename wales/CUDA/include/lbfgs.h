/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/)
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File lbfgs.h: Interface of the minimizer.
 *               This is the core class of the library.
 **/

#ifndef LBFGS_H
#define LBFGS_H

#include "cost_function.h"

class Lbfgs
{
	public:
		Lbfgs(CostFunction &costFunc, Printing &debugPrinting, Timer &timer_total, Timer &timer_updates, Timer &timer_linesearch, 
				Cublas &cublas, size_t numDimensions, double H0, int updates, int maxIter, double gradEps, double maxStep, 
				double maxFkRise, int nRedMax);

		~Lbfgs();

		enum lStatus {
			LBFGS_CONVERGED,
			LBFGS_REACHED_MAX_ITER,
			LBFGS_COLD_FUSION_DIAGNOSED, 
			LBFGS_RMS_IS_NAN
		};

		// Returns a string describing the status indicated by the value of lbfgsStatus.
		static std::string statusToString(lStatus lbfgsStatus);

		// Runs minimization of the cost function cf
		// using the L-BFGS method implemented in CUDA.
		//
		// d_x is the device memory location containing
		// the initial guess as cf.getNumDimensions()
		// consecutive doubles. On output, d_x will
		// contain the solution of argmin_x(cf.f(x)) if
		// minimization succeeded, or the last solution
		// found when minimization was aborted.
		//
		// Returns a status code indicating why minimization
		// has stopped, see also Lbfgs::lStatus and
		// Lbfgs::statusToString.
		lStatus minimize(double *d_x, double *d_fk, double *d_gk, double *outRms, int *itDone, bool resetHistory=true);

		void   setMaxIterations(size_t maxIter); 
		void setEvalPercentEps(double evalPercentEps);
		void setProjectGrad(bool projectGrad);
		void setDevzWork(double *d_zWork, const size_t arraySize);

	private:
		CostFunction &m_costFunction;

		Printing &m_debugPrinting;

		Timer &m_timer_total;
		Timer &m_timer_updates;
		Timer &m_timer_linesearch;

		Cublas &m_cublas;

		size_t       m_maxIter;        // Max. no. of steps allowed in minimization.
		size_t       m_cumulativeIter; // The cumulative no. of iterations done over all calls to lbfgs without resetting history.

		const double m_gradientEps;    // Convergence tolerance for RMS force. 
		const double m_maxStep;        // Max. step size allowed (checked in line search). 
		const double m_maxFkRise;      // Max. energy rise allowed (checked in linesearch).

		double       m_evalPercentEps; // Convergence condition for abs(% change) in eigenvalue for R-R minimization. 
		double       m_H0;             // Initial guess for inverse Hessian diagonal elements. 

		const int    m_nRedMax;        // Maximum number of reductions of step size attempted in linesearch.

		int          m_updates;        // History size for minimization.

		bool         m_projectGrad;    // If true,  project out the components of the gradient along the eigenvector. 


		double      *m_d_H0;           // Current estimate for inverse Hessian diagonal elements, on GPU. 
		double      *m_d_zWork;        // The uphill eigenvector, on GPU.

		// Ring buffers for history. 
		double      *m_d_s;            // History of solution updates, on GPU. 
		double      *m_d_y;            // History of gradient updates, on GPU. 
		double      *m_d_alpha;        // History of alphas (needed for z updates), on GPU. 
		double      *m_d_rho;          // History of rhos (needed for z updates), on GPU. 


		bool linesearch(double *d_x, double *d_z, double *d_fk, double *d_gk, size_t it, Lbfgs::lStatus &lbfgsStatus, 
				double *d_step, double *d_tmp, const double *m_d_zWork, double &evPercent, double *outRms, 
				double *d_gkSave);

};
#endif /* end of include guard: LBFGS_H */
