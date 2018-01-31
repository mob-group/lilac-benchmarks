/* This work is a modification of code written by Jens Wetzl and Oliver Taubamann in 2012. 
 * The original work can be found here: https://github.com/jwetzl/CudaLBFGS (license: http://creativecommons.org/licenses/by/3.0/) 
 * This work is not endorsed by the authors. */

/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File cost_function.h: Cost function base classes.
 *                       Derive your own from one of them
 *                       and feed it to the minimizer.
 *
 **/

#ifndef COST_FUNCTION_H
#define COST_FUNCTION_H

#include "printing.h"
#include "timer.h"
#include "cublas.h"

// GPU based cost functions can directly inherit from this class
// and implement computeEnergyAndGradient(), the passed pointers are supposed
// to reside in device memory.
class CostFunction
{
	public:
		CostFunction(Printing &debugPrinting, Timer &timer_potential, Cublas &cublas, size_t numDimensions, int nDegFreedom, 
				int nRigidBody, int rigidMaxSite, int *nRigidSitesPerBody, int *rigidGroups, double *sitesRigidBody, 
				int *rigidSingles, double *rigidInverse, double *coords, int nSecDiag, bool isAtomisticNotRigid, 
				double aaConvThreshold, double coldFusionLim, bool shouldFreeze, bool *isAtomFrozen, int nFreeze, 
				bool isAaConvergence, int nAddTarget, double *ljAddRep, double *ljAddAtt);

		~CostFunction();

		// Calculates secdiag if R-R minimization and calls rigid body transformations and computeEnergyAndGradient. 
		void basePotential(double *d_x, double *d_f, double *d_gradf);

		// Implement this method computing both function value d_f
		// and gradient d_gradf of your cost function at point d_x.
		virtual void computeEnergyAndGradient(const double *d_x, double *d_f, double *d_gradf) = 0;

		void calculateRms(const double *d_gradf, double *outRms);

		size_t getNumDimensions() const;

		int  getnCalls() const;

		bool getIsColdFusion() const;
		bool getIsRayleighRitz() const;

		void setIsRayleighRitz(bool isRayleighRitz);
		void setDevCoords(double *d_x);
		void setDevEnergy(double *d_energy);
		void setDevGradient(double *d_gradient);
		void setIsBfgsts(bool isBfgsts);

		// Converts rigid body coordinates to atomic coordinates before energy calculation. 
		void transformRigidToC(double *d_x);

		// Converts gradient to rigid body gradient. 
		void transformGrad(double *d_gradf, double *d_x);

		// More accurate calculation of RMS force for rigid body framework (RBF). 
		void aaConvergence(const double *d_gradf, double *outRms);

		// Orthogonalise eigenvector to overall rotations and translations. 
		void orthogopt(double *d_x, const bool shouldNormalise);

		// Freeze specified atoms by setting gradient components to zero. 
		void freeze(double *d_gradf);

	protected:
		Cublas &m_cublas;

		const size_t m_numDimensions; // Number of dimensions of function - usually 3*natoms. 
		bool m_isColdFusion; // True if cold fusion has occurred. 
		const double m_coldFusionLim; // Energy limit below which we say cold fusion has occurred. 
		const int m_nAddTarget; // Target cluster size for addressability. 
                const double *m_ljAddRep; // Repulsive epsilon matrix for addressability. 
                const double *m_ljAddAtt; // Attractive epsilon matrix for addressability. 

                int *m_d_nAddTarget; // Target cluster size for addressability, on GPU.
                double *m_d_ljAddRep; // Repulsive epsilon matrix for addressability, on GPU.
                double *m_d_ljAddAtt; // Attractive epsilon matrix for addressability, on GPU.


	private:
		Printing &m_debugPrinting;

		Timer &m_timer_potential;

		const double  m_aaConvThreshold;      // Threshold under which aaConvergence function is called in RMS force calc.

		const double *m_coords;               // The initial coordinates in Bfgsts on CPU, used for R-R minimization.

		const int     m_nSecDiag;             // Method used to approximate eigenvalue in R-R LBFGS.
		const int     m_nFreeze;              // No. of frozen atoms.

		int           m_nCalls;               // No. of potential calls made during this invocation of the CUDA code.

		const bool    m_shouldFreeze;         // If true, freeze some specified coordinates.

		const bool   *m_isAtomFrozen;         // Logical array specifying frozen atoms.

		bool          m_isBfgsts;             // True if costFunction is being used in a Bfgsts calculation - stops aaConvergence call.
		bool          m_isRayleighRitz;       // If true, curvature of PES calculated along direction vec at point coords (R-R min).

		const bool    m_isAtomisticNotRigid;  // If false, use rigid body coordinates.
		const bool    m_isAaConvergence;      // RBF: if true, use more accurate method of calculating RMS force.

		const double *m_sitesRigidBody;       // RBF: coordinates of the rigid body sites.
		const double *m_rigidInverse;         // RBF: inverse eigenvalues of the unweighted tensor of gyration. 

		const int     m_nDegFreedom;          // RBF: no. of degrees of freedom.
		const int     m_nRigidBody;           // RBF: no. of rigid bodies.
		const int     m_rigidMaxSite;         // RBF: max. no. of sites in a rigid body.

		const int    *m_nRigidSitesPerBody;   // RBF: no. of rigid body sites.
		const int    *m_rigidGroups;          // RBF: list of atoms in rigid bodies.
		const int    *m_rigidSingles;         // RBF: list of atoms not in rigid bodies.

		int           m_nLargeRigidBody;      // RBF: no. of rigid bodies with more than 32 sites.

		int          *m_largeRigidIndices;    // RBF: array containing indices for large bodies in nRigidSitesPerBody array. 


		double       *m_d_coords;             // The current coordinates in Bfgsts on GPU, used in R-R minimization. 
		double       *m_d_energy;             // The energy at the current coordinates in Bfgsts on GPU, used in R-R minimization. 
		double       *m_d_gradient;           // The gradient at the current coordinates in Bfgsts on GPU, used in R-R minimization. 
		double       *m_d_zerosx;             // Array used in cuBLAS operations with ones at all the x positions. 
		double       *m_d_zerosy;             // Array used in cuBLAS operations with ones at all the y positions.
		double       *m_d_zerosz;             // Array used in cuBLAS operations with ones at all the z positions.

		bool         *m_d_isAtomFrozen;       // Logical array specifying frozen atoms, on GPU.

		double       *m_d_sitesRigidBody;     // RBF: coordinates of the rigid body sites, on GPU.
		double       *m_d_rigidInverse;       // RBF: IINVERSE from genrigid, used in aaconvergence subroutine, on GPU.
		double       *m_d_xRigid;             // RBF: rigid body coordinates, on GPU.
		double       *m_d_gkRigid;            // RBF: rigid body gradient, on GPU.
		double       *m_d_gkAtoms;            // RBF: copy of atomistic gradient for use with aaConvergence, on GPU.

		int          *m_d_nRigidSitesPerBody; // RBF: no. of rigid body sites, on GPU. 
		int          *m_d_rigidGroups;        // RBF: list of atoms in rigid bodies, on GPU. 
		int          *m_d_rigidSingles;       // RBF: list of atoms not in rigid bodies, on GPU. 
		int          *m_d_nRigidBody;         // RBF: no. of rigid bodies, on GPU. 
		int          *m_d_rigidMaxSite;       // RBF: maximum number of sites in a rigid body, on GPU. 
		int          *m_d_largeRigidIndices;  // RBF: array containing indices for large bodies in nRigidSitesPerBody array, on GPU. 
		int          *m_d_nLargeRigidBody;    // RBF: no. of rigid bodies with more than 32 sites, on GPU. 

		// Rigid body gradient transformations either side of energy/gradient call.
		void rigidTransforms(double *d_x, double *d_f, double *d_gradf);
};

#endif /* end of include guard: COST_FUNCTION_H */
