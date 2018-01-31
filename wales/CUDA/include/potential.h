/**
 *   ___ _   _ ___   _     _       ___ ___ ___ ___
 *  / __| | | |   \ /_\   | |  ___| _ ) __/ __/ __|
 * | (__| |_| | |) / _ \  | |_|___| _ \ _| (_ \__ \
 *  \___|\___/|___/_/ \_\ |____|  |___/_| \___|___/
 *
 * File potential.h: Potential classes deriving from CostFunction. 
 **/

#ifndef POTENTIAL_H
#define POTENTIAL_H

class LjPotential : public CostFunction
{
	public:
		LjPotential(Printing &debugPrinting, Timer &timer_potential, Cublas &cublas, size_t numDimensions, int nDegFreedom, 
				int nRigidBody, int rigidMaxSite, int *nRigidSitesPerBody, int *rigidGroups, double *sitesRigidBody, 
				int *rigidSingles, double *rigidInverse, double *coords, int nSecDiag, bool isAtomisticNotRigid, 
				double aaConvThreshold, double coldFusionLim, bool shouldFreeze, bool *isAtomFrozen, int nFreeze, 
				bool isAaConvergence, int nAddTarget, double *ljAddRep, double *ljAddAtt);

		void computeEnergyAndGradient(const double *d_x, double *d_f, double *d_gradf);
};

class AmberPotential : public CostFunction
{
	public:
		AmberPotential(Printing &debugPrinting, Timer &timer_potential, Cublas &cublas, size_t numDimensions, int nDegFreedom, 
				int nRigidBody, int rigidMaxSite, int *nRigidSitesPerBody, int *rigidGroups, double *sitesRigidBody, 
				int *rigidSingles, double *rigidInverse, double *x, int nSecDiag, bool isAtomisticNotRigid, 
				double aaConvThreshold, double coldFusionLim, bool shouldFreeze, bool *isAtomFrozen, int nFreeze, 
				bool isAaConvergence, int nAddTarget, double *ljAddRep, double *ljAddAtt);

		void computeEnergyAndGradient(const double *d_x, double *d_f, double *d_gradf);
};

#endif /* end of include guard: POTENTIAL_H */
