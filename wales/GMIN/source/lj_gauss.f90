!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2017 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
! This module provides a routine for calculating the energy and gradients with
! a double well Lennard-Jones + Gaussian potential. The possibility of varying
! and optimising the potential parameneter is included.
! Also included is a framework for a varying triclinic unit cell. This could be
! separated out and generalised to all isotropic pair potentials.
! Questions to jwrm2.
!
MODULE LJ_GAUSS_MOD

! Public parameters
PUBLIC :: LJ_GAUSS_MODE, LJ_GAUSS_RCUT, LJ_GAUSS_EPS, LJ_GAUSS_R0
PUBLIC :: LJ_GAUSS_SIGMASQ, LJ_GAUSS_PARAMS
! Public functions
PUBLIC :: LJ_GAUSS, LJ_GAUSS_INITIALISE, VIEW_LJ_GAUSS, LJ_GAUSS_TAKESTEP
! Private parameters
PRIVATE :: RCUTINV6, RCUTINV12, EXP_RCUT, PREF_RCUT, PI, SCL_FCT
PRIVATE :: IMAGE_CUTOFF, MAX_LENGTH_STEP, MAX_ANGLE_STEP, MAX_EPS_STEP
PRIVATE :: MAX_R0_STEP, MAX_SIGMASQ_STEP, V_EPS, V_SIGMA, EPS_LOWER, EPS_UPPER
PRIVATE :: EPS_REPULSION, R0_LOWER, R0_UPPER, R0_REPULSION, SIGMASQ_LOWER
PRIVATE :: SIGMASQ_UPPER, SIGMASQ_REPULSION
! Private functions
PRIVATE :: PERIODIC_RESET, LJ_GAUSS_PER, CALC_ENERGY, CALC_ENERGY_LJ
PRIVATE :: CALC_ENERGY_GAUSS, CALC_GRAD, CALC_FCTS, CONSTRAIN_PARAMETERS
PRIVATE :: CHECK_ANGLES, REJECT, CONSTRAIN_VOLUME

! EPS: strength of the Gaussian potential. The LJ strength is 1.
! R0: Equilibrium distance of the minimum in the Gaussian potential. The LJ
!     rmin is 1.
! SIGMASQ: variation (width) of the Gaussian potential
! RCUTSQ: cutoff distance at which the energy, gradients and Hessian of
!         both parts of the potential go smoothly to zero.
DOUBLE PRECISION :: LJ_GAUSS_EPS, LJ_GAUSS_R0, LJ_GAUSS_SIGMASQ, LJ_GAUSS_RCUT

! LJ_GAUSS_MODE: 0 standard minimisation
!                1 the final triplet of coordinates will be interpreted as the
!                  unit cell parameters of an orthogonal unit cell, and
!                  optimised
!                2 the final triplet of coordinates are the unit cell length
!                  parameters and the second last triplet are the unit cell
!                  angles, giving a triclinic unit cell, which will be optimised
!                3 as for 2, but the third last triplet are the potential
!                  parameters LJ_GAUSS_EPS, LJ_GAUSS_R0, LJ_GAUSS_SIGMASQ and
!                  the potential parameters will be optimised, depending on the
!                  value of LJ_GAUSS_PARAMS
! LJ_GAUSS_PARAMS: Specifies which of the three potential parameters will be
!                  optimised. An integer between 1 and 7.
!                  1 bit set: LJ_GAUSS_EPS will be optimised
!                  2 bit set: LJ_GAUSS_R0 will be optimised
!                  3 bit set: LJ_GAUSS_SIGMASQ will be optimised
INTEGER :: LJ_GAUSS_MODE, LJ_GAUSS_PARAMS

! RCUTINV6, RCUTINV12: factors for LJ which don't depend on particle distance
DOUBLE PRECISION :: RCUTINV6, RCUTINV12

! EXP_RCUT, PREF_RCUT: factors for Gauss which don't depend on the particle
!                      distance
! SCL_FCT: overall scale factor for the potential, to allow comparison between
!          different parameters
DOUBLE PRECISION :: EXP_RCUT, PREF_RCUT, SCL_FCT

! PI: the usual meaning
DOUBLE PRECISION, PARAMETER :: PI = 4*ATAN(1.D0)

! IMAGE_CUTOFF: if more periodic images than this in any direction fall
!               within the cutoff, the evaluation will not be carried
!               out and the step will be rejected.
DOUBLE PRECISION, PARAMETER :: IMAGE_CUTOFF = 30

! Parameters specifying the largest steps in LJ_GAUSS_TAKESTEP
! MAX_LENGTH_STEP: largest allowed step in the length of a lattice
!                  vector
! MAX_ANGLE_STEP: largest allowed step in the angle of a lattice vector
! MAX_EPS_STEP: largest allowed step in LJ_GAUSS_EPS
! MAX_R0_STEP: largest allowed step in LJ_GAUSS_R0
! MAX_SIGMASQ_STEP: largest allowed step in LJ_GAUSS_SIGMASQ
DOUBLE PRECISION, PARAMETER :: MAX_LENGTH_STEP = 0.3D0
DOUBLE PRECISION, PARAMETER :: MAX_ANGLE_STEP = 0.1D0
DOUBLE PRECISION, PARAMETER :: MAX_EPS_STEP = 0.2D0
DOUBLE PRECISION, PARAMETER :: MAX_R0_STEP = 0.1D0
DOUBLE PRECISION, PARAMETER :: MAX_SIGMASQ_STEP = 0.001D0

! Parameters used for constraining the unit cell volume and the potential
! parameters
! V_EPS: scaling factor for the unit cell volume contraint energy
! V_SIGMA: WCA cutoff distance
! x_LOWER, x_UPPER: bounds on potential parameters
! x_REPULSION: harmonic force constant for each potential parameter
DOUBLE PRECISION, PARAMETER :: V_EPS = 1.D-3
DOUBLE PRECISION, PARAMETER :: V_SIGMA = 3.D-1
DOUBLE PRECISION, PARAMETER :: EPS_LOWER = 0.D0, EPS_UPPER = 5.D0
DOUBLE PRECISION, PARAMETER :: EPS_REPULSION = 1.D4
DOUBLE PRECISION, PARAMETER :: R0_LOWER = 1.D0, R0_UPPER = 2.1D0
DOUBLE PRECISION, PARAMETER :: R0_REPULSION = 1.D4
DOUBLE PRECISION, PARAMETER :: SIGMASQ_LOWER = 0.01D0, SIGMASQ_UPPER = 0.03D0
DOUBLE PRECISION, PARAMETER :: SIGMASQ_REPULSION = 1.D9

CONTAINS

!-------------------------------------------------------------------------------
!
! Main routine for calculating the energy and gradients.
! X: coordinates array
! GRAD: gradients array
! ENERGY: calculated energy
! GTEST: whether gradients should be calculated
!
!-------------------------------------------------------------------------------
    SUBROUTINE LJ_GAUSS(X, GRAD, ENERGY, GTEST)

        ! MYUNIT: file unit for main output file
        ! NATOMS: number of particles
        ! PERIODIC: whether to use periodic boundary conditions
        USE COMMONS, ONLY: MYUNIT, NATOMS, PERIODIC

        ! BOX_DERIV: type for storing box derivative information
        ! CALC_BOX_DERIV: calculates box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV, CALC_BOX_DERIV

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS)
        DOUBLE PRECISION, INTENT(OUT) :: GRAD(3*NATOMS), ENERGY
        LOGICAL, INTENT(IN) :: GTEST

        ! I, J: indices for looping over particles
        ! NMOL: actual number of particles
        INTEGER :: I, J, NMOL

        ! RIJ: Interparticle vector
        ! MODRIJ: distance between a pair of particles
        ! DVDR: derivative of pair potential wrt MODRIJ
        DOUBLE PRECISION :: RIJ(3), MODRIJ, DVDR

        ! Factors for triclinic lattice
        TYPE(BOX_DERIV) :: BD

!        DO I = 1, NATOMS
!            WRITE (MYUNIT, *) X(3*I-2:3*I)
!        END DO

        ! Initialise output variables
        ENERGY = 0.D0
        GRAD(:) = 0.D0

        ! Calculate factors that do depend on parameters, but not on RIJ
        CALL CALC_FCTS(X(NATOMS*3-8:NATOMS*3-6))

        IF (.NOT. PERIODIC) THEN
            IF (LJ_GAUSS_MODE .EQ. 0) THEN
                ! Normal cluster
                NMOL = NATOMS
                DO I = 1, NMOL-1 ! Outer loop over particles
                    DO J = I+1, NMOL ! Inner loop over particles
                        ! Get the particle-particle distance
                        RIJ(:) = X(3*I-2:3*I) - X(3*J-2:3*J)
                        MODRIJ = SQRT(DOT_PRODUCT(RIJ(:), RIJ(:)))

                        ! Add energy and gradients
                        IF (MODRIJ .LT. LJ_GAUSS_RCUT) THEN
                            ENERGY = ENERGY + CALC_ENERGY(MODRIJ)

                            IF (GTEST) THEN
                                DVDR = CALC_GRAD(MODRIJ)
                                GRAD(3*I-2:3*I) = GRAD(3*I-2:3*I) + &
                                                  DVDR*RIJ(:)/MODRIJ
                                GRAD(3*J-2:3*J) = GRAD(3*J-2:3*J) - &
                                                  DVDR*RIJ(:)/MODRIJ
                            END IF ! GTEST
                        END IF ! If less than cutoff distance
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles
            ELSE ! mode not equal to zero
                WRITE (MYUNIT, *) 'LJ_GAUSS> ERROR, PERIODIC must be ', &
                                  'specified if mode is not 0'
                STOP
            END IF
        ELSE ! PERIODIC
            SELECT CASE (LJ_GAUSS_MODE)
            CASE(0)
                ! Periodic, fixed box
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

                NMOL = NATOMS
                ! We need to include self-interactions of particles in different
                ! unit cells
                DO I = 1, NMOL ! Outer loop over particles
                    DO J = I, NMOL ! Inner loop over particles
                        ! Add energy and gradients
                        CALL LJ_GAUSS_PER(I, J, X(3*I-2:3*I), X(3*J-2:3*J), &
                                          GRAD(3*I-2:3*I), GRAD(3*J-2:3*J), &
                                          ENERGY, GTEST)
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

            CASE(1)
                ! Periodic varying lattice
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

                NMOL = NATOMS-1
                ! We need to include self-interactions of particles in different
                ! unit cells
                DO I = 1, NMOL! Outer loop over particles
                    DO J = I, NMOL ! Inner loop over particles
                        ! Add energy and gradients
                        CALL LJ_GAUSS_PER(I, J, X(3*I-2:3*I), X(3*J-2:3*J), &
                                          GRAD(3*I-2:3*I), GRAD(3*J-2:3*J), &
                                          ENERGY, GTEST, &
                                          GRAD(3*NATOMS-2:3*NATOMS))
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

            CASE(2)
                ! Triclinic varying lattice
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

                ! Calculate box derivative factors
                BD = CALC_BOX_DERIV(X(3*NATOMS-5:3*NATOMS-3), &
                                    X(3*NATOMS-2:3*NATOMS  ), LJ_GAUSS_RCUT)

!                WRITE(MYUNIT, *) 'Box lengths: ', X(3*NATOMS-2:3*NATOMS)
!                WRITE(MYUNIT, *) 'Box angles: ', X(3*NATOMS-5:3*NATOMS-3)
!                WRITE(MYUNIT, *) 'CELL_RANGE(:) ', BD%CELL_RANGE(:)
!                WRITE(MYUNIT, *) 'V = ', BD%V

                ! Check whether the unit cell angles are physically possible.
                ! Reject if not.
                IF (.NOT. CHECK_ANGLES(X(3*NATOMS-5:3*NATOMS-3))) THEN
                    CALL REJECT(ENERGY, GRAD)
!                    WRITE(MYUNIT, *) 'Reject for CHECK_ANGLES'
                    RETURN
                END IF

                ! Reject steps where the cell range has got bigger than the 
                ! cutoff. Such unit cells are highly distorted and there are 
                ! probably equivalent ones of a more usual shape. Allowing them
                ! leads to very slow run times.
                IF (.NOT. ALL(BD%CELL_RANGE .LE. IMAGE_CUTOFF)) THEN
                    CALL REJECT(ENERGY, GRAD)
!                    WRITE(MYUNIT, *) 'Reject for CELL_RANGE'
                    RETURN
                END IF

                NMOL = NATOMS-2
                ! We need to include self-interactions of particles in different
                ! unit cells
                DO I = 1, NMOL! Outer loop over particles
                    DO J = I, NMOL ! Inner loop over particles
                        ! Add energy and gradients
                        CALL LJ_GAUSS_PER_TRI(I, J, X(:), GRAD(3*I-2:3*I), &
                            GRAD(3*J-2:3*J), ENERGY, GTEST, &
                            GRAD(3*NATOMS-5:3*NATOMS), BD)
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

                ! Constrain the box volume to be greater than 0 with a WCA-style
                ! repulsion
                CALL CONSTRAIN_VOLUME(ENERGY, GRAD(3*NATOMS-5:3*NATOMS-3), BD, &
                                      GTEST)

            CASE(3)
                ! Triclinic varying lattice, with varying parameters
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

                ! Calculate box derivative factors
                BD = CALC_BOX_DERIV(X(3*NATOMS-5:3*NATOMS-3), &
                                    X(3*NATOMS-2:3*NATOMS  ), LJ_GAUSS_RCUT)

                ! Check whether the unit cell angles are physically possible.
                ! If we have a disallowed unit cell, we set the energy to very
                ! high and the gradients to very small, so the step is
                ! 'converged' but gets rejected.
                IF (.NOT. CHECK_ANGLES(X(3*NATOMS-5:3*NATOMS-3))) THEN
                    CALL REJECT(ENERGY, GRAD)
                END IF

                ! Reject steps where the cell range has got bigger than the 
                ! cutoff. Such unit cells are highly distorted and there are 
                ! probably equivalent ones of a more usual shape. Allowing them
                ! leads to very slow run times.
                IF (.NOT. ALL(BD%CELL_RANGE .LE. IMAGE_CUTOFF)) THEN
                    CALL REJECT(ENERGY, GRAD)
                    RETURN
                END IF

                NMOL = NATOMS-3
                ! We need to include self-interactions of particles in different
                ! unit cells
                DO I = 1, NMOL! Outer loop over particles
                    DO J = I, NMOL ! Inner loop over particles
                        ! Add energy and gradients
                        CALL LJ_GAUSS_PER_TRI(I, J, X(:), GRAD(3*I-2:3*I), &
                            GRAD(3*J-2:3*J), ENERGY, GTEST, &
                            GRAD(3*NATOMS-5:3*NATOMS), BD, &
                            GRAD(3*NATOMS-8:3*NATOMS-6))
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

                ! Constrain the box volume to be greater than 0 with a WCA-style
                ! repulsion
                CALL CONSTRAIN_VOLUME(ENERGY, GRAD(3*NATOMS-5:3*NATOMS-3), BD, &
                                      GTEST)

                ! Constrain the potential parameters with a harmonic repulsion
                CALL CONSTRAIN_PARAMETERS(ENERGY, GRAD(3*NATOMS-8:3*NATOMS-6), &
                                          GTEST)
            END SELECT ! Mode selections
        END IF ! Periodic

        ! Apply the overall scaling
        ENERGY = ENERGY * SCL_FCT
        GRAD(:) = GRAD(:) * SCL_FCT

        ! The gradients wrt the potential parameters already include SCL_FCT
        IF (LJ_GAUSS_MODE .EQ. 3) THEN
            GRAD(3*NATOMS-8:3*NATOMS-6) = GRAD(3*NATOMS-8:3*NATOMS-6) / SCL_FCT
        END IF

    END SUBROUTINE LJ_GAUSS

!-------------------------------------------------------------------------------
!
! Rejects a step by setting the energy very high and the gradients very low.
! This tricks L-BFGS into thinking it's converged, but the high energy of the
! quench will ensure it is rejected.
! ENERGY: system energy
! GRAD: system gradients
!
!-------------------------------------------------------------------------------

    SUBROUTINE REJECT(ENERGY, GRAD)

        ! NATOMS: number of atoms
        USE COMMONS, ONLY: NATOMS

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(OUT) :: ENERGY, GRAD(3*NATOMS)

        ENERGY = 1.D20
        GRAD(:) = 1.D-20

    END SUBROUTINE REJECT

!-------------------------------------------------------------------------------
!
! Calculates the factors that depend on the potential parameters but not on the
! distance between particles.
! PARAMS: the potential parameters when those are to be adjusted, otherwise
!         ignored
!
!-------------------------------------------------------------------------------

    SUBROUTINE CALC_FCTS(PARAMS)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: PARAMS(3)

        ! Adjust the potential parameters, if the mode requires it
        IF (LJ_GAUSS_MODE .EQ. 3) THEN
            ! For each parameter, if it is varying, set it from the coordinate.
            ! If it's not varying, reset the coordinate, which TAKESTEP will
            ! have changed.
            IF (BTEST(LJ_GAUSS_PARAMS, 0)) THEN
                LJ_GAUSS_EPS = PARAMS(1)
            ELSE
                PARAMS(1) = LJ_GAUSS_EPS
            END IF
            IF (BTEST(LJ_GAUSS_PARAMS, 1)) THEN
                LJ_GAUSS_R0 = PARAMS(2)
            ELSE
                PARAMS(2) = LJ_GAUSS_R0
            END IF
            IF (BTEST(LJ_GAUSS_PARAMS, 2)) THEN
                LJ_GAUSS_SIGMASQ = PARAMS(3)
            ELSE
                PARAMS(3) = LJ_GAUSS_SIGMASQ
            END IF
        END IF

        ! Overall scale factor
        SCL_FCT = 1.D0 / (SQRT(2.D0 * PI * LJ_GAUSS_SIGMASQ) * LJ_GAUSS_EPS + &
                          2.D0**(5.D0/6.D0) * 12.D0/55.D0)

        ! Factors for the Gaussian
        EXP_RCUT = LJ_GAUSS_EPS * EXP(- (LJ_GAUSS_RCUT - LJ_GAUSS_R0)**2 / &
                                        (2.D0*LJ_GAUSS_SIGMASQ))
        PREF_RCUT = (LJ_GAUSS_RCUT - LJ_GAUSS_R0) / LJ_GAUSS_SIGMASQ

    END SUBROUTINE CALC_FCTS

!-------------------------------------------------------------------------------
!
! Calculates the energy contribution for a pair of particles.
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_ENERGY(MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_ENERGY = CALC_ENERGY_LJ(MODRIJ) + CALC_ENERGY_GAUSS(MODRIJ)

    END FUNCTION CALC_ENERGY

!-------------------------------------------------------------------------------
!
! Calculates the LJ contribution to the energy
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------

    PURE DOUBLE PRECISION FUNCTION CALC_ENERGY_LJ (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        ! LJ6, LJ12: factors for the LJ potential
        ! ENERGY: temporary energy store
        DOUBLE PRECISION :: LJ6, LJ12, ENERGY

        ! Evaluate LJ energy
        LJ6 = (1.D0/MODRIJ)**6
        LJ12 = LJ6**2
        ENERGY = LJ12 - 2.D0*LJ6
        ! Correction to LJ to make potential go to zero at RCUT
        ENERGY = ENERGY - (RCUTINV12 - 2*RCUTINV6)
        ! Correction to LJ to make derivative go to zero at RCUT
        ENERGY = ENERGY + (12.D0/LJ_GAUSS_RCUT) * &
                          (RCUTINV12 - RCUTINV6)* &
                          (MODRIJ - LJ_GAUSS_RCUT)
        ! Correction to LJ to make second derivatives go to zero at
        ! RCUT
!        ENERGY = ENERGY + (1.D0/LJ_GAUSS_RCUT**2) * &
!                          (42.D0*RCUTINV6 - 78.D0*RCUTINV12) * &
!                          (MODRIJ - LJ_GAUSS_RCUT)**2
        
        ! Set output
        CALC_ENERGY_LJ = ENERGY

    END FUNCTION CALC_ENERGY_LJ

!-------------------------------------------------------------------------------
!
! Calculates the Gaussian contribution to the energy
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------

    PURE DOUBLE PRECISION FUNCTION CALC_ENERGY_GAUSS (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        ! EXP_RIJ: factor for the Gauss potential
        ! ENERGY: temporary energy store
        DOUBLE PRECISION :: EXP_RIJ, ENERGY
        
        ! Set output
        CALC_ENERGY_GAUSS = ENERGY

        ! Evaluate Gaussian energy
        EXP_RIJ = LJ_GAUSS_EPS * &
                  EXP(- (MODRIJ - LJ_GAUSS_R0)**2 / &
                        (2.D0*LJ_GAUSS_SIGMASQ))
        ENERGY = - EXP_RIJ
        ! Correction to Gauss to make potential go to zero at RCUT
        ENERGY = ENERGY + EXP_RCUT
        ! Correction to Gauss to make derivatives go to zero at RCUT
        ENERGY = ENERGY - PREF_RCUT*EXP_RCUT*(MODRIJ-LJ_GAUSS_RCUT)
        ! Correction to Gauss to make second derivatives go to zero at
        ! RCUT
!        ENERGY = ENERGY + EXP_RCUT * (PREF_RCUT**2 - 1.D0/LJ_GAUSS_SIGMASQ) * &
!                          0.5D0*((MODRIJ - LJ_GAUSS_RCUT)**2)

        ! Set output
        CALC_ENERGY_GAUSS = ENERGY

    END FUNCTION CALC_ENERGY_GAUSS

!-------------------------------------------------------------------------------
!
! Calculates the pair potential gradient for a pair of particles.
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_GRAD (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        ! LJ6, LJ12: factors for the LJ potential
        DOUBLE PRECISION :: LJ6, LJ12
        ! EXP_RIJ, PREF_RIJ: factor for the Gauss potential
        DOUBLE PRECISION :: EXP_RIJ, PREF_RIJ
        ! DVDR: temporary store for the result
        DOUBLE PRECISION :: DVDR

        ! Initialise output variables
        DVDR = 0.D0

        ! LJ factors
        LJ6 = (1.D0/MODRIJ)**6
        LJ12 = LJ6**2

        ! Evaluate LJ derivatives
        DVDR = DVDR + (12.D0/MODRIJ) * (LJ6 - LJ12)
        ! Correction to LJ to make derivaitves go to zero at
        ! RCUT
        DVDR = DVDR + (12.D0/LJ_GAUSS_RCUT) * (RCUTINV12 - RCUTINV6)
        ! Correction to LJ to make second derivatives go to zero
        ! at RCUT
!        DVDR = DVDR + (12.D0/LJ_GAUSS_RCUT**2) * &
!                      (7.D0*RCUTINV6 - 13.D0*RCUTINV12) * &
!                      (MODRIJ - LJ_GAUSS_RCUT)

        ! Gauss factors
        PREF_RIJ = (MODRIJ - LJ_GAUSS_R0) / LJ_GAUSS_SIGMASQ
        EXP_RIJ = LJ_GAUSS_EPS * EXP(- (MODRIJ - LJ_GAUSS_R0)**2 / &
                                       (2.D0*LJ_GAUSS_SIGMASQ))

        ! Evaluate Gauss derivatives
        DVDR = DVDR + PREF_RIJ * EXP_RIJ
        ! Correction to Gauss to make derivatives go to zero at
        ! RCUT
        DVDR = DVDR - PREF_RCUT * EXP_RCUT
        ! Correction to Gauss to make second derivatives go to
        ! zero at RCUT
!        DVDR = DVDR - (EXP_RCUT / LJ_GAUSS_SIGMASQ - &
!                       (PREF_RCUT**2) * EXP_RCUT) * &
!                      (MODRIJ - LJ_GAUSS_RCUT)

        CALC_GRAD = DVDR

    END FUNCTION CALC_GRAD

!-------------------------------------------------------------------------------
!
! Initialisation. Calculate some factors that won't be changing
!
!-------------------------------------------------------------------------------
    SUBROUTINE LJ_GAUSS_INITIALISE ()

        ! MYUNIT: file unit for main GMIN output
        USE COMMONS, ONLY: MYUNIT

        IMPLICIT NONE

        ! Check we have sensible value for the mode
        IF (LJ_GAUSS_MODE .NE. 1 .AND. LJ_GAUSS_MODE .NE. 0 .AND. &
            LJ_GAUSS_MODE .NE. 2 .AND. LJ_GAUSS_MODE .NE. 3) THEN
            WRITE (MYUNIT, *) 'LJ_GAUSS> ERROR: mode must be 0, 1, 2 or 3'
            STOP
        END IF

        ! Check we have sensible values for the params
        IF (LJ_GAUSS_MODE .EQ. 3) THEN
            IF (LJ_GAUSS_PARAMS .LT. 1 .OR. LJ_GAUSS_PARAMS .GT. 7) THEN
                WRITE (MYUNIT, *) 'LJ_GAUSS> ERROR: params must be between 1', &
                                  ' and 7'
                STOP
            END IF
        END IF          

        ! Calculate some factors for LJ that don't depend on RIJ or the
        ! potential parameters
        RCUTINV6 = (1.D0/LJ_GAUSS_RCUT)**6
        RCUTINV12 = RCUTINV6**2

    END SUBROUTINE LJ_GAUSS_INITIALISE

!-------------------------------------------------------------------------------
!
! Takes a step, perturbing the Cartesian coordinates, the unit cell parameters
! and the potential parameters, according to the specified mode.
! NP: the index in the main COORDS array to take the step from
!
!-------------------------------------------------------------------------------

    SUBROUTINE LJ_GAUSS_TAKESTEP(NP)

        ! COORDS: full set of Markov chain coordinates
        ! MYUNIT: file unit of the main GMIN output file
        ! NATOMS: number of coordinates
        ! PERIODIC: whetehr periodic boundary conditions are used
        ! PERCOLATET: whether percolation is used to constrain particles
        ! RADIUS: container radius
        ! STEP: maximum step size for this move
        ! TMOVE: specifies which step to do a translational move on
        ! TWOD: whetehr the system is two dimensional
        USE COMMONS, ONLY: COORDS, MYUNIT, NATOMS, PERIODIC, PERCOLATET, & 
                           RADIUS, STEP, TMOVE, TWOD
        
        ! VEC_LEN: function, returns length of a 3D vector
        ! VEC_RANDOM: function, returns a randomly orientated 3D unit vector
        USE VEC3, ONLY: VEC_LEN, VEC_RANDOM

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: NP

        ! X: local copy of the coordinates
        ! RAND_STEP: unit vector in a random step direction
        ! STEP_SIZE: random fraction of the max step size to move
        ! NE_ANGLES: new set of unit cell angles, to be checked
        DOUBLE PRECISION :: X(3*NATOMS), RAND_STEP(3), STEP_SIZE, NEW_ANGLES(3)

        ! I: iteration index
        ! NMOL: actual number of atoms
        INTEGER :: I, NMOL

        ! Set the local copy of the coordinates
        X(:) = COORDS(:, NP)

        ! Set NMOL
        NMOL = NATOMS - LJ_GAUSS_MODE

        ! Random translational steps for each of the particles
        IF (TMOVE(NP)) THEN
            DO I = 1, NMOL
                ! Generate either a 2D or 3D random vector
                ! The magnitude will be uniformly distributed in a circle or 
                ! sphere with radius given by the value of STEP for this step
                IF (TWOD) THEN
                    DO
                        CALL RANDOM_NUMBER(RAND_STEP(1))
                        CALL RANDOM_NUMBER(RAND_STEP(2))
                        RAND_STEP(3) = 0.D0
                        IF(VEC_LEN(RAND_STEP) .LE. 1.D0) EXIT
                    END DO
                    RAND_STEP(:) = RAND_STEP(:) * STEP(NP)
                ELSE ! 3D
                    RAND_STEP(:) = VEC_RANDOM()
                    CALL RANDOM_NUMBER(STEP_SIZE)
                    ! Sqrt for random volume distribution
                    RAND_STEP(:) = RAND_STEP(:) * SQRT(STEP_SIZE * STEP(NP))
                END IF
                ! Add the geenrated ste onto the coordinates
                X(I*3-2:I*3) = X(I*3-2:I*3) + RAND_STEP(:)
            END DO
        END IF

        ! If not periodic and we're using a radius container, bring particles in
        IF (.NOT. PERIODIC .AND. .NOT. PERCOLATET) THEN
            DO I = 1, NMOL
                IF (VEC_LEN(X(I*3-2:I*3)) .GT. RADIUS) THEN
                    WRITE(MYUNIT, *) 'LJ_GAUSS_TAKESTEP> ', &
                        'coord outside container, bringing in'
                    X(I*3-2:I*3) = X(I*3-2:I*3) - &
                        SQRT(RADIUS) * NINT(X(I*3-2:I*3) / SQRT(RADIUS))
                END IF
            END DO
        END IF

        ! Box length steps
        IF (LJ_GAUSS_MODE .GE. 1) THEN
            X(3*NATOMS-2:3*NATOMS) = X(3*NATOMS-2:3*NATOMS) + &
                                     VEC_RANDOM() * MAX_LENGTH_STEP
        END IF

        ! Box angle steps
        IF (LJ_GAUSS_MODE .GE. 2) THEN
            DO
                NEW_ANGLES(:) = X(3*NATOMS-5:3*NATOMS-3) + &
                                VEC_RANDOM() * MAX_ANGLE_STEP

            ! See whether the unit cell we've generated is valid
                IF (CHECK_ANGLES(NEW_ANGLES(:))) THEN
                    X(3*NATOMS-5:3*NATOMS-3) =  NEW_ANGLES(:)
                    EXIT
                END IF
            END DO
        END IF

        ! Parameter steps
        IF (LJ_GAUSS_MODE .EQ. 3) THEN
            IF (BTEST(LJ_GAUSS_PARAMS, 0)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-8) = X(3*NATOMS-8) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_EPS_STEP
            END IF
            IF (BTEST(LJ_GAUSS_PARAMS, 1)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-7) = X(3*NATOMS-7) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_R0_STEP
            END IF
            IF (BTEST(LJ_GAUSS_PARAMS, 2)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-6) = X(3*NATOMS-6) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_SIGMASQ_STEP
            END IF            
        END IF

        COORDS(:, NP) = X(:)

    END SUBROUTINE LJ_GAUSS_TAKESTEP

!-------------------------------------------------------------------------------
!
! Resets all coordinates to within the unit cell
! The unit cell runs from 0 to 1 in all directions. Coordinates are later scaled
! by the box lengths.
! COORDS: coordinates of the particles
!
!-------------------------------------------------------------------------------
    SUBROUTINE PERIODIC_RESET(COORDS)

        ! BOXLX, BOXLY, BOXLZ: dimensions of box
        ! NATOMS: number of particles
        USE COMMONS, ONLY: BOXLX, BOXLY, BOXLZ, NATOMS

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATOMS)

        ! Sort out the unit cell sizes, if they are varying
        IF (LJ_GAUSS_MODE .EQ. 1 .OR. LJ_GAUSS_MODE .EQ. 2) THEN
            BOXLX = COORDS(3*NATOMS - 2)
            BOXLY = COORDS(3*NATOMS - 1)
            BOXLZ = COORDS(3*NATOMS)
        END IF

        ! Reset coordinates to within the box
        SELECT CASE (LJ_GAUSS_MODE)
        CASE(0)
            COORDS(:) = COORDS(:) - FLOOR(COORDS(:))
        CASE(1)
            COORDS(1:3*NATOMS-3) = COORDS(1:3*NATOMS-3) - &
                                   FLOOR(COORDS(1:3*NATOMS-3))
        CASE(2)
            COORDS(1:3*NATOMS-6) = COORDS(1:3*NATOMS-6) - &
                                   FLOOR(COORDS(1:3*NATOMS-6))
            ! Make the lattice angles modulo 2*pi
            ! They will be checked later
            COORDS(3*NATOMS-5:3*NATOMS-3) = COORDS(3*NATOMS-5:3*NATOMS-3) - &
                2.D0*PI*FLOOR(COORDS(3*NATOMS-5:3*NATOMS-3)/(2.D0*PI))
        CASE(3)
           COORDS(1:3*NATOMS-9) = COORDS(1:3*NATOMS-9) - &
                                   FLOOR(COORDS(1:3*NATOMS-9))
            ! Make the lattice angles modulo 2*pi
            ! They will be checked later
            COORDS(3*NATOMS-5:3*NATOMS-3) = COORDS(3*NATOMS-5:3*NATOMS-3) - &
                2.D0*PI*FLOOR(COORDS(3*NATOMS-5:3*NATOMS-3)/(2.D0*PI))

        END SELECT

    END SUBROUTINE PERIODIC_RESET

!-------------------------------------------------------------------------------
!
! Finds the energy and gradients for the periodic system
! IDI, IDJ: id numbers of the particles
! COORDS_I: coordinates of particle I
! COORDS_J: coordinates of particle J
! GRAD_I: gradients for particle I
! GRAD_J: gradients for particle J
! ENERGY: energy contribution for this pair
! RIJ: mininum length vector between the particles
! GTEST: whether to calculate gradients
! GRAD_BOX: gradients for the box
!
!-------------------------------------------------------------------------------

    SUBROUTINE LJ_GAUSS_PER(IDI, IDJ, COORDS_I,COORDS_J,GRAD_I,GRAD_J,ENERGY, &
                            GTEST, GRAD_BOX)

        ! BOXLX, BOXLY, BOXLZ: dimensions of the box
        USE COMMONS, ONLY: BOXLX, BOXLY, BOXLZ

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: IDI, IDJ
        DOUBLE PRECISION, INTENT(IN) :: COORDS_I(3), COORDS_J(3)
        DOUBLE PRECISION, INTENT(INOUT) :: GRAD_I(3), GRAD_J(3), ENERGY
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: GRAD_BOX(3)
        LOGICAL, INTENT(IN) :: GTEST

        ! XJ: copy of j coordinates
        ! RIJ: particle-particle vector
        ! MODRIJ: particle-particle distance
        ! DVDR: derivative of pair potential wrt MODRIJ
        ! BOX_SIZE: convenience copy of the box size
        DOUBLE PRECISION :: XJ(3), RIJ(3), MODRIJ, DVDR, BOX_SIZE(3)

        ! PER_k: integers for looping over cell possibilities
        ! CELL_RANGE: number of cells in each direction required
        INTEGER :: PER_X, PER_Y, PER_Z, CELL_RANGE(3)

!        DOUBLE PRECISION :: START_ENERGY
!        START_ENERGY = ENERGY

        BOX_SIZE(1) = BOXLX
        BOX_SIZE(2) = BOXLY
        BOX_SIZE(3) = BOXLZ

        ! Set the number of cells we need to check, based on the cell size and
        ! the potential cutoff.
        CELL_RANGE(:) = CEILING(LJ_GAUSS_RCUT/BOX_SIZE(:))

        ! Loop over different coordinate possibilities.
        DO PER_X = - CELL_RANGE(1), CELL_RANGE(1)
            DO PER_Y = - CELL_RANGE(2), CELL_RANGE(2)
                DO PER_Z = - CELL_RANGE(3), CELL_RANGE(3)

                    ! Skip the interaction of a particle with itself in the same
                    ! cell and don't double count interactions within the unit
                    ! cell.
                    IF (PER_X .EQ. 0 .AND. PER_Y .EQ. 0 .AND. PER_Z .EQ. 0 &
                        .AND. IDI .EQ. IDJ) CYCLE

                    ! Adjust the j coordinates for the periodic image
                    XJ(1) = COORDS_J(1) + PER_X
                    XJ(2) = COORDS_J(2) + PER_Y
                    XJ(3) = COORDS_J(3) + PER_Z

                    ! Find the distance and scale by the box size
                    RIJ(:) = (COORDS_I(:) - XJ(:))*BOX_SIZE(:)
                    MODRIJ = SQRT(DOT_PRODUCT(RIJ(:), RIJ(:)))

                    IF (MODRIJ .LT. LJ_GAUSS_RCUT) THEN
                        ! Add energy and gradients
                        ! Divide by 2 for images of the same atom, to avoid
                        ! double counting
                        ENERGY = ENERGY + MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                          CALC_ENERGY(MODRIJ)

                        IF (GTEST) THEN
                            ! Divide gradient by 2 for images of the same atom,
                            ! to avoid double counting
                            DVDR = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                   CALC_GRAD(MODRIJ)
                            ! We must undo the box size coordinate scaling
                            GRAD_I(:) = GRAD_I(:) + &
                                        DVDR*RIJ(:)*BOX_SIZE(:)/MODRIJ
                            GRAD_J(:) = GRAD_J(:) - &
                                        DVDR*RIJ(:)*BOX_SIZE(:)/MODRIJ

                            IF (PRESENT(GRAD_BOX)) THEN
                                ! Lattice derivatives
                                GRAD_BOX(:) = GRAD_BOX(:) + &
                                         DVDR*RIJ(:)*RIJ(:)/(BOX_SIZE(:)*MODRIJ)
                            END IF ! lattice derivatievs
                        END IF ! GTEST
                    END IF ! End if less than cutoff
                END DO ! z loop
            END DO ! y loop
        END DO ! x loop

    END SUBROUTINE LJ_GAUSS_PER

!-------------------------------------------------------------------------------
!
! Checks whether the unit cell angles are allowed. It is possible to generate
! sets of angles which do not correspond to a physically possible unit cell.
! This would manifest by an imaginary unit cell volume. See
! https://journals.iucr.org/a/issues/2011/01/00/au5114/au5114.pdf
! 'On the allowed values for the triclinic unit-cell angles',
! J. Foadi and G. Evans, Acta Cryst. (2011) A67, 93–95
! ANGLES: the three unit cell angles
!
!-------------------------------------------------------------------------------

    PURE LOGICAL FUNCTION CHECK_ANGLES(ANGLES)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: ANGLES(3)

        ! SUMS: sums and differences of the angles
        DOUBLE PRECISION :: SUMS(4)

        ! Calculate the necessary sums
        SUMS(1) =   ANGLES(1) + ANGLES(2) + ANGLES(3)
        SUMS(2) = - ANGLES(1) + ANGLES(2) + ANGLES(3)
        SUMS(3) =   ANGLES(1) - ANGLES(2) + ANGLES(3)
        SUMS(4) =   ANGLES(1) + ANGLES(2) - ANGLES(3)

        ! Check all sums 0 < SUM < 2 pi and all angles 0 < ANGLE < pi
        CHECK_ANGLES = ALL(SUMS .GT. 0.D0) .AND. ALL(SUMS .LT. 2*PI) .AND. &
                       ALL(ANGLES .GT. 0.D0) .AND. ALL(ANGLES .LT. PI)

    END FUNCTION CHECK_ANGLES

!-------------------------------------------------------------------------------
!
! Finds the energy and gradients for the periodic system with a triclinic
! varying unit cell. Optionally with varying potential parameters
! IDI, IDJ: id numbers of the particles
! COORDS: coordinates of particles (in crystallographic space)
!         and lattice parameters
! GRAD_I: gradients for particle I
! GRAD_J: gradients for particle J
! ENERGY: energy contribution for this pair
! RIJ: mininum length vector between the particles
! GTEST: whether to calculate gradients
! GRAD_BOX: gradients for the box (angles, then lengths)
! BD: factors for the lattice derivatives
! GRAD_PARAMS: optional, derivatives of potential parameters
!
!-------------------------------------------------------------------------------

    SUBROUTINE LJ_GAUSS_PER_TRI(IDI, IDJ, COORDS, GRAD_I, GRAD_J, ENERGY, &
                                GTEST, GRAD_BOX, BD, GRAD_PARAMS)

        ! NATOMS: number of coordinates (including lattice parameters)
        USE COMMONS, ONLY: NATOMS

        ! BOX_DERIV: type for storing box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: IDI, IDJ
        DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
        DOUBLE PRECISION, INTENT(INOUT) :: GRAD_I(3), GRAD_J(3), ENERGY
        DOUBLE PRECISION, INTENT(INOUT) :: GRAD_BOX(6)
        LOGICAL, INTENT(IN) :: GTEST
        TYPE(BOX_DERIV), INTENT(IN) :: BD
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: GRAD_PARAMS(3)

        ! RJ: particle coordinates in crystallographic space
        ! RIJ: particle-particle vector in crystallographic space
        ! YIJ: particle-particle vector in orthogonal space
        ! MODYIJ: particle-particle distance in orthogonal space
        ! DVDY: derivative of pair potential wrt MODYIJ
        ! TEMP_GRAD: temporary gradient store
        DOUBLE PRECISION :: RJ(3), RIJ(3), YIJ(3), MODYIJ
        DOUBLE PRECISION :: DVDY, TEMP_GRAD(6)

        ! PER_k: integers for looping over cell possibilities
        ! CELL_RANGE: number of cells in each direction required
        ! I, J: loop indices
        INTEGER :: PER_X, PER_Y, PER_Z
        INTEGER :: I, J

        ! Loop over different coordinate possibilities.
        DO PER_X = - BD%CELL_RANGE(1), BD%CELL_RANGE(1)
            DO PER_Y = - BD%CELL_RANGE(2), BD%CELL_RANGE(2)
                DO PER_Z = - BD%CELL_RANGE(3), BD%CELL_RANGE(3)

                    ! Skip the interaction of a particle with itself in the same
                    ! cell and don't double count interactions within the unit
                    ! cell.
                    IF (PER_X .EQ. 0 .AND. PER_Y .EQ. 0 .AND. PER_Z .EQ. 0 &
                        .AND. IDI .EQ. IDJ) CYCLE

                    ! Adjust the j coordinates for the periodic image
                    RJ(1) = COORDS(3*IDJ-2) + PER_X
                    RJ(2) = COORDS(3*IDJ-1) + PER_Y
                    RJ(3) = COORDS(3*IDJ  ) + PER_Z

                    ! Find the distance in crystallographic space
                    RIJ(:) = (COORDS(3*IDI-2:3*IDI) - RJ(:))

                    ! Find the distance in orthogonal space
                    YIJ(:) = MATMUL(BD%ORTHOG(:,:), RIJ(:))
                    MODYIJ = SQRT(DOT_PRODUCT(YIJ(:), YIJ(:)))

                    IF (MODYIJ .LT. LJ_GAUSS_RCUT) THEN
                        ! Add energy and gradients
                        ! Divide by 2 for images of the same atom, to avoid
                        ! double counting
                        ENERGY = ENERGY + MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                          CALC_ENERGY(MODYIJ)

                        IF (GTEST) THEN
                            ! Divide gradient by 2 for images of the same atom,
                            ! to avoid double counting
                            DVDY = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                   CALC_GRAD(MODYIJ)

                            ! Atom position derivatives
                            DO I = 1, 3
                                TEMP_GRAD(I) = DVDY / MODYIJ * &
                                             DOT_PRODUCT(YIJ(:), BD%ORTHOG(:,I))
                            END DO
                            GRAD_I(:) = GRAD_I(:) + TEMP_GRAD(1:3)
                            GRAD_J(:) = GRAD_J(:) - TEMP_GRAD(1:3)

                            ! Lattice parameter derivatives
                            TEMP_GRAD(:) = 0.D0
                            DO I = 1, 6 ! Loop over the different box parameters
                                DO J = 1, 3 ! Loop over x, y, z contributions
                                    TEMP_GRAD(I) = TEMP_GRAD(I) + &
                                    YIJ(J)*DOT_PRODUCT(BD%DERIV(J,:,I), RIJ(:))
                                END DO
                            END DO
                            TEMP_GRAD(:) = TEMP_GRAD(:) * DVDY / MODYIJ
                            GRAD_BOX(:) = GRAD_BOX(:) + TEMP_GRAD(:)

                            ! If necessary, calculate the parameter derivatives
                            ! Divide gradient by 2 for images of the same atom,
                            ! to avoid double counting
                            IF(PRESENT(GRAD_PARAMS) .AND. LJ_GAUSS_MODE .EQ. 3)&
                                GRAD_PARAMS(:) = GRAD_PARAMS(:) + &
                                    MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                    CALC_GRAD_PARAMS(MODYIJ)

                        END IF ! GTEST
                    END IF ! End if less than cutoff
                END DO ! z loop
            END DO ! y loop
        END DO ! x loop

    END SUBROUTINE LJ_GAUSS_PER_TRI

!-------------------------------------------------------------------------------
!
! Calculates the contribution to the derivatives of the potential wrt the
! potential parameters LJ_GAUSS_EPSILON and LJ_GAUSS_R0
! MODYIJ: distance between the particles for this contribution
!
!-------------------------------------------------------------------------------

    PURE FUNCTION CALC_GRAD_PARAMS(MODYIJ)

        IMPLICIT NONE

        DOUBLE PRECISION :: CALC_GRAD_PARAMS(3)
        DOUBLE PRECISION, INTENT(IN) :: MODYIJ

        ! DVDE: derivative of pair potential wrt LJ_GAUSS_EPS
        ! DSFDE: derivative of SCL_FCT wrt LJ_GAUSS_EPS
        ! DVDR0: derivative of pair potential wrt LJ_GAUSS_R0
        ! DVDS: derivative of pair potential wrt LJ_GAUSS_SIGMASQ
        ! DSFDS: derivative of SCL_FCT wrt LJ_GAUSS_SIGMASQ
        ! PREF_RIJ, EXP_RIJ: factors for the Guass potential
        ! V: the pair potential energy
        DOUBLE PRECISION :: DVDE, DSFDE, DVDR0, DVDS, DSFDS
        DOUBLE PRECISION :: PREF_RIJ, EXP_RIJ, V

        ! Initialise output
        CALC_GRAD_PARAMS(:) = 0.D0

        ! Gauss factors
        PREF_RIJ = (MODYIJ - LJ_GAUSS_R0) / LJ_GAUSS_SIGMASQ
        EXP_RIJ = LJ_GAUSS_EPS * EXP(- (MODYIJ - LJ_GAUSS_R0)**2 / &
                                       (2.D0*LJ_GAUSS_SIGMASQ))

        ! Current energy
        V = CALC_ENERGY(MODYIJ)

        IF (BTEST(LJ_GAUSS_PARAMS, 0)) THEN
            ! Derivative of the pair potential wrt LJ_GAUSS_EPS just requires
            ! dividing by LJ_GAUSS_EPS
            DVDE = CALC_ENERGY_GAUSS(MODYIJ) / LJ_GAUSS_EPS

            ! Derivative of the scale factor wrt LJ_GAUSS_EPS
            DSFDE = - SQRT(2.D0 * PI * LJ_GAUSS_SIGMASQ) * SCL_FCT**2

            ! Now we have everything we need for the LJ_GAUSS_EPS derivative
            CALC_GRAD_PARAMS(1) = SCL_FCT * DVDE + DSFDE * V
        END IF

        IF (BTEST(LJ_GAUSS_PARAMS, 1)) THEN
            ! Next we look at the derivative wrt LJ_GAUSS_R0. This is simpler
            ! because SCL_FCT does not depend on LJ_GAUSS_R0, but DVDR0 is more
            ! complicated than DVDE
            ! Gaussian derivative
            DVDR0 = - PREF_RIJ * EXP_RIJ
            ! Correction to make potential got to zero at RCUT
            DVDR0 = DVDR0 + PREF_RCUT * EXP_RCUT
            ! Correction to make the derivate go to zero at RCUT
            DVDR0 = DVDR0 + (1.D0 / LJ_GAUSS_SIGMASQ - PREF_RCUT**2) * &
                             EXP_RCUT * (MODYIJ - LJ_GAUSS_RCUT)
            ! Correction to make the second derivative go to zero at RCUT
!            DVDR0 = DVDR0 - (3.D0 * PREF_RCUT/LJ_GAUSS_SIGMASQ - PREF_RCUT**3)*&
!                             EXP_RCUT * 0.5D0 * (MODYIJ - LJ_GAUSS_RCUT)**2

            ! Now we can do the LJ_GAUSS_R0 derivative
            CALC_GRAD_PARAMS(2) = SCL_FCT * DVDR0
        END IF

        IF (BTEST(LJ_GAUSS_PARAMS, 2)) THEN
            ! Finally we do the derivative wrt LJ_GAUSS_SIGMASQ. This one is
            ! really complicated...
            ! Gaussian derivative
            DVDS = - PREF_RIJ**2 * EXP_RIJ * 0.5D0
            ! Correction to make potential got to zero at RCUT
            DVDS = DVDS + PREF_RCUT**2 * EXP_RCUT * 0.5D0
            ! Correction to make the derivative go to zero at RCUT
            DVDS = DVDS + PREF_RCUT * EXP_RCUT * (MODYIJ - LJ_GAUSS_RCUT) * &
                   (1.D0 / LJ_GAUSS_SIGMASQ - PREF_RCUT**2 * 0.5D0)
            ! Correction to make the second derivative go to zero at RCUT
!            DVDS = DVDS + EXP_RCUT * 0.5D0 * (MODYIJ - LJ_GAUSS_RCUT)**2 * &
!                   (1.D0 / LJ_GAUSS_SIGMASQ**2 + PREF_RCUT**4 * 0.5D0 - &
!                    5.D0 * PREF_RCUT**2 / (2.D0 * LJ_GAUSS_SIGMASQ))

            ! Derivative of SCL_FCT wrt LJ_GAUSS_SIGMASQ
            DSFDS = - SCL_FCT**2 * PI * LJ_GAUSS_EPS / &
                      SQRT(2.D0 * PI * LJ_GAUSS_SIGMASQ)

            ! Now we have all the necessary factors
            CALC_GRAD_PARAMS(3) = SCL_FCT * DVDS + DSFDS * V
        END IF

    END FUNCTION CALC_GRAD_PARAMS

!-------------------------------------------------------------------------------
!
! Repel the unit cell volume from approaching 0, which repels the unit cell from
! adopting physically impossible angle combinations. Uses a WCA potential
! ENERGY: energy of the system
! GRAD_ANGLES: gradients for the box angles
! BD: information on the box and its derivatives
! GTEST: whether to calculate gradients
!
!-------------------------------------------------------------------------------

    SUBROUTINE CONSTRAIN_VOLUME(ENERGY, GRAD_ANGLES, BD, GTEST)

        ! BOX_DERIV: type for storing box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: ENERGY, GRAD_ANGLES(3)
        TYPE(BOX_DERIV), INTENT(IN) :: BD
        LOGICAL, INTENT(IN) :: GTEST

        ! Add a purely repulsive (away from zero) WCA energy term
        IF (BD%V .LT. V_SIGMA * 2.D0**(1.D0/6.D0)) THEN
            ENERGY = ENERGY + 4.D0 * V_EPS * &
                     ((V_SIGMA / BD%V)**(12) - (V_SIGMA / BD%V)**(6)) + V_EPS

            IF (GTEST) THEN
                ! Add the gradients, if necessary
                GRAD_ANGLES(:) = GRAD_ANGLES(:) + 24.D0 * V_EPS / V_SIGMA * &
                    ((V_SIGMA / BD%V)**(7) - 2.D0 * (V_SIGMA / BD%V)**(13)) * &
                    BD%V_DERIV(:)
            END IF
        END IF

    END SUBROUTINE CONSTRAIN_VOLUME

!-------------------------------------------------------------------------------
!
! Impose limits on the allowed range of the potential parameters LJ_GAUSS_R0 and
! LJ_GAUSS_EPSILON, by adding a harmonic repulsion outside the allowed range.
! We divide this energy contribution by SCL_FCT (it will be multiplied later) so
! the constraint is indepenednt of SCL_FCT.
! ENERGY: energy of the system
! GRAD_PARAMS: gradients for the potential parameters
! GTEST: whether to calculate gradients
!
!-------------------------------------------------------------------------------
    SUBROUTINE CONSTRAIN_PARAMETERS(ENERGY, GRAD_PARAMS, GTEST)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: ENERGY, GRAD_PARAMS(3)
        LOGICAL, INTENT(IN) :: GTEST

        ! EN_CON: constraint contribution to the energy
        DOUBLE PRECISION :: EN_CON

        ! Initialise output variable
        EN_CON = 0.D0

        ! Epsilon energy contribution
        IF (BTEST(LJ_GAUSS_PARAMS, 0)) THEN
            EN_CON = EN_CON + 0.5D0 * EPS_REPULSION * &
                (MERGE(LJ_GAUSS_EPS-EPS_UPPER, 0.D0, &
                       LJ_GAUSS_EPS .GT. EPS_UPPER)**2 + &
                 MERGE(LJ_GAUSS_EPS-EPS_LOWER, 0.D0, &
                       LJ_GAUSS_EPS .LT. EPS_LOWER)**2)
            IF (GTEST) THEN
            ! Epsilon gradient contribution
                GRAD_PARAMS(1) = GRAD_PARAMS(1) + EPS_REPULSION * &
                    (MERGE(LJ_GAUSS_EPS-EPS_UPPER, 0.D0, &
                           LJ_GAUSS_EPS .GT. EPS_UPPER) + &
                     MERGE(LJ_GAUSS_EPS-EPS_LOWER, 0.D0, &
                           LJ_GAUSS_EPS .LT. EPS_LOWER))
            END IF
        END IF

        ! R0 energy contribution
        IF (BTEST(LJ_GAUSS_PARAMS, 1)) THEN
            EN_CON = EN_CON + 0.5D0 * R0_REPULSION * &
                (MERGE(LJ_GAUSS_R0 - R0_UPPER, 0.D0, &
                       LJ_GAUSS_R0 .GT. R0_UPPER)**2 + &
                 MERGE(LJ_GAUSS_R0 - R0_LOWER, 0.D0, &
                       LJ_GAUSS_R0 .LT. R0_LOWER)**2)
            IF (GTEST) THEN
                ! R0 gradient contribution
                GRAD_PARAMS(2) = GRAD_PARAMS(2) + R0_REPULSION * &
                    (MERGE(LJ_GAUSS_R0 - R0_UPPER, 0.D0, &
                          LJ_GAUSS_R0 .GT. R0_UPPER) + &
                     MERGE(LJ_GAUSS_R0 - R0_LOWER, 0.D0, &
                           LJ_GAUSS_R0 .LT. R0_LOWER))
            END IF
        END IF

        ! SIGMASQ energy contribution
        IF (BTEST(LJ_GAUSS_PARAMS, 2)) THEN
            EN_CON = EN_CON + 0.5D0 * SIGMASQ_REPULSION * &
                (MERGE(LJ_GAUSS_SIGMASQ - SIGMASQ_UPPER, 0.D0, &
                       LJ_GAUSS_SIGMASQ .GT. SIGMASQ_UPPER)**2 + &
                 MERGE(LJ_GAUSS_SIGMASQ - SIGMASQ_LOWER, 0.D0, &
                       LJ_GAUSS_SIGMASQ .LT. SIGMASQ_LOWER)**2)
            IF (GTEST) THEN
                ! SIGMASQ gradient contribution
                GRAD_PARAMS(3) = GRAD_PARAMS(3) + SIGMASQ_REPULSION * &
                    (MERGE(LJ_GAUSS_SIGMASQ - SIGMASQ_UPPER, 0.D0, &
                           LJ_GAUSS_SIGMASQ .GT. SIGMASQ_UPPER) + &
                     MERGE(LJ_GAUSS_SIGMASQ - SIGMASQ_LOWER, 0.D0, &
                           LJ_GAUSS_SIGMASQ .LT. SIGMASQ_LOWER))
            END IF
        END IF

        ! Scale the contribtion and add it to the energy
        ENERGY = ENERGY + EN_CON / SCL_FCT

    END SUBROUTINE CONSTRAIN_PARAMETERS

!-------------------------------------------------------------------------------
!
! Outputs the coordinates to the file ljgauss.xyz
!
!-------------------------------------------------------------------------------
    SUBROUTINE VIEW_LJ_GAUSS()

        ! NATOMS: number of coordinates
        ! NSAVE: number of minima to write
        ! BOXx: Box dimensions
        ! PERIODIC: whether periodic boundary conditions are in use
        USE COMMONS, ONLY: NATOMS, NSAVE, PERIODIC, BOXLX, BOXLY, BOXLZ

        ! QMIN: energies of the saved minima
        ! QMINP: coordinates of the saved minima
        ! FF: step at which the saved minima was first found
        USE QMODULE, ONLY: QMIN, QMINP, FF

        ! BOX_DERIV: type for storing box derivative information
        ! CALC_BOX_DERIV: calculates box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV, CALC_BOX_DERIV

        IMPLICIT NONE

        ! GETUNIT: function for fiding a free file unit
        ! FUNIT: output file unit
        ! NMOL: actual number of particles
        ! I, J: loop iterators
        INTEGER :: GETUNIT, FUNIT, NMOL, I, J

        ! BOX_SIZE: convenience array for the box sizes
        ! BOX_ANGLES: convenience array for the box angles
        ! COORD: coordinate in orthogonal space
        DOUBLE PRECISION :: BOX_SIZE(3), BOX_ANGLES(3), COORD(3)

        ! Information to transform from crystallographic to orthogonal space
        TYPE(BOX_DERIV) :: BD

        ! Set NMOL and box parameters
        SELECT CASE (LJ_GAUSS_MODE)
        CASE (0)
            NMOL = NATOMS
            IF (PERIODIC) THEN
                BOX_SIZE(1) = BOXLX
                BOX_SIZE(2) = BOXLY
                BOX_SIZE(3) = BOXLZ
            ELSE
                BOX_SIZE(:) = 1.D0
            END IF
            BOX_ANGLES(:) = PI/2.D0
        CASE (1)
            NMOL = NATOMS - 1
            BOX_ANGLES(:) = PI/2.D0
        CASE (2)
            NMOL = NATOMS - 2
        CASE (3)
            NMOL = NATOMS - 3
        END SELECT

        ! Open the file
        FUNIT = GETUNIT()
        OPEN(UNIT=FUNIT, FILE='ljgauss.xyz', STATUS='UNKNOWN')

        ! Loop over the configurations
        DO I = 1, NSAVE

            WRITE(FUNIT,'(I6)') NMOL

            WRITE(FUNIT, '(A,I6,A,F20.10,A,I8,A,F20.10)') 'Energy of minimum ',&
                I, '=', QMIN(I), ' first found at step ', FF(I), &
                ' Energy per particle = ', QMIN(I)/NMOL

            ! Set box parameters for varying box
            IF (LJ_GAUSS_MODE .NE. 0) THEN
                BOX_SIZE(:) = QMINP(I, 3*NATOMS-2:3*NATOMS)
            END IF
            IF (LJ_GAUSS_MODE .EQ. 2 .OR. LJ_GAUSS_MODE .EQ. 3) THEN
                BOX_ANGLES(:) = QMINP(I, 3*NATOMS-5:3*NATOMS-3)
                BD = CALC_BOX_DERIV(BOX_ANGLES(:), BOX_SIZE(:))
            END IF

            DO J = 1, NMOL

                ! If necessary, transform the coordinate from crytallographic
                ! space to orthogonal space
                SELECT CASE (LJ_GAUSS_MODE)
                CASE (0)
                    COORD(:) = QMINP(I, 3*J-2:3*J)
                CASE (1)
                    COORD(:) = QMINP(I, 3*J-2:3*J) * BOX_SIZE(:)
                CASE (2, 3)
                    COORD(:) = MATMUL(BD%ORTHOG(:,:), QMINP(I, 3*J-2:3*J))
                END SELECT

                IF (PERIODIC .AND. J .EQ. 1) THEN
                    WRITE(FUNIT, '(A,3F20.10,A,3(A,F20.10))') 'O ', &
                        COORD(:), ' bbox_xyz', &
                        ' 0.0', BOX_SIZE(1), ' 0.0', BOX_SIZE(2), &
                        ' 0.0', BOX_SIZE(3)
                ELSE
                    WRITE(FUNIT, '(A,3F20.10)') 'O ', COORD(:)
                END IF

            END DO ! Loop over particles

        END DO ! Loop over the configurations

        CLOSE(FUNIT)

    END SUBROUTINE VIEW_LJ_GAUSS

END MODULE LJ_GAUSS_MOD
