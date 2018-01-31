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
! an oscillating pair potential (OPP). The possibility of varying
! and optimising the potential parameters is included.
! Also included is a framework for a varying triclinic unit cell. This could be
! separated out and generalised to all isotropic pair potentials.
! Questions to jwrm2.
!
MODULE OPP_MOD

! Public parameters
PUBLIC :: OPP_MODE, OPP_RCUT, OPP_RON, OPP_K, OPP_PHI, OPP_PARAMS
PUBLIC :: OPP_K_LOWER, OPP_K_UPPER
! Public subroutines
PUBLIC :: OPP, OPP_INITIALISE, OPP_TAKESTEP, OPP_VIEW
! Private parameters
PRIVATE :: V_EPS, V_SIGMA, IMAGE_CUTOFF
PRIVATE :: MAX_LENGTH_STEP, MAX_ANGLE_STEP, MAX_K_STEP, MAX_PHI_STEP
PRIVATE :: K_REPULSION
! Private subroutines
PRIVATE :: REJECT, CALC_FCTS, PERIODIC_RESET, OPP_PER, OPP_PER_TRI, HESS_CONTRIB
PRIVATE :: HESS_CONTRIB_PER, CALC_HESS_PARAMS, CALC_HESS_PARAMS_PER
PRIVATE :: CONSTRAIN_VOLUME, CONSTRAIN_PARAMETERS
! Private functions
PRIVATE :: CALC_ENERGY, CALC_GRAD, CALC_SEC, CALC_V, CALC_DVDR, CALC_D2VDR2
PRIVATE :: CALC_XPLOR, CALC_DXPLORDR, CALC_D2XPLORDR2, CHECK_ANGLES
PRIVATE :: CALC_GRAD_PARAMS, CALC_DVDK, CALC_DVDPHI, CALC_D2VDK2, CALC_D2VDPHI2
PRIVATE :: CALC_D2VDKDPHI, CALC_D2VDRDK, CALC_D2VDRDPHI

! OPP_K: frequency parameter, controls the number of wells
! OPP_PHI: phase parameter, controls the position of the first minimum
! OPP_RCUT: cutoff distance at which the energy and gradients of the potential go 
!       smoothly to zero.
! RON: distance at which the XPLOR smoothing begins
DOUBLE PRECISION :: OPP_K, OPP_PHI, OPP_RCUT, OPP_RON

! OPP_MODE: 0 standard minimisation
!           1 the final triplet of coordinates are the unit cell length
!             parameters and the second last triplet are the unit cell
!             angles, giving a triclinic unit cell, which will be optimised
!           2 as for 1, but the third last triplet are the potential
!             parameters OPP_K, OPP_PHI and 0 and
!             the potential parameters will be optimised, depending on the
!             value of OPP_PARAMS
!           3 for clusters, but with the last set of coordinates as the
!             potential parameters
! OPP_PARAMS: Specifies which of the three potential parameters will be
!             optimised. An integer between 1 and 3.
!           1 OPP_K will be optimised
!           2 OPP_PHI will be optimised
!           3 both parameters will be optimised
INTEGER :: OPP_MODE, OPP_PARAMS

! PI: the usual meaning
DOUBLE PRECISION, PARAMETER :: PI = 4*ATAN(1.D0)

! IMAGE_CUTOFF: if more periodic images than this in any direction fall
!               within the cutoff, the evaluation will not be carried
!               out and the step will be rejected.
DOUBLE PRECISION, PARAMETER :: IMAGE_CUTOFF = 30

! Parameters specifying the largest steps in OPP_TAKESTEP
! MAX_LENGTH_STEP: largest allowed step in the length of a lattice
!                  vector
! MAX_ANGLE_STEP: largest allowed step in the angle of a lattice vector
! MAX_K_STEP: largest allowed step in OPP_K
! MAX_PHI_STEP: largest allowed step in OPP_PHI
DOUBLE PRECISION, PARAMETER :: MAX_LENGTH_STEP = 0.3D0
DOUBLE PRECISION, PARAMETER :: MAX_ANGLE_STEP = 0.1D0
DOUBLE PRECISION, PARAMETER :: MAX_K_STEP = 3.D0
DOUBLE PRECISION, PARAMETER :: MAX_PHI_STEP = 0.5D0

! Parameters used for constraining the unit cell volume and the potential
! parameters
! V_EPS: scaling factor for the unit cell volume contraint energy
! V_SIGMA: WCA cutoff distance
! OPP_K_LOWER, OPP_K_UPPER: bounds on OPP_K
! K_REPULSION: harmonic force constant for OPP_K
! OPP_PHI is periodic over the range 0 to 2*PI, so doesn't need constraints
DOUBLE PRECISION, PARAMETER :: V_EPS = 1.D-3
DOUBLE PRECISION, PARAMETER :: V_SIGMA = 3.D-1
DOUBLE PRECISION :: OPP_K_LOWER, OPP_K_UPPER
DOUBLE PRECISION, PARAMETER :: K_REPULSION = 1.D4

CONTAINS

!-------------------------------------------------------------------------------
!
! Main routine for calculating the energy and gradients.
! X: coordinates array
! GRAD: gradients array
! ENERGY: calculated energy
! GTEST: whether gradients should be calculated
! STEST: whether the Hessian should be calculated
!
!-------------------------------------------------------------------------------
    SUBROUTINE OPP(X, GRAD, ENERGY, GTEST, STEST)

        ! MYUNIT: file unit for main output file
        ! NATOMS: number of particles
        ! PERIODIC: whether to use periodic boundary conditions
        USE COMMONS, ONLY: MYUNIT, NATOMS, PERIODIC

        ! HESS; the Hessian matrix
        USE MODHESS, ONLY: HESS

        ! BOX_DERIV: type for storing box derivative information
        ! CALC_BOX_DERIV: calculates box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV, CALC_BOX_DERIV

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: X(3*NATOMS)
        DOUBLE PRECISION, INTENT(OUT) :: GRAD(3*NATOMS), ENERGY
        LOGICAL, INTENT(IN) :: GTEST, STEST

        ! I, J: indices for looping over particles
        ! NMOL: actual number of particles
        INTEGER :: I, J, NMOL

        ! RIJ: Interparticle vector
        ! MODRIJ: distance between a pair of particles
        ! DVDR: derivative of pair potential wrt MODRIJ
        ! D2VDR2: D(DVDR/MODRIJ)/DMODRIJ * 1/MODRIJ
        DOUBLE PRECISION :: RIJ(3), MODRIJ, DVDR

        ! Factors for triclinic lattice
        TYPE(BOX_DERIV) :: BD

!        WRITE (MYUNIT, *) 'OPP> beginning coords'
!        DO I = 1, NATOMS
!            WRITE (MYUNIT, *) X(3*I-2:3*I)
!        END DO

        ! Initialise output variables
        ENERGY = 0.D0
        GRAD(:) = 0.D0
        IF (STEST) HESS(:,:) = 0.D0

        IF (.NOT. PERIODIC) THEN
            IF (OPP_MODE .EQ. 0) THEN
                ! Normal cluster
                NMOL = NATOMS
                DO I = 1, NMOL-1 ! Outer loop over particles
                    DO J = I+1, NMOL ! Inner loop over particles

                        ! Get the particle-particle distance
                        RIJ(:) = X(3*I-2:3*I) - X(3*J-2:3*J)
                        MODRIJ = SQRT(DOT_PRODUCT(RIJ(:), RIJ(:)))

                        ! Check against cutoff
                        IF (MODRIJ .LT. OPP_RCUT) THEN
                            ! Add energy and gradients
                            ENERGY = ENERGY + CALC_ENERGY(MODRIJ)

                            IF (GTEST) THEN
                                DVDR = CALC_GRAD(MODRIJ) / MODRIJ
                                GRAD(3*I-2:3*I) = GRAD(3*I-2:3*I) + DVDR*RIJ(:)
                                GRAD(3*J-2:3*J) = GRAD(3*J-2:3*J) - DVDR*RIJ(:)
                            END IF ! GTEST

                            IF (STEST) THEN
                                CALL HESS_CONTRIB(RIJ(:), I, J)
                            END IF ! End second derivative
                        END IF ! End less than cutoff
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles
            ELSE IF (OPP_MODE .EQ. 3) THEN
                ! Cluster with potential parameter minimisation
                CALL CALC_FCTS(X(NATOMS*3-2:NATOMS*3))
                ! Reset phi parameter
                CALL PERIODIC_RESET(X(:))
                NMOL = NATOMS - 1
                DO I = 1, NMOL-1 ! Outer loop over particles
                    DO J = I+1, NMOL ! Inner loop over particles

                        ! Get the particle-particle distance
                        RIJ(:) = X(3*I-2:3*I) - X(3*J-2:3*J)
                        MODRIJ = SQRT(DOT_PRODUCT(RIJ(:), RIJ(:)))

                        ! Check against cutoff
                        IF (MODRIJ .LT. OPP_RCUT) THEN
                            ! Add energy and gradients
                            ENERGY = ENERGY + CALC_ENERGY(MODRIJ)

                            IF (GTEST) THEN
                                DVDR = CALC_GRAD(MODRIJ) / MODRIJ
                                GRAD(3*I-2:3*I) = GRAD(3*I-2:3*I) + DVDR*RIJ(:)
                                GRAD(3*J-2:3*J) = GRAD(3*J-2:3*J) - DVDR*RIJ(:)
                                GRAD(3*NATOMS-2:3*NATOMS) = &
                                    GRAD(3*NATOMS-2:3*NATOMS) + &
                                    CALC_GRAD_PARAMS(MODRIJ)
                            END IF ! GTEST

                            IF (STEST) THEN
                                CALL HESS_CONTRIB(RIJ(:), I, J)
                            END IF ! End second derivative
                        END IF ! End less than cutoff
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

                ! Constrain the potential parameters with a harmonic repulsion
                CALL CONSTRAIN_PARAMETERS(ENERGY, GRAD(3*NATOMS-2:3*NATOMS), &
                                          GTEST, STEST)
            ELSE ! mode not equal to zero
                WRITE (MYUNIT, *) 'OPP> ERROR, PERIODIC must be ', &
                                  'specified if mode is not 0'
                STOP
            END IF
        ELSE ! PERIODIC
            SELECT CASE (OPP_MODE)
            CASE(0)
                ! Periodic, fixed box
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

                NMOL = NATOMS
                ! We need to include self-interactions of particles in different
                ! unit cells
                DO I = 1, NMOL ! Outer loop over particles
                    DO J = 1, NMOL ! Inner loop over particles
                        ! Add energy and gradients
                        CALL OPP_PER(I, J, X(3*I-2:3*I), X(3*J-2:3*J), &
                                     GRAD(3*I-2:3*I), ENERGY, GTEST, STEST)
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

            CASE(1)
                ! Triclinic varying lattice
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

!                DO I = 1, NATOMS
!                    WRITE (MYUNIT, *) X(3*I-2:3*I)
!               END DO

                ! Calculate box derivative factors
                BD = CALC_BOX_DERIV(X(3*NATOMS-5:3*NATOMS-3), &
                                    X(3*NATOMS-2:3*NATOMS  ), OPP_RCUT)

                ! Check whether the unit cell angles are physically possible.
                ! Reject if not.
                IF (.NOT. CHECK_ANGLES(X(3*NATOMS-5:3*NATOMS-3))) THEN
                    CALL REJECT(ENERGY, GRAD)
                    RETURN
                END IF

                ! Reject steps where the cell range has got bigger than the 
                ! cutoff. Such unit cells are highly distorted and there are 
                ! probably equivalent ones of a more usual shape. Allowing them
                ! leads to very slow run times.
                IF (.NOT. ALL(BD%CELL_RANGE .LE. IMAGE_CUTOFF)) THEN
                    CALL REJECT(ENERGY, GRAD)
                    RETURN
                END IF

                NMOL = NATOMS-2
                ! We need to include self-interactions of particles in different
                ! unit cells
                DO I = 1, NMOL! Outer loop over particles
                    DO J = I, NMOL ! Inner loop over particles
                        ! Add energy and gradients
                        CALL OPP_PER_TRI(I, J, X(:), GRAD(3*I-2:3*I), &
                            GRAD(3*J-2:3*J), ENERGY, GTEST, STEST, &
                            GRAD(3*NATOMS-5:3*NATOMS), BD)
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

                ! Constrain the box volume to be greater than 0 with a WCA-style
                ! repulsion
                CALL CONSTRAIN_VOLUME(ENERGY, GRAD(3*NATOMS-5:3*NATOMS-3), BD, &
                                      GTEST, STEST)

            CASE(2)
                ! Triclinic varying lattice, with varying parameters
                CALL CALC_FCTS(X(NATOMS*3-8:NATOMS*3-6))
                ! Reset all particles to within the box
                CALL PERIODIC_RESET(X(:))

                ! Calculate box derivative factors
                BD = CALC_BOX_DERIV(X(3*NATOMS-5:3*NATOMS-3), &
                                    X(3*NATOMS-2:3*NATOMS  ), OPP_RCUT)

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
                        CALL OPP_PER_TRI(I, J, X(:), GRAD(3*I-2:3*I), &
                            GRAD(3*J-2:3*J), ENERGY, GTEST, STEST, &
                            GRAD(3*NATOMS-5:3*NATOMS), BD, &
                            GRAD(3*NATOMS-8:3*NATOMS-6))
                    END DO ! Inner loop over particles
                END DO ! Outer loop over particles

                ! Constrain the box volume to be greater than 0 with a WCA-style
                ! repulsion
                CALL CONSTRAIN_VOLUME(ENERGY, GRAD(3*NATOMS-5:3*NATOMS-3), BD, &
                                      GTEST, STEST)

                ! Constrain the potential parameters with a harmonic repulsion
                CALL CONSTRAIN_PARAMETERS(ENERGY, GRAD(3*NATOMS-8:3*NATOMS-6), &
                                          GTEST, STEST)

                ! Fiddle to lock in cubic
!                X(3*NATOMS-1) = X(3*NATOMS-2)
!                X(3*NATOMS  ) = X(3*NATOMS-2)
!                GRAD(3*NATOMS-1) = GRAD(3*NATOMS-2)
!                GRAD(3*NATOMS  ) = GRAD(3*NATOMS-2)
            END SELECT ! Mode selections
        END IF ! Periodic

    END SUBROUTINE OPP

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
! Calculate factors that do depend on parameters, but not on RIJ
! Some book keeping of the potential parameters to deal with different modes
! PARAMS: the potential parameters when those are to be adjusted, otherwise
!         ignored
!
!-------------------------------------------------------------------------------

    SUBROUTINE CALC_FCTS(PARAMS)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: PARAMS(3)

        ! Adjust the potential parameters, if the mode requires it
        ! For each parameter, if it is varying, set it from the coordinate.
        ! If it's not varying, reset the coordinate, which TAKESTEP will
        ! have changed.
        IF (BTEST(OPP_PARAMS, 0)) THEN
            OPP_K = PARAMS(1)
        ELSE
            PARAMS(1) = OPP_K
        END IF
        IF (BTEST(OPP_PARAMS, 1)) THEN
            OPP_PHI = PARAMS(2)
        ELSE
            PARAMS(2) = OPP_PHI
        END IF
        PARAMS(3) = 0.D0

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

        CALC_ENERGY = CALC_V(MODRIJ) * CALC_XPLOR(MODRIJ)

    END FUNCTION CALC_ENERGY

!-------------------------------------------------------------------------------
!
! Calculates the pair potential gradient for a pair of particles.
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_GRAD (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_GRAD = CALC_DVDR(MODRIJ) * CALC_XPLOR(MODRIJ) + &
            CALC_V(MODRIJ) * CALC_DXPLORDR(MODRIJ)

    END FUNCTION CALC_GRAD

!-------------------------------------------------------------------------------
!
! Calculates the pair potential second derivative for a pair of particles.
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_SEC (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_SEC = CALC_D2VDR2(MODRIJ) * CALC_XPLOR(MODRIJ) + &
            2.D0 * CALC_DVDR(MODRIJ) * CALC_DXPLORDR(MODRIJ) + &
            CALC_V(MODRIJ) * CALC_D2XPLORDR2(MODRIJ)

    END FUNCTION CALC_SEC

!-------------------------------------------------------------------------------
!
! Calculates the potential unmodified by XPLOR
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_V (MODRIJ)

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_V = MODRIJ**(-15) + MODRIJ**(-3) * &
            COS(OPP_K * (MODRIJ - 1) - OPP_PHI)

    END FUNCTION CALC_V

!-------------------------------------------------------------------------------
!
! Calculates the gradient of the potential unmodified by XPLOR
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_DVDR (MODRIJ)

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_DVDR = - 15.D0 * MODRIJ**(-16) - &
            OPP_K * MODRIJ**(-3) * SIN(OPP_K * (MODRIJ - 1) - OPP_PHI) - & 
            3.D0 * MODRIJ**(-4) * COS(OPP_K * (MODRIJ - 1) - OPP_PHI)

    END FUNCTION CALC_DVDR

!-------------------------------------------------------------------------------
!
! Calculates the second derivative of the potential unmodified by XPLOR
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2VDR2 (MODRIJ)

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_D2VDR2 = 240.D0 * MODRIJ**(-17) - &
            OPP_K**2 * MODRIJ**(-3) * COS(OPP_K * (MODRIJ - 1) - OPP_PHI) + & 
            6.D0 * OPP_K * MODRIJ**(-4) * SIN(OPP_K * (MODRIJ - 1) - OPP_PHI) +&
            12.D0 * MODRIJ**(-5) * COS(OPP_K * (MODRIJ - 1) - OPP_PHI)

    END FUNCTION CALC_D2VDR2

!-------------------------------------------------------------------------------
!
! Calculates the XPLOR smoothing potential
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_XPLOR (MODRIJ)

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        IF (MODRIJ .LT. OPP_RON) THEN
            CALC_XPLOR = 1.D0
        ELSE IF (MODRIJ .GT. OPP_RCUT) THEN
            CALC_XPLOR = 0.D0
        ELSE
            CALC_XPLOR = (OPP_RCUT**2 - MODRIJ**2)**2 * &
                (OPP_RCUT**2 + 2.D0 * MODRIJ**2 - 3.D0 * OPP_RON**2) / &
                (OPP_RCUT**2 - OPP_RON**2)**3
        END IF

    END FUNCTION CALC_XPLOR

!-------------------------------------------------------------------------------
!
! Calculates the derivative of the XPLOR smoothing potential
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_DXPLORDR (MODRIJ)

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        IF (MODRIJ .LT. OPP_RON .OR. MODRIJ .GT. OPP_RCUT) THEN
            CALC_DXPLORDR = 0.D0
        ELSE
            CALC_DXPLORDR = 12.D0 * MODRIJ * (OPP_RCUT**2 - MODRIJ**2) * &
                (OPP_RON**2 - MODRIJ**2) / (OPP_RCUT**2 - OPP_RON**2)**3
        END IF

    END FUNCTION CALC_DXPLORDR

!-------------------------------------------------------------------------------
!
! Calculates the second derivative of the XPLOR smoothing potential
! MODRIJ: distance beteen the particle pair
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2XPLORDR2 (MODRIJ)

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        IF (MODRIJ .LT. OPP_RON .OR. MODRIJ .GT. OPP_RCUT) THEN
            CALC_D2XPLORDR2 = 0.D0
        ELSE
            CALC_D2XPLORDR2 = 12.D0 * &
                ((OPP_RCUT**2 - MODRIJ**2) * (OPP_RON**2 - MODRIJ**2) - &
                2.D0 * MODRIJ**2 * (OPP_RCUT**2 - MODRIJ**2) - &
                2.D0 * MODRIJ**2 * (OPP_RON**2 - MODRIJ**2)) / &
                (OPP_RCUT**2 - OPP_RON**2)**3
        END IF

    END FUNCTION CALC_D2XPLORDR2

!-------------------------------------------------------------------------------
!
! Initialisation. Check values for OPP_MODE and OPP_PARAMS.
!
!-------------------------------------------------------------------------------
    SUBROUTINE OPP_INITIALISE()

        ! MYUNIT: file unit for main GMIN output
        USE COMMONS, ONLY: MYUNIT

        IMPLICIT NONE

        ! Check we have sensible value for the mode
        IF (OPP_MODE .LT. 0 .OR. OPP_MODE .GT. 3 ) THEN
            WRITE (MYUNIT, *) 'OPP> ERROR: mode must be between 0 and 3'
            STOP
        END IF

        ! Check we have sensible values for the params
        IF (OPP_MODE .EQ. 2 .OR. OPP_MODE .EQ. 3) THEN
            IF (OPP_PARAMS .LT. 1 .OR. OPP_PARAMS .GT. 3) THEN
                WRITE (MYUNIT, *) 'OPP> ERROR: params must be between 1 and 3'
                STOP
            END IF
        END IF

        ! Check the XPLOR cut-off values
        IF (OPP_RCUT .LE. OPP_RON) THEN
            WRITE (MYUNIT, *) 'OPP> ERROR: XPLOR switch on must be less than', &
                              ' the cut-off'
            STOP
        END IF    

    END SUBROUTINE OPP_INITIALISE

!-------------------------------------------------------------------------------
!
! Takes a step, perturbing the Cartesian coordinates, the unit cell parameters
! and the potential parameters, according to the specified mode.
! NP: the index in the main COORDS array to take the step from
!
!-------------------------------------------------------------------------------

    SUBROUTINE OPP_TAKESTEP(NP)

        ! COORDS: full set of Markov chain coordinates
        ! FROZEN: array of atom indices not to move
        ! MYUNIT: file unit of the main GMIN output file
        ! NATOMS: number of coordinates
        ! PERIODIC: whetehr periodic boundary conditions are used
        ! PERCOLATET: whether percolation is used to constrain particles
        ! RADIUS: container radius
        ! STEP: maximum step size for this move
        ! TMOVE: specifies which step to do a translational move on
        ! TWOD: whetehr the system is two dimensional
        USE COMMONS, ONLY: COORDS, FROZEN, MYUNIT, NATOMS, PERIODIC, &
                           PERCOLATET, RADIUS, STEP, TMOVE, TWOD
        
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

!        WRITE (MYUNIT, *) 'OPP_TAKESTEP> beginning'
!        FLUSH(MYUNIT)

!        WRITE (MYUNIT, *) 'Beginning takestep, coordinates are:'
!        DO I = 1, NATOMS
!            WRITE (MYUNIT, *) COORDS(3*I-2:3*I, NP)
!        END DO

        ! Set the local copy of the coordinates
        X(:) = COORDS(:, NP)

        ! Set NMOL
        SELECT CASE(OPP_MODE)
        CASE(0)
            NMOL = NATOMS
        CASE(1)
            NMOL = NATOMS - 2
        CASE(2)
            NMOL = NATOMS - 3
        CASE(3)
            NMOL = NATOMS - 1
        END SELECT

        ! Random translational steps for each of the particles
        IF (TMOVE(NP)) THEN
            DO I = 1, NMOL
                ! Skip if frozen
                IF (FROZEN(I)) CYCLE

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
                ! Add the geenrated step onto the coordinates
                X(I*3-2:I*3) = X(I*3-2:I*3) + RAND_STEP(:)
            END DO
        END IF

        ! If not periodic and we're using a radius container, bring particles in
        IF (.NOT. PERIODIC .AND. .NOT. PERCOLATET) THEN
            DO I = 1, NMOL
                IF (VEC_LEN(X(I*3-2:I*3)) .GT. RADIUS) THEN
                    WRITE(MYUNIT, *) 'OPP_TAKESTEP> ', &
                        'coord outside container, bringing in'

                    ! Skip if frozen
                    IF (FROZEN(I)) CYCLE
                    X(I*3-2:I*3) = X(I*3-2:I*3) - &
                        SQRT(RADIUS) * NINT(X(I*3-2:I*3) / SQRT(RADIUS))
                END IF
            END DO
        END IF

        ! Box length steps
        IF ((OPP_MODE .EQ. 1 .OR. OPP_MODE .EQ. 2) .AND. &
            .NOT. FROZEN(NATOMS)) THEN
            X(3*NATOMS-2:3*NATOMS) = X(3*NATOMS-2:3*NATOMS) + &
                                     VEC_RANDOM() * MAX_LENGTH_STEP
            ! Fiddle to lock in cubic
!            X(3*NATOMS-1) = X(3*NATOMS-2)
!            X(3*NATOMS  ) = X(3*NATOMS-2)
        END IF

        ! Box angle steps
        IF ((OPP_MODE .EQ. 1 .OR. OPP_MODE .EQ. 2) .AND. &
            .NOT. FROZEN(NATOMS-1)) THEN
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
        IF (OPP_MODE .EQ. 2 .AND. .NOT. FROZEN(NATOMS-2)) THEN
            IF (BTEST(OPP_PARAMS, 0)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-8) = X(3*NATOMS-8) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_K_STEP
            END IF
            IF (BTEST(OPP_PARAMS, 1)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-7) = X(3*NATOMS-7) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_PHI_STEP
            END IF
        ELSE IF (OPP_MODE .EQ. 3 .AND. .NOT. FROZEN(NATOMS)) THEN
            IF (BTEST(OPP_PARAMS, 0)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-2) = X(3*NATOMS-2) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_K_STEP
            END IF
            IF (BTEST(OPP_PARAMS, 1)) THEN
                CALL RANDOM_NUMBER(STEP_SIZE)
                X(3*NATOMS-1) = X(3*NATOMS-1) + &
                                (STEP_SIZE - 0.5D0) * 2.D0 * MAX_PHI_STEP
            END IF            
        END IF

        COORDS(:, NP) = X(:)

!        WRITE (MYUNIT, *) 'OPP_TAKESTEP> ending'
!        FLUSH(MYUNIT)
!        WRITE (MYUNIT, *) 'Ending takestep, coordinates are:'
!        DO I = 1, NATOMS
!            WRITE (MYUNIT, *) COORDS(3*I-2:3*I, NP)
!        END DO

    END SUBROUTINE OPP_TAKESTEP

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
        IF (OPP_MODE .EQ. 1 .OR. OPP_MODE .EQ. 2) THEN
            BOXLX = COORDS(3*NATOMS - 2)
            BOXLY = COORDS(3*NATOMS - 1)
            BOXLZ = COORDS(3*NATOMS)
        END IF

        ! Reset coordinates to within the box
        SELECT CASE (OPP_MODE)
        CASE(0)
            COORDS(:) = COORDS(:) - FLOOR(COORDS(:))
        CASE(1)
            COORDS(1:3*NATOMS-6) = COORDS(1:3*NATOMS-6) - &
                                   FLOOR(COORDS(1:3*NATOMS-6))
            ! Make the lattice angles modulo 2*pi
            ! They will be checked later
            COORDS(3*NATOMS-5:3*NATOMS-3) = COORDS(3*NATOMS-5:3*NATOMS-3) - &
                2.D0*PI*FLOOR(COORDS(3*NATOMS-5:3*NATOMS-3)/(2.D0*PI))
        CASE(2)
           COORDS(1:3*NATOMS-9) = COORDS(1:3*NATOMS-9) - &
                                   FLOOR(COORDS(1:3*NATOMS-9))
            ! Make the lattice angles modulo 2*pi
            ! They will be checked later
            COORDS(3*NATOMS-5:3*NATOMS-3) = COORDS(3*NATOMS-5:3*NATOMS-3) - &
                2.D0*PI*FLOOR(COORDS(3*NATOMS-5:3*NATOMS-3)/(2.D0*PI))
            ! Make OPP_PHI modulo 2*pi
            COORDS(3*NATOMS-7) = COORDS(3*NATOMS-7) - &
                2.D0*PI*FLOOR(COORDS(3*NATOMS-7)/(2.D0*PI))
            ! Reset the unused coordinate to zero
            COORDS(3*NATOMS-6) = 0.D0
        CASE(3)
            ! Make OPP_PHI modulo 2*pi
            COORDS(3*NATOMS-1) = COORDS(3*NATOMS-1) - &
                2.D0*PI*FLOOR(COORDS(3*NATOMS-1)/(2.D0*PI))
            ! Reset the unused coordinate to zero
            COORDS(3*NATOMS) = 0.D0
        END SELECT

    END SUBROUTINE PERIODIC_RESET

!-------------------------------------------------------------------------------
!
! Finds the energy and gradients for the othogonal fixed periodic system
! IDI, IDJ: id numbers of the particles
! COORDS_I: coordinates of particle I
! COORDS_J: coordinates of particle J
! GRAD_I: gradients for particle I
! GRAD_J: gradients for particle J
! ENERGY: energy contribution for this pair
! GTEST: whether to calculate gradients
! STEST: whether to calculate Hessian
!
!-------------------------------------------------------------------------------

    SUBROUTINE OPP_PER(IDI, IDJ, COORDS_I, COORDS_J, GRAD_I, ENERGY, &
                       GTEST, STEST)

        ! BOXLX, BOXLY, BOXLZ: dimensions of the box
        USE COMMONS, ONLY: BOXLX, BOXLY, BOXLZ

        ! HESS; the Hessian matrix
        USE MODHESS, ONLY: HESS

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: IDI, IDJ
        DOUBLE PRECISION, INTENT(IN) :: COORDS_I(3), COORDS_J(3)
        DOUBLE PRECISION, INTENT(INOUT) :: GRAD_I(3), ENERGY
        LOGICAL, INTENT(IN) :: GTEST, STEST

        ! XJ: copy of j coordinates
        ! RIJ: particle-particle vector
        ! MODRIJ: particle-particle distance
        ! DVDR: derivative of pair potential wrt MODRIJ
        ! D2VDR2: D(DVDR/MODRIJ)/DMODRIJ * 1/MODRIJ
        ! BOX_SIZE: convenience copy of the box size
        DOUBLE PRECISION :: XJ(3), RIJ(3), MODRIJ, DVDR, D2VDR2, BOX_SIZE(3)

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
        CELL_RANGE(:) = CEILING(OPP_RCUT/BOX_SIZE(:))

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

                    IF (MODRIJ .LT. OPP_RCUT) THEN
                        ! Add energy and gradients
                        ! Divide by 2 for images of the same atom, to avoid
                        ! double counting
                        ! Don't double count energy
                        IF (IDJ .GE. IDI) THEN
                        ENERGY = ENERGY + &
                            MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                  CALC_ENERGY(MODRIJ)
                        END IF

                        IF (GTEST) THEN
                            ! Divide gradient by 2 for images of the same atom,
                            ! to avoid double counting
                            DVDR = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                   CALC_GRAD(MODRIJ) / MODRIJ
                            ! We must undo the box size coordinate scaling
                            GRAD_I(:) = GRAD_I(:) + DVDR * RIJ(:) * BOX_SIZE(:)
                        END IF ! GTEST

                        IF (STEST) THEN
                            ! Hessian terms
                            DVDR = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                CALC_GRAD(MODRIJ) / MODRIJ
                            D2VDR2 = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                (CALC_SEC(MODRIJ) - DVDR) / MODRIJ**2
                            ! Same atom, same coordinate
                            HESS(3*IDI-2,3*IDI-2) = HESS(3*IDI-2,3*IDI-2) + &
                                (D2VDR2 * RIJ(1)**2 + DVDR) * BOX_SIZE(1)**2
                            HESS(3*IDI-1,3*IDI-1) = HESS(3*IDI-1,3*IDI-1) + &
                                (D2VDR2 * RIJ(2)**2 + DVDR) * BOX_SIZE(2)**2
                            HESS(3*IDI  ,3*IDI  ) = HESS(3*IDI  ,3*IDI  ) + &
                                (D2VDR2 * RIJ(3)**2 + DVDR) * BOX_SIZE(3)**2
                            ! Same atom, different coordinate
                            HESS(3*IDI-2,3*IDI-1) = HESS(3*IDI-2,3*IDI-1) + &
                                D2VDR2 * RIJ(1) * RIJ(2) * &
                                BOX_SIZE(1) * BOX_SIZE(2)
                            HESS(3*IDI-1,3*IDI  ) = HESS(3*IDI-1,3*IDI  ) + &
                                D2VDR2 * RIJ(2) * RIJ(3) * &
                                BOX_SIZE(2) * BOX_SIZE(3)
                            HESS(3*IDI  ,3*IDI-2) = HESS(3*IDI  ,3*IDI-2) + &
                                D2VDR2 * RIJ(3) * RIJ(1) * &
                                BOX_SIZE(3) * BOX_SIZE(1)
                            HESS(3*IDI-1,3*IDI-2) = HESS(3*IDI-1,3*IDI-2) + &
                                D2VDR2 * RIJ(2) * RIJ(1) * &
                                BOX_SIZE(2) * BOX_SIZE(1)
                            HESS(3*IDI  ,3*IDI-1) = HESS(3*IDI  ,3*IDI-1) + &
                                D2VDR2 * RIJ(3) * RIJ(2) * &
                                BOX_SIZE(3) * BOX_SIZE(2)
                            HESS(3*IDI-2,3*IDI  ) = HESS(3*IDI-2,3*IDI) + &
                                D2VDR2 * RIJ(1) * RIJ(3) * &
                                BOX_SIZE(1) * BOX_SIZE(3)
                            ! Different atom, same coordinate
                            HESS(3*IDI-2,3*IDJ-2) = HESS(3*IDI-2,3*IDJ-2) + &
                                (- D2VDR2 * RIJ(1)**2 - DVDR) * BOX_SIZE(1)**2
                            HESS(3*IDI-1,3*IDJ-1) = HESS(3*IDI-1,3*IDJ-1) + &
                                (- D2VDR2 * RIJ(2)**2 - DVDR) * BOX_SIZE(2)**2
                            HESS(3*IDI  ,3*IDJ  ) = HESS(3*IDI  ,3*IDJ  ) + &
                                (- D2VDR2 * RIJ(3)**2 - DVDR) * BOX_SIZE(3)**2
                            ! Different atom, different coordinate
                            HESS(3*IDI-2,3*IDJ-1) = HESS(3*IDI-2,3*IDJ-1) - &
                                D2VDR2*RIJ(1)*RIJ(2)*BOX_SIZE(1)*BOX_SIZE(2)
                            HESS(3*IDI-1,3*IDJ  ) = HESS(3*IDI-1,3*IDJ  ) - & 
                                D2VDR2*RIJ(2)*RIJ(3)*BOX_SIZE(2)*BOX_SIZE(3)
                            HESS(3*IDI-2,3*IDJ  ) = HESS(3*IDI-2,3*IDJ  ) - & 
                                D2VDR2*RIJ(3)*RIJ(1)*BOX_SIZE(3)*BOX_SIZE(1)
                            HESS(3*IDJ-1,3*IDI-2) = HESS(3*IDJ-1,3*IDI-2) - &
                                D2VDR2*RIJ(1)*RIJ(2)*BOX_SIZE(2)*BOX_SIZE(1)
                            HESS(3*IDJ  ,3*IDI-1) = HESS(3*IDJ  ,3*IDI-1) - &
                                D2VDR2*RIJ(2)*RIJ(3)*BOX_SIZE(3)*BOX_SIZE(2)
                            HESS(3*IDJ  ,3*IDI-2) = HESS(3*IDJ  ,3*IDI-2) - &
                                D2VDR2*RIJ(3)*RIJ(1)*BOX_SIZE(1)*BOX_SIZE(3)
                        END IF ! End second derivative
                    END IF ! End if less than cutoff
                END DO ! z loop
            END DO ! y loop
        END DO ! x loop

    END SUBROUTINE OPP_PER

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
! STEST: whether to calculate Hessian
! GRAD_BOX: gradients for the box (angles, then lengths)
! BD: factors for the lattice derivatives
! GRAD_PARAMS: optional, derivatives of potential parameters
!
!-------------------------------------------------------------------------------

    SUBROUTINE OPP_PER_TRI(IDI, IDJ, COORDS, GRAD_I, GRAD_J, ENERGY, GTEST, &
                           STEST, GRAD_BOX, BD, GRAD_PARAMS)

        ! NATOMS: number of coordinates (including lattice parameters)
        USE COMMONS, ONLY: NATOMS

        ! BOX_DERIV: type for storing box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: IDI, IDJ
        DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
        DOUBLE PRECISION, INTENT(INOUT) :: GRAD_I(3), GRAD_J(3), ENERGY
        DOUBLE PRECISION, INTENT(INOUT) :: GRAD_BOX(6)
        LOGICAL, INTENT(IN) :: GTEST, STEST
        TYPE(BOX_DERIV), INTENT(IN) :: BD
        DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: GRAD_PARAMS(3)

        ! RJ: particle coordinates in crystallographic space
        ! RIJ: particle-particle vector in crystallographic space
        ! YIJ: particle-particle vector in orthogonal space
        ! MODYIJ: particle-particle distance in orthogonal space
        ! DVDY: derivative of pair potential wrt MODYIJ / MODYIJ
        ! TEMP_GRAD: temporary gradient store
        ! TEMP: temporary calculation store
        DOUBLE PRECISION :: RJ(3), RIJ(3), YIJ(3), MODYIJ
        DOUBLE PRECISION :: DVDY, TEMP_GRAD(3), TEMP

        ! PER_k: integers for looping over cell possibilities
        ! CELL_RANGE: number of cells in each direction required
        ! I: loop index
        INTEGER :: PER_X, PER_Y, PER_Z
        INTEGER :: I

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

                    IF (MODYIJ .LT. OPP_RCUT) THEN
                        ! Add energy and gradients
                        ! Divide by 2 for images of the same atom, to avoid
                        ! double counting
                        ENERGY = ENERGY + MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                  CALC_ENERGY(MODYIJ)

                        ! Gradients
                        IF (GTEST) THEN
                            ! Divide gradient by 2 for images of the same atom,
                            ! to avoid double counting
                            DVDY = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                   CALC_GRAD(MODYIJ) / MODYIJ

                            ! Atom position derivatives
                            TEMP_GRAD(1:3) = DVDY * &
                                MATMUL(TRANSPOSE(BD%ORTHOG(:,:)), YIJ(:))
                            GRAD_I(:) = GRAD_I(:) + TEMP_GRAD(1:3)
                            GRAD_J(:) = GRAD_J(:) - TEMP_GRAD(1:3)

                            ! Lattice parameter derivatives
                            DO I = 1, 6 ! Loop over the different box parameters
                                TEMP = DVDY * &
                                    DOT_PRODUCT(YIJ(:), &
                                    MATMUL(BD%DERIV(:,:,I), RIJ(:)))
                                GRAD_BOX(I) = GRAD_BOX(I) + TEMP
                            END DO

                            ! If necessary, calculate the parameter derivatives
                            ! Divide gradient by 2 for images of the same atom,
                            ! to avoid double counting
                            IF(PRESENT(GRAD_PARAMS) .AND. OPP_MODE .EQ. 2)&
                                GRAD_PARAMS(:) = GRAD_PARAMS(:) + &
                                    MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                                    CALC_GRAD_PARAMS(MODYIJ)

                        END IF ! GTEST

                        IF (STEST) THEN
                            CALL HESS_CONTRIB_PER(RIJ(:), BD, IDI, IDJ)
                        END IF ! End second derivative
                    END IF ! End if less than cutoff
                END DO ! z loop
            END DO ! y loop
        END DO ! x loop

    END SUBROUTINE OPP_PER_TRI

!-------------------------------------------------------------------------------
!
! Calculates the second derivative contributions for a cluster
! RIJ: vector between particles in fractional coordinates
! BD: box derivative information
! IDI, IDJ: indices of the atom pair
!
!-------------------------------------------------------------------------------
    SUBROUTINE HESS_CONTRIB(RIJ, IDI, IDJ)

        ! HESS: the Hessian matrix
        USE MODHESS, ONLY: HESS

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: RIJ(3)
        INTEGER, INTENT(IN) :: IDI, IDJ

        ! MODRIJ: particle-particle distance in Cartesian space
        ! DVDR: derivative of pair potential wrt MODYIJ
        ! D2VDR2: D(DVDR/MODRIJ)/DMODRIJ * 1/MODRIJ
        DOUBLE PRECISION :: MODRIJ, DVDR, D2VDR2

        ! I, J: loop indices
        INTEGER :: I, J

        ! Calculate some factors
        MODRIJ = SQRT(DOT_PRODUCT(RIJ(:), RIJ(:)))
        DVDR = CALC_GRAD(MODRIJ) / MODRIJ
        D2VDR2 = (CALC_SEC(MODRIJ) - DVDR) / MODRIJ**2

        DO I = 1, 3 ! outer loop over x, y, z
            DO J = 1, 3 ! inner loop over x, y, z
                ! Same atom, same or different coordinate
                HESS(3*IDI-3+I,3*IDI-3+J) = HESS(3*IDI-3+I,3*IDI-3+J) + &
                    D2VDR2 * RIJ(I) * RIJ(J)
                HESS(3*IDJ-3+I,3*IDJ-3+J) = HESS(3*IDJ-3+I,3*IDJ-3+J) + &
                    D2VDR2 * RIJ(I) * RIJ(J)
                ! Same atom, same coordinate
                IF (I .EQ. J) THEN
                    HESS(3*IDI-3+I,3*IDI-3+J) = HESS(3*IDI-3+I,3*IDI-3+J) + DVDR
                    HESS(3*IDJ-3+I,3*IDJ-3+J) = HESS(3*IDJ-3+I,3*IDJ-3+J) + DVDR
                END IF
                ! Different atom, same or different coordinate
                HESS(3*IDI-3+I,3*IDJ-3+J) = HESS(3*IDI-3+I,3*IDJ-3+J) - &
                    D2VDR2 * RIJ(I) * RIJ(J)
                HESS(3*IDJ-3+I,3*IDI-3+J) = HESS(3*IDJ-3+I,3*IDI-3+J) - &
                    D2VDR2 * RIJ(I) * RIJ(J)
                ! Same atom, same coordinate
                IF (I .EQ. J) THEN
                    HESS(3*IDI-3+I,3*IDJ-3+J) = HESS(3*IDI-3+I,3*IDJ-3+J) - DVDR
                    HESS(3*IDJ-3+I,3*IDI-3+J) = HESS(3*IDJ-3+I,3*IDI-3+J) - DVDR
                END IF
            END DO ! end inner loop over x, y, z
        END DO ! end outer loop over x, y, z

        IF (OPP_MODE .EQ. 3) THEN
            ! Potential parameter contributions
            CALL CALC_HESS_PARAMS(RIJ, IDI, IDJ)
        END IF

    END SUBROUTINE HESS_CONTRIB

!-------------------------------------------------------------------------------
!
! Calculates the second derivative contributions for a periodic system
! RIJ: vector between particles in fractional coordinates
! BD: box derivative information
! IDI, IDJ: indices of the atom pair
!
!-------------------------------------------------------------------------------
    SUBROUTINE HESS_CONTRIB_PER(RIJ, BD, IDI, IDJ)

        ! Number of atoms, including unit cell parameters
        USE COMMONS, ONLY: NATOMS

        ! HESS: the Hessian matrix
        USE MODHESS, ONLY: HESS

        ! BOX_DERIV: type for storing box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: RIJ(3)
        TYPE(BOX_DERIV), INTENT(IN) :: BD
        INTEGER, INTENT(IN) :: IDI, IDJ

        ! YIJ: particle-particle vector in Cartesian space
        ! MODYIJ: particle-particle distance in Cartesian space
        ! DVDY: derivative of pair potential wrt MODYIJ
        ! D2VDY2: D(DVDY/MODYIJ)/DMODYIJ * 1/MODYIJ
        ! TEMP, TEMP_GRAD: temporary calculation store
        DOUBLE PRECISION :: YIJ(3), MODYIJ, DVDY, D2VDY2, TEMP, TEMP_GRAD(3)

        ! I, J: loop indices
        INTEGER :: I, J

        ! Calculate some factors
        YIJ(:) = MATMUL(BD%ORTHOG(:,:), RIJ(:))
        MODYIJ = SQRT(DOT_PRODUCT(YIJ(:), YIJ(:)))
        DVDY = MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * CALC_GRAD(MODYIJ) / MODYIJ
        D2VDY2 = (MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * CALC_SEC(MODYIJ) - DVDY) &
            / MODYIJ**2

        ! Pure atom translations
        ! Loop over the atom coordinates twice
        DO I = 1, 3
            DO J = 1, 3
                TEMP = D2VDY2 * DOT_PRODUCT(YIJ(:),BD%ORTHOG(:,I)) * &
                    DOT_PRODUCT(YIJ(:),BD%ORTHOG(:,J)) + DVDY * BD%ORTHOGSQ(I,J)
                ! Same atom, same or different coord
                HESS(3*IDI-3+I,3*IDI-3+J) = HESS(3*IDI-3+I,3*IDI-3+J) + TEMP
                HESS(3*IDJ-3+I,3*IDJ-3+J) = HESS(3*IDJ-3+I,3*IDJ-3+J) + TEMP
                ! Different atom, same or different coord
                HESS(3*IDI-3+I,3*IDJ-3+J) = HESS(3*IDI-3+I,3*IDJ-3+J) - TEMP
                HESS(3*IDJ-3+J,3*IDI-3+I) = HESS(3*IDJ-3+J,3*IDI-3+I) - TEMP
            END DO
        END DO

        ! Mixed atom and unit cell parameter
        DO I = 1, 6 ! loop over unit cell parameters
            TEMP_GRAD(:) = D2VDY2 * &
                DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV(:,:,I), RIJ(:))) * &
                MATMUL(TRANSPOSE(BD%ORTHOG(:,:)), YIJ(:)) + &
                DVDY * (MATMUL(TRANSPOSE(BD%DERIV(:,:,I)), YIJ(:)) + &
                MATMUL(TRANSPOSE(BD%ORTHOG(:,:)), &
                MATMUL(BD%DERIV(:,:,I), RIJ(:))))
            HESS(3*IDI-2:3*IDI,3*NATOMS-6+I) = &
                HESS(3*IDI-2:3*IDI,3*NATOMS-6+I) + TEMP_GRAD(:)
            HESS(3*NATOMS-6+I,3*IDI-2:3*IDI) = &
                HESS(3*NATOMS-6+I,3*IDI-2:3*IDI) + TEMP_GRAD(:)
            HESS(3*IDJ-2:3*IDJ,3*NATOMS-6+I) = &
                HESS(3*IDJ-2:3*IDJ,3*NATOMS-6+I) - TEMP_GRAD(:)
            HESS(3*NATOMS-6+I,3*IDJ-2:3*IDJ) = &
                HESS(3*NATOMS-6+I,3*IDJ-2:3*IDJ) - TEMP_GRAD(:)
        END DO ! end loop over unit cell parameters

        ! Pure unit cell parameters
        DO I = 1, 6 ! outer loop over cell parameters
            ! Diagonals
            TEMP = D2VDY2 * &
                DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV(:,:,I), RIJ(:)))**2 + &
                DVDY * (DOT_PRODUCT(MATMUL(BD%DERIV(:,:,I), RIJ(:)), &
                MATMUL(BD%DERIV(:,:,I), RIJ(:))) + &
                DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV2(:,:,I,I), RIJ(:))))
            HESS(3*NATOMS-6+I,3*NATOMS-6+I) = HESS(3*NATOMS-6+I,3*NATOMS-6+I) +&
                TEMP
            DO J = I + 1, 6 ! inner loop over cell parameters
                TEMP = D2VDY2 * &
                    DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV(:,:,I), RIJ(:))) * &
                    DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV(:,:,J), RIJ(:))) + &
                    DVDY * (DOT_PRODUCT(MATMUL(BD%DERIV(:,:,I), RIJ(:)), &
                    MATMUL(BD%DERIV(:,:,J), RIJ(:))) + &
                    DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV2(:,:,I,J), RIJ(:))))
                HESS(3*NATOMS-6+I,3*NATOMS-6+J) = &
                    HESS(3*NATOMS-6+I,3*NATOMS-6+J) + TEMP
                HESS(3*NATOMS-6+J,3*NATOMS-6+I) = &
                    HESS(3*NATOMS-6+J,3*NATOMS-6+I) + TEMP
            END DO ! end inner loop over cell parameters
        END DO ! end outer loop over cell parameters

        IF (OPP_MODE .EQ. 2) THEN
            ! Potential parameter contributions
            CALL CALC_HESS_PARAMS_PER(RIJ, BD, IDI, IDJ)
        END IF

    END SUBROUTINE HESS_CONTRIB_PER

!-------------------------------------------------------------------------------
!
! Calculates the contribution to the derivatives of the potential wrt the
! potential parameters OPP_K and OPP_PHI
! MODRIJ: distance between the particles for this contribution
!
!-------------------------------------------------------------------------------

    PURE FUNCTION CALC_GRAD_PARAMS(MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION :: CALC_GRAD_PARAMS(3)
        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        ! Initialise output
        CALC_GRAD_PARAMS(:) = 0.D0

        IF (BTEST(OPP_PARAMS, 0)) THEN
            CALC_GRAD_PARAMS(1) = CALC_DVDK(MODRIJ) * CALC_XPLOR(MODRIJ)
        END IF
        IF (BTEST(OPP_PARAMS, 1)) THEN
            CALC_GRAD_PARAMS(2) = CALC_DVDPHI(MODRIJ) * CALC_XPLOR(MODRIJ)
        END IF

    END FUNCTION CALC_GRAD_PARAMS

!-------------------------------------------------------------------------------
!
! Calculates the second derivative contributions for the potential parameters
! for a cluster
! RIJ: vector between particles
! IDI, IDJ: indices of the atom pair
!
!-------------------------------------------------------------------------------

    SUBROUTINE CALC_HESS_PARAMS(RIJ, IDI, IDJ)

        ! Number of atoms, including unit cell parameters
        USE COMMONS, ONLY: NATOMS

        ! HESS: the Hessian matrix
        USE MODHESS, ONLY: HESS

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: RIJ(3)
        INTEGER, INTENT(IN) :: IDI, IDJ

        ! MODRIJ: particle-particle distance
        ! TEMP_GRAD: temporary calculation store
        DOUBLE PRECISION :: MODRIJ, TEMP_GRAD(3)

        ! Calculate some factors
        MODRIJ = SQRT(DOT_PRODUCT(RIJ(:), RIJ(:)))

        IF (BTEST(OPP_PARAMS, 0)) THEN
            ! OPP_K pure derivative
            HESS(3*NATOMS-2,3*NATOMS-2) = HESS(3*NATOMS-2,3*NATOMS-2) + &
                CALC_D2VDK2(MODRIJ) * CALC_XPLOR(MODRIJ)

            ! OPP_K and position mixed derivative
            TEMP_GRAD(:) = (CALC_D2VDRDK(MODRIJ) * CALC_XPLOR(MODRIJ) + &
                CALC_DVDK(MODRIJ) * CALC_DXPLORDR(MODRIJ)) * RIJ(:) / MODRIJ
            HESS(3*IDI-2:3*IDI,3*NATOMS-2) = &
                HESS(3*IDI-2:3*IDI,3*NATOMS-2) + TEMP_GRAD(:)
            HESS(3*NATOMS-2,3*IDI-2:3*IDI) = &
                HESS(3*NATOMS-2,3*IDI-2:3*IDI) + TEMP_GRAD(:)
            HESS(3*IDJ-2:3*IDJ,3*NATOMS-2) = &
                HESS(3*IDJ-2:3*IDJ,3*NATOMS-2) - TEMP_GRAD(:)
            HESS(3*NATOMS-2,3*IDJ-2:3*IDJ) = &
                HESS(3*NATOMS-2,3*IDJ-2:3*IDJ) - TEMP_GRAD(:)
        END IF        

        IF (BTEST(OPP_PARAMS, 1)) THEN
            ! OPP_PHI pure derivative
            HESS(3*NATOMS-1,3*NATOMS-1) = HESS(3*NATOMS-1,3*NATOMS-1) + &
                CALC_D2VDPHI2(MODRIJ) * CALC_XPLOR(MODRIJ)

            ! OPP_PHI and position mixed derivative
            TEMP_GRAD(:) = (CALC_D2VDRDPHI(MODRIJ) * CALC_XPLOR(MODRIJ) + &
                CALC_DVDPHI(MODRIJ) * CALC_DXPLORDR(MODRIJ)) * RIJ(:) / MODRIJ
            HESS(3*IDI-2:3*IDI,3*NATOMS-1) = &
                HESS(3*IDI-2:3*IDI,3*NATOMS-1) + TEMP_GRAD(:)
            HESS(3*NATOMS-1,3*IDI-2:3*IDI) = &
                HESS(3*NATOMS-1,3*IDI-2:3*IDI) + TEMP_GRAD(:)
            HESS(3*IDJ-2:3*IDJ,3*NATOMS-1) = &
                HESS(3*IDJ-2:3*IDJ,3*NATOMS-1) - TEMP_GRAD(:)
            HESS(3*NATOMS-1,3*IDJ-2:3*IDJ) = &
                HESS(3*NATOMS-1,3*IDJ-2:3*IDJ) - TEMP_GRAD(:)
        END IF

        IF (BTEST(OPP_PARAMS, 0) .AND. BTEST(OPP_PARAMS, 1)) THEN
            ! OPP_K and OPP_PHI derivative
            HESS(3*NATOMS-2,3*NATOMS-1) = HESS(3*NATOMS-2,3*NATOMS-1) + &
                CALC_D2VDKDPHI(MODRIJ) * CALC_XPLOR(MODRIJ)
            HESS(3*NATOMS-1,3*NATOMS-2) = HESS(3*NATOMS-1,3*NATOMS-2) + &
                CALC_D2VDKDPHI(MODRIJ) * CALC_XPLOR(MODRIJ)
        END IF

    END SUBROUTINE CALC_HESS_PARAMS

!-------------------------------------------------------------------------------
!
! Calculates the second derivative contributions for the potential parameters in
! a periodic system.
! RIJ: vector between particles in fractional coordinates
! BD: box derivative information
! IDI, IDJ: indices of the atom pair
!
!-------------------------------------------------------------------------------

    SUBROUTINE CALC_HESS_PARAMS_PER(RIJ, BD, IDI, IDJ)

        ! Number of atoms, including unit cell parameters
        USE COMMONS, ONLY: NATOMS

        ! HESS: the Hessian matrix
        USE MODHESS, ONLY: HESS

        ! BOX_DERIV: type for storing box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: RIJ(3)
        TYPE(BOX_DERIV), INTENT(IN) :: BD
        INTEGER, INTENT(IN) :: IDI, IDJ

        ! YIJ: particle-particle vector in Cartesian space
        ! MODYIJ: particle-particle distance in Cartesian space
        ! TEMP, TEMP_GRAD: temporary calculation stores
        DOUBLE PRECISION :: YIJ(3), MODYIJ, TEMP, TEMP_GRAD(3)

        ! I: loop index
        INTEGER :: I

        ! Calculate some factors
        YIJ(:) = MATMUL(BD%ORTHOG(:,:), RIJ(:))
        MODYIJ = SQRT(DOT_PRODUCT(YIJ(:), YIJ(:)))

        IF (BTEST(OPP_PARAMS, 0)) THEN
            ! OPP_K pure derivative
            HESS(3*NATOMS-8,3*NATOMS-8) = HESS(3*NATOMS-8,3*NATOMS-8) + &
                CALC_D2VDK2(MODYIJ) * CALC_XPLOR(MODYIJ) * &
                MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ)

            ! OPP_K and position mixed derivative
            TEMP_GRAD(:) = (CALC_D2VDRDK(MODYIJ) * CALC_XPLOR(MODYIJ) + &
                CALC_DVDK(MODYIJ) * CALC_DXPLORDR(MODYIJ)) * &
                MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                MATMUL(TRANSPOSE(BD%ORTHOG(:,:)), YIJ(:)) / MODYIJ
            HESS(3*IDI-2:3*IDI,3*NATOMS-8) = &
                HESS(3*IDI-2:3*IDI,3*NATOMS-8) + TEMP_GRAD(:)
            HESS(3*NATOMS-8,3*IDI-2:3*IDI) = &
                HESS(3*NATOMS-8,3*IDI-2:3*IDI) + TEMP_GRAD(:)
            HESS(3*IDJ-2:3*IDJ,3*NATOMS-8) = &
                HESS(3*IDJ-2:3*IDJ,3*NATOMS-8) - TEMP_GRAD(:)
            HESS(3*NATOMS-8,3*IDJ-2:3*IDJ) = &
                HESS(3*NATOMS-8,3*IDJ-2:3*IDJ) - TEMP_GRAD(:)

            ! OPP_K and unit cell parameter mixed derivative
            DO I = 1, 6 ! loop over unit cell parameters
                TEMP = (CALC_D2VDRDK(MODYIJ) * CALC_XPLOR(MODYIJ) + &
                    CALC_DVDK(MODYIJ) * CALC_DXPLORDR(MODYIJ)) * &
                    MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                    DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV(:,:,I), RIJ(:))) / &
                    MODYIJ
                HESS(3*NATOMS-6+I,3*NATOMS-8) = &
                    HESS(3*NATOMS-6+I,3*NATOMS-8) + TEMP
                HESS(3*NATOMS-8,3*NATOMS-6+I) = &
                    HESS(3*NATOMS-8,3*NATOMS-6+I) + TEMP
            END DO ! end loop over unit cell parameters
        END IF        

        IF (BTEST(OPP_PARAMS, 1)) THEN
            ! OPP_PHI pure derivative
            HESS(3*NATOMS-7,3*NATOMS-7) = HESS(3*NATOMS-7,3*NATOMS-7) + &
                CALC_D2VDPHI2(MODYIJ) * CALC_XPLOR(MODYIJ) * &
                MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ)

            ! OPP_PHI and position mixed derivative
            TEMP_GRAD(:) = (CALC_D2VDRDPHI(MODYIJ) * CALC_XPLOR(MODYIJ) + &
                CALC_DVDPHI(MODYIJ) * CALC_DXPLORDR(MODYIJ)) * &
                MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                MATMUL(TRANSPOSE(BD%ORTHOG(:,:)), YIJ(:)) / MODYIJ
            HESS(3*IDI-2:3*IDI,3*NATOMS-7) = &
                HESS(3*IDI-2:3*IDI,3*NATOMS-7) + TEMP_GRAD(:)
            HESS(3*NATOMS-7,3*IDI-2:3*IDI) = &
                HESS(3*NATOMS-7,3*IDI-2:3*IDI) + TEMP_GRAD(:)
            HESS(3*IDJ-2:3*IDJ,3*NATOMS-7) = &
                HESS(3*IDJ-2:3*IDJ,3*NATOMS-7) - TEMP_GRAD(:)
            HESS(3*NATOMS-7,3*IDJ-2:3*IDJ) = &
                HESS(3*NATOMS-7,3*IDJ-2:3*IDJ) - TEMP_GRAD(:)

            ! OPP_K and unit cell parameter mixed derivative
            DO I = 1, 6 ! loop over unit cell parameters
                TEMP = (CALC_D2VDRDPHI(MODYIJ) * CALC_XPLOR(MODYIJ) + &
                    CALC_DVDPHI(MODYIJ) * CALC_DXPLORDR(MODYIJ)) * &
                    MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ) * &
                    DOT_PRODUCT(YIJ(:), MATMUL(BD%DERIV(:,:,I), RIJ(:))) / &
                    MODYIJ
                HESS(3*NATOMS-6+I,3*NATOMS-7) = &
                    HESS(3*NATOMS-6+I,3*NATOMS-7) + TEMP
                HESS(3*NATOMS-7,3*NATOMS-6+I) = &
                    HESS(3*NATOMS-7,3*NATOMS-6+I) + TEMP
            END DO ! end loop over unit cell parameters
        END IF

        IF (BTEST(OPP_PARAMS, 0) .AND. BTEST(OPP_PARAMS, 1)) THEN
            ! OPP_K and OPP_PHI derivative
            TEMP = CALC_D2VDKDPHI(MODYIJ) * CALC_XPLOR(MODYIJ) * &
                MERGE(0.5D0, 1.0D0, IDI .EQ. IDJ)
            HESS(3*NATOMS-8,3*NATOMS-7) = HESS(3*NATOMS-8,3*NATOMS-7) + TEMP
            HESS(3*NATOMS-7,3*NATOMS-8) = HESS(3*NATOMS-7,3*NATOMS-8) + TEMP
        END IF

    END SUBROUTINE CALC_HESS_PARAMS_PER

!-------------------------------------------------------------------------------
!
! Derivative of potential wrt OPP_K without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_DVDK (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_DVDK = - SIN(OPP_K * (MODRIJ-1) - OPP_PHI) * (MODRIJ-1) / MODRIJ**3

    END FUNCTION CALC_DVDK

!-------------------------------------------------------------------------------
!
! Derivative of potential wrt OPP_PHI without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_DVDPHI (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_DVDPHI = SIN(OPP_K * (MODRIJ - 1) - OPP_PHI) / MODRIJ**3

    END FUNCTION CALC_DVDPHI

!-------------------------------------------------------------------------------
!
! Second derivative of potential wrt OPP_K twice without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2VDK2 (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_D2VDK2 = &
            - COS(OPP_K * (MODRIJ - 1) - OPP_PHI) * (MODRIJ - 1)**2 / MODRIJ**3

    END FUNCTION CALC_D2VDK2

!-------------------------------------------------------------------------------
!
! Second derivative of potential wrt OPP_PHI twice without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2VDPHI2 (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_D2VDPHI2 = - COS(OPP_K * (MODRIJ - 1) - OPP_PHI) / MODRIJ**3

    END FUNCTION CALC_D2VDPHI2

!-------------------------------------------------------------------------------
!
! Second derivative of potential wrt OPP_K and OPP_PHI without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2VDKDPHI (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_D2VDKDPHI = &
            COS(OPP_K * (MODRIJ - 1) - OPP_PHI) * (MODRIJ - 1) / MODRIJ**3

    END FUNCTION CALC_D2VDKDPHI

!-------------------------------------------------------------------------------
!
! Second derivative of potential wrt position and OPP_K without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2VDRDK (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_D2VDRDK = &
            (- COS(OPP_K * (MODRIJ - 1) - OPP_PHI) * OPP_K * (MODRIJ - 1) - &
            SIN(OPP_K * (MODRIJ - 1) - OPP_PHI) + &
            SIN(OPP_K * (MODRIJ - 1) - OPP_PHI) * 3.D0 * (MODRIJ - 1) / &
            MODRIJ) / MODRIJ**3

    END FUNCTION CALC_D2VDRDK

!-------------------------------------------------------------------------------
!
! Second derivative of potential wrt position and OPP_PHI without XPLOR
! MODRIJ: distance between particles
!
!-------------------------------------------------------------------------------
    PURE DOUBLE PRECISION FUNCTION CALC_D2VDRDPHI (MODRIJ)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: MODRIJ

        CALC_D2VDRDPHI = (COS(OPP_K * (MODRIJ - 1) - OPP_PHI) * OPP_K - &
            SIN(OPP_K * (MODRIJ - 1) - OPP_PHI) * 3.D0 / MODRIJ) / MODRIJ**3

    END FUNCTION CALC_D2VDRDPHI

!-------------------------------------------------------------------------------
!
! Repel the unit cell volume from approaching 0, which repels the unit cell from
! adopting physically impossible angle combinations. Uses a WCA potential
! ENERGY: energy of the system
! GRAD_ANGLES: gradients for the box angles
! BD: information on the box and its derivatives
! GTEST: whether to calculate gradients
! STEST: whether to calculate second derivatives
!
!-------------------------------------------------------------------------------

    SUBROUTINE CONSTRAIN_VOLUME(ENERGY, GRAD_ANGLES, BD, GTEST, STEST)

        ! NATOMS: number of atoms, including unit cell parameters
        USE COMMONS, ONLY: NATOMS

        ! HESS: the Hessian matrix
        USE MODHESS, ONLY: HESS

        ! BOX_DERIV: type for storing box derivative information
        USE TRICLINIC_MOD, ONLY: BOX_DERIV

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: ENERGY, GRAD_ANGLES(3)
        TYPE(BOX_DERIV), INTENT(IN) :: BD
        LOGICAL, INTENT(IN) :: GTEST, STEST

        ! DEDV: derivative of energy wrt volume
        ! D2EDV2: second derivative of energy wrt volume
        DOUBLE PRECISION :: DEDV, D2EDV2

        ! I: loop index
        INTEGER :: I

        ! Add a purely repulsive (away from zero) WCA energy term
        IF (BD%V .LT. V_SIGMA * 2.D0**(1.D0/6.D0)) THEN
            ENERGY = ENERGY + 4.D0 * V_EPS * &
                     ((V_SIGMA / BD%V)**(12) - (V_SIGMA / BD%V)**(6)) + V_EPS

            IF (GTEST) THEN
                ! Add the gradients, if necessary
                DEDV = 24.D0 * V_EPS / V_SIGMA * &
                    ((V_SIGMA / BD%V)**(7) - 2.D0 * (V_SIGMA / BD%V)**(13))
                GRAD_ANGLES(:) = GRAD_ANGLES(:) +  DEDV * BD%V_DERIV(:)
            END IF

            IF (STEST) THEN
                ! Add the Hessian contributions, if necessary
                DEDV = 24.D0 * V_EPS / V_SIGMA * &
                    ((V_SIGMA / BD%V)**(7) - 2.D0 * (V_SIGMA / BD%V)**(13))
                D2EDV2 = 24.D0 * V_EPS / V_SIGMA**2 * &
                    (26.D0 * (V_SIGMA/BD%V)**(14) - 7.D0 * (V_SIGMA/BD%V)**(8))

                DO I = 1, 3 ! loop over angles
                    HESS(3*NATOMS-5+I,3*NATOMS-5:3*NATOMS-3) = &
                        HESS(3*NATOMS-5+I,3*NATOMS-5:3*NATOMS-3) + &
                        D2EDV2 * BD%V_DERIV(I) * BD%V_DERIV(:) + &
                        DEDV * BD%V_DERIV2(I, :)
                END DO ! end loop over angles
            END IF
        END IF

    END SUBROUTINE CONSTRAIN_VOLUME

!-------------------------------------------------------------------------------
!
! Impose limits on the allowed range of the potential parameters OPP_K by adding
! a harmonic repulsion outside the allowed range.
! ENERGY: energy of the system
! GRAD_PARAMS: gradients for the potential parameters
! GTEST: whether to calculate gradients
! STEST: whether to calculate second derivatives
!
!-------------------------------------------------------------------------------
    SUBROUTINE CONSTRAIN_PARAMETERS(ENERGY, GRAD_PARAMS, GTEST, STEST)

        ! NATOMS: number of atoms, including unit cell and potential parameters
        USE COMMONS, ONLY: NATOMS

        ! HESS: the Hessian matrix
        USE MODHESS, ONLY: HESS        

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: ENERGY, GRAD_PARAMS(3)
        LOGICAL, INTENT(IN) :: GTEST, STEST

        ! EN_CON: constraint contribution to the energy
        DOUBLE PRECISION :: EN_CON

        ! OFFSET: position of the end of the atom indices of the potential
        !         paramters
        INTEGER :: OFFSET

        ! Initialise output variable
        EN_CON = 0.D0

        IF (OPP_MODE .EQ. 2) THEN
            OFFSET = - 6
        ELSE IF (OPP_MODE .EQ. 3) THEN
            OFFSET = 0
        END IF

        ! OPP_K energy contribution
        IF (BTEST(OPP_PARAMS, 0)) THEN
            EN_CON = EN_CON + 0.5D0 * K_REPULSION * &
                (MERGE(OPP_K - OPP_K_UPPER, 0.D0, &
                       OPP_K .GT. OPP_K_UPPER)**2 + &
                 MERGE(OPP_K - OPP_K_LOWER, 0.D0, &
                       OPP_K .LT. OPP_K_LOWER)**2)

            IF (GTEST) THEN
                ! OPP_K gradient contribution
                GRAD_PARAMS(1) = GRAD_PARAMS(1) + K_REPULSION * &
                    (MERGE(OPP_K - OPP_K_UPPER, 0.D0, &
                           OPP_K .GT. OPP_K_UPPER) + &
                     MERGE(OPP_K - OPP_K_LOWER, 0.D0, &
                           OPP_K .LT. OPP_K_LOWER))
            END IF

            IF (STEST) THEN
                ! OPP_K Hessian contribution
                HESS(3*NATOMS+OFFSET-2,3*NATOMS+OFFSET-2) = &
                    HESS(3*NATOMS+OFFSET-2,3*NATOMS+OFFSET-2) + &
                    MERGE(K_REPULSION, 0.D0, OPP_K .GT. OPP_K_UPPER) + &
                    MERGE(K_REPULSION, 0.D0, OPP_K .LT. OPP_K_LOWER)
            END IF 
        END IF

        ! Add the energy contribution
        ENERGY = ENERGY + EN_CON

    END SUBROUTINE CONSTRAIN_PARAMETERS

!-------------------------------------------------------------------------------
!
! Outputs the coordinates to the file opp.xyz
!
!-------------------------------------------------------------------------------
    SUBROUTINE VIEW_OPP()

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

        ! GETUNIT: function for finding a free file unit
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
        SELECT CASE (OPP_MODE)
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
            NMOL = NATOMS - 2
        CASE (2)
            NMOL = NATOMS - 3
        END SELECT

        ! Open the file
        FUNIT = GETUNIT()
        OPEN(UNIT=FUNIT, FILE='opp.xyz', STATUS='UNKNOWN')

        ! Loop over the configurations
        DO I = 1, NSAVE

            WRITE(FUNIT,'(I6)') NMOL

            WRITE(FUNIT, '(A,I6,A,F20.10,A,I8,A,F20.10)') 'Energy of minimum ',&
                I, '=', QMIN(I), ' first found at step ', FF(I), &
                ' Energy per particle = ', QMIN(I)/NMOL

            ! Set box parameters for varying box
            IF (OPP_MODE .NE. 0) THEN
                BOX_SIZE(:) = QMINP(I, 3*NATOMS-2:3*NATOMS)
            END IF
            IF (OPP_MODE .EQ. 1 .OR. OPP_MODE .EQ. 2) THEN
                BOX_ANGLES(:) = QMINP(I, 3*NATOMS-5:3*NATOMS-3)
                BD = CALC_BOX_DERIV(BOX_ANGLES(:), BOX_SIZE(:))
            END IF

            DO J = 1, NMOL

                ! If necessary, transform the coordinate from crytallographic
                ! space to orthogonal space
                SELECT CASE (OPP_MODE)
                CASE (0)
                    COORD(:) = QMINP(I, 3*J-2:3*J)
                CASE (1, 2)
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

    END SUBROUTINE VIEW_OPP

END MODULE OPP_MOD
