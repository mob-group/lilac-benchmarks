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
! This module provides a routine to calculatye transformation matrices and 
! derivatives for a triclinic unit cell, as well as a type to store them.
! Questions to jwrm2
!
MODULE TRICLINIC_MOD

PUBLIC :: BOX_DERIV, CALC_BOX_DERIV

! BOX_DERIV: stores information for calculating lattice derivatives
TYPE :: BOX_DERIV
    ! ANGLES: unit cell angles
    ! SIZES: unit cell lengths
    DOUBLE PRECISION :: ANGLES(3), SIZES(3)

    ! ORTHOG: the orthogonal tansformation matrix, for changing from
    !         crystallographic space to orthogonal space first index is row,
    !         second is column
    ! ORTHOGSQ: MATMUL(ORTHOG, ORTHOG), needed for second derivatives
    DOUBLE PRECISION :: ORTHOG(3,3), ORTHOGSQ(3, 3)
    ! ORTHOG: the inverse orthogonal tansformation matrix, for changing from
    !             orthogonal space to crystallographic space
    DOUBLE PRECISION :: ORTHOG_INV(3,3)

    ! DERIV: derivatives of the orthogonal matrix wrt lattice parameters,
    !        first index indicates the row, second index the column and
    !        third index the parameter, with 1-3 being the angles and 4-6 being
    !        the lengths
    ! DERIV2: derivatives of DERIV wrt the six unit cell parameters
    DOUBLE PRECISION :: DERIV(3, 3, 6), DERIV2(3, 3, 6, 6)

    ! DERIV_INV: derivatives of the inverse orthogonal matrix wrt lattice 
    !        parameters, first index indicates the row, second index the column
    !        and third index the parameter, with 1-3 being the angles and 4-6 
    !        being the lengths
    DOUBLE PRECISION :: DERIV_INV(3, 3, 6)

    ! RECIP: lengths of the reciprocal lattice vectors
    DOUBLE PRECISION :: RECIP(3)

    ! CELL_RANGE: number of cells in each direction that must be evaluated
    INTEGER :: CELL_RANGE(3)

    ! V: dimensionless volume of the unit cell, ie volume if all lengths were 1
    DOUBLE PRECISION :: V

    ! V_DERIV: derivatives of V wrt the lattice angles
    ! V_DERIV2: sedond derivatives of V wrt the lattice angles
    DOUBLE PRECISION :: V_DERIV(3), V_DERIV2(3, 3)
END TYPE BOX_DERIV

CONTAINS

!-------------------------------------------------------------------------------
!
! Calculates factors for the lattice derivatives that only depend on the cell
! parameters.
! ANGLES: alpha, beta and gamma lattice parameters
! SIZES: a, b and c lattice parameters
! CUTOFF: (optional) provide for working out the number of periodic repeats in
!         each direction within the cutoff sphere
!
!-------------------------------------------------------------------------------
    PURE TYPE(BOX_DERIV) FUNCTION CALC_BOX_DERIV(ANGLES, SIZES, CUTOFF)

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: ANGLES(3), SIZES(3)
        DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: CUTOFF

        ! V: factor, related to (but not quite equal to) the unit cell volume
        ! DVDA: derivatives of V wrt the three lattice angles
        ! C, S: cos and sin of the lattice angles
        ! AL, BL, CL: unit cell lengths
        ! ORTHOG, INV_ORTHOG, DERIV, DERIV2, DERIV_INV: parts of the output
        DOUBLE PRECISION :: V, DVDA(3), C(3), S(3), AL, BL, CL
        DOUBLE PRECISION :: ORTHOG(3, 3), ORTHOG_INV(3, 3)
        DOUBLE PRECISION :: DERIV(3, 3, 6), DERIV2(3,3,6,6), DERIV_INV(3, 3, 6)

        ! Initialise output variables
        ORTHOG(:,:) = 0.D0
        DERIV(:,:,:) = 0.D0
        DERIV2(:,:,:,:) = 0.D0
        ORTHOG_INV(:,:) = 0.D0
        DERIV_INV(:,:,:) = 0.D0

        ! Calculate factors
        AL = SIZES(1)
        BL = SIZES(2)
        CL = SIZES(3)
        C(:) = COS(ANGLES(:))
        S(:) = SIN(ANGLES(:))
        V = SQRT(1 - C(1)**2 - C(2)**2 - C(3)**2 + 2.D0 * C(1) * C(2) * C(3))
        DVDA(1) = S(1) * (C(1) - C(2) * C(3)) / V
        DVDA(2) = S(2) * (C(2) - C(3) * C(1)) / V
        DVDA(3) = S(3) * (C(3) - C(1) * C(2)) / V

        ! Calculate the orthogonal transformation matrix
        ORTHOG(1, 1) = AL
        ORTHOG(1, 2) = BL * C(3)
        ORTHOG(1, 3) = CL * C(2)
        ORTHOG(2, 2) = BL * S(3)
        ORTHOG(2, 3) = CL * (C(1) - C(2) * C(3)) / S(3)
        ORTHOG(3, 3) = CL * V / S(3)

        ! Calculate the derivatives of the othogonal matrix
        ! Derivatives of the top row wrt angles
        DERIV(1, 3, 2) = - CL * S(2)
        DERIV(1, 2, 3) = - BL * S(3)
        ! Derivatives of the top row wrt to lengths
        DERIV(1, 1, 4) = 1.D0
        DERIV(1, 2, 5) = C(3)
        DERIV(1, 3, 6) = C(2)
        ! Derivatives of the middle row wrt angles
        DERIV(2, 3, 1) = - CL * S(1) / S(3)
        DERIV(2, 3, 2) = CL * S(2) * C(3) / S(3)
        DERIV(2, 2, 3) = BL * C(3)
        DERIV(2, 3, 3) = CL * (C(2) - C(1) * C(3)) / S(3)**2
        ! Derivatives of the middle row wrt to lengths
        DERIV(2, 2, 5) = S(3)
        DERIV(2, 3, 6) = (C(1) - C(2) * C(3)) / S(3)
        ! Derivatives of the bottom row wrt angles
        DERIV(3, 3, 1) = CL * DVDA(1) / S(3)
        DERIV(3, 3, 2) = CL * DVDA(2) / S(3)
        DERIV(3, 3, 3) = CL * (DVDA(3) - C(3) * V / S(3)) / S(3)
        ! Derivatives of the bottom row wrt lengths
        DERIV(3, 3, 6) = V / S(3)

        ! Calculate the second derivatives of the orthogonal matrix. FML.
        ! Derivatives of element 12
        DERIV2(1, 2, 3, 3) = - BL * C(3)
        DERIV2(1, 2, 3, 5) = - S(3)
        DERIV2(1, 2, 5, 3) = DERIV2(1, 2, 3, 5)
        ! Derivatives of element 13
        DERIV2(1, 3, 2, 2) = - CL * C(2)
        DERIV2(1, 3, 2, 6) = - S(2)
        DERIV2(1, 3, 6, 2) = DERIV2(1, 3, 2, 6)
        ! Derivatives of element 22
        DERIV2(2, 2, 3, 3) = - BL * S(3)
        DERIV2(2, 2, 3, 5) = C(3)
        DERIV2(2, 2, 5, 3) = DERIV2(2, 2, 3, 5)
        ! Derivatives of element 23
        DERIV2(2, 3, 1, 1) = - CL * C(1) / S(3)
        DERIV2(2, 3, 1, 3) = CL * S(1) * C(3) / S(3)**2
        DERIV2(2, 3, 3, 1) = DERIV2(2, 3, 1, 3)
        DERIV2(2, 3, 1, 6) = - S(1) / S(3)
        DERIV2(2, 3, 6, 1) = DERIV2(2, 3, 1, 6)
        DERIV2(2, 3, 2, 2) = CL * C(2) * C(3) / S(3)
        DERIV2(2, 3, 2, 3) = - CL * S(2) / S(3)**2
        DERIV2(2, 3, 3, 2) = DERIV2(2, 3, 2, 3)
        DERIV2(2, 3, 2, 6) = S(2) * C(3) / S(3)
        DERIV2(2, 3, 6, 2) = DERIV2(2, 3, 2, 6)
        DERIV2(2, 3, 3, 3) = CL*(2.D0 * C(3) * (C(1) * C(3) - C(2)) / S(3)**3 +&
            C(1) / S(3))
        DERIV2(2, 3, 3, 6) = (C(2) - C(1) * C(3)) / S(3)**2
        DERIV2(2, 3, 6, 3) = DERIV2(2, 3, 3, 6)
        ! Derivatives of element 33
        DERIV2(3, 3, 1, 1) = CL * ((C(1)**2 - S(1)**2 - C(1) * C(2) * C(3)) - &
            DVDA(1)**2) / (V * S(3))
        DERIV2(3, 3, 1, 2) = CL * (S(1) * S(2) * C(3) - DVDA(1) * DVDA(2)) / &
            (V * S(3))
        DERIV2(3, 3, 2, 1) = DERIV2(3, 3, 1, 2)
        DERIV2(3, 3, 1, 3) = CL * (S(1) * C(2) / V - &
            DVDA(1) * DVDA(3) / (S(3) * V) - DVDA(1) * C(3) / S(3)**2)
        DERIV2(3, 3, 3, 1) = DERIV2(3, 3, 1, 3)
        DERIV2(3, 3, 1, 6) = DVDA(1) / S(3)
        DERIV2(3, 3, 6, 1) = DERIV2(3, 3, 1, 6)
        DERIV2(3, 3, 2, 2) = CL * ((C(2)**2 - S(2)**2 - C(1) * C(2) * C(3)) - &
            DVDA(2)**2) / (V * S(3))
        DERIV2(3, 3, 2, 3) = CL * (S(2) * C(1) / V - &
            DVDA(2) * DVDA(3) / (S(3) * V) - DVDA(2) * C(3) / S(3)**2)
        DERIV2(3, 3, 3, 2) = DERIV2(3, 3, 2, 3)
        DERIV2(3, 3, 2, 6) = DVDA(2) / S(3)
        DERIV2(3, 3, 6, 2) = DERIV2(3, 3, 2, 6)
        DERIV2(3, 3, 3, 3) = CL * ((C(3)**2 - S(3)**2 - C(1) * C(2) * C(3))/V -&
            DVDA(3)**2 / V - C(3) * DVDA(3) / S(3) + V / S(3)**2)/ S(3) - &
            CL * C(3) * (DVDA(3) - V * C(3) / S(3)) / S(3)**2
        DERIV2(3, 3, 3, 6) = (DVDA(3) - V * C(3) / S(3)) / S(3)
        DERIV2(3, 3, 6, 3) = DERIV2(3, 3, 3, 6)

        ! The inverse of the orthogonal transformation matrix
        ORTHOG_INV(1, 1) = 1.D0 / AL
        ORTHOG_INV(1, 2) = - C(3) / ( AL * S(3))
        ORTHOG_INV(1, 3) = (C(1) * C(3) - C(2)) / (AL * S(3) * V)
        ORTHOG_INV(2, 2) = 1.D0 / (BL * S(3))
        ORTHOG_INV(2, 3) = (C(2) * C(3) - C(1)) / (BL * S(3) * V)
        ORTHOG_INV(3, 3) = S(3) / (CL * V)

        ! Calculate the derivatives of the inverse othogonal matrix
        ! Derivatives of the top row wrt angles
        DERIV_INV(1, 3, 1) = - S(1) * C(3) / (AL * S(3) * V) + &
            (C(1) * C(3) - C(2)) * (C(2) * C(3) - C(1))*S(1)/(AL * S(3) * V**3)
        DERIV_INV(1, 3, 2) = S(2) / (AL * S(3) * V) + &
            (C(1) * C(3) - C(2))**2 * S(2) / (AL * S(3) * V**3)
        DERIV_INV(1, 2, 3) = 1.D0 / (AL * S(3)**2)
        DERIV_INV(1, 3, 3) = - C(1) / (AL * V) + &
            (C(1) * C(3) - C(2)) * ((C(1) * C(2) - C(3)) / &
            (AL * V**3) - C(3) / (AL * S(3)**2 * V))
        ! Derivatives of the top row wrt to lengths
        DERIV_INV(1, 1, 4) = - 1.D0 / AL**2
        DERIV_INV(1, 2, 4) = C(3) / (AL**2 * S(3))
        DERIV_INV(1, 3, 4) = (C(2) - C(1) * C(3)) / (AL**2 * S(3) * V)
        ! Derivatives of the middle row wrt angles
        DERIV_INV(2, 3, 1) = S(1) / (BL * S(3) * V) + &
            (C(2) * C(3) - C(1))**2 * S(1) / (BL * S(3) * V**3)
        DERIV_INV(2, 3, 2) = - S(2) * C(3) / (BL * S(3) * V) + &
            (C(2) * C(3) - C(1)) * (C(1) * C(3) - C(2))*S(2)/(BL * S(3) * V**3)
        DERIV_INV(2, 2, 3) = - C(3) / (BL * S(3)**2)
        DERIV_INV(2, 3, 3) = - C(2) / (BL * V) + &
            (C(2) * C(3) - C(1)) * ((C(1) * C(2) - C(3)) / &
            (BL * V**3) - C(3) / (BL * S(3)**2 * V))
        ! Derivatives of the middle row wrt to lengths
        DERIV_INV(2, 2, 5) = - 1.D0 / (BL**2 * S(3))
        DERIV_INV(2, 3, 5) = (C(1) - C(2) * C(3)) / (BL**2 * S(3) * V)
        ! Derivatives of the bottom row wrt angles
        DERIV_INV(3, 3, 1) = S(1) * S(3) * (C(2) * C(3) - C(1)) / (CL * V**3)
        DERIV_INV(3, 3, 2) = S(2) * S(3) * (C(1) * C(3) - C(2)) / (CL * V**3)
        DERIV_INV(3, 3, 3) = C(3) / (CL * V) + &
            S(3)**2 * (C(1) * C(2) - C(3)) / (CL * V**3)
        ! Derivatives of the bottom row wrt lengths
        DERIV_INV(3, 3, 6) = - S(3) / (CL**2 * V)

        ! Calculate the derivative of volume wrt the angles
        CALC_BOX_DERIV%V_DERIV(1) = S(1) * (C(1) - C(2) * C(3)) / V
        CALC_BOX_DERIV%V_DERIV(2) = S(2) * (C(2) - C(3) * C(1)) / V
        CALC_BOX_DERIV%V_DERIV(3) = S(3) * (C(3) - C(1) * C(2)) / V

        ! Calculate second derivatives of volume wrt the angles
        CALC_BOX_DERIV%V_DERIV2(1,1) = (S(1) * (C(1) - C(2) * C(3))**2 - &
            C(1)**2 + S(1)**2 - C(1) * C(2) * C(3)) / V**3
        CALC_BOX_DERIV%V_DERIV2(2,2) = (S(2) * (C(2) - C(3) * C(1))**2 - &
            C(2)**2 + S(2)**2 - C(2) * C(3) * C(1)) / V**3
        CALC_BOX_DERIV%V_DERIV2(3,3) = (S(3) * (C(3) - C(1) * C(2))**2 - &
            C(3)**2 + S(3)**2 - C(3) * C(1) * C(2)) / V**3
        CALC_BOX_DERIV%V_DERIV2(1,2) = (C(1) * C(2) - &
            C(3)) * S(1) * S(2) * S(3)**2 / V**3
        CALC_BOX_DERIV%V_DERIV2(2,1) = CALC_BOX_DERIV%V_DERIV2(1,2)
        CALC_BOX_DERIV%V_DERIV2(2,3) = (C(2) * C(3) - &
            C(1)) * S(2) * S(3) * S(1)**2 / V**3
        CALC_BOX_DERIV%V_DERIV2(3,2) = CALC_BOX_DERIV%V_DERIV2(2,3)
        CALC_BOX_DERIV%V_DERIV2(3,1) = (C(3) * C(1) - &
            C(2)) * S(3) * S(1) * S(2)**2 / V**3
        CALC_BOX_DERIV%V_DERIV2(1,3) = CALC_BOX_DERIV%V_DERIV2(3,1)

        ! Calculate the lengths of the reciprocal lattice vectors
        CALC_BOX_DERIV%RECIP(:) = S(:) / (SIZES(:) * V)

        ! Calculate the number of cells to check in each direction
        IF (PRESENT(CUTOFF)) THEN
            CALC_BOX_DERIV%CELL_RANGE(:) = CEILING(CUTOFF * &
                CALC_BOX_DERIV%RECIP(:))
        ELSE
            CALC_BOX_DERIV%CELL_RANGE(:) = 0.D0
        END IF

        ! Set the output
        CALC_BOX_DERIV%ANGLES(:) = ANGLES(:)
        CALC_BOX_DERIV%SIZES(:) = SIZES(:)
        CALC_BOX_DERIV%ORTHOG(:,:) = ORTHOG(:,:)
        CALC_BOX_DERIV%ORTHOGSQ(:,:) = &
            MATMUL(TRANSPOSE(ORTHOG(:,:)), ORTHOG(:,:))
        CALC_BOX_DERIV%DERIV(:,:,:) = DERIV(:,:,:)
        CALC_BOX_DERIV%DERIV2(:,:,:,:) = DERIV2(:,:,:,:)
        CALC_BOX_DERIV%ORTHOG_INV(:,:) = ORTHOG_INV(:,:)
        CALC_BOX_DERIV%DERIV_INV(:,:,:) = DERIV_INV(:,:,:)
        CALC_BOX_DERIV%V = V

    END FUNCTION CALC_BOX_DERIV

END MODULE TRICLINIC_MOD
