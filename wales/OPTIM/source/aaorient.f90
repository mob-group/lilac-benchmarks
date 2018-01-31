!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2013 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!----------------------------------------------------------------------------------------------!
!                                                                                              !
! AAORIENT                                                                                     !
! (Angle-Axis ORIENTation)                                                                     !
!                                                                                              !
! Calculates energy, gradient and hessian for the angle-axis orientation potential             !
! Questions to John Morgan (jwrm2)                                                             !
! See:                                                                                         !
! 'Simulations of rigid bodies in an angle axis framework', Dwaipayan Chakrabarti and          !
! David Wales, Phys. Chem. Chem. Phys., 2009, 11, 1970-1976                                    !
! for a description of angle-axis coordinates                                                  !
!                                                                                              !
! X: positions and orientations of the bodies in angle-axis coordinates                        !
! G: gradients with respect to each angle-axis coordinate                                      !
! ENERGY: the calculated energy                                                                !
! GTEST: logical, true if gradients are to be calculated                                       !
! STEST: logical, true if hessian is to be calculated                                          !
!                                                                                              !
!----------------------------------------------------------------------------------------------!

SUBROUTINE AAORIENT (X, G, ENERGY, GTEST, STEST)

! HESS: the hessian matrix
USE MODHESS, ONLY: HESS
! NATOMS: twice the number of rigid bodies
USE COMMONS, ONLY: NATOMS
! KAA: controls the strength of the potential
! SIGMAAA: controls the range of the potential
USE KEY, ONLY: KAA

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)    :: X(3*NATOMS)
DOUBLE PRECISION, INTENT(INOUT) :: G(3*NATOMS), ENERGY
LOGICAL, INTENT(IN)             :: GTEST, STEST

! REALATOMS: the actual number of rigid bodies
! OFFSET: the index of the last translational coordinate
INTEGER :: REALATOMS, OFFSET

! I, IT, IR, I2, J, JT, JR: Counters to iterate over, and indices of coordinates
INTEGER :: I, IT, IR, I2, J, JT, JR

! P: a set of rotational coordinates
! SX, SY, S2: vector sums of the orientations of the bodies
DOUBLE PRECISION, DIMENSION(3) :: P, SX, SY, SZ

! RMI: rotation matrix for a body
! DRMI: first derivative of rotation matrix wrt rotational coordinates (order is 1, 2, 3)
! D2RMI: second derivative of rotation matrix wrt rotational coordinates (order is 11, 22, 33, 12, 23, 31)
! PREFACTOR: appears in all expressions, equal to KAA/REALATOMS**2
DOUBLE PRECISION :: RMI(3,3), DRMI(3,3,3), D2RMI(3,3,6), PREFACTOR

! Dk: derivatives of rotation matrix acting on Cartesian unit vector k
! D2k: second derivatives of rotation matrix acting on Cartesian unit vector k (order is 11, 22, 33, 12, 23, 31)
! X0, Y0, Z0: lab frame cartesian unit vectors
DOUBLE PRECISION, DIMENSION(3, 3*NATOMS/2) :: DX, DY, DZ 
DOUBLE PRECISION, DIMENSION(3, 6*NATOMS/2) :: D2X, D2Y, D2Z
DOUBLE PRECISION, DIMENSION(3) :: X0, Y0, Z0

! Initialise Cartesian unit vectors
X0 = (/ 1.0, 0.0, 0.0 /)
Y0 = (/ 0.0, 1.0, 0.0 /)
Z0 = (/ 0.0, 0.0, 1.0 /)

! Initialise SX, SY, SZ, REALATOMS, OFFSET and PREFACTOR
SX(:) = 0.0D0
SY(:) = 0.0D0
SZ(:) = 0.0D0
REALATOMS = NATOMS/2
OFFSET    = 3*REALATOMS
PREFACTOR = KAA / REALATOMS**2

! I Loop over all bodies  
DO I = 1, REALATOMS

! Set IT as the index of the third translational coordinate of the Ith body, and IR as the index of the third rotational
! coordinate of the Ith body
! Set P to be the rotation vector of body I
  P(:) = X(3*I + OFFSET - 2:3*I + OFFSET)

! Calculate the rotation matrix and derivatives thereof
  CALL RMDFAS(P, RMI, DRMI(:,:,1), DRMI(:,:,2), DRMI(:,:,3), D2RMI(:,:,1), D2RMI(:,:,2), &
    &         D2RMI(:,:,3), D2RMI(:,:,4), D2RMI(:,:,5), D2RMI(:,:,6), GTEST, STEST)

! Now calculate the Sk values
  SX(:) = SX(:) + MATMUL(RMI(:,:), X0)
  SY(:) = SY(:) + MATMUL(RMI(:,:), Y0)
  SZ(:) = SZ(:) + MATMUL(RMI(:,:), Z0)

! Calculate first derivatives of the rotation matrix acting on the Cartesian unit vectors, if required
  IF (GTEST .OR. STEST) THEN
    DO J = 1, 3
      DX(:,3*I-3+J) = MATMUL(DRMI(:,:,J), X0)
      DY(:,3*I-3+J) = MATMUL(DRMI(:,:,J), Y0)
      DZ(:,3*I-3+J) = MATMUL(DRMI(:,:,J), Z0)
    END DO ! J loop over derivative matrixes
  END IF ! end if GTEST or STEST

! Calculate second derivatives of the rotation matrix acting on the Cartesian unit vectors, if required
  IF (STEST) THEN
    DO J = 1, 6
      D2X(:,6*I-6+J) = MATMUL(D2RMI(:,:,J), X0)
      D2Y(:,6*I-6+J) = MATMUL(D2RMI(:,:,J), Y0)
      D2Z(:,6*I-6+J) = MATMUL(D2RMI(:,:,J), Z0)
    END DO ! J loop over coordinate pairs
  END IF ! end if STEST
END DO ! I loop over each body

! Add the energy contribution
ENERGY = ENERGY - PREFACTOR * (DOT_PRODUCT(SX(:),SX(:)) + DOT_PRODUCT(SY(:),SY(:)) + DOT_PRODUCT(SZ(:),SZ(:)))

! Find the gradients (rotational only, translations are zero)
IF (GTEST) THEN
  DO I = 1, OFFSET
    IR = I + OFFSET
    G(IR) = G(IR) - 2.0D0 * PREFACTOR * (DOT_PRODUCT(SX(:),DX(:,I)) + DOT_PRODUCT(SY(:),DY(:,I)) + DOT_PRODUCT(SZ(:),DZ(:,I)))
  END DO ! I loop over rotational coordinates
END IF ! end if gtest

! Find the Hessian (pure rotational only, pure translation and mixed are zero)
IF (STEST) THEN
  DO IR = OFFSET+1, 2*OFFSET
    DO JR = OFFSET+1, 2*OFFSET
! Set I and J to the body index
      IT =  IR-OFFSET
      JT =  JR-OFFSET
      I = (IT+2)/3
      J = (JT+2)/3

      IF (I == J) THEN
! Same body
! Set I2 index for D2k
        IF (IR .EQ. JR) THEN 
          I2 = 6*I - 5 + MOD(IR-1,3)
        ELSE IF ((MOD(IR,3) .EQ. 1 .AND. MOD(JR,3) .EQ. 2) .OR. (MOD(IR,3) .EQ. 2 .AND. MOD(JR,3) .EQ. 1)) THEN 
          I2 = 6*I - 2
        ELSE IF ((MOD(IR,3) .EQ. 2 .AND. MOD(JR,3) .EQ. 0) .OR. (MOD(IR,3) .EQ. 0 .AND. MOD(JR,3) .EQ. 2)) THEN 
          I2 = 6*I - 1
        ELSE IF ((MOD(IR,3) .EQ. 0 .AND. MOD(JR,3) .EQ. 1) .OR. (MOD(IR,3) .EQ. 1 .AND. MOD(JR,3) .EQ. 0)) THEN 
          I2 = 6*I
        ENDIF
!        PRINT *, 'IR = ', IR, '; JR = ', JR, '; I = ', I, '; J = ', J, '; I2 = ', I2 
! Calculate the Hessian
        HESS(IR,JR) = HESS(IR,JR) - 2.D0 * PREFACTOR * (DOT_PRODUCT(SX(:),D2X(:,I2)) + DOT_PRODUCT(SY(:),D2Y(:,I2)) + &
     &              DOT_PRODUCT(SZ(:),D2Z(:,I2)) + DOT_PRODUCT(DX(:,IT),DX(:,JT)) + DOT_PRODUCT(DY(:,IT),DY(:,JT))  + &
     &              DOT_PRODUCT(DZ(:,IT),DZ(:,JT)))
      ELSE
! Different bodies
        HESS(IR,JR) = HESS(IR,JR) - 2.D0 * PREFACTOR * (DOT_PRODUCT(DX(:,IT),DX(:,JT)) + &
     &                DOT_PRODUCT(DY(:,IT),DY(:,JT)) + DOT_PRODUCT(DZ(:,IT),DZ(:,JT)))
      END IF ! end if same body
    END DO ! J loop over rotational coordinates
  END DO ! I loop over rotational coordinates
END IF ! end if STEST
END SUBROUTINE AAORIENT

!----------------------------------------------------------------------------------------------!
!                                                                                              !
! AAORIENTSR                                                                                   !
! (Angle-Axis ORIENTation Short-Range)                                                         !
!                                                                                              !
! Calculates energy, gradient and hessian for the angle-axis orientation potential, including  !
! a gaussian distance dependence so the range can be controlled.                               !
! Questions to John Morgan (jwrm2)                                                             !
! See:                                                                                         !
! 'Simulations of rigid bodies in an angle axis framework', Dwaipayan Chakrabarti and          !
! David Wales, Phys. Chem. Chem. Phys., 2009, 11, 1970-1976                                    !
! for a description of angle-axis coordinates                                                  !
!                                                                                              !
! X: positions and orientations of the bodies in angle-axis coordinates                        !
! G: gradients with respect to each angle-axis coordinate                                      !
! ENERGY: the calculated energy                                                                !
! GTEST: logical, true if gradients are to be calculated                                       !
! STEST: logical, true if hessian is to be calculated                                          !
!                                                                                              !
!----------------------------------------------------------------------------------------------!

SUBROUTINE AAORIENTSR (X, G, ENERGY, GTEST, STEST)

! HESS: the hessian matrix
USE MODHESS, ONLY: HESS
! NATOMS: twice the number of rigid bodies
USE COMMONS, ONLY: NATOMS
! KAA: controls the strength of the potential
! SIGMAAA: controls the range of the potential
USE KEY, ONLY: KAA, SIGMAAA

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)    :: X(3*NATOMS)
DOUBLE PRECISION, INTENT(INOUT) :: G(3*NATOMS), ENERGY
LOGICAL, INTENT(IN)             :: GTEST, STEST

! REALATOMS: the actual number of rigid bodies
! OFFSET: the index of the last translational coordinate
INTEGER :: REALATOMS, OFFSET

! I, IT, IR, I2, J, JT, JR, J2, K: Counters to iterate over, and indices of coordinates
INTEGER :: I, IT, IR, I2, J, JT, JR, J2, K

! P: a set of rotational coordinates
! RIJ: vector from body J to body I
! SX2, SY2, SZ2: squared sums of the orientations of the bodies, including range modulation
! SUMkJEXP, SUMkJEXPRIJ, SUMkJEXPDRIJ: sums needed for derivatives 
DOUBLE PRECISION :: P(3), SX2, SY2, SZ2
DOUBLE PRECISION, DIMENSION(3, NATOMS/2)   :: SUMXJEXP, SUMYJEXP, SUMZJEXP
DOUBLE PRECISION, DIMENSION(3, 3*NATOMS/2) :: SUMXJEXPRIJ, SUMYJEXPRIJ, SUMZJEXPRIJ
DOUBLE PRECISION, DIMENSION(3, 6*NATOMS/2) :: SUMXJEXPDRIJ, SUMYJEXPDRIJ, SUMZJEXPDRIJ

! Various short term stores for calculations
! EXPFACT: gaussian exponential term
! RIJ: vector from body J to body I
! LRIJ2: length of RIJ
! DSk2DR, DSk2DR: derivatives of Sk2 wrt translation or rotation
! SIGMA2: SIGMAAA squared
DOUBLE PRECISION :: EXPFACT, RIJ(3), LRIJ2, DSX2DR, DSY2DR, DSZ2DR, DSX2DP, DSY2DP, DSZ2DP, SIGMA2

! RMI: rotation matrix for a body
! DRMI: first derivative of rotation matrix wrt rotational coordinates (order is 1, 2, 3)
! D2RMI: second derivative of rotation matrix wrt rotational coordinates (order is 11, 22, 33, 12, 23, 31) 
DOUBLE PRECISION :: RMI(3,3), DRMI(3,3,3), D2RMI(3,3,6)

! XB, YB, ZB: rotation matrix acting on Cartesian unit vectors, for each body
! Dk: derivatives of rotation matrix acting on Cartesian unit vector k
! D2k: second derivatives of rotation matrix acting on Cartesian unit vector k (order is x2, y2, z2, xy, yz, zx)
! X0, Y0, Z0: lab frame cartesian unit vectors
DOUBLE PRECISION, DIMENSION(3, NATOMS/2) :: XB, YB, ZB
DOUBLE PRECISION, DIMENSION(3, 3*NATOMS/2) :: DX, DY, DZ 
DOUBLE PRECISION, DIMENSION(3, 6*NATOMS/2) :: D2X, D2Y, D2Z
DOUBLE PRECISION, DIMENSION(3) :: X0, Y0, Z0

! Initialise Cartesian unit vectors
X0 = (/ 1.0, 0.0, 0.0 /)
Y0 = (/ 0.0, 1.0, 0.0 /)
Z0 = (/ 0.0, 0.0, 1.0 /)

! Find REALATOMS as the actual number of bodies, OFFSET as the start of rotational coordinates and SIGMA2 as SIGMAAA squared
REALATOMS = NATOMS/2
OFFSET    = 3*REALATOMS
SIGMA2    = SIGMAAA**2

! I Loop over all bodies  
DO I = 1, REALATOMS

! Set IT as the index of the third translational coordinate of the Ith body, and IR as the index of the third rotational
! coordinate of the Ith body
! Set P to be the rotation vector of body I
  IT = 3*I
  IR = OFFSET + IT
  P(:) = X(IR-2:IR)

! Calculate the rotation matrix and derivatives thereof
  CALL RMDFAS(P, RMI, DRMI(:,:,1), DRMI(:,:,2), DRMI(:,:,3), D2RMI(:,:,1), D2RMI(:,:,2), &
    &         D2RMI(:,:,3), D2RMI(:,:,4), D2RMI(:,:,5), D2RMI(:,:,6), GTEST, STEST)

! Now calculate the rotation matrix acting on the Cartesian unit vectors
  XB(:,I) = MATMUL(RMI(:,:), X0)
  YB(:,I) = MATMUL(RMI(:,:), Y0)
  ZB(:,I) = MATMUL(RMI(:,:), Z0)

! Calculate first derivatives of the rotation matrix acting on the Cartesian unit vectors, if required
  IF (GTEST .OR. STEST) THEN
    DO J = 1, 3
      DX(:,3*I-3+J) = MATMUL(DRMI(:,:,J), X0)
      DY(:,3*I-3+J) = MATMUL(DRMI(:,:,J), Y0)
      DZ(:,3*I-3+J) = MATMUL(DRMI(:,:,J), Z0)
    END DO ! J loop over derivative matrixes
  END IF ! end if GTEST or STEST

! Calculate second derivatives of the rotation matrix acting on the Cartesian unit vectors, if required
  IF (STEST) THEN
    DO J = 1, 6
      D2X(:,6*I-6+J) = MATMUL(D2RMI(:,:,J), X0)
      D2Y(:,6*I-6+J) = MATMUL(D2RMI(:,:,J), Y0)
      D2Z(:,6*I-6+J) = MATMUL(D2RMI(:,:,J), Z0)
    END DO ! J loop over coordinate pairs
  END IF ! end if STEST

END DO ! I loop over each body

! Calculate various sums required for the potential and derivatives
! Start by initialising
SX2 = REALATOMS
SY2 = REALATOMS
SZ2 = REALATOMS
SUMXJEXP(:,:) = 0.0
SUMYJEXP(:,:) = 0.0
SUMZJEXP(:,:) = 0.0
SUMXJEXPRIJ(:,:) = 0.0
SUMYJEXPRIJ(:,:) = 0.0
SUMZJEXPRIJ(:,:) = 0.0
SUMXJEXPDRIJ(:,:) = 0.0
SUMYJEXPDRIJ(:,:) = 0.0
SUMZJEXPDRIJ(:,:) = 0.0

DO I = 1, REALATOMS
! Set IT as the third translational coordinate of I, and IR as the third rotational coordinate
  IT = 3*I
  IR = IT + OFFSET

  DO J = 1, REALATOMS
! Only include cross terms
    IF (I .EQ. J) CYCLE
! Set JT as the third translational coordinate of J
    JT = 3*J

! Find the vector between bodies I and J, and it's length
    RIJ(:) = X(IT-2:IT) - X(JT-2:JT)
    LRIJ2 = DOT_PRODUCT(RIJ(:), RIJ(:))
! Calculate the Gaussian range term
    EXPFACT = EXP(-LRIJ2 / (2 * SIGMA2))

! First find the terms for the potential
    SX2 = SX2 + DOT_PRODUCT(XB(:,I), XB(:,J)) * EXPFACT
    SY2 = SY2 + DOT_PRODUCT(YB(:,I), YB(:,J)) * EXPFACT
    SZ2 = SZ2 + DOT_PRODUCT(ZB(:,I), ZB(:,J)) * EXPFACT

! Sums needed for first and second derivatives
    IF (GTEST .OR. STEST) THEN
      SUMXJEXP(:,I) = SUMXJEXP(:,I) + XB(:,J) * EXPFACT
      SUMYJEXP(:,I) = SUMYJEXP(:,I) + YB(:,J) * EXPFACT
      SUMZJEXP(:,I) = SUMZJEXP(:,I) + ZB(:,J) * EXPFACT
      DO K = 1, 3
        SUMXJEXPRIJ(:,IT-3+K) = SUMXJEXPRIJ(:,IT-3+K) + XB(:,J) * EXPFACT * RIJ(K)
        SUMYJEXPRIJ(:,IT-3+K) = SUMYJEXPRIJ(:,IT-3+K) + YB(:,J) * EXPFACT * RIJ(K)
        SUMZJEXPRIJ(:,IT-3+K) = SUMZJEXPRIJ(:,IT-3+K) + ZB(:,J) * EXPFACT * RIJ(K)
      END DO ! K loop over coordinates
    END IF ! end if GTEST or STEST

! Sums needed for second derivatives
    IF (STEST) THEN
      SUMXJEXPDRIJ(:,6*I-5) = SUMXJEXPDRIJ(:,6*I-5) + XB(:,J) * EXPFACT * (1.0D0 - RIJ(1)**2 / SIGMA2)
      SUMXJEXPDRIJ(:,6*I-4) = SUMXJEXPDRIJ(:,6*I-4) + XB(:,J) * EXPFACT * (1.0D0 - RIJ(2)**2 / SIGMA2)
      SUMXJEXPDRIJ(:,6*I-3) = SUMXJEXPDRIJ(:,6*I-3) + XB(:,J) * EXPFACT * (1.0D0 - RIJ(3)**2 / SIGMA2)
      SUMXJEXPDRIJ(:,6*I-2) = SUMXJEXPDRIJ(:,6*I-2) - XB(:,J) * EXPFACT * RIJ(1) * RIJ(2) / SIGMA2 
      SUMXJEXPDRIJ(:,6*I-1) = SUMXJEXPDRIJ(:,6*I-1) - XB(:,J) * EXPFACT * RIJ(2) * RIJ(3) / SIGMA2
      SUMXJEXPDRIJ(:,6*I  ) = SUMXJEXPDRIJ(:,6*I  ) - XB(:,J) * EXPFACT * RIJ(3) * RIJ(1) / SIGMA2

      SUMYJEXPDRIJ(:,6*I-5) = SUMYJEXPDRIJ(:,6*I-5) + YB(:,J) * EXPFACT * (1.0D0 - RIJ(1)**2 / SIGMA2)
      SUMYJEXPDRIJ(:,6*I-4) = SUMYJEXPDRIJ(:,6*I-4) + YB(:,J) * EXPFACT * (1.0D0 - RIJ(2)**2 / SIGMA2)
      SUMYJEXPDRIJ(:,6*I-3) = SUMYJEXPDRIJ(:,6*I-3) + YB(:,J) * EXPFACT * (1.0D0 - RIJ(3)**2 / SIGMA2)
      SUMYJEXPDRIJ(:,6*I-2) = SUMYJEXPDRIJ(:,6*I-2) - YB(:,J) * EXPFACT * RIJ(1) * RIJ(2) / SIGMA2 
      SUMYJEXPDRIJ(:,6*I-1) = SUMYJEXPDRIJ(:,6*I-1) - YB(:,J) * EXPFACT * RIJ(2) * RIJ(3) / SIGMA2
      SUMYJEXPDRIJ(:,6*I  ) = SUMYJEXPDRIJ(:,6*I  ) - YB(:,J) * EXPFACT * RIJ(3) * RIJ(1) / SIGMA2

      SUMZJEXPDRIJ(:,6*I-5) = SUMZJEXPDRIJ(:,6*I-5) + ZB(:,J) * EXPFACT * (1.0D0 - RIJ(1)**2 / SIGMA2)
      SUMZJEXPDRIJ(:,6*I-4) = SUMZJEXPDRIJ(:,6*I-4) + ZB(:,J) * EXPFACT * (1.0D0 - RIJ(2)**2 / SIGMA2)
      SUMZJEXPDRIJ(:,6*I-3) = SUMZJEXPDRIJ(:,6*I-3) + ZB(:,J) * EXPFACT * (1.0D0 - RIJ(3)**2 / SIGMA2)
      SUMZJEXPDRIJ(:,6*I-2) = SUMZJEXPDRIJ(:,6*I-2) - ZB(:,J) * EXPFACT * RIJ(1) * RIJ(2) / SIGMA2 
      SUMZJEXPDRIJ(:,6*I-1) = SUMZJEXPDRIJ(:,6*I-1) - ZB(:,J) * EXPFACT * RIJ(2) * RIJ(3) / SIGMA2
      SUMZJEXPDRIJ(:,6*I  ) = SUMZJEXPDRIJ(:,6*I  ) - ZB(:,J) * EXPFACT * RIJ(3) * RIJ(1) / SIGMA2
    END IF ! end if STEST
  END DO ! J loop over each body
END DO ! I loop over each body

! Find the energy
ENERGY = ENERGY - KAA * (SX2 + SY2 + SZ2) / REALATOMS**2

! Find the gradients
IF (GTEST) THEN
! Loop over each translational (or rotational) coordinate
  DO IT = 1, 3*REALATOMS
    I = (IT+2)/3
! First the contributions to the translational derivatives
    DSX2DR = - DOT_PRODUCT(XB(:,I),SUMXJEXPRIJ(:,IT))
    DSY2DR = - DOT_PRODUCT(YB(:,I),SUMYJEXPRIJ(:,IT))
    DSZ2DR = - DOT_PRODUCT(ZB(:,I),SUMZJEXPRIJ(:,IT))
    G(IT) = G(IT) - 2 * KAA * (DSX2DR + DSY2DR + DSZ2DR) / (REALATOMS**2 * SIGMA2)
! And then the contributions to the rotational derivatives
    DSX2DP = DOT_PRODUCT(DX(:,IT),SUMXJEXP(:,I))
    DSY2DP = DOT_PRODUCT(DY(:,IT),SUMYJEXP(:,I))
    DSZ2DP = DOT_PRODUCT(DZ(:,IT),SUMZJEXP(:,I))
    G(IT+OFFSET) = G(IT+OFFSET) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / REALATOMS**2
  END DO ! IT loop over each coordinate
END IF ! end if GTEST

! And finally for the Hessian
IF (STEST) THEN
  DO IT = 1, 3*REALATOMS
    DO JT = 1, 3*REALATOMS
! Work out some indices
      I  = (IT+2)/3
      IR = IT+OFFSET
      J  = (JT+2)/3
      JR = JT + OFFSET
! Start with a few necessary factors
      RIJ(:) = X(3*I-2:3*I) - X(3*J-2:3*J)
      LRIJ2 = DOT_PRODUCT(RIJ(:), RIJ(:))
      EXPFACT = EXP(-LRIJ2 / (2 * SIGMA2))

      IF(I .EQ. J) THEN
! Diagonal (same body)
! Set up the I2 index in the same way as SUMkJEXPDRIJ and D2k
        IF (IT .EQ. JT) THEN 
          I2 = 6*I - 5 + MOD(IT-1,3)
        ELSE IF ((MOD(IT,3) .EQ. 1 .AND. MOD(JT,3) .EQ. 2) .OR. (MOD(IT,3) .EQ. 2 .AND. MOD(JT,3) .EQ. 1)) THEN 
          I2 = 6*I - 2
        ELSE IF ((MOD(IT,3) .EQ. 2 .AND. MOD(JT,3) .EQ. 0) .OR. (MOD(IT,3) .EQ. 0 .AND. MOD(JT,3) .EQ. 2)) THEN 
          I2 = 6*I - 1
        ELSE IF ((MOD(IT,3) .EQ. 0 .AND. MOD(JT,3) .EQ. 1) .OR. (MOD(IT,3) .EQ. 1 .AND. MOD(JT,3) .EQ. 0)) THEN 
          I2 = 6*I
        ENDIF 
! Pure translation
        DSX2DR = - DOT_PRODUCT(XB(:,I), SUMXJEXPDRIJ(:,I2))
        DSY2DR = - DOT_PRODUCT(YB(:,I), SUMYJEXPDRIJ(:,I2))
        DSZ2DR = - DOT_PRODUCT(ZB(:,I), SUMZJEXPDRIJ(:,I2))
        HESS(IT,JT) = HESS(IT,JT) - 2 * KAA * (DSX2DR + DSY2DR + DSZ2DR) / (REALATOMS**2 * SIGMA2)

! Mixed rotation and translation
        DSX2DP = - DOT_PRODUCT(DX(:,JT), SUMXJEXPRIJ(:,IT))
        DSY2DP = - DOT_PRODUCT(DY(:,JT), SUMYJEXPRIJ(:,IT))
        DSZ2DP = - DOT_PRODUCT(DZ(:,JT), SUMZJEXPRIJ(:,IT))
        HESS(IT,JR) = HESS(IT,JR) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / (REALATOMS**2 * SIGMA2)
        DSX2DP = - DOT_PRODUCT(DX(:,IT), SUMXJEXPRIJ(:,JT))
        DSY2DP = - DOT_PRODUCT(DY(:,IT), SUMYJEXPRIJ(:,JT))
        DSZ2DP = - DOT_PRODUCT(DZ(:,IT), SUMZJEXPRIJ(:,JT))
        HESS(IR,JT) = HESS(IR,JT) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / (REALATOMS**2 * SIGMA2)

! Pure rotation
        DSX2DP = DOT_PRODUCT(D2X(:,I2), SUMXJEXP(:,I))
        DSY2DP = DOT_PRODUCT(D2Y(:,I2), SUMYJEXP(:,I))
        DSZ2DP = DOT_PRODUCT(D2Z(:,I2), SUMZJEXP(:,I))
        HESS(IR,JR) = HESS(IR,JR) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / REALATOMS**2

      ELSE
! Off diagonal (different body)
! Set up I2 and J2 to be 1 to 3
        I2 = MOD(IT-1,3) + 1
        J2 = MOD(JT-1,3) + 1

! Pure translation
        IF (I2 .EQ. J2) THEN
! Same coordinate
          DSX2DR = DOT_PRODUCT(XB(:,I), XB(:,J)) * EXPFACT * (1.0D0 - RIJ(I2)**2 / SIGMA2)
          DSY2DR = DOT_PRODUCT(YB(:,I), YB(:,J)) * EXPFACT * (1.0D0 - RIJ(I2)**2 / SIGMA2)
          DSZ2DR = DOT_PRODUCT(ZB(:,I), ZB(:,J)) * EXPFACT * (1.0D0 - RIJ(I2)**2 / SIGMA2)
        ELSE
! Different coodinate
          DSX2DR = - DOT_PRODUCT(XB(:,I), XB(:,J)) * EXPFACT * RIJ(I2) * RIJ(J2) / SIGMA2 
          DSY2DR = - DOT_PRODUCT(YB(:,I), YB(:,J)) * EXPFACT * RIJ(I2) * RIJ(J2) / SIGMA2 
          DSZ2DR = - DOT_PRODUCT(ZB(:,I), ZB(:,J)) * EXPFACT * RIJ(I2) * RIJ(J2) / SIGMA2 
        END IF ! end if same coordinate
        HESS(IT,JT) = HESS(IT,JT) - 2 * KAA * (DSX2DR + DSY2DR + DSZ2DR) / (REALATOMS**2 * SIGMA2)

! Mixed translation and rotation
        DSX2DP = - DOT_PRODUCT(XB(:,I), DX(:,JT)) * EXPFACT * RIJ(I2)
        DSY2DP = - DOT_PRODUCT(YB(:,I), DY(:,JT)) * EXPFACT * RIJ(I2)
        DSZ2DP = - DOT_PRODUCT(ZB(:,I), DZ(:,JT)) * EXPFACT * RIJ(I2)
        HESS(IT,JR) = HESS(IT,JR) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / (REALATOMS**2 * SIGMA2)
        DSX2DP = DOT_PRODUCT(XB(:,J), DX(:,IT)) * EXPFACT * RIJ(J2)
        DSY2DP = DOT_PRODUCT(YB(:,J), DY(:,IT)) * EXPFACT * RIJ(J2)
        DSZ2DP = DOT_PRODUCT(ZB(:,J), DZ(:,IT)) * EXPFACT * RIJ(J2)
        HESS(IR,JT) = HESS(IR,JT) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / (REALATOMS**2 * SIGMA2)

! Pure rotation
        DSX2DP = DOT_PRODUCT(DX(:,IT), DX(:,JT)) * EXPFACT
        DSY2DP = DOT_PRODUCT(DY(:,IT), DY(:,JT)) * EXPFACT
        DSZ2DP = DOT_PRODUCT(DZ(:,IT), DZ(:,JT)) * EXPFACT
        HESS(IR,JR) = HESS(IR,JR) - 2 * KAA * (DSX2DP + DSY2DP + DSZ2DP) / REALATOMS**2
      END IF ! end if same body
    END DO ! JT loop over each coordinate
  END DO ! IT loop over each coordinate
END IF ! end if STEST
END SUBROUTINE AAORIENTSR
