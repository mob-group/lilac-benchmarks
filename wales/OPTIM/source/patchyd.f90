!|gd351>

SUBROUTINE PATCHYD (X, G, ENERGY, GTEST, SECT)
  USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, RBSTLA
  USE MODHESS

  IMPLICIT NONE
  
  INTEGER          :: I, I1, I2, I3, J, J1, J2, J3, J4, J5, J6, J7, J8, J3T, J4T, REALNATOMS, OFFSET 
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
  DOUBLE PRECISION :: SIGMASQ, RANGESQ, FACTOR
  DOUBLE PRECISION :: ENERGY, R, RSQ, RCUB, R2, R6, RHALF, ARG1, ARG2, ALPHA1, ALPHA2, VLJ, VEXP1, VEXP2, VEXP, V, DVLJ, A0, A1, A2
  DOUBLE PRECISION :: DARG1(6), DALPHA1(6), DARG2(6), DALPHA2(6), DVDX(3)
  DOUBLE PRECISION :: RI(NATOMS/2,3), RIJ(3), P(NATOMS/2,3), A(3) 
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
  DOUBLE PRECISION :: PATCHPOS(NATOMS/2,NRBSITES,3), DPATCHPOS1(NATOMS/2,NRBSITES,3), &
  &                   DPATCHPOS2(NATOMS/2,NRBSITES,3), DPATCHPOS3(NATOMS/2,NRBSITES,3)
  DOUBLE PRECISION :: TEMPX, TEMPLEFT, TEMPRIGHT, TEMPPATCHPOS(NRBSITES,3)
  LOGICAL          :: GTEST, SECT, LOGICAL
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
  DOUBLE PRECISION :: ALPHAMIN, ALPHA1T, ALPHA2T, VF1, DVF1, VF2, DVF2, DVEXP(1:3), DVF1DX(1:3)
  DOUBLE PRECISION :: S2A(0:3), S1A(0:3), AI(0:3), BI(0:3)


  !SIGMASQ: squared parameter sigma of Lennard-Jones potential
  !RANGESQ: squared range of potential
  !FACTOR: controls patchwidth
  !if SIGMASQ or RANGESQ are changed, parameters of smoothing functions VF1 and VF2 have to be adjusted; if FACTOR is changed, ALPHALEFT, ALPHARIGHT have to be adjusted
  SIGMASQ = 1.0D0
  RANGESQ = 3.61D0
  FACTOR =  (2.D0*PI*0.05)**2

  !parameters for smoothing function VF2 (polynomial of grade 3) - has to fulfill
  !VF2(0.9*range) = VLJ(0.9*range)
  !VF2(range) = 0
  !VF2'(0.9*range) = VLJ'(0.9*range)
  !VF2'(range) = 0
  S2A(0:3) =  (/ 4.3196421660790450D1, -7.2976415249550730D1, 4.0919975890836945D1, -7.6195284520433670D0 /)

  ENERGY = 0.D0
  IF (GTEST) G(:) = 0.D0

  IF (SECT) THEN
    PRINT*,'2nd derivatives not implemented (yet)'
    STOP
  END IF

  REALNATOMS = NATOMS/2
  OFFSET = 3*REALNATOMS
  
  DO J1 = 1, REALNATOMS
    
    J3 = 3*J1
    J5 = OFFSET + J3
    RI(J1,1:3) = X(J3-2:J3)
    P(J1,1:3)  = X(J5-2:J5)

    !calculate actual position of patches
    CALL RMDFAS(P(J1,1:3), RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, SECT)

    DO J2 = 1, NRBSITES
      PATCHPOS(J1,J2,:) = MATMUL(RMI(:,:),RBSITE(J2,:))
      IF (GTEST) THEN
        DPATCHPOS1(J1,J2,:) = MATMUL(DRMI1(:,:),RBSITE(J2,:))
        DPATCHPOS2(J1,J2,:) = MATMUL(DRMI2(:,:),RBSITE(J2,:))
        DPATCHPOS3(J1,J2,:) = MATMUL(DRMI3(:,:),RBSITE(J2,:))
      END IF
    END DO
       
  END DO


  DO J1 = 1, REALNATOMS

    DO J2 = J1+1, REALNATOMS

      RIJ(:) = RI(J1,:) - RI(J2,:)
      RSQ = DOT_PRODUCT(RIJ(:),RIJ(:))

      IF (RSQ.LE.0.81D0*SIGMASQ) THEN
      !core - LJ repulsion
        R2 = 1.D0/RSQ
        R6 = R2*R2*R2

        ENERGY = ENERGY + (R6 - 1.D0) * R6

        IF (GTEST) THEN
          !calculate derivatives
          DVLJ = -6.D0 * (2.D0 * R6 - 1.D0) * R6 * R2 !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVLJ * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVLJ * RIJ(:)
        END IF 

      ELSEIF ((RSQ.GT.0.81D0*SIGMASQ).AND.(RSQ.LE.SIGMASQ)) THEN
      !near edge of core - smoothend repulsion
        ALPHAMIN = 2.D0*PI

        R = SQRT(RSQ)
        RCUB = R*RSQ
        R2 = 1.D0 / RSQ
        RHALF = R / 2.D0

        IF (GTEST) THEN
          A(:) = -RIJ(:) * R2
        END IF

        DO J3T = 1, NRBSITES
          ARG1 = DOT_PRODUCT(PATCHPOS(J1,J3T,:),-RIJ(:)) / RHALF
          IF (ARG1.GE.1.D0) THEN
            ALPHA1T = 0.D0
          ELSEIF (ARG1.LE.-1.D0) THEN
            ALPHA1T = PI
          ELSE
            ALPHA1T = ACOS(ARG1)
          END IF
          DO J4T = 1, NRBSITES
            ARG2 = DOT_PRODUCT(PATCHPOS(J2,J4T,:),RIJ(:)) / RHALF
            IF (ARG2.GE.1.D0) THEN
              ALPHA2T = 0.D0
            ELSEIF (ARG2.LE.-1.D0) THEN
              ALPHA2T = PI
            ELSE
              ALPHA2T = ACOS(ARG2)
            END IF
            IF ((ALPHA1T+ALPHA2T).LE.ALPHAMIN) THEN
              ALPHAMIN=ALPHA1T+ALPHA2T
              ALPHA1=ALPHA1T
              ALPHA2=ALPHA2T
              J3=J3T
              J4=J4T
            END IF
          END DO
        END DO

        VEXP1 = EXP(-ALPHA1**2/FACTOR)              
        VEXP2 = EXP(-ALPHA2**2/FACTOR)
        VEXP = VEXP1 * VEXP2
        IF (GTEST) THEN
          IF ((ALPHA1.EQ.0.D0).OR.(ALPHA1.EQ.PI)) THEN
            DALPHA1(:)=0.D0
          ELSE 
            !derivatives of ARG1 wrt interparticle vector RIJ
            DARG1(1:3) = ( -PATCHPOS(J1,J3,:) + DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) * A(:) ) / RHALF
            !derivatives of ARG1 wrt orientational vector P of particles J1
            DARG1(4) = DOT_PRODUCT(DPATCHPOS1(J1,J3,:),-RIJ(:)) / RHALF
            DARG1(5) = DOT_PRODUCT(DPATCHPOS2(J1,J3,:),-RIJ(:)) / RHALF
            DARG1(6) = DOT_PRODUCT(DPATCHPOS3(J1,J3,:),-RIJ(:)) / RHALF
            DALPHA1(:) = -DARG1(:) / SQRT(1.D0-ARG1**2)
          END IF
          IF ((ALPHA2.EQ.0.D0).OR.(ALPHA2.EQ.PI)) THEN
            DALPHA2(:)=0.D0
          ELSE 
            !derivatives of ARG2 wrt interparticle vector RIJ
            DARG2(1:3) = ( PATCHPOS(J2,J4,:) + DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) * A(:) ) / RHALF
            !derivatives of ARG2 wrt orientational vector P of particle J2
            DARG2(4) = DOT_PRODUCT(DPATCHPOS1(J2,J4,:),RIJ(:)) / RHALF
            DARG2(5) = DOT_PRODUCT(DPATCHPOS2(J2,J4,:),RIJ(:)) / RHALF
            DARG2(6) = DOT_PRODUCT(DPATCHPOS3(J2,J4,:),RIJ(:)) / RHALF
            DALPHA2(:) = -DARG2(:) / SQRT(1.D0-ARG2**2)
          END IF
        END IF

        !parameters for smoothing function VF1 (polynomial of grade 3) - has to fulfill
        !VF1(0.9*sigma) = VLJ(0.9*sigma)
        !VF1(sigma) = VLJ(sigma)*VANG(THETA) = 0
        !VF1'(0.9*sigma) = VLJ'(0.9*sigma)
        !VF1'(sigma) = VLJ'(sigma)*VANG(THETA)
        AI(0:3) = (/ 299.49098473874074D0,-747.4130927079459D0,596.3532311996697D0,-148.4311232304645D0 /)
        BI(0:3) = (/ 485.999999999995D0,-1565.9999999999836D0,1679.9999999999825D0,-599.9999999999939D0 /)
        S1A(0:3) = AI(0:3) + BI(0:3)*VEXP
        VF1 = S1A(0) + S1A(1) * R + S1A(2) * RSQ + S1A(3) * RCUB
        ENERGY = ENERGY + VF1


        IF (GTEST) THEN
          !calculate derivatives
          DVF1 = S1A(1)/R + 2.D0*S1A(2) + 3.D0*S1A(3)*R
          A0 = ((-2.D0*VEXP)/FACTOR)*(BI(0)+BI(1)*R+BI(2)*RSQ+BI(3)*RCUB)
          DVEXP(1:3) = ALPHA1*DALPHA1(1:3) + ALPHA2*DALPHA2(1:3)
          DVF1DX(:) = DVF1 * RIJ(:) + A0 * DVEXP(:)
          !particle 1, derivatives wrt cartesian coordinates RI(J1,j) ( d (RIJ(i)) / d (RI(J1,j)) = delta(ij) )
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVF1DX(1:3)
          !particle 2, derivatives wrt cartesian coordinates RI(J2,j) ( d (RIJ(i)) / d (RI(J2,j)) = -delta(ij) )
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVF1DX(1:3)

          A1 = A0*ALPHA1
          A2 = A0*ALPHA2
          !particle 1, derivatives wrt orientational coordinates
          G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) = G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) + A1 * DALPHA1(4:6)
          !particle 2, derivatives wrt orientational coordinates
          G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) = G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) + A2 * DALPHA2(4:6)
        END IF

      ELSEIF ((RSQ.GT.SIGMASQ).AND.(RSQ.LE.0.81D0*RANGESQ)) THEN
      !patchy region - LJ attraction if patches are aligned
        R = SQRT(RSQ)
        R2 = 1.D0 / RSQ
        R6 = R2*R2*R2
        RHALF = R / 2.D0
        VLJ = (R6 - 1.D0) * R6


        IF (GTEST) THEN
          A(:) = -RIJ(:) * R2
          DVLJ = -6.D0 * (2.D0 * R6 - 1.D0) * R6 * R2 !factor 1/R from dR/dXi = Xi/R already included
        END IF

        DO J3 = 1, NRBSITES
          ARG1 = DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) / RHALF

          IF (ARG1.GE.1.D0) THEN
            ALPHA1 = 0.D0
            IF (GTEST) THEN
              DALPHA1(:) = 0.D0
            END IF 
          ELSEIF (ARG1.LE.-1.D0) THEN
            ALPHA1 = PI
            IF (GTEST) THEN
              DALPHA1(:) = 0.D0
            END IF 
          ELSE
            ALPHA1 = ACOS(ARG1)
            IF (GTEST) THEN
              !derivatives of ARG1 wrt interparticle vector RIJ
              DARG1(1:3) = ( -PATCHPOS(J1,J3,:) + DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) * A(:) ) / RHALF
              !derivatives of ARG1 wrt orientational vector P of particles J1
              DARG1(4) = DOT_PRODUCT(DPATCHPOS1(J1,J3,:),-RIJ(:)) / RHALF
              DARG1(5) = DOT_PRODUCT(DPATCHPOS2(J1,J3,:),-RIJ(:)) / RHALF
              DARG1(6) = DOT_PRODUCT(DPATCHPOS3(J1,J3,:),-RIJ(:)) / RHALF
              DALPHA1(:) = -DARG1(:) / SQRT(1.D0-ARG1**2)
            END IF 
          END IF
          !some kind of angle-cutoff could be added..
          VEXP1 = EXP(-ALPHA1**2/FACTOR)


          DO J4 = 1, NRBSITES
            ARG2 = DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) / RHALF
            IF (ARG2.GE.1.D0) THEN
              ALPHA2 = 0.D0
              IF (GTEST) THEN
                DALPHA2(:) = 0.D0
              END IF 
            ELSEIF (ARG2.LE.-1.D0) THEN
              ALPHA2 = PI
              IF (GTEST) THEN
                DALPHA2(:) = 0.D0
              END IF 
            ELSE
              ALPHA2 = ACOS(ARG2)
              IF (GTEST) THEN
                !derivatives of ARG2 wrt interparticle vector RIJ
                DARG2(1:3) = ( PATCHPOS(J2,J4,:) + DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) * A(:) ) / RHALF
                !derivatives of ARG2 wrt orientational vector P of particle J2
                DARG2(4) = DOT_PRODUCT(DPATCHPOS1(J2,J4,:),RIJ(:)) / RHALF
                DARG2(5) = DOT_PRODUCT(DPATCHPOS2(J2,J4,:),RIJ(:)) / RHALF
                DARG2(6) = DOT_PRODUCT(DPATCHPOS3(J2,J4,:),RIJ(:)) / RHALF
                DALPHA2(:) = -DARG2(:) / SQRT(1.D0-ARG2**2)
              END IF 
            END IF

            !some kind of angle-cutoff could be added..
            VEXP2 = EXP(-ALPHA2**2/FACTOR)
            VEXP = VEXP1 * VEXP2
            V = VLJ * VEXP

            ENERGY = ENERGY + V

            IF (GTEST) THEN

              A0 = -2.D0 * VLJ / FACTOR
              DVDX(1:3) = VEXP * (DVLJ * RIJ(:) + A0 * (ALPHA1*DALPHA1(1:3) + ALPHA2*DALPHA2(1:3)))

              !particle 1, derivatives wrt cartesian coordinates RI(J1,j) ( d (RIJ(i)) / d (RI(J1,j)) = delta(ij) )
              G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVDX(1:3)
              !particle 2, derivatives wrt cartesian coordinates RI(J2,j) ( d (RIJ(i)) / d (RI(J2,j)) = -delta(ij) )
              G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVDX(1:3)

              A1 = -2.D0 * V * ALPHA1 / FACTOR
              A2 = -2.D0 * V * ALPHA2 / FACTOR

              !particle 1, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) = G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) + A1 * DALPHA1(4:6)
              !particle 2, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) = G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) + A2 * DALPHA2(4:6)

            END IF 

          END DO
        END DO

      ELSEIF ((RSQ.GT.0.81*RANGESQ).AND.(RSQ.LE.RANGESQ)) THEN
      !patchy region near cutoff - smoothend attraction if patches are aligned
        R = SQRT(RSQ)
        R2 = 1.D0/RSQ
        RHALF = R / 2.D0
        VF2 = S2A(0) + S2A(1) * R + S2A(2) * RSQ + S2A(3) * R * RSQ


        IF (GTEST) THEN
          A(:) = -RIJ(:) * R2
          DVF2 = S2A(1)/R  + 2.D0*S2A(2)  + 3.D0*S2A(3) * R !factor 1/R from dR/dXi = Xi/R already included
        END IF


        DO J3 = 1, NRBSITES
          ARG1 = DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) / RHALF

          IF (ARG1.GE.1.D0) THEN
            ALPHA1 = 0.D0
            IF (GTEST) THEN
              DALPHA1(:) = 0.D0
            END IF 
          ELSEIF (ARG1.LE.-1.D0) THEN
            ALPHA1 = PI
            IF (GTEST) THEN
              DALPHA1(:) = 0.D0
            END IF 
          ELSE
            ALPHA1 = ACOS(ARG1)
            IF (GTEST) THEN
              !derivatives of ARG1 wrt interparticle vector RIJ
              DARG1(1:3) = ( -PATCHPOS(J1,J3,:) + DOT_PRODUCT(PATCHPOS(J1,J3,:),-RIJ(:)) * A(:) ) / RHALF
              !derivatives of ARG1 wrt orientational vector P of particles J1
              DARG1(4) = DOT_PRODUCT(DPATCHPOS1(J1,J3,:),-RIJ(:)) / RHALF
              DARG1(5) = DOT_PRODUCT(DPATCHPOS2(J1,J3,:),-RIJ(:)) / RHALF
              DARG1(6) = DOT_PRODUCT(DPATCHPOS3(J1,J3,:),-RIJ(:)) / RHALF
              DALPHA1(:) = -DARG1(:) / SQRT(1.D0-ARG1**2)
            END IF 
          END IF
          !some kind of angle-cutoff could be added..
          VEXP1 = EXP(-ALPHA1**2/FACTOR)


          DO J4 = 1, NRBSITES
            ARG2 = DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) / RHALF
            IF (ARG2.GE.1.D0) THEN
              ALPHA2 = 0.D0
              IF (GTEST) THEN
                DALPHA2(:) = 0.D0
              END IF 
            ELSEIF (ARG2.LE.-1.D0) THEN
              ALPHA2 = PI
              IF (GTEST) THEN
                DALPHA2(:) = 0.D0
              END IF 
            ELSE
              ALPHA2 = ACOS(ARG2)
              IF (GTEST) THEN
                !derivatives of ARG2 wrt interparticle vector RIJ
                DARG2(1:3) = ( PATCHPOS(J2,J4,:) + DOT_PRODUCT(PATCHPOS(J2,J4,:),RIJ(:)) * A(:) ) / RHALF
                !derivatives of ARG2 wrt orientational vector P of particle J2
                DARG2(4) = DOT_PRODUCT(DPATCHPOS1(J2,J4,:),RIJ(:)) / RHALF
                DARG2(5) = DOT_PRODUCT(DPATCHPOS2(J2,J4,:),RIJ(:)) / RHALF
                DARG2(6) = DOT_PRODUCT(DPATCHPOS3(J2,J4,:),RIJ(:)) / RHALF
                DALPHA2(:) = -DARG2(:) / SQRT(1.D0-ARG2**2)
              END IF 
            END IF

            !some kind of angle-cutoff could be added..
            VEXP2 = EXP(-ALPHA2**2/FACTOR)
            VEXP = VEXP1 * VEXP2
            
            V = VF2 * VEXP

            ENERGY = ENERGY + V

            IF (GTEST) THEN

              A0 = -2.D0 * VF2 / FACTOR
              DVDX(1:3) = VEXP * (DVF2 * RIJ(:) + A0 * (ALPHA1*DALPHA1(1:3) + ALPHA2*DALPHA2(1:3)))

              !particle 1, derivatives wrt cartesian coordinates RI(J1,j) ( d (RIJ(i)) / d (RI(J1,j)) = delta(ij) )
              G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DVDX(1:3)
              !particle 2, derivatives wrt cartesian coordinates RI(J2,j) ( d (RIJ(i)) / d (RI(J2,j)) = -delta(ij) )
              G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DVDX(1:3)

              A1 = -2.D0 * V * ALPHA1 / FACTOR
              A2 = -2.D0 * V * ALPHA2 / FACTOR

              !particle 1, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) = G(3*(REALNATOMS+J1)-2:3*(REALNATOMS+J1)) + A1 * DALPHA1(4:6)
              !particle 2, derivatives wrt orientational coordinates
              G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) = G(3*(REALNATOMS+J2)-2:3*(REALNATOMS+J2)) + A2 * DALPHA2(4:6)

            END IF 

          END DO
        END DO

      END IF

    END DO
  END DO

END SUBROUTINE PATCHYD



SUBROUTINE DEFPATCHES()
!distribution of patches on particle surface - this is called from keyword.f
!A: geometry parameter

  USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, RBSTLA

  IMPLICIT NONE
  DOUBLE PRECISION :: A
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0


  IF (NRBSITES.EQ.4) THEN

    A = 7.29808138D0
!    A = 7.60
    
!    RBSITE(1,1)=5.D-1*COS(0.D0)*SIN(A*PI/12.D0)
!    RBSITE(1,2)=5.D-1*SIN(0.D0)*SIN(A*PI/12.D0)
!    RBSITE(1,3)=5.D-1*COS(A*PI/12.D0)
!    RBSITE(2,1)=5.D-1*COS(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
!    RBSITE(2,2)=5.D-1*SIN(2.D0*PI/3.D0)*SIN(A*PI/12.D0)
!    RBSITE(2,3)=5.D-1*COS(A*PI/12.D0)
!    RBSITE(3,1)=5.D-1*COS(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
!    RBSITE(3,2)=5.D-1*SIN(4.D0*PI/3.D0)*SIN(A*PI/12.D0)
!    RBSITE(3,3)=5.D-1*COS(A*PI/12.D0)
!    RBSITE(4,1)=5.D-1*COS(0.D0)*SIN(0.D0)
!    RBSITE(4,2)=5.D-1*SIN(0.D0)*SIN(0.D0)
!    RBSITE(4,3)=5.D-1*COS(0.D0)

!direct definition of tetrahedral distribution (equivalent: A~7.29808138D0)
     RBSITE(1,:)= (5.D-1/SQRT(3.D0))*(/  1.D0,  1.D0,  1.D0/)
     RBSITE(2,:)= (5.D-1/SQRT(3.D0))*(/ -1.D0, -1.D0,  1.D0/)
     RBSITE(3,:)= (5.D-1/SQRT(3.D0))*(/ -1.D0,  1.D0, -1.D0/)
     RBSITE(4,:)= (5.D-1/SQRT(3.D0))*(/  1.D0, -1.D0, -1.D0/)

     RBSTLA(:,:)=RBSITE(:,:)

  ELSE
    PRINT*,'Patchnumber ',NRBSITES,' not implemented.'
    STOP
  END IF

END SUBROUTINE DEFPATCHES

!<gd351|
