      SUBROUTINE PTSTST (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE 
      USE KEY, ONLY: PAPEPS, PAPCD, PAPID  

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, JL, JU, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), LJN 
      DOUBLE PRECISION :: ENERGY, RIJ(3), RIJSQ, R2, RSS(3), RSSSQ, FCTR, EPS, SIGMA, SIGSQ, RLJN, R2LJN, LJSIG2
      DOUBLE PRECISION :: RI(3), P(3), R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DVDR(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), D2VDR2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R1(NRBSITES*NATOMS/2,3), D2R2(NRBSITES*NATOMS/2,3), D2R3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R12(NRBSITES*NATOMS/2,3), D2R23(NRBSITES*NATOMS/2,3), D2R31(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DOTI1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTI2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTI3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      DOUBLE PRECISION :: DOTJ1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTJ2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTJ3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DVDRCC(NATOMS/2,NATOMS/2), D2VDR2CC(NATOMS/2,NATOMS/2), DUMMY
      LOGICAL          :: GTEST, STEST

      D2VDR2(:,:) = 0.D0; DVDR(:,:) = 0.D0
      DVDRCC(1:NATOMS/2,1:NATOMS/2) = 0.D0; D2VDR2CC(1:NATOMS/2,1:NATOMS/2) = 0.D0
      DOTI1(:,:) = 0.D0; DOTI2(:,:) = 0.D0; DOTI3(:,:) = 0.D0
      DOTJ1(:,:) = 0.D0; DOTJ2(:,:) = 0.D0; DOTJ3(:,:) = 0.D0

      EPS = PAPEPS; SIGMA = PAPCD
      LJSIG2 = 1.D0; LJN = 6.D0
      SIGSQ = SIGMA*SIGMA

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (STEST) HESS(:,:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
 
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI(:) = X(J3-2:J3)
         P(:)  = X(J5-2:J5)

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, STEST)

         DO J2 = 1, NRBSITES-1
            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,RBSITE(J2,:))
            IF (GTEST .OR. STEST) THEN
               DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))
            ENDIF
            IF (STEST) THEN
               D2R1(J4,:) = MATMUL(D2RMI1,RBSITE(J2,:))
               D2R2(J4,:) = MATMUL(D2RMI2,RBSITE(J2,:))
               D2R3(J4,:) = MATMUL(D2RMI3,RBSITE(J2,:))
                  
               D2R12(J4,:) = MATMUL(D2RMI12,RBSITE(J2,:))
               D2R23(J4,:) = MATMUL(D2RMI23,RBSITE(J2,:))
               D2R31(J4,:) = MATMUL(D2RMI31,RBSITE(J2,:))
            ENDIF
         ENDDO
      ENDDO

      DO J1 = 1, REALNATOMS 
         J3 = 3*J1
         J5 = OFFSET + J3
         DO J2 = J1 + 1, REALNATOMS
            J4 = 3*J2
            J6 = OFFSET + J4

            RIJ(:) = X(J3-2:J3) - X(J4-2:J4)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            R2     = 1.D0/RIJSQ
            RLJN   = (LJSIG2*R2)**(LJN/2.D0)
            R2LJN  = RLJN**2
            ENERGY = ENERGY + 4.D0*(R2LJN-RLJN)

            IF (GTEST .OR. STEST) THEN
               DVDRCC(J1,J2) = 4.D0*LJN*R2*(RLJN - 2.D0*R2LJN)
               DVDRCC(J2,J1) = DVDRCC(J1,J2)
               G(J3-2:J3)  = G(J3-2:J3) + DVDRCC(J1,J2)*RIJ(:)
               G(J4-2:J4)  = G(J4-2:J4) - DVDRCC(J1,J2)*RIJ(:)
               D2VDR2CC(J1,J2) = 4.D0*LJN*R2*R2*(4.D0*(LJN+1.D0)*R2LJN - (LJN+2.D0)*RLJN)
               D2VDR2CC(J2,J1) = D2VDR2CC(J1,J2)
            ENDIF

            DO I = 1, NRBSITES-1 
               J7 = NRBSITES*(J1-1) + I
               IF (I <= (NRBSITES-1)/2) THEN
                  JL = (NRBSITES-1)/2 + 1
                  JU = NRBSITES-1
               ELSE
                  JL = 1
                  JU = (NRBSITES-1)/2
               ENDIF
               DO J = JL, JU 
                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  RSSSQ  = DOT_PRODUCT(RSS(:),RSS(:))
                  FCTR   = EXP(-RSSSQ/(2.D0*SIGSQ))
                  ENERGY = ENERGY - EPS*FCTR

                  IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R
                     DVDR(J7,J8) = EPS*FCTR/SIGSQ
                     DVDR(J8,J7) = DVDR(J7,J8)

                     DOTI1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                     DOTI2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                     DOTI3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J7,:))
                     DOTJ1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J8,:))
                     DOTJ2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J8,:))
                     DOTJ3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J8,:))

                     G(J3-2:J3)  = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
                     G(J4-2:J4)  = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)
                     G(J5-2)     = G(J5-2) + DVDR(J7,J8)*DOTI1(J7,J8)
                     G(J5-1)     = G(J5-1) + DVDR(J7,J8)*DOTI2(J7,J8)
                     G(J5)       = G(J5)   + DVDR(J7,J8)*DOTI3(J7,J8)
                     G(J6-2)     = G(J6-2) - DVDR(J7,J8)*DOTJ1(J7,J8)
                     G(J6-1)     = G(J6-1) - DVDR(J7,J8)*DOTJ2(J7,J8)
                     G(J6)       = G(J6)   - DVDR(J7,J8)*DOTJ3(J7,J8)

                     D2VDR2(J7,J8) =-DVDR(J7,J8)/SIGSQ
                  ENDIF

                  IF (STEST) THEN
                     D2VDR2(J8,J7) = D2VDR2(J7,J8)
                     DOTI1(J8,J7)  =-DOTJ1(J7,J8)
                     DOTI2(J8,J7)  =-DOTJ2(J7,J8)
                     DOTI3(J8,J7)  =-DOTJ3(J7,J8)
                     DOTJ1(J8,J7)  =-DOTI1(J7,J8)
                     DOTJ2(J8,J7)  =-DOTI2(J7,J8)
                     DOTJ3(J8,J7)  =-DOTI3(J7,J8)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!     Centre-Centre Mie contribution
      IF (STEST) THEN
         DO J1 = 1, REALNATOMS
            J3 = 3*J1
            DO J2 = 1, REALNATOMS
               IF (J1 == J2) CYCLE
               J4 = 3*J2
               RIJ(:) = X(J3-2:J3) - X(J4-2:J4)
!     [1] THREE COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES

!     xi,xi
               HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2CC(J1,J2)*RIJ(1)*RIJ(1) + DVDRCC(J1,J2)
!     yi,yi             
               HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2CC(J1,J2)*RIJ(2)*RIJ(2) + DVDRCC(J1,J2)
!     zi,zi
               HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2CC(J1,J2)*RIJ(3)*RIJ(3) + DVDRCC(J1,J2)
!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME MOLECULE, DIFFERENT COORDINATES
!     xi,yi
               DUMMY           = D2VDR2CC(J1,J2)*RIJ(1)*RIJ(2)
               HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
               HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
               DUMMY           = D2VDR2CC(J1,J2)*RIJ(2)*RIJ(3)
               HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
               HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
               DUMMY           = D2VDR2CC(J1,J2)*RIJ(3)*RIJ(1)
               HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
               HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE
!     xi,xj
               HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2CC(J1,J2)*RIJ(1)*RIJ(1) - DVDRCC(J1,J2)
!     yi,yj
               HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2CC(J1,J2)*RIJ(2)*RIJ(2) - DVDRCC(J1,J2)
!     zi,zj
               HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2CC(J1,J2)*RIJ(3)*RIJ(3) - DVDRCC(J1,J2)
!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES
!     xi,yj
               DUMMY           =-D2VDR2CC(J1,J2)*RIJ(1)*RIJ(2)
               HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
               HESS(J4-1,J3-2) = HESS(J4-1,J3-2) + DUMMY
!     yi,zj
               DUMMY           =-D2VDR2CC(J1,J2)*RIJ(2)*RIJ(3)
               HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
               HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
               DUMMY           =-D2VDR2CC(J1,J2)*RIJ(3)*RIJ(1)
               HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
               HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY
            ENDDO
         ENDDO
      ENDIF
!     Site-Site contribution
      IF (STEST) THEN
         DO J1 = 1, REALNATOMS
            J3 = 3*J1
            J5 = OFFSET + J3
            DO J2 = 1, REALNATOMS
               IF (J1 == J2) CYCLE
               J4 = 3*J2
               J6 = OFFSET + J4
               DO I = 1, NRBSITES-1
                  J7 = NRBSITES*(J1 - 1) + I
                  IF (I <= (NRBSITES-1)/2) THEN
                     JL = (NRBSITES-1)/2 + 1
                     JU = NRBSITES-1
                  ELSE
                     JL = 1
                     JU = (NRBSITES-1)/2
                  ENDIF
                  DO J = JL, JU !1, NRBSITES 
                     J8 = NRBSITES*(J2 - 1) + J
                     RSS(:) = R(J7,:) - R(J8,:) 
!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES
!     xi,xi
                     HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2(J7,J8)*RSS(1)*RSS(1) + DVDR(J7,J8)
!     yi,yi             
                     HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2(J7,J8)*RSS(2)*RSS(2) + DVDR(J7,J8)
!     zi,zi
                     HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2(J7,J8)*RSS(3)*RSS(3) + DVDR(J7,J8)
!     pi1,pi1
                     HESS(J5-2,J5-2) = HESS(J5-2,J5-2) + D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTI1(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR1(J7,:),DR1(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R1(J7,:))
!     pi2,pi2
                     HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTI2(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R2(J7,:))
!     pi3,pi3
                     HESS(J5,J5)     = HESS(J5,J5) + D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTI3(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R3(J7,:))
!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME MOLECULE, DIFFERENT COORDINATES
!     xi,yi
                     DUMMY           = D2VDR2(J7,J8)*RSS(1)*RSS(2)
                     HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
                     HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
                     DUMMY           = D2VDR2(J7,J8)*RSS(2)*RSS(3)
                     HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
                     HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
                     DUMMY           = D2VDR2(J7,J8)*RSS(3)*RSS(1)
                     HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
                     HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
!     xi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(1) + DVDR(J7,J8)*DR1(J7,1)
                     HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
                     HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(2) + DVDR(J7,J8)*DR1(J7,2)
                     HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
                     HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(3) + DVDR(J7,J8)*DR1(J7,3)
                     HESS(J3,J5-2)   = HESS(J3,J5-2) + DUMMY
                     HESS(J5-2,J3)   = HESS(J5-2,J3) + DUMMY
!     xi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(1) + DVDR(J7,J8)*DR2(J7,1)
                     HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
                     HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(2) + DVDR(J7,J8)*DR2(J7,2)
                     HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
                     HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(3) + DVDR(J7,J8)*DR2(J7,3)
                     HESS(J3,J5-1)   = HESS(J3,J5-1) + DUMMY
                     HESS(J5-1,J3)   = HESS(J5-1,J3) + DUMMY
!     xi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(1) + DVDR(J7,J8)*DR3(J7,1)
                     HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
                     HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(2) + DVDR(J7,J8)*DR3(J7,2)
                     HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
                     HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(3) + DVDR(J7,J8)*DR3(J7,3)
                     HESS(J3,J5)     = HESS(J3,J5) + DUMMY
                     HESS(J5,J3)     = HESS(J5,J3) + DUMMY
!     pi1,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTI2(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR2(J7,:),DR1(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R12(J7,:))
                     HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
                     HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTI3(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR3(J7,:),DR2(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R23(J7,:))
                     HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
                     HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTI1(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR1(J7,:),DR3(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R31(J7,:))
                     HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
                     HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY
!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE
!     xi,xj
                     HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2(J7,J8)*RSS(1)*RSS(1) - DVDR(J7,J8)
!     yi,yj
                     HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2(J7,J8)*RSS(2)*RSS(2) - DVDR(J7,J8)
!     zi,zj
                     HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2(J7,J8)*RSS(3)*RSS(3) - DVDR(J7,J8)
!     pi1,pj1
                     HESS(J5-2,J6-2) = HESS(J5-2,J6-2) - D2VDR2(J7,J8)*DOTJ1(J7,J8)*DOTI1(J7,J8) &
                                     - DVDR(J7,J8)*DOT_PRODUCT(DR1(J8,:),DR1(J7,:))
!     pi2,pj2
                     HESS(J5-1,J6-1) = HESS(J5-1,J6-1) - D2VDR2(J7,J8)*DOTJ2(J7,J8)*DOTI2(J7,J8) &
                                     - DVDR(J7,J8)*DOT_PRODUCT(DR2(J8,:),DR2(J7,:))
!     pi3,pj3
                     HESS(J5,J6)     = HESS(J5,J6)     - D2VDR2(J7,J8)*DOTJ3(J7,J8)*DOTI3(J7,J8) &
                                    - DVDR(J7,J8)*DOT_PRODUCT(DR3(J8,:),DR3(J7,:))
!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES
!     xi,yj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(1)*RSS(2)
                     HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
                     HESS(J4-1,J3-2) = HESS(J4-1,J3-2) + DUMMY
!     yi,zj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(2)*RSS(3)
                     HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
                     HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(3)*RSS(1)
                     HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
                     HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY
!     xi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(1) - DVDR(J7,J8)*DR1(J8,1)
                     HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
                     HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(2) - DVDR(J7,J8)*DR1(J8,2)
                     HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
                     HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(3) - DVDR(J7,J8)*DR1(J8,3)
                     HESS(J3,J6-2)   = HESS(J3,J6-2) + DUMMY
                     HESS(J6-2,J3)   = HESS(J6-2,J3) + DUMMY
!     xi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(1) - DVDR(J7,J8)*DR2(J8,1)
                     HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
                     HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(2) - DVDR(J7,J8)*DR2(J8,2)
                     HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
                     HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(3) - DVDR(J7,J8)*DR2(J8,3)
                     HESS(J3,J6-1)   = HESS(J3,J6-1) + DUMMY
                     HESS(J6-1,J3)   = HESS(J6-1,J3) + DUMMY
!     xi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(1) - DVDR(J7,J8)*DR3(J8,1)
                     HESS(J3-2,J6)   = HESS(J3-2,J6) + DUMMY
                     HESS(J6,J3-2)   = HESS(J6,J3-2) + DUMMY
!     yi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(2) - DVDR(J7,J8)*DR3(J8,2)
                     HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
                     HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(3) - DVDR(J7,J8)*DR3(J8,3)
                     HESS(J3,J6)     = HESS(J3,J6) + DUMMY
                     HESS(J6,J3)     = HESS(J6,J3) + DUMMY
!     pi1,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTJ2(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR2(J8,:),DR1(J7,:))
                     HESS(J5-2,J6-1) = HESS(J5-2,J6-1) + DUMMY
                     HESS(J6-1,J5-2) = HESS(J6-1,J5-2) + DUMMY
!     pi2,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTJ3(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR3(J8,:),DR2(J7,:))
                     HESS(J5-1,J6)   = HESS(J5-1,J6) + DUMMY
                     HESS(J6,J5-1)   = HESS(J6,J5-1) + DUMMY
!     pi3,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTJ1(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR1(J8,:),DR3(J7,:))
                     HESS(J5,J6-2)   = HESS(J5,J6-2) + DUMMY
                     HESS(J6-2,J5)   = HESS(J6-2,J5) + DUMMY
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE PTSTST 

!     ---------------------------------------------------------------------------------------------

      SUBROUTINE DEFPTSTST()

      USE COMMONS, ONLY: RBSITE
      USE KEY, ONLY: PAPID 

      IMPLICIT NONE

      IF (PAPID == 3) THEN
         RBSITE(1,:) = 0.5D0*(/ -0.096225D0, 0.301232D0, 0.948682D0/)
         RBSITE(2,:) = 0.5D0*(/ -0.096225D0,-0.301232D0,-0.948682D0/)
         RBSITE(3,:) = 0.D0
      ELSE
         PRINT *, 'PAPID should be 3 as of now'
      ENDIF
      END SUBROUTINE DEFPTSTST
