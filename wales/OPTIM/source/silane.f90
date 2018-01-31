      SUBROUTINE SILANE (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
      USE KEY, ONLY: EPS11, EPS22, EPS12, MRHO11, MRHO22, MRHO12, REQ11, REQ22, REQ12  

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), FRQN(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, DSS, FCTR, RHO, REQ, EPS, DUMMY
      DOUBLE PRECISION :: RI(3), RSS(3), P(3), R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DVDR(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), D2VDR2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
      DOUBLE PRECISION :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
      DOUBLE PRECISION :: D2RMI10(3,3), D2RMI20(3,3), D2RMI30(3,3), D2RMI120(3,3), D2RMI230(3,3), D2RMI310(3,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R1(NRBSITES*NATOMS/2,3), D2R2(NRBSITES*NATOMS/2,3), D2R3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R12(NRBSITES*NATOMS/2,3), D2R23(NRBSITES*NATOMS/2,3), D2R31(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DOTI1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTI2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTI3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      DOUBLE PRECISION :: DOTJ1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTJ2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTJ3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      LOGICAL          :: GTEST, STEST

      D2VDR2(:,:) = 0.D0; DVDR(:,:) = 0.D0
      DOTI1(:,:) = 0.D0; DOTI2(:,:) = 0.D0; DOTI3(:,:) = 0.D0
      DOTJ1(:,:) = 0.D0; DOTJ2(:,:) = 0.D0; DOTJ3(:,:) = 0.D0

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (STEST) HESS(:,:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
 
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         P  = X(J5-2:J5)

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, STEST)

         IF (RBAANORMALMODET) THEN
            P = (/0.0D0, 0.0D0, 0.0D0/)
            CALL RMDFAS(P, RMI0, DRMI10, DRMI20, DRMI30, D2RMI10, D2RMI20, D2RMI30, D2RMI120, D2RMI230, D2RMI310, GTEST, STEST)
         ENDIF

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,RBSITE(J2,:))

            IF (GTEST .OR. STEST) THEN

               IF ( RBAANORMALMODET ) THEN
                  DR1(J4,:) = MATMUL(DRMI10,MATMUL(RMI,RBSITE(J2,:)))
                  DR2(J4,:) = MATMUL(DRMI20,MATMUL(RMI,RBSITE(J2,:)))
                  DR3(J4,:) = MATMUL(DRMI30,MATMUL(RMI,RBSITE(J2,:)))
               ELSE
                  DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
                  DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
                  DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))
               ENDIF

            ENDIF

            IF (STEST) THEN

               IF (RBAANORMALMODET) THEN
                  D2R1(J4,:) = MATMUL(D2RMI10,MATMUL(RMI,RBSITE(J2,:)))
                  D2R2(J4,:) = MATMUL(D2RMI20,MATMUL(RMI,RBSITE(J2,:)))
                  D2R3(J4,:) = MATMUL(D2RMI30,MATMUL(RMI,RBSITE(J2,:)))
                  
                  D2R12(J4,:) = MATMUL(D2RMI120,MATMUL(RMI,RBSITE(J2,:)))
                  D2R23(J4,:) = MATMUL(D2RMI230,MATMUL(RMI,RBSITE(J2,:)))
                  D2R31(J4,:) = MATMUL(D2RMI310,MATMUL(RMI,RBSITE(J2,:)))
               ELSE
                  D2R1(J4,:) = MATMUL(D2RMI1,RBSITE(J2,:))
                  D2R2(J4,:) = MATMUL(D2RMI2,RBSITE(J2,:))
                  D2R3(J4,:) = MATMUL(D2RMI3,RBSITE(J2,:))
                  
                  D2R12(J4,:) = MATMUL(D2RMI12,RBSITE(J2,:))
                  D2R23(J4,:) = MATMUL(D2RMI23,RBSITE(J2,:))
                  D2R31(J4,:) = MATMUL(D2RMI31,RBSITE(J2,:))
               ENDIF

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS 
         J3 = 3*J1
         J5 = OFFSET + J3

         DO I = 1, NRBSITES 
            J7 = NRBSITES*(J1-1) + I

            DO J2 = J1 + 1, REALNATOMS
            J4 = 3*J2
            J6 = OFFSET + J4

               DO J = 1, NRBSITES 
                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  DSS    = DSQRT(R2)
                  R2     = 1.D0/R2

                  IF ((I == 1) .AND. (J == 1))  THEN
                     RHO = MRHO11; REQ = REQ11; EPS = EPS11
                  ELSEIF ((I > 1) .AND. (J > 1)) THEN
                     RHO = MRHO22; REQ = REQ22; EPS = EPS22
                  ELSE
                     RHO = MRHO12; REQ = REQ12; EPS = EPS12
                  ENDIF
                  FCTR   = EXP(RHO*(REQ - DSS))
                  ENERGY = ENERGY + EPS*((1.D0 - FCTR)*(1.D0 - FCTR) - 1.D0)

                  IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R
                     DVDR(J7,J8) = EPS*2.D0*RHO*(-EXP(2.D0*RHO*(REQ - DSS)) + FCTR)/DSS
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

                     D2VDR2(J7,J8) =-DVDR(J7,J8)*R2 + 2.D0*EPS*RHO*RHO*R2*FCTR*(2.D0*FCTR -1.D0)
 
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

      IF (STEST) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4

               DO I = 1, NRBSITES

                  J7 = NRBSITES*(J1 - 1) + I

                  DO J = 1, NRBSITES 

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

      END SUBROUTINE SILANE 

!     ---------------------------------------------------------------------------------------------

      SUBROUTINE DEFSILANE()

      USE COMMONS, ONLY: RBSITE
      USE KEY, ONLY: EPS11, EPS22, EPS12, MRHO11, MRHO22, MRHO12, REQ11, REQ22, REQ12  

      IMPLICIT NONE
      DOUBLE PRECISION :: FCTR
!     SiH4 : Si - 1, H - 2

      FCTR      = 0.85391048D0
      RBSITE(1,:) = (/ 0.0D0, 0.0D0, 0.0D0/)
      RBSITE(2,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      RBSITE(3,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      RBSITE(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)
      RBSITE(5,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)

      EPS11  = 1.06D-6
      EPS12  = 0.042D0         ! kcal/mol
      EPS22  = 0.085D0
      MRHO11 = 1.352d0
      MRHO12 = 8.190d0         ! A-1
      MRHO22 = 1.216d0
      REQ11  = 8.743d0
      REQ12  = 2.836d0         ! A
      REQ22  = 3.348d0

      END SUBROUTINE DEFSILANE
