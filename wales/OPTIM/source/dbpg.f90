      SUBROUTINE DUMBBELLP (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
      USE KEY, ONLY: DBEPSBB, DBEPSAB, DBSIGBB, DBSIGAB, DBPMU, EFIELDT, EFIELD

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, DPFCT, DUMMY
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), RSS(3), NR(3), P(3), DU(3), EI(3), EJ(3), R(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DVDR(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), D2VDR2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: E(NATOMS/2,3), DE1(NATOMS/2,3), DE2(NATOMS/2,3), DE3(NATOMS/2,3)
      DOUBLE PRECISION :: D2E1(NATOMS/2,3), D2E2(NATOMS/2,3), D2E3(NATOMS/2,3)
      DOUBLE PRECISION :: D2E12(NATOMS/2,3), D2E23(NATOMS/2,3), D2E31(NATOMS/2,3)   
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3) 
      DOUBLE PRECISION :: D2R1(NRBSITES*NATOMS/2,3), D2R2(NRBSITES*NATOMS/2,3), D2R3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R12(NRBSITES*NATOMS/2,3), D2R23(NRBSITES*NATOMS/2,3), D2R31(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DOTI1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTI2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTI3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTJ1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTJ2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTJ3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3) 
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      DOUBLE PRECISION :: CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS
      DOUBLE PRECISION :: DVRDR, DVRDA, DVRDB, DVRDG, DVADR, DVBDR, DVGDR, FCT1, FCT2
      DOUBLE PRECISION :: DVADB, DVBDA, DADR(3), DBDR(3), D2ADX2(3), D2BDX2(3)
      DOUBLE PRECISION :: D2ADYX, D2ADZY, D2ADXZ, D2BDYX, D2BDZY, D2BDXZ
      DOUBLE PRECISION :: D2API1, D2API2, D2API3, D2BPJ1, D2BPJ2, D2BPJ3
      DOUBLE PRECISION :: D2GPI1, D2GPI2, D2GPI3, D2GPJ1, D2GPJ2, D2GPJ3
      LOGICAL          :: GTEST, STEST

! hk286
!      DOUBLE PRECISION :: k1, k2, DR, alpha
!      k1 = 400.0D0
!      k2 = 400.0D0
!      alpha = 2.75D0
! hk286

      CALL DEFDUM(DU, DPFCT, CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS)

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

         E(J1,:)   = MATMUL(RMI(:,:),DU(:))
         DE1(J1,:) = MATMUL(DRMI1(:,:),DU(:))
         DE2(J1,:) = MATMUL(DRMI2(:,:),DU(:))
         DE3(J1,:) = MATMUL(DRMI3(:,:),DU(:))

         IF (STEST) THEN

            D2E1(J1,:) = MATMUL(D2RMI1(:,:),DU(:))
            D2E2(J1,:) = MATMUL(D2RMI2(:,:),DU(:))
            D2E3(J1,:) = MATMUL(D2RMI3(:,:),DU(:))

            D2E12(J1,:) = MATMUL(D2RMI12(:,:),DU(:))
            D2E23(J1,:) = MATMUL(D2RMI23(:,:),DU(:))
            D2E31(J1,:) = MATMUL(D2RMI31(:,:),DU(:))

         ENDIF

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI(:,:),RBSITE(J2,:))

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

      IF (.NOT. STEST) THEN

         DO J1 = 1, REALNATOMS - 1 

            J3 = 3*J1
            J5 = OFFSET + J3
 
            RI(:)  = X(J3-2:J3)
            EI(:)  = E(J1,:)

            DO J2 = J1 + 1, REALNATOMS

               J4 = 3*J2
               J6 = OFFSET + J4

!     LJ CONTRIBUTION

               DO I = 1, NRBSITES - 1

                  J7 = NRBSITES*(J1-1) + I

                  DO J = 1, NRBSITES - 1

                     J8     = NRBSITES*(J2-1) + J
                     RSS(:) = R(J7,:) - R(J8,:)
                     R2     = DOT_PRODUCT(RSS(:),RSS(:))
                     R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
                     R6     = R2*R2*R2
                     R12    = R6*R6

                     IF (I == 1 .AND. J == 1) THEN
                        ENERGY = ENERGY + (CLJ12BB*R6 - CLJ6BB)*R6
                     ELSEIF (I == 2 .AND. J == 2) THEN
                        ENERGY = ENERGY + (CLJ12SS*R6 - CLJ6SS)*R6
                     ELSE
                        ENERGY = ENERGY + (CLJ12BS*R6 - CLJ6BS)*R6
                     ENDIF

                     IF (GTEST) THEN
!     DVDR = DVDR/R
                        IF (I == 1 .AND. J == 1) THEN
                           DVDR(J7,J8)   = -6.D0*(2.D0*CLJ12BB*R6 - CLJ6BB)*R6*R2
                        ELSEIF (I == 2 .AND. J == 2) THEN
                           DVDR(J7,J8)   = -6.D0*(2.D0*CLJ12SS*R6 - CLJ6SS)*R6*R2
                        ELSE
                           DVDR(J7,J8)   = -6.D0*(2.D0*CLJ12BS*R6 - CLJ6BS)*R6*R2
                        ENDIF

                        G(J3-2:J3)  = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
                        G(J4-2:J4)  = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)

                        G(J5-2)     = G(J5-2) + DVDR(J7,J8)*DOT_PRODUCT(RSS,DR1(J7,:))
                        G(J5-1)     = G(J5-1) + DVDR(J7,J8)*DOT_PRODUCT(RSS,DR2(J7,:))
                        G(J5)       = G(J5)   + DVDR(J7,J8)*DOT_PRODUCT(RSS,DR3(J7,:))

                        G(J6-2)     = G(J6-2) - DVDR(J7,J8)*DOT_PRODUCT(RSS,DR1(J8,:))
                        G(J6-1)     = G(J6-1) - DVDR(J7,J8)*DOT_PRODUCT(RSS,DR2(J8,:))
                        G(J6)       = G(J6)   - DVDR(J7,J8)*DOT_PRODUCT(RSS,DR3(J8,:))

                     ENDIF
                   
                  ENDDO

               ENDDO

!     DIPOLAR CONTRIBUTIONS

               RJ(:)  = X(J4-2:J4)
               RIJ(:) = RI(:) - RJ(:)
               RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
               R2     = 1.D0/RIJSQ
               ABSRIJ = DSQRT(RIJSQ)
               NR(:)  = RIJ(:)/ABSRIJ
               EJ(:)  = E(J2,:)
               ALP    = DOT_PRODUCT(NR(:),EI(:))
               BET    = DOT_PRODUCT(NR(:),EJ(:))
               GAM    = DOT_PRODUCT(EI(:),EJ(:))

               ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

               IF (GTEST) THEN

                  VR  = -DPFCT*R2*R2*(GAM - 3.D0*ALP*BET)
                  VA  = -DPFCT*BET*R2/ABSRIJ
                  VB  = -DPFCT*ALP*R2/ABSRIJ
                  VG  =  DPFCT*R2/(3.D0*ABSRIJ)

                  FIJN   = VR - (VA*ALP+VB*BET)/ABSRIJ
                  FIJEI  = VA/ABSRIJ
                  FIJEJ  = VB/ABSRIJ
                  FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

                  DADPI1 = DOT_PRODUCT(NR(:),DE1(J1,:))
                  DADPI2 = DOT_PRODUCT(NR(:),DE2(J1,:))
                  DADPI3 = DOT_PRODUCT(NR(:),DE3(J1,:))

                  DBDPJ1 = DOT_PRODUCT(NR(:),DE1(J2,:))
                  DBDPJ2 = DOT_PRODUCT(NR(:),DE2(J2,:))
                  DBDPJ3 = DOT_PRODUCT(NR(:),DE3(J2,:))

                  DGDPI1 = DOT_PRODUCT(DE1(J1,:),EJ(:))
                  DGDPI2 = DOT_PRODUCT(DE2(J1,:),EJ(:))
                  DGDPI3 = DOT_PRODUCT(DE3(J1,:),EJ(:))

                  DGDPJ1 = DOT_PRODUCT(EI(:),DE1(J2,:))
                  DGDPJ2 = DOT_PRODUCT(EI(:),DE2(J2,:))
                  DGDPJ3 = DOT_PRODUCT(EI(:),DE3(J2,:))

                  G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
                  G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

                  G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
                  G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
                  G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

                  G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
                  G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
                  G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

               ENDIF

            ENDDO

         ENDDO

      ELSE 

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3
 
            RI(:)  = X(J3-2:J3)
            EI(:)  = E(J1,:)

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4

!     LJ CONTRIBUTION

               DO I = 1, NRBSITES - 1

                  J7 = NRBSITES*(J1-1) + I

                  DO J = 1, NRBSITES - 1

                     J8     = NRBSITES*(J2-1) + J
                     RSS(:) = R(J7,:) - R(J8,:)
                     R2     = DOT_PRODUCT(RSS(:),RSS(:))
                     R2     = 1.D0/DOT_PRODUCT(RSS(:),RSS(:))
                     R6     = R2*R2*R2
                     R12    = R6*R6
!     DVDR = DVDR/R
                     IF (I == 1 .AND. J == 1) THEN
                        ENERGY = ENERGY + (CLJ12BB*R6 - CLJ6BB)*R6
                        DVDR(J7,J8)   = -6.D0*(2.D0*CLJ12BB*R6 - CLJ6BB)*R6*R2
                        D2VDR2(J7,J8) = 168.D0*CLJ12BB*R12*R2*R2 - 48.D0*CLJ6BB*R6*R2*R2 
                     ELSEIF (I == 2 .AND. J == 2) THEN
                        ENERGY = ENERGY + (CLJ12SS*R6 - CLJ6SS)*R6
                        DVDR(J7,J8)   = -6.D0*(2.D0*CLJ12SS*R6 - CLJ6SS)*R6*R2
                        D2VDR2(J7,J8) = 168.D0*CLJ12SS*R12*R2*R2 - 48.D0*CLJ6SS*R6*R2*R2
                     ELSE
                        ENERGY = ENERGY + (CLJ12BS*R6 - CLJ6BS)*R6
                        DVDR(J7,J8)   = -6.D0*(2.D0*CLJ12BS*R6 - CLJ6BS)*R6*R2
                        D2VDR2(J7,J8) = 168.D0*CLJ12BS*R12*R2*R2 - 48.D0*CLJ6BS*R6*R2*R2
                     ENDIF

                     DVDR(J8,J7)  = DVDR(J7,J8)
                     DOTI1(J7,J8) = DOT_PRODUCT(RSS,DR1(J7,:))
                     DOTI2(J7,J8) = DOT_PRODUCT(RSS,DR2(J7,:))
                     DOTI3(J7,J8) = DOT_PRODUCT(RSS,DR3(J7,:))

                     DOTJ1(J7,J8) = DOT_PRODUCT(RSS,DR1(J8,:))
                     DOTJ2(J7,J8) = DOT_PRODUCT(RSS,DR2(J8,:))
                     DOTJ3(J7,J8) = DOT_PRODUCT(RSS,DR3(J8,:))

                     G(J3-2:J3)  = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
                     G(J4-2:J4)  = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)

                     G(J5-2)     = G(J5-2) + DVDR(J7,J8)*DOTI1(J7,J8)
                     G(J5-1)     = G(J5-1) + DVDR(J7,J8)*DOTI2(J7,J8)
                     G(J5)       = G(J5)   + DVDR(J7,J8)*DOTI3(J7,J8)

                     G(J6-2)     = G(J6-2) - DVDR(J7,J8)*DOTJ1(J7,J8)
                     G(J6-1)     = G(J6-1) - DVDR(J7,J8)*DOTJ2(J7,J8)
                     G(J6)       = G(J6)   - DVDR(J7,J8)*DOTJ3(J7,J8)

                     D2VDR2(J8,J7) = D2VDR2(J7,J8)
                     DOTI1(J8,J7)  = -DOTJ1(J7,J8)
                     DOTI2(J8,J7)  = -DOTJ2(J7,J8)
                     DOTI3(J8,J7)  = -DOTJ3(J7,J8)
                     DOTJ1(J8,J7)  = -DOTI1(J7,J8)
                     DOTJ2(J8,J7)  = -DOTI2(J7,J8)
                     DOTJ3(J8,J7)  = -DOTI3(J7,J8)

                  ENDDO

               ENDDO

!     DIPOLAR CONTRIBUTIONS

               RJ(:)  = X(J4-2:J4)
               RIJ(:) = RI(:) - RJ(:)
               RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
               R2     = 1.D0/RIJSQ
               R4     = R2*R2
               ABSRIJ = DSQRT(RIJSQ)
               NR(:)  = RIJ(:)/ABSRIJ
               EJ(:)  = E(J2,:)
               ALP    = DOT_PRODUCT(NR(:),EI(:))
               BET    = DOT_PRODUCT(NR(:),EJ(:))
               GAM    = DOT_PRODUCT(EI(:),EJ(:))

               ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

               VR  = -DPFCT*R2*R2*(GAM - 3.D0*ALP*BET)
               VA  = -DPFCT*BET*R2/ABSRIJ
               VB  = -DPFCT*ALP*R2/ABSRIJ
               VG  =  DPFCT*R2/(3.D0*ABSRIJ)

               FIJN   = VR - (VA*ALP+VB*BET)/ABSRIJ
               FIJEI  = VA/ABSRIJ
               FIJEJ  = VB/ABSRIJ
               FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

               DADPI1 = DOT_PRODUCT(NR(:),DE1(J1,:))
               DADPI2 = DOT_PRODUCT(NR(:),DE2(J1,:))
               DADPI3 = DOT_PRODUCT(NR(:),DE3(J1,:))

               DBDPJ1 = DOT_PRODUCT(NR(:),DE1(J2,:))
               DBDPJ2 = DOT_PRODUCT(NR(:),DE2(J2,:))
               DBDPJ3 = DOT_PRODUCT(NR(:),DE3(J2,:))

               DGDPI1 = DOT_PRODUCT(DE1(J1,:),EJ(:))
               DGDPI2 = DOT_PRODUCT(DE2(J1,:),EJ(:))
               DGDPI3 = DOT_PRODUCT(DE3(J1,:),EJ(:))

               DGDPJ1 = DOT_PRODUCT(EI(:),DE1(J2,:))
               DGDPJ2 = DOT_PRODUCT(EI(:),DE2(J2,:))
               DGDPJ3 = DOT_PRODUCT(EI(:),DE3(J2,:))

               G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
               G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

!     HESSIAN CALCULATION

               DVRDR = 4.D0*DPFCT*R4*(GAM - 3.D0*ALP*BET)/ABSRIJ
               DVADR = 3.D0*DPFCT*BET*R4
               DVBDR = 3.D0*DPFCT*ALP*R4
               DVGDR = -DPFCT*R4
               DVADB = -DPFCT*R2/ABSRIJ
               DVRDA = DVADR
               DVRDB = DVBDR
               DVRDG = DVGDR
               DVBDA = DVADB

               DADR(:) = EI(:)/ABSRIJ - ALP*R2*RIJ(:)
               DBDR(:) = EJ(:)/ABSRIJ - BET*R2*RIJ(:)

               D2ADX2(1) = - 2.D0*R2*RIJ(1)*EI(1)/ABSRIJ + 3.D0*ALP*RIJ(1)*RIJ(1)*R4 - ALP*R2
               D2ADX2(2) = - 2.D0*R2*RIJ(2)*EI(2)/ABSRIJ + 3.D0*ALP*RIJ(2)*RIJ(2)*R4 - ALP*R2
               D2ADX2(3) = - 2.D0*R2*RIJ(3)*EI(3)/ABSRIJ + 3.D0*ALP*RIJ(3)*RIJ(3)*R4 - ALP*R2

               D2BDX2(1) = - 2.D0*R2*RIJ(1)*EJ(1)/ABSRIJ + 3.D0*BET*RIJ(1)*RIJ(1)*R4 - BET*R2
               D2BDX2(2) = - 2.D0*R2*RIJ(2)*EJ(2)/ABSRIJ + 3.D0*BET*RIJ(2)*RIJ(2)*R4 - BET*R2
               D2BDX2(3) = - 2.D0*R2*RIJ(3)*EJ(3)/ABSRIJ + 3.D0*BET*RIJ(3)*RIJ(3)*R4 - BET*R2

               D2ADYX    = - R2*(EI(1)*RIJ(2)+EI(2)*RIJ(1))/ABSRIJ + 3*ALP*RIJ(1)*RIJ(2)*R4
               D2ADZY    = - R2*(EI(3)*RIJ(2)+EI(2)*RIJ(3))/ABSRIJ + 3*ALP*RIJ(3)*RIJ(2)*R4
               D2ADXZ    = - R2*(EI(3)*RIJ(1)+EI(1)*RIJ(3))/ABSRIJ + 3*ALP*RIJ(1)*RIJ(3)*R4

               D2BDYX    = - R2*(EJ(1)*RIJ(2)+EJ(2)*RIJ(1))/ABSRIJ + 3*BET*RIJ(1)*RIJ(2)*R4
               D2BDZY    = - R2*(EJ(3)*RIJ(2)+EJ(2)*RIJ(3))/ABSRIJ + 3*BET*RIJ(3)*RIJ(2)*R4
               D2BDXZ    = - R2*(EJ(3)*RIJ(1)+EJ(1)*RIJ(3))/ABSRIJ + 3*BET*RIJ(1)*RIJ(3)*R4

!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES
!     xi,xi
               HESS(J3-2,J3-2) = HESS(J3-2,J3-2) &
                               + DVRDR*R2*RIJ(1)*RIJ(1) + (DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(1)              &
                               + VR*(1.D0 - R2*RIJ(1)*RIJ(1))/ABSRIJ + (DVADR*NR(1) + DVADB*DBDR(1))*DADR(1) &
                               + VA*D2ADX2(1) + (DVBDR*NR(1) + DVBDA*DADR(1) )*DBDR(1)  + VB*D2BDX2(1)
!     yi,yi
               HESS(J3-1,J3-1) = HESS(J3-1,J3-1) &
                               + DVRDR*R2*RIJ(2)*RIJ(2) + (DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(2)              &
                               + VR*(1.D0 - R2*RIJ(2)*RIJ(2))/ABSRIJ + (DVADR*NR(2) + DVADB*DBDR(2))*DADR(2) &
                               + VA*D2ADX2(2) + (DVBDR*NR(2) + DVBDA*DADR(2))*DBDR(2) + VB*D2BDX2(2)
!     zi,zi
               HESS(J3,J3)     = HESS(J3,J3)     &
                               + DVRDR*R2*RIJ(3)*RIJ(3) + (DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(3)              &
                               + VR*(1.D0 - R2*RIJ(3)*RIJ(3))/ABSRIJ+(DVADR*NR(3) + DVADB*DBDR(3))*DADR(3)   &
                               + VA*D2ADX2(3) + (DVBDR*NR(3) + DVBDA*DADR(3))*DBDR(3) + VB*D2BDX2(3)
!     pi1,pi1
               HESS(J5-2,J5-2) = HESS(J5-2,J5-2) + VA*DOT_PRODUCT(NR(:),D2E1(J1,:)) &
                               + VG*DOT_PRODUCT(D2E1(J1,:),EJ(:))
!     pi2,pi2
               HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + VA*DOT_PRODUCT(NR(:),D2E2(J1,:)) &
                               + VG*DOT_PRODUCT(D2E2(J1,:),EJ(:))
!     pi3,pi3
               HESS(J5,J5)     = HESS(J5,J5)     + VA*DOT_PRODUCT(NR(:),D2E3(J1,:)) &
                               + VG*DOT_PRODUCT(D2E3(J1,:),EJ(:))

!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCK: SAME MOLECULE, DIFFERENT COORDINATES

!     xi,yi
               DUMMY = (DVRDR*NR(2) + DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(1)  &
                     - VR*R2*NR(1)*RIJ(2) + (DVADR*NR(2) + DVADB*DBDR(2))*DADR(1)  &
                     + VA*D2ADYX + (DVBDR*NR(2) + DVBDA*DADR(2))*DBDR(1) + VB*D2BDYX
               HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
               HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
               DUMMY = (DVRDR*NR(3) + DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(2)  &
                     - VR*R2*NR(2)*RIJ(3) + (DVADR*NR(3) + DVADB*DBDR(3))*DADR(2)  &
                     + VA*D2ADZY + (DVBDR*NR(3) + DVBDA*DADR(3))*DBDR(2) + VB*D2BDZY
               HESS(J3-1,J3) = HESS(J3-1,J3) + DUMMY
               HESS(J3,J3-1) = HESS(J3,J3-1) + DUMMY
!     zi,xi
               DUMMY = (DVRDR*NR(1) + DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(3)  &
                     - VR*R2*NR(3)*RIJ(1) + (DVADR*NR(1) + DVADB*DBDR(1))*DADR(3)  &
                     + VA*D2ADXZ + (DVBDR*NR(1) + DVBDA*DADR(1))*DBDR(3) + VB*D2BDXZ
               HESS(J3-2,J3) = HESS(J3-2,J3) + DUMMY
               HESS(J3,J3-2) = HESS(J3,J3-2) + DUMMY

               FCT1  = DVRDA*DADPI1 + DVRDG*DGDPI1 - (VA*DADPI1 + DVBDA*DADPI1*BET)/ABSRIJ
               FCT2  = DVBDA*DADPI1/ABSRIJ
!     xi,pi1
               DUMMY = FCT1*NR(1) + VA*DE1(J1,1)/ABSRIJ + FCT2*EJ(1)
               HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
               HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1
               DUMMY = FCT1*NR(2) + VA*DE1(J1,2)/ABSRIJ + FCT2*EJ(2)
               HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
               HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
               DUMMY = FCT1*NR(3) + VA*DE1(J1,3)/ABSRIJ + FCT2*EJ(3)
               HESS(J3,J5-2)   = HESS(J3,J5-2)   + DUMMY
               HESS(J5-2,J3)   = HESS(J5-2,J3)   + DUMMY

               FCT1  = DVRDA*DADPI2 + DVRDG*DGDPI2 - (VA*DADPI2 + DVBDA*DADPI2*BET)/ABSRIJ
               FCT2  = DVBDA*DADPI2/ABSRIJ
!     xi,pi2
               DUMMY = FCT1*NR(1) + VA*DE2(J1,1)/ABSRIJ + FCT2*EJ(1)
               HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
               HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2
               DUMMY = FCT1*NR(2) + VA*DE2(J1,2)/ABSRIJ + FCT2*EJ(2)
               HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
               HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
               DUMMY = FCT1*NR(3) + VA*DE2(J1,3)/ABSRIJ + FCT2*EJ(3)
               HESS(J3,J5-1)   = HESS(J3,J5-1)   + DUMMY
               HESS(J5-1,J3)   = HESS(J5-1,J3)   + DUMMY

               FCT1  = DVRDA*DADPI3 + DVRDG*DGDPI3 - (VA*DADPI3 + DVBDA*DADPI3*BET)/ABSRIJ
               FCT2  = DVBDA*DADPI3/ABSRIJ
!     xi,pi2
               DUMMY = FCT1*NR(1) + VA*DE3(J1,1)/ABSRIJ + FCT2*EJ(1)
               HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
               HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi2
               DUMMY = FCT1*NR(2) + VA*DE3(J1,2)/ABSRIJ + FCT2*EJ(2)
               HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
               HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi2
               DUMMY = FCT1*NR(3) + VA*DE3(J1,3)/ABSRIJ + FCT2*EJ(3)
               HESS(J3,J5)     = HESS(J3,J5)   + DUMMY
               HESS(J5,J3)   = HESS(J5,J3)   + DUMMY
!     pi1,pi2
               DUMMY = VA*DOT_PRODUCT(NR(:),D2E12(J1,:)) + VG*DOT_PRODUCT(D2E12(J1,:),EJ(:))
               HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
               HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
               DUMMY = VA*DOT_PRODUCT(NR(:),D2E23(J1,:)) + VG*DOT_PRODUCT(D2E23(J1,:),EJ(:))
               HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
               HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
               DUMMY = VA*DOT_PRODUCT(NR(:),D2E31(J1,:)) + VG*DOT_PRODUCT(D2E31(J1,:),EJ(:))
               HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
               HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE

!     xi,xj
               HESS(J3-2,J4-2) = - DVRDR*R2*RIJ(1)*RIJ(1) - (DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(1)      &
                               - VR*(1.D0 - R2*RIJ(1)*RIJ(1))/ABSRIJ - (DVADR*NR(1) + DVADB*DBDR(1)) &
                               *DADR(1) - VA*D2ADX2(1) - (DVBDR*NR(1) + DVBDA*DADR(1))*DBDR(1) - VB*D2BDX2(1)
!     yi,yj
               HESS(J3-1,J4-1) = - DVRDR*R2*RIJ(2)*RIJ(2) - (DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(2)      &
                               - VR*(1.D0 - R2*RIJ(2)*RIJ(2))/ABSRIJ - (DVADR*NR(2) + DVADB*DBDR(2)) &
                               *DADR(2) - VA*D2ADX2(2) - (DVBDR*NR(2) + DVBDA*DADR(2))*DBDR(2) - VB*D2BDX2(2)
!     zi,zj
               HESS(J3,J4)     = - DVRDR*R2*RIJ(3)*RIJ(3) - (DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(3)      &
                               - VR*(1.D0 - R2*RIJ(3)*RIJ(3))/ABSRIJ - (DVADR*NR(3) + DVADB*DBDR(3)) &
                               *DADR(3) - VA*D2ADX2(3) - (DVBDR*NR(3) + DVBDA*DADR(3))*DBDR(3) - VB*D2BDX2(3)
!     pi1,pj1
               HESS(J5-2,J6-2) = HESS(J5-2,J6-2) + DVADB*DBDPJ1*DADPI1 + VG*DOT_PRODUCT(DE1(J1,:),DE1(J2,:))
!     pi2,pj2
               HESS(J5-1,J6-1) = HESS(J5-1,J6-1) + DVADB*DBDPJ2*DADPI2 + VG*DOT_PRODUCT(DE2(J1,:),DE2(J2,:))
!     pi3,pj3
               HESS(J5,J6)     = HESS(J5,J6)     + DVADB*DBDPJ3*DADPI3 + VG*DOT_PRODUCT(DE3(J1,:),DE3(J2,:))

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES

!     xi,yj
               HESS(J3-2,J4-1) = -(DVRDR*NR(2) + DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(1) &
                               + VR*R2*NR(1)*RIJ(2) - (DVADR*NR(2) + DVADB*DBDR(2))*DADR(1)  &
                               - VA*D2ADYX - (DVBDR*NR(2) + DVBDA*DADR(2))*DBDR(1) - VB*D2BDYX
               HESS(J3-1,J4-2) = HESS(J3-2,J4-1)
!     yi,zj
               HESS(J3-1,J4)   = -(DVRDR*NR(3) + DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(2)  &
                               + VR*R2*NR(2)*RIJ(3) - (DVADR*NR(3) + DVADB*DBDR(3))*DADR(2)  &
                               - VA*D2ADZY - (DVBDR*NR(3) + DVBDA*DADR(3))*DBDR(2) - VB*D2BDZY
               HESS(J3,J4-1)   = HESS(J3-1,J4)
!     xi,zj
               HESS(J3-2,J4)   = -(DVRDR*NR(1) + DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(3)  &
                               + VR*R2*NR(3)*RIJ(1) - (DVADR*NR(1) + DVADB*DBDR(1))*DADR(3)  &
                               - VA*D2ADXZ - (DVBDR*NR(1) + DVBDA*DADR(1))*DBDR(3) - VB*D2BDXZ
               HESS(J3,J4-2)   = HESS(J3-2,J4)

               FCT1 = DVRDB*DBDPJ1 + DVRDG*DGDPJ1 - (DVADB*DBDPJ1*ALP + VB*DBDPJ1)/ABSRIJ 
               FCT2 = DVADB*DBDPJ1/ABSRIJ
!     xi,pj1
               DUMMY = FCT1*NR(1) + FCT2*EI(1) + VB*DE1(J2,1)/ABSRIJ
               HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
               HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
               DUMMY = FCT1*NR(2) + FCT2*EI(2) + VB*DE1(J2,2)/ABSRIJ
               HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
               HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
               DUMMY = FCT1*NR(3) + FCT2*EI(3) + VB*DE1(J2,3)/ABSRIJ
               HESS(J3,J6-2)   = HESS(J3,J6-2)   + DUMMY
               HESS(J6-2,J3)   = HESS(J6-2,J3)   + DUMMY

               FCT1  = DVRDB*DBDPJ2 + DVRDG*DGDPJ2 - (DVADB*DBDPJ2*ALP + VB*DBDPJ2)/ABSRIJ
               FCT2  = DVADB*DBDPJ2/ABSRIJ
!     xi,pj2
               DUMMY = FCT1*NR(1) + FCT2*EI(1) + VB*DE2(J2,1)/ABSRIJ
               HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
               HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
               DUMMY = FCT1*NR(2) + FCT2*EI(2) + VB*DE2(J2,2)/ABSRIJ
               HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
               HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
               DUMMY = FCT1*NR(3) + FCT2*EI(3) + VB*DE2(J2,3)/ABSRIJ
               HESS(J3,J6-1)   = HESS(J3,J6-1)   + DUMMY
               HESS(J6-1,J3)   = HESS(J6-1,J3)   + DUMMY

               FCT1  = DVRDB*DBDPJ3 + DVRDG*DGDPJ3 - (DVADB*DBDPJ3*ALP + VB*DBDPJ3)/ABSRIJ
               FCT2  = DVADB*DBDPJ3/ABSRIJ
!     xi,pj3
               DUMMY = FCT1*NR(1) + FCT2*EI(1) + VB*DE3(J2,1)/ABSRIJ
               HESS(J3-2,J6)   = HESS(J3-2,J6)   + DUMMY
               HESS(J6,J3-2)   = HESS(J6,J3-2)   + DUMMY
!     yi,pj3
               DUMMY = FCT1*NR(2) + FCT2*EI(2) + VB*DE3(J2,2)/ABSRIJ
               HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
               HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
               DUMMY = FCT1*NR(3) + FCT2*EI(3) + VB*DE3(J2,3)/ABSRIJ
               HESS(J3,J6)     = HESS(J3,J6)   + DUMMY
               HESS(J6,J3)     = HESS(J6,J3)   + DUMMY
!     pi1,pj2
               HESS(J5-2,J6-1) = HESS(J5-2,J6-1) + DVADB*DBDPJ2*DADPI1 + VG*DOT_PRODUCT(DE1(J1,:),DE2(J2,:))
               HESS(J6-1,J5-2) = HESS(J5-2,J6-1)
!     pi2,pj3
               HESS(J5-1,J6)   = HESS(J5-1,J6) + DVADB*DBDPJ3*DADPI2 + VG*DOT_PRODUCT(DE2(J1,:),DE3(J2,:))
               HESS(J6,J5-1)   = HESS(J5-1,J6)
!     pi3,pj1
               HESS(J5,J6-2)   = HESS(J5,J6-2) + DVADB*DBDPJ1*DADPI3 + VG*DOT_PRODUCT(DE3(J1,:),DE1(J2,:))
               HESS(J6-2,J5)   = HESS(J5,J6-2)

            ENDDO

         ENDDO

      ENERGY = 0.5D0*ENERGY
      G(:)   = 0.5D0*G(:)

      ENDIF


      IF (STEST) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4

               DO I = 1, NRBSITES - 1

                  J7 = NRBSITES*(J1 - 1) + I

                  DO J = 1, NRBSITES - 1

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

      IF (EFIELDT) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3
            EI(:)  = E(J1,:)

            ENERGY = ENERGY - DBPMU*EFIELD*EI(3)

            IF (GTEST) THEN

               G(J5-2) = G(J5-2) - DBPMU*EFIELD*DE1(J1,3)
               G(J5-1) = G(J5-1) - DBPMU*EFIELD*DE2(J1,3)
               G(J5)   = G(J5)   - DBPMU*EFIELD*DE3(J1,3)

            ENDIF

            IF (STEST) THEN
!     pi1,pi1
               HESS(J5-2,J5-2) = HESS(J5-2,J5-2) - DBPMU*EFIELD*D2E1(J1,3)
!     pi2,pi2  
               HESS(J5-1,J5-1) = HESS(J5-1,J5-1) - DBPMU*EFIELD*D2E2(J1,3)
!     pi3,pi3 
               HESS(J5,J5)     = HESS(J5,J5)     - DBPMU*EFIELD*D2E3(J1,3) 
!     pi1,pi2
               DUMMY           = -DBPMU*EFIELD*D2E12(J1,3)
               HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
               HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
               DUMMY           = -DBPMU*EFIELD*D2E23(J1,3)
               HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
               HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
               DUMMY           = -DBPMU*EFIELD*D2E31(J1,3)
               HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
               HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

            ENDIF

         ENDDO

      ENDIF

      
! hk286 - 25 is around the limit, 50 seems to be ok. try also a bit lower and higher - for freezing
!      J1 = 1
!      J2 = 2
!      DR = (X(3*J1-2) - X(3*J2-2))**2 + (X(3*J1-1) - X(3*J2-1))**2 + (X(3*J1) - X(3*J2))**2 
!      DR = DSQRT(DR)
!      IF (GTEST) THEN       
!         G(3*J1)   = G(3*J1)   + (DR - alpha)/DR * (X(3*J1)   - X(3*J2)  ) * k1      
!         G(3*J1-1) = G(3*J1-1) + (DR - alpha)/DR * (X(3*J1-1) - X(3*J2-1)) * k1
!         G(3*J1-2) = G(3*J1-2) + (DR - alpha)/DR * (X(3*J1-2) - X(3*J2-2)) * k1

!         G(3*J2)   = G(3*J2)   - (DR - alpha)/DR * (X(3*J1)   - X(3*J2)  ) * k1      
!         G(3*J2-1) = G(3*J2-1) - (DR - alpha)/DR * (X(3*J1-1) - X(3*J2-1)) * k1
!         G(3*J2-2) = G(3*J2-2) - (DR - alpha)/DR * (X(3*J1-2) - X(3*J2-2)) * k1
!      ENDIF
!      ENERGY = ENERGY + 0.5D0* (DR - alpha)**2 * k1

!      IF (STEST) THEN
!         HESS(3*J1  ,3*J1  ) = HESS(3*J1  ,3*J1  ) + (1.0D0 - alpha/DR + alpha/(DR**3) * (X(3*J1)   - X(3*J2  ))**2 ) * k1
!         HESS(3*J1-1,3*J1-1) = HESS(3*J1-1,3*J1-1) + (1.0D0 - alpha/DR + alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1))**2 ) * k1
!         HESS(3*J1-2,3*J1-2) = HESS(3*J1-2,3*J1-2) + (1.0D0 - alpha/DR + alpha/(DR**3) * (X(3*J1-2) - X(3*J2-2))**2 ) * k1
!         HESS(3*J1-1,3*J1-2) = HESS(3*J1-1,3*J1-2) + (alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1)) * (X(3*J1-2) - X(3*J2-2)) ) * k1
!         HESS(3*J1-2,3*J1-1) = HESS(3*J1-2,3*J1-1) + (alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1)) * (X(3*J1-2) - X(3*J2-2)) ) * k1
!         HESS(3*J1-1,3*J1  ) = HESS(3*J1-1,3*J1  ) + (alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1)) * (X(3*J1  ) - X(3*J2  )) ) * k1
!         HESS(3*J1  ,3*J1-1) = HESS(3*J1  ,3*J1-1) + (alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1)) * (X(3*J1  ) - X(3*J2  )) ) * k1
!         HESS(3*J1  ,3*J1-2) = HESS(3*J1  ,3*J1-2) + (alpha/(DR**3) * (X(3*J1  ) - X(3*J2  )) * (X(3*J1-2) - X(3*J2-2)) ) * k1
!         HESS(3*J1-2,3*J1  ) = HESS(3*J1-2,3*J1  ) + (alpha/(DR**3) * (X(3*J1  ) - X(3*J2  )) * (X(3*J1-2) - X(3*J2-2)) ) * k1

!         HESS(3*J2  ,3*J2  ) = HESS(3*J2  ,3*J2  ) + (1.0D0 - alpha/DR + alpha/(DR**3) * (X(3*J2)   - X(3*J1  ))**2 ) * k1
!         HESS(3*J2-1,3*J2-1) = HESS(3*J2-1,3*J2-1) + (1.0D0 - alpha/DR + alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1))**2 ) * k1
!         HESS(3*J2-2,3*J2-2) = HESS(3*J2-2,3*J2-2) + (1.0D0 - alpha/DR + alpha/(DR**3) * (X(3*J2-2) - X(3*J1-2))**2 ) * k1
!         HESS(3*J2-1,3*J2-2) = HESS(3*J2-1,3*J2-2) + (alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1)) * (X(3*J2-2) - X(3*J1-2)) ) * k1
!         HESS(3*J2-2,3*J2-1) = HESS(3*J2-2,3*J2-1) + (alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1)) * (X(3*J2-2) - X(3*J1-2)) ) * k1
!         HESS(3*J2-1,3*J2  ) = HESS(3*J2-1,3*J2  ) + (alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1)) * (X(3*J2  ) - X(3*J1  )) ) * k1
!         HESS(3*J2  ,3*J2-1) = HESS(3*J2  ,3*J2-1) + (alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1)) * (X(3*J2  ) - X(3*J1  )) ) * k1
!         HESS(3*J2  ,3*J2-2) = HESS(3*J2  ,3*J2-2) + (alpha/(DR**3) * (X(3*J2  ) - X(3*J1  )) * (X(3*J2-2) - X(3*J1-2)) ) * k1
!         HESS(3*J2-2,3*J2  ) = HESS(3*J2-2,3*J2  ) + (alpha/(DR**3) * (X(3*J2  ) - X(3*J1  )) * (X(3*J2-2) - X(3*J1-2)) ) * k1

!         HESS(3*J1  ,3*J2  ) = HESS(3*J1  ,3*J2  ) + (-1.0D0 + alpha/DR - alpha/(DR**3) * (X(3*J1)   - X(3*J2  ))**2 ) * k1
!         HESS(3*J1-1,3*J2-1) = HESS(3*J1-1,3*J2-1) + (-1.0D0 + alpha/DR - alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1))**2 ) * k1
!         HESS(3*J1-2,3*J2-2) = HESS(3*J1-2,3*J2-2) + (-1.0D0 + alpha/DR - alpha/(DR**3) * (X(3*J1-2) - X(3*J2-2))**2 ) * k1
!         HESS(3*J1-1,3*J2-2) = HESS(3*J1-1,3*J2-2) - (alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1)) * (X(3*J1-2) - X(3*J2-2)) ) * k1
!         HESS(3*J1-2,3*J2-1) = HESS(3*J1-2,3*J2-1) - (alpha/(DR**3) * (X(3*J1-2) - X(3*J2-2)) * (X(3*J1-1) - X(3*J2-1)) ) * k1
!         HESS(3*J1-1,3*J2  ) = HESS(3*J1-1,3*J2  ) - (alpha/(DR**3) * (X(3*J1-1) - X(3*J2-1)) * (X(3*J1  ) - X(3*J2  )) ) * k1
!         HESS(3*J1  ,3*J2-1) = HESS(3*J1  ,3*J2-1) - (alpha/(DR**3) * (X(3*J1  ) - X(3*J2  )) * (X(3*J1-1) - X(3*J2-1)) ) * k1
!         HESS(3*J1  ,3*J2-2) = HESS(3*J1  ,3*J2-2) - (alpha/(DR**3) * (X(3*J1  ) - X(3*J2  )) * (X(3*J1-2) - X(3*J2-2)) ) * k1
!         HESS(3*J1-2,3*J2  ) = HESS(3*J1-2,3*J2  ) - (alpha/(DR**3) * (X(3*J1-2) - X(3*J2-2)) * (X(3*J1  ) - X(3*J2  )) ) * k1

!         HESS(3*J2  ,3*J1  ) = HESS(3*J2  ,3*J1  ) + (-1.0D0 + alpha/DR - alpha/(DR**3) * (X(3*J2)   - X(3*J1  ))**2 ) * k1
!         HESS(3*J2-1,3*J1-1) = HESS(3*J2-1,3*J1-1) + (-1.0D0 + alpha/DR - alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1))**2 ) * k1
!         HESS(3*J2-2,3*J1-2) = HESS(3*J2-2,3*J1-2) + (-1.0D0 + alpha/DR - alpha/(DR**3) * (X(3*J2-2) - X(3*J1-2))**2 ) * k1
!         HESS(3*J2-1,3*J1-2) = HESS(3*J2-1,3*J1-2) - (alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1)) * (X(3*J2-2) - X(3*J1-2)) ) * k1
!         HESS(3*J2-2,3*J1-1) = HESS(3*J2-2,3*J1-1) - (alpha/(DR**3) * (X(3*J2-2) - X(3*J1-2)) * (X(3*J2-1) - X(3*J1-1)) ) * k1
!         HESS(3*J2-1,3*J1  ) = HESS(3*J2-1,3*J1  ) - (alpha/(DR**3) * (X(3*J2-1) - X(3*J1-1)) * (X(3*J2  ) - X(3*J1  )) ) * k1
!         HESS(3*J2  ,3*J1-1) = HESS(3*J2  ,3*J1-1) - (alpha/(DR**3) * (X(3*J2  ) - X(3*J1  )) * (X(3*J2-1) - X(3*J1-1)) ) * k1
!         HESS(3*J2  ,3*J1-2) = HESS(3*J2  ,3*J1-2) - (alpha/(DR**3) * (X(3*J2  ) - X(3*J1  )) * (X(3*J2-2) - X(3*J1-2)) ) * k1
!         HESS(3*J2-2,3*J1  ) = HESS(3*J2-2,3*J1  ) - (alpha/(DR**3) * (X(3*J2-2) - X(3*J1-2)) * (X(3*J2  ) - X(3*J1  )) ) * k1

!      ENDIF


      END SUBROUTINE DUMBBELLP 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDUM(DU, DPFCT, CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS)
    
      USE COMMONS, ONLY: RBSITE
      USE KEY, ONLY: DBEPSBB, DBEPSAB, DBSIGBB, DBSIGAB, DBPMU  
     
      IMPLICIT NONE

      DOUBLE PRECISION :: DBSIGAA, DBEPSAA, DPFCT, DU(3)
      DOUBLE PRECISION :: CLJ6BB, CLJ12BB, CLJ6SS, CLJ12SS, CLJ6BS, CLJ12BS

      DBSIGAA  = 1.D0
      DBEPSAA  = 1.D0

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.5D0*DBSIGAA

      RBSITE(2,1) = 0.D0
      RBSITE(2,2) = 0.D0
      RBSITE(2,3) = -0.5D0*DBSIGBB

      RBSITE(3,1) = 0.D0
      RBSITE(3,2) = 0.D0
      RBSITE(3,3) = 0.D0

      DU(:) = (/0.D0, 1.D0, 0.D0/)
      DPFCT = 3.D0*DBPMU*DBPMU*DBSIGAA**3.D0

      CLJ6BB  = 4.D0*DBEPSAA*DBSIGAA**6 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      CLJ12BB = 4.D0*DBEPSAA*DBSIGAA**12
      CLJ6SS  = 4.D0*DBEPSBB*DBSIGBB**6
      CLJ12SS = 4.D0*DBEPSBB*DBSIGBB**12
      CLJ6BS  = 4.D0*DBEPSAB*DBSIGAB**6
      CLJ12BS = 4.D0*DBEPSAB*DBSIGAB**12

      END SUBROUTINE DEFDUM
