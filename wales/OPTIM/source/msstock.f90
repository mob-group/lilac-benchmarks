      SUBROUTINE MSSTOCKGH (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, RBSTLA, DPMU
      USE KEY, ONLY: NTSITES, EFIELDT, EFIELD


      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET, J2INT
!      INTEGER, PARAMETER :: NTSITES = NRBSITES*NATOMS/2
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, DUMMY, DVDR, D2VDR2, FCTR, DPFCT
      DOUBLE PRECISION :: RI(3), RJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3), UVX(3), UVY(3), UVZ(3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG
      DOUBLE PRECISION :: FIJN, FIJEI, FIJEJ, FIJ(3)
      DOUBLE PRECISION :: R(NTSITES,3), E(NTSITES,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
      DOUBLE PRECISION :: DR1(NTSITES,3), DR2(NTSITES,3), DR3(NTSITES,3) 
      DOUBLE PRECISION :: D2R1(NTSITES,3), D2R2(NTSITES,3), D2R3(NTSITES,3) 
      DOUBLE PRECISION :: D2R12(NTSITES,3), D2R23(NTSITES,3), D2R31(NTSITES,3) 
      DOUBLE PRECISION :: DE1(NTSITES,3), DE2(NTSITES,3), DE3(NTSITES,3)
      DOUBLE PRECISION :: D2E1(NTSITES,3), D2E2(NTSITES,3), D2E3(NTSITES,3)
      DOUBLE PRECISION :: D2E12(NTSITES,3), D2E23(NTSITES,3), D2E31(NTSITES,3)
      DOUBLE PRECISION :: DRIJDI1, DRIJDI2, DRIJDI3, DRIJDJ1, DRIJDJ2, DRIJDJ3
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPI1, DBDPI2, DBDPI3
      DOUBLE PRECISION :: DADPJ1, DADPJ2, DADPJ3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      DOUBLE PRECISION :: DVRDR, DVRDA, DVRDB, DVRDG, DVADR, DVBDR, DVGDR, DVADB, DVBDA
      DOUBLE PRECISION :: DADR(3), DBDR(3), D2ADRX(3), D2ADRY(3), D2ADRZ(3), D2BDRX(3), D2BDRY(3), D2BDRZ(3)
      DOUBLE PRECISION :: D2ADRIX(3), D2ADRIY(3), D2ADRIZ(3), D2BDRIX(3), D2BDRIY(3), D2BDRIZ(3)
      DOUBLE PRECISION :: D2ADRJX(3), D2ADRJY(3), D2ADRJZ(3), D2BDRJX(3), D2BDRJY(3), D2BDRJZ(3)
      DOUBLE PRECISION :: D2RDIXIX, D2RDIXJX, D2ADIXIX, D2ADIXJX, D2BDIXIX, D2BDIXJX
      DOUBLE PRECISION :: D2RDIYIY, D2RDIYJY, D2ADIYIY, D2ADIYJY, D2BDIYIY, D2BDIYJY
      DOUBLE PRECISION :: D2RDIZIZ, D2RDIZJZ, D2ADIZIZ, D2ADIZJZ, D2BDIZIZ, D2BDIZJZ
      DOUBLE PRECISION :: D2RDIXIY, D2RDIYIZ, D2RDIZIX, D2RDIXJY, D2RDIYJZ, D2RDIZJX
      DOUBLE PRECISION :: D2ADIXIY, D2ADIYIZ, D2ADIZIX, D2ADIXJY, D2ADIYJZ, D2ADIZJX
      DOUBLE PRECISION :: D2BDIXIY, D2BDIYIZ, D2BDIZIX, D2BDIXJY, D2BDIYJZ, D2BDIZJX
      LOGICAL          :: GTEST, STEST

      UVX = (/1.D0, 0.D0, 0.D0/)
      UVY = (/0.D0, 1.D0, 0.D0/)
      UVZ = (/0.D0, 0.D0, 1.D0/)
      ENERGY  = 0.D0
      J2INT = 1
      IF (GTEST) G(:) = 0.D0
      IF (STEST) THEN
         HESS(:,:) = 0.D0 
         FCTR = 0.5D0
      ELSE
         FCTR = 1.D0
      ENDIF

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         P  = X(J5-2:J5)

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, STEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,RBSITE(J2,:))
            E(J4,:)   = MATMUL(RMI(:,:),RBSTLA(J2,:))

            IF (GTEST .OR. STEST) THEN
 
               DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))
               DE1(J4,:) = MATMUL(DRMI1(:,:),RBSTLA(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBSTLA(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBSTLA(J2,:))

            ENDIF

            IF (STEST) THEN

               D2R1(J4,:)  = MATMUL(D2RMI1(:,:),RBSITE(J2,:))
               D2R2(J4,:)  = MATMUL(D2RMI2(:,:),RBSITE(J2,:))
               D2R3(J4,:)  = MATMUL(D2RMI3(:,:),RBSITE(J2,:))

               D2R12(J4,:) = MATMUL(D2RMI12(:,:),RBSITE(J2,:))
               D2R23(J4,:) = MATMUL(D2RMI23(:,:),RBSITE(J2,:))
               D2R31(J4,:) = MATMUL(D2RMI31(:,:),RBSITE(J2,:))

               D2E1(J4,:) = MATMUL(D2RMI1(:,:),RBSTLA(J2,:))
               D2E2(J4,:) = MATMUL(D2RMI2(:,:),RBSTLA(J2,:))
               D2E3(J4,:) = MATMUL(D2RMI3(:,:),RBSTLA(J2,:))

               D2E12(J4,:) = MATMUL(D2RMI12(:,:),RBSTLA(J2,:))
               D2E23(J4,:) = MATMUL(D2RMI23(:,:),RBSTLA(J2,:))
               D2E31(J4,:) = MATMUL(D2RMI31(:,:),RBSTLA(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS ! - 1  

         J3 = 3*J1
         J5 = OFFSET + J3

         DO I = 1, NRBSITES

            J7 = NRBSITES*(J1-1) + I
            EI(:) = E(J7,:)

            IF (.NOT. STEST) J2INT = J1 + 1

            DO J2 = J2INT, REALNATOMS

               IF (J1 == J2) CYCLE
               J4 = 3*J2
               J6 = OFFSET + J4

               DO J = 1, NRBSITES 

                  J8     = NRBSITES*(J2-1) + J
                  EJ(:)  = E(J8,:)
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
!     LJ CONTRIBUTION
                  R2     = 1.D0/R2
                  R6     = R2*R2*R2
                  R12    = R6*R6
                  ENERGY = ENERGY + 4.D0*(R12 - R6)
!     DIPOLAR CONTRIBUTION
                  R4     = R2*R2
                  ALP    = DOT_PRODUCT(NR(:),EI(:))
                  BET    = DOT_PRODUCT(NR(:),EJ(:))
                  GAM    = DOT_PRODUCT(EI(:),EJ(:))
                  DPFCT  = 3.D0*DPMU(I)*DPMU(J)
                  ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

                  IF (GTEST .OR. STEST) THEN
!     LJ CONTRIBUTION
                     DVDR = 4.D0*(-12.D0*R12 + 6.D0*R6)*R2

                     DRIJDI1 = DOT_PRODUCT(RSS,DR1(J7,:))
                     DRIJDI2 = DOT_PRODUCT(RSS,DR2(J7,:))
                     DRIJDI3 = DOT_PRODUCT(RSS,DR3(J7,:))

                     DRIJDJ1 =-DOT_PRODUCT(RSS,DR1(J8,:))
                     DRIJDJ2 =-DOT_PRODUCT(RSS,DR2(J8,:))
                     DRIJDJ3 =-DOT_PRODUCT(RSS,DR3(J8,:))

                     G(J3-2:J3)  = G(J3-2:J3) + DVDR*RSS(:)
                     G(J4-2:J4)  = G(J4-2:J4) - DVDR*RSS(:)

                     G(J5-2)     = G(J5-2) + DVDR*DRIJDI1
                     G(J5-1)     = G(J5-1) + DVDR*DRIJDI2
                     G(J5)       = G(J5)   + DVDR*DRIJDI3

                     G(J6-2)     = G(J6-2) + DVDR*DRIJDJ1
                     G(J6-1)     = G(J6-1) + DVDR*DRIJDJ2
                     G(J6)       = G(J6)   + DVDR*DRIJDJ3
!     DIPOLAR CONTRIBUTION
                     VR = -DPFCT*R4*(GAM - 3.D0*ALP*BET)
                     VA = -DPFCT*BET*R2/ABSRIJ
                     VB = -DPFCT*ALP*R2/ABSRIJ
                     VG =  DPFCT*R2/(3.D0*ABSRIJ)

                     DADR(:) = EI(:)/ABSRIJ - ALP*R2*RSS(:)
                     DBDR(:) = EJ(:)/ABSRIJ - BET*R2*RSS(:)

                     FIJ(:)     = VR*NR(:) + VA*DADR(:) + VB*DBDR(:)

                     G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
                     G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

                     DADPI1 = DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - ALP*R2*DRIJDI1 + DOT_PRODUCT(NR(:),DE1(J7,:))
                     DADPI2 = DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - ALP*R2*DRIJDI2 + DOT_PRODUCT(NR(:),DE2(J7,:))
                     DADPI3 = DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - ALP*R2*DRIJDI3 + DOT_PRODUCT(NR(:),DE3(J7,:))

                     DADPJ1 =-DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - ALP*R2*DRIJDJ1
                     DADPJ2 =-DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - ALP*R2*DRIJDJ2
                     DADPJ3 =-DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - ALP*R2*DRIJDJ3

                     DBDPI1 = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - BET*R2*DRIJDI1
                     DBDPI2 = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - BET*R2*DRIJDI2
                     DBDPI3 = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - BET*R2*DRIJDI3

                     DBDPJ1 =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ - BET*R2*DRIJDJ1 + DOT_PRODUCT(NR(:),DE1(J8,:))
                     DBDPJ2 =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ - BET*R2*DRIJDJ2 + DOT_PRODUCT(NR(:),DE2(J8,:))
                     DBDPJ3 =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ - BET*R2*DRIJDJ3 + DOT_PRODUCT(NR(:),DE3(J8,:))

                     DGDPI1 = DOT_PRODUCT(DE1(J7,:),EJ(:))
                     DGDPI2 = DOT_PRODUCT(DE2(J7,:),EJ(:))
                     DGDPI3 = DOT_PRODUCT(DE3(J7,:),EJ(:))

                     DGDPJ1 = DOT_PRODUCT(EI(:),DE1(J8,:))
                     DGDPJ2 = DOT_PRODUCT(EI(:),DE2(J8,:))
                     DGDPJ3 = DOT_PRODUCT(EI(:),DE3(J8,:))

                     G(J5-2) = G(J5-2) + VR*DRIJDI1/ABSRIJ + VA*DADPI1 + VB*DBDPI1 + VG*DGDPI1
                     G(J5-1) = G(J5-1) + VR*DRIJDI2/ABSRIJ + VA*DADPI2 + VB*DBDPI2 + VG*DGDPI2
                     G(J5)   = G(J5)   + VR*DRIJDI3/ABSRIJ + VA*DADPI3 + VB*DBDPI3 + VG*DGDPI3

                     G(J6-2) = G(J6-2) + VR*DRIJDJ1/ABSRIJ + VA*DADPJ1 + VB*DBDPJ1 + VG*DGDPJ1
                     G(J6-1) = G(J6-1) + VR*DRIJDJ2/ABSRIJ + VA*DADPJ2 + VB*DBDPJ2 + VG*DGDPJ2
                     G(J6)   = G(J6)   + VR*DRIJDJ3/ABSRIJ + VA*DADPJ3 + VB*DBDPJ3 + VG*DGDPJ3

                  ENDIF

                  IF (STEST) THEN
!     LJ CONTRIBUTION
                     D2VDR2 = 672.D0*R12*R2*R2 - 192.D0*R6*R2*R2
!     DIPOLAR CONTRIBUTION
                     DVRDR = 4.D0*DPFCT*R4*(GAM - 3.D0*ALP*BET)/ABSRIJ
                     DVADR = 3.D0*DPFCT*BET*R4
                     DVBDR = 3.D0*DPFCT*ALP*R4
                     DVGDR =-DPFCT*R4
                     DVADB =-DPFCT*R2/ABSRIJ

                     DVRDA = DVADR
                     DVRDB = DVBDR
                     DVRDG = DVGDR
                     DVBDA = DVADB

                     D2RDIXIX = (DOT_PRODUCT(DR1(J7,:),DR1(J7,:)) + DOT_PRODUCT(RSS(:),D2R1(J7,:)) - R2*DRIJDI1**2)/ABSRIJ
                     D2RDIYIY = (DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) + DOT_PRODUCT(RSS(:),D2R2(J7,:)) - R2*DRIJDI2**2)/ABSRIJ
                     D2RDIZIZ = (DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) + DOT_PRODUCT(RSS(:),D2R3(J7,:)) - R2*DRIJDI3**2)/ABSRIJ

                     D2RDIXIY = (DOT_PRODUCT(DR1(J7,:),DR2(J7,:)) + DOT_PRODUCT(RSS(:),D2R12(J7,:)) - R2*DRIJDI1*DRIJDI2)/ABSRIJ
                     D2RDIYIZ = (DOT_PRODUCT(DR2(J7,:),DR3(J7,:)) + DOT_PRODUCT(RSS(:),D2R23(J7,:)) - R2*DRIJDI2*DRIJDI3)/ABSRIJ
                     D2RDIZIX = (DOT_PRODUCT(DR3(J7,:),DR1(J7,:)) + DOT_PRODUCT(RSS(:),D2R31(J7,:)) - R2*DRIJDI3*DRIJDI1)/ABSRIJ

                     D2RDIXJX = (-DOT_PRODUCT(DR1(J7,:),DR1(J8,:)) - R2*DRIJDI1*DRIJDJ1)/ABSRIJ
                     D2RDIYJY = (-DOT_PRODUCT(DR2(J7,:),DR2(J8,:)) - R2*DRIJDI2*DRIJDJ2)/ABSRIJ
                     D2RDIZJZ = (-DOT_PRODUCT(DR3(J7,:),DR3(J8,:)) - R2*DRIJDI3*DRIJDJ3)/ABSRIJ

                     D2RDIXJY = (-DOT_PRODUCT(DR1(J7,:),DR2(J8,:)) - R2*DRIJDI1*DRIJDJ2)/ABSRIJ
                     D2RDIYJZ = (-DOT_PRODUCT(DR2(J7,:),DR3(J8,:)) - R2*DRIJDI2*DRIJDJ3)/ABSRIJ
                     D2RDIZJX = (-DOT_PRODUCT(DR3(J7,:),DR1(J8,:)) - R2*DRIJDI3*DRIJDJ1)/ABSRIJ

                     D2ADRX(:) =-R2*RSS(1)/ABSRIJ*EI(:) - R2*DADR(1)*RSS(:) - R2*ALP*UVX(:) + 2.D0*R4*ALP*RSS(1)*RSS(:)
                     D2ADRY(:) =-R2*RSS(2)/ABSRIJ*EI(:) - R2*DADR(2)*RSS(:) - R2*ALP*UVY(:) + 2.D0*R4*ALP*RSS(2)*RSS(:)
                     D2ADRZ(:) =-R2*RSS(3)/ABSRIJ*EI(:) - R2*DADR(3)*RSS(:) - R2*ALP*UVZ(:) + 2.D0*R4*ALP*RSS(3)*RSS(:)

                     D2BDRX(:) =-R2*RSS(1)/ABSRIJ*EJ(:) - R2*DBDR(1)*RSS(:) - R2*BET*UVX(:) + 2.D0*R4*BET*RSS(1)*RSS(:)
                     D2BDRY(:) =-R2*RSS(2)/ABSRIJ*EJ(:) - R2*DBDR(2)*RSS(:) - R2*BET*UVY(:) + 2.D0*R4*BET*RSS(2)*RSS(:)
                     D2BDRZ(:) =-R2*RSS(3)/ABSRIJ*EJ(:) - R2*DBDR(3)*RSS(:) - R2*BET*UVZ(:) + 2.D0*R4*BET*RSS(3)*RSS(:)

                     D2ADRIX(:) = DE1(J7,:)/ABSRIJ - R2*DRIJDI1/ABSRIJ*EI(:) - DADPI1*R2*RSS(:)   &
                                + 2.D0*ALP*R2*R2*DRIJDI1*RSS(:) - ALP*R2*DR1(J7,:)
                     D2ADRIY(:) = DE2(J7,:)/ABSRIJ - R2*DRIJDI2/ABSRIJ*EI(:) - DADPI2*R2*RSS(:)   &
                                + 2.D0*ALP*R2*R2*DRIJDI2*RSS(:) - ALP*R2*DR2(J7,:)
                     D2ADRIZ(:) = DE3(J7,:)/ABSRIJ - R2*DRIJDI3/ABSRIJ*EI(:) - DADPI3*R2*RSS(:)   &
                                + 2.D0*ALP*R2*R2*DRIJDI3*RSS(:) - ALP*R2*DR3(J7,:)
                     
                     D2ADRJX(:) =-R2*DRIJDJ1/ABSRIJ*EI(:) - DADPJ1*R2*RSS(:) + 2.D0*ALP*R2*R2*DRIJDJ1*RSS(:) + ALP*R2*DR1(J8,:)
                     D2ADRJY(:) =-R2*DRIJDJ2/ABSRIJ*EI(:) - DADPJ2*R2*RSS(:) + 2.D0*ALP*R2*R2*DRIJDJ2*RSS(:) + ALP*R2*DR2(J8,:)
                     D2ADRJZ(:) =-R2*DRIJDJ3/ABSRIJ*EI(:) - DADPJ3*R2*RSS(:) + 2.D0*ALP*R2*R2*DRIJDJ3*RSS(:) + ALP*R2*DR3(J8,:)

                     D2BDRIX(:) =-R2*DRIJDI1/ABSRIJ*EJ(:) - DBDPI1*R2*RSS(:) + 2.D0*BET*R2*R2*DRIJDI1*RSS(:)      &
                                - BET*R2*DR1(J7,:)
                     D2BDRIY(:) =-R2*DRIJDI2/ABSRIJ*EJ(:) - DBDPI2*R2*RSS(:) + 2.D0*BET*R2*R2*DRIJDI2*RSS(:)      &
                                - BET*R2*DR2(J7,:)
                     D2BDRIZ(:) =-R2*DRIJDI3/ABSRIJ*EJ(:) - DBDPI3*R2*RSS(:) + 2.D0*BET*R2*R2*DRIJDI3*RSS(:)      &
                                - BET*R2*DR3(J7,:)

                     D2BDRJX(:) = DE1(J8,:)/ABSRIJ - R2*DRIJDJ1/ABSRIJ*EJ(:) - DBDPJ1*R2*RSS(:)   &
                                + 2.D0*BET*R2*R2*DRIJDJ1*RSS(:) + BET*R2*DR1(J8,:)
                     D2BDRJY(:) = DE2(J8,:)/ABSRIJ - R2*DRIJDJ2/ABSRIJ*EJ(:) - DBDPJ2*R2*RSS(:)   &
                                + 2.D0*BET*R2*R2*DRIJDJ2*RSS(:) + BET*R2*DR2(J8,:)
                     D2BDRJZ(:) = DE3(J8,:)/ABSRIJ - R2*DRIJDJ3/ABSRIJ*EJ(:) - DBDPJ3*R2*RSS(:)   &
                                + 2.D0*BET*R2*R2*DRIJDJ3*RSS(:) + BET*R2*DR3(J8,:)

                     D2ADIXIX = (DOT_PRODUCT(D2R1(J7,:),EI(:)) - R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDI1      &
                              + DOT_PRODUCT(DR1(J7,:),DE1(J7,:)) - DADPI1*DRIJDI1/ABSRIJ  &
                              + ALP*R2*DRIJDI1**2/ABSRIJ - ALP*D2RDIXIX + DOT_PRODUCT(DR1(J7,:),DE1(J7,:))    &
                              - R2*DRIJDI1*DOT_PRODUCT(RSS(:),DE1(J7,:)) + DOT_PRODUCT(RSS(:),D2E1(J7,:)))/ABSRIJ
                     D2ADIYIY = (DOT_PRODUCT(D2R2(J7,:),EI(:)) - R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDI2      &
                              + DOT_PRODUCT(DR2(J7,:),DE2(J7,:)) - DADPI2*DRIJDI2/ABSRIJ  &
                              + ALP*R2*DRIJDI2**2/ABSRIJ - ALP*D2RDIYIY + DOT_PRODUCT(DR2(J7,:),DE2(J7,:))    &
                              - R2*DRIJDI2*DOT_PRODUCT(RSS(:),DE2(J7,:)) + DOT_PRODUCT(RSS(:),D2E2(J7,:)))/ABSRIJ
                     D2ADIZIZ = (DOT_PRODUCT(D2R3(J7,:),EI(:)) - R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDI3      &
                              + DOT_PRODUCT(DR3(J7,:),DE3(J7,:)) - DADPI3*DRIJDI3/ABSRIJ  &
                              + ALP*R2*DRIJDI3**2/ABSRIJ - ALP*D2RDIZIZ + DOT_PRODUCT(DR3(J7,:),DE3(J7,:))    &
                              - R2*DRIJDI3*DOT_PRODUCT(RSS(:),DE3(J7,:)) + DOT_PRODUCT(RSS(:),D2E3(J7,:)))/ABSRIJ
                   
                     D2BDIXIX = (DOT_PRODUCT(D2R1(J7,:),EJ(:)) - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDI1      &
                              - DBDPI1*DRIJDI1/ABSRIJ + BET*R2*DRIJDI1**2/ABSRIJ - BET*D2RDIXIX)/ABSRIJ
                     D2BDIYIY = (DOT_PRODUCT(D2R2(J7,:),EJ(:)) - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDI2      &
                              - DBDPI2*DRIJDI2/ABSRIJ + BET*R2*DRIJDI2**2/ABSRIJ - BET*D2RDIYIY)/ABSRIJ
                     D2BDIZIZ = (DOT_PRODUCT(D2R3(J7,:),EJ(:)) - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDI3      &
                              - DBDPI3*DRIJDI3/ABSRIJ + BET*R2*DRIJDI3**2/ABSRIJ - BET*D2RDIZIZ)/ABSRIJ

                     D2ADIXJX =-R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDJ1/ABSRIJ - DADPJ1*R2*DRIJDI1 + ALP*R2*R2*DRIJDJ1*DRIJDI1    &
                               - ALP*D2RDIXJX/ABSRIJ - (DOT_PRODUCT(DR1(J8,:),DE1(J7,:))  &
                               + R2*DRIJDJ1*DOT_PRODUCT(RSS(:),DE1(J7,:)))/ABSRIJ 
                     D2ADIYJY =-R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDJ2/ABSRIJ - DADPJ2*R2*DRIJDI2 + ALP*R2*R2*DRIJDJ2*DRIJDI2    &
                               - ALP*D2RDIYJY/ABSRIJ - (DOT_PRODUCT(DR2(J8,:),DE2(J7,:))  &
                               + R2*DRIJDJ2*DOT_PRODUCT(RSS(:),DE2(J7,:)))/ABSRIJ 
                     D2ADIZJZ =-R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDJ3/ABSRIJ - DADPJ3*R2*DRIJDI3 + ALP*R2*R2*DRIJDJ3*DRIJDI3    &
                               - ALP*D2RDIZJZ/ABSRIJ - (DOT_PRODUCT(DR3(J8,:),DE3(J7,:))  &
                               + R2*DRIJDJ3*DOT_PRODUCT(RSS(:),DE3(J7,:)))/ABSRIJ 

                     D2BDIXJX = (DOT_PRODUCT(DR1(J7,:),DE1(J8,:)) - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDJ1   & 
                              - DBDPJ1*DRIJDI1/ABSRIJ + BET*R2*DRIJDJ1*DRIJDI1/ABSRIJ - BET*D2RDIXJX)/ABSRIJ
                     D2BDIYJY = (DOT_PRODUCT(DR2(J7,:),DE2(J8,:)) - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDJ2   & 
                              - DBDPJ2*DRIJDI2/ABSRIJ + BET*R2*DRIJDJ2*DRIJDI2/ABSRIJ - BET*D2RDIYJY)/ABSRIJ
                     D2BDIZJZ = (DOT_PRODUCT(DR3(J7,:),DE3(J8,:)) - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDJ3   & 
                              - DBDPJ3*DRIJDI3/ABSRIJ + BET*R2*DRIJDJ3*DRIJDI3/ABSRIJ - BET*D2RDIZJZ)/ABSRIJ

                     D2ADIXIY = (DOT_PRODUCT(D2R12(J7,:),EI(:)) - R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDI2     &
                              + DOT_PRODUCT(DR1(J7,:),DE2(J7,:)) - DADPI2*DRIJDI1/ABSRIJ  &
                              + ALP*R2*DRIJDI2*DRIJDI1/ABSRIJ - ALP*D2RDIXIY + DOT_PRODUCT(DR2(J7,:),DE1(J7,:))    &
                              - R2*DRIJDI2*DOT_PRODUCT(RSS(:),DE1(J7,:)) + DOT_PRODUCT(RSS(:),D2E12(J7,:)))/ABSRIJ
                     D2ADIYIZ = (DOT_PRODUCT(D2R23(J7,:),EI(:)) - R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDI3     &
                              + DOT_PRODUCT(DR2(J7,:),DE3(J7,:)) - DADPI3*DRIJDI2/ABSRIJ  &
                              + ALP*R2*DRIJDI3*DRIJDI2/ABSRIJ - ALP*D2RDIYIZ + DOT_PRODUCT(DR3(J7,:),DE2(J7,:))    &
                              - R2*DRIJDI3*DOT_PRODUCT(RSS(:),DE2(J7,:)) + DOT_PRODUCT(RSS(:),D2E23(J7,:)))/ABSRIJ
                     D2ADIZIX = (DOT_PRODUCT(D2R31(J7,:),EI(:)) - R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDI1     &
                              + DOT_PRODUCT(DR3(J7,:),DE1(J7,:)) - DADPI1*DRIJDI3/ABSRIJ  &
                              + ALP*R2*DRIJDI1*DRIJDI3/ABSRIJ - ALP*D2RDIZIX + DOT_PRODUCT(DR1(J7,:),DE3(J7,:))    &
                              - R2*DRIJDI1*DOT_PRODUCT(RSS(:),DE3(J7,:)) + DOT_PRODUCT(RSS(:),D2E31(J7,:)))/ABSRIJ

                     D2BDIXIY = (DOT_PRODUCT(D2R12(J7,:),EJ(:)) - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDI2     &
                              - DBDPI2*DRIJDI1/ABSRIJ + BET*R2*DRIJDI2*DRIJDI1/ABSRIJ - BET*D2RDIXIY)/ABSRIJ
                     D2BDIYIZ = (DOT_PRODUCT(D2R23(J7,:),EJ(:)) - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDI3     &
                              - DBDPI3*DRIJDI2/ABSRIJ + BET*R2*DRIJDI3*DRIJDI2/ABSRIJ - BET*D2RDIYIZ)/ABSRIJ
                     D2BDIZIX = (DOT_PRODUCT(D2R31(J7,:),EJ(:)) - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDI1     &
                              - DBDPI1*DRIJDI3/ABSRIJ + BET*R2*DRIJDI1*DRIJDI3/ABSRIJ - BET*D2RDIZIX)/ABSRIJ

                     D2ADIXJY = (-R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDJ2 - DADPJ2*DRIJDI1/ABSRIJ   &
                              + ALP*R2*DRIJDJ2*DRIJDI1/ABSRIJ - ALP*D2RDIXJY - DOT_PRODUCT(DR2(J8,:),DE1(J7,:))         &
                              - R2*DRIJDJ2*DOT_PRODUCT(RSS(:),DE1(J7,:)))/ABSRIJ
                     D2ADIYJZ = (-R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDJ3 - DADPJ3*DRIJDI2/ABSRIJ   &
                              + ALP*R2*DRIJDJ3*DRIJDI2/ABSRIJ - ALP*D2RDIYJZ - DOT_PRODUCT(DR3(J8,:),DE2(J7,:))         &
                              - R2*DRIJDJ3*DOT_PRODUCT(RSS(:),DE2(J7,:)))/ABSRIJ
                     D2ADIZJX = (-R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDJ1 - DADPJ1*DRIJDI3/ABSRIJ   &
                              + ALP*R2*DRIJDJ1*DRIJDI3/ABSRIJ - ALP*D2RDIZJX - DOT_PRODUCT(DR1(J8,:),DE3(J7,:))         &
                              - R2*DRIJDJ1*DOT_PRODUCT(RSS(:),DE3(J7,:)))/ABSRIJ

                     D2BDIXJY = (DOT_PRODUCT(DR1(J7,:),DE2(J8,:)) - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDJ2   &
                              - DBDPJ2*DRIJDI1/ABSRIJ + BET*R2*DRIJDJ2*DRIJDI1/ABSRIJ - BET*D2RDIXJY)/ABSRIJ 
                     D2BDIYJZ = (DOT_PRODUCT(DR2(J7,:),DE3(J8,:)) - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDJ3   &
                              - DBDPJ3*DRIJDI2/ABSRIJ + BET*R2*DRIJDJ3*DRIJDI2/ABSRIJ - BET*D2RDIYJZ)/ABSRIJ 
                     D2BDIZJX = (DOT_PRODUCT(DR3(J7,:),DE1(J8,:)) - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDJ1   &
                              - DBDPJ1*DRIJDI3/ABSRIJ + BET*R2*DRIJDJ1*DRIJDI3/ABSRIJ - BET*D2RDIZJX)/ABSRIJ 

!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES

!     xi,xi
                     HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2*RSS(1)*RSS(1) + DVDR      &         ! LJ
                                     + DVRDR*R2*RSS(1)*RSS(1) + (DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(1)         &
                                     + VR*(1.D0 - R2*RSS(1)*RSS(1))/ABSRIJ  + (DVADR*NR(1) + DVADB*DBDR(1))*DADR(1)     &
                                     + VA*D2ADRX(1) + (DVBDR*NR(1) + DVBDA*DADR(1) )*DBDR(1) + VB*D2BDRX(1)
!     yi,yi             
                     HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2*RSS(2)*RSS(2) + DVDR      &         ! LJ
                                      + DVRDR*R2*RSS(2)*RSS(2) + (DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(2)        &
                                     + VR*(1.D0 - R2*RSS(2)*RSS(2))/ABSRIJ + (DVADR*NR(2) + DVADB*DBDR(2))*DADR(2)      &
                                     + VA*D2ADRY(2) + (DVBDR*NR(2) + DVBDA*DADR(2) )*DBDR(2) + VB*D2BDRY(2)
!     zi,zi
                     HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2*RSS(3)*RSS(3) + DVDR      &         ! LJ
                                     + DVRDR*R2*RSS(3)*RSS(3) + (DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(3)         &
                                     + VR*(1.D0 - R2*RSS(3)*RSS(3))/ABSRIJ + (DVADR*NR(3) + DVADB*DBDR(3))*DADR(3)      &
                                     + VA*D2ADRZ(3) + (DVBDR*NR(3) + DVBDA*DADR(3) )*DBDR(3) + VB*D2BDRZ(3)
!     pi1,pi1
                     HESS(J5-2,J5-2) = HESS(J5-2,J5-2) + D2VDR2*DRIJDI1*DRIJDI1 + DVDR*DOT_PRODUCT(DR1(J7,:),DR1(J7,:)) &
                                     + DVDR*DOT_PRODUCT(RSS,D2R1(J7,:)) + (DVRDR*DRIJDI1/ABSRIJ + DVRDA*DADPI1 + DVRDB*DBDPI1     &
                                     + DVRDG*DGDPI1)*DRIJDI1/ABSRIJ + VR*(DOT_PRODUCT(DR1(J7,:),DR1(J7,:))    &
                                     + DOT_PRODUCT(RSS(:),D2R1(J7,:)))/ABSRIJ - VR*DRIJDI1**2*R2/ABSRIJ + (DVADR*DRIJDI1/ABSRIJ   &
                                     + DVADB*DBDPI1)*DADPI1 + VA*D2ADIXIX + (DVBDR*DRIJDI1/ABSRIJ + DVBDA*DADPI1)*DBDPI1          &
                                     + VB*D2BDIXIX + DVGDR*DRIJDI1*DGDPI1/ABSRIJ + VG*DOT_PRODUCT(D2E1(J7,:),EJ(:))
!     pi2,pi2
                     HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + D2VDR2*DRIJDI2*DRIJDI2 + DVDR*DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) &
                                     + DVDR*DOT_PRODUCT(RSS,D2R2(J7,:)) + (DVRDR*DRIJDI2/ABSRIJ + DVRDA*DADPI2 + DVRDB*DBDPI2     &
                                     + DVRDG*DGDPI2)*DRIJDI2/ABSRIJ + VR*(DOT_PRODUCT(DR2(J7,:),DR2(J7,:))    &
                                     + DOT_PRODUCT(RSS(:),D2R2(J7,:)))/ABSRIJ - VR*DRIJDI2**2*R2/ABSRIJ + (DVADR*DRIJDI2/ABSRIJ   &
                                     + DVADB*DBDPI2)*DADPI2 + VA*D2ADIYIY + (DVBDR*DRIJDI2/ABSRIJ + DVBDA*DADPI2)*DBDPI2          &
                                     + VB*D2BDIYIY + DVGDR*DRIJDI2*DGDPI2/ABSRIJ + VG*DOT_PRODUCT(D2E2(J7,:),EJ(:))
!     pi3,pi3
                     HESS(J5,J5)     = HESS(J5,J5)     + D2VDR2*DRIJDI3*DRIJDI3 + DVDR*DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) &
                                     + DVDR*DOT_PRODUCT(RSS,D2R3(J7,:)) + (DVRDR*DRIJDI3/ABSRIJ + DVRDA*DADPI3 + DVRDB*DBDPI3     &
                                     + DVRDG*DGDPI3)*DRIJDI3/ABSRIJ + VR*(DOT_PRODUCT(DR3(J7,:),DR3(J7,:))    &
                                     + DOT_PRODUCT(RSS(:),D2R3(J7,:)))/ABSRIJ - VR*DRIJDI3**2*R2/ABSRIJ + (DVADR*DRIJDI3/ABSRIJ   &
                                     + DVADB*DBDPI3)*DADPI3 + VA*D2ADIZIZ + (DVBDR*DRIJDI3/ABSRIJ + DVBDA*DADPI3)*DBDPI3          &
                                     + VB*D2BDIZIZ + DVGDR*DRIJDI3*DGDPI3/ABSRIJ + VG*DOT_PRODUCT(D2E3(J7,:),EJ(:))

!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME MOLECULE, DIFFERENT COORDINATES
!     xi,yi
                     DUMMY = D2VDR2*RSS(1)*RSS(2) + DVRDR*R2*RSS(2)*RSS(1) + (DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(1)      &
                           - VR*R2*RSS(2)*RSS(1)/ABSRIJ + (DVADR*NR(2) + DVADB*DBDR(2))*DADR(1) + VA*D2ADRY(1)          &
                           + (DVBDR*NR(2) + DVBDA*DADR(2))*DBDR(1) + VB*D2BDRY(1) 
                     HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
                     HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
                     DUMMY = D2VDR2*RSS(2)*RSS(3) + DVRDR*R2*RSS(3)*RSS(2) + (DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(2)      &
                           - VR*R2*RSS(3)*RSS(2)/ABSRIJ + (DVADR*NR(3) + DVADB*DBDR(3))*DADR(2) + VA*D2ADRZ(2)          &
                           + (DVBDR*NR(3) + DVBDA*DADR(3))*DBDR(2) + VB*D2BDRZ(2) 
                     HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
                     HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
                     DUMMY = D2VDR2*RSS(3)*RSS(1) + DVRDR*R2*RSS(1)*RSS(3) + (DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(3)      &
                           - VR*R2*RSS(1)*RSS(3)/ABSRIJ + (DVADR*NR(1) + DVADB*DBDR(1))*DADR(3) + VA*D2ADRX(3)          &
                           + (DVBDR*NR(1) + DVBDA*DADR(1))*DBDR(3) + VB*D2BDRX(3)
                     HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
                     HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
!     xi,pi1
                     DUMMY = D2VDR2*DRIJDI1*RSS(1) + DVDR*DR1(J7,1) + DVRDR*R2*DRIJDI1*RSS(1) + (DVRDA*DADPI1 + DVRDB*DBDPI1      &
                           + DVRDG*DGDPI1)*NR(1) + VR*DR1(J7,1)/ABSRIJ - VR*R2*DRIJDI1*RSS(1)/ABSRIJ + (DVADR*DRIJDI1/ABSRIJ      &
                           + DVADB*DBDPI1)*DADR(1) + VA*D2ADRIX(1) + (DVBDR*DRIJDI1/ABSRIJ + DVBDA*DADPI1)*DBDR(1) + VB*D2BDRIX(1)
                     HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
                     HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1
                     DUMMY = D2VDR2*DRIJDI1*RSS(2) + DVDR*DR1(J7,2) + DVRDR*R2*DRIJDI1*RSS(2) + (DVRDA*DADPI1 + DVRDB*DBDPI1      &
                           + DVRDG*DGDPI1)*NR(2) + VR*DR1(J7,2)/ABSRIJ - VR*R2*DRIJDI1*RSS(2)/ABSRIJ + (DVADR*DRIJDI1/ABSRIJ      &
                           + DVADB*DBDPI1)*DADR(2) + VA*D2ADRIX(2) + (DVBDR*DRIJDI1/ABSRIJ + DVBDA*DADPI1)*DBDR(2) + VB*D2BDRIX(2)
                     HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
                     HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
                     DUMMY = D2VDR2*DRIJDI1*RSS(3) + DVDR*DR1(J7,3) + DVRDR*R2*DRIJDI1*RSS(3) + (DVRDA*DADPI1 + DVRDB*DBDPI1      &
                           + DVRDG*DGDPI1)*NR(3) + VR*DR1(J7,3)/ABSRIJ - VR*R2*DRIJDI1*RSS(3)/ABSRIJ + (DVADR*DRIJDI1/ABSRIJ      &
                           + DVADB*DBDPI1)*DADR(3) + VA*D2ADRIX(3) + (DVBDR*DRIJDI1/ABSRIJ + DVBDA*DADPI1)*DBDR(3) + VB*D2BDRIX(3)
                     HESS(J3,J5-2)   = HESS(J3,J5-2) + DUMMY
                     HESS(J5-2,J3)   = HESS(J5-2,J3) + DUMMY
!     xi,pi2
                     DUMMY = D2VDR2*DRIJDI2*RSS(1) + DVDR*DR2(J7,1) + DVRDR*R2*DRIJDI2*RSS(1) + (DVRDA*DADPI2 + DVRDB*DBDPI2      &
                           + DVRDG*DGDPI2)*NR(1) + VR*DR2(J7,1)/ABSRIJ - VR*R2*DRIJDI2*RSS(1)/ABSRIJ + (DVADR*DRIJDI2/ABSRIJ      &
                           + DVADB*DBDPI2)*DADR(1) + VA*D2ADRIY(1) + (DVBDR*DRIJDI2/ABSRIJ + DVBDA*DADPI2)*DBDR(1) + VB*D2BDRIY(1)
                     HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
                     HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2
                     DUMMY = D2VDR2*DRIJDI2*RSS(2) + DVDR*DR2(J7,2) + DVRDR*R2*DRIJDI2*RSS(2) + (DVRDA*DADPI2 + DVRDB*DBDPI2      &
                           + DVRDG*DGDPI2)*NR(2) + VR*DR2(J7,2)/ABSRIJ - VR*R2*DRIJDI2*RSS(2)/ABSRIJ + (DVADR*DRIJDI2/ABSRIJ      &
                           + DVADB*DBDPI2)*DADR(2) + VA*D2ADRIY(2) + (DVBDR*DRIJDI2/ABSRIJ + DVBDA*DADPI2)*DBDR(2) + VB*D2BDRIY(2)
                     HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
                     HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
                     DUMMY = D2VDR2*DRIJDI2*RSS(3) + DVDR*DR2(J7,3) + DVRDR*R2*DRIJDI2*RSS(3) + (DVRDA*DADPI2 + DVRDB*DBDPI2      &
                           + DVRDG*DGDPI2)*NR(3) + VR*DR2(J7,3)/ABSRIJ - VR*R2*DRIJDI2*RSS(3)/ABSRIJ + (DVADR*DRIJDI2/ABSRIJ      &
                           + DVADB*DBDPI2)*DADR(3) + VA*D2ADRIY(3) + (DVBDR*DRIJDI2/ABSRIJ + DVBDA*DADPI2)*DBDR(3) + VB*D2BDRIY(3)
                     HESS(J3,J5-1)  = HESS(J3,J5-1) + DUMMY
                     HESS(J5-1,J3)  = HESS(J5-1,J3) + DUMMY
!     xi,pi3
                     DUMMY = D2VDR2*DRIJDI3*RSS(1) + DVDR*DR3(J7,1) + DVRDR*R2*DRIJDI3*RSS(1) + (DVRDA*DADPI3 + DVRDB*DBDPI3      &
                           + DVRDG*DGDPI3)*NR(1) + VR*DR3(J7,1)/ABSRIJ - VR*R2*DRIJDI3*RSS(1)/ABSRIJ + (DVADR*DRIJDI3/ABSRIJ      &
                           + DVADB*DBDPI3)*DADR(1) + VA*D2ADRIZ(1) + (DVBDR*DRIJDI3/ABSRIJ + DVBDA*DADPI3)*DBDR(1) + VB*D2BDRIZ(1)
                     HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
                     HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi3
                     DUMMY = D2VDR2*DRIJDI3*RSS(2) + DVDR*DR3(J7,2) + DVRDR*R2*DRIJDI3*RSS(2) + (DVRDA*DADPI3 + DVRDB*DBDPI3      &
                           + DVRDG*DGDPI3)*NR(2) + VR*DR3(J7,2)/ABSRIJ - VR*R2*DRIJDI3*RSS(2)/ABSRIJ + (DVADR*DRIJDI3/ABSRIJ      &
                           + DVADB*DBDPI3)*DADR(2) + VA*D2ADRIZ(2) + (DVBDR*DRIJDI3/ABSRIJ + DVBDA*DADPI3)*DBDR(2) + VB*D2BDRIZ(2)
                     HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
                     HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi3
                     DUMMY = D2VDR2*DRIJDI3*RSS(3) + DVDR*DR3(J7,3) + DVRDR*R2*DRIJDI3*RSS(3) + (DVRDA*DADPI3 + DVRDB*DBDPI3      &
                           + DVRDG*DGDPI3)*NR(3) + VR*DR3(J7,3)/ABSRIJ - VR*R2*DRIJDI3*RSS(3)/ABSRIJ + (DVADR*DRIJDI3/ABSRIJ      &
                           + DVADB*DBDPI3)*DADR(3) + VA*D2ADRIZ(3) + (DVBDR*DRIJDI3/ABSRIJ + DVBDA*DADPI3)*DBDR(3) + VB*D2BDRIZ(3)
                     HESS(J3,J5)     = HESS(J3,J5) + DUMMY
                     HESS(J5,J3)     = HESS(J5,J3) + DUMMY
!     pi1,pi2
                     DUMMY = D2VDR2*DRIJDI1*DRIJDI2 + DVDR*DOT_PRODUCT(DR2(J7,:),DR1(J7,:)) + DVDR*DOT_PRODUCT(RSS,D2R12(J7,:))   &
                           + (DVRDR*DRIJDI2/ABSRIJ + DVRDA*DADPI2 + DVRDB*DBDPI2 + DVRDG*DGDPI2)*DRIJDI1/ABSRIJ         &
                           + VR*(DOT_PRODUCT(DR1(J7,:),DR2(J7,:)) + DOT_PRODUCT(RSS(:),D2R12(J7,:)))/ABSRIJ   &
                           - VR*DRIJDI2*DRIJDI1*R2/ABSRIJ + (DVADR*DRIJDI2/ABSRIJ + DVADB*DBDPI2)*DADPI1 + VA*D2ADIXIY  &
                           + (DVBDR*DRIJDI2/ABSRIJ + DVBDA*DADPI2)*DBDPI1 + VB*D2BDIXIY + DVGDR*DRIJDI2*DGDPI1/ABSRIJ   &
                           + VG*DOT_PRODUCT(D2E12(J7,:),EJ(:))
                     HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
                     HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
                     DUMMY = D2VDR2*DRIJDI2*DRIJDI3 + DVDR*DOT_PRODUCT(DR3(J7,:),DR2(J7,:)) + DVDR*DOT_PRODUCT(RSS,D2R23(J7,:))   &
                           + (DVRDR*DRIJDI3/ABSRIJ + DVRDA*DADPI3 + DVRDB*DBDPI3 + DVRDG*DGDPI3)*DRIJDI2/ABSRIJ         &
                           + VR*(DOT_PRODUCT(DR2(J7,:),DR3(J7,:)) + DOT_PRODUCT(RSS(:),D2R23(J7,:)))/ABSRIJ   &
                           - VR*DRIJDI3*DRIJDI2*R2/ABSRIJ + (DVADR*DRIJDI3/ABSRIJ + DVADB*DBDPI3)*DADPI2 + VA*D2ADIYIZ  &
                           + (DVBDR*DRIJDI3/ABSRIJ + DVBDA*DADPI3)*DBDPI2 + VB*D2BDIYIZ + DVGDR*DRIJDI3*DGDPI2/ABSRIJ   &
                           + VG*DOT_PRODUCT(D2E23(J7,:),EJ(:))
                     HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
                     HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
                     DUMMY = D2VDR2*DRIJDI3*DRIJDI1 + DVDR*DOT_PRODUCT(DR1(J7,:),DR3(J7,:)) + DVDR*DOT_PRODUCT(RSS,D2R31(J7,:))   &
                           + (DVRDR*DRIJDI1/ABSRIJ + DVRDA*DADPI1 + DVRDB*DBDPI1 + DVRDG*DGDPI1)*DRIJDI3/ABSRIJ         &
                           + VR*(DOT_PRODUCT(DR3(J7,:),DR1(J7,:)) + DOT_PRODUCT(RSS(:),D2R31(J7,:)))/ABSRIJ   &
                           - VR*DRIJDI1*DRIJDI3*R2/ABSRIJ + (DVADR*DRIJDI1/ABSRIJ + DVADB*DBDPI1)*DADPI3 + VA*D2ADIZIX  &
                           + (DVBDR*DRIJDI1/ABSRIJ + DVBDA*DADPI1)*DBDPI3 + VB*D2BDIZIX + DVGDR*DRIJDI1*DGDPI3/ABSRIJ   &
                           + VG*DOT_PRODUCT(D2E31(J7,:),EJ(:))
                     HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
                     HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE

!     xi,xj
                     HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2*RSS(1)*RSS(1) - DVDR      &         ! LJ
                                     - DVRDR*R2*RSS(1)*RSS(1) - (DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(1)         &
                                     - VR*(1.D0 - R2*RSS(1)*RSS(1))/ABSRIJ  - (DVADR*NR(1) + DVADB*DBDR(1))*DADR(1)     &
                                     - VA*D2ADRX(1) - (DVBDR*NR(1) + DVBDA*DADR(1) )*DBDR(1) - VB*D2BDRX(1)
!     yi,yj             
                     HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2*RSS(2)*RSS(2) - DVDR      &         ! LJ
                                     - DVRDR*R2*RSS(2)*RSS(2) - (DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(2)        &
                                     - VR*(1.D0 - R2*RSS(2)*RSS(2))/ABSRIJ -(DVADR*NR(2) + DVADB*DBDR(2))*DADR(2)      &
                                     - VA*D2ADRY(2) - (DVBDR*NR(2) + DVBDA*DADR(2) )*DBDR(2) - VB*D2BDRY(2)
!     zi,zj
                     HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2*RSS(3)*RSS(3) - DVDR      &         ! LJ
                                     - DVRDR*R2*RSS(3)*RSS(3) - (DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(3)         &
                                     - VR*(1.D0 - R2*RSS(3)*RSS(3))/ABSRIJ - (DVADR*NR(3) + DVADB*DBDR(3))*DADR(3)      &
                                     - VA*D2ADRZ(3) - (DVBDR*NR(3) + DVBDA*DADR(3) )*DBDR(3) - VB*D2BDRZ(3)
!     pi1,pi1
                     HESS(J5-2,J6-2) = HESS(J5-2,J6-2) + D2VDR2*DRIJDJ1*DRIJDI1 - DVDR*DOT_PRODUCT(DR1(J7,:),DR1(J8,:)) &
                                     + (DVRDR*DRIJDJ1/ABSRIJ + DVRDA*DADPJ1 + DVRDB*DBDPJ1 + DVRDG*DGDPJ1)*DRIJDI1/ABSRIJ         &
                                     - VR*DOT_PRODUCT(DR1(J7,:),DR1(J8,:))/ABSRIJ - VR*DRIJDJ1*DRIJDI1*R2/ABSRIJ        &
                                     + (DVADR*DRIJDJ1/ABSRIJ + DVADB*DBDPJ1)*DADPI1  + VA*D2ADIXJX + (DVBDR*DRIJDJ1/ABSRIJ        &
                                     + DVBDA*DADPJ1)*DBDPI1 + VB*D2BDIXJX + DVGDR*DRIJDJ1*DGDPI1/ABSRIJ       &
                                     + VG*DOT_PRODUCT(DE1(J7,:),DE1(J8,:))
!     pi2,pi2
                     HESS(J5-1,J6-1) = HESS(J5-1,J6-1) + D2VDR2*DRIJDJ2*DRIJDI2 - DVDR*DOT_PRODUCT(DR2(J7,:),DR2(J8,:)) &
                                     + (DVRDR*DRIJDJ2/ABSRIJ + DVRDA*DADPJ2 + DVRDB*DBDPJ2 + DVRDG*DGDPJ2)*DRIJDI2/ABSRIJ         &
                                     - VR*DOT_PRODUCT(DR2(J7,:),DR2(J8,:))/ABSRIJ - VR*DRIJDJ2*DRIJDI2*R2/ABSRIJ        &
                                     + (DVADR*DRIJDJ2/ABSRIJ + DVADB*DBDPJ2)*DADPI2 + VA*D2ADIYJY + (DVBDR*DRIJDJ2/ABSRIJ         &
                                     + DVBDA*DADPJ2)*DBDPI2 + VB*D2BDIYJY + DVGDR*DRIJDJ2*DGDPI2/ABSRIJ       &
                                     + VG*DOT_PRODUCT(DE2(J7,:),DE2(J8,:))
!     pi3,pi3
                     HESS(J5,J6)     = HESS(J5,J6)     + D2VDR2*DRIJDJ3*DRIJDI3 - DVDR*DOT_PRODUCT(DR3(J7,:),DR3(J8,:)) &
                                     + (DVRDR*DRIJDJ3/ABSRIJ + DVRDA*DADPJ3 + DVRDB*DBDPJ3 + DVRDG*DGDPJ3)*DRIJDI3/ABSRIJ         &
                                     - VR*DOT_PRODUCT(DR3(J7,:),DR3(J8,:))/ABSRIJ - VR*DRIJDJ3*DRIJDI3*R2/ABSRIJ        &
                                     + (DVADR*DRIJDJ3/ABSRIJ + DVADB*DBDPJ3)*DADPI3 + VA*D2ADIZJZ + (DVBDR*DRIJDJ3/ABSRIJ         &
                                     + DVBDA*DADPJ3)*DBDPI3 + VB*D2BDIZJZ + DVGDR*DRIJDJ3*DGDPI3/ABSRIJ       &
                                     + VG*DOT_PRODUCT(DE3(J7,:),DE3(J8,:))

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES

!     xi,yj
                     DUMMY =-D2VDR2*RSS(1)*RSS(2) - DVRDR*R2*RSS(2)*RSS(1) - (DVRDA*DADR(2) + DVRDB*DBDR(2))*NR(1)      &
                           + VR*R2*RSS(2)*RSS(1)/ABSRIJ - (DVADR*NR(2) + DVADB*DBDR(2))*DADR(1) - VA*D2ADRY(1)          &
                           - (DVBDR*NR(2) + DVBDA*DADR(2))*DBDR(1) - VB*D2BDRY(1) 
                     HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
                     HESS(J3-1,J4-2) = HESS(J3-1,J4-2) + DUMMY
!     yi,zj
                     DUMMY =-D2VDR2*RSS(2)*RSS(3) - DVRDR*R2*RSS(3)*RSS(2) - (DVRDA*DADR(3) + DVRDB*DBDR(3))*NR(2)      &
                           + VR*R2*RSS(3)*RSS(2)/ABSRIJ - (DVADR*NR(3) + DVADB*DBDR(3))*DADR(2) - VA*D2ADRZ(2)          &
                           - (DVBDR*NR(3) + DVBDA*DADR(3))*DBDR(2) - VB*D2BDRZ(2) 
                     HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
                     HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
                     DUMMY =-D2VDR2*RSS(3)*RSS(1) - DVRDR*R2*RSS(1)*RSS(3) - (DVRDA*DADR(1) + DVRDB*DBDR(1))*NR(3)      &
                           + VR*R2*RSS(1)*RSS(3)/ABSRIJ - (DVADR*NR(1) + DVADB*DBDR(1))*DADR(3) - VA*D2ADRX(3)          &
                           - (DVBDR*NR(1) + DVBDA*DADR(1))*DBDR(3) - VB*D2BDRX(3)
                     HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
                     HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY
!     xi,pj1
                     DUMMY = D2VDR2*DRIJDJ1*RSS(1) - DVDR*DR1(J8,1) + DVRDR*R2*DRIJDJ1*RSS(1) + (DVRDA*DADPJ1 + DVRDB*DBDPJ1      &
                           + DVRDG*DGDPJ1)*NR(1) - VR*DR1(J8,1)/ABSRIJ - VR*R2*DRIJDJ1*RSS(1)/ABSRIJ + (DVADR*DRIJDJ1/ABSRIJ      &
                           + DVADB*DBDPJ1)*DADR(1) + VA*D2ADRJX(1) + (DVBDR*DRIJDJ1/ABSRIJ + DVBDA*DADPJ1)*DBDR(1) + VB*D2BDRJX(1)
                     HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
                     HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
                     DUMMY = D2VDR2*DRIJDJ1*RSS(2) - DVDR*DR1(J8,2) + DVRDR*R2*DRIJDJ1*RSS(2) + (DVRDA*DADPJ1 + DVRDB*DBDPJ1      &
                           + DVRDG*DGDPJ1)*NR(2) - VR*DR1(J8,2)/ABSRIJ - VR*R2*DRIJDJ1*RSS(2)/ABSRIJ + (DVADR*DRIJDJ1/ABSRIJ      &
                           + DVADB*DBDPJ1)*DADR(2) + VA*D2ADRJX(2) + (DVBDR*DRIJDJ1/ABSRIJ + DVBDA*DADPJ1)*DBDR(2) + VB*D2BDRJX(2)
                     HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
                     HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
                     DUMMY = D2VDR2*DRIJDJ1*RSS(3) - DVDR*DR1(J8,3) + DVRDR*R2*DRIJDJ1*RSS(3) + (DVRDA*DADPJ1 + DVRDB*DBDPJ1      &
                           + DVRDG*DGDPJ1)*NR(3) - VR*DR1(J8,3)/ABSRIJ - VR*R2*DRIJDJ1*RSS(3)/ABSRIJ + (DVADR*DRIJDJ1/ABSRIJ      &
                           + DVADB*DBDPJ1)*DADR(3) + VA*D2ADRJX(3) + (DVBDR*DRIJDJ1/ABSRIJ + DVBDA*DADPJ1)*DBDR(3) + VB*D2BDRJX(3)
                     HESS(J3,J6-2)   = HESS(J3,J6-2) + DUMMY
                     HESS(J6-2,J3)   = HESS(J6-2,J3) + DUMMY
!     xi,pj2
                     DUMMY = D2VDR2*DRIJDJ2*RSS(1) - DVDR*DR2(J8,1) + DVRDR*R2*DRIJDJ2*RSS(1) + (DVRDA*DADPJ2 + DVRDB*DBDPJ2      &
                           + DVRDG*DGDPJ2)*NR(1) - VR*DR2(J8,1)/ABSRIJ - VR*R2*DRIJDJ2*RSS(1)/ABSRIJ + (DVADR*DRIJDJ2/ABSRIJ      &
                           + DVADB*DBDPJ2)*DADR(1) + VA*D2ADRJY(1) + (DVBDR*DRIJDJ2/ABSRIJ + DVBDA*DADPJ2)*DBDR(1) + VB*D2BDRJY(1)
                     HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
                     HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
                     DUMMY = D2VDR2*DRIJDJ2*RSS(2) - DVDR*DR2(J8,2) + DVRDR*R2*DRIJDJ2*RSS(2) + (DVRDA*DADPJ2 + DVRDB*DBDPJ2      &
                           + DVRDG*DGDPJ2)*NR(2) - VR*DR2(J8,2)/ABSRIJ - VR*R2*DRIJDJ2*RSS(2)/ABSRIJ + (DVADR*DRIJDJ2/ABSRIJ      &
                           + DVADB*DBDPJ2)*DADR(2) + VA*D2ADRJY(2) + (DVBDR*DRIJDJ2/ABSRIJ + DVBDA*DADPJ2)*DBDR(2) + VB*D2BDRJY(2)
                     HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
                     HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
                     DUMMY = D2VDR2*DRIJDJ2*RSS(3) - DVDR*DR2(J8,3) + DVRDR*R2*DRIJDJ2*RSS(3) + (DVRDA*DADPJ2 + DVRDB*DBDPJ2      &
                           + DVRDG*DGDPJ2)*NR(3) - VR*DR2(J8,3)/ABSRIJ - VR*R2*DRIJDJ2*RSS(3)/ABSRIJ + (DVADR*DRIJDJ2/ABSRIJ      &
                           + DVADB*DBDPJ2)*DADR(3) + VA*D2ADRJY(3) + (DVBDR*DRIJDJ2/ABSRIJ + DVBDA*DADPJ2)*DBDR(3) + VB*D2BDRJY(3)
                     HESS(J3,J6-1)  = HESS(J3,J6-1) + DUMMY
                     HESS(J6-1,J3)  = HESS(J6-1,J3) + DUMMY
!     xi,pj3
                     DUMMY = D2VDR2*DRIJDJ3*RSS(1) - DVDR*DR3(J8,1) + DVRDR*R2*DRIJDJ3*RSS(1) + (DVRDA*DADPJ3 + DVRDB*DBDPJ3      &
                           + DVRDG*DGDPJ3)*NR(1) - VR*DR3(J8,1)/ABSRIJ - VR*R2*DRIJDJ3*RSS(1)/ABSRIJ + (DVADR*DRIJDJ3/ABSRIJ      &
                           + DVADB*DBDPJ3)*DADR(1) + VA*D2ADRJZ(1) + (DVBDR*DRIJDJ3/ABSRIJ + DVBDA*DADPJ3)*DBDR(1) + VB*D2BDRJZ(1)
                     HESS(J3-2,J6)   = HESS(J3-2,J6) + DUMMY
                     HESS(J6,J3-2)   = HESS(J6,J3-2) + DUMMY
!     yi,pj3
                     DUMMY = D2VDR2*DRIJDJ3*RSS(2) - DVDR*DR3(J8,2) + DVRDR*R2*DRIJDJ3*RSS(2) + (DVRDA*DADPJ3 + DVRDB*DBDPJ3      &
                           + DVRDG*DGDPJ3)*NR(2) - VR*DR3(J8,2)/ABSRIJ - VR*R2*DRIJDJ3*RSS(2)/ABSRIJ + (DVADR*DRIJDJ3/ABSRIJ      &
                           + DVADB*DBDPJ3)*DADR(2) + VA*D2ADRJZ(2) + (DVBDR*DRIJDJ3/ABSRIJ + DVBDA*DADPJ3)*DBDR(2) + VB*D2BDRJZ(2)
                     HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
                     HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
                     DUMMY = D2VDR2*DRIJDJ3*RSS(3) - DVDR*DR3(J8,3) + DVRDR*R2*DRIJDJ3*RSS(3) + (DVRDA*DADPJ3 + DVRDB*DBDPJ3      &
                           + DVRDG*DGDPJ3)*NR(3) - VR*DR3(J8,3)/ABSRIJ - VR*R2*DRIJDJ3*RSS(3)/ABSRIJ + (DVADR*DRIJDJ3/ABSRIJ      &
                           + DVADB*DBDPJ3)*DADR(3) + VA*D2ADRJZ(3) + (DVBDR*DRIJDJ3/ABSRIJ + DVBDA*DADPJ3)*DBDR(3) + VB*D2BDRJZ(3)
                     HESS(J3,J6)     = HESS(J3,J6) + DUMMY
                     HESS(J6,J3)     = HESS(J6,J3) + DUMMY
!     pi1,pj2
                     DUMMY = D2VDR2*DRIJDJ2*DRIJDI1 - DVDR*DOT_PRODUCT(DR2(J8,:),DR1(J7,:))         &
                           + (DVRDR*DRIJDJ2/ABSRIJ + DVRDA*DADPJ2 + DVRDB*DBDPJ2 + DVRDG*DGDPJ2)*DRIJDI1/ABSRIJ         &
                           - (VR*DOT_PRODUCT(DR2(J8,:),DR1(J7,:)) + VR*DRIJDJ2*DRIJDI1*R2)/ABSRIJ   &
                           + (DVADR*DRIJDJ2/ABSRIJ + DVADB*DBDPJ2)*DADPI1 + VA*D2ADIXJY  &
                           + (DVBDR*DRIJDJ2/ABSRIJ + DVBDA*DADPJ2)*DBDPI1 + VB*D2BDIXJY + DVGDR*DRIJDJ2*DGDPI1/ABSRIJ   &
                           + VG*DOT_PRODUCT(DE1(J7,:),DE2(J8,:))
                     HESS(J5-2,J6-1) = HESS(J5-2,J6-1) + DUMMY
                     HESS(J6-1,J5-2) = HESS(J6-1,J5-2) + DUMMY
!     pi2,pj3
                     DUMMY = D2VDR2*DRIJDJ3*DRIJDI2 - DVDR*DOT_PRODUCT(DR3(J8,:),DR2(J7,:))         &
                           + (DVRDR*DRIJDJ3/ABSRIJ + DVRDA*DADPJ3 + DVRDB*DBDPJ3 + DVRDG*DGDPJ3)*DRIJDI2/ABSRIJ         &
                           - (VR*DOT_PRODUCT(DR3(J8,:),DR2(J7,:)) + VR*DRIJDJ3*DRIJDI2*R2)/ABSRIJ   &
                           + (DVADR*DRIJDJ3/ABSRIJ + DVADB*DBDPJ3)*DADPI2 + VA*D2ADIYJZ  &
                           + (DVBDR*DRIJDJ3/ABSRIJ + DVBDA*DADPJ3)*DBDPI2 + VB*D2BDIYJZ + DVGDR*DRIJDJ3*DGDPI2/ABSRIJ   &
                           + VG*DOT_PRODUCT(DE2(J7,:),DE3(J8,:))
                     HESS(J5-1,J6)   = HESS(J5-1,J6) + DUMMY
                     HESS(J6,J5-1)   = HESS(J6,J5-1) + DUMMY
!     pi3,pj1
                     DUMMY = D2VDR2*DRIJDJ1*DRIJDI3 - DVDR*DOT_PRODUCT(DR1(J8,:),DR3(J7,:))         &
                           + (DVRDR*DRIJDJ1/ABSRIJ + DVRDA*DADPJ1 + DVRDB*DBDPJ1 + DVRDG*DGDPJ1)*DRIJDI3/ABSRIJ         &
                           - (VR*DOT_PRODUCT(DR1(J8,:),DR3(J7,:)) + VR*DRIJDJ1*DRIJDI3*R2)/ABSRIJ   &
                           + (DVADR*DRIJDJ1/ABSRIJ + DVADB*DBDPJ1)*DADPI3 + VA*D2ADIZJX  &
                           + (DVBDR*DRIJDJ1/ABSRIJ + DVBDA*DADPJ1)*DBDPI3 + VB*D2BDIZJX + DVGDR*DRIJDJ1*DGDPI3/ABSRIJ   &
                           + VG*DOT_PRODUCT(DE3(J7,:),DE1(J8,:))
                     HESS(J5,J6-2)   = HESS(J5,J6-2) + DUMMY
                     HESS(J6-2,J5)   = HESS(J6-2,J5) + DUMMY

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

      ENDDO

      ENERGY = FCTR*ENERGY
      G(:)   = FCTR*G(:)

      IF (EFIELDT) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3

            DO I = 1, NRBSITES

               J7 = NRBSITES*(J1-1) + I
               EI(:) = E(J7,:)

               ENERGY = ENERGY - DPMU(I)*EFIELD*EI(3)

               IF (GTEST) THEN

                  G(J5-2) = G(J5-2) - DPMU(I)*EFIELD*DE1(J7,3)
                  G(J5-1) = G(J5-1) - DPMU(I)*EFIELD*DE2(J7,3)
                  G(J5)   = G(J5)   - DPMU(I)*EFIELD*DE3(J7,3)

               ENDIF

               IF (STEST) THEN
!     pi1,pi1
                  HESS(J5-2,J5-2) = HESS(J5-2,J5-2) - DPMU(I)*EFIELD*D2E1(J7,3)
!     pi2,pi2  
                  HESS(J5-1,J5-1) = HESS(J5-1,J5-1) - DPMU(I)*EFIELD*D2E2(J7,3)
!     pi3,pi3 
                  HESS(J5,J5)     = HESS(J5,J5)     - DPMU(I)*EFIELD*D2E3(J7,3) 
!     pi1,pi2
                  DUMMY           =-DPMU(I)*EFIELD*D2E12(J7,3)
                  HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
                  HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
                  DUMMY           =-DPMU(I)*EFIELD*D2E23(J7,3)
                  HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
                  HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
                  DUMMY           =-DPMU(I)*EFIELD*D2E31(J7,3)
                  HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
                  HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

               ENDIF

            ENDDO

         ENDDO

      ENDIF

      END SUBROUTINE MSSTOCKGH 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMULTSTOCK()

      USE COMMONS, ONLY: NRBSITES, RBSITE, RBSTLA, DPMU

      IMPLICIT NONE

      IF (NRBSITES == 4) THEN

         RBSITE(1,:) = (/0.D0, 1.D0/DSQRT(3.D0), 0.D0/)
         RBSITE(2,:) = (/0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
         RBSITE(3,:) = (/-0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
         RBSITE(4,:) = (/ 0.D0, 0.D0, 0.D0/)

         RBSTLA(1,:) = (/-1.D0, 0.D0, 0.D0/)
         RBSTLA(2,:) = (/0.5D0, 0.5D0*DSQRT(3.D0), 0.D0/)
         RBSTLA(3,:) = (/0.5D0, -0.5D0*DSQRT(3.D0), 0.D0/)
         RBSTLA(4,:) = (/0.D0, 0.D0, 1.D0/)

      ELSEIF (NRBSITES == 3) THEN

         RBSITE(1,:) = (/0.D0, 1.D0/DSQRT(3.D0), 0.D0/)
         RBSITE(2,:) = (/0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
         RBSITE(3,:) = (/-0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)

         RBSTLA(1,:) = (/-1.D0, 0.D0, 0.D0/)
         RBSTLA(2,:) = (/0.5D0, 0.5D0*DSQRT(3.D0), 0.D0/)
         RBSTLA(3,:) = (/0.5D0, -0.5D0*DSQRT(3.D0), 0.D0/)

      ENDIF

      END SUBROUTINE DEFMULTSTOCK
