      SUBROUTINE PAHAGH (X, G, ENERGY, GTEST, SECT)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, RBSTLA, STCHRG
      USE KEY, ONLY: RHOCC0, RHOCC10, RHOCC20,  RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, RHOCH20,   &
                     ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, CCKJ, NCARBON, NTSITES 

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET, FCT(6), J2INT  
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R4, R6, R12, ABSRIJ, RIJSQ, ENERGY1, ENERGY2, ENERGY3, DUMMY
      DOUBLE PRECISION :: RI(3), RJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3), FRIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: R(NTSITES,3), E(NTSITES,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
      DOUBLE PRECISION :: DR1(NTSITES,3), DR2(NTSITES,3), DR3(NTSITES,3) 
      DOUBLE PRECISION :: D2R1(NTSITES,3), D2R2(NTSITES,3), D2R3(NTSITES,3) 
      DOUBLE PRECISION :: D2R12(NTSITES,3), D2R23(NTSITES,3), D2R31(NTSITES,3) 
      DOUBLE PRECISION :: DE1(NTSITES,3), DE2(NTSITES,3), DE3(NTSITES,3)
      DOUBLE PRECISION :: D2E1(NTSITES,3), D2E2(NTSITES,3), D2E3(NTSITES,3)
      DOUBLE PRECISION :: D2E12(NTSITES,3), D2E23(NTSITES,3), D2E31(NTSITES,3)
      DOUBLE PRECISION :: DCADR(3), DCBDR(3)
      DOUBLE PRECISION :: RHOCC, RHOHH, RHOCH, COSTA, COSTB, DMPFCT, DDMPDR, D2DDR2, EXPFCT 
      DOUBLE PRECISION :: DRIJDPI(3), DRIJDPJ(3), DCADPI(3), DCBDPI(3), DCADPJ(3), DCBDPJ(3)
      DOUBLE PRECISION :: DRHODR(3), DRHODPI(3), DRHODPJ(3) 
      DOUBLE PRECISION :: D2CARX(3), D2CARY(3), D2CARZ(3), D2CBRX(3), D2CBRY(3), D2CBRZ(3)
      DOUBLE PRECISION :: D2RHOX(3), D2RHOY(3), D2RHOZ(3)
      DOUBLE PRECISION :: D2RHODIXIX, D2RHODIXIY, D2RHODIXIZ, D2RHODIYIY, D2RHODIYIZ, D2RHODIZIZ
      DOUBLE PRECISION :: D2RHODIXJX, D2RHODIYJY, D2RHODIZJZ, D2RHODIXJY, D2RHODIYJZ, D2RHODIXJZ
      DOUBLE PRECISION :: D2RHORDIX(3), D2RHORDIY(3), D2RHORDIZ(3), D2RHORDJX(3), D2RHORDJY(3), D2RHORDJZ(3)
      DOUBLE PRECISION :: D2CADIXIX, D2CADIXIY, D2CADIXIZ, D2CADIYIY, D2CADIYIZ, D2CADIZIZ
      DOUBLE PRECISION :: D2CBDIXIX, D2CBDIXIY, D2CBDIXIZ, D2CBDIYIY, D2CBDIYIZ, D2CBDIZIZ
      DOUBLE PRECISION :: D2CADIXJX, D2CADIXJY, D2CADIXJZ, D2CADIYJY, D2CADIYJZ, D2CADIZJZ
      DOUBLE PRECISION :: D2CBDIXJX, D2CBDIXJY, D2CBDIXJZ, D2CBDIYJY, D2CBDIYJZ, D2CBDIZJZ
      DOUBLE PRECISION :: D2CARDIX(3), D2CARDIY(3), D2CARDIZ(3), D2CARDJX(3), D2CARDJY(3), D2CARDJZ(3)
      DOUBLE PRECISION :: D2CBRDIX(3), D2CBRDIY(3), D2CBRDIZ(3), D2CBRDJX(3), D2CBRDJY(3), D2CBRDJZ(3)
      DOUBLE PRECISION :: D2RDIXIX, D2RDIXIY, D2RDIXIZ, D2RDIYIY, D2RDIYIZ, D2RDIZIZ
      DOUBLE PRECISION :: D2RDIXJX, D2RDIYJY, D2RDIZJZ, D2RDIXJY, D2RDIYJZ, D2RDIXJZ
      DOUBLE PRECISION :: D2RDJXJX, D2RDJYJY, D2RDJZJZ
      DOUBLE PRECISION, PARAMETER :: B = 1.6485D0
      LOGICAL          :: GTEST, SECT
      DOUBLE PRECISION :: DVDR, D2VDR2
      DOUBLE PRECISION :: RHO, RHO0, RHO10, RHO01, RHO20, RHO02, ALPHA, DC6, FCTR 
      DOUBLE PRECISION :: FVAL, DF

      FCT(1) = 1; FCT(2) = 2; FCT(3) = 6; FCT(4) = 24; FCT(5) = 120; FCT(6) = 720
      ENERGY = 0.D0; ENERGY1 = 0.D0; ENERGY2 = 0.D0; ENERGY3 = 0.D0
      J2INT = 1

      IF (GTEST) G(:) = 0.D0
      IF (SECT) THEN
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

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, SECT)

         DO J2 = 1, NRBSITES

            J4      = NRBSITES*(J1-1) + J2
            R(J4,:) = RI(:) + MATMUL(RMI(:,:),RBSITE(J2,:))
            E(J4,:) = MATMUL(RMI(:,:),RBSTLA(J2,:))

            IF (GTEST .OR. SECT) THEN

               DR1(J4,:) = MATMUL(DRMI1(:,:),RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),RBSITE(J2,:))

               DE1(J4,:) = MATMUL(DRMI1(:,:),RBSTLA(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBSTLA(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBSTLA(J2,:))

            ENDIF

            IF (SECT) THEN

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

      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3

         RI(:)  = X(J3-2:J3)

         DO I = 1, NRBSITES

            J7    = NRBSITES*(J1-1) + I
            EI(:) = E(J7,:)
               
            IF (.NOT. SECT) J2INT = J1 + 1

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
                  R2     = 1.D0/R2
                  R6     = R2*R2*R2

!     CALCULATE THE DAMPING FACTOR

                  DMPFCT = 1.D0
                  DDMPDR = B
                  D2DDR2 = B

                  DO K = 1, 6

                     DMPFCT = DMPFCT + (B*ABSRIJ)**K/DBLE(FCT(K))
                     IF (K > 1) DDMPDR = DDMPDR + (B**K)*(ABSRIJ)**(K-1)/DBLE(FCT(K-1))
                     IF (K > 2) D2DDR2 = D2DDR2 + (B**K)*(ABSRIJ)**(K-2)/DBLE(FCT(K-2))

                  END DO

                  EXPFCT = DEXP(-B*ABSRIJ)
                  D2DDR2 = (-B*B*EXPFCT*DMPFCT + 2.D0*B*EXPFCT*DDMPDR - EXPFCT*D2DDR2)/ABSRIJ 
                  DDMPDR = (B*EXPFCT*DMPFCT - EXPFCT*DDMPDR)/ABSRIJ
                  DMPFCT = 1.D0 - EXPFCT*DMPFCT

!     NOW CALCULATE RHOAB

                  COSTA      =-DOT_PRODUCT(NR(:),EI(:))
                  COSTB      = DOT_PRODUCT(NR(:),EJ(:))
 
                  IF (I <= NCARBON .AND. J <= NCARBON) THEN

                     RHO0 = RHOCC0; RHO10 = RHOCC10; RHO01 = RHOCC10; RHO20 = RHOCC20; RHO02 = RHOCC20; ALPHA = ALPHACC; DC6 = DC6CC

                  ELSEIF (I > NCARBON .AND. J > NCARBON) THEN

                     RHO0 = RHOHH0; RHO10 = RHOHH10; RHO01 = RHOHH10; RHO20 = RHOHH20; RHO02 = RHOHH20; ALPHA = ALPHAHH; DC6 = DC6HH

                  ELSE IF (I <= NCARBON .AND. J > NCARBON) THEN

                     RHO0 = RHOCH0; RHO10 = RHOC10H; RHO01 = RHOCH10; RHO20 = RHOC20H; RHO02 = RHOCH20; ALPHA = ALPHACH; DC6 = DC6CH

                  ELSE !IF(I > NCARBON .AND. J <= NCARBON) THEN

                     RHO0 = RHOCH0; RHO10 = RHOCH10; RHO01 = RHOC10H; RHO20 = RHOCH20; RHO02 = RHOC20H; ALPHA = ALPHACH; DC6 = DC6CH
                     

                  ENDIF

                  RHO     = RHO0 + RHO10*COSTA + RHO01*COSTB + RHO20*(1.5D0*COSTA*COSTA-0.5D0) + RHO02*(1.5D0*COSTB*COSTB-0.5D0)
                  EXPFCT  = KKJ*DEXP(-ALPHA*(ABSRIJ - RHO))
                  ENERGY1 = ENERGY1 + EXPFCT
                  ENERGY2 = ENERGY2 - DC6*DMPFCT*R6
                  ENERGY3 = ENERGY3 + CCKJ*STCHRG(I)*STCHRG(J)/ABSRIJ

                  IF (GTEST) THEN

                     DVDR   = 6.D0*DC6*R6*R2*DMPFCT - DC6*R6*DDMPDR - CCKJ*STCHRG(I)*STCHRG(J)*R2/ABSRIJ

                     DCADR(:)   =-EI(:)/ABSRIJ - COSTA*R2*RSS(:)
                     DCBDR(:)   = EJ(:)/ABSRIJ - COSTB*R2*RSS(:)

                     DRIJDPI(1) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                     DRIJDPI(2) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                     DRIJDPI(3) = DOT_PRODUCT(RSS(:),DR3(J7,:))

                     DRIJDPJ(1) =-DOT_PRODUCT(RSS(:),DR1(J8,:))
                     DRIJDPJ(2) =-DOT_PRODUCT(RSS(:),DR2(J8,:))
                     DRIJDPJ(3) =-DOT_PRODUCT(RSS(:),DR3(J8,:))

                     DCADPI(1)  =-DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE1(J7,:)) - COSTA*R2*DRIJDPI(1)
                     DCADPI(2)  =-DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE2(J7,:)) - COSTA*R2*DRIJDPI(2)
                     DCADPI(3)  =-DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(NR(:),DE3(J7,:)) - COSTA*R2*DRIJDPI(3)

                     DCBDPI(1)  = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(1)
                     DCBDPI(2)  = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(2)
                     DCBDPI(3)  = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - COSTB*R2*DRIJDPI(3)

                     DCADPJ(1)  = DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(1)
                     DCADPJ(2)  = DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(2)
                     DCADPJ(3)  = DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - COSTA*R2*DRIJDPJ(3)

                     DCBDPJ(1)  =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE1(J8,:)) - COSTB*R2*DRIJDPJ(1)
                     DCBDPJ(2)  =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE2(J8,:)) - COSTB*R2*DRIJDPJ(2)
                     DCBDPJ(3)  =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ + DOT_PRODUCT(NR(:),DE3(J8,:)) - COSTB*R2*DRIJDPJ(3)

                     DRHODR(:)  = (RHO10 + 3.D0*RHO20*COSTA)*DCADR(:) + (RHO01 + 3.D0*RHO02*COSTB)*DCBDR(:)
                     DRHODPI(:) = (RHO10 + 3.D0*RHO20*COSTA)*DCADPI(:) + (RHO01 + 3.D0*RHO02*COSTB)*DCBDPI(:)
                     DRHODPJ(:) = (RHO10 + 3.D0*RHO20*COSTA)*DCADPJ(:) + (RHO01 + 3.D0*RHO02*COSTB)*DCBDPJ(:)

                     FRIJ(:) = ALPHA*EXPFCT*(-NR(:) + DRHODR(:))
                     TIJ(:)  = ALPHA*EXPFCT*(-DRIJDPI(:)/ABSRIJ + DRHODPI(:))
                     TJI(:)  = ALPHA*EXPFCT*(-DRIJDPJ(:)/ABSRIJ + DRHODPJ(:))

                     G(J3-2:J3) = G(J3-2:J3) + DVDR*RSS(:) + FRIJ(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*RSS(:) - FRIJ(:)

                     G(J5-2:J5) = G(J5-2:J5) + DVDR*DRIJDPI(:) + TIJ(:)
                     G(J6-2:J6) = G(J6-2:J6) + DVDR*DRIJDPJ(:) + TJI(:)

!                     FVAL = DRHODPI(2) 

                  ENDIF

                  IF (SECT) THEN

                     D2VDR2 =-48.D0*DC6*R6*R2*R2*DMPFCT + 6.D0*DC6*R6*R2*DDMPDR + 7.D0*DC6*R6*R2*DDMPDR - DC6*R6*D2DDR2/ABSRIJ    &
                            + 3.D0*CCKJ*STCHRG(I)*STCHRG(J)*R2*R2/ABSRIJ

                     D2CARX(:) = R2*RSS(1)/ABSRIJ*EI(:) - R2*DCADR(1)*RSS(:) - R2*COSTA*(/1.D0, 0.D0, 0.D0/)  &
                               + 2.D0*R2*COSTA*NR(1)/ABSRIJ*RSS(:)
                     D2CARY(:) = R2*RSS(2)/ABSRIJ*EI(:) - R2*DCADR(2)*RSS(:) - R2*COSTA*(/0.D0, 1.D0, 0.D0/)  &
                               + 2.D0*R2*COSTA*NR(2)/ABSRIJ*RSS(:)
                     D2CARZ(:) = R2*RSS(3)/ABSRIJ*EI(:) - R2*DCADR(3)*RSS(:) - R2*COSTA*(/0.D0, 0.D0, 1.D0/)  &
                               + 2.D0*R2*COSTA*NR(3)/ABSRIJ*RSS(:)
           
                     D2CARDIX(:) = R2*DRIJDPI(1)/ABSRIJ*EI(:) - DE1(J7,:)/ABSRIJ - R2*DCADPI(1)*RSS(:)        &
                                 - R2*COSTA*DR1(J7,:) + 2.D0*R2*R2*COSTA*DRIJDPI(1)*RSS(:)
                     D2CARDIY(:) = R2*DRIJDPI(2)/ABSRIJ*EI(:) - DE2(J7,:)/ABSRIJ - R2*DCADPI(2)*RSS(:)        &
                                 - R2*COSTA*DR2(J7,:) + 2.D0*R2*R2*COSTA*DRIJDPI(2)*RSS(:)
                     D2CARDIZ(:) = R2*DRIJDPI(3)/ABSRIJ*EI(:) - DE3(J7,:)/ABSRIJ - R2*DCADPI(3)*RSS(:)        &
                                 - R2*COSTA*DR3(J7,:) + 2.D0*R2*R2*COSTA*DRIJDPI(3)*RSS(:)

                     D2CARDJX(:) = R2*DRIJDPJ(1)/ABSRIJ*EI(:) - R2*DCADPJ(1)*RSS(:) + R2*COSTA*DR1(J8,:)      &
                                 + 2.D0*R2*R2*COSTA*DRIJDPJ(1)*RSS(:)
                     D2CARDJY(:) = R2*DRIJDPJ(2)/ABSRIJ*EI(:) - R2*DCADPJ(2)*RSS(:) + R2*COSTA*DR2(J8,:)      &
                                 + 2.D0*R2*R2*COSTA*DRIJDPJ(2)*RSS(:)
                     D2CARDJZ(:) = R2*DRIJDPJ(3)/ABSRIJ*EI(:) - R2*DCADPJ(3)*RSS(:) + R2*COSTA*DR3(J8,:)      &
                                 + 2.D0*R2*R2*COSTA*DRIJDPJ(3)*RSS(:)

                     D2CBRX(:) =-R2*RSS(1)/ABSRIJ*EJ(:) - R2*DCBDR(1)*RSS(:) - R2*COSTB*(/1.D0, 0.D0, 0.D0/)  &
                               + 2.D0*R2*COSTB*NR(1)/ABSRIJ*RSS(:)
                     D2CBRY(:) =-R2*RSS(2)/ABSRIJ*EJ(:) - R2*DCBDR(2)*RSS(:) - R2*COSTB*(/0.D0, 1.D0, 0.D0/)  &
                               + 2.D0*R2*COSTB*NR(2)/ABSRIJ*RSS(:)
                     D2CBRZ(:) =-R2*RSS(3)/ABSRIJ*EJ(:) - R2*DCBDR(3)*RSS(:) - R2*COSTB*(/0.D0, 0.D0, 1.D0/)  &
                               + 2.D0*R2*COSTB*NR(3)/ABSRIJ*RSS(:)

                     D2CBRDIX(:) =-R2*DRIJDPI(1)/ABSRIJ*EJ(:) - R2*DCBDPI(1)*RSS(:) - R2*COSTB*DR1(J7,:) &
                               + 2.D0*R2*R2*COSTB*DRIJDPI(1)*RSS(:)
                     D2CBRDIY(:) =-R2*DRIJDPI(2)/ABSRIJ*EJ(:) - R2*DCBDPI(2)*RSS(:) - R2*COSTB*DR2(J7,:) &
                               + 2.D0*R2*R2*COSTB*DRIJDPI(2)*RSS(:)
                     D2CBRDIZ(:) =-R2*DRIJDPI(3)/ABSRIJ*EJ(:) - R2*DCBDPI(3)*RSS(:) - R2*COSTB*DR3(J7,:) &
                               + 2.D0*R2*R2*COSTB*DRIJDPI(3)*RSS(:)

                     D2CBRDJX(:) =-R2*DRIJDPJ(1)/ABSRIJ*EJ(:) + DE1(J8,:)/ABSRIJ - R2*DCBDPJ(1)*RSS(:)        &
                                 + R2*COSTB*DR1(J8,:) + 2.D0*R2*R2*COSTB*DRIJDPJ(1)*RSS(:)
                     D2CBRDJY(:) =-R2*DRIJDPJ(2)/ABSRIJ*EJ(:) + DE2(J8,:)/ABSRIJ - R2*DCBDPJ(2)*RSS(:)        &
                                 + R2*COSTB*DR2(J8,:) + 2.D0*R2*R2*COSTB*DRIJDPJ(2)*RSS(:)
                     D2CBRDJZ(:) =-R2*DRIJDPJ(3)/ABSRIJ*EJ(:) + DE3(J8,:)/ABSRIJ - R2*DCBDPJ(3)*RSS(:)        &
                                 + R2*COSTB*DR3(J8,:) + 2.D0*R2*R2*COSTB*DRIJDPJ(3)*RSS(:)

                     D2RDIXIX = (DOT_PRODUCT(DR1(J7,:),DR1(J7,:)) + DOT_PRODUCT(RSS(:),D2R1(J7,:)) - R2*DRIJDPI(1)**2)/ABSRIJ
                     D2RDIYIY = (DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) + DOT_PRODUCT(RSS(:),D2R2(J7,:)) - R2*DRIJDPI(2)**2)/ABSRIJ
                     D2RDIZIZ = (DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) + DOT_PRODUCT(RSS(:),D2R3(J7,:)) - R2*DRIJDPI(3)**2)/ABSRIJ

                     D2RDIXIY = (DOT_PRODUCT(DR1(J7,:),DR2(J7,:)) + DOT_PRODUCT(RSS(:),D2R12(J7,:)) &
  &                                                               - R2*DRIJDPI(1)*DRIJDPI(2))/ABSRIJ
                     D2RDIYIZ = (DOT_PRODUCT(DR2(J7,:),DR3(J7,:)) + DOT_PRODUCT(RSS(:),D2R23(J7,:)) &
  &                                                               - R2*DRIJDPI(2)*DRIJDPI(3))/ABSRIJ
                     D2RDIXIZ = (DOT_PRODUCT(DR3(J7,:),DR1(J7,:)) + DOT_PRODUCT(RSS(:),D2R31(J7,:)) &
  &                                                               - R2*DRIJDPI(3)*DRIJDPI(1))/ABSRIJ

                     D2RDJXJX = (DOT_PRODUCT(DR1(J8,:),DR1(J8,:)) - DOT_PRODUCT(RSS(:),D2R1(J8,:)) - R2*DRIJDPJ(1)**2)/ABSRIJ
                     D2RDJYJY = (DOT_PRODUCT(DR2(J8,:),DR2(J8,:)) - DOT_PRODUCT(RSS(:),D2R2(J8,:)) - R2*DRIJDPJ(2)**2)/ABSRIJ
                     D2RDJZJZ = (DOT_PRODUCT(DR3(J8,:),DR3(J8,:)) - DOT_PRODUCT(RSS(:),D2R3(J8,:)) - R2*DRIJDPJ(3)**2)/ABSRIJ

                     D2RDIXJX = (-DOT_PRODUCT(DR1(J7,:),DR1(J8,:)) - R2*DRIJDPI(1)*DRIJDPJ(1))/ABSRIJ
                     D2RDIYJY = (-DOT_PRODUCT(DR2(J7,:),DR2(J8,:)) - R2*DRIJDPI(2)*DRIJDPJ(2))/ABSRIJ
                     D2RDIZJZ = (-DOT_PRODUCT(DR3(J7,:),DR3(J8,:)) - R2*DRIJDPI(3)*DRIJDPJ(3))/ABSRIJ

                     D2RDIXJY = (-DOT_PRODUCT(DR1(J7,:),DR2(J8,:)) - R2*DRIJDPI(1)*DRIJDPJ(2))/ABSRIJ
                     D2RDIYJZ = (-DOT_PRODUCT(DR2(J7,:),DR3(J8,:)) - R2*DRIJDPI(2)*DRIJDPJ(3))/ABSRIJ
                     D2RDIXJZ = (-DOT_PRODUCT(DR3(J7,:),DR1(J8,:)) - R2*DRIJDPI(3)*DRIJDPJ(1))/ABSRIJ
                     
                     D2CADIXIX     =-DOT_PRODUCT(D2R1(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(DR1(J7,:),DE1(J7,:))/ABSRIJ     &
                                   + R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDPI(1)/ABSRIJ - DOT_PRODUCT(DR1(J7,:),DE1(J7,:))/ABSRIJ  &
                                   - DOT_PRODUCT(NR(:),D2E1(J7,:)) + R2*DOT_PRODUCT(RSS(:),DE1(J7,:))*DRIJDPI(1)/ABSRIJ &
                                   - DCADPI(1)*R2*DRIJDPI(1) + COSTA*R2*R2*DRIJDPI(1)*DRIJDPI(1) - COSTA*D2RDIXIX/ABSRIJ

                     D2CADIXIY     =-DOT_PRODUCT(D2R12(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(DR2(J7,:),DE1(J7,:))/ABSRIJ    &
                                   + R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDPI(2)/ABSRIJ - DOT_PRODUCT(DR2(J7,:),DE1(J7,:))/ABSRIJ  &
                                   - DOT_PRODUCT(NR(:),D2E12(J7,:)) + R2*DOT_PRODUCT(RSS(:),DE1(J7,:))*DRIJDPI(2)/ABSRIJ     &
                                   - DCADPI(2)*R2*DRIJDPI(1) + 2.D0*COSTA*R2*R2*DRIJDPI(1)*DRIJDPI(2)         &
                                   - COSTA*R2*DOT_PRODUCT(DR2(J7,:),DR1(J7,:)) - COSTA*R2*DOT_PRODUCT(RSS(:),D2R12(J7,:))
                     
                     D2CADIXIZ     =-DOT_PRODUCT(D2R31(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(DR3(J7,:),DE1(J7,:))/ABSRIJ    &
                                   + R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDPI(3)/ABSRIJ - DOT_PRODUCT(DR3(J7,:),DE1(J7,:))/ABSRIJ  &
                                   - DOT_PRODUCT(NR(:),D2E31(J7,:)) + R2*DOT_PRODUCT(RSS(:),DE1(J7,:))*DRIJDPI(3)/ABSRIJ     &
                                   - DCADPI(3)*R2*DRIJDPI(1) + 2.D0*COSTA*R2*R2*DRIJDPI(1)*DRIJDPI(3)         &
                                   - COSTA*R2*DOT_PRODUCT(DR3(J7,:),DR1(J7,:)) - COSTA*R2*DOT_PRODUCT(RSS(:),D2R31(J7,:))

                     D2CADIYIY     =-DOT_PRODUCT(D2R2(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(DR2(J7,:),DE2(J7,:))/ABSRIJ     &
                                   + R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDPI(2)/ABSRIJ - DOT_PRODUCT(DR2(J7,:),DE2(J7,:))/ABSRIJ  &
                                   - DOT_PRODUCT(NR(:),D2E2(J7,:)) + R2*DOT_PRODUCT(RSS(:),DE2(J7,:))*DRIJDPI(2)/ABSRIJ &
                                   - DCADPI(2)*R2*DRIJDPI(2) + 2.D0*COSTA*R2*R2*DRIJDPI(2)*DRIJDPI(2)         &
                                   - COSTA*R2*DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) - COSTA*R2*DOT_PRODUCT(RSS(:),D2R2(J7,:))

                     D2CADIYIZ     =-DOT_PRODUCT(D2R23(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(DR3(J7,:),DE2(J7,:))/ABSRIJ     &
                                   + R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDPI(2)/ABSRIJ - DOT_PRODUCT(DR3(J7,:),DE2(J7,:))/ABSRIJ  &
                                   - DOT_PRODUCT(NR(:),D2E23(J7,:)) + R2*DOT_PRODUCT(RSS(:),DE3(J7,:))*DRIJDPI(2)/ABSRIJ &
                                   - DCADPI(3)*R2*DRIJDPI(2) + 2.D0*COSTA*R2*R2*DRIJDPI(3)*DRIJDPI(2)         &
                                   - COSTA*R2*DOT_PRODUCT(DR3(J7,:),DR2(J7,:)) - COSTA*R2*DOT_PRODUCT(RSS(:),D2R23(J7,:))

                     D2CADIZIZ     =-DOT_PRODUCT(D2R3(J7,:),EI(:))/ABSRIJ - DOT_PRODUCT(DR3(J7,:),DE3(J7,:))/ABSRIJ     &
                                   + R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDPI(3)/ABSRIJ - DOT_PRODUCT(DR3(J7,:),DE3(J7,:))/ABSRIJ  &
                                   - DOT_PRODUCT(NR(:),D2E3(J7,:)) + R2*DOT_PRODUCT(RSS(:),DE3(J7,:))*DRIJDPI(3)/ABSRIJ &
                                   - DCADPI(3)*R2*DRIJDPI(3) + 2.D0*COSTA*R2*R2*DRIJDPI(3)*DRIJDPI(3)         &
                                   - COSTA*R2*DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) - COSTA*R2*DOT_PRODUCT(RSS(:),D2R3(J7,:))

                     D2CADIXJX     = R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDPJ(1)/ABSRIJ + DOT_PRODUCT(DR1(J8,:),DE1(J7,:))/ABSRIJ  &
                                   + R2*DOT_PRODUCT(RSS(:),DE1(J7,:))*DRIJDPJ(1)/ABSRIJ &
                                   - DCADPJ(1)*R2*DRIJDPI(1) + COSTA*R2*R2*DRIJDPI(1)*DRIJDPJ(1) - COSTA*D2RDIXJX/ABSRIJ
                     D2CADIYJY     = R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDPJ(2)/ABSRIJ + DOT_PRODUCT(DR2(J8,:),DE2(J7,:))/ABSRIJ  &
                                   + R2*DOT_PRODUCT(RSS(:),DE2(J7,:))*DRIJDPJ(2)/ABSRIJ &
                                   - DCADPJ(2)*R2*DRIJDPI(2) + COSTA*R2*R2*DRIJDPI(2)*DRIJDPJ(2) - COSTA*D2RDIYJY/ABSRIJ
                     D2CADIZJZ     = R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDPJ(3)/ABSRIJ + DOT_PRODUCT(DR3(J8,:),DE3(J7,:))/ABSRIJ  &
                                   + R2*DOT_PRODUCT(RSS(:),DE3(J7,:))*DRIJDPJ(3)/ABSRIJ &
                                   - DCADPJ(3)*R2*DRIJDPI(3) + COSTA*R2*R2*DRIJDPI(3)*DRIJDPJ(3) - COSTA*D2RDIZJZ/ABSRIJ

                     D2CADIXJY     = R2*DOT_PRODUCT(DR1(J7,:),EI(:))*DRIJDPJ(2)/ABSRIJ + DOT_PRODUCT(DR2(J8,:),DE1(J7,:))/ABSRIJ  &
                                   + R2*DOT_PRODUCT(RSS(:),DE1(J7,:))*DRIJDPJ(2)/ABSRIJ &
                                   - DCADPJ(2)*R2*DRIJDPI(1) + COSTA*R2*R2*DRIJDPI(1)*DRIJDPJ(2) - COSTA*D2RDIXJY/ABSRIJ
                     D2CADIYJZ     = R2*DOT_PRODUCT(DR2(J7,:),EI(:))*DRIJDPJ(3)/ABSRIJ + DOT_PRODUCT(DR3(J8,:),DE2(J7,:))/ABSRIJ  &
                                   + R2*DOT_PRODUCT(RSS(:),DE2(J7,:))*DRIJDPJ(3)/ABSRIJ &
                                   - DCADPJ(3)*R2*DRIJDPI(2) + COSTA*R2*R2*DRIJDPI(2)*DRIJDPJ(3) - COSTA*D2RDIYJZ/ABSRIJ
                     D2CADIXJZ     = R2*DOT_PRODUCT(DR3(J7,:),EI(:))*DRIJDPJ(1)/ABSRIJ + DOT_PRODUCT(DR1(J8,:),DE3(J7,:))/ABSRIJ  &
                                   + R2*DOT_PRODUCT(RSS(:),DE3(J7,:))*DRIJDPJ(1)/ABSRIJ &
                                   - DCADPJ(1)*R2*DRIJDPI(3) + COSTA*R2*R2*DRIJDPI(3)*DRIJDPJ(1) - COSTA*D2RDIXJZ/ABSRIJ

                     D2CBDIXJX     = DOT_PRODUCT(DR1(J7,:),DE1(J8,:))/ABSRIJ - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDPJ(1)/ABSRIJ  &
                                   - DCBDPJ(1)*R2*DRIJDPI(1) + COSTB*R2*R2*DRIJDPI(1)*DRIJDPJ(1) - COSTB*D2RDIXJX/ABSRIJ
                     D2CBDIYJY     = DOT_PRODUCT(DR2(J7,:),DE2(J8,:))/ABSRIJ - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDPJ(2)/ABSRIJ  &
                                   - DCBDPJ(2)*R2*DRIJDPI(2) + COSTB*R2*R2*DRIJDPI(2)*DRIJDPJ(2) - COSTB*D2RDIYJY/ABSRIJ
                     D2CBDIZJZ     = DOT_PRODUCT(DR3(J7,:),DE3(J8,:))/ABSRIJ - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDPJ(3)/ABSRIJ  &
                                   - DCBDPJ(3)*R2*DRIJDPI(3) + COSTB*R2*R2*DRIJDPI(3)*DRIJDPJ(3) - COSTB*D2RDIZJZ/ABSRIJ

                     D2CBDIXJY     = DOT_PRODUCT(DR1(J7,:),DE2(J8,:))/ABSRIJ - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDPJ(2)/ABSRIJ  &
                                   - DCBDPJ(2)*R2*DRIJDPI(1) + COSTB*R2*R2*DRIJDPI(1)*DRIJDPJ(2) - COSTB*D2RDIXJY/ABSRIJ
                     D2CBDIYJZ     = DOT_PRODUCT(DR2(J7,:),DE3(J8,:))/ABSRIJ - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDPJ(3)/ABSRIJ  &
                                   - DCBDPJ(3)*R2*DRIJDPI(2) + COSTB*R2*R2*DRIJDPI(2)*DRIJDPJ(3) - COSTB*D2RDIYJZ/ABSRIJ
                     D2CBDIXJZ     = DOT_PRODUCT(DR3(J7,:),DE1(J8,:))/ABSRIJ - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDPJ(1)/ABSRIJ  &
                                   - DCBDPJ(1)*R2*DRIJDPI(3) + COSTB*R2*R2*DRIJDPI(3)*DRIJDPJ(1) - COSTB*D2RDIXJZ/ABSRIJ

                     D2CBDIXIX     = DOT_PRODUCT(D2R1(J7,:),EJ(:))/ABSRIJ - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDPI(1)/ABSRIJ     &
                                   - DCBDPI(1)*R2*DRIJDPI(1) + COSTB*R2*R2*DRIJDPI(1)*DRIJDPI(1) - COSTB*D2RDIXIX/ABSRIJ
                     D2CBDIXIY     = DOT_PRODUCT(D2R12(J7,:),EJ(:))/ABSRIJ - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDPI(2)/ABSRIJ     &
                                   - DCBDPI(2)*R2*DRIJDPI(1) + COSTB*R2*R2*DRIJDPI(2)*DRIJDPI(1) - COSTB*D2RDIXIY/ABSRIJ
                     D2CBDIXIZ     = DOT_PRODUCT(D2R31(J7,:),EJ(:))/ABSRIJ - R2*DOT_PRODUCT(DR1(J7,:),EJ(:))*DRIJDPI(3)/ABSRIJ     &
                                   - DCBDPI(3)*R2*DRIJDPI(1) + COSTB*R2*R2*DRIJDPI(3)*DRIJDPI(1) - COSTB*D2RDIXIZ/ABSRIJ
                     D2CBDIYIY     = DOT_PRODUCT(D2R2(J7,:),EJ(:))/ABSRIJ - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDPI(2)/ABSRIJ     &
                                   - DCBDPI(2)*R2*DRIJDPI(2) + COSTB*R2*R2*DRIJDPI(2)*DRIJDPI(2) - COSTB*D2RDIYIY/ABSRIJ
                     D2CBDIYIZ     = DOT_PRODUCT(D2R23(J7,:),EJ(:))/ABSRIJ - R2*DOT_PRODUCT(DR2(J7,:),EJ(:))*DRIJDPI(3)/ABSRIJ     &
                                   - DCBDPI(3)*R2*DRIJDPI(2) + COSTB*R2*R2*DRIJDPI(3)*DRIJDPI(2) - COSTB*D2RDIYIZ/ABSRIJ
                     D2CBDIZIZ     = DOT_PRODUCT(D2R3(J7,:),EJ(:))/ABSRIJ - R2*DOT_PRODUCT(DR3(J7,:),EJ(:))*DRIJDPI(3)/ABSRIJ     &
                                   - DCBDPI(3)*R2*DRIJDPI(3) + COSTB*R2*R2*DRIJDPI(3)*DRIJDPI(3) - COSTB*D2RDIZIZ/ABSRIJ

                     D2RHOX(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARX(:) + 3.D0*RHO20*DCADR(1)*DCADR(:)          &
                               + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRX(:) + 3.D0*RHO02*DCBDR(1)*DCBDR(:)
                     D2RHOY(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARY(:) + 3.D0*RHO20*DCADR(2)*DCADR(:)          &
                               + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRY(:) + 3.D0*RHO02*DCBDR(2)*DCBDR(:)
                     D2RHOZ(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARZ(:) + 3.D0*RHO20*DCADR(3)*DCADR(:)          &
                               + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRZ(:) + 3.D0*RHO02*DCBDR(3)*DCBDR(:)

                     D2RHODIXIX = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIXIX + 3.D0*RHO20*DCADPI(1)*DCADPI(1)       &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIXIX + 3.D0*RHO02*DCBDPI(1)*DCBDPI(1)
                     D2RHODIXIY = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIXIY + 3.D0*RHO20*DCADPI(2)*DCADPI(1)       &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIXIY + 3.D0*RHO02*DCBDPI(2)*DCBDPI(1)
                     D2RHODIXIZ = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIXIZ + 3.D0*RHO20*DCADPI(3)*DCADPI(1)       &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIXIZ + 3.D0*RHO02*DCBDPI(3)*DCBDPI(1)
                     D2RHODIYIY = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIYIY + 3.D0*RHO20*DCADPI(2)*DCADPI(2)       &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIYIY + 3.D0*RHO02*DCBDPI(2)*DCBDPI(2)
                     D2RHODIYIZ = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIYIZ + 3.D0*RHO20*DCADPI(3)*DCADPI(2)       &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIYIZ + 3.D0*RHO02*DCBDPI(3)*DCBDPI(2)
                     D2RHODIZIZ = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIZIZ + 3.D0*RHO20*DCADPI(3)*DCADPI(3)       &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIZIZ + 3.D0*RHO02*DCBDPI(3)*DCBDPI(3)

                     D2RHODIXJX = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIXJX + 3.D0*RHO20*DCADPJ(1)*DCADPI(1)   &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIXJX + 3.D0*RHO02*DCBDPJ(1)*DCBDPI(1)
                     D2RHODIYJY = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIYJY + 3.D0*RHO20*DCADPJ(2)*DCADPI(2)   &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIYJY + 3.D0*RHO02*DCBDPJ(2)*DCBDPI(2)
                     D2RHODIZJZ = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIZJZ + 3.D0*RHO20*DCADPJ(3)*DCADPI(3)   &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIZJZ + 3.D0*RHO02*DCBDPJ(3)*DCBDPI(3)

                     D2RHODIXJY = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIXJY + 3.D0*RHO20*DCADPJ(2)*DCADPI(1)   &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIXJY + 3.D0*RHO02*DCBDPJ(2)*DCBDPI(1)
                     D2RHODIYJZ = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIYJZ + 3.D0*RHO20*DCADPJ(3)*DCADPI(2)   &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIYJZ + 3.D0*RHO02*DCBDPJ(3)*DCBDPI(2)
                     D2RHODIXJZ = (RHO10 + 3.D0*RHO20*COSTA)*D2CADIXJZ + 3.D0*RHO20*DCADPJ(1)*DCADPI(3)   &
                                + (RHO01 + 3.D0*RHO02*COSTB)*D2CBDIXJZ + 3.D0*RHO02*DCBDPJ(1)*DCBDPI(3)

                     D2RHORDIX(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARDIX(:) + 3.D0*RHO20*DCADPI(1)*DCADR(:)    &
                                  + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRDIX(:) + 3.D0*RHO02*DCBDPI(1)*DCBDR(:)
                     D2RHORDIY(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARDIY(:) + 3.D0*RHO20*DCADPI(2)*DCADR(:)    &
                                  + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRDIY(:) + 3.D0*RHO02*DCBDPI(2)*DCBDR(:)
                     D2RHORDIZ(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARDIZ(:) + 3.D0*RHO20*DCADPI(3)*DCADR(:)    &
                                  + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRDIZ(:) + 3.D0*RHO02*DCBDPI(3)*DCBDR(:)

                     D2RHORDJX(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARDJX(:) + 3.D0*RHO20*DCADPJ(1)*DCADR(:)    &
                                  + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRDJX(:) + 3.D0*RHO02*DCBDPJ(1)*DCBDR(:)
                     D2RHORDJY(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARDJY(:) + 3.D0*RHO20*DCADPJ(2)*DCADR(:)    &
                                  + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRDJY(:) + 3.D0*RHO02*DCBDPJ(2)*DCBDR(:)
                     D2RHORDJZ(:) = (RHO10 + 3.D0*RHO20*COSTA)*D2CARDJZ(:) + 3.D0*RHO20*DCADPJ(3)*DCADR(:)    &
                                  + (RHO01 + 3.D0*RHO02*COSTB)*D2CBRDJZ(:) + 3.D0*RHO02*DCBDPJ(3)*DCBDR(:)

!                     DF         = D2RHODIYJZ

                  ENDIF

                  IF (SECT) THEN

!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES

!     xi,xi
                     HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2*RSS(1)*RSS(1) + DVDR - ALPHA*FRIJ(1)*NR(1)    &
                                     + ALPHA*EXPFCT/ABSRIJ*(NR(1)*NR(1) - 1.D0) + ALPHA*FRIJ(1)*DRHODR(1) + ALPHA*EXPFCT*D2RHOX(1)
!     yi,yi
                     HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2*RSS(2)*RSS(2) + DVDR - ALPHA*FRIJ(2)*NR(2)    &
                                     + ALPHA*EXPFCT/ABSRIJ*(NR(2)*NR(2) - 1.D0) + ALPHA*FRIJ(2)*DRHODR(2) + ALPHA*EXPFCT*D2RHOY(2)
!     zi,zi
                     HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2*RSS(3)*RSS(3) + DVDR - ALPHA*FRIJ(3)*NR(3)    &
                                     + ALPHA*EXPFCT/ABSRIJ*(NR(3)*NR(3) - 1.D0) + ALPHA*FRIJ(3)*DRHODR(3) + ALPHA*EXPFCT*D2RHOZ(3)
!     pi1,pi1
                     HESS(J5-2,J5-2) = HESS(J5-2,J5-2) +  D2VDR2*DRIJDPI(1)*DRIJDPI(1) + DVDR*DOT_PRODUCT(DR1(J7,:),DR1(J7,:))    &
                                     + DVDR*DOT_PRODUCT(RSS(:),D2R1(J7,:)) - ALPHA*(TIJ(1)*DRIJDPI(1)/ABSRIJ + EXPFCT*D2RDIXIX    &
                                     - TIJ(1)*DRHODPI(1) - EXPFCT*D2RHODIXIX)
!     pi2,pi2
                     HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + D2VDR2*DRIJDPI(2)*DRIJDPI(2) + DVDR*DOT_PRODUCT(DR2(J7,:),DR2(J7,:))     &
                                     + DVDR*DOT_PRODUCT(RSS(:),D2R2(J7,:)) - ALPHA*(TIJ(2)*DRIJDPI(2)/ABSRIJ + EXPFCT*D2RDIYIY    &
                                     - TIJ(2)*DRHODPI(2) - EXPFCT*D2RHODIYIY)
!     pi3,pi3
                     HESS(J5,J5)     = HESS(J5,J5)     + D2VDR2*DRIJDPI(3)*DRIJDPI(3) + DVDR*DOT_PRODUCT(DR3(J7,:),DR3(J7,:))     &
                                     + DVDR*DOT_PRODUCT(RSS(:),D2R3(J7,:)) - ALPHA*(TIJ(3)*DRIJDPI(3)/ABSRIJ + EXPFCT*D2RDIZIZ    &
                                     - TIJ(3)*DRHODPI(3) - EXPFCT*D2RHODIZIZ) 

!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME MOLECULE, DIFFERENT COORDINATES

!     xi,yi
                     DUMMY = D2VDR2*RSS(1)*RSS(2) - ALPHA*FRIJ(2)*NR(1) + ALPHA*EXPFCT/ABSRIJ*NR(2)*NR(1)     &
                           + ALPHA*FRIJ(2)*DRHODR(1) + ALPHA*EXPFCT*D2RHOY(1)
                     HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
                     HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
                     DUMMY = D2VDR2*RSS(2)*RSS(3) - ALPHA*FRIJ(3)*NR(2) + ALPHA*EXPFCT/ABSRIJ*NR(3)*NR(2)     &
                           + ALPHA*FRIJ(3)*DRHODR(2) + ALPHA*EXPFCT*D2RHOZ(2)
                     HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
                     HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
                     DUMMY = D2VDR2*RSS(3)*RSS(1) - ALPHA*FRIJ(1)*NR(3) + ALPHA*EXPFCT/ABSRIJ*NR(1)*NR(3)     &
                           + ALPHA*FRIJ(1)*DRHODR(3) + ALPHA*EXPFCT*D2RHOX(3)
                     HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
                     HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
!     xi,pi1
                     DUMMY = D2VDR2*DRIJDPI(1)*RSS(1) + DVDR*DR1(J7,1) - ALPHA*TIJ(1)*NR(1) - ALPHA*EXPFCT*DR1(J7,1)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(1)*RSS(1)/ABSRIJ + ALPHA*TIJ(1)*DRHODR(1) + ALPHA*EXPFCT*D2RHORDIX(1)
                     HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
                     HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1
                     DUMMY = D2VDR2*DRIJDPI(1)*RSS(2) + DVDR*DR1(J7,2) - ALPHA*TIJ(1)*NR(2) - ALPHA*EXPFCT*DR1(J7,2)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(1)*RSS(2)/ABSRIJ + ALPHA*TIJ(1)*DRHODR(2) + ALPHA*EXPFCT*D2RHORDIX(2)
                     HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
                     HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
                     DUMMY = D2VDR2*DRIJDPI(1)*RSS(3) + DVDR*DR1(J7,3) - ALPHA*TIJ(1)*NR(3) - ALPHA*EXPFCT*DR1(J7,3)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(1)*RSS(3)/ABSRIJ + ALPHA*TIJ(1)*DRHODR(3) + ALPHA*EXPFCT*D2RHORDIX(3)
                     HESS(J3,J5-2)   = HESS(J3,J5-2) + DUMMY
                     HESS(J5-2,J3)   = HESS(J5-2,J3) + DUMMY
!     xi,pi2
                     DUMMY = D2VDR2*DRIJDPI(2)*RSS(1) + DVDR*DR2(J7,1) - ALPHA*TIJ(2)*NR(1) - ALPHA*EXPFCT*DR2(J7,1)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(2)*RSS(1)/ABSRIJ + ALPHA*TIJ(2)*DRHODR(1) + ALPHA*EXPFCT*D2RHORDIY(1) 
                     HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
                     HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2
                     DUMMY = D2VDR2*DRIJDPI(2)*RSS(2) + DVDR*DR2(J7,2) - ALPHA*TIJ(2)*NR(2) - ALPHA*EXPFCT*DR2(J7,2)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(2)*RSS(2)/ABSRIJ + ALPHA*TIJ(2)*DRHODR(2) + ALPHA*EXPFCT*D2RHORDIY(2) 
                     HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
                     HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
                     DUMMY = D2VDR2*DRIJDPI(2)*RSS(3) + DVDR*DR2(J7,3) - ALPHA*TIJ(2)*NR(3) - ALPHA*EXPFCT*DR2(J7,3)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(2)*RSS(3)/ABSRIJ + ALPHA*TIJ(2)*DRHODR(3) + ALPHA*EXPFCT*D2RHORDIY(3) 
                     HESS(J3,J5-1)   = HESS(J3,J5-1) + DUMMY
                     HESS(J5-1,J3)   = HESS(J5-1,J3) + DUMMY
!     xi,pi3
                     DUMMY = D2VDR2*DRIJDPI(3)*RSS(1) + DVDR*DR3(J7,1) - ALPHA*TIJ(3)*NR(1) - ALPHA*EXPFCT*DR3(J7,1)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(3)*RSS(1)/ABSRIJ + ALPHA*TIJ(3)*DRHODR(1) + ALPHA*EXPFCT*D2RHORDIZ(1)
                     HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
                     HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi3
                     DUMMY = D2VDR2*DRIJDPI(3)*RSS(2) + DVDR*DR3(J7,2) - ALPHA*TIJ(3)*NR(2) - ALPHA*EXPFCT*DR3(J7,2)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(3)*RSS(2)/ABSRIJ + ALPHA*TIJ(3)*DRHODR(2) + ALPHA*EXPFCT*D2RHORDIZ(2)
                     HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
                     HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi3
                     DUMMY = D2VDR2*DRIJDPI(3)*RSS(3) + DVDR*DR3(J7,3) - ALPHA*TIJ(3)*NR(3) - ALPHA*EXPFCT*DR3(J7,3)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPI(3)*RSS(3)/ABSRIJ + ALPHA*TIJ(3)*DRHODR(3) + ALPHA*EXPFCT*D2RHORDIZ(3)
                     HESS(J3,J5)     = HESS(J3,J5) + DUMMY
                     HESS(J5,J3)     = HESS(J5,J3) + DUMMY
!     pi1,pi2
                     DUMMY = D2VDR2*DRIJDPI(1)*DRIJDPI(2) + DVDR*DOT_PRODUCT(DR2(J7,:),DR1(J7,:))   &
                           + DVDR*DOT_PRODUCT(RSS,D2R12(J7,:)) - ALPHA*TIJ(2)*DRIJDPI(1)/ABSRIJ - ALPHA*EXPFCT*D2RDIXIY &
                           + ALPHA*TIJ(2)*DRHODPI(1) + ALPHA*EXPFCT*D2RHODIXIY
                     HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
                     HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
                     DUMMY = D2VDR2*DRIJDPI(2)*DRIJDPI(3) + DVDR*DOT_PRODUCT(DR3(J7,:),DR2(J7,:))   &
                           + DVDR*DOT_PRODUCT(RSS,D2R23(J7,:)) -ALPHA*TIJ(3)*DRIJDPI(2)/ABSRIJ - ALPHA*EXPFCT*D2RDIYIZ  &
                           + ALPHA*TIJ(3)*DRHODPI(2) + ALPHA*EXPFCT*D2RHODIYIZ
                     HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
                     HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
                     DUMMY = D2VDR2*DRIJDPI(3)*DRIJDPI(1) + DVDR*DOT_PRODUCT(DR1(J7,:),DR3(J7,:))   &
                           + DVDR*DOT_PRODUCT(RSS,D2R31(J7,:)) - ALPHA*TIJ(1)*DRIJDPI(3)/ABSRIJ - ALPHA*EXPFCT*D2RDIXIZ &
                           + ALPHA*TIJ(1)*DRHODPI(3) + ALPHA*EXPFCT*D2RHODIXIZ
                     HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
                     HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE

!     xi,xj
                     HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2*RSS(1)*RSS(1) - DVDR + ALPHA*FRIJ(1)*NR(1)    &
                                     - ALPHA*EXPFCT/ABSRIJ*(NR(1)*NR(1) - 1.D0) - ALPHA*FRIJ(1)*DRHODR(1) - ALPHA*EXPFCT*D2RHOX(1)
!     yi,yj
                     HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2*RSS(2)*RSS(2) - DVDR + ALPHA*FRIJ(2)*NR(2)    &
                                     - ALPHA*EXPFCT/ABSRIJ*(NR(2)*NR(2) - 1.D0) - ALPHA*FRIJ(2)*DRHODR(2) - ALPHA*EXPFCT*D2RHOY(2)
!     zi,zj
                     HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2*RSS(3)*RSS(3) - DVDR + ALPHA*FRIJ(3)*NR(3)    &
                                     - ALPHA*EXPFCT/ABSRIJ*(NR(3)*NR(3) - 1.D0) - ALPHA*FRIJ(3)*DRHODR(3) - ALPHA*EXPFCT*D2RHOZ(3)
!     pi1,pj1
                     HESS(J5-2,J6-2) = HESS(J5-2,J6-2) + D2VDR2*DRIJDPJ(1)*DRIJDPI(1) - DVDR*DOT_PRODUCT(DR1(J8,:),DR1(J7,:))     &
                                     - ALPHA*(TJI(1)*DRIJDPI(1)/ABSRIJ + EXPFCT*D2RDIXJX - TJI(1)*DRHODPI(1) - EXPFCT*D2RHODIXJX)
!     pi2,pj2
                     HESS(J5-1,J6-1) = HESS(J5-1,J6-1) + D2VDR2*DRIJDPJ(2)*DRIJDPI(2) - DVDR*DOT_PRODUCT(DR2(J8,:),DR2(J7,:))     &
                                     - ALPHA*(TJI(2)*DRIJDPI(2)/ABSRIJ + EXPFCT*D2RDIYJY - TJI(2)*DRHODPI(2) - EXPFCT*D2RHODIYJY)
!     pi3,pj3
                     HESS(J5,J6)     = HESS(J5,J6)     + D2VDR2*DRIJDPJ(3)*DRIJDPI(3) - DVDR*DOT_PRODUCT(DR3(J8,:),DR3(J7,:))     &
                                     - ALPHA*(TJI(3)*DRIJDPI(3)/ABSRIJ + EXPFCT*D2RDIZJZ - TJI(3)*DRHODPI(3) - EXPFCT*D2RHODIZJZ)

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES

!     xi,yj
                     DUMMY =-D2VDR2*RSS(1)*RSS(2) + ALPHA*FRIJ(2)*NR(1) - ALPHA*EXPFCT/ABSRIJ*NR(2)*NR(1)     &
                           - ALPHA*FRIJ(2)*DRHODR(1) - ALPHA*EXPFCT*D2RHOY(1)
                     HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
                     HESS(J4-1,J3-2) = HESS(J4-1,J3-2) + DUMMY
!     yi,zj
                     DUMMY =-D2VDR2*RSS(2)*RSS(3) + ALPHA*FRIJ(3)*NR(2) - ALPHA*EXPFCT/ABSRIJ*NR(3)*NR(2)     &
                            - ALPHA*FRIJ(3)*DRHODR(2) - ALPHA*EXPFCT*D2RHOZ(2)
                     HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
                     HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
                     DUMMY =-D2VDR2*RSS(3)*RSS(1) + ALPHA*FRIJ(1)*NR(3) - ALPHA*EXPFCT/ABSRIJ*NR(1)*NR(3)     &
                           - ALPHA*FRIJ(1)*DRHODR(3) - ALPHA*EXPFCT*D2RHOX(3)
                     HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
                     HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY
!     xi,pj1
                     DUMMY = D2VDR2*DRIJDPJ(1)*RSS(1) - DVDR*DR1(J8,1) - ALPHA*TJI(1)*NR(1) + ALPHA*EXPFCT*DR1(J8,1)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(1)*RSS(1)/ABSRIJ + ALPHA*TJI(1)*DRHODR(1) + ALPHA*EXPFCT*D2RHORDJX(1) 
                     HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
                     HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
                     DUMMY = D2VDR2*DRIJDPJ(1)*RSS(2) - DVDR*DR1(J8,2) - ALPHA*TJI(1)*NR(2) + ALPHA*EXPFCT*DR1(J8,2)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(1)*RSS(2)/ABSRIJ + ALPHA*TJI(1)*DRHODR(2) + ALPHA*EXPFCT*D2RHORDJX(2) 
                     HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
                     HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
                     DUMMY = D2VDR2*DRIJDPJ(1)*RSS(3) - DVDR*DR1(J8,3) - ALPHA*TJI(1)*NR(3) + ALPHA*EXPFCT*DR1(J8,3)/ABSRIJ      &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(1)*RSS(3)/ABSRIJ + ALPHA*TJI(1)*DRHODR(3) + ALPHA*EXPFCT*D2RHORDJX(3) 
                     HESS(J3,J6-2)   = HESS(J3,J6-2) + DUMMY
                     HESS(J6-2,J3)   = HESS(J6-2,J3) + DUMMY
!     xi,pj2
                     DUMMY = D2VDR2*DRIJDPJ(2)*RSS(1) - DVDR*DR2(J8,1) - ALPHA*TJI(2)*NR(1) + ALPHA*EXPFCT*DR2(J8,1)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(2)*RSS(1)/ABSRIJ + ALPHA*TJI(2)*DRHODR(1) + ALPHA*EXPFCT*D2RHORDJY(1) 
                     HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
                     HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
                     DUMMY = D2VDR2*DRIJDPJ(2)*RSS(2) - DVDR*DR2(J8,2) - ALPHA*TJI(2)*NR(2) + ALPHA*EXPFCT*DR2(J8,2)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(2)*RSS(2)/ABSRIJ + ALPHA*TJI(2)*DRHODR(2) + ALPHA*EXPFCT*D2RHORDJY(2)
                     HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
                     HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
                     DUMMY = D2VDR2*DRIJDPJ(2)*RSS(3) - DVDR*DR2(J8,3) - ALPHA*TJI(2)*NR(3) + ALPHA*EXPFCT*DR2(J8,3)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(2)*RSS(3)/ABSRIJ + ALPHA*TJI(2)*DRHODR(3) + ALPHA*EXPFCT*D2RHORDJY(3)
                     HESS(J3,J6-1)   = HESS(J3,J6-1) + DUMMY
                     HESS(J6-1,J3)   = HESS(J6-1,J3) + DUMMY
!     xi,pj3
                     DUMMY = D2VDR2*DRIJDPJ(3)*RSS(1) - DVDR*DR3(J8,1) - ALPHA*TJI(3)*NR(1) + ALPHA*EXPFCT*DR3(J8,1)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(3)*RSS(1)/ABSRIJ + ALPHA*TJI(3)*DRHODR(1) + ALPHA*EXPFCT*D2RHORDJZ(1) 
                     HESS(J3-2,J6)   = HESS(J3-2,J6) + DUMMY
                     HESS(J6,J3-2)   = HESS(J6,J3-2) + DUMMY
!     yi,pj3
                     DUMMY = D2VDR2*DRIJDPJ(3)*RSS(2) - DVDR*DR3(J8,2) - ALPHA*TJI(3)*NR(2) + ALPHA*EXPFCT*DR3(J8,2)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(3)*RSS(2)/ABSRIJ + ALPHA*TJI(3)*DRHODR(2) + ALPHA*EXPFCT*D2RHORDJZ(2) 
                     HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
                     HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
                     DUMMY = D2VDR2*DRIJDPJ(3)*RSS(3) - DVDR*DR3(J8,3) - ALPHA*TJI(3)*NR(3) + ALPHA*EXPFCT*DR3(J8,3)/ABSRIJ       &
                           + ALPHA*EXPFCT*R2*DRIJDPJ(3)*RSS(3)/ABSRIJ + ALPHA*TJI(3)*DRHODR(3) + ALPHA*EXPFCT*D2RHORDJZ(3) 
                     HESS(J3,J6)     = HESS(J3,J6) + DUMMY
                     HESS(J6,J3)     = HESS(J6,J3) + DUMMY
!     pi1,pj2
                     DUMMY = D2VDR2*DRIJDPI(1)*DRIJDPJ(2) - DVDR*DOT_PRODUCT(DR2(J8,:),DR1(J7,:)) - ALPHA*TJI(2)*DRIJDPI(1)/ABSRIJ &
                           - ALPHA*EXPFCT*D2RDIXJY + ALPHA*TJI(2)*DRHODPI(1) + ALPHA*EXPFCT*D2RHODIXJY
                     HESS(J5-2,J6-1) = HESS(J5-2,J6-1) + DUMMY
                     HESS(J6-1,J5-2) = HESS(J6-1,J5-2) + DUMMY
!     pi2,pj3
                     DUMMY = D2VDR2*DRIJDPI(2)*DRIJDPJ(3) - DVDR*DOT_PRODUCT(DR3(J8,:),DR2(J7,:)) - ALPHA*TJI(3)*DRIJDPI(2)/ABSRIJ &
                           - ALPHA*EXPFCT*D2RDIYJZ + ALPHA*TJI(3)*DRHODPI(2) + ALPHA*EXPFCT*D2RHODIYJZ
                     HESS(J5-1,J6)   = HESS(J5-1,J6) + DUMMY
                     HESS(J6,J5-1)   = HESS(J6,J5-1) + DUMMY
!     pi3,pj1
                     DUMMY = D2VDR2*DRIJDPI(3)*DRIJDPJ(1) - DVDR*DOT_PRODUCT(DR1(J8,:),DR3(J7,:)) - ALPHA*TJI(1)*DRIJDPI(3)/ABSRIJ &
                           - ALPHA*EXPFCT*D2RDIXJZ + ALPHA*TJI(1)*DRHODPI(3) + ALPHA*EXPFCT*D2RHODIXJZ 
                     HESS(J5,J6-2)   = HESS(J5,J6-2) + DUMMY
                     HESS(J6-2,J5)   = HESS(J6-2,J5) + DUMMY

                  ENDIF

               ENDDO

            ENDDO
 
         ENDDO

      ENDDO

      ENERGY = FCTR*2625.499D0*(ENERGY1 + ENERGY2 + ENERGY3) 
      IF (GTEST) G(:) = FCTR*2625.499D0*G(:)
      IF (SECT) HESS(:,:) = 2625.499D0*HESS(:,:)

      END SUBROUTINE PAHAGH

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPAHA()

      USE KEY, ONLY: RHOCC0, RHOCC10, RHOCC20,  RHOHH0, RHOHH10, RHOHH20, RHOCH0, RHOC10H, RHOCH10, RHOC20H, RHOCH20,   &
                     ALPHACC, ALPHAHH, ALPHACH, DC6CC, DC6HH, DC6CH, KKJ, CCKJ

      IMPLICIT NONE

      ALPHACC = 1.861500D0
      ALPHAHH = 1.431200D0
      ALPHACH = 1.775600D0

      DC6CC    = 30.469D0
      DC6HH    = 5.359D0
      DC6CH    = 12.840D0

      RHOCC0  = 5.814700D0
      RHOCC10 = 0.021700D0
      RHOCC20 =-0.220800D0

      RHOHH0  = 4.486200D0
      RHOHH10 =-0.271800D0
      RHOHH20 = 0.0D0

      RHOCH0  = 5.150500D0
      RHOC10H = 0.021700D0
      RHOCH10 =-0.271800D0
      RHOC20H =-0.220800D0
      RHOCH20 = 0.0D0

      KKJ     = 1.D-03
      CCKJ    = 1.D0   !1389.354848D0

      END SUBROUTINE DEFPAHA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBENZENE()

      USE COMMONS, ONLY: NRBSITES, RBSITE, RBSTLA, STCHRG

      IMPLICIT NONE
 
      INTEGER :: J1

!     C6H6

!      RBSITE(1,:)  = (/-2.639231D0,    0.D0,         0.D0/)   !C1
!      RBSITE(2,:)  = (/ 2.63923100D0,  0.D0,         0.D0/)   !C2
!      RBSITE(3,:)  = (/-1.31961456D0, -2.28564437D0, 0.D0/)   !C3
!      RBSITE(4,:)  = (/ 1.31961456D0, -2.28564437D0, 0.D0/)   !C4
!      RBSITE(5,:)  = (/-1.31961456D0,  2.28564437D0, 0.D0/)   !C5
!      RBSITE(6,:)  = (/ 1.31961456D0,  2.28564437D0, 0.D0/)   !C6
!      RBSITE(7,:)  = (/-4.69338583D0,  0.D0,         0.D0/)   !H1
!      RBSITE(8,:)  = (/ 4.69338583D0,  0.D0,         0.D0/)   !H2
!      RBSITE(9,:)  = (/ 2.34669008D0,  4.06459651D0, 0.D0/)   !H3
!      RBSITE(10,:) = (/-2.34669008D0,  4.06459651D0, 0.D0/)   !H4
!      RBSITE(11,:) = (/ 2.34669008D0, -4.06459651D0, 0.D0/)   !H5
!      RBSITE(12,:) = (/-2.34669008D0, -4.06459651D0, 0.D0/)   !H6

      RBSITE(1,:)  = (/ 2.63923430843701,   0.00000000000000,   0.00000000000000/)
      RBSITE(2,:)  = (/ 1.31961715421850,  -2.28564395764590,   0.00000000000000/)
      RBSITE(3,:)  = (/-1.31961715421850,  -2.28564395764590,   0.00000000000000/)
      RBSITE(4,:)  = (/-2.63923430843701,   0.00000000000000,   0.00000000000000/)
      RBSITE(5,:)  = (/-1.31961715421850,   2.28564395764590,   0.00000000000000/)
      RBSITE(6,:)  = (/ 1.31961715421850,   2.28564395764590,   0.00000000000000/)
      RBSITE(7,:)  = (/ 4.69338981379532,   0.00000000000000,   0.00000000000000/)
      RBSITE(8,:)  = (/ 2.34669490689766,  -4.06459480860986,   0.00000000000000/)
      RBSITE(9,:)  = (/-2.34669490689766,  -4.06459480860986,   0.00000000000000/)
      RBSITE(10,:) = (/-4.69338981379532,   0.00000000000000,   0.00000000000000/)
      RBSITE(11,:) = (/-2.34669490689766,   4.06459480860986,   0.00000000000000/)
      RBSITE(12,:) = (/ 2.34669490689766,   4.06459480860986,   0.00000000000000/)

!      RBSITE(:,:) = RBSITE(:,:)*0.529177D0              ! au to angstrom

!      RBSTLA(1,:)  = RBSITE(7,:)  - RBSITE(1,:)                 ! Z FROM C1 TO H1
!      RBSTLA(2,:)  = RBSITE(8,:)  - RBSITE(2,:)                 ! Z FROM C2 TO H2
!      RBSTLA(3,:)  = RBSITE(12,:) - RBSITE(3,:)                 ! Z FROM C3 TO H6
!      RBSTLA(4,:)  = RBSITE(11,:) - RBSITE(4,:)                 ! Z FROM C4 TO H5
!      RBSTLA(5,:)  = RBSITE(10,:) - RBSITE(5,:)                 ! Z FROM C5 TO H4
!      RBSTLA(6,:)  = RBSITE(9,:)  - RBSITE(6,:)                 ! Z FROM C6 TO H3
!      RBSTLA(7,:)  = RBSITE(7,:)  - RBSITE(1,:)                 ! Z FROM C1 TO H1
!      RBSTLA(8,:)  = RBSITE(8,:)  - RBSITE(2,:)                 ! Z FROM C2 TO H2
!      RBSTLA(9,:)  = RBSITE(9,:)  - RBSITE(6,:)                 ! Z FROM C6 TO H3
!      RBSTLA(10,:) = RBSITE(10,:) - RBSITE(5,:)                 ! Z FROM C5 TO H4
!      RBSTLA(11,:) = RBSITE(11,:) - RBSITE(4,:)                 ! Z FROM C4 TO H5
!      RBSTLA(12,:) = RBSITE(12,:) - RBSITE(3,:)                 ! Z FROM C3 TO H6

      RBSTLA(1,:)  = RBSITE(7,:)  - RBSITE(1,:)                 ! Z FROM C1 TO H1
      RBSTLA(2,:)  = RBSITE(8,:)  - RBSITE(2,:)                 ! Z FROM C2 TO H2
      RBSTLA(3,:)  = RBSITE(9,:)  - RBSITE(3,:)                 ! Z FROM C3 TO H3
      RBSTLA(4,:)  = RBSITE(10,:) - RBSITE(4,:)                 ! Z FROM C4 TO H4
      RBSTLA(5,:)  = RBSITE(11,:) - RBSITE(5,:)                 ! Z FROM C5 TO H5
      RBSTLA(6,:)  = RBSITE(12,:) - RBSITE(6,:)                 ! Z FROM C6 TO H6
      RBSTLA(7,:)  = RBSITE(7,:)  - RBSITE(1,:)                 ! Z FROM C1 TO H1!
      RBSTLA(8,:)  = RBSITE(8,:)  - RBSITE(2,:)                 ! Z FROM C2 TO H2!
      RBSTLA(9,:)  = RBSITE(9,:) -  RBSITE(3,:)                 ! Z FROM C3 TO H3!
      RBSTLA(10,:) = RBSITE(10,:) - RBSITE(4,:)                 ! Z FROM C4 TO H4!
      RBSTLA(11,:) = RBSITE(11,:) - RBSITE(5,:)                 ! Z FROM C5 TO H5!
      RBSTLA(12,:) = RBSITE(12,:) - RBSITE(6,:)                 ! Z FROM C6 TO H6!
      
      DO J1 = 1, NRBSITES
 
         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      STCHRG(1:6)  =-0.11114D0
      STCHRG(7:12) = 0.11114D0

      END SUBROUTINE DEFBENZENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFNAPHTHALENE()

      USE COMMONS, ONLY: NRBSITES, RBSITE, RBSTLA, STCHRG

      IMPLICIT NONE
 
      INTEGER :: J1

!     C10H8

      RBSITE(1,:)  = (/-1.33862D0, -4.59918D0, 0.D0/)   ! C1
      RBSITE(2,:)  = (/-2.65019D0, -2.35249D0, 0.D0/)   ! C2
      RBSITE(3,:)  = (/-1.35523D0, 0.D0, 0.D0/)         ! C3
      RBSITE(4,:)  = (/ 1.35523D0, 0.D0, 0.D0/)         ! C4
      RBSITE(5,:)  = (/ 2.65019D0,-2.35249D0, 0.D0/)    ! C5
      RBSITE(6,:)  = (/ 1.33862D0,-4.59918D0, 0.D0/)    ! C6
      RBSITE(7,:)  = (/-2.65019D0, 2.35249D0, 0.D0/)    ! C7
      RBSITE(8,:)  = (/ 2.65019D0, 2.35249D0, 0.D0/)    ! C8
      RBSITE(9,:)  = (/ 1.33862D0, 4.59918D0, 0.D0/)    ! C9
      RBSITE(10,:) = (/-1.33862D0, 4.59918D0, 0.D0/)    ! C10
      RBSITE(11,:) = (/-4.70575D0, 2.34799D0, 0.D0/)    ! H1
      RBSITE(12,:) = (/-2.35493D0,-6.38388D0, 0.D0/)    ! H2
      RBSITE(13,:) = (/-4.70575D0,-2.34799D0, 0.D0/)    ! H3
      RBSITE(14,:) = (/ 4.70575D0,-2.34799D0, 0.D0/)    ! H4
      RBSITE(15,:) = (/ 2.35493D0,-6.38388D0, 0.D0/)    ! H5
      RBSITE(16,:) = (/ 4.70575D0, 2.34799D0, 0.D0/)    ! H6
      RBSITE(17,:) = (/ 2.35493D0, 6.38388D0, 0.D0/)    ! H7
      RBSITE(18,:) = (/-2.35493D0, 6.38388D0, 0.D0/)    ! H8

!      RBSITE(:,:) = RBSITE(:,:)*0.529177D0              ! au to angstrom

      STCHRG(1)  =-0.10048D0 
      STCHRG(2)  =-0.29796D0
      STCHRG(3)  =-0.24018D0
      STCHRG(4)  = 0.24018D0
      STCHRG(5)  =-0.29796D0
      STCHRG(6)  =-0.10048D0
      STCHRG(7)  =-0.29796D0
      STCHRG(8)  =-0.29796D0
      STCHRG(9)  =-0.10048D0 
      STCHRG(10) =-0.10048D0
      STCHRG(11) = 0.15530D0
      STCHRG(12) = 0.12304D0
      STCHRG(13) = 0.15530D0
      STCHRG(14) = 0.15530D0
      STCHRG(15) = 0.12304D0
      STCHRG(16) = 0.15530D0
      STCHRG(17) = 0.12304D0
      STCHRG(18) = 0.12304D0

      RBSTLA(1,:)  = RBSITE(12,:) - RBSITE(1,:)                 ! Z FROM C1 TO H2
      RBSTLA(2,:)  = RBSITE(13,:) - RBSITE(2,:)                 ! Z FROM C2 TO H3
      RBSTLA(3,:)  = RBSITE(4,:)  - RBSITE(3,:)                 ! Z FROM C3 TO C4
      RBSTLA(4,:)  = RBSITE(3,:)  - RBSITE(4,:)                 ! Z FROM C4 TO C3
      RBSTLA(5,:)  = RBSITE(14,:) - RBSITE(5,:)                 ! Z FROM C5 TO H4
      RBSTLA(6,:)  = RBSITE(15,:) - RBSITE(6,:)                 ! Z FROM C6 TO H5
      RBSTLA(7,:)  = RBSITE(11,:) - RBSITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(8,:)  = RBSITE(16,:) - RBSITE(8,:)                 ! Z FROM C8 TO H6
      RBSTLA(9,:)  = RBSITE(17,:) - RBSITE(9,:)                 ! Z FROM C9 TO H7
      RBSTLA(10,:) = RBSITE(18,:) - RBSITE(10,:)                ! Z FROM C10 TO H8
      RBSTLA(11,:) = RBSITE(11,:) - RBSITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(12,:) = RBSITE(12,:) - RBSITE(1,:)                 ! Z FROM C1 TO H2
      RBSTLA(13,:) = RBSITE(13,:) - RBSITE(2,:)                 ! Z FROM C2 TO H3
      RBSTLA(14,:) = RBSITE(14,:) - RBSITE(5,:)                 ! Z FROM C5 TO H4
      RBSTLA(15,:) = RBSITE(15,:) - RBSITE(6,:)                 ! Z FROM C6 TO H5
      RBSTLA(16,:) = RBSITE(16,:) - RBSITE(8,:)                 ! Z FROM C8 TO H6
      RBSTLA(17,:) = RBSITE(17,:) - RBSITE(9,:)                 ! Z FROM C9 TO H7
      RBSTLA(18,:) = RBSITE(18,:) - RBSITE(10,:)                ! Z FROM C10 TO H8
      
      DO J1 = 1, NRBSITES
 
         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFNAPHTHALENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFANTHRACENE()

      USE COMMONS, ONLY: NRBSITES, RBSITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C14H10

      RBSITE(1,:)  = (/ 1.36540D0, 2.31298D0, 0.D0/)    ! C1
      RBSITE(2,:)  = (/-1.36540D0, 2.31298D0, 0.D0/)    ! C2     
      RBSITE(3,:)  = (/ 1.36540D0,-2.31298D0, 0.D0/)    ! C3
      RBSITE(4,:)  = (/-1.36540D0,-2.31298D0, 0.D0/)    ! C4
      RBSITE(5,:)  = (/-2.65253D0, 0.D0, 0.D0/)         ! C5
      RBSITE(6,:)  = (/ 2.65253D0, 0.D0, 0.D0/)         ! C6
      RBSITE(7,:)  = (/ 2.65927D0, 4.68538D0, 0.D0/)    ! C7
      RBSITE(8,:)  = (/-2.65927D0, 4.68538D0, 0.D0/)    ! C8
      RBSITE(9,:)  = (/ 2.65927D0,-4.68538D0, 0.D0/)    ! C9
      RBSITE(10,:) = (/-2.65927D0,-4.68538D0, 0.D0/)    ! C10
      RBSITE(11,:) = (/ 1.34762D0,-6.91760D0, 0.D0/)    ! C11
      RBSITE(12,:) = (/-1.34762D0,-6.91760D0, 0.D0/)    ! C12
      RBSITE(13,:) = (/ 1.34762D0, 6.91760D0, 0.D0/)    ! C13
      RBSITE(14,:) = (/-1.34762D0, 6.91760D0, 0.D0/)    ! C14
      RBSITE(15,:) = (/ 4.71450D0,-4.67888D0, 0.D0/)    ! H1
      RBSITE(16,:) = (/ 2.35428D0,-8.70751D0, 0.D0/)    ! H2
      RBSITE(17,:) = (/-2.35428D0,-8.70751D0, 0.D0/)    ! H3
      RBSITE(18,:) = (/-4.71450D0,-4.67888D0, 0.D0/)    ! H4
      RBSITE(19,:) = (/ 4.71450D0, 4.67888D0, 0.D0/)    ! H5
      RBSITE(20,:) = (/ 2.35428D0, 8.70751D0, 0.D0/)    ! H6
      RBSITE(21,:) = (/-2.35428D0, 8.70751D0, 0.D0/)    ! H7
      RBSITE(22,:) = (/-4.71450D0, 4.67888D0, 0.D0/)    ! H8
      RBSITE(23,:) = (/-4.70918D0, 0.D0, 0.D0/)         ! H9
      RBSITE(24,:) = (/ 4.70918D0, 0.D0, 0.D0/)         ! H10

!      RBSITE(:,:) = RBSITE(:,:)*0.529177D0              ! au to angstrom

      STCHRG(1)  = 0.23448D0 
      STCHRG(2)  = 0.23448D0
      STCHRG(3)  = 0.23448D0
      STCHRG(4)  = 0.23448D0
      STCHRG(5)  =-0.47174D0
      STCHRG(6)  =-0.47174D0
      STCHRG(7)  =-0.25252D0
      STCHRG(8)  =-0.25252D0
      STCHRG(9)  =-0.25252D0
      STCHRG(10) =-0.25252D0
      STCHRG(11) =-0.11389D0
      STCHRG(12) =-0.11389D0
      STCHRG(13) =-0.11389D0
      STCHRG(14) =-0.11389D0
      STCHRG(15) = 0.14291D0
      STCHRG(16) = 0.12531D0
      STCHRG(17) = 0.12531D0
      STCHRG(18) = 0.14291D0
      STCHRG(19) = 0.14291D0
      STCHRG(20) = 0.12531D0
      STCHRG(21) = 0.12531D0
      STCHRG(22) = 0.14291D0
      STCHRG(23) = 0.19915D0
      STCHRG(24) = 0.19915D0

      RBSTLA(1,:)  = RBSITE(2,:)  - RBSITE(1,:)                 ! Z FROM C1 TO C2
      RBSTLA(2,:)  = RBSITE(1,:)  - RBSITE(2,:)                 ! Z FROM C2 TO C1
      RBSTLA(3,:)  = RBSITE(4,:)  - RBSITE(3,:)                 ! Z FROM C3 TO C4
      RBSTLA(4,:)  = RBSITE(3,:)  - RBSITE(4,:)                 ! Z FROM C4 TO C3
      RBSTLA(5,:)  = RBSITE(23,:) - RBSITE(5,:)                 ! Z FROM C5 TO H9
      RBSTLA(6,:)  = RBSITE(24,:) - RBSITE(6,:)                 ! Z FROM C6 TO H10
      RBSTLA(7,:)  = RBSITE(19,:) - RBSITE(7,:)                 ! Z FROM C7 TO H5
      RBSTLA(8,:)  = RBSITE(22,:) - RBSITE(8,:)                 ! Z FROM C8 TO H8
      RBSTLA(9,:)  = RBSITE(15,:) - RBSITE(9,:)                 ! Z FROM C9 TO H1
      RBSTLA(10,:) = RBSITE(18,:) - RBSITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(11,:) = RBSITE(16,:) - RBSITE(11,:)                ! Z FROM C11 TO H2
      RBSTLA(12,:) = RBSITE(17,:) - RBSITE(12,:)                ! Z FROM C12 TO H3
      RBSTLA(13,:) = RBSITE(20,:) - RBSITE(13,:)                ! Z FROM C13 TO H6
      RBSTLA(14,:) = RBSITE(21,:) - RBSITE(14,:)                ! Z FROM C14 TO H7
      RBSTLA(15,:) = RBSITE(15,:) - RBSITE(9,:)                 ! Z FROM C9 TO H1
      RBSTLA(16,:) = RBSITE(16,:) - RBSITE(11,:)                ! Z FROM C11 TO H2
      RBSTLA(17,:) = RBSITE(17,:) - RBSITE(12,:)                ! Z FROM C12 TO H3
      RBSTLA(18,:) = RBSITE(18,:) - RBSITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(19,:) = RBSITE(19,:) - RBSITE(7,:)                 ! Z FROM C7 TO H5
      RBSTLA(20,:) = RBSITE(20,:) - RBSITE(13,:)                ! Z FROM C13 TO H6
      RBSTLA(21,:) = RBSITE(21,:) - RBSITE(14,:)                ! Z FROM C14 TO H7
      RBSTLA(22,:) = RBSITE(22,:) - RBSITE(8,:)                 ! Z FROM C8 TO H8
      RBSTLA(23,:) = RBSITE(23,:) - RBSITE(5,:)                 ! Z FROM C5 TO H9
      RBSTLA(24,:) = RBSITE(24,:) - RBSITE(6,:)                 ! Z FROM C6 TO H10

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFANTHRACENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPYRENE()

      USE COMMONS, ONLY: NRBSITES, RBSITE, RBSTLA, STCHRG

      IMPLICIT NONE

      INTEGER :: J1

!     C16H10

      RBSITE(1,:)  = (/-1.34794D0, 0.D0, 0.D0/)         ! C1
      RBSITE(2,:)  = (/ 1.34794D0, 0.D0, 0.D0/)         ! C2
      RBSITE(3,:)  = (/ 2.70059D0, 2.33625D0, 0.D0/)    ! C3
      RBSITE(4,:)  = (/ 2.70059d0,-2.33625D0, 0.D0/)    ! C4
      RBSITE(5,:)  = (/-2.70059D0,-2.33625D0, 0.D0/)    ! C5
      RBSITE(6,:)  = (/-2.70059D0, 2.33625D0, 0.D0/)    ! C6
      RBSITE(7,:)  = (/ 1.28651D0, 4.65603D0, 0.D0/)    ! C7
      RBSITE(8,:)  = (/ 5.35355D0, 2.28771D0, 0.D0/)    ! C8
      RBSITE(9,:)  = (/ 1.28651D0,-4.65603D0, 0.D0/)    ! C9
      RBSITE(10,:) = (/ 5.35355D0,-2.28771D0, 0.D0/)    ! C10
      RBSITE(11,:) = (/-1.28651D0,-4.65603D0, 0.D0/)    ! C11
      RBSITE(12,:) = (/-5.35355D0,-2.28771D0, 0.D0/)    ! C12
      RBSITE(13,:) = (/-1.28651D0, 4.65603D0, 0.D0/)    ! C13
      RBSITE(14,:) = (/-5.35355D0, 2.28771D0, 0.D0/)    ! C14
      RBSITE(15,:) = (/ 6.65929D0, 0.D0, 0.D0/)         ! C15
      RBSITE(16,:) = (/-6.65929D0, 0.D0, 0.D0/)         ! C16
      RBSITE(17,:) = (/ 2.32543D0, 6.42907D0, 0.D0/)    ! H1
      RBSITE(18,:) = (/ 6.38694D0, 4.06382D0, 0.D0/)    ! H2
      RBSITE(19,:) = (/ 2.32543D0,-6.42907D0, 0.D0/)    ! H3
      RBSITE(20,:) = (/ 6.38694D0,-4.06382D0, 0.D0/)    ! H4
      RBSITE(21,:) = (/-2.32543D0,-6.42907D0, 0.D0/)    ! H5
      RBSITE(22,:) = (/-6.38694D0,-4.06382D0, 0.D0/)    ! H6
      RBSITE(23,:) = (/-2.32543D0, 6.42907D0, 0.D0/)    ! H7
      RBSITE(24,:) = (/-6.38694D0, 4.06382D0, 0.D0/)    ! H8
      RBSITE(25,:) = (/ 8.71284D0, 0.D0, 0.D0/)         ! H9
      RBSITE(26,:) = (/-8.71284D0, 0.D0, 0.D0/)         ! H10

!      RBSITE(:,:) = RBSITE(:,:)*0.529177D0              ! au to angstrom

      STCHRG(1)  =-0.04275D0 
      STCHRG(2)  =-0.04275D0
      STCHRG(3)  = 0.22339D0
      STCHRG(4)  = 0.22339D0
      STCHRG(5)  = 0.22339D0
      STCHRG(6)  = 0.22339D0
      STCHRG(7)  =-0.24782D0
      STCHRG(8)  =-0.29542D0
      STCHRG(9)  =-0.24782D0
      STCHRG(10) =-0.29542D0
      STCHRG(11) =-0.24782D0
      STCHRG(12) =-0.29542D0
      STCHRG(13) =-0.24782D0
      STCHRG(14) =-0.29542D0
      STCHRG(15) =-0.05466D0
      STCHRG(16) =-0.05466D0
      STCHRG(17) = 0.15533D0
      STCHRG(18) = 0.15109D0
      STCHRG(19) = 0.15533D0
      STCHRG(20) = 0.15109D0
      STCHRG(21) = 0.15533D0
      STCHRG(22) = 0.15109D0
      STCHRG(23) = 0.15533D0
      STCHRG(24) = 0.15109D0
      STCHRG(25) = 0.12425D0
      STCHRG(26) = 0.12425D0

      RBSTLA(1,:)  = RBSITE(2,:)  - RBSITE(1,:)                 ! Z FROM C1 TO C2
      RBSTLA(2,:)  = RBSITE(1,:)  - RBSITE(2,:)                 ! Z FROM C2 TO C1
      RBSTLA(3,:)  = RBSITE(2,:)  - RBSITE(3,:)                 ! Z FROM C3 TO C2
      RBSTLA(4,:)  = RBSITE(2,:)  - RBSITE(4,:)                 ! Z FROM C4 TO C2
      RBSTLA(5,:)  = RBSITE(1,:) - RBSITE(5,:)                  ! Z FROM C5 TO C1
      RBSTLA(6,:)  = RBSITE(1,:) - RBSITE(6,:)                  ! Z FROM C6 TO C1
      RBSTLA(7,:)  = RBSITE(17,:) - RBSITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(8,:)  = RBSITE(18,:) - RBSITE(8,:)                 ! Z FROM C8 TO H2
      RBSTLA(9,:)  = RBSITE(19,:) - RBSITE(9,:)                 ! Z FROM C9 TO H3
      RBSTLA(10,:) = RBSITE(20,:) - RBSITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(11,:) = RBSITE(21,:) - RBSITE(11,:)                ! Z FROM C11 TO H5
      RBSTLA(12,:) = RBSITE(22,:) - RBSITE(12,:)                ! Z FROM C12 TO H6
      RBSTLA(13,:) = RBSITE(23,:) - RBSITE(13,:)                ! Z FROM C13 TO H7
      RBSTLA(14,:) = RBSITE(24,:) - RBSITE(14,:)                ! Z FROM C14 TO H8
      RBSTLA(15,:) = RBSITE(25,:) - RBSITE(15,:)                ! Z FROM C15 TO H9
      RBSTLA(16,:) = RBSITE(26,:) - RBSITE(16,:)                ! Z FROM C16 TO H10
      RBSTLA(17,:) = RBSITE(17,:) - RBSITE(7,:)                 ! Z FROM C7 TO H1
      RBSTLA(18,:) = RBSITE(18,:) - RBSITE(8,:)                 ! Z FROM C8 TO H2
      RBSTLA(19,:) = RBSITE(19,:) - RBSITE(9,:)                 ! Z FROM C9 TO H3
      RBSTLA(20,:) = RBSITE(20,:) - RBSITE(10,:)                ! Z FROM C10 TO H4
      RBSTLA(21,:) = RBSITE(21,:) - RBSITE(11,:)                ! Z FROM C11 TO H5
      RBSTLA(22,:) = RBSITE(22,:) - RBSITE(12,:)                ! Z FROM C12 TO H6
      RBSTLA(23,:) = RBSITE(23,:) - RBSITE(13,:)                ! Z FROM C13 TO H7
      RBSTLA(24,:) = RBSITE(24,:) - RBSITE(14,:)                ! Z FROM C14 TO H8
      RBSTLA(25,:) = RBSITE(25,:) - RBSITE(15,:)                ! Z FROM C15 TO H9
      RBSTLA(26,:) = RBSITE(26,:) - RBSITE(16,:)                ! Z FROM C16 TO H10

      DO J1 = 1, NRBSITES

         RBSTLA(J1,:)   = RBSTLA(J1,:)/DSQRT(DOT_PRODUCT(RBSTLA(J1,:),RBSTLA(J1,:)))

      ENDDO

      END SUBROUTINE DEFPYRENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEW(X)

      USE COMMONS

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J5
      DOUBLE PRECISION :: X(3*NATOMS), P1(3), P2(3), Q(3*NRBSITES), RMI(3,3), DRMI(3,3)
      LOGICAL          :: GTEST
      CHARACTER (LEN = 2), DIMENSION(NRBSITES) :: AL   

      AL(1:6)   = 'C'
      AL(7:12)  = 'H' 

      GTEST = .FALSE.

      OPEN(UNIT = 11, FILE = 'config1.xyz', STATUS = 'UNKNOWN')

      WRITE(11, *) (NATOMS/2)*NRBSITES
      WRITE(11, *)

      DO J1 = 1, NATOMS/2

         J3    = 3*J1
         J5    = 3*NATOMS/2 + J3
         P1(:) = X(J5-2:J5)

         CALL RMDRVT(P1(:), RMI(:,:), DRMI(:,:), DRMI(:,:), DRMI(:,:), GTEST)

         DO J2 = 1, NRBSITES

            J5 = 3*J2
            Q(3*J2-2:3*J2) = X(J3-2:J3) + MATMUL(RMI(:,:),RBSITE(J2,:))
            WRITE(11,'(A2,3F20.10,2X,A12,2X,3F20.10,2X,A12,2X,3F20.10,2X,A12,2X,3F20.10)')    &
                 AL(J2), Q(3*J2-2), Q(3*J2-1), Q(3*J2),                         &
                 'atom_vector', RBSTLA(J5-2,1), RBSTLA(J5-2,2), RBSTLA(J5-2,3), &
                 'atom_vector', RBSTLA(J5-1,1), RBSTLA(J5-1,2), RBSTLA(J5-1,3), &
                 'atom_vector', RBSTLA(J5,1), RBSTLA(J5,2), RBSTLA(J5,3)  
 
         ENDDO

      ENDDO

      CLOSE(11)

      END SUBROUTINE VIEW

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DAMP(ABSRIJ, DMPFCT, DDMPDR)

      IMPLICIT NONE

      INTEGER :: K, FCT(6)
      DOUBLE PRECISION, PARAMETER :: B = 1.648520924D0  !3.1152D0
      DOUBLE PRECISION :: ABSRIJ, EXPFCT, DMPFCT, DDMPDR 
      LOGICAL          :: GTEST

      DMPFCT = 1.D0
      DDMPDR = B

      FCT(1) = 1; FCT(2) = 2; FCT(3) = 6; FCT(4) = 24; FCT(5) = 120; FCT(6) = 720

      DO K = 1, 6

         DMPFCT = DMPFCT + (B*ABSRIJ)**K/DBLE(FCT(K))
         IF (K > 1) DDMPDR = DDMPDR + (B**K)*(ABSRIJ)**(K-1)/DBLE(FCT(K-1))

      END DO 

      EXPFCT = DEXP(-B*ABSRIJ)
      DDMPDR = (B*EXPFCT*DMPFCT - EXPFCT*DDMPDR)/ABSRIJ
      DMPFCT = 1.D0 - EXPFCT*DMPFCT

      END SUBROUTINE DAMP     
