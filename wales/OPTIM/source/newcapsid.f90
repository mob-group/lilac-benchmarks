      SUBROUTINE NEWCAPSID (X, G, ENERGY, GTEST, SECT)

      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
 
      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), HESS(3*NATOMS,3*NATOMS)
      DOUBLE PRECISION :: RISITE(3), RJSITE(3), RI(3), RJ(3), DSS(3), PI(3)
      DOUBLE PRECISION :: DSS2, RSS, DVDR, ENERGY, SIGMASQ, VSS
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: CAPRHO, CAPEPS2, CAPRAD, CAPHEIGHT
      LOGICAL          :: GTEST, SECT
      COMMON / CAPS /  CAPRHO, CAPEPS2, CAPRAD, CAPHEIGHT

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      
      OFFSET  = 3*NATOMS/2
      SIGMASQ = (1.D0 + CAPRAD*SQRT(0.5D0*(5.D0 + SQRT(5.D0))))**2
     
      CALL DEFCAPSID(CAPRAD,CAPHEIGHT)

      DO J1 = 1, NATOMS/2

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         PI = X(J5-2:J5)

         CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)
 
         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,RBSITE(J2,:))

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, NATOMS/2 - 1 

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, NATOMS/2

            J4 = 3*J2
            J6 = OFFSET + J4

!     REPULSIVE SITE CONTRIBUTION

            J7     = NRBSITES*(J1-1) + 6
            J8     = NRBSITES*(J2-1) + 6
            DSS(:) = R(J7,:) - R(J8,:)
            DSS2   = DOT_PRODUCT(DSS,DSS)
              
            VSS    = CAPEPS2*(SIGMASQ/DSS2)**6
            ENERGY = ENERGY + VSS

            IF (GTEST) THEN

               DVDR   = -12.D0*VSS/DSS2
               G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

               G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
               G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
               G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

               G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
               G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
               G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

            ENDIF

!     SUM OVER MORSE SITES

            DO I = 1, NRBSITES - 1

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES - 1

                  J8 = NRBSITES*(J2-1) + J
                  DSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DOT_PRODUCT(DSS(:),DSS(:))
                  RSS    = SQRT(DSS2)
                  VSS    = EXP(CAPRHO*(1.D0 - RSS))
                  ENERGY = ENERGY + (1.D0 - VSS)*(1.D0 - VSS) - 1.D0 

                  IF (GTEST) THEN
 
                     DVDR       = 2.D0*CAPRHO*(-VSS*VSS + VSS)/RSS
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

      ENDDO

      END SUBROUTINE NEWCAPSID 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFCAPSID(CAPRAD,CAPHEIGHT)
    
      USE COMMONS, ONLY: RBSITE  
      
      IMPLICIT NONE
  
      DOUBLE PRECISION :: CAPRAD, CAPHEIGHT

      RBSITE(1,1) = 1.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.D0
      RBSITE(2,1) = ((-1.D0 + SQRT(5.D0)))/4.D0
      RBSITE(2,2) = (SQRT((5.D0 + SQRT(5.D0))/2.D0))/2.D0
      RBSITE(2,3) = 0.D0
      RBSITE(3,1) = ((-1.D0 - SQRT(5.D0)))/4.0D0
      RBSITE(3,2) = (SQRT((5.D0 - SQRT(5.D0))/2.D0))/2.D0
      RBSITE(3,3) = 0.0D0
      RBSITE(4,1) = ((-1 - SQRT(5.D0)))/4.D0
      RBSITE(4,2) = -(SQRT((5.D0 - SQRT(5.D0))/2.D0))/2.D0
      RBSITE(4,3) = 0.D0
      RBSITE(5,1) = ((-1 + SQRT(5.D0)))/4.D0
      RBSITE(5,2) = -(SQRT((5.D0 + SQRT(5.D0))/2.D0))/2.D0
      RBSITE(5,3) = 0.0D0
      RBSITE(6,1) = 0.0D0
      RBSITE(6,2) = 0.0D0
      RBSITE(6,3) = CAPHEIGHT

!      RBSITE(1:5,:) = CAPRAD*RBSITE(1:5,:)
      RBSITE(:,:)   = CAPRAD*RBSITE(:,:)

      END SUBROUTINE DEFCAPSID
