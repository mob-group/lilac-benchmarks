      SUBROUTINE LWOTPBOX (X, G, ENERGY, GTEST, SECT)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, PARAM1, PARAM2, PARAM4
      USE KEY

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R6, R12, DUMMY, DSS2, D2CMCM
      DOUBLE PRECISION :: LWCNSTA, LWCNSTB, LWRCUT, LWRCUTSQ, LWRCUT2SQ, BOXL
      DOUBLE PRECISION :: RI(3), RJ(3), DR(3), RSS(3), P(3), R(NRBSITES*NATOMS/2,3)
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
      LOGICAL          :: GTEST, SECT

      !js850> a cubic box is assumed, so ignore PARAM2 and PARAM3
      BOXL   = PARAM1
      LWRCUT = PARAM4

      !js850> 12/04/2012 the PARAMS keyword for bulk OTP was changed from 
      ! PARAMS boxl cutoff
      !to
      ! PARAMS boxlx boxly boxlz cutoff
      !add some check to catch odata files that use the old format
      IF ( ABS(PARAM1 - PARAM2) .GT. 1E-5 ) THEN
         WRITE(*,*) "lwotpbox> WARNING: periodic box is not cubic, but this is not implimented."
         write(*,*) "lwotpbox>          Perhaps you're using the old style PARAMS keyword with only 2 arguments?"
      ENDIF
      IF ( LWRCUT .LT. 0.5 .OR. LWRCUT .GT. 10. ) THEN
         WRITE(*,*) "lwotpbox> WARNING: cutoff seems too large or too small."
         write(*,*) "lwotpbox>          Perhaps you're using the old style PARAMS keyword with only 2 arguments?"
      ENDIF

      CALL DEFLWOTPBX(LWCNSTA, LWCNSTB, LWRCUT, LWRCUTSQ, LWRCUT2SQ)

      ENERGY  = 0.D0

      IF (GTEST) THEN
         G(:) = 0.D0; DVDR(:,:) = 0.D0 
         DOTI1(:,:) = 0.D0; DOTI2(:,:) = 0.D0; DOTI3(:,:) = 0.D0
         DOTJ1(:,:) = 0.D0; DOTJ1(:,:) = 0.D0; DOTJ3(:,:) = 0.D0
      ENDIF

      IF (SECT) THEN
         HESS(:,:) = 0.D0 
         D2VDR2(:,:) = 0.D0
      ENDIF

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, SECT)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = MATMUL(RMI,RBSITE(J2,:))

            IF (GTEST .OR. SECT) THEN
 
               DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))

            ENDIF

            IF (SECT) THEN

               D2R1(J4,:) = MATMUL(D2RMI1,RBSITE(J2,:))
               D2R2(J4,:) = MATMUL(D2RMI2,RBSITE(J2,:))
               D2R3(J4,:) = MATMUL(D2RMI3,RBSITE(J2,:))

               D2R12(J4,:) = MATMUL(D2RMI12,RBSITE(J2,:))
               D2R23(J4,:) = MATMUL(D2RMI23,RBSITE(J2,:))
               D2R31(J4,:) = MATMUL(D2RMI31,RBSITE(J2,:))

            ENDIF

         ENDDO

      ENDDO


      DO J1 = 1, REALNATOMS - 1 

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4
            RJ = X(J4-2:J4)

            DR(:)  = RI(:) - RJ(:) - ANINT((RI(:) - RJ(:))/BOXL)*BOXL
            D2CMCM = DOT_PRODUCT(DR(:),DR(:))

            IF (D2CMCM < LWRCUT2SQ) THEN 

               DO I = 1, NRBSITES 

                  J7 = NRBSITES*(J1-1) + I

                  DO J = 1, NRBSITES 

                     J8     = NRBSITES*(J2-1) + J
                     RSS(:) = DR(:) + R(J7,:) - R(J8,:)
                     DSS2   = DOT_PRODUCT(RSS(:),RSS(:))

                     IF (DSS2 < LWRCUTSQ) THEN

                        R2     = 1.D0/DSS2
                        R6     = R2*R2*R2
                        R12    = R6*R6
                        ENERGY = ENERGY + R12 - R6 + LWCNSTA*DSS2 + LWCNSTB

                        IF (GTEST .OR. SECT) THEN
!     DVDR = DVDR/R
                           DVDR(J7,J8) = -(6.D0*R12 - 3.D0*R6)*R2 + LWCNSTA
                           DVDR(J8,J7) = DVDR(J7,J8)

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

                        ENDIF

                        IF (SECT) THEN

                           D2VDR2(J7,J8) = 84.D0*R12*R2*R2 - 24.D0*R6*R2*R2
                           D2VDR2(J8,J7) = D2VDR2(J7,J8)
                           DOTI1(J8,J7)  = -DOTJ1(J7,J8)
                           DOTI2(J8,J7)  = -DOTJ2(J7,J8)
                           DOTI3(J8,J7)  = -DOTJ3(J7,J8)
                           DOTJ1(J8,J7)  = -DOTI1(J7,J8)
                           DOTJ2(J8,J7)  = -DOTI2(J7,J8)
                           DOTJ3(J8,J7)  = -DOTI3(J7,J8)

                        ENDIF

                     ENDIF

                  ENDDO

               ENDDO

            ENDIF

         ENDDO

      ENDDO

      IF (SECT) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3
            RI(:) = X(J3-2:J3)

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4
               RJ(:) = X(J4-2:J4)

               DR(:)  = RI(:) - RJ(:) - ANINT((RI(:) - RJ(:))/BOXL)*BOXL
               D2CMCM = DOT_PRODUCT(DR(:),DR(:))

               IF (D2CMCM < LWRCUT2SQ) THEN

               DO I = 1, NRBSITES

                  J7 = NRBSITES*(J1 - 1) + I

                  DO J = 1, NRBSITES 

                     J8 = NRBSITES*(J2 - 1) + J

                     RSS(:) = DR(:) + R(J7,:) - R(J8,:) 
                     DSS2   = DOT_PRODUCT(RSS(:),RSS(:))

                     IF (DSS2 < LWRCUTSQ) THEN
                 
       
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

                  ENDIF

                  ENDDO

               ENDDO

               ENDIF

            ENDDO

         ENDDO

      ENDIF

      ENERGY    = 4.D0*ENERGY
      IF (GTEST .OR. SECT) G(:) = 8.D0*G(:)
      IF (SECT) HESS(:,:) = 8.D0*HESS(:,:)

      END SUBROUTINE LWOTPBOX 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFLWOTPBX(LWCNSTA, LWCNSTB, LWRCUT, LWRCUTSQ, LWRCUT2SQ)

      USE COMMONS, ONLY: NRBSITES, RBSITE

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: LWCNSTA, LWCNSTB, LWRCUT, LWRCUTSQ, LWRCUT2SQ
      DOUBLE PRECISION :: PI, RCUT2, RCUT6, RCUT12, DELRC, DIST(3)

      LWRCUTSQ = LWRCUT*LWRCUT
      RCUT6    = (1.D0/LWRCUT)**6
      RCUT12   = RCUT6*RCUT6
      LWCNSTA  = (6.D0*RCUT12 - 3.D0*RCUT6)/LWRCUTSQ
      LWCNSTB  = 4.D0*RCUT6 - 7.D0*RCUT12

      PI       = 4.D0*DATAN(1.D0)

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = - 2.D0 * DSIN(7.D0*PI/24.D0) / 3.D0
      RBSITE(1,3) = 0.D0

      RBSITE(2,1) = DCOS(7.D0*PI/24.D0)
      RBSITE(2,2) = DSIN(7.D0*PI/24.D0) / 3.D0
      RBSITE(2,3) = 0.D0

      RBSITE(3,1) = - DCOS(7.D0*PI/24.D0)
      RBSITE(3,2) = DSIN(7.D0*PI/24.D0) / 3.D0
      RBSITE(3,3) = 0.D0

      DO I = 1, NRBSITES 

         DIST(I) = DSQRT(DOT_PRODUCT(RBSITE(I,:),RBSITE(I,:)))

      ENDDO

      DELRC      = 2.D0*MAXVAL(DIST)
      RCUT2      = LWRCUT + DELRC
      LWRCUT2SQ  = RCUT2*RCUT2

      END SUBROUTINE DEFLWOTPBX
