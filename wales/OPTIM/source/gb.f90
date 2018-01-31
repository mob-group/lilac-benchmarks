      SUBROUTINE GB (NATOMS, X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE KEY, ONLY: GBSIGNOT, GBEPSNOT, GBMU, GBNU, GBCHI, GBCHIPRM 

      IMPLICIT NONE

      INTEGER          :: NATOMS, I, J, J1, J2, J3, J4, J5, J6, K, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: SCSIG, SCSIG3, ENERGY
      DOUBLE PRECISION :: EPS1, EPS2, EPS
      DOUBLE PRECISION :: FCT1, FCT2, FCT3, FCT4, FCT5, FCT6, FCT7, FCT8
      DOUBLE PRECISION :: FCT9, FCT10, FCT11, FCT12, FCT13, FCT14, FCT15
      DOUBLE PRECISION :: FCT16, FCT17, FCT18, FCT19, FCT20, FCT21
      DOUBLE PRECISION :: FCT22, FCT23, FCT24
      DOUBLE PRECISION :: FCT3P4, FCT3M4, FCT7P8, FCT7M8
      DOUBLE PRECISION :: ALP, BET, GAM, APB, AMB, DSIGDA, DSIGDB
      DOUBLE PRECISION :: RIJSQ, R2, ABSRIJ, INVR, SCR 
      DOUBLE PRECISION :: SRM1, SRM2, SRM6, SRM7, SRM12, SRM13, SR12M6
      DOUBLE PRECISION :: VR, VA, VB, VG, FIJ(3), FIJN, FIJEI, FIJEJ 
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NRSS(3), P(3), OST(3), EI(3), EJ(3), RMI(3,3), RMJ(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DEDPI1(3), DEDPI2(3), DEDPI3(3), DEDPJ1(3), DEDPJ2(3), DEDPJ3(3)
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      DOUBLE PRECISION :: DVRDR, DVRDA, DVRDB, DVRDG, DVADR, DVBDR, D2VDA2, D2VDB2
      DOUBLE PRECISION :: DVADB, DVBDA, DADR(3), DBDR(3), D2ADX2(3), D2BDX2(3)
      DOUBLE PRECISION :: D2ADYX, D2ADZY, D2ADXZ, D2BDYX, D2BDZY, D2BDXZ
      DOUBLE PRECISION :: D2VDG2, DVADG, DVGDA, DVGDB, DVBDG
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3)
      DOUBLE PRECISION :: D2RMJ1(3,3), D2RMJ2(3,3), D2RMJ3(3,3)
      DOUBLE PRECISION :: D2RI12(3,3), D2RI23(3,3), D2RI31(3,3)
      DOUBLE PRECISION :: D2RJ12(3,3), D2RJ23(3,3), D2RJ31(3,3)
      DOUBLE PRECISION :: D2API1, D2API2, D2API3, D2BPJ1, D2BPJ2, D2BPJ3
      DOUBLE PRECISION :: D2GPI1, D2GPI2, D2GPI3, D2GPJ1, D2GPJ2, D2GPJ3
      DOUBLE PRECISION :: DUMMY
      LOGICAL          :: GTEST, STEST

      OST = (/0.D0, 0.D0, 1.D0/)
      REALNATOMS = NATOMS/2
      OFFSET = 3*REALNATOMS

      IF (GTEST .AND. .NOT. STEST) THEN

         ENERGY    = 0.D0
         G(:)      = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               IF (J1 == J2) CYCLE
               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4)
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDRVT (P, RMJ, DRMJ1, DRMJ2, DRMJ3, GTEST)

!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NRSS   = RIJ / ABSRIJ
               SCR    = ABSRIJ/GBSIGNOT
               INVR   = 1.D0/ABSRIJ
               R2     = 1.D0/RIJSQ

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

               EI  = MATMUL(RMI,OST)
               EJ  = MATMUL(RMJ,OST)
               ALP = DOT_PRODUCT(NRSS,EI)
               BET = DOT_PRODUCT(NRSS,EJ)
               GAM = DOT_PRODUCT(EI,EJ)

!     CALCULATE USEFUL QUANTITIES

               APB    = ALP+BET
               AMB    = ALP-BET

               FCT1   = 1.D0/(1.D0+GBCHI*GAM)
               FCT2   = 1.D0/(1.D0-GBCHI*GAM)
               FCT3   = APB*FCT1
               FCT4   = AMB*FCT2
               FCT3P4 = FCT3+FCT4
               FCT3M4 = FCT3-FCT4

               FCT5   = 1.D0/(1.D0+GBCHIPRM*GAM)
               FCT6   = 1.D0/(1.D0-GBCHIPRM*GAM)
               FCT7   = (ALP+BET)*FCT5
               FCT8   = (ALP-BET)*FCT6
               FCT7P8 = FCT7+FCT8
               FCT7M8 = FCT7-FCT8

!     CALCULATE $\epsilon$

               EPS1   = DSQRT(FCT1*FCT2)
               EPS2   = 1.D0-0.5D0*GBCHIPRM*(APB*FCT7+AMB*FCT8)
               EPS    = GBEPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\sigma/\sigma_{0})^3$

               SCSIG  = 1.d0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
               SCSIG3 = SCSIG*SCSIG*SCSIG

!     CALCULATE DEL(V)/DEL(R)

               SRM1   = 1.D0/(SCR-SCSIG+1.D0)
               SRM2   = SRM1*SRM1
               SRM6   = SRM2*SRM2*SRM2
               SRM7   = SRM6*SRM1
               SRM12  = SRM6*SRM6
               SRM13  = SRM12*SRM1
               SR12M6 = SRM12-SRM6
               FCT9   = 2.D0*SRM13-SRM7
               VR     = -(24.D0/GBSIGNOT)*EPS*FCT9

!     CALCULATE ENERGY

               ENERGY = ENERGY + EPS*SR12M6

!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA)

               FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(GBSIGNOT*EPS2)
               FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
               FCT12  = 12.D0*EPS*GBCHI*SCSIG3

               VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
               VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)

               FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
               FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2
               FCT15  = FCT13 + FCT14
               FCT16  = FCT3*FCT3 - FCT4*FCT4
               FCT17  = 4.D0*EPS*FCT15
               FCT18  = 6.D0*EPS*GBCHI*GBCHI*SCSIG3*FCT16
               FCT19  = FCT17*SR12M6
               FCT20  = FCT18*FCT9
               FCT21  = 26.D0*SRM12*SRM2-7.D0*SRM6*SRM2

               VG     = FCT19 - FCT20

               DEDPI1 = MATMUL(DRMI1,OST)
               DEDPI2 = MATMUL(DRMI2,OST)
               DEDPI3 = MATMUL(DRMI3,OST)

               DEDPJ1 = MATMUL(DRMJ1,OST)
               DEDPJ2 = MATMUL(DRMJ2,OST)
               DEDPJ3 = MATMUL(DRMJ3,OST)

               DADPI1 = DOT_PRODUCT(NRSS,DEDPI1)
               DADPI2 = DOT_PRODUCT(NRSS,DEDPI2)
               DADPI3 = DOT_PRODUCT(NRSS,DEDPI3)

               DBDPJ1 = DOT_PRODUCT(NRSS,DEDPJ1)
               DBDPJ2 = DOT_PRODUCT(NRSS,DEDPJ2)
               DBDPJ3 = DOT_PRODUCT(NRSS,DEDPJ3)

               DGDPI1 = DOT_PRODUCT(DEDPI1,EJ)
               DGDPI2 = DOT_PRODUCT(DEDPI2,EJ)
               DGDPI3 = DOT_PRODUCT(DEDPI3,EJ)

               DGDPJ1 = DOT_PRODUCT(EI,DEDPJ1)
               DGDPJ2 = DOT_PRODUCT(EI,DEDPJ2)
               DGDPJ3 = DOT_PRODUCT(EI,DEDPJ3)

!     CALCULATE CONTRIBUTION TO FORCES

               FIJN   = VR - (VA*ALP+VB*BET)*INVR
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR

               FIJ  = FIJN*NRSS + FIJEI*EI + FIJEJ*EJ

               G(J3-2:J3) = G(J3-2:J3) + FIJ
               G(J4-2:J4) = G(J4-2:J4) - FIJ

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

            ENDDO

         ENDDO

!     CALCULATE THE GAY-BERNE POTENTIAL

         ENERGY = 4.D0 * ENERGY

      ELSE IF (GTEST .AND. STEST) THEN

         ENERGY    = 0.D0
         G(:)      = 0.D0
         HESS(:,:) = 0.D0

         DO J1 = 1, REALNATOMS 

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RI12, D2RI23, D2RI31, GTEST, STEST) 

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE 
               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDFAS(P, RMJ, DRMJ1, DRMJ2, DRMJ3, D2RMJ1, D2RMJ2, D2RMJ3, D2RJ12, D2RJ23, D2RJ31, GTEST, STEST) 
               
!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NRSS   = RIJ / ABSRIJ
               SCR    = ABSRIJ/GBSIGNOT
               INVR   = 1.D0/ABSRIJ
               R2     = 1.D0/RIJSQ

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

               EI  = MATMUL(RMI,OST)
               EJ  = MATMUL(RMJ,OST)
               ALP = DOT_PRODUCT(NRSS,EI) 
               BET = DOT_PRODUCT(NRSS,EJ)
               GAM = DOT_PRODUCT(EI,EJ)

!     CALCULATE USEFUL QUANTITIES

               APB    = ALP+BET
               AMB    = ALP-BET

               FCT1   = 1.D0/(1.D0+GBCHI*GAM)
               FCT2   = 1.D0/(1.D0-GBCHI*GAM)
               FCT3   = APB*FCT1
               FCT4   = AMB*FCT2
               FCT3P4 = FCT3+FCT4
               FCT3M4 = FCT3-FCT4

               FCT5   = 1.D0/(1.D0+GBCHIPRM*GAM)
               FCT6   = 1.D0/(1.D0-GBCHIPRM*GAM)
               FCT7   = (ALP+BET)*FCT5
               FCT8   = (ALP-BET)*FCT6
               FCT7P8 = FCT7+FCT8
               FCT7M8 = FCT7-FCT8

!     CALCULATE $\epsilon$

               EPS1   = DSQRT(FCT1*FCT2)
               EPS2   = 1.D0-0.5D0*GBCHIPRM*(APB*FCT7+AMB*FCT8)
               EPS    = GBEPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\sigma/\sigma_{0})^3$

               SCSIG  = 1.d0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
               SCSIG3 = SCSIG*SCSIG*SCSIG

!     CALCULATE DEL(V)/DEL(R)

               SRM1   = 1.D0/(SCR-SCSIG+1.D0)
               SRM2   = SRM1*SRM1
               SRM6   = SRM2*SRM2*SRM2
               SRM7   = SRM6*SRM1
               SRM12  = SRM6*SRM6
               SRM13  = SRM12*SRM1
               SR12M6 = SRM12-SRM6
               FCT9   = 2.D0*SRM13-SRM7
               VR     = -(24.D0/GBSIGNOT)*EPS*FCT9

!     CALCULATE ENERGY

               ENERGY = ENERGY + EPS*SR12M6

!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA)

               FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(GBSIGNOT*EPS2)
               FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
               FCT12  = 12.D0*EPS*GBCHI*SCSIG3

               VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
               VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)

               FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
               FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2
               FCT15  = FCT13 + FCT14
               FCT16  = FCT3*FCT3 - FCT4*FCT4
               FCT17  = 4.D0*EPS*FCT15
               FCT18  = 6.D0*EPS*GBCHI*GBCHI*SCSIG3*FCT16
               FCT19  = FCT17*SR12M6
               FCT20  = FCT18*FCT9
               FCT21  = 26.D0*SRM12*SRM2-7.D0*SRM6*SRM2

               VG     = FCT19 - FCT20

               DEDPI1 = MATMUL(DRMI1,OST)
               DEDPI2 = MATMUL(DRMI2,OST)
               DEDPI3 = MATMUL(DRMI3,OST)

               DEDPJ1 = MATMUL(DRMJ1,OST)
               DEDPJ2 = MATMUL(DRMJ2,OST)
               DEDPJ3 = MATMUL(DRMJ3,OST)

               DADPI1 = DOT_PRODUCT(NRSS,DEDPI1)
               DADPI2 = DOT_PRODUCT(NRSS,DEDPI2)
               DADPI3 = DOT_PRODUCT(NRSS,DEDPI3)

               DBDPJ1 = DOT_PRODUCT(NRSS,DEDPJ1)
               DBDPJ2 = DOT_PRODUCT(NRSS,DEDPJ2)
               DBDPJ3 = DOT_PRODUCT(NRSS,DEDPJ3)

               DGDPI1 = DOT_PRODUCT(DEDPI1,EJ)
               DGDPI2 = DOT_PRODUCT(DEDPI2,EJ)
               DGDPI3 = DOT_PRODUCT(DEDPI3,EJ)

               DGDPJ1 = DOT_PRODUCT(EI,DEDPJ1)
               DGDPJ2 = DOT_PRODUCT(EI,DEDPJ2)
               DGDPJ3 = DOT_PRODUCT(EI,DEDPJ3)

!     CALCULATE CONTRIBUTION TO FORCES

               FIJN   = VR - (VA*ALP+VB*BET)*INVR
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR

               FIJ  = FIJN*NRSS + FIJEI*EI + FIJEJ*EJ

               G(J3-2:J3) = G(J3-2:J3) + FIJ
               G(J4-2:J4) = G(J4-2:J4) - FIJ

               G(J5-2) = G(J5-2) + VA*DADPI1 + VG*DGDPI1
               G(J5-1) = G(J5-1) + VA*DADPI2 + VG*DGDPI2
               G(J5)   = G(J5)   + VA*DADPI3 + VG*DGDPI3

               G(J6-2) = G(J6-2) + VB*DBDPJ1 + VG*DGDPJ1
               G(J6-1) = G(J6-1) + VB*DBDPJ2 + VG*DGDPJ2
               G(J6)   = G(J6)   + VB*DBDPJ3 + VG*DGDPJ3

!     HESSIAN CALCULATION

               DSIGDA = 0.5D0*GBSIGNOT*GBCHI*SCSIG3*FCT3P4
               DSIGDB = 0.5D0*GBSIGNOT*GBCHI*SCSIG3*FCT3M4

               DADR(:) = INVR*EI(:) - ALP*INVR*NRSS(:)
               DBDR(:) = INVR*EJ(:) - BET*INVR*NRSS(:)

               DVRDR  = (624.D0*SRM12*SRM2-168.D0*SRM6*SRM2)*EPS / (GBSIGNOT*GBSIGNOT)
               DVRDA  = FCT10*FCT9*FCT7P8-DVRDR*DSIGDA
               DVRDB  = FCT10*FCT9*FCT7M8-DVRDR*DSIGDB


               DVADR  = -FCT12*FCT3P4*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2) & 
                      /GBSIGNOT+6.D0*FCT11*FCT9*FCT7P8/GBSIGNOT

               DVBDR  = -FCT12*FCT3M4*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2) &
                      /GBSIGNOT+6.D0*FCT11*FCT9*FCT7M8/GBSIGNOT

               D2VDA2 = -FCT11*FCT12*FCT9*FCT3P4*FCT7P8/(4.D0*EPS) &
                      +FCT11*GBCHIPRM*(GBMU-1.D0)*SR12M6*FCT7P8**2.D0/EPS2 &
                      +(FCT12*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2)*FCT3P4 &
                      -6.D0*FCT11*FCT9*FCT7P8+3.D0*FCT12*FCT9*FCT3P4 &
                      /SCSIG)*(DSIGDA/GBSIGNOT)-FCT11*SR12M6*(FCT5+FCT6) &
                      +FCT12*FCT9*(FCT1+FCT2)


               D2VDB2 = -FCT11*FCT12*FCT9*FCT3M4*FCT7M8/(4.D0*EPS) &
                      +FCT11*GBCHIPRM*(GBMU-1.D0)*SR12M6*FCT7M8**2.D0/EPS2 &
                      +(FCT12*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2)*FCT3M4 &
                      -6.D0*FCT11*FCT9*FCT7M8+3.D0*FCT12*FCT9*FCT3M4 &
                      /SCSIG)*(DSIGDB/GBSIGNOT)-FCT11*SR12M6*(FCT5+FCT6) &
                      +FCT12*FCT9*(FCT1+FCT2)

               DVADB  = -FCT11*FCT12*FCT9*FCT3P4*FCT7M8/(4.D0*EPS)+ &
                      FCT11*GBCHIPRM*(GBMU-1.D0)*SR12M6*FCT7P8*FCT7M8/EPS2 &
                      +(FCT12*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2)*FCT3P4 &
                      -6.D0*FCT11*FCT9*FCT7P8+3.D0*FCT12*FCT9*FCT3P4 &
                      /SCSIG)*(DSIGDB/GBSIGNOT)+FCT12*FCT9*(FCT1-FCT2) &
                      -FCT11*SR12M6*(FCT5-FCT6)

               DVBDA  = DVADB

               D2ADX2(1) = - 2.D0*R2*INVR*RIJ(1)*EI(1) + 3.D0*ALP*RIJ(1)*RIJ(1)*R2*R2 - ALP*R2
               D2ADX2(2) = - 2.D0*R2*INVR*RIJ(2)*EI(2) + 3.D0*ALP*RIJ(2)*RIJ(2)*R2*R2 - ALP*R2
               D2ADX2(3) = - 2.D0*R2*INVR*RIJ(3)*EI(3)   + 3.D0*ALP*RIJ(3)*RIJ(3)*R2*R2 - ALP*R2

               D2BDX2(1) = - 2.D0*R2*INVR*RIJ(1)*EJ(1) + 3.D0*BET*RIJ(1)*RIJ(1)*R2*R2 - BET*R2
               D2BDX2(2) = - 2.D0*R2*INVR*RIJ(2)*EJ(2) + 3.D0*BET*RIJ(2)*RIJ(2)*R2*R2 - BET*R2
               D2BDX2(3) = - 2.D0*R2*INVR*RIJ(3)*EJ(3) + 3.D0*BET*RIJ(3)*RIJ(3)*R2*R2 - BET*R2

               D2ADYX    = - INVR*R2*(EI(1)*RIJ(2)+EI(2)*RIJ(1)) + 3*ALP*RIJ(1)*RIJ(2)*R2*R2
               D2ADZY    = - INVR*R2*(EI(3)*RIJ(2)+EI(2)*RIJ(3)) + 3*ALP*RIJ(3)*RIJ(2)*R2*R2
               D2ADXZ    = - INVR*R2*(EI(3)*RIJ(1)+EI(1)*RIJ(3)) + 3*ALP*RIJ(1)*RIJ(3)*R2*R2

               D2BDYX    = - INVR*R2*(EJ(1)*RIJ(2)+EJ(2)*RIJ(1)) + 3*BET*RIJ(1)*RIJ(2)*R2*R2
               D2BDZY    = - INVR*R2*(EJ(3)*RIJ(2)+EJ(2)*RIJ(3))  + 3*BET*RIJ(3)*RIJ(2)*R2*R2
               D2BDXZ    = - INVR*R2*(EJ(3)*RIJ(1)+EJ(1)*RIJ(3))  + 3*BET*RIJ(1)*RIJ(3)*R2*R2


               D2VDG2 = FCT19*FCT15-FCT20*FCT15+4.D0*EPS*SR12M6* &
                      (FCT13/GAM+2.D0*FCT13*FCT13/GBNU-0.5D0*FCT14 &
                      *GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2 &
                      -GBMU*GBCHIPRM**3.D0*(FCT7*FCT7*FCT5+FCT8*FCT8*FCT6) &
                      /EPS2)-FCT20*FCT15+FCT18*FCT16*GBCHI*GBCHI*SCSIG3* &
                      (26.D0*SRM12*SRM2-7.D0*SRM6*SRM2)/4.D0+FCT20* &
                      FCT16*SCSIG*SCSIG*3.D0*GBCHI*GBCHI/4.D0 + 12.D0*EPS* &
                      GBCHI*GBCHI*GBCHI*SCSIG3*FCT9*(FCT3*FCT3*FCT1+FCT4*FCT4*FCT2)


               DVGDA  =-FCT19*GBCHIPRM*GBMU*FCT7P8/EPS2+3.D0*GBCHI*SCSIG3 &
                      *FCT3P4*FCT17*FCT9+4.D0*EPS*SR12M6*(GBCHIPRM*FCT14 &
                      *FCT7P8/EPS2+GBMU*GBCHIPRM**2*(FCT7*FCT5-FCT8*FCT6)/ &
                      EPS2)+FCT20*GBCHIPRM*GBMU*FCT7P8/EPS2-0.5D0*GBCHI &
                      *SCSIG3*FCT18*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2)* &
                      FCT3P4-1.5D0*GBCHI*SCSIG*SCSIG*FCT20*FCT3P4- 12.D0* &
                      EPS*GBCHI*GBCHI*SCSIG3*FCT9*(FCT3*FCT1-FCT4*FCT2)

               DVADG  = DVGDA

               DVGDB  = -FCT19*GBCHIPRM*GBMU*FCT7M8/EPS2+3.D0*GBCHI*SCSIG3 &
                      *FCT3M4*FCT17*FCT9+4.D0*EPS*SR12M6*(GBCHIPRM*FCT14 &
                      *FCT7M8/EPS2+GBMU*GBCHIPRM**2*(FCT7*FCT5+FCT8*FCT6)/ &
                      EPS2)+FCT20*GBCHIPRM*GBMU*FCT7M8/EPS2-0.5D0*GBCHI &
                      *SCSIG3*FCT18*(26.D0*SRM12*SRM2-7.D0*SRM6*SRM2)* &
                      FCT3M4-1.5D0*GBCHI*SCSIG*SCSIG*FCT20*FCT3M4- 12.D0* &
                      EPS*GBCHI*GBCHI*SCSIG3*FCT9*(FCT3*FCT1+FCT4*FCT2)

               DVBDG  = DVGDB

               DVRDG  = (-24.D0*EPS*FCT9*FCT15 + FCT21*FCT18)/GBSIGNOT

               DEDPI1 = MATMUL(DRMI1,OST)
               DEDPI2 = MATMUL(DRMI2,OST)
               DEDPI3 = MATMUL(DRMI3,OST)

               DEDPJ1 = MATMUL(DRMJ1,OST)
               DEDPJ2 = MATMUL(DRMJ2,OST)
               DEDPJ3 = MATMUL(DRMJ3,OST)

               DADPI1 = DOT_PRODUCT(NRSS,DEDPI1)
               DADPI2 = DOT_PRODUCT(NRSS,DEDPI2)
               DADPI3 = DOT_PRODUCT(NRSS,DEDPI3)

               DBDPJ1 = DOT_PRODUCT(NRSS,DEDPJ1)
               DBDPJ2 = DOT_PRODUCT(NRSS,DEDPJ2)
               DBDPJ3 = DOT_PRODUCT(NRSS,DEDPJ3)

               DGDPI1 = DOT_PRODUCT(DEDPI1,EJ)
               DGDPI2 = DOT_PRODUCT(DEDPI2,EJ)
               DGDPI3 = DOT_PRODUCT(DEDPI3,EJ)

               DGDPJ1 = DOT_PRODUCT(EI,DEDPJ1)
               DGDPJ2 = DOT_PRODUCT(EI,DEDPJ2)
               DGDPJ3 = DOT_PRODUCT(EI,DEDPJ3)

               D2API1 = DOT_PRODUCT(NRSS,MATMUL(D2RMI1,OST))
               D2API2 = DOT_PRODUCT(NRSS,MATMUL(D2RMI2,OST))
               D2API3 = DOT_PRODUCT(NRSS,MATMUL(D2RMI3,OST))

               D2GPI1 = DOT_PRODUCT(MATMUL(D2RMI1,OST),EJ)
               D2GPI2 = DOT_PRODUCT(MATMUL(D2RMI2,OST),EJ)
               D2GPI3 = DOT_PRODUCT(MATMUL(D2RMI3,OST),EJ)

!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES
!     xi,xi
               HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + DVRDR*R2*RIJ(1)*RIJ(1) &
                               + (DVRDA*DADR(1)+DVRDB*DBDR(1))*INVR*RIJ(1) &
                               + VR*INVR*(1.D0-R2*RIJ(1)*RIJ(1))+(DVADR*INVR*RIJ(1) &
                               + D2VDA2*DADR(1)+DVADB*DBDR(1))*DADR(1) &
                               + VA*D2ADX2(1)+(DVBDR*INVR*RIJ(1) + DVBDA*DADR(1)+D2VDB2*DBDR(1)) &
                               * DBDR(1)  + VB*D2BDX2(1)
!     yi,yi
               HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + DVRDR*R2*RIJ(2)*RIJ(2) &
                               + (DVRDA*DADR(2)+DVRDB*DBDR(2))*INVR*RIJ(2) &
                               + VR*INVR*(1.D0-R2*RIJ(2)*RIJ(2))+(DVADR*INVR*RIJ(2) &
                               + D2VDA2*DADR(2)+DVADB*DBDR(2))*DADR(2) &
                               + VA*D2ADX2(2)+(DVBDR*INVR*RIJ(2) + DVBDA*DADR(2)+D2VDB2*DBDR(2)) &
                               *DBDR(2) + VB*D2BDX2(2)
!     zi,zi
               HESS(J3,J3)     = HESS(J3,J3) + DVRDR*R2*RIJ(3)*RIJ(3) &
                               + (DVRDA*DADR(3)+DVRDB*DBDR(3))*INVR*RIJ(3) &
                               + VR*INVR*(1.D0-R2*RIJ(3)*RIJ(3))+(DVADR*INVR*RIJ(3) &
                               + D2VDA2*DADR(3)+DVADB*DBDR(3))*DADR(3) &
                               + VA*D2ADX2(3)+(DVBDR*INVR*RIJ(3) + DVBDA*DADR(3)+D2VDB2*DBDR(3)) & 
                               * DBDR(3) + VB*D2BDX2(3)
!     pi1,pi1
               HESS(J5-2,J5-2) = HESS(J5-2,J5-2) + D2VDA2*DADPI1*DADPI1 + 2.D0*DVADG*DGDPI1*DADPI1 &
                               + VA*D2API1 + D2VDG2*DGDPI1*DGDPI1 + VG*D2GPI1
!     pi2,pi2
               HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + D2VDA2*DADPI2*DADPI2 + 2.D0*DVADG*DGDPI2*DADPI2 &
                               + VA*D2API2 + D2VDG2*DGDPI2*DGDPI2 + VG*D2GPI2
!     pi3,pi3
               HESS(J5,J5)     = HESS(J5,J5)     + D2VDA2*DADPI3*DADPI3 + 2.D0*DVADG*DGDPI3*DADPI3 &
                               + VA*D2API3 + D2VDG2*DGDPI3*DGDPI3 + VG*D2GPI3

!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCK: SAME MOLECULE, DIFFERENT COORDINATES

!     xi,yi
               DUMMY = (DVRDR*INVR*RIJ(2) + DVRDA*DADR(2) + DVRDB*DBDR(2))*INVR*RIJ(1) &
                     - VR*INVR*R2*RIJ(1)*RIJ(2) + (DVADR*INVR*RIJ(2) + D2VDA2*DADR(2) &
                     + DVADB*DBDR(2))*DADR(1) + VA*D2ADYX + (DVBDR*INVR*RIJ(2) + DVBDA*DADR(2) &
                     + D2VDB2*DBDR(2))*DBDR(1) + VB*D2BDYX
               HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
               HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
               DUMMY = (DVRDR*INVR*RIJ(3) + DVRDA*DADR(3) + DVRDB*DBDR(3))*INVR*RIJ(2) & 
                     - VR*INVR*R2*RIJ(2)*RIJ(3) + (DVADR*INVR*RIJ(3) + D2VDA2*DADR(3) & 
                     + DVADB*DBDR(3))*DADR(2) + VA*D2ADZY + (DVBDR*INVR*RIJ(3) + DVBDA*DADR(3) & 
                     + D2VDB2*DBDR(3))*DBDR(2) + VB*D2BDZY
               HESS(J3-1,J3) = HESS(J3-1,J3) + DUMMY
               HESS(J3,J3-1) = HESS(J3,J3-1) + DUMMY
!     zi,xi
               DUMMY = (DVRDR*INVR*RIJ(1) + DVRDA*DADR(1) + DVRDB*DBDR(1))*INVR*RIJ(3) & 
                     - VR*INVR*R2*RIJ(3)*RIJ(1) + (DVADR*INVR*RIJ(1) + D2VDA2*DADR(1) &
                     + DVADB*DBDR(1))*DADR(3) + VA*D2ADXZ + (DVBDR*INVR*RIJ(1) + DVBDA*DADR(1) &
                     + D2VDB2*DBDR(1))*DBDR(3) + VB*D2BDXZ
               HESS(J3-2,J3) = HESS(J3-2,J3) + DUMMY
               HESS(J3,J3-2) = HESS(J3,J3-2) + DUMMY
               
               FCT22 = DVRDA*DADPI1 + DVRDG*DGDPI1 - ((D2VDA2*DADPI1 + DVADG*DGDPI1)*ALP &
                     + VA*DADPI1 + (DVBDA*DADPI1 + DVBDG*DGDPI1)*BET)*INVR 
               FCT23 = (D2VDA2*DADPI1 + DVADG*DGDPI1)*INVR
               FCT24 = (DVBDA*DADPI1 + DVBDG*DGDPI1)*INVR
!     xi,pi1
               DUMMY = FCT22*NRSS(1) + FCT23*EI(1) + VA*INVR*DEDPI1(1) + FCT24*EJ(1)
               HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
               HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1   
               DUMMY = FCT22*NRSS(2) + FCT23*EI(2) + VA*INVR*DEDPI1(2) + FCT24*EJ(2)
               HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
               HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
               DUMMY = FCT22*NRSS(3) + FCT23*EI(3) + VA*INVR*DEDPI1(3) + FCT24*EJ(3)
               HESS(J3,J5-2)   = HESS(J3,J5-2)   + DUMMY
               HESS(J5-2,J3)   = HESS(J5-2,J3)   + DUMMY
               
               FCT22 = DVRDA*DADPI2 + DVRDG*DGDPI2 - ((D2VDA2*DADPI2 + DVADG*DGDPI2)*ALP &
                     + VA*DADPI2 + (DVBDA*DADPI2 + DVBDG*DGDPI2)*BET)*INVR 
               FCT23 = (D2VDA2*DADPI2 + DVADG*DGDPI2)*INVR
               FCT24 = (DVBDA*DADPI2 + DVBDG*DGDPI2)*INVR
!     xi,pi2
               DUMMY = FCT22*NRSS(1) + FCT23*EI(1) + VA*INVR*DEDPI2(1) + FCT24*EJ(1)
               HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
               HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2   
               DUMMY = FCT22*NRSS(2) + FCT23*EI(2) + VA*INVR*DEDPI2(2) + FCT24*EJ(2)
               HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
               HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
               DUMMY = FCT22*NRSS(3) + FCT23*EI(3) + VA*INVR*DEDPI2(3) + FCT24*EJ(3)
               HESS(J3,J5-1)   = HESS(J3,J5-1)   + DUMMY
               HESS(J5-1,J3)   = HESS(J5-1,J3)   + DUMMY
               
               FCT22 = DVRDA*DADPI3 + DVRDG*DGDPI3 - ((D2VDA2*DADPI3 + DVADG*DGDPI3)*ALP &
                     + VA*DADPI3 + (DVBDA*DADPI3 + DVBDG*DGDPI3)*BET)*INVR 
               FCT23 = (D2VDA2*DADPI3 + DVADG*DGDPI3)*INVR
               FCT24 = (DVBDA*DADPI3 + DVBDG*DGDPI3)*INVR
!     xi,pi2
               DUMMY = FCT22*NRSS(1) + FCT23*EI(1) + VA*INVR*DEDPI3(1) + FCT24*EJ(1)
               HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
               HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi2   
               DUMMY = FCT22*NRSS(2) + FCT23*EI(2) + VA*INVR*DEDPI3(2) + FCT24*EJ(2)
               HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
               HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi2
               DUMMY = FCT22*NRSS(3) + FCT23*EI(3) + VA*INVR*DEDPI3(3) + FCT24*EJ(3)
               HESS(J3,J5)     = HESS(J3,J5)   + DUMMY
               HESS(J5,J3)   = HESS(J5,J3)   + DUMMY
!     pi1,pi2
               DUMMY = D2VDA2*DADPI2*DADPI1 + DVADG*DGDPI2*DADPI1 &
                     + VA*DOT_PRODUCT(NRSS,MATMUL(D2RI12,OST)) + D2VDG2*DGDPI2*DGDPI1 &
                     + DVGDA*DADPI2*DGDPI1 + VG*DOT_PRODUCT(MATMUL(D2RI12,OST),EJ)
               HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
               HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
               DUMMY = D2VDA2*DADPI3*DADPI2 + DVADG*DGDPI3*DADPI2 &
                     + VA*DOT_PRODUCT(NRSS,MATMUL(D2RI23,OST)) + D2VDG2*DGDPI3*DGDPI2 &
                     + DVGDA*DADPI3*DGDPI2 + VG*DOT_PRODUCT(MATMUL(D2RI23,OST),EJ)
               HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
               HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
               DUMMY = D2VDA2*DADPI1*DADPI3 + DVADG*DGDPI1*DADPI3 &
                     + VA*DOT_PRODUCT(NRSS,MATMUL(D2RI31,OST)) + D2VDG2*DGDPI1*DGDPI3 &
                     + DVGDA*DADPI1*DGDPI3 + VG*DOT_PRODUCT(MATMUL(D2RI31,OST),EJ)
               HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
               HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE

!     xi,xj
               HESS(J3-2,J4-2) = -DVRDR*R2*RIJ(1)*RIJ(1) - (DVRDA*DADR(1) &
                               + DVRDB*DBDR(1))*INVR*RIJ(1) - VR*INVR*(1.D0 - R2*RIJ(1)*RIJ(1)) &
                               - (DVADR*INVR*RIJ(1) + D2VDA2*DADR(1) + DVADB*DBDR(1))*DADR(1) &
                               - VA*D2ADX2(1) - (DVBDR*INVR*RIJ(1) + DVBDA*DADR(1) &
                               + D2VDB2*DBDR(1))*DBDR(1) - VB*D2BDX2(1)
!     yi,yj
               HESS(J3-1,J4-1) = - DVRDR*R2*RIJ(2)*RIJ(2) - (DVRDA*DADR(2) &
                               + DVRDB*DBDR(2))*INVR*RIJ(2) - VR*INVR*(1.D0 - R2*RIJ(2)*RIJ(2)) &
                               - (DVADR*INVR*RIJ(2) + D2VDA2*DADR(2) + DVADB*DBDR(2))*DADR(2) &
                               - VA*D2ADX2(2) - (DVBDR*INVR*RIJ(2) + DVBDA*DADR(2) &
                               + D2VDB2*DBDR(2))*DBDR(2) - VB*D2BDX2(2) 
!     zi,zj
               HESS(J3,J4)     = - DVRDR*R2*RIJ(3)*RIJ(3) - (DVRDA*DADR(3) &
                               + DVRDB*DBDR(3))*INVR*RIJ(3) - VR*INVR*(1.D0 - R2*RIJ(3)*RIJ(3)) &
                               - (DVADR*INVR*RIJ(3) + D2VDA2*DADR(3) + DVADB*DBDR(3))*DADR(3) &
                               - VA*D2ADX2(3) - (DVBDR*INVR*RIJ(3) + DVBDA*DADR(3) &
                               + D2VDB2*DBDR(3))*DBDR(3) - VB*D2BDX2(3) 
!     pi1,pj1
               HESS(J5-2,J6-2) = HESS(J5-2,J6-2) + (DVADB*DBDPJ1 + DVADG*DGDPJ1)*DADPI1 &
                               + (DVGDB*DBDPJ1 + D2VDG2*DGDPJ1)*DGDPI1 & 
                               + VG*DOT_PRODUCT(DEDPI1,DEDPJ1)
!     pi2,pj2
               HESS(J5-1,J6-1) = HESS(J5-1,J6-1) + (DVADB*DBDPJ2 + DVADG*DGDPJ2)*DADPI2 &
                               + (DVGDB*DBDPJ2 + D2VDG2*DGDPJ2)*DGDPI2 &
                               + VG*DOT_PRODUCT(DEDPI2,DEDPJ2)
!     pi3,pj3                
               HESS(J5,J6)     = HESS(J5,J6)     + (DVADB*DBDPJ3 + DVADG*DGDPJ3)*DADPI3 &
                               + (DVGDB*DBDPJ3 + D2VDG2*DGDPJ3)*DGDPI3 &
                               + VG*DOT_PRODUCT(DEDPI3,DEDPJ3)

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES

!     xi,yj 
               HESS(J3-2,J4-1) = -(DVRDR*INVR*RIJ(2) + DVRDA*DADR(2) + DVRDB*DBDR(2))*INVR*RIJ(1) &
                               + VR*INVR*R2*RIJ(1)*RIJ(2) - (DVADR*INVR*RIJ(2) + D2VDA2*DADR(2) &
                               + DVADB*DBDR(2))*DADR(1) - VA*D2ADYX - (DVBDR*INVR*RIJ(2) &
                               + DVBDA*DADR(2) + D2VDB2*DBDR(2))*DBDR(1) - VB*D2BDYX
               HESS(J3-1,J4-2) = HESS(J3-2,J4-1)
!     yi,zj 
               HESS(J3-1,J4)   = -(DVRDR*INVR*RIJ(3) + DVRDA*DADR(3) + DVRDB*DBDR(3))*INVR*RIJ(2) &
                               + VR*INVR*R2*RIJ(2)*RIJ(3) - (DVADR*INVR*RIJ(3) + D2VDA2*DADR(3) &
                               + DVADB*DBDR(3))*DADR(2) - VA*D2ADZY - (DVBDR*INVR*RIJ(3) &
                               + DVBDA*DADR(3) + D2VDB2*DBDR(3))*DBDR(2) - VB*D2BDZY
               HESS(J3,J4-1)   = HESS(J3-1,J4)
!     xi,zj 
               HESS(J3-2,J4)   = -(DVRDR*INVR*RIJ(1) + DVRDA*DADR(1) + DVRDB*DBDR(1))*INVR*RIJ(3) &
                               + VR*INVR*R2*RIJ(3)*RIJ(1) - (DVADR*INVR*RIJ(1) + D2VDA2*DADR(1) &
                               + DVADB*DBDR(1))*DADR(3) - VA*D2ADXZ - (DVBDR*INVR*RIJ(1) &
                               + DVBDA*DADR(1) + D2VDB2*DBDR(1))*DBDR(3) - VB*D2BDXZ  
               HESS(J3,J4-2)   = HESS(J3-2,J4)

               FCT22 = DVRDB*DBDPJ1 + DVRDG*DGDPJ1 - ((DVADB*DBDPJ1 + DVADG*DGDPJ1)*ALP &
                     + (D2VDB2*DBDPJ1 + DVBDG*DGDPJ1)*BET + VB*DBDPJ1)*INVR 
               FCT23 = (DVADB*DBDPJ1 + DVADG*DGDPJ1)*INVR
               FCT24 = (D2VDB2*DBDPJ1 + DVBDG*DGDPJ1)*INVR
!     xi,pj1
               DUMMY = FCT22*NRSS(1) + FCT23*EI(1) + FCT24*EJ(1) + VB*INVR*DEDPJ1(1)
               HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
               HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
               DUMMY = FCT22*NRSS(2) + FCT23*EI(2) + FCT24*EJ(2) + VB*INVR*DEDPJ1(2)
               HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
               HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
               DUMMY = FCT22*NRSS(3) + FCT23*EI(3) + FCT24*EJ(3) + VB*INVR*DEDPJ1(3)
               HESS(J3,J6-2)   = HESS(J3,J6-2)   + DUMMY
               HESS(J6-2,J3)   = HESS(J6-2,J3)   + DUMMY

               FCT22 = DVRDB*DBDPJ2 + DVRDG*DGDPJ2 - ((DVADB*DBDPJ2 + DVADG*DGDPJ2)*ALP &
                     + (D2VDB2*DBDPJ2 + DVBDG*DGDPJ2)*BET + VB*DBDPJ2)*INVR 
               FCT23 = (DVADB*DBDPJ2 + DVADG*DGDPJ2)*INVR
               FCT24 = (D2VDB2*DBDPJ2 + DVBDG*DGDPJ2)*INVR
!     xi,pj2
               DUMMY = FCT22*NRSS(1) + FCT23*EI(1) + FCT24*EJ(1) + VB*INVR*DEDPJ2(1)
               HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
               HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
               DUMMY = FCT22*NRSS(2) + FCT23*EI(2) + FCT24*EJ(2) + VB*INVR*DEDPJ2(2)
               HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
               HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
               DUMMY = FCT22*NRSS(3) + FCT23*EI(3) + FCT24*EJ(3) + VB*INVR*DEDPJ2(3)
               HESS(J3,J6-1)   = HESS(J3,J6-1)   + DUMMY
               HESS(J6-1,J3)   = HESS(J6-1,J3)   + DUMMY

               FCT22 = DVRDB*DBDPJ3 + DVRDG*DGDPJ3 - ((DVADB*DBDPJ3 + DVADG*DGDPJ3)*ALP &
                     + (D2VDB2*DBDPJ3 + DVBDG*DGDPJ3)*BET + VB*DBDPJ3)*INVR 
               FCT23 = (DVADB*DBDPJ3 + DVADG*DGDPJ3)*INVR
               FCT24 = (D2VDB2*DBDPJ3 + DVBDG*DGDPJ3)*INVR
!     xi,pj3
               DUMMY = FCT22*NRSS(1) + FCT23*EI(1) + FCT24*EJ(1) + VB*INVR*DEDPJ3(1)
               HESS(J3-2,J6)   = HESS(J3-2,J6)   + DUMMY
               HESS(J6,J3-2)   = HESS(J6,J3-2)   + DUMMY
!     yi,pj3 
               DUMMY = FCT22*NRSS(2) + FCT23*EI(2) + FCT24*EJ(2) + VB*INVR*DEDPJ3(2)
               HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
               HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
               DUMMY = FCT22*NRSS(3) + FCT23*EI(3) + FCT24*EJ(3) + VB*INVR*DEDPJ3(3)
               HESS(J3,J6)     = HESS(J3,J6)     + DUMMY
               HESS(J6,J3)     = HESS(J6,J3)   + DUMMY
!     pi1,pj2
               HESS(J5-2,J6-1) = (DVADB*DBDPJ2 + DVADG*DGDPJ2)*DADPI1 + (DVGDB*DBDPJ2 &
                               + D2VDG2*DGDPJ2)*DGDPI1 + VG*DOT_PRODUCT(DEDPI1,DEDPJ2) 
               HESS(J6-1,J5-2) = HESS(J5-2,J6-1)
!     pi2,pj3
               HESS(J5-1,J6)   = (DVADB*DBDPJ3 + DVADG*DGDPJ3)*DADPI2 + (DVGDB*DBDPJ3 &
                               + D2VDG2*DGDPJ3)*DGDPI2 + VG*DOT_PRODUCT(DEDPI2,DEDPJ3) 
               HESS(J6,J5-1)   = HESS(J5-1,J6)
!     pi3,pj1
               HESS(J5,J6-2)   = (DVADB*DBDPJ1 + DVADG*DGDPJ1)*DADPI3 + (DVGDB*DBDPJ1 &
                               + D2VDG2*DGDPJ1)*DGDPI3 + VG*DOT_PRODUCT(DEDPI3,DEDPJ1) 
               HESS(J6-2,J5)   = HESS(J5,J6-2)
               
            ENDDO

         ENDDO

!     CALCULATE THE GAY-BERNE POTENTIAL: EACH PAIR HAS BEEN CONSIDERED TWICE

         ENERGY = 2.D0 * ENERGY
         G(:)   = 0.5D0*G(:)

!         DO I = 1, 3*NATOMS

!            DO J = 1, 3*NATOMS
  
!               WRITE(*,*) I, J, HESS(I,J) 

!            ENDDO

!         ENDDO

      ENDIF

      END SUBROUTINE GB		 
