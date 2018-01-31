! Major commentary by Scott Olesen swo24 (July 2011)
! Last updates are on (at earliest) 11 July 2011
! 
! Subroutines in this file:
!     MULTISITEPY - called by MSPYG; I don't think it works
!     DEFMSPY - defines the structure of the multisite molecule
!     SITEBF - 
!     DEFINEPYMULTISITES - reads from pysites.xyz to get PST, etc
!     DEFINELJMULTISITES - reads from ljsites.xyz to get LJGSITECOORDS, etc
!     PYSITEORIENTATIONS -
!     MULTISITEPY2 - called by MULTISITEPY; the main subroutine
!     AAtoSites - 
!     TAKESTEPMULTISITEPY - 

!-------------------------------------------------------!
! This is the older version of the multisite subroutine !
!-------------------------------------------------------!

      SUBROUTINE MULTISITEPY (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS
      USE KEY, ONLY :  NPYSITE, PYA1, PYA2, PYSIGNOT, PYEPSNOT, RADIFT, FROZEN, EFIELDT 
      USE PYMODULE, ONLY : pi, twopi
      IMPLICIT NONE

      INTEGER          :: I, I1, J, J1, J2, J3, J4, J5, J6, K1, K2, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), energy_contrib2, grad_contrib2
      DOUBLE PRECISION :: PST(NATOMS,3), EZRI1(3,3), EZRI2(3,3), EZRJ1(3,3), EZRJ2(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: I3(3,3), AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: RSTI(3), RSTJ(3), FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DV1DF1, DV2DF2, DV1DR, DV2DR
      DOUBLE PRECISION :: DRDPI1, DRDPI2, DRDPI3, DRDPJ1, DRDPJ2, DRDPJ3 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      double precision :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3), angle

      LOGICAL          :: GTEST

      I3(:,:)      = 0.D0
       NPYSITE = 100

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      DO K1 = 1, 3

         I3(K1,K1)      = 1.D0

      ENDDO


      CALL DEFMSPY(PST)

      IF (GTEST) THEN

         ENERGY = 0.D0
         G(:)   = 0.D0


         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

            angle=sqrt(dot_product(P,P))
            if(angle>twopi) then
! normalise angle-axis coordinates
                X(J5-2:J5)=X(J5-2:J5)/angle
                do
                  angle=angle-twopi
                  if(angle<2*pi) exit
                end do
! multiply with new angle
                X(J5-2:J5)=X(J5-2:J5)*angle
            end if


!     ROTATION MATRIX

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDRVT(P, RMJ, DRMJ1, DRMJ2, DRMJ3, GTEST)

               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))
                  
                  CALL SITEBF (I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF
               
                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))

                     CALL SITEBF (J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN
   
                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF
                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)

                     FCNT1   = - FMIN
 
                     SRTFI1  = 1.D0 / DSQRT(FCNT1)
                     APB     = LAMDAC * AE1 + (1.D0 - LAMDAC) * BE1

                     CALL MTRXIN (APB, APBINV)

                     ARIBRJ =  LAMDAC * MATMUL(AE1,(RI+RSTI)) + (1.D0 - LAMDAC) * MATMUL(BE1,(RJ+RSTJ))
                     XC     =  MATMUL(APBINV, ARIBRJ)
                     XCMRI  = XC - RI - RSTI
                     XCMRJ  = XC - RJ - RSTJ
                     DF1DR  = - 2.D0 * LAMDAC * MATMUL(AE1,XCMRI)

                     D1ABEZ = MATMUL(DRMI1,EZRI1(:,:))
                     D1ABE  = MATMUL(D1ABEZ,TRANSPOSE(RMI)) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                     D2ABEZ = MATMUL(DRMI2,EZRI1(:,:))
                     D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                     D3ABEZ = MATMUL(DRMI3,EZRI1(:,:))
                     D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                     DF1PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI1,PST(I,:))))) 
                     DF1PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI2,PST(I,:)))))
                     DF1PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI3,PST(I,:)))))

                     D1ABEZ = MATMUL(DRMJ1,EZRJ1(:,:))
                     D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                     D2ABEZ = MATMUL(DRMJ2,EZRJ1(:,:))
                     D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                     D3ABEZ = MATMUL(DRMJ3,EZRJ1(:,:))
                     D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
                     DF1PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ1,PST(J,:)))))
                     DF1PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ2,PST(J,:)))))
                     DF1PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ3,PST(J,:))))) 

                     RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
                     RHO1SQ = RHO1*RHO1
                     RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                     RHO112 = RHO16 * RHO16

                     FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
                     DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
                     DV1DF1 = -2.D0*RHO112*RHO1*FCTR1
                     DV1DR  = -2.D0*RHO112*RHO1*(1.D0-SRTFI1)/PYSIGNOT

                     IF (RADIFT) THEN

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC, FMIN)

                        FCNT2   = - FMIN
                        SRTFI2  = 1.D0 / DSQRT(FCNT2)
                        APB     = LAMDAC * AE2 + (1.D0 - LAMDAC) * BE2

                        CALL MTRXIN (APB, APBINV)

                        ARIBRJ = LAMDAC * MATMUL(AE2,(RI+RSTI)) + (1.D0 - LAMDAC) * MATMUL(BE2,(RJ+RSTJ))
                        XC     = MATMUL(APBINV, ARIBRJ)
                        XCMRI  = XC - RI - RSTI
                        XCMRJ  = XC - RJ - RSTJ
                        DF2DR  = - 2.D0 * LAMDAC * MATMUL(AE2,XCMRI)

                        RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                        FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                        DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                        DV2DF2 = RHO26*RHO2*FCTR2
                        DV2DR  = RHO26*RHO2*(1.D0-SRTFI2)/PYSIGNOT

                        D1ABEZ = MATMUL(DRMI1,EZRI2(:,:))
                        D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                        D2ABEZ = MATMUL(DRMI2,EZRI2(:,:))
                        D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                        D3ABEZ = MATMUL(DRMI3,EZRI2(:,:))
                        D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                        DF2PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI1,PST(I,:)))))
                        DF2PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI2,PST(I,:)))))
                        DF2PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI3,PST(I,:)))))
                        
                        D1ABEZ = MATMUL(DRMJ1,EZRJ2(:,:))
                        D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                        D2ABEZ = MATMUL(DRMJ2,EZRJ2(:,:))
                        D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                        D3ABEZ = MATMUL(DRMJ3,EZRJ2(:,:))
                        D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                        DF2PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ1,PST(J,:)))))
                        DF2PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ2,PST(J,:)))))
                        DF2PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ3,PST(J,:)))))

                     ELSE

                        RHO2   = RHO1
                        RHO26  = RHO16
                        DG2DR  = DG1DR
                        DV2DF2 = RHO26*RHO2*FCTR1
                        DV2DR  = RHO26*RHO2*(1.D0-SRTFI1)/PYSIGNOT
                        DF2PI1 = DF1PI1
                        DF2PI2 = DF1PI2
                        DF2PI3 = DF1PI3
                        DF2PJ1 = DF1PJ1
                        DF2PJ2 = DF1PJ2
                        DF2PJ3 = DF1PJ3

                     ENDIF
             
!     CALCULATE PY POTENTIAL ENERGY

!     CALCULATE GRADIENT

                     FIJ = 2.D0*RHO112*RHO1*DG1DR - RHO26*RHO2*DG2DR
                    
                     DRDPI1 = DOT_PRODUCT(NR,MATMUL(DRMI1,PST(I,:)))
                     DRDPI2 = DOT_PRODUCT(NR,MATMUL(DRMI2,PST(I,:)))
                     DRDPI3 = DOT_PRODUCT(NR,MATMUL(DRMI3,PST(I,:)))

                     TIJ(1) = DV1DF1*DF1PI1 + DV1DR*DRDPI1 + DV2DF2*DF2PI1 + DV2DR*DRDPI1
                     TIJ(2) = DV1DF1*DF1PI2 + DV1DR*DRDPI2 + DV2DF2*DF2PI2 + DV2DR*DRDPI2
                     TIJ(3) = DV1DF1*DF1PI3 + DV1DR*DRDPI3 + DV2DF2*DF2PI3 + DV2DR*DRDPI3

                     DRDPJ1 = DOT_PRODUCT(NR,MATMUL(DRMJ1,-PST(J,:)))
                     DRDPJ2 = DOT_PRODUCT(NR,MATMUL(DRMJ2,-PST(J,:)))
                     DRDPJ3 = DOT_PRODUCT(NR,MATMUL(DRMJ3,-PST(J,:)))

                     TJI(1) = DV1DF1*DF1PJ1 + DV1DR*DRDPJ1 + DV2DF2*DF2PJ1 + DV2DR*DRDPJ1
                     TJI(2) = DV1DF1*DF1PJ2 + DV1DR*DRDPJ2 + DV2DF2*DF2PJ2 + DV2DR*DRDPJ2
                     TJI(3) = DV1DF1*DF1PJ3 + DV1DR*DRDPJ3 + DV2DF2*DF2PJ3 + DV2DR*DRDPJ3

                     G(J3-2:J3) = G(J3-2:J3) - FIJ
                     G(J4-2:J4) = G(J4-2:J4) + FIJ
                     G(J5-2:J5) = G(J5-2:J5) + TIJ
                     G(J6-2:J6) = G(J6-2:J6) + TJI

                  ENDDO

               ENDDO

!     END INNER LOOP OVER PARTICLES

            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO

         ENERGY = 4.D0*PYEPSNOT*ENERGY
         G      = 24.D0*PYEPSNOT*G

    do i1 = 1, natoms/2
!sf344> gravity field with constant gradient, and an r^-12 repulsion around z=0
      if(efieldt) then
        energy_contrib2=0.0
        grad_contrib2=0.0 ! this is a scalar, we have only contribution along the z direction
        call gravity(x(3*i1),energy_contrib2,grad_contrib2)
        energy = energy + energy_contrib2
        g(3*i1)=g(3*i1)+grad_contrib2
      end if
    end do

      ELSE

         ENERGY = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            THETA  = DSQRT(DOT_PRODUCT(P,P))

            IF (THETA == 0.D0) THEN

               RMI = I3

            ELSE
            
               THETA2 = THETA * THETA
               CT      = COS(THETA)
               ST      = SIN(THETA)
               E(:,:)  = 0.D0
               E(1,2)  = -P(3)
               E(1,3)  =  P(2)
               E(2,3)  = -P(1)
               E(2,1)  = -E(1,2)
               E(3,1)  = -E(1,3)
               E(3,2)  = -E(2,3)
               E       = E/THETA

               RMI     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4)
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               THETA  = DSQRT(DOT_PRODUCT(P,P))
              
               IF (THETA == 0.D0) THEN

                  RMJ = I3

               ELSE

                  THETA2 = THETA * THETA
                  CT      = COS(THETA)
                  ST      = SIN(THETA)
                  E(:,:)  = 0.D0
                  E(1,2)  = -P(3)
                  E(1,3)  =  P(2)
                  E(2,3)  = -P(1)
                  E(2,1)  = -E(1,2)
                  E(3,1)  = -E(1,3)
                  E(3,2)  = -E(2,3)
                  E       = E/THETA

                  RMJ     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

               ENDIF
 
               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))

                  CALL SITEBF(I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF

                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))
 
                     CALL SITEBF (J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN

                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF

                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)

                     FCNT1   = - FMIN
                     SRTFI1  = 1.D0 / DSQRT(FCNT1)

                     RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
                     RHO1SQ = RHO1*RHO1
                     RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                     RHO112 = RHO16 * RHO16

                     IF (RADIFT) THEN

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC, FMIN)

                        FCNT2   = - FMIN
                        SRTFI2  = 1.D0 / DSQRT(FCNT2)

                        RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ

                     ELSE

                        RHO2   = RHO1
                        RHO26  = RHO16

                     ENDIF

                     ENERGY = ENERGY + RHO112 - RHO26

                  ENDDO

               ENDDO

!     CALCULATE PY POTENTIAL ENERGY

            ENDDO

         ENDDO

         ENERGY = 4.D0*PYEPSNOT*ENERGY
    do i1 = 1, natoms/2
!sf344> gravity field with constant gradient, and an r^-12 repulsion around z=0
      if(efieldt) then
        energy_contrib2=0.0
        grad_contrib2=0.0 ! this is a scalar, we have only contribution along the z direction
        call gravity(x(3*i1),energy_contrib2,grad_contrib2)
        energy = energy + energy_contrib2
        g(3*i1)=g(3*i1)+grad_contrib2
      end if
    end do


      ENDIF
      END SUBROUTINE MULTISITEPY 


!---------------------------------!
! DEFMSPY - Called by MULTISITEPY !
!---------------------------------!
! This subroutine will define where the multiple PY sites are supposed to go.
! Because this is all hardcoded and is for MULTISITEPY, I think this is obsolete.
      SUBROUTINE DEFMSPY(PST)

!        PYA1 ~ components a_1k as of eqns between 25 & 26
!        NPYSITE ~ number of PY sites in each molecule
      USE KEY, ONLY: PYA1, NPYSITE  

      IMPLICIT NONE

!        LENGTH ~ length scale used in dimensions of the molecule
!        PST ~ ? Position of SiTes ?. PST(1,i) ~ i-th coordinate of 1st site
      DOUBLE PRECISION :: LENGTH, PST(NPYSITE,3)

!      LENGTH   = 2.D0*PYA1(1)

!        Looks like, for the time being, the sites are right on top of each other.
      LENGTH = 0.D0
      PST(1,1) = - 0.25D0*LENGTH
      PST(1,2) = 0.D0
      PST(1,3) = 0.D0

      PST(2,1) = 0.25D0*LENGTH
      PST(2,2) = 0.D0
      PST(2,3) = 0.D0

      END SUBROUTINE DEFMSPY

!     --------------------------------------------------------------------------

      SUBROUTINE SITEBF (K, EZR1, EZR2)

      USE KEY, ONLY: PYA1, PYA2
    
      IMPLICIT NONE
      
      INTEGER          :: K 
      DOUBLE PRECISION :: EZR1(3,3), EZR2(3,3)

      EZR1(:,:) = 0.D0
      EZR2(:,:) = 0.D0

      IF (K == 1) THEN

         EZR1(1,1) = 1.D0/(PYA1(1)*PYA1(1))
         EZR1(2,2) = 1.D0/(PYA1(2)*PYA1(2))
         EZR1(3,3) = 1.D0/(PYA1(3)*PYA1(3))

         EZR2(1,1) = 1.D0/(PYA2(1)*PYA2(1))
         EZR2(2,2) = 1.D0/(PYA2(2)*PYA2(2))
         EZR2(3,3) = 1.D0/(PYA2(3)*PYA2(3))

      ELSE

         EZR1(1,1) = 1.D0/(PYA1(2)*PYA1(2))
         EZR1(2,2) = 1.D0/(PYA1(1)*PYA1(1))
         EZR1(3,3) = 1.D0/(PYA1(3)*PYA1(3))

         EZR2(1,1) = 1.D0/(PYA2(2)*PYA2(2))
         EZR2(2,2) = 1.D0/(PYA2(1)*PYA2(1))
         EZR2(3,3) = 1.D0/(PYA2(3)*PYA2(3))

      ENDIF

      END SUBROUTINE SITEBF

!----------------------------------------!
! BEGIN SUBROUTINE DEFINEPYMULTISITES #P !
!----------------------------------------!
! Called by io1.f when using keyword MULTISITEPY. This subroutine reads off
! the site positions and orientations from a file pysites.xyz. It checks which
! sites have equal repulsive and attractive parts and notes which of the semiaxes
! is longest.
!
! Syntax of pysites.xyz files:
! [number of sites]
! [blank line]
! site [x coords] [y] [z] shapes [repsulvie coefficient] [attr coeff] [repsulvie x semiaxis] [y] [z] ...
!   [attr x] [y] [z] orient [x component of p vector] [y] [z]

      SUBROUTINE DEFINEPYMULTISITES

      ! NPYSITE = number of sites per molecule
      ! NATOMS = number of lines in coords = twice the number of rigid bodies ("molecules")
      ! MYUNIT = file unit for GMIN_out
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
      USE KEY, ONLY: NPYSITE, NTSITES

      ! PST[site][xyz] = position of ellipsoid center relative to molecule position (ie, in body frame)
      ! OST[site][xyz] = p vector describing site rotation in body frame
      ! ELLST1[site][xyz] = semiaxes of repulsive ellipsoid in site frame
      ! ELLST2[site][xyz] = semiaxes of attractive ellipsoid in site frame
      ! RADIFTARRAY[site] = logical: if rep & attr shape matrices are different for the site
      ! pyrepp[site] = coefficent that adjusts contributions to repulsive potential energy for that site
      ! pyattrp[site] = coefficent that adjusts contributions to attractive potential energy for that site
      ! LONGRSAX[site] = LONGest Repulsive SemiAXis [PBC]
      ! LONGASAX[site] = longest attractive semiaxis [PBC]
      ! SMRBF[site][][] = Shape Matrix Repulsive in Body Frame = ELLST1 rotated according to OST
      ! SMABF[site][][] = shape matrix attractive in body frame = ELLST2 rotated according to OST
      ! SMRDET[site] = DETerminant of the Shape Matrix for the Repulsive part
      ! SMADET[site] = DETerminant of the Shape Matrix for the Attractive part
      ! SITECOORDS[site][xyz] = ?, maybe holder for final output coordinates
      ! ELLMAT[site][][] = ? used by AAtoSites
      USE PYMODULE, ONLY : PST, OST, ELLST1, ELLST2, RADIFTARRAY, & 
            & pyattrp, pyrepp, LONGRSAX, LONGASAX, SMRBF, SMABF, &
            & SMRDET, SMADET, ELLMAT, SITECOORDS

      IMPLICIT NONE

      ! J1, J2 = counters
      ! LABEL, etc = placeholders for labelling strings in pysites.xyz ("site", "shapes", "orient")
      INTEGER          :: J1, J2, J3
      CHARACTER(LEN=5) :: LABEL
      CHARACTER(LEN=7) :: DUMMYLABEL1
      CHARACTER(LEN=11):: DUMMYLABEL2

      ! RMAT[][] = rotation matrix computed from angle-axis OST
      ! DRMAT1, etc = placeholders since RMDRVT is called without computing derivatives
      DOUBLE PRECISION :: RMAT(3,3),DRMAT1(3,3),DRMAT2(3,3),DRMAT3(3,3),vecsbf(3)
      DOUBLE PRECISION, ALLOCATABLE :: XEND(:,:,:)

      ! Open the file. Read 1st line as NPYSITE. Skip a line.
      OPEN(UNIT=299,FILE="pysites.xyz",STATUS="old")
      READ(299,*) NPYSITE 
      READ(299,*) 
      NRBSITES=NPYSITE*2
      NTSITES = NATOMS*NRBSITES/2
      ! Allocate arrays that use NPYSITE as a bound. Assume shapes at
      ! a site are different until proven equal.
      IF(.NOT.ALLOCATED(PST)) ALLOCATE(PST(NPYSITE,3))
      IF(.NOT.ALLOCATED(OST)) ALLOCATE(OST(NPYSITE,3))
      IF(.NOT.ALLOCATED(ELLST1)) ALLOCATE(ELLST1(NPYSITE,3))
      IF(.NOT.ALLOCATED(ELLST2)) ALLOCATE(ELLST2(NPYSITE,3))
      IF(.NOT.ALLOCATED(RADIFTARRAY)) ALLOCATE(RADIFTARRAY(NPYSITE))
      IF(.NOT.ALLOCATED(pyattrp)) ALLOCATE(pyattrp(NPYSITE))
      IF(.NOT.ALLOCATED(pyrepp)) ALLOCATE(pyrepp(NPYSITE))
      IF(.NOT.ALLOCATED(LONGRSAX)) ALLOCATE(LONGRSAX(NPYSITE))
      IF(.NOT.ALLOCATED(LONGASAX)) ALLOCATE(LONGASAX(NPYSITE))
      IF(.NOT.ALLOCATED(SMRBF)) ALLOCATE(SMRBF(NPYSITE,3,3))
      IF(.NOT.ALLOCATED(SMABF)) ALLOCATE(SMABF(NPYSITE,3,3))
      IF(.NOT.ALLOCATED(SMRDET)) ALLOCATE(SMRDET(NPYSITE))
      IF(.NOT.ALLOCATED(SMADET)) ALLOCATE(SMADET(NPYSITE))
      IF(.NOT.ALLOCATED(ELLMAT)) ALLOCATE(ELLMAT(NPYSITE,3,3))
      IF(.NOT.ALLOCATED(SITECOORDS)) ALLOCATE(SITECOORDS(NPYSITE,3))
      IF(.NOT.ALLOCATED(XEND)) ALLOCATE(XEND(2,NPYSITE,3))

      IF(.NOT.ALLOCATED(RBSITE)) ALLOCATE(RBSITE(NRBSITES,3))
      RADIFTARRAY(:) = .TRUE.
          vecsbf(1)=0.0D0
          vecsbf(2)=0.0D0
          vecsbf(3)=1.0D0
!    IF(LJSITECOORDST) THEN
!        vecsbf(:)=LJSITECOORDS(:)
!        WRITE(*,'(A,3F8.3)') ' repulsive LJ site coordinates will be ', LJSITECOORDS(:)
!    END IF

      ! Begin loop lines in input. Each line corresponds to a site in the molecule.
      DO J1=1,NPYSITE
            ! Read line as per syntax above.
            READ(299,*) LABEL, PST(J1,1), PST(J1,2), PST(J1,3), &
             & DUMMYLABEL1, pyrepp(J1), pyattrp(J1), &
             & ELLST1(J1,1),ELLST1(J1,2), ELLST1(J1,3), & 
             & ELLST2(J1,1),ELLST2(J1,2), ELLST2(J1,3), &
             & DUMMYLABEL2, OST(J1,1), OST(J1,2), OST(J1,3)
      END DO ! loop over all sites

      ! Close pysites.xyz.
      CLOSE(299)
      DO J1=1,NPYSITE
            ! Decide which sites have rep & attr parts equal.
            IF(      ELLST2(J1,1)==ELLST1(J1,1) &
             & .AND. ELLST2(J1,2)==ELLST1(J1,2) & 
             & .AND. ELLST2(J1,3)==ELLST1(J1,3)) THEN
                  WRITE(*,*) 'attr and rep semiaxes are same for &
                   &site ', J1 
                  RADIFTARRAY(J1) = .FALSE.
            END IF

            ! Note which semiaxis is the longest for each ellipsoid @ site J1
            LONGRSAX(J1)=MAXVAL(ELLST1(J1,:))
            LONGASAX(J1)=MAXVAL(ELLST2(J1,:))

            ! Calculate the shape matrices for the ellipsoids in the body frame:
            ! Initialize matrices. Compute rotation matrix RMAT from angle-axis coords OST.
            ! Calculate shape matrices with no rotation as per eq 23.
            SMRBF(J1,:,:) = 0.D0
            SMABF(J1,:,:) = 0.D0
            DO J3=1, 3
                  SMRBF(J1,J3,J3) = 1.D0/(ELLST1(J1,J3)**2)
                  SMABF(J1,J3,J3) = 1.D0/(ELLST2(J1,J3)**2)
            ENDDO
            
            ! Compute the inverse determinants of the shape matrices. They do not change as the
            ! matrices are rotated.
            SMRDET(J1) = 1.D0 / (SMRBF(J1,1,1)*SMRBF(J1,2,2)*SMRBF(J1,3,3))
            SMADET(J1) = 1.D0 / (SMABF(J1,1,1)*SMABF(J1,2,2)*SMABF(J1,3,3))

            ! Rewrite the shape matrices as they appear in the body frame
            CALL RMDRVT(OST(J1,:), RMAT, DRMAT1, DRMAT2, DRMAT3, .FALSE.)
            SMRBF(J1,:,:) = MATMUL(RMAT,(MATMUL(SMRBF(J1,:,:),TRANSPOSE(RMAT))))
            SMABF(J1,:,:) = MATMUL(RMAT,(MATMUL(SMABF(J1,:,:),TRANSPOSE(RMAT))))

! generate coordinates for the end of the largest ellipsoid axis
        ! sf344> reset RBSITES to the centre of mass for each ellipsoid. Forget about RBSITE for now.
! sf344> try to define sites for axially symmetric particles...
            xend(1,J1,:)=PST(J1,:)+MATMUL(RMAT,vecsbf) !+ELLST1(J1,2)*MATMUL(RMAT,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid
            xend(2,J1,:)=PST(J1,:)-MATMUL(RMAT,vecsbf) !-ELLST1(J1,2)*MATMUL(RMAT,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid

      END DO ! loop over sites

      J1=0
      J2=0
      DO 
        J1=J1+2
        J2=J2+1
         RBSITE(J1-1,:)=XEND(1,J2,:)
         RBSITE(J1,:)=XEND(2,J2,:)
        IF(J1==NRBSITES) EXIT
      END DO
        WRITE(*,*) 'PY sites defined as: '
        
      DO J1=1,NRBSITES
        WRITE(*,*) RBSITE(J1,:)
      END DO
      END SUBROUTINE DEFINEPYMULTISITES


!-------------------------------------!
! BEGIN SUBROUTINE DEFINELJMULTISITES !
!-------------------------------------!
! Called by io1.f. This subroutine reads off
! the positions of LJ sites from a file ljsites.xyz
!
! Syntax for ljsites.xyz:
! [number of sites]
! [blank line]
! coord [x] [y] [z] coeffs [rep coeff] [attr coeff]

      SUBROUTINE DEFINELJMULTISITES

      ! See DEFINEPYMULTISITES for description of these variables
      USE COMMONS, ONLY: NATOMS

      ! LJGSITECOORDS[site][xyz] = position of site in body frame
      ! ljattrp[site] = coefficient to adjust attractive contributions to energy
      ! ljrepp[site] = as above for repulsive
      USE PYMODULE, ONLY : LJGSITECOORDS, NLJSITE, ljattrp, ljrepp

      IMPLICIT NONE

      ! J1 = counter
      ! LABEL, DUMMY = comments in ljsites.xyz
      INTEGER          :: J1
      CHARACTER(LEN=5) :: LABEL
      CHARACTER(LEN=7) :: DUMMYLABEL1

      ! Open the file. First line should read as number of LJ sites in 
      ! molecule. Then read a blank line. 
      OPEN(UNIT=300,FILE="ljsites.xyz",STATUS="old")
      READ(300,*) NLJSITE 
      READ(300,*) 

      ! Allocate arrays
      IF(.NOT.ALLOCATED(LJGSITECOORDS)) ALLOCATE(LJGSITECOORDS(NLJSITE,3))
      IF(.NOT.ALLOCATED(ljattrp)) ALLOCATE(ljattrp(NLJSITE))
      IF(.NOT.ALLOCATED(ljrepp)) ALLOCATE(ljrepp(NLJSITE))

      ! Begin loop over the LJ sites in the model molecule
      DO J1=1,NLJSITE

        ! Read the input from ljsites.xyz as per above syntax
        READ(300,*) LABEL, LJGSITECOORDS(J1,1), LJGSITECOORDS(J1,2), &
             & LJGSITECOORDS(J1,3), DUMMYLABEL1, ljrepp(J1), ljattrp(J1)
        WRITE(*,*) &
           & 'defineljmultisites> rep and attr params for body ', J1, &
           & ' are: ' , ljrepp(J1), ' & ', ljattrp(J1)

      END DO ! loop over all sites

      ! Close ljsites.xyz
      CLOSE(300)

      END SUBROUTINE DEFINELJMULTISITES


!-------------------------------------!
! BEGIN SUBROUTINE PYSITEORIENTATIONS !
!-------------------------------------!
! Obsolete. Previously used by MULTISITEPY2. Code was merged into
! DEFINEPYMULTISITES.

      SUBROUTINE PYSITEORIENTATIONS (K, EZR1R, EZR2R)

      USE KEY, ONLY: NPYSITE
      USE PYMODULE, ONLY : PST, OST, ELLST1, ELLST2
   
      IMPLICIT NONE
      
      INTEGER          :: K 
      DOUBLE PRECISION :: EZR1(3,3), EZR2(3,3), EZR1R(3,3), EZR2R(3,3)
      DOUBLE PRECISION :: RMAT(3,3),DRMAT1(3,3),DRMAT2(3,3),DRMAT3(3,3)

      ! Initialize arrays
      EZR1(:,:) = 0.D0
      EZR2(:,:) = 0.D0

      ! Compute rotation matrix RMAT from angle-axis coords OST specified
      ! in pysites.xyz
      CALL RMDRVT(OST(K,:), RMAT, DRMAT1, DRMAT2, DRMAT3, .FALSE.)

      ! Calculate shape matrices with no rotation as per eq 23
      EZR1(1,1) = 1.D0/(ELLST1(K,1)*ELLST1(K,1))
      EZR1(2,2) = 1.D0/(ELLST1(K,2)*ELLST1(K,2))
      EZR1(3,3) = 1.D0/(ELLST1(K,3)*ELLST1(K,3))

      EZR2(1,1) = 1.D0/(ELLST2(K,1)*ELLST2(K,1))
      EZR2(2,2) = 1.D0/(ELLST2(K,2)*ELLST2(K,2))
      EZR2(3,3) = 1.D0/(ELLST2(K,3)*ELLST2(K,3))

      ! Rewrite matrices for the body frame, ie rotated by angle-axes coords
      ! specified in pysites.xyz
      EZR1R = MATMUL(RMAT,(MATMUL(EZR1(:,:),(TRANSPOSE(RMAT)))))
      EZR2R = MATMUL(RMAT,(MATMUL(EZR2(:,:),(TRANSPOSE(RMAT)))))

      END SUBROUTINE PYSITEORIENTATIONS




!------------------------------------!
!------------------------------------!
! BEGIN MAIN MULTISITE SUBROUTINE #A !
!------------------------------------!
!------------------------------------!
!
! Basic overview of the subroutine follows. Unless otherwise specified, equation
! references are to Chakrabarti & Wales, Phys Chem Chem Phys 11 1970-6 (2009).
! References labeled 'PY' refer to Paramonov & Yaliraki, JCP 123 194111ff (2005).
! Cutoff scheme is borrowed from Stoddard & Ford, PRA 8(3) 1504ff (1973).
! In derivatives, I use 'd' to mean partial, not full, derivative. It's easier to type.
! Square brackets refer to the three optional conditions for running the
! subroutine: cutoff [CUTOFF], periodic boundary conditions [PBC], and the use
! of general LJ sites [LJGSITE]. The 'G' is for general(issimo).
!
! Initialize variables, etc (#B)
! Outer molecule loop (#C)
! Inner molecule loop (#D)
!     Loop over sites in outer molecule (#E)
!     Loop over sites in inner molecule (#F)
!           Compute repulsive ECF (#G)
!           Calculate rep variables (#H)
!           Compute attractive ECF (#I)
!           Compute attr vars (#J)
!           Add repulsive contributions to energy and grad (#K)
!           Add attr contribs (#L)
!     [LJGSITE] (#M)
! Finalize (#N) 

      SUBROUTINE MULTISITEPY2 (X, G, ENERGY, GTEST)

      ! NATOMS, NPYSITE, MYUNIT = cf DEFINEPYMULTISITES above
      ! PYSIGNOT = sigma_0 (eq 25)
      ! PYEPSNOT = epsilon_0 (eq 24)
      ! FROZEN[body] = logical if body is fixed through the run
      ! LJGSITET  = logical if LJ sites are specified in data/ljsites.xyz
      ! LJGSITEEPS, LJGSITESIGMA = cf DEFINELJMULTISITES above
      ! [CUTOFF]
      ! PARAMONOVPBCX = logical if there is a PBC in the x direction
      ! BOXLX = size of PBC box in x direction in absolute units
      ! PARAMONOVCUTOFF = logical if cutoff is used
      ! PCUTOFF = cutoff distance in absolute units
      USE COMMONS, ONLY: NATOMS 
      USE KEY, ONLY: NPYSITE, PYSIGNOT, PYEPSNOT, FROZEN, &
            & LJGSITET, LJGSITEEPS, LJGSITESIGMA, &
            & PARAMONOVPBCX, PARAMONOVPBCY, PARAMONOVPBCZ, &
            & PARAMONOVCUTOFF, PCUTOFF, EFIELDT, BULKT
      USE PYMODULE, ONLY : BOXLX, BOXLY, BOXLZ

      ! PST, ..., SMABF = cf DEFINEPYMULTISITES above
      ! NLJSITE, ..., ljattrp = cf DEFINELJMULTISITES above
      ! LONGRSAX, LONGASAX = cf DEFINEPYMULTISITES above
      USE PYMODULE, ONLY : PST, OST, RADIFTARRAY, pi, twopi, &
            & pyrepp, pyattrp, LONGRSAX, LONGASAX, SMRBF, SMABF, &
            & NLJSITE, LJGSITECOORDS, ljrepp, ljattrp, &
            & SMRDET, SMADET
 
      IMPLICIT NONE

      ! J1 = index for outer loop over all molecules (cf sum over I in eq 9)
      ! J2 = index for inner loop over molecules with J2 > J1 (cf sum over J<I in eq 9)
      ! J3 = index of z coord of J1-th molecule
      ! J4 = index of z coord of J2-th molecule
      ! J5 = index of z component of p vector for J1-th molecule
      ! J6 = index of z component of p vector for J2-th molecule
      ! REALNATOMS = number of molecules
      ! OFFSET = indexical offset between translation coordinates and rotational components of p in array X
      INTEGER :: I, J, J1, J2, J3, J4, J5, J6, REALNATOMS, OFFSET

      ! X[mol][xyz] = translational coords and components of p for molecule
      ! G[mol][xyz] = gradient with respect to the variables in X
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)

      ! RI[xyz] = translational coordinates of molecule in outer loop (sum over I in eq 9)
      ! RJ[xyz] = translational coordinates of molecule in inner loop (sum over J<I in eq 9)
      ! RIJ[xyz] = separation between sites (r_ij in eq 13)
      ! NR[xyz] = unit vector parallel to r_ij
      ! RIJSQ[xyz] = |r_ij|^2
      ! ABSRIJ = |r_ij|
      ! P[xyz] = p (cf eq 2), used for both inner and outer loops
      ! THETA = used for unwrapping extra 2*pi's from P
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA

      ! I3[][] = 3x3 identity matrix (eq 1, etc)
      ! REI[xyz][] = Rep Ellipsoid for outer loop site (aka site I) in lab
      !  frame (ie, rotated SMRBF) (A before eq 33)
      ! REJ[xyz][] = Rep Ellipsoid for inner loop site (aka site J) in lab frame
      ! AEI, AEJ as above for Attractive Ellipsoids (using SMABF)
      ! RE[mol][site][][] = array for storing the Repulsive Ellipsoids (ie, the REI and REJ)
      ! AE[mol][site][][] = array for storing the Attractive Ellipsoids (ie, the AEI and AEJ)
      ! APB[][] = lambda*A . (1-lambda)*B (first term inside braces in PY eq 5)
      ! APBINV[][] = inverse of APB (first term in PY eq 5)
      ! RSTI[xyz] = position of ellipsoids in outer mol in lab frame relative to molecule position
      ! RSTJ[xyz] = position of ell in inner loop mol in lab frame rel to mol pos
      ! RST[mol][site][xyz] = array for storing the RSTI's
      ! FCNT1 = F_1(A_1,B_1) (eq 21 and after eq 25)
      ! SRTFI1 = F_1^-1/2
      ! FMIN = output from minimization problem in eq 21
      ! LAMDAC = lambda_c (after eq 22)
      ! ENERGY = energy due to PY interactions (and LJ interactions too if they exist)
      ! DUMMY = dummy variable for holding energy contributions
      ! DUMMYMAT =  dummy matrix for doing ARIBRJ calculations
      DOUBLE PRECISION :: I3(3,3), REI(3,3), REJ(3,3), AEI(3,3), AEJ(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: RE(NATOMS/2,NPYSITE,3,3), AE(NATOMS/2,NPYSITE,3,3)
      DOUBLE PRECISION :: REINV(NATOMS/2,NPYSITE,3,3), AEINV(NATOMS/2,NPYSITE,3,3)
      DOUBLE PRECISION :: RSTI(3), RSTJ(3), RST(NATOMS/2,NPYSITE,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC, ENERGY, DUMMY, DUMMYMAT(3,3)

      ! RHO1 = g_1^-1 (eq 25 & 26)
      ! RHO1SQ = g_1^-2; RHO16 ~ g_1^-6; RHO112 ~ g_1^-12
      ! RHO2, etc are as above but for g_2
      ! FCTR1 = 1/2*r_ij*F_1^(-3/2)/sig_0 = dg_1/dF_1
      ! DV1DF1 = dU^{PY}/dF_1
      ! DV1DR[xyz] = dU_repulsive/dr{vec}
      ! DRDPI1 = d|r_ij|/dp_1^I
      ! DF1PI1 = dF_1/dp_1^I (eq 32)
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DV1DF1, DV2DF2, DV1DR, DV2DR
      DOUBLE PRECISION :: DRDPI1, DRDPI2, DRDPI3, DRDPJ1, DRDPJ2, DRDPJ3 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3

      ! RMI[][] = matrix for molecule in outer loop from body to lab frame
      ! RMJ for inner loop
      ! RM[mol][][] = array for storing the RMI's
      ! DRMI1[][] = x derivative of RMI
      ! DRMJ3[][] = z deriv of RMJ
      ! DRM1[mol][][] = array that stores the DRMI1's
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), RM(NATOMS/2,3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DRM1(NATOMS/2,3,3), DRM2(NATOMS/2,3,3), DRM3(NATOMS/2,3,3)

      ! DF1DR[] = dF_1/dr_i{vec} (eq 28)
      ! DG1DR[] = dg_1/dr_i{vec} (cf first partial on RHS of eq 26)
      ! DG1DRIJ = dg_1/d|r_ij|
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3), DG1DRIJ, DG2DRIJ

      ! ARIBRJ[xyz] = lambda*A . r_i + (1-lambda)*B . r_j (2nd term in PY eq 5)
      ! ARIBRJMAT[A/B][mol][site][xyz] = stores two the A and B components in the above
      ! XC = contact point x_c (after eq 22 and PY eq 5)
      ! D1ABE = dA/dp_1^i (or j) (eq 33)
      ! D1ABEMAT[rep/att][mol][site][][] = array for storing D1ABE etc
      DOUBLE PRECISION :: ARIBRJ(3), ARIBRJMAT(2,NATOMS/2,NPYSITE,3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: D1ABEMAT(2,NATOMS/2,NPYSITE,3,3), D1ABE(3,3)
      DOUBLE PRECISION :: D2ABEMAT(2,NATOMS/2,NPYSITE,3,3), D2ABE(3,3)
      DOUBLE PRECISION :: D3ABEMAT(2,NATOMS/2,NPYSITE,3,3), D3ABE(3,3)
      LOGICAL          :: GTEST

      ! For computing the polynomial ECF:
      ! RCOF[mol][site][3x3] = cofactor matrix of the repulsive inverse lab frame shape matrices
      DOUBLE PRECISION :: RCOF(NATOMS/2,NPYSITE,3,3), ACOF(NATOMS/2,NPYSITE,3,3) !!!!

      ! [LJGSITE]
      ! R6 = (sigma_0 / |r_ij|)^6; similarly for R12
      ! LJENERGY = for storing LJ contributions to energy; similarly for LJG
      ! DFDR = f_ij'(r_ij) (eq 10)
      ! LJVT[mol] = LJ contributions to binding energy of molecule
      DOUBLE PRECISION :: R6, R12, LJENERGY, LJG(3*NATOMS), DFDR, LJVT(NATOMS/2)

      ! [CUTOFF]
      ! PCUTOFF = value of cutoff
      ! PCUT2 = PCUTOFF^-2, etc
      ! [LJGSITE] assigns PCUT2 = (pcut/lj sig_0)^-2
      ! GCUT = g value corresponding to d_r as cutoff of PCUTOFF (cf original PY paper)
      ! GCUT2 = GCUT^2, etc
      DOUBLE PRECISION :: PCUT2, PCUT6, PCUT12, GCUT, GCUT2, GCUT6, GCUT12

      ! [PBC]
      ! PBC = indicates if PBCs are active
      ! K = index for looping over images
      ! KOUT = marker for closest image
      ! PBCCOUNT = number of images to test for closeness
      ! ROUT = holder for putative min for g_1 (or max of RHO1,2)
      ! RJST[image][xyz] = "r_ij" sotre holding image locations
      ! LOUT = lambda out for lambda_c val from putative g_1 min
      LOGICAL :: PBC
      INTEGER :: K, KOUT, PBCCOUNT
      DOUBLE PRECISION :: ROUT, RJST(8,3), LOUT

      INTEGER :: I1
      DOUBLE PRECISION :: energy_contrib2,grad_contrib2

      !-------------------------!
      ! Initialize variables #B !
      !-------------------------!

      ! NATOMS is the number of lines in the coords file. Because one line specifies
      ! the mol's position and the next specifies the components of p for a molecule,
      ! the real number of atoms is NATOMS/2. Set up the offset in indices.
      ! Initialize PBC check.
      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
      PBC        = PARAMONOVPBCX.OR.PARAMONOVPBCY.OR.PARAMONOVPBCZ

      ! Set up 3x3 identity matrix I3
      I3(:,:)      = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0; 

      ! [PBC]
      ! Move molecules inside the box by subtracting multiples of
      ! the box length from the appropriate components
      IF(PBC) THEN
            DO J1 = 1, REALNATOMS
                  J3=J1*3 ! X(J3) is z comp of position of mol J1
                  IF(PARAMONOVPBCX) X(J3-2)=X(J3-2)-BOXLX*NINT(X(J3-2)/BOXLX)
                  IF(PARAMONOVPBCY) X(J3-1)=X(J3-1)-BOXLY*NINT(X(J3-1)/BOXLY)
                  IF(PARAMONOVPBCZ) X(J3)=X(J3)-BOXLZ*NINT(X(J3)/BOXLZ)
            ENDDO
      ENDIF ! [PBC]

      ! [CUTOFF]
      ! Initialize variables
      IF(PARAMONOVCUTOFF) THEN
            PCUT2 = (1.D0/PCUTOFF)**2
            PCUT6 = PCUT2**3
            PCUT12 = PCUT6**2
            GCUT = PYSIGNOT/(PCUTOFF+PYSIGNOT)
            GCUT2 = GCUT*GCUT
            GCUT6 = GCUT2**3
            GCUT12 = GCUT6**2
      ENDIF ! [CUTOFF]

      ! Initialize energies (and gradient if required)
      ENERGY = 0.D0
      IF(GTEST) G(:)   = 0.D0


      ! [LJGSITE]
      IF(LJGSITET) THEN
          LJENERGY = 0.D0
          IF(GTEST) LJG(:) = 0.D0
      ENDIF ! [LJGSITE]

    ! There are a number of calculations that depend only on the positions, orientations,
    ! and parameters of each site, so it is more efficient to do all those once at the 
    ! beginning rather than repeat them inside the site loops.
    ! Loop over all mols
    DO J1=1, REALNATOMS
        ! Use indices are described above. Then find the rotation matrix and its derivatives.
        ! Store these values
        J3 = 3*J1
        J5 = OFFSET + J3
        
        ! Get |p| to within (-2*pi,2*pi).
        P = X(J5-2:J5)
        THETA=SQRT(DOT_PRODUCT(P,P))
        IF(THETA>twopi) then
            X(J5-2:J5)=X(J5-2:J5)/THETA
            DO
              THETA=THETA-twopi
              IF(THETA<2*pi) EXIT
            ENDDO
            X(J5-2:J5)=X(J5-2:J5)*THETA
        ENDIF
        P = X(J5-2:J5)

        CALL RMDRVT(P, RM(J1,:,:), DRM1(J1,:,:), DRM2(J1,:,:), DRM3(J1,:,:), GTEST)
    
        ! Loop over all sites
        DO I=1, NPYSITE

            ! Calculate the position of a site relative to the molecule center in the
            ! lab frame (RST), and calculate the repulsive ellipsoid's shape in the lab
            ! frame (RE). 
            RST(J1,I,:)  = MATMUL(RM(J1,:,:),PST(I,:))
            RE(J1,I,:,:) = MATMUL(RM(J1,:,:),MATMUL(SMRBF(I,:,:),TRANSPOSE(RM(J1,:,:))))
            CALL MTRXINDET(RE(J1,I,:,:),SMRDET(I),REINV(J1,I,:,:))
    
            ! Compute the matrix of cofactors of REINV
            CALL COFACTORS(REINV(J1,I,:,:),RCOF(J1,I,:,:))

            ! Calculate and store the components that will be used to computer ARIBRJ later.
            ARIBRJMAT(1,J1,I,:) = MATMUL(RE(J1,I,:,:),(X(J3-2:J3)+RST(J1,I,:)))
    
            ! Calculate and store the values for D1ABE, etc that will be used later.
            DUMMYMAT = MATMUL(DRM1(J1,:,:),SMRBF(I,:,:)) ! dummy (previously D1ABEZ) = R_1^i . A^0
            D1ABEMAT(1,J1,I,:,:) = MATMUL(DUMMYMAT,TRANSPOSE(RM(J1,:,:))) + MATMUL(RM(J1,:,:),(TRANSPOSE(DUMMYMAT)))
            DUMMYMAT = MATMUL(DRM2(J1,:,:),SMRBF(I,:,:))
            D2ABEMAT(1,J1,I,:,:) = MATMUL(DUMMYMAT,TRANSPOSE(RM(J1,:,:))) + MATMUL(RM(J1,:,:),(TRANSPOSE(DUMMYMAT)))
            DUMMYMAT = MATMUL(DRM3(J1,:,:),SMRBF(I,:,:))
            D3ABEMAT(1,J1,I,:,:) = MATMUL(DUMMYMAT,TRANSPOSE(RM(J1,:,:))) + MATMUL(RM(J1,:,:),(TRANSPOSE(DUMMYMAT)))

            ! If the parameters of the repulsive and attractive ellipsoids are the same,
            ! don't bother calculating the attractive one, just copy the repulsive result.
            IF(RADIFTARRAY(I)) THEN
                AE(J1,I,:,:) = MATMUL(RM(J1,:,:),MATMUL(SMABF(I,:,:),TRANSPOSE(RM(J1,:,:))))
                CALL MTRXINDET(AE(J1,I,:,:),SMADET(I),AEINV(J1,I,:,:))

                CALL COFACTORS(AEINV(J1,I,:,:),ACOF(J1,I,:,:))

                ARIBRJMAT(2,J1,I,:) = MATMUL(AE(J1,I,:,:),(X(J3-2:J3)+RST(J1,I,:)))

                DUMMYMAT = MATMUL(DRM1(J1,:,:),SMABF(I,:,:))
                D1ABEMAT(2,J1,I,:,:) = MATMUL(DUMMYMAT,TRANSPOSE(RM(J1,:,:))) + MATMUL(RM(J1,:,:),(TRANSPOSE(DUMMYMAT)))
                DUMMYMAT = MATMUL(DRM2(J1,:,:),SMABF(I,:,:))
                D2ABEMAT(2,J1,I,:,:) = MATMUL(DUMMYMAT,TRANSPOSE(RM(J1,:,:))) + MATMUL(RM(J1,:,:),(TRANSPOSE(DUMMYMAT)))
                DUMMYMAT = MATMUL(DRM3(J1,:,:),SMABF(I,:,:)) 
                D3ABEMAT(2,J1,I,:,:) = MATMUL(DUMMYMAT,TRANSPOSE(RM(J1,:,:))) + MATMUL(RM(J1,:,:),(TRANSPOSE(DUMMYMAT)))
            ELSE
                AE(J1,I,:,:)    = RE(J1,I,:,:)
                AEINV(J1,I,:,:) = REINV(J1,I,:,:)
                ACOF(J1,I,:,:)  = RCOF(J1,I,:,:)
                ARIBRJMAT(2,J1,I,:) = ARIBRJMAT(1,J1,I,:)
                D1ABEMAT(2,J1,I,:,:) = D1ABEMAT(1,J1,I,:,:)
                D2ABEMAT(2,J1,I,:,:) = D2ABEMAT(1,J1,I,:,:)
                D3ABEMAT(2,J1,I,:,:) = D3ABEMAT(1,J1,I,:,:)
            ENDIF

        ENDDO
    ENDDO

      !------------------------------------!
      ! Begin outer loop over molecules #C !
      !------------------------------------!

      ! Loop over all molecules (cf sum over I in eq 9)
      DO J1 = 1, REALNATOMS 

      ! Set up indices J3 and J5.
      J3 = 3*J1
      J5 = OFFSET + J3

      ! Export trans coords of this molecule to RI and P respectively.
      RI = X(J3-2:J3)

      ! Pull the rotation matrix and its derivatives for this molecule
      RMI(:,:) = RM(J1,:,:)
      DRMI1(:,:) = DRM1(J1,:,:)
      DRMI2(:,:) = DRM2(J1,:,:)
      DRMI3(:,:) = DRM3(J1,:,:)

      !------------------------------------!
      ! Begin inner loop over molecules #D !
      !------------------------------------!

      ! Loop over molecules (cf sum over J<I in eq 9)
      DO J2 = J1 + 1, REALNATOMS

      ! Define indices and store trans and rot data as per the outer loop
      J4 = 3*J2
      J6 = OFFSET + J4
      RJ = X(J4-2:J4) 

      RMJ(:,:) = RM(J2,:,:)
      DRMJ1(:,:) = DRM1(J2,:,:)
      DRMJ2(:,:) = DRM2(J2,:,:)
      DRMJ3(:,:) = DRM3(J2,:,:)

      !--------------------------------!
      ! Begin outer loop over sites #E !
      !--------------------------------!

      ! Outer loop over sites in outer molecule (cf sum over i in I in eq 9)
      DO I = 1, NPYSITE

      ! Pull center of ellipsoid relative to molecule position in lab frame and
      ! the shape matrices for this site in the lab frame.
      RSTI(:) = RST(J1,I,:)
      REI(:,:) = RE(J1,I,:,:)
      AEI(:,:) = AE(J1,I,:,:)

      !--------------------------------!
      ! Begin inner loop over sites #F !
      !--------------------------------!
   
      ! Loop over sites in inner molecule (cf sum over j in J in eq 9)
      DO J = 1, NPYSITE

            RSTJ(:) = RST(J2,J,:)
            REJ(:,:) = RE(J2,J,:,:)
            AEJ(:,:) = AE(J2,J,:,:)

            ! Compute r_ij (eq 13)
            RIJ    = RI - RJ + RSTI - RSTJ
            
            !------------------------------------!
            ! Compute repulsive ellipsoid ECF #G !
            !------------------------------------!

            ! [PBC]
            ! Figure out which images will need to be tested, then
            ! calculate and compare the resulting g values.
            IF(PBC) THEN
                  ! Initially assume that the closest image is the one that is inside
                  ! the PBC box. Store original RJ.
                  PBCCOUNT = 1
                  RJST(1,:) = RJ(:)

                  ! If the both ellipsoids fit entirely inside half a box length, then
                  ! the closest image will be the one already in the box. However,
                  ! if Lx/2 < |r_x|+a+b, where a and b are the longest semiaxes of the
                  ! two ellipsoids, we might need to use an x-shifted image.
                  IF( PARAMONOVPBCX.AND. (BOXLX/2 < ABS(RIJ(1))+LONGRSAX(I)+LONGRSAX(J) ) ) THEN
                        ! Put the shifted image location in the 2nd slot. If x comp
                        ! of RIJ is negative, move it +x; etc. Increase counter.
                        RJST(2,:)=RJ(:)
                        RJST(2,1)=RJST(2,1)-SIGN(BOXLX,RJ(1))
                        PBCCOUNT=PBCCOUNT*2
                  ENDIF
                  IF (PARAMONOVPBCY.AND.(BOXLY/2<ABS(RIJ(2))+LONGRSAX(I)+LONGRSAX(J))) THEN
                        DO K=1,PBCCOUNT
                              ! If no X testing, then we are copying into 2nd slot. If there was an X test,
                              ! copy into 3rd and 4th slots.
                              RJST(K+PBCCOUNT,:)=RJST(PBCCOUNT,:)
                              ! So if no X, then RJST[1] is original, [2] is y. If there was X, RJST[1]
                              ! is orig, [2] = x, [3] = y,[4] = x and y shifted.
                              RJST(K+PBCCOUNT,2)=RJST(K+PBCCOUNT,2)-SIGN(BOXLY,RJST(K+PBCCOUNT,2))
                        ENDDO
                        PBCCOUNT=PBCCOUNT*2
                  ENDIF
                  IF (PARAMONOVPBCZ.AND.(BOXLZ/2<ABS(RIJ(3))+LONGRSAX(I)+LONGRSAX(J))) THEN
                        DO K=1,PBCCOUNT
                              RJST(K+PBCCOUNT,:)=RJST(PBCCOUNT,:)
                              ! Now we have one of these: orig, z; orig, x, xz;
                              ! orig,x,y,z,xy,xz,yz,xyz; etc, depending
                              RJST(K+PBCCOUNT,3)=RJST(K+PBCCOUNT,3)-SIGN(BOXLZ,RJST(K+PBCCOUNT,3))
                        ENDDO
                        PBCCOUNT=PBCCOUNT*2
                  ENDIF

                  ! Compute inter-ellipsoid distances. Look for the minimum.
                  DO K=1, PBCCOUNT
                        RIJ = RI - RJST(K,:) + RSTI - RSTJ
                        RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                        ABSRIJ = DSQRT(RIJSQ)

                        ! [CUTOFF]
                        ! The contact distance is at least |r_ij|-a-b, so if this is over
                        ! the cutoff then we know we don't want to use this image.
                        IF(PARAMONOVCUTOFF.AND.ABSRIJ-LONGRSAX(I)-LONGRSAX(J)>PCUTOFF) THEN
                              RHO1 = 0.D0 ! Note that any physical RHO1 > 0
                              KOUT = 1
                        ELSE
                              ! Calculate ECF
                              CALL POLYECF(REINV(J1,I,:,:), REINV(J2,J,:,:), RCOF(J1,I,:,:), RCOF(J2,J,:,:), &
                                   & SMRDET(I), SMRDET(J), RIJ, LAMDAC, FCNT1)
                        
                              NR     = RIJ / ABSRIJ
                              SRTFI1  = 1.D0 / DSQRT(FCNT1) ! ~ F_1(A_1,B_1)^-1/2 (eq 25)
                              RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT) ! ~ g_1^-1 (eq 25)
                        ENDIF ! [CUTOFF] and dist>cut

                        IF(K==1) THEN
                              ! This is the first go. Let these results be putative minimum.
                              ROUT = RHO1
                              KOUT = 1
                              LOUT = LAMDAC
                        ELSE
                              ! This is not the first go. Check for new min.
                              ! Since RHO1 ~ g_1^-1, to minimize g_1 we max RHO1
                              IF(RHO1 > ROUT) THEN
                                    ROUT = RHO1
                                    KOUT = K
                                    LOUT = LAMDAC
                              ENDIF
                        ENDIF
                  ENDDO ! loop for computing potential ECFs

                  ! Use good results
                  RHO1 = ROUT
                  RJ(:) = RJST(KOUT,:)
                  LAMDAC = LOUT
                  RIJ = RI - RJ + RSTI - RSTJ
                  RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                  ABSRIJ = DSQRT(RIJSQ)
                  
            ! No PBCs: Compute the ECF just once. r_ij was computed before
            ! the PBC block, so keep that value.
            ELSE
                  RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                  ABSRIJ = DSQRT(RIJSQ)
                  
                  ! [CUTOFF]: If distance is over cutoff, immediately set RHO1=0 as
                  ! signal that the interaction should vanish
                  IF(PARAMONOVCUTOFF.AND.ABSRIJ-LONGRSAX(I)-LONGRSAX(J)>PCUTOFF) THEN
                        RHO1 = 0.D0
                  ELSE
                        CALL POLYECF(REINV(J1,I,:,:), REINV(J2,J,:,:), RCOF(J1,I,:,:), RCOF(J2,J,:,:), &
                             & SMRDET(I), SMRDET(J), RIJ, LAMDAC, FCNT1)

                        NR     = RIJ / ABSRIJ
                        SRTFI1  = 1.D0 / DSQRT(FCNT1) ! ~ F_1(A_1,B_1)^-1/2 cf eq 25
                        RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT) ! ~ g_1^-1 cf eq 25
                  ENDIF ! [CUTOFF] and rij>cut
            ENDIF ! [PBC]

            !---------------------------------!
            ! Do repulsive ellipsoid calcs #H !
            !---------------------------------!

            ! [CUTOFF]: If we have RHO1=0, then we are over the cutoff, 
            ! so we can skip all the calculations below.
            ! Check also if g_1=1/RHO1 > g val corresponding to d_r of PCUTOFF = 1/GCUT.
            ! If it is, we are past the cutoff.
            IF(PARAMONOVCUTOFF.AND.(RHO1.EQ.0.D0.OR.RHO1<GCUT)) THEN
                  RHO1   = 0.D0
            ELSE
                  ! Calculate values needed for computing the energy
                  RHO1SQ = RHO1*RHO1
                  RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                  RHO112 = RHO16 * RHO16

                  ! Calculate values needed for computing the gradient
                  IF(GTEST) THEN
                        APB     = LAMDAC * REI + (1.D0 - LAMDAC) * REJ
                        CALL MTRXIN (APB, APBINV)

                        ARIBRJ = LAMDAC * ARIBRJMAT(1,J1,I,:) + (1.D0 - LAMDAC) * ARIBRJMAT(1,J2,J,:)
                        XC     =  MATMUL(APBINV, ARIBRJ)
                        XCMRI  = XC - RI - RSTI
                        XCMRJ  = XC - RJ - RSTJ
                        DF1DR  = - 2.D0 * LAMDAC * MATMUL(REI,XCMRI)

                        D1ABE = D1ABEMAT(1,J1,I,:,:)
                        D2ABE = D2ABEMAT(1,J1,I,:,:)
                        D3ABE = D3ABEMAT(1,J1,I,:,:)

                        DF1PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                              - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(REI,MATMUL(DRMI1,PST(I,:))))) 
                        DF1PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                              - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(REI,MATMUL(DRMI2,PST(I,:)))))
                        DF1PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                              - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(REI,MATMUL(DRMI3,PST(I,:)))))

                        D1ABE = D1ABEMAT(1,J2,J,:,:)
                        D2ABE = D2ABEMAT(1,J2,J,:,:)
                        D3ABE = D3ABEMAT(1,J2,J,:,:)
         
                        DF1PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                              - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(REJ,MATMUL(DRMJ1,PST(J,:)))))
                        DF1PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                              - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(REJ,MATMUL(DRMJ2,PST(J,:)))))
                        DF1PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                              - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(REJ,MATMUL(DRMJ3,PST(J,:))))) 

                        FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT) ! compare last term eq 27
                        DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR ! cf eq 27
                        DV1DF1 = -2.D0*pyrepp(I)*pyrepp(J)*RHO112*RHO1*FCTR1 ! cf eq 30
                        DV1DR  = -2.D0*pyrepp(I)*pyrepp(J)*RHO112*RHO1*(1.D0-SRTFI1)/PYSIGNOT

                        IF(PARAMONOVCUTOFF) DG1DRIJ = (1.D0 - SRTFI1)/PYSIGNOT
                  ENDIF ! gradients are required
            ENDIF ! [CUTOFF] exceeded or not

            !-------------------------!
            ! Compute attr ell ECF #I !
            !-------------------------!

            ! If the two shapes at site I are the same and the two shapes at site J are the same,
            ! then the calculations for the attractive ellipsoids will be the same as for the
            ! repulsive. If either of these pairs differ, we'll need to do the attr calcs
            ! separately.
            IF (RADIFTARRAY(I).OR.RADIFTARRAY(J)) THEN

                  IF(PBC) THEN
                        PBCCOUNT = 1
                        ! No need to reset RJST(1,:), since that is still at the original value
                        ! But we do need to reset RJ, RIJ for the image checks
                        RJ(:)  = RJST(1,:)
                        RIJ(:) = RI - RJ + RSTI - RSTJ

                        IF( PARAMONOVPBCX.AND. (BOXLX/2 < ABS(RIJ(1))+LONGASAX(I)+LONGASAX(J) ) ) THEN
                              RJST(2,:)=RJ(:)
                              RJST(2,1)=RJST(2,1)-SIGN(BOXLX,RJ(1))
                              PBCCOUNT=PBCCOUNT*2
                        ENDIF
                        IF (PARAMONOVPBCY.AND.(BOXLY/2<ABS(RIJ(2))+LONGASAX(I)+LONGASAX(J))) THEN
                              DO K=1,PBCCOUNT
                                    RJST(K+PBCCOUNT,:)=RJST(PBCCOUNT,:)
                                    RJST(K+PBCCOUNT,2)=RJST(K+PBCCOUNT,2)-SIGN(BOXLY,RJST(K+PBCCOUNT,2))
                              ENDDO
                              PBCCOUNT=PBCCOUNT*2
                        ENDIF
                        IF (PARAMONOVPBCZ.AND.(BOXLZ/2<ABS(RIJ(3))+LONGASAX(I)+LONGASAX(J))) THEN
                              DO K=1,PBCCOUNT
                                    RJST(K+PBCCOUNT,:)=RJST(PBCCOUNT,:)
                                    RJST(K+PBCCOUNT,3)=RJST(K+PBCCOUNT,3)-SIGN(BOXLZ,RJST(K+PBCCOUNT,3))
                              ENDDO
                              PBCCOUNT=PBCCOUNT*2
                        ENDIF

                        DO K=1, PBCCOUNT
                              RIJ = RI - RJST(K,:) + RSTI - RSTJ
                              RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                              ABSRIJ = DSQRT(RIJSQ)

                              ! [CUTOFF]
                              IF(PARAMONOVCUTOFF.AND.ABSRIJ-LONGASAX(I)-LONGASAX(J)>PCUTOFF) THEN
                                    RHO2 = 0.D0
                                    KOUT = 1
                              ELSE
                                    CALL POLYECF(AEINV(J1,I,:,:), AEINV(J2,J,:,:), ACOF(J1,I,:,:), ACOF(J2,J,:,:), &
                                         & SMADET(I), SMADET(J), RIJ, LAMDAC, FCNT2)

                                    NR     = RIJ / ABSRIJ
                                    SRTFI2 = 1.D0 / DSQRT(FCNT2) ! ~ F_1(A_1,B_1)^-1/2 cf eq 25
                                    RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT) ! ~ g_1^-1 cf eq 25
                              ENDIF ! [CUTOFF]

                              IF(K==1) THEN
                                    ROUT = RHO2
                                    KOUT = 1
                                    LOUT = LAMDAC
                              ELSE
                                    IF(RHO2 > ROUT) THEN
                                          ROUT = RHO2
                                          KOUT = K
                                          LOUT = LAMDAC
                                    ENDIF
                              ENDIF
                        ENDDO ! loop for computing potential ECFs

                        RHO2 = ROUT
                        RJ(:) = RJST(KOUT,:)
                        LAMDAC = LOUT
                        RIJ = RI - RJ + RSTI - RSTJ
                        RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                        ABSRIJ = DSQRT(RIJSQ)
                        
                  ELSE
                        ! No need to reset RIJSQ, etc since they haven't changed
                        IF(PARAMONOVCUTOFF.AND.ABSRIJ-LONGASAX(I)-LONGASAX(J)>PCUTOFF) THEN
                              RHO2 = 0.D0
                        ELSE
                              CALL POLYECF(AEINV(J1,I,:,:), AEINV(J2,J,:,:), ACOF(J1,I,:,:), ACOF(J2,J,:,:), &
                                   & SMADET(I), SMADET(J), RIJ, LAMDAC, FCNT2)
                              SRTFI2  = 1.D0 / DSQRT(FCNT2) ! ~ F_1(A_1,B_1)^-1/2 cf eq 25
                              RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT) ! ~ g_1^-1 cf eq 25
                        ENDIF
                  ENDIF ! [PBC]

                  !----------------------!
                  ! Do attr ell calcs #J !
                  !----------------------!

                  IF(PARAMONOVCUTOFF.AND.(RHO2.EQ.0.D0.OR.RHO2<GCUT)) THEN
                        RHO2   = 0.D0
                  ELSE
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
                        
                        IF(GTEST) THEN
                              APB     = LAMDAC * AEI + (1.D0 - LAMDAC) * AEJ

                              CALL MTRXIN (APB, APBINV)

                              ARIBRJ = LAMDAC * ARIBRJMAT(2,J1,I,:) + (1.D0 - LAMDAC) * ARIBRJMAT(2,J2,J,:)
                              XC     = MATMUL(APBINV, ARIBRJ)
                              XCMRI  = XC - RI - RSTI
                              XCMRJ  = XC - RJ - RSTJ
                              DF2DR  = - 2.D0 * LAMDAC * MATMUL(AEI,XCMRI)

                              FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                              DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                              DV2DF2 = pyattrp(I)*pyattrp(J)*RHO26*RHO2*FCTR2
                              DV2DR  = pyattrp(I)*pyattrp(J)*RHO26*RHO2*(1.D0-SRTFI2)/PYSIGNOT

                              D1ABE = D1ABEMAT(2,J1,I,:,:)
                              D2ABE = D2ABEMAT(2,J1,I,:,:)
                              D3ABE = D3ABEMAT(2,J1,I,:,:)

                              DF2PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                                    - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AEI,MATMUL(DRMI1,PST(I,:)))))
                              DF2PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                                    - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AEI,MATMUL(DRMI2,PST(I,:)))))
                              DF2PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                                    - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AEI,MATMUL(DRMI3,PST(I,:)))))
                  
                              D1ABE = D1ABEMAT(2,J2,J,:,:)
                              D2ABE = D2ABEMAT(2,J2,J,:,:)
                              D3ABE = D3ABEMAT(2,J2,J,:,:)

                              DF2PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                                    - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(AEJ,MATMUL(DRMJ1,PST(J,:)))))
                              DF2PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                                    - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(AEJ,MATMUL(DRMJ2,PST(J,:)))))
                              DF2PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                                    - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(AEJ,MATMUL(DRMJ3,PST(J,:)))))
                        
                              IF(PARAMONOVCUTOFF) DG2DRIJ = (1.D0 - SRTFI2)/PYSIGNOT 
                        ENDIF ! gradients are required
                  ENDIF ! [CUTOFF]

            ! If the attr ellipsoids require calculating the same ECF, just copy
            ! the results from the repulsive ones.
            ELSE
                  RHO2   = RHO1
                  RHO2SQ = RHO1SQ
                  RHO26  = RHO16
                  IF(GTEST) THEN
                        DG2DR  = DG1DR
                        DV2DF2 = pyattrp(I)*pyattrp(J)*RHO26*RHO2*FCTR1
                        DV2DR  = pyattrp(I)*pyattrp(J)*RHO26*RHO2*(1.D0-SRTFI1)/PYSIGNOT
                        DF2PI1 = DF1PI1
                        DF2PI2 = DF1PI2
                        DF2PI3 = DF1PI3
                        DF2PJ1 = DF1PJ1
                        DF2PJ2 = DF1PJ2
                        DF2PJ3 = DF1PJ3
                        IF(PARAMONOVCUTOFF) THEN
                              DG2DRIJ = DG1DRIJ
                              FCTR2   = FCTR1
                              SRTFI2  = SRTFI1
                        ENDIF
                  ENDIF ! gradients are required
            ENDIF ! check if the two sites have the same shape matrices
 
            !---------------------------------------------------!
            ! Calculate rep contributions to energy and grad #K !
            !---------------------------------------------------!

            ! If RHO1 ne 0, then the repulsive contact was not 
            ! over the cutoff and we need to add its contributions.
            IF(RHO1.NE.0.D0) THEN
                  ! Calculate the contribution of these to sites to the repuslive part
                  ! of the energy. Modulate by parameters that adjust the repulsive
                  ! contribution (cf eq 24).
                  DUMMY = pyrepp(I)*pyrepp(J)*RHO112

                  ! [CUTOFF]
                  ! Add spline terms to energy contribution
                  IF(PARAMONOVCUTOFF) THEN
                        DUMMY = DUMMY + pyrepp(I)*pyrepp(J)*GCUT12* &
                        & (6.D0*GCUT2/RHO1SQ - 7.D0) 
                  ENDIF

                  ! Add this contribution to the total energy and the pairwise
                  ! energies
                  ENERGY = ENERGY + DUMMY

                  ! Calculate the gradient. 
                  IF(GTEST) THEN
                        ! Calculate dU/dr_i (eq 26)
                        FIJ = -2.D0*pyrepp(I)*pyrepp(J)*RHO112*RHO1*DG1DR

                        DRDPI1 = DOT_PRODUCT(NR,MATMUL(DRMI1,PST(I,:)))
                        DRDPI2 = DOT_PRODUCT(NR,MATMUL(DRMI2,PST(I,:)))
                        DRDPI3 = DOT_PRODUCT(NR,MATMUL(DRMI3,PST(I,:)))

                        TIJ(1) = DV1DF1*DF1PI1 + DV1DR*DRDPI1
                        TIJ(2) = DV1DF1*DF1PI2 + DV1DR*DRDPI2
                        TIJ(3) = DV1DF1*DF1PI3 + DV1DR*DRDPI3

                        DRDPJ1 = DOT_PRODUCT(NR,MATMUL(DRMJ1,-PST(J,:)))
                        DRDPJ2 = DOT_PRODUCT(NR,MATMUL(DRMJ2,-PST(J,:)))
                        DRDPJ3 = DOT_PRODUCT(NR,MATMUL(DRMJ3,-PST(J,:)))

                        TJI(1) = DV1DF1*DF1PJ1 + DV1DR*DRDPJ1
                        TJI(2) = DV1DF1*DF1PJ2 + DV1DR*DRDPJ2
                        TJI(3) = DV1DF1*DF1PJ3 + DV1DR*DRDPJ3

                        ! [CUTOFF]: Add spline terms
                        IF(PARAMONOVCUTOFF) THEN
                              FIJ = FIJ + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1*DG1DR
                              TIJ(1) = TIJ(1) + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1* &
                               & (FCTR1*DF1PI1 + DG1DRIJ*DRDPI1)
                              TIJ(2) = TIJ(2) + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1* &
                               & (FCTR1*DF1PI2 + DG1DRIJ*DRDPI2)
                              TIJ(3) = TIJ(3) + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1* &
                               & (FCTR1*DF1PI3 + DG1DRIJ*DRDPI3)
                              TJI(1) = TJI(1) + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1* &
                               & (FCTR1*DF1PJ1 + DG1DRIJ*DRDPJ1)
                              TJI(2) = TJI(2) + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1* &
                               & (FCTR1*DF1PJ2 + DG1DRIJ*DRDPJ2)
                              TJI(3) = TJI(3) + pyrepp(I)*pyrepp(J)*GCUT12*2.D0*GCUT2/RHO1* &
                               & (FCTR1*DF1PJ3 + DG1DRIJ*DRDPJ3)
                        ENDIF ! [CUTOFF]

                        ! Add these contributions to the whole, which is the sum
                        ! over all these interacting pairs.
                        G(J3-2:J3) = G(J3-2:J3) + FIJ
                        G(J4-2:J4) = G(J4-2:J4) - FIJ
                        G(J5-2:J5) = G(J5-2:J5) + TIJ
                        G(J6-2:J6) = G(J6-2:J6) + TJI

                  ENDIF ! gradients are required
            ENDIF ! [CUTOFF] was not tripped

            !---------------------------------------------------!
            ! Calculate att contributions to energy and grad #L !
            !---------------------------------------------------!

            IF(RHO2.NE.0.D0) THEN
                  DUMMY = -pyattrp(I)*pyattrp(J)*RHO26 

                  IF(PARAMONOVCUTOFF) THEN
                        DUMMY = DUMMY + pyattrp(I)*pyattrp(J)*GCUT6* &
                        & (-3.D0*GCUT2/RHO2SQ + 4.D0) 
                  ENDIF

                  ENERGY = ENERGY + DUMMY
                  IF(GTEST) THEN
                        FIJ = pyattrp(I)*pyattrp(J)*RHO26*RHO2*DG2DR
                     
                        DRDPI1 = DOT_PRODUCT(NR,MATMUL(DRMI1,PST(I,:)))
                        DRDPI2 = DOT_PRODUCT(NR,MATMUL(DRMI2,PST(I,:)))
                        DRDPI3 = DOT_PRODUCT(NR,MATMUL(DRMI3,PST(I,:)))

                        TIJ(1) = DV2DF2*DF2PI1 + DV2DR*DRDPI1
                        TIJ(2) = DV2DF2*DF2PI2 + DV2DR*DRDPI2
                        TIJ(3) = DV2DF2*DF2PI3 + DV2DR*DRDPI3

                        DRDPJ1 = DOT_PRODUCT(NR,MATMUL(DRMJ1,-PST(J,:)))
                        DRDPJ2 = DOT_PRODUCT(NR,MATMUL(DRMJ2,-PST(J,:)))
                        DRDPJ3 = DOT_PRODUCT(NR,MATMUL(DRMJ3,-PST(J,:)))

                        TJI(1) = DV2DF2*DF2PJ1 + DV2DR*DRDPJ1
                        TJI(2) = DV2DF2*DF2PJ2 + DV2DR*DRDPJ2
                        TJI(3) = DV2DF2*DF2PJ3 + DV2DR*DRDPJ3

                        IF(PARAMONOVCUTOFF) THEN
                              FIJ = FIJ - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2*DG2DR
                              TIJ(1) = TIJ(1) - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2* &
                               & (FCTR2*DF2PI1 + DG2DRIJ*DRDPI1)
                              TIJ(2) = TIJ(2) - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2* &
                               & (FCTR2*DF2PI2 + DG2DRIJ*DRDPI2)
                              TIJ(3) = TIJ(3) - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2* &
                               & (FCTR2*DF2PI3 + DG2DRIJ*DRDPI3)
                              TJI(1) = TJI(1) - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2* &
                               & (FCTR2*DF2PJ1 + DG2DRIJ*DRDPJ1)
                              TJI(2) = TJI(2) - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2* &
                               & (FCTR2*DF2PJ2 + DG2DRIJ*DRDPJ2)
                              TJI(3) = TJI(3) - pyattrp(I)*pyattrp(J)*GCUT6*GCUT2/RHO2* &
                               & (FCTR2*DF2PJ3 + DG2DRIJ*DRDPJ3)
                        ENDIF

                        G(J3-2:J3) = G(J3-2:J3) + FIJ
                        G(J4-2:J4) = G(J4-2:J4) - FIJ
                        G(J5-2:J5) = G(J5-2:J5) + TIJ
                        G(J6-2:J6) = G(J6-2:J6) + TJI
        
                  ENDIF ! gradients are required
            ENDIF ! [CUTOFF]

      ENDDO ! inner loop over sites
      ENDDO ! outer loop over sites

      !--------------------------!
      ! LJ site contributions #M !
      !--------------------------!

      IF(LJGSITET) THEN
            ! [CUTOFF]
            IF(PARAMONOVCUTOFF) THEN
                  ! Rescale cutoff values
                  PCUT2 = (LJGSITESIGMA/PCUTOFF)**2
                  PCUT6 = PCUT2**3
                  PCUT12 = PCUT6**2
            ENDIF ! [CUTOFF]
      
            ! Begin outer loop over LJ sites
            DO I = 1, NLJSITE

            ! Calculate sep between mol position and site pos in lab frame
            RSTI(:) = MATMUL(RMI,LJGSITECOORDS(I,:))

            ! Begin inner loop over LJ sites
            DO J = 1, NLJSITE
                  RSTJ(:) = MATMUL(RMJ,LJGSITECOORDS(J,:))

                  ! Calculate separation between LJ sites
                  RIJ    = RI - RJ + RSTI - RSTJ

                  ! [PBC]
                  ! Minimize the particle separation. Unlike with ECFs, I can directly
                  ! change RIJ because the gradient calculation does not use RI and RJ separately.
                  IF(PBC) THEN
                        IF(PARAMONOVPBCX.AND.BOXLX/2<ABS(RIJ(1))) RIJ(1) = RIJ(1) - SIGN(BOXLX,RIJ(1))
                        IF(PARAMONOVPBCY.AND.BOXLY/2<ABS(RIJ(2))) RIJ(2) = RIJ(2) - SIGN(BOXLY,RIJ(2))
                        IF(PARAMONOVPBCZ.AND.BOXLZ/2<ABS(RIJ(3))) RIJ(3) = RIJ(3) - SIGN(BOXLZ,RIJ(3))
                  ENDIF

                  ! Compute |r_ij|. This will be used to check against cutoffs.
                  RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                  ABSRIJ = DSQRT(RIJSQ)

                  ! Initially assume that there is no energy contribution
                  DUMMY = 0.D0

                  ! [CUTOFF]
                  ! If there is no cutoff, or if the cutoff is on but the distance is below
                  ! cutoff, then compute basic LJ values. If over the cutoff, all these blocks
                  ! will be skipped and there will be no contribution to energy and G.
                  IF(.NOT.(PARAMONOVCUTOFF.AND.ABSRIJ>PCUTOFF)) THEN
                        NR     = RIJ / ABSRIJ

                        ! For ease of calculation, reassign RIJSQ ~ (sigma_0 / |r_ij|)^2
                        RIJSQ  = LJGSITESIGMA**2 / RIJSQ
                        R6     = RIJSQ**3
                        R12    = R6*R6

                        ! Calculate energy. Attractive and repulsive parts are modulated by
                        !  ljrepp and ljattrp
                        DUMMY = DUMMY + & 
                         & ljrepp(I)*ljrepp(J)*R12 - ljattrp(I)*ljattrp(J)* R6

                        IF(GTEST) DFDR   = (-12.D0*R12 + 6.D0*R6)/ABSRIJ ! factor of 4eps_0 will be added later
                  ENDIF ! [CUTOFF] off or dist<cut

                  ! [CUTOFF]
                  ! If cutoff is on and the distance is below cutoff, add the spline terms. See
                  ! Stoddard & Ford paper.
                  IF(PARAMONOVCUTOFF.AND.ABSRIJ<PCUTOFF) THEN
                        DUMMY = DUMMY + ljrepp(I)*ljrepp(J)* &
                         & (6.D0*PCUT2/RIJSQ - 7.D0)*PCUT12 + &
                         & ljattrp(I)*ljattrp(J)* ( &
                         & (-3.D0*PCUT2/RIJSQ + 4.D0)*PCUT6 )

                        IF(GTEST) DFDR = DFDR + (6.D0*PCUT12 - 3.D0*PCUT6)*2.0D0*ABSRIJ/PCUTOFF**2
                  ENDIF ! [CUTOFF] on and dist<cut

                  ! Add LJ contributions to energies
                  LJENERGY = LJENERGY + DUMMY

                  IF(GTEST.AND..NOT.(PARAMONOVCUTOFF.AND.ABSRIJ>PCUTOFF)) THEN
                        ! Calculate contributions to LJ gradient. Start by doing the first term in
                        ! the sum in eq 10, namely DFDR ~ f_ij'. For translation components,
                        ! (ie FIJcutoff and dist<cutoff) cf eq 10 & 11; for rotational, cf 10 & 12
                        FIJ    = NR*DFDR ! cf eq 10 & 11 
                        TIJ(1) = DOT_PRODUCT(NR,MATMUL(DRMI1,LJGSITECOORDS(I,:)))*DFDR
                        TIJ(2) = DOT_PRODUCT(NR,MATMUL(DRMI2,LJGSITECOORDS(I,:)))*DFDR
                        TIJ(3) = DOT_PRODUCT(NR,MATMUL(DRMI3,LJGSITECOORDS(I,:)))*DFDR
                        TJI(1) = DOT_PRODUCT(NR,MATMUL(DRMJ1,LJGSITECOORDS(J,:)))*DFDR
                        TJI(2) = DOT_PRODUCT(NR,MATMUL(DRMJ2,LJGSITECOORDS(J,:)))*DFDR
                        TJI(3) = DOT_PRODUCT(NR,MATMUL(DRMJ3,LJGSITECOORDS(J,:)))*DFDR

                        ! Add these new contributions to the whole.
                        LJG(J3-2:J3) = LJG(J3-2:J3) + FIJ
                        LJG(J4-2:J4) = LJG(J4-2:J4) - FIJ
                        LJG(J5-2:J5) = LJG(J5-2:J5) + TIJ
                        LJG(J6-2:J6) = LJG(J6-2:J6) - TJI
                  ENDIF ! gradients


            ENDDO ! inner loop over LJ sites
            ENDDO ! outer loop over LJ sites

         ENDIF ! [LJGSITE]

      ENDDO ! inner loop over molecules

      ENDDO ! outer loop over molecules

      !-------------!
      ! Finalize #N !
      !-------------!
      
      ! Tack on the final factors to the energy and gradient as per eq 24 & 30
      ENERGY = 4.D0*PYEPSNOT*ENERGY
      IF(GTEST) G = 24.D0*PYEPSNOT*G

     ![GRAVITY]
     do i1=1,realnatoms
!sf344> gravity field with constant gradient, and an r^-12 repulsion around z=-10 (for backwards compatibility with existing databases)
      if(efieldt) then
        energy_contrib2=0.0
        grad_contrib2=0.0 ! this is a scalar, we have only contribution along the z direction
        call gravity(x(3*i1),energy_contrib2,grad_contrib2)
        energy = energy + energy_contrib2
        g(3*i1)=g(3*i1)+grad_contrib2
      end if
     end do

      ! [LJGSITE]
      ! Add the LJ contributions to the PY contributions
      IF(LJGSITET) THEN
            ! Tack on final factors to energy and gradient
            LJENERGY = 4.D0*LJGSITEEPS*LJENERGY
            IF(GTEST) LJG = 4.D0*LJGSITEEPS*LJG 

            ENERGY = ENERGY + LJENERGY
            IF(GTEST) G = G + LJG
      ENDIF ! [LJGSITE]

      ! To freeze a molecule, make all the components of its gradient vanish.
      DO J1=1,NATOMS
            J2=3*J1
            IF(FROZEN(J1).AND.GTEST) THEN
                  G(J2-2)=0.0D0
                  G(J2-1)=0.0D0
                  G(J2  )=0.0D0
            END IF
      END DO

      END SUBROUTINE MULTISITEPY2 


!     --------------------------------------------------------------------------
! Called by finalio.f

      SUBROUTINE AAtoSites(X,P,XS,RMATBEST)
      use key, only : NPYSITE
      use pymodule, only : PST,OST,ELLST1,ELLMAT
      implicit none
      
      DOUBLE PRECISION :: X(3),P(3),XS(NPYSITE,3),P1(3),I3(3,3),phi,theta,psi, RMATBEST(3,3)
      DOUBLE PRECISION :: RMAT(3,3), RMAT1(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3), EZR1(3,3), EZR1R(3,3)
      INTEGER          :: J1
        
        CALL RMDRVT(P, RMAT, DRM1, DRM2, DRM3, .FALSE.)
      I3(:,:)=0.0D0
      DO J1=1,3
        I3(J1,J1)=1.0D0
      END DO
      IF(RMATBEST(1,1)==0.0D0.AND.RMATBEST(1,2)==0.0D0.AND.RMATBEST(1,3)==0.0D0.AND.&
         RMATBEST(2,1)==0.0D0.AND.RMATBEST(2,2)==0.0D0.AND.RMATBEST(2,3)==0.0D0.AND.&
         RMATBEST(3,1)==0.0D0.AND.RMATBEST(3,2)==0.0D0.AND.RMATBEST(3,3)==0.0D0) THEN
        RMATBEST(:,:)=I3(:,:)
      END IF
      DO J1=1,NPYSITE
        CALL RMDRVT(OST(J1,:), RMAT1, DRM1, DRM2, DRM3, .FALSE.)

        XS(J1,:)=MATMUL(RMATBEST,MATMUL(RMAT,PST(J1,:)))
        XS(J1,1)=XS(J1,1) + X(1)
        XS(J1,2)=XS(J1,2) + X(2)
        XS(J1,3)=XS(J1,3) + X(3)

!        AE1 = MATMUL(RMAT,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMAT)))))
! original matrix for angle-axis part

         EZR1(1,1) = 1.D0/(ELLST1(J1,1)*ELLST1(J1,1))
         EZR1(2,2) = 1.D0/(ELLST1(J1,2)*ELLST1(J1,2))
         EZR1(3,3) = 1.D0/(ELLST1(J1,3)*ELLST1(J1,3))

! now rotate

        EZR1R = MATMUL(RMAT1,(MATMUL(I3,(TRANSPOSE(RMAT1)))))

! now rotate again

!        ELLMAT(J1,:,:) = MATMUL(RMAT,(MATMUL(EZR1R,(TRANSPOSE(RMAT)))))
        ELLMAT(J1,:,:) = MATMUL(RMATBEST,MATMUL(RMAT,RMAT1))
!        WRITE(*,*) 'body ', J1
!        WRITE(*,*) 'original matrix: '
!        WRITE(*,*) EZR1(:,:)

!        WRITE(*,*) 'rotation matrix: '
!        WRITE(*,*) ELLMAT(J1,:,:)


!   XS: cartesian coordinates of each site
!   P: angle-axis coordinates of the whole body
!        PS(J1,:)=MATMUL(RMAT,OST(J1,:))
      END DO

      END SUBROUTINE AAtoSites

!     --------------------------------------------------------------------------

SUBROUTINE MULTISITEPYSECDER(OLDX,STEST)
use commons
!use keyword
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST

GTEST=.TRUE.
ksi=0.00001D0
X(:)=OLDX(:)

!WRITE(*,*) 'sf344> in Paramonovsecder', NATOMS, X(1)
!X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO i=1,3*NATOMS

!       IF ((i.GT.3*NATOMS/2).and.(MOD(i,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(i)=X(i)-ksi
 
!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL MULTISITEPY2(X,V,EGB,.TRUE.)
!         WRITE(*,*) 'EGB1',EGB
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL MULTISITEPY2(X,V,EGB,.TRUE.)
!         WRITE(*,*) 'EGB2',EGB

         VTEMP(2,:)=V(:)

                DO j=i,3*NATOMS
                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                        HESS(j,i)=HESS(i,j)
                END DO
!        END IF
END DO
!WRITE(*,*) 'sf344> exiting Paramonovsecder', NATOMS, X(1)
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'(12F10.3)') HESS(:,:)
END SUBROUTINE MULTISITEPYSECDER


!-------------------------------------------!
! BEGIN SUBROUTINE TAKESTEPMULTISITEPY (#T) !
!-------------------------------------------!

! This subroutine takes a random MC step for the multisite PY potential.
! First, each molecule is moved in the xyz coordinates and rotated.
! Second, the most weakly bound molecule is identified and an extra angular
! step is taken for that molecule. The algorithm for finding the angular
! step is borrowed from gay-berne.f90's TAKESTEPELLIPSOIDS with some
! modifications.

SUBROUTINE TAKESTEPMULTISITEPY (COORDS)

! NATOMS, ..., MYUNIT = cf above
! COORDS[xyzp][step] = configuration of molecules in this step of run
! PYOVERLAPTHRESH = min val allowed for ECF before overlap is assumed
! RADIUS = radius of spherical container for molecules
! TMOVE[step] = logical if this step has a translational move
! OMOVE[step] = if this step has orientational (ie p vector) move
! STEP[step] = size of translational move
! OSTEP[step] = size of orientational move
! PERCOLATET = logical if percolate keyword is used
! PERCCUT = length for percolation cutoff, if used
! VAT[mol] = binding energy of the molecule
! CENT = logical if CENTRE keyword is used
USE COMMONS, ONLY: NATOMS
USE KEY, ONLY: NPYSITE, PYLOCALSTEP, PYOVERLAPTHRESH

! PST, ... = cf DEFINEPYMULTISITES above
USE PYMODULE, ONLY : PST, OST, SMRBF, SMABF, PI, LONGRSAX, SMRDET

IMPLICIT NONE

! I, ..., OFFSET = cf MULTISITEPY2 above
! NP = this step, passed from arguments
! J1MAX = index for molecule that is most weakly bound
INTEGER :: I, J, J1, J2, J3, J4, J5, J6, K1, K3, REALNATOMS, OFFSET
INTEGER :: NP, J1MAX

! X, ..., D3ABE = cf above
! DRM1, etc are just dummies because we don't need them
DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), COORDS(3*NATOMS)
DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), RIJSQ, ABSRIJ, P(3)
DOUBLE PRECISION :: I3(3,3), REI(3,3), REJ(3,3)
DOUBLE PRECISION :: RSTI(3), RSTJ(3), FCNT1, FMIN, LAMDAC
DOUBLE PRECISION :: RMI(3,3), RMJ(3,3)
DOUBLE PRECISION :: DRM1(3,3), DRM2(3,3), DRM3(3,3)
DOUBLE PRECISION :: RST(NATOMS/2,NPYSITE,3), RE(NATOMS/2,NPYSITE,3,3), RM(NATOMS/2,3,3)
DOUBLE PRECISION :: REINV(NATOMS/2,NPYSITE,3,3), RCOF(NATOMS/2,NPYSITE,3,3)

! RANDOM = random real between -1 and +1
! DPRAND = double precision random real between 0 and 1 (cf dprand.f)
! LOCALSTEP = local dummy for step sizes
! SAVECOORDS[xyzp] = local backup for COORDS
! CM[xyz] = center of mass of cluster
! CMMAX = max{distance of mol from center of mass}
! DUMMY = dummy variable for holding radii, etc for comparisons
! THETA, PHI = angles for angular step
! VMAX = max{binding energies}, so weakest binding
! OVERLAPT = logical flag for overlap
! RADIUSGOOD = logical flag for whether the angular step was OK
DOUBLE PRECISION :: RANDOM, DPRAND, LOCALSTEP, SAVECOORDS(3*NATOMS)
DOUBLE PRECISION :: CM(3), CMMAX, DUMMY, THETA, PHI, VMAX, STEPARRAY(NATOMS)
LOGICAL          :: OVERLAPT, RADIUSGOOD, OVERLAPARRAY(NATOMS)

REALNATOMS = NATOMS/2
OFFSET     = 3*REALNATOMS

! Define identity matrix
I3(:,:) = 0.D0
I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0


! Calculate rotation matrices. cf main subroutine for an explanation
DO J1=1, REALNATOMS
  J3 = 3*J1
  J5 = OFFSET + J3
  P  = COORDS(J5-2:J5)
  CALL RMDRVT(P, RM(J1,:,:), DRM1, DRM2, DRM3, .FALSE.)
  DO I=1, NPYSITE
      RST(J1,I,:)  = MATMUL(RM(J1,:,:),PST(I,:))
      RE(J1,I,:,:) = MATMUL(RM(J1,:,:),MATMUL(SMRBF(I,:,:),TRANSPOSE(RM(J1,:,:))))
      CALL MTRXINDET(RE(J1,I,:,:),SMRDET(I),REINV(J1,I,:,:))
      CALL COFACTORS(REINV(J1,I,:,:),RCOF(J1,I,:,:))
  ENDDO
ENDDO

LOCALSTEP = 0.0D0
STEPARRAY(:)=0.0D0
OVERLAPARRAY(:)=.FALSE.

! Loop over all molecules
DO J1=1, REALNATOMS-1 

      J3 = 3*J1
      J5 = OFFSET + J3

      ! If the overlap is over the threshold, then return to this line
95    CONTINUE

      !-------------------!
      ! Take initial step !
      !-------------------!
      
      ! Initially set translation step to zero. Then choose a new translation step
      ! from the STEP array.
      LOCALSTEP = 0.0D0
      STEPARRAY(:)=0.0D0
      IF (OVERLAPARRAY(J1)) THEN
              OVERLAPARRAY(J1)=.FALSE.
              ! Create random numbers and use them to shift the original translational coords
              RANDOM          = (DPRAND() - 0.5D0)*2.0D0
              COORDS(J3-2) = COORDS(J3-2) + PYLOCALSTEP(1)*RANDOM
              RANDOM          = (DPRAND() - 0.5D0)*2.0D0
              COORDS(J3-1) = COORDS(J3-1) + PYLOCALSTEP(1)*RANDOM
              RANDOM          = (DPRAND() - 0.5D0)*2.0D0
              COORDS(J3)   = COORDS(J3) + PYLOCALSTEP(1)*RANDOM
        
              ! Repeat above for orientational coordinates (ie, p vector).
              RANDOM          = (DPRAND() - 0.5D0)*2.0D0
              COORDS(J5-2) = COORDS(J5-2) + PYLOCALSTEP(2)*RANDOM
              RANDOM          = (DPRAND() - 0.5D0)*2.0D0
              COORDS(J5-1) = COORDS(J5-1) + PYLOCALSTEP(2)*RANDOM
              RANDOM          = (DPRAND() - 0.5D0)*2.0D0
              COORDS(J5)   = COORDS(J5) + PYLOCALSTEP(2)*RANDOM
              ! Recompute the rotation matrices for the ellipsoid that has been rotated
              P  = COORDS(J5-2:J5)
              CALL RMDRVT(P, RM(J1,:,:), DRM1, DRM2, DRM3, .FALSE.)
              DO I=1, NPYSITE
                  RST(J1,I,:)  = MATMUL(RM(J1,:,:),PST(I,:))
                  RE(J1,I,:,:) = MATMUL(RM(J1,:,:),MATMUL(SMRBF(I,:,:),TRANSPOSE(RM(J1,:,:))))
                  CALL MTRXINDET(RE(J1,I,:,:),SMRDET(I),REINV(J1,I,:,:))
                  CALL COFACTORS(REINV(J1,I,:,:),RCOF(J1,I,:,:))
              ENDDO
      END IF
      ! Pull the values for the current positions from the altered COORDS array
      RI = COORDS(J3-2:J3)

      !-----------------------!
      ! Now check for overlap !
      !-----------------------!

      ! Pull rotation matrix for outer molecule
      RMI(:,:) = RM(J1,:,:)

      ! Begin inner loop over molecules
      DO J2 = 1, REALNATOMS

            ! Do not check a molecule against itself
            IF(J2==J1) CYCLE

            ! Grab indices and coords for inner loop molecule from arrays
            J4 = 3*J2
            J6 = OFFSET + J4
            RJ = COORDS(J4-2:J4) 

            ! Pull rotation matrix for inner loop molecule
            RMJ(:,:) = RM(J2,:,:)

            ! Loop over all sites in the outer molecule
            DO I = 1, NPYSITE

                  ! Pull center of ell relative to molecule position in lab frame
                  ! Pull repulsive shape matrix for site in lab frame
                  RSTI(:) = RST(J1,I,:)
                  REI(:,:) = RE(J1,I,:,:)
               
                  ! Loop over sites in inner molecule
                  DO J = 1, NPYSITE

                     ! Find relative center of mass and shape mat
                     RSTJ(:) = RST(J2,J,:)
                     REJ(:,:) = RE(J2,J,:,:)

                     ! Calculate separation between the two sites
                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)

                     ! If the distance between the two centers is greater than the sum of
                     ! the longest semiaxes of the two ellipsoids, then there will be
                     ! no overlap, so we can just skip on to the next pair. Otherwise we
                     ! need to calculate the ECF and check if the overlap is there.
                     IF(ABSRIJ<LONGRSAX(I)+LONGRSAX(J)) THEN
                         ! Compute the ECF for these two sites
                         CALL POLYECF(REINV(J1,I,:,:), REINV(J2,J,:,:), RCOF(J1,I,:,:), RCOF(J2,J,:,:), &
                              & SMRDET(I), SMRDET(J), RIJ, LAMDAC, FCNT1)

                         ! Allow an overlap threshold, but if the ECF is too small, then
                         ! flag for overlap and go back to the beginning of the while loop.
                         IF (FCNT1 < PYOVERLAPTHRESH) THEN
                            OVERLAPARRAY(J1) = .TRUE.
!                            WRITE(*,*) 'takestepmultisitepy> bodies overlapping', J1, J2
                            GO TO 95
                         ENDIF
                     ENDIF

                  ENDDO  ! end loop over sites in the second body

            ENDDO ! end loop over sites in the first body

      ENDDO ! end inner mol loop

ENDDO ! end outer mol loop


END SUBROUTINE TAKESTEPMULTISITEPY 




!-------------------------------!
! BEGIN SUBROUTINE POLYECF (#Z) !
!-------------------------------!

! This subroutine computes the ECF using the polynomial method (cf
! Scott Olesen's "computing the ECF using polynomials", probably in my
! MPhil thesis by the time I leave). The subroutine takes the two matrices,
! their cofactor matrices, and their determinants and computes the
! coefficients of the polynomial representations of S(lambda). Then
! it uses Newton's method to find the zero of the derivative of the
! numerator of S (with a modificiation to deal with the fact that the
! root of S' must lie within [0,1].

SUBROUTINE POLYECF(A, B, AT, BT, DETA, DETB, X, X1, SMAX)
implicit none
! A[][] = inverse shape matrix of one ellipsoid
! AT[][] = cofactor matrix of A
! B[][] = inverse shape matrix of other ellipsoid
! BT[][] = cofactor matrix of B
! DETA = determinant of A
! X[] = separation between ellipsoid centers
! X1 = x_n in Newton's method, and output value for lambda
! SMAX = S(lambda_c)
DOUBLE PRECISION :: A(3,3), AT(3,3), B(3,3), BT(3,3)
DOUBLE PRECISION :: DETA, DETB, X(3), X1, SMAX

! AS = A*, etc
DOUBLE PRECISION :: AS, BS, CS, DE, ES, DS

! F1 = f_1, coefficient of numerator of S(lambda)
! G0 = g_0, coefficient of denominator of S(l)
! H0 = h_0, coefficient of numerator of S'(l)
DOUBLE PRECISION :: F1, F2, F3, F4, G0, G1, G2, G3
DOUBLE PRECISION :: H0, H1, H2, H3, H4, H5, H6

! X0 = x_(n-1) in Newton's method, test values for lambda
! XX = placeholder for powers of X0 (aka lambda)
! P  = value of polynomial h at lambda=X0
! DP = value of derivative of h
! N = numerator of S
! D = denominator of S
DOUBLE PRECISION :: X0, XX, P, DP, N, D

! I = counter for Newton's method iterations
! CONVERGED = if X1-X0 got sufficiently small
INTEGER :: I
LOGICAL :: CONVERGED

! Compute A* = Sum_ij X_i X_j AT_ij (noting that AT is symmetric)
AS =      X(1)*X(1)*AT(1,1) +      X(2)*X(2)*AT(2,2) +      X(3)*X(3)*AT(3,3) + &
   & 2.D0*X(1)*X(2)*AT(1,2) + 2.D0*X(1)*X(3)*AT(1,3) + 2.D0*X(2)*X(3)*AT(2,3) 

! Compute B* = Sum_ij X_i X_j BT_ij (noting that BT is symmetric)
BS =      X(1)*X(1)*BT(1,1) +      X(2)*X(2)*BT(2,2) +      X(3)*X(3)*BT(3,3) + &
   & 2.D0*X(1)*X(2)*BT(1,2) + 2.D0*X(1)*X(3)*BT(1,3) + 2.D0*X(2)*X(3)*BT(2,3) 

! Compute C* = Sum_ij X_i X_j CT_ij (noting that CT is symmetric)
CS =      X(1)*X(1)*(A(2,2)*B(3,3) - A(2,3)*B(3,2) + B(2,2)*A(3,3) - B(2,3)*A(3,2)) + &
   &      X(2)*X(2)*(A(3,3)*B(1,1) - A(3,1)*B(1,3) + B(3,3)*A(1,1) - B(3,1)*A(1,3)) + &
   &      X(3)*X(3)*(A(1,1)*B(2,2) - A(1,2)*B(2,1) + B(1,1)*A(2,2) - B(1,2)*A(2,1)) + &
   & 2.D0*X(1)*X(2)*(A(2,3)*B(3,1) - A(2,1)*B(3,3) + B(2,3)*A(3,1) - B(2,1)*A(3,3)) + &
   & 2.D0*X(1)*X(3)*(A(2,1)*B(3,2) - A(2,2)*B(3,1) + B(2,1)*A(3,2) - B(2,2)*A(3,1)) + &
   & 2.D0*X(2)*X(3)*(A(3,1)*B(1,2) - A(3,2)*B(1,1) + B(3,1)*A(1,2) - B(3,2)*A(1,1))

! Compute D* = Sum_ij A_ij BT_ij (noting symmetry)
DS =      A(1,1)*BT(1,1) +      A(2,2)*BT(2,2) +      A(3,3)*BT(3,3) + &
   & 2.D0*A(1,2)*BT(1,2) + 2.D0*A(1,3)*BT(1,3) + 2.D0*A(2,3)*BT(2,3)

! Compute E* = Sum_ij B_ij AT_ij (noting symmetry)
ES =      B(1,1)*AT(1,1) +      B(2,2)*AT(2,2) +      B(3,3)*AT(3,3) + &
   & 2.D0*B(1,2)*AT(1,2) + 2.D0*B(1,3)*AT(1,3) + 2.D0*B(2,3)*AT(2,3)

! Compute coefficients of f = numerator of S
F1 = AS
F2 = -3.D0*AS + CS
F3 = 3.D0*AS + BS - 2.D0*CS
F4 = - AS - BS + CS

! Compute coefficients of g = denominator of S
G0 = DETA
G1 = -3.D0*DETA + ES
G2 = 3.D0*DETA + DS - 2.D0*ES
G3 = - DETA + DETB - DS + ES

! Compute coefficients of h = numerator of S'
H0 = F1*G0
H1 = 2.D0*F2*G0
H2 = 3.D0*F3*G0 + F2*G1 - F1*G2
H3 = 4.D0*F4*G0 + 2.D0*F3*G1 - 2.D0*F1*G3
H4 = 3.D0*F4*G1 + F3*G2 - F2*G3
H5 = 2.D0*F4*G2
H6 = F4*G3

! Begin iteration loop. Guess initial position.
X0 = 0.51D0
CONVERGED = .FALSE.
DO I = 1, 25
    ! Compute polynomial and its derivative: First, add the constant parts. Each
    ! step adds the next power onto both, and XX increases to a new power of X0
    ! each time.
    P  = H0
    DP = H1
    XX = X0
    P  = P  + H1*XX
    DP = DP + 2.D0*H2*XX
    XX = XX*X0
    P  = P  + H2*XX
    DP = DP + 3.D0*H3*XX
    XX = XX*X0
    P  = P  + H3*XX
    DP = DP + 4.D0*H4*XX
    XX = XX*X0
    P  = P  + H4*XX
    DP = DP + 5.D0*H5*XX
    XX = XX*X0
    P  = P  + H5*XX
    DP = DP + 6.D0*H6*XX
    XX = XX*X0
    P  = P  + H6*XX

    ! Compute new position. If DP=0, we are at a critical point and need to 
    ! do some extra work.
    IF(DP .NE. 0.D0) THEN
        ! Take the normal Newton-Rapheson step
        X1 = X0 - P / DP

        ! If X1 ended up outside [0,1], then pretend that X0 was the first point
        ! computed in a bisection method. So if P>0, go right, and if P<0 go left.
        IF (X1 .LT. 0.D0 .OR. X1 .GT. 1.D0) THEN
            IF(P .GT. 0.D0) X1 = (X0 + 1.D0)/2.D0
            IF(P .LT. 0.D0) X1 =  X0/2.D0

        ! If the new position is sufficiently close to the old one, we say we found
        ! the root.
        ELSEIF (ABS(X1 - X0) .LT. 1.D-8) THEN
            CONVERGED = .TRUE.
            EXIT ! leave DO loop
        ENDIF

    ! If DP=0, take a bisection step.
    ELSE
        IF(P .GT. 0.D0) X1 = (X0 + 1.D0)/2.D0
        IF(P .LT. 0.D0) X1 =  X0/2.D0
    ENDIF

    ! Update old position
    X0 = X1
ENDDO

! Compute numerator and denominator of S(lambda_c)
D  = G0
XX = X1
N  = F1*XX
D  = D + G1*XX
XX = XX*X1
N  = N + F2*XX
D  = D + G2*XX
XX = XX*X1
N  = N + F3*XX
D  = D + G3*XX
XX = XX*X1
N  = N + F4*XX

! Compute S(lambda_c), which gets returned as a parameter
SMAX = N / D

! If we did not find a root, notify the user. If this happens,
! something has gone very wrong!
IF (.NOT. CONVERGED) THEN
   STOP 'polyecf> Newton did not converge'
END IF

RETURN

END SUBROUTINE POLYECF


SUBROUTINE COFACTORS(A, COFA)

! This subroutine takes a 3x3 matrix A and returns the matrix
! of its cofactors COFA

DOUBLE PRECISION :: A(3,3), COFA(3,3)

COFA(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
COFA(1,2) = A(2,3)*A(3,1) - A(3,3)*A(2,1)
COFA(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
COFA(2,1) = COFA(1,2)
COFA(2,2) = A(3,3)*A(1,1) - A(1,3)*A(3,1)
COFA(2,3) = A(3,1)*A(1,2) - A(1,1)*A(3,2)
COFA(3,1) = COFA(1,3)
COFA(3,2) = COFA(2,3)
COFA(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

RETURN

END SUBROUTINE COFACTORS


SUBROUTINE INERTIAPY(RMI, KBLOCK, TMASS)
! sf344> first implementation of the inertia routine for multisite PY particles
!        applicable only for those centred around the origin
!        copied over from INERTIAPAP

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN)  :: RMI(3,3)
        DOUBLE PRECISION, INTENT(OUT) :: TMASS, KBLOCK(3,3)
! LFCT: half the equilibrium pair distance squared
! MOMENT: scalar moment of inertia of a sphere about any axis
        DOUBLE PRECISION :: LFCT, MOMENT

! Set the total mass as unity
        TMASS  = 1.D0
! Set the value of LFCT to half the equilibrium pair distance squared
        LFCT = 1.44D0
! Set the value of MOMENT using the formula I= 0.2*M*R*R
        MOMENT = 0.2D0*TMASS*LFCT

! Set the moment of inertia tensor to zero.
! Off diagonal elements will remain zero
        KBLOCK(:,:) = 0.D0

! Set the diagonal elements to MOMENT
        KBLOCK(1,1) = MOMENT
        KBLOCK(2,2) = MOMENT
        KBLOCK(3,3) = MOMENT

      END SUBROUTINE INERTIAPY

