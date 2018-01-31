      SUBROUTINE PYG (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
      USE KEY

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5, J6, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: AEZR1(3,3), AEZR2(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC1, LAMDAC2, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DVDF1, DVDF2 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: DEZ(3,3), D1E(3,3), D2E(3,3), D3E(3,3) 
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: SE1(3*NATOMS/2,3), SD1E1(3*NATOMS/2,3), SD2E1(3*NATOMS/2,3), SD3E1(3*NATOMS/2,3)
      DOUBLE PRECISION :: SE2(3*NATOMS/2,3), SD1E2(3*NATOMS/2,3), SD2E2(3*NATOMS/2,3), SD3E2(3*NATOMS/2,3)
      LOGICAL          :: GTEST, STEST

      ENERGY = 0.D0
      IF (GTEST) G(:)   = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      AEZR1(:,:) = 0.D0
      AEZR2(:,:) = 0.D0

      DO J1 = 1, 3

         AEZR1(J1,J1) = 1.D0/(PYA1(J1)*PYA1(J1))
         AEZR2(J1,J1) = 1.D0/(PYA2(J1)*PYA2(J1))

      ENDDO

      DO J1 = 1, REALNATOMS

         J3      = 3*J1
         J5      = OFFSET + J3
         RI      = X(J3-2:J3)
         P       = X(J5-2:J5)

!     ROTATION MATRIX AND ITS DERIVATIVES

         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         SE1(J3-2:J3,:)   = MATMUL(RMI(:,:),(MATMUL(AEZR1(:,:),(TRANSPOSE(RMI(:,:))))))

         IF (RADIFT) THEN

            SE2(J3-2:J3,:)   = MATMUL(RMI(:,:),(MATMUL(AEZR2(:,:),(TRANSPOSE(RMI(:,:))))))

         ENDIF

         IF (GTEST) THEN

            DEZ(:,:)         = MATMUL(DRMI1(:,:),AEZR1(:,:))
            SD1E1(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI2(:,:),AEZR1(:,:))
            SD2E1(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI3(:,:),AEZR1(:,:))
            SD3E1(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))

            IF (RADIFT) THEN

            DEZ(:,:)         = MATMUL(DRMI1(:,:),AEZR2(:,:))
            SD1E2(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI2(:,:),AEZR2(:,:))
            SD2E2(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
            DEZ(:,:)         = MATMUL(DRMI3(:,:),AEZR2(:,:))
            SD3E2(J3-2:J3,:) = MATMUL(DEZ(:,:),(TRANSPOSE(RMI(:,:)))) + MATMUL(RMI(:,:),(TRANSPOSE(DEZ(:,:))))
 
            ENDIF

         ENDIF

      ENDDO

      DO J1 = 1, REALNATOMS - 1

         J3      = 3*J1
         J5      = OFFSET + J3
         RI      = X(J3-2:J3)

         AE1(:,:) = SE1(J3-2:J3,:)

         IF (RADIFT) THEN

            AE2(:,:) = SE2(J3-2:J3,:)         

         ENDIF

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4
            RJ = X(J4-2:J4) 

            BE1(:,:) = SE1(J4-2:J4,:)

            IF (RADIFT) THEN
   
               BE2(:,:) = SE2(J4-2:J4,:)

            ENDIF

            RIJ(:) = RI(:) - RJ(:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            ABSRIJ = DSQRT(RIJSQ)
            NR(:)  = RIJ(:)/ABSRIJ

!     CALCULATE ECF

            CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC1, FMIN)

            FCNT1   = -FMIN
            SRTFI1  = 1.D0 / DSQRT(FCNT1)

            RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
            RHO1SQ = RHO1*RHO1
            RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
            RHO112 = RHO16 * RHO16

            IF (RADIFT) THEN

               CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC2, FMIN)

               FCNT2    = -FMIN
               SRTFI2   = 1.D0/DSQRT(FCNT2)
               RHO2     = PYSIGNOT/(ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
               RHO2SQ   = RHO2*RHO2
               RHO26    = RHO2SQ*RHO2SQ*RHO2SQ

            ELSE

               RHO2   = RHO1
               RHO26  = RHO16

            ENDIF

            ENERGY = ENERGY + RHO112 - RHO26

            IF (GTEST) THEN

               APB(:,:) = LAMDAC1*AE1(:,:) + (1.D0-LAMDAC1)*BE1(:,:)

               CALL MTRXIN (APB(:,:), APBINV(:,:))

               ARIBRJ(:) = LAMDAC1*MATMUL(AE1(:,:),RI(:)) + (1.D0-LAMDAC1)*MATMUL(BE1(:,:),RJ(:))
               XC(:)     = MATMUL(APBINV(:,:), ARIBRJ(:))
               XCMRI(:)  = XC(:) - RI(:)
               XCMRJ(:)  = XC(:) - RJ(:)
               DF1DR(:)  = -2.D0*LAMDAC1*MATMUL(AE1(:,:),XCMRI(:))
               FCTR1     = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
               DG1DR(:)  = (1.D0-SRTFI1)*NR(:)/PYSIGNOT + FCTR1*DF1DR(:)
               DVDF1     = -2.D0*RHO112*RHO1*FCTR1
                
               D1E(:,:)  = SD1E1(J3-2:J3,:)
               D2E(:,:)  = SD2E1(J3-2:J3,:)
               D3E(:,:)  = SD3E1(J3-2:J3,:)

               DF1PI1    = LAMDAC1*DOT_PRODUCT(XCMRI(:),MATMUL(D1E(:,:),XCMRI(:)))
               DF1PI2    = LAMDAC1*DOT_PRODUCT(XCMRI(:),MATMUL(D2E(:,:),XCMRI(:)))
               DF1PI3    = LAMDAC1*DOT_PRODUCT(XCMRI(:),MATMUL(D3E(:,:),XCMRI(:)))

               D1E(:,:)  = SD1E1(J4-2:J4,:)
               D2E(:,:)  = SD2E1(J4-2:J4,:)
               D3E(:,:)  = SD3E1(J4-2:J4,:)
 
               DF1PJ1    = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ(:),MATMUL(D1E(:,:),XCMRJ(:)))
               DF1PJ2    = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ(:),MATMUL(D2E(:,:),XCMRJ(:)))
               DF1PJ3    = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ(:),MATMUL(D3E(:,:),XCMRJ(:)))

               IF (RADIFT) THEN

                  APB(:,:) = LAMDAC2*AE2(:,:) + (1.D0-LAMDAC2)*BE2(:,:)

                  CALL MTRXIN (APB(:,:), APBINV(:,:))

                  ARIBRJ(:) = LAMDAC2*MATMUL(AE2(:,:),RI(:)) + (1.D0-LAMDAC2)*MATMUL(BE2(:,:),RJ(:))
                  XC(:)     = MATMUL(APBINV(:,:), ARIBRJ(:))
                  XCMRI(:)  = XC(:) - RI(:)
                  XCMRJ(:)  = XC(:) - RJ(:)
                  DF2DR(:)  = -2.D0*LAMDAC2*MATMUL(AE2(:,:),XCMRI(:))
                 
                  FCTR2     = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                  DG2DR(:)  = (1.D0-SRTFI2)*NR(:)/PYSIGNOT + FCTR2*DF2DR(:)
                  DVDF2     = RHO26*RHO2*FCTR2

                  D1E(:,:)  = SD1E2(J3-2:J3,:)
                  D2E(:,:)  = SD2E2(J3-2:J3,:)
                  D3E(:,:)  = SD3E2(J3-2:J3,:)

                  DF2PI1    = LAMDAC2*DOT_PRODUCT(XCMRI(:),MATMUL(D1E(:,:),XCMRI(:)))
                  DF2PI2    = LAMDAC2*DOT_PRODUCT(XCMRI(:),MATMUL(D2E(:,:),XCMRI(:)))
                  DF2PI3    = LAMDAC2*DOT_PRODUCT(XCMRI(:),MATMUL(D3E(:,:),XCMRI(:)))

                  D1E(:,:)  = SD1E2(J4-2:J4,:)
                  D2E(:,:)  = SD2E2(J4-2:J4,:)
                  D3E(:,:)  = SD3E2(J4-2:J4,:)

                  DF2PJ1    = (1.D0-LAMDAC2)*DOT_PRODUCT(XCMRJ(:),MATMUL(D1E(:,:),XCMRJ(:)))
                  DF2PJ2    = (1.D0-LAMDAC2)*DOT_PRODUCT(XCMRJ(:),MATMUL(D2E(:,:),XCMRJ(:)))
                  DF2PJ3    = (1.D0-LAMDAC2)*DOT_PRODUCT(XCMRJ(:),MATMUL(D3E(:,:),XCMRJ(:)))

               ELSE

                  DG2DR(:)  = DG1DR(:)
                  DVDF2     = RHO26*RHO2*FCTR1
                  DF2PI1    = DF1PI1
                  DF2PI2    = DF1PI2
                  DF2PI3    = DF1PI3
                  DF2PJ1    = DF1PJ1
                  DF2PJ2    = DF1PJ2
                  DF2PJ3    = DF1PJ3

               ENDIF

!     CALCULATE GRADIENT

               FIJ(:) = 2.D0*RHO112*RHO1*DG1DR(:) - RHO26*RHO2*DG2DR(:)
               TIJ(1) = DVDF1*DF1PI1 + DVDF2*DF2PI1
               TIJ(2) = DVDF1*DF1PI2 + DVDF2*DF2PI2
               TIJ(3) = DVDF1*DF1PI3 + DVDF2*DF2PI3
               TJI(1) = DVDF1*DF1PJ1 + DVDF2*DF2PJ1
               TJI(2) = DVDF1*DF1PJ2 + DVDF2*DF2PJ2
               TJI(3) = DVDF1*DF1PJ3 + DVDF2*DF2PJ3

               G(J3-2:J3) = G(J3-2:J3) - FIJ(:)
               G(J4-2:J4) = G(J4-2:J4) + FIJ(:)
               G(J5-2:J5) = G(J5-2:J5) + TIJ(:)
               G(J6-2:J6) = G(J6-2:J6) + TJI(:)

            ENDIF
          
         ENDDO

      ENDDO

      ENERGY = 4.D0*PYEPSNOT*ENERGY
      IF (GTEST) G(:) = 24.D0*PYEPSNOT*G(:)

      END SUBROUTINE PYG 

!     --------------------------------------------------------------------------

      SUBROUTINE MTRXIN(S, SINV)

      IMPLICIT NONE

      INTEGER            :: I, J
      INTEGER, PARAMETER :: M = 3
      DOUBLE PRECISION   :: S(M,M), SINV(M,M), A(M,M), DET, INVDET

!     ADJOINT

      A(1,1) = S(2,2) * S(3,3) - S(2,3) * S(2,3)
      A(2,2) = S(1,1) * S(3,3) - S(1,3) * S(1,3)
      A(3,3) = S(1,1) * S(2,2) - S(1,2) * S(1,2)
      A(1,2) = S(1,3) * S(2,3) - S(1,2) * S(3,3)
      A(1,3) = S(1,2) * S(2,3) - S(1,3) * S(2,2)
      A(2,3) = S(1,2) * S(1,3) - S(1,1) * S(2,3)

!     DETERMINANT

      DET    =  S(1,1)* A(1,1) + S(1,2) * A(1,2) + S(1,3) * A(1,3)
      INVDET = 1.D0 / DET

      DO I = 1, 3

         DO J = 1, 3

            IF (I > J) A(I,J) = A(J,I)
            SINV(I,J) = A(I,J) * INVDET

         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE MTRXIN

!     --------------------------------------------------------------------------

      SUBROUTINE OBJCTF(AE, BE, RIJ, LAMDA, SLMD)

      IMPLICIT NONE

      INTEGER          :: I, J
      DOUBLE PRECISION :: AE(3,3), BE(3,3), AEINV(3,3), BEINV(3,3), RIJ(3)
      DOUBLE PRECISION :: MG(3,3), MGINV(3,3), MGINVR(3)
      DOUBLE PRECISION :: LAMDA, SLMD

      CALL MTRXIN(AE, AEINV)
      CALL MTRXIN(BE, BEINV)

      MG = (1.D0 - LAMDA) * AEINV + LAMDA * BEINV

      CALL MTRXIN(MG, MGINV)

      MGINVR  =  MATMUL(MGINV, RIJ) 

      SLMD =  - LAMDA * (1.D0 - LAMDA) * DOT_PRODUCT(RIJ,MGINVR)

      RETURN
      END SUBROUTINE OBJCTF

!     --------------------------------------------------------------------------

      SUBROUTINE BRENTMIN (AX, BX, CX, AE, BE, RIJ, XMIN, FMIN)

      IMPLICIT NONE

      INTEGER            :: ITR
      INTEGER, PARAMETER :: ITRMX = 2000
      DOUBLE PRECISION   :: MA(3,3), MB(3,3)
      DOUBLE PRECISION   :: AX, BX, CX, A, B, D, E, P, Q, R, U, V, W, X, XM
      DOUBLE PRECISION   :: XMIN, FX, FU, FV, FW, FMIN, F, ETMP
      DOUBLE PRECISION   :: AE(3,3), BE(3,3), RIJ(3)
      DOUBLE PRECISION   :: TOL1, TOL2, CGOLD
      DOUBLE PRECISION, PARAMETER :: TOL = 1.D-11, ZEPS = 1.D-12

      CGOLD = 0.5D0 * (3 - DSQRT(5.D0))
      A     = MIN(AX,CX)
      B     = MAX(AX,CX)
      V     = BX
      W     = V
      X     = V
      E     = 0.D0
!     FX    = F(X)

      CALL OBJCTF (AE, BE, RIJ, X, FX)

      FV    = FX
      FW    = FX

      DO 10 ITR = 1, ITRMX

         XM   = 0.5D0 * (A + B)
         TOL1 = TOL * ABS(X) + ZEPS
         TOL2 = 2.D0 * TOL1
         IF (ABS(X - XM) <= (TOL2 - 0.5D0 * (B - A))) GOTO 3
         IF (ABS(E) > TOL1) THEN
            R    = (X - W) * (FX - FV)
            Q    = (X - V) * (FX - FW)
            P    = (X - V) * Q - (X - W) * R
            Q    = 2.D0 * (Q - R)
            IF (Q > 0.D0) P = -P
            Q    = ABS(Q)
            ETMP = E
            E    = D
            IF (ABS(P) >= ABS(0.5D0*Q*ETMP) .OR. P <= Q*(A-X) .OR. P >= Q*(B-X)) GOTO 1
!     The above conditions determine the acceptability of the parabolic fit. Here it is o.k.
            D = P / Q ! Take the parabolic step.
            U = X + D
            IF(U-A < TOL2 .OR. B-U < TOL2) D=SIGN(TOL1,XM-X)
            GOTO 2 !Skip over the golden section step.

         ENDIF

1        IF (X >= XM) THEN
!      We arrive here for a golden section step, which we take
            E = A - X !into the larger of the two segments.
         ELSE
            E = B - X
         ENDIF

         D = CGOLD * E !Take the golden section step.
2        IF (ABS(D) >= TOL1) THEN ! Arrive here with d computed either from parabolic fit, or
            U = X + D !else from golden section.
         ELSE
            U = X + SIGN(TOL1,D)
         ENDIF

!         FU = F(U) !This is the one function evaluation per iteration,

         CALL OBJCTF (AE, BE, RIJ, U, FU)

         IF (FU <= FX) THEN !and now we have to decide what to do with our function
            IF (U >= X) THEN !evaluation. Housekeeping follows:
               A = X
            ELSE
               B = X
            ENDIF
            V  = W
            FV = FW
            W  = X
            FW = FX
            X  = U
            FX = FU
         ELSE
            IF (U < X) THEN
               A = U
            ELSE
               B = U
            ENDIF

            IF (FU <= FW .OR. W == X) THEN
               V  = W
               FV = FW
               W  = U
               FW = FU
            ELSE IF (FU <= FV .OR. V == X .OR. V == W) THEN
               V = U
               FV = FU
            ENDIF

         ENDIF ! Done with housekeeping. Back for another iteration.

10    CONTINUE
      PRINT*, 'brent exceed maximum iterations'
3     XMIN = X ! Arrive here ready to exit with best values.
      FMIN = FX
!      WRITE(*,*) 'exiting brentmin, ITR=',ITR
      RETURN
      END SUBROUTINE BRENTMIN

      SUBROUTINE MTRXINDET(S, INVDET, SINV)

      ! This subroutine inverts a symmetric matrix whose determinant is known.
      ! Works by finding the upper-left part of the adjugate matrix, then
      ! divides by the determinant to get the inverse.

      ! Note that input is the inverse of the determinant!

      IMPLICIT NONE

      INTEGER, PARAMETER :: M = 3
      DOUBLE PRECISION   :: S(M,M), SINV(M,M), A(M,M), INVDET

      ! Calculate the upper half of the adjugate matrix
      A(1,1) = S(2,2) * S(3,3) - S(2,3) * S(2,3)
      A(2,2) = S(1,1) * S(3,3) - S(1,3) * S(1,3)
      A(3,3) = S(1,1) * S(2,2) - S(1,2) * S(1,2)
      A(1,2) = S(1,3) * S(2,3) - S(1,2) * S(3,3)
      A(1,3) = S(1,2) * S(2,3) - S(1,3) * S(2,2)
      A(2,3) = S(1,2) * S(1,3) - S(1,1) * S(2,3)

      ! Assign elements to matrix. The adjugate is symmetric,
      ! so assign half the values that way.
      SINV(1,1) = INVDET * A(1,1)
      SINV(1,2) = INVDET * A(1,2)
      SINV(1,3) = INVDET * A(1,3)
      SINV(2,1) = SINV(1,2)
      SINV(2,2) = INVDET * A(2,2)
      SINV(2,3) = INVDET * A(2,3)
      SINV(3,1) = SINV(1,3)
      SINV(3,2) = SINV(2,3)
      SINV(3,3) = INVDET * A(3,3)

      RETURN
      END SUBROUTINE MTRXINDET

