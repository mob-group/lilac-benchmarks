! **************************************************************************
     
      SUBROUTINE AMB_GETKDNAT(KD)
!get width of diagonal for GCS matrix, for natural internals
!for now, assume the dihedrals or bonds give the maximum width
      USE MODAMBER9


      INTEGER KD,ATOMWIDTH,DIFF,IPHI
      INTEGER I,J,K,L
!
      ATOMWIDTH=0
      DIFF=0
      DO IPHI=1,NPHIA
        I=ix(i50+IPHI)/3 + 1
        J=ix(i52+IPHI)/3 + 1
        K=IABS(ix(i54+IPHI))/3 + 1
        L=IABS(ix(i56+IPHI))/3 + 1
!     print *,'pdihe I J K L',I,J,K,L
        DIFF=ABS(I-J)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(I-K)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(I-L)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(J-K)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(J-L)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(K-L)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
      ENDDO
      DO IPHI=1,NPHIH
        I=ix(i40+IPHI)/3 + 1
        J=ix(i42+IPHI)/3 + 1
        K=IABS(ix(i44+IPHI))/3 + 1
        L=IABS(ix(i46+IPHI))/3 + 1
!     print *,'pdihe I J K L',I,J,K,L
        DIFF=ABS(I-J)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(I-K)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(I-L)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(J-K)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(J-L)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
        DIFF=ABS(K-L)
        IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
      ENDDO

      IF(.NOT.ALLOCATED(IBIB)) THEN
         PRINT*, "ib going wrong!"
         ALLOCATE(IBIB(nbona+nbonh))
         DO I =1,nbona
            IBIB(I) = ix(iiba+I-1)/3+1
         ENDDO
         DO I=1, nbonh
            IBIB(nbona+I) = ix(iibh+I-1)/3+1
         ENDDO
      ENDIF

      IF(.NOT.ALLOCATED(JBJB)) THEN
          ALLOCATE(JBJB(nbona+nbonh))
          DO I =1,nbona
             JBJB(I) = ix(ijba+I-1)/3+1
          ENDDO
          DO I=1, nbonh
             JBJB(nbona+I) = ix(ijbh+I-1)/3+1
          ENDDO
      ENDIF

      DO IPHI = 1,nbona+nbonh
         I = IBIB(IPHI)
         J = JBJB(IPHI)
         DIFF=ABS(I-J)
         IF (DIFF.GT.ATOMWIDTH) ATOMWIDTH=DIFF
      ENDDO
!
      KD=3*ATOMWIDTH+2
!     print *,'amb_getkdnat> ',KD
!
      RETURN
      END


! ***************************************************************************

      SUBROUTINE AMBGETANGLE(I, J, K, THETA, dTdI, dTdJ, dTdK, NOCOOR,NODERV)
! given the coordinates for three atoms I, J, K
! find the angle with J as the central atom; return it in THETA
! Also return all derivatives: dT/dIx, dT/dIy, dT/dIz in triplet dTdI
!same with dTdJ and dTdK

      LOGICAL NOCOOR,NODERV
      DOUBLE PRECISION I(3), J(3), K(3), THETA
      DOUBLE PRECISION dTdI(3), dTdJ(3), dTdK(3)
      DOUBLE PRECISION DXI, DYI, DZI, DXJ, DYJ, DZJ, RI2, RJ2, RI, RJ, RIR, RJR
      DOUBLE PRECISION DXIR, DYIR, DZIR, DXJR, DYJR, DZJR, CST, ISNT

      DXI=I(1)-J(1)
      DYI=I(2)-J(2)
      DZI=I(3)-J(3)
      DXJ=K(1)-J(1)
      DYJ=K(2)-J(2)
      DZJ=K(3)-J(3)
      RI2=DXI*DXI+DYI*DYI+DZI*DZI
      RJ2=DXJ*DXJ+DYJ*DYJ+DZJ*DZJ
      RI=SQRT(RI2)
      RJ=SQRT(RJ2)
      RIR=1/RI
      RJR=1/RJ
      DXIR=DXI*RIR
      DYIR=DYI*RIR
      DZIR=DZI*RIR
      DXJR=DXJ*RJR
      DYJR=DYJ*RJR
      DZJR=DZJ*RJR
      
      CST=DXIR*DXJR+DYIR*DYJR+DZIR*DZJR
      IF (.NOT.NOCOOR) THETA = ACOS(CST)

      IF (.NOT.NODERV) THEN
      ISNT=1/SQRT(1-CST*CST)
      dTdI(1)=(DXIR*CST-DXJR)*RIR*ISNT
      dTdI(2)=(DYIR*CST-DYJR)*RIR*ISNT
      dTdI(3)=(DZIR*CST-DZJR)*RIR*ISNT
      dTdJ(1)=(DXIR*(RI-RJ*CST)+DXJR*(RJ-RI*CST))*RIR*RJR*ISNT
      dTdJ(2)=(DYIR*(RI-RJ*CST)+DYJR*(RJ-RI*CST))*RIR*RJR*ISNT
      dTdJ(3)=(DZIR*(RI-RJ*CST)+DZJR*(RJ-RI*CST))*RIR*RJR*ISNT
      dTdK(1)=(DXJR*CST-DXIR)*RJR*ISNT
      dTdK(2)=(DYJR*CST-DYIR)*RJR*ISNT
      dTdK(3)=(DZJR*CST-DZIR)*RJR*ISNT
      ENDIF

      RETURN
      END SUBROUTINE

      SUBROUTINE AMBGETOUTOFPLANE(I, J, K, L, PHI, dPdI, dPdJ, dPdK, dPdL, NOCOOR,NODERV)
!EFK:
!given the coordinates for four atoms I, J, K, L with I the central atom
!find the angle between J and the plane defined by I,K,L; return it in PHI
!Also return all derivatives: dP/dIx, dP/dIy, dP/dIz in triplet dPdI
!same with dPdJ, dPdK, and dPdL
!Calculations from Wilson, Delcius, and Cross (WDC) Ch.4
!
      LOGICAL NOCOOR, NODERV
      DOUBLE PRECISION I(0:2), J(0:2), K(0:2), L(0:2), PHI
      DOUBLE PRECISION dPdI(0:2), dPdJ(0:2), dPdK(0:2), dPdL(0:2)
      DOUBLE PRECISION JIX, JIY, JIZ, KIX, KIY, KIZ, LIX, LIY, LIZ
      DOUBLE PRECISION RJI2, RKI2, RLI2, RJI, RKI, RLI, RJIR, RKIR, RLIR
      DOUBLE PRECISION AX, AY, AZ, BX, BY, BZ, CX, CY, CZ
      DOUBLE PRECISION RA2, RA, RB2, RB, RC2, RC, RCR
      DOUBLE PRECISION SNP, CSP2, CSP, CSPR, TNP, DUMMY
      INTEGER X

!get vector differences
      JIX=J(0)-I(0) ! r41 in WDC
      JIY=J(1)-I(1)
      JIZ=J(2)-I(2)
      KIX=K(0)-I(0) ! r42 in WDC
      KIY=K(1)-I(1)
      KIZ=K(2)-I(2)
      LIX=L(0)-I(0) ! r43 in WDC
      LIY=L(1)-I(1)
      LIZ=L(2)-I(2)
!Get vector norms
      RJI2 = JIX*JIX+JIY*JIY+JIZ*JIZ
      RKI2 = KIX*KIX+KIY*KIY+KIZ*KIZ
      RLI2 = LIX*LIX+LIY*LIY+LIZ*LIZ
      RJI = SQRT(RJI2)
      RKI = SQRT(RKI2)
      RLI = SQRT(RLI2)
      RJIR = 1/RJI
      RKIR = 1/RKI
      RLIR = 1/RLI
!Get A = r41 x r42
      AX = JIY*KIZ-KIY*JIZ
      AY = -JIX*KIZ+KIX*JIZ
      AZ = JIX*KIY-KIX*JIY
!Get B = r43 x r41
      BX = LIY*JIZ-JIY*LIZ
      BY = -LIX*JIZ +JIX*LIZ
      BZ = LIX*JIY-JIX*LIY
!Get C = r42 x r43
      CX = KIY*LIZ-LIY*KIZ
      CY = -KIX*LIZ+LIX*KIZ
      CZ = KIX*LIY-LIX*KIY
      RC2 = CX*CX+CY*CY+CZ*CZ
      RC = SQRT(RC2)
      RCR = 1/RC
!sin(phi1) = |r42 x r43| / (|r42| *|r43|)
!DUMMY = (r42 x r43)*r41
      DUMMY = CX*JIX+CY*JIY+CZ*JIZ
!sin(PHI) = DUMMY / (sin(phi1) * |r42||r43||r41|)
      SNP = DUMMY * RCR*RJIR 
!NOTE: POSSIBLE BUG; may run into problems if angle over 90 degrees
!should fix this later
      IF (.NOT.NOCOOR) PHI = ASIN(SNP)

      IF (NODERV) RETURN

!Calculate derivatives for B matrix
      CSP2 = 1-SNP*SNP
      CSP = SQRT(CSP2) !cos(PHI)
      CSPR = 1/CSP !sec(PHI)
      TNP = SNP*CSPR !tan(PHI)
!POTENTIAL BUG: singularity at PHI->pi/2
!
      dPdJ(0) = RJIR*(CX*CSPR*RCR - JIX*TNP*RJIR)
      dPdJ(1) = RJIR*(CY*CSPR*RCR - JIY*TNP*RJIR)
      dPdJ(2) = RJIR*(CZ*CSPR*RCR - JIZ*TNP*RJIR)
!
!CONTINUE CHECKING WITH 2nd TERM HERE....
      DUMMY = LIX*KIX+LIY*KIY+LIZ*KIZ !r43*r42
      dPdK(0) = BX*CSPR*RCR*RJIR - KIX*TNP*RLI2*RCR*RCR &
     &     + LIX*TNP*DUMMY*RCR*RCR
      dPdK(1) = BY*CSPR*RCR*RJIR - KIY*TNP*RLI2*RCR*RCR &
     &     + LIY*TNP*DUMMY*RCR*RCR
      dPdK(2) = BZ*CSPR*RCR*RJIR - KIZ*TNP*RLI2*RCR*RCR &
     &     + LIZ*TNP*DUMMY*RCR*RCR
!    
      dPdL(0) = AX*CSPR*RCR*RJIR - LIX*TNP*RKI2*RCR*RCR &
     &     + KIX*TNP*DUMMY*RCR*RCR
      dPdL(1) = AY*CSPR*RCR*RJIR - LIY*TNP*RKI2*RCR*RCR &
     &     + KIY*TNP*DUMMY*RCR*RCR
      dPdL(2) = AZ*CSPR*RCR*RJIR - LIZ*TNP*RKI2*RCR*RCR &
     &     + KIZ*TNP*DUMMY*RCR*RCR
!
      DO X = 0,2
         dPdI(X) = -dPdJ(X)-dPdK(X)-dPdL(X)
      ENDDO
!
      RETURN
      END

      SUBROUTINE AMBGETTORSION(I, J, K, L, PHI, dPdI, dPdJ, dPdK, dPdL, NOCOOR, NODERV)
!given the coordinates for four atoms I, J, K, L, bound in order, 
!find the dihedral torsion angle; return it in PHI
!Also return all derivatives: dP/dIx, dP/dIy, dP/dIz in triplet dPdI
!same with dPdJ, dPdK, and dPdL

      LOGICAL NOCOOR, NODERV
      DOUBLE PRECISION I(3), J(3), K(3), L(3), PHI
      DOUBLE PRECISION dPdI(3), dPdJ(3), dPdK(3), dPdL(3)
      DOUBLE PRECISION FX, FY, FZ, GX, GY, GZ, HX, HY, HZ
      DOUBLE PRECISION AX, AY, AZ, BX, BY, Bz
      DOUBLE PRECISION RF, RG, RH, RF2, RG2, RH2, RFR, RGR, RHR
      DOUBLE PRECISION CSTTWO, SNTTWO2, CSTTHREE, SNTTHREE2, SNTTWO2R, SNTTHREE2R
      DOUBLE PRECISION RA2, RB2, RA2R, RB2R, RABR, CP
      DOUBLE PRECISION MYTX, MYTY, MYTZ, MYSCALAR
      DOUBLE PRECISION DUMMY, DUMMY2

      FX=I(1)-J(1)
      FY=I(2)-J(2)
      FZ=I(3)-J(3)
      GX=J(1)-K(1)
      GY=J(2)-K(2)
      GZ=J(3)-K(3)
      HX=L(1)-K(1)
      HY=L(2)-K(2)
      HZ=L(3)-K(3)
!A=F^G, B=H^G
      AX=FY*GZ-FZ*GY
      AY=FZ*GX-FX*GZ
      AZ=FX*GY-FY*GX
      BX=HY*GZ-HZ*GY
      BY=HZ*GX-HX*GZ
      BZ=HX*GY-HY*GX
!RG=|G|, RGR=1/|G|
      RG2=GX*GX+GY*GY+GZ*GZ
      RG=SQRT(RG2)
      RGR=1/RG
!dae for use in evaluating B-matrix
      RF2=FX*FX+FY*FY+FZ*FZ
      RF=SQRT(RF2)
      RFR=1/RF
      RH2=HX*HX+HY*HY+HZ*HZ
      RH=SQRT(RH2)
      RHR=1/RH
!    jmc        CSTTWO=(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
      CSTTWO=-(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
      SNTTWO2=1-CSTTWO*CSTTWO
      SNTTWO2R=1/SNTTWO2
      CSTTHREE=(HX*GX+HY*GY+HZ*GZ)*RHR*RGR
      SNTTHREE2=1-CSTTHREE*CSTTHREE
      SNTTHREE2R=1/SNTTHREE2      
!
      IF (.NOT.NOCOOR) THEN
         RA2=AX*AX+AY*AY+AZ*AZ
         RB2=BX*BX+BY*BY+BZ*BZ
         RA2R=1/RA2
         RB2R=1/RB2
         RABR=SQRT(RA2R*RB2R)
!    CP=cos(phi)
         CP=(AX*BX+AY*BY+AZ*BZ)*RABR
         IF (CP.GT.1.0D0) THEN
            PHI=0.0D0
         ELSEIF (CP.LT.-1.0D0) THEN
            PHI=PI
         ELSE
            PHI=ACOS(CP)
         ENDIF

!    jmc Test to determine whether the dihedral angle should be > 0 or < 0, from unres (intcor.f)
!    > 0 if cw rotation from 2-> 1 to 3-> 4.
         MYTX=AY*BZ-BY*AZ
         MYTY=-AX*BZ+BX*AZ
         MYTZ=AX*BY-BX*AY
         MYSCALAR=MYTX*GX+MYTY*GY+MYTZ*GZ
         IF (MYSCALAR.GT.0.0D0) PHI = -PHI
      ENDIF

      IF (NODERV) RETURN

!dae calculate B-matrix elements
      DUMMY=RFR*RFR*RGR*SNTTWO2R
      dPdI = (/-AX*DUMMY, -AY*DUMMY, -AZ*DUMMY/)

      DUMMY=RFR*RFR*RGR*RGR*SNTTWO2R*(RG-RF*CSTTWO)
      DUMMY2=RHR*RGR*RGR*SNTTHREE2R*CSTTHREE
      dPdJ = (/AX*DUMMY-BX*DUMMY2, AY*DUMMY-BY*DUMMY2, AZ*DUMMY-BZ*DUMMY2/)

      DUMMY=RHR*RHR*RGR*RGR*SNTTHREE2R*(RG-RH*CSTTHREE)
      DUMMY2=RFR*RGR*RGR*SNTTWO2R*CSTTWO
      dPdK = (/-BX*DUMMY+AX*DUMMY2, -BY*DUMMY+AY*DUMMY2,-BZ*DUMMY+AZ*DUMMY2/)

      DUMMY=RHR*RHR*RGR*SNTTHREE2R
      dPdL = (/BX*DUMMY,BY*DUMMY,BZ*DUMMY/)

      RETURN
      END SUBROUTINE


!************************************************************************************************

      SUBROUTINE AMBALIGNDIH(PHI,DIHC,RELPHI,PREVRELPHI,N)
!aligning dihedral angles properly
!if aligndir is set, and it's not the first dihedral in the group (N.NE.1) then 
!adjust angle to minimise change in relative dihedral to the first in the group
!otherwise adjust angle to minimise change from previous value of same angle

      USE INTCOMMONS, ONLY : PREVDIH, ALIGNDIR

      IMPLICIT NONE

      DOUBLE PRECISION PHI, RELPHI, PREVRELPHI, DRP, PHIORIG
      INTEGER DIHC, N
      DOUBLE PRECISION PI, TWOPI
      PARAMETER(PI=3.141592653589793D0)

      TWOPI = 2.0D0*PI
      PHIORIG = PHI

      IF (ALIGNDIR.AND.N.NE.1) THEN
         DRP = RELPHI-PREVRELPHI
         IF (DRP.GT.PI) THEN
            PHI = PHI - TWOPI
         ELSE IF (DRP.LT.-PI) THEN
            PHI = PHI + TWOPI
         ENDIF         
      ELSE
         IF (PHI-PREVDIH(DIHC).GT.PI) THEN
            PHI = PHI - TWOPI
         ELSE IF (PHI-PREVDIH(DIHC).LT.-PI) THEN
            PHI = PHI + TWOPI
         ENDIF
      ENDIF

      PREVDIH(DIHC) = PHI

      RETURN
      END SUBROUTINE


!*********************************************************************
      SUBROUTINE AMBGETNNZNAT(NNZ)
!    get number of nonzero elements in B matrix, for natural internals      
      USE COMMONS, ONLY : NATOMS
      USE INTCOMMONS, ONLY : CARTATMSTART, ATOMxPOSinC, CENTERS, NBDS,&
     &     NIMP, NFRG, NCRT, NCNT, LINDIH, RINGS, NRNG, NLDH

      IMPLICIT NONE

      INTEGER NNZ, I, NANGC, CNTTYPE, NATMSC

      NNZ = 6*NBDS+12*NIMP+18*NFRG + 3*NCRT
      DO I=1,NCNT
         CNTTYPE = CENTERS(I,0)
         IF (CNTTYPE.LT.3) THEN 
            NNZ = NNZ + ATOMxPOSinC(CNTTYPE,5,5) + 3
         ELSE IF (CNTTYPE.EQ.3) THEN
            NNZ = NNZ + ATOMxPOSinC(CNTTYPE,5,4) + 3
         ELSE IF (CNTTYPE.EQ.4.OR.CNTTYPE.EQ.5) THEN
            NNZ = NNZ + ATOMxPOSinC(CNTTYPE,2,4) + 3
         ELSE IF (CNTTYPE.EQ.6) THEN
            NNZ = NNZ + ATOMxPOSinC(CNTTYPE,1,3) + 3
         ELSE IF (CNTTYPE.EQ.7) THEN
            NNZ = NNZ + ATOMxPOSinC(CNTTYPE,1,4) + 3
         ELSE IF (CNTTYPE.EQ.8) THEN
            NNZ = NNZ +ATOMxPOSinC(CNTTYPE,4,5) +3
!    ELSE IF (CNTTYPE.EQ.8) THEN
!    NNZ = NNZ + ATOMxPOSinC(3,4,4) + 3
!    ELSE IF (CNTTYPE.EQ.9) THEN
!    NNZ = NNZ + ATOMxPOSinC(2,4,5) + 3
         ENDIF
      ENDDO
      DO I=1,NRNG
         IF(RINGS(I,6).EQ.-1) THEN
            NNZ = NNZ + 3*5*4
         ELSE
            NNZ = NNZ + 3*6*6
         ENDIF
      ENDDO
      DO I=1,NLDH
         NNZ = NNZ + 6 + 3*(LINDIH(I,-1)+LINDIH(I,0))
      ENDDO
!     print*, 'getnnznat>> nonzeros found is:', NNZ
      RETURN
      END


!***********************************************************************
!
!subroutine to do internal coordinate manipulations using B-matrix
!calls GETBEE to calcualte internal coordinates and B-matrix
!then transforms cartesian gradients to internal gradients
!Nemeth et al. JCP 113,5598(2000)
!

      SUBROUTINE AMBTRANSFORM(XCART,GCART,XINT,GINT,NINTC,NCART,NNZ,NOCOOR,NODERV,KD,INTEPSILON)
      USE modmxatms, ONLY: MXATMS
 
       USE SPFUNCTS, ONLY : DUMPCOORDS

       IMPLICIT NONE
!SAT: two implicit none - line below commented!
!#INCLUDE '~/charmmcode/fcm/dimens.fcm'
!#INCLUDE '~/charmmcode/fcm/psf.fcm'
!NNZ is number of non-zero elements in B-matrix
!NINTC is number of internal coordinates
!NCART is number of cartesians
!
       INTEGER NNZ,NINTC,NCART,KD,K
!
       DOUBLE PRECISION XCART(3*MXATMS),GCART(3*MXATMS),XINT(NINTC),GINT(NINTC)
       DOUBLE PRECISION VAL(NNZ)
       INTEGER INDX(NNZ),JNDX(NNZ)
       INTEGER TRANSA,M,N,DESCRA(9),I1,J1,K1,LDGCS,LDA,LDB,INFO,NRHS,LWORK,ATOMWIDTH
       DOUBLE PRECISION ALPHA,BETA,WORK(NINTC),GCS(KD+1,NCART),GCS2(2*KD+1,NCART)
       DOUBLE PRECISION INTEPSILON,GX(NCART),GY(NINTC),SUM,SUMGY2,MEANSUMGY2,GCSREAL
       CHARACTER*1 UPLO
       LOGICAL NOCOOR,LCONVERGE, NODERV
       DOUBLE PRECISION :: GCARTORIG(NCART), GNORM, GCARTINPUT(NCART)
       INTEGER :: NITS

!---------------------------------

!get B-matrix and internal coordiantes
!

       CALL AMB_GETBEENAT(XCART,XINT,VAL,INDX,JNDX,NINTC,NCART, NNZ,GCS2,NOCOOR,NODERV,KD)
       IF (NODERV) RETURN

       GINT(1:NINTC) = 0.0D0
!     print *,'tf here1'
!
!VAL is column of non-zero BEE elements indexed by INDX and JNDX
!GCS2 = BEEtBEE as calculated in GETBEE, stored in sparse format
!need to convert double sparse format of GCS2 to lapack upperdiagonal only GCS
!KD is number of superdiagonals of GC
!see dpbtrf.f http://www.netlib.org/lapack/double/dpbtrf.f
!     
!move to array of dimension KD+1
       DO J1=1,KD
          DO I1=1,J1
             GCS(KD+1+I1-J1,J1)=GCS2(KD+1+I1-J1,J1)
          ENDDO
       ENDDO
       DO J1=KD+1,NCART
          DO I1=J1-KD,J1
             GCS(KD+1+I1-J1,J1)=GCS2(KD+1+I1-J1,J1)
          ENDDO
       ENDDO

!
!XINT now contains internal coordinates
!now need to transform gradients
!
!form GC = GC + eI where e is small
!
       DO I1=1,NCART
          GCS(KD+1,I1)=GCS(KD+1,I1)+INTEPSILON
       ENDDO

!
!Cholesky decompose GCS
!use lapack routine DPBTRF
!LDGCS is leading dimension of GCS
!
          UPLO='U'
          N=NCART
          LDGCS=KD+1

          CALL DPBTRF(UPLO,N,KD,GCS,LDGCS,INFO)

          IF (INFO.NE.0) PRINT*,'WARNING - after DPBTRF INFO=',INFO
!     print *,'tf here2'
!
!GCS is now the upper triangular matrix after decomposition
!SOLVE G*GCART2=GCART FOR GCART2, WHERE G IS ORIGINAL GC(=BEETBEE)+INTEPSILON*I
!done using lapack routine DPBTRS with input of the upper matrix from
!cholesky decomposition
!
          NRHS=1
          LDB=NCART
          CALL DPBTRS(UPLO,N,KD,NRHS,GCS,LDGCS,GCART,LDB,INFO)
          IF (INFO.NE.0) PRINT*,'WARNING - after DPBTRS INFO=',INFO

!    print *,'tf here3'
!
!GCART is now GCART2
!GINT = BEE*GCART as first approximation
!use DCOOMM to do this sparse multiplication
!sparse matrix multiplication routine from BLAS
!http://math.nist.gov/~KRemington/fspblas/
!
       TRANSA=0
       M=NINTC
       N=1
       K=NCART
       ALPHA=1.0D0
       BETA=0.0D0
       DO I1=1,9
         DESCRA(I1)=0
       ENDDO
         DESCRA(4)=1
       LWORK=NINTC
!     print *,'tf here3.5'
       CALL DCOOMM(TRANSA,M,N,K,ALPHA,DESCRA,VAL,INDX,JNDX,NNZ,GCART,NCART,BETA,GINT,NINTC,WORK,LWORK)
       IF (INFO.NE.0) PRINT*,'WARNING - after DCOOMM INFO=',INFO

!      IF(ITERATEG) PREITGRAD(:) = GINT(:)
!     print *,'tf here4'
!
!can skip below iteration and use this first approximation
!
!print out GINT
!     DO I1=1,NINTC
!      print *,'GINT',GINT(I1)
!     ENDDO
!C
!now iterate GINT
!GINT(K+1)=GINT(K)+At(GRAD-BEEt*GINT(K))
!multiplication At*GX is the same as solving G*GCART2=GX
!C
!     print *,'transform here4'
!      LCONVERGE = .FALSE.
!      NITS = 0
!      IF (ITERATEG) THEN
!         DO WHILE (.NOT.LCONVERGE)
!    GX=BEEt*GINT
!GX=GCART-GX
!            TRANSA=1
!            M=NINTC; N = 1; K = NCART
!            ALPHA=1.0D0
!            BETA=0.0D0
!            DO I1=1,9
!               DESCRA(I1)=0
!            ENDDO
!            DESCRA(4)=1
!            LWORK=NINTC
!            GX(:) = 0.0D0
!            CALL DCOOMM(TRANSA,M,N,K,ALPHA,DESCRA,VAL,INDX,JNDX,NNZ,
!    $            GINT,NINTC,BETA,GX,NCART,WORK,LWORK)
!            DO I1=1,NCART
!               GX(I1)=GCARTORIG(I1)-GX(I1)
!            ENDDO
!C!             print*, 'TESTXC: ', SQRT(DOT_PRODUCT(GX,GX)/NCART)
!C!GY=At(GRAD-BEEt*GINT(K))=At*GX
!            NRHS=1
!            LDB=NCART
!            UPLO='U'
!            N=NCART
!            LDGCS=KD+1
!            CALL DPBTRS(UPLO,N,KD,NRHS,GCS,LDGCS,GX,LDB,INFO)
!            IF (INFO.NE.0) PRINT*,'WARNING - after DPBTRS 2 INFO=',INFO
!C!       CALL DPOTRS(UPLO,N,NRHS,GC,LDA,GX,LDB,INFO)
!            TRANSA=0
!            M=NINTC; N = 1; K=NCART
!            GY(:) = 0.0D0
!            CALL DCOOMM(TRANSA,M,N,K,ALPHA,DESCRA,VAL,INDX,JNDX,NNZ,GX,
!    $            NCART,BETA,GY,NINTC,WORK,LWORK) 

!             print*, 'TESTXI: ', SQRT(DOT_PRODUCT(GY,GY)/NINTC)
!update GINT 
!            SUMGY2=0.0
!            DO I1=1,NINTC
!               GINT(I1)=GINT(I1)+GY(I1)
!               SUMGY2=SUMGY2+GY(I1)*GY(I1)
!            ENDDO
!            MEANSUMGY2=SQRT(SUMGY2/NINTC)
!            IF (MEANSUMGY2.LT.ITGRADCUTOFF) THEN
!               LCONVERGE=.TRUE.
!            ENDIF
!            
!            NITS = NITS+1
!            IF (NITS.GT.1000) THEN
!               print*, 'ERROR: Grad iteration failed to converge!'
!               CALL DUMPCOORDS(XCART,'badcoords.xyz', .FALSE.)
!               CALL DUMPCOORDS(GCARTORIG, 'badgrad.xyz', .FALSE.)
!               STOP
!            ENDIF
!      print *,'transform MEANSUMGY2',MEANSUMGY2
!         ENDDO
!      ENDIF
!C
!move final internal gradients to GRAD
!C
!     DO I1=1,NINTC
!       GRAD(I1)=GINT(I1)
!     ENDDO
!C
!zero GCS2 for next call to getbee
!     DO J1=1,KD
!       DO I1=1,J1
!         GCS2(KD+1+I1-J1,J1)=0.0
!         GCS2(KD+1+J1-I1,I1)=0.0
!       ENDDO
!     ENDDO
!     DO J1=KD+1,NCART
!       DO I1=J1-KD,J1
!         GCS2(KD+1+I1-J1,J1)=0.0
!         GCS2(KD+1+J1-I1,I1)=0.0
!       ENDDO
!     ENDDO
!C
       RETURN
       END
!


! *************************************************************************

!     
!     subroutine to do internal coordinate manipulations using B-matrix
!     calls GETBEE to calcualte B-matrix
!     then transforms a step in internals to a step in cartesians 
!     Nemeth et al. JCP 113,5598(2000)
!     
      SUBROUTINE AMB_TRANSBACKDELTA(X,CARTX,COORDS,NINTC,NCART,NNZ,KD,FAILED,PTEST2,INTEPSILON)
      USE COMMONS, only: NATOMS
      USE MODAMBER9, only: nbona, nbonh, ntheth, ntheta
      USE INTCOMMONS, only : INTNEWT, NATINT, BACKTCUTOFF, BACKTPRINTERR
      USE SPFUNCTS, only : DUMPCOORDS

      IMPLICIT NONE
!     

!     REAL*8 PREVDIH
!     COMMON /DIHINFO/ PREVDIH(MAXP)
!     
!     NNZ is number of non-zero elements in B-matrix
!     NINTC is number of internal coordinates
!     NCART is number of cartesians
!    E 
      INTEGER NNZ,NINTC,NCART
!     
!     COORDS is cartesian coordinates
!     OLDQ is old internal coordinates
!     X is eigenvector in internals, CARTX is ev in cartesians which we need to calculate
!     need to zero CARTX
!     
      DOUBLE PRECISION COORDS(3*NATOMS),DUMCOORDS(3*NATOMS)
      DOUBLE PRECISION VAL(NNZ),WORK(NINTC)
      INTEGER INDX(NNZ),JNDX(NNZ),NITS
      DOUBLE PRECISION ALPHA,BETA
      INTEGER TRANSA,M,N,DESCRA(9),I1,LDGCS,LDB,INFO,NRHS,LWORK,ATOMWIDTH,KD
      DOUBLE PRECISION INTEPSILON,CARTX(NCART),X(NINTC),NEWX(NINTC),XNEW2(NINTC),DUMX(NINTC),SUMDC2,DUMMYQ(NINTC)
      DOUBLE PRECISION MEANDC2,GCS2(2*KD+1,NCART),GCS(KD+1,NCART)
      DOUBLE PRECISION CARTDOTX,CARTDOTY,CARTDOTZ
      CHARACTER*1 UPLO
      LOGICAL NOCOOR,LCONVERGE,FAILED,PTEST, NODERV
      INTEGER J1, J2, K
!    EFK
      INTEGER DUMINDX(NNZ), DUMJNDX(NNZ)
      DOUBLE PRECISION DUMVAL(NNZ), DUMGCS2(2*KD+1,NCART), DUMXINT1(NINTC), DUMXINT2(NINTC)
      DOUBLE PRECISION QORIG(NINTC), TESTQ(NINTC), QTARGET(NINTC)
      LOGICAL PTEST2


!      INTEGER CENTERS, RINGS, FRINGS, LINDIH, IMPDIH, NCNT, NRNG, NFRG,
!     $     NLDH, NIMP, NBDS, ATOMxPOSINC
!      REAL*8 COEFF
!      COMMON /NATINTERN/ COEFF(30,6,6),
!     $     CENTERS(MAXA, 0:5), RINGS(MAXA, 6),
!     $     FRINGS(MAXA, 6), LINDIH(MAXA, -1:8), IMPDIH(MAXA, 4), NCNT,
!     $     NRNG, NFRG, NLDH, NIMP, NBDS, ATOMxPOSINC(10,5,5)
      
      CHARACTER*30 FNAME
!msb50
      INTEGER NBOND, NTHETS, LUNIT, GETUNIT
      DOUBLE PRECISION PI, TWOPI
      PARAMETER(PI=3.141592653589793D0)


      TWOPI = 2.0D0*PI
      NBOND = nbona+nbonh 
      NTHETS = ntheta +ntheth

!      PTEST=.TRUE.
      PTEST = .FALSE.

!     make dummy array NEWX that can change (error in estimate of X)
!     and DUMX which is cumulative estimate of X
!     and dummy array DUMCOORDS that can change
!     
      DO I1=1,NINTC
         NEWX(I1)=X(I1)
         DUMX(I1)=0.0D0
      ENDDO
      DO I1=1,NCART
         CARTX(I1)=0.0D0
         DUMCOORDS(I1)=COORDS(I1)
      ENDDO
!    
!    start iterative loop to do conversion of internal coordinates to cartesians
!    have X, need to satisfy CARTX=A*X and X=B*CARTX
!    
      LCONVERGE=.FALSE.
      FAILED=.TRUE.
      NITS=0
!    
!    Attempt to prevent NaN on Solaris
!    
      DUMMYQ(1:NINTC) = 0.0D0
      XNEW2(1:NINTC)=0.0D0

      DO WHILE (.NOT.LCONVERGE)
!        IF (PTEST) print '(A,I6,G20.10)',
!    $        'transbackdelta>> internals error', NITS,
!    $        SQRT(DOT_PRODUCT(NEWX(1:NINTC),NEWX(1:NINTC)))
!    
!    get B-matrix without changing cartesians to internals
!    added DUMMYQ argument 13/2/04 DJW
!    
         IF (PTEST) THEN
            PRINT*,'DUMCOORDS BEFORE GETBEE, INTEPSILON=',INTEPSILON
            WRITE(*,'(6G18.10)') (DUMCOORDS(J1),J1=1,NCART)
            PRINT*,'DUMMYQ before GETBEE'
            WRITE(*,'(6G18.10)') (DUMMYQ(J1),J1=1,NINTC)
         ENDIF

         IF (NITS.EQ.0) THEN ! store original and target internals
            NOCOOR = .FALSE.; NODERV = .FALSE.
            CALL AMB_GETBEENAT(DUMCOORDS,QORIG,VAL,INDX,JNDX,NINTC,NCART,NNZ,GCS2,NOCOOR,NODERV,KD)
            QTARGET(:) = QORIG(:) + X(:)
         ELSE
            NOCOOR = .TRUE.; NODERV = .FALSE.
            CALL AMB_GETBEENAT(DUMCOORDS,DUMMYQ,VAL,INDX,JNDX,NINTC,NCART,NNZ,GCS2,NOCOOR,NODERV,KD)
         ENDIF

        IF (PTEST) THEN
           PRINT*,'GCS2 AFTER GETBEE:',INTEPSILON
           WRITE(*,'(6G18.10)') ((GCS2(J1,J2),J1=1,2*KD+1),J2=1,NCART)
        ENDIF

!    
!    want to find CARTX=A*X
!    A = G^-1*BEEt
!    same as solving G*CARTX=BEEt*X
!    G=GC+INTEPSILON*I
!    GC=BEEtBEE (as in TRANSFORM), stored as sparse GCS
!    KD is number of superdiagonals of GC
!    see dpbtrf.f http://www.netlib.org/lapack/double/dpbtrf.f
!    
!    KD=3*ATOMWIDTH+2       
!    need to move GCS2 into GCS

         GCS(:,:) = 0.0D0
         DO J1=1,KD
            DO I1=1,J1
               GCS(KD+1+I1-J1,J1)=GCS2(KD+1+I1-J1,J1)
            ENDDO
         ENDDO
         DO J1=KD+1,NCART
            DO I1=J1-KD,J1
               GCS(KD+1+I1-J1,J1)=GCS2(KD+1+I1-J1,J1)
               !print*, GCS(KD+1+I1-J1,J1)
            ENDDO
         ENDDO
           
!    
!    form G = GC + eI where e is small
!    
         DO I1=1,NCART
            GCS(KD+1,I1)=GCS(KD+1,I1)+INTEPSILON
         ENDDO
!    
!    NEWX is error in internals at end of last cycle
!    multiply NEWX by BEEt to give CARTX
!    Uses sparse matrix multiplication routine from BLAS
!    http://math.nist.gov/~KRemington/fspblas/
!    
         TRANSA=1
         M=NINTC
         N=1
         K=NCART
         ALPHA=1.0D0
         BETA=0.0D0
         DO I1=1,9
            DESCRA(I1)=0
         ENDDO
         DESCRA(4)=1
         LWORK=NINTC
         IF (PTEST) THEN
            PRINT*,'before dcoomm 1 NEWX:'
            WRITE(*,'(6G20.10)') (NEWX(J1),J1=1,NINTC)
            PRINT*,'before dcoomm 1 CARTX:'
            WRITE(*,'(6G20.10)') (CARTX(J1),J1=1,NCART)
         ENDIF
         CALL DCOOMM(TRANSA,M,N,K,ALPHA,DESCRA,VAL,INDX,JNDX,NNZ, &
     &        NEWX,NINTC,BETA,CARTX,NCART,WORK,LWORK)
         IF (PTEST) THEN
            PRINT*,'after dcoomm 1 NEWX:'
            WRITE(*,'(6G20.10)') (NEWX(J1),J1=1,NINTC)
            PRINT*,'after dcoomm 1 CARTX:'
            WRITE(*,'(6G20.10)') (CARTX(J1),J1=1,NCART)
         ENDIF         
!    
!    Cholesky decompose GCS
!    use lapack routine DPBTRF
!    LDGCS is leading dimension of GCS
!    
         UPLO='U'
         N=NCART
         LDGCS=KD+1
        IF (PTEST) THEN
           PRINT*,'GCS before DPBTRF:'
           WRITE(*,'(6G18.10)') ((GCS(J1,J2),J1=1,KD+1),J2=1,NCART)
        ENDIF
!    CALL MYCPU_TIME(TESTTIME1, .FALSE.)
         CALL DPBTRF(UPLO,N,KD,GCS,LDGCS,INFO)
!    CALL MYCPU_TIME(TESTTIME2, .FALSE.)
!    TOTCHOLTIME = TOTCHOLTIME + TESTTIME2 - TESTTIME1
         IF (INFO.NE.0) PRINT*,'after DPBTRF INFO=',INFO

!msb50
!        IF (PTEST) THEN
!           PRINT*,'GCS after DPBTRF:'
!           WRITE(*,'(6G18.10)') ((GCS(J1,J2),J1=1,KD+1),J2=1,NCART)
!        ENDIF
!    
!    GCS is now the upper triangular matrix after decomposition
!    SOLVE G*DELTAX=CARTX FOR DELTAX, WHERE G IS ORIGINAL BEETBEE+INTEPSILON*I
!    done using lapack routine DPBTRS with input of the upper matrix from
!    cholesky decomposition
!    
         NRHS=1
         LDB=NCART
         CALL DPBTRS(UPLO,N,KD,NRHS,GCS,LDGCS,CARTX,LDB,INFO)
         IF (INFO.NE.0) PRINT*,'after DPBTRS INFO=',INFO
        IF (PTEST) THEN
           PRINT*,'GCS after DPBTRS:'
           WRITE(*,'(6G18.10)') ((GCS(J1,J2),J1=1,KD+1),J2=1,NCART)
           PRINT*,'after DPBTRS CARTX:'
           WRITE(*,'(6G20.10)') (CARTX(J1),J1=1,NCART)
        ENDIF

!    
!    CARTX is now solution of CARTX=A*X
!    evaluate XNEW2=B*CARTX
!    

         IF (INTNEWT) THEN
            ! evaluate XNEW2 directly from new cartesian coords
            NOCOOR = .FALSE.; NODERV = .TRUE.
            CALL AMB_GETBEENAT(DUMCOORDS+CARTX,DUMXINT2,DUMVAL,DUMINDX, &
     &           DUMJNDX,NINTC,NCART,NNZ,DUMGCS2,NOCOOR,NODERV,KD)

            ! align dihedral angles if using primitive internals
            IF (.NOT.NATINT) THEN 
               !msb 50 - probably bug for amber
               PRINT*, "Warning - AMB_TRANSBACKDELTA has natint .FALSE.&
               &- are you sure? -> uncomment stop in ambnatinern_extras"
               STOP 
               DO I1 = NBOND+NTHETS+1, NINTC
                  IF (DUMXINT2(I1) - QTARGET(I1).GT.PI) THEN
                     DUMXINT2(I1) = DUMXINT2(I1) - TWOPI
                  ELSE IF (DUMXINT2(I1) - QTARGET(I1).LT.-PI) THEN
                     DUMXINT2(I1) = DUMXINT2(I1) + TWOPI
                  ENDIF                  
               ENDDO
            ENDIF
            NEWX = QTARGET - DUMXINT2

         ELSE ! solve linear equation only; iterate away effect of epsilon
            TRANSA=0
            N=1
            CALL DCOOMM(TRANSA,M,N,K,ALPHA,DESCRA,VAL,INDX,JNDX,NNZ, &
     &           CARTX,NCART,BETA,XNEW2,NINTC,WORK,LWORK)

            IF (PTEST) THEN
               PRINT*,'XNEW2 after DCOOMM 2:'
               WRITE(*,'(6G18.10)') (XNEW2(J1),J1=1,NINTC)
               PRINT*,'after DCOOMM 2 CARTX:'
               WRITE(*,'(6G20.10)') (CARTX(J1),J1=1,NCART)
            ENDIF
            DO I1=1,NINTC
               DUMX(I1)=DUMX(I1)+XNEW2(I1) ! sum of f(x_i) instead of f(sum(x_i))
               NEWX(I1)=X(I1)-DUMX(I1)
            ENDDO
         ENDIF

!    project translation out of CARTX
         CARTDOTX=0.0D0
         CARTDOTY=0.0D0
         CARTDOTZ=0.0D0
         DO I1=1,NATOMS
            CARTDOTX=CARTDOTX+CARTX(3*(I1-1)+1)
            CARTDOTY=CARTDOTY+CARTX(3*(I1-1)+2)
            CARTDOTZ=CARTDOTZ+CARTX(3*(I1-1)+3)
         ENDDO
         DO I1=1,NATOMS
            CARTX(3*(I1-1)+1)=CARTX(3*(I1-1)+1)-CARTDOTX/NATOMS
            CARTX(3*(I1-1)+2)=CARTX(3*(I1-1)+2)-CARTDOTY/NATOMS
            CARTX(3*(I1-1)+3)=CARTX(3*(I1-1)+3)-CARTDOTZ/NATOMS
         ENDDO
!    update DUMCOORDS for next call to GETBEE; check for convergence
         SUMDC2=0.0D0
         DO I1=1,NCART
            DUMCOORDS(I1)=DUMCOORDS(I1)+CARTX(I1)
            SUMDC2=SUMDC2+CARTX(I1)*CARTX(I1)
         ENDDO
         MEANDC2=SUMDC2/NINTC

         IF (MEANDC2.LT.BACKTCUTOFF) THEN                  
            LCONVERGE=.TRUE.
            FAILED=.FALSE.
         ENDIF
         NITS=NITS+1

         IF (PTEST2) PRINT '(A,I5,2G20.10)',&
             'transbackdelta>> NITS,MEANDC2, int erro =',NITS, MEANDC2,&
             SQRT(DOT_PRODUCT(NEWX(1:NINTC),NEWX(1:NINTC)))

         IF (SQRT(DOT_PRODUCT(NEWX(1:NINTC),NEWX(1:NINTC))).GT.100) THEN
            print*, 'ERROR! deltatransform diverged.'
            CALL DUMPCOORDS(COORDS,'diverged.xyz', .FALSE.)
            LUNIT=GETUNIT()
            OPEN(unit=LUNIT,FILE='badstp', STATUS = 'UNKNOWN')
            WRITE(LUNIT,'(G25.15)') (X(J1), J1 = 1,NINTC)
            CLOSE(LUNIT)
            LCONVERGE = .TRUE.
            FAILED = .TRUE.
            CARTX(:) = 0.0D0
            RETURN
         ENDIF

         IF (NITS.GT.100) THEN
            PRINT*,'WARNING internal coordinate transformation&
               & did not converge, MEANDC2=',MEANDC2
            FAILED = .TRUE.
            LCONVERGE=.TRUE.
         ENDIF
!    
      ENDDO

!    update CARTX so that it is the difference between final DUMCOORDS and COORDS
!    print *,'tdb here1'

      DO I1=1,3*NATOMS
         CARTX(I1)=DUMCOORDS(I1)-COORDS(I1)
      ENDDO
!    
!    project translation out of CARTX
      CARTDOTX=0.0D0
      CARTDOTY=0.0D0
      CARTDOTZ=0.0D0
      DO I1=1,NATOMS
         CARTDOTX=CARTDOTX+CARTX(3*(I1-1)+1)
         CARTDOTY=CARTDOTY+CARTX(3*(I1-1)+2)
         CARTDOTZ=CARTDOTZ+CARTX(3*(I1-1)+3)
      ENDDO
      DO I1=1,NATOMS
         CARTX(3*(I1-1)+1)=CARTX(3*(I1-1)+1)-CARTDOTX/NATOMS
         CARTX(3*(I1-1)+2)=CARTX(3*(I1-1)+2)-CARTDOTY/NATOMS
         CARTX(3*(I1-1)+3)=CARTX(3*(I1-1)+3)-CARTDOTZ/NATOMS
      ENDDO

      IF(BACKTPRINTERR) print*, 'transbackdelta>>> error = ', SQRT(DOT_PRODUCT(NEWX(1:NINTC),NEWX(1:NINTC)))

      RETURN
      END


      SUBROUTINE GETDIHONLY(XCART)
!subroutine from AMBGETBEENAT -- only dihedral parts to provide easier access
      USE COMMONS, ONLY: NATOMS
      USE INTCOMMONS
      USE modamber9, only: ih, m04

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN):: XCART(3*NATOMS)     
      INTEGER DIHC, RATMS, NANGC, NDIHC
      INTEGER XIND(8),INFO(-1:8)
      DOUBLE PRECISION PHI, dPdI(0:2), dPdJ(0:2), dPdK(0:2), dPdL(0:2)  
      DOUBLE PRECISION NANGCR
      INTEGER RNG, FRG,LDH, I1, TCOUNT, J1, K1,L1
      INTEGER CURPHI
      INTEGER Ix, Jx, Kx, Lx
      DOUBLE PRECISION PHI2, dPdI2(0:2), dPdJ2(0:2), dPdK2(0:2), dPdL2(0:2)
      LOGICAL NOCOOR
      DOUBLE PRECISION PHI1, PREVPHI1
      INTEGER ATOMS2PERM

      DIHC=0 
      NOCOOR = .FALSE.
      DO 25 RNG=1,NRNG

         IF(RINGS(RNG,6).EQ.-1) THEN
            RATMS = 5 ! 5 atoms in ring
            NANGC = 2; NDIHC = 2 ! 2 angle & 2 dihedral coords
         ELSE
            RATMS = 6 ! 6 atoms in ring
            NANGC = 3; NDIHC = 3
         ENDIF

         DO I1=1,RATMS
            XIND(I1) = 3*(RINGS(RNG,I1) - 1) + 1
         ENDDO

         TCOUNT = 1
         DO I1 = 1,RATMS
            J1 = MOD(I1,RATMS)+1
            K1 = MOD(I1+1,RATMS)+1
            L1 = MOD(I1+2,RATMS)+1

            CALL AMBGETTORSION(XCART(XIND(I1):XIND(I1)+2),                 &
     &           XCART(XIND(J1):XIND(J1)+2), XCART(XIND(K1):XIND(K1)+2),&
     &           XCART(XIND(L1):XIND(L1)+2), PHI, dPdI, dPdJ, dPdK, dPdL,&
     &           .FALSE.,.TRUE.)

            DIHC = DIHC + 1

            IF (TCOUNT.EQ.1) PREVPHI1 = PREVDIH(DIHC)
            IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,TCOUNT)
            IF (TCOUNT.EQ.1) PHI1 = PHI
            TCOUNT = TCOUNT + 1
!           PRINT '(a5,i4,4a4, f11.6,f11.6)',"RING", DIHC,ih(m04+RINGS(RNG,I1)-1), ih(m04+RINGS(RNG,J1)-1), &
! &                  ih(m04+RINGS(RNG,K1)-1), ih(m04+RINGS(RNG,L1)-1), PREVDIH(DIHC),(PREVDIH(DIHC)/3.14159265358979)*180.0
          ENDDO

 25   CONTINUE

      DO FRG = 1,NFRG
          DO I1=1,6
             XIND(I1) = 3*(FRINGS(FRG,I1) - 1)+1
          ENDDO
    
          CALL AMBGETTORSION(XCART(XIND(3):XIND(3)+2),&
     &         XCART(XIND(1):XIND(1)+2), XCART(XIND(2):XIND(2)+2), &
     &         XCART(XIND(6):XIND(6)+2), PHI, dPdI, dPdJ, dPdK, dPdL,&
     &         NOCOOR,.TRUE.)   
          DIHC = DIHC + 1
          PREVPHI1 = PREVDIH(DIHC)
          IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, 0.0D0, 0.0D0, 1) 
!         PRINT '(a5,i4,4a4,f11.6,f11.6)',"FRNG", DIHC,ih(m04+RINGS(RNG,3)-1), ih(m04+RINGS(RNG,1)-1), 
!                ih(m04+RINGS(RNG,2)-1), ih(m04+RINGS(RNG,6)-1), PREVDIH(DIHC),  (PREVDIH(DIHC)/3.14159265358979)*180.0
 
          CALL AMBGETTORSION(XCART(XIND(5):XIND(5)+2),&
     &         XCART(XIND(1):XIND(1)+2), XCART(XIND(2):XIND(2)+2),  &
     &         XCART(XIND(4):XIND(4)+2), PHI2, dPdI2, dPdJ2, dPdK2, &
     &         dPdL2, NOCOOR,.TRUE.)
         DIHC = DIHC + 1

         IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI2, DIHC, PHI2-PHI, PREVDIH(DIHC)-PREVPHI1, 2)         
!           PRINT '(a5,i4,4a4,f11.6,f11.6)',"FRNG", DIHC,ih(m04+RINGS(RNG,5)-1), ih(m04+RINGS(RNG,1)-1), 
!                   ih(m04+RINGS(RNG,2)-1), ih(m04+RINGS(RNG,4)-1), PREVDIH(DIHC), (PREVDIH(DIHC)/3.14159265358979)*180.0
         TCOUNT = TCOUNT + 1
      ENDDO 

      DO 30 LDH=1,NLDH
         INFO(-1:8) = LINDIH(LDH,-1:8)
         !INFO(-1) and INFO(0) denote the number of possible I and L atoms
         DO I1=1,INFO(0)+INFO(-1)+2
            XIND(I1) = 3*(INFO(I1) - 1) + 1
         ENDDO
         NANGC = INFO(-1)*INFO(0) ! number of torsions involved

         NANGCR = 1.0D0/NANGC

         CURPHI = 0
         DO I1=3,INFO(-1)+2
            DO L1=INFO(-1)+3,INFO(-1)+INFO(0)+2
               CURPHI = CURPHI + 1

               Ix = XIND(I1)
               Jx = XIND(1)
               Kx = XIND(2)
               Lx = XIND(L1)
               
!               PRINT*, "getdihonly phi1", PHI1
!               IF (CURPHI.EQ.1) THEN
!                  DO I2 = 1, NPERMGROUP
!                     IF (PERMNEIGHBOURS(I2,4)==LDH) THEN
!                        PERMNEIGHBOURS(I2, 8)==PHI1
!                        I3= PERMNEIGHBOURS(I2,1)
!                        PERMNEIGHBOURS(I3,
!                      IF (PERMNEIGHBOURS(I2,5)==LDH) THEN
               CALL AMBGETTORSION(XCART(Ix:Ix+2), XCART(Jx:Jx+2),  &
     &              XCART(Kx:Kx+2), XCART(Lx:Lx+2), PHI, dPdI, dPdJ,   &
     &              dPdK, dPdL, NOCOOR,.TRUE.)                              

               DIHC = DIHC + 1
               IF (CURPHI.EQ.1) PREVPHI1 = PREVDIH(DIHC)

               IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,CURPHI)
!              PRINT '(a5,i4,x,4a4,f11.6, f11.6)',"LDH", DIHC,ih(m04+INFO(I1)-1), ih(m04+INFO(1)-1), 
!    &                 ih(m04+INFO(2)-1), ih(m04+INFO(L1)-1), PREVDIH(DIHC) , (PREVDIH(DIHC)/3.14159265358979)*180.0

               IF (CURPHI.EQ.1) PHI1 = PHI
            ENDDO
         ENDDO
 30   CONTINUE 

      END SUBROUTINE GETDIHONLY
! *****************************************
      SUBROUTINE GETLDH(LDH, STDIHC,PHI1,XCART)
         use intcommons, only: LINDIH, PREVDIH
         use modamber9, only: ih, m04
         use commons, only: NATOMS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: LDH, STDIHC
         DOUBLE PRECISION, INTENT(INOUT):: PHI1 
         DOUBLE PRECISION, INTENT(IN):: XCART(3*NATOMS)
         INTEGER :: INFO(-1:8), XIND(8), NANGC, CURPHI, Ix,Jx,Kx,Lx,DIHC
         DOUBLE PRECISION :: NANGCR, PHI, dPdI(3), dPdJ(3), dPdK(3), dPdL(3)
         DOUBLE PRECISION :: PREVPHI1
         INTEGER :: I1,L1
         DIHC = STDIHC -1 
!find ph1
         INFO(-1:8) = LINDIH(LDH,-1:8)
         INFO(-1:8) = LINDIH(LDH,-1:8)
         !INFO(-1) and INFO(0) denote the number of possible I and L atoms
         DO I1=1,INFO(0)+INFO(-1)+2
            XIND(I1) = 3*(INFO(I1) - 1) + 1
         ENDDO
         NANGC = INFO(-1)*INFO(0) ! number of torsions involved

         NANGCR = 1.0D0/NANGC

         CURPHI = 0
         DO I1=3,INFO(-1)+2
            DO L1=INFO(-1)+3,INFO(-1)+INFO(0)+2
               CURPHI = CURPHI + 1

               Ix = XIND(I1)
               Jx = XIND(1)
               Kx = XIND(2)
               Lx = XIND(L1)

               CALL AMBGETTORSION(XCART(Ix:Ix+2), XCART(Jx:Jx+2),  &
     &              XCART(Kx:Kx+2), XCART(Lx:Lx+2), PHI, dPdI, dPdJ,   &
     &              dPdK, dPdL, .FALSE., .TRUE.)                      
               DIHC = DIHC + 1
               IF (CURPHI.EQ.1) PREVPHI1 = PREVDIH(DIHC)

!               print '(a15,3i4,4(i4,a4),f8.2)','GETBEE>dih', & 
!     &       LDH,DIHC, INFO(I1), ih(m04+INFO(I1)-1),INFO(1),ih(m04+INFO(1)-1),&
!              INFO(2), ih(m04+INFO(2)-1),INFO(L1),ih(m04+INFO(L1)-1), PHI
               CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,CURPHI)

               IF (CURPHI.EQ.1) PHI1 = PHI
            ENDDO
         ENDDO
      END SUBROUTINE GETLDH
        



! ******************************************
      SUBROUTINE FINDPERMDIH
!msb50 - find dihedral angles which involve atoms from different permgroups
!        s.t. for INTMINPERM the minimum dihedral can be found
      use key, only: NPERMGROUP
      use intcommons
      use modamber9, only: ih, m04, nres, ix,i02,m02
      use commons, only: NATOMS,DEBUG     

      IMPLICIT NONE
      INTEGER i1, A1, A2, B, k1, j1,l1, B1P,h2,hh, Bh,B1, NCTa
      INTEGER :: PGH, PGH2, PGBh
      CHARACTER(LEN=4) :: RESNAME,  B1PN,hname,BhN
      CHARACTER(LEN=4) :: permac(10,3)
      INTEGER :: RESNUM
      INTEGER ATOMS2PERM(NATOMS) !includes pg of some (important atoms, rest 0)
      INTEGER ATCENTRE2PERM(NATOMS) !pg for certain centres (i.e. the atoms around it!!)
      INTEGER thisroundfound
      LOGICAL critical, arginine
      INTEGER pg_centre, pg_centre2, pg_centreAC
      INTEGER alreadypg(NATOMS)
      data permac(1,:) /'LEU ', 'NLEU','CLEU'/
      data permac(2,:) /'ILE ', 'NILE', 'CILE'/
      data permac(3,:) /'LYS ', 'NLYS', 'CLYS'/
      data permac(4,:) /'PHE ', 'NPHE', 'CPHE'/
      data permac(5,:) /'TYR ', 'NTYR', 'CTYR'/ 
      data permac(6,:) /'GLU ', 'NGLU', 'CGLU'/
      data permac(7,:) /'ASP ', 'NASP', 'CASP'/
      data permac(8,:) /'MET', 'NMET', 'CMET'/
      data permac(9,:) /'GLN', 'NGLN', 'CGLN'/
      data permac(10,:) /'ARG', 'NARG','CARG'/
      CHARACTER(LEN=2) :: bondedHs(4)
      data bondedHs /'HB','HE', 'HD', 'HG'/
     

      IF (.NOT.ALLOCATED(MBADJACENT)) THEN
          CALL AMB_NATINTSETUP
          CALL AMBGETNATINTERN
      ENDIF
      IF (.NOT.ALLOCATED(PERMNANGLE)) ALLOCATE(PERMNANGLE(NPERMGROUP,2))
      ATOMS2PERM(:)=0
      PGH = 0; PGH2 = 0; PGBh = 0
      IF (.NOT.ALLOCATED(PERMNEIGHBOURS)) ALLOCATE(PERMNEIGHBOURS(NPERMGROUP,7))
      PERMNEIGHBOURS(:,:) = 0
      alreadypg(:)=0
      !alreadypg(ATOM)=ATOM2 if this centre has already been found and is in permneighbours
      IF (DEBUG) PRINT*, "ncnt2", NCNT2
      !this scheme only works if of the two adjacent permgroups at least one is centre2 - 
      ! which it is for the amino acids i have
      arginine = .FALSE.
      DO I1 = 1, NCNT2
         PGH = 0; PGH2 = 0; PGBh = 0
         hh = 0; h2 = 0; Bh=0
         pg_centre2=0;pg_centreAC=0;pg_centre=0
         thisroundfound=0
         NCTa = CENTER2(I1)
         A1 = CENTERS(NCTa, 1)
         CALL GETRESID(A1, resname, resnum)
         critical = .FALSE.
         DO k1=1,10
            IF (permac(k1,1).EQ.resname.OR.permac(k1,2).EQ.resname &
     &        .OR.permac(k1,3).EQ.resname) THEN
               critical = .TRUE.
               IF (k1==10) arginine=.TRUE.
               EXIT
            ENDIF
         ENDDO
         !PRINT*,I1, "findminpermn, neigbours of", A1, "critical", critical, resname
         IF (.NOT.critical) CYCLE !no two groups after one another
         !need to find hydrogens/permgroup atoms to find permgroup
         IF (ih(m04+A1-1)=="CA") CYCLE !happen to be center2s as well
         DO j1 = 2, 5 !loop over atoms bonded to A1
            B1 = CENTERS(NCTa,j1)
            DO k1= 1, NCNT2 !center 2
               IF(CENTERS(CENTER2(k1),1)==B1.AND.alreadypg(A1).NE.B1) THEN !found 2-2
                  DO l1 = 5, 2,-1
                    pg_centre2=B1
                    hname=ih(m04+CENTERS(CENTER2(k1),l1)-1)
                    IF ((hname.EQ.('HB2').OR.hname.EQ.('HG2')).OR. &
     &              (hname.EQ.('HD2').OR.hname.EQ.('HE23')).OR.hname.EQ.('HE2 ')) THEN
                      h2 = CENTERS(CENTER2(k1),l1)
                      thisroundfound = thisroundfound + 1
                      alreadypg(B1)=A1
                      !PRINT*, "neighbour", B1, "h", h2
                      !PRINT*, "alreadypg", alreadypg(A1), "b1", b1,alreadypg(B1)
                      EXIT
                    ENDIF
                  ENDDO !once a particular neighbour with centertype 2 is found 
                  EXIT  !exit loop over 2centres in that case
               ENDIF
            ENDDO 
            DO k1 = 1,3  !all other centers
               IF ((resname.EQ.permac(1,k1).AND.ih(m04+B1-1).EQ.'CG  ').OR.&
     &             (ih(m04+B1-1).EQ.'CG  '.AND.resname.EQ.permac(4,k1)).OR.&
     &             (ih(m04+B1-1).EQ.'CG  '.AND.resname.EQ.permac(5,k1)))THEN
                  pg_centreAC=B1
                  DO l1 = 1, MBBONDNUM(B1)
                   BhN = ih(m04+MBADJACENT(B1,l1)-1)
                   IF (BhN.EQ.'CD1 '.OR.BhN.EQ.'CD2 ') THEN
!                     CALL FINDPERMGROUP(Bh, PGBh)
                      Bh=MBADJACENT(B1,l1)
                      thisroundfound = thisroundfound +1
                      EXIT                     
                    ENDIF 
                  ENDDO   
                  EXIT               
               ELSEIF (ih(m04+B1-1).EQ.'CD1 '.AND.resname.EQ.permac(2,k1))THEN
                   pg_centreAC=B1
                  DO l1 = 1, MBBONDNUM(B1)    
                   BhN = ih(m04+MBADJACENT(B1,l1)-1)
                   IF (BhN.EQ.'HD11 '.OR.BhN.EQ.'HD22'.OR.&
     &                     BhN.EQ.'HD33') THEN
                      Bh=MBADJACENT(B1,l1)
                      thisroundfound = thisroundfound +1
                      EXIT
                    ENDIF
                  ENDDO                  
                  EXIT
               ELSEIF (ih(m04+B1-1).EQ.'NZ  '.AND.resname.EQ.permac(3,k1))THEN
                  pg_centreAC=B1
                  DO l1 = 1, MBBONDNUM(B1)
                   BhN = ih(m04+MBADJACENT(B1,l1)-1)
                   IF (BhN.EQ.'HZ1 '.OR.BhN.EQ.'HZ2 '.OR.&
     &                     BhN.EQ.'HZ3 ') THEN
                      Bh=MBADJACENT(B1,l1)
                      thisroundfound = thisroundfound +1
                      EXIT
                    ENDIF
                  ENDDO                  
                  EXIT
!               ELSEIF (ih(m04+B1-1).EQ.'CG  '.AND.resname.EQ.permac(4,k1)).OR. &
!    &               (ih(m04+B1-1).EQ.'CG  '.AND.resname.EQ.permac(5,k1)) THEN

               ELSEIF(((ih(m04+B1-1).EQ.'CD  '.AND.resname.EQ.permac(6,k1))).OR. &
     &              (ih(m04+B1-1).EQ.'CG  '.AND.resname.EQ.permac(7,k1))) THEN
                  pg_centreAC=B1
                  DO l1 = 1, MBBONDNUM(B1)
                   BhN = ih(m04+MBADJACENT(B1,l1)-1)
                   IF (BhN.EQ.'OE1 '.OR.BhN.EQ.'OE2 ' .OR.&
     &                  BhN.EQ.'OD1 '.OR.BhN.EQ.'OD2 ') THEN
                      thisroundfound = thisroundfound +1
                      Bh=MBADJACENT(B1,l1)
                      EXIT
                    ENDIF
                  ENDDO
                  EXIT
               ENDIF
            ENDDO

            IF ((ih(m04+b1-1).EQ.('HB2 ').OR.ih(m04+b1-1).EQ.('HG2 ')).OR. &
     &        (ih(m04+b1-1).EQ.('HD2 ').OR.ih(m04+b1-1).EQ.('HE23').OR.    &
     &         ih(m04+b1-1).EQ.('HE2 ').OR.ih(m04+b1-1).EQ.('HG12'))) THEN
               hh = b1
               pg_centre=A1
!               CALL FINDPERMGROUP(j1, PGH)  
               thisroundfound = thisroundfound +1
               EXIT !exit whole loop - H will be found last
            ENDIF
         ENDDO 

         IF ((thisroundfound.GT.3).OR.(thisroundfound.LE.1)) THEN
            IF ((alreadypg(A1).NE.0).AND.(thisroundfound.EQ.1)) THEN !could be if we've already included this
               PRINT*, "cycle"
               CYCLE
            ELSE
            PRINT*, "impossible numbers of neighbours found for",A1, thisroundfound
            PRINT*, "ERROR in findpermdih"
            STOP
            ENDIF
         ENDIF
         IF (j1==0) THEN
           PRINT*, "no hydrogen found"; STOP
         ELSE
            CALL FINDPERMGROUP(hh,h2,Bh,PGH,PGH2,PGBh)
            !PRINT*, "after findpg hh,h2,Bh,Pgh,pgh2,pgbh", hh,h2,Bh,PGH,PGH2,PGBh
            IF (PGBh.NE.0) THEN
               PERMNEIGHBOURS(PGH,1) = PERMNEIGHBOURS(PGH,1)+1
               PERMNEIGHBOURS(PGBh,1) = PERMNEIGHBOURS(PGBh,1)+1
               IF (PERMNEIGHBOURS(PGH, PERMNEIGHBOURS(PGH,1)+1).NE.0 .OR. &
     &              PERMNEIGHBOURS(PGBh,PERMNEIGHBOURS(PGBh,1)+1).NE.0) THEN
                   PRINT*,i1, "ERROR", PGH, PERMNEIGHBOURS(PGH,1),PGBh
                   STOP
               ENDIF
               PERMNEIGHBOURS(PGH,PERMNEIGHBOURS(PGH,1)+1) = PGBh
               PERMNEIGHBOURS(PGBh,PERMNEIGHBOURS(PGBh,1)+1) = PGH
               ATOMS2PERM(hh)=PGH
               ATOMS2PERM(Bh)=PGBH
               ATCENTRE2PERM(pg_centre) = PGH
               ATCENTRE2PERM(pg_centreAC)= PGBH
            ENDIF
            IF (PGH2.NE.0) THEN !not bonded to anything new
               PERMNEIGHBOURS(PGH,1) = PERMNEIGHBOURS(PGH,1)+1
               PERMNEIGHBOURS(PGH2,1) = PERMNEIGHBOURS(PGH2,1)+1
               IF( PERMNEIGHBOURS(PGH, PERMNEIGHBOURS(PGH,1)+1).NE.0 .OR. &
     &              PERMNEIGHBOURS(PGH2,PERMNEIGHBOURS(PGH2,1)+1).NE.0) THEN
               PRINT*,i1, "ERROR", PGH, PERMNEIGHBOURS(PGH,1),PGH2
               STOP
               ENDIF
               PERMNEIGHBOURS(PGH, PERMNEIGHBOURS(PGH,1)+1) = PGH2
               PERMNEIGHBOURS(PGH2,PERMNEIGHBOURS(PGH2,1)+1) = PGH
               ATOMS2PERM(hh)=PGH
               ATOMS2PERM(h2)=PGH2
               ATCENTRE2PERM(pg_centre) = PGH
               ATCENTRE2PERM(pg_centre2)= PGH2
            ENDIF
         ENDIF

      ENDDO


      !one special permgroup for arginine, not involving a centre 2
      IF (arginine) THEN
         DO I1=1,NRES
               h2=0;Bh=0;hh=0; pg_centre2=0;pg_centre=0; pg_centreAC=0
               IF (ih(m02+I1-1)=='ARG'.OR.ih(m02+I1-1)=='NARG' &
     &                 .OR.ih(m02+I1-1)=='CARG') THEN
                   DO k1=ix(i02+I1-1),ix(i02+I1)-1
                      SELECT CASE (ih(m04+k1-1))
                         CASE ('NH1 ')
                              pg_centre2=k1
                         CASE ('CZ  ')
                              pg_centre=k1
                         CASE ('NH2 ')
                              pg_centreAC=k1; hh=k1
                         CASE ('HH11')
                              h2=k1
                         CASE ('HH21')
                              Bh=k1
                      END SELECT 
                   ENDDO
                   CALL FINDPG_NOSWAPS(hh,h2,Bh,PGH,PGH2,PGBh)
            !PRINT*, "after findpg hh,h2,Bh,Pgh,pgh2,pgbh", hh,h2,Bh,PGH,PGH2,PGBh
                   IF (PERMNEIGHBOURS(PGH,1).NE.0) THEN
                      PRINT*, "error, only 1 permgroup about", PGH
                      STOP
                   ENDIF
                   PERMNEIGHBOURS(PGH,1)= 2
                   PERMNEIGHBOURS(PGH2,1)= 1
                   PERMNEIGHBOURS(PGBh,1)= 1
                   PERMNEIGHBOURS(PGH,2) = PGBh;PERMNEIGHBOURS(PGH,3)=PGH2
                   PERMNEIGHBOURS(PGBh,2) = PGH  
                   PERMNEIGHBOURS(PGH2,2)= PGH
                   ATOMS2PERM(hh)=PGH;ATOMS2PERM(h2)=PGH2;ATOMS2PERM(Bh)=PGBh
                   ATCENTRE2PERM(pg_centre) = PGH
                   ATCENTRE2PERM(pg_centre2)= PGH2 
                   ATCENTRE2PERM(pg_centreAC)=PGBh
               ENDIF
         ENDDO
      ENDIF
         
      IF (DEBUG) THEN
         DO I1=1, NPERMGROUP
            PRINT '(a15,i2,2i4)', "permneighbours", PERMNEIGHBOURS(I1,1), PERMNEIGHBOURS(I1,2), PERMNEIGHBOURS(I1,3)
         ENDDO
      ENDIF
      DO I1=1, NPERMGROUP
         !IF (I1.EQ.8.OR.I1.EQ.10.OR.I1.EQ.20.OR.I1.EQ.30.OR.I1.EQ.124.OR.I1.EQ.29.OR.I1.EQ.41.OR.I1.EQ.43.OR.I1.EQ.38) CYCLE
         CALL FINDPGCHAIN(I1)
      ENDDO
      IF (DEBUG) THEN
         DO I1=1, NPERMGROUP
            PRINT*, PERMCHAIN(I1,:)
         ENDDO
      ENDIF

      !this would be if I just wanted to calculate the relevant dihedrals
      !CALL LINDH2PERMG(ATOMS2PERM, ATCENTRE2PERM)
 
      DEALLOCATE(MBBONDNUM)
      DEALLOCATE(MBADJACENT)
      DEALLOCATE(CENTER2)


      END SUBROUTINE FINDPERMDIH
! *************************************************
      SUBROUTINE LINDH2PERMG(ATOMS2PERM,ATCENTRE2PERM)
!msb50 see amb_natinters.f90
!finiding number of dihedral (LDH) for particular permgroup
!also need DIHC to find PHI1 (from previous ldh first prevdih in group?)
!               and for alignment
      use intcommons, only: startlindh, NLDH, LINDIH, PERMNEIGHBOURS
      use commons, only: NATOMS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ATOMS2PERM(NATOMS)
      INTEGER, INTENT(IN) ::ATCENTRE2PERM(NATOMS)
      INTEGER:: LDH, I1, DIHC, CURPHI, PG
      INTEGER:: INFO(-1:8)
      INTEGER:: first(NATOMS), second(NATOMS)
      INTEGER:: ATOMS(2)

      first(NATOMS)=0; second(NATOMS) = 0
      DO LDH =1, NLDH
         IF (LINDIH(LDH,-1).NE.1.AND.LINDIH(LDH,0).NE.1) THEN
            ATOMS(1) = LINDIH(LDH,1)
            ATOMS(2) = LINDIH(LDH,2)
            DO I1 = 1,2
            IF (ATCENTRE2PERM(ATOMS(I1)).NE.0) THEN
!               PG = ATOMS2PERM(ATOMS(I1)) !dunno?
                PG = ATCENTRE2PERM(ATOMS(I1))
               IF (first(ATOMS(I1))==0) THEN
                  first(ATOMS(I1))=1
                  PERMNEIGHBOURS(PG, 3+PERMNEIGHBOURS(PG,1))=LDH 
               !second dihedral involved
! every group has got two dihedrals that change though, so this is not quite right yet. Store the other one too!
               ELSEIF (first(ATOMS(I1))==1.AND.second(ATOMS(I1))==0) THEN
                  second(ATOMS(I1))=1
                  PERMNEIGHBOURS(PG, 5+PERMNEIGHBOURS(PG,1))=LDH 
               ELSEIF (first(ATOMS(I1))==1.AND.second(ATOMS(I1))==1) THEN
                  PRINT*, "error in LINDH2PERM"; STOP
               ENDIF
            ENDIF
            ENDDO
            !maybe do the above for all PGs, could facilitate things
        ENDIF
      ENDDO
      END SUBROUTINE LINDH2PERMG
            
!*****************************************************
      SUBROUTINE FINDPERMGROUP(i1,i2, i3,PG1, PG2, PG3)
!msb50 - for atoms in three diff permgroups find PERMGROUP (see keywords.f, PERMDIST)
      use key
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i1, i2, i3
      INTEGER, INTENT(OUT) :: PG1, PG2, PG3
      INTEGER :: start, III, total, sw,II,JJ
      INTEGER :: pgATOMS(8), PGS(3), atoms(3)
      LOGICAL :: found(3)

      ATOMS(:)=0; PGS(:)=0; found(:)=.FALSE.
      ATOMS(1)=i1; ATOMS(2)=i2; ATOMS(3)=i3
      PG1=0;PG2=0;PG3=0
!      PG1(1)=PG1; PGS(2)=PG2; PGS(3)=PG3
      IF (i2.EQ.0) found(2)=.TRUE.
      IF (i3.EQ.0) found(3)=.TRUE.
      start = 1
!       WRITE(*,*) 'in findpermgroup, PGS(:)', PGS(:)
!       STOP
      DO III = 1, NPERMGROUP !loop over all pg's
         pgATOMS(:) = 0 !pgAtoms are all atoms in this permgroup
         IF (NPERMSIZE(III) == 2) THEN
            pgATOMS(1) = PERMGROUP(start);pgATOMS(2)= PERMGROUP(start+1)
            DO sw = 1, NSETS(III)
!              pgATOMS(2+2*(sw-1)+1) = SWAP1(III, sw)
!              pgATOMS(2+2*sw)=SWAP2(III, sw)
               pgATOMS(2+2*(sw-1)+1) = SETS(pgATOMS(1), sw)
               pgATOMS(2+2*sw)=SETS(pgATOMS(2), sw)
            ENDDO
            total = 2 + 2*NSETS(III)
         ELSEIF (NPERMSIZE(III) == 3) THEN
            pgATOMS(1) = PERMGROUP(start);pgATOMS(2)= PERMGROUP(start+1)
            pgATOMS(3) = PERMGROUP(start+2)
            total = 3
         ELSE 
            STOP
         ENDIF

         start = start + NPERMSIZE(III)
         DO II = 1, total
            DO JJ = 1,3 !check whether one of your atoms is in this permgroup
               IF (ATOMS(JJ).EQ.pgATOMS(II)) THEN
                  PGS(JJ)=III   
                  found(JJ)=.TRUE. 
               ENDIF
            ENDDO
         ENDDO
         IF (found(1).AND.found(2).AND.found(3)) THEN
            PG1=PGS(1); PG2=PGS(2);PG3=PGS(3)
            EXIT
         ENDIF
      ENDDO
     
      END SUBROUTINE FINDPERMGROUP
! **************************************************************
      SUBROUTINE FINDPG_NOSWAPS(i1,i2,i3,PG1,PG2,PG3)
!msb50 - for atoms in three diff permgroups find PERMGROUP (see keywords.f, PERMDIST)
!here don't include swaps, so if i2 and i3 are in PG1, but also turn out in PGS independently
      use key
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i1, i2, i3
      INTEGER, INTENT(OUT) :: PG1, PG2, PG3
      INTEGER :: start, III, total, II,JJ
      INTEGER :: pgATOMS(3), PGS(3), atoms(3)
      LOGICAL :: found(3) 

      ATOMS(:)=0; PGS(:)=0; found(:)=.FALSE.
      ATOMS(1)=i1; ATOMS(2)=i2; ATOMS(3)=i3
      PGS(1)=PG1; PGS(2)=PG2; PGS(3)=PG3
      IF (i2.EQ.0) found(2)=.TRUE.
      IF (i3.EQ.0) found(3)=.TRUE.
      start = 1
      DO III = 1, NPERMGROUP !loop over all pg's
         pgATOMS(:) = 0 !pgAtoms are all atoms in this permgroup
         IF (NPERMSIZE(III) == 2) THEN
            pgATOMS(1) = PERMGROUP(start);pgATOMS(2)= PERMGROUP(start+1)
            total = 2 
         ELSEIF (NPERMSIZE(III) == 3) THEN
            pgATOMS(1) = PERMGROUP(start);pgATOMS(2)= PERMGROUP(start+1)
            pgATOMS(3) = PERMGROUP(start+2)
            total = 3
         ELSE
            STOP
         ENDIF

         start = start + NPERMSIZE(III)
         DO II = 1, total
            DO JJ = 1,3 !check whether one of your atoms is in this permgroup
               IF (ATOMS(JJ).EQ.pgATOMS(II)) THEN
                  PGS(JJ)=III
                  found(JJ)=.TRUE.
               ENDIF
            ENDDO
         ENDDO
         IF (found(1).AND.found(2).AND.found(3)) THEN
            PG1=PGS(1); PG2=PGS(2);PG3=PGS(3)
            EXIT
         ENDIF
      ENDDO
      END SUBROUTINE FINDPG_NOSWAPS

! ***************************************************
      SUBROUTINE LOOPCENTERS(CN1, CT2, CN2, CN3,i1)
!msb50 - finds atom of certain centertype bonded to centertypes 2 or 4
        use intcommons
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: CN1 ! first center - number, not atom
        INTEGER, INTENT(IN) :: CT2! centertype of atom adjacent to it -not atom
        INTEGER, INTENT(OUT) :: CN2, CN3 ! number of above center (maximally 2)
        INTEGER, INTENT(IN)  :: i1 !where to loop from, in case you know
        INTEGER :: l1, A2, m1,A3, CT1
     
        CN2 = 0
        CN3 = 0
        CT1 = CENTERS(CN1, 0)
        IF (CT1 == 2) THEN
           A2 = CENTERS(CN1,1) !adjacent to it -can do this as non-terminal at beginning
           A3 = CENTERS(CN1,2)
           DO m1 = i1, NCNT
              IF ((CENTERS(m1,1)==A2 .AND. CENTERS(m1,0)==CT2)) THEN
                 CN2 = m1
              ELSEIF ((CENTERS(m1,1)==A3 .AND. CENTERS(m1,0)==CT2)) THEN
                 CN3 = m1
              ENDIF
              IF (CN2 .NE. 0 .AND. CN3 .NE. 0 ) EXIT 
           ENDDO
        ELSEIF (CT1 == 4) THEN
           A2 = CENTERS(CN1, 1)
           DO m1 = i1, NCNT
             IF ((CENTERS(m1,1)==A2 .AND. CENTERS(m1,0)==CT2)) THEN
                 CN2 = m1
                 EXIT
             ENDIF
           ENDDO
        ENDIF
      END SUBROUTINE LOOPCENTERS

! ****************************************************
      SUBROUTINE FINDPGCHAIN(I)
!msb50 - stores per PERMGROUP which PERMNEIGHBOURS(PG,1) .GT. 1 (i.e. coupled
! permgroups) how large the set is it is part of - i.e. for lysine it would be
! five, which is currently maximum
         use INTCOMMONS, only: PERMNEIGHBOURS, PERMCHAIN
         use key, only: NPERMGROUP
         IMPLICIT NONE
         INTEGER, INTENT(IN) ::I
         INTEGER :: PGCHAIN(5)
         INTEGER :: PG2, PG3, PG4, PG5, II, LEN_CHAIN,START_PG

         IF (.NOT.ALLOCATED(PERMCHAIN)) THEN
             ALLOCATE(PERMCHAIN(NPERMGROUP,6))
             PERMCHAIN(:,:)=0
         ENDIF

         LEN_CHAIN=0; PGCHAIN(:)=0
         IF (PERMNEIGHBOURS(I,1).EQ.0) THEN
            RETURN
         ELSEIF (PERMNEIGHBOURS(I,1).EQ.1) THEN
            START_PG=I
         ELSEIF (PERMNEIGHBOURS(I,1).EQ.2) THEN
            PG2=PERMNEIGHBOURS(I,2); PG3= PERMNEIGHBOURS(I,3)
            IF (PERMNEIGHBOURS(PG2,1)==1) THEN
              START_PG= PG2
            ELSEIF (PERMNEIGHBOURS(PG3,1)==1) THEN
              START_PG= PG3
            ELSE
              IF (PERMNEIGHBOURS(PG2,2).NE.I) THEN
                 IF (PERMNEIGHBOURS(PERMNEIGHBOURS(PG2,2),1).NE.1) THEN
                     PRINT*, "findpgchain: error! chain too long"
                 ELSE
                     START_PG=PERMNEIGHBOURS(PG2,2)
                 ENDIF
              ELSE
                 IF (PERMNEIGHBOURS(PERMNEIGHBOURS(PG2,1),1).NE.1) THEN
                     PRINT*, "findpgchain: error! chain too long"
                 ELSE
                     START_PG=PERMNEIGHBOURS(PG2,1)
                 ENDIF
              ENDIF
            ENDIF
         ENDIF

         PG2=0;PG3=0;PG4=0;PG5=0; PGCHAIN(:)=0
         LEN_CHAIN = LEN_CHAIN+1
         PGCHAIN(1) = START_PG
         PG2=PERMNEIGHBOURS(START_PG,2) !has only one neighbour
         LEN_CHAIN = LEN_CHAIN+1; PGCHAIN(2) = PG2
         IF (PERMNEIGHBOURS(PG2,1).EQ.2) THEN !another neighbour than start_pg
            IF (PERMNEIGHBOURS(PG2,2).NE.START_PG) THEN
                PG3=PERMNEIGHBOURS(PG2,2)
                LEN_CHAIN=LEN_CHAIN+1; PGCHAIN(3)=PG3
            ELSE
                PG3=PERMNEIGHBOURS(PG2,3)
                LEN_CHAIN=LEN_CHAIN+1; PGCHAIN(3)=PG3
            ENDIF
            IF (PG3.NE.0.AND.PERMNEIGHBOURS(PG3,1).EQ.2) THEN !another neighbour-other side 
                IF (PERMNEIGHBOURS(PG3,2).NE.PG2) THEN
                   PG4=PERMNEIGHBOURS(PG3,2)
                   LEN_CHAIN=LEN_CHAIN+1; PGCHAIN(4)=PG4
                ELSE
                   PG4=PERMNEIGHBOURS(PG3,3)
                   LEN_CHAIN=LEN_CHAIN+1; PGCHAIN(4)=PG4
                ENDIF
                IF (PG4.NE.0.AND.(PERMNEIGHBOURS(PG4,1).EQ.2)) THEN
                   IF (PERMNEIGHBOURS(PG4,2).NE.PG3) THEN
                      PG5=PERMNEIGHBOURS(PG4,2)
                      LEN_CHAIN=LEN_CHAIN+1; PGCHAIN(5)=PG5
                   ELSE
                      PG5=PERMNEIGHBOURS(PG4,3)
                      LEN_CHAIN=LEN_CHAIN+1; PGCHAIN(5)=PG5
                   ENDIF
                   IF (PG5.NE.0.AND.(PERMNEIGHBOURS(PG5,1).EQ.2)) THEN
                      PRINT*,"findpgchain>chain longer than atm allowed"
                      !as only lysine - 5 coupled permgroups
                      STOP
                   ENDIF
                ENDIF
             ENDIF
         ENDIF

         IF (LEN_CHAIN.LT.2.OR.LEN_CHAIN.GT.5) THEN
            PRINT*, "findpgchain - sth gone wrong"
            STOP
         ELSE
            PERMCHAIN(START_PG,1) = LEN_CHAIN; PERMCHAIN(PG2,1) = LEN_CHAIN
            PERMCHAIN(START_PG,2:)=PGCHAIN(:); PERMCHAIN(PG2,2:)=PGCHAIN(:)
            IF (LEN_CHAIN.GE.3) THEN
               PERMCHAIN(PG3,1)=LEN_CHAIN
               PERMCHAIN(PG3,2:)=PGCHAIN(:)
               IF (LEN_CHAIN.GE.4) THEN
                  PERMCHAIN(PG4,1)=LEN_CHAIN
                  PERMCHAIN(PG4,2:)=PGCHAIN(:)
                  IF (LEN_CHAIN.GE.5) THEN
                     PERMCHAIN(PG5,1)=LEN_CHAIN
                     PERMCHAIN(PG5,2:)=PGCHAIN(:)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      END SUBROUTINE FINDPGCHAIN



! *****************************************************
!msb50 - give it atom index and it returns name of residue and residue number
      SUBROUTINE GETRESID(ATNUM, RESNAME, RESNUM)
      use modamber9, only: ix, i02, nres, m02,ih
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ATNUM
      CHARACTER(LEN=4), INTENT(OUT):: RESNAME
      INTEGER, INTENT(OUT) :: RESNUM
      INTEGER II

      DO II=1, NRES 
          IF (ix(i02+II-1).LE.ATNUM .AND. ix(i02+II).GT.ATNUM) THEN
             resnum = II
             resname = ih(m02+II-1)
          EXIT
          ENDIF
      ENDDO

      END SUBROUTINE GETRESID

!*********************************************************
      SUBROUTINE PROL_PERMUTE(A1,A2,PTEST, RS,RF, DISTANCE)
      use modamber9, only:  nres, m02, ih, m04
      use commons, only: NATOMS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: A1, A2
      LOGICAL, INTENT(IN)::PTEST
      DOUBLE PRECISION, INTENT(IN) ::RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT)::RF(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT):: DISTANCE
      DOUBLE PRECISION:: RFTMP(3*NATOMS)
      DOUBLE PRECISION:: DIST_orig

      RFTMP(:)=RF(:)
      DIST_orig=DOT_PRODUCT(RF-RS, RF-RS)
      IF (((ih(m04+A1-1).EQ.'HB2 '.OR.ih(m04+A1-1).EQ.'HB3 ').AND.&
     &  (ih(m04+A2-1).EQ.'HB2 '.OR.ih(m04+A2-1).EQ.'HB3 ')).OR.   &
     & ((ih(m04+A1-1).EQ.'HG2 '.OR.ih(m04+A1-1).EQ.'HG3 ').AND.&
     &  (ih(m04+A2-1).EQ.'HG2 '.OR.ih(m04+A2-1).EQ.'HG3 ')).OR.   &
     & ((ih(m04+A1-1).EQ.'HD2 '.OR.ih(m04+A1-1).EQ.'HD3 ').AND.&
     &  (ih(m04+A2-1).EQ.'HD2 '.OR.ih(m04+A2-1).EQ.'HD3 '))) THEN
         
          IF (PTEST) PRINT*, "prol found"
          RF(3*(A1-1)+1:3*A1)=RFTMP(3*(A2-1)+1:3*A2)
          RF(3*(A2-1)+1:3*A2)=RFTMP(3*(A1-1)+1:3*A1)
          DISTANCE=DOT_PRODUCT(RF-RS, RF-RS)
          IF (PTEST) PRINT*, "prol_permute> distance old, new", DIST_orig, DISTANCE
          IF (DISTANCE.LT.DIST_orig) THEN
             IF (PTEST) PRINT*, "prol_permute> keep permutation", & 
     &          DIST_orig,A1, A2, distance
          ELSE
             RF(:) =RFTMP(:)
          ENDIF
      ELSE
         PRINT*, "prol_permute> error! not proline -check intminperm"
         STOP
      ENDIF
      RETURN
 
      END SUBROUTINE PROL_PERMUTE 


