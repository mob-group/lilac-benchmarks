C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C*******************************************************************
C
C subroutine MIEL for bulk periodic Mie potential calculates the
C lattice constant with minimum energy:
C
C*******************************************************************
C
      SUBROUTINE MIEL(N,P,EPS,SIG,NEXP,MEXP,BOXLX,BOXLY,BOXLZ,CUTOFF,Q)
      IMPLICIT NONE
      INTEGER I, J, J1, J2, N,NEXP,MEXP
      DOUBLE PRECISION P(3*N), Q(3*N), VEC(N,N,3),RSQ1(N,N)
      DOUBLE PRECISION RNF, RMF, RNP, RMP, AX, BX, CX, X0, X3, X1, X2, F1, F2, XMIN,
     1                 EPS, SIG, BOXLX,BOXLY,BOXLZ,CUTOFF,SN,SM,RO,C,TOL
      PARAMETER (RO=0.61803399D0,C=0.38196602D0,TOL=1.0D-10) 
C     COMMON /WORK/ VEC(N,N,3),RSQ1(N,N),DUMMY(5*N*N)
C
      SN=0.0D00
      SM=0.0D00
      RNF=FLOAT(MEXP)/FLOAT(NEXP-MEXP)
      RMF=FLOAT(NEXP)/FLOAT(NEXP-MEXP)
      RNP=FLOAT(NEXP)/2.0D0
      RMP=FLOAT(MEXP)/2.0D0
C
      DO 211 I=1,N
         DO 10 J=I+1,N
            IF(DSQRT(RSQ1(I,J)).LT.CUTOFF) THEN
               SN=SN + RNF/RSQ1(I,J)**RNP
               SM=SM + RMF/RSQ1(I,J)**RMP
            ENDIF
10       CONTINUE
211   CONTINUE
      AX=0.9*SIG
      BX=SIG
      CX=1.1*SIG
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      F1=EPS*X1**NEXP*SN - EPS*X1**MEXP*SM
      F2=EPS*X2**NEXP*SN - EPS*X2**MEXP*SM
1     IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))THEN
        IF(F2.LT.F1)THEN
          X0=X1
          X1=X2
          X2=RO*X1+C*X3
          F1=F2
          F2=EPS*X2**NEXP*SN - EPS*X2**MEXP*SM
        ELSE
          X3=X2
          X2=X1
          X1=RO*X2+C*X0
          F2=F1
          F1=EPS*X1**NEXP*SN - EPS*X1**MEXP*SM
        ENDIF
      GOTO 1
      ENDIF
      IF(F1.LT.F2)THEN
        XMIN=SIG/X1
      ELSE
        XMIN=SIG/X2
      ENDIF
      BOXLX=BOXLX*XMIN
      BOXLY=BOXLY*XMIN
      BOXLZ=BOXLZ*XMIN
      CUTOFF=CUTOFF*XMIN
      DO 29 J1=1,N
         DO 28 J2=1,N
            RSQ1(J2,J1)=RSQ1(J2,J1)*XMIN*XMIN
            VEC(J2,J1,1)=VEC(J2,J1,1)*XMIN
            VEC(J2,J1,2)=VEC(J2,J1,2)*XMIN
            VEC(J2,J1,3)=VEC(J2,J1,3)*XMIN
28       CONTINUE
29    CONTINUE
      DO 30 J1=1,3*N
         P(J1)=P(J1)*XMIN
         Q(J1)=Q(J1)*XMIN
30    CONTINUE
      RETURN
      END
