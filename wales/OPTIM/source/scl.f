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
C*************************************************************************
C
C  Subroutine SCL calculates the optimum lattice constant for the SC 
C  potential
C
C*************************************************************************
C
      SUBROUTINE SCL(N,P,EPS,C,SIG,BOXLX,BOXLY,BOXLZ,CUTOFF)
      IMPLICIT NONE
      INTEGER N, J1, J2, I, NN, MM, J
      DOUBLE PRECISION R, CC, TOL
      PARAMETER (R=.61803399,CC=.38196602,TOL=.1D-05)
      COMMON /POWERS/ NN, MM
      DOUBLE PRECISION RM(N,N), RN(N,N),VEC(N,N,3)
      DOUBLE PRECISION RHO(N), P(3*N), XMIN, SN, SM, AX, BX, CX, X0, X3, X1, X2, F1, F2,
     1                 EPS, C, SIG, BOXLX, BOXLY, BOXLZ, CUTOFF
C
      DO 11 I=1,N
         RHO(I)=0.0D00
         DO 122 J=1,N
            RHO(I)=RHO(I) + RM(I,J)
122      CONTINUE
11    CONTINUE
C
C First calculate the sums:
C
      SN=0.0D0
      SM=0.0D0
      DO 13 I=1,N
         DO 14 J=1,N
            SN=SN + 0.50D00*EPS*RN(I,J)
14       CONTINUE
         SM=SM + EPS*DSQRT(RHO(I))*C
13    CONTINUE
C
C now optimise using golden sectioning:
C
      AX=0.9*SIG
      BX=SIG
      CX=1.1*SIG
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+CC*(CX-BX)
      ELSE
        X2=BX
        X1=BX-CC*(BX-AX)
      ENDIF
      F1=X1**NN*SN - X1**(MM/2.0D0)*SM
      F2=X2**NN*SN - X2**(MM/2.0D0)*SM
1     IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))THEN
        IF(F2.LT.F1)THEN
          X0=X1
          X1=X2
          X2=R*X1+CC*X3
          F1=F2
          F2=X2**NN*SN - X2**(MM/2.0D0)*SM
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+CC*X0
          F2=F1
          F1=X1**NN*SN - X1**(MM/2.0D0)*SM
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
      DO 28 J1=1,N
         DO 29 J2=1,N
            RM(J2,J1)=RM(J2,J1)*XMIN**(-FLOAT(MM))
            RN(J2,J1)=RN(J2,J1)*XMIN**(-FLOAT(NN))
            VEC(J2,J1,1)=VEC(J2,J1,1)*XMIN
            VEC(J2,J1,2)=VEC(J2,J1,2)*XMIN
            VEC(J2,J1,3)=VEC(J2,J1,3)*XMIN
29       CONTINUE
28    CONTINUE
      DO 30 J1=1,3*N
         P(J1)=P(J1)*XMIN
30    CONTINUE
      RETURN
      END

