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
C Energy subroutine for bulk periodic Mie potential
C
C*******************************************************************
C
      SUBROUTINE EMIE(N,P,POTEL1,RHO,RE2,BOXLX,BOXLY,BOXLZ,CUTOFF,PRESSURE)   
      IMPLICIT NONE
      INTEGER J, I, MX, MY, MZ, J1, J2, N, NEXP, MEXP
      DOUBLE PRECISION P2, P(3*N),RHO,RE2
      LOGICAL PRESSURE
      DOUBLE PRECISION VEC(N,N,3),RSQ1(N,N),POT,Q,RMP,RNP,RMF,RNF,CUTOFF,
     1                  POTEL1, BOXLX, BOXLZ, BOXLY, SIG, EPS
      PARAMETER (SIG=1.0D0, EPS=1.00D0) 

      NEXP=DINT(RHO)
      MEXP=DINT(RE2)

      POTEL1=0.0D00
      P2=0.0D00
      RNF=FLOAT(MEXP)/FLOAT(NEXP-MEXP)*EPS*SIG**NEXP
      RMF=FLOAT(NEXP)/FLOAT(NEXP-MEXP)*EPS*SIG**MEXP
      RNP=FLOAT(NEXP)/2.0D0
      RMP=FLOAT(MEXP)/2.0D0
      DO 21 J1=1,N
         P(3*(J1-1)+1)=P(3*(J1-1)+1) - BOXLX *
     1                  DNINT(P(3*(J1-1)+1)/BOXLX)
         P(3*(J1-1)+2)=P(3*(J1-1)+2) -  BOXLY *
     1                  DNINT(P(3*(J1-1)+2)/BOXLY)
         P(3*(J1-1)+3)=P(3*(J1-1)+3) -  BOXLZ *
     1                  DNINT(P(3*(J1-1)+3)/BOXLZ)
21    CONTINUE
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=P(3*(J2-1)+1)-P(3*(J1-1)+1)
            VEC(J2,J1,2)=P(3*(J2-1)+2)-P(3*(J1-1)+2)
            VEC(J2,J1,3)=P(3*(J2-1)+3)-P(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            MZ=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(MX)
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(MY)
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(MZ)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
            RSQ1(J1,J2)=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1                  VEC(J1,J2,3)**2
            RSQ1(J2,J1)=RSQ1(J1,J2)
15       CONTINUE
25    CONTINUE
C
C call miel for lattice constant optimisation if required:
C
      IF (PRESSURE) THEN
         CALL MIEL(N,P,EPS,SIG,NEXP,MEXP,BOXLX,BOXLY,BOXLZ,CUTOFF,Q) 
         PRINT*,'Energy minimised with respect to lattice constants' 
      ENDIF
      DO 211 I=1,N
         DO 10 J=I+1,N
            IF(DSQRT(RSQ1(I,J)).LT.CUTOFF) THEN
               POT=RNF/RSQ1(I,J)**RNP-RMF/RSQ1(I,J)**RMP
               P2=P2+POT
            ENDIF
10       CONTINUE
211   CONTINUE
      POTEL1=P2
      RETURN
      END
