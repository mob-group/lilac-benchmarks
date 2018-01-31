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
C  Calculate the analytic first and second derivatives for the bulk
C  periodic Mie potential.
C
C*************************************************************************
C
      SUBROUTINE MIED(N,X,V,RHO,RE2,ALPHA,XN,AWELL,RE1) 
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), RNF, RMF,RN,RM
      DOUBLE PRECISION VEC(N,N,3),RSQ1(N,N),
     1              R8(N,N),R10(N,N),R14(N,N),
     2              R16(N,N)
C 
      DOUBLE PRECISION RHO, RE2, ALPHA, XN, AWELL, RE1, BOXLX,BOXLY,BOXLZ,CUTOFF
      INTEGER NEXP,MEXP
      NEXP=DINT(RHO)
      MEXP=DINT(RE2)
      BOXLX=ALPHA
      BOXLY=XN
      BOXLZ=AWELL
      CUTOFF=RE1
      
      RN=FLOAT(NEXP)
      RM=FLOAT(MEXP)
      RNF=RM/(RN-RM)
      RMF=RN/(RN-RM)
      PRINT*, 'N, M, N/(N-M), M/(N-M)=', IDNINT(RN), IDNINT(RM), RNF, RMF
C     PRINT*,'incoords'
C     PRINT*, (X(I),X(I+1),X(I+2),I=1,3*N,3)
C
C  Store distance matrices.
C
      PRINT*,'WARNING - VEC was not set! Not tested!'

      DO 20 J1=1,N
         R8(J1,J1)=0.0D0
         R10(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         R16(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 10 J2=1,J1-1
            IF (DSQRT(RSQ1(J2,J1)).LT.CUTOFF) THEN
               R8(J2,J1)=RSQ1(J2,J1)**(-(RM/2.0D0+1.0D0))
               R10(J2,J1)=RSQ1(J2,J1)**(-(RM/2.0D0+2.0D0))
               R14(J2,J1)=RSQ1(J2,J1)**(-(RN/2.0D0+1.0D0))
               R16(J2,J1)=RSQ1(J2,J1)**(-(RN/2.0D0+2.0D0))
               VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
               VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
               VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            ELSE
               R8(J2,J1)=0.0D0
               R10(J2,J1)=0.0D0
               R14(J2,J1)=0.0D0
               R16(J2,J1)=0.0D0
               VEC(J2,J1,1)=0.0D0
               VEC(J2,J1,2)=0.0D0
               VEC(J2,J1,3)=0.0D0
            ENDIF
            R8(J1,J2)=R8(J2,J1)
            R10(J1,J2)=R10(J2,J1)
            R14(J1,J2)=R14(J2,J1)
            R16(J1,J2)=R16(J2,J1)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
10       CONTINUE
20    CONTINUE
C
C  First calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
               V(J3)=V(J3)-(RN*RNF*R14(J1,J4)-RM*RMF*R8(J1,J4))
     1              *VEC(J1,J4,J2)
30          CONTINUE
40       CONTINUE
50    CONTINUE
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
               HESS(J3,J3)=HESS(J3,J3)+(RN*(RN+2.0D0)*RNF*R16(J4,J1)
     1                  -RM*(RM+2.0D0)*RMF*R10(J4,J1))
     2                  *VEC(J1,J4,J2)**2
     3                  -RN*RNF*R14(J4,J1)+RM*RMF*R8(J4,J1)
60          CONTINUE
70       CONTINUE
80    CONTINUE
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=1,J2-1
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 90 J5=1,N
                  HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +(RN*(RN+2.0D0)*RNF*R16(J1,J5)
     2           -RM*(RM+2.0D0)*RMF*R10(J1,J5))
     3           *VEC(J1,J5,J4)*VEC(J1,J5,J2)
90             CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=1,J1-1
               HESS(J3,3*(J4-1)+J2)=-(RN*(RN+2.0D0)*RNF*R16(J1,J4)
     1                            -RM*(RM+2.0D0)*RMF*R10(J1,J4)) 
     2                            *VEC(J1,J4,J2)**2
     3                            +RN*RNF*R14(J1,J4)
     4                            -RM*RMF*R8(J1,J4)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
130         CONTINUE
140      CONTINUE
150   CONTINUE
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=1,J1-1
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=-(RN*(RN+2.0D0)*RNF*R16(J1,J4)
     1                     -RM*(RM+2.0D0)*RMF*R10(J1,J4))
     2                     *VEC(J1,J4,J2)*VEC(J1,J4,J5)
                  HESS(J6,J3)=HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=-(RN*(RN+2)*RNF*R16(J1,J4)
     1                     -RM*(RM+2)*RMF*R10(J1,J4))
     2                     *VEC(J1,J4,J2)*VEC(J1,J4,J5)
                  HESS(J6,J3)=HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
