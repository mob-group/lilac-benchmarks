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
C  Subroutine DTRAP calculates the cartesian gradient and second
C  derivative matrix analytically.
C
C*************************************************************************
C
      SUBROUTINE DTRAP(N, X, V, C1, C2, C3, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION X(3*N), V(3*N),DIST,K
      DOUBLE PRECISION R3(N,N), R5(N,N), CENTRE(3), C1, C2, C3 
      LOGICAL STEST, GTEST
C 
      PARAMETER (K=1.0D0)

      IF (.NOT.GTEST) RETURN
C
C  Assume the fixed origin for the harmonic potential is at the origin, not the
C  centre of mass, otherwise the derivatives must include terms for the variation
C  of the centre of mass with respect to coordinates.
C
C     CENTRE(1)=C1
C     CENTRE(2)=C2
C     CENTRE(3)=C3
C
C  Store distance matrices.
C
      DO 10 J1=1,N
         R3(J1,J1)=0.0D0
         R5(J1,J1)=0.0D0
         DO 20 J2=J1+1,N
            DIST=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1              +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2              +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R3(J1,J2)=1.0D0/DIST**(1.5D0)
            R5(J1,J2)=1.0D0/DIST**(2.5D0)
            R3(J2,J1)=R3(J1,J2)
            R5(J2,J1)=R5(J1,J2)
20       CONTINUE
10    CONTINUE
C
C  First calculate the gradient analytically.
C
      DO 30 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 50 J4=1,N
               V(J3)=V(J3)+(X(3*(J4-1)+J2)-X(J3))*R3(J4,J1)
50          CONTINUE
C           V(J3)=V(J3)+K*(X(J3)-CENTRE(J2))
            V(J3)=V(J3)+K*X(J3)
40       CONTINUE
30    CONTINUE
      IF (.NOT.STEST) RETURN
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 60 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 80 J4=1,N
               HESS(J3,J3)=HESS(J3,J3)+3.0D0*(X(J3)-X(3*(J4-1)+J2))**2*R5(J4,J1)-R3(J4,J1)
80          CONTINUE
            HESS(J3,J3)=HESS(J3,J3)+K
70       CONTINUE
60    CONTINUE
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO 90 J1=1,N
         DO 100 J2=1,3
            J3=3*(J1-1)+J2
            DO 110 J4=J2+1,3
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 120 J5=1,N
                  HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))*3.0D0*R5(J5,J1)
120            CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
110         CONTINUE
100      CONTINUE
90    CONTINUE
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO 130 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 150 J4=J1+1,N
               HESS(J3,3*(J4-1)+J2)=-3.0D0*R5(J4,J1)*(X(J3)-X(3*(J4-1)+J2))**2+R3(J4,J1)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
150         CONTINUE
140      CONTINUE
130   CONTINUE
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=J1+1,N
               DO 155 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=-R5(J4,J1)*3.0D0
     1                       *(X(J3)-X(3*(J4-1)+J2))
     2                       *(X(3*(J1-1)+J5)-X(J6))
                  HESS(J6,J3)=HESS(J3,J6)
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J3,J6)
                  HESS(3*(J1-1)+J5,3*(J4-1)+J2)=HESS(J3,J6)
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
