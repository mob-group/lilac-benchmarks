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
C  Subroutine MDIFF calculates the cartesian gradient and second
C  derivative matrix analytically for the Morse potential.
C
C*************************************************************************
C
      SUBROUTINE MDIFF(N, X, V, RHO) 
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6 
      DOUBLE PRECISION X(3*N), RHO,
     1                 V(3*N), DIST, R(N,N), RR(N,N)
C 
C  Store distance matrices.  
C 
      DO 20 J1=1,N 
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1          +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2          +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R(J2,J1)=RHO*DSQRT(DIST)
            RR(J2,J1)=1.0D0/R(J2,J1)
            R(J1,J2)=R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
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
               V(J3)=V(J3)-2.0*RHO*RHO*(X(J3)-X(3*(J4-1)+J2))*
     1               DEXP(RHO-R(J4,J1))* 
     2              (DEXP(RHO-R(J4,J1))-1.0D0)*RR(J4,J1)
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
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
               HESS(J3,J3)=HESS(J3,J3)+2.0*RHO*RHO*DEXP(RHO-R(J4,J1))
     1         *(  (DEXP(RHO-R(J4,J1))-1.0D0) 
     2         * ( (RHO * (X(J3)-X(3*(J4-1)+J2) )*RR(J4,J1) )**2-1.0D0)  
     4         + ( RHO*(X(J3)-X(3*(J4-1)+J2)) )**2*RR(J4,J1)*
     5             (2.0*DEXP(RHO-R(J4,J1))-1.0D0) ) * RR(J4,J1)
60          CONTINUE
C           PRINT*,'HESS(',J3,',',J3,')=',HESS(J3,J3)
70       CONTINUE
80    CONTINUE
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=J2+1,3
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 90 J5=1,N
                 HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +2.0*RHO**4*DEXP(RHO-R(J5,J1))
     2           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     3           *( (DEXP(RHO-R(J5,J1))-1.0D0)*RR(J5,J1)+
     4              2.0*DEXP(RHO-R(J5,J1))-1.0D0) * RR(J5,J1)**2
90             CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
C              PRINT*,'HESS(',J3,',',3*(J1-1)+J4,')=',HESS(J3,3*(J1-1)+J4)
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=J1+1,N
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*DEXP(RHO-R(J4,J1))
     1      *(  (-RHO**2*(X(J3)-X(3*(J4-1)+J2))**2*RR(J4,J1)**2+1.0) 
     2           *   (DEXP(RHO-R(J4,J1))-1.0D0)
     3           -   (X(J3)-X(3*(J4-1)+J2))**2 * RR(J4,J1) * RHO**2
     4           *   (2.0*DEXP(RHO-R(J4,J1))-1.0D0)   )  * RR(J4,J1)
C              PRINT*,'HESS(',J3,',',3*(J4-1)+J2,')=',HESS(J3,3*(J4-1)+J2)
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
            DO 160 J4=J1+1,N
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1           -2.0*RHO**4*DEXP(RHO-R(J4,J1))
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( (DEXP(RHO-R(J4,J1))-1.0D0)*RR(J4,J1)+
     4              2.0*DEXP(RHO-R(J4,J1))-1.0D0) * RR(J4,J1)**2

                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1           -2.0*RHO**4*DEXP(RHO-R(J4,J1))
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( (DEXP(RHO-R(J4,J1))-1.0D0)*RR(J4,J1)+
     4              2.0*DEXP(RHO-R(J4,J1))-1.0D0) * RR(J4,J1)**2
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
