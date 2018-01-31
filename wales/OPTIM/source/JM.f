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
C  Two body term of Murrell Potential
C  Subroutine V2DIFF calculates the cartesian gradient and second
C  derivative matrix analytically. 
C
C*************************************************************************
C
      SUBROUTINE JM2(N, X, V)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION RE, D, AN
      PARAMETER(D=2.918D0)
      PARAMETER(AN=6.50D0)
      PARAMETER(RE=2.389D0)
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX, 
     2                 RAB(N,N), RRAB(N,N),
     3                 RHOAB(N,N)
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         RRAB(J1,J1)=0.0D0
         RAB(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            RAB(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RAB(J2,J1)=DSQRT(RAB(J2,J1))
            RAB(J1,J2)=RAB(J2,J1)
            RRAB(J2,J1)=1.0D0/RAB(J2,J1)
            RRAB(J1,J2)=RRAB(J2,J1)
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
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
               ABX=X(3*(J4-1)+J2)-X(J3)
               V(J3)=V(J3)
     1          +AN*ABX*D*(RRAB(J4,J1)-(1+AN*RHOAB(J4,J1))*RRAB(J4,J1)) 
     2           /(DEXP(AN*RHOAB(J4,J1))*RE) 
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
               ABX=X(3*(J4-1)+J2)-X(J3)
               HESS(J3,J3)=HESS(J3,J3)+
     1   D*(AN**2*ABX**2*(1-AN*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2   AN*(AN*RHOAB(J4,J1)*RRAB(J4,J1)+ABX**2*(RRAB(J4,J1)**3-
     3   (1+AN*RHOAB(J4,J1))*RRAB(J4,J1)**3))/RE)/DEXP(AN*RHOAB(J4,J1)) 
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
            DO 100 J4=J2+1,3
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO 90 J5=1,N
                  ABX=-X(J3)+X(3*(J5-1)+J2)
                  ABY=-X(3*(J1-1)+J4)+X(3*(J5-1)+J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABX*ABY*D*(AN**2*(1-AN*RHOAB(J5,J1))*RRAB(J5,J1)**2/RE**2+  
     2    AN*(RRAB(J5,J1)**3-(1+AN*RHOAB(J5,J1))*RRAB(J5,J1)**3)/RE)/ 
     3    DEXP((AN*RHOAB(J5,J1)))
90             CONTINUE
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
              ABX=-X(J3)+X(3*(J4-1)+J2)
                HESS(3*(J4-1)+J2,J3)=
     1     D*(AN**2*ABX**2*(-1+AN*RHOAB(J4,J1))/(RAB(J4,J1)**2*RE**2)+   
     2     AN*(RRAB(J4,J1)-(1+AN*RHOAB(J4,J1))*RRAB(J4,J1)+
     3     AN*ABX**2*RHOAB(J4,J1)*RRAB(J4,J1)**3)/RE)/
     4     DEXP((AN*RHOAB(J4,J1)))
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
                  ABX=-X(J3)+X(3*(J4-1)+J2)
                  ABY=-X(3*(J1-1)+J5)+X(J6)
                  HESS(J6,J3)=
     1          AN**2*ABX*ABY*D*((-1+AN*RHOAB(J4,J1))/
     2          (RAB(J4,J1)**2*RE**2)+
     3          RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/DEXP(AN*RHOAB(J4,J1)) 
                HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
C
C  Symmetrise Hessian
C
      DO 200 J1=1,3*N
         DO 190 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Two body term of Murrell Potential - analytic derivatives with
C  periodic boundary conditions
C
C*************************************************************************
C
      SUBROUTINE JM2C (N, X, V, CUTOFF, RAB, RRAB, VEC, RHOAB, M, NDUM, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      INTEGER M(N,N,3)
      LOGICAL STEST
      DOUBLE PRECISION RAB(N,N), RRAB(N,N), 
     1              VEC(N,N,3), RHOAB(N,N), NDUM(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF, SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX
C 
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         DO 10 J2=J1+1,N
            IF (RAB(J2,J1).GE.CUTOFF) THEN
               RRAB(J2,J1)=0.0D0
               RRAB(J1,J2)=0.0D0
            ENDIF
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
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
                  ABX=VEC(J4,J1,J2)
                  V(J3)=V(J3)
     1              +AN2*ABX*D*(RRAB(J4,J1)-
     2              (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)) 
     3              /(DEXP(AN2*RHOAB(J4,J1))*RE) 
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
40       CONTINUE
50    CONTINUE
      IF (.NOT.STEST) RETURN
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
               ABX=VEC(J4,J1,J2)
               HESS(J3,J3)=HESS(J3,J3)+
     1 D*(AN2**2*ABX**2*(1-AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2 AN2*(AN2*RHOAB(J4,J1)*RRAB(J4,J1)+ABX**2*(RRAB(J4,J1)**3-
     3 (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**3))/RE)
     4 /DEXP(AN2*RHOAB(J4,J1)) 
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
            DO 100 J4=J2+1,3
                HESS(3*(J1-1)+J4,J3)=0.0D0
                DO 90 J5=1,N
                  ABX=VEC(J5,J1,J2)
                  ABY=VEC(J5,J1,J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABX*ABY*D*(AN2**2*(1-AN2*RHOAB(J5,J1))*RRAB(J5,J1)**2/RE**2+  
     2    AN2*(RRAB(J5,J1)**3-(1+AN2*RHOAB(J5,J1))*RRAB(J5,J1)**3)/RE)/ 
     3    DEXP((AN2*RHOAB(J5,J1)))
90             CONTINUE
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
              ABX=VEC(J4,J1,J2)
                HESS(3*(J4-1)+J2,J3)=
     1     D*(AN2**2*ABX**2*(-1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2     AN2*(RRAB(J4,J1)-(1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)+
     3     AN2*ABX**2*RHOAB(J4,J1)*RRAB(J4,J1)**3)/RE)/
     4     DEXP((AN2*RHOAB(J4,J1)))
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
                  ABX=VEC(J4,J1,J2)
                  ABY=VEC(J4,J1,J5)
                  HESS(J6,J3)=AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     1                     RRAB(J4,J1)**2/RE**2+RHOAB(J4,J1)
     2                    *RRAB(J4,J1)**3/RE)/DEXP(AN2*RHOAB(J4,J1)) 

                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=
     1          + AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     2          RRAB(J4,J1)**2/RE**2+RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/
     3          DEXP(AN2*RHOAB(J4,J1)) 
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
C
C  Symmetrise Hessian
C
      DO 200 J1=1,3*N
         DO 190 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Two body term of Murrell Potential - analytic derivatives
C
C*************************************************************************
C
      SUBROUTINE JM2CC(N, X, V,RAB,RRAB,VEC,RHOAB,M,NDUM)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      INTEGER M(N,N,3)
      DOUBLE PRECISION RAB(N,N), RRAB(N,N),
     1              VEC(N,N,3), RHOAB(N,N), NDUM(N,N)
      DOUBLE PRECISION CUTOFF
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX, ABIG, B, C
      PARAMETER (ABIG=13.29D0, B=0.574D0, C=1.0D0, CUTOFF=2.8D0)
C   
C      INITIALISE THE FIRST AND SECOND DERIVATIVES
C
      DO 21 J1=1,3*N
         V(J1)=0.0D0
         DO 22 J2=1,3*N
            HESS(J2,J1)=0.0D0
22       CONTINUE
21    CONTINUE
C
C  First calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            DO 30 J4=1,N
            IF ((RAB(J4,J1).LT.CUTOFF).AND.(J4.NE.J1)) THEN
               ABX=VEC(J4,J1,J2)
               V(J3)=V(J3)+
     1         ABIG*ABX*DEXP(C/RHOAB(J4,J1))*(4.0D0*B/RAB(J4,J1)**6 
     2      + C*(-1.0D0 + B/RAB(J4,J1)**4)*RRAB(J4,J1)/RHOAB(J4,J1)**2)
            ENDIF
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
            DO 60 J4=1,N
            IF ((RAB(J4,J1).GE.CUTOFF).OR.(J4.EQ.J1)) GOTO 60
               ABX=VEC(J4,J1,J2)
               HESS(J3,J3)=HESS(J3,J3)+
     1   ABIG*DEXP(C/RHOAB(J4,J1))*(B*(-4/RAB(J4,J1)**6 + 
     2   ABX**2*(24/RAB(J4,J1)**8+8*C/(RAB(J4,J1)**7*RHOAB(J4,J1)**2))) 
     3   + (-1 + B/RAB(J4,J1)**4)*(-(C*RRAB(J4,J1)/RHOAB(J4,J1)**2) + 
     4   ABX**2*((C**2/RHOAB(J4,J1)**4 +
     5   2*C/RHOAB(J4,J1)**3)/RAB(J4,J1)**2 + 
     6   C*RRAB(J4,J1)**3/RHOAB(J4,J1)**2)))
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
            DO 100 J4=J2+1,3
                HESS(3*(J1-1)+J4,J3)=0.0D0
                DO 90 J5=1,N
                  IF ((RAB(J5,J1).GE.CUTOFF).OR.(J5.EQ.J1)) GOTO 90 
                  ABX=VEC(J5,J1,J2)
                  ABY=VEC(J5,J1,J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABIG*ABX*ABY*DEXP(C/RHOAB(J5,J1))*(B*(24/RAB(J5,J1)**8 +
     2    8*C/(RAB(J5,J1)**7*RHOAB(J5,J1)**2)) +
     3    (-1.0D0+B/RAB(J5,J1)**4)*((C**2/RHOAB(J5,J1)**4 + 
     4    2*C/RHOAB(J5,J1)**3)/RAB(J5,J1)**2 + 
     5    C*RRAB(J5,J1)**3/RHOAB(J5,J1)**2))
90             CONTINUE
            HESS(J3,3*(J1-1)+J4)=HESS(3*(J1-1)+J4,J3)
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
                IF (RAB(J4,J1).GE.CUTOFF) GOTO 130 
                ABX=VEC(J4,J1,J2)
                HESS(3*(J4-1)+J2,J3)=
     1    ABIG*DEXP(C/RHOAB(J4,J1))*(B*(4/RAB(J4,J1)**6 +
     2    ABX**2*(-24/RAB(J4,J1)**8 
     3    - 8*C/(RAB(J4,J1)**7*RHOAB(J4,J1)**2))) + 
     4    (-1.0D0 + B/RAB(J4,J1)**4)*(C*RRAB(J4,J1)/RHOAB(J4,J1)**2 + 
     5    ABX**2*((-C**2/RHOAB(J4,J1)**4 -
     6    2*C/RHOAB(J4,J1)**3)/RAB(J4,J1)**2 - 
     7    C*RRAB(J4,J1)**3/RHOAB(J4,J1)**2)))
            HESS(J3,3*(J4-1)+J2)=HESS(3*(J4-1)+J2,J3)
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
               IF (RAB(J4,J1).GE.CUTOFF) GOTO 160 
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  ABX=VEC(J4,J1,J2)
                  ABY=VEC(J4,J1,J5)
                  HESS(J6,J3)=
     1    ABIG*ABX*ABY*DEXP(C/RHOAB(J4,J1))*(B*(-24/RAB(J4,J1)**8 -  
     2    8*C/(RAB(J4,J1)**7*RHOAB(J4,J1)**2)) + 
     3    (-1.0D0+B/RAB(J4,J1)**4)*((-C**2/RHOAB(J4,J1)**4 -
     4    2.0D0*C/RHOAB(J4,J1)**3)/RAB(J4,J1)**2 - 
     5    C*RRAB(J4,J1)**3/RHOAB(J4,J1)**2))
               HESS(J3,J6)=HESS(J6,J3)
               HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               HESS(3*(J1-1)+J5,3*(J4-1)+J2)=HESS(J6,J3)
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
C
C  Symmetrise Hessian
C
      DO 200 J1=1,3*N
         DO 190 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Two body term of Murrell Potential - analytic derivatives with
C  periodic boundary conditions
C
C*************************************************************************
C
      SUBROUTINE JM2P (N, X, V, BOXLX, BOXLY, BOXLZ, CUTOFF, RAB, VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION RAB(N,N), RRAB(N,N),
     1              VEC(N,N,3), RHOAB(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF
      DOUBLE PRECISION X(3*N), V(3*N), ABY, ABX,
     1                 SR2, SR3, SR6, C0, C1, C2, C3, C4, BOXLX, BOXLY, BOXLZ,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
C 
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         DO 10 J2=J1+1,N
            IF (RAB(J2,J1).GE.CUTOFF) THEN
               RRAB(J2,J1)=0.0D0
               RRAB(J1,J2)=0.0D0
            ENDIF
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
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
                  ABX=VEC(J4,J1,J2)
                  V(J3)=V(J3)
     1              +AN2*ABX*D*(RRAB(J4,J1)-
     2              (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)) 
     3              /(DEXP(AN2*RHOAB(J4,J1))*RE) 
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
C           HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
               ABX=VEC(J4,J1,J2)
               HESS(J3,J3)=HESS(J3,J3)+
     1 D*(AN2**2*ABX**2*(1-AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2 AN2*(AN2*RHOAB(J4,J1)*RRAB(J4,J1)+ABX**2*(RRAB(J4,J1)**3-
     3 (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**3))/RE)
     4 /DEXP(AN2*RHOAB(J4,J1)) 
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
            DO 100 J4=J2+1,3
C               HESS(3*(J1-1)+J4,J3)=0.0D0
                DO 90 J5=1,N
                  ABX=VEC(J5,J1,J2)
                  ABY=VEC(J5,J1,J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABX*ABY*D*(AN2**2*(1-AN2*RHOAB(J5,J1))*RRAB(J5,J1)**2/RE**2+  
     2    AN2*(RRAB(J5,J1)**3-(1+AN2*RHOAB(J5,J1))*RRAB(J5,J1)**3)/RE)/ 
     3    DEXP((AN2*RHOAB(J5,J1)))
90             CONTINUE
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
              ABX=VEC(J4,J1,J2)
                HESS(3*(J4-1)+J2,J3)=HESS(3*(J4-1)+J2,J3) +
     1     D*(AN2**2*ABX**2*(-1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+   
     2     AN2*(RRAB(J4,J1)-(1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)+
     3     AN2*ABX**2*RHOAB(J4,J1)*RRAB(J4,J1)**3)/RE)/
     4     DEXP((AN2*RHOAB(J4,J1)))
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
                  ABX=VEC(J4,J1,J2)
                  ABY=VEC(J4,J1,J5)
                  HESS(J6,J3)=
     1          AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     2          RRAB(J4,J1)**2/RE**2+
     3          RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/DEXP(AN2*RHOAB(J4,J1)) 
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=
     1          + AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     2          RRAB(J4,J1)**2/RE**2+
     3          RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/DEXP(AN2*RHOAB(J4,J1)) 
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
C
C  Symmetrise Hessian
C
      DO 200 J1=1,3*N
         DO 190 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Derivatives of the JM three-body terms - to save space these are
C  added directly to the supplied gradient and second derivative
C  matrix.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3(N, X, V)
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5
      DOUBLE PRECISION SR2, SR3, SR6, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, RE, D, AN
      PARAMETER(SR2=1.414213562D0)
      PARAMETER(SR3=1.732050808D0)
      PARAMETER(SR6=2.449489743D0)
      PARAMETER(C0=3.598D0)
      PARAMETER(C1=-11.609D0)
      PARAMETER(C2=13.486D0)
      PARAMETER(C3=-18.174D0)
      PARAMETER(C4=-5.570D0)
      PARAMETER(C5=79.210D0)
      PARAMETER(C6=-6.458D0)
      PARAMETER(C7=23.383D0)
      PARAMETER(C8=-111.890D0)
      PARAMETER(C9=9.705D0)
      PARAMETER(C10=38.297D0)
      PARAMETER(RE=2.389D0)
      PARAMETER(D=2.918D0)
      PARAMETER(AN=6.50D0)
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX,
     1                 R2(N,N), V(3*N), 
     2                 VEC(N,N,3), TEMP, RR2(N,N), 
     3                 RAB, RRAB, RAC, RRAC, RBC, RRBC, 
     4                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     5                 QQ3, QB, QC
      DO 20 J1=1,N
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 10 J2=J1+1,N
            R2(J2,J1)=
     1                ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2                ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3                ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            R2(J2,J1)=DSQRT(R2(J2,J1))
            RR2(J2,J1)=1.0/R2(J2,J1)
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            RR2(J1,J2)=RR2(J2,J1) 
            R2(J1,J2)=R2(J2,J1) 
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
10       CONTINUE 
20    CONTINUE 
C
C  First the gradient.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/re)/
     9           DEXP(AN*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=TEMP
     1                     + V(3*(J1-1)+J2)
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
C
C  Diagonal bits of the Hessian.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 130 J4=J3+1,N
                  IF (J4.NE.J1) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=(RBC-RE)/RE
                  RHOAC=(RAC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ2=(RHOBC-RHOAC)/SR2
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QB=(QQ2**2)+(QQ3**2)
                  QC=(QQ3**3)-3*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1      D*((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     2      AN**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     3      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))*
     4      (-(abx*RRAB/re)-acx*RRAC/re)**2/SR3**2+
     5      ((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*(c2+c8*QB)+
     6      c10*QC-AN*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     7      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))*
     8      (RRAB-abx**2*RRAB**3+RRAC-acx**2*RRAC**3)/re+
     9      (-(abx*RRAB/re)-acx*RRAC/re)*
     a      (-2*AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     b      2*QA*(c2+c8*QB)+c10*QC)*
     c      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     d      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     e      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     f      +(c6+c10*QA)*
     g      (-6*acx*QQ2*QQ3*RRAC/SR2+
     h      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     i      )/re)+2*
     j      ((2*C5+4*c8*QA)*
     k      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     l      c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     m      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/
     n      re))/SR3+2*c9*
     o      (QB*(-2*QQ2*RRAC/(re*SR2)+
     p      2*(acx**2*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     q      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)+
     r      QQ3*(2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/
     s      (re*SR6)))+
     t      4*(acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)**2/
     u      re**2)+(C3+C5*QA+c8*QA**2)*
     v      (-2*QQ2*RRAC/(re*SR2)+
     w      2*(acx**2*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+ 
     x      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)+
     y      QQ3*(2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/
     z      (re*SR6)))+(c6+c10*QA)*
     A      (QQ3*(-6*acx**2*(1/(RAC**2*re**2*SR2**2)+
     B      QQ2*RRAC**3/(re*SR2))+ 
     C      6*(QQ2*RRAC/(re*SR2)+ 
     D      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)))+
     E      ((-3*QQ2**2+3*QQ3**2)*
     F      (2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/re-
     G      12*acx*QQ2*RRAC*(-2*abx*RRAB+acx*RRAC)/(re**2*SR2))/SR6))/ 
     H      DEXP(AN*QA)
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=TEMP
     1                                 + HESS(3*(J1-1)+J2,3*(J1-1)+J2)
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  Same atom, different component.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                IF (J4.NE.J1) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2) 
                QC=QQ3**3-3*QQ3*QQ2**2
               TEMP=TEMP+ 
     1 D*(((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+2*QA*(c2+c8*QB)+c10*QC- 
     2           AN*(c0+c4*QA**3+c7*QA**4+c3*QB+c9*QB**2+
     3     QA**2*(c2+c8*QB)+c6*QC+QA*(c1+c5*QB+c10*QC)))*
     4         (-abx*aby*RRAB**3/re-acx*acy*RRAC**3/re)+
     5        (-aby*RRAB/re-acy*RRAC/re)*
     6         (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+
     7        2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9    (2*(c3+c5*QA+c8*QA**2+2*c9*QB)*
     A     (acx*QQ2*RRAC/SR2+
     B       QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C       (c6+c10*QA)*
     D     (-6*acx*QQ2*QQ3*RRAC/SR2+
     E       (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/
     F        SR6))/re))+
     G     (2*(c5+2*c8*QA)*
     H      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     I     c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     J     (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/
     K      re))/SR3+(-(abx*RRAB/re)-acx*RRAC/re)*
     L   ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     M      AN**2*(c0+c4*QA**3+c7*QA**4+c3*QB+c9*QB**2+
     N      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+c5*QB+c10*QC)))* 
     O    (-(aby*RRAB/re)-acy*RRAC/re)/SR3**2+
     P     (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+
     Q     2*QA*(c2+c8*QB)+c10*QC)*
     R   (-(aby*RRAB/re)-acy*RRAC/re)/SR3+
     S    (2*(c3+c5*QA+c8*QA**2+2*c9*QB)*
     T     (acy*QQ2*RRAC/SR2+
     U    QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)+
     V    (c6+c10*QA)*
     W     (-6*acy*QQ2*QQ3*RRAC/SR2+
     X    (-3*QQ2**2+3*QQ3**2)*(-2*aby*RRAB+acy*RRAC)/
     Y     SR6))/re))+
     Z     (2*(c5+2*c8*QA)*
     a   (acy*QQ2*RRAC/SR2+QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)+
     b  c10*(-6*acy*QQ2*QQ3*RRAC/SR2+
     c  (-3*QQ2**2+3*QQ3**2)*(-2*aby*RRAB+acy*RRAC)/SR6))/
     d      re)/SR3)+2*(c9*
     e   (4*(acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     f    (acy*QQ2*RRAC/SR2+QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)/re**2 
     g   +2*QB*(acx*acy*
     h   (1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     i  (-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     j   (re**2*SR6**2)+
     k  QQ3*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/(re*SR6)))+
     l  (c3+c5*QA+c8*QA**2)*
     m   (acx*acy*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     n     (-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     o   (re**2*SR6**2)+
     p  QQ3*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/(re*SR6)))+
     q     (c6+c10*QA)*(-6*(acx*
     r   (acy*QQ3*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     s  QQ2*RRAC*(-2*aby*RRAB+acy*RRAC)/(re**2*SR2*SR6))+
     t  acy*QQ2*RRAC*(-2*abx*RRAB+acx*RRAC)/(re**2*SR2*SR6))+
     u  6*QQ3*(-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     v   (re**2*SR6**2)+
     w  (-3*QQ2**2+3*QQ3**2)*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/ 
     x   (re*SR6)))*DEXP(-AN*QA)
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=TEMP
     1                                 + HESS(3*(J1-1)+J5,3*(J1-1)+J2)
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  Different atoms, same component.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              IF ((J4.NE.J3).AND.(J4.NE.J1)) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=(RAB-RE)/RE
              RHOBC=(RBC-RE)/RE
              RHOAC=(RAC-RE)/RE
              QA=(RHOAB+RHOBC+RHOAC)/SR3
              QQ2=(RHOBC-RHOAC)/SR2
              QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
              QB=(QQ2**2)+(QQ3**2)
              QC=QQ3**3-3*QQ3*QQ2**2
              TEMP=TEMP+
     1        D*((C3+C5*QA+c8*QA**2+2*c9*QB)*
     2        (-2*acx*bcx*RRAC*RRBC/(re**2*SR2**2)+
     3        2*((-2*abx*RRAB+acx*RRAC)*(2*abx*RRAB+bcx*RRBC)/
     4        (re**2*SR6**2)+QQ3*(-2*RRAB+2*abx**2*RRAB**3)/(re*SR6)))
     5        +(((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*(c2+c8*QB)+
     6        c10*QC-AN*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     7        QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     8        *(-RRAB+abx**2*RRAB**3)+
     9        (abx*RRAB-bcx*RRBC)*
     A        (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     B        2*QA*(c2+c8*QB)+c10*QC)*
     C        (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     D        (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     E        (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     F        (c6+c10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     G        (-3*QQ2**2+3*QQ3**2)*
     H        (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     I        (2*(C5+2*c8*QA)*
     J        (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     K        c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     L        (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     M        )/re))/SR3+(-(abx*RRAB/re)-acx*RRAC/re)*
     N        ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     O        AN**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     P        QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     Q        )*(abx*RRAB-bcx*RRBC)/SR3**2+
     R        (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     S        2*QA*(c2+c8*QB)+c10*QC)*(abx*RRAB-bcx*RRBC)
     T        /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     U        (-2*bcx*QQ2*RRBC/SR2+
     V        2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)+
     W        (c6+c10*QA)*(6*bcx*QQ2*QQ3*RRBC/SR2+
     X        (-3*QQ2**2+3*QQ3**2)*(2*abx*RRAB+bcx*RRBC)/SR6
     Y        )))+(C5+2*c8*QA)*
     Z        (-2*bcx*QQ2*RRBC/SR2+2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)
     a        +c10*(6*bcx*QQ2*QQ3*RRBC/SR2+
     b        (-3*QQ2**2+3*QQ3**2)*(2*abx*RRAB+bcx*RRBC)/SR6))/
     c        SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     d        QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     e        (-2*bcx*QQ2*RRBC/SR2+2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)/re)/
     f        re+(c6+c10*QA)*(6*
     g        (QQ3*(acx*bcx*RRAC*RRBC/SR2**2+
     h        (-2*abx*RRAB+acx*RRAC)*(2*abx*RRAB+bcx*RRBC)/SR6**2)
     i        +bcx*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(SR2*SR6))/re**2
     j        +((-3*QQ2**2+3*QQ3**2)*(-2*RRAB+2*abx**2*RRAB**3)/re-
     k        6*acx*QQ2*RRAC*(2*abx*RRAB+bcx*RRBC)/(re**2*SR2))/SR6))/
     l        DEXP(AN*QA)

            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=TEMP
     1                                 + HESS(3*(J3-1)+J2,3*(J1-1)+J2)
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  Different atoms and different components
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((abx*aby*(c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*
     2      (c2+c8*QB)+
     3      c10*QC-AN*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     4      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     5      *RRAB**3+(aby*RRAB-bcy*RRBC)*
     6      (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     7      2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     A      (acx*QQ2*RRAC/SR2+
     B      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C      (c6+c10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     D      (-3*QQ2**2+3*QQ3**2)*
     E      (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     F      (2*(C5+2*c8*QA)*
     G      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     H      +c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     I      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     J      )/re))/SR3+
     K      (-(abx*RRAB/re)-acx*RRAC/re)*
     L      ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     M      AN**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     N      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     O      )*(aby*RRAB-bcy*RRBC)/SR3**2+
     P      (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     Q      2*QA*(c2+c8*QB)+c10*QC)*(aby*RRAB-bcy*RRBC)
     R      /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     S      (-2*bcy*QQ2*RRBC/SR2+
     T      2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)+
     U      (c6+c10*QA)*(6*bcy*QQ2*QQ3*RRBC/SR2+
     V      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6
     W      )))+(C5+2*c8*QA)*
     X      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)
     Y      +c10*(6*bcy*QQ2*QQ3*RRBC/SR2+
     Z      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6))/
     a      SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     b      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     c      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)/re)/ 
     d      re+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     e      ((-2*acx*bcy*RRAC*RRBC/SR2**2+
     f      2*(-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     g      re**2+4*abx*aby*QQ3*RRAB**3/(re*SR6))+
     h      (c6+c10*QA)*(6*(QQ3*
     i      (acx*bcy*RRAC*RRBC/SR2**2+
     j      (-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     k      re**2+(abx*aby*QQ3**2*RRAB**3+
     l      bcy*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(re*SR2))/(re*SR6))
     m      -6*(abx*aby*QQ2**2*RRAB**3+
     n      acx*QQ2*RRAC*(2*aby*RRAB+bcy*RRBC)/(re*SR2))/(re*SR6)))/ 
     o      DEXP(AN*QA)
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=TEMP
     1                                   +HESS(3*(J3-1)+J5,3*(J1-1)+J2)
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOAC+RHOBC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOAC-RHOBC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((abx*aby*(c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*
     2      (c2+c8*QB)+
     3      c10*QC-AN*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     4      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     5      *RRAB**3+(aby*RRAB-bcy*RRBC)*
     6      (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     7      2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     A      (acx*QQ2*RRAC/SR2+
     B      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C      (c6+c10*QA)*
     D      (-6*acx*QQ2*QQ3*RRAC/SR2+
     E      (-3*QQ2**2+3*QQ3**2)*
     F      (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     G      (2*(C5+2*c8*QA)*
     H      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     I      +c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     J      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     K      )/re))/SR3+
     L      (-(abx*RRAB/re)-acx*RRAC/re)*
     M      ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     O      AN**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     P      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     Q      )*(aby*RRAB-bcy*RRBC)/SR3**2+
     R      (-(AN*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     S      2*QA*(c2+c8*QB)+c10*QC)*(aby*RRAB-bcy*RRBC)
     T      /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     U      (-2*bcy*QQ2*RRBC/SR2+
     V      2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)+
     W      (c6+c10*QA)*(6*bcy*QQ2*QQ3*RRBC/SR2+
     X      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6
     Y      )))+(C5+2*c8*QA)*
     Z      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)
     a      +c10*(6*bcy*QQ2*QQ3*RRBC/SR2+
     b      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6))/
     c      SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     d      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     e      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)/re)/
     f      re+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     g      ((-2*acx*bcy*RRAC*RRBC/SR2**2+
     h      2*(-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     i      re**2+4*abx*aby*QQ3*RRAB**3/(re*SR6))+
     j      (c6+c10*QA)*(6*(QQ3*
     k      (acx*bcy*RRAC*RRBC/SR2**2+
     l      (-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     m      re**2+(abx*aby*QQ3**2*RRAB**3+
     n      bcy*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(re*SR2))/(re*SR6))
     o      -6*(abx*aby*QQ2**2*RRAB**3+
     p      acx*QQ2*RRAC*(2*aby*RRAB+bcy*RRBC)/(re*SR2))/(re*SR6)))/
     q      DEXP(AN*QA)
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=TEMP
     1                                   + HESS(3*(J3-1)+J5,3*(J1-1)+J2)
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  Symmetrise
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Here we calculate the analytic gradient and second derivatives
C  for the three-body term with periodic boundary conditions.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3C (N, X, V, CUTOFF, R2, RR2, VEC, RH, M, NDUM, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5
      INTEGER M(N,N,3)
      LOGICAL STEST
      DOUBLE PRECISION R2(N,N), RR2(N,N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     4                 QQ3, QB, QC, 
     5                 SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF
C
C
C  First the gradient.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/re)/
     9           DEXP(AN3*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
      IF (.NOT.STEST) RETURN
C
C  Diagonal bits of the Hessian.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 130 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=(RBC-RE)/RE
                  RHOAC=(RAC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ2=(RHOBC-RHOAC)/SR2
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QB=(QQ2**2)+(QQ3**2)
                  QC=(QQ3**3)-3*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1      D*((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     2      AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     3      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))*
     4      (-(abx*RRAB/re)-acx*RRAC/re)**2/SR3**2+
     5      ((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*(c2+c8*QB)+
     6      c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     7      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))*
     8      (RRAB-abx**2*RRAB**3+RRAC-acx**2*RRAC**3)/re+
     9      (-(abx*RRAB/re)-acx*RRAC/re)*
     a      (-2*AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     b      2*QA*(c2+c8*QB)+c10*QC)*
     c      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     d      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     e      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     f      +(c6+c10*QA)*
     g      (-6*acx*QQ2*QQ3*RRAC/SR2+
     h      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     i      )/re)+2*
     j      ((2*C5+4*c8*QA)*
     k      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     l      c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     m      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/
     n      re))/SR3+2*c9*
     o      (QB*(-2*QQ2*RRAC/(re*SR2)+
     p      2*(acx**2*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     q      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)+
     r      QQ3*(2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/
     s      (re*SR6)))+
     t      4*(acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)**2/
     u      re**2)+(C3+C5*QA+c8*QA**2)*
     v      (-2*QQ2*RRAC/(re*SR2)+
     w      2*(acx**2*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+ 
     x      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)+
     y      QQ3*(2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/
     z      (re*SR6)))+(c6+c10*QA)*
     A      (QQ3*(-6*acx**2*(1/(RAC**2*re**2*SR2**2)+
     B      QQ2*RRAC**3/(re*SR2))+ 
     C      6*(QQ2*RRAC/(re*SR2)+ 
     D      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)))+
     E      ((-3*QQ2**2+3*QQ3**2)*
     F      (2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/re-
     G      12*acx*QQ2*RRAC*(-2*abx*RRAB+acx*RRAC)/(re**2*SR2))/SR6))/ 
     H      DEXP(AN3*QA)
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=HESS(3*(J1-1)+J2,3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  Same atom, different component.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2) 
                QC=QQ3**3-3*QQ3*QQ2**2
               TEMP=TEMP+ 
     1 D*(((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+2*QA*(c2+c8*QB)+c10*QC- 
     2           AN3*(c0+c4*QA**3+c7*QA**4+c3*QB+c9*QB**2+
     3     QA**2*(c2+c8*QB)+c6*QC+QA*(c1+c5*QB+c10*QC)))*
     4         (-abx*aby*RRAB**3/re-acx*acy*RRAC**3/re)+
     5        (-aby*RRAB/re-acy*RRAC/re)*
     6         (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+
     7        2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9    (2*(c3+c5*QA+c8*QA**2+2*c9*QB)*
     A     (acx*QQ2*RRAC/SR2+
     B       QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C       (c6+c10*QA)*
     D     (-6*acx*QQ2*QQ3*RRAC/SR2+
     E       (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/
     F        SR6))/re))+
     G     (2*(c5+2*c8*QA)*
     H      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     I     c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     J     (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/
     K      re))/SR3+(-(abx*RRAB/re)-acx*RRAC/re)*
     L   ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     M      AN3**2*(c0+c4*QA**3+c7*QA**4+c3*QB+c9*QB**2+
     N      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+c5*QB+c10*QC)))* 
     O    (-(aby*RRAB/re)-acy*RRAC/re)/SR3**2+
     P     (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+
     Q     2*QA*(c2+c8*QB)+c10*QC)*
     R   (-(aby*RRAB/re)-acy*RRAC/re)/SR3+
     S    (2*(c3+c5*QA+c8*QA**2+2*c9*QB)*
     T     (acy*QQ2*RRAC/SR2+
     U    QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)+
     V    (c6+c10*QA)*
     W     (-6*acy*QQ2*QQ3*RRAC/SR2+
     X    (-3*QQ2**2+3*QQ3**2)*(-2*aby*RRAB+acy*RRAC)/
     Y     SR6))/re))+
     Z     (2*(c5+2*c8*QA)*
     a   (acy*QQ2*RRAC/SR2+QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)+
     b  c10*(-6*acy*QQ2*QQ3*RRAC/SR2+
     c  (-3*QQ2**2+3*QQ3**2)*(-2*aby*RRAB+acy*RRAC)/SR6))/
     d      re)/SR3)+2*(c9*
     e   (4*(acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     f    (acy*QQ2*RRAC/SR2+QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)/re**2 
     g   +2*QB*(acx*acy*
     h   (1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     i  (-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     j   (re**2*SR6**2)+
     k  QQ3*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/(re*SR6)))+
     l  (c3+c5*QA+c8*QA**2)*
     m   (acx*acy*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     n     (-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     o   (re**2*SR6**2)+
     p  QQ3*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/(re*SR6)))+
     q     (c6+c10*QA)*(-6*(acx*
     r   (acy*QQ3*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     s  QQ2*RRAC*(-2*aby*RRAB+acy*RRAC)/(re**2*SR2*SR6))+
     t  acy*QQ2*RRAC*(-2*abx*RRAB+acx*RRAC)/(re**2*SR2*SR6))+
     u  6*QQ3*(-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     v   (re**2*SR6**2)+
     w  (-3*QQ2**2+3*QQ3**2)*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/ 
     x   (re*SR6)))*DEXP(-AN3*QA)
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=HESS(3*(J1-1)+J5,3*(J1-1)+J2)+TEMP
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  Different atoms, same component.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1            R2(J4,J3).LT.CUTOFF.AND.
     2            R2(J4,J1).LT.CUTOFF.AND.
     3            R2(J3,J1).LT.CUTOFF) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=(RAB-RE)/RE
              RHOBC=(RBC-RE)/RE
              RHOAC=(RAC-RE)/RE
              QA=(RHOAB+RHOBC+RHOAC)/SR3
              QQ2=(RHOBC-RHOAC)/SR2
              QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
              QB=(QQ2**2)+(QQ3**2)
              QC=QQ3**3-3*QQ3*QQ2**2
              TEMP=TEMP+
     1        D*((C3+C5*QA+c8*QA**2+2*c9*QB)*
     2        (-2*acx*bcx*RRAC*RRBC/(re**2*SR2**2)+
     3        2*((-2*abx*RRAB+acx*RRAC)*(2*abx*RRAB+bcx*RRBC)/
     4        (re**2*SR6**2)+QQ3*(-2*RRAB+2*abx**2*RRAB**3)/(re*SR6)))
     5        +(((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*(c2+c8*QB)+
     6        c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     7        QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     8        *(-RRAB+abx**2*RRAB**3)+
     9        (abx*RRAB-bcx*RRBC)*
     A        (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     B        2*QA*(c2+c8*QB)+c10*QC)*
     C        (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     D        (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     E        (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     F        (c6+c10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     G        (-3*QQ2**2+3*QQ3**2)*
     H        (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     I        (2*(C5+2*c8*QA)*
     J        (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     K        c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     L        (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     M        )/re))/SR3+(-(abx*RRAB/re)-acx*RRAC/re)*
     N        ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     O        AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     P        QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     Q        )*(abx*RRAB-bcx*RRBC)/SR3**2+
     R        (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     S        2*QA*(c2+c8*QB)+c10*QC)*(abx*RRAB-bcx*RRBC)
     T        /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     U        (-2*bcx*QQ2*RRBC/SR2+
     V        2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)+
     W        (c6+c10*QA)*(6*bcx*QQ2*QQ3*RRBC/SR2+
     X        (-3*QQ2**2+3*QQ3**2)*(2*abx*RRAB+bcx*RRBC)/SR6
     Y        )))+(C5+2*c8*QA)*
     Z        (-2*bcx*QQ2*RRBC/SR2+2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)
     a        +c10*(6*bcx*QQ2*QQ3*RRBC/SR2+
     b        (-3*QQ2**2+3*QQ3**2)*(2*abx*RRAB+bcx*RRBC)/SR6))/
     c        SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     d        QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     e        (-2*bcx*QQ2*RRBC/SR2+2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)/re)/
     f        re+(c6+c10*QA)*(6*
     g        (QQ3*(acx*bcx*RRAC*RRBC/SR2**2+
     h        (-2*abx*RRAB+acx*RRAC)*(2*abx*RRAB+bcx*RRBC)/SR6**2)
     i        +bcx*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(SR2*SR6))/re**2
     j        +((-3*QQ2**2+3*QQ3**2)*(-2*RRAB+2*abx**2*RRAB**3)/re-
     k        6*acx*QQ2*RRAC*(2*abx*RRAB+bcx*RRBC)/(re**2*SR2))/SR6))/
     l        DEXP(AN3*QA)

            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=HESS(3*(J3-1)+J2,3*(J1-1)+J2)+TEMP
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  Different atoms and different components
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((abx*aby*(c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*
     2      (c2+c8*QB)+
     3      c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     4      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     5      *RRAB**3+(aby*RRAB-bcy*RRBC)*
     6      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     7      2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     A      (acx*QQ2*RRAC/SR2+
     B      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C      (c6+c10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     D      (-3*QQ2**2+3*QQ3**2)*
     E      (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     F      (2*(C5+2*c8*QA)*
     G      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     H      +c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     I      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     J      )/re))/SR3+
     K      (-(abx*RRAB/re)-acx*RRAC/re)*
     L      ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     M      AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     N      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     O      )*(aby*RRAB-bcy*RRBC)/SR3**2+
     P      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     Q      2*QA*(c2+c8*QB)+c10*QC)*(aby*RRAB-bcy*RRBC)
     R      /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     S      (-2*bcy*QQ2*RRBC/SR2+
     T      2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)+
     U      (c6+c10*QA)*(6*bcy*QQ2*QQ3*RRBC/SR2+
     V      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6
     W      )))+(C5+2*c8*QA)*
     X      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)
     Y      +c10*(6*bcy*QQ2*QQ3*RRBC/SR2+
     Z      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6))/
     a      SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     b      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     c      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)/re)/ 
     d      re+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     e      ((-2*acx*bcy*RRAC*RRBC/SR2**2+
     f      2*(-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     g      re**2+4*abx*aby*QQ3*RRAB**3/(re*SR6))+
     h      (c6+c10*QA)*(6*(QQ3*
     i      (acx*bcy*RRAC*RRBC/SR2**2+
     j      (-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     k      re**2+(abx*aby*QQ3**2*RRAB**3+
     l      bcy*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(re*SR2))/(re*SR6))
     m      -6*(abx*aby*QQ2**2*RRAB**3+
     n      acx*QQ2*RRAC*(2*aby*RRAB+bcy*RRBC)/(re*SR2))/(re*SR6)))/ 
     o      DEXP(AN3*QA)
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOAC+RHOBC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOAC-RHOBC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((abx*aby*(c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*
     2      (c2+c8*QB)+
     3      c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     4      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     5      *RRAB**3+(aby*RRAB-bcy*RRBC)*
     6      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     7      2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     A      (acx*QQ2*RRAC/SR2+
     B      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C      (c6+c10*QA)*
     D      (-6*acx*QQ2*QQ3*RRAC/SR2+
     E      (-3*QQ2**2+3*QQ3**2)*
     F      (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     G      (2*(C5+2*c8*QA)*
     H      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     I      +c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     J      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     K      )/re))/SR3+
     L      (-(abx*RRAB/re)-acx*RRAC/re)*
     M      ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     O      AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     P      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     Q      )*(aby*RRAB-bcy*RRBC)/SR3**2+
     R      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     S      2*QA*(c2+c8*QB)+c10*QC)*(aby*RRAB-bcy*RRBC)
     T      /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     U      (-2*bcy*QQ2*RRBC/SR2+
     V      2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)+
     W      (c6+c10*QA)*(6*bcy*QQ2*QQ3*RRBC/SR2+
     X      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6
     Y      )))+(C5+2*c8*QA)*
     Z      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)
     a      +c10*(6*bcy*QQ2*QQ3*RRBC/SR2+
     b      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6))/
     c      SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     d      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     e      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)/re)/
     f      re+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     g      ((-2*acx*bcy*RRAC*RRBC/SR2**2+
     h      2*(-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     i      re**2+4*abx*aby*QQ3*RRAB**3/(re*SR6))+
     j      (c6+c10*QA)*(6*(QQ3*
     k      (acx*bcy*RRAC*RRBC/SR2**2+
     l      (-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     m      re**2+(abx*aby*QQ3**2*RRAB**3+
     n      bcy*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(re*SR2))/(re*SR6))
     o      -6*(abx*aby*QQ2**2*RRAB**3+
     p      acx*QQ2*RRAC*(2*aby*RRAB+bcy*RRBC)/(re*SR2))/(re*SR6)))/
     q      DEXP(AN3*QA)
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  Symmetrise
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END

C*************************************************************************
C
C  Here we calculate the analytic gradient and second derivatives
C  for the three-body term 
C                                        
C*************************************************************************
C
      SUBROUTINE JM3CC(N, X, V,R2,RR2,VEC,RH,M,NDUM)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5
      INTEGER M(N,N,3)
      DOUBLE PRECISION R2(N,N), RR2(N,N), C,
     1              VEC(N,N,3), RH(N,N), 
     2              NDUM(N,N)
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC, V0, LAMBDA, SR2, SR3, SR6,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC, QQ2,
     4                 QQ3, QC
      DOUBLE PRECISION CUTOFF
      PARAMETER (SR2=1.414213562D0, SR3=1.732050808D0,SR6=2.449489743D0,
     1           C=1.0D0, 
     2           V0=60.0D0, LAMBDA=-4.0D0, CUTOFF=2.8D0)
C
C
C  First the gradient.
C 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=RH(J3,J1)
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(J4.NE.J3).AND.(R2(J4,J3).LT.CUTOFF)) THEN 
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=RH(J4,J1)
                  RHOBC=RH(J4,J3)
                  QQ3=(2*RAB-RBC-RAC)/SR6
                  QQ2=(RBC-RAC)/SR2
                  QC=QQ3**3-3.0D0*QQ3*QQ2**2
                  TEMP=TEMP+
     1   DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2   (C*(1.0D0+LAMBDA*QC)
     3    *(ABX*RRAB/RHOAB**2+ACX*RRAC/RHOAC**2) + 
     4   LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     5   (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6))*V0
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
C
C  Diagonal bits of the Hessian.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=RH(J3,J1)
               DO 130 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(J4.NE.J3).AND.(R2(J4,J3).LT.CUTOFF)) THEN 
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=RH(J4,J3)
                  RHOAC=RH(J4,J1)
                  QQ2=(RBC-RAC)/SR2
                  QQ3=(2*RAB-RBC-RAC)/SR6
                  QC=(QQ3**3)-3.0D0*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)**2 + 
     4  C*(2*(ABX**2/(RAB**2*RHOAB**3) + ACX**2/(RAC**2*RHOAC**3))+  
     5  (-RRAB + ABX**2*RRAB**3)/RHOAB**2 + 
     6  (-RRAC + ACX**2*RRAC**3)/RHOAC**2)) + 
     7  2*C*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     8  (LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     a  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     f  LAMBDA*(QQ3*(-6*ACX**2*(1/(RAC**2*SR2**2) + QQ2*RRAC**3/SR2) +   
     g  6*(QQ2*RRAC/SR2 + (-2*ABX*RRAB + ACX*RRAC)**2/SR6**2)) + 
     h  ((-3*QQ2**2 + 3*QQ3**2)*
     i  (2*RRAB - 2*ABX**2*RRAB**3 - RRAC + ACX**2*RRAC**3) - 
     j  12*ACX*QQ2*RRAC*(-2*ABX*RRAB + ACX*RRAC)/SR2)/SR6))*V0
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=HESS(3*(J1-1)+J2,3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  Same atom, different component.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(J4.NE.J3).AND.(R2(J4,J3).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=RH(J3,J1)
                RHOBC=RH(J4,J3)
                RHOAC=RH(J4,J1)
                QQ2=(RBC-RAC)/SR2
                QQ3=(2*RAB-RBC-RAC)/SR6
                QC=QQ3**3-3.0D0*QQ3*QQ2**2
               TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (c**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     4  (ABY*RRAB/RHOAB**2 + ACY*RRAC/RHOAC**2) + 
     5  C*(2*(ABX*ABY/(RAB**2*RHOAB**3) + ACX*ACY/(RAC**2*RHOAC**3)) + 
     6  ABX*ABY*RRAB**3/RHOAB**2 + ACX*ACY*RRAC**3/RHOAC**2)) + 
     7  C*((ABY*RRAB/RHOAB**2 + ACY*RRAC/RHOAC**2)*
     8  (  LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     a  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     b  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     d  (LAMBDA*(-6*ACY*QQ2*QQ3*RRAC/SR2 + 
     e  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABY*RRAB + ACY*RRAC)/SR6))) + 
     i  LAMBDA*(-6*(ACX*(ACY*QQ3*(1/(RAC**2*SR2**2) + QQ2*RRAC**3/SR2)+ 
     j  QQ2*RRAC*(-2*ABY*RRAB + ACY*RRAC)/(SR2*SR6)) + 
     k  ACY*QQ2*RRAC*(-2*ABX*RRAB + ACX*RRAC)/(SR2*SR6)) + 
     l  6*QQ3*(-2*ABX*RRAB + ACX*RRAC)*(-2*ABY*RRAB + ACY*RRAC)/SR6**2 + 
     m  (-3*QQ2**2+3*QQ3**2)*(-2*ABX*ABY*RRAB**3 + ACX*ACY*RRAC**3)/SR6 
     n  ))*V0
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=HESS(3*(J1-1)+J5,3*(J1-1)+J2)+TEMP
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  Different atoms, same component.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1            R2(J4,J3).LT.CUTOFF.AND.
     2            R2(J4,J1).LT.CUTOFF.AND.
     3            R2(J3,J1).LT.CUTOFF) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=RH(J3,J1)
              RHOBC=RH(J4,J3)
              RHOAC=RH(J4,J1)
              QQ2=(RBC-RAC)/SR2
              QQ3=(2*RAB-RBC-RAC)/SR6
              QC=QQ3**3-3.0D0*QQ3*QQ2**2
              TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (C*(RRAB/RHOAB**2 + abx**2*
     4  (-2/(RAB**2*RHOAB**3) - RRAB**3/RHOAB**2)) + 
     5  C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     6  (-(abx*RRAB/RHOAB**2) + BCX*RRBC/RHOBC**2)) + 
     7  C*((-(ABX*RRAB/RHOAB**2) + bcx*RRBC/RHOBC**2)*
     8  (  LAMBDA*(-6*acx*QQ2*QQ3*RRAC/SR2 + 
     A  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     B  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     D  (LAMBDA*(6*BCX*QQ2*QQ3*RRBC/SR2 + 
     E  (-3*QQ2**2 + 3*QQ3**2)*(2*abx*RRAB + BCX*RRBC)/SR6))) + 
     I  LAMBDA*(6*(QQ3*(ACX*BCX*RRAC*RRBC/SR2**2 + 
     J  (-2*ABX*RRAB + ACX*RRAC)*(2*ABX*RRAB + BCX*RRBC)/SR6**2) +  
     K  BCX*QQ2*(-2*ABX*RRAB + acx*RRAC)*RRBC/(SR2*SR6)) + 
     L  ((-3*QQ2**2 + 3*QQ3**2)*(-2*RRAB + 2*abx**2*RRAB**3) - 
     M  6*ACX*QQ2*RRAC*(2*ABX*RRAB + BCX*RRBC)/SR2)/SR6))*V0
            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=HESS(3*(J3-1)+J2,3*(J1-1)+J2)+TEMP
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  Different atoms and different components
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1
              TEMP=0.0D0
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 270 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             (R2(J3,J1).LT.CUTOFF).AND.
     2             (R2(J4,J3).LT.CUTOFF).AND.
     3             (R2(J4,J1).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=RH(J3,J1)
                RHOBC=RH(J4,J3)
                RHOAC=RH(J4,J1)
                QQ2=(RBC-RAC)/SR2
                QQ3=(2*RAB-RBC-RAC)/SR6
                QC=QQ3**3-3.0D0*QQ3*QQ2**2
                TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (ABX*ABY*C*(-2/(RAB**2*RHOAB**3) - RRAB**3/RHOAB**2) + 
     4  C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     5  (-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)) + 
     6  C*((-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)*
     8  (LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     9  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     a  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     c  (LAMBDA*(6*BCY*QQ2*QQ3*RRBC/SR2 + 
     d  (-3*QQ2**2 + 3*QQ3**2)*(2*ABY*RRAB + BCY*RRBC)/SR6))) + 
     h  LAMBDA*(6*(QQ3*(ACX*BCY*RRAC*RRBC/SR2**2 + 
     i  (-2*ABX*RRAB + ACX*RRAC)*(2*ABY*RRAB + BCY*RRBC)/SR6**2) + 
     j  (ABX*ABY*QQ3**2*RRAB**3 + 
     k  BCY*QQ2*(-2*ABX*RRAB + ACX*RRAC)*RRBC/SR2)/SR6) - 
     l  6*(ABX*ABY*QQ2**2*RRAB**3 + 
     m  ACX*QQ2*RRAC*(2*ABY*RRAB + BCY*RRBC)/SR2)/SR6))*V0
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP 
290         CONTINUE
            DO 295 J5=J2+1,3
              TEMP=0.0D0
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 275 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             (R2(J3,J1).LT.CUTOFF).AND.
     2             (R2(J4,J3).LT.CUTOFF).AND.
     3             (R2(J4,J1).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=RH(J3,J1)
                RHOBC=RH(J4,J3)
                RHOAC=RH(J4,J1)
                QQ2=(RBC-RAC)/SR2
                QQ3=(2*RAB-RAC-RBC)/SR6
                QC=QQ3**3-3.0D0*QQ3*QQ2**2
                TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (ABX*ABY*C*(-2/(RAB**2*RHOAB**3) - RRAB**3/RHOAB**2) + 
     4  C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     5  (-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)) + 
     6  C*((-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)*
     8  (LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     9  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     a  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     c  (LAMBDA*(6*BCY*QQ2*QQ3*RRBC/SR2 + 
     d  (-3*QQ2**2 + 3*QQ3**2)*(2*ABY*RRAB + BCY*RRBC)/SR6))) + 
     h  LAMBDA*(6*(QQ3*(ACX*BCY*RRAC*RRBC/SR2**2 + 
     i  (-2*ABX*RRAB + ACX*RRAC)*(2*ABY*RRAB + BCY*RRBC)/SR6**2) + 
     j  (ABX*ABY*QQ3**2*RRAB**3 + 
     k  BCY*QQ2*(-2*ABX*RRAB + ACX*RRAC)*RRBC/SR2)/SR6) - 
     l  6*(ABX*ABY*QQ2**2*RRAB**3 + 
     m  ACX*QQ2*RRAC*(2*ABY*RRAB + BCY*RRBC)/SR2)/SR6))*V0
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  Symmetrise
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Here we calculate the analytic gradient and second derivatives
C  for the three-body term with periodic boundary conditions.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3P (N, X, V, BOXLX, BOXLY, BOXLZ, CUTOFF, R2, RR2, VEC, M)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, NMX, NMY, NMZ, J3, J4, J5
      INTEGER M(N,N,3)
      DOUBLE PRECISION R2(N,N), RR2(N,N),
     1              VEC(N,N,3)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     4                 QQ3, QB, QC, 
     5                 SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3, BOXLX, BOXLY, BOXLZ
      DOUBLE PRECISION CUTOFF
C
C
C  First the gradient.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  NMX=M(J4,J3,1)+M(J3,J1,1)
                  NMY=M(J4,J3,2)+M(J3,J1,2)
                  NMZ=M(J4,J3,3)+M(J3,J1,3)
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF).AND.
     3                NMZ.EQ.M(J4,J1,3).AND.
     2                NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/re)/
     9           DEXP(AN3*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
C
C  Diagonal bits of the Hessian.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 130 J4=J3+1,N
                  NMX=M(J4,J3,1)+M(J3,J1,1)
                  NMY=M(J4,J3,2)+M(J3,J1,2)
                  NMZ=M(J4,J3,3)+M(J3,J1,3)
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF).AND.
     3                NMZ.EQ.M(J4,J1,3).AND.
     2                NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=(RBC-RE)/RE
                  RHOAC=(RAC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ2=(RHOBC-RHOAC)/SR2
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QB=(QQ2**2)+(QQ3**2)
                  QC=(QQ3**3)-3*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1      D*((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     2      AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     3      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))*
     4      (-(abx*RRAB/re)-acx*RRAC/re)**2/SR3**2+
     5      ((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*(c2+c8*QB)+
     6      c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     7      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))*
     8      (RRAB-abx**2*RRAB**3+RRAC-acx**2*RRAC**3)/re+
     9      (-(abx*RRAB/re)-acx*RRAC/re)*
     a      (-2*AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     b      2*QA*(c2+c8*QB)+c10*QC)*
     c      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     d      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     e      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     f      +(c6+c10*QA)*
     g      (-6*acx*QQ2*QQ3*RRAC/SR2+
     h      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     i      )/re)+2*
     j      ((2*C5+4*c8*QA)*
     k      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     l      c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     m      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/
     n      re))/SR3+2*c9*
     o      (QB*(-2*QQ2*RRAC/(re*SR2)+
     p      2*(acx**2*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     q      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)+
     r      QQ3*(2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/
     s      (re*SR6)))+
     t      4*(acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)**2/
     u      re**2)+(C3+C5*QA+c8*QA**2)*
     v      (-2*QQ2*RRAC/(re*SR2)+
     w      2*(acx**2*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+ 
     x      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)+
     y      QQ3*(2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/
     z      (re*SR6)))+(c6+c10*QA)*
     A      (QQ3*(-6*acx**2*(1/(RAC**2*re**2*SR2**2)+
     B      QQ2*RRAC**3/(re*SR2))+ 
     C      6*(QQ2*RRAC/(re*SR2)+ 
     D      (-2*abx*RRAB+acx*RRAC)**2/(re**2*SR6**2)))+
     E      ((-3*QQ2**2+3*QQ3**2)*
     F      (2*RRAB-2*abx**2*RRAB**3-RRAC+acx**2*RRAC**3)/re-
     G      12*acx*QQ2*RRAC*(-2*abx*RRAB+acx*RRAC)/(re**2*SR2))/SR6))/ 
     H      DEXP(AN3*QA)
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=HESS(3*(J1-1)+J2,3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  Same atom, different component.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                  NMX=M(J4,J3,1)+M(J3,J1,1)
                  NMY=M(J4,J3,2)+M(J3,J1,2)
                  NMZ=M(J4,J3,3)+M(J3,J1,3)
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF).AND.
     3                NMZ.EQ.M(J4,J1,3).AND.
     2                NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2) 
                QC=QQ3**3-3*QQ3*QQ2**2
               TEMP=TEMP+ 
     1 D*(((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+2*QA*(c2+c8*QB)+c10*QC- 
     2           AN3*(c0+c4*QA**3+c7*QA**4+c3*QB+c9*QB**2+
     3     QA**2*(c2+c8*QB)+c6*QC+QA*(c1+c5*QB+c10*QC)))*
     4         (-abx*aby*RRAB**3/re-acx*acy*RRAC**3/re)+
     5        (-aby*RRAB/re-acy*RRAC/re)*
     6         (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+
     7        2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9    (2*(c3+c5*QA+c8*QA**2+2*c9*QB)*
     A     (acx*QQ2*RRAC/SR2+
     B       QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C       (c6+c10*QA)*
     D     (-6*acx*QQ2*QQ3*RRAC/SR2+
     E       (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/
     F        SR6))/re))+
     G     (2*(c5+2*c8*QA)*
     H      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     I     c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     J     (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6))/
     K      re))/SR3+(-(abx*RRAB/re)-acx*RRAC/re)*
     L   ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     M      AN3**2*(c0+c4*QA**3+c7*QA**4+c3*QB+c9*QB**2+
     N      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+c5*QB+c10*QC)))* 
     O    (-(aby*RRAB/re)-acy*RRAC/re)/SR3**2+
     P     (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+c5*QB+
     Q     2*QA*(c2+c8*QB)+c10*QC)*
     R   (-(aby*RRAB/re)-acy*RRAC/re)/SR3+
     S    (2*(c3+c5*QA+c8*QA**2+2*c9*QB)*
     T     (acy*QQ2*RRAC/SR2+
     U    QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)+
     V    (c6+c10*QA)*
     W     (-6*acy*QQ2*QQ3*RRAC/SR2+
     X    (-3*QQ2**2+3*QQ3**2)*(-2*aby*RRAB+acy*RRAC)/
     Y     SR6))/re))+
     Z     (2*(c5+2*c8*QA)*
     a   (acy*QQ2*RRAC/SR2+QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)+
     b  c10*(-6*acy*QQ2*QQ3*RRAC/SR2+
     c  (-3*QQ2**2+3*QQ3**2)*(-2*aby*RRAB+acy*RRAC)/SR6))/
     d      re)/SR3)+2*(c9*
     e   (4*(acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     f    (acy*QQ2*RRAC/SR2+QQ3*(-2*aby*RRAB+acy*RRAC)/SR6)/re**2 
     g   +2*QB*(acx*acy*
     h   (1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     i  (-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     j   (re**2*SR6**2)+
     k  QQ3*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/(re*SR6)))+
     l  (c3+c5*QA+c8*QA**2)*
     m   (acx*acy*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     n     (-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     o   (re**2*SR6**2)+
     p  QQ3*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/(re*SR6)))+
     q     (c6+c10*QA)*(-6*(acx*
     r   (acy*QQ3*(1/(RAC**2*re**2*SR2**2)+QQ2*RRAC**3/(re*SR2))+
     s  QQ2*RRAC*(-2*aby*RRAB+acy*RRAC)/(re**2*SR2*SR6))+
     t  acy*QQ2*RRAC*(-2*abx*RRAB+acx*RRAC)/(re**2*SR2*SR6))+
     u  6*QQ3*(-2*abx*RRAB+acx*RRAC)*(-2*aby*RRAB+acy*RRAC)/
     v   (re**2*SR6**2)+
     w  (-3*QQ2**2+3*QQ3**2)*(-2*abx*aby*RRAB**3+acx*acy*RRAC**3)/ 
     x   (re*SR6)))*DEXP(-AN3*QA)
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=HESS(3*(J1-1)+J5,3*(J1-1)+J2)+TEMP
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  Different atoms, same component.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
               NMX=M(J4,J3,1)+M(J3,J1,1)
               NMY=M(J4,J3,2)+M(J3,J1,2)
               NMZ=M(J4,J3,3)+M(J3,J1,3)
              IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1            R2(J4,J3).LT.CUTOFF.AND.
     2            R2(J4,J1).LT.CUTOFF.AND.
     3            R2(J3,J1).LT.CUTOFF.AND.
     3            NMZ.EQ.M(J4,J1,3).AND.
     4            NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=(RAB-RE)/RE
              RHOBC=(RBC-RE)/RE
              RHOAC=(RAC-RE)/RE
              QA=(RHOAB+RHOBC+RHOAC)/SR3
              QQ2=(RHOBC-RHOAC)/SR2
              QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
              QB=(QQ2**2)+(QQ3**2)
              QC=QQ3**3-3*QQ3*QQ2**2
              TEMP=TEMP+
     1        D*((C3+C5*QA+c8*QA**2+2*c9*QB)*
     2        (-2*acx*bcx*RRAC*RRBC/(re**2*SR2**2)+
     3        2*((-2*abx*RRAB+acx*RRAC)*(2*abx*RRAB+bcx*RRBC)/
     4        (re**2*SR6**2)+QQ3*(-2*RRAB+2*abx**2*RRAB**3)/(re*SR6)))
     5        +(((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*(c2+c8*QB)+
     6        c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     7        QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     8        *(-RRAB+abx**2*RRAB**3)+
     9        (abx*RRAB-bcx*RRBC)*
     A        (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     B        2*QA*(c2+c8*QB)+c10*QC)*
     C        (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     D        (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     E        (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     F        (c6+c10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     G        (-3*QQ2**2+3*QQ3**2)*
     H        (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     I        (2*(C5+2*c8*QA)*
     J        (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     K        c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     L        (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     M        )/re))/SR3+(-(abx*RRAB/re)-acx*RRAC/re)*
     N        ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     O        AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     P        QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     Q        )*(abx*RRAB-bcx*RRBC)/SR3**2+
     R        (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     S        2*QA*(c2+c8*QB)+c10*QC)*(abx*RRAB-bcx*RRBC)
     T        /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     U        (-2*bcx*QQ2*RRBC/SR2+
     V        2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)+
     W        (c6+c10*QA)*(6*bcx*QQ2*QQ3*RRBC/SR2+
     X        (-3*QQ2**2+3*QQ3**2)*(2*abx*RRAB+bcx*RRBC)/SR6
     Y        )))+(C5+2*c8*QA)*
     Z        (-2*bcx*QQ2*RRBC/SR2+2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)
     a        +c10*(6*bcx*QQ2*QQ3*RRBC/SR2+
     b        (-3*QQ2**2+3*QQ3**2)*(2*abx*RRAB+bcx*RRBC)/SR6))/
     c        SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     d        QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     e        (-2*bcx*QQ2*RRBC/SR2+2*QQ3*(2*abx*RRAB+bcx*RRBC)/SR6)/re)/
     f        re+(c6+c10*QA)*(6*
     g        (QQ3*(acx*bcx*RRAC*RRBC/SR2**2+
     h        (-2*abx*RRAB+acx*RRAC)*(2*abx*RRAB+bcx*RRBC)/SR6**2)
     i        +bcx*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(SR2*SR6))/re**2
     j        +((-3*QQ2**2+3*QQ3**2)*(-2*RRAB+2*abx**2*RRAB**3)/re-
     k        6*acx*QQ2*RRAC*(2*abx*RRAB+bcx*RRBC)/(re**2*SR2))/SR6))/
     l        DEXP(AN3*QA)

            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=HESS(3*(J3-1)+J2,3*(J1-1)+J2)+TEMP
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  Different atoms and different components
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                NMX=M(J4,J3,1)+M(J3,J1,1)
                NMY=M(J4,J3,2)+M(J3,J1,2)
                NMZ=M(J4,J3,3)+M(J3,J1,3)
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF.AND.
     3             NMZ.EQ.M(J4,J1,3).AND.
     4             NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((abx*aby*(c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*
     2      (c2+c8*QB)+
     3      c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     4      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     5      *RRAB**3+(aby*RRAB-bcy*RRBC)*
     6      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     7      2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     A      (acx*QQ2*RRAC/SR2+
     B      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C      (c6+c10*QA)*(-6*acx*QQ2*QQ3*RRAC/SR2+
     D      (-3*QQ2**2+3*QQ3**2)*
     E      (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     F      (2*(C5+2*c8*QA)*
     G      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     H      +c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     I      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     J      )/re))/SR3+
     K      (-(abx*RRAB/re)-acx*RRAC/re)*
     L      ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     M      AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     N      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     O      )*(aby*RRAB-bcy*RRBC)/SR3**2+
     P      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     Q      2*QA*(c2+c8*QB)+c10*QC)*(aby*RRAB-bcy*RRBC)
     R      /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     S      (-2*bcy*QQ2*RRBC/SR2+
     T      2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)+
     U      (c6+c10*QA)*(6*bcy*QQ2*QQ3*RRBC/SR2+
     V      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6
     W      )))+(C5+2*c8*QA)*
     X      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)
     Y      +c10*(6*bcy*QQ2*QQ3*RRBC/SR2+
     Z      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6))/
     a      SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     b      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     c      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)/re)/ 
     d      re+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     e      ((-2*acx*bcy*RRAC*RRBC/SR2**2+
     f      2*(-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     g      re**2+4*abx*aby*QQ3*RRAB**3/(re*SR6))+
     h      (c6+c10*QA)*(6*(QQ3*
     i      (acx*bcy*RRAC*RRBC/SR2**2+
     j      (-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     k      re**2+(abx*aby*QQ3**2*RRAB**3+
     l      bcy*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(re*SR2))/(re*SR6))
     m      -6*(abx*aby*QQ2**2*RRAB**3+
     n      acx*QQ2*RRAC*(2*aby*RRAB+bcy*RRBC)/(re*SR2))/(re*SR6)))/ 
     o      DEXP(AN3*QA)
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                NMX=M(J4,J3,1)+M(J3,J1,1)
                NMY=M(J4,J3,2)+M(J3,J1,2)
                NMZ=M(J4,J3,3)+M(J3,J1,3)
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF.AND.
     3             NMZ.EQ.M(J4,J1,3).AND.
     4             NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOAC+RHOBC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOAC-RHOBC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((abx*aby*(c1+3*c4*QA**2+4*c7*QA**3+C5*QB+2*QA*
     2      (c2+c8*QB)+
     3      c10*QC-AN3*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     4      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC)))
     5      *RRAB**3+(aby*RRAB-bcy*RRBC)*
     6      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     7      2*QA*(c2+c8*QB)+c10*QC)*
     8      (-(abx*RRAB/re)-acx*RRAC/re)/SR3+
     9      (2*(C3+C5*QA+c8*QA**2+2*c9*QB)*
     A      (acx*QQ2*RRAC/SR2+
     B      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)+
     C      (c6+c10*QA)*
     D      (-6*acx*QQ2*QQ3*RRAC/SR2+
     E      (-3*QQ2**2+3*QQ3**2)*
     F      (-2*abx*RRAB+acx*RRAC)/SR6))/re))+
     G      (2*(C5+2*c8*QA)*
     H      (acx*QQ2*RRAC/SR2+QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)
     I      +c10*(-6*acx*QQ2*QQ3*RRAC/SR2+
     J      (-3*QQ2**2+3*QQ3**2)*(-2*abx*RRAB+acx*RRAC)/SR6)
     K      )/re))/SR3+
     L      (-(abx*RRAB/re)-acx*RRAC/re)*
     M      ((6*c4*QA+12*c7*QA**2+2*(c2+c8*QB)+
     O      AN3**2*(c0+c4*QA**3+c7*QA**4+C3*QB+c9*QB**2+
     P      QA**2*(c2+c8*QB)+c6*QC+QA*(c1+C5*QB+c10*QC))
     Q      )*(aby*RRAB-bcy*RRBC)/SR3**2+
     R      (-(AN3*((c1+3*c4*QA**2+4*c7*QA**3+C5*QB+
     S      2*QA*(c2+c8*QB)+c10*QC)*(aby*RRAB-bcy*RRBC)
     T      /SR3+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     U      (-2*bcy*QQ2*RRBC/SR2+
     V      2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)+
     W      (c6+c10*QA)*(6*bcy*QQ2*QQ3*RRBC/SR2+
     X      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6
     Y      )))+(C5+2*c8*QA)*
     Z      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)
     a      +c10*(6*bcy*QQ2*QQ3*RRBC/SR2+
     b      (-3*QQ2**2+3*QQ3**2)*(2*aby*RRAB+bcy*RRBC)/SR6))/
     c      SR3)+4*c9*(acx*QQ2*RRAC/SR2+
     d      QQ3*(-2*abx*RRAB+acx*RRAC)/SR6)*
     e      (-2*bcy*QQ2*RRBC/SR2+2*QQ3*(2*aby*RRAB+bcy*RRBC)/SR6)/re)/
     f      re+(C3+C5*QA+c8*QA**2+2*c9*QB)*
     g      ((-2*acx*bcy*RRAC*RRBC/SR2**2+
     h      2*(-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     i      re**2+4*abx*aby*QQ3*RRAB**3/(re*SR6))+
     j      (c6+c10*QA)*(6*(QQ3*
     k      (acx*bcy*RRAC*RRBC/SR2**2+
     l      (-2*abx*RRAB+acx*RRAC)*(2*aby*RRAB+bcy*RRBC)/SR6**2)/
     m      re**2+(abx*aby*QQ3**2*RRAB**3+
     n      bcy*QQ2*(-2*abx*RRAB+acx*RRAC)*RRBC/(re*SR2))/(re*SR6))
     o      -6*(abx*aby*QQ2**2*RRAB**3+
     p      acx*QQ2*RRAC*(2*aby*RRAB+bcy*RRBC)/(re*SR2))/(re*SR6)))/
     q      DEXP(AN3*QA)
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  Symmetrise
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  Here we calculate the potential using the two and three-body
C  terms of the Murrell potential
C                                        
C*************************************************************************
C
      SUBROUTINE JME(N, X, P2, P3, POTEL)
      IMPLICIT NONE 
      INTEGER N, J1, J2, I, J, IJ
      DOUBLE PRECISION SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN
      PARAMETER(SR2=1.414213562D0)
      PARAMETER(SR3=1.732050808D0)
      PARAMETER(SR6=2.449489743D0)
      PARAMETER(C0=3.598D0)
      PARAMETER(C1=-11.609D0)
      PARAMETER(C2=13.486D0)
      PARAMETER(C3=-18.174D0)
      PARAMETER(C4=-5.570D0)
      PARAMETER(C5=79.210D0)
      PARAMETER(C6=-6.458D0)
      PARAMETER(C7=23.383D0)
      PARAMETER(C8=-111.890D0)
      PARAMETER(C9=9.705D0)
      PARAMETER(C10=38.297D0)
      PARAMETER(RE=2.389D0)
      PARAMETER(D=2.918D0)
      PARAMETER(AN=6.50D0)
      DOUBLE PRECISION X(3*N), POTEL, P2, P3, RAB, RHOAB,
     1                 DIST(N,N), SDIST(N,N), RBC, 
     2                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3
C
      P2=0.0D0
      P3=0.0D0
C
C  It is vital to take as much as possible out of the loops,
C  especially the inner loop.
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            SDIST(J1,J2)=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1                   ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2                   ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            DIST(J1,J2)=DSQRT(SDIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  Calculate the energy
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J) THEN
               RAB=DIST(I,J)
               RHOAB=(RAB-RE)/RE
               P2=P2-D*(1+AN*RHOAB)*DEXP(-AN*RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J)) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      RETURN
      END
C
C*************************************************************************
C
C  Two- and three-body terms in the energy of the JM potential with
C  periodic boundary conditions. This must be called before jm2p
C  or jm3p as it does some setting up for them.
C                                        
C*************************************************************************
C
      SUBROUTINE JMEC(N,X,P2,P3,VNEW,POTEL,CUTOFF,GTEST,STEST)
      use porfuncs
      IMPLICIT NONE
      INTEGER N, J1, J2, I, J, IJ
      LOGICAL YESNO, GTEST,STEST
      INTEGER M(N,N,3)
      DOUBLE PRECISION DIST(N,N), SDIST(N,N), VNEW(3*N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N),
     2              SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB, RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3
C
      SR2=1.414213562D0
      SR3=1.732050808D0
      SR6=2.449489743D0
      P2=0.0D0
      P3=0.0D0
      CUTOFF=1.0D10
      INQUIRE(FILE='JMparams',EXIST=YESNO)
      IF (.NOT.YESNO) THEN
         PRINT*,'Data file JMparams not found - quit'
         STOP
      ELSE
         OPEN(UNIT=33,FILE='JMparams',STATUS='OLD')
         READ(33,*) C0
         READ(33,*) C1
         READ(33,*) C2
         READ(33,*) C3
         READ(33,*) C4
         READ(33,*) C5
         READ(33,*) C6
         READ(33,*) C7
         READ(33,*) C8
         READ(33,*) C9
         READ(33,*) C10
         READ(33,*) RE
         READ(33,*) D
         READ(33,*) AN2
         READ(33,*) AN3
         READ(33,*,END=666) CUTOFF
666      CLOSE(33)
      ENDIF
C
C  Calculation of connecting vectors
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  Calculation of distances:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  Calculate the energy
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB=DIST(I,J)
               RHOAB=(RAB-RE)/RE
               P2=P2-D*(1+AN2*RHOAB)*DEXP(-AN2*RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN3*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      IF (GTEST.OR.STEST) THEN
         CALL JM2C(N, X, VNEW, CUTOFF, DIST, SDIST, VEC, RH, M, NDUM,STEST)
         CALL JM3C(N, X, VNEW, CUTOFF, DIST, SDIST, VEC, RH, M, NDUM,STEST)
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  Two- and three-body terms in the energy of the JM potential with
C  periodic boundary conditions. This MUst be called before jm2p
C  or jm3p as it does some setting up for them.
C                                        
C*************************************************************************
C
      SUBROUTINE JMECC(N,X,P2,P3,VNEW,POTEL)
      IMPLICIT NONE
      INTEGER N, J1, J2, I, J, IJ
      INTEGER M(N,N,3)
      DOUBLE PRECISION DIST(N,N), SDIST(N,N), VNEW(3*N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N)
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB, RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QC, QQ3, V0, 
     2                 LAMBDA, SR2, SR6, ABIG, B, C, ASMALL
      PARAMETER (SR2=1.414213562D0, SR6=2.449489743D0,  
     1           ABIG=13.29D0, B=0.574D0, C=1.0D0, ASMALL=2.8D0,
     2           V0=60.0D0, LAMBDA=-4.0D0, CUTOFF=ASMALL)
C
      P2=0.0D0
      P3=0.0D0
C
C  Calculation of connecting vectors
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  Calculation of distances:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
            RH(J1,J2)=(DIST(J1,J2)-ASMALL)
            RH(J2,J1)=RH(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  Calculate the energy
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB=DIST(I,J)
               RHOAB=RH(I,J)
               P2=P2+ABIG*((B/RAB**4)-1.0D0)*DEXP(C/RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=RH(IJ,J)
                     RHOAC=RH(IJ,I)
                     QQ2=(RBC-RAC)/SR2
                     QQ3=(2*RAB-RBC-RAC)/SR6
                     QC=QQ3**3-3.0D0*QQ3*QQ2**2
                     P3=P3+(1.0D0+LAMBDA*QC)*V0
     1                 *DEXP(C*(1.0D0/RHOAB+1.0D0/RHOAC+1.0D0/RHOBC))
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3

      CALL JM2CC(N, X, VNEW,DIST,SDIST,VEC,RH,M,NDUM)
      CALL JM3CC(N, X, VNEW,DIST,SDIST,VEC,RH,M,NDUM)
      RETURN
      END
C
C*************************************************************************
C
C  Two- and three-body terms in the energy of the JM potential with
C  periodic boundary conditions. This must be called before jm2p
C  or jm3p as it does some setting up for them.
C                                        
C*************************************************************************
C

      SUBROUTINE JMEP(N,X,P2,P3,VNEW,POTEL,BOXLX,BOXLY,BOXLZ,CUTOFF)
      use porfuncs
      IMPLICIT NONE
      INTEGER N, J1, J2, I, J, IJ, NMX, NMY, NMZ
      LOGICAL YESNO
      INTEGER M(N,N,3)
      DOUBLE PRECISION DIST(N,N), SDIST(N,N),
     1              VEC(N,N,3), R2(N,N), RR2(N,N),
     2              SR2, SR3, SR6, C0, C1, C2, C3, C4, VNEW(3*N),
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB(N,N), RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3, BOXLX, BOXLY, BOXLZ
C
      P2=0.0D0
      P3=0.0D0
      SR2=1.414213562D0
      SR3=1.732050808D0
      SR6=2.449489743D0
      PRINT*,'WARNING - R2, RAB and RR2 were not set - not tested!'
      INQUIRE(FILE='JMparams',EXIST=YESNO)
      IF (.NOT.YESNO) THEN
         PRINT*,'Data file JMparams not found - quit'
         STOP
      ELSE
         OPEN(UNIT=33,FILE='JMparams',STATUS='OLD')
         READ(33,*) C0
         READ(33,*) C1
         READ(33,*) C2
         READ(33,*) C3
         READ(33,*) C4
         READ(33,*) C5
         READ(33,*) C6
         READ(33,*) C7
         READ(33,*) C8
         READ(33,*) C9
         READ(33,*) C10
         READ(33,*) RE
         READ(33,*) D
         READ(33,*) AN2
         READ(33,*) AN3
         CLOSE(33)
      ENDIF
C
C  atoms leaving the box in x or y direction enter again on the
C  other side.
C
      DO 21 J1=1,N
         X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX *
     1                  DNINT(X(3*(J1-1)+1)/BOXLX)
         X(3*(J1-1)+2)=X(3*(J1-1)+2) -  BOXLY *
     1                  DNINT(X(3*(J1-1)+2)/BOXLY)
         X(3*(J1-1)+3)=X(3*(J1-1)+3) -  BOXLZ *
     1                  DNINT(X(3*(J1-1)+3)/BOXLZ)
21    CONTINUE
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is
C  used:
C
      DO 25 J1=1,N
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         M(J1,J1,1)=0
         M(J1,J1,2)=0
         M(J1,J1,3)=0
         DO 15 J2=J1+1,N
            R2(J2,J1)=
     1                ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2                ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3                ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            R2(J2,J1)=DSQRT(R2(J2,J1))
            RR2(J2,J1)=1.0/R2(J2,J1)
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            M(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
            M(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
            M(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(M(J2,J1,1))
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(M(J2,J1,2))
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(M(J2,J1,3))
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
            M(J1,J2,1)=-M(J2,J1,1)
            M(J1,J2,2)=-M(J2,J1,2)
            M(J1,J2,3)=-M(J2,J1,3)
            RR2(J1,J2)=RR2(J2,J1) 
            R2(J1,J2)=R2(J2,J1) 
15       CONTINUE
25    CONTINUE
C
C  Calculation of distances:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  Calculate the energy
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB(I,J)=DIST(I,J)
               RHOAB=(RAB(I,J)-RE)/RE
               P2=P2-D*(1+AN2*RHOAB)*DEXP(-AN2*RHOAB)
               DO 24 IJ=1,N
                  NMX=M(I,J,1)+M(J,IJ,1)
                  NMY=M(I,J,2)+M(J,IJ,2)
                  NMZ=M(I,J,3)+M(J,IJ,3)
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF.AND.
     3                NMZ.EQ.M(I,IJ,3).AND.
     3                NMX.EQ.M(I,IJ,1).AND.NMY.EQ.M(I,IJ,2)) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN3*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      CALL JM2P(N,X,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF,RAB,VEC)
      CALL JM3P(N,X,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF, R2, RR2, VEC, M)
      RETURN
      END
