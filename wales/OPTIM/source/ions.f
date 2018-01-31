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
C  Shielded Born-Mayer potential. Ions must be enerted alternating
C  plus, minus, plus and there must be equal numbers of each sign.
C
C*************************************************************************
C
      SUBROUTINE IONS(N, X, V, ENERGY, GAMMA, CHARGE, RHO, APP, AMM, APM, NTEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTEST
!
! should really make the 2D arrays allocatble
!
      DOUBLE PRECISION X(3*N), ENERGY, GAMMA, CHARGE,
     1                 V(3*N), R(N,N), EXG(N,N), RHO, APP, 
     2                 R3(N,N), G(N,N), AMM, APM, RIJ, 
     3                 R2(N,N), EXR(N,N), F(N,N)

      IF (CHARGE.EQ.0.0D0) CHARGE=1.0D0
      IF (RHO.EQ.0.0D0) RHO=0.636848D0
      IF (APP.EQ.0.0D0) APP=57.153D0
      IF (AMM.EQ.0.0D0) AMM=70.735D0
      IF (APM.EQ.0.0D0) APM=65.667D0
C 
C  Store distance matrices.
C
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         R2(J1,J1)=0.0D0
         R3(J1,J1)=0.0D0
         EXR(J1,J1)=1.0D0
         EXG(J1,J1)=1.0D0
         DO 10 J2=J1+1,N
            R(J2,J1)=DSQRT( (X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                     +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                     +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2 )
            EXR(J2,J1)=DEXP(-R(J2,J1)/RHO)
            EXG(J2,J1)=DEXP(-R(J2,J1)*GAMMA)
            R2(J2,J1)=1.0D0/R(J2,J1)**2
            R3(J2,J1)=R2(J2,J1)/R(J2,J1)
            EXR(J1,J2)=EXR(J2,J1)
            EXG(J1,J2)=EXG(J2,J1)
            R(J1,J2)=R(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R3(J1,J2)=R3(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  Calculate the g and f tensors and the energy.
C
C  plus - plus
C
      ENERGY=0.0D0

      DO J1=1,N-1,2
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=J1+2,N-1,2 
            RIJ=R(J2,J1)
            G(J2,J1)=-APP*EXR(J2,J1)/(RHO*RIJ) - CHARGE**2*(1.0D0 + 
     1                GAMMA*RIJ)*EXG(J2,J1)*R3(J2,J1)
            F(J2,J1)=APP*EXR(J2,J1)/RHO**2 + CHARGE**2*(2.0D0 + 
     1                2.0D0*GAMMA*RIJ + (GAMMA*RIJ)**2)*EXG(J2,J1)
     2               *R3(J2,J1)-G(J2,J1)            
            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
            ENERGY=ENERGY+CHARGE**2*EXG(J2,J1)/RIJ + APP*EXR(J2,J1)
         ENDDO
      ENDDO
C
C  minus - minus
C
      DO J1=2,N,2
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=J1+2,N,2 
            RIJ=R(J2,J1)
            G(J2,J1)=-(AMM*EXR(J2,J1)/(RHO*RIJ)) - CHARGE**2*(1.0D0 
     1               + GAMMA*RIJ)*EXG(J2,J1)*R3(J2,J1)
            F(J2,J1)=AMM*EXR(J2,J1)/RHO**2 + CHARGE**2*(2.0D0 
     1 + 2.0D0*GAMMA*RIJ + (GAMMA*RIJ)**2)*EXG(J2,J1)*R3(J2,J1)-G(J2,J1)
            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
            ENERGY=ENERGY+CHARGE**2*EXG(J2,J1)/RIJ + AMM*EXR(J2,J1)
         ENDDO
      ENDDO
C
C  plus - minus
C
      DO J1=1,N-1,2
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=2,N,2 
            RIJ=R(J2,J1)
            G(J2,J1)=-(APM*EXR(J2,J1)/(RHO*RIJ)) + CHARGE**2*(1.0D0 
     1               + GAMMA*RIJ)*EXG(J2,J1)*R3(J2,J1)
            F(J2,J1)=APM*EXR(J2,J1)/RHO**2 - CHARGE**2*(2.0D0 + 
     1   2.0D0*GAMMA*RIJ + (GAMMA*RIJ)**2)*EXG(J2,J1)*R3(J2,J1)-G(J2,J1)
            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
            ENERGY=ENERGY-CHARGE**2*EXG(J2,J1)/RIJ + APM*EXR(J2,J1)
         ENDDO
      ENDDO

      IF (NTEST.EQ.0) RETURN
C
C  From here on down the code is system-independent!
C  First calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
               V(J3)=V(J3)+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
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
               HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
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
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3) + 
     1            F(J5,J1)*R2(J5,J1)* 
     2           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
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
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
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
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
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
