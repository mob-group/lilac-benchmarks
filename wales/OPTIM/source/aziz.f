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
C***********************************************************************
C
C  Subroutine aziz calculates the cartesian gradient and second
C  derivative matrix analytically for the accurate Ar-Ar pair potential
C  of Aziz (JCP, 99, 4518, 1993). Epsilon and R_e are set to unity.
C  Input coordinates are assumed to be in sigma=R_e/2**(1/6).
C
C***********************************************************************
C
      SUBROUTINE AZIZ(N, X, V, ENERGY, NTEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTEST
      DOUBLE PRECISION X(3*N), CONV, SR2, SR3, SR6, 
     1                 SR10, SR14,
     1              V(3*N), R2(N,N), ENERGY, R, R6(N,N),
     2              R8(N,N), G(N,N), R1(N,N), R10(N,N),
     3              R14(N,N), F(N,N), R12(N,N)
C     PARAMETER (CONV=2.0D0**(1.0D0/6.0D0), SR14=DSQRT(14.0D0), 
C    1           SR3=DSQRT(3.0D0))
C     PARAMETER (SR10=DSQRT(10.0D0), SR2=DSQRT(2.0D0), SR6=DSQRT(6.0D0))
      PARAMETER (SR10=3.16227766D0, SR2=1.414213562D0, SR6=2.44948743D0)
      PARAMETER (CONV=1.122462048D0, SR14=3.741657387D0, 
     1           SR3=1.732050808D0)
      DOUBLE PRECISION AA, ALPHA, BETA, C6, C8, C10, C12, C14, RHO, RM,
     1                 atob, z1, z2, z3, rc1, rc2, rc3, rc4, rc5
      PARAMETER (AA= 8.73933927D4, alpha=9.0328328D0, 
     1           beta=-2.37132823D0,
     2           c6=1.09309955D0, c8=0.51568309D0, c10=0.32521242D0, 
     3           c12=0.27818156D0, c14=0.31111959D0,
     4           rho=1.107D0, rm=3.757D0, atob=1.889726164D0, 
     5           z1=2.1D0, z2=0.109D0, z3=0.78D0)
C
C  Convert distances to R_e rather than sigma.
C
      DO J1=1,3*N
         X(J1)=X(J1)/CONV
      ENDDO
C 
C  Store distance matrices.
C
      DO 20 J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
         R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R1(J2,J1)=DSQRT(R2(J2,J1))
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R8(J2,J1)=R2(J2,J1)**4
            R6(J2,J1)=R8(J2,J1)/R2(J2,J1)
            R10(J2,J1)=R8(J2,J1)*R2(J2,J1)
            R12(J2,J1)=R10(J2,J1)*R2(J2,J1)
            R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)

            R1(J1,J2)=R1(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R6(J1,J2)=R6(J2,J1)
            R8(J1,J2)=R8(J2,J1)
            R10(J1,J2)=R10(J2,J1)
            R12(J1,J2)=R12(J2,J1)
            R14(J1,J2)=R14(J2,J1)
C
10       CONTINUE
20    CONTINUE
C
C  Calculate the energy in epsilon.
C
      ENERGY=0.0D0
      DO J1=1,N
         DO J2=J1+1,N
            R=R1(J2,J1)
            rc1=atob*R*rho*rm*z1
            rc2=atob**2*R**2*rho**2*rm**2*z2
            rc3=atob*R*rho*rm
            ENERGY=ENERGY + AA*DEXP(-(alpha*R) + beta*R**2) - 
     1             (1 - rc3**1.68/DEXP(rc3*z3))*
     2   (c10*(1 - DEXP(-rc1/10 - rc2/SR10))**10*R10(J2,J1) + 
     3   c12*(1 - DEXP(-rc1/12 - rc2/(2*SR3)))**12*R12(J2,J1) + 
     4    c14*(1 - DEXP(-rc1/14 - rc2/SR14))**14*R14(J2,J1) 
     5   + c6*(1 - DEXP(-rc1/6 - rc2/SR6))**6*R6(J2,J1) + 
     6     c8*(1 - DEXP(-rc1/8 - rc2/(2*SR2)))**8*R8(J2,J1))
         ENDDO
      ENDDO

      IF (NTEST.EQ.0) GOTO 11
C
C  Calculate the g and f tensors.
C
      DO 21 J1=1,N
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO 22 J2=J1+1,N 
            R=R1(J2,J1)
            rc1=atob*R*rho*rm*z1
            rc2=atob**2*R**2*rho**2*rm**2*z2
            rc3=atob*R*rho*rm
            rc4=atob*rho*rm*z1
            rc5=atob**2*R*rho**2*rm**2

            G(J2,J1)= (AA*DEXP(-(alpha*R) + beta*R**2)*(-alpha 
     %       + 2*beta*R) - 
     %                (1 - rc3**1.68/DEXP(rc3*z3))*
     %     (-14*c14*(1 - DEXP(-rc1/14 - rc2/SR14))**14/R**15 - 
     %       12*c12*(1 - DEXP(-rc1/12 - rc2/(2*SR3)))**12/R**13 - 
     %       10*c10*(1 - DEXP(-rc1/10 - rc2/SR10))**10/R**11 - 
     %       8*c8*(1 - DEXP(-rc1/8 - rc2/(2*SR2)))**8/R**9 - 
     %       6*c6*(1 - DEXP(-rc1/6 - rc2/SR6))**6/R**7 - 
     %       14*c14*DEXP(-rc1/14 - rc2/SR14)*(1 - DEXP(-rc1/14 
     %       - rc2/SR14))**13*
     %        R14(J2,J1)*(-rc4/14 - DSqrt(2.0D0/7.0D0)*rc5*z2) - 
     %       10*c10*DEXP(-rc1/10 - rc2/SR10)*(1 - DEXP(-rc1/10 
     %       - rc2/SR10))**9*
     %        R10(J2,J1)*(-rc4/10 - DSqrt(2.0D0/5.0D0)*rc5*z2) - 
     %       6*c6*DEXP(-rc1/6 - rc2/SR6)*(1 - DEXP(-rc1/6 
     %       - rc2/SR6))**5*R6(J2,J1)*
     %        (-rc4/6 - DSqrt(2.0D0/3.0D0)*rc5*z2) - 
     %       8*c8*DEXP(-rc1/8 - rc2/(2*SR2))*
     %        (1 - DEXP(-rc1/8 - rc2/(2*SR2)))**7*R8(J2,J1)*(-rc4/8 
     %        - rc5*z2/SR2) - 
     %       12*c12*DEXP(-rc1/12 - rc2/(2*SR3))*
     %        (1 - DEXP(-rc1/12 - rc2/(2*SR3)))**11*R12(J2,J1)
     %       *(-rc4/12 - rc5*z2/SR3))
     %      - (c10*(1 - DEXP(-rc1/10 - rc2/SR10))**10*R10(J2,J1) + 
     %       c12*(1 - DEXP(-rc1/12 - rc2/(2*SR3)))**12*R12(J2,J1) + 
     %       c14*(1 - DEXP(-rc1/14 - rc2/SR14))**14*R14(J2,J1) + 
     %       c6*(1 - DEXP(-rc1/6 - rc2/SR6))**6*R6(J2,J1) + 
     %       c8*(1 - DEXP(-rc1/8 - rc2/(2*SR2)))**8*R8(J2,J1))*
     %     (-1.68*atob*rc3**0.68*rho*rm/DEXP(rc3*z3) + 
     %       atob*rc3**1.68*rho*rm*z3/DEXP(rc3*z3)))/R

            F(J2,J1)= 2*AA*beta*DEXP(-(alpha*R) + beta*R**2) + 
     -  AA*DEXP(-(alpha*R) + beta*R**2)*(-alpha + 2*beta*R)**2 - 
     -  (1.0D0 - rc3**1.68/DEXP(rc3*z3))*
     -   (210*c14*(1.0D0 - DEXP(-rc1/14.0D0 - rc2/SR14))**14/R**16 + 
     -   72*c8*(1.0D0 - DEXP(-rc1/8.0D0 - rc2/(2*SR2)))**8*R10(J2,J1) + 
     -   110*c10*(1.0D0 - DEXP(-rc1/10.0D0 - rc2/SR10))**10*R12(J2,J1) +
     -  156*c12*(1.0D0-DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**12*R14(J2,J1) +
     -     42*c6*(1.0D0 - DEXP(-rc1/6.0D0 - rc2/SR6))**6*R8(J2,J1) + 
     -     2*SR10*atob**2*c10*DEXP(-rc1/10.0D0 - rc2/SR10)*
     -      (1.0D0 - DEXP(-rc1/10.0D0 - rc2/SR10))**9
     -    *rho**2*rm**2*R10(J2,J1)*z2 + 
     -     4*SR3*atob**2*c12*DEXP(-rc1/12.0D0 - rc2/(2*SR3))*
     -      (1.0D0 - DEXP(-rc1/12.0D0 - rc2/(2.0D0*SR3)))**11
     -   *rho**2*rm**2*R12(J2,J1)*z2 + 
     -     2*SR14*atob**2*c14*DEXP(-rc1/14.0D0 - rc2/SR14)*
     -      (1.0D0 - DEXP(-rc1/14.0D0 - rc2/SR14))**13*rho**2
     -    *rm**2*R14(J2,J1)*z2 + 
     -     2*SR6*atob**2*c6*DEXP(-rc1/6.0D0 - rc2/SR6)*
     -      (1.0D0 - DEXP(-rc1/6.0D0 - rc2/SR6))**5
     -    *rho**2*rm**2*R6(J2,J1)*z2 + 
     -     4*SR2*atob**2*c8*DEXP(-rc1/8.0D0 - rc2/(2*SR2))*
     -      (1.0D0 - DEXP(-rc1/8.0D0 - rc2/(2*SR2)))**7*rho**2
     -     *rm**2*R8(J2,J1)*z2 + 
     -     392*c14*DEXP(-rc1/14 - rc2/SR14)*(1 - DEXP(-rc1/14 - 
     -     rc2/SR14))**13*
     -       (-rc4/14.0D0 - DSqrt(2.0D0/7.0D0)*rc5*z2)/R**15 + 
     -     182*c14*DEXP(-rc1/7.0D0 - DSqrt(2.0D0/7.0D0)*rc2)*
     -      (1.0D0 - DEXP(-rc1/14.0D0 - rc2/SR14))**12*R14(J2,J1)*
     -      (-rc4/14.0D0 - DSqrt(2.0D0/7.0D0)*rc5*z2)**2 - 
     -     14*c14*DEXP(-rc1/14.0D0 - rc2/SR14)*(1.0D0 - DEXP(-rc1/14 
     -     - rc2/SR14))**13*
     -      R14(J2,J1)*(-rc4/14.0D0 - DSQRT(2.0D0/7.0D0)*rc5*z2)**2 + 
     -     200*c10*DEXP(-rc1/10.0D0 - rc2/SR10)*(1.0D0 - 
     -      DEXP(-rc1/10.0D0 - rc2/SR10))**9*
     -       (-rc4/10.0D0 - DSQRT(2.0D0/5.0D0)*rc5*z2)/R**11 + 
     -     90*c10*DEXP(-rc1/5.0D0 - DSQRT(2.0D0/5.0D0)*rc2)*
     -      (1.0D0 - DEXP(-rc1/10.0D0 - rc2/SR10))**8*R10(J2,J1)*
     -      (-rc4/10.0D0 - DSQRT(2.0D0/5.0D0)*rc5*z2)**2 - 
     -     10*c10*DEXP(-rc1/10.0D0 - rc2/SR10)*(1.0D0 - 
     -     DEXP(-rc1/10.0D0 - rc2/SR10))**9*
     -      R10(J2,J1)*(-rc4/10.0D0 - DSQRT(2.0D0/5.0D0)*rc5*z2)**2 + 
     -     72*c6*DEXP(-rc1/6.0D0 - rc2/SR6)*(1.0D0 - DEXP(-rc1/6.0D0 
     -      - rc2/SR6))**5*
     -       (-rc4/6.0D0 - DSQRT(2.0D0/3.0D0)*rc5*z2)/R**7 + 
     -     30*c6*DEXP(-rc1/3.0D0 - DSQRT(2.0D0/3.0D0)*rc2)*(1.0D0 
     -     - DEXP(-rc1/6.0D0 - rc2/SR6))**4*
     -      R6(J2,J1)*(-rc4/6.0D0 - DSQRT(2.0D0/3.0D0)*rc5*z2)**2 - 
     -     6*c6*DEXP(-rc1/6.0D0 - rc2/SR6)*(1.0D0 - DEXP(-rc1/6 
     -     - rc2/SR6))**5*R6(J2,J1)*
     -      (-rc4/6.0D0 - DSQRT(2.0D0/3.0D0)*rc5*z2)**2 + 
     -     128*c8*DEXP(-rc1/8.0D0 - rc2/(2*SR2))*
     -       (1.0D0 - DEXP(-rc1/8.0D0 - rc2/(2*SR2)))**7*(-rc4/8 
     -      - rc5*z2/SR2)/R**9 + 
     -     56*c8*DEXP(-rc1/4.0D0 - rc2/SR2)*(1.0D0 - DEXP(-rc1/8 
     -    - rc2/(2*SR2)))**6*
     -      R8(J2,J1)*(-rc4/8.0D0 - rc5*z2/SR2)**2 - 
     -     8*c8*DEXP(-rc1/8.0D0 - rc2/(2*SR2))*(1.0D0 - DEXP(-rc1/8 
     -     - rc2/(2*SR2)))**7*
     -      R8(J2,J1)*(-rc4/8.0D0 - rc5*z2/SR2)**2 + 
     -     288*c12*DEXP(-rc1/12.0D0 - rc2/(2*SR3))*
     -       (1.0D0 - DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**11*(-rc4/12 
     -     - rc5*z2/SR3)/R**13
     -       + 132*c12*DEXP(-rc1/6.0D0 - rc2/SR3)*
     -      (1.0D0 - DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**10*R12(J2,J1)*
     -      (-rc4/12.0D0 - rc5*z2/SR3)**2 - 
     -     12*c12*DEXP(-rc1/12.0D0 - rc2/(2*SR3))*
     -      (1.0D0 - DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**11*R12(J2,J1)*
     -      (-rc4/12.0D0 - rc5*z2/SR3)**2) - 
     -  2*(-14*c14*(1.0D0 - DEXP(-rc1/14.0D0 - rc2/SR14))**14/R**15 - 
     -     12*c12*(1.0D0 - DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**12/R**13 - 
     -     10*c10*(1.0D0 - DEXP(-rc1/10.0D0 - rc2/SR10))**10/R**11 - 
     -     8*c8*(1.0D0 - DEXP(-rc1/8.0D0 - rc2/(2*SR2)))**8/R**9 - 
     -     6*c6*(1.0D0 - DEXP(-rc1/6.0D0 - rc2/SR6))**6/R**7 - 
     -     14*c14*DEXP(-rc1/14.0D0 - rc2/SR14)*(1.0D0 - DEXP(-rc1/14 
     -      - rc2/SR14))**13*
     -      R14(J2,J1)*(-rc4/14.0D0 - DSQRT(2.0D0/7.0D0)*rc5*z2) - 
     -     10*c10*DEXP(-rc1/10.0D0 - rc2/SR10)*(1.0D0 - DEXP(-rc1/10 
     -     - rc2/SR10))**9*
     -      R10(J2,J1)*(-rc4/10.0D0 - DSQRT(2.0D0/5.0D0)*rc5*z2) - 
     -     6*c6*DEXP(-rc1/6.0D0 - rc2/SR6)*(1.0D0 - DEXP(-rc1/6 - 
     -     rc2/SR6))**5*R6(J2,J1)*
     -      (-rc4/6.0D0 - DSQRT(2.0D0/3.0D0)*rc5*z2) - 
     -     8*c8*DEXP(-rc1/8.0D0 - rc2/(2*SR2))*(1.0D0 - DEXP(-rc1/8 
     -     - rc2/(2*SR2)))**7*
     -      R8(J2,J1)*(-rc4/8.0D0 - rc5*z2/SR2) - 
     -     12*c12*DEXP(-rc1/12.0D0 - rc2/(2*SR3))*
     -      (1.0D0 - DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**11*R12(J2,J1)*
     -      (-rc4/12.0D0 - rc5*z2/SR3))*
     -   (-1.68*atob*rc3**0.68*rho*rm/DEXP(rc3*z3) + 
     -     atob*rc3**1.68*rho*rm*z3/DEXP(rc3*z3)) - 
     -  (c10*(1.0D0 - DEXP(-rc1/10.0D0 - rc2/SR10))**10*R10(J2,J1) + 
     -     c12*(1.0D0-DEXP(-rc1/12.0D0 - rc2/(2*SR3)))**12*R12(J2,J1) + 
     -     c14*(1.0D0 - DEXP(-rc1/14.0D0 - rc2/SR14))**14*R14(J2,J1) + 
     -     c6*(1.0D0 - DEXP(-rc1/6.0D0 - rc2/SR6))**6*R6(J2,J1) + 
     -     c8*(1.0D0 - DEXP(-rc1/8.0D0 - rc2/(2*SR2)))**8*R8(J2,J1))*
     -   (-1.1424*atob**2*rho**2*rm**2/
     -      (DEXP(rc3*z3)*rc3**0.32) + 
     -     3.36*atob**2*rc3**0.68*rho**2*rm**2*z3/DEXP(rc3*z3) - 
     -     atob**2*rc3**1.68*rho**2*rm**2*z3**2/DEXP(rc3*z3)) 

            F(J2,J1)=F(J2,J1)-G(J2,J1)

            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
22       CONTINUE
21    CONTINUE
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
         V(J1)=V(J1)/CONV
         HESS(J1,J1)=HESS(J1,J1)/(CONV**2)
         DO 190 J2=J1+1,3*N
            HESS(J2,J1)=HESS(J2,J1)/(CONV**2)
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
C
C  Convert distances back to sigma.
C
11    DO J1=1,3*N
          X(J1)=X(J1)*CONV
      ENDDO
      RETURN
      END
