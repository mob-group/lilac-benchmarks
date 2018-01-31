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
C**************************************************************************
C
C  Subroutine ROTD calculates analytically the cartesian gradient and
C  second derivative matrix due to the ROTATIONAL TERM ONLY and increments
C  the existing matrices containing the potential contributions 
C  accordingly.
C
C**************************************************************************
C
      SUBROUTINE ROTD(N, X, V, MASS, ZJ, FLAG, ROTE)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION X(3*N), ROTE,
     1                 V(3*N), MASS, ZJ, IZINV, XA, XB, YA, YB, ZA,ZB,  
     2                 A, B, C, D, E, F, TOP, BOTTOM, XM, YM, ZM, XOM
      LOGICAL FLAG
      XOM=1.0D0-1.0D0/FLOAT(N)
C
C  Calculate centre of mass components XM, YM, ZM
C
      XM=0.0D0
      YM=0.0D0
      ZM=0.0D0
      DO J1=1,N
         XM=XM+X(3*(J1-1)+1)
         YM=YM+X(3*(J1-1)+2)
         ZM=ZM+X(3*(J1-1)+3)
      ENDDO
      XM=XM/N
      YM=YM/N
      ZM=ZM/N
C
C  Calculate constants A, B, C, D, E, F
C
      A=0.0D0
      B=0.0D0
      C=0.0D0
      D=0.0D0
      E=0.0D0
      F=0.0D0
      DO J1=1,N
         J3=3*(J1-1)
         A=A+(X(J3+2)-YM)**2
         B=B+(X(J3+3)-ZM)**2
         E=E+(X(J3+1)-XM)**2
         C=C+(X(J3+2)-YM)*(X(J3+1)-XM)
         D=D+(X(J3+3)-ZM)*(X(J3+1)-XM)
         F=F+(X(J3+3)-ZM)*(X(J3+2)-YM)
      ENDDO
      TOP=A*B+B*B-C*C+A*E+B*E
      BOTTOM=A*A*B+A*B*B-A*C*C-B*D*D+A*A*E+2.0D0*A*B*E+B*B*E-C*C*E-D*D*E
     1      +A*E*E+B*E*E-2.0D0*C*D*F-A*F*F-B*F*F
      IF (DABS(BOTTOM).LT.1.0D-6) THEN
         PRINT*,'PANIC CALLED, BOTTOM=',BOTTOM
         CALL EPANIC(X,N,ZJ,ROTE,FLAG,V,MASS)
         RETURN
      ELSE
         IZINV=TOP/(MASS*BOTTOM)
         ROTE=ZJ**2*IZINV/2.0D0
      ENDIF

      IF (FLAG) RETURN
C
C  Analytic gradient:
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         V(J3+1)=V(J3+1)+ZJ**2*(2*(-(a*XM)-b*XM+a*xa+b*xa+c*YM-c*ya)
     &    /bottom +
     &    2*top*(a**2*XM+2*a*b*XM+b**2*XM-c**2*XM-d**2*XM+2*a*e*XM+  
     &    2*b*e*XM-a**2*xa-2*a*b*xa-b**2*xa+c**2*xa+d**2*xa -
     &    2*a*e*xa-2*b*e*xa-a*c*YM-c*e*YM-d*f*YM+a*c*ya+c*e*ya +
     &    d*f*ya-b*d*ZM-d*e*ZM-c*f*ZM+b*d*za+d*e*za+c*f*za)/
     &    bottom**2) / (2.0D0*MASS) 

         V(J3+2)=V(J3+2)+ZJ**2*(
     &    2*(c*XM-c*xa-b*YM-e*YM+b*ya+e*ya)/bottom -
     &    2*top*(a*c*XM+c*e*XM+d*f*XM-a*c*xa-c*e*xa-d*f*xa-2*a*b*YM- 
     &    b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-e**2*YM+f**2*YM +
     &    2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+e**2*ya -
     &    f**2*ya+c*d*ZM+a*f*ZM+b*f*ZM-c*d*za-a*f*za-b*f*za)/
     &    bottom**2 ) / (2.0D0*MASS)

         V(J3+3)=V(J3+3)+ZJ**2*(
     &    2*(a+2*b+e)*(-ZM+za)/bottom -
     &    2*top*(b*d*XM+d*e*XM+c*f*XM-b*d*xa-d*e*xa-c*f*xa+c*d*YM +
     &    a*f*YM+b*f*YM-c*d*ya-a*f*ya-b*f*ya-a**2*ZM-2*a*b*ZM +
     &    d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM+a**2*za +
     &    2*a*b*za-d**2*za+2*a*e*za+2*b*e*za+e**2*za-f**2*za)/
     &    bottom**2) / (2.0D0*MASS)

      ENDDO
C
C  Case I: diagonal terms.
C
      DO J1 = 1, N
         J3 = 3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         HESS(J3+1, J3+1) = HESS(J3+1, J3+1)+ZJ**2 * (
     &  (2*a+2*b-2*a/N-2*b/N-2*ya**2+4*ya*YM-2*YM**2)/bottom- 
     &  top*(2*a**2*XOM+4*a*b*XOM+
     &      2*b**2*XOM-2*c**2*XOM-
     &      2*d**2*XOM+4*a*e*XOM+
     &      4*b*e*XOM+8*a*(xa-XM)**2+8*b*(xa-XM)**2-
     &      8*c*(xa-XM)*(ya-YM)-2*a*(ya-YM)**2-2*e*(ya-YM)**2-
     &      8*d*(xa-XM)*(za-ZM)-4*f*(ya-YM)*(za-ZM)-
     &      2*b*(za-ZM)**2-2*e*(za-ZM)**2)/bottom**2-
     &  8*(a*xa+b*xa-a*XM-b*XM-c*ya+c*YM)*
     &    (a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     &      2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     &      2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     &      c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     &      c*f*ZM)/bottom**2+
     &  8*top*(a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     &       2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     &       2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     &       c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     &       c*f*ZM)**2/bottom**3)/(2.0D0*MASS)

         HESS(J3+2, J3+2) = HESS(J3+2, J3+2)+ZJ**2 * (
     &   (2*b+2*e-2*b/N-2*e/N-2*xa**2+4*xa*XM-2*XM**2)/bottom-   
     &  top*(4*a*b*XOM+2*b**2*XOM-
     &      2*c**2*XOM+4*a*e*XOM+
     &      4*b*e*XOM+2*e**2*XOM-
     &      2*f**2*XOM-2*a*(xa-XM)**2-2*e*(xa-XM)**2-
     &      8*c*(xa-XM)*(ya-YM)+8*b*(ya-YM)**2+8*e*(ya-YM)**2-
     &      4*d*(xa-XM)*(za-ZM)-8*f*(ya-YM)*(za-ZM)-
     &      2*a*(za-ZM)**2-2*b*(za-ZM)**2)/bottom**2+
     &  8*(c*xa-c*XM-b*ya-e*ya+b*YM+e*YM)*
     &    (-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+e**2*ya-
     &      f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za+c*d*ZM+a*f*ZM+
     &      b*f*ZM)/bottom**2+
     &  8*top*(-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM+d*f*XM+
     &       2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+e**2*ya-
     &       f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &       e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za+c*d*ZM+a*f*ZM+
     &       b*f*ZM)**2/bottom**3)/(2.0D0*MASS)

        HESS(J3+3, J3+3) = HESS(J3+3, J3+3)+ZJ**2 * (
     &   -(top*(2*a**2*XOM+4*a*b*XOM-2*d**2*XOM+
     &        4*a*e*XOM+4*b*e*XOM+
     &        2*e**2*XOM-2*f**2*XOM-
     &        2*b*(xa-XM)**2-2*e*(xa-XM)**2-4*c*(xa-XM)*(ya-YM)-
     &        2*a*(ya-YM)**2-2*b*(ya-YM)**2-8*d*(xa-XM)*(za-ZM)-
     &        8*f*(ya-YM)*(za-ZM)+8*a*(za-ZM)**2+8*e*(za-ZM)**2)/
     &     bottom**2)+8*top*(b*d*xa+d*e*xa+c*f*xa-b*d*XM-d*e*XM-
     &       c*f*XM+c*d*ya+a*f*ya+b*f*ya-c*d*YM-a*f*YM-b*f*YM-
     &       a**2*za-2*a*b*za+d**2*za-2*a*e*za-2*b*e*za-e**2*za+
     &       f**2*za+a**2*ZM+2*a*b*ZM-d**2*ZM+2*a*e*ZM+2*b*e*ZM+
     &       e**2*ZM-f**2*ZM)**2/bottom**3-
     &  8*(a+2*b+e)*(za-ZM)*
     &    (-(b*d*xa)-d*e*xa-c*f*xa+b*d*XM+d*e*XM+c*f*XM-c*d*ya-
     &      a*f*ya-b*f*ya+c*d*YM+a*f*YM+b*f*YM+a**2*za+2*a*b*za-
     &      d**2*za+2*a*e*za+2*b*e*za+e**2*za-f**2*za-a**2*ZM-
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**2+(2*a+4*b+2*e-2*a/N-4*b/N-2*e/N+
     &     8*za**2-16*za*ZM+8*ZM**2)/bottom)/(2.0D0*MASS)

      ENDDO
C
C  Case II: same Cartesian coordinate, different atom.
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         DO J2=J1+1,N
            J4=3*(J2-1)
            XB=X(J4+1)
            YB=X(J4+2)
            ZB=X(J4+3)
            HESS(J3+1,J4+1)=HESS(J3+1,J4+1)+ZJ**2*(
     1      (-2*a/N-2*b/N-2*(ya-YM)*(yb-YM))/bottom-
     1     top*(-2*a**2/N-4*a*b/N-2*b**2/N+2*c**2/N+
     1     2*d**2/N-4*a*e/N-4*b*e/N+
     1     8*a*(xa-XM)*(xb-XM)+8*b*(xa-XM)*(xb-XM)-
     1     4*c*(xb-XM)*(ya-YM)-4*c*(xa-XM)*(yb-YM)-
     1     2*a*(ya-YM)*(yb-YM)-2*e*(ya-YM)*(yb-YM)-
     1     4*d*(xb-XM)*(za-ZM)-2*f*(yb-YM)*(za-ZM)-
     1     4*d*(xa-XM)*(zb-ZM)-2*f*(ya-YM)*(zb-ZM)-
     1     2*b*(za-ZM)*(zb-ZM)-2*e*(za-ZM)*(zb-ZM))/bottom**2-
     1     4*(a*xb+b*xb-a*XM-b*XM-c*yb+c*YM)*
     1     (a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     1     2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     1     2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     1     c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     1     c*f*ZM)/bottom**2-
     1     4*(a*xa+b*xa-a*XM-b*XM-c*ya+c*YM)*
     1     (a**2*xb+2*a*b*xb+b**2*xb-c**2*xb-d**2*xb+2*a*e*xb+
     1     2*b*e*xb-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     1     2*a*e*XM-2*b*e*XM-a*c*yb-c*e*yb-d*f*yb+a*c*YM+
     1     c*e*YM+d*f*YM-b*d*zb-d*e*zb-c*f*zb+b*d*ZM+d*e*ZM+
     1     c*f*ZM)/bottom**2+
     1     8*top*(a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     1     2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     1     2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     1     c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     1     c*f*ZM)*(a**2*xb+2*a*b*xb+b**2*xb-c**2*xb-d**2*xb+
     1     2*a*e*xb+2*b*e*xb-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+
     1     d**2*XM-2*a*e*XM-2*b*e*XM-a*c*yb-c*e*yb-d*f*yb+
     1     a*c*YM+c*e*YM+d*f*YM-b*d*zb-d*e*zb-c*f*zb+b*d*ZM+
     1     d*e*ZM+c*f*ZM)/bottom**3)/(2.0D0*MASS)
            HESS(J4+1,J3+1)=HESS(J3+1,J4+1)

            HESS(J3+2,J4+2)=HESS(J3+2,J4+2)+ZJ**2*(
     &      (-2*b/N-2*e/N-2*(xa-XM)*(xb-XM))/bottom-
     &  top*(-4*a*b/N-2*b**2/N+2*c**2/N-4*a*e/N-
     &      4*b*e/N-2*e**2/N+2*f**2/N-
     &      2*a*(xa-XM)*(xb-XM)-2*e*(xa-XM)*(xb-XM)-
     &      4*c*(xb-XM)*(ya-YM)-4*c*(xa-XM)*(yb-YM)+
     &      8*b*(ya-YM)*(yb-YM)+8*e*(ya-YM)*(yb-YM)-
     &      2*d*(xb-XM)*(za-ZM)-4*f*(yb-YM)*(za-ZM)-
     &      2*d*(xa-XM)*(zb-ZM)-4*f*(ya-YM)*(zb-ZM)-
     &      2*a*(za-ZM)*(zb-ZM)-2*b*(za-ZM)*(zb-ZM))/bottom**2+
     &  8*top*(a*c*xa+c*e*xa+d*f*xa-a*c*XM-c*e*XM-d*f*XM-
     &      2*a*b*ya-b**2*ya+c**2*ya-2*a*e*ya-2*b*e*ya-e**2*ya+
     &      f**2*ya+2*a*b*YM+b**2*YM-c**2*YM+2*a*e*YM+2*b*e*YM+
     &      e**2*YM-f**2*YM+c*d*za+a*f*za+b*f*za-c*d*ZM-a*f*ZM-
     &      b*f*ZM)*(a*c*xb+c*e*xb+d*f*xb-a*c*XM-c*e*XM-d*f*XM-
     &      2*a*b*yb-b**2*yb+c**2*yb-2*a*e*yb-2*b*e*yb-e**2*yb+
     &      f**2*yb+2*a*b*YM+b**2*YM-c**2*YM+2*a*e*YM+2*b*e*YM+
     &      e**2*YM-f**2*YM+c*d*zb+a*f*zb+b*f*zb-c*d*ZM-a*f*ZM-
     &      b*f*ZM)/bottom**3+
     &  4*(c*xb-c*XM-b*yb-e*yb+b*YM+e*YM)*
     &    (-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+e**2*ya-
     &      f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za+c*d*ZM+a*f*ZM+
     &      b*f*ZM)/bottom**2+
     &  4*(c*xa-c*XM-b*ya-e*ya+b*YM+e*YM)*
     &    (-(a*c*xb)-c*e*xb-d*f*xb+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*yb+b**2*yb-c**2*yb+2*a*e*yb+2*b*e*yb+e**2*yb-
     &      f**2*yb-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*zb-a*f*zb-b*f*zb+c*d*ZM+a*f*ZM+
     &      b*f*ZM)/bottom**2)/(2.0D0*MASS)
            HESS(J4+2,J3+2)=HESS(J3+2,J4+2)

            HESS(J3+3,J4+3)=HESS(J3+3,J4+3)+ZJ**2*(
     &   (-2*a/N-4*b/N-2*e/N+8*(za-ZM)*(zb-ZM))/bottom-  
     &  top*(-2*a**2/N-4*a*b/N+2*d**2/N-4*a*e/N-
     &      4*b*e/N-2*e**2/N+2*f**2/N-
     &      2*b*(xa-XM)*(xb-XM)-2*e*(xa-XM)*(xb-XM)-
     &      2*c*(xb-XM)*(ya-YM)-2*c*(xa-XM)*(yb-YM)-
     &      2*a*(ya-YM)*(yb-YM)-2*b*(ya-YM)*(yb-YM)-
     &      4*d*(xb-XM)*(za-ZM)-4*f*(yb-YM)*(za-ZM)-
     &      4*d*(xa-XM)*(zb-ZM)-4*f*(ya-YM)*(zb-ZM)+
     &      8*a*(za-ZM)*(zb-ZM)+8*e*(za-ZM)*(zb-ZM))/bottom**2-
     &  4*(a+2*b+e)*(zb-ZM)*
     &    (-(b*d*xa)-d*e*xa-c*f*xa+b*d*XM+d*e*XM+c*f*XM-c*d*ya-
     &      a*f*ya-b*f*ya+c*d*YM+a*f*YM+b*f*YM+a**2*za+2*a*b*za-
     &      d**2*za+2*a*e*za+2*b*e*za+e**2*za-f**2*za-a**2*ZM-
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**2-4*(a+2*b+e)*(za-ZM)*
     &    (-(b*d*xb)-d*e*xb-c*f*xb+b*d*XM+d*e*XM+c*f*XM-c*d*yb-
     &      a*f*yb-b*f*yb+c*d*YM+a*f*YM+b*f*YM+a**2*zb+2*a*b*zb-
     &      d**2*zb+2*a*e*zb+2*b*e*zb+e**2*zb-f**2*zb-a**2*ZM-
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**2+8*top*(-(b*d*xa)-d*e*xa-c*f*xa+b*d*XM+d*e*XM+
     &      c*f*XM-c*d*ya-a*f*ya-b*f*ya+c*d*YM+a*f*YM+b*f*YM+
     &      a**2*za+2*a*b*za-d**2*za+2*a*e*za+2*b*e*za+e**2*za-
     &      f**2*za-a**2*ZM-2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-
     &      e**2*ZM+f**2*ZM)*(-(b*d*xb)-d*e*xb-c*f*xb+b*d*XM+
     &      d*e*XM+c*f*XM-c*d*yb-a*f*yb-b*f*yb+c*d*YM+a*f*YM+
     &      b*f*YM+a**2*zb+2*a*b*zb-d**2*zb+2*a*e*zb+2*b*e*zb+
     &      e**2*zb-f**2*zb-a**2*ZM-2*a*b*ZM+d**2*ZM-2*a*e*ZM-
     &      2*b*e*ZM-e**2*ZM+f**2*ZM)/bottom**3)/(2.0D0*MASS)
            HESS(J4+3,J3+3)=HESS(J3+3,J4+3)
         ENDDO
      ENDDO
C
C  Case III: same atom, different Cartesian coordinate
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         HESS(J3+1,J3+2)=HESS(J3+1,J3+2)+ZJ**2*(
     &  (-2*c+2*c/N+2*xa*ya-2*XM*ya-2*xa*YM+2*XM*YM)/bottom-
     &  top*(-2*a*c*XOM-2*c*e*XOM-
     &      2*d*f*XOM-4*c*(xa-XM)**2+
     &      6*a*(xa-XM)*(ya-YM)+8*b*(xa-XM)*(ya-YM)+
     &      6*e*(xa-XM)*(ya-YM)-4*c*(ya-YM)**2-
     &      2*f*(xa-XM)*(za-ZM)-2*d*(ya-YM)*(za-ZM)-
     &      2*c*(za-ZM)**2)/bottom**2-
     &  4*(a*xa+b*xa-a*XM-b*XM-c*ya+c*YM)*
     &    (-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+e**2*ya-
     &      f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za+c*d*ZM+a*f*ZM+
     &      b*f*ZM)/bottom**2+
     &  4*(c*xa-c*XM-b*ya-e*ya+b*YM+e*YM)*
     &    (a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     &      2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     &      2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     &      c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     &      c*f*ZM)/bottom**2+
     &  8*top*(-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+e**2*ya-
     &      f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za+c*d*ZM+a*f*ZM+
     &      b*f*ZM)*(a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+
     &      2*a*e*xa+2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+
     &      d**2*XM-2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+
     &      a*c*YM+c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+
     &      d*e*ZM+c*f*ZM)/bottom**3)/(2.0D0*MASS)
         HESS(J3+2,J3+1)=HESS(J3+1,J3+2)

         HESS(J3+1,J3+3)=HESS(J3+1,J3+3)+ZJ**2*(
     & -(top*(-2*b*d*XOM-2*d*e*XOM-2*c*f*XOM-
     &        4*d*(xa-XM)**2-2*f*(xa-XM)*(ya-YM)-2*d*(ya-YM)**2+
     &        8*a*(xa-XM)*(za-ZM)+6*b*(xa-XM)*(za-ZM)+
     &        6*e*(xa-XM)*(za-ZM)-2*c*(ya-YM)*(za-ZM)-
     &        4*d*(za-ZM)**2)/bottom**2)+4*(xa-XM)*(za-ZM)/bottom-
     &  4*(a+2*b+e)*(za-ZM)*
     &    (a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     &      2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     &      2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     &      c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     &      c*f*ZM)/bottom**2+
     &  4*(a*xa+b*xa-a*XM-b*XM-c*ya+c*YM)*
     &    (b*d*xa+d*e*xa+c*f*xa-b*d*XM-d*e*XM-c*f*XM+c*d*ya+
     &      a*f*ya+b*f*ya-c*d*YM-a*f*YM-b*f*YM-a**2*za-2*a*b*za+
     &      d**2*za-2*a*e*za-2*b*e*za-e**2*za+f**2*za+a**2*ZM+
     &      2*a*b*ZM-d**2*ZM+2*a*e*ZM+2*b*e*ZM+e**2*ZM-f**2*ZM)/
     &   bottom**2+8*top*(a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-
     &      d**2*xa+2*a*e*xa+2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+
     &      c**2*XM+d**2*XM-2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-
     &      d*f*ya+a*c*YM+c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+
     &      b*d*ZM+d*e*ZM+c*f*ZM)*
     &    (-(b*d*xa)-d*e*xa-c*f*xa+b*d*XM+d*e*XM+c*f*XM-c*d*ya-
     &      a*f*ya-b*f*ya+c*d*YM+a*f*YM+b*f*YM+a**2*za+2*a*b*za-
     &      d**2*za+2*a*e*za+2*b*e*za+e**2*za-f**2*za-a**2*ZM-
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**3 )/(2.0D0*MASS)
         HESS(J3+3,J3+1)=HESS(J3+1,J3+3)

         HESS(J3+2,J3+3)=HESS(J3+2,J3+3)+ZJ**2*(
     &  -(top*(-2*c*d*XOM-2*a*f*XOM-2*b*f*XOM-
     &        2*f*(xa-XM)**2-2*d*(xa-XM)*(ya-YM)-4*f*(ya-YM)**2-
     &        2*c*(xa-XM)*(za-ZM)+6*a*(ya-YM)*(za-ZM)+
     &        6*b*(ya-YM)*(za-ZM)+8*e*(ya-YM)*(za-ZM)-
     &        4*f*(za-ZM)**2)/bottom**2)+4*(ya-YM)*(za-ZM)/bottom+   
     &  4*(a+2*b+e)*(za-ZM)*
     &    (a*c*xa+c*e*xa+d*f*xa-a*c*XM-c*e*XM-d*f*XM-2*a*b*ya-
     &      b**2*ya+c**2*ya-2*a*e*ya-2*b*e*ya-e**2*ya+f**2*ya+
     &      2*a*b*YM+b**2*YM-c**2*YM+2*a*e*YM+2*b*e*YM+e**2*YM-
     &      f**2*YM+c*d*za+a*f*za+b*f*za-c*d*ZM-a*f*ZM-b*f*ZM)/
     &   bottom**2+4*(-(c*xa)+c*XM+b*ya+e*ya-b*YM-e*YM)*
     &    (b*d*xa+d*e*xa+c*f*xa-b*d*XM-d*e*XM-c*f*XM+c*d*ya+
     &      a*f*ya+b*f*ya-c*d*YM-a*f*YM-b*f*YM-a**2*za-2*a*b*za+
     &      d**2*za-2*a*e*za-2*b*e*za-e**2*za+f**2*za+a**2*ZM+
     &      2*a*b*ZM-d**2*ZM+2*a*e*ZM+2*b*e*ZM+e**2*ZM-f**2*ZM)/
     &   bottom**2+8*top*(-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM+
     &      d*f*XM+2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya+
     &      e**2*ya-f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-
     &      2*b*e*YM-e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za+
     &      c*d*ZM+a*f*ZM+b*f*ZM)*
     &    (-(b*d*xa)-d*e*xa-c*f*xa+b*d*XM+d*e*XM+c*f*XM-c*d*ya-
     &      a*f*ya-b*f*ya+c*d*YM+a*f*YM+b*f*YM+a**2*za+2*a*b*za-
     &      d**2*za+2*a*e*za+2*b*e*za+e**2*za-f**2*za-a**2*ZM-
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**3 )/(2.0D0*MASS)
         HESS(J3+3,J3+2)=HESS(J3+2,J3+3)
      ENDDO
C
C  Case IV: Different atoms and different Cartesian coordinates.
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         DO J2=1,N
            IF (J2.NE.J1) THEN
            J4=3*(J2-1)
            XB=X(J4+1)
            YB=X(J4+2)
            ZB=X(J4+3)
            HESS(J3+1,J4+2)=HESS(J3+1,J4+2)+ZJ**2*(
     & (2*c/N-2*(xb-XM)*(ya-YM)+4*(xa-XM)*(yb-YM))/bottom- 
     &  top*(2*a*c/N+2*c*e/N+2*d*f/N-
     &      4*c*(xa-XM)*(xb-XM)-2*a*(xb-XM)*(ya-YM)-
     &      2*e*(xb-XM)*(ya-YM)+8*a*(xa-XM)*(yb-YM)+
     &      8*b*(xa-XM)*(yb-YM)+8*e*(xa-XM)*(yb-YM)-
     &      4*c*(ya-YM)*(yb-YM)-2*f*(xb-XM)*(za-ZM)-
     &      2*d*(ya-YM)*(zb-ZM)-2*c*(za-ZM)*(zb-ZM))/bottom**2-
     &  4*(a*xa+b*xa-a*XM-b*XM-c*ya+c*YM)*
     &    (-(a*c*xb)-c*e*xb-d*f*xb+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*yb+b**2*yb-c**2*yb+2*a*e*yb+2*b*e*yb+e**2*yb-
     &      f**2*yb-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*zb-a*f*zb-b*f*zb+c*d*ZM+a*f*ZM+
     &      b*f*ZM)/bottom**2+
     &  4*(c*xb-c*XM-b*yb-e*yb+b*YM+e*YM)*
     &    (a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa+
     &      2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM-
     &      2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM+
     &      c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM+
     &      c*f*ZM)/bottom**2+
     &  8*top*(-(a*c*xb)-c*e*xb-d*f*xb+a*c*XM+c*e*XM+d*f*XM+
     &      2*a*b*yb+b**2*yb-c**2*yb+2*a*e*yb+2*b*e*yb+e**2*yb-
     &      f**2*yb-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM-2*b*e*YM-
     &      e**2*YM+f**2*YM-c*d*zb-a*f*zb-b*f*zb+c*d*ZM+a*f*ZM+
     &      b*f*ZM)*(a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+
     &      2*a*e*xa+2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+
     &      d**2*XM-2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+
     &      a*c*YM+c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+
     &      d*e*ZM+c*f*ZM)/bottom**3 )/(2.0D0*MASS)
            HESS(J4+2,J3+1)=HESS(J3+1,J4+2)

            HESS(J3+1,J4+3)=HESS(J3+1,J4+3)+ZJ**2*(
     &  -(top*(2*b*d/N+2*d*e/N+2*c*f/N-4*d*(xa-XM)*(xb-XM) -
     &        2*f*(xb-XM)*(ya-YM)-2*d*(ya-YM)*(yb-YM) -
     &        2*b*(xb-XM)*(za-ZM)-2*e*(xb-XM)*(za-ZM) -
     &        2*c*(yb-YM)*(za-ZM)+8*a*(xa-XM)*(zb-ZM) +
     &        8*b*(xa-XM)*(zb-ZM)+8*e*(xa-XM)*(zb-ZM) -
     &        4*d*(za-ZM)*(zb-ZM))/bottom**2) +
     &  4*(xa-XM)*(zb-ZM)/bottom -
     &  4*(a+2*b+e)*(zb-ZM)*
     &    (a**2*xa+2*a*b*xa+b**2*xa-c**2*xa-d**2*xa+2*a*e*xa +
     &      2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM+c**2*XM+d**2*XM -
     &      2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya-d*f*ya+a*c*YM +
     &      c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za+b*d*ZM+d*e*ZM +
     &      c*f*ZM)/bottom**2 +
     &  4*(a*xa+b*xa-a*XM-b*XM-c*ya+c*YM)*
     &    (b*d*xb+d*e*xb+c*f*xb-b*d*XM-d*e*XM-c*f*XM+c*d*yb +
     &      a*f*yb+b*f*yb-c*d*YM-a*f*YM-b*f*YM-a**2*zb-2*a*b*zb +
     &      d**2*zb-2*a*e*zb-2*b*e*zb-e**2*zb+f**2*zb+a**2*ZM +
     &      2*a*b*ZM-d**2*ZM+2*a*e*ZM+2*b*e*ZM+e**2*ZM-f**2*ZM)/
     &   bottom**2+8*top*(a**2*xa+2*a*b*xa+b**2*xa-c**2*xa -
     &      d**2*xa+2*a*e*xa+2*b*e*xa-a**2*XM-2*a*b*XM-b**2*XM +
     &      c**2*XM+d**2*XM-2*a*e*XM-2*b*e*XM-a*c*ya-c*e*ya -
     &      d*f*ya+a*c*YM+c*e*YM+d*f*YM-b*d*za-d*e*za-c*f*za +
     &      b*d*ZM+d*e*ZM+c*f*ZM)*
     &    (-(b*d*xb)-d*e*xb-c*f*xb+b*d*XM+d*e*XM+c*f*XM-c*d*yb -
     &      a*f*yb-b*f*yb+c*d*YM+a*f*YM+b*f*YM+a**2*zb+2*a*b*zb -
     &      d**2*zb+2*a*e*zb+2*b*e*zb+e**2*zb-f**2*zb-a**2*ZM -
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**3 )/(2.0D0*MASS)
            HESS(J4+3,J3+1)=HESS(J3+1,J4+3)

            HESS(J3+2,J4+3)=HESS(J3+2,J4+3)+ZJ**2*(
     &   -(top*(2*c*d/N+2*a*f/N+2*b*f/N-2*f*(xa-XM)*(xb-XM) -
     &        2*d*(xa-XM)*(yb-YM)-4*f*(ya-YM)*(yb-YM) -
     &        2*c*(xb-XM)*(za-ZM)-2*a*(yb-YM)*(za-ZM) -
     &        2*b*(yb-YM)*(za-ZM)+8*a*(ya-YM)*(zb-ZM) +
     &        8*b*(ya-YM)*(zb-ZM)+8*e*(ya-YM)*(zb-ZM) -
     &        4*f*(za-ZM)*(zb-ZM))/bottom**2) +
     &  4*(ya-YM)*(zb-ZM)/bottom +
     &  4*(a+2*b+e)*(zb-ZM)*
     &    (a*c*xa+c*e*xa+d*f*xa-a*c*XM-c*e*XM-d*f*XM-2*a*b*ya -
     &      b**2*ya+c**2*ya-2*a*e*ya-2*b*e*ya-e**2*ya+f**2*ya +
     &      2*a*b*YM+b**2*YM-c**2*YM+2*a*e*YM+2*b*e*YM+e**2*YM -
     &      f**2*YM+c*d*za+a*f*za+b*f*za-c*d*ZM-a*f*ZM-b*f*ZM)/
     &   bottom**2+4*(-(c*xa)+c*XM+b*ya+e*ya-b*YM-e*YM)*
     &    (b*d*xb+d*e*xb+c*f*xb-b*d*XM-d*e*XM-c*f*XM+c*d*yb +
     &      a*f*yb+b*f*yb-c*d*YM-a*f*YM-b*f*YM-a**2*zb-2*a*b*zb +
     &      d**2*zb-2*a*e*zb-2*b*e*zb-e**2*zb+f**2*zb+a**2*ZM +
     &      2*a*b*ZM-d**2*ZM+2*a*e*ZM+2*b*e*ZM+e**2*ZM-f**2*ZM)/
     &   bottom**2+8*top*(-(a*c*xa)-c*e*xa-d*f*xa+a*c*XM+c*e*XM +
     &      d*f*XM+2*a*b*ya+b**2*ya-c**2*ya+2*a*e*ya+2*b*e*ya +
     &      e**2*ya-f**2*ya-2*a*b*YM-b**2*YM+c**2*YM-2*a*e*YM -
     &      2*b*e*YM-e**2*YM+f**2*YM-c*d*za-a*f*za-b*f*za +
     &      c*d*ZM+a*f*ZM+b*f*ZM)*
     &    (-(b*d*xb)-d*e*xb-c*f*xb+b*d*XM+d*e*XM+c*f*XM-c*d*yb -
     &      a*f*yb-b*f*yb+c*d*YM+a*f*YM+b*f*YM+a**2*zb+2*a*b*zb -
     &      d**2*zb+2*a*e*zb+2*b*e*zb+e**2*zb-f**2*zb-a**2*ZM -
     &      2*a*b*ZM+d**2*ZM-2*a*e*ZM-2*b*e*ZM-e**2*ZM+f**2*ZM)/
     &   bottom**3 )/(2.0D0*MASS)
            HESS(J4+3,J3+2)=HESS(J3+2,J4+3)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
C
C  Energy for linear system
C
      SUBROUTINE EPANIC(X,N,ZJ,ROTE,FLAG,V,MASS)
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, N, J2
      LOGICAL FLAG
      DOUBLE PRECISION X(3*N), ROTE, ZJ, 
     1                 V(3*N), TEMP2, TEMP1,DIF,
     4                 V1, V2, V3, V4, MASS
      CALL NROTE(X,N,ROTE,MASS,ZJ)
      IF (FLAG) RETURN
      DIF=1.0D-4
      DO 20 J1=1,3*N
         TEMP1=X(J1)
         X(J1)=X(J1)+DIF
C        PRINT*,'J1,X=',J1,X(J1)
         CALL NROTE(X,N,V1,MASS,ZJ)
         X(J1)=X(J1)-2.0D0*DIF
C        PRINT*,'J1,X=',J1,X(J1)
         CALL NROTE(X,N,V2,MASS,ZJ)
         V(J1)=(V1-V2)/(2.0D0*DIF)
         X(J1)=TEMP1
C        PRINT*,'J1,GRAD=',J1,V(J1)
         DO 10 J2=J1,3*N
            TEMP1=X(J1)
            TEMP2=X(J2)
            X(J1)=X(J1)+DIF
            X(J2)=X(J2)+DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V1,MASS,ZJ)
            X(J1)=X(J1)-2.0D0*DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V2,MASS,ZJ)
            X(J2)=X(J2)-2.0D0*DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V3,MASS,ZJ)
            X(J1)=X(J1)+2.0D0*DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V4,MASS,ZJ)
            X(J1)=TEMP1
            X(J2)=TEMP2
            HESS(J1,J2)=(V1-V2-V4+V3)/(4.0D0*DIF*DIF)
            HESS(J2,J1)=HESS(J1,J2)
C           PRINT*,'V1,V2,V3,V4=',V1,V2,V3,V4
C           PRINT*,'J1,J2,HESS=',J1,J2,HESS(J1,J2)
10       CONTINUE
20    CONTINUE

      RETURN
      END

      SUBROUTINE NROTE(X,N,ROTE,MASS,ZJ)
      IMPLICIT NONE
      INTEGER J1, N, J3
      DOUBLE PRECISION MASS, X(3*N), ROTE, ZJ,
     1                 XM, YM, ZM, TOP, BOTTOM, 
     2                 A, B, C, D, E, F, XT(N)
C
C  Calculate centre of mass components XM, YM, ZM
C
      XM=0.0D0
      YM=0.0D0
      ZM=0.0D0
      DO J1=1,N
         XM=XM+X(3*(J1-1)+1)
         YM=YM+X(3*(J1-1)+2)
         ZM=ZM+X(3*(J1-1)+3)
      ENDDO
      XM=XM/N
      YM=YM/N
      ZM=ZM/N
C
C  Calculate constants A, B, C, D, E, F
C
      A=0.0D0
      B=0.0D0
      C=0.0D0
      D=0.0D0
      E=0.0D0
      F=0.0D0
      DO J1=1,N
         J3=3*(J1-1)
         A=A+(X(J3+2)-YM)**2
         B=B+(X(J3+3)-ZM)**2
         E=E+(X(J3+1)-XM)**2
         C=C+(X(J3+2)-YM)*(X(J3+1)-XM)
         D=D+(X(J3+3)-ZM)*(X(J3+1)-XM)
         F=F+(X(J3+3)-ZM)*(X(J3+2)-YM)
      ENDDO
      TOP=A*B+B*B-C*C+A*E+B*E
      BOTTOM=A*A*B+A*B*B-A*C*C-B*D*D+A*A*E+2.0D0*A*B*E+B*B*E-C*C*E-D*D*E
     1      +A*E*E+B*E*E-2.0D0*C*D*F-A*F*F-B*F*F
C     IF (DABS(BOTTOM).LT.1.0D-10) THEN
         DO J1=1,N
            IF (X(3*(J1-1)+1).NE.0.0D0) THEN
               XT(J1)=DSQRT(X(3*(J1-1)+1)**2
     1         +X(3*(J1-1)+2)**2)*X(3*(J1-1)+1)/DABS(X(3*(J1-1)+1))
            ELSE
               XT(J1)=0.0D0
            ENDIF
C           PRINT*,'J1,XT=',J1,XT(J1)
         ENDDO
         XM=0.0D0
         DO J1=1,N
            XM=XM+XT(J1)
         ENDDO
         XM=XM/N
         E=0.0D0
         DO J1=1,N
            E=E+(XT(J1)-XM)**2
         ENDDO
         ROTE=ZJ**2/(2.0D0*E*MASS)
C        PRINT*,'PANIC ROTE, E, BOTTOM=',ROTE,E, BOTTOM
C     ELSE
C        IZINV=TOP/(MASS*BOTTOM)
C        ROTE=ZJ**2*IZINV/2.0D0
C        PRINT*,'normal ROTE, BOTTOM=',ROTE, BOTTOM
C     ENDIF
      RETURN
      END
C
C**************************************************************************
C
C  Subroutine ROTDERIV calculates analytically the cartesian gradient and
C  second derivative matrix due to the ROTATIONAL TERM ONLY and increments
C  the existing matrices containing the potential contributions
C  accordingly.
C
C**************************************************************************
C
      SUBROUTINE ROTDERIV(N, X, V, MASS, AVEL, IZ)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), MASS, AVEL, IZ
C
C  Analytic gradient:
C
      DO J1 = 1, N
         DO J2 = 1, 2
            J3 = 3*(J1-1)+J2
            V(J3) = V(J3) + (AVEL**2) * X(J3) * MASS
         ENDDO
      ENDDO
C
C  Case I: entirely diagonal terms.
C
      DO J1 = 1, N
         DO J2 = 1, 2
            J3 = 3*(J1-1)+J2
            HESS(J3, J3) = HESS(J3, J3) + (AVEL**2) * MASS
         ENDDO
      ENDDO
      RETURN
      END
C
C*************************************************************************
C
C  THIS SUBROUTINE CALCULATES THE MOMENT OF INERTIA IZ ABOUT THE Z AXIS
C  AND HENCE THE ROTATIONAL ENERGY
C
C  ENERGY = omega^2 * I_{zz} / 2 
C
C*************************************************************************
C
      SUBROUTINE ROTENERGY(N, X, AVEL, MASS, IZ, ROT)
      IMPLICIT NONE
      INTEGER N, I
      DOUBLE PRECISION X(3*N), AVEL, ROT, IZ, MASS
      IZ = 0.0D0
      DO 10 I = 0, N-1
         IZ = IZ + MASS * ( X(I*3+1)**2 + X(I*3+2)**2 )
10    CONTINUE
      ROT = 0.5D0 * IZ * AVEL**2
      RETURN
      END
