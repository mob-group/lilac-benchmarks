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
C********************************************************************
C
C Periodic bulk C60
C
C********************************************************************
C
      SUBROUTINE C60P(N,X,V,ENERGY,BOXLX,BOXLY,BOXLZ,CUTOFF)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, MX, MY, MZ, J3, J4, J5, J6
      DOUBLE PRECISION X(3*N), R2(N,N),
     1                 V(3*N), F(N,N), ENERGY, TEMP,
     2                 R(N,N), G(N,N), 
     3                 T2, T4, T3, T9, D, DSQ, D2, VEC(N,N,3),
     4                 CM3, CP3, CM4, CP4, CM9, CP9, CM10, CP10,  
     5                 D2M, D2P, CM5, CP5, CM11, CP11, T6, T12, DIST,
     6                 BOXLX, BOXLY, BOXLZ, CUTOFF
      PARAMETER(D=1.023349668D0, DSQ=D*D, D2=2.0D0*D)
C
C  Deal with atoms leaving the box:
C
      DO 41 J1=1,N
         X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX*DNINT(X(3*(J1-1)+1)/BOXLX)  
         X(3*(J1-1)+2)=X(3*(J1-1)+2) - BOXLY*DNINT(X(3*(J1-1)+2)/BOXLY)
         X(3*(J1-1)+3)=X(3*(J1-1)+3) - BOXLZ*DNINT(X(3*(J1-1)+3)/BOXLZ)
41    CONTINUE
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is
C  used:
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            MZ=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(MX)
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(MY)
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(MZ)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         R2(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
            R2(J2,J1)=1.0D0/DIST
            R(J2,J1)=DSQRT(DIST)
            R(J1,J2)=R(J2,J1)
            R2(J1,J2)=R2(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  Calculate the energy and the g and f tensors.
C
      ENERGY=0.0D0
      DO 21 J1=1,N
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO 22 J2=J1+1,N 
            TEMP=R(J2,J1)
            IF (TEMP.LT.CUTOFF) THEN
               T2=1.0D0/(TEMP*TEMP)
               T3=T2/TEMP
               T4=T2*T2
               T6=T4*T2
               T9=T6*T3
               T12=T9*T3
               D2M=1.0D0/(-D2+TEMP)
               D2P=1.0D0/(D2+TEMP)
               CM3=D2M**3
               CP3=D2P**3
               CM4=CM3*D2M
               CP4=CP3*D2P
               CM5=CM4*D2M
               CP5=CP4*D2P
               CM9=CM5*CM4
               CP9=CP5*CP4
               CM10=CM9*D2M
               CP10=CP9*D2P
               CM11=CM10*D2M
               CP11=CP10*D2P
               G(J2,J1)=(4.0D0*((-10.0D0*(-2.0D0*T9+CM9+CP9)
     1          +75.0D0*(-2.0D0*T3+CM3+CP3))*T2
     2         +(10.0D0*(18.0D0*T9/TEMP-9.0D0*(CM10+CP10))
     3          -75.0D0*(6.0D0*T4-3.0D0*(CM4+CP4)))/TEMP))/(DSQ*TEMP) 
               G(J1,J2)=G(J2,J1)

               F(J2,J1)=(-9600*T12+14400*T6
     1                 +3600*(CM11+CP11-CM5-CP5)/TEMP
     3                 +(1080*(CM10+CP10)-2700*(CM4+CP4))*T2
     5                 +(120*(CM9+CP9)-900*(CM3+CP3))*T3)/DSQ
               F(J1,J2)=F(J2,J1)

               ENERGY=ENERGY+(10.0D0*(-2.0D0*T9+CM9+CP9)
     1                     -75.0D0*(-2.0D0*T3+CM3+CP3))/TEMP  
            ENDIF
22       CONTINUE
21    CONTINUE
      ENERGY=ENERGY*4.0D0/DSQ
      ENERGY=ENERGY/96.794820624738D0
C
C  From here on down the code is system-independent!
C
C  First the gradient.
C
      DO 13 J1=1,N
         DO 14 J2=1,3
            J3=3*(J1-1)+J2
            TEMP=0.0D0
            DO 17 J4=1,N
               IF (R(J4,J1).LT.CUTOFF) TEMP=TEMP+G(J4,J1)*VEC(J1,J4,J2)  
17          CONTINUE
            V(J3)=TEMP/96.794820624738D0
14       CONTINUE
13    CONTINUE
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            TEMP=0.0D0
            DO 60 J4=1,N
               IF (R(J4,J1).LT.CUTOFF) THEN
                  TEMP=TEMP+F(J4,J1)*R2(J4,J1)*VEC(J1,J4,J2)**2+G(J4,J1)     
               ENDIF
60          CONTINUE
            HESS(J3,J3)=TEMP/96.794820624738D0
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
               TEMP=0.0D0
               DO 90 J5=1,N
                  IF (R(J5,J1).LT.CUTOFF) THEN
                     TEMP=TEMP+
     1               F(J5,J1)*R2(J5,J1)*VEC(J1,J5,J2)*VEC(J1,J5,J4) 
                  ENDIF
90             CONTINUE
              HESS(3*(J1-1)+J4,J3)=TEMP/96.794820624738D0
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
               IF (R(J4,J1).LT.CUTOFF) THEN
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                              VEC(J1,J4,J2)**2-G(J4,J1)
                  HESS(3*(J4-1)+J2,J3)=HESS(3*(J4-1)+J2,J3)/96.794820624738D0
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
               ENDIF
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
                  IF (R(J4,J1).LT.CUTOFF) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                       *VEC(J1,J4,J2)*VEC(J1,J4,J5)
                     HESS(J6,J3)=HESS(J6,J3)/96.794820624738D0
                     HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
                  ELSE
                     HESS(J6,J3)=0.0D0
                     HESS(3*(J4-1)+J2,3*(J1-1)+J5)=0.0D0
                  ENDIF
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
