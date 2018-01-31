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
      SUBROUTINE FD(NATOMS,P,GRAD,ENERGY)
      USE MODHESS
      USE KEY
      IMPLICIT NONE
      INTEGER J1, J2, NATOMS
      DOUBLE PRECISION P(3*NATOMS), ND5H, NTD, NOH, NIH, PI, ENERGY, GRAD(3*NATOMS), DUMMY1,
     1                 DUMMY2, X, Y, Z, R, DUMMY3, XMASS, YMASS, ZMASS
     2                 
      PARAMETER (PI=3.141592654D0)
C     PARAMETER (ND5H=16.0D0*DSQRT(1.0D0/7.0D0)/(3.0D0*PI), 
C    1           NTD=8.0D0/PI,
C    2           NOH=64.0D0*DSQRT(2.0D0/4449.0D0)/PI,
C    3           NIH=256.0D0*DSQRT(2.0D0/105614964576169.0D0)/PI)
      PARAMETER (ND5H=0.641652418D0,
     1           NTD=2.546479089D0,
     2           NOH=0.431930524D0,
     3           NIH=3.546029643D-3)

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
C     DO 10 I=1,NATOMS
C        XMASS=XMASS+P(3*I-2)
C        YMASS=YMASS+P(3*I-1)
C        ZMASS=ZMASS+P(3*I)
C10    CONTINUE
C     XMASS=XMASS/NATOMS
C     YMASS=YMASS/NATOMS
C     ZMASS=ZMASS/NATOMS

      IF (D5HT) THEN
         DUMMY1=0.0D0
         DUMMY2=FD5H*ND5H
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=DSQRT(X*X+Y*Y+Z*Z)
            DUMMY3=(X**5 - 10.0D0*X**3*Y**2 + 5.0D0*X*Y**4)/R**5
            DUMMY1=DUMMY1 + DUMMY3
            GRAD(J2-2)=GRAD(J2-2) + (5.0D0*(R**2*X**4 - X**6 - 6.0D0*R**2*X**2*Y**2 + 10.0D0*X**4*Y**2 + R**2*Y**4 -
     1                               5.0D0*X**2*Y**4)/R**7)*DUMMY2
            GRAD(J2-1)=GRAD(J2-1) + (5.0D0*X*Y*(-4.0D0*R**2*X**2 - X**4 + 4.0D0*R**2*Y**2 + 10.0D0*X**2*Y**2 - 
     1                               5.0D0*Y**4)/R**7)*DUMMY2
            GRAD(J2)=GRAD(J2)     + (-5.0D0*(X**5 - 10.0D0*X**3*Y**2 + 5.0D0*X*Y**4)*Z/R**7)*DUMMY2
            HESS(J2-2,J2-2)=HESS(J2-2,J2-2)+(5.0D0*X*(4.0D0*R**4*X**2 - 11.0D0*R**2*X**4 + 7.0D0*X**6 
     1                                      - 12.0D0*R**4*Y**2 + 70*R**2*X**2*Y**2 -  
     1                                       70.0D0*X**4*Y**2 - 15.0D0*R**2*Y**4 + 35.0D0*X**2*Y**4)/R**9)*DUMMY2
            HESS(J2-2,J2-1)=HESS(J2-2,J2-1)+(5.0D0*Y*(-12.0D0*R**4*X**2 + 15.0D0*R**2*X**4 + 7.0D0*X**6 + 4.0D0*R**4*Y**2 +
     1                                       10.0D0*R**2*X**2*Y**2 - 70.0D0*X**4*Y**2 - 5.0D0*R**2*Y**4 
     1                                       + 35.0D0*X**2*Y**4)/R**9)*DUMMY2
            HESS(J2-1,J2-2)=HESS(J2-2,J2-1)
            HESS(J2-2,J2)=HESS(J2-2,J2)+(5.0D0*(-5.0D0*R**2*X**4 + 7.0D0*X**6 + 30.0D0*R**2*X**2*Y**2 
     1                                   - 70.0D0*X**4*Y**2 - 5.0D0*R**2*Y**4 +
     1                                   35.0D0*X**2*Y**4)*Z/R**9)*DUMMY2
            HESS(J2,J2-2)=HESS(J2-2,J2)
            HESS(J2-1,J2-1)=HESS(J2-1,J2-1)+(5.0D0*X*(-4.0D0*R**4*X**2 - R**2*X**4 + 12.0D0*R**4*Y**2 
     1                                       + 50.0D0*R**2*X**2*Y**2 +
     1                                       7.0D0*X**4*Y**2 - 45.0D0*R**2*Y**4 - 70.0D0*X**2*Y**4 
     1                                       + 35.0D0*Y**6)/R**9)*DUMMY2
            HESS(J2-1,J2)=HESS(J2-1,J2)+(5.0D0*X*Y*(20.0D0*R**2*X**2 + 7.0D0*X**4 - 20.0D0*R**2*Y**2 
     1                                   - 70.0D0*X**2*Y**2 + 35.0D0*Y**4)*Z/R**9)*DUMMY2  
            HESS(J2,J2-1)=HESS(J2-1,J2)
            HESS(J2,J2)=HESS(J2,J2)+(5.0D0*X*(X**4 - 10.0D0*X**2*Y**2 + 5.0D0*Y**4)*(-R**2 + 7.0D0*Z**2)/R**9)*DUMMY2
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF

      IF (TDT) THEN
         DUMMY1=0.0D0
         DUMMY2=FTD*NTD
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=DSQRT(X*X+Y*Y+Z*Z)
            DUMMY3=X*Y*Z/R**3
            DUMMY1=DUMMY1 + DUMMY3
            GRAD(J2-2)=GRAD(J2-2) + (Y*Z/R**3-3.0D0*DUMMY3*X/R**2)*DUMMY2
            GRAD(J2-1)=GRAD(J2-1) + (X*Z/R**3-3.0D0*DUMMY3*Y/R**2)*DUMMY2
            GRAD(J2)=GRAD(J2)     + (X*Y/R**3-3.0D0*DUMMY3*Z/R**2)*DUMMY2
            HESS(J2-2,J2-2)=HESS(J2-2,J2-2)+(3*X*(-3*R**2 + 5*X**2)*Y*Z/R**7)*DUMMY2
            HESS(J2-2,J2-1)=HESS(J2-2,J2-1)+((R**4 - 3*R**2*X**2 - 3*R**2*Y**2 + 15*X**2*Y**2)*Z/R**7)*DUMMY2
            HESS(J2-2,J2)=HESS(J2-2,J2)+(Y*(R**4 - 3*R**2*X**2 - 3*R**2*Z**2 + 15*X**2*Z**2)/R**7)*DUMMY2
            HESS(J2-1,J2-1)=HESS(J2-1,J2-1)+(3*X*Y*(-3*R**2 + 5*Y**2)*Z/R**7)*DUMMY2
            HESS(J2-1,J2)=HESS(J2-1,J2)+(X*(R**4 - 3*R**2*Y**2 - 3*R**2*Z**2 + 15*Y**2*Z**2)/R**7)*DUMMY2
            HESS(J2,J2)=HESS(J2,J2)+(3*X*Y*Z*(-3*R**2 + 5*Z**2)/R**7)*DUMMY2
            HESS(J2-1,J2-2)=HESS(J2-2,J2-1)
            HESS(J2,J2-2)=HESS(J2-2,J2)
            HESS(J2,J2-1)=HESS(J2-1,J2)
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF
      IF (OHT) THEN
         DUMMY1=0.0D0
         DUMMY2=FOH*NOH
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=DSQRT(X*X+Y*Y+Z*Z)
            DUMMY3=(X**4 + Y**4 + Z**4 - 3.0D0*(X**2*Y**2 + X**2*Z**2 + Y**2*Z**2))/R**4
            DUMMY1=DUMMY1 + DUMMY3
            GRAD(J2-2)=GRAD(J2-2) + ((4.0D0*X**3 - 6.0D0*(X*Y**2 + X*Z**2))/R**4
     1                                -4.0D0*DUMMY3*X/R**2)*DUMMY2
            GRAD(J2-1)=GRAD(J2-1) + ((4.0D0*Y**3 - 6.0D0*(X**2*Y + Y*Z**2))/R**4
     1                                -4.0D0*DUMMY3*Y/R**2)*DUMMY2
            GRAD(J2)=GRAD(J2)     + ((4.0D0*Z**3 - 6.0D0*(X**2*Z + Y**2*Z))/R**4
     1                                -4.0D0*DUMMY3*Z/R**2)*DUMMY2
            HESS(J2-2,J2-2)=HESS(J2-2,J2-2)+(2*(6*R**4*X**2 - 18*R**2*X**4 + 12*X**6 - 3*R**4*Y**2 + 
     1      30*R**2*X**2*Y**2 -
     1      36*X**4*Y**2 - 2*R**2*Y**4 + 12*X**2*Y**4 - 3*R**4*Z**2 +
     1      30*R**2*X**2*Z**2 - 36*X**4*Z**2 + 6*R**2*Y**2*Z**2 - 
     1      36*X**2*Y**2*Z**2 - 2*R**2*Z**4 + 12*X**2*Z**4)/R**8)*DUMMY2
            HESS(J2-2,J2-1)=HESS(J2-2,J2-1)+(4*X*Y*(-3*R**4 + 2*R**2*X**2 + 6*X**4 + 2*R**2*Y**2 - 18*X**2*Y**2 + 
     1      6*Y**4 + 
     1      12*R**2*Z**2 - 18*X**2*Z**2 - 18*Y**2*Z**2 + 6*Z**4)/R**8)*DUMMY2
            HESS(J2-2,J2)=HESS(J2-2,J2)+(4*X*Z*(-3*R**4 + 2*R**2*X**2 + 6*X**4 + 12*R**2*Y**2 - 18*X**2*Y**2 + 
     1      6*Y**4 + 2*R**2*Z**2 - 18*X**2*Z**2 - 18*Y**2*Z**2 + 6*Z**4)/R**8)*DUMMY2
            HESS(J2-1,J2-1)=HESS(J2-1,J2-1)+(2*(-3*R**4*X**2 - 2*R**2*X**4 + 6*R**4*Y**2 + 30*R**2*X**2*Y**2 + 
     1      12*X**4*Y**2 - 18*R**2*Y**4 - 36*X**2*Y**4 + 12*Y**6 -
     1      3*R**4*Z**2 + 6*R**2*X**2*Z**2 + 30*R**2*Y**2*Z**2 -
     1      36*X**2*Y**2*Z**2 - 36*Y**4*Z**2 - 2*R**2*Z**4 + 12*Y**2*Z**4)/R**8)*DUMMY2
            HESS(J2-1,J2)=HESS(J2-1,J2)+(4*Y*Z*(-3*R**4 + 12*R**2*X**2 + 6*X**4 + 2*R**2*Y**2 - 18*X**2*Y**2 +
     1      6*Y**4 + 2*R**2*Z**2 - 18*X**2*Z**2 - 18*Y**2*Z**2 + 6*Z**4)/R**8)*DUMMY2
            HESS(J2,J2)=HESS(J2,J2)+(2*(-3*R**4*X**2 - 2*R**2*X**4 - 3*R**4*Y**2 + 6*R**2*X**2*Y**2 -
     1      2*R**2*Y**4 + 6*R**4*Z**2 + 30*R**2*X**2*Z**2 + 12*X**4*Z**2 +
     1      30*R**2*Y**2*Z**2 - 36*X**2*Y**2*Z**2 + 12*Y**4*Z**2 -
     1      18*R**2*Z**4 - 36*X**2*Z**4 - 36*Y**2*Z**4 + 12*Z**6)/R**8)*DUMMY2
            HESS(J2-1,J2-2)=HESS(J2-2,J2-1)
            HESS(J2,J2-2)=HESS(J2-2,J2)
            HESS(J2,J2-1)=HESS(J2-1,J2)
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF
      IF (IHT) THEN
         DUMMY1=0.0D0
         DUMMY2=FIH*NIH
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=DSQRT(X*X+Y*Y+Z*Z)
            DUMMY3=(742500.0D0*X**2*Y**2*Z**2 - 25574.0D0*(X**4*Y**2 + Y**4*Z**2 + X**2*Z**4) +
     1                     131824.0D0*(X**2*Y**4 + X**4*Z**2 + Y**2*Z**4) + 8250.0D0*(X**6 + Y**6 + Z**6))/R**6
            DUMMY1=DUMMY1 + DUMMY3
            GRAD(J2-2)=GRAD(J2-2) + ((49500.0D0*X**5 + 1485000.0D0*X*Y**2*Z**2 + 131824.0D0*(2*X*Y**4 + 4*X**3*Z**2) -
     1                            25574.0D0*(4*X**3*Y**2 + 2*X*Z**4))/R**6
     2                              -6.0D0*DUMMY3*X/R**2)*DUMMY2
            GRAD(J2-1)=GRAD(J2-1) + ((49500.D0*Y**5 + 1485000.D0*X**2*Y*Z**2 - 25574.D0*(2*X**4*Y + 4*Y**3*Z**2) +
     1                            131824.D0*(4*X**2*Y**3 + 2*Y*Z**4))/R**6
     2                              -6.0D0*DUMMY3*Y/R**2)*DUMMY2
            GRAD(J2)=GRAD(J2)     + ((1485000.0D0*X**2*Y**2*Z + 49500.0D0*Z**5 - 25574.0D0*(2*Y**4*Z + 4*X**2*Z**3) +
     1                            131824.0D0*(2*X**4*Z + 4*Y**2*Z**3))/R**6
     2                              -6.0D0*DUMMY3*Z/R**2)*DUMMY2
            HESS(J2-2,J2-2)=HESS(J2-2,J2-2)+(4*(61875*R**4*X**4 - 160875*R**2*X**6 + 99000*X**8 - 76722*R**4*X**2*Y**2 +
     1     345249*R**2*X**4*Y**2 - 306888*X**6*Y**2 + 65912*R**4*Y**4 -
     1     988680*R**2*X**2*Y**4 + 1581888*X**4*Y**4 - 12375*R**2*Y**6 +
     1     99000*X**2*Y**6 + 395472*R**4*X**2*Z**2 - 1779624*R**2*X**4*Z**2 +
     1     1581888*X**6*Z**2 + 371250*R**4*Y**2*Z**2 -
     1     5568750*R**2*X**2*Y**2*Z**2 + 8910000*X**4*Y**2*Z**2 +
     1     38361*R**2*Y**4*Z**2 - 306888*X**2*Y**4*Z**2 - 12787*R**4*Z**4 +
     1     191805*R**2*X**2*Z**4 - 306888*X**4*Z**4 - 197736*R**2*Y**2*Z**4 +
     1     1581888*X**2*Y**2*Z**4 - 12375*R**2*Z**6 + 99000*X**2*Z**6)/R**10)*DUMMY2
            HESS(J2-2,J2-1)=HESS(J2-2,J2-1)+(8*X*Y*(-25574*R**4*X**2 + 1236*R**2*X**4 + 49500*X**6 + 131824*R**4*Y**2 -
     1     318750*R**2*X**2*Y**2 - 153444*X**4*Y**2 - 234861*R**2*Y**4 +
     1     790944*X**2*Y**4 + 49500*Y**6 + 371250*R**4*Z**2 -
     1     1509222*R**2*X**2*Z**2 + 790944*X**4*Z**2 -
     1     1037028*R**2*Y**2*Z**2 + 4455000*X**2*Y**2*Z**2 -
     1     153444*Y**4*Z**2 - 159375*R**2*Z**4 - 153444*X**2*Z**4 +
     1     790944*Y**2*Z**4 + 49500*Z**6)/R**10)*DUMMY2
            HESS(J2-2,J2)=HESS(J2-2,J2)+(8*X*Z*(131824*R**4*X**2 - 234861*R**2*X**4 + 49500*X**6 + 371250*R**4*Y**2 -
     1     1037028*R**2*X**2*Y**2 - 153444*X**4*Y**2 - 159375*R**2*Y**4 +
     1     790944*X**2*Y**4 + 49500*Y**6 - 25574*R**4*Z**2 -
     1     318750*R**2*X**2*Z**2 + 790944*X**4*Z**2 -
     1     1509222*R**2*Y**2*Z**2 + 4455000*X**2*Y**2*Z**2 -
     1     153444*Y**4*Z**2 + 1236*R**2*Z**4 - 153444*X**2*Z**4 +
     1     790944*Y**2*Z**4 + 49500*Z**6)/R**10)*DUMMY2
            HESS(J2-1,J2-1)=HESS(J2-1,J2-1)+(4*(-12787*R**4*X**4 - 12375*R**2*X**6 + 395472*R**4*X**2*Y**2 +
     1     191805*R**2*X**4*Y**2 + 99000*X**6*Y**2 + 61875*R**4*Y**4 -
     1     1779624*R**2*X**2*Y**4 - 306888*X**4*Y**4 - 160875*R**2*Y**6 +
     1     1581888*X**2*Y**6 + 99000*Y**8 + 371250*R**4*X**2*Z**2 -
     1     197736*R**2*X**4*Z**2 - 76722*R**4*Y**2*Z**2 -
     1     5568750*R**2*X**2*Y**2*Z**2 + 1581888*X**4*Y**2*Z**2 +
     1     345249*R**2*Y**4*Z**2 + 8910000*X**2*Y**4*Z**2 -
     1     306888*Y**6*Z**2 + 65912*R**4*Z**4 + 38361*R**2*X**2*Z**4 -
     1     988680*R**2*Y**2*Z**4 - 306888*X**2*Y**2*Z**4 +
     1     1581888*Y**4*Z**4 - 12375*R**2*Z**6 + 99000*Y**2*Z**6)/R**10)*DUMMY2
            HESS(J2-1,J2)=HESS(J2-1,J2)+(8*Y*Z*(371250*R**4*X**2 - 159375*R**2*X**4 + 49500*X**6 - 25574*R**4*Y**2 -
     1     1509222*R**2*X**2*Y**2 - 153444*X**4*Y**2 + 1236*R**2*Y**4 +
     1     790944*X**2*Y**4 + 49500*Y**6 + 131824*R**4*Z**2 -
     1     1037028*R**2*X**2*Z**2 + 790944*X**4*Z**2 -
     1     318750*R**2*Y**2*Z**2 + 4455000*X**2*Y**2*Z**2 -
     1     153444*Y**4*Z**2 - 234861*R**2*Z**4 - 153444*X**2*Z**4 +
     1     790944*Y**2*Z**4 + 49500*Z**6)/R**10)*DUMMY2
            HESS(J2,J2)=HESS(J2,J2)+(4*(65912*R**4*X**4 - 12375*R**2*X**6 + 371250*R**4*X**2*Y**2 +
     1     38361*R**2*X**4*Y**2 - 12787*R**4*Y**4 - 197736*R**2*X**2*Y**4 -
     1     12375*R**2*Y**6 - 76722*R**4*X**2*Z**2 - 988680*R**2*X**4*Z**2 +
     1     99000*X**6*Z**2 + 395472*R**4*Y**2*Z**2 -
     1     5568750*R**2*X**2*Y**2*Z**2 - 306888*X**4*Y**2*Z**2 +
     1     191805*R**2*Y**4*Z**2 + 1581888*X**2*Y**4*Z**2 + 99000*Y**6*Z**2 +
     1     61875*R**4*Z**4 + 345249*R**2*X**2*Z**4 + 1581888*X**4*Z**4 -
     1     1779624*R**2*Y**2*Z**4 + 8910000*X**2*Y**2*Z**4 -
     1     306888*Y**4*Z**4 - 160875*R**2*Z**6 - 306888*X**2*Z**6 +
     1     1581888*Y**2*Z**6 + 99000*Z**8)/R**10)*DUMMY2
            HESS(J2-1,J2-2)=HESS(J2-2,J2-1)
            HESS(J2,J2-2)=HESS(J2-2,J2)
            HESS(J2,J2-1)=HESS(J2-1,J2)

         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF

      RETURN
      END

