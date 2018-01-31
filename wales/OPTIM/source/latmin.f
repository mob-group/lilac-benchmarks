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
      SUBROUTINE LATMIN(N,XC,BOXLX,CUTOFF)
      IMPLICIT NONE
      INTEGER N, NCOUNT, J1
      DOUBLE PRECISION XC(3*N),CUTOFF,F1,F2,F3,GRAD,BOXLX,DIF,TEMP1,
     1                 SECOND, XSAVE(3*N), BSAVE, CSAVE

C
C  Value of DIF is the order of magnitude to which the lattice
C  constant can be optimised. Setting it smaller than 10^(-7)
C  causes numerical problems on the DEC.
C
      DIF=1.0D-7
      BSAVE=BOXLX
      CSAVE=CUTOFF
      NCOUNT=1
      DO 20 J1=1,3*N
         XSAVE(J1)=XC(J1)
20    CONTINUE

10    TEMP1=BOXLX
      BOXLX=TEMP1+DIF
      CUTOFF=CSAVE*BOXLX/BSAVE
      CALL C60PE(N,XC,F1,BOXLX,BOXLX,BOXLX,CUTOFF)

      BOXLX=TEMP1-DIF
      CUTOFF=CSAVE*BOXLX/BSAVE
      DO 30 J1=1,3*N
         XC(J1)=XSAVE(J1)
30    CONTINUE
      CALL C60PE(N,XC,F2,BOXLX,BOXLX,BOXLX,CUTOFF)

      GRAD=(F1-F2)/(2.0D0*DIF)

      BOXLX=TEMP1
      CUTOFF=CSAVE*BOXLX/BSAVE
      DO 40 J1=1,3*N
         XC(J1)=XSAVE(J1)
40    CONTINUE
      CALL C60PE(N,XC,F1,BOXLX,BOXLX,BOXLX,CUTOFF)

      BOXLX=TEMP1+2.0D0*DIF
      CUTOFF=CSAVE*BOXLX/BSAVE
      DO 50 J1=1,3*N
         XC(J1)=XSAVE(J1)
50    CONTINUE
      CALL C60PE(N,XC,F2,BOXLX,BOXLX,BOXLX,CUTOFF)

      BOXLX=TEMP1-2.0D0*DIF
      CUTOFF=CSAVE*BOXLX/BSAVE
      DO 60 J1=1,3*N
         XC(J1)=XSAVE(J1)
60    CONTINUE
      CALL C60PE(N,XC,F3,BOXLX,BOXLX,BOXLX,CUTOFF)

      SECOND=(F3+F2-2.0D0*F1)/(4.0D0*DIF*DIF)
      PRINT*,'Energy for lattice cycle ',NCOUNT,' is ',F1
      PRINT*,'Gradient wrt box length=',GRAD
      PRINT*,'Second derivative wrt box length=',SECOND
      NCOUNT=NCOUNT+1
      IF (DABS(GRAD/SECOND).GT.1.0D-4) THEN
         BOXLX=BOXLX-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         PRINT*,'Step=',-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         GOTO 10
      ELSE
         BOXLX=BOXLX-GRAD/SECOND
         PRINT*,'Step=',-GRAD/SECOND
         IF (DABS(GRAD/SECOND).GT.1.0D-6) GOTO 10
      ENDIF
      CUTOFF=CSAVE*BOXLX/BSAVE
      RETURN
      END
C
C********************************************************************
C
C Periodic bulk C60 energy only
C
C********************************************************************
C
      SUBROUTINE C60PE(N,X,ENERGY,BOXLX,BOXLY,BOXLZ,CUTOFF)
      IMPLICIT NONE
      INTEGER N, J1, J2, MX, MY, MZ
      DOUBLE PRECISION D
      PARAMETER(D=1.023349668D0)
      DOUBLE PRECISION X(3*N), ENERGY, TEMP, R(N,N), VEC(N,N,3),
     1                 BOXLX,BOXLY,BOXLZ,CUTOFF,DIST
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
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
            R(J2,J1)=DSQRT(DIST)
            R(J1,J2)=R(J2,J1)
10       CONTINUE
20    CONTINUE
C
C calculate the potential energy:
C
      ENERGY=0.0D0
      DO 12 J1=1,N
         DO 11 J2=J1+1,N
            TEMP=R(J2,J1)
            IF (TEMP.LT.CUTOFF) THEN
               ENERGY=ENERGY+4*(10*(-2/TEMP**9+(-2*D+TEMP)**(-9)+
     1        (2*D+TEMP)**(-9))/(D**2*TEMP)-75*(-2/TEMP**3
     2        +(-2*D+TEMP)**(-3)+(2*D+TEMP)**(-3))/(D**2*TEMP))
            ENDIF
11       CONTINUE
12    CONTINUE
      ENERGY=ENERGY/96.794820624738D0
C     PRINT*,'Energy for box length ',BOXLX,' is ',ENERGY

      RETURN
      END
