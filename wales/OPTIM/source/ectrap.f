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
C Pair potential for a 2D trapped ion cluster
C                                        
      SUBROUTINE ECTRAP (N, X, POTEL, C1, C2, C3)
      USE KEY
      IMPLICIT NONE 
      INTEGER N, I, J, J1, J2
      DOUBLE PRECISION X(3*N), K, POTEL, R, RO, C1, C2, C3
      PARAMETER (K=20.0D0)

C      C1=0.0D0
C      C2=0.0D0
C      DO 27 J1=1,N
C         C2=C2+X(3*J1-1)
C         C1=C1+X(3*J1-2)
C27    CONTINUE
C      C2=C2/N
C      C1=C1/N

C     DO I=1,N
C        X(3*(J1-1)+1)=X(3*(J1-1)+1)-C1
C        X(3*(J1-1)+2)=X(3*(J1-1)+2)-C2
C     ENDDO
C     C2=0.0D0
C     C1=0.0D0

      POTEL=0.0D0
      DO I=1,N
         J1=3*(I-1)+1
C        RO=(X(J1)-C1)**2 + (X(J1+1)-C2)**2
         RO=X(J1)**2 + X(J1+1)**2
         POTEL=POTEL+0.5D0*K*RO
         DO J=I+1,N
            J2=3*(J-1)+1
            R=(X(J1)-X(J2))**2+(X(J1+1)-X(J2+1))**2
            POTEL=POTEL+CHARGE(J)*CHARGE(I)/DSQRT(R)
         ENDDO
      ENDDO
      RETURN
      END
