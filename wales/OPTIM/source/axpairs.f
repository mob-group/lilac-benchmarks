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
C  Here we calculate the potential using the two and three-body
C  terms of Lennard-Jones plus Axilrod-Teller.
C                                        
C*************************************************************************
C
      SUBROUTINE AXPAIRS (N, X, P2, P3, POTEL, ZSTAR)
      IMPLICIT NONE 
      INTEGER N, J1, J2, I, J, IJ
      DOUBLE PRECISION X(3*N), POTEL, P2, P3, ZSTAR,
     1                 DIST(N,N), SDIST(N,N), DOT(N,N)
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
      DO 25 J1=1,N
         DO 35 J2=J1,N
            DOT(J2,J1)=X(3*(J2-1)+1)*X(3*(J1-1)+1)
     1                +X(3*(J2-1)+2)*X(3*(J1-1)+2)
     2                +X(3*(J2-1)+3)*X(3*(J1-1)+3)
            DOT(J1,J2)=DOT(J2,J1)
35       CONTINUE
25    CONTINUE
C
C  Calculate the energy
C  Try to make Sun Fortran compile the inner loop properly!
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J) THEN
               P2=P2+4.0D0/SDIST(I,J)**6-4.0D0/SDIST(I,J)**3
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J)) THEN
                     P3=P3+(1.0-3.0*
     1                   (DOT(J,I)+DOT(IJ,J)-DOT(IJ,I)-DOT(J,J))
     2                  *(DOT(J,I)+DOT(IJ,IJ)-DOT(IJ,J)-DOT(IJ,I))
     3                  *(DOT(I,I)+DOT(IJ,J)-DOT(J,I)-DOT(IJ,I))
     4                  / (DIST(I,J)*DIST(I,IJ)*DIST(J,IJ))**2 ) 
     5                  / (DIST(I,J)*DIST(I,IJ)*DIST(J,IJ))**3
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3*ZSTAR/6.0
      P2=P2/2.0D0
      POTEL=P2+P3
      RETURN
      END
