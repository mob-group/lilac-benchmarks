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
      SUBROUTINE ROTCON(IT,IPRNT,IERR,PTEST)
C
C  COMPUTES ROTATIONAL CONSTANTS (IN CM-1), OPTIONALLY
C  PRINTS THEM AND RETURNS AN ERROR MESSAGE IF INERTIA
C  TENSOR PASSED IS NOT DIAGONAL.
C
      IMPLICIT NONE
      DOUBLE PRECISION IT(3,3),RC(3),PI,FACTOR,Z
      INTEGER NORD(3),IPRNT,IERR,IBOT,I,J
      DATA FACTOR /189.123013D0/
      DATA PI /3.141592654D0/
      LOGICAL PTEST

      
      Z=0.D0
      IBOT=0
      IERR=0
      DO I=1,3
         DO J=I+1,3
            Z=DABS(IT(I,J))+DABS(IT(J,I))+Z
         ENDDO
      ENDDO
      IF(Z.GT.1.D-9)THEN
       WRITE(6,100)
 100   FORMAT(' ***PROGRAM ERROR***, Inertia tensor not diagonal in ROTCON.')
       PRINT*,'Z=',Z
       IERR=1
       IBOT=1
      ELSE
       IBOT=1
      ENDIF
      DO 20 I=1,3
       IF(IT(I,I).GT.0)THEN
        RC(I)=FACTOR/IT(I,I)
        RC(I)=RC(I)/PI
       ELSE
        IBOT=2
       ENDIF
 20   CONTINUE
      CALL PIKSR2(3,RC,NORD)
      IF(IPRNT.NE.0)THEN
       IF (PTEST) WRITE(6,200)
       IF (PTEST) WRITE(6,300) (RC(J),J=IBOT,3),IT(1,1),IT(2,2),IT(3,3)
 200   FORMAT(' rotcon> Rotational constants (in cm-1) and principal moments of inertia: ')
 300   FORMAT((6(F15.5,1X)))
      ENDIF
      RETURN
      END        
