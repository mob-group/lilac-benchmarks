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
      SUBROUTINE ADM(Q)
      use porfuncs
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      INTEGER IEXIT, JBOT, NTIMES, ICOUNT, I, ICN, J
      DOUBLE PRECISION DIS, TOL, DIST, Q(3*NATOMS)
C
C  PRINTS LOWER DIAGONAL ATOMIC DISTANCE MATRIX 
C
      PARAMETER (TOL = 0.5D0)

       WRITE(6,1092)
1092   FORMAT(' Interatomic distance matrix ')
       JBOT=1
       NTIMES=1+(NATOMS-1)/5
       DO 20 ICOUNT=1,NTIMES
          WRITE(6,*)
          WRITE(6,142)(ZSYM(ICN),ICN=JBOT,MIN(NATOMS,JBOT+4))
142       FORMAT(18X,A3,4(10X,A3))
          WRITE(6,144)(ICN,ICN=JBOT,MIN(NATOMS,JBOT+4))
144       FORMAT(16X,:'[',I3,']',4(8X,:'[',I3,']'))
          DO 10 I=JBOT,NATOMS
             WRITE(6,143)ZSYM(I),I,(DIST(Q(3*i-2),Q(3*J-2)),J=JBOT
     &      ,MIN(I,JBOT+4))
143          FORMAT(T3,A3,'[',I3,']',5(3X,F10.5))
10     CONTINUE
       JBOT=JBOT+5
20     CONTINUE
C
C CHECK DISTANCES TO SEE IF ANY ARE TOO SHORT.  CALL EXIT IF BELOW
C  A CERTAIN TOLERANCE.
C
       IF (DCHECK) THEN
          DO I=1,NATOMS
             DO J=I+1,NATOMS
                DIS=DIST(Q(3*I-2),Q(3*J-2))
                IF(DIS.LT.TOL)THEN
C
C IF ONE IS A DUMMY, GO ON.
C
                    IF(ATMASS(I).LT.1.D-3.OR.ATMASS(J).LT.1.D-3) GOTO 12
                    WRITE(6,231)I,J,DIS,TOL
231                 FORMAT(T3,' Atoms ',i2,' and ',i2,
     1               ' are rather close.',
     &               /,T3,' Distance of ',f6.4,' is below threshold of '
     &                        ,f6.4,'.')
                    IEXIT=0
                ENDIF
12              CONTINUE
             ENDDO
          ENDDO
       ENDIF
       IF(IEXIT .eq. 1) then
          Write (6,*) ' *Inspect distance matrix and correct ',
     $       'Z-matrix.'
          STOP
       ENDIF
      RETURN
      END

