!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!
!  Subroutine to provide candidate pairs of minima based on distance analysis
!  for pairs of minima.
!
SUBROUTINE CONNECTDIST(NAVAIL)
USE COMMONS, ONLY: NATOMS, UMIN, CONNECTMIN1, CONNECTMIN2, CONNECTDIST, DMIN1, DMIN2, TCONNECTDIST, TOPPOINTER, POINTERP, &
  &               POINTERM, NMIN, PLUS, MINUS, BULKT, TWOD, ZSYM, DMINMAX, DEBUG,INTERPCOSTFUNCTION
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, J1, J2, J3
INTEGER, ALLOCATABLE :: VINT(:)
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), ELAPSED, TNEW, DISTANCE

CALL CPU_TIME(ELAPSED)
ALLOCATE(VINT(DMINMAX))

NAVAIL=0
DO J1=CONNECTMIN1,NMIN
   READ(UMIN,REC=J1) SPOINTS(1:3*NATOMS)
   loop2: DO J2=MAX(CONNECTMIN2,J1+1),NMIN
      J3=TOPPOINTER(J2)  !  sets J3 to the ts connected to minimum J2 with the highest value
      IF (J3.LE.0) GOTO 555 
11    IF (PLUS(J3).EQ.J2) THEN
         IF (MINUS(J3).EQ.J1) CYCLE loop2 ! they are already connected
         J3=POINTERP(J3)
      ELSE
         IF (PLUS(J3).EQ.J1) CYCLE loop2 ! they are already connected
         J3=POINTERM(J3)
      ENDIF
      IF (J3.GT.0) GOTO 11 ! check the next ts connected to minimum J2
555   CONTINUE
      READ(UMIN,REC=J2) FPOINTS(1:3*NATOMS)
!     CALL MIND(SPOINTS,FPOINTS,NATOMS,DISTANCE,BULKT,TWOD,ZSYM)
      CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
      IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2, &
  &                                            RIGIDBODY,RMAT,INTERPCOSTFUNCTION)

      IF (DISTANCE.LE.CONNECTDIST) THEN
         NAVAIL=NAVAIL+1
         IF (NAVAIL.GT.DMINMAX) THEN
            VINT(1:DMINMAX)=DMIN1(1:DMINMAX)
            DEALLOCATE(DMIN1)
            ALLOCATE(DMIN1(2*DMINMAX))
            DMIN1(1:DMINMAX)=VINT(1:DMINMAX)
      
            VINT(1:DMINMAX)=DMIN2(1:DMINMAX)
            DEALLOCATE(DMIN2)
            ALLOCATE(DMIN2(2*DMINMAX))
            DMIN2(1:DMINMAX)=VINT(1:DMINMAX)
            DMINMAX=2*DMINMAX
            DEALLOCATE(VINT)
            ALLOCATE(VINT(DMINMAX))
         ENDIF
         DMIN1(NAVAIL)=J1
         DMIN2(NAVAIL)=J2
         IF (DEBUG) PRINT '(3(A,I8),A,G20.10)','connectdist> connection ',NAVAIL,' pair ',J1,' and ',J2,' distance=',DISTANCE
         IF (DMINMAX.GT.10000) THEN ! set a maximum number of candidates for each call to CONNECTDIST
            CONNECTMIN1=J1
            CONNECTMIN2=J2+1
            GOTO 777
         ENDIF
      ENDIF
   ENDDO loop2
   CONNECTMIN2=0 ! this is only to allow an efficient restart
ENDDO
CONNECTMIN1=NMIN+1
777 CONTINUE
PRINT '(A,I8)','connectdist> Number of available connection pairs to try=',NAVAIL

CALL CPU_TIME(TNEW)
TCONNECTDIST=TCONNECTDIST+TNEW-ELAPSED
DEALLOCATE(VINT)

END SUBROUTINE CONNECTDIST
