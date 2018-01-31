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
!  Subroutine to provide candidate pairs of minima based on the difference in Pfold values
!  and their distances. 
!
SUBROUTINE GETPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: UMIN, NATOMS, DMIN1, DMIN2, DIJINITT, NATTEMPT, NCPU, PSCALE, DSCALE, BULKT, TWOD, ZSYM, DEBUG, GPFOLD, NMIN, &
  &               NPAIRFRQ, PAIR1, PAIR2, NPAIRFRQ, NPAIRDONE, MAXPAIRS, PERMDIST, BOXLX, BOXLY, BOXLZ, RIGIDBODY, ANGLEAXIS, &
  &               INTERPCOSTFUNCTION, NOPT
USE PORFUNCS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, J1, NSTEPS, MAXSTEPS, NDUMMY, J2, J3, J4
DOUBLE PRECISION, ALLOCATABLE :: DISTLIST(:), PDIFFLIST(:)
DOUBLE PRECISION SPOINTS(NOPT), FPOINTS(NOPT), DISTANCE, RMAT(3,3), DIST2

10 CONTINUE
IF (NAVAIL.EQ.0) THEN
   PAIRSTODO=NCPU*NPAIRFRQ
   IF (NPAIRFRQ.LT.1) PAIRSTODO=NATTEMPT*NCPU ! just one set of pairs unless we run out
   IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
   ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO),DISTLIST(PAIRSTODO),PDIFFLIST(PAIRSTODO))
   DISTLIST(1:PAIRSTODO)=1.0D100
   DO J1=1,NMIN
!     IF (GPFOLD(J1).EQ.0.0D0) CYCLE
      READ(UMIN,REC=J1) (SPOINTS(J2),J2=1,NOPT)
      min2: DO J2=J1+1,NMIN
!        IF (GPFOLD(J2).EQ.0.0D0) CYCLE
!        IF (ABS(GPFOLD(J1)-GPFOLD(J2)).LT.PSCALE) CYCLE
         DO J3=1,NPAIRDONE
            IF ((PAIR1(J3).EQ.J1).AND.(PAIR2(J3).EQ.J2)) CYCLE min2 ! do not repeat searches
            IF ((PAIR1(J3).EQ.J2).AND.(PAIR2(J3).EQ.J1)) CYCLE min2 ! do not repeat searches
         ENDDO
         READ(UMIN,REC=J2) (FPOINTS(J3),J3=1,NOPT)
         CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
  &                                               DIST2,RIGIDBODY, &
  &                       RMAT,INTERPCOSTFUNCTION)
         NAVAIL=NAVAIL+1
         sortloop: DO J3=1,MIN(NAVAIL,PAIRSTODO) ! sort the shortest PAIRSTODO values
            IF (DISTANCE.LT.DISTLIST(J3)) THEN
               DO J4=MIN(NAVAIL,PAIRSTODO),J3+1,-1
                  DMIN1(J4)=DMIN1(J4-1)
                  DMIN2(J4)=DMIN2(J4-1)
                  DISTLIST(J4)=DISTLIST(J4-1)
                  PDIFFLIST(J4)=PDIFFLIST(J4-1)
               ENDDO
               DMIN1(J3)=J1
               DMIN2(J3)=J2
               DISTLIST(J3)=DISTANCE
               PDIFFLIST(J3)=ABS(GPFOLD(J1)-GPFOLD(J2))
               EXIT sortloop
            ENDIF
         ENDDO sortloop
         IF (DEBUG) PRINT '(3(A,I8),2(A,G20.10))','getpair> connection ',NAVAIL,' pair ',J1,' and ',J2,' diff=', &
  &                                            ABS(GPFOLD(J2)-GPFOLD(J1)),' distance=',DISTANCE
      ENDDO min2
   ENDDO
   NAVAIL=MIN(NAVAIL,PAIRSTODO) 
   PRINT '(A,I8,A)','getpair> sorted list of ',NAVAIL,' pairs'
   PRINT '(2I8,2G20.10)',(DMIN1(J1),DMIN2(J1),PDIFFLIST(J1),DISTLIST(J1),J1=1,NAVAIL)
   DEALLOCATE(DISTLIST,PDIFFLIST)
   IF (NAVAIL.LE.0) THEN
      PRINT '(A)','getpair> No more candidate pairs of minima in getpair - quit'
!     STOP
   ENDIF
   NUSED=0
ENDIF
NUSED=NUSED+1
NAVAIL=NAVAIL-1
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)
WRITE(*,'(5(A,I8))') 'getpair> connecting minima ',MINS,' and ',MINF, ' pairs used=',  &
  &    NUSED,' remaining=',NAVAIL,' total pairs=',NPAIRDONE
CALL FLUSH(6)
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
READ(UMIN,REC=MINS) SPOINTS(1:NOPT)
READ(UMIN,REC=MINF) FPOINTS(1:NOPT)

END SUBROUTINE GETPAIR
