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
!  Subroutine to provide candidate pairs of minima based on the largest
!  barriers between regrouped free energy minima.
!
SUBROUTINE GETFREEPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: UMIN, NATOMS, DMIN1, DMIN2, NATTEMPT, NCPU, NMIN, &
  &               NPAIRFRQ, PAIR1, PAIR2, NPAIRFRQ, NPAIRDONE, MAXPAIRS, DMINMAX, DEBUG, LOCATIONA, LOCATIONB, BULKT, &
  &               ZSYM, TWOD, DIRECTION, PLUS, MINUS, NMINA, NMINB, EMIN, NTS, ETS, EUNTRAPTHRESH, EINC, DEBUG, ANGLEAXIS, &
  &               TOPPOINTER, POINTERP, POINTERM
USE PORFUNCS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, J1, J3, J4, CLOSEST(NMIN), MINVAL, OLDBASIN(NMIN), BASIN(NMIN), NBASIN, J2, NDONE
INTEGER MAXNEIGHBOURS, NP, NNEIGH, FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1)
INTEGER, ALLOCATABLE :: NEIGHBOURS(:), ITEMP(:)
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), BARRIER(NMIN), POINTS1(3*NATOMS), POINTS2(3*NATOMS), DISTANCE, RMAT(3,3), &
  &              HIGHESTTS, LOWESTTARG, DUMMY, ETHRESH
INTEGER, ALLOCATABLE :: VINT(:)
LOGICAL ISA(NMIN), ISB(NMIN), BASINT(NMIN), CHANGED, DONE(NMIN), MATCHED, OLDBASINT(NMIN)

ALLOCATE(VINT(DMINMAX))

10 CONTINUE
IF (NAVAIL.EQ.0) THEN
   PAIRSTODO=NCPU*NPAIRFRQ
   IF (NPAIRFRQ.LT.1) PAIRSTODO=NATTEMPT*NCPU ! just one set of pairs unless we run out
   IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1)
   IF (ALLOCATED(DMIN2)) DEALLOCATE(DMIN2)
   ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO))
   CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
   CALL REGROUPFREE2(.TRUE.,PAIRSTODO,FREEMINLIST,FREEMINPOINT,NAVAIL)
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getfpair> No more candidate pairs of minima in getfpair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF
NUSED=NUSED+1
NAVAIL=NAVAIL-1
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)
WRITE(*,'(5(A,I8))') 'getfpair> connecting minima ',MINS,' and ',MINF, ' pairs used=',  &
  &    NUSED,' remaining=',NAVAIL,' total pairs=',NPAIRDONE
CALL FLUSH(6)
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
READ(UMIN,REC=MINS) SPOINTS(1:3*NATOMS)
READ(UMIN,REC=MINF) FPOINTS(1:3*NATOMS)

END SUBROUTINE GETFREEPAIR

