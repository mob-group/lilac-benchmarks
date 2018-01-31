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
!  Subroutine to provide candidate pairs of minima based on the list
!  of NUSEPAIRS in array USEPAIRSMIN
!
SUBROUTINE GETCONNECTPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: NCONNECTPAIRS, CONNECTPAIRSMIN, UMIN, NATOMS, DMIN1, DMIN2, NATTEMPT, NCPU, MAXBARRIER,  &
  &               DEBUG, NPAIRFRQ, PAIR1, PAIR2, NPAIRFRQ, NPAIRDONE, MAXPAIRS, LOCATIONA, LOCATIONB, NCONNMAX, &
                  NTS, NMIN, NMINA, NMINB, DIRECTION, PLUS, MINUS, KPLUS, KMINUS, NCONN, &
  &               ETS, EMIN, SKIPPAIRST, NOPT
USE PORFUNCS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, J1, J2, J3, NDIFF
DOUBLE PRECISION SPOINTS(NOPT), FPOINTS(NOPT)
DOUBLE PRECISION DMATMC(NCONNMAX,NMIN), KSUM(NMIN)
INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN
INTEGER :: NDISTSTART(NMIN), NUNCONSTART ! sn402
LOGICAL DEADTS(NTS), ISA(NMIN), ISB(NMIN), CHANGED, CHECKCONN
INTEGER DMAX, NUNCONA, NUNCONB
DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0

IF (NAVAIL.EQ.0) THEN
!
! If called a second time we won't get any more candidate pairs because the routine
! will detect each pair has been tried before. 
!
   PAIRSTODO=NCONNECTPAIRS
   IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
   ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO))

   minloop: DO J1=1,NCONNECTPAIRS
      DO J3=1,NPAIRDONE
         IF ((PAIR1(J3).EQ.CONNECTPAIRSMIN(J1,1)).AND.(PAIR2(J3).EQ.CONNECTPAIRSMIN(J1,2))) CYCLE minloop ! do not repeat searches
         IF ((PAIR1(J3).EQ.CONNECTPAIRSMIN(J1,2)).AND.(PAIR2(J3).EQ.CONNECTPAIRSMIN(J1,1))) CYCLE minloop ! do not repeat searches
      ENDDO
      NAVAIL=NAVAIL+1
      DMIN1(NAVAIL)=CONNECTPAIRSMIN(J1,1)
      DMIN2(NAVAIL)=CONNECTPAIRSMIN(J1,2)

      IF (DEBUG) PRINT '(3(A,I8))','getconnectpair> connection ',NAVAIL,' pair ',CONNECTPAIRSMIN(J1,1),' and ',CONNECTPAIRSMIN(J1,2)
   ENDDO minloop

   NAVAIL=MIN(NAVAIL,PAIRSTODO) 
   PRINT '(A,I8,A)','getconnectpair> list of ',NAVAIL,' pairs'
   PRINT '(2I8)',(DMIN1(J1),DMIN2(J1),J1=1,NAVAIL)
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getconnectpair> No more candidate pairs of minima in getconnectpair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF

NUSED=NUSED+1
NAVAIL=NAVAIL-1
IF (NAVAIL.LT.0) THEN
   PRINT '(A)','getconnectpair> No more candidate pairs of minima in getusepair - quit'
    STOP
ENDIF

MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)

WRITE(*,'(5(A,I8))') 'getconnectpair> connecting minima ',MINS,' and ',MINF, ' pairs used=',  &
  &  NUSED,' remaining=',NAVAIL,' total pairs=',NPAIRDONE
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
CALL FLUSH(6)
READ(UMIN,REC=MINS) SPOINTS(1:NOPT)
READ(UMIN,REC=MINF) FPOINTS(1:NOPT)

END SUBROUTINE GETCONNECTPAIR
