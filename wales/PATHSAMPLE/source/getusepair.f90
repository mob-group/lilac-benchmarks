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
SUBROUTINE GETUSEPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: NUSEPAIRS, USEPAIRSMIN, UMIN, NATOMS, DMIN1, DMIN2, NATTEMPT, NCPU, MAXBARRIER,  &
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
   NDIFF=1
   PAIRSTODO=NCPU*NPAIRFRQ
   IF (NPAIRFRQ.LT.1) PAIRSTODO=NATTEMPT*NCPU ! just one set of pairs unless we run out
   IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
   ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO))

10 CONTINUE
   minloop: DO J1=1,NUSEPAIRS-NDIFF
      J2=J1+NDIFF
      DO J3=1,NPAIRDONE
         IF ((PAIR1(J3).EQ.USEPAIRSMIN(J1)).AND.(PAIR2(J3).EQ.USEPAIRSMIN(J2))) CYCLE minloop ! do not repeat searches
         IF ((PAIR1(J3).EQ.USEPAIRSMIN(J2)).AND.(PAIR2(J3).EQ.USEPAIRSMIN(J1))) CYCLE minloop ! do not repeat searches
      ENDDO
      NAVAIL=NAVAIL+1
      DMIN1(NAVAIL)=USEPAIRSMIN(J1)
      DMIN2(NAVAIL)=USEPAIRSMIN(J2)
      IF (NAVAIL.EQ.PAIRSTODO) EXIT minloop

      IF (DEBUG) PRINT '(3(A,I8))','getusepair> connection ',NAVAIL,' pair ',USEPAIRSMIN(J1),' and ',USEPAIRSMIN(J2)
   ENDDO minloop

   IF (NAVAIL.LT.PAIRSTODO) THEN
      NDIFF=NDIFF+1
      IF (NDIFF.LE.NUSEPAIRS-1) GOTO 10
   ENDIF

   NAVAIL=MIN(NAVAIL,PAIRSTODO) 
   PRINT '(A,I8,A)','getusepair> sorted list of ',NAVAIL,' pairs'
   PRINT '(2I8)',(DMIN1(J1),DMIN2(J1),J1=1,NAVAIL)
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getusepair> No more candidate pairs of minima in getusepair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF
!
!  Check whether we already have a connection, and stop if we do.
!
KSUM(1:NMIN)=0.0D0
DO J1=1,NTS
   CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,DEADTS(J1))
   IF ((.NOT. DEADTS(J1)) .AND. (PLUS(J1).NE.MINUS(J1))) THEN
      KSUM(PLUS(J1))=KSUM(PLUS(J1))+EXP(KPLUS(J1))
      KSUM(MINUS(J1))=KSUM(MINUS(J1))+EXP(KMINUS(J1))
   ENDIF
ENDDO
DO J1=1,NMIN
   IF (KSUM(J1).GT.0.0D0) THEN
      KSUM(J1)=LOG(KSUM(J1))
   ENDIF
ENDDO

CALL MAKED(DMATMC,NCOL,NVAL,DEADTS,.TRUE.,ISA,ISB,KSUM)
NDISTA(1:NMIN)=1000000
DO J1=1,NMINA
   NDISTA(LOCATIONA(J1))=0
ENDDO
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN=100000
DMAX=0
NUNCONA=0
DO J1=1,NMIN
   IF (NDISTA(J1).EQ.0) CYCLE ! A minimum
   DO J2=1,NCOL(J1)
      IF (NDISTA(NVAL(J2,J1))+1.LT.NDISTA(J1)) THEN
         CHANGED=.TRUE.
         NDISTA(J1)=NDISTA(NVAL(J2,J1))+1
      ENDIF
   ENDDO
   IF ((NDISTA(J1).GT.DMAX).AND.(NDISTA(J1).NE.1000000)) DMAX=NDISTA(J1)
   IF (NDISTA(J1).LT.DMIN) DMIN=NDISTA(J1)
   IF (NDISTA(J1).EQ.1000000) NUNCONA=NUNCONA+1
ENDDO
IF (CHANGED) GOTO 5
PRINT '(3(A,I8))','getusepair> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONA
!
!  Calculate minimum number of steps of each minimum from the B set.
!
NDISTB(1:NMIN)=1000000
DO J1=1,NMINB
   NDISTB(LOCATIONB(J1))=0
ENDDO
NCYCLE=0
51    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN=100000
DMAX=0
NUNCONB=0
DO J1=1,NMIN
   IF (NDISTB(J1).EQ.0) CYCLE ! B minimum
   DO J2=1,NCOL(J1)
      IF (NDISTB(NVAL(J2,J1))+1.LT.NDISTB(J1)) THEN
         CHANGED=.TRUE.
         NDISTB(J1)=NDISTB(NVAL(J2,J1))+1
      ENDIF
   ENDDO
   IF ((NDISTB(J1).GT.DMAX).AND.(NDISTB(J1).NE.1000000)) DMAX=NDISTB(J1)
   IF (NDISTB(J1).LT.DMIN) DMIN=NDISTB(J1)
   IF (NDISTB(J1).EQ.1000000) NUNCONB=NUNCONB+1
ENDDO
IF (CHANGED) GOTO 51
PRINT '(3(A,I8))','getusepair> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONB
CHECKCONN=.FALSE.
IF (DIRECTION.EQ.'AB') THEN
   DO J1=1,NMINB
      IF (NDISTA(LOCATIONB(J1)).LT.1000000) THEN
         CHECKCONN=.TRUE.
         EXIT
      ENDIF
   ENDDO
ELSE
   DO J1=1,NMINA
      IF (NDISTB(LOCATIONA(J1)).LT.1000000) THEN
         CHECKCONN=.TRUE.
         EXIT
      ENDIF
   ENDDO
ENDIF
IF (.NOT.CHECKCONN) THEN
   PRINT '(A)','getusepair> There is no connection between the A and B regions'
ELSE
   PRINT '(A)','getusepair> Connection found between A and B regions'
   STOP
ENDIF
!
!  End of connection check.
!

7 CONTINUE ! Jump back to here if the previously selected pair of minima were already connected
NUSED=NUSED+1
NAVAIL=NAVAIL-1
IF (NAVAIL.EQ.0) THEN
   PRINT '(A)','getusepair> No more candidate pairs of minima in getusepair - quit'
 !  IF(SKIPPAIRST) STOP
 !  RETURN
    STOP
ENDIF
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)

! Check whether this particular pair is already connected. SKIPPAIRS allows us to skip over the pair if it is.
IF(SKIPPAIRST) THEN
    NDISTSTART(1:NMIN)=1000000
    NDISTSTART(MINS) = 0

    NCYCLE=0
52   CHANGED=.FALSE.
    NCYCLE=NCYCLE+1
    DMIN=100000
    DMAX=0
    NUNCONSTART=0
    DO J1=1,NMIN
        IF (NDISTSTART(J1).EQ.0) CYCLE ! Starting minimum
        DO J2=1,NCOL(J1)
            IF (NDISTSTART(NVAL(J2,J1))+1.LT.NDISTSTART(J1)) THEN
                CHANGED=.TRUE.
                NDISTSTART(J1)=NDISTSTART(NVAL(J2,J1))+1
            ENDIF
        ENDDO
        IF ((NDISTSTART(J1).GT.DMAX).AND.(NDISTSTART(J1).NE.1000000)) DMAX=NDISTSTART(J1)
        IF (NDISTSTART(J1).LT.DMIN) DMIN=NDISTSTART(J1)
        IF (NDISTSTART(J1).EQ.1000000) NUNCONSTART=NUNCONSTART+1
    ENDDO

    IF(NDISTSTART(MINF).LT.DMAX+1) THEN
        PRINT '(3(A,I8))','getusepair> Connectivity of start minimum converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONB       
        WRITE(*,*) "Minima ", MINS, "and", MINF, "are already connected. Skipping."
        GOTO 7
    ENDIF
    IF (CHANGED) GOTO 52
ENDIF


WRITE(*,'(5(A,I8))') 'getusepair> connecting minima ',MINS,' and ',MINF, ' pairs used=',  &
  &  NUSED,' remaining=',NAVAIL,' total pairs=',NPAIRDONE
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
CALL FLUSH(6)
READ(UMIN,REC=MINS) SPOINTS(1:NOPT)
READ(UMIN,REC=MINF) FPOINTS(1:NOPT)

END SUBROUTINE GETUSEPAIR
