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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  DIJKSTRA analysis subroutine
!
SUBROUTINE DIJKSTRA(NWORST,GETCANDIDATES,PAIRSTODO,NMINSAVE,MINMAP)
USE PORFUNCS
USE COMMONS
IMPLICIT NONE

LOGICAL, INTENT(IN) :: GETCANDIDATES
INTEGER J1, J2, J3, J4, PARENT(NMIN), JMINW, NPERM, OTHERMIN, J5, LJ1, LJ2, BESTEND, VERYBESTSTART, ENDMIN
INTEGER VERYBESTEND, BESTPARENT(NMIN), J6, NWORST, JTSWORST, JWORST, PAIRSTODO, NDUMMY, NDUMMY2
INTEGER, PARAMETER :: MAXPATH=10000 ! longest path length allowed - surely 10000 is enough?!
INTEGER VERYBESTPARENT(NMIN), NMINSTART, NMINEND, NDISTEND(NMIN), NSTEPS, TSPATHID(MAXPATH), MINPATHID(MAXPATH), NTSMIN, ISTAT
INTEGER, ALLOCATABLE :: LOCATIONSTART(:), LOCATIONEND(:)
LOGICAL DEADTS(NTS), PERMANENT(NMIN), ISA(NMIN), ISB(NMIN), CHANGED, ISSTART(NMIN), CHECKCONN, NOTNEW, SAVEGROUP(NMIN)
INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NDEAD, NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN
INTEGER DMAX, NUNCONA, NUNCONB, WORSTTS(NMIN), MINPATHID1, MINPATHIDNSP1
DOUBLE PRECISION MINWEIGHT, TMPWEIGHT, EMKSUM(NMIN), WEIGHT(NMIN), VERYBEST, BEST, ETSMIN, RMAT(3,3)
DOUBLE PRECISION TSWEIGHT(NMIN), TNEW, ELAPSED, EWORST, EWORSTANY, EREF(NMIN)
DOUBLE PRECISION DMATMC(NCONNMAX,NMIN), DUMMY, PFTOTALSTART, HUGESAVE, LOCALPOINTS(3*NATOMS), KSUM(NMIN)
DOUBLE PRECISION STARTPOINTS(3*NATOMS), FINISHPOINTS(3*NATOMS), DISTANCES
DOUBLE PRECISION EREFTMP, TSWEIGHTTMP, ALIGNPOINTS(3*NATOMS), DIST2, BARRIER(NMIN), RAT1, RAT2, KNPLUS, RAT3, RAT4, KMEAN
INTEGER DMIN1TMP, DMIN2TMP, WORSTTSTMP, CLOSEST(NMIN)
CHARACTER(LEN=4) SEGID(NATOMS), TYPE(NATOMS), RES(NATOMS)
INTEGER IRES(NATOMS), NMINSAVE, MINMAP(NMIN)
DOUBLE PRECISION, ALLOCATABLE :: ETA(:), THETA(:)
DOUBLE PRECISION WAITAB, WAITBA, THETAS, ETAS, THETASM1, ETASM1, LOWESTTARG
LOGICAL, ALLOCATABLE :: INCLUDEMIN(:)
LOGICAL FORWARDT, TRIEDSHIFT(NTS)
INTEGER, ALLOCATABLE :: NEWINDEX(:), MAXDOWNHILLTS(:)
CHARACTER(LEN=130) DUMMYSTRING
DOUBLE PRECISION DETS,DFVIBTS,DIXTS,DIYTS,DIZTS, CUT_UNDERFLOW
DOUBLE PRECISION, ALLOCATABLE :: MAXDOWNHILLB(:)
INTEGER DHORDERTS, DPLUS, DMINUS, NNMINA, NNMINB

! sf344> additions for dumping of pdb
CHARACTER(len=80) DUMMYLINE
CHARACTER(len=8) PRMFORMAT
CHARACTER(len=4) ATOMLABELS(NATOMS)
INTEGER PRMLINES
! end sf344

! csw34> added STATUS integer to allow calls to MYSYSTEM
INTEGER :: STATUS, INPCRDCOL

CALL CPU_TIME(ELAPSED)
NTRYING=0
TRIEDSHIFT(1:NTS)=.FALSE.
IF (ALLOCATED(SHIFTABLE)) DEALLOCATE(SHIFTABLE)
ALLOCATE(SHIFTABLE(NTS))
SHIFTABLE(1:NTS)=.FALSE.
975 CONTINUE ! branch for testing whether ts corresponding to large downhill barriers can be removed
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO
WORSTTS(1:PAIRSTODO)=-1

!  Original KSUM values are restored before return so that DIJKSTRA can be
!  called more than once.

CUT_UNDERFLOW=-300.0D0
CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.FALSE.,CUT_UNDERFLOW)

DO J1=1,NMIN
   EMKSUM(J1)=EXP(-KSUM(J1))
ENDDO

VERYBEST=-1.0D0
IF (DIJKSTRAWAITT) VERYBEST=HUGE(1.0D0)

!!!!!!!!!!!!!!!!!!!   Dijkstra calculation    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Find the single paths from a to b or b to a with the largest product of
!  branching probabilities. 
!
!  We need to allow transitions to all A and B minima - DMATMC is different from P^fold part
!  DMATMC((J2,J4) contains the branching probability from J4 to J2 - not used if DIJKSTRAWAITT is true.
!
CALL MAKED(DMATMC,NCOL,NVAL,DEADTS,.TRUE.,ISA,ISB,KSUM)
!
!  Check that the stationary point database is actually connected, 
!  and remove minima that lie in disjoint graphs.
!  Calculate minimum number of steps of each minimum from the A set.
!
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
IF (.NOT.SHANNONT) PRINT '(3(A,I8))','Dijkstra> steps to A region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONA
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
IF (.NOT.SHANNONT) PRINT '(3(A,I8))','Dijkstra> steps to B region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!
!  This could happen if disconnected minima lie in the A or B region
!
IF (NUNCONB.NE.NUNCONA) PRINT '(A)', &
&                   'Dijkstra> WARNING - number of disconnected minima from A and B is different'
!
!  Check that we actually have a connection between the A and B regions.
!  If not, STOP.
!
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
   PRINT '(A)','Dijkstra> There is no connection between the A and B regions'
   IF (SHANNONT) THEN
      NWORST=10000000
      RETURN
   ENDIF
!  OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
!  WRITE(1,'(I8)') TSATTEMPT(1:NTS)
!  CLOSE(1)
   IF (NTRYING.NE.0) THEN
      SHIFTABLE(NTRYING)=.FALSE.
      TRIEDSHIFT(NTRYING)=.TRUE.
      PRINT '(A,I8,A)','Dijkstra> transition state ',NTRYING,' cannot be shifted - retry'
      NTRYING=0
      GOTO 975
   ENDIF
   STOP
ENDIF
IF (NTRYING.NE.0) THEN
   SHIFTABLE(NTRYING)=.TRUE.
   TRIEDSHIFT(NTRYING)=.TRUE.
   PRINT '(A,I8,A)','Dijkstra> transition state ',NTRYING,' can be shifted - continue'
   NTRYING=0
ENDIF

IF (DIRECTION.EQ.'AB') THEN
   NMINSTART=NMINB
   NMINEND=NMINA
   IF (ALLOCATED(LOCATIONSTART)) DEALLOCATE(LOCATIONSTART)
   IF (ALLOCATED(LOCATIONEND)) DEALLOCATE(LOCATIONEND)
   ALLOCATE(LOCATIONSTART(NMINB),LOCATIONEND(NMINA))
   LOCATIONSTART(1:NMINB)=LOCATIONB(1:NMINB)
   LOCATIONEND(1:NMINA)=LOCATIONA(1:NMINA)
   NDISTEND(1:NMIN)=NDISTA(1:NMIN)
   ISSTART(1:NMIN)=ISB(1:NMIN)
   PFTOTALSTART=PFTOTALB
ELSEIF (DIRECTION.EQ.'BA') THEN
   NMINSTART=NMINA
   NMINEND=NMINB
   IF (ALLOCATED(LOCATIONSTART)) DEALLOCATE(LOCATIONSTART)
   IF (ALLOCATED(LOCATIONEND)) DEALLOCATE(LOCATIONEND)
   ALLOCATE(LOCATIONSTART(NMINA),LOCATIONEND(NMINB))
   LOCATIONSTART(1:NMINA)=LOCATIONA(1:NMINA)
   LOCATIONEND(1:NMINB)=LOCATIONB(1:NMINB)
   NDISTEND(1:NMIN)=NDISTB(1:NMIN)
   ISSTART(1:NMIN)=ISA(1:NMIN)
   PFTOTALSTART=PFTOTALA
ENDIF
!
!  Find largest weight for each B(A) minimum to all A(B) minima.
!
EREF(1:PAIRSTODO)=1.0D100
ENDMIN=NMINSTART
IF (GETCANDIDATES) THEN
   ENDMIN=NMIN
!
!  Use estimate of largest barrier on lowest energy path to minima in the
!  target region to speed up DIJPAIR by skipping over minima that have
!  barriers or ts energies that would not get onto the current sorted list,
!  thus avoiding a Dijkstra search for them.
!
   IF (UNTRAPT.AND.BARRIERSORT) CALL GETBARRIER(BARRIER,CLOSEST,LOWESTTARG)
ENDIF
NWORST=0
loopstart: DO J1=1,ENDMIN
   BEST=-1.0D0
   IF (DIJKSTRAWAITT) BEST=HUGE(1.0D0)
   IF (GETCANDIDATES) THEN
      LJ1=J1
   ELSE
      LJ1=LOCATIONSTART(J1)
   ENDIF
   IF (GETCANDIDATES.AND.(NWORST.GT.PAIRSTODO)) THEN ! filter out minima that are not expected to get on the list
      IF (.NOT.BARRIERSORT) THEN
         IF (EMIN(LJ1).GT.EREF(PAIRSTODO)) THEN
            PRINT '(A,I6,2(A,F15.5))','Dijkstra> Skipping minimum ',LJ1,' based on energy of minimum: ',EMIN(LJ1), &
   &                            ' and current highest value on list: ',EREF(PAIRSTODO)
            CYCLE loopstart 
         ENDIF
      ELSE
         IF (UNTRAPT) THEN
            IF (-BARRIER(LJ1).GT.EREF(PAIRSTODO)) THEN
               PRINT '(A,I6,2(A,F15.5))','Dijkstra> Skipping minimum ',LJ1,' based on estimated barrier ',BARRIER(LJ1), &
   &                            ' and current smallest value on list ',-EREF(PAIRSTODO)
               CYCLE loopstart 
            ENDIF
         ENDIF
      ENDIF
   ENDIF

   IF (NCOL(LJ1).LE.0) CYCLE loopstart
   IF (NDISTEND(LJ1).EQ.1000000) CYCLE loopstart
   IF (NDISTEND(LJ1).EQ.0) CYCLE loopstart
   WEIGHT(1:NMIN)=HUGE(1.0D0)
   HUGESAVE=WEIGHT(1)
   WEIGHT(LJ1)=0.0D0
   PERMANENT(1:NMIN)=.FALSE.
   PERMANENT(LJ1)=.TRUE.
   NPERM=1
   DO J2=1,NMIN
      IF (NCOL(J2).EQ.0) THEN
         IF (.NOT.PERMANENT(J2)) THEN
            PERMANENT(J2)=.TRUE.
            NPERM=NPERM+1
         ENDIF
      ENDIF
   ENDDO
   PARENT(1:NMIN)=0 ! parent is initially undefined
   J4=LJ1

   dijkstraloop: DO

      DO J2=1,NCOL(J4)
         OTHERMIN=NVAL(J2,J4)
         IF (.NOT.PERMANENT(OTHERMIN)) THEN
!
! This needs to be the branching probability out of the current minimum J4 into OTHERMIN
! via connection number J2.
! If DIJKSTRAWAITT is true then let the edge weight be the waiting time in the minimum
! we are coming from, which is stored in EMKSUM(J4)
!
            IF (DIJKSTRAWAITT) THEN
               TMPWEIGHT = EMKSUM(J4)
            ELSE
               IF (DMATMC(J2,J4).LE.0.0D0) THEN
                  TMPWEIGHT=HUGE(1.0D0)
               ELSE
                  TMPWEIGHT = -LOG(DMATMC(J2,J4))
               ENDIF
            ENDIF
            IF (TMPWEIGHT.LT.-0.001D0) THEN
               PRINT*,'error - tmpweight<0'
               STOP
            ENDIF
            IF (TMPWEIGHT+WEIGHT(J4).LT.WEIGHT(OTHERMIN)) THEN ! relax OTHERMIN
               WEIGHT(OTHERMIN)=TMPWEIGHT+WEIGHT(J4)
               PARENT(OTHERMIN)=J4
            ENDIF
         ENDIF
      ENDDO

      MINWEIGHT=HUGE(1.0D0)
      DO J2=1,NMIN
         IF (.NOT.PERMANENT(J2)) THEN
            IF (WEIGHT(J2).LT.MINWEIGHT) THEN
               MINWEIGHT=WEIGHT(J2)
               JMINW=J2
            ENDIF
         ENDIF
      ENDDO

      J4=JMINW
      PERMANENT(J4)=.TRUE.
      NPERM=NPERM+1
!     IF (DEBUG) PRINT '(A,I8)','NPERM=',NPERM
      PRINT '(A,2I8,G20.10)','permanent minimum J4,NPERM,WEIGHT=',J4,NPERM,WEIGHT(J4)

      IF (NPERM.EQ.NMIN) EXIT DIJKSTRALOOP

   ENDDO dijkstraloop
!
!  The next block should be in a separate subroutine for finding pairs of minima
!  to connect.
!
!  We only want one candidiate ts from a given non-endpoint minimum to the endpoint set,
!  otherwise dimensions will get out of hand. EWORSTANY helps us to identify this one.
!
   EWORSTANY=1.0D100
   DO J2=1,NMINEND
      LJ2=LOCATIONEND(J2)
!     PRINT '(A,I8,L5,I8,G20.10)','LJ2,ISSTART,NCOL,WEIGHT=',LJ2,ISSTART(LJ2),NCOL(LJ2),WEIGHT(LJ2)
      IF (ISSTART(LJ2)) CYCLE
      IF (NCOL(LJ2).EQ.0) CYCLE
      IF (WEIGHT(LJ2).EQ.HUGESAVE) CYCLE ! end min LJ2 is disconnected from start min LJ1, 
                                         ! hence from other end minima.
      DUMMY=1.0D0
      IF (DIJKSTRAWAITT) DUMMY=0.0D0
      J5=LJ2
      NSTEPS=0
      DO 
!
!  First identify the branching probability to calculate the weight of
!  this path.
!
!  Need DMATMC entry for branching probability from PARENT to J5.
!  This is simply a waiting time if DIJKSTRAWAITT is true
!
         IF (DIJKSTRAWAITT) THEN
            DUMMY=DUMMY+EMKSUM(PARENT(J5))
         ELSE
            DO J3=1,NCOL(PARENT(J5))
               IF (NVAL(J3,PARENT(J5)).EQ.J5) THEN
                  DUMMY=DUMMY*DMATMC(J3,PARENT(J5))
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         IF (DEBUG) PRINT '(A,2I8)','J5,parent(J5)=',J5,parent(J5)
!
!  Now identify the lowest transition state between the two minima.
!
         NSTEPS=NSTEPS+1
         ETSMIN=1.0D100
         J6=TOPPOINTER(J5) ! sets J6 to the transition state connected to J5 with the
                           ! highest number that isn't DEADTS
         pointa: DO
            IF (PLUS(J6).EQ.J5) THEN
               IF ((MINUS(J6).EQ.PARENT(J5)).AND.(.NOT.DEADTS(J6))) THEN
                  IF (ETS(J6).LT.ETSMIN) THEN
                     ETSMIN=ETS(J6)
                     NTSMIN=J6
                  ENDIF
               ENDIF
               J6=POINTERP(J6)
            ELSE IF (MINUS(J6).EQ.J5) THEN
               IF ((PLUS(J6).EQ.PARENT(J5)).AND.(.NOT.DEADTS(J6))) THEN
                  IF (ETS(J6).LT.ETSMIN) THEN
                     ETSMIN=ETS(J6)
                     NTSMIN=J6
                  ENDIF
               ENDIF
               J6=POINTERM(J6)
            ENDIF
            IF (J6.LT.0) EXIT pointa
         ENDDO pointa

         IF (NSTEPS.GT.MAXPATH) THEN
            PRINT '(A,I8)','ERROR in Dijkstra, path is longer than MAXPATH=',MAXPATH
            STOP
         ENDIF

         TSPATHID(NSTEPS)=NTSMIN
         IF (DEADTS(NTSMIN)) THEN
            PRINT '(A,2I8,L5)','NTSMIN,TSATTEMPT,DEADTS=',NTSMIN,TSATTEMPT(NTSMIN),DEADTS(NTSMIN)
            PRINT '(A)','Dijkstra> ERROR - this should never happen'
            STOP
         ENDIF
         MINPATHID(NSTEPS)=J5
!
!  End of TS identification.
!
         J5=PARENT(J5)
         IF (J5.EQ.LJ1) EXIT
         IF (J5.EQ.0) EXIT
      ENDDO

      MINPATHID(NSTEPS+1)=LJ1
      IF (DEBUG) THEN
         IF (DIJKSTRAWAITT) THEN
            IF (.NOT.SHANNONT) PRINT '(A,I6,A,A1,A,I6,2(A,G20.10),A,I8)', &
   &           'Dijkstra> Best path for minimum ',LJ2,' and ', DIRECTION(2:2),' minimum ',LJ1,' sum tau=', &
   &            WEIGHT(LJ2)+EMKSUM(LJ1),' sum*prob=',(WEIGHT(LJ2)+EMKSUM(LJ1))*EXP(PFMIN(LJ1)-PFTOTALSTART),' steps=',NSTEPS
           ELSE
            IF (.NOT.SHANNONT) PRINT '(A,I6,A,A1,A,I6,A,G20.10,A,I8)', &
 &           'Dijkstra> Best path for minimum ',LJ2,' and ', DIRECTION(2:2),' minimum ',LJ1,' k^SS=', &
 &            EXP(-WEIGHT(LJ2))*EXP(PFMIN(LJ1)-PFTOTALSTART)/EMKSUM(LJ1),' steps=',NSTEPS
           ENDIF
      ENDIF
!
!  Summarise the best path for given LJ1 and LJ2
!
      J5=LJ2
      EWORST=-1.0D100
      DO J6=1,NSTEPS
         IF (DEBUG) PRINT '(3(I8,F20.10))',J5,EMIN(J5),TSPATHID(J6),ETS(TSPATHID(J6)),PARENT(J5),EMIN(PARENT(J5))
!
! Highest ts condition:
!
         IF (ETS(TSPATHID(J6)).GT.EWORST) THEN
            IF (TSATTEMPT(TSPATHID(J6)).LT.MAXATTEMPT) THEN
               JWORST=J6
               JTSWORST=TSPATHID(J6)
               EWORST=ETS(JTSWORST)
            ENDIF
         ENDIF
!
! Barrier height condition:
!
!        IF (ETS(TSPATHID(J6))-EMIN(PARENT(J5)).GT.EWORST) THEN
!           IF (TSATTEMPT(TSPATHID(J6)).LT.MAXATTEMPT) THEN
!              JWORST=J6
!              JTSWORST=TSPATHID(J6)
!              EWORST=ETS(TSPATHID(J6))-EMIN(PARENT(J5))
!           ENDIF
!        ENDIF
!
!        PRINT '(5(A,I6),3(A,F15.5))','step=',J6,' min1=',J5,' ts=',TSPATHID(J6),' min2=',PARENT(J5),' JWORST=',JWORST, &
!    &              'eworst=',EWORST,' ETS=',ETS(TSPATHID(J6)),' bar=',ETS(TSPATHID(J6))-EMIN(PARENT(J5))
         J5=PARENT(J5)
      ENDDO
      IF (EWORST.GT.-1.0D100) THEN ! we might have tried all of transition states MAXATTEMPT times
         IF (EWORST.LT.EWORSTANY) THEN ! select the lowest of the worst transition states for the given
                                       ! starting minimum with any of the final state minima.
            EWORSTANY=EWORST
!
!  Choose the pair of minima TSATTEMPT(JTSWORST) steps away from the worst
!  transition state, without going beyond the ends of the path.
!
            DMIN1TMP=MINPATHID(MAX(1,JWORST-TSATTEMPT(JTSWORST)))
            DMIN2TMP=MINPATHID(MIN(NSTEPS+1,JWORST+TSATTEMPT(JTSWORST)+1))
            EREFTMP=EMIN(LJ1)
            IF (BARRIERSORT) EREFTMP=-(ETS(JTSWORST)-EMIN(LJ1)) ! sort is done lowest to highest, so change sign here
                                                                ! to get biggest barriers
            WORSTTSTMP=JTSWORST
            TSWEIGHTTMP=EXP(-WEIGHT(LJ2))*EXP(PFMIN(LJ1)-PFTOTALB)/EMKSUM(LJ1)
            MINPATHID1=MINPATHID(1)
            MINPATHIDNSP1=MINPATHID(NSTEPS+1)
         ENDIF
      ENDIF
!
!  End of summary
!
      IF (DUMMY.NE.0.0D0) THEN
         IF (DIJKSTRAWAITT) THEN
            IF (DABS((WEIGHT(LJ2)-DUMMY)/DUMMY).GT.1.0D-2) THEN
               PRINT '(A,G20.10)','WARNING - alternative sum of MFPTs=',DUMMY
            ENDIF
         ELSE
            IF (DABS((EXP(-WEIGHT(LJ2))-DUMMY)/DUMMY).GT.1.0D-2) THEN
               PRINT '(A,G20.10)','WARNING - alternative product of branching probabilities=',DUMMY
            ENDIF
         ENDIF
      ENDIF
!
!  Note that VERYBEST contains the conditional probability and is divided by the waiting time.
!  We could get the contribution to NSS rate constants if we did some short KMC runs to get
!  the appropriate waiting time.
!
      IF (DIJKSTRAWAITT) THEN
         IF ((DUMMY+EMKSUM(LJ1))*EXP(PFMIN(LJ1)-PFTOTALSTART).LT.VERYBEST) THEN
            VERYBEST=(DUMMY+EMKSUM(LJ1))*EXP(PFMIN(LJ1)-PFTOTALSTART)
            VERYBESTEND=LJ2
            VERYBESTSTART=LJ1
            VERYBESTPARENT(1:NMIN)=PARENT(1:NMIN)
         ENDIF
         IF (DUMMY.LT.BEST) THEN
            BEST=DUMMY
            BESTEND=LJ2
            BESTPARENT(1:NMIN)=PARENT(1:NMIN)
         ENDIF
      ELSE
         IF (DUMMY*EXP(PFMIN(LJ1)-PFTOTALSTART)/EMKSUM(LJ1).GT.VERYBEST) THEN
            VERYBEST=DUMMY*EXP(PFMIN(LJ1)-PFTOTALSTART)/EMKSUM(LJ1)
            VERYBESTEND=LJ2
            VERYBESTSTART=LJ1
            VERYBESTPARENT(1:NMIN)=PARENT(1:NMIN)
         ENDIF
         IF (DUMMY.GT.BEST) THEN
            BEST=DUMMY
            BESTEND=LJ2
            BESTPARENT(1:NMIN)=PARENT(1:NMIN)
         ENDIF
      ENDIF
   ENDDO
   IF (GETCANDIDATES.AND.(EWORSTANY.LT.1.0D100)) THEN
      NWORST=NWORST+1
!
!  Don;t retry for this transition state any more if we reached the ends of the path in question.
!
      IF ((DMIN1TMP.EQ.MINPATHID1).AND. (DMIN2TMP.EQ.MINPATHIDNSP1)) TSATTEMPT(WORSTTSTMP)=MAXATTEMPT
      PRINT '(5(A,I8))','Dijkstra> Candidate end points for minimum ',LJ1,' and product minimum ',MINPATHID1, &
   &                    ' are ',DMIN1TMP,' and ',DMIN2TMP,' for highest ts ',WORSTTSTMP
!
!  Put the new pair of endpoints into the list sorted according to increasing energy.
!  We are sorting so that we seek connections for the lowest minima or highest barrier first.
!
!  Check that we don't already have this pair. Use the lowest value of EREF for sorting (note
!  that this is the negative barrier if sorting is done on barrier heights, so that we can
!  still sort from smallest to largest). We may also need to move this pair up the list
!  according to EREF.
!
      NOTNEW=.FALSE.
      DO J4=1,MIN(NWORST,PAIRSTODO)
         IF ((DMIN1TMP.EQ.DMIN1(J4)).AND.(DMIN2TMP.EQ.DMIN2(J4))) THEN
            NOTNEW=.TRUE.
            IF (EREFTMP.LT.EREF(J4)) EREF(J4)=EREFTMP
            NDUMMY=J4
            DO
               IF ((EREF(NDUMMY).GT.EREF(NDUMMY-1)).OR.(NDUMMY.EQ.1)) EXIT
               NDUMMY2=DMIN1(NDUMMY-1)
               DMIN1(NDUMMY-1)=DMIN1(NDUMMY)
               DMIN1(NDUMMY)=NDUMMY2
               NDUMMY2=DMIN2(NDUMMY-1)
               DMIN2(NDUMMY-1)=DMIN2(NDUMMY)
               DMIN2(NDUMMY)=NDUMMY2
               DUMMY=EREF(NDUMMY-1)
               EREF(NDUMMY-1)=EREF(NDUMMY)
               EREF(NDUMMY)=DUMMY
               NDUMMY2=WORSTTS(NDUMMY-1)
               WORSTTS(NDUMMY-1)=WORSTTS(NDUMMY)
               WORSTTS(NDUMMY)=NDUMMY2
               DUMMY=TSWEIGHT(NDUMMY-1)
               TSWEIGHT(NDUMMY-1)=TSWEIGHT(NDUMMY)
               TSWEIGHT(NDUMMY)=DUMMY
               NDUMMY=NDUMMY-1
            ENDDO
            EXIT
         ELSE IF ((DMIN1TMP.EQ.DMIN2(J4)).AND.(DMIN2TMP.EQ.DMIN1(J4))) THEN
            NOTNEW=.TRUE.
            IF (EREFTMP.LT.EREF(J4)) EREF(J4)=EREFTMP
            NDUMMY=J4
            DO
               IF ((EREF(NDUMMY).GT.EREF(NDUMMY-1)).OR.(NDUMMY.EQ.1)) EXIT
               NDUMMY2=DMIN1(NDUMMY-1)
               DMIN1(NDUMMY-1)=DMIN1(NDUMMY)
               DMIN1(NDUMMY)=NDUMMY2
               NDUMMY2=DMIN2(NDUMMY-1)
               DMIN2(NDUMMY-1)=DMIN2(NDUMMY)
               DMIN2(NDUMMY)=NDUMMY2
               DUMMY=EREF(NDUMMY-1)
               EREF(NDUMMY-1)=EREF(NDUMMY)
               EREF(NDUMMY)=DUMMY
               NDUMMY2=WORSTTS(NDUMMY-1)
               WORSTTS(NDUMMY-1)=WORSTTS(NDUMMY)
               WORSTTS(NDUMMY)=NDUMMY2
               DUMMY=TSWEIGHT(NDUMMY-1)
               TSWEIGHT(NDUMMY-1)=TSWEIGHT(NDUMMY)
               TSWEIGHT(NDUMMY)=DUMMY
               NDUMMY=NDUMMY-1
            ENDDO
            EXIT
         ENDIF
      ENDDO
      IF (.NOT.NOTNEW) THEN 
         sortloop: DO J3=1,MIN(NWORST,PAIRSTODO)
            IF (EREFTMP.LT.EREF(J3)) THEN
!              PRINT '(A,I6,2G20.10)','J3,EREF(J3),EREFTMP=',J3,EREF(J3),EREFTMP
               IF (WORSTTS(MIN(NWORST,PAIRSTODO)).GT.0) THEN
!                 PRINT *,'MIN(NWORST,PAIRSTODO)=',MIN(NWORST,PAIRSTODO)
!                 PRINT *,'WORSTTS(MIN(NWORST,PAIRSTODO))=',WORSTTS(MIN(NWORST,PAIRSTODO))
                  TSATTEMPT(WORSTTS(MIN(NWORST,PAIRSTODO)))=TSATTEMPT(WORSTTS(MIN(NWORST,PAIRSTODO)))-1
               ENDIF
               DO J4=MIN(NWORST,PAIRSTODO),J3+1,-1
                  DMIN1(J4)=DMIN1(J4-1)
                  DMIN2(J4)=DMIN2(J4-1)
                  EREF(J4)=EREF(J4-1)
                  WORSTTS(J4)=WORSTTS(J4-1)
                  TSWEIGHT(J4)=TSWEIGHT(J4-1)
               ENDDO
               DMIN1(J3)=DMIN1TMP
               DMIN2(J3)=DMIN2TMP
               EREF(J3)=EREFTMP
               WORSTTS(J3)=WORSTTSTMP
               TSWEIGHT(J3)=TSWEIGHTTMP
               TSATTEMPT(WORSTTSTMP)=TSATTEMPT(WORSTTSTMP)+1
!              PRINT '(A)','summary of current list, min1, min2, ts, ts attempts, eref:'
!              PRINT '(4I6,F12.2)',(DMIN1(J4),DMIN2(J4),WORSTTS(J4),TSATTEMPT(WORSTTS(J4)),EREF(J4),J4=1,MIN(NWORST,PAIRSTODO))
               EXIT sortloop
            ENDIF
         ENDDO sortloop
!        PRINT '(A,32G20.10)','EREF: ',(EREF(J3),J3=1,MIN(NWORST,PAIRSTODO))
      ENDIF
   ENDIF
!  PRINT '(A)','summary of current list, min1, min2, ts, ts attempts, eref:'
!  PRINT '(4I6,F12.2)',(DMIN1(J4),DMIN2(J4),WORSTTS(J4),TSATTEMPT(WORSTTS(J4)),EREF(J4),J4=1,MIN(NWORST,PAIRSTODO))
! 
!  Summarise the best path for any A(B) and this particular starting minimum (LJ1)
!
   IF (DIJKSTRAWAITT) THEN
      IF (.NOT.SHANNONT) &
  &       PRINT '(A,I8,A,A1,A,A1,A,A1,G20.10,A,G20.10)','Dijkstra> Best path for min ',LJ1,' and any ',DIRECTION(1:1), &
  &         ' minimum, sum tau ', &
  &         DIRECTION(1:1),'<-',DIRECTION(2:2),BEST+EMKSUM(LJ1),' sum*prob=',(BEST+EMKSUM(LJ1))*EXP(PFMIN(LJ1)-PFTOTALSTART)
   ELSE
      IF (.NOT.SHANNONT) &
  &     PRINT '(A,I8,A,A1,A,A1,A,A1,G20.10)','Dijkstra> Best path for min ',LJ1,' and any ',DIRECTION(1:1),' minimum, k^SS ', &
  &         DIRECTION(1:1),'<-',DIRECTION(2:2),BEST*EXP(PFMIN(LJ1)-PFTOTALSTART)/EMKSUM(LJ1)
   ENDIF

   J5=BESTEND
   IF (.NOT.SHANNONT) WRITE(*,'(I8)',ADVANCE='NO') J5
   DO
     J5=BESTPARENT(J5)
     IF (.NOT.SHANNONT) WRITE(*,'(I8)',ADVANCE='NO') J5
     IF (J5.EQ.LJ1) EXIT
     IF (J5.EQ.0) EXIT
   ENDDO
   IF (.NOT.SHANNONT) PRINT*

ENDDO loopstart
! 
!  Summarise the best path for any A(B) and any B(A)
!
IF (.NOT.SHANNONT) THEN
   IF (.NOT.GETCANDIDATES) THEN
      PRINT '(2(A,A1),A)','Dijkstra> Best path between any ',DIRECTION(2:2),' minimum and any ',DIRECTION(1:1),' minimum:'
   ELSE
      PRINT '(A,A1,A)','Dijkstra> Best path between any ',DIRECTION(1:1),' product minimum and any other minimum:'
   ENDIF
ENDIF
J5=VERYBESTEND
IF (.NOT.SHANNONT) WRITE(*,'(I8)',ADVANCE='NO') J5
CALL FLUSH(6)
NSTEPS=0
DO
  NSTEPS=NSTEPS+1
!
!  Find the transition state with the lowest barrier between minima
!  J5 and VERYBESTPARENT(J5).
!
  J1=TOPPOINTER(J5) ! sets J1 to the transition state connected to J5 with the
                       ! highest number
  ETSMIN=1.0D100
  point: DO
     IF (PLUS(J1).EQ.J5) THEN
        IF (MINUS(J1).EQ.VERYBESTPARENT(J5)) THEN
           IF (ETS(J1).LT.ETSMIN) THEN
              ETSMIN=ETS(J1)
              NTSMIN=J1
           ENDIF
        ENDIF
        J1=POINTERP(J1)
     ELSE IF (MINUS(J1).EQ.J5) THEN
        IF (PLUS(J1).EQ.VERYBESTPARENT(J5)) THEN
           IF (ETS(J1).LT.ETSMIN) THEN
              ETSMIN=ETS(J1)
              NTSMIN=J1
           ENDIF
        ENDIF
        J1=POINTERM(J1)
     ENDIF
    IF (J1.LT.0) EXIT point
  ENDDO point
  TSPATHID(NSTEPS)=NTSMIN

  J5=VERYBESTPARENT(J5)
  IF (.NOT.SHANNONT) WRITE(*,'(I8)',ADVANCE='NO') J5
  IF (J5.EQ.VERYBESTSTART) EXIT
  IF (J5.EQ.0) EXIT
ENDDO
IF (.NOT.SHANNONT) PRINT*
IF (.NOT.GETCANDIDATES) THEN
   IF (.NOT.DIJKSTRAWAITT) THEN
      IF (.NOT.SHANNONT) PRINT '(2(A,A1),A,G20.10,A,I6,A)','Dijkstra> Largest contribution to SS rate constant ', &
  &     DIRECTION(1:1),'<-',DIRECTION(2:2), &
  &          ' for any A and B is ',VERYBEST,' for ',NSTEPS,' transition states:'
   ELSE
      IF (.NOT.SHANNONT) PRINT '(2(A,A1),A,G20.10,A,I6,A)','Dijkstra> Smallest sum tau * prob ', &
  &     DIRECTION(1:1),'<-',DIRECTION(2:2), &
  &          ' for any A and B is ',VERYBEST,' for ',NSTEPS,' transition states:'
   ENDIF
ELSE
   IF (DIJKSTRAWAITT) THEN
      IF (.NOT.SHANNONT) PRINT '(2(A,A1),A,G20.10,A,I6,A)','Dijkstra> Smallest sum tau * prob ', &
  &     DIRECTION(1:1),'<-',DIRECTION(2:2), &
  &          ' for any A and B is ',VERYBEST,' for ',NSTEPS,' transition states:'
   ELSE
      IF (.NOT.SHANNONT) PRINT '(A,A1,A,G20.10,A,I6,A)','Dijkstra> Largest contribution to SS rate constant ', &
  &     DIRECTION(1:1),'<-any other minimum is ',VERYBEST,' for ',NSTEPS,' transition states:'
   ENDIF
ENDIF

!!!!!!!!!!!!!!!!!!!   end Dijkstra calculation    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Now summarise the path corresponding to the steps for the VERYBEST path
!  Note that this path is printed backwards!
!
IF (.NOT.SHANNONT) THEN
   IF (DIRECTION.EQ.'AB') PRINT '(A)','Dijkstra> Note that path is printed backwards starting with A, ending with B'
   IF (DIRECTION.EQ.'BA') PRINT '(A)','Dijkstra> Note that path is printed backwards starting with B, ending with A'
   PRINT '(A)','                    E+                          Ets                         E-'
ENDIF
OPEN(UNIT=1,FILE='Epath',STATUS='UNKNOWN')

J5=VERYBESTEND
IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.REGROUPPERSISTT.OR.NOPOINTS)) &
   & READ(UMIN,REC=VERYBESTEND) (STARTPOINTS(J2),J2=1,3*NATOMS)
DO J1=1,NSTEPS
   IF (.NOT.SHANNONT) THEN
      PRINT '(3(I8,F20.10))',J5,EMIN(J5),TSPATHID(J1),ETS(TSPATHID(J1)),VERYBESTPARENT(J5), &
   &                      EMIN(VERYBESTPARENT(J5))
      WRITE(1,'(I8,G20.10,I8)') 2*J1-1,EMIN(J5),J5
      WRITE(1,'(I8,G20.10,I8)') 2*J1,ETS(TSPATHID(J1)),TSPATHID(J1)
      IF (J1.EQ.NSTEPS) WRITE(1,'(I8,G20.10,I8)') 2*NSTEPS+1,EMIN(VERYBESTPARENT(J5)),VERYBESTPARENT(J5)
      IF ((.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.REGROUPPERSISTT.OR.NOPOINTS)).AND.(J1.EQ.NSTEPS)) &
  &                 READ(UMIN,REC=VERYBESTPARENT(J5)) (FINISHPOINTS(J2),J2=1,3*NATOMS)
   ENDIF

   J5=VERYBESTPARENT(J5)
ENDDO
CLOSE(1)
!
! Identify the highest downhill barriers.
!
FORWARDT=.TRUE.
IF (EMIN(VERYBESTEND).LT.EMIN(J5)) FORWARDT=.FALSE.
IF (ALLOCATED(MAXDOWNHILLB)) DEALLOCATE(MAXDOWNHILLB)
ALLOCATE(MAXDOWNHILLB(NSTEPS))
IF (ALLOCATED(MAXDOWNHILLTS)) DEALLOCATE(MAXDOWNHILLTS)
ALLOCATE(MAXDOWNHILLTS(NSTEPS))
MAXDOWNHILLB(1:NSTEPS)=-1.0D100
MAXDOWNHILLTS(1:NSTEPS)=1
J5=VERYBESTEND
DO J1=1,NSTEPS
   IF (FORWARDT) THEN
      DUMMY=ETS(TSPATHID(J1))-EMIN(J5)
   ELSE
      DUMMY=ETS(TSPATHID(J1))-EMIN(VERYBESTPARENT(J5))
   ENDIF
   DO J2=1,NSTEPS
      IF (DUMMY.GT.MAXDOWNHILLB(J2)) THEN
         DO J3=NSTEPS,J2+1,-1
            MAXDOWNHILLB(J3)=MAXDOWNHILLB(J3-1)
            MAXDOWNHILLTS(J3)=MAXDOWNHILLTS(J3-1)
         ENDDO
         MAXDOWNHILLB(J2)=DUMMY
         MAXDOWNHILLTS(J2)=TSPATHID(J1)
         EXIT
      ENDIF
   ENDDO
   J5=VERYBESTPARENT(J5)
ENDDO
IF (.NOT.SHANNONT) PRINT '(A)','Dijkstra> Ordered downhill barriers,    ts        barrier'
DO J1=1,NSTEPS
   IF (MAXDOWNHILLB(J1).LT.-1.0D99) EXIT
   IF (.NOT.SHANNONT) PRINT '(A,I8,G20.10)','                                   ',MAXDOWNHILLTS(J1),MAXDOWNHILLB(J1)
ENDDO
DO J1=1,NSTEPS
   IF ((MAXDOWNHILLB(J1).GT.MAXDOWNBARRIER).AND.(.NOT.TRIEDSHIFT(MAXDOWNHILLTS(J1)))) THEN
      NTRYING=MAXDOWNHILLTS(J1)
      PRINT '(A,I8,A)','Dijkstra> downhill barrier for transition state ',NTRYING,' is over threshold - trying to shift'
      GOTO 975
   ENDIF
ENDDO
DO J1=1,NTS
   IF (TRIEDSHIFT(J1).AND.SHIFTABLE(J1)) PRINT '(A,I8,A)','Dijkstra> transition state ',J1,' was shifted successfully'
ENDDO
IF (.NOT.DIJPAIRT) THEN
   IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST.OR.REGROUPPERSISTT.OR.NOPOINTS.OR.SHANNONT)) PRINT '(A,I6)', &
  &          'Dijkstra> dumping minima and transition state coordinates to file redopoints, steps=',NSTEPS

   IF (CHARMMT) THEN
!
! Figure out how many columns are present to determine if the SEGID is specified
!
      CALL MYSYSTEM(STATUS,DEBUG,"sed -n '2p' input.crd | awk '{print NF}' > ncol")
      OPEN(UNIT=4,FILE='ncol',STATUS='OLD')
      READ(4,*) INPCRDCOL
      CLOSE(4)
      CALL MYSYSTEM(STATUS,DEBUG,"rm ncol") 
!
! Open input.crd to read system information
!
      OPEN(UNIT=4,FILE='input.crd',STATUS='OLD')
      READ(4,*) NDUMMY
!
! There will be a problem here if input.crd does not have the last column corresponding to SEGID
!
      DO J1=1,NATOMS
         IF(INPCRDCOL.GT.7) THEN
! If the input.crd file has more than 7 columns, we assume that the SEGID is present
            READ(4,*) NDUMMY,IRES(J1),RES(J1),TYPE(J1),DUMMY,DUMMY,DUMMY,SEGID(J1)
         ELSE
! If the SEGID is missing, we just set it to 1 for all atoms
            READ(4,*) NDUMMY,IRES(J1),RES(J1),TYPE(J1),DUMMY,DUMMY,DUMMY
            SEGID(J1)="1"
         ENDIF
      ENDDO
      CLOSE(4)
      IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.REGROUPPERSISTT.OR.DUMMYTST.OR.NOPOINTS)) &
  &            OPEN(UNIT=3,FILE='stationary.points.pdb',STATUS='UNKNOWN')
!
!  Use our free format .crd file to get the necessary data to make a fixed format pdb file
!
      J5=VERYBESTEND
      DO J1=1,NSTEPS
         IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.REGROUPPERSISTT.OR.DUMMYTST.OR.NOPOINTS)) THEN
            READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            IF (J1.GT.1) THEN
               CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
   &                            RMAT,.FALSE.)
            ENDIF
            DO J2=1,NATOMS
!              WRITE(3,'(a4,2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
               IF (.NOT.SHANNONT) WRITE(3,'(a4,2x,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
  &               'ATOM',J2,TYPE(J2),RES(J2),IRES(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3), &
  &                1.00,0.00,SEGID(J2)
            ENDDO
            WRITE(3,'(A)') 'END'
            ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
            READ(UTS,REC=TSPATHID(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
   &                         RMAT,.FALSE.)
            PRINT *,'J1, B, DISTANCES=',J1,DISTANCES
            DO J2=1,NATOMS
!              WRITE(3,'(a4,2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
               WRITE(3,'(a4,2x,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
  &               'ATOM',J2,TYPE(J2),RES(J2),IRES(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3), &
  &                1.00,0.00,SEGID(J2)
            ENDDO
            WRITE(3,'(A)') 'END'
            ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
            IF (J1.EQ.NSTEPS) THEN
               READ(UMIN,REC=VERYBESTPARENT(J5)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                             RMAT,.FALSE.)
               PRINT *,'J1, C, DISTANCES=',J1,DISTANCES
               DO J2=1,NATOMS
                  WRITE(3,'(a4,2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
  &                  'ATOM',J2,TYPE(J2),RES(J2),IRES(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3), &
  &                   1.00,0.00,SEGID(J2)
               ENDDO
               WRITE(3,'(A)') 'END'
            ENDIF
         ENDIF
         J5=VERYBESTPARENT(J5)
      ENDDO
   ELSE IF (AMBERT) THEN
! get the atom labels
       PRMLINES=(NATOMS-MOD(NATOMS,20))/20
       WRITE(PRMFORMAT,'(A,I2,A)') "(",MAX(20,MOD(NATOMS,20)),"A4)"
       OPEN(UNIT=4,FILE='coords.prmtop',STATUS='OLD')
       DO
        READ(4,'(A)',IOSTAT=ISTAT) DUMMYLINE
        IF(TRIM(ADJUSTL(DUMMYLINE))=='%FLAG ATOM_NAME') THEN
           READ(4,'(A)') DUMMYLINE
           EXIT
        END IF
        IF(ISTAT<0) EXIT
       END DO
       DO J1=1,PRMLINES
        READ(4,'(20A4)') ATOMLABELS(20*(J1-1)+1:20*J1)
       END DO
       IF(MOD(NATOMS,20)/=0) THEN
       WRITE(PRMFORMAT,'(A,I2,A)') "(",MOD(NATOMS,20),"A4)"
       READ(4,TRIM(ADJUSTL(PRMFORMAT))) ATOMLABELS(PRMLINES*20+1:PRMLINES*20+MOD(NATOMS,20))
        write(*,*) ATOMLABELS(PRMLINES+1:PRMLINES+MOD(NATOMS,20))
        WRITE(*,*) 'prmformat=',prmformat
!       READ(4,TRIM(ADJUSTL(PRMFORMAT))) ATOMLABELS(PRMLINES+1:PRMLINES+MOD(NATOMS,20))
       END IF
       WRITE(*,*) 'ATOMLABELS='
        WRITE(*,*) ATOMLABELS(1:NATOMS)
!! get the residue pointers (to be added later, now we just dump xyz files)
!       REWIND(4)
!       DO
!        READ(4,'(A)',STATUS=ISTAT) DUMMYLINE
!        IF(TRIM(ADJUSTL(DUMMYLINE))=='%FLAG RESIDUE_POINTER') THEN
!           READ(4,'(A)') DUMMYLINE
!           EXIT
!        END IF
!        IF(ISTAT<0) EXIT
!       END DO
!       DO J1=1,PRMLINES
!        READ(4,'(20A4)') ATOMLABELS(20*(J1-1)+1:20*J1)
!       END DO
!       IF(MOD(NATOMS,20)/=0) THEN
!        READ(4,TRIM(ADJUSTL(PRMFORMAT))) ATOMLABELS(PRMLINES+1:PRMLINES+MOD(NATOMS,20))
!       END IF
!       WRITE(*,*) 'ATOMLABELS=',ATOMLABELS(:)
!
       CLOSE(4)

      IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.REGROUPPERSISTT.OR.NOPOINTS)) &
  &                   OPEN(UNIT=3,FILE='stationary.points.xyz',STATUS='UNKNOWN')
      
      J5=VERYBESTEND
      DO J1=1,NSTEPS
         IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.REGROUPPERSISTT.OR.NOPOINTS)) THEN
            READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            IF (J1.GT.1) THEN
               CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                             RMAT,.FALSE.)
            ENDIF
            WRITE(3,'(I6)') NATOMS
            WRITE(3,*) 
            DO J2=1,NATOMS
               WRITE(3,'(a4,1x,3f8.3)') ATOMLABELS(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3)
            ENDDO
            ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
            READ(UTS,REC=TSPATHID(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                          RMAT,.FALSE.)
            WRITE(3,'(I6)') NATOMS
            WRITE(3,*) 
            DO J2=1,NATOMS
               WRITE(3,'(a4,1x,3f8.3)') ATOMLABELS(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3)
            ENDDO
            ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
            IF (J1.EQ.NSTEPS) THEN
               READ(UMIN,REC=VERYBESTPARENT(J5)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
   &                            RMAT,.FALSE.)
               WRITE(3,'(I6)') NATOMS
               WRITE(3,*) 
               DO J2=1,NATOMS
                 WRITE(3,'(a4,1x,3f8.3)') ATOMLABELS(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3)
               ENDDO
            ENDIF
         ENDIF
         J5=VERYBESTPARENT(J5)
      ENDDO
   ENDIF
   IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST.OR.REGROUPPERSISTT.OR.NOPOINTS)) CLOSE(3)

   IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST.OR.REGROUPPERSISTT.OR.NOPOINTS)) &
  &                   OPEN(UNIT=1,FILE='redopoints',STATUS='UNKNOWN')
   IF (ALLOCATED(BESTPATH)) DEALLOCATE(BESTPATH)
   ALLOCATE(BESTPATH(2*NSTEPS+1))
   J5=VERYBESTEND
   BESTPATHLENGTH=2*NSTEPS+1
   DO J1=1,NSTEPS
      IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST.OR.REGROUPPERSISTT.OR.NOPOINTS)) THEN
         READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         IF (J1.GT.1) THEN
          IF(.NOT.RBAAT)  CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
   &                         RMAT,.FALSE.)
         ENDIF
!
! We MUST use fixed format here in case the coordinates are read back into CHARMM. 
! CHARMM cannot read exponential format, even if the free keyword is specified!!
!
         WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
         READ(UTS,REC=TSPATHID(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
          IF(.NOT.RBAAT)  CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
   &                      RMAT,.FALSE.)
         WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
         IF (J1.EQ.NSTEPS) THEN ! last minimum
            READ(UMIN,REC=VERYBESTPARENT(J5)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
          IF(.NOT.RBAAT)  CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
   &                         RMAT,.FALSE.)
            WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDIF
      ENDIF
      IF (J1.EQ.NSTEPS) BESTPATH(2*J1+1)=VERYBESTPARENT(J5) ! last min
      BESTPATH(2*J1-1)=J5         ! min
      BESTPATH(2*J1)=TSPATHID(J1) ! ts
      J5=VERYBESTPARENT(J5)
   ENDDO
   IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST.OR.REGROUPPERSISTT.OR.NOPOINTS)) CLOSE(1)

   GOTO 555
!
! MFPT calculations for the best path IN ISOLATION from all the other stationary points.
! Uses the algorith of Weiss, Adv. Chem. Phys., 13, 1-18, 1967.
! To maintain precision we have to calculate ratios that are close to one and then subtract one.
!
! Now scrapped - GT approach is much better!
!
   ALLOCATE(ETA(NSTEPS+1), THETA(NSTEPS+1))

   ETA(1)=0.0D0
   THETA(1)=1.0D0
   ETAS=0.0D0
   THETAS=0.0D0
   ETASM1=0.0D0
   THETASM1=0.0D0
!
!  Direction opposite from current DIRECTION directive.
!
   DO J1=1,NSTEPS
      IF (PLUS(BESTPATH(2*J1)).EQ.BESTPATH(2*J1-1)) THEN
         ETA(J1+1)=(EXP(KPLUS(BESTPATH(2*J1)))*ETA(J1)+1.0D0)*EXP(-KMINUS(BESTPATH(2*J1)))
         THETA(J1+1)=THETA(J1)*EXP(KPLUS(BESTPATH(2*J1))-KMINUS(BESTPATH(2*J1)))
         IF (J1.EQ.NSTEPS) KNPLUS=EXP(KPLUS(BESTPATH(2*J1)))
!        PRINT '(A,G20.10)','case 1 ratio=', &
!  &                EXP(PFMIN(BESTPATH(2*J1-1))+KPLUS(BESTPATH(2*J1))-PFMIN(BESTPATH(2*J1+1))-KMINUS(BESTPATH(2*J1)))
      ELSEIF (MINUS(BESTPATH(2*J1)).EQ.BESTPATH(2*J1-1)) THEN
         ETA(J1+1)=(EXP(KMINUS(BESTPATH(2*J1)))*ETA(J1)+1.0D0)*EXP(-KPLUS(BESTPATH(2*J1)))
         THETA(J1+1)=THETA(J1)*EXP(KMINUS(BESTPATH(2*J1))-KPLUS(BESTPATH(2*J1)))
         IF (J1.EQ.NSTEPS) KNPLUS=EXP(KMINUS(BESTPATH(2*J1)))
!        PRINT '(A,G20.10)','case 2 ratio=', &
!  &                EXP(PFMIN(BESTPATH(2*J1-1))+KMINUS(BESTPATH(2*J1))-PFMIN(BESTPATH(2*J1+1))-KPLUS(BESTPATH(2*J1)))
      ELSE
         PRINT '(A)','Dijkstra> ERROR'
         STOP
      ENDIF
      ETAS=ETAS+ETA(J1)
      THETAS=THETAS+THETA(J1) 
      IF (J1.NE.NSTEPS) ETASM1=ETASM1+ETA(J1)
      IF (J1.NE.NSTEPS) THETASM1=THETASM1+THETA(J1) 
!     PRINT '(A,I6,4G20.10)','J1,ETA,THETA,ETAS,THETAS=',J1,ETA(J1),THETA(J1),ETAS,THETAS
   ENDDO
   IF (DIRECTION.EQ.'AB') THEN
      IF (NSTEPS.GT.1) THEN
         RAT1=(ETASM1/ETA(NSTEPS))
         RAT2=(THETASM1/THETA(NSTEPS))
         RAT3=(ETAS/ETA(NSTEPS))
         RAT4=(THETAS/THETA(NSTEPS))
         WAITBA=(RAT2-RAT1)/RAT3
         WAITBA=WAITBA*ETAS+RAT4/KNPLUS
         PRINT '(A,2G20.10)','Dijkstra> WAITBA new and old=',WAITBA,ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ELSE
         WAITBA=ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ENDIF
   ELSE
      IF (NSTEPS.GT.1) THEN
         RAT2=(THETASM1/THETA(NSTEPS))
         RAT3=(ETAS/ETA(NSTEPS))
         RAT4=(THETAS/THETA(NSTEPS))
         WAITAB=(RAT2-RAT1)/RAT3
         WAITAB=WAITAB*ETAS+RAT4/KNPLUS
         PRINT '(A,2G20.10)','Dijkstra> WAITAB new and old=',WAITAB,ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ELSE
         WAITAB=ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ENDIF
   ENDIF
!
!  Direction same as current DIRECTION directive.
!
   ETAS=0.0D0
   THETAS=0.0D0
   ETASM1=0.0D0
   THETASM1=0.0D0
   DO J1=1,NSTEPS
      J2=NSTEPS+1-J1
      IF (PLUS(BESTPATH(2*J2)).EQ.BESTPATH(2*J2+1)) THEN
         ETA(J1+1)=(EXP(KPLUS(BESTPATH(2*J2)))*ETA(J1)+1.0D0)*EXP(-KMINUS(BESTPATH(2*J2)))
         THETA(J1+1)=THETA(J1)*EXP(KPLUS(BESTPATH(2*J2))-KMINUS(BESTPATH(2*J2)))
         IF (J1.EQ.NSTEPS) KNPLUS=EXP(KPLUS(BESTPATH(2*J2)))
!        PRINT '(A,G20.10)','case 3 ratio=', &
!  &                EXP(PFMIN(BESTPATH(2*J2+1))+KPLUS(BESTPATH(2*J2))-PFMIN(BESTPATH(2*J2-1))-KMINUS(BESTPATH(2*J2)))
      ELSEIF (MINUS(BESTPATH(2*J2)).EQ.BESTPATH(2*J2+1)) THEN
         ETA(J1+1)=(EXP(KMINUS(BESTPATH(2*J2)))*ETA(J1)+1.0D0)*EXP(-KPLUS(BESTPATH(2*J2)))
         THETA(J1+1)=THETA(J1)*EXP(KMINUS(BESTPATH(2*J2))-KPLUS(BESTPATH(2*J2)))
         IF (J1.EQ.NSTEPS) KNPLUS=EXP(KMINUS(BESTPATH(2*J2)))
!        PRINT '(A,G20.10)','case 4 ratio=', &
!  &                EXP(PFMIN(BESTPATH(2*J2+1))+KMINUS(BESTPATH(2*J2))-PFMIN(BESTPATH(2*J2-1))-KPLUS(BESTPATH(2*J2)))
      ELSE
         PRINT '(A)','Dijkstra> ERROR'
         STOP
      ENDIF
      ETAS=ETAS+ETA(J1)
      THETAS=THETAS+THETA(J1) 
      IF (J1.NE.NSTEPS) ETASM1=ETASM1+ETA(J1)
      IF (J1.NE.NSTEPS) THETASM1=THETASM1+THETA(J1) 
!     PRINT '(A,I6,4G20.10)','J1,ETA,THETA,ETAS,THETAS=',J1,ETA(J1),THETA(J1),ETAS,THETAS
   ENDDO
!
!  Division by zero for zero MFPT's should now be avoided.
!
   IF (DIRECTION.EQ.'AB') THEN
      IF (NSTEPS.GT.1) THEN
         RAT1=(ETASM1/ETA(NSTEPS))
         RAT2=(THETASM1/THETA(NSTEPS))
         RAT3=(ETAS/ETA(NSTEPS))
         RAT4=(THETAS/THETA(NSTEPS))
         WAITAB=(RAT2-RAT1)/RAT3
         WAITAB=WAITAB*ETAS+RAT4/KNPLUS
         PRINT '(A,2G20.10)','Dijkstra> WAITAB new and old=',WAITAB,ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ELSE
         WAITAB=ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ENDIF
      IF (WAITAB.LE.0.0D0) THEN
         PRINT '(A)','Dijkstra> WARNING - Precision lost for WAITAB'
         WAITAB=-1.0D0
      ENDIF
      IF (WAITBA.LE.0.0D0) THEN
         PRINT '(A,G20.10)','Dijkstra> WARNING - Precision lost for WAITBA, RAT2-RAT1=',RAT2-RAT1
         WAITBA=-1.0D0
      ENDIF
      IF (.NOT.SHANNONT) THEN
         WRITE(*,'(A,G20.10)') 'Dijkstra> detailed balance, ratio should be one if SS applies: ', &
   &                         WAITAB*EXP(PFMIN(BESTPATH(1))-PFMIN(BESTPATH(2*NSTEPS+1)))/WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> MFPT for fastest path: A<-B value=',WAITAB,' B<-A value=',WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> rates without conditional probability: A<-B value=',1.0D0/WAITAB, &
   &                           ' B<-A value=',1.0D0/WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> rates with    conditional probability: A<-B value=', &
   &                        EXP(PFMIN(BESTPATH(2*NSTEPS+1))-PFTOTALB)/WAITAB, &
   &         ' B<-A value=',EXP(PFMIN(BESTPATH(1))-PFTOTALA)/WAITBA
      ENDIF
   ELSE
      IF (NSTEPS.GT.1) THEN
         RAT1=(ETASM1/ETA(NSTEPS))
         RAT2=(THETASM1/THETA(NSTEPS))
         RAT3=(ETAS/ETA(NSTEPS))
         RAT4=(THETAS/THETA(NSTEPS))
         WAITBA=(RAT2-RAT1)/RAT3
         WAITBA=WAITBA*ETAS+RAT4/KNPLUS
         IF (.NOT.SHANNONT) PRINT '(A,2G20.10)','Dijkstra> WAITBA new and old=',WAITBA,ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ELSE
         WAITBA=ETA(NSTEPS+1)*(THETAS/THETA(NSTEPS+1))-ETAS
      ENDIF
      IF (WAITAB.LE.0.0D0) THEN
         PRINT '(A)','Dijkstra> WARNING - Precision lost for WAITAB'
         WAITAB=-1.0D0
      ENDIF
      IF (WAITBA.LE.0.0D0) THEN
         PRINT '(A,G20.10)','Dijkstra> WARNING - Precision lost for WAITBA, RAT2-RAT1=',RAT2-RAT1
         WAITBA=-1.0D0
      ENDIF
      IF (.NOT.SHANNONT) WRITE(*,'(A,G20.10)') 'Dijkstra> detailed balance, ratio should be one if SS applies: ', &
   &                         WAITBA*EXP(PFMIN(BESTPATH(1))-PFMIN(BESTPATH(2*NSTEPS+1)))/WAITAB
      IF (.NOT.SHANNONT) PRINT '(2(A,G20.10))','Dijkstra> MFPT for fastest path: A<-B value=',WAITAB,' B<-A value=',WAITBA
!     PRINT '(A,2G20.10)','WAITAB,WAITBA=',WAITAB,WAITBA
      IF (.NOT.SHANNONT) PRINT '(2(A,G20.10))','Dijkstra> rates without conditional probability: A<-B value=',1.0D0/WAITAB, &
   &                           ' B<-A value=',1.0D0/WAITBA
      IF (.NOT.SHANNONT) PRINT '(2(A,G20.10))','Dijkstra> rates with    conditional probability: A<-B value=', &
   &                        EXP(PFMIN(BESTPATH(1))-PFTOTALB)/WAITAB, &
   &         ' B<-A value=',EXP(PFMIN(BESTPATH(2*NSTEPS+1))-PFTOTALA)/WAITBA
   ENDIF

   DEALLOCATE(ETA,THETA)
555 CONTINUE

!
!  Dump potential energy stationary points corresponding to fastest path. Must allow for
!  regrouping and sorting of free energy groups using the MINGROUP and MINMAP arrays.
!  MINGROUP(J1) is the original free energy group to which potential energy minimum J1
!  belongs. MINMAP(J2) tells us which original free energy group the sorted free energy
!  group J2 corresponds to. This is enough information!
!

   IF ((.NOT.NOPOINTS).AND.(REGROUPFREET.OR.REGROUPFREEABT.OR.REGROUPRATET.OR.REGROUPPERSISTT.OR.REGROUPPET)) THEN ! deal with free energy groups
      OPEN(UNIT=2,FILE='min.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='points.min.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      ALLOCATE(INCLUDEMIN(NMINSAVE),NEWINDEX(NMINSAVE))
      INCLUDEMIN(1:NMINSAVE)=.FALSE.
      NEWINDEX(1:NMINSAVE)=0
      SAVEGROUP(1:NMIN)=.FALSE.
      DO J1=1,NSTEPS+1 ! these are the free energy minima
         SAVEGROUP(MINMAP(BESTPATH(2*J1-1)))=.TRUE.
      ENDDO
      DO J2=1,NMINSAVE
         IF (MINGROUP(J2).LE.0) CYCLE ! PE minima with insufficient connections etc. can have MINGROUP zero.
         IF (SAVEGROUP(MINGROUP(J2))) THEN
            IF (.NOT.INCLUDEMIN(J2)) THEN
               INCLUDEMIN(J2)=.TRUE.
            ENDIF
            IF (DEBUG) PRINT '(4(A,I8))','Dijkstra> Including PE minimum ',J2,' from group ',MINGROUP(J2), &
  &                          ' originally ',MINMAP(MINGROUP(J2))
         ENDIF
      ENDDO
!     OPEN(UNIT=7,FILE='min.data',STATUS='OLD')
      IF (CLOSEFILEST) THEN
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN')
      ELSE
         REWIND(UMINDATA) ! is already open in setup
      ENDIF
      NDUMMY=0
      DO J2=1,NMINSAVE
         READ(UMINDATA,'(A)') DUMMYSTRING
         IF (.NOT.INCLUDEMIN(J2)) CYCLE
         NDUMMY=NDUMMY+1
         NEWINDEX(J2)=NDUMMY
         WRITE(2,'(A,2X,I8)') TRIM(ADJUSTL(DUMMYSTRING)),J2
         READ(UMIN,REC=J2) (LOCALPOINTS(J3),J3=1,3*NATOMS)
         WRITE(4,REC=NDUMMY) (LOCALPOINTS(J3),J3=1,3*NATOMS)
      ENDDO
      CLOSE(2); CLOSE(4)  ! ; CLOSE(7)
      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
!     OPEN(UNIT=7,FILE='ts.data',STATUS='OLD')
      IF (CLOSEFILEST) THEN
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD')
      ELSE
         REWIND(UTSDATA) ! is already open in setup
      ENDIF
      OPEN(UNIT=3,FILE='ts.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='points.ts.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      NDUMMY=0
      NDUMMY2=0
      DO 
         NDUMMY=NDUMMY+1
         READ(UTSDATA,*,END=40) DETS,DFVIBTS,DHORDERTS,DPLUS,DMINUS,DIXTS,DIYTS,DIZTS
         IF (NOFRQS) DFVIBTS=1.0D0 ! for consistency with setup etc.
         IF (INCLUDEMIN(DPLUS).AND.INCLUDEMIN(DMINUS)) THEN
            IF (DEBUG) PRINT '(3(A,I8))','Dijkstra> Including PE ts ',NDUMMY
            NDUMMY2=NDUMMY2+1
            WRITE(3,'(2F20.10,3I8,3F20.10,2X,I8)')  &
  &                DETS,DFVIBTS,DHORDERTS,NEWINDEX(DPLUS),NEWINDEX(DMINUS),DIXTS,DIYTS,DIZTS,NDUMMY
            IF (.NOT.(DUMMYTST.OR.NOPOINTS)) THEN
               READ(UTS,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(5,REC=NDUMMY2) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDIF
         ENDIF
      ENDDO
40    CONTINUE
      CLOSE(3); CLOSE(5) ! ; CLOSE(7)
      IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
      OPEN(UNIT=2,FILE='min.A.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=3,FILE='min.A',STATUS='OLD')
      READ(3,*) NNMINA
      NDUMMY2=0
      DO J1=1,NNMINA
         READ(3,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) NDUMMY2=NDUMMY2+1
      ENDDO
      WRITE(2,*) NDUMMY2
      REWIND(3)
      READ(3,*) NNMINA
      DO J1=1,NNMINA
         READ(3,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) WRITE(2,*) NEWINDEX(NDUMMY)
      ENDDO
      CLOSE(2); CLOSE(3)

      OPEN(UNIT=4,FILE='min.B.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='min.B',STATUS='OLD')
      READ(5,*) NNMINB
      NDUMMY2=0
      DO J1=1,NNMINB
         READ(5,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) NDUMMY2=NDUMMY2+1
      ENDDO
      WRITE(4,*) NDUMMY2
      REWIND(5)
      READ(5,*) NNMINB
      DO J1=1,NNMINB
         READ(5,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) WRITE(4,*) NEWINDEX(NDUMMY)
      ENDDO
      CLOSE(4); CLOSE(5)
      DEALLOCATE(INCLUDEMIN,NEWINDEX)
   ELSEIF (.NOT.NOPOINTS) THEN ! the PE minima were not regrouped
      OPEN(UNIT=2,FILE='min.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='points.min.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      DO J1=1,NSTEPS+1
         J5=BESTPATH(2*J1-1)
         WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J5), FVIBMIN(J5), HORDERMIN(J5), IXMIN(J5), IYMIN(J5), IZMIN(J5)
         READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         WRITE(4,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      ENDDO
      CLOSE(2); CLOSE(4)
      OPEN(UNIT=3,FILE='ts.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='points.ts.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      DO J1=1,NSTEPS
         J5=TSPATHID(J1)
         IF (IMFRQT) THEN
            WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J5),FVIBTS(J5),HORDERTS(J5),J1,J1+1,IXTS(J5),IYTS(J5),IZTS(J5),NEGEIG(J5)
         ELSE
            WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J5),FVIBTS(J5),HORDERTS(J5),J1,J1+1,IXTS(J5),IYTS(J5),IZTS(J5)
         ENDIF
         IF (.NOT.(DUMMYTST.OR.NOPOINTS)) THEN
            READ(UTS,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(5,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDIF
      ENDDO
      CLOSE(3); CLOSE(5)
      OPEN(UNIT=2,FILE='min.A.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='min.B.fastest',STATUS='UNKNOWN')
      WRITE(2,'(I6)') 1
      WRITE(4,'(I6)') 1
      IF (DIRECTION.EQ.'AB') THEN
!
! jmc swapped units - check!!!
!
         WRITE(2,'(I6)') 1
         WRITE(4,'(I6)') NSTEPS+1
      ELSE
         WRITE(4,'(I6)') 1
         WRITE(2,'(I6)') NSTEPS+1
      ENDIF
      CLOSE(2); CLOSE(4)
   ENDIF
!
! call MFPT subroutine instead of the Weiss MFPT calculation skipped using label 555.
!
   CALL MFPT(NSTEPS, TSPATHID, WAITAB, WAITBA)

   IF (DIRECTION.EQ.'AB') THEN
      IF (.NOT.SHANNONT) THEN
         WRITE(*,'(A,G20.10)') 'Dijkstra> detailed balance, ratio should be one if SS applies: ', &
   &                         WAITAB*EXP(PFMIN(BESTPATH(1))-PFMIN(BESTPATH(2*NSTEPS+1)))/WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> MFPT for fastest path: A<-B value=',WAITAB,' B<-A value=',WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> rates without conditional probability: A<-B value=',1.0D0/WAITAB, &
   &                           ' B<-A value=',1.0D0/WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> rates with    conditional probability: A<-B value=', &
   &                        EXP(PFMIN(BESTPATH(2*NSTEPS+1))-PFTOTALB)/WAITAB, &
   &         ' B<-A value=',EXP(PFMIN(BESTPATH(1))-PFTOTALA)/WAITBA
      ENDIF
   ELSE
      IF (.NOT.SHANNONT) THEN
         WRITE(*,'(A,G20.10)') 'Dijkstra> detailed balance, ratio should be one if SS applies: ', &
   &                         WAITBA*EXP(PFMIN(BESTPATH(1))-PFMIN(BESTPATH(2*NSTEPS+1)))/WAITAB
         PRINT '(2(A,G20.10))','Dijkstra> MFPT for fastest path: A<-B value=',WAITAB,' B<-A value=',WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> rates without conditional probability: A<-B value=',1.0D0/WAITAB, &
   &                           ' B<-A value=',1.0D0/WAITBA
         PRINT '(2(A,G20.10))','Dijkstra> rates with    conditional probability: A<-B value=', &
   &                        EXP(PFMIN(BESTPATH(1))-PFTOTALB)/WAITAB, &
   &         ' B<-A value=',EXP(PFMIN(BESTPATH(2*NSTEPS+1))-PFTOTALA)/WAITBA
      ENDIF
   ENDIF

ENDIF
!
!  Summarise pairs of minima for the worst transition states, i.e. the highest ones on
!  each path, if we did this calculation.
!
NWORST=MIN(NWORST,PAIRSTODO)
IF (NWORST.GT.0) THEN
   IF (.NOT.SHANNONT) PRINT '(A)','Dijkstra> Highest transition states on the best A<->B paths and neighbouring minima:'
   IF (.NOT.SHANNONT) PRINT '(A)','    min1    min2      ts         weight     attempts'
ENDIF
DO J1=1,NWORST
   IF (.NOT.SHANNONT) PRINT '(3I8,G20.10,I8)',DMIN1(J1),DMIN2(J1),WORSTTS(J1),TSWEIGHT(J1),TSATTEMPT(WORSTTS(J1))
ENDDO

IF (ALLOCATED(MAXDOWNHILLB)) DEALLOCATE(MAXDOWNHILLB)
IF (ALLOCATED(MAXDOWNHILLTS)) DEALLOCATE(MAXDOWNHILLTS)
IF (ALLOCATED(SHIFTABLE)) DEALLOCATE(SHIFTABLE)
CALL CPU_TIME(TNEW)
TDIJKSTRA=TDIJKSTRA+TNEW-ELAPSED
IF (SHANNONT) NWORST=NSTEPS ! return the number of steps in this case (number of ts).

RETURN
END

!
!     This subprogram performs a sort on the input data and
!     arranges it from smallest to biggest. The exchange-sort
!     algorithm is used.
!
      SUBROUTINE SORT3(NWORST,EREF,DMIN1,DMIN2)
      IMPLICIT NONE
      INTEGER J1, L, NWORST, J2, DMIN1(*), DMIN2(*), NTEMP
      DOUBLE PRECISION EREF(*), TEMP
!
      DO 20 J1=1,NWORST-1
         L=J1
         DO 10 J2=J1+1,NWORST
            IF (EREF(L).GT.EREF(J2)) L=J2
10       CONTINUE
         TEMP=EREF(L)
         EREF(L)=EREF(J1)
         EREF(J1)=TEMP
         NTEMP=DMIN1(L)
         DMIN1(L)=DMIN1(J1)
         DMIN1(J1)=NTEMP
         NTEMP=DMIN2(L)
         DMIN2(L)=DMIN2(J1)
         DMIN2(J1)=NTEMP
20    CONTINUE
      RETURN
      END

! MFPT from a GT calculation for the fastest SS path. 
! Hopefully won`t suffer from the same numerical problems!
! All intervening minima have two connections, and the end
! minima only one each.
! The minima and transition states in question are saved in integer vector 
! BESTPATH from BESTPATH(1) to BESTPATH(2*NSTEPS+1).
! Each minimum needs a plus/minus branching probability and a waiting time,
! both based on the rate constants, and subsequently renormalised. 
! If we remove minima from 2*NSTEPS-1 to 3 then only the waiting time of
! minimum NSTEPS+1 changes - its branching probability is unity to the
! next minimum in the chain.
!
SUBROUTINE MFPT(NSTEPS, TSPATHID, WAITAB, WAITBA)
USE PORFUNCS
USE COMMONS
IMPLICIT NONE

INTEGER, PARAMETER :: MAXPATH=10000 ! longest path length allowed - surely 10000 is enough?!
INTEGER :: J1, NSTEPS, TSPATHID(MAXPATH)
DOUBLE PRECISION :: WAITAB, WAITBA, FACTOR
DOUBLE PRECISION, ALLOCATABLE :: TAU(:), BRANCHP(:), BRANCHM(:)

   ALLOCATE(TAU(NSTEPS+1), BRANCHP(NSTEPS+1), BRANCHM(NSTEPS+1))
   TAU(1:NSTEPS+1)=0.0D0
   BRANCHM(1:NSTEPS+1)=0.0D0
   BRANCHP(1:NSTEPS+1)=0.0D0
!
! Initialise waiting times and +/- branching probabilities.
!
   DO J1=1,NSTEPS
!     PRINT *,'J1,TSPATHID,kplus,kminus,BESTPATH-,BESTPATH+: ',J1,TSPATHID(J1),PLUS(TSPATHID(J1)),MINUS(TSPATHID(J1)), &
! &              BESTPATH(2*J1-1),BESTPATH(2*J1+1)
      IF (PLUS(TSPATHID(J1)).EQ.BESTPATH(2*J1+1)) THEN
         IF (MINUS(TSPATHID(J1)).NE.BESTPATH(2*J1-1)) THEN
            PRINT *,'ERROR ts ',TSPATHID(J1),' plus and minus minima are ',PLUS(TSPATHID(J1)),MINUS(TSPATHID(J1))
            PRINT *,'ERROR BESTPATH minima are: ',BESTPATH(2*J1-1),BESTPATH(2*J1+1)
         ENDIF
         BRANCHM(J1+1)=EXP(KPLUS(TSPATHID(J1)))
         BRANCHP(J1)  =EXP(KMINUS(TSPATHID(J1)))
      ELSE
         IF (MINUS(TSPATHID(J1)).NE.BESTPATH(2*J1+1)) THEN
            PRINT *,'ERROR ts ',TSPATHID(J1),' plus and minus minima are ',PLUS(TSPATHID(J1)),MINUS(TSPATHID(J1))
            PRINT *,'ERROR BESTPATH minima are: ',BESTPATH(2*J1-1),BESTPATH(2*J1+1)
         ENDIF
         BRANCHM(J1+1)=EXP(KMINUS(TSPATHID(J1)))
         BRANCHP(J1)  =EXP(KPLUS(TSPATHID(J1)))
      ENDIF
   ENDDO
   DO J1=1,NSTEPS+1
!     PRINT *,'J1,BRANCHM,BRANCHP=',J1,BRANCHM(J1),BRANCHP(J1)
      TAU(J1)=1.0D0/(BRANCHM(J1)+BRANCHP(J1))
      BRANCHM(J1)=BRANCHM(J1)*TAU(J1)
      BRANCHP(J1)=BRANCHP(J1)*TAU(J1)
   ENDDO
!
!  Renormalise NSTEPS-1 times for the NSTEPS-1 intervening minima.
!
   DO J1=NSTEPS,2,-1
!
! Change waiting times for adjacent minima J1-1 and NSTEPS+1
! Change branching probabilities for adjacent minimum J1-1 
! (BRANCHM(NSTEPS+1)=1 throughout).
!
      FACTOR=BRANCHM(J1-1)+BRANCHP(J1)-BRANCHM(J1-1)*BRANCHP(J1)
      FACTOR=1.0D0/FACTOR
      TAU(NSTEPS+1)=(TAU(NSTEPS+1)+TAU(J1))/BRANCHM(J1)
      TAU(J1-1)=(TAU(J1-1)+TAU(J1)*BRANCHP(J1-1))*FACTOR
      BRANCHP(J1-1)=BRANCHP(J1-1)*BRANCHP(J1)*FACTOR
      BRANCHM(J1-1)=BRANCHM(J1-1)*FACTOR
      IF (DEBUG) PRINT *,'MFPT> J1-1,tau,b+,b-,sum=',J1-1,TAU(J1-1),BRANCHP(J1-1),BRANCHM(J1-1),BRANCHP(J1-1)+BRANCHM(J1-1)
!     PRINT *,'MFPT> J1-1,    tau,b+,b-,sum=',J1-1,TAU(J1-1),BRANCHP(J1-1),BRANCHM(J1-1),BRANCHP(J1-1)+BRANCHM(J1-1)
!     PRINT *,'MFPT> NSTEPS+1,tau,b+,b-,sum=',NSTEPS+1,TAU(NSTEPS+1),BRANCHP(NSTEPS+1),BRANCHM(NSTEPS+1), &
! &                       BRANCHP(NSTEPS+1)+BRANCHM(NSTEPS+1)
!     DO J2=1,J1-1
!        PRINT *,'after removing minimum ',J1,' minimum,TAU,p+,p-=',J2,TAU(J2),BRANCHP(J2),BRANCHM(J2)
!     ENDDO
!     PRINT *,'after removing minimum ',J1,' minimum,TAU,p+,p-=',NSTEPS+1,TAU(NSTEPS+1),BRANCHP(NSTEPS+1),BRANCHM(NSTEPS+1)
   ENDDO

   IF (DIRECTION.EQ.'AB') THEN
      WAITAB=TAU(NSTEPS+1)
      WAITBA=TAU(1)
   ELSE
      WAITAB=TAU(1)
      WAITBA=TAU(NSTEPS+1)
   ENDIF

   DEALLOCATE(TAU,BRANCHP,BRANCHM)

   RETURN

END SUBROUTINE MFPT
