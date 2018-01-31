!
!  Persistence analysis using getbarrier2 routine.
!  Persistance diagram is a plot of local minimum energy
!  versus transition state energy where that minimum joins a
!  superbasin with a lower minimum in it. This ts energy is
!  called the culminance of the minimum in question.
!
!  Usual superbasin analysis at fine grained total energies.
!  Could use the sorted list of ts energies?
!  Check if minimum has a threshold assigned.
!  If not, check if it is the lowest minimum in the superbasin.
!  If not, assign that threshold as its culminance, and mark it done.
!
SUBROUTINE GETBARRIER2(BARRIER,MINTARG,GMIN)
USE COMMONS,ONLY : NMIN, NTS, ETS, EMIN, PLUS, MINUS, KPLUS, KMINUS, EINC, NCONNMIN, NCONN, EUNTRAPTHRESH, SHANNONT
IMPLICIT NONE
DOUBLE PRECISION HIGHESTTS, LOWESTTARG, ETHRESH, BARRIER(NMIN), BASINE(NMIN), DUMMY, GMIN
DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0 ! set here but NOT used in the current coding of checkTS
INTEGER J1, BASIN(NMIN), NBASIN, MINTARG, NDEAD, NDISTA(NMIN), NDISTB(NMIN), NCYCLE
INTEGER NUNCONA, NLEFTMIN, NLEFTTS, DJWBASIN(NMIN)
LOGICAL CHANGED, DEADTS(NTS)

CALL GETNCONN
DEADTS(1:NTS)=.FALSE.
IF (NCONNMIN.GE.0) THEN
   NDEAD=0
   DO J1=1,NMIN
      IF (NCONN(J1).LE.NCONNMIN) THEN
         NDEAD=NDEAD+1
      ENDIF 
   ENDDO
   PRINT '(3(I8,A))',NDEAD,' minima with ',NCONNMIN,' connections or fewer will not be considered'
ENDIF
!
! Determine position of global minimum
!
GMIN=1.0D100
DO J1=1,NMIN
   IF (EMIN(J1).LT.GMIN .AND. (NCONN(J1).GT.NCONNMIN)) THEN
      GMIN=EMIN(J1)
      MINTARG=J1
   ENDIF
ENDDO
!
!  Flag transition states to underconnected minima as DEAD.
!  NCONN only counts non-degenerate rearrangements as connections.
!
NLEFTTS=0
DO J1=1,NTS
   CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,DEADTS(J1))
   IF (.NOT.DEADTS(J1)) NLEFTTS=NLEFTTS+1
ENDDO
!
!  Check that the stationary point database is actually connected, and remove
!  minima that lie in disjoint graphs.
!  Calculate minimum number of steps of each minimum from the global minimum.
!
NDISTA(1:NMIN)=1000000
NDISTA(MINTARG)=0
NCYCLE=0
5  CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DO J1=1,NTS
   IF (DEADTS(J1)) CYCLE
   IF (NDISTA(MINUS(J1))+1.LT.NDISTA(PLUS(J1))) THEN
      CHANGED=.TRUE.
      NDISTA(PLUS(J1))=NDISTA(MINUS(J1))+1
   ENDIF
   IF (NDISTA(PLUS(J1))+1.LT.NDISTA(MINUS(J1))) THEN
      CHANGED=.TRUE.
      NDISTA(MINUS(J1))=NDISTA(PLUS(J1))+1
   ENDIF
ENDDO
IF (CHANGED) GOTO 5
NUNCONA=0
NLEFTMIN=0
DO J1=1,NMIN
   IF (NDISTA(J1).EQ.1000000) THEN
      NUNCONA=NUNCONA+1
      NCONN(J1)=0
      NDISTA(MINTARG)=0
   ELSEIF (NCONN(J1).GT.NCONNMIN) THEN
      NLEFTMIN=NLEFTMIN+1
   ENDIF
ENDDO
PRINT '(2(A,I8))','Steps to global minimum converged in ',NCYCLE-1,' cycles; disconnected=',NUNCONA
PRINT '(A,2I8)','Number of remaining minima and transition states=',NLEFTMIN,NLEFTTS

!  Assign the minima to their basins at each energy.

WRITE (6, '(A)') 'Assigning minima to basins.'
LOWESTTARG=EMIN(MINTARG) ! energy of global minimum, which has no culminance
HIGHESTTS=-1.0D100
DO J1=1,NTS
   IF ((ETS(J1).GT.HIGHESTTS).AND.(.NOT.DEADTS(J1))) THEN
      HIGHESTTS=ETS(J1)
   ENDIF
ENDDO
PRINT '(A,G20.10)','getbarrier2> highest transition state lies at ',HIGHESTTS

ETHRESH=LOWESTTARG ! start at global minimum.
!  BASIN=0
BARRIER(1:NMIN)=-1.0D0

DO 

  DJWBASIN(1:NMIN)=0
  BASINE(1:NMIN)=1.0D200 ! initialise energy of lowest minimum in superbasin
  NBASIN=0
  DO 
     CHANGED=.FALSE.
     DO J1=1,NTS
        IF (DEADTS(J1)) CYCLE
        IF (ETS(J1).LT.ETHRESH) THEN
           IF ((DJWBASIN(PLUS(J1)).EQ.0).AND.(DJWBASIN(MINUS(J1)).EQ.0)) THEN
              CHANGED=.TRUE.
              NBASIN=NBASIN+1
              DJWBASIN(PLUS(J1))=NBASIN
              DJWBASIN(MINUS(J1))=NBASIN
              BASINE(NBASIN)=MIN(EMIN(PLUS(J1)),EMIN(MINUS(J1)))
           ELSEIF (DJWBASIN(PLUS(J1)).NE.DJWBASIN(MINUS(J1))) THEN
              CHANGED=.TRUE.
              IF (DJWBASIN(PLUS(J1)).EQ.0) THEN
                 BASINE(DJWBASIN(MINUS(J1)))=MIN(EMIN(PLUS(J1)),BASINE(DJWBASIN(MINUS(J1))))
                 DJWBASIN(PLUS(J1))=DJWBASIN(MINUS(J1))
              ELSEIF (DJWBASIN(MINUS(J1)).EQ.0) THEN
                 BASINE(DJWBASIN(PLUS(J1)))=MIN(EMIN(MINUS(J1)),BASINE(DJWBASIN(PLUS(J1))))
                 DJWBASIN(MINUS(J1))=DJWBASIN(PLUS(J1))
              ELSE
                 DUMMY=MIN(BASINE(DJWBASIN(MINUS(J1))),BASINE(DJWBASIN(PLUS(J1))))
                 DJWBASIN(PLUS(J1))=MIN(DJWBASIN(PLUS(J1)),DJWBASIN(MINUS(J1)))
                 DJWBASIN(MINUS(J1))=DJWBASIN(PLUS(J1))
                 BASINE(DJWBASIN(MINUS(J1)))=DUMMY
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     IF (.NOT.CHANGED) EXIT
  ENDDO 
!
!  No need to remove empty groups and close up gaps. 
!  Disconnected minima should all belong to group 0.
!  Minima that are not connected at this energy level but lie below it
!  and will be connected should be assigned to a non-zero basin. DJW
!
  DO J1=1,NMIN
     IF ((DJWBASIN(J1).EQ.0).AND.(NCONN(J1).GT.NCONNMIN).AND.(EMIN(J1).LT.ETHRESH)) THEN
        NBASIN=NBASIN+1
        DJWBASIN(J1)=NBASIN
        BASINE(NBASIN)=EMIN(J1)
     ENDIF
  ENDDO 

  IF (NBASIN == 1) THEN
     IF (.NOT.SHANNONT) WRITE (6, '(I6, A, F18.10)') 1, ' basin at energy threshold ', ETHRESH
  ELSE
     IF (.NOT.SHANNONT) WRITE (6, '(I6, A, F18.10)') NBASIN, ' basins at energy threshold ', ETHRESH
  ENDIF
   DO J1=1,NMIN
      IF (EMIN(J1).GE.ETHRESH) CYCLE
      IF (NCONN(J1).LE.NCONNMIN) CYCLE
      IF (J1.EQ.MINTARG) CYCLE
!     PRINT '(A,I6,2G20.10)','J1,EMIN(J1),DJWBASIN(J1)=',J1,EMIN(J1),DJWBASIN(J1)
!     PRINT '(A,I6,G20.10,I6,G20.10)','J1,BARRIER,DJWBASIN,BASINE(DJWBASIN(J1))=',J1,BARRIER(J1),DJWBASIN(J1),BASINE(DJWBASIN(J1))
      IF (BARRIER(J1).GT.0.0D0) CYCLE
      IF (DJWBASIN(J1).EQ.0) CYCLE ! this should not be possible?
      IF (EMIN(J1).GT.BASINE(DJWBASIN(J1))) THEN
         BARRIER(J1)=ETHRESH-EMIN(J1) ! defined positive as barrier not ts energy
!        PRINT '(A,I6,F20.10)','j1,barrier=',J1,BARRIER(J1)
      ENDIF
   ENDDO
   ETHRESH=ETHRESH+EINC
   IF ((ETHRESH.GT.HIGHESTTS+EINC) .OR. (ETHRESH-LOWESTTARG.GT.EUNTRAPTHRESH)) EXIT 
ENDDO
WRITE (6, '(A, /)') 'Done.'

END SUBROUTINE GETBARRIER2
