!
!  Persistence analysis using exact TS energy, rather 
!  than an estimate from a superbasin analysis.
!  Persistence diagram is a plot of local minimum energy
!  versus transition state energy where that minimum joins a
!  superbasin with a lower minimum in it.
!
!  uses a sorted list of ts energies
!
SUBROUTINE GETEXACTBARRIER(BARRIER,MINTARG,GMIN)
USE COMMONS,ONLY : NMIN, NTS, ETS, EMIN, PLUS, MINUS, KPLUS, KMINUS, NCONNMIN, NCONN, EUNTRAPTHRESH, ALLCOMPONENTST
IMPLICIT NONE
DOUBLE PRECISION LOWESTTARG, BARRIER(NMIN), DUMMY, GMIN, TEMP, NEWETS(NTS), E1, E2
DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0 ! set here but NOT used in the current coding of checkTS
INTEGER NEWPLUS(NTS), NEWMINUS(NTS), L, J2, NTEMP
INTEGER J1, MINTARG, NDEAD, NDISTA(NMIN), NDISTB(NMIN), NCYCLE
INTEGER NUNCONA, NLEFTMIN, NLEFTTS, COMPONENT(NMIN)
LOGICAL CHANGED, DEADTS(NTS), PERSIST(NMIN), LTEMP, NEWDEADTS(NTS)

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

LOWESTTARG=EMIN(MINTARG) ! energy of global minimum, which has no culminance

BARRIER(1:NMIN)=-1.0D0

DO J1=1,NMIN
   COMPONENT(J1)=J1 ! which connected component J1 belongs to. Initially all min
                    ! are in their own component; later the label for the
                    ! component is determined by the lowest minimum in that set.
   PERSIST(J1)=.FALSE.
ENDDO

! sort TSs lowest to highest
NEWPLUS=PLUS
NEWMINUS=MINUS
NEWETS=ETS
NEWDEADTS=DEADTS
DO J1=1,NTS-1
   L=J1
   DO J2=J1+1,NTS
      IF (NEWETS(L).GT.NEWETS(J2)) L=J2
   ENDDO
   TEMP=NEWETS(L)
   NEWETS(L)=NEWETS(J1)
   NEWETS(J1)=TEMP
   NTEMP=NEWPLUS(L)
   NEWPLUS(L)=NEWPLUS(J1)
   NEWPLUS(J1)=NTEMP
   NTEMP=NEWMINUS(L)
   NEWMINUS(L)=NEWMINUS(J1)
   NEWMINUS(J1)=NTEMP
   LTEMP=NEWDEADTS(L)
   NEWDEADTS(L)=NEWDEADTS(J1)
   NEWDEADTS(J1)=LTEMP
ENDDO

DO J1=1,NTS
   IF (NEWDEADTS(J1)) CYCLE ! ignore this ts
! do nothing if the two minima connected by the TS are already in the same component
   IF (COMPONENT(NEWPLUS(J1)).EQ.COMPONENT(NEWMINUS(J1))) CYCLE
   IF ((EMIN(NEWPLUS(J1))-LOWESTTARG).GT.EUNTRAPTHRESH) CYCLE
   IF ((EMIN(NEWMINUS(J1))-LOWESTTARG).GT.EUNTRAPTHRESH) CYCLE
   IF (EMIN(NEWPLUS(J1)).GT.EMIN(NEWMINUS(J1))) THEN
      IF (.NOT.PERSIST(NEWPLUS(J1))) THEN
         BARRIER(NEWPLUS(J1))=NEWETS(J1)-EMIN(NEWPLUS(J1))
         PERSIST(NEWPLUS(J1))=.TRUE.
! set the component of all min already connected to NEWPLUS(J1) (and NEWPLUS(J1) itself) 
! to that of the lower-energy minimum
         WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
      ELSE
! we dont reset the persistence of NEWPLUS(J1), but because this TS now links two
! separate components we set the barrier for the higher of the two min that are
! lowest in their respective components (NB this barrier should not previously have
! been set). We also put all minima in this new component into the component of the 
! lowest member of this set.
         E1=EMIN(COMPONENT(NEWPLUS(J1)))
         E2=EMIN(COMPONENT(NEWMINUS(J1)))
         IF (E1.GT.E2) THEN
            BARRIER(COMPONENT(NEWPLUS(J1)))=NEWETS(J1)-EMIN(COMPONENT(NEWPLUS(J1)))
            PERSIST(COMPONENT(NEWPLUS(J1)))=.TRUE.
            WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
         ELSEIF (E2.GT.E1) THEN
            BARRIER(COMPONENT(NEWMINUS(J1)))=NEWETS(J1)-EMIN(COMPONENT(NEWMINUS(J1)))
            PERSIST(COMPONENT(NEWMINUS(J1)))=.TRUE.
            WHERE (COMPONENT.EQ.COMPONENT(NEWMINUS(J1))) COMPONENT=COMPONENT(NEWPLUS(J1))
         ELSE
            PRINT *,'degeneracy here: ',J1, E1, E2
            WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
         ENDIF
      ENDIF
! and if the two connected minima are the other way round in energy...
   ELSEIF (EMIN(NEWMINUS(J1)).GT.EMIN(NEWPLUS(J1))) THEN
      IF (.NOT.PERSIST(NEWMINUS(J1))) THEN
         BARRIER(NEWMINUS(J1))=NEWETS(J1)-EMIN(NEWMINUS(J1))
         PERSIST(NEWMINUS(J1))=.TRUE.
         WHERE (COMPONENT.EQ.COMPONENT(NEWMINUS(J1))) COMPONENT=COMPONENT(NEWPLUS(J1))
      ELSE
         E1=EMIN(COMPONENT(NEWPLUS(J1)))
         E2=EMIN(COMPONENT(NEWMINUS(J1)))
         IF (E1.GT.E2) THEN
            BARRIER(COMPONENT(NEWPLUS(J1)))=NEWETS(J1)-EMIN(COMPONENT(NEWPLUS(J1)))
            PERSIST(COMPONENT(NEWPLUS(J1)))=.TRUE.
            WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
         ELSEIF (E2.GT.E1) THEN
            BARRIER(COMPONENT(NEWMINUS(J1)))=NEWETS(J1)-EMIN(COMPONENT(NEWMINUS(J1)))
            PERSIST(COMPONENT(NEWMINUS(J1)))=.TRUE.
            WHERE (COMPONENT.EQ.COMPONENT(NEWMINUS(J1))) COMPONENT=COMPONENT(NEWPLUS(J1))
         ELSE
            PRINT *,'degeneracy here ',J1, E1, E2
            WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
         ENDIF
      ENDIF
! if the two connected minima are degenerate, then dont set their own persistence but see 
! what they are connected to and at least update their connectivity.
   ELSE
      print *,'degeneracy here ',J1,NEWETS(J1),EMIN(NEWMINUS(J1)),EMIN(NEWPLUS(J1))
      E1=EMIN(COMPONENT(NEWPLUS(J1)))
      E2=EMIN(COMPONENT(NEWMINUS(J1)))
      IF (E1.GT.E2) THEN
         BARRIER(COMPONENT(NEWPLUS(J1)))=NEWETS(J1)-EMIN(COMPONENT(NEWPLUS(J1)))
         PERSIST(COMPONENT(NEWPLUS(J1)))=.TRUE.
         WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
      ELSEIF (E2.GT.E1) THEN
         BARRIER(COMPONENT(NEWMINUS(J1)))=NEWETS(J1)-EMIN(COMPONENT(NEWMINUS(J1)))
         PERSIST(COMPONENT(NEWMINUS(J1)))=.TRUE.
         WHERE (COMPONENT.EQ.COMPONENT(NEWMINUS(J1))) COMPONENT=COMPONENT(NEWPLUS(J1))
      ELSE
         PRINT *,'degeneracy here ',J1, E1, E2
         WHERE (COMPONENT.EQ.COMPONENT(NEWPLUS(J1))) COMPONENT=COMPONENT(NEWMINUS(J1))
      ENDIF
   ENDIF
ENDDO

IF (.NOT.ALLCOMPONENTST) THEN
   DO J1=1,NMIN
      IF (NCONN(J1).EQ.0) BARRIER(J1)=-1.0D0
! as NCONN is set to zero earlier in this routine if this minimum is not connected to the global min.
! NB in persistence.f90, a min is only processed if it has a positive barrier.
   ENDDO
ENDIF

END SUBROUTINE GETEXACTBARRIER
