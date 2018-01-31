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
!  Subroutine to provide candidate pairs of minima based on equilibrium occupation probability
!  times waiting time. We try to make a connection to minima on the fastest path for which
!  no connection exists in order of increasing distance. For a minimum that is itself on the
!  fastest path we should only try connections to minima that are fewer steps away from the
!  product.
!
SUBROUTINE GETPPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: UMIN, NATOMS, DMIN1, DMIN2, NATTEMPT, NCPU, MINSEP, BULKT, TWOD, ZSYM, DEBUG, BESTPATHLENGTH, ETS, &
  &               NPAIRFRQ, PAIR1, PAIR2, NPAIRFRQ, NPAIRDONE, MAXPAIRS, PERMDIST, BOXLX, BOXLY, BOXLZ, RIGIDBODY, BESTPATH, &
  &               BARRIERSHORT, EMIN, PLUS, MINUS, KPLUS, KMINUS, RATESHORT, ANGLEAXIS, NCONN, PFMIN, &
  &               NCONNMIN, TSTHRESH, NMIN, NTS, LOCATIONA, LOCATIONB, TOPPOINTER, POINTERM, POINTERP, NCONNMAX, NMINA, NMINB, &
  &               TEMPERATURE, REGROUPRATET, REGROUPPET, MINGROUP, MAXBARRIER
USE SAVESTATE
USE PORFUNCS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, J1, NSTEPS, J2, J3, J4, N1, N2, NLEFT, NUNCONA, NUNCONB, M1, M2, NCYCLE, DMIN, DMAX
INTEGER NMINSAVE, MINMAP(NMIN)
INTEGER, ALLOCATABLE :: MINLIST(:), POSITION(:)
INTEGER NCOL(NMIN), NDISTA(NMIN), NDISTB(NMIN), MINMAP(NMIN), MIN1, MIN2, NDEAD
DOUBLE PRECISION, ALLOCATABLE :: DISTLIST(:), BHEIGHT(:)
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), DISTANCE, RMAT(3,3), DIST2, PTIMESTAU(NMIN), LKSUM(NMIN), DINDEX(NMIN), &
  &              DUMMY, P1, P2, PBRANCHMIN1MIN2, PPROD
LOGICAL DEADTS(NTS), MATCHED, CHANGED
INTEGER, ALLOCATABLE :: NVAL(:,:),  NVTEMP(:)
DOUBLE PRECISION, ALLOCATABLE :: PBRANCH(:,:), PBTEMP(:)

10 CONTINUE
IF (NAVAIL.EQ.0) THEN
   PAIRSTODO=NCPU*NPAIRFRQ
   IF (NPAIRFRQ.LT.1) PAIRSTODO=NATTEMPT*NCPU ! just one set of pairs unless we run out
   CALL GETNCONN
!
!  NMINSAVE and MINMAP are just dummies here.
!
   NMINSAVE=NMIN
   DO J1=1,NMIN
      MINMAP(J1)=J1
   ENDDO
   CALL DIJKSTRA(NAVAIL,.FALSE.,0,NMINSAVE,MINMAP)
   NAVAIL=0
   NSTEPS=(BESTPATHLENGTH+1)/2
   ALLOCATE(MINLIST(NSTEPS))
   DO J1=1,NSTEPS
      MINLIST(J1)=BESTPATH(2*J1-1)
   ENDDO
50 CONTINUE
   PRINT '(A,I8,A)','getspair> best path contains ',NSTEPS,' minima'
   IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
   ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO),DISTLIST(PAIRSTODO))

!
!  Save state.
!
   ALLOCATE(EMINSAVE(NMIN),PFMINSAVE(NMIN),ETSSAVE(NTS),KPLUSSAVE(NTS),KMINUSSAVE(NTS),TOPPOINTERSAVE(NMIN), &
     &         PLUSSAVE(NTS),MINUSSAVE(NTS),POINTERMSAVE(NTS),POINTERPSAVE(NTS),MINGROUP(NMIN), &
     &         LOCATIONASAVE(NMINA),LOCATIONBSAVE(NMINB))
   NMINASAVE=NMINA; NMINBSAVE=NMINB; NMINSAVE=NMIN; NTSSAVE=NTS; LOCATIONASAVE(1:NMINA)=LOCATIONA(1:NMINA)
   LOCATIONBSAVE(1:NMINB)=LOCATIONB(1:NMINB); EMINSAVE(1:NMIN)=EMIN(1:NMIN); PFMINSAVE(1:NMIN)=PFMIN(1:NMIN)
   ETSSAVE(1:NTS)=ETS(1:NTS); KPLUSSAVE(1:NTS)=KPLUS(1:NTS); KMINUSSAVE(1:NTS)=KMINUS(1:NTS)
   TOPPOINTERSAVE(1:NMIN)=TOPPOINTER(1:NMIN); PLUSSAVE(1:NTS)=PLUS(1:NTS); MINUSSAVE(1:NTS)=MINUS(1:NTS)
   POINTERMSAVE(1:NTS)=POINTERM(1:NTS); POINTERPSAVE(1:NTS)=POINTERP(1:NTS)

   IF (REGROUPRATET.OR.REGROUPPET) THEN
      CALL REGROUPFREE
      TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                           ! before calling GETNCONN
      MAXBARRIER=HUGE(1.0D0)
   ENDIF
   CALL GETNCONN
!
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!
   CALL REGROUP(MINMAP)

   CALL RATECONST_SETUP(LKSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  PBRANCH(J1,M2) = KMC-type probability of taking connection J1 from minimum M2 to minimum NVAL(J1,M2)
!  Degenerate rearrangements are excluded.
!
   ALLOCATE(NVAL(NCONNMAX,NMIN),PBRANCH(NCONNMAX,NMIN),PBTEMP(NCONNMAX),NVTEMP(NCONNMAX))
   IF (NCONNMAX*1.0D0*NMIN*1.0D0*8*1.0D0.GT.1.0D9) PRINT '(A)', &
   &    'regroupfree> WARNING - about to try and allocate more than a Gb of RAM'
   IF (DEBUG) PRINT '(A,2I8)','regroupfree> allocating PBRANCH and NVAL: NCONNMAX, MIN1=',NCONNMAX, NMIN
   NCOL(1:NMIN)=0

   FROMLOOP: DO M2=1,NMIN
      J1=TOPPOINTER(M2)  !  sets J1 to the TS connected to minimum M2 with the highest id
      IF (J1.LE.0) CYCLE FROMLOOP
      DO WHILE (J1.GT.0) 
         IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
            MATCHED=.FALSE.
            IF (PLUS(J1).EQ.M2) THEN  !  M2 M1
               MATCHCOLP: DO M1=1,NCOL(M2)
                  IF (NVAL(M1,M2).EQ.MINUS(J1)) THEN ! A previous TS also links this pair
                     PBRANCH(M1,M2)=MIN(PBRANCH(M1,M2)+EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
                     MATCHED=.TRUE.
                     EXIT MATCHCOLP
                  ENDIF
               ENDDO MATCHCOLP
               IF (.NOT.MATCHED) THEN ! This minimum has not been connected to from M1 before
                  NCOL(M2)=NCOL(M2)+1 
                  NVAL(NCOL(M2),M2)=MINUS(J1)
                  PBRANCH(NCOL(M2),M2)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
               ENDIF
            ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2
               MATCHCOLM: DO M1=1,NCOL(M2)  
                  IF (NVAL(M1,M2).EQ.PLUS(J1)) THEN ! A PREVIOUS TS ALSO LINKS THIS PAIR
                     PBRANCH(M1,M2)=MIN(PBRANCH(M1,M2)+EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
                     MATCHED=.TRUE.
                     EXIT MATCHCOLM
                  ENDIF
               ENDDO MATCHCOLM
               IF (.NOT.MATCHED) THEN ! THIS MINIMUM HAS NOT BEEN CONNECTED TO FROM M1 BEFORE
                  NCOL(M2)=NCOL(M2)+1
                  NVAL(NCOL(M2),M2)=PLUS(J1)
                  PBRANCH(NCOL(M2),M2)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
               ENDIF
            ENDIF
         ENDIF
         IF (PLUS(J1).EQ.M2) THEN
            J1=POINTERP(J1)
         ELSE IF (MINUS(J1).EQ.M2) THEN
            J1=POINTERM(J1)
         ENDIF
      ENDDO
   ENDDO FROMLOOP
!
!  Check row normalisation.
!  
   IF (.TRUE.) THEN
      DO J1=1,NMIN
         DUMMY=0.0D0
         DO J2=1,NCOL(J1)
            DUMMY=DUMMY+PBRANCH(J2,J1)
            IF (DEBUG) WRITE(*,'(A,3I6,3G20.10)') 'J1,J2,NVAL,PBRANCH,sum=',J1,J2,NVAL(J2,J1),PBRANCH(J2,J1),DUMMY
         ENDDO 
         IF (DEBUG.AND.(NCOL(J1).GT.0)) WRITE(*,'(A,2I6,3G20.10)') 'J1,ncol,sum=',J1,NCOL(J1),DUMMY
         IF ((NCOL(J1).GT.0).AND.(ABS(DUMMY-1.0D0).GT.1.0D-10)) THEN
            PRINT*,'ERROR - J1,NCOL(J1),DUMMY=',J1,NCOL(J1),DUMMY
            STOP
         ENDIF
      ENDDO
   ENDIF
!
!  Check that the stationary point database is actually connected, and remove
!  minima that lie in disjoint graphs.
!  Calculate minimum number of steps of each minimum from the A set.
!
   DO J1=1,NMIN
      NDISTA(J1)=1000000
   ENDDO
   DO J1=1,NMINA
   NDISTA(LOCATIONA(J1))=0
   ENDDO
   NCYCLE=0
5  CHANGED=.FALSE.
   NCYCLE=NCYCLE+1
   DMIN=100000
   DMAX=0
   NUNCONA=0
   DO J1=1,NMIN
      IF (NDISTA(J1).EQ.0) CYCLE ! A MINIMUM
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
   PRINT '(3(A,I8))','regroupfree> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONA
!
!  Calculate minimum number of steps of each minimum from the B set.
!
   NDISTB(1:NMIN)=1000000
   DO J1=1,NMINB
      NDISTB(LOCATIONB(J1))=0
   ENDDO
   NCYCLE=0
51 CHANGED=.FALSE.
   NCYCLE=NCYCLE+1
   DMIN=100000
   DMAX=0
   NUNCONB=0
   DO J1=1,NMIN
      IF (NDISTB(J1).EQ.0) CYCLE ! B MINIMUM
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
   PRINT '(3(A,I8))','regroupfree> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!  This could happen if disconnected minima lie in the A or B region
   IF (NUNCONB.NE.NUNCONA) PRINT '(A)','regroupfree> WARNING - number of disconnected minima from A and B is different'
!
!  Remove disconnected minima from consideration.
!
   NLEFT=0
   DO J1=1,NMIN
      IF ((NDISTA(J1).EQ.1000000).OR.(NDISTB(J1).EQ.1000000)) THEN
         NCONN(J1)=0
         IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the A or B region'
      ENDIF
      IF (NCONN(J1).GT.NCONNMIN) NLEFT=NLEFT+1
   ENDDO
   PRINT '(A,I8)','regroupfree> Number of connected minima remaining with sufficient neighbours=',NLEFT

   DO MIN1=1,NMIN
      IF (NCONN(MIN1).GT.NCONNMIN) THEN
         PTIMESTAU(MIN1)=PFMIN(MIN1) ! PFMIN is ln(Z), so this is proportional to ln(p^eq*tau) 
!        PTIMESTAU(MIN1)=PFMIN(MIN1)-LKSUM(MIN1) ! PFMIN is ln(Z), so this is proportional to ln(p^eq*tau) 
!        PTIMESTAU(MIN1)=0.0D0
         GOTO 777
         DUMMY=EXP(PFMIN(MIN1))
         DO J2=1,NCOL(MIN1)
            MIN2=NVAL(J2,MIN1)
!
!  MIN1 must be in the neighbour list of MIN2.
!  We can use detailed balance to get the product of branching probabilities instead?
!  But it doesn't yield a helpful formula in the interesting case when 1-PPROD is small.
!  
            MATCHED=.FALSE.
            DO J3=1,NCOL(MIN2)
               IF (NVAL(J3,MIN2).EQ.MIN1) THEN
                  PBRANCHMIN1MIN2=PBRANCH(J3,MIN2)
                  PPROD=PBRANCH(J2,MIN1)*PBRANCH(J3,MIN2)
                  P1=0.0D0
                  DO J4=1,NCOL(MIN1)
                     IF (NVAL(J4,MIN1).EQ.MIN2) CYCLE
                     P1=P1+PBRANCH(J4,MIN1)
                  ENDDO
                  P2=0.0D0
                  DO J4=1,NCOL(MIN2)
                     IF (NVAL(J4,MIN2).EQ.MIN1) CYCLE
                     P2=P2+PBRANCH(J4,MIN2)
                  ENDDO
                  IF (DEBUG) PRINT '(A,2I8,A,2G20.10)','MIN1,MIN2=',MIN1,MIN2,' original 1-PPROD, alternative: ', &
  &                                              1.0D0-PPROD,P1+P2-P1*P2
                  PPROD=P1+P2-P1*P2
                  IF (PPROD.LE.0.0D0) THEN
                     PRINT '(A,G20.10,A,2I8)','ERROR in GT, 1-branching probability product is',PPROD,' MIN1, MIN2=',MIN1,MIN2
                     PRINT '(A,3G20.10)','             P1,P2,PPROD=',P1,P2,PPROD
                     STOP
                  ELSE
                     PTIMESTAU(MIN1)=PTIMESTAU(MIN1)+(EXP(-LKSUM(MIN2))+EXP(-LKSUM(MIN1))*PBRANCHMIN1MIN2)/PPROD
                     DUMMY=DUMMY+EXP(PFMIN(MIN2))
                     PRINT '(2I6,6F20.10)',MIN1,MIN2,PFMIN(MIN1),-LKSUM(MIN1),PFMIN(MIN2),-LKSUM(MIN2),PPROD,PBRANCHMIN1MIN2
                  ENDIF
                  MATCHED=.TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT.MATCHED) THEN
               PRINT '(A)','ERROR - failed match in GT - this should never happen'
               STOP
            ENDIF
         ENDDO
         PTIMESTAU(MIN1)=PTIMESTAU(MIN1)*DUMMY
777      CONTINUE
      ELSE
         PTIMESTAU(MIN1)=-HUGE(1.0D0)
      ENDIF
      DINDEX(MIN1)=MIN1*1.0D0
   ENDDO

   CALL DSORT(PTIMESTAU, DINDEX, NMIN, -2)
   DO J1=1,NLEFT
      PRINT '(I6,G20.10,F10.1,2F20.10)',J1,PTIMESTAU(J1),DINDEX(J1),-TEMPERATURE*PFMIN(INT(DINDEX(J1))),-LKSUM(INT(DINDEX(J1)))
      DO J2=1,NMINSAVE
!        IF (MINGROUP(J2).EQ.MINMAP(INT(DINDEX(J1)))) PRINT '(3(A,I6))','getppair> pe minimum ',J2,' original group ', &
! &              MINMAP(INT(DINDEX(J1))),' renumbered group ',INT(DINDEX(J1))
         IF (MINGROUP(J2).EQ.MINMAP(INT(DINDEX(J1)))) PRINT '(A,I6,3(A,G20.10))','getppair> pe minimum ',J2, &
  &              ' Z=',PFMIN(INT(DINDEX(J1))),' tau=',-LKSUM(INT(DINDEX(J1))), &
  &              ' Z*tau=',PFMIN(INT(DINDEX(J1)))-LKSUM(INT(DINDEX(J1)))
      ENDDO
!
!  We want to print the pe minima corresponding to free energy minima here. 
!  However, the original NMIN has been overwritten, and we have regrouped and
!  renumbered twice! 
!  Need to save the original stationary point information.
!  If we want to run this regrouping more than once then we need to
!  restore it after running GT/GT2/GETPPAIRS etc.
!  From REGROUP: MINMAP(J1), the location of minimum J1 in the original scheme.
!  If we have REGROUPRATET or REGROUPPE then the groups referred to by MINMAP
!  are free energy groups, not pe minima!
!  MINGROUP(J1) is the index of the group containing minimum J1
!
   ENDDO

   NMINA=NMINASAVE; NMINB=NMINBSAVE; NMIN=NMINSAVE; NTS=NTSSAVE; LOCATIONA(1:NMINA)=LOCATIONASAVE(1:NMINA)
   LOCATIONB(1:NMINB)=LOCATIONBSAVE(1:NMINB); EMIN(1:NMIN)=EMINSAVE(1:NMIN); PFMIN(1:NMIN)=PFMINSAVE(1:NMIN)
   ETS(1:NTS)=ETSSAVE(1:NTS); KPLUS(1:NTS)=KPLUSSAVE(1:NTS); KMINUS(1:NTS)=KMINUSSAVE(1:NTS)
   TOPPOINTER(1:NMIN)=TOPPOINTERSAVE(1:NMIN); PLUS(1:NTS)=PLUSSAVE(1:NTS); MINUS(1:NTS)=MINUSSAVE(1:NTS)
   POINTERM(1:NTS)=POINTERMSAVE(1:NTS); POINTERP(1:NTS)=POINTERPSAVE(1:NTS) 
   DEALLOCATE(EMINSAVE)
   DEALLOCATE(PFMINSAVE)
   DEALLOCATE(ETSSAVE)
   DEALLOCATE(KPLUSSAVE)
   DEALLOCATE(KMINUSSAVE)
   DEALLOCATE(TOPPOINTERSAVE)
   DEALLOCATE(PLUSSAVE)
   DEALLOCATE(MINUSSAVE)
   DEALLOCATE(POINTERMSAVE)
   DEALLOCATE(POINTERPSAVE)
   DEALLOCATE(MINGROUP)
   DEALLOCATE(LOCATIONASAVE)
   DEALLOCATE(LOCATIONBSAVE)


!    IF (BARRIERSHORT.OR.RATESHORT) THEN
! !
! ! Try connection pairs on either side of the highest barriers on the path, separated
! ! by MINSEP steps or fewer.
! !
!       ALLOCATE(BHEIGHT(NSTEPS-1),POSITION(NSTEPS-1))
!       DO J1=1,NSTEPS-1
!          IF (BARRIERSHORT) THEN
!             BHEIGHT(J1)=ETS(BESTPATH(2*J1))-EMIN(BESTPATH(2*J1+1))
!          ELSEIF (RATESHORT) THEN
!             IF (PLUS(BESTPATH(2*J1)).EQ.BESTPATH(2*J1+1)) THEN
!                BHEIGHT(J1)=-KPLUS(BESTPATH(2*J1))
!             ELSEIF (MINUS(BESTPATH(2*J1)).EQ.BESTPATH(2*J1+1)) THEN
!                BHEIGHT(J1)=-KMINUS(BESTPATH(2*J1))
!             ELSE
!                PRINT '(3(A,I6))','getspair> ERROR - ts=',BESTPATH(2*J1),' plus=',PLUS(BESTPATH(2*J1)),' minus=',MINUS(BESTPATH(2*J1))
!                PRINT '(A,I6)','getspair> ERROR - minimum id = ',BESTPATH(2*J1+1)
!             ENDIF
!          ENDIF
!          POSITION(J1)=2*J1 ! position in BESTPATH list
!       ENDDO
!       CALL SORT(NSTEPS-1,NSTEPS-1,BHEIGHT,POSITION)
!       IF (BARRIERSHORT) THEN
!          PRINT '(A)','sorted barriers on best path labelled according to the transition state:'
!          PRINT '(I6,G20.10)',(BESTPATH(POSITION(J1)),BHEIGHT(J1),J1=1,NSTEPS-1)
!       ELSE IF (RATESHORT) THEN
!          PRINT '(A)','sorted rates on best path labelled according to the transition state:'
!          PRINT '(I6,G20.10)',(BESTPATH(POSITION(J1)),-BHEIGHT(J1),J1=1,NSTEPS-1)
!       ENDIF
!       DO J1=1,NSTEPS-1 ! loop over barriers/rates
! !
! ! Do TS number POSITION(J1) in BESTPATH first etc.
! ! Must check that we aren't off the end of the best path, and that we
! ! haven't done this pair before.
! ! Exit when we have PAIRSTODO pairs.
! !
!          loop1: DO J2=1,MINSEP ! allow up to MINSEP steps away from the ts in question
!             N1=POSITION(J1)-(2*J2-1)
!             N2=POSITION(J1)+(2*J2-1)
! !           IF ((N1.LT.1).OR.(N2.GT.BESTPATHLENGTH)) EXIT      ! do not go off the end of the path!
! !           N1=BESTPATH(N1)
! !           N2=BESTPATH(N2)
!             IF ((N1.LT.1).AND.(N2.GT.BESTPATHLENGTH)) EXIT
!             N1=BESTPATH(MAX(N1,1))
!             N2=BESTPATH(MIN(N2,BESTPATHLENGTH))
!             DO J3=1,NPAIRDONE
!                IF ((PAIR1(J3).EQ.N1).AND.(PAIR2(J3).EQ.N2)) CYCLE loop1 ! do not repeat searches!
!                IF ((PAIR1(J3).EQ.N2).AND.(PAIR2(J3).EQ.N1)) CYCLE loop1 ! do not repeat searches!
!             ENDDO
!             NAVAIL=NAVAIL+1
!             DMIN1(NAVAIL)=N1
!             DMIN2(NAVAIL)=N2
!             DISTLIST(NAVAIL)=BHEIGHT(J1)
!             IF (NAVAIL.GE.PAIRSTODO) EXIT
!          ENDDO loop1
!          IF (NAVAIL.GE.PAIRSTODO) EXIT
!       ENDDO
!       DEALLOCATE(BHEIGHT,POSITION)
!    ELSE
!       DISTLIST(1:PAIRSTODO)=1.0D100
!       DO J1=1,NSTEPS
!          READ(UMIN,REC=MINLIST(J1)) (SPOINTS(J2),J2=1,3*NATOMS)
!          min2: DO J2=J1+1,NSTEPS
!             DO J3=1,NPAIRDONE
!                IF ((PAIR1(J3).EQ.MINLIST(J1)).AND.(PAIR2(J3).EQ.MINLIST(J2))) CYCLE min2 ! do not repeat searches
!                IF ((PAIR1(J3).EQ.MINLIST(J2)).AND.(PAIR2(J3).EQ.MINLIST(J1))) CYCLE min2 ! do not repeat searches
!             ENDDO
!             IF (J2-J1.GE.MINSEP) THEN ! find distance if separation is >= MINSEP
!                READ(UMIN,REC=MINLIST(J2)) (FPOINTS(J3),J3=1,3*NATOMS)
!                CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
!  ^                              RMAT,.FALSE.)
!                IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
!                                                         DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
!                NAVAIL=NAVAIL+1
!                sortloop: DO J3=1,MIN(NAVAIL,PAIRSTODO) ! sort the shortest PAIRSTODO values
!                   IF (DISTANCE.LT.DISTLIST(J3)) THEN
!                      DO J4=MIN(NAVAIL,PAIRSTODO),J3+1,-1
!                         DMIN1(J4)=DMIN1(J4-1)
!                         DMIN2(J4)=DMIN2(J4-1)
!                         DISTLIST(J4)=DISTLIST(J4-1)
!                      ENDDO
!                      DMIN1(J3)=MINLIST(J1)
!                      DMIN2(J3)=MINLIST(J2)
!                      DISTLIST(J3)=DISTANCE
!                      EXIT sortloop
!                   ENDIF
!                ENDDO sortloop
!                IF (DEBUG) PRINT '(3(A,I8),A,G20.10)','getspair> connection ',NAVAIL,' pair ',MINLIST(J1),  &
!      &                                               ' and ',MINLIST(J2),' distance=',DISTANCE
!             ENDIF
!          ENDDO min2
!       ENDDO
!    ENDIF
!    NAVAIL=MIN(NAVAIL,PAIRSTODO) 
!    PRINT '(A,I8,A)','getspair> sorted list of ',NAVAIL,' pairs'
!    PRINT '(2I8,F15.5)',(DMIN1(J1),DMIN2(J1),DISTLIST(J1),J1=1,NAVAIL)
!    DEALLOCATE(MINLIST,DISTLIST)
!    IF (NAVAIL.EQ.0) THEN
!       PRINT '(A)','getspair> No more candidate pairs of minima in getspair - quit'
!       STOP
!    ENDIF
!    NUSED=0

ENDIF

STOP

NUSED=NUSED+1
NAVAIL=NAVAIL-1
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)
WRITE(*,'(5(A,I8))') 'getspair> connecting minima ',MINS,' and ',MINF, ' pairs used=',  &
  &  NUSED,' remaining=',NAVAIL,' total pairs=',NPAIRDONE
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
CALL FLUSH(6)
READ(UMIN,REC=MINS) SPOINTS(1:3*NATOMS)
READ(UMIN,REC=MINF) FPOINTS(1:3*NATOMS)

END SUBROUTINE GETPPAIR
