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

SUBROUTINE PFOLD
   USE COMMONS, ONLY: NMIN,NTS,MAXMIN,MAXTS,PLUS,MINUS,GPFOLD,OMEGA,DEBUG,KPLUS,KMINUS,MAXBARRIER, &
  &                  NPFOLD,TPFOLD,DIRECTION,NMINA,NMINB,LOCATIONA,LOCATIONB,NCONNMIN,ETS,EMIN,PFOLDCONV
   IMPLICIT NONE
   INTEGER J1, J2, J3, NDMAX, JMAX, NZERO, NCONNECTED, PNCONNECTED
   INTEGER LNCONN(MAXMIN), NDIST(NMIN), NCYCLE, DMIN, DMAX, NUNCON
   INTEGER NCOL(NMIN)
   LOGICAL CONNECTED(MAXMIN), DEADTS(MAXTS), CHANGED, LDEADTS
   DOUBLE PRECISION LDUMMY, NEWPFOLD(NMIN), ELAPSED, TNEW, LKSUM(NMIN), GPDIFF(NMIN)
   DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0
   INTEGER NCOLPREV, ROW_PTR(NMIN+1), NCOUNT, NONZERO
   INTEGER, ALLOCATABLE :: COL_IND(:), NVAL(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: DVEC(:), DMATMC(:,:)

   CALL CPU_TIME(ELAPSED)
!
!  Record the number of connections for each minimum in LNCONN.
!
   DO J1=1,NMIN
      CONNECTED(J1)=.TRUE.
   ENDDO
   NCONNECTED=0
11 DO J1=1,NMIN
      LNCONN(J1)=0
   ENDDO
   PNCONNECTED=NCONNECTED
   DO J1=1,NTS
! JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's in place of NCONN() for the plus and minus minima.
      CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),HUGE(1),HUGE(1), &
                   PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,LDEADTS)
      IF ((.NOT.LDEADTS) .AND. (PLUS(J1).NE.MINUS(J1))) THEN
         IF (CONNECTED(MINUS(J1))) LNCONN(PLUS(J1))=LNCONN(PLUS(J1))+1
         IF (CONNECTED(PLUS(J1)))  LNCONN(MINUS(J1))=LNCONN(MINUS(J1))+1
      ENDIF
   ENDDO
   NCONNECTED=0
   DO J1=1,NMIN
      CONNECTED(J1)=.FALSE.
      IF (LNCONN(J1).GT.NCONNMIN) THEN
         CONNECTED(J1)=.TRUE.
         NCONNECTED=NCONNECTED+1
      ENDIF
   ENDDO
!  IF (DEBUG) PRINT*,'NCONNECTED,PNCONNECTED=',NCONNECTED,PNCONNECTED
   IF (NCONNECTED.NE.PNCONNECTED) GOTO 11

   NDMAX=LNCONN(1)
   NZERO=0
   IF (LNCONN(1).EQ.0) NZERO=1
   JMAX=1
   DO J1=2,NMIN
      IF (LNCONN(J1).EQ.0) NZERO=NZERO+1
      IF (LNCONN(J1).GT.NDMAX) THEN
         NDMAX=LNCONN(J1)
         JMAX=J1
      ENDIF
   ENDDO
   WRITE(*,'(3(A,I6))') 'Pfold> largest connectivity is ',NDMAX,' for minimum ',JMAX, &
  &                     ' number with zero connections=',NZERO
   DEADTS(1:NTS)=.FALSE.
!  PRINT '(A,2I8,A,I8)','in Pfold, about to allocate NVAL, DMATMC: NDMAX, NMIN=',NDMAX, NMIN,' prod=',NDMAX*NMIN
!  PRINT '(A,2L5)','ALLOCATED(NVAL),ALLOCATED(DMATMC)=',ALLOCATED(NVAL),ALLOCATED(DMATMC)
   ALLOCATE(NVAL(NDMAX,NMIN),DMATMC(NDMAX,NMIN))

!!!!!!!!!!!!!!!!!!!   P^fold calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  IF DIRECTION is AB then we want P^fold->A and A minima are sinks
!  IF DIRECTION is BA then we want P^fold->B and B minima are sinks
!

   DO J1=1,NMIN
      LKSUM(J1)=0.0D0
   ENDDO
   DO J1=1,NTS
      CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),LNCONN(PLUS(J1)),LNCONN(MINUS(J1)), &
                   PLUS(J1),MINUS(J1),.TRUE.,CUT_UNDERFLOW,DEADTS(J1))
      IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
!        LKSUM(PLUS(J1))=LKSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
!        LKSUM(MINUS(J1))=LKSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
         LKSUM(PLUS(J1))=LKSUM(PLUS(J1))+EXP(KPLUS(J1))
         LKSUM(MINUS(J1))=LKSUM(MINUS(J1))+EXP(KMINUS(J1))
      ENDIF
   ENDDO
   DO J1=1,NMIN
      IF (LKSUM(J1).GT.0.0D0) THEN
!        LKSUM(J1)=LOG(LKSUM(J1))+KMEAN
         LKSUM(J1)=LOG(LKSUM(J1))
      ENDIF
   ENDDO

   CALL MAKED2(DMATMC,NCOL,NDMAX,NVAL,DEADTS,LKSUM,LNCONN)
!  
!  Check that the stationary point database is actually connected, and remove minima that lie in disjoint graphs.
!  Calculate minimum number of steps of each minimum from the A set.
!
   IF (DIRECTION.EQ.'AB') THEN
      NDIST(1:NMIN)=1000000
      DO J1=1,NMINA
         NDIST(LOCATIONA(J1))=0
      ENDDO 
      NCYCLE=0
5     CHANGED=.FALSE.
      NCYCLE=NCYCLE+1
      DMIN=100000
      DMAX=0
      NUNCON=0
      DO J1=1,NMIN
         IF (NDIST(J1).EQ.0) CYCLE ! A minimum
         DO J2=1,NCOL(J1)
            IF (NDIST(NVAL(J2,J1))+1.LT.NDIST(J1)) THEN 
               CHANGED=.TRUE.
               NDIST(J1)=NDIST(NVAL(J2,J1))+1
            ENDIF
         ENDDO
         IF ((NDIST(J1).GT.DMAX).AND.(NDIST(J1).NE.1000000)) DMAX=NDIST(J1)
         IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
         IF (NDIST(J1).EQ.1000000) NUNCON=NUNCON+1
      ENDDO
      IF (CHANGED) GOTO 5
      PRINT '(3(A,I8))','pfold> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCON
!
!  Might as well exclude disconnected minima.
!
      DO J1=1,NMIN
         IF (NDIST(J1).EQ.1000000) NCOL(J1)=0
      ENDDO
!
!  Calculate minimum number of steps of each minimum from the B set.
!
   ELSE
      NDIST(1:NMIN)=1000000
      DO J1=1,NMINB
         NDIST(LOCATIONB(J1))=0
      ENDDO
      NCYCLE=0
51    CHANGED=.FALSE.
      NCYCLE=NCYCLE+1
      DMIN=100000
      DMAX=0
      NUNCON=0
      DO J1=1,NMIN
         IF (NDIST(J1).EQ.0) CYCLE ! B minimum
         DO J2=1,NCOL(J1)
            IF (NDIST(NVAL(J2,J1))+1.LT.NDIST(J1)) THEN
               CHANGED=.TRUE.
               NDIST(J1)=NDIST(NVAL(J2,J1))+1
            ENDIF
         ENDDO
         IF ((NDIST(J1).GT.DMAX).AND.(NDIST(J1).NE.1000000)) DMAX=NDIST(J1)
         IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
         IF (NDIST(J1).EQ.1000000) NUNCON=NUNCON+1      
      ENDDO
      IF (CHANGED) GOTO 51
      PRINT '(3(A,I8))','pfold> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCON
!
!  Might as well exclude disconnected minima.
!
      DO J1=1,NMIN
         IF (NDIST(J1).EQ.1000000) NCOL(J1)=0
      ENDDO
   ENDIF
!
!  Now iterate GPFOLD's for a fixed number of cycles.
!  For DIRECTION AB calculate PFA directly: initial values for A minima are 1, the rest are 0.
!  For DIRECTION BA calculate PFB directly: initial values for B minima are 1, the rest are 0.
!  Gauss-Seidel iteration if OMEGA=1: successive over-relaxation if 1<OMEGA<2.
!  We get the same answer whether we calculate PFA or PFB, but it may converge
!  faster one way.
!
   NEWPFOLD(1:NMIN)=GPFOLD(1:NMIN)
   GPDIFF=0.0D0
!  
!  Make compressed row storage for DMAT.
!
   NONZERO=0
   DO J1=1,NMIN
      NONZERO=NONZERO+NCOL(J1)
   ENDDO
   ALLOCATE(DVEC(NONZERO),COL_IND(NONZERO))
   NCOUNT=0
   ROW_PTR(1)=1
   DO J1=1,NMIN
      IF (J1.GT.1) ROW_PTR(J1)=ROW_PTR(J1-1)+NCOLPREV
      NCOLPREV=NCOL(J1)
      DO J2=1,NCOL(J1)
         NCOUNT=NCOUNT+1
         DVEC(NCOUNT)=DMATMC(J2,J1)
         COL_IND(NCOUNT)=NVAL(J2,J1)
      ENDDO
   ENDDO
   ROW_PTR(NMIN+1)=NONZERO+1
!
!  Main P^fold loop.
!  OMEGA is the damping factor for successive overrelaxation method (SOR)
!  OMEGA=1 is pure Gauss-Seidel. OMEGA should be < 2
!
   IF (OMEGA .NE. 1.0) THEN
     STOP 1
   ENDIF

   DO J1=1,NPFOLD
     CALL SPMV_HARNESS(GPFOLD, NEWPFOLD, DVEC, ROW_PTR, COL_IND, NMIN-1)
     ! CALL SPMV(GPFOLD, DVEC, NEWPFOLD, ROW_PTR, COL_IND, NMIN, NONZERO, NCOL)
     NEWPFOLD(1:NMIN) = GPFOLD(1:NMIN)
   ENDDO

   DEALLOCATE(DVEC,COL_IND,NVAL,DMATMC)

   CALL CPU_TIME(TNEW)
   TPFOLD=TPFOLD+TNEW-ELAPSED

   RETURN

END SUBROUTINE PFOLD

!
! Calculate the mean wating time for a transition to any product minimum using an
! iterative first step type analysis as for Pfold.
!
SUBROUTINE TFOLD
USE COMMONS
USE PORFUNCS
IMPLICIT NONE
INTEGER J1, J2, J3, NDMAX, JMAX, NZERO, NCONNECTED, PNCONNECTED, J4, M1, M2, ISTAT, STEPMIN,NAVAIL
INTEGER LNCONN(MAXMIN), NDIST(NMIN), NCYCLE, DMIN, DMAX, NUNCON
INTEGER NCOL(NMIN), FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1)
LOGICAL CONNECTED(MAXMIN), DEADTS(MAXTS), CHANGED, MATCHED
DOUBLE PRECISION LDUMMY, NEWTFOLD(NMIN), ELAPSED, TNEW, LKSUM(NMIN), EMKSUM(NMIN), KAB, KBA, DEVIATION, DUMMY, XJ1
DOUBLE PRECISION KABOLD, KBAOLD
INTEGER NDEAD, NDISTA(NMIN), NDISTB(NMIN), NUNCONA, NUNCONB, NLEFT, MINMAP(NMIN)
INTEGER NCOLPREV, ROW_PTR(NMIN+1), NCOUNT, NONZERO
INTEGER, ALLOCATABLE :: COL_IND(:), NVAL(:,:)
DOUBLE PRECISION, ALLOCATABLE :: DVEC(:), PBRANCH(:,:)

CALL CPU_TIME(ELAPSED)
!
!  REGROUP and REGROUPFREE change NMINA, NMINB, LOCATIONA, LOCATIONB, so save the values and reset
!  to call GT more than once.
!  In fact, REGROUPFREE changes NMIN, NTS, etc. so we must stop after such a run or
!  do a complete reset somehow. This is done explicitly in routines like getppair
!  using the SAVESTATE module. Not needed here because we assume that NGT cannot be
!  called more than once!
!
IF (REGROUPFREET.OR.REGROUPFREEABT) THEN
   CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
   CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
   TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                        ! before calling GETNCONN
   MAXBARRIER=HUGE(1.0D0)
ENDIF

IF (REGROUPPERSISTT) THEN
   ! CALL GETNCONN done in getbarrier2 called from persistence
   CALL PERSISTENCE
   TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                        ! before calling GETNCONN
   MAXBARRIER=HUGE(1.0D0)
ENDIF


IF (REGROUPRATET.OR.REGROUPPET) THEN
   CALL REGROUPFREE
   TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                        ! before calling GETNCONN
   MAXBARRIER=HUGE(1.0D0)
ENDIF

CALL GETNCONN
!
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!  REGROUP reorders the minima A, B, I.
!
CALL REGROUP(MINMAP)

NZERO=0
DO J1=1,NMIN
   IF (NCONN(J1).EQ.0) NZERO=NZERO+1
ENDDO
WRITE(*,'(2(A,I6))') 'Pfold> largest connectivity is ',NCONNMAX,' number with zero connections=',NZERO

CALL RATECONST_SETUP(LKSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)

DO J1=1,NMIN
   EMKSUM(J1)=EXP(-LKSUM(J1)) ! exponent minus local KSUM is the waiting time
!  PRINT '(A,I6,A,G20.10)','Tfold> min ',J1,' initial waiting time=',EMKSUM(J1)
ENDDO

!!!!!!!!!!!!!!!!!!!   MFPT calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  PBRANCH(J1,M2) = KMC-type probability of taking connection J1 from minimum M2 to minimum NVAL(J1,M2)
!  Degenerate rearrangements are excluded.
!
!  NCONN may have been set to zero for minima in a disconnected region. DEADTS could still
!  be false for these, so exclude them explicitly.
!
PRINT '(A,F4.1,A,I8,A,I6)','Tfold> about to try and allocate ',NCONNMAX*1.0D0*NMIN*1.0D0*8*2.0/1.0D9,'Gb of RAM for ',NMIN, &
     &                     ' minima, maximum connectivity ',NCONNMAX
CALL FLUSH(6)
ALLOCATE(NVAL(NCONNMAX,NMIN),PBRANCH(NCONNMAX,NMIN))
NCOL(1:NMIN)=0
FROMLOOP: DO M2=1,NMIN
   IF (NCONN(M2).LE.NCONNMIN) CYCLE FROMLOOP
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
!
! Put it in the neighbour list, maintaining a sorted list.
!
               NCOL(M2)=NCOL(M2)+1
               IF (NCOL(M2).EQ.1) THEN
                  NVAL(NCOL(M2),M2)=MINUS(J1)
                  PBRANCH(NCOL(M2),M2)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
               ELSEIF (MINUS(J1).GT.NVAL(NCOL(M2)-1,M2)) THEN
                  NVAL(NCOL(M2),M2)=MINUS(J1)
                  PBRANCH(NCOL(M2),M2)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
               ELSE
                  j3loop: DO J3=1,NCOL(M2)-1
                     IF (MINUS(J1).LT.NVAL(J3,M2)) THEN
!
! Move the rest up.
!
                        DO J4=NCOL(M2),J3+1,-1
                           NVAL(J4,M2)=NVAL(J4-1,M2)
                           PBRANCH(J4,M2)=PBRANCH(J4-1,M2)
                        ENDDO
                        NVAL(J3,M2)=MINUS(J1)
                        PBRANCH(J3,M2)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
                        EXIT j3loop
                     ENDIF
                  ENDDO j3loop
               ENDIF
            ENDIF
         ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2
            MATCHCOLM: DO M1=1,NCOL(M2)
               IF (NVAL(M1,M2).EQ.PLUS(J1)) THEN ! A PREVIOUS TS ALSO LINKS THIS PAIR
                  PBRANCH(M1,M2)=MIN(PBRANCH(M1,M2)+EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
                  MATCHED=.TRUE.
                  EXIT MATCHCOLM
               ENDIF
            ENDDO MATCHCOLM
            IF (.NOT.MATCHED) THEN ! This minimum has not been connected to from M1 before
!
! Put it in the neighbour list, maintaining a sorted list.
! 
               NCOL(M2)=NCOL(M2)+1
               IF (NCOL(M2).EQ.1) THEN
                  NVAL(NCOL(M2),M2)=PLUS(J1)
                  PBRANCH(NCOL(M2),M2)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
               ELSEIF (PLUS(J1).GT.NVAL(NCOL(M2)-1,M2)) THEN
                  NVAL(NCOL(M2),M2)=PLUS(J1)
                  PBRANCH(NCOL(M2),M2)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
               ELSE
                  j3loop2: DO J3=1,NCOL(M2)-1
                     IF (PLUS(J1).LT.NVAL(J3,M2)) THEN
!
! Move the rest up.
! 
                        DO J4=NCOL(M2),J3+1,-1
                           NVAL(J4,M2)=NVAL(J4-1,M2)
                           PBRANCH(J4,M2)=PBRANCH(J4-1,M2)
                        ENDDO 
                        NVAL(J3,M2)=PLUS(J1)
                        PBRANCH(J3,M2)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
                        EXIT j3loop2
                     ENDIF
                  ENDDO j3loop2
               ENDIF
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
         WRITE (*,'(A,2I8,G20.10)') 'Tfold> ERROR - J1,NCOL(J1),DUMMY=',J1,NCOL(J1),DUMMY
         ! STOP
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
5     CHANGED=.FALSE.
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
PRINT '(3(A,I8))','Tfold> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONA
STEPMIN=DMAX
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
PRINT '(3(A,I8))','Tfold> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONB
STEPMIN=MAX(DMAX,STEPMIN)
!  This could happen if disconnected minima lie in the A or B region
IF (NUNCONB.NE.NUNCONA) PRINT '(A)','Tfold> WARNING - number of disconnected minima from A and B is different'
!
!  Remove disconnected minima from consideration.
!  
NLEFT=0
DO J1=1,NMIN
   IF ((NDISTA(J1).EQ.1000000).OR.(NDISTB(J1).EQ.1000000)) THEN
      NCONN(J1)=0
      NCOL(J1)=0
      IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the A or B region'
   ENDIF 
   IF (NCONN(J1).GT.NCONNMIN) NLEFT=NLEFT+1
ENDDO 
PRINT '(A,I8)','Tfold> Number of connected minima remaining with sufficient neighbours=',NLEFT
!
!  IF DIRECTION is AB then we want MFPT->A and A minima are sinks
!  IF DIRECTION is BA then we want MFPT->B and B minima are sinks
!  Iterate waiting times for a maximum number of cycles.
!  Gauss-Seidel iteration if TOMEGA=1: successive over-relaxation if 1<TOMEGA<2.
!  We get the same answer whether we calculate PFA or PFB, but it may converge
!  faster one way.
!  
!  Make compressed row storage for PBRANCH.
!
NONZERO=0
DO J1=1,NMIN
   NONZERO=NONZERO+NCOL(J1)
ENDDO
ALLOCATE(DVEC(NONZERO),COL_IND(NONZERO))
NCOUNT=0
ROW_PTR(1)=1
DO J1=1,NMIN
   IF (J1.GT.1) ROW_PTR(J1)=ROW_PTR(J1-1)+NCOLPREV
   NCOLPREV=NCOL(J1)
   DO J2=1,NCOL(J1)
      NCOUNT=NCOUNT+1
      DVEC(NCOUNT)=PBRANCH(J2,J1)
      COL_IND(NCOUNT)=NVAL(J2,J1)
   ENDDO
ENDDO
ROW_PTR(NMIN+1)=NONZERO+1
!
!  Main MFPT loop.
!  TOMEGA is the damping factor for successive overrelaxation method (SOR)
!  TOMEGA=1 is pure Gauss-Seidel. TOMEGA should be < 2
!
!  The waiting time for sinks should be zero. We might as well skip them in the 
!  loop over minima.
!
IF ((TOMEGA.LE.0.0D0).OR.(TOMEGA.GE.2.0D0)) THEN
   PRINT '(A,G20.10,A)','Tfold> ERROR, value of TOMEGA is out of range: ',TOMEGA,' resetting to one'
   TOMEGA=1.0D0
ENDIF

NEWTFOLD(1:NMINA)=EMKSUM(1:NMINA)
NEWTFOLD(NMINA+1:NMIN)=EMKSUM(NMINA+1:NMIN)
!
!  First do direction A<-B, so skip A minima (sinks).
!
KABOLD=0.0D0
DO XJ1=1.0D0,NTFOLD
   DEVIATION=0.0D0
   KAB=0.0D0
   DO J3=NMINA+1,NMIN
      IF (NCOL(J3).EQ.0) CYCLE
      LDUMMY=0.0D0
      DO J2=ROW_PTR(J3),ROW_PTR(J3+1)-1
         LDUMMY=LDUMMY-DVEC(J2)*NEWTFOLD(COL_IND(J2))    ! Gauss-Seidel
      ENDDO
!     IF (NEWTFOLD(J3).NE.0.0D0) DEVIATION=DEVIATION+((NEWTFOLD(J3)-EMKSUM(J3)+LDUMMY)/NEWTFOLD(J3))**2
      NEWTFOLD(J3)=TOMEGA*(EMKSUM(J3)-LDUMMY)+(1.0D0-TOMEGA)*NEWTFOLD(J3) !SOR
   ENDDO
!  IF (DEVIATION.NE.0.0D0) THEN
   IF (XJ1.GT.STEPMIN) THEN
      DO J3=NMINA+1,NMINA+NMINB
         KAB=KAB+EXP(PFMIN(J3)-PFTOTALB)/NEWTFOLD(J3)
      ENDDO
      IF (KAB.GT.0.0D0) THEN
         IF (ABS(KAB-KABOLD)/KAB.LT.TFOLDTHRESH) THEN
            PRINT '(A,F15.1,2(A,G20.10))','Tfold> deviation after ',XJ1,' steps=',ABS(KAB-KABOLD)/KAB,' kAB=',KAB
            EXIT
         ENDIF
      ENDIF
             
!     IF (DEVIATION.LT.TFOLDTHRESH) THEN
!        PRINT '(A,F15.1,2(A,G20.10))','Tfold> deviation after ',XJ1,' steps=',DEVIATION,' kAB=',KAB
!        EXIT
!     ENDIF

      IF (MOD(XJ1,1000.0D0).EQ.0) PRINT '(A,F15.1,2(A,G20.10))','Tfold> deviation after ',XJ1,' steps=', &
  &                                     ABS(KAB-KABOLD)/KAB,' kAB=',KAB
      KABOLD=KAB
   ENDIF
ENDDO 
IF (DEBUG) PRINT '(A)','final MFPT values in TFOLD:'
IF (DEBUG) PRINT '(6G20.10)',NEWTFOLD(1:NMIN)
!
!  Now direction B<-A, so skip B minima (sinks).
!
NEWTFOLD(1:NMINA)=EMKSUM(1:NMINA)
NEWTFOLD(NMINA+1:NMINA+NMINB)=0.0D0
NEWTFOLD(NMINB+1:NMIN)=EMKSUM(NMINB+1:NMIN)
NEWTFOLD(1:NMIN)=0.0D0
KBAOLD=0.0D0
DO XJ1=1,NTFOLD
   DEVIATION=0.0D0
   KBA=0.0D0
   DO J3=1,NMIN
      IF ((J3.GT.NMINA).AND.(J3.LE.NMINA+NMINB)) CYCLE
      IF (NCOL(J3).EQ.0) CYCLE
      LDUMMY=0.0D0
      DO J2=ROW_PTR(J3),ROW_PTR(J3+1)-1
         LDUMMY=LDUMMY-DVEC(J2)*NEWTFOLD(COL_IND(J2))    ! Gauss-Seidel
      ENDDO
!     IF (XJ1.GT.STEPMIN) DEVIATION=DEVIATION+((NEWTFOLD(J3)-EMKSUM(J3)+LDUMMY)/NEWTFOLD(J3))**2
      IF (NEWTFOLD(J3).GT.0.0D0) DEVIATION=DEVIATION+((NEWTFOLD(J3)-EMKSUM(J3)+LDUMMY)/NEWTFOLD(J3))**2
      NEWTFOLD(J3)=TOMEGA*(EMKSUM(J3)-LDUMMY)+(1.0D0-TOMEGA)*NEWTFOLD(J3) !SOR
   ENDDO
   IF (XJ1.GT.STEPMIN) THEN
!  IF (DEVIATION.NE.0.0D0) THEN
      DO J3=1,NMINA
         KBA=KBA+EXP(PFMIN(J3)-PFTOTALA)/NEWTFOLD(J3)
      ENDDO
      IF (KBA.GT.0.0D0) THEN
         IF (ABS(KBA-KBAOLD)/KBA.LT.TFOLDTHRESH) THEN
            PRINT '(A,F15.1,2(A,G20.10))','Tfold> deviation after ',XJ1,' steps=',ABS(KBA-KBAOLD)/KBA,' kBA=',KBA
            EXIT
         ENDIF
      ENDIF

!     IF (DEVIATION.LT.TFOLDTHRESH) THEN
!        PRINT '(A,F15.1,2(A,G20.10))','Tfold> deviation after ',XJ1,' steps=',DEVIATION,' kBA=',KBA
!        EXIT
!     ENDIF
      IF (MOD(XJ1,1000.0D0).EQ.0) PRINT '(A,F15.1,2(A,G20.10))','Tfold> deviation after ',XJ1,' steps=', &
  &                               ABS(KBA-KBAOLD)/KBA,' kBA=',KBA
      KBAOLD=KBA
   ENDIF
ENDDO 
IF (DEBUG) PRINT '(A)','Tfold> final MFPT values in TFOLD:'
IF (DEBUG) PRINT '(6G20.10)',NEWTFOLD(1:NMIN)
WRITE(*,'(4(A,G20.10))') 'Tfold> kKMC(A<-B)=',KAB,'   kKMC(B<-A)=',KBA
WRITE(*,'(A,G20.10)') 'Tfold> detailed balance for kKMC, ratio should be one if SS applies: ',KAB*EXP(PFTOTALB-PFTOTALA)/KBA

DEALLOCATE(DVEC,COL_IND,NVAL,PBRANCH)

CALL CPU_TIME(TNEW)
TTFOLD=TTFOLD+TNEW-ELAPSED

RETURN

END SUBROUTINE TFOLD

SUBROUTINE MAKED2(DMAT0,NCOL,NDMAX,NVAL,DEADTS,KSUM,LNCONN)
!
!  Construct DMATMC, which contains branching probabilities with minima in
!  regions A and B acting as sinks, i.e. no escape probability, except for
!  the B minimum NSOURCE.
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  DMAT0(J1,M2) = KMC-type probability of taking connection J1 from minimum 
!                 M2 to minimum NVAL(J1,M2)
!
!  Degenerate rearrangements are excluded.
!
   USE COMMONS
   IMPLICIT NONE
   INTEGER NCOL(NMIN), NDMAX, J1, M1, M2, NVAL(NDMAX,NMIN), OTHER, LNCONN(MAXMIN)
   LOGICAL DEADTS(NTS), MATCHED, ALLOWED
   DOUBLE PRECISION :: DMAT0(NDMAX,NMIN), KSUM(NMIN)
   LOGICAL ISA(NMIN), ISB(NMIN)
!
!  We are setting up DMAT0(thing,M2)
!  Cycle over connected transition states using pointers for minimum M2.
!
   NCOL(1:NMIN)=0
   ISA(1:NMIN)=.FALSE.
   DO J1=1,NMINA
      ISA(LOCATIONA(J1))=.TRUE.
   ENDDO
   ISB(1:NMIN)=.FALSE.
   DO J1=1,NMINB
      ISB(LOCATIONB(J1))=.TRUE.
   ENDDO

   fromloop: DO M2=1,NMIN   
      IF (LNCONN(M2).LE.NCONNMIN) CYCLE fromloop
      IF ((DIRECTION.EQ.'AB').AND.ISA(M2)) THEN
         CYCLE fromloop ! no escape from A
      ELSEIF ((DIRECTION.EQ.'BA').AND.ISB(M2)) THEN
         CYCLE fromloop ! no escape from B
      ENDIF
      J1=TOPPOINTER(M2)  !  sets J1 to the ts connected to minimum M2 with the highest value
!     PRINT '(A,2I8)','M2,TOPPOINTER(M2)=',M2,TOPPOINTER(M2)
      IF (J1.LE.0) CYCLE fromloop
11    CONTINUE
      IF ((PLUS(J1).NE.M2).AND.(MINUS(J1).NE.M2)) PRINT '(A,4I8)','ERROR, J1,M2,plus,minus=',J1,M2,PLUS(J1),MINUS(J1)
      IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
         IF (PLUS(J1).EQ.M2) THEN
            OTHER=MINUS(J1)
         ELSE 
            OTHER=PLUS(J1)
         ENDIF
         MATCHED=.FALSE.
         ALLOWED=.TRUE.
         IF ((DIRECTION.EQ.'BA').AND.ISA(OTHER)) THEN
            ALLOWED=.FALSE. ! no contribution from transitions into A minima
         ELSEIF ((DIRECTION.EQ.'AB').AND.ISB(OTHER)) THEN
            ALLOWED=.FALSE. ! no contribution from transitions into B minima
         ENDIF
         IF (ALLOWED) THEN
            IF (PLUS(J1).EQ.M2) THEN  !  M2 M1  
               matchcolp: DO M1=1,NCOL(M2)
                  IF (NVAL(M1,M2).EQ.MINUS(J1)) THEN ! a previous ts also links this pair
                     DMAT0(M1,M2)=MIN(DMAT0(M1,M2)+EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                     MATCHED=.TRUE.
                     EXIT matchcolp
                  ENDIF
               ENDDO matchcolp
               IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                  NCOL(M2)=NCOL(M2)+1
                  NVAL(NCOL(M2),M2)=MINUS(J1)
                  DMAT0(NCOL(M2),M2)=MIN(EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
               ENDIF
!              IF (DEBUG) PRINT '(A,3I6,G20.10)','Pfold> M2,NCOL,NVAL,DMAT0=', &
!    &                                            M2,NCOL(M2),NVAL(NCOL(M2),M2),DMAT0(NCOL(M2),M2)
            ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2 
               matchcolm: DO M1=1,NCOL(M2)
                  IF (NVAL(M1,M2).EQ.PLUS(J1)) THEN ! a previous ts also links this pair
                     DMAT0(M1,M2)=MIN(DMAT0(M1,M2)+EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                     MATCHED=.TRUE.
                     EXIT matchcolm
                  ENDIF
               ENDDO matchcolm
               IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                  NCOL(M2)=NCOL(M2)+1
                  NVAL(NCOL(M2),M2)=PLUS(J1)
                  DMAT0(NCOL(M2),M2)=MIN(EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
               ENDIF
!              IF (DEBUG) PRINT '(A,3I6,G20.10)','Pfold> M2,NCOL,NVAL,DMAT0=', &
!    &                                            M2,NCOL(M2),NVAL(NCOL(M2),M2),DMAT0(NCOL(M2),M2)
            ENDIF
         ENDIF
      ENDIF
      IF (PLUS(J1).EQ.M2) THEN
         J1=POINTERP(J1)
      ELSE IF (MINUS(J1).EQ.M2) THEN
         J1=POINTERM(J1)
      ENDIF
      IF (J1.GT.0) GOTO 11
   ENDDO fromloop

   RETURN

END SUBROUTINE MAKED2
