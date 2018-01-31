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
!  Calculate rate constant by branching probability analysis combined with a waiting time,
!  but this time calculated by removing all the I minima and successively renormalising
!  the branching probabilities and waiting times.
!
      SUBROUTINE GT
      USE PORFUNCS
      USE SAVESTATE
      USE COMMONS
      USE GRAPHTRANSFORMATIONMODULE,ONLY:SPARSEGRAPHTRANSFORMATION,DenseGraphTransformation
      USE INPUTMODULE
      USE FREEMEMORYMODULE
      IMPLICIT NONE
      INTEGER M1, M2, J1, J2, J3, J4, MINMAP(NMIN), MIN1, MIN2, MIN3, NEWNCOL, NVALMIN1, NLEFT, J6, LASTMIN
      INTEGER NCOL(NMIN), NDEAD, NCONNMAXO, NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN, DMAX, NUNCONA, NUNCONB
      INTEGER NCOLOURS, NCONNDUMM(NMIN), INDEX(NMIN), NCOUNT, INVMAP(NMIN), NEWNDISTA(NMIN), NEWNDISTB(NMIN), NEWNCONN(NMIN)
      LOGICAL DEADTS(NTS), MATCHED, CHANGED
      INTEGER NC1, LASTONE, J5, NC2, NAVAIL, ISTAT, FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1)
      DOUBLE PRECISION EMKSUM(NMIN), PPROD, PBRANCHMIN1MIN2, COMMIT, KBA, KAB, NORM, P1, P2, DUMMY, LKSUM(NMIN), GBMAX
      DOUBLE PRECISION TNEW, ELAPSED, PEMKSUM(NMIN), NEWPFMIN(NMIN), NEWEMKSUM(NMIN)
      DOUBLE PRECISION, ALLOCATABLE :: PBRANCH(:,:), PBRANCHTEMP(:,:), PBTEMP(:), PBTEMPTEMP(:)
      DOUBLE PRECISION, ALLOCATABLE ::  PB1(:), PB1TEMP(:), PBRANCHTMP(:), PBTMP(:)
      INTEGER, ALLOCATABLE :: NVAL(:,:), NVALTEMP(:,:), NVTEMP(:), NVTEMPTEMP(:)
      INTEGER, ALLOCATABLE ::  NV1(:), NV1TEMP(:), MERGEDNV(:)
      CHARACTER(LEN=10) COLOURSTRING
      INTEGER, ALLOCATABLE :: NEWNVAL(:,:), NNEWNCOL(:), NVALTMP(:), NVTMP(:)
      DOUBLE PRECISION, ALLOCATABLE :: NEWPBRANCH(:,:), MERGEDPB(:)
      DOUBLE PRECISION DUMMYA, DUMMYB
C     INTEGER, PARAMETER :: CHANGETOFULL=1000
C     DOUBLE PRECISION PBRANCHFULL(CHANGETOFULL,CHANGETOFULL)

      CALL CPU_TIME(ELAPSED)
!
!  REGROUP and REGROUPFREE change NMINA, NMINB, LOCATIONA, LOCATIONB, so save the values and reset 
!  to call GT more than once.
!  In fact, REGROUPFREE changes NMIN, NTS, etc. so we must stop after such a run or
!  do a complete reset somehow.
!
!  Save state.
!
      ALLOCATE(EMINSAVE(NMIN),PFMINSAVE(NMIN),ETSSAVE(NTS),KPLUSSAVE(NTS),KMINUSSAVE(NTS),TOPPOINTERSAVE(NMIN), 
     &         PLUSSAVE(NTS),MINUSSAVE(NTS),POINTERMSAVE(NTS),POINTERPSAVE(NTS),LOCATIONASAVE(NMINA),LOCATIONBSAVE(NMINB))
      NMINASAVE=NMINA; NMINBSAVE=NMINB; NMINSAVE=NMIN; NTSSAVE=NTS; LOCATIONASAVE(1:NMINA)=LOCATIONA(1:NMINA)
      LOCATIONBSAVE(1:NMINB)=LOCATIONB(1:NMINB); EMINSAVE(1:NMIN)=EMIN(1:NMIN); PFMINSAVE(1:NMIN)=PFMIN(1:NMIN)
      ETSSAVE(1:NTS)=ETS(1:NTS); KPLUSSAVE(1:NTS)=KPLUS(1:NTS); KMINUSSAVE(1:NTS)=KMINUS(1:NTS)
      TOPPOINTERSAVE(1:NMIN)=TOPPOINTER(1:NMIN); PLUSSAVE(1:NTS)=PLUS(1:NTS); MINUSSAVE(1:NTS)=MINUS(1:NTS)
      POINTERMSAVE(1:NTS)=POINTERM(1:NTS); POINTERPSAVE(1:NTS)=POINTERP(1:NTS)

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
!
      CALL REGROUP(MINMAP)

      CALL RATECONST_SETUP(LKSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)

      DO J1=1,NMIN
         EMKSUM(J1)=EXP(-LKSUM(J1)) ! exponent minus local KSUM is the waiting time
!        PRINT '(A,I6,A,G20.10)','GT> min ',J1,' initial waiting time=',EMKSUM(J1)
      ENDDO
      PEMKSUM(1:NMIN)=EMKSUM(1:NMIN)

!!!!!!!!!!!!!!!!!!!   GT calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  PBRANCH(J1,M2) = KMC-type probability of taking connection J1 from minimum M2 to minimum NVAL(J1,M2)
!  Degenerate rearrangements are excluded.
!
!  NCONN may have been set to zero for minima in a disconnected region. DEADTS could still
!  be false for these, so exclude them explicitly.
!
      PRINT '(A,F4.1,A,I8,A,I6)','GT> about to try and allocate ',NCONNMAX*1.0D0*NMIN*1.0D0*8*2.0/1.0D9,'Gb of RAM for ',NMIN,
     &                           ' minima, maximum connectivity ',NCONNMAX
      CALL FLUSH(6)
      GBMAX=NCONNMAX*1.0D0*NMIN*1.0D0*8*2.0/1.0D9
      ALLOCATE(NVAL(NCONNMAX,NMIN),PBRANCH(NCONNMAX,NMIN),PBTEMP(NCONNMAX),NVTEMP(NCONNMAX))
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
C
C  Check row normalisation.
C
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
C
C  Reorder PBRANCH and NVAL so that the connected minima for each minimum
C  appear in ascending order. 
C
      ALLOCATE(NVALTMP(NCONNMAX),PBRANCHTMP(NCONNMAX))
      DO J1=1,NMIN
         IF (NCOL(J1).GT.1) THEN
            NVALTMP(1:NCOL(J1))=NVAL(1:NCOL(J1),J1)
            PBRANCHTMP(1:NCOL(J1))=PBRANCH(1:NCOL(J1),J1)
            CALL SORT4(NCOL(J1),NCONNMAX,PBRANCHTMP,NVALTMP)
            NVAL(1:NCOL(J1),J1)=NVALTMP(1:NCOL(J1))
            PBRANCH(1:NCOL(J1),J1)=PBRANCHTMP(1:NCOL(J1))
         ENDIF
      ENDDO
      DEALLOCATE(NVALTMP,PBRANCHTMP)
C
C  Check that the stationary point database is actually connected, and remove
C  minima that lie in disjoint graphs.
C  Calculate minimum number of steps of each minimum from the A set.
C
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
      PRINT '(3(A,I8))','GT> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONA
C
C  Calculate minimum number of steps of each minimum from the B set.
C
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
      PRINT '(3(A,I8))','GT> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!  This could happen if disconnected minima lie in the A or B region
      IF (NUNCONB.NE.NUNCONA) PRINT '(A)','GT> WARNING - number of disconnected minima from A and B is different'
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
      PRINT '(A,I8)','GT> Number of connected minima remaining with sufficient neighbours=',NLEFT

!     DO MIN2=1,NMIN
!            DO J3=1,NCOL(MIN2) ! debug printing
!               PRINT '(A,I6,A,I6,A,G20.10)','initial branching probability from ',MIN2,' to ',NVAL(J3,MIN2),' is ',PBRANCH(J3,MIN2)
!            ENDDO
!     ENDDO

      IF (GT2T) THEN ! SAT 
!         PRINT *,'GT> size(nconn)=',size(nconn)
          CALL READINPUT(NDISTA,NDISTB,NCOL,NVAL,PBRANCH,EMKSUM)
          IF (GT2SPARSE) THEN
               CALL SPARSEGRAPHTRANSFORMATION
          ELSE
               CALL DenseGraphTransformation
          ENDIF
          CALL FREEMEMORY
      ELSE ! SAT: Continue with David's version of GT

!
!  Remove I minima from the bottom up and renormalise the waiting times and
!  branching probabilities out of adjacent minima.
!  When removing minimum i the waiting times change for all minima, j \in J, adjacent
!  to i. The branching probabilities out of these minima also change.
!  For each j \in J new connections are created to minima in J that were not
!  originally connected to j. We finally end up with only A and B minima.
!  Must change the information in the following arrays:
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = id of minimum involved in connection J1 from minimum M2
!  PBRANCH(J1,M2) = branching probability of taking connection J1 from minimum M2 to minimum NVAL(J1,M2)
!  We have to DEALLOCATE and ALLOCATE NVAL and PBRANCH whenever necessary.
!  Waiting time is stored in EMKSUM. 
!  NCONNMAX is the maximum number of connections.
!
!  For NGT we use a different renormalisation and allow steps back to the starting minimum.
!  The initial return probability is zero.
!
      ALLOCATE(PB1(NCONNMAX),NV1(NCONNMAX))
      ALLOCATE(NVTMP(NCONNMAX),PBTMP(NCONNMAX))
      ALLOCATE(MERGEDNV(NCONNMAX),MERGEDPB(NCONNMAX))
      LASTMIN=NMIN-NMINA-NMINB
      DO MIN1=NMIN,1,-1
!        PRINT '(A)','GT> PBRANCH summary:'
!        DO J1=1,MIN1
!           PRINT '(A,I8)','PBRANCH for minimum ',J1 
!           PRINT '(6(I3,F15.5))',(NVAL(J2,J1),PBRANCH(J2,J1),J2=1,NCOL(J1))
!        ENDDO
         IF (NCONN(MIN1).LE.NCONNMIN) CYCLE ! minima with <= nconnmin connections are ignored
         IF ((NDISTA(MIN1).EQ.0).OR.(NDISTB(MIN1).EQ.0)) CYCLE
         NC1=NCOL(MIN1)
         PB1(1:NC1)=PBRANCH(1:NC1,MIN1)
         NV1(1:NC1)=NVAL(1:NC1,MIN1)
!
!  Cycle over adjacent minima of MIN1 and calculate their new waiting times
!  and branching probabilities.
!
         DO J2=1,NC1
            MIN2=NV1(J2)
!
!  MIN1 must be in the neighbour list of MIN2 unless MIN2 is a sink. 
!  MIN1 must be entry NCOL(MIN2) if we maintain sorted lists!
!  Could we use detailed balance to get the product of branching probabilities instead?
!
             J3=NCOL(MIN2)
             IF (NVAL(J3,MIN2).NE.MIN1) THEN
                PRINT *,'ERROR, MIN1 is not the last minimum in the list for MIN2'
                STOP
             ENDIF
             PBRANCHMIN1MIN2=PBRANCH(J3,MIN2)
             PPROD=PB1(J2)*PBRANCH(J3,MIN2)
             DUMMY=PPROD
C
C There seems to be only a small impact on execution time for doing this, so it could
C be the default. The loss of precision is a serious issue!
C
             IF (1.0D0-PPROD.LT.1.0D-1) THEN
                CALL PROBACC(NC1,NCONNMAX,NV1,MIN1,MIN2,PB1,PBRANCH,PPROD,NVAL,DEBUG,NCOL)
             ELSE
                PPROD=1.0D0-PPROD
             ENDIF
             NVALMIN1=J3
             IF (PPROD.LE.0.0D0) THEN
                PRINT '(A,G20.10,A,2I8)','ERROR in GT, 1-branching probability product is',PPROD,
     &                                     ' MIN1, MIN2=',MIN1,MIN2
                PRINT '(A,3G20.10)','             P1,P2,PPROD=',P1,P2,PPROD
                STOP
             ELSE
                PPROD=1.0D0/PPROD
             ENDIF
!            IF (MIN2.EQ.1) PRINT '(A,2I6,4G15.5)','MIN2,MIN1,EMKSUM(M1),EMKSUM(M1),P,new EMKSUM=',
!    &                      MIN2,MIN1,EMKSUM(MIN1),EMKSUM(MIN2),PBRANCHMIN1MIN2,
!    &                      (EMKSUM(MIN2)+EMKSUM(MIN1)*PBRANCHMIN1MIN2)*PPROD
!            PRINT '(A,I6,A,G20.10)','GT> initial waiting time for min ',MIN2,' is ',EMKSUM(MIN2)
!            PRINT '(A,I6,A,G20.10)','GT> new waiting time for min ',MIN2,' is ',(EMKSUM(MIN2)+EMKSUM(MIN1)*PBRANCHMIN1MIN2)*PPROD
             EMKSUM(MIN2)=(EMKSUM(MIN2)+EMKSUM(MIN1)*PBRANCHMIN1MIN2)*PPROD
!            PRINT '(A,3G20.10)','GT> PBRANCHMIN1MIN2,PPROD,EMKSUM(MIN1)=',PBRANCHMIN1MIN2,PPROD,EMKSUM(MIN1)
!
!  New branching probabilities for MIN2 to all its old neighbours and any neighbours
!  of MIN1 that it was not previously adjacent to.
!
             DO J3=1,NCOL(MIN2) ! old connections
                PBRANCH(J3,MIN2)=PBRANCH(J3,MIN2)*PPROD 
             ENDDO

             CALL NEWCONN(NC1,NV1,NCONNMAX,MIN1,MIN2,NCOL,NVAL,PBRANCH,PBRANCHMIN1MIN2,PB1,
     &                   NVTEMP,PBTEMP,PPROD,NEWNCOL)

!            DO J3=1,NCOL(MIN2) ! debug printing
!               PRINT '(A,I6,A,I6,A,G20.10)','new branching probability from ',MIN2,' to ',NVAL(J3,MIN2),' is ',PBRANCH(J3,MIN2)
!            ENDDO
!
!  Remove the MIN1 connection from MIN2 entries
!  It is simply the last entry for sorted lists!
!
             NCOL(MIN2)=NCOL(MIN2)-1
! 
!  Allocate more space if required. We don't need PBRANCH(x,y) and NVAL(x,y) beyond y=MIN1
!  so long as regrouping has been performed to put the A minima first, followed by the B minima,
!  followed by the I minima.
!  This reordering is done in the REGROUP subroutine.
!
             IF (NEWNCOL+NCOL(MIN2).GT.NCONNMAX) THEN
                PRINT '(A,F15.5,2(A,I8))','GT> allocating temporary array size ',8.0D0*NCONNMAX*MIN1*1.0D-6,
     &                                           ' maximum connections=',NCONNMAX,' remaining nodes=',MIN1
                ALLOCATE(PBRANCHTEMP(NCONNMAX,MIN1),NVALTEMP(NCONNMAX,MIN1),NVTEMPTEMP(NCONNMAX), 
     &                   PBTEMPTEMP(NCONNMAX),PB1TEMP(NCONNMAX),NV1TEMP(NCONNMAX))
                PBRANCHTEMP(1:NCONNMAX,1:MIN1)=PBRANCH(1:NCONNMAX,1:MIN1)
                NVALTEMP(1:NCONNMAX,1:MIN1)=NVAL(1:NCONNMAX,1:MIN1)
                PBTEMPTEMP(1:NCONNMAX)=PBTEMP(1:NCONNMAX)
                NVTEMPTEMP(1:NCONNMAX)=NVTEMP(1:NCONNMAX)
                PB1TEMP(1:NCONNMAX)=PB1(1:NCONNMAX)
                NV1TEMP(1:NCONNMAX)=NV1(1:NCONNMAX)
                DEALLOCATE(PBRANCH,NVAL,PBTEMP,NVTEMP,PB1,NV1,PBTMP,NVTMP,MERGEDNV,MERGEDPB)
                NCONNMAXO=NCONNMAX
                NCONNMAX=MAX(NCONNMAX,NEWNCOL+NCOL(MIN2))
                IF (8.0D0*NCONNMAX*MIN1*1.0D-9.GT.GBMAX) GBMAX=8.0D0*NCONNMAX*MIN1*1.0D-9
                PRINT '(A,F15.5,2(A,I8))','GT> allocating new array size       ',8.0D0*NCONNMAX*MIN1*1.0D-6,
     &                                           ' maximum connections=',NCONNMAX,' remaining nodes=',MIN1
                ALLOCATE(PBRANCH(NCONNMAX,MIN1),NVAL(NCONNMAX,MIN1),PBTEMP(NCONNMAX),NVTEMP(NCONNMAX))
                ALLOCATE(PB1(NCONNMAX),NV1(NCONNMAX))
                ALLOCATE(PBTMP(NCONNMAX),NVTMP(NCONNMAX))
                ALLOCATE(MERGEDNV(NCONNMAX),MERGEDPB(NCONNMAX))
                PBRANCH=0.0; NVAL=0; PBTEMP=0.0D0; NVTEMP=0; PB1=0.0D0; NV1=0
                PBRANCH(1:NCONNMAXO,1:MIN1)=PBRANCHTEMP(1:NCONNMAXO,1:MIN1)
                NVAL(1:NCONNMAXO,1:MIN1)=NVALTEMP(1:NCONNMAXO,1:MIN1)
                PBTEMP(1:NCONNMAXO)=PBTEMPTEMP(1:NCONNMAXO)
                NVTEMP(1:NCONNMAXO)=NVTEMPTEMP(1:NCONNMAXO)
                NV1(1:NCONNMAXO)=NV1TEMP(1:NCONNMAXO)
                PB1(1:NCONNMAXO)=PB1TEMP(1:NCONNMAXO)
                DEALLOCATE(PBRANCHTEMP,NVALTEMP,NVTEMPTEMP,PBTEMPTEMP,NV1TEMP,PB1TEMP)
                CALL FLUSH(6)
             ENDIF
!            IF (DEBUG) PRINT '(2(A,I8))','connectivity for minimum ',MIN2,' is now ',NCOL(MIN2)+NEWNCOL
!
! Maintain a sorted list by merging the two already sorted lists into one.
!
             IF (NEWNCOL.GT.0) CALL MYMERGE(NCOL(MIN2),NEWNCOL,NVAL,PBRANCH,NVTEMP,PBTEMP,NVTMP,
     &                         MERGEDNV,MERGEDPB,PBTMP,NCONNMAX,MIN1,MIN2,NCOL)

             NORM=0.0D0
             DO J3=1,NCOL(MIN2)
                NORM=NORM+PBRANCH(J3,MIN2)
             ENDDO
             IF (ABS(NORM-1.0D0).GT.1.0D-6) THEN
                PRINT '(A,I8,A,G20.10)','GT> WARNING - sum of branching probabilities out of minimum ',MIN2,
     &                                  ' deviates from unity: ',NORM
                IF (ABS(NORM-1.0D0).GT.1.0D-1) THEN
                   NORM=0.0D0
                   DO J3=1,NCOL(MIN2)
                      NORM=NORM+PBRANCH(J3,MIN2)
                      PRINT '(A,2I8,2G20.10)','GT> J3,NVAL,PBRANCH,SUM=',J3,NVAL(J3,MIN2),PBRANCH(J3,MIN2),NORM
                   ENDDO
                   PRINT '(A)','GT> Deviation is beyond tolerance'
                   STOP
                ENDIF
             ENDIF
         ENDDO
         IF (DEBUG) PRINT '(A,I8,A)','GT> minimum ',MIN1,' has been disconnected'
      ENDDO

      IF (DEBUG) THEN
         PRINT '(A)','GT> remaining waiting times and branching ratios:'
         DO J1=1,NMINA
            IF (NCOL(LOCATIONA(J1)).EQ.0) CYCLE
            PRINT '(A,I8,A,G20.10)','GT> A minimum ',LOCATIONA(J1),' tau =',EMKSUM(LOCATIONA(J1))
            DO J2=1,NCOL(LOCATIONA(J1))
               PRINT '(A,I8,A,G20.10)', 'GT> branch to minimum ',
     &               NVAL(J2,LOCATIONA(J1)),' probability=',PBRANCH(J2,LOCATIONA(J1))
            ENDDO
         ENDDO
         DO J1=1,NMINB
            IF (NCOL(LOCATIONB(J1)).EQ.0) CYCLE
            PRINT '(A,I8,A,G20.10)','GT> B minimum ',LOCATIONB(J1),' tau =',EMKSUM(LOCATIONB(J1))
            DO J2=1,NCOL(LOCATIONB(J1))
               PRINT '(A,I8,A,G20.10)', 'GT> branch to minimum ',
     &               NVAL(J2,LOCATIONB(J1)),' probability=',PBRANCH(J2,LOCATIONB(J1))
            ENDDO
         ENDDO
      ENDIF

      PRINT '(A)','remaining waiting times and branching ratios:'
      DO J1=1,NMINA
         IF (NCOL(LOCATIONA(J1)).EQ.0) CYCLE
         DUMMYA=0.0D0
         DUMMYB=0.0D0
         DO J2=1,NCOL(LOCATIONA(J1))
            DO J3=1,NMINA
               IF (NVAL(J2,LOCATIONA(J1)).EQ.LOCATIONA(J3)) THEN
                  DUMMYA=DUMMYA+PBRANCH(J2,LOCATIONA(J1))
                  GOTO 444
               ENDIF
            ENDDO
            DO J3=1,NMINB
               IF (NVAL(J2,LOCATIONA(J1)).EQ.LOCATIONB(J3)) THEN
                  DUMMYB=DUMMYB+PBRANCH(J2,LOCATIONA(J1))
                  GOTO 444
               ENDIF
            ENDDO
444         CONTINUE
         ENDDO
         PRINT '(A,I8,A,G20.10,3(A,G20.10))','GT> A min ',LOCATIONA(J1),' tau =',EMKSUM(LOCATIONA(J1)),
     &                           ' P_Aa=',DUMMYA,' P_Ba=',DUMMYB,' kNSS(B<-A) ratio=',DUMMYB/EMKSUM(LOCATIONA(J1))
      ENDDO
      DO J1=1,NMINB
         IF (NCOL(LOCATIONB(J1)).EQ.0) CYCLE
         DUMMYA=0.0D0
         DUMMYB=0.0D0
         DO J2=1,NCOL(LOCATIONB(J1))
            DO J3=1,NMINA
               IF (NVAL(J2,LOCATIONB(J1)).EQ.LOCATIONA(J3)) THEN
                  DUMMYA=DUMMYA+PBRANCH(J2,LOCATIONB(J1))
                  GOTO 555
               ENDIF
            ENDDO
            DO J3=1,NMINB
               IF (NVAL(J2,LOCATIONB(J1)).EQ.LOCATIONB(J3)) THEN
                  DUMMYB=DUMMYB+PBRANCH(J2,LOCATIONB(J1))
                  GOTO 555
               ENDIF
            ENDDO
555         CONTINUE
         ENDDO
         PRINT '(A,I8,A,G20.10,3(A,G20.10))','GT> B min ',LOCATIONB(J1),' tau =',EMKSUM(LOCATIONB(J1)),
     &                           ' P_Ab=',DUMMYA,' P_Bb=',DUMMYB,' kNSS(A<-B) ratio=',DUMMYA/EMKSUM(LOCATIONB(J1))
      ENDDO

      KBA=0.0D0
      DO J3=1,NMINA
         COMMIT=0.0D0
         DO J4=1,NCOL(LOCATIONA(J3)) ! sum over renormalised branching probabilities to B minima
            IF (NDISTB(NVAL(J4,LOCATIONA(J3))).EQ.0) THEN
               COMMIT=COMMIT+PBRANCH(J4,LOCATIONA(J3))
!              PRINT *,'J3,J4,PBRANCH,COMMIT=',J3,J4,PBRANCH(J4,LOCATIONA(J3)),COMMIT
            ENDIF
         ENDDO
         KBA=KBA+COMMIT*EXP(PFMIN(LOCATIONA(J3))-PFTOTALA)/EMKSUM(LOCATIONA(J3))
!        PRINT *,'J3,LOCATIONA,PFMIN,EMKSUM,KBA=',J3,LOCATIONA(J3),PFMIN(LOCATIONA(J3)),EMKSUM(LOCATIONA(J3)),KBA
      ENDDO

      KAB=0.0D0
      DO J3=1,NMINB
         COMMIT=0.0D0
         DO J4=1,NCOL(LOCATIONB(J3)) ! sum over renormalised branching probabilities to A minima
            IF (NDISTA(NVAL(J4,LOCATIONB(J3))).EQ.0) THEN
               COMMIT=COMMIT+PBRANCH(J4,LOCATIONB(J3))
!              PRINT *,'J3,J4,PBRANCH,COMMIT=',J3,J4,PBRANCH(J4,LOCATIONB(J3)),COMMIT
            ENDIF
         ENDDO
         KAB=KAB+COMMIT*EXP(PFMIN(LOCATIONB(J3))-PFTOTALB)/EMKSUM(LOCATIONB(J3))
!        PRINT *,'J3,LOCATIONB,PFMIN,EMKSUM,KAB=',J3,LOCATIONB(J3),PFMIN(LOCATIONB(J3)),EMKSUM(LOCATIONB(J3)),KAB
      ENDDO

      WRITE(*,'(A,F15.5,A)') 'GT> Peak array size for branching ratios=',GBMAX,' Gb'
      WRITE(*,'(4(A,G20.10))') 'GT> direction A<-B:  kNSS(A<-B)=',KAB,' kNSS(B<-A) (from detailed balance)=',
     &                              KAB*EXP(PFTOTALB-PFTOTALA),' total=',KAB*(1.0D0+EXP(PFTOTALB-PFTOTALA))
      WRITE(*,'(4(A,G20.10))') 'GT> direction B<-A:  kNSS(A<-B) (from detailed balance)=',KBA*EXP(PFTOTALA-PFTOTALB),
     &                            ' kNSS(B<-A)=',KBA,' total=',KBA*(1.0D0+EXP(PFTOTALA-PFTOTALB))
      WRITE(*,'(A,G20.10)') 'GT> detailed balance, ratio should be one if SS applies: ',KAB*EXP(PFTOTALB-PFTOTALA)/KBA
      WRITE(*,'(4(A,G20.10))') 'GT> kNSS(A<-B)=',KAB,' kNSS(B<-A) =',KBA,' total=',KAB+KBA

      ENDIF 

!     SAT Thu Mar 15 09:47:16 GMT 2007 -- arrays not correctly deallocated for GT2 runs -- checks added.
      if (allocated(PBRANCH)) deallocate(PBRANCH)
      if (allocated(NVAL)) deallocate(NVAL)
      if (allocated(PBTEMP)) deallocate(PBTEMP)
      if (allocated(NVTEMP)) deallocate(NVTEMP)
      if (allocated(PB1)) deallocate(PB1)
      if (allocated(NV1)) deallocate(NV1)
      if (allocated(NVTMP)) deallocate(NVTMP)
      if (allocated(PBTMP)) deallocate(PBTMP)
      if (allocated(MERGEDNV)) deallocate(MERGEDNV)
      if (allocated(MERGEDPB)) deallocate(MERGEDPB)

      CALL CPU_TIME(TNEW)
      TGT=TGT+TNEW-ELAPSED

      RETURN
      END
