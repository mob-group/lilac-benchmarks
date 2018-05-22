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
!  the branching probabilities and waiting times. In this routine we allow return to the
!  starting minimum.
!
      SUBROUTINE NGT
      USE PORFUNCS
      USE COMMONS
      USE NGTMEM
      IMPLICIT NONE
      INTEGER M1, M2, J1, J2, J3, J4, MINMAP(NMIN), NLEFT, J5, NDUMMY, NLABEL(NMIN), NCONNDUM
      INTEGER NCOL(NMIN), NDEAD, NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN, DMAX, NUNCONA, NUNCONB
      LOGICAL DEADTS(NTS), MATCHED, CHANGED
      INTEGER NAVAIL, ISTAT, NCONNMAXSAVE, NMINSAVE, NMINABSAVE, LNCONN(MAXMIN), FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1)
      DOUBLE PRECISION EMKSUM(NMIN), COMMIT, KBA, KAB, DUMMY, LKSUM(NMIN), GBMAX, SELF(NMIN)
      DOUBLE PRECISION TNEW, ELAPSED, PEMKSUM(NMIN)
      DOUBLE PRECISION DUMMYA, DUMMYB, KSSAB, KSSBA, EDUM
      INTEGER NMINASAVE, NMINBSAVE, NSUMAB, NTOP, NCDUM, NMAXDIM
      DOUBLE PRECISION, ALLOCATABLE :: PBRANCHSAVE(:,:), EMKSUMSAVE(:), BBRANCH(:)
      DOUBLE PRECISION, ALLOCATABLE :: PBRANCHSAVE2(:,:), EMKSUMSAVE2(:)
      DOUBLE PRECISION, ALLOCATABLE :: LPFOLDAB(:), LPFOLDBA(:)
      INTEGER, ALLOCATABLE :: NVALSAVE(:,:), NCOLSAVE(:), NCONNSAVE(:)
      INTEGER, ALLOCATABLE :: NVALSAVE2(:,:), NCOLSAVE2(:), NCONNSAVE2(:)
      LOGICAL, ALLOCATABLE :: BCON(:)
      INTEGER NONZERO, NCOUNT, NCOLPREV, ORIGNMINA, ORIGNMINB
      DOUBLE PRECISION LDUMMY, NEWPFOLD(NMIN), GPDIFF(NMIN)

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
         PRINT '(A)','NGT> calling persistence from NGT'
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

      ALLOCATE(LPFOLDAB(NMIN),LPFOLDBA(NMIN))
      LPFOLDAB=0.0D0
      LPFOLDBA=0.0D0

      DUMMY=0.0D0
      DO J1=1,NMIN
         EMKSUM(J1)=EXP(-LKSUM(J1)) ! exponent minus local KSUM is the waiting time
!        PRINT '(A,I6,A,G20.10)','NGT> min ',J1,' initial waiting time=',EMKSUM(J1)
         DUMMY=DUMMY+EMKSUM(J1)
      ENDDO
      DUMMY=DUMMY/NMIN
      PEMKSUM(1:NMIN)=EMKSUM(1:NMIN)
      DO J1=1,NMIN
         NLABEL(J1)=J1
      ENDDO
      CALL SORT(NMIN,NMIN,PEMKSUM,NLABEL)
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,G20.10,A)','NGT> Mean waiting time=',DUMMY,
     &                    ' Ten shortest waiting times for connected minima:'
      J1=NMIN; NDUMMY=0
      DO WHILE (NDUMMY.LE.MIN(10,NMIN))
         J1=J1-1
         IF (J1.LE.0) EXIT
         IF (NCONN(NLABEL(J1)).GT.NCONNMIN) THEN
            IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(I8,G20.10)',NLABEL(J1),PEMKSUM(J1)
            NDUMMY=NDUMMY+1
         ENDIF
      ENDDO
      PEMKSUM(1:NMIN)=EMKSUM(1:NMIN)

!!!!!!!!!!!!!!!!!!!   NGT calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  PBRANCH(J1,M2) = KMC-type probability of taking connection J1 from minimum M2 to minimum NVAL(J1,M2)
!  Degenerate rearrangements are excluded.
!
!  NCONN may have been set to zero for minima in a disconnected region. DEADTS could still
!  be false for these, so exclude them explicitly.
!  Initial self-connection branching probability is zero.
!
!  For NGT add an initially zero self-branching probability.
!
      NCONNMAXSAVE=NCONNMAX
      NMINSAVE=NMIN
      NMINABSAVE=NMINA+NMINB
      ORIGNMINA=NMINA
      ORIGNMINB=NMINB
      NCONNMAX=NCONNMAX+1
      GBMAX=NCONNMAX*1.0D0*NMIN*8.0D0/1.0D9

      IF (GBMAX*1.5D0 .GT. NGTCRSWITCH) THEN
         WRITE(*,*) 'NGT> memory required for full, rectangular storage = ',GBMAX*1.5D0,' Gb'
         WRITE(*,*) 'NGT> switching to compressed row storage instead'
         CALL NGT_CRSTORAGE(GBMAX,DEADTS,PEMKSUM,EMKSUM,LKSUM,NCOL,KBA,KAB,MINMAP,LPFOLDAB,LPFOLDBA,LNCONN)
         NCONNMAX=NMINA+NMINB
         GO TO 558 ! to do the disconnection of other sources and sinks using rectangular storage
      END IF

      IF (.NOT.SHANNONT) PRINT '(A,F4.1,A,I8,A,I6)','NGT> about to try and allocate ',NCONNMAX*1.0D0*NMIN*1.0D0*8*1.5D0/1.0D9, 
     &                            'Gb of RAM for ',NMIN,
     &                           ' minima, maximum connectivity ',NCONNMAX
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
C
C  Initial self-branching probability, but only for minima with enough connections!
C  Must maintain sorted lists!
C
C The number of connections in NCONN counts connections between the same minima through
C different transition states separately. Hence NCOL can end up less than NCONN.
C
C
      DO M2=1,NMIN
         IF (NCONN(M2).GT.NCONNMIN) THEN
            NCOL(M2)=NCOL(M2)+1
            IF (M2.GT.NVAL(NCOL(M2)-1,M2)) THEN
               NVAL(NCOL(M2),M2)=M2
               PBRANCH(NCOL(M2),M2)=0.0D0
            ELSE
               j3loop3: DO J3=1,NCOL(M2)-1
                  IF (M2.LT.NVAL(J3,M2)) THEN
!
! Move the rest up.
!
                     DO J4=NCOL(M2),J3+1,-1
                        NVAL(J4,M2)=NVAL(J4-1,M2)
                        PBRANCH(J4,M2)=PBRANCH(J4-1,M2)
                     ENDDO
                     NVAL(J3,M2)=M2
                     PBRANCH(J3,M2)=0.0D0
                     EXIT j3loop3
                  ENDIF
               ENDDO j3loop3
            ENDIF
!
!           NVAL(NCOL(M2),M2)=M2
!           PBRANCH(NCOL(M2),M2)=0.0D0
         ENDIF
      ENDDO
C
C  Check row normalisation.
C
      IF (.TRUE.) THEN
         DO J1=1,NMIN
            DUMMY=0.0D0
            DO J2=1,NCOL(J1)
               DUMMY=DUMMY+PBRANCH(J2,J1)
               IF (DEBUG) WRITE(*,'(A,3I6,3G20.10)') 'NGT> J1,J2,NVAL,PBRANCH,sum=',J1,J2,NVAL(J2,J1),PBRANCH(J2,J1),DUMMY
            ENDDO
            IF (DEBUG.AND.(NCOL(J1).GT.0)) WRITE(*,'(A,2I6,3G20.10)') 'NGT> J1,ncol,sum=',J1,NCOL(J1),DUMMY
            IF ((NCOL(J1).GT.0).AND.(ABS(DUMMY-1.0D0).GT.1.0D-10)) THEN
               WRITE (*,'(A,2I8,G20.10)') 'NGT> ERROR - J1,NCOL(J1),DUMMY=',J1,NCOL(J1),DUMMY
               STOP
            ENDIF
         ENDDO
      ENDIF
C
C  Reorder PBRANCH and NVAL so that the connected minima for each minimum
C  appear in ascending order. 
C  No longer necessary with sorted lists made from scratch.
C
!     ALLOCATE(NVALTMP(NCONNMAX),PBRANCHTMP(NCONNMAX))
!     DO J1=1,NMIN
!        IF (NCOL(J1).GT.1) THEN
!           NVALTMP(1:NCOL(J1))=NVAL(1:NCOL(J1),J1)
!           PBRANCHTMP(1:NCOL(J1))=PBRANCH(1:NCOL(J1),J1)
!           CALL SORT4(NCOL(J1),NCONNMAX,PBRANCHTMP,NVALTMP)
!           NVAL(1:NCOL(J1),J1)=NVALTMP(1:NCOL(J1))
!           PBRANCH(1:NCOL(J1),J1)=PBRANCHTMP(1:NCOL(J1))
!        ENDIF
!     ENDDO
!     DEALLOCATE(NVALTMP,PBRANCHTMP)
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
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(3(A,I8))','NGT> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',
     &                 DMAX,' disconnected=',NUNCONA
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
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(3(A,I8))','NGT> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',
     &                           DMAX,' disconnected=',NUNCONB
!  This could happen if disconnected minima lie in the A or B region
      IF (NUNCONB.NE.NUNCONA) PRINT '(A)','NGT> WARNING - number of disconnected minima from A and B is different'
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
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8)','NGT> Number of connected minima remaining with sufficient neighbours=',
     &              NLEFT
      LNCONN(1:NMIN)=NCONN(1:NMIN) ! save the nconn values for use in the pfold calculation at the end.
!
!  Remove I minima from the bottom up, i.e. from NMIN down to NMINA+NMINB+1.
!
      CALL NGTREMOVEI(NMIN,NMINA,NMINB,NCONNMAX,NCONNMIN,DEBUG,NCOL,GBMAX,EMKSUM,.TRUE.)
!
!  Having removed all the I minima we now have committor probabilities and we can calculate
!  k^SS and k^NSS.
!
      SELF(1:NMIN)=0.0D0
      IF (DEBUG) THEN
         IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A)','NGT> remaining waiting times and branching ratios:'
         DO J1=1,NMINA
            IF (NCOL(LOCATIONA(J1)).EQ.0) CYCLE
            IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,A,G20.10)','NGT> A minimum ',LOCATIONA(J1),' tau =', 
     &                                            EMKSUM(LOCATIONA(J1))
            DO J2=1,NCOL(LOCATIONA(J1))
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,A,G20.10)', 'NGT> branch to minimum ',
     &               NVAL(J2,LOCATIONA(J1)),' probability=',PBRANCH(J2,LOCATIONA(J1))
               IF (NVAL(J2,LOCATIONA(J1)).EQ.LOCATIONA(J1)) THEN
                  SELF(LOCATIONA(J1))=PBRANCH(J2,LOCATIONA(J1))
                  IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,G20.10)','NGT> P_aa=',SELF(LOCATIONA(J1))
               ENDIF
            ENDDO
         ENDDO
         DO J1=1,NMINB
            IF (NCOL(LOCATIONB(J1)).EQ.0) CYCLE
            IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,A,G20.10)','NGT> B minimum ',LOCATIONB(J1),' tau =',
     &                                            EMKSUM(LOCATIONB(J1))
            IF (.NOT.(RATESCYCLET)) PRINT '(A,I8,A,G20.10)','NGT> B minimum ',LOCATIONB(J1),' tau =',
     &                                            EMKSUM(LOCATIONB(J1))
            DO J2=1,NCOL(LOCATIONB(J1))
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,A,G20.10)', 'NGT> branch to minimum ',
     &               NVAL(J2,LOCATIONB(J1)),' probability=',PBRANCH(J2,LOCATIONB(J1))
               IF (NVAL(J2,LOCATIONB(J1)).EQ.LOCATIONB(J1)) THEN
                  SELF(LOCATIONB(J1))=PBRANCH(J2,LOCATIONB(J1))
                  IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,G20.10)','NGT> P_bb=',SELF(LOCATIONB(J1))
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A)',' '
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) 
     &   PRINT '(A)','NGT> remaining waiting times and branching ratios (these are committor probabilities):'
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A)',' '
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
         IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I7,A,G18.10,3(A,G18.10))','NGT> A min ',
     &      MINMAP(LOCATIONA(J1)),
     &                           ' tau =',EMKSUM(LOCATIONA(J1)),
     &                           ' P_Aa=',DUMMYA,' P_Ba=',DUMMYB,' kNSS(B<-A) ratio=',DUMMYB/EMKSUM(LOCATIONA(J1))
         IF (DEBUG.AND.(1.0D0-SELF(LOCATIONA(J1)).GT.0.0D0)) THEN
            PRINT '(A,I8,A,G18.10,3(A,G18.10))','NGT> A min ',LOCATIONA(J1),' tau/(1-P_aa)=',
     &                                           EMKSUM(LOCATIONA(J1))/(1.0D0-SELF(LOCATIONA(J1))),
     &                                           ' P_Ba/(1-P_aa)=',DUMMYB/(1.0D0-SELF(LOCATIONA(J1)))
         ENDIF
         LPFOLDAB(MINMAP(LOCATIONA(J1)))=DUMMYA
         LPFOLDBA(MINMAP(LOCATIONA(J1)))=DUMMYB
         IF (DIRECTION.EQ.'AB') THEN
            GPFOLD(LOCATIONA(J1))=DUMMYA
         ELSE
            GPFOLD(LOCATIONA(J1))=DUMMYB
         END IF
      ENDDO
      IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A)',' '
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
         IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I7,A,G18.10,3(A,G18.10))','NGT> B min ',MINMAP(LOCATIONB(J1)), 
     &                           ' tau =',EMKSUM(LOCATIONB(J1)),
     &                           ' P_Ab=',DUMMYA,' P_Bb=',DUMMYB,' kNSS(A<-B) ratio=',DUMMYA/EMKSUM(LOCATIONB(J1))
         IF (DEBUG.AND.(1.0D0-SELF(LOCATIONB(J1)).GT.0.0D0)) THEN
            PRINT '(A,I8,A,G18.10,3(A,G18.10))','NGT> B min ',LOCATIONB(J1),' tau/(1-P_bb)=',
     &                                           EMKSUM(LOCATIONB(J1))/(1.0D0-SELF(LOCATIONB(J1))),
     &                                           ' P_Ab/(1-P_bb)=',DUMMYA/(1.0D0-SELF(LOCATIONB(J1)))
         ENDIF
         LPFOLDAB(MINMAP(LOCATIONB(J1)))=DUMMYA
         LPFOLDBA(MINMAP(LOCATIONB(J1)))=DUMMYB
         IF (DIRECTION.EQ.'AB') THEN
            GPFOLD(LOCATIONB(J1))=DUMMYA
         ELSE
            GPFOLD(LOCATIONB(J1))=DUMMYB
         END IF
      ENDDO
      IF (.NOT.SHANNONT) PRINT '(A)',' '

      KBA=0.0D0
      KSSBA=0.0D0
      DO J3=1,NMINA
         COMMIT=0.0D0
         DO J4=1,NCOL(LOCATIONA(J3)) ! sum over renormalised branching probabilities to B minima
            IF (NDISTB(NVAL(J4,LOCATIONA(J3))).EQ.0) THEN
               COMMIT=COMMIT+PBRANCH(J4,LOCATIONA(J3))
            ENDIF
         ENDDO
         KBA  =KBA  +COMMIT*EXP(PFMIN(LOCATIONA(J3))-PFTOTALA)/EMKSUM(LOCATIONA(J3))
         KSSBA=KSSBA+COMMIT*EXP(PFMIN(LOCATIONA(J3))-PFTOTALA)/PEMKSUM(LOCATIONA(J3))
      ENDDO

      KAB=0.0D0
      KSSAB=0.0D0
      DO J3=1,NMINB
         COMMIT=0.0D0
         DO J4=1,NCOL(LOCATIONB(J3)) ! sum over renormalised branching probabilities to A minima
            IF (NDISTA(NVAL(J4,LOCATIONB(J3))).EQ.0) THEN
               COMMIT=COMMIT+PBRANCH(J4,LOCATIONB(J3))
            ENDIF
         ENDDO
         KAB  =KAB  +COMMIT*EXP(PFMIN(LOCATIONB(J3))-PFTOTALB)/EMKSUM(LOCATIONB(J3))
         KSSAB=KSSAB+COMMIT*EXP(PFMIN(LOCATIONB(J3))-PFTOTALB)/PEMKSUM(LOCATIONB(J3))
      ENDDO

      IF (.NOT.(SHANNONT)) THEN
      WRITE(*,'(A,F15.5,A)') 'NGT> Peak array size for branching ratios=',GBMAX,' Gb'
      WRITE(*,'(4(A,G20.10))') 'NGT> For direction A<-B:  kNSS(A<-B)=',KAB,' kNSS(B<-A) (detailed balance)=',
     &                              KAB*EXP(PFTOTALB-PFTOTALA) ! ,' total=',KAB*(1.0D0+EXP(PFTOTALB-PFTOTALA))
      WRITE(*,'(4(A,G20.10))') 'NGT> For direction B<-A:  kNSS(A<-B) (detailed balance)=',KBA*EXP(PFTOTALA-PFTOTALB),
     &                            ' kNSS(B<-A)=',KBA ! ,' total=',KBA*(1.0D0+EXP(PFTOTALA-PFTOTALB))
      WRITE(*,'(A,G20.10)') 'NGT> detailed balance for kSS,  ratio should be one if SS applies: ',KSSAB*EXP(PFTOTALB-PFTOTALA)/KSSBA
      WRITE(*,'(A,G20.10)') 'NGT> detailed balance for kNSS, ratio should be one if SS applies: ',KAB*EXP(PFTOTALB-PFTOTALA)/KBA
      WRITE(*,'(4(A,G20.10))') 'NGT> kSS (A<-B)=',KSSAB,'   kSS (B<-A)=',KSSBA
      WRITE(*,'(4(A,G20.10))') 'NGT> kNSS(A<-B)=',KAB,'   kNSS(B<-A)=',KBA
      ENDIF
      RATEAB=KAB
      RATEBA=KBA ! for communication with shannon subroutine

558   CONTINUE ! JMC jump to here after disconnection of I minima using the compressed-row storage scheme

! Write the committor probabilities for the end point minima to commit.ngt.AB for A as products and commit.ngt.BA for B as products.
      OPEN(UNIT=1,FILE='commit.ngt.AB',STATUS='UNKNOWN')
      WRITE(1,'(G20.10)') LPFOLDAB(1:NMIN)
      CLOSE(1)
      OPEN(UNIT=1,FILE='commit.ngt.BA',STATUS='UNKNOWN')
      WRITE(1,'(G20.10)') LPFOLDBA(1:NMIN)
      CLOSE(1)
      DEALLOCATE(LPFOLDAB,LPFOLDBA)

      IF (NGTDISCONNECTALL) THEN
!
!  For each remaining A or B minimum disconnect the other minima in the same set. 
!  This enables us to calculate k^KMC (now renamed simply k as in recent NGT paper) equivalent rate constants.
!  Can do this just by pretending that the other sources are type I and reordering the
!  lists appropriately. Need to save the current values of NMINA, NMINB, PBRANCH, NVAL and NCOL etc.
!
         NMINASAVE=NMINA; NMINBSAVE=NMINB
         NSUMAB=NMINA+NMINB
         ALLOCATE(BCON(NMINASAVE),BBRANCH(NMINASAVE))
         ALLOCATE(PBRANCHSAVE(NSUMAB,NSUMAB),NVALSAVE(NSUMAB,NSUMAB),NCOLSAVE(NSUMAB),EMKSUMSAVE(NSUMAB),NCONNSAVE(NSUMAB))
         ALLOCATE(PBRANCHSAVE2(NMINASAVE+1,NMINASAVE+1),NVALSAVE2(NMINASAVE+1,NMINASAVE+1),NCOLSAVE2(NMINASAVE+1), 
     &            EMKSUMSAVE2(NMINASAVE+1),NCONNSAVE2(NMINASAVE+1))
!
!  Without this initialisation we may access uninitialised memory.
!  The first elements of these arrays are not set directly, and
!  should be irrelevant, but best to set them!
!
         NCONNSAVE2(1:NMINASAVE+1)=0
         EMKSUMSAVE2(1:NMINASAVE+1)=0.0D0
         NCOLSAVE2(1:NMINASAVE+1)=0
         NVALSAVE2(1:NMINASAVE+1,1:NMINASAVE+1)=0
         PBRANCHSAVE2(1:NMINASAVE+1,1:NMINASAVE+1)=0.0D0
!
!  The first dimension of PBRANCH, NVAL should be the maximum number of connections.
!  This may be < NSUMAB, so need to get the dimension right!
!
         IF (DEBUG) PRINT '(A,2I6)','NGT> Dimensions of branching array are ',SIZE(PBRANCH,1),SIZE(PBRANCH,2)
         NMAXDIM=MIN(SIZE(PBRANCH,1),NSUMAB)
         PBRANCHSAVE(1:NMAXDIM,1:NSUMAB)=PBRANCH(1:NMAXDIM,1:NSUMAB)
         NVALSAVE(1:NMAXDIM,1:NSUMAB)=NVAL(1:NMAXDIM,1:NSUMAB)
         NCOLSAVE(1:NSUMAB)=NCOL(1:NSUMAB)
         EMKSUMSAVE(1:NSUMAB)=EMKSUM(1:NSUMAB)
         NCONNSAVE(1:NSUMAB)=NCONN(1:NSUMAB)
         NMIN=NMINASAVE+1
         NMINA=1
         NMINB=1
         IF (NMINASAVE.EQ.1) GOTO 557 ! preserve KBA
         PRINT '(A)','NGT> Sum branching probability to all sinks for each A minimum'
!
!  Sum branching probability to all sinks for each A minimum and save in 
!  BBRANCH(J4). Need a logical array BCON(J1) to remember if there is actually
!  a connection from A minimum J4 to a B minimum at all.
!
         DO J4=1,NMINASAVE
            BCON(J4)=.FALSE.
            BBRANCH(J4)=0.0D0
            DO J5=NCOLSAVE(J4),1,-1
               IF (NVAL(J5,J4).LE.NMINASAVE) EXIT
               BCON(J4)=.TRUE.
               BBRANCH(J4)=BBRANCH(J4)+PBRANCH(J5,J4)
            ENDDO
         ENDDO
!
!  Renumber the minima so that B sink = 1, then the A minima from 2 to NMINA+1. 
!  All the A minima move up one index. Create an entry for B with
!  zero branching probability if necessary to simplify later bookkeeping.
!  This could change the maximum connectivity, so need to reset NCONNMAX!!
!
         NCONNMAX=NCONNMAX+1
         DO J4=NMINASAVE,1,-1
            NCOLSAVE2(J4+1)=1
            EMKSUMSAVE2(J4+1)=EMKSUMSAVE(J4)
            NCONNSAVE2(J4+1)=NCONNSAVE(J4)
            NVALSAVE2(1,J4+1)=1 ! for connection to the B sink at position 1.
            PBRANCHSAVE2(1,J4+1)=BBRANCH(J4)
            DO J5=1,NCOLSAVE(J4)
               IF (NVALSAVE(J5,J4).GT.NMINASAVE) EXIT
               NVALSAVE2(J5+1,J4+1)=NVALSAVE(J5,J4)+1
               PBRANCHSAVE2(J5+1,J4+1)=PBRANCHSAVE(J5,J4)
               NCOLSAVE2(J4+1)=NCOLSAVE2(J4+1)+1
            ENDDO
         ENDDO
         NCOLSAVE2(1)=0 ! no connections out of sinks needed here.

         KBA=0.0D0
         ALLOCATE(NVALTMP(NMINASAVE+1),PBRANCHTMP(NMINASAVE+1))
!
!  For each A minimum remove all the other A sources and add this contribution to the 
!  overall rate constant.
!  Renumber the A minima so that on pass J3 through the loop the original minimum J3
!  becomes number 2. Subtract one from all the A indices and move 2 to NMINASAVE+2
!  NVAL and PBRANCH can be redimensioned in NGTREALLOC, called from 
!  NGTRENORM via NGTREMOVEID, so we need to reallocate if NMINASAVE > 1.
!
         PRINT '(A)','NGT> For each A minimum remove all the other A sources'
         DO J3=1,NMINASAVE
            IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I6)','NGT> doing minimum ',J3
!
!  Put the values for the renumbered minima into EMKSUM, NCOL, NVAL, PBRANCH
!  for renormalisation by NGTREMOVEI, which removes minima 3 to NMINASAVE+1
!
            IF (ALLOCATED(NVAL)) DEALLOCATE(NVAL)
            IF (ALLOCATED(PBRANCH)) DEALLOCATE(PBRANCH)
            ALLOCATE(NVAL(NCONNMAX,NMINASAVE+1),PBRANCH(NCONNMAX,NMINASAVE+1))
!
            DO J4=1,NMINASAVE+1
               NVAL(1:NCOLSAVE2(J4),J4)=NVALSAVE2(1:NCOLSAVE2(J4),J4)
               PBRANCH(1:NCOLSAVE2(J4),J4)=PBRANCHSAVE2(1:NCOLSAVE2(J4),J4)
            ENDDO
            EMKSUM(1:NMINASAVE+1)=EMKSUMSAVE2(1:NMINASAVE+1)
            NCONN(1:NMINASAVE+1)=NCONNSAVE2(1:NMINASAVE+1)
            NCOL(1:NMINASAVE+1)=NCOLSAVE2(1:NMINASAVE+1)

            IF (NCOL(2).GT.1) THEN
               CALL NGTREMOVEI(NMIN,NMINA,NMINB,NCONNMAX,NCONNMIN,DEBUG,NCOL,GBMAX,EMKSUM,.FALSE.)
!
!  There should now only be one source at position 2 and one sink at position 1.
!  PBRANCH(1,2) is not unity, in contrast to the BKL-type original GT scheme.
!
!              IF (DEBUG) PRINT '(A,I8,2(A,G20.10))','NGT> for A minimum ',J3,' P_{Ba}=',PBRANCH(1,2),' and time=',EMKSUM(2) 
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,2(A,G20.10))','NGT> for A minimum ',J3,' P_{Ba}=',PBRANCH(1,2),
     &                                                          ' and time=',EMKSUM(2) 
               KBA=KBA+PBRANCH(1,2)*EXP(PFMIN(LOCATIONA(J3))-PFTOTALA)/EMKSUM(2)
            ELSE
               PRINT '(A,I8,2(A,G20.10))','NGT> A minimum ',J3,' has no connections'
               PBRANCH(1,2)=0.0D0
               EMKSUM(2)=HUGE(1.0D0)
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,2(A,G20.10))','NGT> for A minimum ',J3,' P_{Ba}=',
     &                                     PBRANCH(1,2),' and time= infinity'
            ENDIF
!
!  Renumber A minima J1 -> J1-1 etc. with 2 -> NMINASAVE+1
!  NCOLSAVE2 and EMKSUMSAVE2 and the second index of NVALSAVE2 and PBRANCHSAVE2 
!  refer to absolute positions, so these need to change.
!
            NCDUM=NCOLSAVE2(2)
            EDUM=EMKSUMSAVE2(2)
            NCONNDUM=NCONNSAVE2(2)
            NVALTMP(1:NCDUM)=NVALSAVE2(1:NCDUM,2)
            PBRANCHTMP(1:NCDUM)=PBRANCHSAVE2(1:NCDUM,2)

            DO J4=2,NMINASAVE
               NVALSAVE2(2:NCOLSAVE2(J4+1),J4)=NVALSAVE2(2:NCOLSAVE2(J4+1),J4+1)-1
               PBRANCHSAVE2(1:NCOLSAVE2(J4+1),J4)=PBRANCHSAVE2(1:NCOLSAVE2(J4+1),J4+1)
               EMKSUMSAVE2(J4)=EMKSUMSAVE2(J4+1)
               NCONNSAVE2(J4)=NCONNSAVE2(J4+1)
               NCOLSAVE2(J4)=NCOLSAVE2(J4+1)
            ENDDO

            NVALSAVE2(2:NCDUM,NMINASAVE+1)=NVALTMP(2:NCDUM)-1
            PBRANCHSAVE2(1:NCDUM,NMINASAVE+1)=PBRANCHTMP(1:NCDUM)
            EMKSUMSAVE2(NMINASAVE+1)=EDUM
            NCONNSAVE2(NMINASAVE+1)=NCONNDUM
            NCOLSAVE2(NMINASAVE+1)=NCDUM
!
!  Need to reorder NVALSAVE2 and PBRANCHSAVE2 if they contain old minimum 2 in the list.
!  All A minima have a connection to the B sink in position 1, which must be the first
!  entry in NVALSAVE2, so we can check the second entry. If there is a connection to old minimum
!  2 it would have to occur at position 2 in NVALSAVE2 and PBRANCHSAVE2, and the index will have changed
!  to 1 above.
!
           DO J4=2,NMINASAVE+1
              IF (NVALSAVE2(2,J4).NE.1) CYCLE
              DUMMY=PBRANCHSAVE2(2,J4)
              DO J5=2,NCOLSAVE2(J4)-1
                 NVALSAVE2(J5,J4)=NVALSAVE2(J5+1,J4)
                 PBRANCHSAVE2(J5,J4)=PBRANCHSAVE2(J5+1,J4)
              ENDDO
              PBRANCHSAVE2(NCOLSAVE2(J4),J4)=DUMMY
              NVALSAVE2(NCOLSAVE2(J4),J4)=NMINASAVE+1
           ENDDO

         ENDDO
         DEALLOCATE(NVALTMP,PBRANCHTMP)
557      CONTINUE
         IF (NMINBSAVE.EQ.1) GOTO 556
!
!  Now do the B minima as sources.
!  Swap the A and B minima round then continue as for the A minima above.
!  This makes the bookkeeping much easier!
!
         IF (ALLOCATED(NVAL)) DEALLOCATE(NVAL)
         IF (ALLOCATED(PBRANCH)) DEALLOCATE(PBRANCH)
         ALLOCATE(NVAL(NCONNMAX,NMINASAVE+NMINBSAVE),PBRANCH(NCONNMAX,NMINASAVE+NMINBSAVE))
         DO J3=1,NMINBSAVE
            NCOL(J3)=NCOLSAVE(NMINASAVE+J3)
            EMKSUM(J3)=EMKSUMSAVE(NMINASAVE+J3)
            NCONN(J3)=NCONNSAVE(NMINASAVE+J3)
            DO J4=1,NCOL(J3)
               NVAL(J4,J3)=NVALSAVE(J4,NMINASAVE+J3)
               PBRANCH(J4,J3)=PBRANCHSAVE(J4,NMINASAVE+J3)
            ENDDO
         ENDDO
         DO J3=1,NMINASAVE
            NCOL(NMINBSAVE+J3)=NCOLSAVE(J3)
            EMKSUM(NMINBSAVE+J3)=EMKSUMSAVE(J3)
            NCONN(NMINBSAVE+J3)=NCONNSAVE(J3)
            DO J4=1,NCOL(NMINBSAVE+J3)
               NVAL(J4,NMINBSAVE+J3)=NVALSAVE(J4,J3)
               PBRANCH(J4,NMINBSAVE+J3)=PBRANCHSAVE(J4,J3)
            ENDDO
         ENDDO
         DO J3=1,NSUMAB
            DO J4=1,NCOL(J3)
               IF (NVAL(J4,J3).LE.NMINASAVE) THEN
                  NVAL(J4,J3)=NMINBSAVE+NVAL(J4,J3)
               ELSE
                  NVAL(J4,J3)=NVAL(J4,J3)-NMINASAVE
               ENDIF
            ENDDO
         ENDDO
!
! Resort the connections. They will be arranged original A then original B.
! Can speed things up using this information.
!
         ALLOCATE(NVALTMP(NSUMAB),PBRANCHTMP(NSUMAB))
         DO J3=1,NSUMAB
            NTOP=0 ! in case we don't execute the next loop!
            inloop: DO J4=2,NCOL(J3)
               NTOP=0
               IF (NVAL(J4-1,J3).GT.NVAL(J4,J3)) THEN
                  NTOP=J4-1
                  EXIT inloop
               ENDIF
            ENDDO inloop
            IF (NTOP.EQ.0) CYCLE
!
! Move entries 1 to NTOP after NTOP+1 to NCOL(J3) in both NVAL(:,J3) and PBRANCH(:,J3)
!
            NVALTMP(1:NCOL(J3))=NVAL(1:NCOL(J3),J3)
            PBRANCHTMP(1:NCOL(J3))=PBRANCH(1:NCOL(J3),J3)
            NVAL(1:NCOL(J3)-NTOP,J3)=NVALTMP(NTOP+1:NCOL(J3))
            NVAL(NCOL(J3)-NTOP+1:NCOL(J3),J3)=NVALTMP(1:NTOP)
            PBRANCH(1:NCOL(J3)-NTOP,J3)=PBRANCHTMP(NTOP+1:NCOL(J3))
            PBRANCH(NCOL(J3)-NTOP+1:NCOL(J3),J3)=PBRANCHTMP(1:NTOP)
         ENDDO
         DEALLOCATE(NVALTMP,PBRANCHTMP)
!
! Put these reordered values in the SAVE arrays.
!
         PBRANCHSAVE(1:NMAXDIM,1:NSUMAB)=PBRANCH(1:NMAXDIM,1:NSUMAB)
         NVALSAVE(1:NMAXDIM,1:NSUMAB)=NVAL(1:NMAXDIM,1:NSUMAB)
         NCOLSAVE(1:NSUMAB)=NCOL(1:NSUMAB)
         EMKSUMSAVE(1:NSUMAB)=EMKSUM(1:NSUMAB)
         NCONNSAVE(1:NSUMAB)=NCONN(1:NSUMAB)
         NDUMMY=NMINASAVE
         NMINASAVE=NMINBSAVE
         NMINBSAVE=NDUMMY
         NMIN=NMINASAVE+1
!
!  New dimension NMINASAVE could be larger, since it is really the B dimension!
!
         DEALLOCATE(BCON,BBRANCH,PBRANCHSAVE2,NVALSAVE2,NCOLSAVE2,EMKSUMSAVE2,NCONNSAVE2)
         ALLOCATE(BCON(NMINASAVE),BBRANCH(NMINASAVE))
         ALLOCATE(PBRANCHSAVE2(NMINASAVE+1,NMINASAVE+1),NVALSAVE2(NMINASAVE+1,NMINASAVE+1),NCOLSAVE2(NMINASAVE+1), 
     &            EMKSUMSAVE2(NMINASAVE+1),NCONNSAVE2(NMINASAVE+1))
!
!  Without these initialisations we may access uninitialised memory.
!
         NCONNSAVE2(1:NMINASAVE+1)=0
         EMKSUMSAVE2(1:NMINASAVE+1)=0.0D0
         NCOLSAVE2(1:NMINASAVE+1)=0
         NVALSAVE2(1:NMINASAVE+1,1:NMINASAVE+1)=0
         PBRANCHSAVE2(1:NMINASAVE+1,1:NMINASAVE+1)=0.0D0
!  
!  Sum branching probability to all sinks for each A (really B) minimum and save in  
!  BBRANCH(J1). Need a logical array BCON(J1) to remember if there is actually
!  a connection from A (really B) minimum J1 to a B (really A) minimum at all.
!
         PRINT '(A)','NGT> Sum branching probability to all sinks for each A (really B) minimum'
         DO J4=1,NMINASAVE
            BCON(J4)=.FALSE.
            BBRANCH(J4)=0.0D0
            DO J5=NCOLSAVE(J4),1,-1
               IF (NVAL(J5,J4).LE.NMINASAVE) EXIT
               BCON(J4)=.TRUE.
               BBRANCH(J4)=BBRANCH(J4)+PBRANCH(J5,J4)
            ENDDO
         ENDDO
!
!  Renumber the minima so that B (really A) i sink = 1, then the a (really b) minima.
!  All the A (really B) minima move up one index. Create an entry for B (really A) with
!  zero branching probability if necessary to simplify later bookkeeping.
!  
         DO J4=NMINASAVE,1,-1 
            NCOLSAVE2(J4+1)=1
            EMKSUMSAVE2(J4+1)=EMKSUMSAVE(J4)
            NCONNSAVE2(J4+1)=NCONNSAVE(J4)
            NVALSAVE2(1,J4+1)=1
            PBRANCHSAVE2(1,J4+1)=BBRANCH(J4)
            DO J5=1,NCOLSAVE(J4)
               IF (NVALSAVE(J5,J4).GT.NMINASAVE) EXIT
               NVALSAVE2(J5+1,J4+1)=NVALSAVE(J5,J4)+1
               PBRANCHSAVE2(J5+1,J4+1)=PBRANCHSAVE(J5,J4)
               NCOLSAVE2(J4+1)=NCOLSAVE2(J4+1)+1
            ENDDO
         ENDDO
         NCOLSAVE2(1)=0 ! no connections out of sinks needed here.

         ALLOCATE(NVALTMP(NMINASAVE+1),PBRANCHTMP(NMINASAVE+1))
         KAB=0.0D0
         DO J3=1,NMINASAVE
!  
!  Put the values for the renumbered minima into EMKSUM, NCOL, NVAL, PBRANCH
!  for renormalisation by NGTREMOVEI, which removes minima 3 to NMINASAVE+1
!  NVAL and PBRANCH can be redimensioned in NGTREALLOC, called from 
!  NGTRENORM via NGTREMOVEID, so we need to reallocate here.
!
            IF (ALLOCATED(NVAL)) DEALLOCATE(NVAL)
            IF (ALLOCATED(PBRANCH)) DEALLOCATE(PBRANCH)
            ALLOCATE(NVAL(NCONNMAX,NMINASAVE+1),PBRANCH(NCONNMAX,NMINASAVE+1))
            DO J4=1,NMINASAVE+1
               NVAL(1:NCOLSAVE2(J4),J4)=NVALSAVE2(1:NCOLSAVE2(J4),J4)
               PBRANCH(1:NCOLSAVE2(J4),J4)=PBRANCHSAVE2(1:NCOLSAVE2(J4),J4)
            ENDDO
            EMKSUM(1:NMINASAVE+1)=EMKSUMSAVE2(1:NMINASAVE+1)
            NCONN(1:NMINASAVE+1)=NCONNSAVE2(1:NMINASAVE+1)
            NCOL(1:NMINASAVE+1)=NCOLSAVE2(1:NMINASAVE+1)
            IF (NCOL(2).GT.1) THEN
               CALL NGTREMOVEI(NMIN,NMINA,NMINB,NCONNMAX,NCONNMIN,DEBUG,NCOL,GBMAX,EMKSUM,.FALSE.)
!  
!  There should now only be one source at position 2 and one sink at position 1.
!  PBRANCH(1,2) is not unity, in contrast to the BKL-type original GT scheme.
!
!              IF (DEBUG) PRINT '(A,I8,2(A,G20.10))','NGT> for B minimum ',J3,' P_{Ab}=',PBRANCH(1,2),' and time=',EMKSUM(2)
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,2(A,G20.10))','NGT> for B minimum ',J3,
     &                                     ' P_{Ab}=',PBRANCH(1,2),' and time=',EMKSUM(2)
!
!  Note that we have to use the relevant B partition function here, not A!
!
               KAB=KAB+PBRANCH(1,2)*EXP(PFMIN(LOCATIONB(J3))-PFTOTALB)/EMKSUM(2)
            ELSE
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,2(A,G20.10))','NGT> B minimum ',J3,' has no connections'
               PBRANCH(1,2)=0.0D0
               EMKSUM(2)=HUGE(1.0D0)
               IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,I8,2(A,G20.10))','NGT> for B minimum ',J3,
     &                                     ' P_{Ba}=',PBRANCH(1,2),' and time= infinity'
            ENDIF
!  
!  Renumber A (really B) minima J1 -> J1-1 etc. with 2 -> NMINASAVE+1
!  NCOLSAVE2 and EMKSUMSAVE2 and the second index of NVALSAVE2 and PBRANCHSAVE2
!  refer to absolute positions, so these need to change.
!
            NCDUM=NCOLSAVE2(2)
            EDUM=EMKSUMSAVE2(2)
            NCONNDUM=NCONNSAVE2(2)
            NVALTMP(1:NCDUM)=NVALSAVE2(1:NCDUM,2)
            PBRANCHTMP(1:NCDUM)=PBRANCHSAVE2(1:NCDUM,2)

            DO J4=2,NMINASAVE
               NVALSAVE2(2:NCOLSAVE2(J4+1),J4)=NVALSAVE2(2:NCOLSAVE2(J4+1),J4+1)-1
               PBRANCHSAVE2(1:NCOLSAVE2(J4+1),J4)=PBRANCHSAVE2(1:NCOLSAVE2(J4+1),J4+1)
               EMKSUMSAVE2(J4)=EMKSUMSAVE2(J4+1)
               NCONNSAVE2(J4)=NCONNSAVE2(J4+1)
               NCOLSAVE2(J4)=NCOLSAVE2(J4+1)
            ENDDO

            NVALSAVE2(2:NCDUM,NMINASAVE+1)=NVALTMP(2:NCDUM)-1
            PBRANCHSAVE2(1:NCDUM,NMINASAVE+1)=PBRANCHTMP(1:NCDUM)
            EMKSUMSAVE2(NMINASAVE+1)=EDUM
            NCONNSAVE2(NMINASAVE+1)=NCONNDUM
            NCOLSAVE2(NMINASAVE+1)=NCDUM
!  
!  Need to reorder NVALSAVE2 and PBRANCHSAVE2 if they contain old minimum 2 in the list.
!  All A minima have a connection to the B sink in position 1, which must be the first
!  entry in NVALSAVE2, so we can check the second entry. If there is a connection to old minimum 
!  2 it would have to occur at position 2 in NVALSAVE2 and PBRANCHSAVE2, and the index will have changed
!  to 1 above.
!
            DO J4=2,NMINASAVE+1
               IF (NVALSAVE2(2,J4).NE.1) CYCLE
               DUMMY=PBRANCHSAVE2(2,J4)
               DO J5=2,NCOLSAVE2(J4)-1
                  NVALSAVE2(J5,J4)=NVALSAVE2(J5+1,J4)
                  PBRANCHSAVE2(J5,J4)=PBRANCHSAVE2(J5+1,J4)
               ENDDO
               PBRANCHSAVE2(NCOLSAVE2(J4),J4)=DUMMY
               NVALSAVE2(NCOLSAVE2(J4),J4)=NMINASAVE+1
            ENDDO
         ENDDO
         DEALLOCATE(NVALTMP,PBRANCHTMP)

556      CONTINUE
         DEALLOCATE(PBRANCHSAVE,NVALSAVE,NCOLSAVE,EMKSUMSAVE,NCONNSAVE)
         DEALLOCATE(PBRANCHSAVE2,NVALSAVE2,NCOLSAVE2,EMKSUMSAVE2,NCONNSAVE2)
         DEALLOCATE(BCON,BBRANCH)
         WRITE(*,'(4(A,G20.10))') 'NGT> k(A<-B)=',KAB,'   k(B<-A)=',KBA
         WRITE(*,'(A,G20.10)') 'NGT> detailed balance for k, ratio should be one if SS applies: ',KAB*EXP(PFTOTALB-PFTOTALA)/KBA
         IF (RATESCYCLET) THEN
            WRITE(RATESUNIT,'(3G20.10)') TEMPERATURE,KAB,KBA
            IF (RATETARGETT) THEN
               IF ((KAB.GE.RATETARGETFRAC*RATETARGETAB).AND.(KBA.GE.RATETARGETFRAC*RATETARGETBA)) THEN
                  PRINT '(A)','NGT> Target rates achieved within 10% tolerance'
                  TARGETHIT=.TRUE.
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF (ALLOCATED(PBRANCH)) DEALLOCATE(PBRANCH)
      IF (ALLOCATED(NVAL)) DEALLOCATE(NVAL)

      CALL CPU_TIME(TNEW)
      TGT=TGT+TNEW-ELAPSED

      IF (NPFOLD.LE.0) RETURN ! Pfold calculation has not been requested.

! JMC Pfold calculation seeded with the committor probabilites that we have just calculated for the end points

      WRITE(*,'(A)') 'NGT> Pfold calculation seeded with the known committor probabilites for the end points'
      IF (NGTDISCONNECTALL) THEN ! need to reset these values
         NMIN=NMINSAVE
         NMINA=ORIGNMINA
         NMINB=ORIGNMINB
      END IF
      ALLOCATE(PBRANCH(NCONNMAXSAVE,NMIN),NVAL(NCONNMAXSAVE,NMIN))

! Subroutine MAKED2 is in Pfold.f90 and sets up the arrays as appropriate for the Pfold calculation (different from those 
! used in NGT above).
      CALL MAKED2(PBRANCH,NCOL,NCONNMAXSAVE,NVAL,DEADTS,LKSUM,LNCONN)

!  Now iterate GPFOLD's for a fixed number of cycles.
!  Gauss-Seidel iteration if OMEGA=1: successive over-relaxation if 1<OMEGA<2.
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
      ALLOCATE(DVEC(NONZERO),COL_IND(NONZERO),ROW_PTR(NMIN+1))
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
!  Main P^fold loop.
!  OMEGA is the damping factor for successive overrelaxation method (SOR)
!  OMEGA=1 is pure Gauss-Seidel. OMEGA should be < 2
!
      CALL CPU_TIME(TNEW)
      DO J1=1,NPFOLD
        ! CALL SPMV(GPFOLD, DVEC, NEWPFOLD, ROW_PTR, COL_IND, NMIN, NONZERO, NCOL)
        CALL SPMV_HARNESS(GPFOLD, NEWPFOLD, DVEC, ROW_PTR, COL_IND, NMIN-1)
      ENDDO

      CALL CPU_TIME(ELAPSED)
      WRITE(*,*) 'NGT> CPU time spent iterating committor probability =',ELAPSED-TNEW,' s'

! Now put the Pfold values into the ordering of the original database for writing to commit.data
      DO J3=1,NMIN
         NEWPFOLD(MINMAP(J3))=GPFOLD(J3)
      END DO
      WRITE(*,'(A)') 'NGT> Writing Pfold values to commit.data for minima in their original, possibly grouped, order (NOT RESORTED)'
      OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
      DO J3=1,NMIN
         WRITE(1,'(G20.10)') NEWPFOLD(J3)
      END DO
      CLOSE(1)

      DEALLOCATE(DVEC,COL_IND,ROW_PTR,NVAL,PBRANCH)

      RETURN
      END
