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
! JMC This routine is called from NGT, as an alternative to the main code there, if the memory requirements 
! are such that we have chosen to use the compressed row storage scheme for PBRANCH and NVAL
! (called DVEC and COL_IND, respectively, in the CR scheme)
!
      SUBROUTINE NGT_CRSTORAGE(GBMAX,DEADTS,PEMKSUM,EMKSUM,LKSUM,NCOL,KBA,KAB,MINMAP,LPFOLDAB,LPFOLDBA,LNCONN)
      USE PORFUNCS
      USE COMMONS
      USE NGTMEM
      IMPLICIT NONE
      INTEGER M1, M2, J1, J2, J3, J4, NLEFT, ISTAT
      INTEGER NCOL(NMIN), NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN, DMAX, NUNCONA, NUNCONB, MINMAP(NMIN)
      LOGICAL DEADTS(NTS), MATCHED, CHANGED
      DOUBLE PRECISION EMKSUM(NMIN), COMMIT, KBA, KAB, DUMMY, LKSUM(NMIN), GBMAX, SELF(NMIN)
      DOUBLE PRECISION PEMKSUM(NMIN), MEANNCONN
      DOUBLE PRECISION DUMMYA, DUMMYB, KSSAB, KSSBA
      DOUBLE PRECISION LPFOLDAB(NMIN), LPFOLDBA(NMIN)
      INTEGER NONZERO, NCOLPREV, NMIN_CONNECTED, NCONNTOT, LNCONN(MAXMIN)

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
! JMC set ROW_PTR initially using NCONN and including the self-connection branching.
! Note that this is an overestimate as NCONN treats multiple TSs between the same pair 
! of minima as extra connections but NCOL below does not.  But an overestimate is OK!

      NONZERO=0
      NCONNTOT=0
      NMIN_CONNECTED=0
      DO J1=1,NMIN
         IF (NCONN(J1).GT.NCONNMIN .AND. NCONN(J1).GT.1) THEN ! the 2nd criterion is because ncol = 2 (incl the self branching) is a special case in NGTrenorm_crstorage.f
            NCONNTOT=NCONNTOT+(NCONN(J1)+1)**2 ! the +1 is for the self-connection
            NMIN_CONNECTED=NMIN_CONNECTED+1
         END IF
      ENDDO
      MEANNCONN=DBLE(NCONNTOT) / DBLE(NMIN_CONNECTED)
      DO J1=1,NMIN
         IF (NCONN(J1).GT.NCONNMIN) NONZERO=NONZERO+NCONN(J1)+1+INT(MEANNCONN) ! the +1 is for the self-connection
      ENDDO
      WRITE(*,*)'NGT_crstorage> about to try and allocate ',NONZERO*8.0D0*1.5D0/1.0D9,'Gb of RAM for ',NMIN,
     &                           ' minima, maximum connectivity ',NCONNMAX,' average squared connectivity ',INT(MEANNCONN)
      CALL FLUSH(6)
      ALLOCATE(DVEC(NONZERO),COL_IND(NONZERO),ROW_PTR(NMIN))
      GBMAX=NONZERO*8.0D0/1.0D9
      ROW_PTR(1)=1
      IF (NCONN(1).GT.NCONNMIN) THEN
         NCOLPREV=NCONN(1)+1+INT(MEANNCONN)
      ELSE
         NCOLPREV=0
      END IF
      DO J1=2,NMIN
         ROW_PTR(J1)=ROW_PTR(J1-1)+NCOLPREV
         IF (NCONN(J1).GT.NCONNMIN) THEN
            NCOLPREV=NCONN(J1)+1+INT(MEANNCONN)
         ELSE
            NCOLPREV=0
         END IF
      ENDDO 

! JMC end of the compressed row storage additional setup

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
                     IF (COL_IND(M1+ROW_PTR(M2)-1).EQ.MINUS(J1)) THEN ! A previous TS also links this pair
                        DVEC(M1+ROW_PTR(M2)-1)=MIN(DVEC(M1+ROW_PTR(M2)-1)+EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
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
                        COL_IND(NCOL(M2)+ROW_PTR(M2)-1)=MINUS(J1)
                        DVEC(NCOL(M2)+ROW_PTR(M2)-1)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
                     ELSEIF (MINUS(J1).GT.COL_IND(NCOL(M2)-1+ROW_PTR(M2)-1)) THEN
                        COL_IND(NCOL(M2)+ROW_PTR(M2)-1)=MINUS(J1)
                        DVEC(NCOL(M2)+ROW_PTR(M2)-1)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
                     ELSE
                        j3loop: DO J3=1,NCOL(M2)-1
                           IF (MINUS(J1).LT.COL_IND(J3+ROW_PTR(M2)-1)) THEN
!
! Move the rest up.
!
                              DO J4=NCOL(M2),J3+1,-1
                                 COL_IND(J4+ROW_PTR(M2)-1)=COL_IND(J4-1+ROW_PTR(M2)-1)
                                 DVEC(J4+ROW_PTR(M2)-1)=DVEC(J4-1+ROW_PTR(M2)-1)
                              ENDDO 
                              COL_IND(J3+ROW_PTR(M2)-1)=MINUS(J1)
                              DVEC(J3+ROW_PTR(M2)-1)=MIN(EXP(KPLUS(J1)-LKSUM(PLUS(J1))),1.0D0)
                              EXIT j3loop
                           ENDIF
                        ENDDO j3loop
                     ENDIF
                  ENDIF
               ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2 
                  MATCHCOLM: DO M1=1,NCOL(M2)
                     IF (COL_IND(M1+ROW_PTR(M2)-1).EQ.PLUS(J1)) THEN ! A PREVIOUS TS ALSO LINKS THIS PAIR
                        DVEC(M1+ROW_PTR(M2)-1)=MIN(DVEC(M1+ROW_PTR(M2)-1)+EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
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
                        COL_IND(NCOL(M2)+ROW_PTR(M2)-1)=PLUS(J1)
                        DVEC(NCOL(M2)+ROW_PTR(M2)-1)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
                     ELSEIF (PLUS(J1).GT.COL_IND(NCOL(M2)-1+ROW_PTR(M2)-1)) THEN
                        COL_IND(NCOL(M2)+ROW_PTR(M2)-1)=PLUS(J1)
                        DVEC(NCOL(M2)+ROW_PTR(M2)-1)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
                     ELSE
                        j3loop2: DO J3=1,NCOL(M2)-1
                           IF (PLUS(J1).LT.COL_IND(J3+ROW_PTR(M2)-1)) THEN
!
! Move the rest up.
!
                              DO J4=NCOL(M2),J3+1,-1
                                 COL_IND(J4+ROW_PTR(M2)-1)=COL_IND(J4-1+ROW_PTR(M2)-1)
                                 DVEC(J4+ROW_PTR(M2)-1)=DVEC(J4-1+ROW_PTR(M2)-1)
                              ENDDO 
                              COL_IND(J3+ROW_PTR(M2)-1)=PLUS(J1)
                              DVEC(J3+ROW_PTR(M2)-1)=MIN(EXP(KMINUS(J1)-LKSUM(MINUS(J1))),1.0D0)
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
            IF (NCOL(M2).GT.NCONN(M2)) THEN ! JMC adding this extra test for now, for extra security
               WRITE(*,*)'NGT> Error, for minimum ',M2,' NCOL, NCONN=',NCOL(M2),NCONN(M2)
               STOP
            END IF
            NCOL(M2)=NCOL(M2)+1
            IF (M2.GT.COL_IND(NCOL(M2)-1+ROW_PTR(M2)-1)) THEN
               COL_IND(NCOL(M2)+ROW_PTR(M2)-1)=M2
               DVEC(NCOL(M2)+ROW_PTR(M2)-1)=0.0D0
            ELSE
               j3loop3: DO J3=1,NCOL(M2)-1
                  IF (M2.LT.COL_IND(J3+ROW_PTR(M2)-1)) THEN
!
! Move the rest up.
!
                     DO J4=NCOL(M2),J3+1,-1
                        COL_IND(J4+ROW_PTR(M2)-1)=COL_IND(J4-1+ROW_PTR(M2)-1)
                        DVEC(J4+ROW_PTR(M2)-1)=DVEC(J4-1+ROW_PTR(M2)-1)
                     ENDDO
                     COL_IND(J3+ROW_PTR(M2)-1)=M2
                     DVEC(J3+ROW_PTR(M2)-1)=0.0D0
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
               DUMMY=DUMMY+DVEC(J2+ROW_PTR(J1)-1)
               IF (DEBUG) WRITE(*,'(A,3I6,3G20.10)') 'NGT_crstorage> J1,J2,NVAL,PBRANCH,sum=',J1,J2,COL_IND(J2+ROW_PTR(J1)-1),
     &                                                                         DVEC(J2+ROW_PTR(J1)-1),DUMMY
            ENDDO
            IF (DEBUG.AND.(NCOL(J1).GT.0)) WRITE(*,'(A,2I6,3G20.10)') 'NGT_crstorage> J1,ncol,sum=',J1,NCOL(J1),DUMMY
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
            IF (NDISTA(COL_IND(J2+ROW_PTR(J1)-1))+1.LT.NDISTA(J1)) THEN
               CHANGED=.TRUE.
               NDISTA(J1)=NDISTA(COL_IND(J2+ROW_PTR(J1)-1))+1
            ENDIF
         ENDDO
         IF ((NDISTA(J1).GT.DMAX).AND.(NDISTA(J1).NE.1000000)) DMAX=NDISTA(J1)
         IF (NDISTA(J1).LT.DMIN) DMIN=NDISTA(J1)
         IF (NDISTA(J1).EQ.1000000) NUNCONA=NUNCONA+1
      ENDDO
      IF (CHANGED) GOTO 5
      PRINT '(3(A,I8))','NGT> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONA
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
            IF (NDISTB(COL_IND(J2+ROW_PTR(J1)-1))+1.LT.NDISTB(J1)) THEN
               CHANGED=.TRUE.
               NDISTB(J1)=NDISTB(COL_IND(J2+ROW_PTR(J1)-1))+1
            ENDIF
         ENDDO
         IF ((NDISTB(J1).GT.DMAX).AND.(NDISTB(J1).NE.1000000)) DMAX=NDISTB(J1)
         IF (NDISTB(J1).LT.DMIN) DMIN=NDISTB(J1)
         IF (NDISTB(J1).EQ.1000000) NUNCONB=NUNCONB+1
      ENDDO
      IF (CHANGED) GOTO 51
      PRINT '(3(A,I8))','NGT> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCONB
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
      PRINT '(A,I8)','NGT> Number of connected minima remaining with sufficient neighbours=',NLEFT
      LNCONN(1:NMIN)=NCONN(1:NMIN) ! save the nconn values for use in the pfold calculation at the end.
!
!  Remove I minima from the bottom up, i.e. from NMIN down to NMINA+NMINB+1.
!
      CALL NGTREMOVEI_CRSTORAGE(NMIN,NMINA,NMINB,NCONNMAX,NCONNMIN,DEBUG,NCOL,GBMAX,EMKSUM,.TRUE.)
!
!  Having removed all the I minima we now have committor probabilities and we can calculate
!  k^SS and k^NSS.
!
      SELF(1:NMIN)=0.0D0
      IF (DEBUG) THEN
         PRINT '(A)','NGT> remaining waiting times and branching ratios:'
         DO J1=1,NMINA
            IF (NCOL(LOCATIONA(J1)).EQ.0) CYCLE
            PRINT '(A,I8,A,G20.10)','NGT> A minimum ',LOCATIONA(J1),' tau =',EMKSUM(LOCATIONA(J1))
            DO J2=1,NCOL(LOCATIONA(J1))
               PRINT '(A,I8,A,G20.10)', 'NGT> branch to minimum ',
     &               COL_IND(J2+ROW_PTR(LOCATIONA(J1))-1),' probability=',DVEC(J2+ROW_PTR(LOCATIONA(J1))-1)
               IF (COL_IND(J2+ROW_PTR(LOCATIONA(J1))-1).EQ.LOCATIONA(J1)) THEN
                  SELF(LOCATIONA(J1))=DVEC(J2+ROW_PTR(LOCATIONA(J1))-1)
                  PRINT '(A,G20.10)','NGT> P_aa=',SELF(LOCATIONA(J1))
               ENDIF
            ENDDO
         ENDDO
         DO J1=1,NMINB
            IF (NCOL(LOCATIONB(J1)).EQ.0) CYCLE
            PRINT '(A,I8,A,G20.10)','NGT> B minimum ',LOCATIONB(J1),' tau =',EMKSUM(LOCATIONB(J1))
            DO J2=1,NCOL(LOCATIONB(J1))
               PRINT '(A,I8,A,G20.10)', 'NGT> branch to minimum ',
     &               COL_IND(J2+ROW_PTR(LOCATIONB(J1))-1),' probability=',DVEC(J2+ROW_PTR(LOCATIONB(J1))-1)
               IF (COL_IND(J2+ROW_PTR(LOCATIONB(J1))-1).EQ.LOCATIONB(J1)) THEN
                  SELF(LOCATIONB(J1))=DVEC(J2+ROW_PTR(LOCATIONB(J1))-1)
                  PRINT '(A,G20.10)','NGT> P_bb=',SELF(LOCATIONB(J1))
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      PRINT '(A)',' '
      PRINT '(A)','NGT> remaining waiting times and branching ratios (these are committor probabilities):'
      PRINT '(A)',' '
      DO J1=1,NMINA
         IF (NCOL(LOCATIONA(J1)).EQ.0) CYCLE
         DUMMYA=0.0D0
         DUMMYB=0.0D0
         DO J2=1,NCOL(LOCATIONA(J1))
            DO J3=1,NMINA
               IF (COL_IND(J2+ROW_PTR(LOCATIONA(J1))-1).EQ.LOCATIONA(J3)) THEN
                  DUMMYA=DUMMYA+DVEC(J2+ROW_PTR(LOCATIONA(J1))-1)
                  GOTO 444
               ENDIF
            ENDDO
            DO J3=1,NMINB
               IF (COL_IND(J2+ROW_PTR(LOCATIONA(J1))-1).EQ.LOCATIONB(J3)) THEN
                  DUMMYB=DUMMYB+DVEC(J2+ROW_PTR(LOCATIONA(J1))-1)
                  GOTO 444
               ENDIF
            ENDDO
444         CONTINUE
         ENDDO
         PRINT '(A,I7,A,G18.10,3(A,G18.10))','NGT_crstorage> A min ',MINMAP(LOCATIONA(J1)),
     &                           ' tau =',EMKSUM(LOCATIONA(J1)),
     &                           ' P_Aa=',DUMMYA,' P_Ba=',DUMMYB,' kNSS(B<-A) ratio=',DUMMYB/EMKSUM(LOCATIONA(J1))
         IF (DEBUG.AND.(1.0D0-SELF(LOCATIONA(J1)).GT.0.0D0)) THEN
            PRINT '(A,I8,A,G18.10,3(A,G18.10))','NGT_crstorage> A min ',LOCATIONA(J1),' tau/(1-P_aa)=',
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
      PRINT '(A)',' '
      DO J1=1,NMINB
         IF (NCOL(LOCATIONB(J1)).EQ.0) CYCLE
         DUMMYA=0.0D0
         DUMMYB=0.0D0
         DO J2=1,NCOL(LOCATIONB(J1))
            DO J3=1,NMINA
               IF (COL_IND(J2+ROW_PTR(LOCATIONB(J1))-1).EQ.LOCATIONA(J3)) THEN
                  DUMMYA=DUMMYA+DVEC(J2+ROW_PTR(LOCATIONB(J1))-1)
                  GOTO 555
               ENDIF
            ENDDO
            DO J3=1,NMINB
               IF (COL_IND(J2+ROW_PTR(LOCATIONB(J1))-1).EQ.LOCATIONB(J3)) THEN
                  DUMMYB=DUMMYB+DVEC(J2+ROW_PTR(LOCATIONB(J1))-1)
                  GOTO 555
               ENDIF
            ENDDO
555         CONTINUE
         ENDDO
         PRINT '(A,I7,A,G18.10,3(A,G18.10))','NGT_crstorage> B min ',MINMAP(LOCATIONB(J1)),
     &                           ' tau =',EMKSUM(LOCATIONB(J1)),
     &                           ' P_Ab=',DUMMYA,' P_Bb=',DUMMYB,' kNSS(A<-B) ratio=',DUMMYA/EMKSUM(LOCATIONB(J1))
         IF (DEBUG.AND.(1.0D0-SELF(LOCATIONB(J1)).GT.0.0D0)) THEN
            PRINT '(A,I8,A,G18.10,3(A,G18.10))','NGT_crstorage> B min ',LOCATIONB(J1),' tau/(1-P_bb)=',
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
      PRINT '(A)',' '

      KBA=0.0D0
      KSSBA=0.0D0
      DO J3=1,NMINA
         COMMIT=0.0D0
         DO J4=1,NCOL(LOCATIONA(J3)) ! sum over renormalised branching probabilities to B minima
            IF (NDISTB(COL_IND(J4+ROW_PTR(LOCATIONA(J3))-1)).EQ.0) THEN
               COMMIT=COMMIT+DVEC(J4+ROW_PTR(LOCATIONA(J3))-1)
!              PRINT *,'J3,J4,PBRANCH,COMMIT=',J3,J4,DVEC(J4+ROW_PTR(LOCATIONA(J3))-1),COMMIT
            ENDIF
         ENDDO
         KBA  =KBA  +COMMIT*EXP(PFMIN(LOCATIONA(J3))-PFTOTALA)/EMKSUM(LOCATIONA(J3))
         KSSBA=KSSBA+COMMIT*EXP(PFMIN(LOCATIONA(J3))-PFTOTALA)/PEMKSUM(LOCATIONA(J3))
!        PRINT *,'J3,LOCATIONA,PFMIN,EMKSUM,KBA=',J3,LOCATIONA(J3),PFMIN(LOCATIONA(J3)),EMKSUM(LOCATIONA(J3)),KBA
      ENDDO

      KAB=0.0D0
      KSSAB=0.0D0
      DO J3=1,NMINB
         COMMIT=0.0D0
         DO J4=1,NCOL(LOCATIONB(J3)) ! sum over renormalised branching probabilities to A minima
            IF (NDISTA(COL_IND(J4+ROW_PTR(LOCATIONB(J3))-1)).EQ.0) THEN
               COMMIT=COMMIT+DVEC(J4+ROW_PTR(LOCATIONB(J3))-1)
!              PRINT *,'J3,J4,PBRANCH,COMMIT=',J3,J4,DVEC(J4+ROW_PTR(LOCATIONB(J3))-1),COMMIT
            ENDIF
         ENDDO
         KAB  =KAB  +COMMIT*EXP(PFMIN(LOCATIONB(J3))-PFTOTALB)/EMKSUM(LOCATIONB(J3))
         KSSAB=KSSAB+COMMIT*EXP(PFMIN(LOCATIONB(J3))-PFTOTALB)/PEMKSUM(LOCATIONB(J3))
!        PRINT *,'J3,LOCATIONB,PFMIN,EMKSUM,KAB=',J3,LOCATIONB(J3),PFMIN(LOCATIONB(J3)),EMKSUM(LOCATIONB(J3)),KAB
      ENDDO

      WRITE(*,'(A,F15.5,A)') 'NGT> Peak array size for branching ratios=',GBMAX,' Gb'
      WRITE(*,'(4(A,G20.10))') 'NGT> For direction A<-B:  kNSS(A<-B)=',KAB,' kNSS(B<-A) (from detailed balance)=',
     &                              KAB*EXP(PFTOTALB-PFTOTALA) ! ,' total=',KAB*(1.0D0+EXP(PFTOTALB-PFTOTALA))
      WRITE(*,'(4(A,G20.10))') 'NGT> For direction B<-A:  kNSS(A<-B) (from detailed balance)=',KBA*EXP(PFTOTALA-PFTOTALB),
     &                            ' kNSS(B<-A)=',KBA ! ,' total=',KBA*(1.0D0+EXP(PFTOTALA-PFTOTALB))
      WRITE(*,'(A,G20.10)') 'NGT> detailed balance for kSS,  ratio should be one if SS applies: ',KSSAB*EXP(PFTOTALB-PFTOTALA)/KSSBA
      WRITE(*,'(A,G20.10)') 'NGT> detailed balance for kNSS, ratio should be one if SS applies: ',KAB*EXP(PFTOTALB-PFTOTALA)/KBA
      WRITE(*,'(4(A,G20.10))') 'NGT> kSS (A<-B)=',KSSAB,'   kSS (B<-A)=',KSSBA
      WRITE(*,'(4(A,G20.10))') 'NGT> kNSS(A<-B)=',KAB,'   kNSS(B<-A)=',KBA

      IF (NGTDISCONNECTALL) THEN
! Reset to rectangular storage, as we are assuming for now that the CR storage will not be necessary when we're just considering 
! A and B minima

! Note that this allocation may be an overestimate in the 1st dimension; NMINA+NMINB is the maximum possible value for each connectivity.
         ALLOCATE(PBRANCH(NMINA+NMINB,NMINA+NMINB),NVAL(NMINA+NMINB,NMINA+NMINB))

         NVAL=0 ! whole array 
         PBRANCH=0.0D0 ! whole array
         DO J3=1,NMINA+NMINB
            DO J4=1,NCOL(J3)
               PBRANCH(J4,J3)=DVEC(ROW_PTR(J3)-1+J4)
               NVAL(J4,J3)=COL_IND(ROW_PTR(J3)-1+J4)
            END DO
         END DO
      ENDIF

      IF (ALLOCATED(DVEC)) DEALLOCATE(DVEC)
      IF (ALLOCATED(COL_IND)) DEALLOCATE(COL_IND)
      IF (ALLOCATED(ROW_PTR)) DEALLOCATE(ROW_PTR)

      RETURN
      END
