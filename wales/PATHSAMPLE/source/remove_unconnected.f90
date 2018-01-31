!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999- David J. Wales
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
!  Remove stationary points that are not connected to A, B or either A or B, including
!  minima with fewer than NCONMINN connections. Write new database files with renumbered
!  stationary points.
!
      SUBROUTINE REMOVE_UNCONNECTED
      USE PORFUNCS
      USE COMMONS
      IMPLICIT NONE
      INTEGER M1, M2, J1, J2, J3, J4, MINMAP(NMIN), NLEFT, NDUMMY
      INTEGER NCOL(NMIN), NDEAD, NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN, DMAX, NUNCONA, NUNCONB
      LOGICAL DEADTS(NTS), MATCHED, CHANGED  
      INTEGER ISTAT, MINPREV(NMIN)
      DOUBLE PRECISION LKSUM(NMIN)
      DOUBLE PRECISION TNEW, ELAPSED
      INTEGER, ALLOCATABLE :: NVAL(:,:)
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS)

      CALL CPU_TIME(ELAPSED)

      CALL GETNCONN
!
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!
      CALL REGROUP(MINMAP)

      CALL RATECONST_SETUP(LKSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)

      DO J1=1,NMIN
         MINPREV(J1)=0
      ENDDO
!
!  This setup is copied from NGT. DJW 9/3/11
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  Degenerate rearrangements are excluded.
!
!  NCONN may have been set to zero for minima in a disconnected region. DEADTS could still
!  be false for these, so exclude them explicitly.
!  Initial self-connection branching probability is zero.
!
      NCONNMAX=NCONNMAX+1
      PRINT '(A,F4.1,A,I8,A,I6)','remove_unconnected> about to try and allocate ', &
     &     NCONNMAX*1.0D0*NMIN*1.0D0*8*1.0D0/1.0D9,'Gb of RAM for ',NMIN, &
     &                           ' minima, maximum connectivity ',NCONNMAX
      CALL FLUSH(6)
      PRINT *,'allocating NVAL dimensions ',NCONNMAX,NMIN
      ALLOCATE(NVAL(NCONNMAX,NMIN))
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
                     ELSEIF (MINUS(J1).GT.NVAL(NCOL(M2)-1,M2)) THEN
                        NVAL(NCOL(M2),M2)=MINUS(J1)
                     ELSE
                        j3loop: DO J3=1,NCOL(M2)-1
                           IF (MINUS(J1).LT.NVAL(J3,M2)) THEN
!
! Move the rest up.
!
                              DO J4=NCOL(M2),J3+1,-1
                                 NVAL(J4,M2)=NVAL(J4-1,M2)
                              ENDDO 
                              NVAL(J3,M2)=MINUS(J1)
                              EXIT j3loop
                           ENDIF
                        ENDDO j3loop
                     ENDIF
                  ENDIF
               ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2 
                  MATCHCOLM: DO M1=1,NCOL(M2)
                     IF (NVAL(M1,M2).EQ.PLUS(J1)) THEN ! A PREVIOUS TS ALSO LINKS THIS PAIR
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
                     ELSEIF (PLUS(J1).GT.NVAL(NCOL(M2)-1,M2)) THEN
                        NVAL(NCOL(M2),M2)=PLUS(J1)
                     ELSE
                        j3loop2: DO J3=1,NCOL(M2)-1
                           IF (PLUS(J1).LT.NVAL(J3,M2)) THEN
!
! Move the rest up.
!
                              DO J4=NCOL(M2),J3+1,-1
                                 NVAL(J4,M2)=NVAL(J4-1,M2)
                              ENDDO 
                              NVAL(J3,M2)=PLUS(J1)
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
!  Initial self-branching probability, but only for minima with enough connections!
!  Must maintain sorted lists!
!
!  The number of connections in NCONN counts connections between the same minima through
!  different transition states separately. Hence NCOL can end up less than NCONN.
!
!
      DO M2=1,NMIN
         IF (NCONN(M2).GT.NCONNMIN) THEN
            NCOL(M2)=NCOL(M2)+1
            IF (M2.GT.NVAL(NCOL(M2)-1,M2)) THEN
               NVAL(NCOL(M2),M2)=M2
            ELSE
               j3loop3: DO J3=1,NCOL(M2)-1
                  IF (M2.LT.NVAL(J3,M2)) THEN
!
! Move the rest up.
!
                     DO J4=NCOL(M2),J3+1,-1
                        NVAL(J4,M2)=NVAL(J4-1,M2)
                     ENDDO
                     NVAL(J3,M2)=M2
                     EXIT j3loop3
                  ENDIF
               ENDDO j3loop3
            ENDIF
         ENDIF
      ENDDO
!
!  Check connectivity of the stationary point database.
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
      PRINT '(3(A,I8))','remove_unconnected> steps to A region converged in ',NCYCLE-1, &
     &            ' cycles; maximum=',DMAX,' disconnected=',NUNCONA
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
      PRINT '(3(A,I8))','remove_unconnected> steps to B region converged in ',NCYCLE-1, &
     &          ' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!
!  This could happen if disconnected minima lie in the A or B region
!
      IF (NUNCONB.NE.NUNCONA) PRINT '(A)','remove_unconnected> WARNING - number of disconnected minima from A and B is different'
!
!  Remove disconnected minima from consideration.
!
      NLEFT=0
      DO J1=1,NMIN
         IF (TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'A') THEN
            IF (NDISTA(J1).EQ.1000000) THEN
               NCONN(J1)=0
               NCOL(J1)=0
               IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the A region'
            ENDIF
         ELSEIF (TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'B') THEN
            IF (NDISTB(J1).EQ.1000000) THEN
               NCONN(J1)=0
               NCOL(J1)=0
               IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the B region'
            ENDIF
         ELSEIF (TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'AB') THEN
            IF ((NDISTA(J1).EQ.1000000).OR.(NDISTB(J1).EQ.1000000)) THEN
               NCONN(J1)=0
               NCOL(J1)=0
               IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the A or B regions'
            ENDIF
         ELSEIF (TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'AandB'.OR.TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'AANDB') THEN
            IF ((NDISTA(J1).EQ.1000000).AND.(NDISTB(J1).EQ.1000000)) THEN
               NCONN(J1)=0
               NCOL(J1)=0
               IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the A and B regions'
            ENDIF
         ELSEIF (TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'AorB'.OR.TRIM(ADJUSTL(UNCONNECTEDS)).EQ.'AORB') THEN
            IF ((NDISTA(J1).EQ.1000000).OR.(NDISTB(J1).EQ.1000000)) THEN
               NCONN(J1)=0
               NCOL(J1)=0
               IF (DEBUG) PRINT '(A,I8,A)','minimum ',J1,' is disconnected from the A or B regions'
            ENDIF
         ENDIF
         IF (NCONN(J1).GT.NCONNMIN) NLEFT=NLEFT+1
      ENDDO
      PRINT '(A,I8)','remove_unconnected> Number of connected minima remaining=',NLEFT

      NDUMMY=0
      OPEN(UNIT=2,FILE='min.data.removed',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='points.min.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      minloop: DO J1=1,NMIN
         IF (NCONN(J1).LE.NCONNMIN) THEN
            IF (DEBUG) PRINT '(A,I8)','remove_unconnected> removing minimum ',J1
            CYCLE minloop
         ENDIF
         IF (DEBUG) PRINT '(A,I8)','remove_unconnected> not removing minimum ',J1
         NDUMMY=NDUMMY+1
         MINPREV(J1)=NDUMMY
         WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
!
! csw34> Because the minima have been reordered by REGROUP earlier, the
! mapping in MINMAP must be used when extracting coordinates from
! points.min which has NOT been reordered -> MINMAP(J1)
!
         IF (.NOT. NOPOINTS) THEN
            READ(UMIN,REC=MINMAP(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(4,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDIF
      ENDDO minloop
      PRINT '(3(A,I8))','remove_unconnected> ',NDUMMY,' of ',NMIN,' minima left after removing ',NMIN-NDUMMY
      CLOSE(2)
      CLOSE(4)
!
! rewrite min.A and min.B in min.A.removed and min.B.removed since numbering may change.
!
      NDUMMY=0
      Aloop: DO J1=1,NMINA
          IF (NCONN(LOCATIONA(J1)).LE.NCONNMIN) THEN
             PRINT '(A,I8)','remove_unconnected> removing A minimum ',LOCATIONA(J1)
             CYCLE Aloop
          ENDIF
          NDUMMY=NDUMMY+1
      ENDDO Aloop
      OPEN(UNIT=2,FILE='min.A.removed',STATUS='UNKNOWN')
      WRITE(2,'(I8)') NDUMMY
      IF (NDUMMY.EQ.0) THEN
         PRINT '(A)','remove_unconnected> ERROR - all A minima removed'
         STOP
      ENDIF
      DO J1=1,NMINA
         IF (NCONN(LOCATIONA(J1)).GT.NCONNMIN) THEN
            WRITE(2,'(I8)') MINPREV(LOCATIONA(J1))
         ENDIF
      ENDDO 
      CLOSE(2)

      NDUMMY=0
      Bloop: DO J1=1,NMINB
         IF (NCONN(LOCATIONB(J1)).LE.NCONNMIN) THEN
            PRINT '(A,I8)','remove_unconnected> removing B minimum ',LOCATIONB(J1)
            CYCLE Bloop
         ENDIF
         NDUMMY=NDUMMY+1
      ENDDO Bloop
      OPEN(UNIT=2,FILE='min.B.removed',STATUS='UNKNOWN')
      WRITE(2,'(I8)') NDUMMY
      IF (NDUMMY.EQ.0) THEN
         PRINT '(A)','remove_unconnected> ERROR - all B minima removed'
         STOP
      ENDIF
      DO J1=1,NMINB
         IF (NCONN(LOCATIONB(J1)).GT.NCONNMIN) THEN
            WRITE(2,'(I8)') MINPREV(LOCATIONB(J1))
         ENDIF
      ENDDO 
      CLOSE(2)

      OPEN(UNIT=3,FILE='ts.data.removed',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='points.ts.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      NDUMMY=0
      tsloop: DO J1=1,NTS
         IF (MINPREV(PLUS(J1)).EQ.0) THEN
            IF (DEBUG) PRINT '(A,I8,A,I8,A)','remove_unconnected> transition state ',J1,' links minimum ',PLUS(J1),  &
     &                            ' which has been removed - removing this ts'
            CYCLE tsloop
         ENDIF
         IF (MINPREV(MINUS(J1)).EQ.0) THEN
            IF (DEBUG) PRINT '(A,I8,A,I8,A)','remove_unconnected> transition state ',J1,' links minimum ',MINUS(J1),  &
     &                               ' which has been removed - removing this ts'
            CYCLE tsloop
         ENDIF
         NDUMMY=NDUMMY+1
         IF (IMFRQT) THEN
            WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)), &
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
         ELSE
            WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)), &
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
         ENDIF
         IF (.NOT. NOPOINTS) THEN
            READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(5,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDIF
      ENDDO tsloop
      CLOSE(3); CLOSE(5)
      PRINT '(3(A,I8))','remove_unconnected> ',NDUMMY,' of ',NTS,' ts left after removing ',NTS-NDUMMY

      IF (ALLOCATED(NVAL)) DEALLOCATE(NVAL)

      CALL CPU_TIME(TNEW)
      TGT=TGT+TNEW-ELAPSED
      STOP

END SUBROUTINE REMOVE_UNCONNECTED
