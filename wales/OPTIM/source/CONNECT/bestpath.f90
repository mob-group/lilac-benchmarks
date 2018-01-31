!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE BESTPATHDIJKSTRA
     IMPLICIT NONE

     CONTAINS

     SUBROUTINE BESTPATH
          USE PORFUNCS
          USE CONNECTDATA
          USE CONNECTUTILS,ONLY: GETDISTANCE,SETDISTANCE,GETINTERP,NREDO,MAKEMINPATHINFO,MAKETSPATHINFO
          USE KEY, ONLY : INTCONSTRAINTT, INTLJT
          USE KEYDECIDE
          USE KEYCONNECT, ONLY : IMAGEMAX
          IMPLICIT NONE

          ! use Dijkstra's algorithm to find the shortest path between the original two endpoints.
          ! The edge weights are zero for minima connected by a single ts, infinity for pairs
          ! where we have already tried the maximum allowed connections, and the minimum distance from
          ! mind (possibly squared) otherwise.
          ! DIJKSTRA should return the number of connections to try from this shortest path, and the
          ! corresponding pairs of minima.
          ! Implementation of Dijkstra's algorithm for finding shortest paths in an UNconnected database
          ! Weights edges of graph according to mi(i)%data%D(j) or mi(i)%data%INTERP(j) (i>j)
          ! Weight=0 if connected by a single ts
          ! Weight=infinity if the maximum tryconnect attempts have been done for this pair
          ! Weight=mind(i,j) otherwise
INTEGER :: I, NDUMMY, BEGIN, END, NSTEP, J, ISTAT
TYPE MINIMUM
   INTEGER :: PREVMIN                             ! PREVIOUS MIN/TS IN SHORTEST PATH
   DOUBLE PRECISION :: D                          ! DISTANCE FROM STARTING POINT
   LOGICAL :: PERMANENT                           ! ONCE IN PERMANENT SET, THE MINIMUM's
                                                  ! label cannot be changed
END TYPE MINIMUM
TYPE(MINIMUM), DIMENSION(:), ALLOCATABLE :: DIJMIN
TYPE PATHWAY
   INTEGER :: DIJMIN, DIJTS
END TYPE PATHWAY
DOUBLE PRECISION :: WEIGHT, MINWEIGHT, DUMMYDJW, DIST
INTEGER JS, JF, J1, J2, J3, J4, NSP, NPATHTS
INTEGER, DIMENSION(:), ALLOCATABLE :: PATHMINS, PATHTS
LOGICAL, DIMENSION(:), ALLOCATABLE :: PATHSP

975 CONTINUE

! Initialise connections
   ALLOCATE(DIJMIN(NMIN))
   DIJMIN(1:NMIN)%PREVMIN=0
   DIJMIN(1:NMIN)%D=HUGE(1.0D0)
   DIJMIN(1:NMIN)%PERMANENT=.FALSE.

! Start algorithm at 'begin' minimum

   I=1
   DIJMIN(I)%D=0.0D0
   DIJMIN(I)%PERMANENT=.TRUE.
      
   MAIN: DO
      DO J=2, NMIN    ! LOOP OVER ALL MINIMA 
         IF (DIJMIN(J)%PERMANENT) CYCLE          ! ONLY CONSIDER TEMPORARY MINIMA
         IF (GETDISTANCE(I,J)>HUGE(1.0D0)/2) THEN ! DON't raise a HUGE number to any power > 1!!
            WEIGHT=GETDISTANCE(I,J)
         ELSE
            IF (INDEXCOSTFUNCTION) THEN 
               WEIGHT=GETDISTANCE(I,J)
               IF (WEIGHT.GT.0.0D0) THEN
                  WEIGHT=ABS(I-J)
                  IF (I.EQ.2) WEIGHT=ABS((NREDO+1)/2-J)
                  IF (J.EQ.2) WEIGHT=ABS((NREDO+1)/2-I)
               ENDIF
            ELSEIF (INTERPCOSTFUNCTION) THEN 
                 DIST=GETDISTANCE(I,J)
                 WEIGHT=GETINTERP(I,J)/1.0D3+DIST
                 IF (DIST.GT.DIJKSTRADMAX) THEN
                    IF (EXPCOSTFUNCTION) THEN
                       WEIGHT=EXP(MIN(WEIGHT,300.0D0))
                    ELSE
                       WEIGHT=COSTFUNCTIONPOWER*LOG(MAX(WEIGHT,1.0D-10)) 
                       WEIGHT=EXP(MIN(WEIGHT,300.0D0))
                    ENDIF
                 ENDIF
            ELSEIF (EXPCOSTFUNCTION) THEN ! SAVES MEMORY AND CPU WHEN ENDPOINT SEPARATION IS VERY LARGE SAT
                 DIST=GETDISTANCE(I,J)
                 WEIGHT=DIST
                 IF (DIST.GT.DIJKSTRADMAX) WEIGHT=EXP(DIST)
            ELSE ! Compare squares to favour more small jumps over big ones DJW
                 DIST=GETDISTANCE(I,J)
                 WEIGHT=DIST
                 IF (DIST.GT.DIJKSTRADMAX) WEIGHT=DIST**COSTFUNCTIONPOWER
            ENDIF
         ENDIF
         JS=MAX(I,J)
         JF=MIN(I,J)

         IF ((DIJMIN(I)%D + WEIGHT) < DIJMIN(J)%D) THEN    ! UPDATE D FOR EACH CONNECTED MINIMUM
!           PRINT '(A,2I8,3G20.10)','decide> I,J,DIJMIN(J)%D,DIJMIN(I)%D,WEIGHT=',I,J,DIJMIN(J)%D,DIJMIN(I)%D,WEIGHT
            DIJMIN(J)%D=DIJMIN(I)%D + WEIGHT
            DIJMIN(J)%PREVMIN=I
         ENDIF
!        PRINT '(A,2I8,3G20.10)','decide> I,J,DIJMIN(J)%D,DIJMIN(I)%D,WEIGHT=',I,J,DIJMIN(J)%D,DIJMIN(I)%D,WEIGHT
      ENDDO
            
      MINWEIGHT=HUGE(MINWEIGHT)
      NDUMMY=I ! SAT: IF THIS IS THE FIRST CYCLE OF THE MAIN DO LOOP AND ALL MINIMA HAVE
      ! HUGE weight ndummy will never be set which is a bit unfortunate given
      ! that we use it immediately afterwards in the IF test ...
      DO J=2,NMIN                           ! FIND TEMPORARY MINIMUM WITH LOWEST D
!        PRINT '(A,I8,L5,2G20.10,I6)','decide> J,DIJMIN(J)%PERMANENT,DIJMIN(J)%D,MINWEIGHT,NDUMMY=', &
! &            J,DIJMIN(J)%PERMANENT,DIJMIN(J)%D,MINWEIGHT,NDUMMY
         IF (DIJMIN(J)%PERMANENT) CYCLE     ! ONLY CONSIDER TEMPORARY MINIMA
         IF (DIJMIN(J)%D < MINWEIGHT) THEN  ! FIND MINIMUM WEIGHT 
            MINWEIGHT=DIJMIN(J)%D
            NDUMMY=J
         ENDIF
      ENDDO

! If no temporary minima with total weight less than HUGE are present, we will never exit the outer do loop!
      IF (NDUMMY==I) THEN
         WRITE(*,'(A)') ' bestpath> There appear to be no connections left to try in DIJKSTRA'
         IMAGEMAX=MAX(1.1D0*IMAGEMAX,DBLE(IMAGEMAX+1))
         IF (IMAGEMAX.GT.10000) STOP
         WRITE(*,'(A,I6,A)') ' bestpath> Increasing maximum number of images to ',IMAGEMAX,' and clearing weights'
         DO J1=1,NMIN-1
            DO J2=J1+1,NMIN
               IF (MI(J2)%DATA%D(J1).GT.HUGE(1.0D0)/1.0D2) THEN
!
! Perhaps we should calculate the metric again? On the other hand,
! we just need a weight that enables us to continue!
!
                  IF (INTERPCOSTFUNCTION) MI(J2)%DATA%INTERP(J1)=10.0D0
                  MI(J2)%DATA%D(J1)=10.0D0
                  MI(J2)%DATA%INTNTRIES(J1)=1
               ENDIF
            ENDDO
         ENDDO
!        CALL TSUMMARY
         DEALLOCATE(DIJMIN)
         CALL FLUSH(6)
         GOTO 975
!        STOP
      ENDIF

      I=NDUMMY
      DIJMIN(I)%PERMANENT=.TRUE.
      IF (I==2) EXIT   
   ENDDO MAIN

   IF (MINWEIGHT.GT.HUGE(MINWEIGHT)/2) THEN
      WRITE(*,'(A)') ' bestpath> There appear to be no connections left to try in DIJKSTRA'
      IMAGEMAX=MAX(1.1D0*IMAGEMAX,DBLE(IMAGEMAX+1))
      IF (IMAGEMAX.GT.10000) STOP
      WRITE(*,'(A,I6,A)') ' bestpath> Increasing maximum number of images to ',IMAGEMAX,' and clearing weights'
      DO J1=1,NMIN-1
         DO J2=J1+1,NMIN
            IF (MI(J2)%DATA%D(J1).GT.HUGE(1.0D0)/1.0D2) THEN
!
! Perhaps we should calculate the metric again? On the other hand,
! we just need a weight that enables us to continue!
!
               IF (INTERPCOSTFUNCTION) MI(J2)%DATA%INTERP(J1)=10.0D0
               MI(J2)%DATA%D(J1)=10.0D0
               MI(J2)%DATA%INTNTRIES(J1)=1
            ENDIF
         ENDDO
      ENDDO
!     CALL TSUMMARY
      DEALLOCATE(DIJMIN)
      CALL FLUSH(6)
      GOTO 975
!     STOP
   ENDIF

   IF (ABS(MINWEIGHT).LT.1.0D-100) THEN
      WRITE(*,'(A)') ' WARNING - zero weight path already exists in DIJKSTRA'
   ENDIF

! Should have found shortest path now. Get it, backwards from end

  BEGIN=1
  END=2
  I=END
  NSTEP=0
  NDIJPAIRS=0
  DO 
     IF (I==BEGIN) EXIT
     NSTEP=NSTEP + 1
     IF (GETDISTANCE(I,DIJMIN(I)%PREVMIN)>0.0D0) THEN
        NDIJPAIRS=NDIJPAIRS+1
     ENDIF
     I=DIJMIN(I)%PREVMIN
  ENDDO

     ALLOCATE(DIJPAIR(NDIJPAIRS,2))
     ALLOCATE(DIJPAIRDIST(NDIJPAIRS))
!    IF (INDEXCOSTFUNCTION) THEN 
!       MINWEIGHT=MINWEIGHT
!    ELSEIF (EXPCOSTFUNCTION) THEN ! SAT
!         MINWEIGHT=LOG(MINWEIGHT)
!    ELSEIF (INTERPCOSTFUNCTION) THEN 
!         MINWEIGHT=MINWEIGHT
!    ELSE
!         MINWEIGHT=MINWEIGHT**(1.0D0/COSTFUNCTIONPOWER)
!    ENDIF
     WRITE(*,'(A,I6,A,I6,A,G15.5)') ' bestpath> Shortest path in Dijkstra has ',nstep,' steps with ',NDIJPAIRS, &
                                 ' missing connections, weight=',MINWEIGHT
     WRITE(*,'(A)') ' bestpath> The unconnected minima in the chain and their distances are:'
  I=END
  NSTEP=0
  NDIJPAIRS=0
  ALLOCATE(PATHMINS(NMIN))
  ALLOCATE(PATHTS(NTS))
  ALLOCATE(PATHSP(NMIN+NTS))
  NSP=0
  NPATHTS=0
  DO 
     IF (I==BEGIN) EXIT
     NSTEP=NSTEP + 1
     NSP=NSP + 1
! add min to pathlist
     PATHMINS(NSTEP)=I
     PATHSP(NSP)=.TRUE.
     DUMMYDJW=GETDISTANCE(I,DIJMIN(I)%PREVMIN)
!    WRITE(*,'(I6,F12.2)',ADVANCE="NO") I,DUMMYDJW
     IF (DUMMYDJW>0.0D0) THEN
!        PRINT *, I,DUMMYDJW,DIJMIN(I)%PREVMIN
        WRITE(*,'(I6,F12.2,I6)',ADVANCE="NO") I,DUMMYDJW,DIJMIN(I)%PREVMIN
        NDIJPAIRS=NDIJPAIRS+1
        DIJPAIR(NDIJPAIRS,1)=I
        DIJPAIR(NDIJPAIRS,2)=DIJMIN(I)%PREVMIN
        DIJPAIRDIST(NDIJPAIRS)=DUMMYDJW
     ELSE
       NPATHTS=NPATHTS + 1
       NSP=NSP + 1
       PATHSP(NSP)=.FALSE.
       DO J4=1,NTS
          IF (((TS(J4)%data%P).EQ.I).AND.((TS(J4)%data%M).EQ.DIJMIN(I)%PREVMIN)) THEN
             PATHTS(NPATHTS)=J4
          ELSEIF (((TS(J4)%data%M).EQ.I).AND.((TS(J4)%data%P).EQ.DIJMIN(I)%PREVMIN)) THEN
             PATHTS(NPATHTS)=J4
          ENDIF
       ENDDO
     ENDIF
     I=DIJMIN(I)%PREVMIN
  ENDDO
  PATHMINS(NSTEP+1)=I
  PATHSP(NSP+1)=.TRUE.
  NSTEP=NSTEP +1
  NSP=NSP + 1
  NPATHTS=0
  DO J3=1,NSP
     IF (.NOT.(PATHSP(J3))) THEN
        NPATHTS=NPATHTS + 1
        CALL MAKEMINPATHINFO(TS(PATHTS(NPATHTS))%data%P)
        CALL MAKETSPATHINFO(PATHTS(NPATHTS))
        CALL MAKEMINPATHINFO(TS(PATHTS(NPATHTS))%data%M)
     ENDIF
  ENDDO
! WRITE(*,'(I6)') BEGIN
  WRITE(*,'(A)') ' '

  DEALLOCATE(PATHMINS)
  DEALLOCATE(PATHTS)
  DEALLOCATE(PATHSP)
  DEALLOCATE(DIJMIN)
  CALL FLUSH(88)
  CLOSE(88)
     END SUBROUTINE BESTPATH

END MODULE BESTPATHDIJKSTRA 
