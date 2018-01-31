SUBROUTINE GETMETRIC(NSTART,NFINISH) ! (NEWPOINTSMIN)
USE PORFUNCS
USE COMMONS, ONLY : NMIN, PAIRDIST, PAIRLIST, PAIRDISTMAX, UMIN, NATOMS, DEBUG, BOXLX, BOXLY, BOXLZ, &
  &                BULKT, TWOD, RIGIDBODY, INTERPCOSTFUNCTION, ETS, PLUS, MINUS, NTS, PAIR1, PAIR2, &
  &                NPAIRDONE, INDEXCOSTFUNCTION, RANDOMMETRICT, NRANDOMMETRIC, GEOMDIFFTOL, ALLPAIRS, &
  &                INITIALDIST, DISBOUND
IMPLICIT NONE
INTEGER J1, J2, J3, J4, J5, ISTAT, BASIN(NMIN), NBASIN, NSTART, NFINISH, J6, NDONE
INTEGER JM, JN, NPOSITION
DOUBLE PRECISION LOCALPOINTS(3*NATOMS), NEWPOINTSMIN(3*NATOMS), DISTANCE, DIST2, RMAT(3,3), DPRAND
DOUBLE PRECISION HIGHESTTS, ETHRESH
LOGICAL CHANGED


!
! Find highest transition state.
!
! HIGHESTTS=-1.0D100
! DO J1=1,NTS
!    IF (ETS(J1).GT.HIGHESTTS) THEN
!       HIGHESTTS=ETS(J1)
!    ENDIF
! ENDDO 
! PRINT '(A,G20.10)','getmetric> highest transition state lies at ',HIGHESTTS
! ETHRESH=HIGHESTTS+1.0D0

! BASIN(1:NMIN)=0
! NBASIN=0
! DO
!    CHANGED=.FALSE.
!    DO J1=1,NTS
!       IF (ETS(J1).LT.ETHRESH) THEN
!          IF ((BASIN(PLUS(J1)).EQ.0).AND.(BASIN(MINUS(J1)).EQ.0)) THEN
!             CHANGED=.TRUE.
!             NBASIN=NBASIN+1
!             BASIN(PLUS(J1))=NBASIN
!             BASIN(MINUS(J1))=NBASIN
!          ELSEIF (BASIN(PLUS(J1)).NE.BASIN(MINUS(J1))) THEN
!             CHANGED=.TRUE.
!             IF (BASIN(PLUS(J1)).EQ.0) THEN
!                BASIN(PLUS(J1))=BASIN(MINUS(J1))
!             ELSEIF (BASIN(MINUS(J1)).EQ.0) THEN
!                BASIN(MINUS(J1))=BASIN(PLUS(J1))
!             ELSE
!                BASIN(PLUS(J1))=MIN(BASIN(PLUS(J1)),BASIN(MINUS(J1)))
!                BASIN(MINUS(J1))=BASIN(PLUS(J1))
!             ENDIF
!          ENDIF
!       ENDIF
!    ENDDO
!    IF (.NOT.CHANGED) EXIT
! ENDDO
! PRINT '(A,I8)','getmetric> Number of superbasins=',NBASIN

DO J6=NSTART,NFINISH
   READ(UMIN,REC=J6) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)

   min2: DO J3=1,NMIN
      DISTANCE=1.0D100
      DO J4=1,NTS
         IF ((PLUS(J4).EQ.J6).AND.(MINUS(J4).EQ.J3)) DISTANCE=0.0D0
         IF ((PLUS(J4).EQ.J3).AND.(MINUS(J4).EQ.J6)) DISTANCE=0.0D0
!        PRINT '(A,5I6,G20.10)','J6,J3,J4,plus,minus,distance=',J6,J3,J4,plus(J4),minus(J4),distance
         IF (DISTANCE.LT.1.0D-10) EXIT
      ENDDO 
!
! This line is setting the metric to 1.0D100 if two minima are connected by a
! discrete path of any length. Change to using the actual metric if they
! are not directly connected.
!
! If they are in the same superbasin, there is a path between them!
!     IF ((BASIN(J3).EQ.BASIN(J6)).AND.(DISTANCE.GT.1.0D-10)) CYCLE 
      IF ((J3.LE.J6).AND.(J3.GE.NSTART)) CYCLE ! already done by symmetry
!
! Set the pairs for which connections have already been tried to infinite distance,
! so they are not tried again. Don't overwrite zero distance settings for connections
! that have actually been found!
! 
      IF ((DISTANCE.GT.1.0D-10).AND.(.NOT.RANDOMMETRICT)) THEN
         DO J4=1,NPAIRDONE
            IF ((PAIR1(J4).EQ.J6).AND.(PAIR2(J4).EQ.J3)) CYCLE min2
            IF ((PAIR1(J4).EQ.J3).AND.(PAIR2(J4).EQ.J6)) CYCLE min2
         ENDDO 
         IF (INDEXCOSTFUNCTION) THEN
            DISTANCE=ABS(J6-J3) 
         ELSEIF (INITIALDIST) THEN
            IF (DISBOUND.LE.0.0D0) THEN ! do all distances
               READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
  &                             DIST2,RIGIDBODY,RMAT,.FALSE.)
            ELSE
               IF (DISTANCE.GT.1.0D-10) DISTANCE=-DISBOUND
            ENDIF
         ELSE
            READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
  &                          DIST2,RIGIDBODY,RMAT,.FALSE.)
            IF (INTERPCOSTFUNCTION) THEN
               CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
  &                             DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
            ENDIF
!    
! The interpolation metric could be zero for minima that are not actually connected.
! Set a minimum non-zero value to avoid Dijinit thinking that they are connected.
! 
            
            DISTANCE=MAX(DISTANCE,GEOMDIFFTOL/100.0)
            !DISTANCE=MAX(DISTANCE,0.1D0)
         ENDIF
      ENDIF
      IF (INITIALDIST.AND.(J6.NE.J3)) THEN
         JM=MIN(J6,J3)
         JN=MAX(J6,J3)
         NPOSITION=((JN-2)*(JN-1))/2+JM 
         ALLPAIRS(NPOSITION)=DISTANCE
      ELSE
!
! Maintain sorted list of nearest nodes according to the chosen interpolation metric.
! 
         sortloop: DO J4=1,PAIRDISTMAX
            IF (DISTANCE.LT.PAIRDIST(J6,J4)) THEN
               DO J5=PAIRDISTMAX,J4+1,-1
                  PAIRDIST(J6,J5)=PAIRDIST(J6,J5-1)
                  PAIRLIST(J6,J5)=PAIRLIST(J6,J5-1)
               ENDDO
               PAIRDIST(J6,J4)=DISTANCE
               PAIRLIST(J6,J4)=J3
               EXIT sortloop
            ENDIF
         ENDDO sortloop
         sortloop2: DO J4=1,PAIRDISTMAX
            IF (DISTANCE.LT.PAIRDIST(J3,J4)) THEN
               DO J5=PAIRDISTMAX,J4+1,-1
                  PAIRDIST(J3,J5)=PAIRDIST(J3,J5-1)
                  PAIRLIST(J3,J5)=PAIRLIST(J3,J5-1)
               ENDDO
               PAIRDIST(J3,J4)=DISTANCE
               PAIRLIST(J3,J4)=J6
               EXIT sortloop2
            ENDIF
         ENDDO sortloop2
      ENDIF
   ENDDO min2
   PRINT '(A,I8)','getmetric> Finished metric calculation for minimum ',J6
   IF (DEBUG) THEN
      IF (INITIALDIST) THEN
!        IF (J6.GT.1) PRINT '(10G13.5)',ALLPAIRS(((J6-1)*(J6-2))/2+1:(J6*(J6-1)/2))
      ELSE
         PRINT '(10G13.5)',PAIRDIST(J6,1:PAIRDISTMAX)
         PRINT '(10I13)',PAIRLIST(J6,1:PAIRDISTMAX)
      ENDIF
   ENDIF
   CALL FLUSH(6)
ENDDO

IF (INITIALDIST) RETURN
IF (.NOT.RANDOMMETRICT) RETURN
!
! We have assigned all the direct connections in the above loop. Now calculate
! NRANDOMMETRIC values for each specified minimum.
!
DO J6=NSTART,NFINISH
   READ(UMIN,REC=J6) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)

   NDONE=0
   ranmin2: DO WHILE (NDONE.LE.NRANDOMMETRIC)
      J3=MIN(INT(NMIN*DPRAND())+1,NMIN)
      DISTANCE=1.0D100
      DO J4=1,NTS
         IF ((PLUS(J4).EQ.J6).AND.(MINUS(J4).EQ.J3)) DISTANCE=0.0D0
         IF ((PLUS(J4).EQ.J3).AND.(MINUS(J4).EQ.J6)) DISTANCE=0.0D0
         IF (DISTANCE.LT.1.0D-10) CYCLE ranmin2 ! already in the list
      ENDDO 
!
! Change to set the metric if there is a >1 step path between these minima.
!
!     IF ((BASIN(J3).EQ.BASIN(J6))) CYCLE ranmin2 ! if they are in the same superbasin, there is a path between them!
!
! Don;t retry the pairs for which connections have already been attempted,
! or repeat calculations for exiting entries in pairlist.
! 
      DO J4=1,NPAIRDONE
         IF ((PAIR1(J4).EQ.J6).AND.(PAIR2(J4).EQ.J3)) CYCLE ranmin2
         IF ((PAIR1(J4).EQ.J3).AND.(PAIR2(J4).EQ.J6)) CYCLE ranmin2
      ENDDO 
      DO J4=1,PAIRDISTMAX
         IF (PAIRLIST(J6,J4).EQ.J3) CYCLE ranmin2
      ENDDO 

      NDONE=NDONE+1

      IF (INDEXCOSTFUNCTION) THEN
         DISTANCE=ABS(J6-J3) 
      ELSE
         READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
  &                       DIST2,RIGIDBODY,RMAT,.FALSE.)
         IF (INTERPCOSTFUNCTION) THEN
            CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
  &                          DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
         ENDIF
!    
! The interpolation metric could be zero for minima that are not actually connected.
! Set a minimum non-zero value to avoid Dijinit thinking that they are connected.
! 

         DISTANCE=MAX(DISTANCE,GEOMDIFFTOL/100.0)
         !DISTANCE=MAX(DISTANCE,0.1D0)
      ENDIF
!
! Maintain sorted list of nearest nodes according to the chosen interpolation metric.
! 
      ransortloop: DO J4=1,PAIRDISTMAX
         IF (DISTANCE.LT.PAIRDIST(J6,J4)) THEN
            DO J5=PAIRDISTMAX,J4+1,-1
               PAIRDIST(J6,J5)=PAIRDIST(J6,J5-1)
               PAIRLIST(J6,J5)=PAIRLIST(J6,J5-1)
            ENDDO
            PAIRDIST(J6,J4)=DISTANCE
            PAIRLIST(J6,J4)=J3
            EXIT ransortloop
         ENDIF
      ENDDO ransortloop
      ransortloop2: DO J4=1,PAIRDISTMAX
         IF (DISTANCE.LT.PAIRDIST(J3,J4)) THEN
            DO J5=PAIRDISTMAX,J4+1,-1
               PAIRDIST(J3,J5)=PAIRDIST(J3,J5-1)
               PAIRLIST(J3,J5)=PAIRLIST(J3,J5-1)
            ENDDO
            PAIRDIST(J3,J4)=DISTANCE
            PAIRLIST(J3,J4)=J6
            EXIT ransortloop2
         ENDIF
      ENDDO ransortloop2
   ENDDO ranmin2
   PRINT '(A,I8)','getmetric> Finished metric calculation with random pair selection for minimum ',J6
   IF (DEBUG) THEN
      PRINT '(10G13.2)',PAIRDIST(J6,1:PAIRDISTMAX)
      PRINT '(10I13)',PAIRLIST(J6,1:PAIRDISTMAX)
   ENDIF
   CALL FLUSH(6)
ENDDO

END SUBROUTINE GETMETRIC
