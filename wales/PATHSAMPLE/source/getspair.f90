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
!  Subroutine to provide candidate pairs of minima based on minimum distances
!  for pairs separated by a certain number of steps in the best path.
!
SUBROUTINE GETSPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: UMIN, NATOMS, DMIN1, DMIN2, NATTEMPT, NCPU, MINSEP, BULKT, TWOD, ZSYM, DEBUG, BESTPATHLENGTH, ETS, &
  &               NPAIRFRQ, PAIR1, PAIR2, NPAIRFRQ, NPAIRDONE, MAXPAIRS, PERMDIST, BOXLX, BOXLY, BOXLZ, RIGIDBODY, BESTPATH, &
  &               BARRIERSHORT, EMIN, PLUS, MINUS, KPLUS, KMINUS, RATESHORT, ANGLEAXIS, NMIN, INTERPCOSTFUNCTION, NOPT
USE PORFUNCS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, J1, NSTEPS, J2, J3, J4, N1, N2
INTEGER NMINSAVE, MINMAP(NMIN)
INTEGER, ALLOCATABLE :: MINLIST(:), POSITION(:)
DOUBLE PRECISION, ALLOCATABLE :: DISTLIST(:), BHEIGHT(:)
DOUBLE PRECISION SPOINTS(NOPT), FPOINTS(NOPT), DISTANCE, RMAT(3,3), DIST2

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
   IF (BARRIERSHORT.OR.RATESHORT) THEN
!
! Try connection pairs on either side of the highest barriers on the path, separated
! by MINSEP steps or fewer.
!
      ALLOCATE(BHEIGHT(NSTEPS-1),POSITION(NSTEPS-1))
      DO J1=1,NSTEPS-1
         IF (BARRIERSHORT) THEN
            BHEIGHT(J1)=ETS(BESTPATH(2*J1))-EMIN(BESTPATH(2*J1+1))
         ELSEIF (RATESHORT) THEN
            IF (PLUS(BESTPATH(2*J1)).EQ.BESTPATH(2*J1+1)) THEN
               BHEIGHT(J1)=-KPLUS(BESTPATH(2*J1))
            ELSEIF (MINUS(BESTPATH(2*J1)).EQ.BESTPATH(2*J1+1)) THEN
               BHEIGHT(J1)=-KMINUS(BESTPATH(2*J1))
            ELSE
               PRINT '(3(A,I6))','getspair> ERROR - ts=',BESTPATH(2*J1),' plus=',PLUS(BESTPATH(2*J1)), &
  &                              ' minus=',MINUS(BESTPATH(2*J1))
               PRINT '(A,I6)','getspair> ERROR - minimum id = ',BESTPATH(2*J1+1)
            ENDIF
         ENDIF
         POSITION(J1)=2*J1 ! position in BESTPATH list
      ENDDO
      CALL SORT(NSTEPS-1,NSTEPS-1,BHEIGHT,POSITION)
      IF (BARRIERSHORT) THEN
         PRINT '(A)','sorted barriers on best path labelled according to the transition state:'
         PRINT '(I6,G20.10)',(BESTPATH(POSITION(J1)),BHEIGHT(J1),J1=1,NSTEPS-1)
      ELSE IF (RATESHORT) THEN
         PRINT '(A)','sorted rates on best path labelled according to the transition state:'
         PRINT '(I6,G20.10)',(BESTPATH(POSITION(J1)),-BHEIGHT(J1),J1=1,NSTEPS-1)
      ENDIF
      DO J1=1,NSTEPS-1 ! loop over barriers/rates
!
! Do TS number POSITION(J1) in BESTPATH first etc.
! Must check that we aren't off the end of the best path, and that we
! haven't done this pair before.
! Exit when we have PAIRSTODO pairs.
!
         loop1: DO J2=1,MINSEP ! allow up to MINSEP steps away from the ts in question
            N1=POSITION(J1)-(2*J2-1)
            N2=POSITION(J1)+(2*J2-1)
!           IF ((N1.LT.1).OR.(N2.GT.BESTPATHLENGTH)) EXIT      ! do not go off the end of the path!
!           N1=BESTPATH(N1)
!           N2=BESTPATH(N2)
            IF ((N1.LT.1).AND.(N2.GT.BESTPATHLENGTH)) EXIT
            N1=BESTPATH(MAX(N1,1))
            N2=BESTPATH(MIN(N2,BESTPATHLENGTH))
            DO J3=1,NPAIRDONE
               IF ((PAIR1(J3).EQ.N1).AND.(PAIR2(J3).EQ.N2)) CYCLE loop1 ! do not repeat searches!
               IF ((PAIR1(J3).EQ.N2).AND.(PAIR2(J3).EQ.N1)) CYCLE loop1 ! do not repeat searches!
            ENDDO
            NAVAIL=NAVAIL+1
            DMIN1(NAVAIL)=N1
            DMIN2(NAVAIL)=N2
            DISTLIST(NAVAIL)=BHEIGHT(J1)
            IF (NAVAIL.GE.PAIRSTODO) EXIT
         ENDDO loop1
         IF (NAVAIL.GE.PAIRSTODO) EXIT
      ENDDO
      DEALLOCATE(BHEIGHT,POSITION)
   ELSE
      DISTLIST(1:PAIRSTODO)=1.0D100
      DO J1=1,NSTEPS
         READ(UMIN,REC=MINLIST(J1)) (SPOINTS(J2),J2=1,NOPT)
         min2: DO J2=J1+1,NSTEPS
            DO J3=1,NPAIRDONE
               IF ((PAIR1(J3).EQ.MINLIST(J1)).AND.(PAIR2(J3).EQ.MINLIST(J2))) CYCLE min2 ! do not repeat searches
               IF ((PAIR1(J3).EQ.MINLIST(J2)).AND.(PAIR2(J3).EQ.MINLIST(J1))) CYCLE min2 ! do not repeat searches
            ENDDO
            IF (J2-J1.GE.MINSEP) THEN ! find distance if separation is >= MINSEP
               READ(UMIN,REC=MINLIST(J2)) (FPOINTS(J3),J3=1,NOPT)
               CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                             RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
  &                             DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               NAVAIL=NAVAIL+1
               sortloop: DO J3=1,MIN(NAVAIL,PAIRSTODO) ! sort the shortest PAIRSTODO values
                  IF (DISTANCE.LT.DISTLIST(J3)) THEN
                     DO J4=MIN(NAVAIL,PAIRSTODO),J3+1,-1
                        DMIN1(J4)=DMIN1(J4-1)
                        DMIN2(J4)=DMIN2(J4-1)
                        DISTLIST(J4)=DISTLIST(J4-1)
                     ENDDO
                     DMIN1(J3)=MINLIST(J1)
                     DMIN2(J3)=MINLIST(J2)
                     DISTLIST(J3)=DISTANCE
                     EXIT sortloop
                  ENDIF
               ENDDO sortloop
               IF (DEBUG) PRINT '(3(A,I8),A,G20.10)','getspair> connection ',NAVAIL,' pair ',MINLIST(J1),  &
     &                                               ' and ',MINLIST(J2),' distance=',DISTANCE
            ENDIF
         ENDDO min2
      ENDDO
   ENDIF
   NAVAIL=MIN(NAVAIL,PAIRSTODO) 
   PRINT '(A,I8,A)','getspair> sorted list of ',NAVAIL,' pairs with corresponding barrier height or minimum distance'
   PRINT '(2I8,F15.5)',(DMIN1(J1),DMIN2(J1),DISTLIST(J1),J1=1,NAVAIL)
   DEALLOCATE(MINLIST,DISTLIST)
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getspair> No more candidate pairs of minima in getspair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF
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
READ(UMIN,REC=MINS) SPOINTS(1:NOPT)
READ(UMIN,REC=MINF) FPOINTS(1:NOPT)

END SUBROUTINE GETSPAIR
