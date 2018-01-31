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
!  Subroutine to provide candidate pairs of minima based on distance analysis
!  for pairs of minima.
!
SUBROUTINE CONNECTD(NAVAIL)
USE COMMONS, ONLY: NATOMS, UMIN, CONNECTMIN1, CONNECTMIN2, CONNECTDIST, DMIN1, DMIN2, TCONNECTDIST, TOPPOINTER, POINTERP, &
  &               POINTERM, NMIN, PLUS, MINUS, BULKT, TWOD, DMINMAX, DEBUG, PERMDIST, BOXLX, BOXLY, BOXLZ, RIGIDBODY, &
  &               ANGLEAXIS, PAIR1, PAIR2, NPAIRDONE, INTERPCOSTFUNCTION, CONNECTREGIONT, NOPT, CONNECTMIN2F, NCPU, &
  &               CONNECTMIN2SAVE
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, J1, J2, J3, LOOP2END
INTEGER, ALLOCATABLE :: VINT(:)
DOUBLE PRECISION SPOINTS(NOPT), FPOINTS(NOPT), ELAPSED, TNEW, DISTANCE, RMAT(3,3), DIST2

CALL CPU_TIME(ELAPSED)
ALLOCATE(VINT(DMINMAX))

NAVAIL=0
LOOP2END=CONNECTMIN2F
IF (CONNECTMIN2F.EQ.-1) LOOP2END=NMIN
DO J1=CONNECTMIN1,NMIN
   READ(UMIN,REC=J1) SPOINTS(1:NOPT)
   IF (J1.GE.LOOP2END) LOOP2END=NMIN
   loop2: DO J2=MAX(CONNECTMIN2,J1+1),MIN(LOOP2END,NMIN)

!     PRINT'(A,7I10)','connectd> J1,J2,CONNECTMIN1,CONNECTMIN2,CONNECTMIN2F,LOOP2END,NMIN=', &
!  &                   J1,J2,CONNECTMIN1,CONNECTMIN2,CONNECTMIN2F,LOOP2END,NMIN

!
! For a connectregion run allow connection attempts for local minima
! that already have a direct connection. Otherwise we might miss
! a lower path.
!
!     IF (.NOT.CONNECTREGIONT) THEN
         J3=TOPPOINTER(J2)  !  sets J3 to the ts connected to minimum J2 with the highest value
         IF (J3.LE.0) GOTO 555 
11       IF (PLUS(J3).EQ.J2) THEN
            IF (MINUS(J3).EQ.J1) CYCLE loop2 ! they are already connected
            J3=POINTERP(J3)
         ELSE
            IF (PLUS(J3).EQ.J1) CYCLE loop2 ! they are already connected
            J3=POINTERM(J3)
         ENDIF
         IF (J3.GT.0) GOTO 11 ! check the next ts connected to minimum J2
555      CONTINUE
         DO J3=1,NPAIRDONE
            IF ((PAIR1(J3).EQ.J1).AND.(PAIR2(J3).EQ.J2)) THEN
               IF (DEBUG) PRINT '(3(A,I8),A,G20.10)','connectd> pair ',J1,' and ',J2,' already connected'
               CYCLE loop2 ! do not repeat searches!
            ENDIF
            IF ((PAIR1(J3).EQ.J2).AND.(PAIR2(J3).EQ.J1)) THEN
               IF (DEBUG) PRINT '(3(A,I8),A,G20.10)','connectd> pair ',J1,' and ',J2,' already connected'
               CYCLE loop2 ! do not repeat searches!
            ENDIF
         ENDDO
!     ENDIF
   
         READ(UMIN,REC=J2) FPOINTS(1:NOPT)
     
         CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
      IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2, &
  &                                            RIGIDBODY,RMAT,INTERPCOSTFUNCTION)

!     PRINT '(A,2G20.10)', 'DISTANCE,CONNECTDIST=',DISTANCE,CONNECTDIST
      IF (DISTANCE.LE.CONNECTDIST) THEN

         NAVAIL=NAVAIL+1
         IF (NAVAIL.GT.DMINMAX) THEN
            VINT(1:DMINMAX)=DMIN1(1:DMINMAX)
            DEALLOCATE(DMIN1)
            ALLOCATE(DMIN1(2*DMINMAX))
            DMIN1(1:DMINMAX)=VINT(1:DMINMAX)
      
            VINT(1:DMINMAX)=DMIN2(1:DMINMAX)
            DEALLOCATE(DMIN2)
            ALLOCATE(DMIN2(2*DMINMAX))
            DMIN2(1:DMINMAX)=VINT(1:DMINMAX)
            DMINMAX=2*DMINMAX
            DEALLOCATE(VINT)
            ALLOCATE(VINT(DMINMAX))
         ENDIF
         DMIN1(NAVAIL)=J1
         DMIN2(NAVAIL)=J2
         IF (DEBUG) PRINT '(3(A,I8),A,G20.10)','connectd> connection ',NAVAIL,' pair ',J1,' and ',J2,' distance=',DISTANCE
!        IF (DMINMAX.GT.NCPU) THEN ! set a maximum number of candidates for each call to CONNECTDIST
         IF (NAVAIL.GE.NCPU) THEN ! set a maximum number of candidates for each call to CONNECTDIST
            CONNECTMIN1=J1
            CONNECTMIN2=J2+1
            IF (CONNECTMIN2.GT.MIN(LOOP2END,NMIN)) THEN
               CONNECTMIN1=J1+1
               IF (CONNECTMIN1.EQ.LOOP2END) LOOP2END=NMIN
               CONNECTMIN2=MAX(CONNECTMIN2SAVE,CONNECTMIN1+1)
            ENDIF
            GOTO 777
         ENDIF
      ENDIF
   ENDDO loop2
!  CONNECTMIN2=0 ! this is only to allow an efficient restart
   PRINT '(A,I8,A,I8)','setting CONNECTMIN2 ',CONNECTMIN2,' to CONNECTMIN2SAVE ',CONNECTMIN2SAVE
   CONNECTMIN2=CONNECTMIN2SAVE
ENDDO
! JMC commenting out the line below and the CONNECTMIN2=0 after loop2 (above) so that ALL possible new pairs will be 
! considered next time this double loop is executed. pairs.data will make sure that each pair is only tried once. 
! CONNECTMIN1=NMIN+1
777 CONTINUE
PRINT '(A,I8)','connectd> Number of available connection pairs to try=',NAVAIL

CALL CPU_TIME(TNEW)
TCONNECTDIST=TCONNECTDIST+TNEW-ELAPSED
DEALLOCATE(VINT)

END SUBROUTINE CONNECTD
