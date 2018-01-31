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
!  Subroutine to provide candidate pairs of minima based on Dijkstra analysis
!  of the current database.
!
SUBROUTINE GETDPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: UMIN, NATOMS, DMIN1, DMIN2, DIJINITT, NCPU, NPAIRFRQ, NATTEMPT, NMIN, DIJINITFLYT, &
  &  PAIR1, PAIR2, NPAIRDONE, MAXPAIRS, NOPT
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, MINMAP(NMIN), NMINSAVE, J1
DOUBLE PRECISION SPOINTS(NOPT), FPOINTS(NOPT)

10 CONTINUE
IF (NAVAIL.EQ.0) THEN
   IF (DIJINITT) THEN 
      CALL DIJINIT(NAVAIL)
   ELSE IF (DIJINITFLYT) THEN 
      CALL GETNCONN ! must call this first to set NCONNMAX 
      CALL DIJINITFLY(NAVAIL)
   ELSE
      PAIRSTODO=NCPU*NPAIRFRQ
      IF (NPAIRFRQ.LT.1) PAIRSTODO=NATTEMPT*NCPU ! just one set of pairs unless we run out
      IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
      ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO))
      CALL GETNCONN
!
!  NMINSAVE and MINMAP are just dummies here.
!
      NMINSAVE=NMIN
      DO J1=1,NMIN
         MINMAP(J1)=J1
      ENDDO
      CALL DIJKSTRA(NAVAIL,.TRUE.,PAIRSTODO,NMINSAVE,MINMAP)
   ENDIF
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getdpair> No more candidate pairs of minima in getdpair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF
NUSED=NUSED+1
NAVAIL=NAVAIL-1
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)
WRITE(*,'(4(A,I8))') 'getdpair> connecting minima ',MINS,' and ',MINF, ' pairs used=',NUSED,' remaining=',NAVAIL
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
READ(UMIN,REC=MINS) SPOINTS(1:NOPT)
READ(UMIN,REC=MINF) FPOINTS(1:NOPT)

END SUBROUTINE GETDPAIR
