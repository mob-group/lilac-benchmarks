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

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Connect minimum MINS to minimum MINF from the current path, call GETNEWPATH to
C  read in the path.
C
      SUBROUTINE CONNECT(OK,MINS,MINF)
      USE COMMONS
      IMPLICIT NONE
      INTEGER J2, STATUS, NDUMMY, MINS, MINF
      LOGICAL OK
      DOUBLE PRECISION DUMMY
      DOUBLE PRECISION, ALLOCATABLE :: VTEMP(:)

      CALL CPU_TIME(ELAPSED)

      CALL MYSYSTEM(STATUS,DEBUG,'sed -e "1,/Connected path found/d" -e "/Rotational/,/#/d" -e "/=/d" -e "/essian/d"
     1             -e "/printing/d" -e "/Number of/d" -e "/STOP/d" OPTIM.connect > pathout')
      OPEN(1,FILE='pathout',STATUS='OLD')
      READ(1,*,END=666)
      WRITE(*,'(A)') '          E+                  Ets                 E-'
      J2=0
      DO 
         J2=J2+1
         IF (J2+1.GT.MAXSTEPS) THEN
            ALLOCATE(VTEMP(J2))
            VTEMP(1:J2)=NEWEMIN(1:J2)
            DEALLOCATE(NEWEMIN)
            ALLOCATE(NEWEMIN(2*MAXSTEPS))
            NEWEMIN(1:J2)=VTEMP(1:J2)
            VTEMP(1:J2-1)=NEWETS(1:J2-1)
            DEALLOCATE(NEWETS)
            ALLOCATE(NEWETS(2*MAXSTEPS))
            NEWETS(1:J2-1)=VTEMP(1:J2-1)
            MAXSTEPS=2*MAXSTEPS 
            DEALLOCATE(VTEMP)
            PRINT '(A,I8)','Maximum number of steps in a path increased to ',MAXSTEPS
         ENDIF
         READ(1,*,END=50) NDUMMY, NEWEMIN(J2), DUMMY, NEWETS(J2), DUMMY, NEWEMIN(J2+1)
         WRITE(*,'(3F20.10)') NEWEMIN(J2), NEWETS(J2), NEWEMIN(J2+1)
      ENDDO
50    NEWLENGTH=J2-1

      IF (DEBUG) PRINT '(A,3I8)','connect> MINS,MINF,NEWLENGTH=',MINS,MINF,NEWLENGTH
      WRITE(*,'(A,I6,A)') 'connect> new segment has length ',NEWLENGTH,' steps'

      CALL CPU_TIME(TNEW)
      TCONNECT=TCONNECT+TNEW-ELAPSED

      CALL GETNEWPATH(MINS,MINF)
      CALL CPU_TIME(ELAPSED)
      OK=.TRUE.
      CALL CPU_TIME(TNEW)
      TCONNECT=TCONNECT+TNEW-ELAPSED
      RETURN

666   CLOSE(1)

      WRITE(*,'(A)') 'connect> unknown connection failure - this should never happen!'
      STOP

      RETURN
      END
