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

! MYWAIT manages forked processes 
! NCPU OPTIM jobs should now be running in forked child processes
!
SUBROUTINE MYWAIT(NCPU,NOFFSET,PID,NOJOB,KILLED,DEBUG)
USE COMMONS, ONLY : SLEEPTIME1, SLEEPTIME2
USE PORFUNCS
IMPLICIT NONE
INTEGER NCPU, NOFFSET, PID(NCPU+NOFFSET), NDONE, J2, PIDDONE, STATUS, NRUNNING,ISTAT
DOUBLE PRECISION TIME0, TIME1
LOGICAL DEBUG, NOJOB(NCPU+NOFFSET), KILLED(NCPU+NOFFSET)
CHARACTER(LEN=80) SLEEPSTRING

WRITE(SLEEPSTRING,'(A,F20.10)') 'sleep ',SLEEPTIME2

KILLED(1:NCPU+NOFFSET)=.FALSE.
CALL CPU_TIME(TIME0)
NDONE=0
NRUNNING=0
DO J2=1+NOFFSET,NCPU+NOFFSET
   IF (.NOT.NOJOB(J2)) NRUNNING=NRUNNING+1
ENDDO
IF (NRUNNING.EQ.0) RETURN
waitloop: DO
!sf344>   piddone = 0
! csw34> Added call to wait here to allow processes to finish on sinister
!  CALL MYSYSTEM(STATUS,DEBUG,' sleep 1')
   IF (SLEEPTIME2.GT.0.0D0) CALL SYSTEM(TRIM(ADJUSTL(SLEEPSTRING)))
   CALL WAIT_SUBR(PIDDONE,STATUS)
11 CONTINUE
!  PRINT '(A,4I10)','PIDDONE,STATUS,NDONE,PID=',PIDDONE,STATUS,NDONE,PID(1+NOFFSET:NCPU+NOFFSET)
   CALL FLUSH(6)
   IF (PIDDONE.GT.0) THEN
      CALL CPU_TIME(TIME1)
      IF (DEBUG) PRINT '(A,I8,A,I6,A,F15.2)','mywait> PID ',PIDDONE,' has finished with exit status ',STATUS, & 
                                             ' local CPU time taken=',TIME1-TIME0
!     WRITE(*,*) 'NDONE, NRUNNING', NDONE, NRUNNING
      CALL FLUSH(6)
      NDONE=NDONE+1
      DO J2=1+NOFFSET,NCPU+NOFFSET
         IF (PIDDONE.EQ.PID(J2)) THEN
            IF (STATUS.NE.0) KILLED(J2)=.TRUE. ! this is probably useless !
            IF (NDONE.EQ.NRUNNING) EXIT waitloop
            CYCLE waitloop
         ENDIF
      ENDDO
      PRINT*,'ERROR - PID of completed child process not recognised: ',PIDDONE
      STOP
   ELSE
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      PRINT '(A,I20)','mywait> WARNING - WAIT returned system error code ',-PIDDONE
!
! Try calling wait again to see if this fixes things.
! For very short OPTIM jobs WAIT may have trouble keeping up!
!
      CALL MYSYSTEM(STATUS,DEBUG,' sleep 1')
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I8))','mywait> on calling wait again pid=',PIDDONE,' status=',STATUS
      PRINT *,'mywait> after second call to WAIT_SUBR, PIDDONE=',PIDDONE
      IF (PIDDONE.GT.0) GOTO 11
      STOP
   ENDIF
ENDDO waitloop

RETURN
END SUBROUTINE MYWAIT
