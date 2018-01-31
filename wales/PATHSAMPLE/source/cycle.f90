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

SUBROUTINE CYCLE
USE COMMONS
USE PORFUNCS
IMPLICIT NONE
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), RANDOM, DISTANCE, DPRAND, RMAT(3,3), DIST2
INTEGER J3, PID(NCPU), MINS(NCPU), MINF(NCPU), STATUS, NCYCLES, NMINOLD, NTSOLD, ISTAT, NREJ, J1
LOGICAL KILLED(NCPU), NOJOB(NCPU), STOPCYCLE
CHARACTER(LEN=10) CONNSTR

NMINOLD=NMIN; NTSOLD=NTS
CALL PFOLD
OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
CLOSE(1)

cycles: DO NCYCLES=1,NATTEMPT
   IF (DEBUG) PRINT '(A)','cycle> removing previous OPTIM files'
!  CALL MYSYSTEM(STATUS,DEBUG, &
!&   'rm odata.[0-9]* *.xyz.[0-9]* finish.[0-9]* EofS.*[0-9]* OPTIM.*.[0-9]* energies.[0-9]* points.[0-9]*') 
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f odata.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f convcheck.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f path.[0-9]*.xyz.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f EofS.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f OPTIM.connect.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f OPTIM.checksp.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f path.info.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f min.data.info.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f input.crd.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f finish.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f points1.inp.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f points2.inp.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f vector.dump.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f coords.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f path.xyz.*[0-9] > /dev/null ')
   CALL MYSYSTEM(STATUS,DEBUG,'rm -f points.final.*[0-9] > /dev/null ')
!
!  Create NCPU start/finish pairs, which require NCPU OPTIM connect runs.
!  Make the NCONNRUNS odata.connect files, then run them, then process the results.
!
   DO J3=1,NCPU
      NREJ=-1
! 20    CALL RANDOM_NUMBER(RANDOM)
20    RANDOM=DPRAND()
      NREJ=NREJ+1
      MINS(J3)=NINT(0.5D0+NMIN*RANDOM)     ! starting min - change variable name
! 10    CALL RANDOM_NUMBER(RANDOM)
10    RANDOM=DPRAND()
      MINF(J3)=NINT(0.5D0+NMIN*RANDOM)     ! finishing min - change variable name
      IF (MINF(J3).EQ.MINS(J3)) GOTO 10
!
!  Do the check based on GPFOLD difference first so that we save the MIND
!  call if it fails.
!
!  Any difference > PSCALE will be accepted.
!
!     CALL RANDOM_NUMBER(RANDOM)
      RANDOM=DPRAND()
      IF (0.0D0+ABS(GPFOLD(MINS(J3))-GPFOLD(MINF(J3)))/PSCALE.LT.RANDOM) THEN
         IF (NREJ.LT.1000000) THEN
!           IF (DEBUG) PRINT '(A,2I6,A,G20.10)','rejected pair ',MINS(J3),MINF(J3),' difference=', &
! &                     ABS(GPFOLD(MINS(J3))-GPFOLD(MINF(J3)))/PSCALE
            GOTO 20
         ENDIF
      ENDIF
      READ(UMIN,REC=MINS(J3)) SPOINTS(1:3*NATOMS)
      READ(UMIN,REC=MINF(J3)) FPOINTS(1:3*NATOMS)
      CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
!     CALL RANDOM_NUMBER(RANDOM)
      RANDOM=DPRAND()
!
!  Calculate metric for connection attempt probability. We probably want this
!  to depend upon the minimised distance and the difference in committor probabilities.
!
!  Any DISTANCE < DSCALE is accepted; then the probability decreases exponentially.
!
      IF (EXP(-(DISTANCE-DSCALE)/DSCALE).LT.RANDOM) THEN
         IF (NREJ.LT.1000000) THEN
!           IF (DEBUG) PRINT '(A,2I6,A,G20.10,A,G20.10)','rejected pair ',MINS(J3),MINF(J3),' distance=', &
! &                     DISTANCE,' exponent=',EXP(-(DISTANCE-DSCALE)/DSCALE)
            GOTO 20
         ENDIF
      ENDIF

      CALL CONNECTODATA(J3,SPOINTS,FPOINTS)
!
! Printing PFOLD makes no sense, sincw PFOLD isn't being calculated.
!
!     WRITE(*,'(2(A,I6),A,F12.1,A,F12.3,A,I8)') 'cycle> connecting between minima ',MINS(J3),' and ', &
! &      MINF(J3), ' distance=',DISTANCE,' |Pfold diff|=',ABS(GPFOLD(MINS(J3))-GPFOLD(MINF(J3))),' rejects=',NREJ
      WRITE(*,'(2(A,I6),A,F12.1,A,I8)') 'cycle> connecting between minima ',MINS(J3),' and ', &
  &      MINF(J3), ' distance=',DISTANCE,' rejects=',NREJ
      NOJOB(J3)=.FALSE.
!
!  Don;t bother to check whether we already have a single ts that does this connection.
!
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      call fork_subr(PID(J3))
!     IF (DEBUG.AND.(PID(J3).NE.0)) WRITE(*,'(A,I8)') 'forked connect run process id=',PID(J3)
!     IF (PID(J3).NE.0) WRITE(*,'(A,I8)') 'forked connect run process id=',PID(J3)
!     IF (PID(J3).EQ.0) WRITE(*,'(A,I8)') 'I am the child! PID=',PID(J3)
      IF (PID(J3).EQ.0) CALL SUBMITOPTIMJOB(J3,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.connect.')
   ENDDO

   CALL MYWAIT(NCPU,0,PID,NOJOB,KILLED,DEBUG) ! manage the forked processes
   WRITE(*,'(A)') 'cycle> all forked connect runs are completed or killed'

!
!  It is important to identify OPTIM jobs that did not terminate with exit code 0.
!  In such cases KILLED should be .TRUE.
!
   analyse_connect: DO J3=1,NCPU
      WRITE(*,'(A)') ' '
      WRITE(*,'(A,I6,A,I8)') 'cycle> analysing result of search ',J3,' for process id ',PID(J3)
      IF (KILLED(J3)) THEN
         WRITE(*,'(A,I6,A)') 'cycle> connection ',J3,' was unsuccessful'
!
!  Nevertheless, there could be a viable path.info file if DUMPALLPATHS is set in
!  OPTIM. Give it a try?
!
!        IF (TRIPLES) THEN
            PRINT '(A)','cycle> attempting to analyse the path.info file nevertheless'
!        ELSE
!           CYCLE analyse_connect
!        ENDIF
      ENDIF
      WRITE(CONNSTR,'(I10)') PID(J3)
!     CALL MYSYSTEM(STATUS,DEBUG,'cp OPTIM.connect.'//TRIM(ADJUSTL(CONNSTR))//' OPTIM.connect') ! not needed
      CALL MYSYSTEM(STATUS,DEBUG,'cp path.info.'//TRIM(ADJUSTL(CONNSTR))//' path.info')
!     IF (TRIPLES) THEN
         CALL GETALLPATHS
!     ELSE
!        CALL GETNEWPATH(MINS(J3),MINF(J3))
!     ENDIF
   ENDDO analyse_connect

   WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                 '--------------------------------------------------'
   WRITE(*,'(5(A,I8))') 'cycle> end of cycle ',NCYCLES,' new min=',NMIN-NMINOLD,' new ts=',NTS-NTSOLD, &
  &                       ' total min=',NMIN,' total ts=',NTS
   WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                 '--------------------------------------------------'
   NMINOLD=NMIN; NTSOLD=NTS

   IF ((NPAIRDONE.GT.0).AND.(.NOT.DUMMYRUNT)) THEN
      OPEN(UNIT=1,FILE='pairs.data',STATUS='UNKNOWN')
      WRITE(1,'(2I8)') (PAIR1(J1),PAIR2(J1),J1=1,NPAIRDONE)
      CLOSE(1)
   ENDIF
   IF (NMINDONE.GT.0) THEN
      OPEN(UNIT=1,FILE='min.done',STATUS='UNKNOWN')
      WRITE(1,'(2I8)') (MINDONE(J1),J1=1,NMINDONE)
      CLOSE(1)
   ENDIF
   IF (PFOLDINT.NE.0) THEN
      IF (MOD(NCYCLES,PFOLDINT).EQ.0) THEN
         CALL PFOLD
         OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
         WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
         CLOSE(1)
      ENDIF
   ENDIF
   INQUIRE(FILE='STOP',EXIST=STOPCYCLE)
   IF (STOPCYCLE) THEN
      PRINT '(A)','File STOP detected - exit'
      EXIT
   ENDIF

ENDDO cycles

RETURN

END SUBROUTINE CYCLE
