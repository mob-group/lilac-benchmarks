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

SUBROUTINE CYCLE2
USE COMMONS
USE PORFUNCS
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), DPRAND, MINFRQ
INTEGER J1, J2, J3, PID(NCPU), MINS(NCPU), MINF(NCPU), STATUS, NCYCLES, NMINOLD, NTSOLD, ISTAT, &
  & NEWJOB, PIDDONE, MINID(NCPU), LUNIT, NSPCHECK, XYZUNIT, &
  & ADDMINXYZCHECK, CONNDOING, CONNATTEMPTS, L1, L2, CONNECTPAIRDONE
INTEGER NAVAIL, NUSED, NDUMMY, NCONNTOTAL
LOGICAL KILLED, STOPCYCLE, CHANGEPD, OPTEST, NOMORE, LTEST1, LTEST2, MINFRQSDUMP, TSFRQSDUMP
CHARACTER(LEN=10) CONNSTR
CHARACTER(LEN=5) SDUMMY
CHARACTER(LEN=80) SLEEPSTRING1, SLEEPSTRING2

IF (.NOT.ALLOCATED(FROZEN)) ALLOCATE(FROZEN(NATOMS))
WRITE(SLEEPSTRING1,'(A,F20.10)') 'sleep ',SLEEPTIME1
WRITE(SLEEPSTRING2,'(A,F20.10)') 'sleep ',SLEEPTIME2
NAVAIL=0
IF (CHECKSPT) NSPCHECK=CHECKSPS-1
IF (ADDMINXYZT) THEN
   XYZUNIT=GETUNIT()
   OPEN(XYZUNIT,FILE=TRIM(ADJUSTL(ADDMINXYZNAME)),STATUS='OLD')
ENDIF
ADDMINXYZCHECK=0
CONNECTPAIRDONE=0
NMINOLD=NMIN; NTSOLD=NTS
IF ((NPFOLD.GT.0).AND.(.NOT.DIJPAIRT)) CALL PFOLD
OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
CLOSE(1)
IF (NEWCONNECTIONST) THEN
   CONNDOING=CONNMINSTART
   CONNATTEMPTS=0
   DO J1=1,NMIN
      NCONNTOTAL=0
      DO L1=1,NTS
         LTEST1=PLUS(L1).EQ.J1
         LTEST2=MINUS(L1).EQ.J1
         IF ((LTEST1.OR.LTEST2).AND.(PLUS(L1).NE.MINUS(L1))) THEN
            NCONNTOTAL=NCONNTOTAL+1
         ENDIF
      ENDDO
!
! MINCONN is allocated in main.F and doubled by mindouble routine if NEWCONNECTIONST is true.
!
      MINCONN(J1)=NCONNTOTAL
   ENDDO
ENDIF
NOMORE=.FALSE.

IF (DEBUG) PRINT '(A)','cycle2> removing previous OPTIM files'
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
! Calculate initial rates if required.
!
IF (RATESCYCLET) CALL GETRATESCYCLE

! NUSED is otherwise uninitialised
NUSED=0
IF ((NMIN.EQ.2).AND.(DIJINITT.OR.DIJINITFLYT)) THEN ! no point running more than one initial search in this case
   CALL GETDPAIR(NAVAIL,NUSED,MINS(1),MINF(1),SPOINTS,FPOINTS)
   CALL CONNECTODATA(1,SPOINTS,FPOINTS)
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   CALL FORK_SUBR(PID(1)) ! PID is zero in the child, non-zero in the parent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IF (PID(1).NE.0) WRITE(*,'(A,I10)') 'cycle2> forked connect run process id=',PID(1)
!  PRINT *,'PID(1)=',PID(1)
!  IF (PID(1).EQ.0) THEN
!     WRITE(*,'(A,I8)') 'cycle2> I am the child! PID=',PID(1)
!     CALL GETPID_SUBR(CHILDPID)
!     PRINT *,'CHILDPID=',CHILDPID
!  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (PID(1).EQ.0) CALL SUBMITOPTIMJOB(1,CHARMMT,UNRST,1,EXEC,DEBUG,'OPTIM.connect.')
   CALL WAIT_SUBR(PIDDONE,STATUS)
!
! Failure to return the right pid here is not fatal if we have the path.info file.
! It shouldn;t happen, though! Calling wait a second time seems to allow the system
! to catch up!
!
   IF (PIDDONE.NE.PID(1)) THEN
      PRINT '(2(A,I10))','cycle2> ERROR - initial path WAIT returned process id',PIDDONE,' instead of ',PID(1)
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I10))','cycle2> on calling wait again pid=',PIDDONE,' status=',STATUS
      IF (PIDDONE.NE.PID(1)) STOP
   ENDIF
   WRITE(CONNSTR,'(I10)') PID(1)
   WRITE(*,'(A,I8)') 'cycle2> analysing result of initial search for process id ',PID(1)
   IF (BHINTERPT.OR.BISECTT) THEN
      CALL MYSYSTEM(STATUS,DEBUG,'cp min.data.info.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info')
   ELSE
      CALL MYSYSTEM(STATUS,DEBUG,'cp path.info.'//TRIM(ADJUSTL(CONNSTR))//' path.info')
   ENDIF
   IF (STATUS.EQ.0) THEN ! The file exists, so try to analyse it
      IF (BHINTERPT.OR.BISECTT) THEN
         CALL GETALLMIN(MINS(1),MINF(1))
      ELSE ! IF (TRIPLES) THEN
         CALL GETALLPATHS
!     ELSE
!        CALL GETNEWPATH(MINS(1),MINF(1))
      ENDIF
   ELSE
      PRINT '(A)','cycle2> ERROR - no path.info file generated by initial search'
      STOP
   ENDIF
ENDIF

DO J3=1,NCPU
   IF (SLEEPTIME1.GT.0.0D0) CALL SYSTEM(TRIM(ADJUSTL(SLEEPSTRING1))) ! to prevent us running out of source ports. Needed for M13.
   IF (DIJPAIRT) THEN
      CALL GETDPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (CONNECTREGIONT) THEN
      CALL GETRPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (CONUNCONT) THEN
      CALL GETCUPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
!  ELSEIF (PTAUT) THEN
!     CALL GETPPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (SHORTCUTT) THEN
      CALL GETSPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (UNTRAPT) THEN
      CALL GETUPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (FREEPAIRT) THEN
      CALL GETFREEPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (CONNECTPAIRST) THEN
      CONNECTPAIRDONE=CONNECTPAIRDONE+1
      IF (CONNECTPAIRDONE.LE.NCONNECTPAIRS) THEN
         IF (CONNECTPAIRDONE.GE.NCONNECTPAIRS) THEN
            NOMORE=.TRUE.
            PRINT '(A)','cycle2> Completed sweep through all available pairs of minima'
         ELSE
            CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in GETCONNECTPAIR
            CALL GETCONNECTPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
         ENDIF
      ENDIF
   ELSEIF (USEPAIRST) THEN
      CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in GETUSEPAIR
      CALL GETUSEPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (ADDMINXYZT) THEN
      ADDMINXYZCHECK=ADDMINXYZCHECK+1
      IF (ADDMINXYZCHECK.LE.NMINADDXYZ) THEN
         READ(XYZUNIT,'(I6)') NDUMMY
         READ(XYZUNIT,'(A)') SDUMMY
         READ(XYZUNIT,*) (SDUMMY,SPOINTS(3*(J1-1)+1:3*(J1-1)+3),J1=1,NATOMS)
      ENDIF
      IF (DEBUG) PRINT '(A,I6)','cycle2> points read for xyz minmum ',ADDMINXYZCHECK
   ELSEIF (NEWCONNECTIONST) THEN
864   IF (CONNDOING.LE.NMIN) THEN
         IF ((MINCONN(CONNDOING).GE.CONNECTIONS)) THEN
            CONNDOING=CONNDOING+1
            CONNATTEMPTS=0
            IF (CONNDOING.GT.NMIN) THEN
               PRINT '(A)','cycle2> Completed sweep through all available minima'
               NOMORE=.TRUE.
            ELSE
               GOTO 864
            ENDIF
         ENDIF
         IF (.NOT.NOMORE) THEN
            MINS(J3)=CONNDOING
            IF (MLP3T.OR.MLPB3T) THEN
               READ(UMIN,REC=CONNDOING) (SPOINTS(L2),L2=1,NATOMS)
               PRINT '(4(A,I6))','cycle2> Perturbation ',CONNATTEMPTS,' for minimum ',CONNDOING,' with ', &
  &                            MINCONN(CONNDOING),' connections'
               DO L2=1,NATOMS
                  IF (.NOT.FROZEN(L2)) THEN
                     SPOINTS(L2)=SPOINTS(L2)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                  ENDIF
               ENDDO
            ELSE
               READ(UMIN,REC=CONNDOING) (SPOINTS(L2),L2=1,3*NATOMS)
               PRINT '(4(A,I6))','cycle2> Perturbation ',CONNATTEMPTS,' for minimum ',CONNDOING,' with ', &
  &                            MINCONN(CONNDOING),' connections'
               DO L2=1,NATOMS
                  IF (.NOT.FROZEN(L2)) THEN
                     SPOINTS(3*(L2-1)+1)=SPOINTS(3*(L2-1)+1)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                     SPOINTS(3*(L2-1)+2)=SPOINTS(3*(L2-1)+2)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                     SPOINTS(3*(L2-1)+3)=SPOINTS(3*(L2-1)+3)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                  ENDIF
               ENDDO
            ENDIF
            CONNATTEMPTS=CONNATTEMPTS+1
            IF (CONNATTEMPTS.EQ.MAXTSATTEMPTS) THEN
               CONNDOING=CONNDOING+1
               CONNATTEMPTS=0
               IF (CONNDOING.GT.NMIN) THEN
                  PRINT '(A)','cycle2> Completed sweep through all available minima'
                  NOMORE=.TRUE.
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   ELSEIF (CHECKSPT) THEN
      NSPCHECK=NSPCHECK+1
      IF (NSPCHECK.LE.CHECKSPF) THEN
         IF (CHECKMINT) THEN
            READ(UMIN,REC=NSPCHECK) SPOINTS(1:NOPT)
         ELSE
            READ(UTS,REC=NSPCHECK) SPOINTS(1:NOPT)
         ENDIF
      ENDIF
   ELSE
      CALL GETPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ENDIF
   IF (ADDMINXYZT) THEN 
      IF (ADDMINXYZCHECK.LE.NMINADDXYZ) CALL ADDMINXYZODATA(J3,SPOINTS) 
   ELSEIF (CHECKSPT) THEN 
      IF (NSPCHECK.LE.CHECKSPF) CALL CHECKSPODATA(J3,SPOINTS) 
   ELSEIF (NEWCONNECTIONST) THEN 
      IF (.NOT.NOMORE) CALL NEWCONNODATA(J3,SPOINTS) 
   ELSEIF (CONNECTPAIRST) THEN 
      IF (CONNECTPAIRDONE.LE.NCONNECTPAIRS) CALL CONNECTODATA(J3,SPOINTS,FPOINTS)
   ELSE
      CALL CONNECTODATA(J3,SPOINTS,FPOINTS)
   ENDIF
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   IF (.NOT.(NOMORE.OR.(CHECKSPT.AND.(NSPCHECK.GT.CHECKSPF)).OR.(ADDMINXYZT.AND.(ADDMINXYZCHECK.GT.NMINADDXYZ)) &
 &                 .OR.(CONNECTPAIRST.AND.(CONNECTPAIRDONE.GT.NCONNECTPAIRS)))) THEN
      IF (CHECKSPT) MINID(J3)=NSPCHECK
      IF (ADDMINXYZT) MINID(J3)=ADDMINXYZCHECK
      IF (CONNECTPAIRST) MINID(J3)=CONNECTPAIRDONE
      CALL FORK_SUBR(PID(J3)) ! PID is zero in the child, non-zero in the parent
      IF (DEBUG.AND.(PID(J3).NE.0)) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(J3)
      CALL FLUSH(6) 
!     IF (PID(J3).NE.0) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(J3)
!     IF (PID(J3).EQ.0) THEN
!        WRITE(*,'(A,I8)') 'cycle2> I am the child! PID=',PID(J3)
!        CALL GETPID_SUBR(CHILDPID)
!        PRINT *,'CHILDPID=',CHILDPID
!     ENDIF
!     PRINT '(A,2I8)','cycle2> J3,PID=',J3,PID(J3)
      IF (PID(J3).EQ.0) THEN
         IF (ADDMINXYZT) THEN
            CALL SUBMITOPTIMJOB(J3,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.addminxyz.')
         ELSEIF (NEWCONNECTIONST) THEN
            CALL SUBMITOPTIMJOB(J3,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.tssearch.')
         ELSEIF (CHECKSPT) THEN
            CALL SUBMITOPTIMJOB(J3,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.checksp.')
         ELSE
            CALL SUBMITOPTIMJOB(J3,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.connect.')
         ENDIF
      ENDIF
   ELSE
      MINID(J3)=-1
      PID(J3)=-1
   ENDIF
ENDDO

cycles: DO NCYCLES=1,NATTEMPT*NCPU ! the last NCPU steps do not submit new jobs

   IF (NAVAIL.EQ.0) THEN
      IF (SLEEPTIME1.GT.0.0D0) CALL SYSTEM(TRIM(ADJUSTL(SLEEPSTRING1))) ! to allow jobs to catch up for demos
   ENDIF
   KILLED=.FALSE.
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   CALL FLUSH(6)
   IF (CONNECTPAIRST.OR.CHECKSPT.OR.ADDMINXYZT) THEN
      DO J2=1,NCPU
         IF (MINID(J2).GT.0) GOTO 113
         PRINT '(A,I6,I10)','J2,MINID=',J2,MINID(J2)
      ENDDO

      PRINT '(A)','cycle2> No more candidates to run'
      STOP
   ENDIF
113 CONTINUE
   CALL WAIT_SUBR(PIDDONE,STATUS)
11 CONTINUE
!  PRINT '(A,5I8)','cycle2> PIDDONE,STATUS,PID=',PIDDONE,STATUS,PID(1:NCPU)
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   IF (PIDDONE.GT.0) THEN
      IF (DEBUG) PRINT '(A,I8,A,I6)','cycle2> PID ',PIDDONE,' has finished with exit status ',STATUS
      PRINT '(A,I8,A,I6)','cycle2> PID ',PIDDONE,' has finished with exit status ',STATUS
      DO J2=1,NCPU
         IF (PIDDONE.EQ.PID(J2)) THEN
            IF (STATUS.NE.0) KILLED=.TRUE. ! INCOMPLETE OPTIM JOBS WOULD IDEALLY RETURN A NON-ZERO EXIT CODE 
            NEWJOB=J2
            IF (DEBUG) PRINT '(2(A,I8))','cycle2> PID ',PIDDONE,' has finished on cpu ',J2
            PRINT '(2(A,I8))','cycle2> PID ',PIDDONE,' has finished on cpu ',J2
            GOTO 10
         ENDIF
      ENDDO
      PRINT*,'ERROR - PID of completed child process not recognised: ',PIDDONE
      CALL MYSYSTEM(STATUS,DEBUG,' ls -lrt path.info*')
      GOTO 10   
!     STOP
   ELSE
112   CALL FLUSH(6) ! the child process may duplicate output without this line!
      PRINT '(A,I20)','cycle2> WARNING - WAIT returned system error code ',-PIDDONE
!
! Try calling wait again to see if this fixes things. 
! For very short OPTIM jobs WAIT may have trouble keeping up!
!
      CALL MYSYSTEM(STATUS,DEBUG,' sleep 1')
      IF (CHECKSPT.OR.ADDMINXYZT.OR.CONNECTPAIRST) THEN
         DO J2=1,NCPU
            IF (MINID(J2).GT.0) GOTO 114
         ENDDO
         PRINT '(A)','cycle2> No more candidates to check'
         CLOSE(513)
         STOP
      ENDIF
114   CONTINUE
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I8))','cycle2> on calling wait again pid=',PIDDONE,' status=',STATUS
      IF (PIDDONE.LE.0) PRINT '(2(A,I8))','cycle2> WARNING *** continuing for non-positive process id'
      IF (STATUS.EQ.0) GOTO 10
!     IF (PIDDONE.GT.0) GOTO 11
!     CALL MYSYSTEM(STATUS,DEBUG,' sleep 1')
!     GOTO 112
      STOP
   ENDIF
10 CONTINUE
!  WRITE(*,'(3(A,I8))') 'cycle2> forked connect run ',NCYCLES,' on CPU ',NEWJOB,' process id ',PID(NEWJOB)
!
!  It is important to identify OPTIM jobs that did not terminate with exit code 0 
!  In such cases KILLED should be .TRUE.
!
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   WRITE(*,'(3(A,I8))') 'cycle2> analysing result of search ',NCYCLES,' on CPU ',NEWJOB,' for process id ',PID(NEWJOB)
   WRITE(CONNSTR,'(I10)') PID(NEWJOB)
   IF (KILLED) WRITE(*,'(3(A,I8))') 'cycle2> connection ',NCYCLES,' on CPU ',NEWJOB,' was unsuccessful'
!
!  If KILLED is .TRUE. there could be a viable path.info file if DUMPALLPATHS is set in
!  OPTIM. 
!  It would be nice if we had the OPTIM exit status when running in a distributed
!  environment - but we don;t! Instead we have the exit status of the attempt to copy
!  back the path.info file for distributed environments, and the OPTIM exit status
!  for SMP. We can go ahead and try to analyse a path.info file so long as it exists!
!
   IF (BHINTERPT.OR.BISECTT) THEN
      CALL MYSYSTEM(STATUS,DEBUG,'cp min.data.info.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info')
   ELSEIF (CHECKSPT) THEN
      INQUIRE(FILE='frqs.min',EXIST=MINFRQSDUMP)
      IF ((GETMINFRQST).AND.(.NOT.(MINFRQSDUMP))) THEN
         OPEN(513,FILE='frqs.min',POSITION='APPEND',ACTION='WRITE',STATUS='UNKNOWN')
      ENDIF
      INQUIRE(FILE='frqs.ts',EXIST=TSFRQSDUMP)
      IF ((GETTSFRQST).AND.(.NOT.(TSFRQSDUMP))) THEN
         OPEN(513,FILE='frqs.ts',POSITION='APPEND',ACTION='WRITE',STATUS='UNKNOWN')
      ENDIF
      INQUIRE(FILE='OPTIM.checksp.'//TRIM(ADJUSTL(CONNSTR)),EXIST=OPTEST)
      IF (OPTEST) THEN
         CALL MYSYSTEM(STATUS,DEBUG,'grep -c CONVERGED OPTIM.checksp.'//TRIM(ADJUSTL(CONNSTR)) &
  &                    // ' > convcheck.'//TRIM(ADJUSTL(CONNSTR)))
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE='convcheck.'//TRIM(ADJUSTL(CONNSTR)))
         READ(LUNIT,*) NDUMMY 
         CLOSE(LUNIT)
         IF (NDUMMY.EQ.1) THEN
            PRINT '(3(A,I8))','cycle2> stationary point ',MINID(NEWJOB), &
  &                           ' process id ',PID(NEWJOB),' converged on cpu ',NEWJOB
            !kr366> get frequency from frqs.dump
            IF (GETMINFRQST) THEN
               OPEN(514,FILE='frqs.dump.'//TRIM(ADJUSTL(CONNSTR)),ACTION='READ',STATUS='OLD')
               READ(514,*) MINFRQ 
               WRITE(513,'(I8,F20.10)') MINID(NEWJOB) , MINFRQ
               CLOSE(514)
               PRINT '(A,I8)','cycle2> frequency calculated for minimum: ',MINID(NEWJOB)
               IF(.NOT. DEBUG) THEN
                   CALL MYSYSTEM(STATUS,DEBUG,'rm frqs.dump.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info.' &
  &                    //TRIM(ADJUSTL(CONNSTR))//' OPTIM.checksp.'//TRIM(ADJUSTL(CONNSTR))//' convcheck.' &
  &                    //TRIM(ADJUSTL(CONNSTR))//' start.'//TRIM(ADJUSTL(CONNSTR))//' odata.' &
  &                    //TRIM(ADJUSTL(CONNSTR)))
               ENDIF
            ENDIF
            IF (GETTSFRQST) THEN
               OPEN(514,FILE='frqs.dump.'//TRIM(ADJUSTL(CONNSTR)),ACTION='READ',STATUS='OLD')
               READ(514,*) MINFRQ
               WRITE(513,'(I8,F20.10)') MINID(NEWJOB) , MINFRQ
               CLOSE(514)
               PRINT '(A,I8)','cycle2> frequency calculated for minimum: ',MINID(NEWJOB)
               IF(.NOT. DEBUG) THEN
                   CALL MYSYSTEM(STATUS,DEBUG,'rm frqs.dump.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info.' &
  &                    //TRIM(ADJUSTL(CONNSTR))//' OPTIM.checksp.'//TRIM(ADJUSTL(CONNSTR))//' convcheck.' &
  &                    //TRIM(ADJUSTL(CONNSTR))//' start.'//TRIM(ADJUSTL(CONNSTR))//' odata.' &
  &                    //TRIM(ADJUSTL(CONNSTR)))
               ENDIF
            ENDIF
         ELSE
            PRINT '(3(A,I8))','cycle2> WARNING *** stationary point ',MINID(NEWJOB), &
  &                           ' process id ',PID(NEWJOB),' no CONVERGED message ',NEWJOB
         ENDIF
      ELSE
         PRINT '(3(A,I8))','cycle2> ERROR *** no OPTIM output file for stationary point ',MINID(NEWJOB), &
  &                        ' process id ',PID(NEWJOB),' cpu ',NEWJOB
      ENDIF
      STATUS=-1
   ELSEIF (ADDMINXYZT) THEN
      CALL MYSYSTEM(STATUS,DEBUG,'cp min.data.info.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info')
      IF (STATUS.EQ.0) THEN ! the file exists, so try to analyse it!
         CALL GETALLMIN(MINS(NEWJOB),MINF(NEWJOB))
      ELSE
         PRINT '(2(A,I6))','cycle2> WARNING *** no min.data.info file found for xyz run ',MINID(NEWJOB),' process id ',PID(NEWJOB)
      ENDIF
      STATUS=-1
   ELSE
      CALL MYSYSTEM(STATUS,DEBUG,'cp path.info.'//TRIM(ADJUSTL(CONNSTR))//' path.info')
   ENDIF
   IF (STATUS.EQ.0) THEN ! the file exists, so try to analyse it!
      IF (BHINTERPT.OR.BISECTT) THEN
         CALL GETALLMIN(MINS(NEWJOB),MINF(NEWJOB))
      ELSE ! IF (TRIPLES) THEN
         CALL GETALLPATHS
!     ELSE
!        CALL GETNEWPATH(MINS(NEWJOB),MINF(NEWJOB))
      ENDIF
   ELSE IF (DEBUG) THEN
      WRITE(*,*) "cycle2> warning: path.info file does not exist. Can't analyse results."
   ENDIF

   IF (NCYCLES.LE.NATTEMPT*NCPU-NCPU) THEN ! submit replacement job
      IF ((.NOT.DEBUG).AND.(.NOT.COPYOPTIMT)) &
   &  CALL MYSYSTEM(STATUS,DEBUG,'rm -f *' // TRIM(ADJUSTL(CONNSTR)) // ' > /dev/null' ) ! remove old output
      IF (SLEEPTIME2.GT.0.0D0) CALL SYSTEM(TRIM(ADJUSTL(SLEEPSTRING2))) ! to prevent us running out of source ports. Needed for M13.
      IF (DIJPAIRT) THEN
         CALL GETDPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (CONNECTREGIONT) THEN
         CALL GETRPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (CONUNCONT) THEN
         CALL GETCUPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
!     ELSEIF (PTAUT) THEN
!        CALL GETSPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (SHORTCUTT) THEN
         CALL GETSPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (UNTRAPT) THEN
         CALL GETUPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (NEWCONNECTIONST) THEN
862      IF (CONNDOING.LE.NMIN) THEN
            IF ((MINCONN(CONNDOING).GE.CONNECTIONS)) THEN
               CONNDOING=CONNDOING+1
               CONNATTEMPTS=0
               IF (CONNDOING.GT.NMIN) THEN
                  PRINT '(A)','cycle2> Completed sweep through all available minima'
                  NOMORE=.TRUE.
               ELSE
                  GOTO 862
               ENDIF
            ENDIF
            PRINT '(4(A,I6))','cycle2> Perturbation ',CONNATTEMPTS,' for minimum ', &
                CONNDOING,' with ',MINCONN(CONNDOING),' connections'
            MINS(NEWJOB)=CONNDOING
            IF (MLP3T.OR.MLPB3T) THEN
               READ(UMIN,REC=CONNDOING) (SPOINTS(L2),L2=1,NATOMS)
               PRINT '(4(A,I6))','cycle2> Perturbation ',CONNATTEMPTS,' for minimum ',CONNDOING,' with ', &
  &                            MINCONN(CONNDOING),' connections'
               DO L2=1,NATOMS
                  IF (.NOT.FROZEN(L2)) THEN
                     SPOINTS(L2)=SPOINTS(L2)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                  ENDIF
               ENDDO
            ELSE
               READ(UMIN,REC=CONNDOING) (SPOINTS(L2),L2=1,3*NATOMS)
               DO L2=1,NATOMS
                  IF (.NOT.FROZEN(L2)) THEN
                     SPOINTS(3*(L2-1)+1)=SPOINTS(3*(L2-1)+1)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                     SPOINTS(3*(L2-1)+2)=SPOINTS(3*(L2-1)+2)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                     SPOINTS(3*(L2-1)+3)=SPOINTS(3*(L2-1)+3)+(DPRAND()-0.5D0)*2.0D0*PERTVALUE
                  ENDIF
               ENDDO
            ENDIF
            CONNATTEMPTS=CONNATTEMPTS+1
            IF (CONNATTEMPTS.EQ.MAXTSATTEMPTS) THEN
               CONNDOING=CONNDOING+1
               CONNATTEMPTS=0
               IF (CONNDOING.GT.NMIN) THEN
                  PRINT '(A)','cycle2> Completed sweep through all available minima'
                  NOMORE=.TRUE.
               ENDIF
            ENDIF
         ENDIF
      ELSEIF (FREEPAIRT) THEN
         CALL GETFREEPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (CONNECTPAIRST) THEN
         CONNECTPAIRDONE=CONNECTPAIRDONE+1
         IF (CONNECTPAIRDONE.LE.NCONNECTPAIRS) THEN
            IF (CONNECTPAIRDONE.GE.NCONNECTPAIRS) THEN
               NOMORE=.TRUE.
               PRINT '(A)','cycle2> Completed sweep through all available pairs of minima'
            ELSE
               CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in GETCONNECTPAIR
               CALL GETCONNECTPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
            ENDIF
         ENDIF
      ELSEIF (USEPAIRST) THEN
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in GETUSEPAIR
         CALL GETUSEPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (ADDMINXYZT) THEN
         ADDMINXYZCHECK=ADDMINXYZCHECK+1
         IF (ADDMINXYZCHECK.LE.NMINADDXYZ) THEN
            READ(XYZUNIT,'(I6)') NDUMMY
            READ(XYZUNIT,'(A)') SDUMMY
            READ(XYZUNIT,*) (SDUMMY,SPOINTS(3*(J1-1)+1:3*(J1-1)+3),J1=1,NATOMS)
         ENDIF
      ELSEIF (CHECKSPT) THEN
         NSPCHECK=NSPCHECK+1
         IF (NSPCHECK.LE.CHECKSPF) THEN
            IF (CHECKMINT) THEN
               READ(UMIN,REC=NSPCHECK) SPOINTS(1:NOPT)
            ELSE
               READ(UTS,REC=NSPCHECK) SPOINTS(1:NOPT)
            ENDIF
         ENDIF
      ELSE
         CALL GETPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ENDIF

!     WRITE(*,'(2(A,I6),A,F12.1,A,F12.3,A,I8)') 'cycle2> connecting minima ',MINS(NEWJOB),' and ', &
! &      MINF(NEWJOB), ' distance=',DISTANCE,' |Pfold diff|=',ABS(GPFOLD(MINS(NEWJOB))-GPFOLD(MINF(NEWJOB))),' rejects=',NREJ
      IF (CHECKSPT) THEN 
         IF (NSPCHECK.LE.CHECKSPF) CALL CHECKSPODATA(NEWJOB,SPOINTS) 
      ELSEIF (NEWCONNECTIONST) THEN 
         IF (.NOT.NOMORE) CALL NEWCONNODATA(NEWJOB,SPOINTS) 
      ELSEIF (CONNECTPAIRST) THEN 
         IF (CONNECTPAIRDONE.LE.NCONNECTPAIRS) CALL CONNECTODATA(NEWJOB,SPOINTS,FPOINTS)
      ELSEIF (ADDMINXYZT) THEN 
         IF (ADDMINXYZCHECK.LE.NMINADDXYZ) CALL ADDMINXYZODATA(NEWJOB,SPOINTS) 
      ELSE
         CALL CONNECTODATA(NEWJOB,SPOINTS,FPOINTS)
      ENDIF
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      IF (.NOT.(NOMORE.OR.(CHECKSPT.AND.(NSPCHECK.GT.CHECKSPF)).OR.(ADDMINXYZT.AND.(ADDMINXYZCHECK.GT.NMINADDXYZ)))) THEN
         IF (CHECKSPT) MINID(NEWJOB)=NSPCHECK
         IF (ADDMINXYZT) MINID(NEWJOB)=ADDMINXYZCHECK
         CALL FORK_SUBR(PID(NEWJOB))
         IF (DEBUG.AND.(PID(NEWJOB).NE.0)) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(NEWJOB)
         CALL FLUSH(6) ! the child process may duplicate output without this line!
!        IF (DEBUG.AND.(PID(NEWJOB).NE.0)) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(NEWJOB)
!        IF (PID(NEWJOB).NE.0) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(NEWJOB)
!        IF (PID(NEWJOB).EQ.0) WRITE(*,'(A,I8)') 'cycle2> I am the child! PID=',PID(NEWJOB)
         IF (PID(NEWJOB).EQ.0) THEN
            IF (ADDMINXYZT) THEN
               CALL SUBMITOPTIMJOB(NEWJOB,CHARMMT,UNRST,NEWJOB,EXEC,DEBUG,'OPTIM.addminxyz.')
            ELSEIF (CHECKSPT) THEN
               CALL SUBMITOPTIMJOB(NEWJOB,CHARMMT,UNRST,NEWJOB,EXEC,DEBUG,'OPTIM.checksp.')
            ELSEIF (NEWCONNECTIONST) THEN
               CALL SUBMITOPTIMJOB(NEWJOB,CHARMMT,UNRST,NEWJOB,EXEC,DEBUG,'OPTIM.tssearch.')
            ELSE
               CALL SUBMITOPTIMJOB(NEWJOB,CHARMMT,UNRST,NEWJOB,EXEC,DEBUG,'OPTIM.connect.')
            ENDIF
         ENDIF
      ELSE
         MINID(NEWJOB)=-1
         PID(NEWJOB)=-1
      ENDIF
   ENDIF

   IF (MOD(NCYCLES,NCPU).EQ.0) THEN 
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
      WRITE(*,'(5(A,I8))') 'cycle2> end of cycle ',NCYCLES/NCPU,' new min=',NMIN-NMINOLD,' new ts=',NTS-NTSOLD, &
  &                          ' total min=',NMIN,' total ts=',NTS
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
      NMINOLD=NMIN; NTSOLD=NTS
      IF ((NPAIRDONE.GT.0).AND.(.NOT.DUMMYRUNT)) THEN
         OPEN(UNIT=1,FILE='pairs.data',STATUS='UNKNOWN')
         WRITE(1,'(2I8)') (PAIR1(J1),PAIR2(J1),J1=1,NPAIRDONE)
         CLOSE(1)
      ENDIF
      IF (NMINDONE.GT.0) THEN
         OPEN(UNIT=1,FILE='min.done',STATUS='UNKNOWN')
         WRITE(1,'(I8)') (MINDONE(J1),J1=1,NMINDONE)
         CLOSE(1)
      ENDIF
      IF (PFOLDINT.NE.0) THEN
         IF (MOD(NCYCLES,PFOLDINT*NCPU).EQ.0) THEN
            IF ((NPFOLD.GT.0).AND.(.NOT.DIJPAIRT)) CALL PFOLD
            OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
            WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
            CLOSE(1)
         ENDIF
      ENDIF
      IF (RATESCYCLET) CALL GETRATESCYCLE
      IF (TARGETHIT) EXIT
   ENDIF

   INQUIRE(FILE='STOP',EXIST=STOPCYCLE)
   IF (STOPCYCLE) THEN
      PRINT '(A)','cycle2> File STOP detected - exit'
      EXIT
   ENDIF
   INQUIRE(FILE='pathdata.change',EXIST=CHANGEPD)
   IF (CHANGEPD) THEN
      PRINT '(A)','cycle2> rereading parameter file'
      CALL MYSYSTEM(STATUS,DEBUG,'mv pathdata pathdata.orig')
      CALL MYSYSTEM(STATUS,DEBUG,'mv pathdata.change pathdata')
      CALL KEYWORDS
   ENDIF

ENDDO CYCLES

RETURN

END SUBROUTINE CYCLE2

SUBROUTINE GETRATESCYCLE
USE SAVESTATE
USE COMMONS
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE
DOUBLE PRECISION FRICTIONFAC, TEMPSAVE, PFNORM1, PFNORM2
INTEGER J1, J2, STATUS
LOGICAL LDEBUG

TEMPSAVE=TEMPERATURE
RATESUNIT=GETUNIT()
OPEN(RATESUNIT,FILE='NGT.rates',STATUS='UNKNOWN')
!
!  Save state.
!
IF (ALLOCATED(MINGROUP)) DEALLOCATE(MINGROUP)
IF (.NOT.SHANNONT) THEN
ALLOCATE(EMINSAVE(NMIN),PFMINSAVE(NMIN),ETSSAVE(NTS),KPLUSSAVE(NTS),KMINUSSAVE(NTS),TOPPOINTERSAVE(NMIN), &
  &         PLUSSAVE(NTS),MINUSSAVE(NTS),POINTERMSAVE(NTS),POINTERPSAVE(NTS),MINGROUP(NMIN), &
  &         LOCATIONASAVE(NMINA),LOCATIONBSAVE(NMINB))
ALLOCATE(FVIBMINSAVE(NMIN),HORDERMINSAVE(NMIN),IXMINSAVE(NMIN),IYMINSAVE(NMIN),IZMINSAVE(NMIN), &
         NCONNSAVE(NMIN),GPFOLDSAVE(NMIN))
ENDIF
NMINASAVE=NMINA; NMINBSAVE=NMINB; NMINSAVE=NMIN; NTSSAVE=NTS; LOCATIONASAVE(1:NMINA)=LOCATIONA(1:NMINA)
LOCATIONBSAVE(1:NMINB)=LOCATIONB(1:NMINB); EMINSAVE(1:NMIN)=EMIN(1:NMIN); PFMINSAVE(1:NMIN)=PFMIN(1:NMIN)
ETSSAVE(1:NTS)=ETS(1:NTS); KPLUSSAVE(1:NTS)=KPLUS(1:NTS); KMINUSSAVE(1:NTS)=KMINUS(1:NTS)
TOPPOINTERSAVE(1:NMIN)=TOPPOINTER(1:NMIN); PLUSSAVE(1:NTS)=PLUS(1:NTS); MINUSSAVE(1:NTS)=MINUS(1:NTS)
POINTERMSAVE(1:NTS)=POINTERM(1:NTS); POINTERPSAVE(1:NTS)=POINTERP(1:NTS)
PFMEANSAVE=PFMEAN; PFTOTALASAVE=PFTOTALA; PFTOTALBSAVE=PFTOTALB
FVIBMINSAVE(1:NMIN)=FVIBMIN(1:NMIN); HORDERMINSAVE(1:NMIN)=HORDERMIN(1:NMIN)
IXMINSAVE(1:NMIN)=IXMIN(1:NMIN); IYMINSAVE(1:NMIN)=IYMIN(1:NMIN); IZMINSAVE(1:NMIN)=IZMIN(1:NMIN)
GPFOLDSAVE(1:NMIN)=GPFOLD(1:NMIN)

DO J1=1,NRATESCYCLETEMPS
   TEMPERATURE=RATESCYCLETEMPS(J1)
!
!  Calculate partition functions for minima as in setup.
!
      PFMEAN=-HUGE(1.0D0)
      PFNORM1=0.0D0 ! use this to calculate ratios without the pe factor
      PFNORM2=0.0D0 ! use this to calculate ratios with the pe factor

      IF (ENSEMBLE.EQ.'T') THEN
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A,G20.10)','getratescycle> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
         DO J2 = 1,NMIN
            PFMIN(J2) = -EMIN(J2)/TEMPERATURE - FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2))
            IF (PFMIN(J2).GT.PFMEAN) PFMEAN=PFMIN(J2)
         ENDDO
         DO J2 = 1,NMIN
            PFMIN(J2) = -EMIN(J2)/TEMPERATURE - FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2)) - PFMEAN
            PFNORM1=PFNORM1+EXP(- FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2))+ FVIBMIN(1)/2.0D0 + LOG(1.0D0*HORDERMIN(1)))
            PFNORM2=PFNORM2+EXP(PFMIN(J2)-PFMIN(1))
         ENDDO
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         DO J2 = 1,NMIN
            IF (TOTALE.GT.EMIN(J2)) THEN
               PFMIN(J2) = (KAPPA-1)*LOG(TOTALE-EMIN(J2)) - FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2))
               PFNORM1=PFNORM1+EXP(- FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2)) + FVIBMIN(1)/2.0D0 + LOG(1.0D0*HORDERMIN(1)))
               PFNORM2=PFNORM2+EXP(PFMIN(J2)-PFMIN(1))
               IF (PFMIN(J2).GT.PFMEAN) PFMEAN=PFMIN(J2)
            ELSE
               PFMIN(J2) = -1.0D250
            ENDIF
         ENDDO
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF
      IF (DEBUG) THEN
         WRITE(*,'(A,G20.10)') 'getratescycle> mean ln Z=',PFMEAN
      ENDIF
!     DO J2=1,NMIN
!        PFMIN(J2) = PFMIN(J2) - PFMEAN
!     ENDDO

      PFTOTALB=0.0D0
      DO J2=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J2))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0.0D0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J2=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J2))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0.0D0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))
!
!  Calculate rate constants for this temperature. The original values
!  have been saved in PLUSSAVE and MINUSSAVE above.
!
   IF (ENSEMBLE.EQ.'T') THEN
      DO J2=1,NTS
         KPLUS(J2)  = LOG(1.0D0 * HORDERMIN(PLUS(J2))  / (2.0D0 * PI*HORDERTS(J2))) + &
  &             (FVIBMIN(PLUS(J2))  - FVIBTS(J2)) / 2.0D0 - (ETS(J2) - EMIN(PLUS(J2)) )/TEMPERATURE
         IF (FRICTIONT) KPLUS(J2)=KPLUS(J2)+LOG(FRICTIONFAC(NEGEIG(J2)))
         KMINUS(J2) = LOG(1.0D0 * HORDERMIN(MINUS(J2)) / (2.0D0 * PI*HORDERTS(J2))) + &
  &             (FVIBMIN(MINUS(J2)) - FVIBTS(J2)) / 2.0D0 - (ETS(J2) - EMIN(MINUS(J2)))/TEMPERATURE
         IF (FRICTIONT) KMINUS(J2)=KMINUS(J2)+LOG(FRICTIONFAC(NEGEIG(J2)))
         IF (ZSYM(1:2).EQ.'CA') KPLUS(J2)=KPLUS(J2)+30.66356D0
         IF (ZSYM(1:2).EQ.'CA') KMINUS(J2)=KMINUS(J2)+30.66356D0
         IF (PLUS(J2).EQ.MINUS(J2)) KPLUS(J2)=KPLUS(J2)+LOG(2.0D0)
            IF (PLUS(J2).EQ.MINUS(J2)) KMINUS(J2)=KMINUS(J2)+LOG(2.0D0)
         PRINT*,'J2,KPLUS,KMINUS=',J2,KPLUS(J2),KMINUS(J2)
      ENDDO
   ELSE
      DO J2=1,NTS
         IF (TOTALE.GT.ETS(J2)) THEN
            KPLUS(J2)  = LOG(1.0D0 * HORDERMIN(PLUS(J2))  / (2*PI*HORDERTS(J2))) + &
  &            (FVIBMIN(PLUS(J2))  - FVIBTS(J2))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J2))/(TOTALE-EMIN(PLUS(J2))))
            KMINUS(J2) = LOG(1.0D0 * HORDERMIN(MINUS(J2)) / (2*PI*HORDERTS(J2))) + &
  &           (FVIBMIN(MINUS(J2)) - FVIBTS(J2))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J2))/(TOTALE-EMIN(MINUS(J2))))
            IF (ZSYM(1:2).EQ.'CA') KPLUS(J2)=KPLUS(J2)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(J2)=KMINUS(J2)+30.66356D0
            IF (PLUS(J2).EQ.MINUS(J2)) KPLUS(J2)=KPLUS(J2)+LOG(2.0D0)
            IF (PLUS(J2).EQ.MINUS(J2)) KMINUS(J2)=KMINUS(J2)+LOG(2.0D0)
         ELSE
            KPLUS(J2)=-1.0D250
            KMINUS(J2)=-1.0D250
         ENDIF
      ENDDO
   ENDIF
   IF (DEBUG) PRINT '(A,F15.5)','getratescycle> Calling NGT for T=',TEMPERATURE
   CALL NGT
   IF (J1.LT.NRATESCYCLETEMPS) TARGETHIT=.FALSE.
!
! Reset everything.
!
   NMINA=NMINASAVE; NMINB=NMINBSAVE; NMIN=NMINSAVE; NTS=NTSSAVE; LOCATIONA(1:NMINA)=LOCATIONASAVE(1:NMINA)
   LOCATIONB(1:NMINB)=LOCATIONBSAVE(1:NMINB); EMIN(1:NMIN)=EMINSAVE(1:NMIN); PFMIN(1:NMIN)=PFMINSAVE(1:NMIN)
   ETS(1:NTS)=ETSSAVE(1:NTS); KPLUS(1:NTS)=KPLUSSAVE(1:NTS); KMINUS(1:NTS)=KMINUSSAVE(1:NTS)
   TOPPOINTER(1:NMIN)=TOPPOINTERSAVE(1:NMIN); PLUS(1:NTS)=PLUSSAVE(1:NTS); MINUS(1:NTS)=MINUSSAVE(1:NTS)
   POINTERM(1:NTS)=POINTERMSAVE(1:NTS); POINTERP(1:NTS)=POINTERPSAVE(1:NTS)
   PFMEAN=PFMEANSAVE; PFTOTALA=PFTOTALASAVE; PFTOTALB=PFTOTALBSAVE
   FVIBMIN(1:NMIN)=FVIBMINSAVE(1:NMIN); HORDERMIN(1:NMIN)=HORDERMINSAVE(1:NMIN)
   IXMIN(1:NMIN)=IXMINSAVE(1:NMIN); IYMIN(1:NMIN)=IYMINSAVE(1:NMIN); IZMIN(1:NMIN)=IZMINSAVE(1:NMIN)
   GPFOLD(1:NMIN)=GPFOLDSAVE(1:NMIN)
ENDDO
CLOSE(RATESUNIT)
TEMPERATURE=TEMPSAVE
IF (.NOT.SHANNONT) THEN
   DEALLOCATE(EMINSAVE)
   DEALLOCATE(PFMINSAVE)
   DEALLOCATE(ETSSAVE)
   DEALLOCATE(KPLUSSAVE)
   DEALLOCATE(KMINUSSAVE)
   DEALLOCATE(TOPPOINTERSAVE)
   DEALLOCATE(PLUSSAVE)
   DEALLOCATE(MINUSSAVE)
   DEALLOCATE(POINTERMSAVE)
   DEALLOCATE(POINTERPSAVE)
   DEALLOCATE(MINGROUP)
   DEALLOCATE(LOCATIONASAVE)
   DEALLOCATE(LOCATIONBSAVE)
   DEALLOCATE(NCONNSAVE)
   DEALLOCATE(FVIBMINSAVE,HORDERMINSAVE,IXMINSAVE,IYMINSAVE,IZMINSAVE,GPFOLDSAVE)
   CALL MYSYSTEM(STATUS,LDEBUG,'mv NGT.rates NGT.rates.new') 
ENDIF
LDEBUG=DEBUG

END SUBROUTINE GETRATESCYCLE

!
! Do regroupfree2 at multiple temperatures
!
SUBROUTINE RFCYCLE
USE SAVESTATE
USE COMMONS
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE
DOUBLE PRECISION FRICTIONFAC, TEMPSAVE, LNZ, PEGROUP(NMIN), SGROUP(NMIN), DLNZ, FGROUP(NMIN), CGROUP(NMIN)
DOUBLE PRECISION TLOCAL
INTEGER J1, J2, NAVAIL, J3, FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1)
INTEGER EUNIT, SUNIT, CUNIT
INTEGER, ALLOCATABLE :: NMINT(:), MINGROUPT(:,:)
DOUBLE PRECISION, ALLOCATABLE :: FOFT(:,:), EOFT(:,:), SOFT(:,:), FOFT2(:,:), SOFT2(:,:), COFT(:,:), COFT2(:,:)
LOGICAL LDEBUG

TEMPSAVE=TEMPERATURE
RFUNIT=GETUNIT()
PRINT '(A,I12)','rfcycle> Allocating F of T array dimension ',NMIN*(RFMULTIN+1)
ALLOCATE(FOFT(NMIN,RFMULTIN+1),EOFT(NMIN,RFMULTIN+1),SOFT(NMIN,RFMULTIN+1), &
  & SOFT2(NMIN,RFMULTIN+1),NMINT(RFMULTIN+1),FOFT2(NMIN,RFMULTIN+1),COFT(NMIN,RFMULTIN+1),COFT2(NMIN,RFMULTIN+1))
ALLOCATE(MINGROUPT(NMIN,RFMULTIN+1))
FOFT(1:NMIN,1:RFMULTIN+1)=0.0D0
EOFT(1:NMIN,1:RFMULTIN+1)=0.0D0
SOFT(1:NMIN,1:RFMULTIN+1)=0.0D0
SOFT2(1:NMIN,1:RFMULTIN+1)=0.0D0
COFT(1:NMIN,1:RFMULTIN+1)=0.0D0
COFT2(1:NMIN,1:RFMULTIN+1)=0.0D0
FOFT2(1:NMIN,1:RFMULTIN+1)=0.0D0
!
!  Save state.
!
IF (ALLOCATED(MINGROUP)) DEALLOCATE(MINGROUP)
ALLOCATE(EMINSAVE(NMIN),PFMINSAVE(NMIN),ETSSAVE(NTS),KPLUSSAVE(NTS),KMINUSSAVE(NTS),TOPPOINTERSAVE(NMIN), &
  &         PLUSSAVE(NTS),MINUSSAVE(NTS),POINTERMSAVE(NTS),POINTERPSAVE(NTS),MINGROUP(NMIN), &
  &         LOCATIONASAVE(NMINA),LOCATIONBSAVE(NMINB))
ALLOCATE(FVIBMINSAVE(NMIN),HORDERMINSAVE(NMIN),IXMINSAVE(NMIN),IYMINSAVE(NMIN),IZMINSAVE(NMIN), &
         NCONNSAVE(NMIN),GPFOLDSAVE(NMIN))
NMINASAVE=NMINA; NMINBSAVE=NMINB; NMINSAVE=NMIN; NTSSAVE=NTS; LOCATIONASAVE(1:NMINA)=LOCATIONA(1:NMINA)
LOCATIONBSAVE(1:NMINB)=LOCATIONB(1:NMINB); EMINSAVE(1:NMIN)=EMIN(1:NMIN); PFMINSAVE(1:NMIN)=PFMIN(1:NMIN)
ETSSAVE(1:NTS)=ETS(1:NTS); KPLUSSAVE(1:NTS)=KPLUS(1:NTS); KMINUSSAVE(1:NTS)=KMINUS(1:NTS)
TOPPOINTERSAVE(1:NMIN)=TOPPOINTER(1:NMIN); PLUSSAVE(1:NTS)=PLUS(1:NTS); MINUSSAVE(1:NTS)=MINUS(1:NTS)
POINTERMSAVE(1:NTS)=POINTERM(1:NTS); POINTERPSAVE(1:NTS)=POINTERP(1:NTS)
PFMEANSAVE=PFMEAN; PFTOTALASAVE=PFTOTALA; PFTOTALBSAVE=PFTOTALB
FVIBMINSAVE(1:NMIN)=FVIBMIN(1:NMIN); HORDERMINSAVE(1:NMIN)=HORDERMIN(1:NMIN)
IXMINSAVE(1:NMIN)=IXMIN(1:NMIN); IYMINSAVE(1:NMIN)=IYMIN(1:NMIN); IZMINSAVE(1:NMIN)=IZMIN(1:NMIN)
GPFOLDSAVE(1:NMIN)=GPFOLD(1:NMIN)
REGROUPFREETHRESHSAVE=REGROUPFREETHRESH

DO J1=1,RFMULTIN+1
   TEMPERATURE=RFMULTITLOW+(J1-1)*RFMULTITINC
   REGROUPFREETHRESH=TEMPERATURE*LOG(TIMESCALE*TEMPERATURE/PLANCK)
   PRINT '(2(A,G15.5))','rfcycle> Calling regroupfree2 with T=',TEMPERATURE,' threshold=',REGROUPFREETHRESH
!
!  Calculate partition functions for minima as in setup.
!
      PFMEAN=-HUGE(1.0D0)

      IF (ENSEMBLE.EQ.'T') THEN
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A,G20.10)','getratescycle> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
         DO J2 = 1,NMIN
            PFMIN(J2) = -EMIN(J2)/TEMPERATURE - FVIBMIN(J2)/2.0D0 &
  &          +KAPPA*LOG(2.0D0*3.14159265358979D0)- LOG(1.0D0*HORDERMIN(J2)) + KAPPA*LOG(TEMPERATURE)
            IF (PFMIN(J2).GT.PFMEAN) PFMEAN=PFMIN(J2)
         ENDDO
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         DO J2 = 1,NMIN
            IF (TOTALE.GT.EMIN(J2)) THEN
               PFMIN(J2) = (KAPPA-1)*LOG(TOTALE-EMIN(J2)) - FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2))
               IF (PFMIN(J2).GT.PFMEAN) PFMEAN=PFMIN(J2)
            ELSE
               PFMIN(J2) = -1.0D250
            ENDIF
         ENDDO
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF
      IF (DEBUG) THEN
         WRITE(*,'(A,G20.10)') 'getratescycle> mean ln Z=',PFMEAN
      ENDIF
      PFMEAN=PFSHIFT
!
! Don't take out PFMEAN (actually the largest value, not the mean at the moment)
! We want to see how the complete free energy changes with T.
! However, a constant shift may be needed to prevent under/over flow!
!
      DO J2=1,NMIN
         PFMIN(J2)=PFMIN(J2)-PFSHIFT
      ENDDO

      PFTOTALB=0.0D0
      DO J2=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J2))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0.0D0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J2=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J2))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0.0D0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))
!
!  Calculate rate constants for this temperature. The original values
!  have been saved in PLUSSAVE and MINUSSAVE above.
!
   IF (ENSEMBLE.EQ.'T') THEN
      DO J2=1,NTS
         KPLUS(J2)  = LOG(1.0D0 * HORDERMIN(PLUS(J2))  / (2.0D0 * PI*HORDERTS(J2))) + &
  &             (FVIBMIN(PLUS(J2))  - FVIBTS(J2)) / 2.0D0 - (ETS(J2) - EMIN(PLUS(J2)) )/TEMPERATURE
         IF (FRICTIONT) KPLUS(J2)=KPLUS(J2)+LOG(FRICTIONFAC(NEGEIG(J2)))
         KMINUS(J2) = LOG(1.0D0 * HORDERMIN(MINUS(J2)) / (2.0D0 * PI*HORDERTS(J2))) + &
  &             (FVIBMIN(MINUS(J2)) - FVIBTS(J2)) / 2.0D0 - (ETS(J2) - EMIN(MINUS(J2)))/TEMPERATURE
         IF (FRICTIONT) KMINUS(J2)=KMINUS(J2)+LOG(FRICTIONFAC(NEGEIG(J2)))
         IF (ZSYM(1:2).EQ.'CA') KPLUS(J2)=KPLUS(J2)+30.66356D0
         IF (ZSYM(1:2).EQ.'CA') KMINUS(J2)=KMINUS(J2)+30.66356D0
         IF (PLUS(J2).EQ.MINUS(J2)) KPLUS(J2)=KPLUS(J2)+LOG(2.0D0)
            IF (PLUS(J2).EQ.MINUS(J2)) KMINUS(J2)=KMINUS(J2)+LOG(2.0D0)
      ENDDO
   ELSE
      DO J2=1,NTS
         IF (TOTALE.GT.ETS(J2)) THEN
            KPLUS(J2)  = LOG(1.0D0 * HORDERMIN(PLUS(J2))  / (2*PI*HORDERTS(J2))) + &
  &            (FVIBMIN(PLUS(J2))  - FVIBTS(J2))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J2))/(TOTALE-EMIN(PLUS(J2))))
            KMINUS(J2) = LOG(1.0D0 * HORDERMIN(MINUS(J2)) / (2*PI*HORDERTS(J2))) + &
  &           (FVIBMIN(MINUS(J2)) - FVIBTS(J2))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J2))/(TOTALE-EMIN(MINUS(J2))))
            IF (ZSYM(1:2).EQ.'CA') KPLUS(J2)=KPLUS(J2)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(J2)=KMINUS(J2)+30.66356D0
            IF (PLUS(J2).EQ.MINUS(J2)) KPLUS(J2)=KPLUS(J2)+LOG(2.0D0)
            IF (PLUS(J2).EQ.MINUS(J2)) KMINUS(J2)=KMINUS(J2)+LOG(2.0D0)
         ELSE
            KPLUS(J2)=-1.0D250
            KMINUS(J2)=-1.0D250
         ENDIF
      ENDDO
   ENDIF
   IF (DEBUG) PRINT '(A,F15.5)','rfcycle> Calling REGROUPFREE2 for T=',TEMPERATURE
   CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
   MINGROUPT(1:NMINSAVE,J1)=MINGROUP(1:NMINSAVE)
!  PRINT '(A,I6)','mingroupt(1:nmin,j1) for j1=',j1
!  PRINT '(10I5)',MINGROUPT(1:NMINSAVE,J1)
   FGROUP(1:NMIN)=0.0D0
   PEGROUP(1:NMIN)=0.0D0
   SGROUP(1:NMIN)=0.0D0
   CGROUP(1:NMIN)=0.0D0
   DO J2=1,NMINSAVE
      IF (ENSEMBLE.EQ.'T') THEN
         LNZ=-EMINSAVE(J2)/TEMPERATURE-FVIBMINSAVE(J2)/2.0D0 &
  &       +KAPPA*LOG(2.0D0*3.14159265358979D0)-LOG(1.0D0*HORDERMINSAVE(J2)) + KAPPA*LOG(TEMPERATURE) - PFMEAN
        DLNZ= EMINSAVE(J2)/TEMPERATURE**2 + KAPPA/TEMPERATURE
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         IF (TOTALE.GT.EMINSAVE(J2)) THEN
            LNZ=(KAPPA-1)*LOG(TOTALE-EMINSAVE(J2))-FVIBMINSAVE(J2)/2.0D0 - LOG(1.0D0*HORDERMINSAVE(J2)) - PFMEAN
         ENDIF
      ENDIF
      FGROUP(MINGROUP(J2))=FGROUP(MINGROUP(J2))+EXP(LNZ)
      PEGROUP(MINGROUP(J2))=PEGROUP(MINGROUP(J2))+EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)*DLNZ*TEMPERATURE**2
      SGROUP(MINGROUP(J2))=SGROUP(MINGROUP(J2)) + EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)* &
  &                ((LNZ+TEMPERATURE*DLNZ)-(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)) 
      CGROUP(MINGROUP(J2))=CGROUP(MINGROUP(J2)) + EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)* &
  &                        (EMINSAVE(J2) + KAPPA*TEMPERATURE)**2 
   ENDDO
   NMINT(J1)=NMIN
   DO J2=1,NMIN
      FOFT(J2,J1)=EMIN(J2)
      FOFT2(J2,J1)=-TEMPERATURE*LOG(FGROUP(J2)) ! EMIN(J2)
      EOFT(J2,J1)=PEGROUP(J2)
      SOFT(J2,J1)=SGROUP(J2)
      SOFT2(J2,J1)=PEGROUP(J2)/TEMPERATURE-EMIN(J2)/TEMPERATURE
      COFT(J2,J1)=KAPPA+CGROUP(J2)/TEMPERATURE**2
   ENDDO
!
! Now check the alternative heat capacity formula.
!
   CGROUP(1:NMIN)=0.0D0
   DO J2=1,NMINSAVE
      IF (ENSEMBLE.EQ.'T') THEN
         LNZ=-EMINSAVE(J2)/TEMPERATURE-FVIBMINSAVE(J2)/2.0D0 &
  &       +KAPPA*LOG(2.0D0*3.14159265358979D0)-LOG(1.0D0*HORDERMINSAVE(J2)) + KAPPA*LOG(TEMPERATURE) - PFMEAN
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         IF (TOTALE.GT.EMINSAVE(J2)) THEN
            LNZ=(KAPPA-1)*LOG(TOTALE-EMINSAVE(J2))-FVIBMINSAVE(J2)/2.0D0 - LOG(1.0D0*HORDERMINSAVE(J2)) - PFMEAN
         ENDIF
      ENDIF
      CGROUP(MINGROUP(J2))=CGROUP(MINGROUP(J2)) + EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)* &
  &                        (EMINSAVE(J2) + KAPPA*TEMPERATURE-PEGROUP(MINGROUP(J2)))**2
   ENDDO
   DO J2=1,NMIN
      COFT2(J2,J1)=KAPPA+CGROUP(J2)/TEMPERATURE**2
   ENDDO

!
! Reset everything.
!
   NMINA=NMINASAVE; NMINB=NMINBSAVE; NMIN=NMINSAVE; NTS=NTSSAVE; LOCATIONA(1:NMINA)=LOCATIONASAVE(1:NMINA)
   LOCATIONB(1:NMINB)=LOCATIONBSAVE(1:NMINB); EMIN(1:NMIN)=EMINSAVE(1:NMIN); PFMIN(1:NMIN)=PFMINSAVE(1:NMIN)
   ETS(1:NTS)=ETSSAVE(1:NTS); KPLUS(1:NTS)=KPLUSSAVE(1:NTS); KMINUS(1:NTS)=KMINUSSAVE(1:NTS)
   TOPPOINTER(1:NMIN)=TOPPOINTERSAVE(1:NMIN); PLUS(1:NTS)=PLUSSAVE(1:NTS); MINUS(1:NTS)=MINUSSAVE(1:NTS)
   POINTERM(1:NTS)=POINTERMSAVE(1:NTS); POINTERP(1:NTS)=POINTERPSAVE(1:NTS)
   PFMEAN=PFMEANSAVE; PFTOTALA=PFTOTALASAVE; PFTOTALB=PFTOTALBSAVE
   FVIBMIN(1:NMIN)=FVIBMINSAVE(1:NMIN); HORDERMIN(1:NMIN)=HORDERMINSAVE(1:NMIN)
   IXMIN(1:NMIN)=IXMINSAVE(1:NMIN); IYMIN(1:NMIN)=IYMINSAVE(1:NMIN); IZMIN(1:NMIN)=IZMINSAVE(1:NMIN)
   GPFOLD(1:NMIN)=GPFOLDSAVE(1:NMIN)
   REGROUPFREETHRESH=REGROUPFREETHRESHSAVE
ENDDO
! OPEN(RFUNIT,FILE='Fmulti',STATUS='UNKNOWN')
!  DO J2=1,NMIN
!     DO J1=1,RFMULTIN+1
!        IF (J2.GT.NMINT(J1)) CYCLE
!        WRITE(RFUNIT,'(2G15.4)') RFMULTITLOW+(J1-1)*RFMULTITINC, FOFT(J2,J1)
!        WRITE(*,'(2G15.4)') RFMULTITLOW+(J1-1)*RFMULTITINC, FOFT(J2,J1)
!     ENDDO
!     WRITE(RFUNIT,'(A)') ' '
!     WRITE(*,'(A)') ' '
!  ENDDO
! CLOSE(RFUNIT)
RFUNIT=GETUNIT()
OPEN(RFUNIT,FILE='Fmulti',STATUS='UNKNOWN')
EUNIT=GETUNIT()
OPEN(EUNIT,FILE='Emulti',STATUS='UNKNOWN')
SUNIT=GETUNIT()
OPEN(SUNIT,FILE='Smulti',STATUS='UNKNOWN')
CUNIT=GETUNIT()
OPEN(CUNIT,FILE='Cmulti',STATUS='UNKNOWN')
DO J2=1,NMIN
   DO J3=1,J2-1
      DO J1=1,RFMULTIN+1
         IF (ABS(FOFT(MINGROUPT(J2,J1),J1)-FOFT(MINGROUPT(J3,J1),J1)).GT.EDIFFTOL) GOTO 642
      ENDDO
      PRINT '(A,I6,A,I6)','F(T) is the same for PE minimum ',J2,' and PE minimum ',J3
      GOTO 753
642   CONTINUE
   ENDDO
   DO J1=1,RFMULTIN+1
      WRITE(RFUNIT,'(3G20.10)') RFMULTITLOW+(J1-1)*RFMULTITINC, FOFT(MINGROUPT(J2,J1),J1),FOFT2(MINGROUPT(J2,J1),J1)
   ENDDO
   WRITE(RFUNIT,'(A)') ' '
   WRITE(*,'(A)') ' '
   DO J1=1,RFMULTIN+1
      WRITE(EUNIT,'(2G20.10)') RFMULTITLOW+(J1-1)*RFMULTITINC, EOFT(MINGROUPT(J2,J1),J1)
   ENDDO
   WRITE(EUNIT,'(A)') ' '
   WRITE(*,'(A)') ' '
   DO J1=1,RFMULTIN+1
      TLOCAL=RFMULTITLOW+(J1-1)*RFMULTITINC
      WRITE(SUNIT,'(5G20.10)') TLOCAL, &
  &                           (EOFT(MINGROUPT(J2,J1),J1)-FOFT(MINGROUPT(J2,J1),J1))/TLOCAL, &
  &                            EOFT(MINGROUPT(J2,J1),J1)-FOFT(MINGROUPT(J2,J1),J1),&
  &                            SOFT(MINGROUPT(J2,J1),J1),SOFT2(MINGROUPT(J2,J1),J1)
   ENDDO
   WRITE(SUNIT,'(A)') ' '
   WRITE(*,'(A)') ' '
   DO J1=1,RFMULTIN+1
      TLOCAL=RFMULTITLOW+(J1-1)*RFMULTITINC
      WRITE(CUNIT,'(3G20.10)') TLOCAL,COFT(MINGROUPT(J2,J1),J1)-(EOFT(MINGROUPT(J2,J1),J1)/TLOCAL)**2,COFT2(MINGROUPT(J2,J1),J1)
   ENDDO
   WRITE(CUNIT,'(A)') ' '
   WRITE(*,'(A)') ' '
753 CONTINUE
ENDDO
CLOSE(RFUNIT)
CLOSE(EUNIT)
CLOSE(SUNIT)
CLOSE(CUNIT)
TEMPERATURE=TEMPSAVE
DEALLOCATE(EMINSAVE)
DEALLOCATE(PFMINSAVE)
DEALLOCATE(ETSSAVE)
DEALLOCATE(KPLUSSAVE)
DEALLOCATE(KMINUSSAVE)
DEALLOCATE(TOPPOINTERSAVE)
DEALLOCATE(PLUSSAVE)
DEALLOCATE(MINUSSAVE)
DEALLOCATE(POINTERMSAVE)
DEALLOCATE(POINTERPSAVE)
DEALLOCATE(MINGROUP)
DEALLOCATE(LOCATIONASAVE)
DEALLOCATE(LOCATIONBSAVE)
DEALLOCATE(NCONNSAVE)
DEALLOCATE(FVIBMINSAVE,HORDERMINSAVE,IXMINSAVE,IYMINSAVE,IZMINSAVE,GPFOLDSAVE,FOFT,EOFT,SOFT)
LDEBUG=DEBUG

END SUBROUTINE RFCYCLE
!
! Do regroupfree2 over a temperature range and averge thermodynamics
! over groups chosen by KMC criteria.
! Identify group from previous cycle by tracking one potential energy minimum that
! lies in it: NCURRENTMIN
!
SUBROUTINE RFKMC
USE SAVESTATE
USE COMMONS
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE
DOUBLE PRECISION FRICTIONFAC, TEMPSAVE, LNZ, PEGROUP(NMIN), DLNZ, CGROUP(NMIN), DPRAND, RANDOM, DELTATIME
DOUBLE PRECISION TLOCAL, DUMMY, ZPROB, RATE
INTEGER J1, J2, NAVAIL, FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1), J3, K1, NCURRENTMIN, NGROUPCONN, NGROUP
INTEGER NINGROUP, CGMEMBER(NMIN)
INTEGER EUNIT, SUNIT, CUNIT
INTEGER, ALLOCATABLE :: NMINT(:), MINGROUPT(:,:)
DOUBLE PRECISION, ALLOCATABLE :: FOFT(:), EOFT(:), SOFT2(:), COFT(:), COFT2(:)
DOUBLE PRECISION, ALLOCATABLE :: FOFTSUM(:), EOFTSUM(:), SOFT2SUM(:), COFTSUM(:), COFT2SUM(:)
LOGICAL LDEBUG

TEMPSAVE=TEMPERATURE
RFUNIT=GETUNIT()
PRINT '(A,I12)','rfkmc> Allocating F of T array dimension ',NMIN*(RFKMCN+1)
ALLOCATE(FOFT(RFKMCN+1),EOFT(RFKMCN+1),SOFT2(RFKMCN+1),NMINT(RFKMCN+1),COFT(RFKMCN+1),COFT2(RFKMCN+1))
ALLOCATE(FOFTSUM(RFKMCN+1),EOFTSUM(RFKMCN+1),SOFT2SUM(RFKMCN+1),COFTSUM(RFKMCN+1),COFT2SUM(RFKMCN+1))
ALLOCATE(MINGROUPT(NMIN,RFKMCN+1))
FOFTSUM(1:RFKMCN+1)=0.0D0
EOFTSUM(1:RFKMCN+1)=0.0D0
SOFT2SUM(1:RFKMCN+1)=0.0D0
COFTSUM(1:RFKMCN+1)=0.0D0
COFT2SUM(1:RFKMCN+1)=0.0D0

!
!  Save state.
!
IF (ALLOCATED(MINGROUP)) DEALLOCATE(MINGROUP)
ALLOCATE(EMINSAVE(NMIN),PFMINSAVE(NMIN),ETSSAVE(NTS),KPLUSSAVE(NTS),KMINUSSAVE(NTS),TOPPOINTERSAVE(NMIN), &
  &         PLUSSAVE(NTS),MINUSSAVE(NTS),POINTERMSAVE(NTS),POINTERPSAVE(NTS),MINGROUP(NMIN), &
  &         LOCATIONASAVE(NMINA),LOCATIONBSAVE(NMINB))
ALLOCATE(FVIBMINSAVE(NMIN),HORDERMINSAVE(NMIN),IXMINSAVE(NMIN),IYMINSAVE(NMIN),IZMINSAVE(NMIN), &
         NCONNSAVE(NMIN),GPFOLDSAVE(NMIN))
NMINASAVE=NMINA; NMINBSAVE=NMINB; NMINSAVE=NMIN; NTSSAVE=NTS; LOCATIONASAVE(1:NMINA)=LOCATIONA(1:NMINA)
LOCATIONBSAVE(1:NMINB)=LOCATIONB(1:NMINB); EMINSAVE(1:NMIN)=EMIN(1:NMIN); PFMINSAVE(1:NMIN)=PFMIN(1:NMIN)
ETSSAVE(1:NTS)=ETS(1:NTS); KPLUSSAVE(1:NTS)=KPLUS(1:NTS); KMINUSSAVE(1:NTS)=KMINUS(1:NTS)
TOPPOINTERSAVE(1:NMIN)=TOPPOINTER(1:NMIN); PLUSSAVE(1:NTS)=PLUS(1:NTS); MINUSSAVE(1:NTS)=MINUS(1:NTS)
POINTERMSAVE(1:NTS)=POINTERM(1:NTS); POINTERPSAVE(1:NTS)=POINTERP(1:NTS)
PFMEANSAVE=PFMEAN; PFTOTALASAVE=PFTOTALA; PFTOTALBSAVE=PFTOTALB
FVIBMINSAVE(1:NMIN)=FVIBMIN(1:NMIN); HORDERMINSAVE(1:NMIN)=HORDERMIN(1:NMIN)
IXMINSAVE(1:NMIN)=IXMIN(1:NMIN); IYMINSAVE(1:NMIN)=IYMIN(1:NMIN); IZMINSAVE(1:NMIN)=IZMIN(1:NMIN)
GPFOLDSAVE(1:NMIN)=GPFOLD(1:NMIN)
REGROUPFREETHRESHSAVE=REGROUPFREETHRESH

PRINT '(A,G20.10)','rfkmc> Setting observation time scale to T increment over rate of change=',RFKMCTINC/RFKMCTRATE

RFUNIT=GETUNIT()
OPEN(RFUNIT,FILE='FKMC',STATUS='UNKNOWN')
EUNIT=GETUNIT()
OPEN(EUNIT,FILE='EKMC',STATUS='UNKNOWN')
SUNIT=GETUNIT()
OPEN(SUNIT,FILE='SKMC',STATUS='UNKNOWN')
CUNIT=GETUNIT()
OPEN(CUNIT,FILE='CKMC',STATUS='UNKNOWN')

DO K1=1,RFKMCSTEPS

PRINT '(A,I6)','rfkmc> KMC step number ',K1

DO J1=1,RFKMCN+1
   TEMPERATURE=RFKMCTSTART+(J1-1)*RFKMCTINC
!
! The temperature changes by RFKMCTINC per step, which corresponds to a time
! step of RFKMCTINC/RFKMCTRATE, where RFKMCTRATE is the time rate of change.
!
   DELTATIME=ABS(RFKMCTINC/RFKMCTRATE)
   REGROUPFREETHRESH=TEMPERATURE*LOG(TIMESCALE*TEMPERATURE/PLANCK)
   IF (REGROUPKMCT) THEN
      PRINT '(2(A,G15.5))','rfcycle> Calling regroupfree2 with T=',TEMPERATURE,' for stochastic regrouping'
   ELSE
      PRINT '(2(A,G15.5))','rfcycle> Calling regroupfree2 with T=',TEMPERATURE,' threshold=',REGROUPFREETHRESH
   ENDIF
!
!  Calculate partition functions for minima as in setup.
!
      PFMEAN=-HUGE(1.0D0)

      IF (ENSEMBLE.EQ.'T') THEN
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A,G20.10)','getratescycle> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
         DO J2 = 1,NMIN
            PFMIN(J2) = -EMIN(J2)/TEMPERATURE - FVIBMIN(J2)/2.0D0 &
  &          +KAPPA*LOG(2.0D0*3.14159265358979D0)- LOG(1.0D0*HORDERMIN(J2)) + KAPPA*LOG(TEMPERATURE)
            IF (PFMIN(J2).GT.PFMEAN) PFMEAN=PFMIN(J2)
         ENDDO
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         DO J2 = 1,NMIN
            IF (TOTALE.GT.EMIN(J2)) THEN
               PFMIN(J2) = (KAPPA-1)*LOG(TOTALE-EMIN(J2)) - FVIBMIN(J2)/2.0D0 - LOG(1.0D0*HORDERMIN(J2))
               IF (PFMIN(J2).GT.PFMEAN) PFMEAN=PFMIN(J2)
            ELSE
               PFMIN(J2) = -1.0D250
            ENDIF
         ENDDO
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF
      IF (DEBUG) THEN
         WRITE(*,'(A,G20.10)') 'getratescycle> mean ln Z=',PFMEAN
      ENDIF
      PFMEAN=PFSHIFT
!
! Don't take out PFMEAN (actually the largest value, not the mean at the moment)
! We want to see how the complete free energy changes with T.
! However, a constant shift may be needed to prevent under/over flow!
!
      DO J2=1,NMIN
         PFMIN(J2)=PFMIN(J2)-PFSHIFT
      ENDDO

      PFTOTALB=0.0D0
      DO J2=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J2))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0.0D0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J2=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J2))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0.0D0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))
!
!  Calculate rate constants for this temperature. The original values
!  have been saved in PLUSSAVE and MINUSSAVE above.
!
   IF (ENSEMBLE.EQ.'T') THEN
      DO J2=1,NTS
         KPLUS(J2)  = LOG(1.0D0 * HORDERMIN(PLUS(J2))  / (2.0D0 * PI*HORDERTS(J2))) + &
  &             (FVIBMIN(PLUS(J2))  - FVIBTS(J2)) / 2.0D0 - (ETS(J2) - EMIN(PLUS(J2)) )/TEMPERATURE
         IF (FRICTIONT) KPLUS(J2)=KPLUS(J2)+LOG(FRICTIONFAC(NEGEIG(J2)))
         KMINUS(J2) = LOG(1.0D0 * HORDERMIN(MINUS(J2)) / (2.0D0 * PI*HORDERTS(J2))) + &
  &             (FVIBMIN(MINUS(J2)) - FVIBTS(J2)) / 2.0D0 - (ETS(J2) - EMIN(MINUS(J2)))/TEMPERATURE
         IF (FRICTIONT) KMINUS(J2)=KMINUS(J2)+LOG(FRICTIONFAC(NEGEIG(J2)))
         IF (ZSYM(1:2).EQ.'CA') KPLUS(J2)=KPLUS(J2)+30.66356D0
         IF (ZSYM(1:2).EQ.'CA') KMINUS(J2)=KMINUS(J2)+30.66356D0
         IF (PLUS(J2).EQ.MINUS(J2)) KPLUS(J2)=KPLUS(J2)+LOG(2.0D0)
            IF (PLUS(J2).EQ.MINUS(J2)) KMINUS(J2)=KMINUS(J2)+LOG(2.0D0)
      ENDDO
   ELSE
      DO J2=1,NTS
         IF (TOTALE.GT.ETS(J2)) THEN
            KPLUS(J2)  = LOG(1.0D0 * HORDERMIN(PLUS(J2))  / (2*PI*HORDERTS(J2))) + &
  &            (FVIBMIN(PLUS(J2))  - FVIBTS(J2))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J2))/(TOTALE-EMIN(PLUS(J2))))
            KMINUS(J2) = LOG(1.0D0 * HORDERMIN(MINUS(J2)) / (2*PI*HORDERTS(J2))) + &
  &           (FVIBMIN(MINUS(J2)) - FVIBTS(J2))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J2))/(TOTALE-EMIN(MINUS(J2))))
            IF (ZSYM(1:2).EQ.'CA') KPLUS(J2)=KPLUS(J2)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(J2)=KMINUS(J2)+30.66356D0
            IF (PLUS(J2).EQ.MINUS(J2)) KPLUS(J2)=KPLUS(J2)+LOG(2.0D0)
            IF (PLUS(J2).EQ.MINUS(J2)) KMINUS(J2)=KMINUS(J2)+LOG(2.0D0)
         ELSE
            KPLUS(J2)=-1.0D250
            KMINUS(J2)=-1.0D250
         ENDIF
      ENDDO
   ENDIF
   IF (DEBUG) PRINT '(A,F15.5)','rfkmc> Calling REGROUPFREE2 for T=',TEMPERATURE
   CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
   MINGROUPT(1:NMINSAVE,J1)=MINGROUP(1:NMINSAVE)
!
! Set current free energy group as one of the lowest if this is the first
! temperature. If we make a transition to a different group then NCURRENTMIN
! will change.
!
   DUMMY=HUGE(1.0D0)
   IF (J1.EQ.1) THEN
      DO J2=1,NMIN
         IF (EMIN(J2).LT.DUMMY) THEN
            DUMMY=EMIN(J2)
            minloop : DO J3=1,NMINSAVE
               IF (MINGROUP(J3).EQ.J2) THEN
                  NCURRENTMIN=J3
                  EXIT minloop
               ENDIF
            ENDDO minloop
         ENDIF
      ENDDO
      PRINT '(2(A,I6))','rfkmc> Assigning initial reference pe minimum ',NCURRENTMIN,' from group ',MINGROUP(NCURRENTMIN)
!
! Reset memory of pe minima in current free energy group.
!
      NINGROUP=0
      DO J2=1,NMINSAVE
         IF (MINGROUP(J2).EQ.MINGROUP(NCURRENTMIN)) THEN
            NINGROUP=NINGROUP+1
            CGMEMBER(NINGROUP)=J2
         ENDIF
      ENDDO
      PRINT '(3(A,I6))','rfkmc> reference pe minimum ',NCURRENTMIN,' is in group ',MINGROUP(NCURRENTMIN),' size=',NINGROUP
   ELSE
!
! To set the current free energy group we should choose the lowest free energy group
! containing any of the PE minima from the previous temperature. Otherwise, with stochastic
! regrouping (REGROUPKMC) the free energy can jump if we happen to choose the wrong minimum.
! Need to save the pe minima in the current free energy group after each selection of NGROUP.
!
      DUMMY=1.0D100
      DO J2=1,NINGROUP ! loop over PE minima from previous current group
!
! If free energy of new group containing pe minimum is lower than best so far, switch to that one
!
         IF (EMIN(MINGROUP(CGMEMBER(J2))).LT.DUMMY) THEN
            NCURRENTMIN=CGMEMBER(J2)
            NGROUP=MINGROUP(NCURRENTMIN)
            PRINT '(A,2I6,G20.10,I6)','rfkmc> changing current pe minimum and free energy group to ',NCURRENTMIN,NGROUP
            DUMMY=EMIN(MINGROUP(CGMEMBER(J2)))
         ENDIF
         PRINT '(A,2I6,G20.10,I6)','rfkmc> pe minimum, free energy group, free energy, current minimum: ', &
   &                CGMEMBER(J2),MINGROUP(CGMEMBER(J2)),EMIN(MINGROUP(CGMEMBER(J2))),NCURRENTMIN
      ENDDO
!
! Reset memory of pe minima in current free energy group.
!
      NINGROUP=0
      DO J2=1,NMINSAVE
         IF (MINGROUP(J2).EQ.MINGROUP(NCURRENTMIN)) THEN
            NINGROUP=NINGROUP+1
            CGMEMBER(NINGROUP)=J2
         ENDIF
      ENDDO
      PRINT '(2(A,I6))','rfkmc> reference pe minimum ',NCURRENTMIN,' is in group ',MINGROUP(NCURRENTMIN),' size=',NINGROUP
   ENDIF
   PEGROUP(1:NMIN)=0.0D0
   CGROUP(1:NMIN)=0.0D0
   DO J2=1,NMINSAVE
      IF (ENSEMBLE.EQ.'T') THEN
         LNZ=-EMINSAVE(J2)/TEMPERATURE-FVIBMINSAVE(J2)/2.0D0 &
  &       +KAPPA*LOG(2.0D0*3.14159265358979D0)-LOG(1.0D0*HORDERMINSAVE(J2)) + KAPPA*LOG(TEMPERATURE) - PFMEAN
        DLNZ= EMINSAVE(J2)/TEMPERATURE**2 + KAPPA/TEMPERATURE
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         IF (TOTALE.GT.EMINSAVE(J2)) THEN
            LNZ=(KAPPA-1)*LOG(TOTALE-EMINSAVE(J2))-FVIBMINSAVE(J2)/2.0D0 - LOG(1.0D0*HORDERMINSAVE(J2)) - PFMEAN
         ENDIF
      ENDIF
      PEGROUP(MINGROUP(J2))=PEGROUP(MINGROUP(J2))+EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)*DLNZ*TEMPERATURE**2
      CGROUP(MINGROUP(J2))=CGROUP(MINGROUP(J2)) + EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)* &
  &                        (EMINSAVE(J2) + KAPPA*TEMPERATURE)**2 
   ENDDO
   NMINT(J1)=NMIN
!
! Accumulate averages for thermodynamic quantities using the value for the group
! corresponding to the current NCURRENTMIN pe minimum.
!
   FOFT(J1)=EMIN(MINGROUP(NCURRENTMIN))
   EOFT(J1)=PEGROUP(MINGROUP(NCURRENTMIN))
   SOFT2(J1)=PEGROUP(MINGROUP(NCURRENTMIN))/TEMPERATURE-EMIN(MINGROUP(NCURRENTMIN))/TEMPERATURE
   COFT(J1)=KAPPA+CGROUP(MINGROUP(NCURRENTMIN))/TEMPERATURE**2-(PEGROUP(MINGROUP(NCURRENTMIN))/TEMPERATURE)**2
   FOFTSUM(J1)=FOFTSUM(J1)+FOFT(J1)
   EOFTSUM(J1)=EOFTSUM(J1)+EOFT(J1)
   SOFT2SUM(J1)=SOFT2SUM(J1)+SOFT2(J1)
   COFTSUM(J1)=COFTSUM(J1)+COFT(J1)
!
! Now check the alternative heat capacity formula.
!
   CGROUP(1:NMIN)=0.0D0
   DO J2=1,NMINSAVE
      IF (ENSEMBLE.EQ.'T') THEN
         LNZ=-EMINSAVE(J2)/TEMPERATURE-FVIBMINSAVE(J2)/2.0D0 &
  &       +KAPPA*LOG(2.0D0*3.14159265358979D0)-LOG(1.0D0*HORDERMINSAVE(J2)) + KAPPA*LOG(TEMPERATURE) - PFMEAN
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         IF (TOTALE.GT.EMINSAVE(J2)) THEN
            LNZ=(KAPPA-1)*LOG(TOTALE-EMINSAVE(J2))-FVIBMINSAVE(J2)/2.0D0 - LOG(1.0D0*HORDERMINSAVE(J2)) - PFMEAN
         ENDIF
      ENDIF
      CGROUP(MINGROUP(J2))=CGROUP(MINGROUP(J2)) + EXP(LNZ+EMIN(MINGROUP(J2))/TEMPERATURE)* &
  &                        (EMINSAVE(J2) + KAPPA*TEMPERATURE-PEGROUP(MINGROUP(J2)))**2
   ENDDO
   COFT2(J1)=KAPPA+CGROUP(MINGROUP(NCURRENTMIN))/TEMPERATURE**2
   COFT2SUM(J1)=COFT2SUM(J1)+COFT2(J1)
!
! Do we change groups? If we do, then we ought to use a mean waiting time.
! However, this messes up the regularity of the temperature increments, so
! for convenience try assuming the same time interval and temperature step,
! but allow the change in group.
! First calculate the probability of zero transitions for the given timescale,
! which is delta T, the temperature increment, divided by | the rate of change of T |.
!
   RATE=0.0D0
   NGROUPCONN=0
   DO J2=1,NTS
      IF (PLUS(J2).EQ.MINGROUP(NCURRENTMIN)) THEN
         IF (MINUS(J2).NE.MINGROUP(NCURRENTMIN)) THEN
            RATE=RATE+EXP(KPLUS(J2))
            NGROUPCONN=NGROUPCONN+1
            IF (DEBUG) PRINT '(A,3I6,3G20.10)','ts,+,-,k+,k-,rate after adding + =',J2,PLUS(J2),MINUS(J2),KPLUS(J2),KMINUS(J2),RATE
         ENDIF
      ENDIF
      IF (MINUS(J2).EQ.MINGROUP(NCURRENTMIN)) THEN
         IF (PLUS(J2).NE.MINGROUP(NCURRENTMIN)) THEN
            RATE=RATE+EXP(KMINUS(J2))
            NGROUPCONN=NGROUPCONN+1
            IF (DEBUG) PRINT '(A,3I6,3G20.10)','ts,+,-,k+,k-,rate after adding - =',J2,PLUS(J2),MINUS(J2),KPLUS(J2),KMINUS(J2),RATE
         ENDIF
      ENDIF
   ENDDO
   PRINT '(A,I6,A,G20.10)','rfkmc> Sum of rate constants out of current group for ',NGROUPCONN,' connections is ',RATE
   ZPROB=EXP(-RATE*DELTATIME)
   PRINT '(A,G20.10,A,G20.10)','rfkmc> Probability of zero transitions from current group in time ',DELTATIME,' is ',ZPROB
   RANDOM=DPRAND()
   IF (ZPROB.LT.RANDOM) THEN
      RANDOM=DPRAND()*RATE
      RATE=0.0D0
      choosegroup : DO J2=1,NTS
         IF (PLUS(J2).EQ.MINGROUP(NCURRENTMIN)) THEN
            IF (MINUS(J2).NE.MINGROUP(NCURRENTMIN)) THEN
               RATE=RATE+EXP(KPLUS(J2))
               PRINT '(A,3I6,2G20.10)','rfkmc> j2,k+,k-,rate,random=',J2,PLUS(J2),MINUS(J2),RATE,RANDOM
               IF (RATE.GT.RANDOM) THEN
                  NGROUP=MINUS(J2)
                  EXIT choosegroup
               ENDIF
            ENDIF
         ENDIF
         IF (MINUS(J2).EQ.MINGROUP(NCURRENTMIN)) THEN
            IF (PLUS(J2).NE.MINGROUP(NCURRENTMIN)) THEN
               RATE=RATE+EXP(KMINUS(J2))
               PRINT '(A,3I6,2G20.10)','rfkmc> j2,k+,k-,rate,random=',J2,PLUS(J2),MINUS(J2),RATE,RANDOM
               IF (RATE.GT.RANDOM) THEN
                  NGROUP=PLUS(J2)
                  EXIT choosegroup
               ENDIF
            ENDIF
         ENDIF
      ENDDO choosegroup
      PRINT '(A,I6,A,G20.10,A,I6,A,G20.10,A,F20.10)','rfkmc> Making a transition from group ',MINGROUP(NCURRENTMIN),' F=', &
   &                       EMIN(MINGROUP(NCURRENTMIN)),' to group ',NGROUP,' F=',EMIN(NGROUP),' at T=',TEMPERATURE
      minloop2 : DO J3=1,NMINSAVE
         IF (MINGROUP(J3).EQ.NGROUP) THEN
            NCURRENTMIN=J3
            EXIT minloop2
         ENDIF
      ENDDO minloop2
      PRINT '(2(A,I6))','rfkmc> Assigning new reference pe minimum ',NCURRENTMIN,' from group ',MINGROUP(NCURRENTMIN)
   ELSE
      PRINT '(A)','rfkmc> Staying in the current group'
   ENDIF
   PRINT '(A)',' '

!
! Reset everything for next temperature.
!
   NMINA=NMINASAVE; NMINB=NMINBSAVE; NMIN=NMINSAVE; NTS=NTSSAVE; LOCATIONA(1:NMINA)=LOCATIONASAVE(1:NMINA)
   LOCATIONB(1:NMINB)=LOCATIONBSAVE(1:NMINB); EMIN(1:NMIN)=EMINSAVE(1:NMIN); PFMIN(1:NMIN)=PFMINSAVE(1:NMIN)
   ETS(1:NTS)=ETSSAVE(1:NTS); KPLUS(1:NTS)=KPLUSSAVE(1:NTS); KMINUS(1:NTS)=KMINUSSAVE(1:NTS)
   TOPPOINTER(1:NMIN)=TOPPOINTERSAVE(1:NMIN); PLUS(1:NTS)=PLUSSAVE(1:NTS); MINUS(1:NTS)=MINUSSAVE(1:NTS)
   POINTERM(1:NTS)=POINTERMSAVE(1:NTS); POINTERP(1:NTS)=POINTERPSAVE(1:NTS)
   PFMEAN=PFMEANSAVE; PFTOTALA=PFTOTALASAVE; PFTOTALB=PFTOTALBSAVE
   FVIBMIN(1:NMIN)=FVIBMINSAVE(1:NMIN); HORDERMIN(1:NMIN)=HORDERMINSAVE(1:NMIN)
   IXMIN(1:NMIN)=IXMINSAVE(1:NMIN); IYMIN(1:NMIN)=IYMINSAVE(1:NMIN); IZMIN(1:NMIN)=IZMINSAVE(1:NMIN)
   GPFOLD(1:NMIN)=GPFOLDSAVE(1:NMIN)
   REGROUPFREETHRESH=REGROUPFREETHRESHSAVE
ENDDO

DO J2=1,RFKMCN+1
   TLOCAL=RFKMCTSTART+(J2-1)*RFKMCTINC
   WRITE(RFUNIT,'(2G20.10)') TLOCAL, FOFT(J2)
   WRITE(EUNIT,'(2G20.10)') TLOCAL, EOFT(J2)
   WRITE(SUNIT,'(2G20.10)') TLOCAL, SOFT2(J2)
   WRITE(CUNIT,'(3G20.10)') TLOCAL, COFT(J2), COFT2(J2)
ENDDO
WRITE(RFUNIT,'(A)') ' '
WRITE(EUNIT,'(A)') ' '
WRITE(SUNIT,'(A)') ' '
WRITE(CUNIT,'(A)') ' '

ENDDO ! this is the end of the K1 loop over KMC steps
CLOSE(RFUNIT)
CLOSE(EUNIT)
CLOSE(SUNIT)
CLOSE(CUNIT)
OPEN(RFUNIT,FILE='FKMC.final',STATUS='UNKNOWN')
OPEN(EUNIT,FILE='EKMC.final',STATUS='UNKNOWN')
OPEN(SUNIT,FILE='SKMC.final',STATUS='UNKNOWN')
OPEN(CUNIT,FILE='CKMC.final',STATUS='UNKNOWN')
DO J1=1,RFKMCN+1
   TLOCAL=RFKMCTSTART+(J1-1)*RFKMCTINC
   WRITE(RFUNIT,'(2G20.10)') TLOCAL, FOFTSUM(J1)/RFKMCSTEPS
   WRITE(EUNIT,'(2G20.10)') TLOCAL, EOFTSUM(J1)/RFKMCSTEPS
   WRITE(SUNIT,'(2G20.10)') TLOCAL, SOFT2SUM(J1)/RFKMCSTEPS
   WRITE(CUNIT,'(3G20.10)') TLOCAL, COFTSUM(J1)/RFKMCSTEPS, COFT2SUM(J1)/RFKMCSTEPS
ENDDO
CLOSE(RFUNIT)
CLOSE(EUNIT)
CLOSE(SUNIT)
CLOSE(CUNIT)

TEMPERATURE=TEMPSAVE
DEALLOCATE(EMINSAVE)
DEALLOCATE(PFMINSAVE)
DEALLOCATE(ETSSAVE)
DEALLOCATE(KPLUSSAVE)
DEALLOCATE(KMINUSSAVE)
DEALLOCATE(TOPPOINTERSAVE)
DEALLOCATE(PLUSSAVE)
DEALLOCATE(MINUSSAVE)
DEALLOCATE(POINTERMSAVE)
DEALLOCATE(POINTERPSAVE)
DEALLOCATE(MINGROUP)
DEALLOCATE(LOCATIONASAVE)
DEALLOCATE(LOCATIONBSAVE)
DEALLOCATE(NCONNSAVE)
DEALLOCATE(FVIBMINSAVE,HORDERMINSAVE,IXMINSAVE,IYMINSAVE,IZMINSAVE,GPFOLDSAVE,FOFT,EOFT)
LDEBUG=DEBUG

END SUBROUTINE RFKMC
