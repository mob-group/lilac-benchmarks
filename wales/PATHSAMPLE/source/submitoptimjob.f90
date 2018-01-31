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

SUBROUTINE SUBMITOPTIMJOB(ICPU,CHARMMT,UNRST,CONNID,EXEC,LDEBUG,JOBSTRING)
USE PORFUNCS
USE NODES, ONLY: SSHSUBMIT, SSHPARALLEL, SLURMSUBMIT
USE COMMONS, ONLY: DUMMYRUNT, AMBERT, OPEPT, AMHT, MACHINE, SLURMT
IMPLICIT NONE
INTEGER, INTENT(IN) :: ICPU
CHARACTER(LEN=10) CONNSTR1, CONNSTR2
CHARACTER(LEN=*) JOBSTRING
CHARACTER(LEN=*) EXEC
CHARACTER(LEN=256) MYJOBSTRING
INTEGER :: CHILDPID, CONNID, STATUS
LOGICAL LDEBUG, CHARMMT, UNRST

CALL GETPID_SUBR(CHILDPID)
! PRINT '(A,2I8)','in SUBMITOPTIMJOB, CHILDPID,CONNID=',CHILDPID,CONNID
WRITE(CONNSTR1,'(I10)') CHILDPID
WRITE(CONNSTR2,'(I10)') CONNID
CALL MYSYSTEM(STATUS,LDEBUG,'mv odata.' // TRIM(ADJUSTL(CONNSTR2)) // ' odata.' // TRIM(ADJUSTL(CONNSTR1)) )

IF (MACHINE) THEN
     CALL MYSYSTEM(STATUS,LDEBUG,'mv points1.inp.' // TRIM(ADJUSTL(CONNSTR2)) // ' points1.inp.' // TRIM(ADJUSTL(CONNSTR1)) )
     IF (VERIFY('connect',JOBSTRING).EQ.0) &
     & call mysystem(STATUS,LDEBUG,'mv points2.inp.' // TRIM(ADJUSTL(CONNSTR2)) // ' points2.inp.' // TRIM(ADJUSTL(CONNSTR1)) )
else ! preprocess odata.XXX to include correct name instead of input.crd
     CALL SYSTEM('sed -e "s/input.crd/input.crd.' // TRIM(ADJUSTL(CONNSTR1)) // '/"  odata.' // TRIM(ADJUSTL(CONNSTR1)) &
          &                 // ' > temp.'//TRIM(ADJUSTL(CONNSTR1)))
     CALL SYSTEM('mv temp.'//TRIM(ADJUSTL(CONNSTR1)) // ' odata.' // TRIM(ADJUSTL(CONNSTR1)) )

     IF (VERIFY('connect',JOBSTRING).EQ.0) &
      CALL MYSYSTEM(STATUS,LDEBUG,'mv finish.' // TRIM(ADJUSTL(CONNSTR2)) // ' finish.' // TRIM(ADJUSTL(CONNSTR1)) )
     IF (CHARMMT) CALL MYSYSTEM( &
       STATUS,LDEBUG,'mv input.crd.' // TRIM(ADJUSTL(CONNSTR2)) // ' input.crd.' // TRIM(ADJUSTL(CONNSTR1)) )
     IF (AMBERT.OR.OPEPT) CALL MYSYSTEM(STATUS,LDEBUG,'mv start.' // TRIM(ADJUSTL(CONNSTR2)) // ' start.' // TRIM(ADJUSTL(CONNSTR1)) )
     IF (AMHT) CALL MYSYSTEM(STATUS,LDEBUG,'mv start.' // TRIM(ADJUSTL(CONNSTR2)) // ' start.' // TRIM(ADJUSTL(CONNSTR1)) )
endif

IF ((VERIFY('path',JOBSTRING).EQ.0).AND.(VERIFY('sloppy',JOBSTRING).NE.0)) CALL MYSYSTEM( &
  STATUS,LDEBUG,'mv vector.dump.' // TRIM(ADJUSTL(CONNSTR2)) // ' vector.dump.' // TRIM(ADJUSTL(CONNSTR1)) )
IF (UNRST) CALL MYSYSTEM( &
  STATUS,LDEBUG,'mv coords.' // TRIM(ADJUSTL(CONNSTR2)) // ' coords.' // TRIM(ADJUSTL(CONNSTR1)) )

MYJOBSTRING=TRIM(ADJUSTL(EXEC))//' '//TRIM(ADJUSTL(CONNSTR1))//' > '//TRIM(ADJUSTL(JOBSTRING))//TRIM(ADJUSTL(CONNSTR1))
IF (DUMMYRUNT) MYJOBSTRING='sleep 10'
IF (LDEBUG) PRINT '(2A)','submitoptimjob> myjobstring=',TRIM(ADJUSTL(MYJOBSTRING))
IF (SSHPARALLEL) then
     CALL SSHSUBMIT(icpu,status,trim(adjustl(myjobstring)),trim(adjustl(CONNSTR1)),LDEBUG)
ELSE IF (SLURMT) THEN
     CALL SLURMSUBMIT(status,trim(adjustl(myjobstring)),trim(adjustl(CONNSTR1)),LDEBUG)
ELSE
     CALL MYSYSTEM(STATUS,LDEBUG,trim(adjustl(myjobstring)))
ENDIF
IF (STATUS.NE.0) PRINT '(A,I8)','submitoptimjob> WARNING - '//trim(adjustl(MYJOBSTRING))//' exit status=',STATUS
CALL EXIT(STATUS)
STOP

RETURN
END
