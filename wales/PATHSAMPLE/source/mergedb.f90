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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  This subroutine merges the database information in directory PATHNAME.
!
SUBROUTINE MERGEDB
USE COMMONS,ONLY : NATOMS, IYTS, IZTS, UTSDATA, UTS, HORDERMIN, TOPPOINTER, HORDERTS, PLUS, MINUS, GPFOLD, &
   &              MAXMIN, MAXTS, FVIBTS, EMIN, FVIBMIN, IXMIN, IYMIN, IZMIN, NEGEIG, PATHNAME, ETS, DEBUG, &
   &              NMIN, UNRST, CHARMMT, IDIFFTOL, EDIFFTOL, UMINDATA, UMIN, NTS, IXTS, NMINA, NMINB, &
   &              LOCATIONA, LOCATIONB, ANGLEAXIS, PERMDIST, BOXLX, BOXLY, BOXLZ, GEOMDIFFTOL, TWOD, &
   &              RIGIDBODY, BULKT, ZSYM, PERMISOMER, IMFRQT, CLOSEFILEST, AMHT, PAIRDISTMAX, DIJINITT, &
   &              DIJINITCONTT, ALLTST, ADDPT2, ADDPT3, NOPT, LPERMDIST, INITIALDIST
USE UTILS,ONLY : GETUNIT

USE PORFUNCS
IMPLICIT NONE

INTEGER J1, J2, ISTAT, NMINOLD, NTSOLD, NMINDB, NDUMMY, J3, LUNIT, J4, LPLUNIT, LPDUNIT, PLUNIT, PDUNIT
DOUBLE PRECISION LOCALPOINTS(NOPT), NEWEMIN, NEWETS
DOUBLE PRECISION NEWFVIBMIN, NEWFVIBTS, NEWPOINTSMIN(NOPT), NEWNEGEIG, LPLUS(NOPT), LMINUS(NOPT), &
  &  NEWPOINTSTS(NOPT), NEWIXMIN,  NEWIYMIN, NEWIZMIN, LPAIRDIST(PAIRDISTMAX), &
  &  NEWIXTS,  NEWIYTS, NEWIZTS, DISTANCE, DIST2, LOCALPOINTS2(NOPT), RMAT(3,3)
INTEGER NEWHORDERMIN, NEWHORDERTS, NEWMIN, NEWTS, INDEX, LPAIRLIST(PAIRDISTMAX)
INTEGER, ALLOCATABLE :: MINMAP(:), IDUM(:), TSMAP(:), LOCATIONDB(:)
CHARACTER(LEN=130) FNAME, PDFILE, PLFILE, LPDFILE, LPLFILE
LOGICAL YESNO, NEWMINT
INTEGER MINDB, TSDB
INTEGER :: NMINDBMAX=10
INTEGER :: NTSDBMAX=10
LOGICAL I_OPENED
CHARACTER(LEN=80) I_NAME, I_ACTION

!
! First check for new minima in <pathname>/min.data and <pathname>/points.min
!
NMINOLD=NMIN
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/min.data'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/points.min'
OPEN(UNIT=2,FILE=TRIM(ADJUSTL(FNAME)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*NOPT)
NEWMIN=0
MINDB=0
ALLOCATE(MINMAP(NMINDBMAX))
!
! check for pairlist and pairdist files
!
PLUNIT=-1
PDUNIT=-1
IF (DIJINITT.OR.DIJINITCONTT) THEN
   INQUIRE(FILE='pairlist',EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(A)','mergedb> ERROR *** no pairlist file in current directory'
      STOP
   ENDIF
   INQUIRE(FILE='pairdist',EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(A)','mergedb> ERROR *** no pairdist file in current directory'
      STOP
   ENDIF
   LPLUNIT=GETUNIT()
   OPEN(UNIT=LPLUNIT,FILE='pairlist',STATUS='OLD',POSITION='APPEND')
   LPDUNIT=GETUNIT()
   OPEN(UNIT=LPDUNIT,FILE='pairdist',STATUS='OLD',POSITION='APPEND')
   PDFILE=TRIM(ADJUSTL(PATHNAME)) // '/pairdist'
   PLFILE=TRIM(ADJUSTL(PATHNAME)) // '/pairlist'
   INQUIRE(FILE=PLFILE,EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(A)','mergedb> ERROR *** no pairlist file in merge directory'
      STOP
   ENDIF
   INQUIRE(FILE=PDFILE,EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(A)','mergedb> ERROR *** no pairdist file in merge directory'
      STOP
   ENDIF
   PLUNIT=GETUNIT()
   OPEN(UNIT=PLUNIT,FILE=PLFILE,STATUS='OLD')
   PDUNIT=GETUNIT()
   OPEN(UNIT=PDUNIT,FILE=PDFILE,STATUS='OLD')
ENDIF

DO
   MINDB=MINDB+1
   IF (MINDB.GT.NMINDBMAX) THEN
      ALLOCATE(IDUM(NMINDBMAX))
      IDUM(1:NMINDBMAX)=MINMAP(1:NMINDBMAX)
      DEALLOCATE(MINMAP)
      ALLOCATE(MINMAP(2*NMINDBMAX))
      MINMAP(1:NMINDBMAX)=IDUM(1:NMINDBMAX)
      NMINDBMAX=2*NMINDBMAX
      DEALLOCATE(IDUM)
   ENDIF
   INDEX=NMIN+NEWMIN+1
   IF (INDEX.GT.MAXMIN) CALL MINDOUBLE
   IF (IMFRQT) THEN
      READ(1,*,END=30) EMIN(INDEX),FVIBMIN(INDEX),HORDERMIN(INDEX),IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX),NEGEIG(INDEX)
   ELSE
      READ(1,*,END=30) EMIN(INDEX),FVIBMIN(INDEX),HORDERMIN(INDEX),IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
   END IF
   NEWEMIN=EMIN(INDEX)
   NEWFVIBMIN=FVIBMIN(INDEX)
   NEWHORDERMIN=HORDERMIN(INDEX)
   IF (DIJINITT.OR.DIJINITCONTT) THEN
      IF (INITIALDIST) THEN
         PRINT '(A)','mergedb> ERROR *** mergedb has not been changed to handle INITIALDIST yet'
         STOP
      ENDIF
      READ(PLUNIT,'(10I10)') (LPAIRLIST(J4),J4=1,PAIRDISTMAX)
!     PRINT '(A)','mergedb> read pairlist values:'
!     WRITE(*,'(10I10)') (LPAIRLIST(J4),J4=1,PAIRDISTMAX)
      READ(PDUNIT,'(10G20.10)') (LPAIRDIST(J4),J4=1,PAIRDISTMAX)
   ENDIF
!
!  Read in points and check for agreement with moments of inertia as in setup
!
   READ(2,REC=MINDB) (NEWPOINTSMIN(J2),J2=1,NOPT)  
   LOCALPOINTS(1:NOPT)=NEWPOINTSMIN(1:NOPT)
   IF (AMHT) THEN
      WRITE(*,*)'mergedb> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
   ELSE
      CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,NEWIXMIN,NEWIYMIN,NEWIZMIN)
      IF ((ABS(NEWIXMIN-IXMIN(INDEX)).GT.IDIFFTOL).OR. &
  &    (ABS(NEWIYMIN-IYMIN(INDEX)).GT.IDIFFTOL).OR. &
  &    (ABS(NEWIZMIN-IZMIN(INDEX)).GT.IDIFFTOL)) THEN
         WRITE(*,'(A)') 'mergedb> possible error - principal moments of inertia do not agree with input'
         WRITE(*,'(A,3F20.10)') 'mergedb> values from coordinates: ',NEWIXMIN,NEWIYMIN,NEWIZMIN
         WRITE(*,'(A,3F20.10)') 'mergedb> values from ' // TRIM(ADJUSTL(PATHNAME)) // '/min.data', &
   &                                     IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
                     PRINT '(A,3G20.10)',' NEWEMIN=',NEWEMIN
                     PRINT '(A,2I10)',' MINDB=',MINDB
                      PRINT '(A)','LOCALPOINTS:'
                      PRINT '(3F20.10)',LOCALPOINTS(1:NOPT)
        STOP
      ENDIF
   ENDIF
!
!  Is it a new minimum? Set up MINMAP.
!  Testing the new minima against themselves allows us to remove duplicates from 
!  the database we are merging! 
!
!  DO J2=1,NMIN
   DO J2=1,NMIN+NEWMIN
      IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
         DISTANCE=1.0D100
         READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NOPT)
         CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         IF (.NOT.(ADDPT2.OR.ADDPT3)) THEN
            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.GT.GEOMDIFFTOL)) THEN
               IF (PERMDIST.OR.LPERMDIST) THEN
                  IF(DEBUG) THEN
                     PRINT '(A)',' likely error?'
                     PRINT '(A,3G20.10)',' NEWEMIN,EMIN(J2),DISTANCE=',NEWEMIN,EMIN(J2),DISTANCE
                     PRINT '(A,2I10)',' J2,MINDB=',J2,MINDB
                      PRINT '(A)','LOCALPOINTS:'
                      PRINT '(3F20.10)',LOCALPOINTS(1:NOPT)
                      PRINT '(A)','LOCALPOINTS2:'
                      PRINT '(3F20.10)',LOCALPOINTS2(1:NOPT)
                  ENDIF
!                 STOP !!! DJW
               ENDIF
!
! csw34> Added CHARMM structure dumping for debugging purposes
!
               IF (CHARMMT) THEN
                  PRINT '(A)','Dumping both sets of coordinates to .crd files'
                  CALL CHARMMDUMP(LOCALPOINTS,'localpoints.crd')
                  CALL CHARMMDUMP(LOCALPOINTS2,'localpoints2.crd')
               ENDIF 
            ENDIF
         ENDIF
         
         IF (DISTANCE.LT.GEOMDIFFTOL) THEN
            IF (DEBUG) PRINT '(2(A,I10))','mergedb> minimum ',MINDB,' is database minimum ',J2
            IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-3) THEN
               WRITE(*,'(A,F15.5,A,F15.5)') 'mergedb> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
            ENDIF
            IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
               WRITE(*,'(A,I6,A,I6)') 'mergedb> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
               NEWHORDERMIN=HORDERMIN(J2)
               WRITE(*,'(A,I6)') 'mergedb> using existing value: ',HORDERMIN(J2)
            ENDIF
            MINMAP(MINDB)=J2
            GOTO 130
         ENDIF
      ENDIF
   ENDDO
!
!  If we reach this point we have a new minimum
!
   NEWMIN=NEWMIN+1
   MINMAP(MINDB)=NMIN+NEWMIN
   GPFOLD(NMIN+NEWMIN)=0.0D0
   TOPPOINTER(NMIN+NEWMIN)=0
   IF (DEBUG) PRINT '(2(A,I10))','mergedb> new minimum number ',NMIN+NEWMIN,' number of new minima=',NEWMIN
   IF (DEBUG) WRITE(*,'(A,I10,A)') 'mergedb> new minimum ',NMIN+NEWMIN,&
  &                ' writing parameters to file min.data and points to points.min'
   IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
   WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(INDEX),FVIBMIN(INDEX),HORDERMIN(INDEX),IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
   CALL FLUSH(UMINDATA)
   IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
   WRITE(UMIN,REC=INDEX) (NEWPOINTSMIN(J2),J2=1,NOPT)
!
! pairlist and pairdist if applicable
!
   IF (DIJINITT.OR.DIJINITCONTT) THEN
      WRITE(LPLUNIT,'(10I10)') (LPAIRLIST(J4),J4=1,PAIRDISTMAX)
      WRITE(LPDUNIT,'(10G20.10)') (LPAIRDIST(J4),J4=1,PAIRDISTMAX)
   ENDIF
   
130 CONTINUE
ENDDO
30 NMIN=NMINOLD+NEWMIN
CLOSE(1)
CLOSE(2)
IF (DIJINITT.OR.DIJINITCONTT) THEN
   CLOSE(LPDUNIT)
   CLOSE(LPLUNIT)
ENDIF
IF (PDUNIT.GT.0) CLOSE(PDUNIT)
IF (PLUNIT.GT.0) CLOSE(PLUNIT)
IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
MINDB=MINDB-1
PRINT '(A,I10,A,I10,A)','mergedb> merged ',MINDB,' minima - ',NEWMIN,' are new'
!
!  Now for the transition states.
!
NTSOLD=NTS
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/ts.data'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/points.ts'
OPEN(UNIT=2,FILE=TRIM(ADJUSTL(FNAME)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*NOPT)
NEWTS=0
TSDB=0
NDUMMY=0
ALLOCATE(TSMAP(NTSDBMAX))
DO
   TSDB=TSDB+1
   IF (TSDB.GT.NTSDBMAX) THEN
      ALLOCATE(IDUM(NTSDBMAX))
      IDUM(1:NTSDBMAX)=TSMAP(1:NTSDBMAX)
      DEALLOCATE(TSMAP)
      ALLOCATE(TSMAP(2*NTSDBMAX))
      TSMAP(1:NTSDBMAX)=IDUM(1:NTSDBMAX)
      NTSDBMAX=2*NTSDBMAX
      DEALLOCATE(IDUM)
   ENDIF
   INDEX=NTS+NEWTS+1
   IF (INDEX.GT.MAXTS) CALL TSDOUBLE
   IF (IMFRQT) THEN
      READ(1,*,END=40) &
         ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),PLUS(INDEX),MINUS(INDEX),IXTS(INDEX),IYTS(INDEX),IZTS(INDEX),NEGEIG(INDEX)
   ELSE
      READ(1,*,END=40) ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),PLUS(INDEX),MINUS(INDEX),IXTS(INDEX),IYTS(INDEX),IZTS(INDEX)
   ENDIF
   NEWETS=ETS(INDEX)
   NEWFVIBTS=FVIBTS(INDEX)
   NEWHORDERTS=HORDERTS(INDEX)
   NEWNEGEIG=NEGEIG(INDEX)
!
!  Read in points and check for agreement with moments of inertia as in setup
!
   READ(2,REC=TSDB) (NEWPOINTSTS(J2),J2=1,NOPT)  
   LOCALPOINTS(1:NOPT)=NEWPOINTSTS(1:NOPT)
   IF (AMHT) THEN
      WRITE(*,*)'mergedb> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
   ELSE
      CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,NEWIXTS,NEWIYTS,NEWIZTS)
      IF ((ABS(NEWIXTS-IXTS(INDEX)).GT.IDIFFTOL).OR. &
  &       (ABS(NEWIYTS-IYTS(INDEX)).GT.IDIFFTOL).OR. &
  &       (ABS(NEWIZTS-IZTS(INDEX)).GT.IDIFFTOL)) THEN
         WRITE(*,'(A)') 'mergedb> possible error - principal moments of inertia do not agree with input'
         WRITE(*,'(A,3F20.10)') 'mergedb> values from coordinates: ',NEWIXTS,NEWIYTS,NEWIZTS
         WRITE(*,'(A,3F20.10)') 'mergedb> values from ' // TRIM(ADJUSTL(PATHNAME)) // '/ts.data', &
   &                                        IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
         IF (.NOT.BULKT) STOP
      ENDIF
   ENDIF
!
!  Is it a new ts? Set up TSMAP.
!
!  DO J2=1,NTS
   DO J2=1,NTS+NEWTS
      IF (ABS(NEWETS-ETS(J2)).LT.EDIFFTOL) THEN
         DISTANCE=1.0D100
         READ(UTS,REC=J2) (LOCALPOINTS2(J3),J3=1,NOPT)
         CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         IF (DISTANCE.LT.GEOMDIFFTOL) THEN

            IF (DEBUG) PRINT '(2(A,I10),A,2G20.10)','mergedb> ts ',TSDB,' is database ts ',J2,' energies: ',NEWETS,ETS(J2)
            IF (ABS(NEWFVIBTS-FVIBTS(J2))/FVIBTS(J2).GT.1.0D-3) THEN
               WRITE(*,'(A,F15.5,A,F15.5)') 'mergedb> WARNING, NEWFVIBTS=',NEWFVIBTS,' should be ',FVIBTS(J2)
            ENDIF
            IF (NEWHORDERTS.NE.HORDERTS(J2)) THEN
            WRITE(*,'(A,I6,A,I6)') 'mergedb> ERROR, NEWHORDERTS=',NEWHORDERTS,' should be ',HORDERTS(J2)
                  NEWHORDERTS=HORDERTS(J2)
               WRITE(*,'(A,I6)') 'mergedb> using existing value: ',HORDERTS(J2)
            ENDIF
!
!  Check end minima are consistent as well.
!
            NEWMINT=.FALSE.
            IF (.NOT.(((MINMAP(PLUS(INDEX)).EQ.PLUS(J2)).AND.(MINMAP(MINUS(INDEX)).EQ.MINUS(J2))).OR. &
  &                   ((MINMAP(PLUS(INDEX)).EQ.MINUS(J2)).AND.(MINMAP(MINUS(INDEX)).EQ.PLUS(J2)))) ) THEN
               PRINT '(A,I10)', 'mergedb> Inconsistent minima for old transition state ',J2
               PRINT '(A,2I10)','mergedb> Previous minima: ',PLUS(J2),MINUS(J2)
               PRINT '(A,2I10)','mergedb> Minima for database '//TRIM(ADJUSTL(PATHNAME))//'/ts.data: ', &
  &                             MINMAP(PLUS(INDEX)),MINMAP(MINUS(INDEX))
               NEWMINT=.TRUE.
            ENDIF
            IF ((.NOT.ALLTST).OR.(.NOT.NEWMINT)) THEN
               TSMAP(TSDB)=J2
               NDUMMY=NDUMMY+1
               PRINT '(A)', 'mergedb> Not recounting this transition state'
               GOTO 140
            ELSE
               PRINT '(A)', 'mergedb> ALLTS set - treating this as a new transition state'
            ENDIF
         ENDIF
      ENDIF
   ENDDO
!
!  If we reach this point we have a new ts
!
   NEWTS=NEWTS+1
   TSMAP(TSDB)=NTS+NEWTS
   IF (DEBUG) PRINT '(A,I10,A,I10)','mergedb> new ts number ',NTS+NEWTS,' number of new ts=',NEWTS
   IF (DEBUG) WRITE(*,'(A,I10,A)') 'mergedb> new ts ',NEWTS,' writing parameters to file ts.data and points to points.ts'
   IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
!  PRINT *,'CLOSEFILEST=',CLOSEFILEST
!  INQUIRE (UTSDATA, OPENED=I_OPENED, NAME=I_NAME, ACTION=I_ACTION) 
!  PRINT *,'I_OPENED=',I_OPENED
!  PRINT *,'I_NAME=',I_NAME
!  PRINT *,'I_ACTION=',I_ACTION

   IF (IMFRQT) THEN
      WRITE(UTSDATA,'(2F25.15,3I10,4F20.10)') ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),MINMAP(PLUS(INDEX)),MINMAP(MINUS(INDEX)), &
        &                                          IXTS(INDEX),IYTS(INDEX),IZTS(INDEX),NEGEIG(INDEX)
   ELSE
      WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),MINMAP(PLUS(INDEX)),MINMAP(MINUS(INDEX)), &
        &                                          IXTS(INDEX),IYTS(INDEX),IZTS(INDEX)
   ENDIF
   CALL FLUSH(UTSDATA)
   IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
   WRITE(UTS,REC=INDEX) (NEWPOINTSTS(J2),J2=1,NOPT)
   
140 CONTINUE
ENDDO
40  CONTINUE
NTS=NTSOLD+NEWTS

IF (ADDPT2.OR.ADDPT3) THEN
   DO J1=1,NEWTS
      READ(UTS,REC=NTSOLD+NEWTS-NEWTS+J1) (NEWPOINTSTS(J2),J2=1,NOPT)
      READ(UMIN,REC=PLUS(NTSOLD+NEWTS-NEWTS+J1)) (LPLUS(J2),J2=1,NOPT)
      READ(UMIN,REC=MINUS(NTSOLD+NEWTS-NEWTS+J1)) (LMINUS(J2),J2=1,NOPT)
      IF (ADDPT2) CALL ADDPERM2(NTSOLD+NEWTS-NEWTS+J1,NEWPOINTSTS,LPLUS,LMINUS)
      IF (ADDPT3) CALL ADDPERM3(NTSOLD+NEWTS-NEWTS+J1,NEWPOINTSTS,LPLUS,LMINUS)
   ENDDO
ENDIF

CLOSE(1)
CLOSE(2)
IF (NTS.GT.MAXTS) CALL TSDOUBLE
TSDB=TSDB-1
PRINT '(A,I10,A,I10,A)','mergedb> merged ',TSDB,' ts - ',NEWTS,' are new'
!
!  Finally, check min.A and min.B
!

FNAME=TRIM(ADJUSTL(PATHNAME)) // '/min.A'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
READ(1,*) NMINDB
ALLOCATE(LOCATIONDB(NMINDB)) 
READ(1,*) LOCATIONDB(1:NMINDB) 
CLOSE(1)
NEWMIN=0
OPEN(UNIT=3,FILE='min.A.new',STATUS='UNKNOWN')
DO J1=1,NMINDB
   DO J2=1,NMINA
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONA(J1)) GOTO 55 ! it was already a A minimum
   ENDDO
   NEWMIN=NEWMIN+1
55 CONTINUE
ENDDO
PRINT '(A,I10)','mergedb> number of new A minima: ',NEWMIN
WRITE(3,*) NMINA+NEWMIN
IF (NMINA.GT.0) WRITE(3,'(I10)') LOCATIONA(1:NMINA)
DO J1=1,NMINDB
   PRINT '(3(A,I10))','mergedb> A minimum ',J1,' was number ',LOCATIONDB(J1),' maps to number ',MINMAP(LOCATIONDB(J1))
   DO J2=1,NMINA
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONA(J1)) GOTO 50 ! it was already an A minimum
   ENDDO
   WRITE(3,'(I10)') MINMAP(LOCATIONDB(J1))
50 CONTINUE
ENDDO
CLOSE(3)
DEALLOCATE(LOCATIONDB)

FNAME=TRIM(ADJUSTL(PATHNAME)) // '/min.B'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
READ(1,*) NMINDB
ALLOCATE(LOCATIONDB(NMINDB)) 
READ(1,*) LOCATIONDB(1:NMINDB) 
CLOSE(1)
NEWMIN=0
OPEN(UNIT=3,FILE='min.B.new',STATUS='UNKNOWN')
DO J1=1,NMINDB
   DO J2=1,NMINB
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONB(J1)) GOTO 65 ! it was already a B minimum
   ENDDO
   NEWMIN=NEWMIN+1
65 CONTINUE
ENDDO
WRITE(3,*) NMINB+NEWMIN
PRINT '(A,I10)','mergedb> number of new B minima: ',NEWMIN
IF (NMINB.GT.0) WRITE(3,'(I10)') LOCATIONB(1:NMINB)
DO J1=1,NMINDB
   PRINT '(3(A,I10))','mergedb> B minimum ',J1,' was number ',LOCATIONDB(J1),' maps to number ',MINMAP(LOCATIONDB(J1))
   DO J2=1,NMINB
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONB(J1)) GOTO 60 ! it was already a B minimum
   ENDDO
   WRITE(3,'(I10)') MINMAP(LOCATIONDB(J1))
60 CONTINUE
ENDDO
CLOSE(3)
DEALLOCATE(LOCATIONDB)

STOP
END
