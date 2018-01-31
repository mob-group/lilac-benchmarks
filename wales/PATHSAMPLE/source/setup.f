!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling 
!   and perform kinetic analysis
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Read in databases of A and B minima, calculate partition functions, rate constants
C  and sums of rate constants for all the transition states in the database.
C
C  We need pre-existing databases to specify which minima are A and which are B.
C
      SUBROUTINE SETUP
      USE PORFUNCS
      USE UTILS
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2, STATUS, J3, NDUMMY, NRANDOM, NCOUNT, NMINREMOVE, NTSREMOVE, NMINRETAIN, 
     &        J4, NTSRETAIN, ISTAT, LUNIT, LUNIT2
      DOUBLE PRECISION LOCALPOINTS(NOPT), IXM, IYM, IZM, LOCALPOINTS2(NOPT), DISTANCE, RMAT(3,3), DIST2, DPRAND
      DOUBLE PRECISION PFNORM1, PFNORM2
      DOUBLE PRECISION, ALLOCATABLE :: NEWPFMIN(:)
      INTEGER, ALLOCATABLE :: CANDIDATES(:), MINPREV(:), MINREMOVE(:), TSREMOVE(:), MINRETAIN(:), TSRETAIN(:)
      DOUBLE PRECISION NEWEMIN, NEWIXMIN, NEWIYMIN, NEWIZMIN, NEWFVIBMIN, TSEGUESS, NEWMINCURVE, NEWMINFRQ2,
     &                 TSFVIBGUESS, DUMMY, FRICTIONFAC, GLOBALMIN
      DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0
      LOGICAL DEADTS
      CHARACTER(LEN=80) S1, S2, FNAME
      CHARACTER :: CDUMMY
      INTEGER NEWHORDERMIN
      LOGICAL DUPLICATES_CHECK, TESTOP, TESTNAME
      INTEGER :: LINECOUNT = 1 ! For writing pairfiles with MAKEPAIRS
      INTEGER :: LASTMIN = -1  ! For writing pairfiles with MAKEPAIRS

      IF (CHARMMT.and..not.machine) CALL READREF('input.crd')

      OPEN(UNIT=UMIN,FILE='points.min',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT) 
      OPEN(UNIT=UTS,FILE='points.ts',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT) 

!     INQUIRE(UNIT=UMIN,OPENED=TESTOP,NAMED=TESTNAME)
!     PRINT *,'UNIT UMIN is ',UMIN
!     PRINT *,'opened is ',TESTOP
!     PRINT *,'named is ',TESTNAME

!     INQUIRE(UNIT=UTS,OPENED=TESTOP,NAMED=TESTNAME)
!     PRINT *,'UNIT UTS is ',UTS
!     PRINT *,'opened is ',TESTOP
!     PRINT *,'named is ',TESTNAME

      IF (FROMLOWESTT) THEN
         IF (.NOT.NOFRQS) THEN
            WRITE(*,'(A)')'setup> FROMLOWEST does not calculate frequencies; only use with NOFRQS'
            STOP
         ENDIF
         WRITE(*,'(A)')'setup> writing points.min file and min.data from lowest file'
         OPEN(UNIT=1,FILE='lowest',STATUS='UNKNOWN')
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE='points.min.new',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT)
         LUNIT2=GETUNIT()
         OPEN(UNIT=LUNIT2,FILE='min.data.new',STATUS='UNKNOWN')
         NDUMMY=0
!         FVIBMIN(1)=4.6757541330
         FVIBMIN(1)=1.000000000000000
         HORDERMIN(1)=1
         DO
            NDUMMY=NDUMMY+1
            READ(1,*,END=776)
! csw34> The GMIN lowest file usually includes atom names. If for any reason it does not,
!        the 'NOLABELS' option can be used to read just the coordinates
            IF (NOLABELST) THEN
               READ(1,'(A25,F20.10)') CDUMMY,EMIN(1)
               READ(1,*) (LOCALPOINTS(J2),J2=1,NOPT)
            ELSE
               READ(1,'(A25,F20.10,A71)') CDUMMY,EMIN(1),CDUMMY
               DO J1=1,NATOMS
                  READ(1,*) CDUMMY,LOCALPOINTS(3*J1-2),LOCALPOINTS(3*J1-1),LOCALPOINTS(3*J1)
               ENDDO
            ENDIF
            CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,IXM,IYM,IZM)
            WRITE(LUNIT,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
            WRITE(LUNIT2,'(2F25.15,I6,3F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXM,IYM,IZM
         ENDDO
776      NDUMMY=NDUMMY-1
         WRITE(*,'(A,I10)')'setup> number of minima written=',NDUMMY
         CLOSE(LUNIT)
         CLOSE(LUNIT2)
         CLOSE(1)
         STOP
      ENDIF

      IF (EXTRACTMINT) THEN
         OPEN(UNIT=1,FILE='extractedmin', STATUS='UNKNOWN')
         IF (EXTRACTMINFILET) THEN
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='extractmin',STATUS='OLD')
            DO
              READ(LUNIT,*,END=177) WHICHMIN
              READ(UMIN,REC=WHICHMIN) (LOCALPOINTS(J2),J2=1,NOPT)
              IF (MLP3T.OR.MLPB3T) THEN
                  WRITE(1,'(F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ELSE
                  WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ENDIF
            ENDDO
177         CLOSE(LUNIT)
         ELSEIF (WHICHMIN.EQ.-123) THEN
            NDUMMY=0
            PRINT '(A)', 'setup> rewriting binary points.min file from extractedmin file'
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='points.min.new',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT) 
            DO 
               NDUMMY=NDUMMY+1
               READ(1,*,END=777) (LOCALPOINTS(J2),J2=1,NOPT)
               WRITE(LUNIT,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
            ENDDO
777         NDUMMY=NDUMMY-1
            PRINT '(A,I10)','setup> number of minima extracted=',NDUMMY
            CLOSE(LUNIT)
         ELSEIF (WHICHMIN.LE.0) THEN
            NDUMMY=0
            PRINT '(A)', 'setup> extracting all minima '
            DO 
               NDUMMY=NDUMMY+1
               READ(UMIN,REC=NDUMMY,ERR=877) (LOCALPOINTS(J2),J2=1,NOPT)
               IF (MLP3T.OR.MLPB3T) THEN
                  WRITE(1,'(F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ELSE
                  WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ENDIF
            ENDDO
877         NDUMMY=NDUMMY-1
            PRINT '(A,I10)','setup> number of minima extracted=',NDUMMY
         ELSE
            PRINT '(A,I10)', 'setup> extracting minimum ',WHICHMIN
            READ(UMIN,REC=WHICHMIN) (LOCALPOINTS(J2),J2=1,NOPT)

           IF (AMHT) THEN
              CALL AMHDUMP(LOCALPOINTS,'amhmin.pdb')
  
              IF (AMHQT)THEN
                 CALL AMHQ(WHICHMIN)
              ENDIF
  
              IF (AMHQENGMINT)THEN
                 CALL AMHQENGMIN(WHICHMIN)
              ENDIF

              IF (AMHQCONTT)THEN
                 CALL AMHQCONT(WHICHMIN,QCONTCUT)
              ENDIF
  
              IF (AMHRMSDT)THEN
                 CALL AMHRMSD(WHICHMIN)
             ENDIF
  
              IF (AMHRELQT)THEN
                 CALL AMHRELQ(QRELONE, QRELTWO)
              ENDIF
  
              IF (AMH_RELCOT)THEN
                CALL AMH_RELCO(WHICHMIN, RELCOCUT)
              ENDIF
  
              IF (AMHALLATOMMINT)THEN
                 CALL AMHALLATOMMIN
             ENDIF
           ELSE
!
! Write the extractedmin file to UNIT=1
!
              IF (MLP3T.OR.MLPB3T) THEN
                 WRITE(1,'(F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
              ELSE
                 WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
              ENDIF
!
! Write the extractedmin.crd file if using CHARMM
!
               IF (CHARMMT) CALL CHARMMDUMP(LOCALPOINTS,'extractedmin.crd')
!
! csw34> Write the extractedmin.rst file if using AMBER/NAB
!
               IF (AMBERT) THEN
                  OPEN(UNIT=10,FILE='extractedmin.rst', STATUS='UNKNOWN')
                  WRITE(10,'(20A4)') 'MOL'
! Formats used come from the AMBER routine minrit
                  IF (NATOMS.GT.100000) THEN
                     WRITE(10,'(I10)') NATOMS
                  ELSE
                     WRITE(10,'(I5)') NATOMS
                  ENDIF
                  WRITE(10,'(6f12.7)') (LOCALPOINTS(J2),J2=1,NOPT)
                  CLOSE(10)
               ENDIF

           ENDIF

         ENDIF
         CLOSE(1)
         STOP
      ENDIF

      IF (EXTRACTTST) THEN
         OPEN(UNIT=1,FILE='extractedts', STATUS='UNKNOWN')
         IF (EXTRACTTSFILET) THEN
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='extractts',STATUS='OLD')
            DO
              READ(LUNIT,*,END=166) WHICHTS
              READ(UTS,REC=WHICHTS) (LOCALPOINTS(J2),J2=1,NOPT)
              IF (MLP3T.OR.MLPB3T) THEN
                  WRITE(1,'(F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ELSE
                  WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ENDIF
            ENDDO
166         CLOSE(LUNIT)

         ELSEIF (WHICHTS.EQ.-123) THEN
            NDUMMY=0
            PRINT '(A)', 'setup> rewriting binary points.ts file from extractedts file'
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='points.ts.new',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT) 
            DO
               NDUMMY=NDUMMY+1
               READ(1,*,END=778) (LOCALPOINTS(J2),J2=1,NOPT)
               WRITE(LUNIT,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
            ENDDO
778         NDUMMY=NDUMMY-1
            CLOSE(LUNIT)
            PRINT '(A,I10)','setup> number of ts extracted=',NDUMMY
         ELSEIF (WHICHTS.LE.0) THEN
            NDUMMY=0
            PRINT '(A)', 'setup> extracting all ts '
            DO
               NDUMMY=NDUMMY+1
               READ(UTS,REC=NDUMMY,ERR=878) (LOCALPOINTS(J2),J2=1,NOPT)
               IF (MLP3T.OR.MLPB3T) THEN
                  WRITE(1,'(F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ELSE
                  WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
               ENDIF
            ENDDO
878         NDUMMY=NDUMMY-1
            PRINT '(A,I10)','setup> number of ts extracted=',NDUMMY
         ELSE
            PRINT '(A,I10)', 'setup> extracting ts ',WHICHTS
            READ(UTS,REC=WHICHTS) (LOCALPOINTS(J2),J2=1,NOPT)

!           IF (AMHT) THEN
!              CALL AMHDUMP(LOCALPOINTS,'amhts.pdb')
!           ELSE
!
! Write the extractedmin file to UNIT=1
!
               WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,NOPT)
!
! Write the extractedmin.crd file when using CHARMM
!
               IF (CHARMMT) CALL CHARMMDUMP(LOCALPOINTS,'extractedts.crd')
!
! csw34> Write the extractedts.rst file if using AMBER/NAB
!
               IF (AMBERT) THEN
                  OPEN(UNIT=10,FILE='extractedts.rst', STATUS='UNKNOWN')
                  WRITE(10,'(20A4)') 'MOL'
! Formats used come from the AMBER routine minrit
                  IF (NATOMS.GT.100000) THEN
                     WRITE(10,'(I10)') NATOMS
                  ELSE
                     WRITE(10,'(I5)') NATOMS
                  ENDIF
                  WRITE(10,'(6f12.7)') (LOCALPOINTS(J2),J2=1,NOPT)
                  CLOSE(10)
               ENDIF
!          ENDIF
         ENDIF
         CLOSE(1)
         STOP
      ENDIF
!
!  In this case we try to connect two minima, whose coordinates are in
!  files odata.start and odata.finish for a DIJINITSTARTT run. We set up min.data as 
!  well as min.A and min.B, with entries 1 1 and 1 2.
!
!  DIJINITCONTT specifies a continuation of an initial path run.
!
      IF (DIJINITT.OR.DIJINITFLYT) THEN
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (DIJINITCONTT) THEN
            IF (.NOT.YESNO) THEN
               PRINT '(A)','setup> ERROR - min.data must exist for a DIJINIT continuation run'
               STOP
            ENDIF
         ELSEIF (DIJINITSTARTT) THEN
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.data min.data.save')
            PRINT '(A)','setup> Creating new min.data file'
            OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN')
            INQUIRE(FILE='pairs.data',EXIST=YESNO)
            IF (YESNO) THEN
               PRINT '(A)','setup> Moving old pairs.data file'
               CALL MYSYSTEM(STATUS,DEBUG,'mv pairs.data pairs.data.old')
            ENDIF
            INQUIRE(FILE='commit.data',EXIST=YESNO)
            IF (YESNO) THEN
               PRINT '(A)','setup> Moving old commit.data file'
               CALL MYSYSTEM(STATUS,DEBUG,'mv commit.data commit.data.old')
            ENDIF
            INQUIRE(FILE='min.data.info.1',EXIST=YESNO)
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'rm min.data.info.1')
            CALL MYSYSTEM(STATUS,DEBUG,'cp odata.start odata.1')
            IF (CHARMMT) CALL MYSYSTEM(STATUS,DEBUG,'cp start.crd input.crd.1')
            IF (AMHT) CALL MYSYSTEM(STATUS,DEBUG,'cp AMHstart start.1')
            IF (DEBUG .AND. (EXEC.EQ.'UNDEFINED')) THEN
               WRITE(*,*) "setup> ERROR: Executable is undefined"
               STOP
            ENDIF
            CALL MYSYSTEM(STATUS,DEBUG,TRIM(ADJUSTL(EXEC))//' 1 > output.start')
            INQUIRE(FILE='min.data.info.1',EXIST=YESNO)
            IF (YESNO) THEN
               OPEN(UNIT=1,FILE='min.data.info.1',STATUS='OLD')
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*) EMIN(1),FVIBMIN(1),HORDERMIN(1),IXMIN(1),IYMIN(1),IZMIN(1),MINCURVE(1),MINFRQ2(1)
                  WRITE(UMINDATA,'(2F25.15,I6,5F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXMIN(1),
     &                                                      IYMIN(1),IZMIN(1),MINCURVE(1),MINFRQ2(1)
               ELSE
                  READ(1,*) EMIN(1),FVIBMIN(1),HORDERMIN(1),IXMIN(1),IYMIN(1),IZMIN(1)
                  WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXMIN(1), IYMIN(1), IZMIN(1)
               ENDIF
               READ(1,*) (LOCALPOINTS(J2),J2=1,NOPT)
               WRITE(UMIN,REC=1) (LOCALPOINTS(J2),J2=1,NOPT)
               CLOSE(1)
            ELSE
               PRINT *, 'setup> ERROR - no min.data.info.1 found - check OPTIM output in output.start'
               STOP        
            ENDIF
            INQUIRE(FILE='min.data.info.2',EXIST=YESNO)
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'rm min.data.info.2')
            CALL MYSYSTEM(STATUS,DEBUG,'cp odata.finish odata.2')
            IF (CHARMMT) CALL MYSYSTEM(STATUS,DEBUG,'cp finish.crd input.crd.2')
            IF (AMHT) CALL MYSYSTEM(STATUS,DEBUG,'cp AMHfinish start.2')
            CALL MYSYSTEM(STATUS,DEBUG,TRIM(ADJUSTL(EXEC))//' 2 > output.finish')
            INQUIRE(FILE='min.data.info.2',EXIST=YESNO)
            IF (YESNO) THEN
               OPEN(UNIT=1,FILE='min.data.info.2',STATUS='OLD')
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*) EMIN(2),FVIBMIN(2),HORDERMIN(2),IXMIN(2),IYMIN(2),IZMIN(2),MINCURVE(2),MINFRQ2(2)
                  WRITE(UMINDATA,'(2F25.15,I6,5F20.10)') EMIN(2), FVIBMIN(2), HORDERMIN(2),IXMIN(2), 
     &                                                   IYMIN(2),IZMIN(2),MINCURVE(2),MINFRQ2(2)
               ELSE
                  READ(1,*) EMIN(2),FVIBMIN(2),HORDERMIN(2),IXMIN(2),IYMIN(2),IZMIN(2)
                  WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(2), FVIBMIN(2), HORDERMIN(2),IXMIN(2), IYMIN(2), IZMIN(2)
               ENDIF
               CLOSE(UMINDATA)
               READ(1,*) (LOCALPOINTS(J2),J2=1,NOPT)
               WRITE(UMIN,REC=2) (LOCALPOINTS(J2),J2=1,NOPT)
               CLOSE(1)
            ELSE
               PRINT *, 'setup> ERROR - no min.data.info.2 found - check OPTIM output in output.finish'
               STOP
            ENDIF
            INQUIRE(FILE='min.A',EXIST=YESNO)
            IF (DIJINITSTARTT) THEN
               IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.A min.A.save')
               PRINT '(A)','setup> Creating new min.A file'
            ELSEIF (DIJINITCONTT) THEN
               IF (.NOT.YESNO) THEN
                  PRINT '(A)','setup> ERROR - min.A must exist for a DIJINIT continuation run'
                  STOP
               ENDIF
            ENDIF
            OPEN(1,FILE='min.A',STATUS='UNKNOWN')
            WRITE(1,'(I10)') 1
            WRITE(1,'(I10)') 1
            CLOSE(1)
            INQUIRE(FILE='min.B',EXIST=YESNO)
            IF (DIJINITSTARTT) THEN
               IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.B min.B.save')
               PRINT '(A)','setup> Creating new min.B file'
            ELSEIF (DIJINITCONTT) THEN
               IF (.NOT.YESNO) THEN
                  PRINT '(A)','setup> ERROR - min.B must exist for a DIJINIT continuation run'
                  STOP
               ENDIF
            ENDIF
            OPEN(1,FILE='min.B',STATUS='UNKNOWN')
            WRITE(1,'(I10)') 1
            WRITE(1,'(I10)') 2
            CLOSE(1)
            PRINT '(A)','setup> initial OPTIM jobs run for odata.start and odata.finish'
            IF (DUMMYTST) THEN
               READ(UMIN,REC=1) (LOCALPOINTS(J2),J2=1,NOPT)
               READ(UMIN,REC=2) (LOCALPOINTS2(J3),J3=1,NOPT)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                                                  DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               MINDISTMIN(1)=DISTANCE
               MINDISTMIN(2)=DISTANCE
C
C  Must create an entry in ts.data in this case.
C  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
C
               INQUIRE(FILE='ts.data',EXIST=YESNO)
               IF (YESNO) THEN
                  PRINT '(A)','setup> - file ts.data already exists. Copying to ts.data.save'
                  CALL MYSYSTEM(STATUS,DEBUG,'mv ts.data ts.data.save')
               ENDIF
               OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW') 

               IF (IMFRQT) THEN
                  PRINT '(A)',"setup> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                  STOP
               ENDIF
               WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') TSEGUESS(EMIN(1),EMIN(2),MINCURVE(1),MINCURVE(2),DISTANCE), 
     &               TSFVIBGUESS(EMIN(1),EMIN(2),FVIBMIN(1),FVIBMIN(2),MINFRQ2(1),MINFRQ2(2),NATOMS),1,1,2,1.0D0,1.0D0,1.0D0

               CLOSE(UTSDATA)
            ENDIF
         ENDIF
      ENDIF

      IF (STARTFROMPATH) THEN
C
C  Get all the necessary information about the A and B minima from the <PATHNAME> file.
C  The assumption is that we have just one A minimum and one B in this case.
C  Use GETNEWPATH or GETALLPATHS to do the bookkeeping.
C  Detect presence of existing min.data file and STOP if found to prevent overwriting.
C
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.data already exists. Will not overwrite.'
            STOP
         ENDIF
         CALL MYSYSTEM(STATUS,DEBUG,'cp ' // TRIM(ADJUSTL(PATHNAME)) // ' path.info')
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='NEW')
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file ts.data already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW')
         IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
!        IF (STARTTRIPLES) THEN
            CALL GETALLPATHS
!        ELSE
!           CALL GETNEWPATH(0,0)
!        ENDIF
         INQUIRE(FILE='min.A',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.A already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(1,FILE='min.A',STATUS='NEW') 
         WRITE(1,'(I10)') 1
         WRITE(1,'(I10)') STARTMINA
         CLOSE(1)
         INQUIRE(FILE='min.B',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.B already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(1,FILE='min.B',STATUS='NEW') 
         WRITE(1,'(I10)') 1
         WRITE(1,'(I10)') STARTMINB
         CLOSE(1)
         CLOSE(UMINDATA)
         CLOSE(UTSDATA)
         PRINT '(A,I5,A,I5,A)','setup> The unique A and B minima are ',STARTMINA,' and ',STARTMINB,' respectively'
         IF (NATTEMPT.LE.0) STOP
      ENDIF
   
      ! here we read in minima from min.A and min.B
      INQUIRE(FILE='min.A',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(1,FILE='min.A',STATUS='OLD')
         READ(1,*) NMINA
         ALLOCATE(LOCATIONA(NMINA)) ! not deallocated
         READ(1,*) LOCATIONA(1:NMINA)
         CLOSE(1)
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT.OR.ADDMINXYZT)) THEN
         PRINT '(A)','setup> ERROR - no min.A file'
         STOP
      ELSE
         NMINA=0
         ALLOCATE(LOCATIONA(NMINA))
      ENDIF
      IF (DEBUG) THEN
         PRINT '(A,I10,A)','setup> there are ',NMINA,' A minima at locations:'
         PRINT '(10I10)',LOCATIONA(1:NMINA)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I10,A)') 'setup> locations read for ',NMINA,' min of type A'
      ! check there are not duplicates in min.A
      IF (DUPLICATES_CHECK(LOCATIONA, NMINA)) THEN
          write(*,'(A)') "setup> ERROR - there are duplicate values in min.A"
          stop
      ENDIF

      INQUIRE(FILE='min.B',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(1,FILE='min.B',STATUS='OLD')
         READ(1,*) NMINB
         ALLOCATE(LOCATIONB(NMINB)) ! not deallocated
         READ(1,*) LOCATIONB(1:NMINB)
         CLOSE(1)
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT.OR.ADDMINXYZT)) THEN
         PRINT '(A)','setup> ERROR - no min.B file'
         STOP
      ELSE
         NMINB=0
         ALLOCATE(LOCATIONB(NMINB)) ! not deallocated
      ENDIF
      IF (DEBUG) THEN
         PRINT '(A,I10,A)','setup> there are ',NMINB,' B minima at locations:'
         PRINT '(10I6)',LOCATIONB(1:NMINB)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I10,A)') 'setup> locations read for ',NMINB,' min of type B'
      IF (DUPLICATES_CHECK(LOCATIONB, NMINB)) THEN
         write(*,'(A)') "setup> ERROR - there are duplicate values in min.B"
         stop
      ENDIF
C
C  Load the minima.
C
      INQUIRE(FILE='min.data',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='OLD')
         NMIN=0
         J1=0
         DO 
            J1=J1+1
            IF (J1.GT.MAXMIN) CALL MINDOUBLE
            IF (DUMMYTST.AND.LOWESTFRQT) THEN
               READ(UMINDATA,*,END=30) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1),MINCURVE(J1),MINFRQ2(J1)
            ELSE
               READ(UMINDATA,*,END=30) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
            ENDIF
            NMIN=NMIN+1
         ENDDO
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT.OR.ADDMINXYZT)) THEN
         PRINT '(A)','setup> ERROR - no min.data file'
         STOP
      ELSE
         NMIN=0
      ENDIF

30    IF (YESNO) CLOSE(UMINDATA) ! SAT need to reopen this file
      IF (PRINTT) WRITE(*,'(A,I10,A)') 'setup> parameters read for ',NMIN,' min'
!     IF (DEBUG) WRITE(*,'(I6,2F17.7,I6,3F15.5)') (J1,EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),
!    1                                         IXMIN(J1),IYMIN(J1),IZMIN(J1),J1=1,NMIN)
      DO J1=1,NMINA
         IF ((LOCATIONA(J1).GT.NMIN).AND.(.NOT.STARTFROMPATH)) THEN
            PRINT '(3(A,I8))','setup> ERROR - A minimum ',J1,' is number ',LOCATIONA(J1),' but total minima=',NMIN
            STOP
         ENDIF
      ENDDO
      DO J1=1,NMINB
         IF ((LOCATIONB(J1).GT.NMIN).AND.(.NOT.STARTFROMPATH)) THEN
            PRINT '(3(A,I8))','setup> ERROR - B minimum ',J1,' is number ',LOCATIONB(J1),' but total minima=',NMIN
            STOP
         ENDIF
      ENDDO
      IF (CALCORDERT) THEN
         CALL CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
         STOP
      ENDIF
      IF (NMIN.GT.0) THEN
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE") ! read used in Dijkstra
      ENDIF
      IF ((NTS.EQ.0).AND.MERGEDBT) THEN
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='UNKNOWN') ! for MERGEDB startup
      ENDIF
      IF ((NMIN.EQ.0).AND.(READMINT.OR.MERGEDBT.OR.READMINT.OR.ADDMINXYZT)) THEN
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION="APPEND")
      ENDIF
C
C  Check that the necessary coordinates are in fact present. 
C
      IF ((.NOT.NOPOINTS).AND.(NATTEMPT.GT.0)) THEN
          IF (AMHT) THEN
            WRITE(*,*) 'setup> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
          ELSE
           DO J1=1,NMIN
            READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,NOPT)
            CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,IXM,IYM,IZM)
!           IF (DEBUG) WRITE(*,'(2F20.10,I6,3F20.10)') EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXM,IYM,IZM
            IF ((ABS(IXM-IXMIN(J1)).GT.IDIFFTOL).OR.
     1          (ABS(IYM-IYMIN(J1)).GT.IDIFFTOL).OR.
     1          (ABS(IZM-IZMIN(J1)).GT.IDIFFTOL)) THEN
               WRITE(*,'(A,I10)') 'setup> WARNING - principal moments of inertia do not agree with input for minimum ',J1
               WRITE(*,'(A,3F20.10)') 'setup> values from coordinates: ',IXM,IYM,IZM
               WRITE(*,'(A,2F20.10,I6,3F20.10)') 'setup> values from min.data: ', 
     &                     EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
               IXMIN(J1)=IXM
               IYMIN(J1)=IYM
               IZMIN(J1)=IZM
!              STOP  
            ENDIF
           ENDDO
          ENDIF
         IF (PRINTT) WRITE(*,'(A,I10,A)') 'setup> points for ',NMIN,' minima read from file points.min'
      ENDIF
!
! Read data for minima from min.data.info style file. Multiple minima are allowed.
!
      IF (READMINT) THEN
         INQUIRE(FILE=TRIM(ADJUSTL(MINNAME)),EXIST=YESNO)
         IF (YESNO) THEN
            OPEN(UNIT=1,FILE=TRIM(ADJUSTL(MINNAME)),STATUS='OLD')
            IF(MAKEPAIRS) OPEN(UNIT=25,FILE=TRIM(ADJUSTL(MAKEPAIRSFILE)),STATUS='NEW')
            DO 
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*,END=130) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN,NEWMINCURVE,NEWMINFRQ2
               ELSE
                  READ(1,*,END=130) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN
               ENDIF
               READ(1,*) (LOCALPOINTS(J2),J2=1,NOPT)
!
! Must check it is not an old minimum!
!
               DO J2=1,NMIN
                  DISTANCE=1.0D100
                  IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
                     READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NOPT)
                     CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                                DIST2,RIGIDBODY,RMAT,.FALSE.)
                  ENDIF
      
                  IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
                     PRINT '(A,I10)','setup> minimum is database minimum ',J2
                     IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
                        WRITE(*,'(A,F15.5,A,F15.5)') 'setup> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
                     ENDIF
                     IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                        WRITE(*,'(A,I10,A,I10)') 'setup> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                        NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                        WRITE(*,'(A,I10)') 'setup> using maximum value: ',NEWHORDERMIN
                     ENDIF
                     IF (MAKEPAIRS) THEN
                        IF(J2 .NE. LASTMIN) THEN  ! Don't want to try to connect a minimum to itself
                            WRITE(25,*) LINECOUNT, NEWEMIN, J2
                            WRITE(25,*) LINECOUNT+1, "0.0 0.0"
                            LINECOUNT = LINECOUNT+2
                        ENDIF
                        LASTMIN = J2
                     ENDIF
                     GOTO 140
                  ENDIF
               ENDDO
               NMIN=NMIN+1
               IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
               EMIN(NMIN)=NEWEMIN
               FVIBMIN(NMIN)=NEWFVIBMIN
               HORDERMIN(NMIN)=NEWHORDERMIN
               IXMIN(NMIN)=NEWIXMIN
               IYMIN(NMIN)=NEWIYMIN
               IZMIN(NMIN)=NEWIZMIN
               IF (DUMMYTST.AND.LOWESTFRQT) MINCURVE(NMIN)=NEWMINCURVE
               IF (DUMMYTST.AND.LOWESTFRQT) MINFRQ2(NMIN)=NEWMINFRQ2
               WRITE(*,'(A,I10,A)') 'setup> new minimum ',NMIN,
     &                ' writing parameters to file min.data and points to points.min'
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  WRITE(UMINDATA,'(2F25.15,I6,5F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN), MINCURVE(NMIN),MINFRQ2(NMIN)
               ELSE
                  WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
               ENDIF
               IF (MAKEPAIRS) THEN
                   WRITE(25,*) LINECOUNT, EMIN(NMIN), NMIN
                   WRITE(25,*) LINECOUNT+1, "0.0 0.0"
                   LINECOUNT = LINECOUNT+2
                   LASTMIN = NMIN
               ENDIF
               CALL FLUSH(UMINDATA)
               WRITE(UMIN,REC=NMIN) (LOCALPOINTS(J2),J2=1,NOPT)
140            CONTINUE
            ENDDO
130         CLOSE(1)
         ELSE
            PRINT '(A)','setup> ERROR - no file ',TRIM(ADJUSTL(MINNAME))
         ENDIF
      ENDIF

      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
      IF (MAKEPAIRS) CLOSE(25)
C
C  Set all FVIBMIN to ln(2 pi) if NOFRQS is true for consistency.
C  May be useful for running REGROUPFREE on the basis of potential energy only.
C  Won;t work with FREEPAIRS unless the OPTIM runs are also run with NOFRQS.
C
!      IF (NOFRQS) FVIBMIN(1:NMIN)=4.675754133D0 ! 2 ln(2pi) +1
      IF (NOFRQS) FVIBMIN(1:NMIN) = 1.D0
C
C  Calculate partition functions for minima. Note that the total partition function
C  is not needed, just the totals for A and B. Since A and B sets are fixed here
C  we don;t need to change the totals.
C
      PFMEAN=0.0D0
      PFMEAN=-HUGE(1.0D0)
      PFNORM1=0.0D0 ! use this to calculate ratios without the pe factor
      PFNORM2=0.0D0 ! use this to calculate ratios with the pe factor
      
      GLOBALMIN=0.0D0
      IF (RELATIVEET) THEN
         GLOBALMIN=HUGE(1.0D0)
         DO J1=1,NMIN
            IF (EMIN(J1).LT.GLOBALMIN) GLOBALMIN=EMIN(J1)
         ENDDO
         PRINT '(A,G20.10)','setup> Lowest PE minimum=',GLOBALMIN
         PRINT '(A,G20.10)','setup> Shifting all potential energies relative to zero for the lowest minimum'
         DO J1=1,NMIN
            EMIN(J1)=EMIN(J1)-GLOBALMIN
         ENDDO
         PRINT '(A,G20.10)','setup> Shifting transition state threshold for rejection relative to the global minimum'
         TSTHRESH=TSTHRESH-GLOBALMIN
      ENDIF

      IF (ENSEMBLE.EQ.'T') THEN
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A,G20.10)','setup> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
!
! Note that common factors are normally omitted from the canonical partition
! function. For example, the temperature^kappa
!
         DO J1 = 1,NMIN
            PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
            IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
         ENDDO
         DO J1 = 1,NMIN
            PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)) - PFMEAN
            PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))+ FVIBMIN(1)/2.0D0 + LOG(1.0D0*HORDERMIN(1)))
            PFNORM2=PFNORM2+EXP(PFMIN(J1)-PFMIN(1))
         ENDDO
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         DO J1 = 1,NMIN
            IF (TOTALE.GT.EMIN(J1)) THEN
               PFMIN(J1) = (KAPPA-1)*LOG(TOTALE-EMIN(J1)) - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
               PFMEAN=PFMEAN+PFMIN(J1)
               PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)) + FVIBMIN(1)/2.0D0 + LOG(1.0D0*HORDERMIN(1)))
               PFNORM2=PFNORM2+EXP(PFMIN(J1)-PFMIN(1))
               IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
            ELSE
               PFMIN(J1) = -1.0D250
            ENDIF
         ENDDO
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF
!     IF (DEBUG) THEN
!        WRITE(*,'(A,3G20.10)') 'setup> largest ln Z,PFNORM1,PFNORM2=',PFMEAN,PFNORM1,PFNORM2
!        WRITE(*,'(A)') '        V-V_min    pg order     high T/E prob       Peq'
!        DO J1=1,NMIN
!           WRITE(*,'(F20.10,I6,2G20.10)') EMIN(J1),HORDERMIN(J1), 
!    &                    EXP(-FVIBMIN(J1)/2.0D0-LOG(1.0D0*HORDERMIN(J1))-LOG(PFNORM1)+ FVIBMIN(1)/2.0D0 + LOG(1.0D0*HORDERMIN(1))), 
!    &                    EXP(PFMIN(J1)-PFMIN(1)-LOG(PFNORM2))
!        ENDDO
!     ENDIF
!     DO J1=1,NMIN
!        PFMIN(J1) = PFMIN(J1) - PFMEAN
!     ENDDO

      PFTOTALB=0.0D0
      DO J1=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J1))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0.0D0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J1=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J1))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0.0D0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))
!
!     Optional change of reactant minima set via reweighting.
!
      IF (REWEIGHTT) THEN
         ALLOCATE(CANDIDATES(NMIN))
         IF (DIRECTION.EQ.'AB') THEN
            ALLOCATE(NEWPFMIN(NMINB))
            NEWPFMIN(1:NMINB)=0.0D0
!
!  Select NRWREACTANT minima from the B set according to the required weights in RWPROB
!
            PFTOTALB=0.0D0
            DO J1=1,NRWBINS
               IF (NINT(NRWREACTANT*RWPROB(J1)).EQ.0) CYCLE
               NCOUNT=0
               DO J2=1,NMINB ! identify minima in the required energy range
                  IF ((EMIN(LOCATIONB(J2)).GE.RWEMIN+(J1-1)*RWBINWIDTH).AND.
     &                (EMIN(LOCATIONB(J2)).LE.RWEMIN+J1*RWBINWIDTH)) THEN
                      NCOUNT=NCOUNT+1
                      CANDIDATES(NCOUNT)=J2
                  ENDIF
               ENDDO
               PRINT '(3(A,I8),A,G20.10)','setup> number of B minima in energy bin ',J1,' is ',NCOUNT,' number needed=',
     &                           NINT(NRWREACTANT*RWPROB(J1)),' probability=',RWPROB(J1)
               IF (NCOUNT.EQ.0) STOP
               DO J2=1,NINT(NRWREACTANT*RWPROB(J1))
                  NRANDOM=NCOUNT*DPRAND()+1
                  PRINT '(3(A,I8))','setup> selecting B minimum number ',CANDIDATES(NRANDOM),
     &                              ' location ',LOCATIONB(CANDIDATES(NRANDOM)),' for the product set'
                  NEWPFMIN(CANDIDATES(NRANDOM))=NEWPFMIN(CANDIDATES(NRANDOM))+1.0D0
               ENDDO
               PFTOTALB=PFTOTALB+NINT(NRWREACTANT*RWPROB(J1))
            ENDDO
            PFTOTALB=LOG(PFTOTALB) ! partition functions are stored as log's
            NCOUNT=0
            DO J1=1,NMINB
               IF (NEWPFMIN(J1).GT.0.0D0) THEN
                  NCOUNT=NCOUNT+1
                  LOCATIONB(NCOUNT)=LOCATIONB(J1)
                  PFMIN(LOCATIONB(NCOUNT))=LOG(NEWPFMIN(J1))
                  PRINT '(A,I8,A,G20.10)','setup> relative weight for reactant minimum ',LOCATIONB(NCOUNT),' is ',
     &                        EXP(PFMIN(LOCATIONB(NCOUNT))-PFTOTALB)
               ENDIF
            ENDDO
            NMINB=NCOUNT
            PRINT '(A,I8,A)','setup> there are now ',NMINB,' minima of type B'
         ELSE
            ALLOCATE(NEWPFMIN(NMINA))
            NEWPFMIN(1:NMINA)=0.0D0
!
!  Select NRWREACTANT minima from the A set according to the required weights in RWPROB
!
            PFTOTALA=0.0D0
            DO J1=1,NRWBINS
               IF (NINT(NRWREACTANT*RWPROB(J1)).EQ.0) CYCLE
               NCOUNT=0
               DO J2=1,NMINA ! identify minima in the required energy range
                  IF ((EMIN(LOCATIONA(J2)).GE.RWEMIN+(J1-1)*RWBINWIDTH).AND.
     &                (EMIN(LOCATIONA(J2)).LE.RWEMIN+J1*RWBINWIDTH)) THEN
                      NCOUNT=NCOUNT+1
                      CANDIDATES(NCOUNT)=J2
                  ENDIF
               ENDDO
               PRINT '(3(A,I8),A,G20.10)','setup> number of A minima in energy bin ',J1,' is ',NCOUNT,' number needed=',
     &                           NINT(NRWREACTANT*RWPROB(J1)),' probability=',RWPROB(J1)
               IF (NCOUNT.EQ.0) STOP
               DO J2=1,NINT(NRWREACTANT*RWPROB(J1))
                  NRANDOM=NCOUNT*DPRAND()+1
                  PRINT '(3(A,I8))','setup> selecting A minimum number ',CANDIDATES(NRANDOM),
     &                              ' location ',LOCATIONA(CANDIDATES(NRANDOM)),' for the product set'
                  NEWPFMIN(CANDIDATES(NRANDOM))=NEWPFMIN(CANDIDATES(NRANDOM))+1.0D0
               ENDDO
               PFTOTALA=PFTOTALA+NINT(NRWREACTANT*RWPROB(J1))
            ENDDO
            PFTOTALA=LOG(PFTOTALA) ! partition functions are stored as log's
            NCOUNT=0
            DO J1=1,NMINA
               IF (NEWPFMIN(J1).GT.0.0D0) THEN
                  NCOUNT=NCOUNT+1
                  LOCATIONA(NCOUNT)=LOCATIONA(J1)
                  PFMIN(LOCATIONA(NCOUNT))=LOG(NEWPFMIN(J1))
                  PRINT '(A,I8,A,G20.10)','setup> relative weight for reactant minimum ',LOCATIONA(NCOUNT),' is ',
     &                        EXP(PFMIN(LOCATIONA(NCOUNT))-PFTOTALA)
               ENDIF
            ENDDO
            NMINA=NCOUNT
            PRINT '(A,I8,A)','setup> there are now ',NMINA,' minima of type A'
         ENDIF
         DEALLOCATE(NEWPFMIN,CANDIDATES)
      ENDIF
C
C  Load transition states.
C
      DO J1=1,NMIN
         TOPPOINTER(J1)=-1
      ENDDO
      INQUIRE(FILE='ts.data',EXIST=YESNO)
      IF (YESNO.AND.DIJINITSTARTT) THEN
         IF (.NOT.DUMMYTST) THEN
            CALL MYSYSTEM(STATUS,DEBUG,'mv ts.data ts.data.save')
            PRINT '(A)','setup> Removing old ts.data file'
            YESNO=.FALSE.
         ENDIF
      ENDIF
      
      IF (YESNO) THEN
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD')
         J1=0
         DO 
            J1=J1+1
            IF (J1.GT.MAXTS) CALL TSDOUBLE
            IF (IMFRQT) THEN
               READ(UTSDATA,*,END=40) ETS(J1),FVIBTS(J1),HORDERTS(J1),PLUS(J1),MINUS(J1),IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
            ELSE
               READ(UTSDATA,*,END=40) ETS(J1),FVIBTS(J1),HORDERTS(J1),PLUS(J1),MINUS(J1),IXTS(J1),IYTS(J1),IZTS(J1)
            ENDIF
            IF ((PLUS(J1).GT.NMIN).OR.(MINUS(J1).GT.NMIN)) THEN
               PRINT '(A,I10,A)','setup> ERROR - minima specified for ts ',J1,' lie beyond those specified in min.data'
               PRINT '(A,2I10)','setup> plus and minus minima are ',PLUS(J1),MINUS(J1)
               STOP
            ENDIF

            IF (DUMMYTST.AND.(.NOT.NOPOINTS).AND.(NATTEMPT.GT.0)) THEN
               READ(UMIN,REC=PLUS(J1)) (LOCALPOINTS(J2),J2=1,NOPT)
               READ(UMIN,REC=MINUS(J1)) (LOCALPOINTS2(J3),J3=1,NOPT)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                                                     DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               IF (DISTANCE.LT.MINDISTMIN(PLUS(J1))) MINDISTMIN(PLUS(J1))=DISTANCE
               IF (DISTANCE.LT.MINDISTMIN(MINUS(J1))) MINDISTMIN(MINUS(J1))=DISTANCE
            ENDIF
         ENDDO
40       CLOSE(UTSDATA) ! SAT need to reopen this file
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE",FORM="FORMATTED") ! read used in Dijkstra
         NTS=J1-1
         IF (PRINTT) WRITE(*,'(A,I10,A)') 'setup> parameters read for ',NTS,' ts'
         IF (RELATIVEET) PRINT '(A,G20.10)','setup> Shifting TS potential energies relative to zero for the lowest minimum'
C        N.B. GLOBALMIN is zero if RELATIVEET is .false.
         DO J1=1,NTS
            ETS(J1)=ETS(J1)-GLOBALMIN
         ENDDO

         IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN
            INQUIRE(FILE='ts.attempts',EXIST=YESNO)
            TSATTEMPT(1:NTS)=0
            IF (YESNO) THEN
               PRINT '(A)','setup> Reading the number of searches for existing transition states from ts.attempts'
               OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
               J2=0
               DO J1=1,NTS
                  READ(1,'(I8)',END=51) TSATTEMPT(J1)
                  J2=J2+1
               ENDDO
51             CLOSE(1)
               IF (J2.LT.NTS) PRINT '(A)','setup> WARNING - end of file ts.attempts, remaining attempts set to zero'
            ENDIF
         ENDIF

         PRINT '(A)','setup> Setting up ts pointers'
         DO J1=1,NTS
            POINTERP(J1)=-1
            POINTERM(J1)=-1
         ENDDO
         DO J1=NTS,1,-1
            IF (J1.GT.TOPPOINTER(PLUS(J1)))  TOPPOINTER(PLUS(J1))=J1
            IF (J1.GT.TOPPOINTER(MINUS(J1))) TOPPOINTER(MINUS(J1))=J1
            DO J2=J1-1,1,-1
               IF (PLUS(J2).EQ.PLUS(J1)) THEN
                  POINTERP(J1)=J2
                  GOTO 41
               ELSE IF (MINUS(J2).EQ.PLUS(J1)) THEN
                  POINTERP(J1)=J2
                  GOTO 41
               ENDIF
            ENDDO
41          CONTINUE
            DO J2=J1-1,1,-1
               IF (PLUS(J2).EQ.MINUS(J1)) THEN
                  POINTERM(J1)=J2
                  GOTO 42
               ELSE IF (MINUS(J2).EQ.MINUS(J1)) THEN
                  POINTERM(J1)=J2
                  GOTO 42
               ENDIF
            ENDDO
42          CONTINUE
         ENDDO
!        IF (DEBUG) THEN
!           DO J1=1,NMIN
!              WRITE(*,'(A,I10,A,I10)') 'setup> for minimum ',J1,' last occurrence is for ts number ',TOPPOINTER(J1)
!           ENDDO
!           DO J1=1,NTS
!              WRITE(*,'(A,I10,A,4I10)') 'setup> for ts ',J1,' +,-,p+,p-:',PLUS(J1),MINUS(J1),POINTERP(J1),POINTERM(J1)
!           ENDDO
!        ENDIF
    
         IF ((.NOT.NOPOINTS).AND.(NATTEMPT.GT.0).AND.(.NOT.DUMMYTST)) THEN
            PRINT '(A)','setup> Checking ts moments of inertia'
            IF (AMHT) THEN
              WRITE(*,*) 'setup> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
            ELSE
               DO J1=1,NTS
                  READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,NOPT)
                  CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,IXM,IYM,IZM)
!                 IF (DEBUG) WRITE(*,'(A,I10,2F17.7,3I6,3F15.5)') 'setup> ',J1,ETS(J1),FVIBTS(J1),HORDERTS(J1),
!    1                                                            PLUS(J1),MINUS(J1),IXM,IYM,IZM
                  IF ((ABS(IXM-IXTS(J1)).GT.IDIFFTOL).OR.
     1                (ABS(IYM-IYTS(J1)).GT.IDIFFTOL).OR.
     1                (ABS(IZM-IZTS(J1)).GT.IDIFFTOL)) THEN
                     WRITE(*,'(A,I10)') 'setup> WARNING - principal moments of inertia do not agree with input for ts ',J1
                     WRITE(*,'(A)') 'values in ts.data:'
                     WRITE(*,'(3F20.10)') IXTS(J1),IYTS(J1),IZTS(J1)
                     WRITE(*,'(A)') 'values recalculated from points.ts:'
                     WRITE(*,'(3F20.10)') IXM,IYM,IZM
!                    STOP  
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSE
         WRITE(*,'(A)') 'setup> no transition states found'
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file ts.data already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW')
         NTS=0
      ENDIF
C
C  Create a ts entry for DUMMYTS runs if there are minima that seem to lack such entries.
C
      IF (DUMMYTST) THEN
         PRINT '(A)',' setup> shortest distances for local minima:'
         DO J1=1,NMIN
            IF (MINDISTMIN(J1).GT.HUGE(1.0D0)/1.0D1) THEN
               PRINT '(A,I8,A,G20.10)',' setup> in setup, minimum ',J1,' shortest distance unassigned'
               READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,NOPT)
               DO J3=1,NMIN
                  IF (J3.EQ.J1) CYCLE
                  READ(UMIN,REC=J3) (LOCALPOINTS2(J2),J2=1,NOPT)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF ((DISTANCE.LT.MINDISTMIN(J1)).OR.(DISTANCE.LT.MINDISTMIN(J3))) THEN
                     IF (DISTANCE.LT.MINDISTMIN(J1)) MINDISTMIN(J1)=DISTANCE
                     IF (DISTANCE.LT.MINDISTMIN(J3)) MINDISTMIN(J3)=DISTANCE
C
C  Must create an entry in ts.data in this case.
C  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
C
                     IF (IMFRQT) THEN
                        PRINT '()', "setup> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                        STOP
                     ENDIF

                     WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') TSEGUESS(EMIN(J1),EMIN(J3),MINCURVE(J1),MINCURVE(J3),DISTANCE),
     &                           TSFVIBGUESS(EMIN(J1),EMIN(J3),FVIBMIN(J1),FVIBMIN(J3),MINFRQ2(J1),MINFRQ2(J3),NATOMS),
     &                           1,J3,J1,1.0D0,1.0D0,1.0D0
                     CALL FLUSH(UTSDATA)
                     NTS=NTS+1
                     IF (NTS.GT.MAXTS) CALL TSDOUBLE
                     ETS(NTS)=TSEGUESS(EMIN(J1),EMIN(J3),MINCURVE(J1),MINCURVE(J3),DISTANCE)
                     FVIBTS(NTS)=TSFVIBGUESS(EMIN(J1),EMIN(J3),FVIBMIN(J1),FVIBMIN(J3),MINFRQ2(J1),MINFRQ2(J3),NATOMS)
                     HORDERTS(NTS)=1
                     IXTS(NTS)=1.0D0
                     IYTS(NTS)=1.0D0
                     IZTS(NTS)=1.0D0
                     PLUS(NTS)=J3
                     MINUS(NTS)=J1
                     IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
                     IF (DEBUG) WRITE(*,'(A,I10,A)') 'setup> dummy ts ',NTS,' writing parameters to file ts.data'
C
C  Update ts pointers.
C
                     POINTERP(NTS)=-1
                     POINTERM(NTS)=-1
                     IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
                     IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
                     TOPPOINTER(PLUS(NTS))=NTS
                     TOPPOINTER(MINUS(NTS))=NTS
                  ENDIF
               ENDDO
            ENDIF
            PRINT '(A,I8,A,G20.10)',' setup> shortest distance for minimum ',J1,' is ',MINDISTMIN(J1)
         ENDDO
      ENDIF
C
C  Set transition state vibrational product to unity for consistency if NOFRQS is true.
C  Won;t work with FREEPAIRS unless the OPTIM runs are also run with NOFRQS.
C
      IF (NOFRQS) THEN
         FVIBTS(1:NTS)=1.0D0
         NEGEIG(1:NTS)=-1.0D0
      ENDIF
      IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
!
!  Procedure to remove selected stationary points specified by min.remove and ts.remove.
!  First line of each file gives the numbers of structures to remove.
!
      IF (REMOVESP) THEN
         OPEN(UNIT=1,FILE='min.remove',STATUS='OLD')
         READ(1,*) NMINREMOVE
         ALLOCATE(MINPREV(NMIN))
         IF (NMINREMOVE .GT. 0) THEN
            ALLOCATE(MINREMOVE(NMINREMOVE))
            READ(1,*) MINREMOVE(1:NMINREMOVE)
            PRINT '(A)','setup> removing the following minima:'
            PRINT '(10I8)',MINREMOVE(1:NMINREMOVE)
         ENDIF
         CLOSE(1)
         OPEN(UNIT=2,FILE='min.data.removed',STATUS='UNKNOWN')
         OPEN(UNIT=4,FILE='points.min.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT)
         NDUMMY=0
         MINPREV(1:NMIN)=0
         minloop: DO J1=1,NMIN
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.J1) THEN
                  PRINT '(A,I8)','setup> removing minimum ',J1
                  CYCLE minloop
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> not removing minimum ',J1
            NDUMMY=NDUMMY+1
            MINPREV(J1)=NDUMMY
            WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
            READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,NOPT)
            WRITE(4,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
         ENDDO minloop
         CLOSE(2)
         CLOSE(4)
!
! rewrite min.A and min.B in min.A.removed and min.B.removed since numbering may change.
!
         NDUMMY=0
         Aloop: DO J1=1,NMINA
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.LOCATIONA(J1)) THEN
                  PRINT '(A,I8)','setup> removing A minimum ',LOCATIONA(J1)
                  CYCLE Aloop
               ENDIF
            ENDDO
            NDUMMY=NDUMMY+1
         ENDDO Aloop
         OPEN(UNIT=2,FILE='min.A.removed',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all A minima removed'
            STOP
         ENDIF
         DO J1=1,NMINA
            IF (MINPREV(LOCATIONA(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONA(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         NDUMMY=0
         Bloop: DO J1=1,NMINB
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.LOCATIONB(J1)) THEN
                  PRINT '(A,I8)','setup> removing B minimum ',LOCATIONB(J1)
                  CYCLE Bloop
               ENDIF
            ENDDO
            NDUMMY=NDUMMY+1
         ENDDO Bloop
         OPEN(UNIT=2,FILE='min.B.removed',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all B minima removed'
            STOP
         ENDIF
         DO J1=1,NMINB
            IF (MINPREV(LOCATIONB(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONB(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         OPEN(UNIT=1,FILE='ts.remove',STATUS='OLD')
         READ(1,*) NTSREMOVE
         IF (NTSREMOVE.LE.0) GOTO 444
         ALLOCATE(TSREMOVE(NTSREMOVE))
         READ(1,*) TSREMOVE(1:NTSREMOVE)
         CLOSE(1)
         PRINT '(A)','setup> removing the following transition states:'
         PRINT '(10I8)',TSREMOVE(1:NTSREMOVE)
444      CONTINUE
         OPEN(UNIT=3,FILE='ts.data.removed',STATUS='UNKNOWN')
         OPEN(UNIT=5,FILE='points.ts.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT)
         NDUMMY=0
         tsloop: DO J1=1,NTS
            DO J2=1,NTSREMOVE
               IF (TSREMOVE(J2).EQ.J1) CYCLE tsloop
            ENDDO
            IF (MINPREV(PLUS(J1)).EQ.0) THEN
               PRINT '(A,I8,A,I8,A)','setup> possible ERROR - transition state ',J1,' links minimum ',PLUS(J1), 
     &                               ' which has been removed - removing this ts'
               CYCLE tsloop
            ENDIF
            IF (MINPREV(MINUS(J1)).EQ.0) THEN
               PRINT '(A,I8,A,I8,A)','setup> possible ERROR - transition state ',J1,' links minimum ',MINUS(J1), 
     &                               ' which has been removed - removing this ts'
               CYCLE tsloop
            ENDIF
            NDUMMY=NDUMMY+1
            IF (IMFRQT) THEN
               WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
            ELSE
               WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
            ENDIF
            READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,NOPT)
            WRITE(5,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
         ENDDO tsloop
         CLOSE(3); CLOSE(5)
         STOP
      ENDIF
!
!  Procedure to retain selected stationary points specified by min.retain.
!  First line of each file gives the numbers of structures to retain.
!  All ts linking minima in the retain list are themselves retained.
!
      IF (RETAINSP) THEN
         OPEN(UNIT=1,FILE='min.retain',STATUS='OLD')
         READ(1,*) NMINRETAIN
         ALLOCATE(MINRETAIN(NMINRETAIN),MINPREV(NMIN))
         READ(1,*) MINRETAIN(1:NMINRETAIN)
         CLOSE(1)
         PRINT '(A)','setup> retaining the following minima:'
         PRINT '(10I8)',MINRETAIN(1:NMINRETAIN)
         OPEN(UNIT=2,FILE='min.data.retained',STATUS='UNKNOWN')
         OPEN(UNIT=4,FILE='points.min.retained',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT)
         NDUMMY=0
         MINPREV(1:NMIN)=0
         minloop2: DO J1=1,NMIN
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.J1) THEN
                  PRINT '(A,I8)','setup> retaining minimum ',J1
                  NDUMMY=NDUMMY+1
                  MINPREV(J1)=NDUMMY
                  WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
                  IF (.NOT. NOPOINTS) THEN
                     READ(UMIN,REC=J1) (LOCALPOINTS(J3),J3=1,NOPT)
                     WRITE(4,REC=NDUMMY) (LOCALPOINTS(J3),J3=1,NOPT)
                  ENDIF
                  CYCLE minloop2
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> removing minimum ',J1
         ENDDO minloop2
         CLOSE(2)
         CLOSE(4)
!
! rewrite min.A and min.B in min.A.retained and min.B.retained since numbering may change.
!
         NDUMMY=0
         Aloop2: DO J1=1,NMINA
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.LOCATIONA(J1)) THEN
                  NDUMMY=NDUMMY+1
                  PRINT '(A,I8)','setup> retaining A minimum ',LOCATIONA(J1)
                  CYCLE Aloop2
               ENDIF
            ENDDO
         ENDDO Aloop2
         OPEN(UNIT=2,FILE='min.A.retained',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all A minima removed'
            STOP
         ENDIF
         DO J1=1,NMINA
            IF (MINPREV(LOCATIONA(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONA(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         NDUMMY=0
         Bloop2: DO J1=1,NMINB
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.LOCATIONB(J1)) THEN
                  NDUMMY=NDUMMY+1
                  PRINT '(A,I8)','setup> retaining B minimum ',LOCATIONB(J1)
                  CYCLE Bloop2
               ENDIF
            ENDDO
         ENDDO Bloop2
         OPEN(UNIT=2,FILE='min.B.retained',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all B minima removed'
            STOP
         ENDIF
         DO J1=1,NMINB
            IF (MINPREV(LOCATIONB(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONB(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         OPEN(UNIT=3,FILE='ts.data.retained',STATUS='UNKNOWN')
         OPEN(UNIT=5,FILE='points.ts.retained',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*NOPT)
         NDUMMY=0
         tsloop2: DO J1=1,NTS
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.PLUS(J1)) THEN
                  DO J3=1,NMINRETAIN
                     IF (MINRETAIN(J3).EQ.MINUS(J1)) THEN
                        NDUMMY=NDUMMY+1
                        PRINT '(A,I8,A,2I8)','setup> retaining ts ',J1,' connected to retained minima ',PLUS(J1),MINUS(J1)
                        IF (IMFRQT) THEN
                           WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
                        ELSE
                           WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
                        ENDIF
                        IF (.NOT. NOPOINTS) THEN
                           READ(UTS,REC=J1) (LOCALPOINTS(J4),J4=1,NOPT)
                           WRITE(5,REC=NDUMMY) (LOCALPOINTS(J4),J4=1,NOPT)
                        ENDIF
                        CYCLE tsloop2
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> removing ts ',J1
         ENDDO tsloop2
         CLOSE(3); CLOSE(5)
         PRINT '(A,I8,A)','setup> ',NDUMMY,' transition states retained'
         STOP
      ENDIF
C
C  Rate constants.
C
!     KMEAN=0.0D0
      PRINT '(A)','setup> Calculating rate constants'
      IF (ENSEMBLE.EQ.'T') THEN
         DO J1=1,NTS
            KPLUS(J1)  = LOG(1.0D0 * HORDERMIN(PLUS(J1))  / (2.0D0 * PI*HORDERTS(J1))) +
     1             (FVIBMIN(PLUS(J1))  - FVIBTS(J1)) / 2.0D0 - (ETS(J1) - EMIN(PLUS(J1)) )/TEMPERATURE
            IF (FRICTIONT) KPLUS(J1)=KPLUS(J1)+LOG(FRICTIONFAC(NEGEIG(J1)))
            KMINUS(J1) = LOG(1.0D0 * HORDERMIN(MINUS(J1)) / (2.0D0 * PI*HORDERTS(J1))) +
     1             (FVIBMIN(MINUS(J1)) - FVIBTS(J1)) / 2.0D0 - (ETS(J1) - EMIN(MINUS(J1)))/TEMPERATURE
            IF (FRICTIONT) KMINUS(J1)=KMINUS(J1)+LOG(FRICTIONFAC(NEGEIG(J1)))
            IF (ZSYM(1:2).EQ.'CA') KPLUS(J1)=KPLUS(J1)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(J1)=KMINUS(J1)+30.66356D0
            IF (PLUS(J1).EQ.MINUS(J1)) KPLUS(J1)=KPLUS(J1)+LOG(2.0D0)
            IF (PLUS(J1).EQ.MINUS(J1)) KMINUS(J1)=KMINUS(J1)+LOG(2.0D0)
!           KMEAN=KMEAN+KPLUS(J1)+KMINUS(J1)
!           IF (DEBUG) WRITE(*,'(A,3I10,5G20.10,3G20.10)') 'setup> J1,PLUS,MINUS,Ets,E+,E-,k+,k-=',J1,PLUS(J1),MINUS(J1),
!    1                                            ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1)
         ENDDO
      ELSE
         DO J1=1,NTS
            IF (TOTALE.GT.ETS(J1)) THEN
               KPLUS(J1)  = LOG(1.0D0 * HORDERMIN(PLUS(J1))  / (2*PI*HORDERTS(J1))) +
     1                   (FVIBMIN(PLUS(J1))  - FVIBTS(J1))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J1))/(TOTALE-EMIN(PLUS(J1))))
               KMINUS(J1) = LOG(1.0D0 * HORDERMIN(MINUS(J1)) / (2*PI*HORDERTS(J1))) +
     1                   (FVIBMIN(MINUS(J1)) - FVIBTS(J1))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J1))/(TOTALE-EMIN(MINUS(J1))))
               IF (ZSYM(1:2).EQ.'CA') KPLUS(J1)=KPLUS(J1)+30.66356D0
               IF (ZSYM(1:2).EQ.'CA') KMINUS(J1)=KMINUS(J1)+30.66356D0
               IF (PLUS(J1).EQ.MINUS(J1)) KPLUS(J1)=KPLUS(J1)+LOG(2.0D0)
               IF (PLUS(J1).EQ.MINUS(J1)) KMINUS(J1)=KMINUS(J1)+LOG(2.0D0)
!              KMEAN=KMEAN+KPLUS(J1)+KMINUS(J1)
            ELSE
               KPLUS(J1)=-1.0D250
               KMINUS(J1)=-1.0D250
            ENDIF
         ENDDO
      ENDIF
!     IF (NTS.GT.0) KMEAN=KMEAN/(2.0D0*NTS)
!     PRINT '(A,G20.10)', 'setup> Mean log rate constant=', KMEAN
C
C  Sums of rates out of the intermediate minima
C
!       DO J1=1,NMIN
!          KSUM(J1)=0.0D0
!       ENDDO
!       DO J1=1,NTS
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(PLUS(J1))=KSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(MINUS(J1))=KSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
!       ENDDO
!       DO J1=1,NMIN
!          IF (KSUM(J1).GT.0.0D0) THEN
!             KSUM(J1)=LOG(KSUM(J1))+KMEAN
! !           IF (DEBUG) WRITE(*,'(A,I10,2E20.10)') 'setup> J1,KSUM=',J1,KSUM(J1)
!          ENDIF
!       ENDDO
!       DO J1=1,NTS
! !        IF (DEBUG) WRITE(*,'(A,I10,2E20.10)') 'setup> J1,k+,k-=',J1,KPLUS(J1),KMINUS(J1)
!       ENDDO
!
!  Procedure to remove stationary points that are unconnected from A or B sets (or both)
!  according to the prevailing NCONNMIN value.
!
      IF (REMOVEUNCONNECTEDT) THEN
         CALL REMOVE_UNCONNECTED
         STOP
      ENDIF
!
      IF (MERGEDBT) THEN
         CALL MERGEDB
         STOP
      ENDIF

      IF (NPFOLD.GT.0) THEN
         INQUIRE(FILE='commit.data',EXIST=YESNO)
         GPFOLD(1:NMIN)=0.0D0
         IF (YESNO) THEN
            PRINT '(A)','setup> Reading initial committor probabilities read from commit.data'
            OPEN(UNIT=1,FILE='commit.data',STATUS='OLD')
            J2=0
            DO J1=1,NMIN
               READ(1,*,END=110) GPFOLD(J1)
               J2=J2+1
            ENDDO 
110         CLOSE(1)
            IF (J2.LT.NMIN) PRINT '(A)','setup> WARNING - end of file commit.data, remaining probabilities set to zero'
         ELSE
            IF (DIRECTION.EQ.'AB') THEN ! GPFOLD is PFA
               DO J1=1,NMINA
                  GPFOLD(LOCATIONA(J1))=1.0D0
               ENDDO
            ELSE ! GPFOLD is PFB
               DO J1=1,NMINB
                  GPFOLD(LOCATIONB(J1))=1.0D0
               ENDDO
            ENDIF
            PRINT '(A)','setup> Initial committor probabilities set to 0 or 1'
!           PRINT '(6G20.10)',GPFOLD(1:NMIN)
         ENDIF
      ENDIF
C
C  Read in the pairs of minima previously searched in pairs.data exists.
C
      ALLOCATE(PAIR1(MAXPAIRS),PAIR2(MAXPAIRS))
      INQUIRE(FILE='pairs.data',EXIST=YESNO)
      NPAIRDONE=0
      IF (YESNO) THEN
         PRINT '(A)','setup> Reading pairs.data'
         OPEN(UNIT=1,FILE='pairs.data',STATUS='OLD')
         DO 
            NPAIRDONE=NPAIRDONE+1
            IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
            READ(1,*,END=120) PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched pair number ',
     &                                      NPAIRDONE,' is ',PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF ((PAIR1(NPAIRDONE).GT.NMIN).OR.(PAIR2(NPAIRDONE).GT.NMIN)) THEN
               PRINT '(A)','setup> ERROR *** minima specified in pairs.data do not exist in min.data'
               STOP
            ENDIF
         ENDDO
120      CLOSE(1)
         NPAIRDONE=NPAIRDONE-1
         PRINT '(A,I8,A)','setup> ',NPAIRDONE,' pairs of minima already searched read from pairs.data'
      ENDIF
C
C  Read in the minima previously searched in UNTRAP runs.
C
      ALLOCATE(MINDONE(MAXDONE))
      INQUIRE(FILE='min.done',EXIST=YESNO)
      NMINDONE=0
      IF (YESNO) THEN
         PRINT '(A)','setup> Reading min.done'
         OPEN(UNIT=1,FILE='min.done',STATUS='OLD')
         DO 
            NMINDONE=NMINDONE+1
            IF (NMINDONE.GT.MAXDONE) CALL DONEDOUBLE
            READ(1,*,END=121) MINDONE(NMINDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched minimum ',
     &                                      NMINDONE,' is ',MINDONE(NMINDONE)
         ENDDO
121      CLOSE(1)
         NMINDONE=NMINDONE-1
         PRINT '(A,I8,A)','setup> ',NMINDONE,' minima already searched read from min.done'
      ENDIF
C
C Optimised distance calculation for target minimum and any range for the other minima.
C
      IF (DISTANCET) THEN
         READ(UMIN,REC=DISTANCETO) (LOCALPOINTS(J2),J2=1,NOPT)
         DO J1=DISTANCETO1,DISTANCETO2
            READ(UMIN,REC=J1) (LOCALPOINTS2(J2),J2=1,NOPT)

            CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                       DIST2,RIGIDBODY,RMAT,.FALSE.)
            PRINT '(I10,G20.5)',J1,DISTANCE
!           PRINT '(I6,G20.10)',NATOMS
!           PRINT '(A)','reference geometry'
!           DO J2=1,NATOMS
!              PRINT '(A2,1X,3G20.10)','AX',LOCALPOINTS(3*(J2-1)+1:3*(J2-1)+3)
!           ENDDO
!           PRINT '(I6,G20.10)',NATOMS
!           PRINT '(A)','aligned geometry'
!           DO J2=1,NATOMS
!              PRINT '(A2,1X,3G20.10)','AX',LOCALPOINTS2(3*(J2-1)+1:3*(J2-1)+3)
!           ENDDO
         ENDDO
         STOP
      ENDIF
C
C  Initialise PAIRDIST array for use in making an intial connection.
C  PAIRDIST should contain zero if the two minima are linked by a transition state.
C  PAIRLIST contains the index of the other minimum.
C
      IF (DIJINITT) THEN
         DO J1=1,NTS
!
! JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's 
! in place of NCONN() for the plus and minus minima.
!
            CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),HUGE(1),HUGE(1), 
     &                   PLUS(J1),MINUS(J1),.TRUE.,CUT_UNDERFLOW,DEADTS)
         ENDDO
         IF (INITIALDIST) THEN
            ALLPAIRS(1:(NMIN*(NMIN-1))/2)=-DISBOUND
            PRINT '(A,I12)','nmin*(nmin-1)/2=',(NMIN*(NMIN-1))/2
            INQUIRE(FILE='allpairs',EXIST=YESNO)
            IF (YESNO) THEN
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE='allpairs',STATUS='OLD')
               DO J1=1,(NMIN*(NMIN-1))/2
                  READ(LUNIT,*) ALLPAIRS(J1)
               ENDDO
               CLOSE(LUNIT)
               PRINT '(A,I8)','setup> Pair distance values read'
            ELSE
               CALL GETMETRIC(1,NMIN)
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE='allpairs',STATUS='UNKNOWN')
               WRITE(LUNIT,'(G20.10)') ALLPAIRS(1:(NMIN*(NMIN-1))/2)
               CLOSE(LUNIT)
            ENDIF
         ELSE
            PAIRDIST(1:NMIN,1:PAIRDISTMAX)=1.0D100
            PAIRLIST(1:NMIN,1:PAIRDISTMAX)=-1
            INQUIRE(FILE='pairdist',EXIST=YESNO)
            IF (PAIRDIST1.NE.0) YESNO=.FALSE. ! so we can write new entries for READMIN etc.
            IF (YESNO) THEN
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE='pairdist',STATUS='OLD')
               DO J1=1,NMIN
                  READ(LUNIT,*) (PAIRDIST(J1,J2),J2=1,PAIRDISTMAX)
               ENDDO
               CLOSE(LUNIT)
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE='pairlist',STATUS='OLD')
               DO J1=1,NMIN
                  READ(LUNIT,*) (PAIRLIST(J1,J2),J2=1,PAIRDISTMAX)
               ENDDO
               CLOSE(LUNIT)
               PRINT '(A,I8)','setup> Pair distance metric values read'
            ELSE
               IF (PAIRDIST1.EQ.0) PAIRDIST1=1
               IF (PAIRDIST2.EQ.0) PAIRDIST2=NMIN
               IF (PAIRDIST1.GT.NMIN) STOP
               PAIRDIST2=MIN(PAIRDIST2,NMIN)
               CALL GETMETRIC(PAIRDIST1,PAIRDIST2)
               IF ((PAIRDIST1.EQ.1).AND.(PAIRDIST2.EQ.NMIN)) THEN
                  OPEN(UNIT=1,FILE='pairdist',STATUS='UNKNOWN')
                  DO J3=1,NMIN
                     WRITE(1,'(10G20.10)') (PAIRDIST(J3,J4),J4=1,PAIRDISTMAX)
                  ENDDO
                  CLOSE(1)
                  OPEN(UNIT=1,FILE='pairlist',STATUS='UNKNOWN')
                  DO J3=1,NMIN
                     WRITE(1,'(10I10)') (PAIRLIST(J3,J4),J4=1,PAIRDISTMAX)
                  ENDDO
                  CLOSE(1)
               ELSE
                  WRITE(S1,'(I10)') PAIRDIST1
                  WRITE(S2,'(I10)') PAIRDIST2
                  WRITE(FNAME,'(A)') 'pairdist.' // TRIM(ADJUSTL(S1)) // '.' // TRIM(ADJUSTL(S2))
                  OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='UNKNOWN')
                  DO J3=PAIRDIST1,PAIRDIST2
                     WRITE(1,'(10G20.10)') (PAIRDIST(J3,J4),J4=1,PAIRDISTMAX)
                  ENDDO
                  CLOSE(1)
                  WRITE(FNAME,'(A)') 'pairlist.' // TRIM(ADJUSTL(S1)) // '.' // TRIM(ADJUSTL(S2))
                  OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='UNKNOWN')
                  DO J3=PAIRDIST1,PAIRDIST2
                     WRITE(1,'(10I10)') (PAIRLIST(J3,J4),J4=1,PAIRDISTMAX)
                  ENDDO
                  CLOSE(1)
                  STOP
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
C  Add transition states and minima from the <PATHNAME> file.
C  Use GETNEWPATH to do the bookkeeping.
C
      IF (ADDPATH) THEN
         CALL MYSYSTEM(STATUS,DEBUG,'cp ' // TRIM(ADJUSTL(PATHNAME)) // ' path.info')
!        IF (ADDTRIPLES) THEN
            CALL GETALLPATHS
!        ELSE
!           CALL GETNEWPATH(0,0)
!        ENDIF
      ENDIF
C
C IF NOPOINTGROUPT is true then set all point group orders to unity. This is designed for
C addressable potentials, where permutational isomers are enumerated explicitly. We can 
C still lump enantiomers, since there will always be two of these for each stationary point,
C except for planar gemetries.
C
      IF (NOPOINTGROUPT) THEN
         DO J3=1,NMIN
            HORDERMIN(J3)=1
         ENDDO
         DO J3=1,NTS
            HORDERTS(J3)=1
         ENDDO
      ENDIF
C
C If CONNECTPAIRST is true then read the pairs of minima to connect from file CONNECTPAIRSFILE
C
      IF (CONNECTPAIRST) THEN
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(CONNECTPAIRSFILE)),STATUS='OLD')
         NCONNECTPAIRS=0
         DO
            READ(LUNIT,*,END=119) NDUMMY, NDUMMY
            NCONNECTPAIRS=NCONNECTPAIRS+1
         ENDDO
119      REWIND(LUNIT)
         PRINT '(A,A,A,I8)','setup> Number of minima pairs in file ',TRIM(ADJUSTL(CONNECTPAIRSFILE)),' is ',NCONNECTPAIRS
         ALLOCATE(CONNECTPAIRSMIN(NCONNECTPAIRS,2))
         PRINT '(A)','setup> Connection pairs of local minima:'
         DO J1=1,NCONNECTPAIRS
            READ(LUNIT,*) CONNECTPAIRSMIN(J1,1),CONNECTPAIRSMIN(J1,2)
            WRITE(*,'(2I10)') CONNECTPAIRSMIN(J1,1),CONNECTPAIRSMIN(J1,2)
         ENDDO
         CLOSE(LUNIT)
      ENDIF
C
C If USEPAIRST is true then read the sequence of minima from file USEPAIRSFILE
C USEPAIRSFILE must be formatted as a single Epath file
C
      IF (USEPAIRST) THEN
         OPEN(UNIT=1,FILE=TRIM(ADJUSTL(USEPAIRSFILE)),STATUS='OLD')
         NUSEPAIRS=0
         DO
            READ(1,*,END=111) NDUMMY, DUMMY, NDUMMY
            NUSEPAIRS=NUSEPAIRS+1
         ENDDO
111      REWIND(1)
         PRINT '(A,A,A,I8,A,I8,A)','setup> Number of lines in file ',TRIM(ADJUSTL(USEPAIRSFILE)),' is ',NUSEPAIRS,' i.e. ',
     &           (NUSEPAIRS+1)/2,' minima'
         NUSEPAIRS=(NUSEPAIRS+1)/2
         ALLOCATE(USEPAIRSMIN(NUSEPAIRS))
         DO J1=1,NUSEPAIRS
            READ(1,*) NDUMMY, DUMMY, USEPAIRSMIN(J1)
            IF (J1.EQ.NUSEPAIRS) EXIT
            READ(1,*) NDUMMY, DUMMY, NDUMMY
         ENDDO
         CLOSE(1)
         PRINT '(A)','setup> Sequence of local minima:'
         PRINT '(15I8)',USEPAIRSMIN(1:NUSEPAIRS)
      ENDIF

      IF (DOST) CALL DOS
      IF (CVMINIMAT) THEN
         CALL CVMINIMA
         STOP
      ELSEIF (CVT) THEN
         CALL CV
         STOP
      ENDIF
      IF (SHANNONT) CALL SHANNON
      IF (MICROTHERMT) CALL MICROTHERM

      IF ((CONNECTIONS.GT.1).AND.(CHECKCONNECTIONST)) THEN
         WRITE(*,'(A,I6,A)') 'setup> checking for at least ',CONNECTIONS,' connections per minimum'
         WRITE(*,'(A,I6,A)') 'setup> WARNING *** use the NEWCONNECTIONS keyword to use more than one core'
         DO J1=1,NMIN
            CALL TSSEARCH(J1,0)
         ENDDO
      ENDIF

!     INQUIRE(UNIT=UMIN,OPENED=TESTOP,NAMED=TESTNAME)
!     PRINT *,'UNIT UMIN is ',UMIN
!     PRINT *,'opened is ',TESTOP
!     PRINT *,'named is ',TESTNAME

!     INQUIRE(UNIT=UTS,OPENED=TESTOP,NAMED=TESTNAME)
!     PRINT *,'UNIT UTS is ',UTS
!     PRINT *,'opened is ',TESTOP
!     PRINT *,'named is ',TESTNAME
      

      RETURN 
      END

      DOUBLE PRECISION FUNCTION TSEGUESS(E1,E2,C1,C2,DISTANCE)
      IMPLICIT NONE
      DOUBLE PRECISION E1, E2, DISTANCE, C1, C2, ARGUMENT
C
C     ARGUMENT=c1*c2*distance**2 - 6*(c1 - c2)*(e1 - e2)
C     IF (ARGUMENT.LT.0.0D0) ARGUMENT=0.0D0
C     IF (C1.EQ.C2) THEN
C        TSEGUESS=((c1*distance**2 + 6*e1 - 6*e2)**2/(4.*c1*distance**2) + 6*e2)/6.
C     ELSE
C        TSEGUESS=(c1*c2*(c1 + c2)*distance**2 - 2*c1*c2*distance*Sqrt(ARGUMENT)
C    &             + 6*(c1 - c2)*(-(c2*e1) + c1*e2))/(6.*(c1 - c2)**2)
C     ENDIF
C     IF (TSEGUESS.LT.MAX(E1,E2)) TSEGUESS=MAX(E1,E2)+ABS(E1-E2)
C
C  Double quadratic formulation.
C      
C     ARGUMENT=c1*c2*distance**2 - 2*(c1 - c2)*(e1 - e2)
C     IF (ARGUMENT.LT.0.0D0) ARGUMENT=0.0D0
C     IF (C1.EQ.C2) THEN
C        TSEGUESS=((c1*distance**2)/4. + e1 + (e1 - e2)**2/(c1*distance**2) + e2)/2.
C     ELSE
C        TSEGUESS=(c1*c2*(c1 + c2)*distance**2 - 2*c1*c2*distance*
C    &     Sqrt(ARGUMENT) + 2*(c1 - c2)*(-(c2*e1) + c1*e2))/(2.*(c1 - c2)**2)
C     ENDIF

      TSEGUESS=MAX(E1,E2)+DISTANCE
      
      END FUNCTION TSEGUESS

      DOUBLE PRECISION FUNCTION TSFVIBGUESS(E1,E2,FVIB1,FVIB2,MINF1,MINF2,NATOMS)
      USE COMMONS,ONLY : NOPT
      IMPLICIT NONE
      DOUBLE PRECISION E1, E2, FVIB1,FVIB2, MINF1, MINF2
      INTEGER NATOMS
C
C  The conversion factor for CHARMM and AMBER is included in the MINFRQ2 values read from min.data.info
C  The MINFRQ2 values are read as the ln from min.data.info
C
C     IF (E1.GT.E2) THEN
C        TSFVIBGUESS=FVIB1-MINF1
C     ELSE
C        TSFVIBGUESS=FVIB2-MINF2
C     ENDIF
      IF (E1.GT.E2) THEN
         TSFVIBGUESS=FVIB1*(NOPT-7)/(NOPT-6)
      ELSE
         TSFVIBGUESS=FVIB2*(NOPT-7)/(NOPT-6)
      ENDIF

      
      END FUNCTION TSFVIBGUESS

      FUNCTION DUPLICATES_CHECK(ILIST, N)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(IN) :: ILIST(N)
         LOGICAL DUPLICATES_CHECK
         INTEGER SORTED_LIST(N), TEMPLIST(N), I, VPREV
         DUPLICATES_CHECK = .FALSE.
         SORTED_LIST(:) = ILIST(:)

         ! sn402: added following block because of segfaults when trying to build a new database using READMIN
         IF (N==0) THEN
            DUPLICATES_CHECK = .FALSE.
            RETURN
         ENDIF

         DO I=1,N
            TEMPLIST(I) = I
         ENDDO
         !CALL SORT(N, N, TEMPLIST, SORTED_LIST)
         CALL SORT2(N, N, SORTED_LIST, TEMPLIST)

         VPREV = SORTED_LIST(1)
         DO I=2,N
            IF (VPREV .EQ. SORTED_LIST(I) ) THEN
               DUPLICATES_CHECK = .TRUE.
               RETURN
            ENDIF
            !WRITE(*,*) SORTED_LIST(I), VPREV, TEMPLIST(I), ILIST(I)
            VPREV = SORTED_LIST(I)
         ENDDO
      END FUNCTION DUPLICATES_CHECK
