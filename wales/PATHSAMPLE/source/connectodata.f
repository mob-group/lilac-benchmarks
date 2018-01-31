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
C  Make an odata connect file for the minima in LOCALPOINTS1 and LOCALPOINTS2.
C
      SUBROUTINE CONNECTODATA(CONNID,LOCALPOINTS1,LOCALPOINTS2)
      USE COMMONS
      IMPLICIT NONE
      INTEGER J2, CONNID, STATUS
      DOUBLE PRECISION DISTANCE, LOCALPOINTS1(NOPT), LOCALPOINTS2(NOPT), SAVEPOINTS(NOPT), 
     &                 RMAT(3,3), DIST2
      DOUBLE PRECISION DUMMY
      CHARACTER(LEN=20) UNSTRING
      CHARACTER(LEN=10) CONNSTR
      CHARACTER(LEN=80) FPOO
      CHARACTER(LEN=80) BHSTRING1,BHSTRING2

      WRITE(CONNSTR,'(I10)') CONNID
!
! DISTANCE is only used with BHINTERPT - don't really need to calculate it here otherwise?
!
      IF (.NOT.(MLP3T.OR.MLPB3T)) THEN
         IF ((ABS(LOCALPOINTS2(1)-LOCALPOINTS2(4)).LT.1.0D-6).AND. 
     &       (ABS(LOCALPOINTS2(2)-LOCALPOINTS2(5)).LT.1.0D-6).AND. 
     &       (ABS(LOCALPOINTS2(3)-LOCALPOINTS2(6)).LT.1.0D-6)) THEN
            PRINT '(A)','connectodata> here A ERROR *** atoms 1 and 2 are in the same place'
            STOP
         ENDIF
      ENDIF
      CALL MINPERMDIST(LOCALPOINTS1,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
      IF (.NOT.(MLP3T.OR.MLPB3T)) THEN
         IF ((ABS(LOCALPOINTS2(1)-LOCALPOINTS2(4)).LT.1.0D-6).AND. 
     &       (ABS(LOCALPOINTS2(2)-LOCALPOINTS2(5)).LT.1.0D-6).AND. 
     &       (ABS(LOCALPOINTS2(3)-LOCALPOINTS2(6)).LT.1.0D-6)) THEN
            PRINT '(A)','connectodata> here B ERROR *** atoms 1 and 2 are in the same place'
            STOP
         ENDIF
      ENDIF
C
C  Other DOALLMIN might need an alternative to this STOPDIST for BLJ.
C      
      IF (ZSYM.EQ.'LS') THEN
!        FPOO='odata.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
!        OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
!        IF (DIRECTION.EQ.'BA') WRITE(1,'(A,F20.10)') 'STOPDISP   ',ORDERPARAM
!        IF (DIRECTION.EQ.'AB') WRITE(1,'(A,F20.10)') 'STOPDISP   ',-ORDERPARAM
!        CLOSE(1)
         IF (BHINTERPT) THEN
            DUMMY=MAX(DISTANCE/2.0D0,BHDISTTHRESH) ! aim to halve the separation with each call
            IF (ICINTERPT) THEN
               WRITE(BHSTRING1,'(A9,F12.3,1X,F12.3,1X,I6,1X,G12.4,A)') 
     &                  'BHINTERP ',DUMMY,BHMAXENERGY,BHSTEPS,BHCONV,' +++'
               WRITE(BHSTRING2,'(3X,3(F12.4,1X),G12.4,1X,G12.4,A)')
     &                  BHTEMP,BHSTEPSIZE,BHACCREJ,BHK,BHSFRAC,' ICINTERP'
            ELSE
               WRITE(BHSTRING1,'(A9,F12.3,1X,F12.3,1X,I6,1X,G12.4,A)') 
     &                  'BHINTERP ',DUMMY,BHMAXENERGY,BHSTEPS,BHCONV,' +++'
               WRITE(BHSTRING2,'(3X,3(F12.4,1X),G12.4,1X,G12.4)')
     &                  BHTEMP,BHSTEPSIZE,BHACCREJ,BHK,BHSFRAC
            ENDIF
            CALL MYSYSTEM(STATUS,DEBUG,'echo "'//BHSTRING1//'" > odata.'//TRIM(ADJUSTL(CONNSTR)))
            CALL MYSYSTEM(STATUS,DEBUG,'echo "'//BHSTRING2//'" >> odata.'//TRIM(ADJUSTL(CONNSTR)))
            CALL MYSYSTEM(STATUS,DEBUG,'cat odata.bhinterp >> odata.' // TRIM(ADJUSTL(CONNSTR)))
         ELSE IF (BISECTT) THEN
            IF (ICINTERPT) THEN
               WRITE(BHSTRING1,'(A7,F12.3,1X,F12.3,1X,I6,1X,I6,A)') 
     &                  'BISECT ',BISECTMINDIST,BISECTMAXENERGY,BISECTSTEPS,BISECTMAXATTEMPTS,' ICINTERP'
            ELSE
               WRITE(BHSTRING1,'(A7,F12.3,1X,F12.3,1X,I6,1X,I6)')
     &                  'BISECT ',BISECTMINDIST,BISECTMAXENERGY,BISECTSTEPS,BISECTMAXATTEMPTS
            ENDIF
            CALL MYSYSTEM(STATUS,DEBUG,'echo "'//BHSTRING1//'" > odata.'//TRIM(ADJUSTL(CONNSTR)))
            CALL MYSYSTEM(STATUS,DEBUG,'cat odata.bisect >> odata.' // TRIM(ADJUSTL(CONNSTR)))
         ELSE
            CALL MYSYSTEM(STATUS,DEBUG,'cat odata.connect > odata.' // TRIM(ADJUSTL(CONNSTR)))
         ENDIF
      ELSE
         IF (BHINTERPT) THEN
            DUMMY=MAX(DISTANCE/2.0D0,BHDISTTHRESH) ! aim to halve the separation with each call
            IF (ICINTERPT) THEN
               WRITE(BHSTRING1,'(A9,F12.3,1X,F12.3,1X,I6,1X,G12.4,A)') 
     &                  'BHINTERP ',DUMMY,BHMAXENERGY,BHSTEPS,BHCONV,' +++'
               WRITE(BHSTRING2,'(3X,3(F12.4,1X),G12.4,1X,G12.4,A)')
     &                  BHTEMP,BHSTEPSIZE,BHACCREJ,BHK,BHSFRAC,' ICINTERP'
            ELSE
               WRITE(BHSTRING1,'(A9,F12.3,1X,F12.3,1X,I6,1X,G12.4,A)') 
     &                  'BHINTERP ',DUMMY,BHMAXENERGY,BHSTEPS,BHCONV,' +++'
               WRITE(BHSTRING2,'(3X,3(F12.4,1X),G12.4,1X,G12.4)')
     &                  BHTEMP,BHSTEPSIZE,BHACCREJ,BHK,BHSFRAC
            ENDIF
            CALL MYSYSTEM(STATUS,DEBUG,'echo "'//BHSTRING1//'" > odata.'//TRIM(ADJUSTL(CONNSTR)))
            CALL MYSYSTEM(STATUS,DEBUG,'echo "'//BHSTRING2//'" >> odata.'//TRIM(ADJUSTL(CONNSTR)))
            CALL MYSYSTEM(STATUS,DEBUG,'cat odata.bhinterp >> odata.' // TRIM(ADJUSTL(CONNSTR)))
         ELSE IF (BISECTT) THEN
            IF (ICINTERPT) THEN
               WRITE(BHSTRING1,'(A7,F12.3,1X,F12.3,1X,I6,1X,I6,A)')
     &                  'BISECT ',BISECTMINDIST,BISECTMAXENERGY,BISECTSTEPS,BISECTMAXATTEMPTS,' ICINTERP'
            ELSE
               WRITE(BHSTRING1,'(A7,F12.3,1X,F12.3,1X,I6,1X,I6)')
     &                  'BISECT ',BISECTMINDIST,BISECTMAXENERGY,BISECTSTEPS,BISECTMAXATTEMPTS
            ENDIF
            CALL MYSYSTEM(STATUS,DEBUG,'echo "'//BHSTRING1//'" > odata.'//TRIM(ADJUSTL(CONNSTR)))
            CALL MYSYSTEM(STATUS,DEBUG,'cat odata.bisect >> odata.' // TRIM(ADJUSTL(CONNSTR)))
         ELSE
            CALL MYSYSTEM(STATUS,DEBUG,'cat odata.connect > odata.'//TRIM(ADJUSTL(CONNSTR)))
         ENDIF
      ENDIF

      IF (CHARMMT) THEN
         if (machine) then ! SAT
              DO J2=1,3*NATOMS
                 SAVEPOINTS(J2)=LOCALPOINTS1(J2)
              ENDDO
              CALL CHARMMDUMP(SAVEPOINTS,'points1.inp.'//TRIM(ADJUSTL(CONNSTR)))
              DO J2=1,3*NATOMS
                 SAVEPOINTS(J2)=LOCALPOINTS2(J2)
              ENDDO
              CALL CHARMMDUMP(SAVEPOINTS,'points2.inp.'//TRIM(ADJUSTL(CONNSTR)))
         else
              DO J2=1,3*NATOMS
                 SAVEPOINTS(J2)=LOCALPOINTS1(J2)
              ENDDO
              CALL CHARMMDUMP(SAVEPOINTS,'input.crd.'//TRIM(ADJUSTL(CONNSTR)))
              DO J2=1,3*NATOMS
                 SAVEPOINTS(J2)=LOCALPOINTS2(J2)
              ENDDO
              CALL CHARMMDUMP(SAVEPOINTS,'finish.'//TRIM(ADJUSTL(CONNSTR)))
         endif
      ELSE IF (UNRST) THEN
         DO J2=1,3*NATOMS
            SAVEPOINTS(J2)=LOCALPOINTS1(J2)
         ENDDO
         WRITE(UNSTRING,'(A)') 'coords.'//TRIM(ADJUSTL(CONNSTR))
         CALL MYUNRESDUMP(SAVEPOINTS,UNSTRING)
         DO J2=1,3*NATOMS
            SAVEPOINTS(J2)=LOCALPOINTS2(J2)
         ENDDO
         CALL MYUNRESDUMP(SAVEPOINTS,'finish.'//TRIM(ADJUSTL(CONNSTR)))
      ELSE IF (AMBERT.OR.OPEPT) THEN
         FPOO='start.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
         OPEN(2,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(2,'(3F20.10)') (LOCALPOINTS1(3*(J2-1)+1),LOCALPOINTS1(3*(J2-1)+2),
     1                               LOCALPOINTS1(3*(J2-1)+3),J2=1,NATOMS)
         CLOSE(2)
         FPOO='finish.'//TRIM(ADJUSTL(CONNSTR))
         OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(1,'(3F20.10)') (LOCALPOINTS2(3*(J2-1)+1),LOCALPOINTS2(3*(J2-1)+2),LOCALPOINTS2(3*(J2-1)+3),J2=1,NATOMS)
         CLOSE(1)
      ELSE IF (AMHT) THEN

         FPOO='start.'//TRIM(ADJUSTL(CONNSTR)) ! 
         OPEN(2,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(2,'(3G25.15)') LOCALPOINTS1(1:3*NATOMS)
         CLOSE(2)

         FPOO='finish.'//TRIM(ADJUSTL(CONNSTR))
         OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(1,'(3G25.15)') LOCALPOINTS2(1:3*NATOMS)
         CLOSE(1)
      ELSE IF (MLP3T.OR.MLPB3T) THEN
         FPOO='odata.'//TRIM(ADJUSTL(CONNSTR)) 
         OPEN(2,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(2,'(A2,2X,F20.10)') (ZSYMBOL(J2),LOCALPOINTS1(J2),J2=1,NOPT)
         CLOSE(2)
         FPOO='finish.'//TRIM(ADJUSTL(CONNSTR))
         OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(1,'(F20.10)') (LOCALPOINTS2(J2),J2=1,NOPT)
         CLOSE(1)

      ELSE
         FPOO='odata.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
         OPEN(2,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(2,'(A2,2X,3F20.10)') (ZSYMBOL(J2),LOCALPOINTS1(3*(J2-1)+1),LOCALPOINTS1(3*(J2-1)+2),
     1                               LOCALPOINTS1(3*(J2-1)+3),J2=1,NATOMS)
         CLOSE(2)
         FPOO='finish.'//TRIM(ADJUSTL(CONNSTR))
         OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(1,'(3F20.10)') (LOCALPOINTS2(3*(J2-1)+1),LOCALPOINTS2(3*(J2-1)+2),LOCALPOINTS2(3*(J2-1)+3),J2=1,NATOMS)
         CLOSE(1)
         IF (.NOT.(MLP3T.OR.MLPB3T)) THEN
            IF ((ABS(LOCALPOINTS2(1)-LOCALPOINTS2(4)).LT.1.0D-6).AND. 
     &          (ABS(LOCALPOINTS2(2)-LOCALPOINTS2(5)).LT.1.0D-6).AND. 
     &          (ABS(LOCALPOINTS2(3)-LOCALPOINTS2(6)).LT.1.0D-6)) THEN
               PRINT '(A)','connectodata> here Z ERROR *** atoms 1 and 2 are in the same place'
               STOP
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END

