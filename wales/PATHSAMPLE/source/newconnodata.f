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
C  Make an odata tspath file for the minimum in LOCALPOINTS1
C  Note that the last line of odata.tspath should be POINTS, like odata.connect.
C  For historical reasons the POINTS line should not be present for the old style
C  odata.tsseach file!
C
      SUBROUTINE NEWCONNODATA(CONNID,LOCALPOINTS1)
      USE COMMONS
      USE UTILS,ONLY : GETUNIT
      IMPLICIT NONE
      INTEGER J2, CONNID, STATUS, LUNIT
      DOUBLE PRECISION LOCALPOINTS1(3*NATOMS), SAVEPOINTS(3*NATOMS), RANDOM, DPRAND
      CHARACTER(LEN=20) UNSTRING
      CHARACTER(LEN=10) CONNSTR
      CHARACTER(LEN=80) FPOO

      WRITE(CONNSTR,'(I10)') CONNID
      CALL MYSYSTEM(STATUS,DEBUG,'cat odata.tspath > odata.'//TRIM(ADJUSTL(CONNSTR)))

!
! perturbation is done in cycle2
!
!     IF (MLP3T.OR.MLPB3T) THEN
!        DO J2=1,NOPT
!           IF (FROZEN(J2)) THEN
!              LOCALPOINTS1(J2)=LOCALPOINTS1(J2)
!           ELSE
!              RANDOM=DPRAND()
!              LOCALPOINTS1(J2)=LOCALPOINTS1(J2)+(RANDOM-0.5D0)*2.0D0*PERTVALUE
!           ENDIF
!        ENDDO
!     ELSE
!        DO J2=1,3*NATOMS
!           RANDOM=DPRAND()
!           LOCALPOINTS1(J2)=LOCALPOINTS1(J2)+(RANDOM-0.5D0)*2.0D0*PERTVALUE
!        ENDDO
!     ENDIF

              PRINT '(A)','LOCALPOINTS1 in newconnodata:'
              PRINT '(A,2I10)','NATOMS,NOPT=',NATOMS,NOPT
               PRINT '(6G20.10)',LOCALPOINTS1(1:NOPT)

      IF (CHARMMT) THEN
         IF (MACHINE) THEN ! SAT
            DO J2=1,3*NATOMS
               SAVEPOINTS(J2)=LOCALPOINTS1(J2)
            ENDDO
            CALL CHARMMDUMP(SAVEPOINTS,'points1.inp.'//TRIM(ADJUSTL(CONNSTR)))
         ELSE
            DO J2=1,3*NATOMS
               SAVEPOINTS(J2)=LOCALPOINTS1(J2)
            ENDDO
            CALL CHARMMDUMP(SAVEPOINTS,'input.crd.'//TRIM(ADJUSTL(CONNSTR)))
         ENDIF
      ELSE IF (UNRST) THEN
         DO J2=1,3*NATOMS
            SAVEPOINTS(J2)=LOCALPOINTS1(J2)
         ENDDO
         WRITE(UNSTRING,'(A)') 'coords.'//TRIM(ADJUSTL(CONNSTR))
         CALL MYUNRESDUMP(SAVEPOINTS,UNSTRING)
      ELSE IF (AMBERT.OR.OPEPT) THEN
         FPOO='start.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(LUNIT,'(3F20.10)') (LOCALPOINTS1(3*(J2-1)+1),LOCALPOINTS1(3*(J2-1)+2),
     &                               LOCALPOINTS1(3*(J2-1)+3),J2=1,NATOMS)
         CLOSE(LUNIT)
      ELSE IF (AMHT) THEN
         FPOO='start.'//TRIM(ADJUSTL(CONNSTR)) ! 
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
         WRITE(LUNIT,'(3G25.15)') LOCALPOINTS1(1:3*NATOMS)
         CLOSE(LUNIT)
      ELSE IF (MLP3T.OR.MLPB3T) THEN
         FPOO='odata.'//TRIM(ADJUSTL(CONNSTR)) 
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(LUNIT,'(A2,2X,F20.10)') ('  ',LOCALPOINTS1(J2),J2=1,NOPT)
         CLOSE(LUNIT)
      ELSE IF (PHI4MODT) THEN
         FPOO='odata.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(LUNIT,'(A2,2X,F20.10)') ('  ',LOCALPOINTS1(J2),J2=1,NATOMS)
         CLOSE(LUNIT)
      ELSE
         FPOO='odata.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(LUNIT,'(A2,2X,3F20.10)') (ZSYMBOL(J2),LOCALPOINTS1(3*(J2-1)+1),LOCALPOINTS1(3*(J2-1)+2),
     &                               LOCALPOINTS1(3*(J2-1)+3),J2=1,NATOMS)
         CLOSE(LUNIT)
      ENDIF

      RETURN
      END

