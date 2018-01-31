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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Reads a database of minima and saves it in individual files
!
      SUBROUTINE CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
      USE COMMONS,ONLY : RIGIDBODY,TWOD,BULKT,BOXLZ,BOXLY,BOXLX
      IMPLICIT NONE
      INTEGER J1, J2, NMIN,NATOMS,CNT,K,NTS,UMIN,UTS
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS),NATIVEPOINTS(3*NATOMS),HEL,DIST,HELRMSD,RMSD,RMAT(3,3),DISTANCE,DIST2
      CHARACTER*10 STRING
      LOGICAL DEBUG

      WRITE(*,*) 'opening the p53 database file' 
      OPEN(UNIT=10,FILE='points.min',ACCESS='DIRECT', &
     & FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      OPEN(UNIT=91,FILE='UNFOLDED.DBASE',STATUS='UNKNOWN')
      OPEN(UNIT=92,FILE='FOLDED.DBASE',STATUS='UNKNOWN')
      OPEN(UNIT=100,FILE='native.crd',STATUS='UNKNOWN')
      READ(100,*) (NATIVEPOINTS(J2),J2=1,3*NATOMS)
      DO J1=1,NMIN
        RMSD=0.0D0
        HEL=0.0D0
        HELRMSD=0.0D0
        CNT=0	
        READ(10,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
        CALL MINPERMDIST(NATIVEPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
        DO J2=63,187
          HELRMSD=HELRMSD+(LOCALPOINTS(3*(J2-1)+1)- &
     &	  NATIVEPOINTS(3*(J2-1)+1))**2  &
     &    +(LOCALPOINTS(3*(J2-1)+2)-NATIVEPOINTS(3*(J2-1)+2))**2 &
     &    +(LOCALPOINTS(3*(J2-1)+3)-NATIVEPOINTS(3*(J2-1)+3))**2
        ENDDO
        DO J2=1,187
          RMSD=RMSD+(LOCALPOINTS(3*(J2-1)+1)- &
     &    NATIVEPOINTS(3*(J2-1)+1))**2 &
     &    +(LOCALPOINTS(3*(J2-1)+2)-NATIVEPOINTS(3*(J2-1)+2))**2 &
     &    +(LOCALPOINTS(3*(J2-1)+3)-NATIVEPOINTS(3*(J2-1)+3))**2
        ENDDO

        RMSD=DSQRT(RMSD/187.0D0)
        HELRMSD=DSQRT(HELRMSD/125.0D0)

        WRITE (STRING, '(I10)') J1

        DIST=SQRT((LOCALPOINTS(3*(55-1)+1) &
     &	-LOCALPOINTS(3*(85-1)+1))**2+  &
     &  (LOCALPOINTS(3*(55-1)+2)-LOCALPOINTS(3*(85-1)+2))**2+  &
     &  (LOCALPOINTS(3*(55-1)+3)-LOCALPOINTS(3*(85-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        DIST=SQRT((LOCALPOINTS(3*(62-1)+1)- &
     &	LOCALPOINTS(3*(94-1)+1))**2+  &
     &  (LOCALPOINTS(3*(62-1)+2)-LOCALPOINTS(3*(94-1)+2))**2+  &
     &  (LOCALPOINTS(3*(62-1)+3)-LOCALPOINTS(3*(94-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        DIST=SQRT((LOCALPOINTS(3*(67-1)+1)- &
     &	LOCALPOINTS(3*(111-1)+1))**2+  &
     &  (LOCALPOINTS(3*(67-1)+2)-LOCALPOINTS(3*(111-1)+2))**2+  &
     &  (LOCALPOINTS(3*(67-1)+3)-LOCALPOINTS(3*(111-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        DIST=SQRT((LOCALPOINTS(3*(84-1)+1)- &
     &	LOCALPOINTS(3*(128-1)+1))**2+  &
     &  (LOCALPOINTS(3*(84-1)+2)-LOCALPOINTS(3*(128-1)+2))**2+  &
     &  (LOCALPOINTS(3*(84-1)+3)-LOCALPOINTS(3*(128-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        DIST=SQRT((LOCALPOINTS(3*(93-1)+1) &
     &	-LOCALPOINTS(3*(137-1)+1))**2+  &
     &  (LOCALPOINTS(3*(93-1)+2)-LOCALPOINTS(3*(137-1)+2))**2+  &
     &  (LOCALPOINTS(3*(93-1)+3)-LOCALPOINTS(3*(137-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        DIST=SQRT((LOCALPOINTS(3*(110-1)+1) &
     &	-LOCALPOINTS(3*(147-1)+1))**2+  &
     &  (LOCALPOINTS(3*(110-1)+2)-LOCALPOINTS(3*(147-1)+2))**2+  &
     &  (LOCALPOINTS(3*(110-1)+3)-LOCALPOINTS(3*(147-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        DIST=SQRT((LOCALPOINTS(3*(127-1)+1) &
     &	-LOCALPOINTS(3*(157-1)+1))**2+  &
     &  (LOCALPOINTS(3*(127-1)+2)-LOCALPOINTS(3*(157-1)+2))**2+  &
     &  (LOCALPOINTS(3*(127-1)+3)-LOCALPOINTS(3*(157-1)+3))**2)
        IF (DIST.LT.3.6) CNT=CNT+1
        HEL=CNT/7.0D0
        WRITE(*,'(a10,3f15.5)') STRING,HEL,HELRMSD,RMSD

        IF (HEL.LE.0.2D0.AND.HELRMSD.GE.7.5D0) THEN
           WRITE(91,'(a10,4f15.5)')STRING,HEL,HELRMSD,RMSD,DISTANCE
!       ELSEIF ((HEL.GE.0.85D0).AND.(HELRMSD.LT.3.0D0).AND.(RMSD.LT.3.5D0)) THEN 
        ELSEIF ((HEL.GE.0.57D0).AND.(HELRMSD.LT.4.00D0)) THEN 
           WRITE(92,'(a10,4f15.5)') STRING,HEL,HELRMSD, RMSD,DISTANCE
        ENDIF

      ENDDO
      CLOSE(10)
      CLOSE(100)
      CLOSE(91)
      CLOSE(92)
      RETURN 
      END
