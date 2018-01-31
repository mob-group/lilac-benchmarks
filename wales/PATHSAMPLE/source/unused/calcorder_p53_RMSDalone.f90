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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Reads a database of minima for p53 and calculates order parameters based on helix RMSD
!
SUBROUTINE CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
USE COMMONS,ONLY : RIGIDBODY,TWOD,BULKT,BOXLZ,BOXLY,BOXLX
IMPLICIT NONE
INTEGER J1, J2, NMIN,NATOMS,CNT,K,NTS,UMIN,UTS
DOUBLE PRECISION LOCALPOINTS(3*NATOMS),NATIVEPOINTS(3*NATOMS),RMSD,RMAT(3,3),DISTANCE,DIST2
LOGICAL DEBUG

WRITE(*,*) 'calcorder> Opening the database file for p53'
OPEN(UNIT=UMIN,FILE='points.min',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
OPEN(UNIT=91,FILE='FOLDED.DBASE',STATUS='UNKNOWN')
OPEN(UNIT=92,FILE='UNFOLDED.DBASE',STATUS='UNKNOWN')
OPEN(UNIT=93,FILE='NEITHER.DBASE',STATUS='UNKNOWN')

!The native structure
OPEN(UNIT=100,FILE='FINISH.CRD_POINTS',STATUS='UNKNOWN')
READ(100,*) (NATIVEPOINTS(J2),J2=1,3*NATOMS)
close(100)
DO J1=1,NMIN
   RMSD=0.0D0
   CNT=0	
   READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
   CALL MINPERMDIST(NATIVEPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
   DO J2=63,187       !This spans the helical segment of the fragment only
     RMSD=RMSD+(LOCALPOINTS(3*(J2-1)+1)-NATIVEPOINTS(3*(J2-1)+1))**2 &
&    +(LOCALPOINTS(3*(J2-1)+2)-NATIVEPOINTS(3*(J2-1)+2))**2 &
&    +(LOCALPOINTS(3*(J2-1)+3)-NATIVEPOINTS(3*(J2-1)+3))**2
     CNT=CNT+1 
   ENDDO
   RMSD=DSQRT(RMSD/CNT)
   IF (RMSD.LE.3.5D0) THEN               
      WRITE(91,*) J1
   ELSEIF (RMSD.GT.4.5D0) THEN
      WRITE(92,*) J1
   ELSE 
       WRITE(93,*) J1
   ENDIF
ENDDO
CLOSE(100)
CLOSE(91)
CLOSE(92)
CLOSE(93)

RETURN 
END
