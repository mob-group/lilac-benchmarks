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
!  Calculate order parameter chi for protein L according to
!  PNAS, 102, 16685, 2005
!
SUBROUTINE CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
IMPLICIT NONE
LOGICAL DEBUG
INTEGER NATOMS, J1, J2, J3, NDUMMY, NMIN, NTS, UMIN, UTS
DOUBLE PRECISION LOCALPOINTS(3*NATOMS), THISORDER, DIST, CMX, CMY, CMZ, RADG
DOUBLE PRECISION, PARAMETER :: NORM=1378.0D0
DOUBLE PRECISION DISTNATIVE(1378)

OPEN(UNIT=1,FILE='native',STATUS='UNKNOWN')
READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
CLOSE(1)
NDUMMY=0
DO J1=1,52
   DO J2=J1+4,56
      DIST=SQRT( (LOCALPOINTS(3*(J1-1)+1)-LOCALPOINTS(3*(J2-1)+1))**2 + &
  &              (LOCALPOINTS(3*(J1-1)+2)-LOCALPOINTS(3*(J2-1)+2))**2 + &
  &              (LOCALPOINTS(3*(J1-1)+3)-LOCALPOINTS(3*(J2-1)+3))**2 )
      NDUMMY=NDUMMY+1
      DISTNATIVE(NDUMMY)=DIST
   ENDDO
ENDDO

OPEN(UNIT=1,FILE='min.order.native',STATUS='UNKNOWN')
OPEN(UNIT=2,FILE='min.order.unfolded',STATUS='UNKNOWN')
DO J3=1,NMIN
   READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
!
!  First get the centre of mass to calculate radius of gyration.
!
   CMX=0.0D0; CMY=0.0D0; CMZ=0.0D0
   DO J1=1,56
      CMX=CMX+LOCALPOINTS(3*(J1-1)+1)
      CMY=CMY+LOCALPOINTS(3*(J1-1)+2)
      CMZ=CMZ+LOCALPOINTS(3*(J1-1)+3)
   ENDDO
   CMX=CMX/NATOMS
   CMY=CMY/NATOMS
   CMZ=CMZ/NATOMS
   DIST=0.0D0
   DO J1=1,56
      DIST=DIST+(LOCALPOINTS(3*(J1-1)+1)-CMX)**2 + (LOCALPOINTS(3*(J1-1)+2)-CMY)**2 + (LOCALPOINTS(3*(J1-1)+3)-CMZ)**2
   ENDDO
   RADG=SQRT(DIST/56.0D0)
   NDUMMY=0
   THISORDER=0.0D0
   DO J1=1,52
      DO J2=J1+4,56
         DIST=SQRT( (LOCALPOINTS(3*(J1-1)+1)-LOCALPOINTS(3*(J2-1)+1))**2 + &
  &                 (LOCALPOINTS(3*(J1-1)+2)-LOCALPOINTS(3*(J2-1)+2))**2 + &
  &                 (LOCALPOINTS(3*(J1-1)+3)-LOCALPOINTS(3*(J2-1)+3))**2 )
         NDUMMY=NDUMMY+1
         IF (0.2D0-ABS(DIST-DISTNATIVE(NDUMMY)).GT.0.0D0) THISORDER=THISORDER+1.0D0
      ENDDO
   ENDDO
   THISORDER=THISORDER/NORM
   IF (THISORDER.GT.0.4D0) THEN
      WRITE(1,'(I8,2G20.10)') J3, THISORDER, RADG
   ELSEIF ((THISORDER.LT.0.4D0).AND.(RADG.GT.2.5D0).AND.(RADG.LT.4.5D0)) THEN
      WRITE(2,'(I8,2G20.10)') J3, THISORDER, RADG
   ENDIF
   IF (DEBUG) WRITE(*,'(I8,2G20.10)') J3, THISORDER, RADG
ENDDO
CLOSE(1)
CLOSE(2)
PRINT '(A)','calcorder> Order parameters written to min.order.native and min.order.unfolded'

END SUBROUTINE CALCORDER
