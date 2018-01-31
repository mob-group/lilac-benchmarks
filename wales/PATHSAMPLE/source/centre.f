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
C  Move the centre of mass to the origin.
C
      SUBROUTINE CENTRE(T1)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION CMX, CMY, CMZ, T1(3*NATOMS)
      INTEGER I
C
C  Move centre of mass to the origin.
C
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+T1(3*(I-1)+1)
         CMY=CMY+T1(3*(I-1)+2)
         CMZ=CMZ+T1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
      DO I=1,NATOMS
         T1(3*(I-1)+1)=T1(3*(I-1)+1)-CMX
         T1(3*(I-1)+2)=T1(3*(I-1)+2)-CMY
         T1(3*(I-1)+3)=T1(3*(I-1)+3)-CMZ
      ENDDO

      RETURN
      END
