C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE HESSOUT(IR,IC,IDEC,NSYM)
      USE MODHESS
      IMPLICIT NONE
      INTEGER IR, IC, IDEC, NSYM, JX, ITOP, I, J
C
      JX=1
 1    ITOP=MIN(JX+4,IR)
      DO 10 I=1,IR
         IF (NSYM .EQ. 1) WRITE (*,43) (HESS(I,J),J=JX,MIN(IC,ITOP))
         IF (NSYM .EQ. 0) WRITE (*,43) (HESS(I,J),J=JX,MIN(IC,ITOP))
 10   CONTINUE
      JX=JX+5
      IF (JX .GT. IR) RETURN
      WRITE(*,*)
      GOTO 1
 43   FORMAT(5X,T12,F12.6,T25,F12.6,T38,F12.6,T51,F12.6,T64,F12.6)
      END
