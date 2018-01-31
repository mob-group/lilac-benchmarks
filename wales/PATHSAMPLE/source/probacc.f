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

      SUBROUTINE PROBACC(NC1,NCONNMAX,NV1,MIN1,MIN2,PB1,PBRANCH,PPROD,NVAL,DEBUG,NCOL)
      IMPLICIT NONE
      LOGICAL DEBUG
      INTEGER NC1,NCONNMAX,MIN2,J4,MIN1
      INTEGER NV1(NCONNMAX),NVAL(NCONNMAX,MIN1),NCOL(MIN1)
      DOUBLE PRECISION P1,PB1(NCONNMAX),P2,PBRANCH(NCONNMAX,MIN1),PPROD

      P1=0.0D0
      DO J4=1,NC1
         IF (NV1(J4).EQ.MIN2) CYCLE
         P1=P1+PB1(J4)
      ENDDO
      P2=0.0D0
      DO J4=1,NCOL(MIN2)
         IF (NVAL(J4,MIN2).EQ.MIN1) CYCLE
         P2=P2+PBRANCH(J4,MIN2)
      ENDDO
      IF (DEBUG) PRINT '(A,2I8,A,2G20.10)','MIN1,MIN2=',MIN1,MIN2,
     &    ' original 1-PPROD, alternative: ',1.0D0-PPROD,P1+P2-P1*P2
      PPROD=P1+P2-P1*P2

      RETURN
      END
