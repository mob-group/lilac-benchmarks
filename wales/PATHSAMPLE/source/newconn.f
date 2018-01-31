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

      SUBROUTINE NEWCONN(NC1,NV1,NCONNMAX,MIN1,MIN2,NCOL,NVAL,PBRANCH,PBRANCHMIN1MIN2,PB1,
     &                   NVTEMP,PBTEMP,PPROD,NEWNCOL)
      IMPLICIT NONE
      INTEGER LASTONE,J3,NC1,MIN3,MIN2,NEWNCOL,MIN1,NCONNMAX,J4
      INTEGER NV1(NCONNMAX),NCOL(MIN1),NVAL(NCONNMAX,MIN1),NVTEMP(NCONNMAX)
      DOUBLE PRECISION PBRANCH(NCONNMAX,MIN1),PBRANCHMIN1MIN2,PB1(NCONNMAX),PPROD,PBTEMP(NCONNMAX)

      NEWNCOL=0
!
!  Find new connections for MIN2
!
      LASTONE=1
      DO J3=1,NC1 ! potential new connections
         MIN3=NV1(J3)
         IF (MIN3.EQ.MIN2) CYCLE
666      CONTINUE
C        PRINT *,'J3,MIN3,LASTONE,NVAL(LASTONE,MIN2)=',J3,MIN3,LASTONE,NVAL(LASTONE,MIN2)
         IF (MIN3.LT.NVAL(LASTONE,MIN2)) THEN
            NEWNCOL=NEWNCOL+1   ! New connection from MIN2 to MIN3
            NVTEMP(NEWNCOL)=MIN3
            PBTEMP(NEWNCOL)=PBRANCHMIN1MIN2*PB1(J3)*PPROD
         ELSEIF (MIN3.EQ.NVAL(LASTONE,MIN2)) THEN ! existing connection
            PBRANCH(LASTONE,MIN2)=PBRANCH(LASTONE,MIN2)+PBRANCHMIN1MIN2*PB1(J3)*PPROD
888         LASTONE=LASTONE+1
            IF (LASTONE-1.EQ.NCOL(MIN2)) THEN ! all remaining J3+1:NC1 are new connections
               NVTEMP(NEWNCOL+1:NEWNCOL+NC1-J3)=NV1(J3+1:NC1)
               PBTEMP(NEWNCOL+1:NEWNCOL+NC1-J3)=PBRANCHMIN1MIN2*PB1(J3+1:NC1)*PPROD
               NEWNCOL=NEWNCOL+NC1-J3
               RETURN
            ENDIF
            IF (NVAL(LASTONE,MIN2).GE.MIN3) CYCLE
            GOTO 888
         ELSE
777         LASTONE=LASTONE+1
            IF (LASTONE-1.EQ.NCOL(MIN2)) THEN ! all remaining J3+1:NC1 are new connections
               NVTEMP(NEWNCOL+1:NEWNCOL+NC1-J3)=NV1(J3+1:NC1)
               PBTEMP(NEWNCOL+1:NEWNCOL+NC1-J3)=PBRANCHMIN1MIN2*PB1(J3+1:NC1)*PPROD
               NEWNCOL=NEWNCOL+NC1-J3
               RETURN
            ENDIF
            IF (NVAL(LASTONE,MIN2).GE.MIN3) GOTO 666
            GOTO 777
         ENDIF
      ENDDO 
      
      RETURN
      END
