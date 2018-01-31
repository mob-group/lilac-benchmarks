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

      SUBROUTINE MERGE(NCOLMIN2,NEWNCOL,NVAL,PBRANCH,NVTEMP,PBTEMP,NVTMP,
     &                 MERGEDNV,MERGEDPB,PBTMP,NCONNMAX,MIN1,MIN2,NCOL)
      IMPLICIT NONE
      INTEGER NCOLMIN2,NC2,NEWNCOL,MIN1,NCONNMAX,J3,J4,J5,MIN2
      INTEGER NVAL(NCONNMAX,MIN1),NVTEMP(NCONNMAX),NVTMP(NCONNMAX),MERGEDNV(NCONNMAX),NCOL(MIN1)
      DOUBLE PRECISION PBRANCH(NCONNMAX,MIN1),PBTEMP(NCONNMAX),PBTMP(NCONNMAX),MERGEDPB(NCONNMAX)

      NC2=NCOL(MIN2)
      IF (NC2.EQ.0) THEN
         NCOL(MIN2)=NEWNCOL
         NVAL(1:NCOL(MIN2),MIN2)=NVTEMP(1:NCOL(MIN2))
         PBRANCH(1:NCOL(MIN2),MIN2)=PBTEMP(1:NCOL(MIN2))
      ELSE
         NVTMP(1:NC2)=NVAL(1:NC2,MIN2)
         PBTMP(1:NC2)=PBRANCH(1:NC2,MIN2)
         J3=1
         J4=1
         DO J5=1,NC2+NEWNCOL
            IF (NVTMP(J3).LT.NVTEMP(J4)) THEN 
               MERGEDNV(J5)=NVTMP(J3)
               MERGEDPB(J5)=PBTMP(J3)
               IF (J3.EQ.NC2) THEN
                  NVAL(1:J5,MIN2)=MERGEDNV(1:J5)
                  PBRANCH(1:J5,MIN2)=MERGEDPB(1:J5)
                  NVAL(J5+1:NC2+NEWNCOL,MIN2)=NVTEMP(J4:NEWNCOL)
                  PBRANCH(J5+1:NC2+NEWNCOL,MIN2)=PBTEMP(J4:NEWNCOL)
                  EXIT
               ENDIF
               J3=J3+1
            ELSE
               MERGEDNV(J5)=NVTEMP(J4)
               MERGEDPB(J5)=PBTEMP(J4)
                IF (J4.EQ.NEWNCOL) THEN
                  NVAL(1:J5,MIN2)=MERGEDNV(1:J5)
                  PBRANCH(1:J5,MIN2)=MERGEDPB(1:J5)
                  NVAL(J5+1:NC2+NEWNCOL,MIN2)=NVTMP(J3:NC2)
                  PBRANCH(J5+1:NC2+NEWNCOL,MIN2)=PBTMP(J3:NC2)
                  EXIT
                ENDIF
               J4=J4+1
            ENDIF
         ENDDO
         NCOL(MIN2)=NC2+NEWNCOL
      ENDIF
    
      RETURN
      END
