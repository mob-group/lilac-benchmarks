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
!  Remove I minima from the bottom up and renormalise the waiting times and
!  branching probabilities out of adjacent minima.
!  When removing minimum i the waiting times change for all minima, j \in J, adjacent
!  to i (and including i for NGT). The branching probabilities out of these minima also change.
!  For each j \in J new connections are created to minima in J that were not
!  originally connected to j. We finally end up with only A and B minima.
!  Must change the information in the following arrays:
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = id of minimum involved in connection J1 from minimum M2
!  PBRANCH(J1,M2) = branching probability of taking connection J1 from minimum M2 to minimum NVAL(J1,M2)
!  We have to DEALLOCATE and ALLOCATE NVAL and PBRANCH whenever necessary.
!  Waiting time is stored in EMKSUM. 
!  NCONNMAX is the maximum number of connections.
!
SUBROUTINE NGTREMOVEI(NMIN,NMINA,NMINB,NCONNMAX,NCONNMIN,DEBUG,NCOL,GBMAX,EMKSUM,REALLOCMEM)
USE COMMONS,ONLY : NCONN, NGTSWITCH, NGTSIZE, RATESCYCLET, SHANNONT
USE PORFUNCS
USE NGTMEM
IMPLICIT NONE
INTEGER NMIN, NMINA, NMINB, NCONNMAX, NCONNMIN, NC1, J2, MIN1, MINCONN, J1, ISTAT
INTEGER NCOL(NMIN), NV1(NMIN)
DOUBLE PRECISION PII, PPROD, GBMAX, EMKSUM(NMIN), PB1(NMIN)
LOGICAL DEBUG, REALLOCMEM

DO MIN1=NMIN,NMINA+NMINB+1,-1
   IF (NCONN(MIN1).LE.NCONNMIN) CYCLE ! minima with <= nconnmin connections are ignored
   MINCONN=0
   DO J1=1,MIN1
      MINCONN=MINCONN+NCOL(J1)
   ENDDO

   IF ((MINCONN*1.0D0)/(MIN1*1.0D0)**2.GE.NGTSWITCH) THEN
      IF (MIN1.LE.NGTSIZE) THEN 
         IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,F12.3,A)', &
         'NGTremovei> about to switch to dense scheme and allocate storage of ',MIN1*MIN1*8.0D-6,' Mb'
         IF (.NOT.(RATESCYCLET.OR.SHANNONT)) PRINT '(A,F12.6,A,I6)', &
         'NGTremovei> density=',(MINCONN*1.0D0)/(MIN1*MIN1*1.0D0),' remaining minima=',MIN1
         CALL FLUSH(6)
         CALL NGTREMOVEID(NMIN,NMINA,NMINB,NCONNMIN,DEBUG,NCOL,EMKSUM,MIN1,REALLOCMEM) 
         RETURN
      ENDIF
   ENDIF

   NC1=NCOL(MIN1)
   PB1(1:NC1)=PBRANCH(1:NC1,MIN1)
   NV1(1:NC1)=NVAL(1:NC1,MIN1)
   IF (NVAL(NC1,MIN1).NE.MIN1) THEN
      PRINT *,'NGTremovei> ERROR, MIN1 is not the last minimum in the list for MIN1'
      PRINT *,'MIN1=',MIN1
      PRINT *,'NCOL(MIN1)=',NCOL(MIN1)
      PRINT *,'NVAL(NCOL(MIN1),MIN1)=',NVAL(NCOL(MIN1),MIN1)
      PRINT *,'NVAL entries for MIN1:'
      PRINT *,NVAL(1:NCOL(MIN1),MIN1)
      STOP
   ENDIF
   PII=PBRANCH(NC1,MIN1)
!
! Calculate 1-Pii directly by summing over the other branching probabilities if Pii gets
! close to one. Avoids numerical precision loss for calculating 1-Pii.
!
   IF (PII.GT.0.99D0) THEN
      PPROD=0.0D0
      DO J2=1,NC1-1
         PPROD=PPROD+PBRANCH(J2,MIN1)
      ENDDO
   ELSE
      PPROD=1.0D0-PII
   ENDIF
   IF (PPROD.LE.0.0D0) THEN
      PRINT '(A,G20.10,A,I8)','ERROR in NGT, 1-branching probability product is',PPROD,' MIN1=',MIN1
      STOP
   ELSE
      PPROD=1.0D0/PPROD
   ENDIF
!
!  Cycle over adjacent minima of MIN1 and calculate their new waiting times
!  and branching probabilities.
!
   CALL NGTRENORM(MIN1,NC1,NMIN,EMKSUM,PPROD,NCOL,DEBUG,NCONNMAX,GBMAX,PB1,NV1)
   IF (DEBUG) PRINT '(A,I8,A)','NGTremovei> minimum ',MIN1,' has been disconnected'
ENDDO

END SUBROUTINE NGTREMOVEI
