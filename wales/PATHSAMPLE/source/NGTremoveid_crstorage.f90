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
!  This is the dense version, where all the branching probabilities are stored
!  in a square matrix PMAT.
!
SUBROUTINE NGTREMOVEID_CRSTORAGE(NMIN,NMINA,NMINB,NCONNMIN,DEBUG,NCOL,EMKSUM,MIN1,REALLOCMEM)
USE COMMONS,ONLY : NCONN
USE NGTMEM
IMPLICIT NONE
INTEGER NMIN, NMINA, NMINB, NCONNMIN, J2, MIN1, J1, LMIN1, LM1
INTEGER NCOL(NMIN)
DOUBLE PRECISION PII, PPROD, EMKSUM(NMIN), PMAT(MIN1,MIN1), PB(MIN1), PD(MIN1), DUMMY
LOGICAL DEBUG, REALLOCMEM

PMAT(1:MIN1,1:MIN1)=0.0D0
DO J1=1,MIN1
   DO J2=1,NCOL(J1)
      PMAT(COL_IND(J2+ROW_PTR(J1)-1),J1)=DVEC(J2+ROW_PTR(J1)-1)
   ENDDO
ENDDO
IF (REALLOCMEM) DEALLOCATE(DVEC,COL_IND)

DO LMIN1=MIN1,NMINA+NMINB+1,-1
   IF (NCONN(LMIN1).LE.NCONNMIN) CYCLE ! minima with <= nconnmin connections are ignored
   LM1=LMIN1-1
   PD(1:LM1)=PMAT(1:LM1,LMIN1)
   PII=PMAT(LMIN1,LMIN1)
!
! Calculate 1-Pii directly by summing over the other branching probabilities if Pii gets
! close to one. Avoids numerical precision loss for calculating 1-Pii.
!
   IF (PII.GT.0.99D0) THEN
      PPROD=0.0D0
      DO J2=1,LM1
         PPROD=PPROD+PD(J2)
      ENDDO
   ELSE
      PPROD=1.0D0-PII
   ENDIF
   IF (PPROD.LE.0.0D0) THEN
      PRINT '(A,G20.10,A,I8)','NGTremoveid_crstorage> ERROR in NGT, 1-branching probability product is',PPROD,' LMIN1=',LMIN1
      PRINT *,'PD=',PD(1:LMIN1)
      STOP
   ELSE
      PPROD=1.0D0/PPROD
   ENDIF
   PB(1:LM1)=PMAT(LMIN1,1:LM1)*PPROD
   DUMMY=EMKSUM(LMIN1)
   EMKSUM(1:LM1)=EMKSUM(1:LM1)+PB(1:LM1)*DUMMY
!
!  BLAS level 2 routine dger would do this operation.
!  The test for zero PB(J1) was inspired by dger.
!
!  DO J1=1,LM1
!     PMAT(1:LM1,J1)=PMAT(1:LM1,J1)+PD(1:LM1)*PB(J1)
!  ENDDO

   DO J1=1,LM1
      DUMMY=PB(J1)
!     IF (DUMMY.NE.0.0D0) PMAT(1:LM1,J1)=PMAT(1:LM1,J1)+PD(1:LM1)*DUMMY
      IF (DUMMY.EQ.0.0D0) CYCLE
      DO J2=1,LM1
         PMAT(J2,J1)=PMAT(J2,J1)+PD(J2)*DUMMY
      ENDDO
   ENDDO

   IF (DEBUG) PRINT '(A,I8,A)','NGTremoveid_crstorage> minimum ',LMIN1,' has been disconnected'
ENDDO
!
! Now we should only have A then B minima left.
! Revert to sparse storage.
!
IF (REALLOCMEM) ALLOCATE(DVEC((NMINA+NMINB)*(NMINA+NMINB)),COL_IND((NMINA+NMINB)*(NMINA+NMINB)))
ROW_PTR(1)=1
DO J1=1,NMINA+NMINB
   NCOL(J1)=0
   IF (J1.GT.1) ROW_PTR(J1)=ROW_PTR(J1-1)+NMINA+NMINB 
   DO J2=1,NMINA+NMINB
      IF ((PMAT(J2,J1).GT.0.0D0).OR.(J2.EQ.J1)) THEN
         NCOL(J1)=NCOL(J1)+1
         COL_IND(NCOL(J1)+ROW_PTR(J1)-1)=J2
         DVEC(NCOL(J1)+ROW_PTR(J1)-1)=PMAT(J2,J1)
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE NGTREMOVEID_CRSTORAGE
