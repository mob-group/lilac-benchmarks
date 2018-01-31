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
!
!  Cycle over adjacent minima of MIN1 and calculate their new waiting times
!  and branching probabilities.
!
SUBROUTINE NGTRENORM_CRSTORAGE(MIN1,NC1,NMIN,EMKSUM,PPROD,NCOL,DEBUG,NCONNMAX,GBMAX,PB1,NV1)
USE NGTMEM
IMPLICIT NONE
INTEGER J1, J2, NC1, MIN2, J3, MIN1, NMIN, NEWNCOL, LASTONE1, LASTONE2, NCONNMAX, NC2, LNV1, LNV2
INTEGER NCOL(NMIN), NVTEMP(NMIN), NV1(NMIN), NV2(NMIN)
DOUBLE PRECISION EMKSUM(NMIN), PPROD, NORM, GBMAX, PB1LOCAL(NMIN), PBTEMP(NMIN), PB1(NMIN), PB2(NMIN), DUMMY
LOGICAL DEBUG
!
!  MIN1 must be the last entry in its own neighbour list for NGT!
!  Remove it from consideration.
!
NC1=NC1-1
!
! Case where MIN1 only has one entry is common; treat it separately. 
! Here we are renormalising the entries for MIN2.
!
IF (NC1.EQ.1) THEN
   MIN2=NV1(1)
   J3=NCOL(MIN2)
   IF (J3.EQ.0) THEN
      PRINT '(A,I6,A,I6)','ngtrenorm_crstorage> WARNING *** no columns for minimum ',MIN2,' the only neighbour of ',MIN1
      RETURN
   ENDIF
   DUMMY=DVEC(J3+ROW_PTR(MIN2)-1)*PPROD
   DO J1=1,J3-1 ! find the entry for MIN2 to MIN2 and renormalise
      IF (COL_IND(J1+ROW_PTR(MIN2)-1).EQ.MIN2) THEN
         DVEC(J1+ROW_PTR(MIN2)-1)=DVEC(J1+ROW_PTR(MIN2)-1)+PB1(1)*DUMMY
         EMKSUM(MIN2)=EMKSUM(MIN2)+EMKSUM(MIN1)*DUMMY
         EXIT
      ENDIF
   ENDDO
   NCOL(MIN2)=J3-1
   RETURN
ENDIF

DO J2=1,NC1
   MIN2=NV1(J2)
!
!  MIN1 must be in the neighbour list of MIN2 unless MIN2 is a sink. 
!  MIN1 must be entry NCOL(MIN2) if we maintain sorted lists!
!
   J3=NCOL(MIN2)
   IF (J3.EQ.0) CYCLE
   IF (COL_IND(J3+ROW_PTR(MIN2)-1).NE.MIN1) THEN
      PRINT *,'NGTrenorm_crstorage> ERROR, MIN1 is not the last minimum in the list for MIN2'
      STOP
   ENDIF
   PB1LOCAL(1:NC1)=PB1(1:NC1)*DVEC(J3+ROW_PTR(MIN2)-1)*PPROD
!
!  Waiting time for min2 (beta).
!
   EMKSUM(MIN2)=EMKSUM(MIN2)+EMKSUM(MIN1)*DVEC(J3+ROW_PTR(MIN2)-1)*PPROD
!
!
!  New branching probabilities for MIN2 to all its old neighbours and any neighbours
!  of MIN1 that it was not previously adjacent to.
!  For NGT we include a renormalised branching probability to the same minimum.
!  Note that the connection lists for MIN1 and MIN2 are sorted!
!
   NEWNCOL=0
   LASTONE2=1
   LASTONE1=1
!
!  Remove the MIN1 connection from MIN2 entries
!  It is simply the last entry for sorted lists!
!
   NCOL(MIN2)=NCOL(MIN2)-1
   NC2=NCOL(MIN2)
   NV2(1:NC2)=COL_IND(1+ROW_PTR(MIN2)-1:NC2+ROW_PTR(MIN2)-1)
   PB2(1:NC2)=DVEC(1+ROW_PTR(MIN2)-1:NC2+ROW_PTR(MIN2)-1)
!
!  Merge the two sorted lists of remaining connections.
!
   DO
      IF ((LASTONE1.GT.NC1).OR.(LASTONE2.GT.NC2)) EXIT 
      NEWNCOL=NEWNCOL+1
      LNV2=NV2(LASTONE2)
      LNV1=NV1(LASTONE1)
      IF (LNV2.LT.LNV1) THEN
         NVTEMP(NEWNCOL)=LNV2
         PBTEMP(NEWNCOL)=PB2(LASTONE2)
         LASTONE2=LASTONE2+1
      ELSEIF (LNV1.LT.LNV2) THEN
         NVTEMP(NEWNCOL)=LNV1
         PBTEMP(NEWNCOL)=PB1LOCAL(LASTONE1) 
         LASTONE1=LASTONE1+1
      ELSE
         PBTEMP(NEWNCOL)=PB2(LASTONE2)+PB1LOCAL(LASTONE1)
         NVTEMP(NEWNCOL)=LNV1
         LASTONE1=LASTONE1+1
         LASTONE2=LASTONE2+1 
      ENDIF
   ENDDO
   DO 
      IF (LASTONE1.GT.NC1) EXIT
      NEWNCOL=NEWNCOL+1
      NVTEMP(NEWNCOL)=NV1(LASTONE1)
      PBTEMP(NEWNCOL)=PB1LOCAL(LASTONE1) 
      LASTONE1=LASTONE1+1
   ENDDO
   DO 
      IF (LASTONE2.GT.NC2) EXIT
      NEWNCOL=NEWNCOL+1
      NVTEMP(NEWNCOL)=NV2(LASTONE2)
      PBTEMP(NEWNCOL)=PB2(LASTONE2)
      LASTONE2=LASTONE2+1
   ENDDO
! 
!  Allocate more space if required. We don't need PBRANCH(x,y) and NVAL(x,y) beyond y=MIN1
!  so long as regrouping has been performed to put the A minima first, followed by the B minima,
!  followed by the I minima.
!  This reordering is done in the REGROUP subroutine.
!
! JMC need to call NGTREALLOC_CRSTORAGE every time to move the elements within DVEC and COL_IND, 
! even if we don't need to reallocate memory.
   CALL NGTREALLOC_CRSTORAGE(MIN1,NEWNCOL,MIN2,GBMAX,NMIN,NCOL)
   IF (DEBUG) PRINT '(2(A,I8))','NGTrenorm_crstorage> connectivity for minimum ',MIN2,' is now ',NEWNCOL

   NCOL(MIN2)=NEWNCOL
   COL_IND(1+ROW_PTR(MIN2)-1:NEWNCOL+ROW_PTR(MIN2)-1)=NVTEMP(1:NEWNCOL)
   DVEC(1+ROW_PTR(MIN2)-1:NEWNCOL+ROW_PTR(MIN2)-1)=PBTEMP(1:NEWNCOL)
!
!  Checking normalisation here adds 16% to run time for a test case.
!
!  NORM=0.0D0
!  DO J3=1,NCOL(MIN2)
!     NORM=NORM+DVEC(J3+ROW_PTR(MIN2)-1)
!  ENDDO
!  IF (ABS(NORM-1.0D0).GT.1.0D-6) THEN
!     PRINT '(A,I8,A,G20.10)','NGTrenorm_crstorage> WARNING - sum of branching probabilities out of minimum ',MIN2, &
!    &                        ' deviates from unity: ',NORM
!     IF (ABS(NORM-1.0D0).GT.1.0D-1) THEN
!        NORM=0.0D0
!        DO J3=1,NCOL(MIN2)
!           NORM=NORM+DVEC(J3+ROW_PTR(MIN2)-1)
!           PRINT '(A,2I8,2G20.10)','NGTrenorm_crstorage> J3,NVAL,PBRANCH,SUM=',J3,COL_IND(J3+ROW_PTR(MIN2)-1),
!    &                                                    DVEC(J3+ROW_PTR(MIN2)-1),NORM
!        ENDDO
!        PRINT '(A)','NGTrenorm_crstorage> Deviation is beyond tolerance'
!        STOP
!     ENDIF
!  ENDIF
ENDDO

END SUBROUTINE NGTRENORM_CRSTORAGE
