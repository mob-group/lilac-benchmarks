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

SUBROUTINE SWAPNODE(FIRST,SECOND,PBRANCH,NVAL,NCOL,NCONNMAX,NMIN,EMKSUM)
IMPLICIT NONE
INTEGER FIRST, SECOND, NCONNMAX, NMIN, J1, J2
INTEGER NVAL(NCONNMAX, NMIN), NCOL(NMIN), NDUMMY(NCONNMAX), NCDUMMY
DOUBLE PRECISION PBRANCH(NCONNMAX, NMIN), PDUMMY(NCONNMAX), EMKSUM(NMIN), EDUMMY
LOGICAL CHANGED
!
! Need to exchange the entries at FIRST and SECOND and then change them
! over in the NVAL list.
!
! PRINT '(A,I8,A,I8)','swapnode> swapping minima ',FIRST,' and ', SECOND
IF (FIRST.EQ.SECOND) RETURN

NCDUMMY=NCOL(FIRST)
PDUMMY(1:NCOL(FIRST))=PBRANCH(1:NCOL(FIRST),FIRST)
NDUMMY(1:NCOL(FIRST))=NVAL(1:NCOL(FIRST),FIRST)
EDUMMY=EMKSUM(FIRST)

NCOL(FIRST)=NCOL(SECOND)
PBRANCH(1:NCOL(SECOND),FIRST)=PBRANCH(1:NCOL(SECOND),SECOND)
NVAL(1:NCOL(SECOND),FIRST)=NVAL(1:NCOL(SECOND),SECOND)
EMKSUM(FIRST)=EMKSUM(SECOND)

NCOL(SECOND)=NCDUMMY
PBRANCH(1:NCDUMMY,SECOND)=PDUMMY(1:NCDUMMY)
NVAL(1:NCDUMMY,SECOND)=NDUMMY(1:NCDUMMY)
EMKSUM(SECOND)=EDUMMY

DO J1=1,NMIN
   CHANGED=.FALSE.
   DO J2=1,NCOL(J1)
      IF (NVAL(J2,J1).EQ.FIRST) THEN
         NVAL(J2,J1)=SECOND
         CHANGED=.TRUE.
      ELSEIF (NVAL(J2,J1).EQ.SECOND) THEN
         NVAL(J2,J1)=FIRST
         CHANGED=.TRUE.
      ENDIF
   ENDDO
!
!  We rely elsewhere on ordered lists of connections, so resort if necessary.
!
   IF (CHANGED) THEN
      NDUMMY(1:NCOL(J1))=NVAL(1:NCOL(J1),J1)
      PDUMMY(1:NCOL(J1))=PBRANCH(1:NCOL(J1),J1)
      CALL SORT4(NCOL(J1),NCONNMAX,PDUMMY,NDUMMY)
      NVAL(1:NCOL(J1),J1)=NDUMMY(1:NCOL(J1))
      PBRANCH(1:NCOL(J1),J1)=PDUMMY(1:NCOL(J1))
   ENDIF
ENDDO

END SUBROUTINE SWAPNODE
