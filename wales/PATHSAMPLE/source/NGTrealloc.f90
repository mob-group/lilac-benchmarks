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
!  Reallocate NGT memory if required. We don't need PBRANCH(x,y) and NVAL(x,y) beyond y=MIN1 
!  so long as regrouping has been performed to put the A minima first, followed by the B minima,
!  followed by the I minima.
!  This reordering is done in the REGROUP subroutine.
!
!
SUBROUTINE NGTREALLOC(NCONNMAX,MIN1,NEWNCOL,MIN2,GBMAX,NMIN,NCOL)
USE NGTMEM
USE COMMONS, ONLY : SHANNONT
IMPLICIT NONE
INTEGER NCONNMAXO, NCONNMAX, MIN1, NEWNCOL, MIN2, NMIN
INTEGER NCOL(NMIN)
DOUBLE PRECISION GBMAX

   IF (.NOT.SHANNONT) PRINT '(A,F15.5,2(A,I8))','NGTrealloc> allocating temporary array size ',8.0D0*NCONNMAX*MIN1*1.0D-6, &
     &                       ' Mb; maximum connections=',NCONNMAX,' remaining nodes=',MIN1
   ALLOCATE(PBRANCHTEMP(NCONNMAX,MIN1))
   PBRANCHTEMP(1:NCONNMAX,1:MIN1)=PBRANCH(1:NCONNMAX,1:MIN1)
   DEALLOCATE(PBRANCH)
   ALLOCATE(NVALTEMP(NCONNMAX,MIN1))
   NVALTEMP(1:NCONNMAX,1:MIN1)=NVAL(1:NCONNMAX,1:MIN1)
   DEALLOCATE(NVAL)
   NCONNMAXO=NCONNMAX
   NCONNMAX=MAX(NCONNMAX,NEWNCOL+NCOL(MIN2)+20) ! the +20 is to reduce the number of calls to NGTREALLOC
   IF (8.0D0*NCONNMAX*MIN1*1.0D-9.GT.GBMAX) GBMAX=8.0D0*NCONNMAX*MIN1*1.0D-9
   IF (.NOT.SHANNONT) THEN
      PRINT '(A,F15.5,2(A,I8))','NGTrealloc> allocating new array size       ',8.0D0*NCONNMAX*MIN1*1.0D-6, &
     &                       ' Mb; maximum connections=',NCONNMAX,' remaining nodes=',MIN1
   ENDIF
   ALLOCATE(PBRANCH(NCONNMAX,MIN1))
   PBRANCH=0.0D0
   PBRANCH(1:NCONNMAXO,1:MIN1)=PBRANCHTEMP(1:NCONNMAXO,1:MIN1)
   DEALLOCATE(PBRANCHTEMP)
   ALLOCATE(NVAL(NCONNMAX,MIN1))
   NVAL=0
   NVAL(1:NCONNMAXO,1:MIN1)=NVALTEMP(1:NCONNMAXO,1:MIN1)
   DEALLOCATE(NVALTEMP)

END SUBROUTINE NGTREALLOC
