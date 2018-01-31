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
SUBROUTINE NGTREALLOC_CRSTORAGE(MIN1,NEWNCOL,MIN2,GBMAX,NMIN,NCOL)
USE NGTMEM
USE COMMONS, ONLY : NCONN, NCONNMIN, DEBUG
IMPLICIT NONE
INTEGER MIN1, NEWNCOL, MIN2, NMIN, NONZERO, J3, OLDROWPTR1, OLDROWPTR2, OFFSET, SIZEVEC
INTEGER NCOL(NMIN)
DOUBLE PRECISION GBMAX

! JMC work out the new size of the dvec and col_ind arrays
   OFFSET=NEWNCOL - ROW_PTR(MIN2+1) + ROW_PTR(MIN2)
   NONZERO=ROW_PTR(MIN1) + NCOL(MIN1) -1 + OFFSET
   
! JMC the line below can improve efficiency if we have allowed more space than strictly necessary for each minimum in 
! dvec and col_ind, as we don't have to move around array sections.  On the other hand, it may mean 
! reallocation happens more often as we don't "bunch up" 
! the elements in the array when they have too much space.
   IF(OFFSET.LE.0)RETURN

   OLDROWPTR2=ROW_PTR(MIN2+1)
   OLDROWPTR1=ROW_PTR(MIN1)+NCOL(MIN1)-1

! need to reset row_ptr here
   DO J3=MIN2+1,NMIN
      ROW_PTR(J3)=ROW_PTR(J3) + OFFSET
   END DO

   SIZEVEC = SIZE(DVEC)

! Do we need to reallocate DVEC and COL_IND?
   IF (NONZERO.GT.SIZEVEC) THEN
!      PRINT '(A,F15.5,2(A,I8))','NGTrealloc_crstorage> allocating temporary array size ',8.0D0*SIZEVEC*1.0D-6, &
!     &                          ' Mb; remaining nodes=',MIN1
      ALLOCATE(DVECTEMP(SIZEVEC),COL_INDTEMP(SIZEVEC))
      DO J3=1,SIZEVEC ! JMC not using whole array operations for efficiency here, as this routine is a considerable bottleneck!!!
         DVECTEMP(J3)=DVEC(J3)
         COL_INDTEMP(J3)=COL_IND(J3)
      END DO
      DEALLOCATE(DVEC,COL_IND)
      IF (8.0D0*NONZERO*1.0D-9.GT.GBMAX) GBMAX=8.0D0*NONZERO*1.0D-9
      IF (DEBUG) PRINT '(A,F15.5,2(A,I8))','NGTrealloc_crstorage> allocating new array size       ',8.0D0*NONZERO*1.0D-6, &
        &                       ' Mb; remaining nodes=',MIN1
      ALLOCATE(DVEC(NONZERO),COL_IND(NONZERO))
      DO J3=1,ROW_PTR(MIN2)-1
         DVEC(J3)=DVECTEMP(J3)
         COL_IND(J3)=COL_INDTEMP(J3)
      END DO
   ELSE
!      PRINT '(A,F15.5,2(A,I8))','NGTrealloc_crstorage> allocating temporary array size ',8.0D0*SIZEVEC*1.0D-6, &
!     &                          ' Mb; remaining nodes=',MIN1
      ALLOCATE(DVECTEMP(SIZEVEC),COL_INDTEMP(SIZEVEC))
! JMC not using whole array operations for efficiency here, as this routine is a considerable bottleneck!!!
      DO J3=OLDROWPTR2,OLDROWPTR1 
         DVECTEMP(J3)=DVEC(J3)
         COL_INDTEMP(J3)=COL_IND(J3)
      END DO
! JMC commented out the zero-ing below for efficiency, and trusting to the fact that we should not need to access 
! elements > ROW_PTR(MIN1)-1+NCOL(MIN1) anyway
      !DO J3=ROW_PTR(MIN1)+NCOL(MIN1),SIZEVEC
      !   DVEC(J3)=0.0D0
      !   COL_IND(J3)=0
      !END DO
   END IF

   DVEC(ROW_PTR(MIN2+1):ROW_PTR(MIN1)+NCOL(MIN1)-1)=DVECTEMP(OLDROWPTR2:OLDROWPTR1)
   COL_IND(ROW_PTR(MIN2+1):ROW_PTR(MIN1)+NCOL(MIN1)-1)=COL_INDTEMP(OLDROWPTR2:OLDROWPTR1)

   DEALLOCATE(DVECTEMP,COL_INDTEMP)

END SUBROUTINE NGTREALLOC_CRSTORAGE
