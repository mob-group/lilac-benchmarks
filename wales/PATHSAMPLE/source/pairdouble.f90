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

SUBROUTINE PAIRDOUBLE
   USE COMMONS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE :: VINT(:)

   ALLOCATE(VINT(MAXPAIRS))

   VINT(1:MAXPAIRS)=PAIR1(1:MAXPAIRS)
   DEALLOCATE(PAIR1)
   ALLOCATE(PAIR1(2*MAXPAIRS))
   PAIR1(1:MAXPAIRS)=VINT(1:MAXPAIRS)

   VINT(1:MAXPAIRS)=PAIR2(1:MAXPAIRS)
   DEALLOCATE(PAIR2)
   ALLOCATE(PAIR2(2*MAXPAIRS))
   PAIR2(1:MAXPAIRS)=VINT(1:MAXPAIRS)

   MAXPAIRS=2*MAXPAIRS
   PRINT '(A,I8)','pairsdouble> Maximum number of pairs tried increased to ',MAXPAIRS
   DEALLOCATE(VINT)

END SUBROUTINE PAIRDOUBLE
