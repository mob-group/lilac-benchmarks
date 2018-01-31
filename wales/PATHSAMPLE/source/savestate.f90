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

MODULE SAVESTATE
      IMPLICIT NONE
      SAVE

      INTEGER NMINASAVE, NMINBSAVE, NMINSAVE, NTSSAVE
      INTEGER, ALLOCATABLE :: LOCATIONASAVE(:), LOCATIONBSAVE(:)

! NMIN
      DOUBLE PRECISION, ALLOCATABLE :: EMINSAVE(:), PFMINSAVE(:), FVIBMINSAVE(:), GPFOLDSAVE(:)
      DOUBLE PRECISION, ALLOCATABLE :: IXMINSAVE(:), IYMINSAVE(:), IZMINSAVE(:)
! NTS
      DOUBLE PRECISION, ALLOCATABLE :: ETSSAVE(:), KPLUSSAVE(:), KMINUSSAVE(:)

      DOUBLE PRECISION PFTOTALASAVE, PFTOTALBSAVE, PFMEANSAVE, REGROUPFREETHRESHSAVE

! NMIN
      INTEGER, ALLOCATABLE :: TOPPOINTERSAVE(:), HORDERMINSAVE(:), NCONNSAVE(:)
! NTS
      INTEGER, ALLOCATABLE :: PLUSSAVE(:), MINUSSAVE(:), POINTERMSAVE(:), POINTERPSAVE(:)

END MODULE SAVESTATE
