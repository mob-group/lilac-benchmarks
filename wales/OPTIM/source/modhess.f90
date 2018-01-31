!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE MODHESS
      IMPLICIT NONE
      SAVE

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS  !  3*MXATMS,3*MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: NONFROZENHESS
      INTEGER         , ALLOCATABLE, DIMENSION(:)   :: NONFREEZELIST
      ! hk286 - toggle between moving and stationary frame
      LOGICAL :: RBAANORMALMODET

CONTAINS

SUBROUTINE HESSOUT(NSYM)
   IMPLICIT NONE
! Arguments
   INTEGER, INTENT(IN)  :: NSYM

   IF (ALLOCATED(HESS) .AND. (NSYM == 0 .OR. NSYM == 1)) THEN
      WRITE(*, '(12X,5F12.6)') HESS(:, :)
   END IF

END SUBROUTINE HESSOUT

END MODULE MODHESS
