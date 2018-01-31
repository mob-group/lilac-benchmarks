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
MODULE COMMONS
      IMPLICIT NONE

      INTEGER :: NATOMS, ISTATUS, NUNIQUE, NOPT, IPRNT, INR, IVEC, IVEC2, ISTCRT, NMOL, NINTS, NRBSITES
      DOUBLE PRECISION :: ZSTAR, RHO, EVDISTTHRESH, NEWRES_TEMP, PERCCUT
      DOUBLE PRECISION :: PARAM1=0.0D0, PARAM2=0.0D0, PARAM3=0.0D0, PARAM4=0.0D0, PARAM5=0.0D0, PARAM6=0.0D0, PARAM7=0.0D0
      LOGICAL REDOPATH, REDOPATHXYZ, RATIOS, REDOPATHNEB, MPIT, DEBUG

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NR, IATNUM      !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ATMASS !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: STPMAX ! 3*MXATMS
      CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:) :: ZSYM   !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE :: TAGFAC(:)
      DOUBLE PRECISION, ALLOCATABLE :: RBSITE(:,:), RBSTLA(:,:), STCHRG(:), DPMU(:)
      INTEGER, ALLOCATABLE :: TAGNUM(:)
      INTEGER NTAG
      LOGICAL TAGT, CHANGE_TEMP, PERCT
      DOUBLE PRECISION, PARAMETER :: numberPI=3.141592653589793115997963468544185162d0 !DACOS(-1.D0)
      integer, parameter :: rk = kind(1.d0) ! read kind

END MODULE COMMONS
