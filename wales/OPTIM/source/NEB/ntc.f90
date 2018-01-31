!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
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
MODULE NEBTOCONNECT
     IMPLICIT NONE
     SAVE
     
     INTEGER,PARAMETER :: NTSMAX = 3000 !MAXIMAL NUMBER OF TS TO BE CHECKED IN ONE DNEB RUN
     
     TYPE TSFOUNDTYPE
          DOUBLE PRECISION,POINTER :: E
          DOUBLE PRECISION,POINTER :: EVALMIN
          DOUBLE PRECISION,POINTER :: COORD(:)
          DOUBLE PRECISION,POINTER :: VECS(:)
     END TYPE TSFOUNDTYPE

     TYPE (TSFOUNDTYPE) :: TSFOUND(NTSMAX)

     INTEGER :: NTSFOUND=0

     INTEGER,PARAMETER :: NMINMAX = 3000 !MAXIMAL NUMBER OF MIN TO BE CHECKED IN ONE DNEB RUN
     
     TYPE MINFOUNDTYPE
          DOUBLE PRECISION,POINTER :: E
          DOUBLE PRECISION,POINTER :: COORD(:)
     END TYPE MINFOUNDTYPE

     TYPE (MINFOUNDTYPE) :: MINFOUND(NMINMAX)

     INTEGER :: NMINFOUND=0
END MODULE NEBTOCONNECT
