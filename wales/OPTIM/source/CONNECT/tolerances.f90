!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
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
MODULE TOLERANCES
     IMPLICIT NONE
     SAVE
     DOUBLE PRECISION :: GEOMDIFFTOL = 1.0D-3 ! SPECIFIES THE DISTANCE CRITERION USED TO IDENTIFY IDENTICAL ISOMERS IN CONNECT RUNS
     DOUBLE PRECISION :: EDIFFTOL    = 1.0D-9 ! SPECIFIES THE ENERGY DIFFERENCE USED TO IDENTIFY PERMUTATIONAL ISOMERS IN CONNECT RUNS
END MODULE TOLERANCES
