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
MODULE MODGUESS
      IMPLICIT NONE
      SAVE

      INTEGER :: GSTEPS, MAXGCYCLES ! number of steps in sequential coordinates - maximum number of sweeps
                                    ! GSTEPS could be one, then the changes are simply done one at a time
      DOUBLE PRECISION :: GTHRESHOLD, MAXINTE ! THRESHOLD CHANGE IN COORDINATES FOR APPLYING SEQUENTIAL CHANGE
                                             ! and maximum energy of an intermediate geometry
      LOGICAL :: OPTTOTAL=.FALSE., OPTMAXE=.TRUE. ! default is to optimise based on the highest energy 
                                                  ! interpolated structure.
      LOGICAL :: GUESSPATHT=.FALSE.
      INTEGER :: NINTERP  ! number of interpolated geometries
      
END MODULE MODGUESS
