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
MODULE msevb_interfaces
  
  INTERFACE
     SUBROUTINE msevb(ncoords, systemCoords, assignVBstates, energy, onlyAssignStates)
       USE msevb_common
       INTEGER, INTENT(IN) :: ncoords  
       DOUBLE PRECISION, INTENT(IN) :: SYSTEMCOORDS(NCOORDS)
       LOGICAL, INTENT(IN) :: assignVBstates
       DOUBLE PRECISION, INTENT(OUT) :: ENERGY
       LOGICAL, OPTIONAL, INTENT(IN) :: onlyAssignStates
     END SUBROUTINE msevb

     SUBROUTINE assignZundelSpecies(firstVBstate, secondVBstate, exchangeH, checkedHbond, minHbondAngleRad)
       USE msevb_common     
       INTEGER, INTENT(IN) :: firstVBstate, secondVBstate, exchangeH
       LOGICAL, INTENT(IN) :: checkedHbond
       DOUBLE PRECISION, INTENT(IN) :: MINHBONDANGLERAD

     END SUBROUTINE assignZundelSpecies

  END INTERFACE

END MODULE msevb_interfaces
