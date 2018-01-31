C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE TWODFUNC(COORDS,GRAD,ENERGY,GTEST,STEST)
      USE COMMONS
      USE MODHESS
      IMPLICIT NONE
      DOUBLE PRECISION COORDS(3*NATOMS),GRAD(3*NATOMS),X,Y,ENERGY
      LOGICAL GTEST, STEST

      PRINT '(A)','2Dfunc is currently a dummy routine'
      STOP

      RETURN
      END
