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
      SUBROUTINE CLATMIN(COORDS,V)
      USE COMMONS
      IMPLICIT NONE
      INTEGER NELEMENTS, NTYPE(105), NCOUNT
      DOUBLE PRECISION F1,F2,F3,GRAD,SECOND, V(3*NATOMS),
     1                 AMAT(3,3), AINV(3,3), COORDS(3*NATOMS), DIF, TEMP1, RMS
      COMMON /CAS/ AMAT, AINV, NELEMENTS, NTYPE
C
C  Value of DIF is the order of magnitude to which the lattice constant can be optimised. 
C
      DIF=1.0D-3
      NCOUNT=1

10    TEMP1=AMAT(1,1)
      CALL POTENTIAL(COORDS,F3,V,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)

      AMAT(1,1)=AMAT(1,1)+DIF
      AMAT(2,2)=AMAT(2,2)+DIF
      AMAT(3,3)=AMAT(3,3)+DIF

      CALL POTENTIAL(COORDS,F1,V,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)

      AMAT(1,1)=AMAT(1,1)-2*DIF
      AMAT(2,2)=AMAT(2,2)-2*DIF
      AMAT(3,3)=AMAT(3,3)-2*DIF
      CALL POTENTIAL(COORDS,F2,V,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)

      GRAD=(F1-F2)/(2.0D0*DIF)

      SECOND=(F1+F2-2.0D0*F3)/(DIF*DIF)
      PRINT*,'Energy for lattice cycle ',NCOUNT,' is ',F3,' length=',TEMP1
      PRINT*,'Gradient wrt box length=',GRAD
      PRINT*,'Second derivative wrt box length=',SECOND
      WRITE(*,'(A,3F20.10)') 'F1,F2,F3=',F1,F2,F3
      NCOUNT=NCOUNT+1
      IF (ABS(GRAD/SECOND).GT.5.0D-2) THEN
         AMAT(1,1)=TEMP1-5.0D-2*GRAD/DABS(GRAD)
         AMAT(2,2)=AMAT(1,1)
         AMAT(3,3)=AMAT(1,1)
         PRINT*,'Step=',-5.0D-2*GRAD/DABS(GRAD)
         GOTO 10
      ELSE
         AMAT(1,1)=TEMP1-GRAD/ABS(SECOND)
         AMAT(2,2)=AMAT(1,1)
         AMAT(3,3)=AMAT(1,1)
         PRINT*,'Step=',-GRAD/ABS(SECOND)
         IF (DABS(GRAD/SECOND).GT.DIF/10) GOTO 10
      ENDIF
      RETURN
      END
