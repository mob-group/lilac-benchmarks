!   GMIN: A program for finding global minima
!   Copyright (C) 1999- David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE PAIRCOLOUR(NATOMS,P,ZSTRING)
USE COMMONS, ONLY : ZSYM
IMPLICIT NONE

INTEGER J1, J2, J3, NATOMS, J4
DOUBLE PRECISION EREAL, P(3*NATOMS)
DOUBLE PRECISION VTINT, VTMIN, VTMAX, VT(NATOMS), R6, DUMMY, DIST
CHARACTER(LEN=2) ZSTRING(NATOMS)

ZSTRING(1:NATOMS)='X1'
!
! Calculate LJ pair energies.
!
VT(1:NATOMS)=0.0D0
IF (ZSYM(1).EQ.'AX') THEN
   DO J1=1,NATOMS
      J3=3*J1
      DO J2=J1+1,NATOMS
         J4=3*J2
         DIST=(P(J3-2)-P(J4-2))**2+(P(J3-1)-P(J4-1))**2+(P(J3)-P(J4))**2
         DIST=1.0D0/DIST
         R6=DIST**3
         DUMMY=R6*(R6-1.0D0)
         VT(J1)=VT(J1)+DUMMY
         VT(J2)=VT(J2)+DUMMY
      ENDDO
   ENDDO
ELSE
   PRINT '(A)',' paircolour> WARNING *** pair energy colouring not yet coded for potential ' // ZSYM
   RETURN
ENDIF
!
! Sort the atoms according to pair energy.
!
CALL SORT3(NATOMS,NATOMS,VT,P)

VTMIN=1.0D100
VTMAX=-1.0D100
DO J1=1,NATOMS
   IF (VT(J1).GT.VTMAX) VTMAX=VT(J1)
   IF (VT(J1).LT.VTMIN) VTMIN=VT(J1)
ENDDO
VTINT=VTMAX-VTMIN
IF (VTINT.GT.0.0D0) THEN
   DO J1=1,NATOMS
      IF (VT(J1).GT.VTMIN+2.0D0*VTINT/3.0D0) THEN
         ZSTRING(J1)='X3'
      ELSE IF (VT(J1).GT.VTMIN+VTINT/3.0D0) THEN
         ZSTRING(J1)='X2'
      ENDIF
   ENDDO
ENDIF

RETURN
END SUBROUTINE PAIRCOLOUR
!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!     This subprogram performs a sort on the input data and
!     arranges it from smallest to biggest. The exchange-sort
!     algorithm is used.
!
      SUBROUTINE SORT3(N,J3,A,B)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2
      DOUBLE PRECISION TEMP, A(J3), B(3*J3)
!
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).LT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         DO J2=0,2
            TEMP=B(3*L-J2)
            B(3*L-J2)=B(3*J1-J2)
            B(3*J1-J2)=TEMP
         ENDDO
20    CONTINUE
      RETURN
      END
