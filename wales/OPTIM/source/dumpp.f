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
C
C***********************************************************************
C
      SUBROUTINE DUMPP(QL,ENERGY)
      USE COMMONS
      USE KEY
      USE modcharmm
      use porfuncs
      IMPLICIT NONE
      DOUBLE PRECISION QL(3*NATOMS), S, ENERGY, OVEC(3), H1VEC(3), H2VEC(3), XS(3*NTSITES)
      INTEGER J1, J2, J3, ISTAT

      IF (.NOT.PRINTPTS) RETURN
C
C  Write out the Cartesian coordinates to file points and the energy to
C  file energies
C
      WRITE(2,'(G30.20)') ENERGY
      S=1.0D0
      IF ((ZSYM(NATOMS).EQ.'AR').AND.(ZSYM(1).NE.'CA')) S=3.4D0
C     IF (ZSYM(1)(1:1).EQ.'W') THEN
C        WRITE(1,'(I5)') 3*NATOMS
C        WRITE(1,'(F20.10)') ENERGY
C     ENDIF
      DO J1=1,NATOMS
         IF (VARIABLES) THEN
            WRITE(1,'(G20.10)') QL(J1)
         ELSE IF (ZSYM(1)(1:1).EQ.'W') THEN
            IF (J1.LE.NATOMS/2) THEN  !  WCOMMENT new line
            CALL CONVERT(QL(3*(J1-1)+1),QL(3*(J1-1)+2),QL(3*(J1-1)+3),
     1            QL(3*(NATOMS/2+J1-1)+1),QL(3*(NATOMS/2+J1-1)+2),QL(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
C  ! WCOMMENT
C    1            QL(3*(NATOMS+J1-1)+1),QL(3*(NATOMS+J1-1)+2),QL(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
C           WRITE(1,'(A3,3G20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
C           WRITE(1,'(A3,3G20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(1,'(A3,3G20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
            WRITE(1,'(2X,3G20.10)') OVEC(1),OVEC(2),OVEC(3)
            WRITE(1,'(2X,3G20.10)') H1VEC(1),H1VEC(2),H1VEC(3)
            WRITE(1,'(2X,3G20.10)') H2VEC(1),H2VEC(2),H2VEC(3)
            ENDIF
         ELSE IF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
            WRITE(1,'(G20.10)') QL(J1)
         ELSE IF ((IATNUM(J1).NE.0).OR.CHRMMT.OR.UNRST) THEN
            J2=3*(J1-1)
            WRITE(1,'(1X,3G20.10)') QL(J2+1)*S,QL(J2+2)*S,QL(J2+3)*S
         ELSE
            J2 = 3*(J1-1)
            WRITE(1,'(1X,3G20.10)') QL(J2+1), QL(J2+2), QL(J2+3)
         ENDIF
      ENDDO
      CALL FLUSH(1)
      CALL FLUSH(2)
      RETURN
      END
