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
C dae 
C get coordinates from external file 'coords'
C which is CHARMM format but use standard fortran
C reading commands, unlike CHARMM which has a limit
C on total line length and does involved procedure with strings
C
      SUBROUTINE UNEWREAD(X,Y,Z,NATOMS,FILTH2,FILTHSTR)
      IMPLICIT NONE
  
      CHARACTER*4 A1
      CHARACTER(LEN=*) FILTHSTR
      CHARACTER(LEN=20) CFNAME
      INTEGER I,I1,NATOMS,FILTH2
      REAL*8 X(NATOMS),Y(NATOMS),Z(NATOMS)
 
      IF (FILTH2.EQ.0) THEN
         WRITE(CFNAME, '(A)') 'coords'
      ELSE
         WRITE(CFNAME, '(A)') 'coords.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

      OPEN(UNIT=19,FILE=CFNAME,STATUS='UNKNOWN')
      READ (19,*)
      DO I=1,NATOMS
         READ (19,'(A,3F20.10)') A1,X(I),Y(I),Z(I)
      ENDDO
  
      CLOSE(19)

      RETURN
      END


      SUBROUTINE UNRESDUMP2(COORDS,IUNIT)
      USE COMMONS
      IMPLICIT NONE

      INTEGER IUNIT,I,J,K
      REAL*8 COORDS(3*NATOMS)
      REAL*8 PEPCOORDS(3*NATOMS)

C jmc requires an open file
C also writes coordinates for dummy peptide atoms (O [representing C=O],N), for
C visualisation purposes.
      DO J=1,(NATOMS/2)-1
         DO K=1,3
            PEPCOORDS(6*(J-1)+K)=(2.0D0*COORDS(6*(J-1)+K)+COORDS(6*J+K))/3.0D0
            PEPCOORDS(6*(J-1)+K+3)=(COORDS(6*(J-1)+K)+2.0D0*COORDS(6*J+K))/3.0D0
         END DO
      END DO

      WRITE(IUNIT,'(I5)') 2*NATOMS-2
      WRITE(IUNIT,'(A)') ' ' 

      DO J=1,NATOMS/2
C Calpha
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=1,3)
C side chain
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=4,6)
C peptide atoms
         IF (J.LT.(NATOMS/2)) THEN
             WRITE(IUNIT,'(A2,3F20.10)') 'O ',(PEPCOORDS(6*(J-1)+K),K=1,3)
             WRITE(IUNIT,'(A2,3F20.10)') 'N ',(PEPCOORDS(6*(J-1)+K),K=4,6)
         END IF
      ENDDO

      RETURN

      END

      SUBROUTINE MYUNRESDUMP(COORDS,FNAMEF)
      USE COMMONS
      IMPLICIT NONE

      CHARACTER(LEN=*)  FNAMEF 
      INTEGER IUNIT,J,K
      REAL*8 COORDS(3*NATOMS)

      IUNIT=19
      OPEN(UNIT=IUNIT,FILE=FNAMEF,STATUS='UNKNOWN')
      WRITE(IUNIT,'(I5)') NATOMS

      DO J=1,NATOMS/2
C Calpha
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=1,3)
C side chain
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=4,6)
      END DO

      CLOSE(IUNIT)

      RETURN

      END

      SUBROUTINE UNRESDUMP3(COORDS,FNAMEF)
      USE COMMONS
      IMPLICIT NONE

      CHARACTER(LEN=*)  FNAMEF
      INTEGER IUNIT,I,J,K
      REAL*8 COORDS(3*NATOMS)
      REAL*8 PEPCOORDS(3*NATOMS)

      IUNIT=19
      OPEN(UNIT=IUNIT,FILE=FNAMEF,STATUS='UNKNOWN')
C also writes coordinates for dummy peptide atoms (O [representing C=O],N), for
C visualisation purposes.
      DO J=1,(NATOMS/2)-1
         DO K=1,3
            PEPCOORDS(6*(J-1)+K)=(2.0D0*COORDS(6*(J-1)+K)+COORDS(6*J+K))/3.0D0
            PEPCOORDS(6*(J-1)+K+3)=(COORDS(6*(J-1)+K)+2.0D0*COORDS(6*J+K))/3.0D0
         END DO
      END DO

      WRITE(IUNIT,'(I5)') 2*NATOMS-2
      WRITE(IUNIT,'(A)') ' '

      DO J=1,NATOMS/2
C Calpha
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=1,3)
C side chain
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=4,6)
C peptide atoms
         IF (J.LT.(NATOMS/2)) THEN
             WRITE(IUNIT,'(A2,3F20.10)') 'O ',(PEPCOORDS(6*(J-1)+K),K=1,3)
             WRITE(IUNIT,'(A2,3F20.10)') 'N ',(PEPCOORDS(6*(J-1)+K),K=4,6)
         END IF
      ENDDO

      CLOSE(IUNIT)
      RETURN

      END
