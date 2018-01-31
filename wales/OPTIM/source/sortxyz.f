C
C GPL License Info {{{
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
C }}}
C
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR AND ATOMIC MASS VECTOR (ATMAS
C
      SUBROUTINE SORTXYZ(XX,Y,NORD,TOL)
! Doxygen {{{
!> \name SORTXYZ
!> \brief Sort vector of nuclear coordinates 
!> \param[in] XX DP(3*NATOMS)  - input vector to be sorted
!> \param Y  DP(3*NATOMS)      - output sorted vector
!> \param NORD INT(2*NATOMS+1) - vector containing permutations 
!> \param TOL DP - tolerance factor
! }}}
! Declarations {{{ 
      USE COMMONS
      IMPLICIT NONE
! subroutine parameters  
      DOUBLE PRECISION XX(3*NATOMS),Y(3*NATOMS),TOL
      INTEGER NORD(2*NATOMS+1)
! local parameters 
!      DOUBLE PRECISION X(3*NATOMS),XX(*),Y(*),TOL
!      INTEGER I, J, JK, NORD(NATOMS), ILINE
      DOUBLE PRECISION X(3*NATOMS)
      INTEGER I, J, JK, ILINE
! }}}
! Subroutine body {{{
C
      ILINE(J)=1+J/3
C
C     Sort on the X - if two X's are equivalent, sort on Y and so on.
C     If the coordinate < the tolerance we should ignore it! However,
C     if the tolerance is sloppy that can lead to the sorting ignoring
C     genuine small differences between coordinates. Sigh. 
C
      DO I=1,3*NATOMS
         X(I)=XX(I)
      ENDDO
C
C     FIRST GIVE DUMMY ATOMS RIDICULOUS SCRATCH COORDINATES - ENSURES
C     THAT THEY WILL WIND UP AT THE BOTTOM OF THE LIST
C
      DO I=1,3*NATOMS-2,3
         IF (DABS(ATMASS(ILINE(I))).LT.1D-3) THEN
            DO J=0,2
               X(J+I) = -99995.0D0
            ENDDO
         ENDIF
      ENDDO
      JK=1
40    J=1
      DO I=1,3*NATOMS-2,3
         IF (X(I)-X(J).GT.TOL) J=I
         IF (DABS(X(I)-X(J)).LT.TOL) THEN
            IF (X(I+1)-X(J+1).GT.TOL) J=I
            IF (DABS(X(I+1)-X(J+1)).LT.TOL) THEN
               IF (X(I+2)-X(J+2).GT.TOL) J=I
            ENDIF
         ENDIF
      ENDDO
C
C     MASS-WEIGHT SORTED VECTOR - WILL ZERO ELEMENTS CORRESPONDING
C     TO DUMMY ATOMS SO THEY DONT MUCK UP THE SYMMETRY CHECKING.
C
      DO I=0,2
         Y(3*JK-2+I)=X(J+I)*NINT(ATMASS(ILINE(J)))
         X(J)=-99999.D0
      ENDDO
      NORD(JK)=(J+2)/3
      JK=JK+1
      IF (JK.EQ.NATOMS+1) GOTO 70
      GOTO 40
70    CONTINUE

! }}}
      RETURN
      END
