C 
C  GPL License Info {{{
C
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
      SUBROUTINE ROTM(IDIR,ANG,IX,RT)
! Doxygen {{{
!>
!> \name ROTM
!> \brief Forms three dimensional rotation matrix
!>  (IX=0 if ANG is expressed in radians, IX=1 if deg, add 10 to ix to get transpose)                 
!>
!> \param[in] IDIR - axis around which the angle ANG is specified
!> \param[in] ANG - angle; can be expressed either in radians (IX=0) or in
!>              degrees (IX=1).
!> \param[in] IX - specify the input format for the angle ANG. Values:
!>              0 - radians
!>              1 - degrees
!>              Adding +10 to IX results in the matrix being transposed. 
!> \param[out] RT - the 3x3 rotation matrix 
!>
! }}}
C
C     FORMS THREE DIMENSIONAL ROTATION MATRIX (IX=0 IF ANG IN RADS,
c     IX=1 IF DEG,  ADD 10 TO IX TO GET TRANSPOSE)
C
      IMPLICIT NONE
! subroutine parameters 
      INTEGER IDIR,IX
      DOUBLE PRECISION RT(3,3),ANG

! local parameters 
      INTEGER I,J,INK,IEO,IBTAND
      DOUBLE PRECISION ATOR,TANG
! Subroutine body {{{
      IBTAND(I,J) = IAND(I,J)
      IEO(I)=2*IBTAND(I,1)-1
C
      ATOR=DACOS(-1.D0)/180.D0
      INK=IX/10
      RT(1:3,1:3)=0.0D0
      TANG=ANG
      IF(IX.EQ.1)TANG=ANG*ATOR
      DO J=1,3
         RT(J,J)=DCOS(TANG)
         RT(IDIR,J)=IBTAND(IDIR,J)/MAX(IDIR,J)
      ENDDO
      DO I=1,3
         DO J=1,I
            IF(I.NE.J.AND.I.NE.IDIR.AND.J.NE.IDIR)THEN
               RT(I,J)=IEO(I*J)*DSIN(TANG)
               RT(J,I)=-RT(I,J)
            ENDIF
         ENDDO
      ENDDO
      IF(INK.NE.0)CALL MTRANSOPT(RT,RT,3,3,3,3)
      ! }}}
      RETURN
      END
